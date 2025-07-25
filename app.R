library(shiny)
library(bio3d)
library(ggplot2)
library(plotly)
library(shinybusy)
library(r3dmol)
library(bslib)
library(shinydashboard)
library(waiter)

source("theme.R")  # Cargá tu tema

ui <- fluidPage(
  use_waiter(),
  theme = my_theme,
  add_busy_spinner(spin = "fading-circle"),
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  
  tags$div(
    style = "padding: 10px 20px;",
    tags$h1("Rkmeans", style = "margin-bottom: 5px; font-weight: bold; color: black;"),
    tags$h4("Structural k-means clustering (.cif / .pdb files)", style = "margin-top: 0; color: #555555;")
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload .cif or .pdb files", multiple = TRUE,
                accept = c(".cif", ".pdb")),
      sliderInput("k", "Number of clusters (k)", min = 2, max = 10, value = 3),
      checkboxInput("use3D", "3D MDS Plot", FALSE),
      div(
        style = "margin-bottom: 15px;",
        actionButton("run", "Run Clustering", class = "btn-primary")
      ),
      
      div(
        style = "margin-bottom: 15px;",
        downloadButton("downloadData", "Download Cluster Representatives (.zip)")
      ),
      
      div(
        style = "margin-bottom: 15px;",
        downloadButton("downloadTable", "Download Cluster Summary (.csv)")
      ),
      div(
        style = "margin-bottom: 15px;",
        uiOutput("selectClusterUI")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("MDS Plot", plotlyOutput("mdsPlot", height = "600px")),
        tabPanel("Cluster Summary", tableOutput("clusterTable")),
        tabPanel("Cluster Sizes", plotlyOutput("barPlot")),
        tabPanel("3D Structure Viewer", r3dmolOutput("structureViewer", height = "600px")),
        tabPanel("Cif/Pdb Converter",
                 tags$div(style = "padding: 20px; max-width: 600px;",
                          tags$h6("Upload Files for Conversion"),
                          fileInput("convert_files", NULL, 
                                    buttonLabel = "Browse...", 
                                    placeholder = "No file selected",
                                    multiple = TRUE,
                                    accept = c(".cif", ".pdb")),
                          
                          tags$div(style = "margin-top: 20px; margin-bottom: 10px;",
                                   tags$label("Conversion Direction:", style = "font-weight: 500; color: #2c3e50;"),
                                   radioButtons("convert_direction", NULL,
                                                choices = c("CIF → PDB" = "cif2pdb", "PDB → CIF" = "pdb2cif"),
                                                inline = TRUE
                                   )
                          ),
                          
                          div(
                            style = "margin-top: 10px; margin-bottom: 20px;",
                            actionButton("convert_button", "Convert Files", class = "btn btn-primary")
                          ),
                          
                          uiOutput("converted_files_ui"),
                          
                          div(
                            style = "margin-top: 20px;",
                            downloadButton("downloadConverted", "Download Converted Files (.zip)", class = "btn btn-primary")
                          )
                 )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  results <- reactiveValues()
  
  observeEvent(input$run, {
    req(input$files)
    
    # waiter
    waiter_show(html = tagList(
      spin_fading_circles(), 
      tags$br(),
      tags$h4("Processing and clustering...", style = "color: #333;")
    ), color = "#ffffffcc")
    
    on.exit(waiter_hide(), add = TRUE)
    #
    
    files <- input$files
    paths <- setNames(files$datapath, files$name)
    
    structs <- lapply(paths, function(p) {
      ext <- tools::file_ext(p)
      if (tolower(ext) == "cif") read.cif(p) else read.pdb(p)
    })
    
    coords <- lapply(structs, function(p) p$xyz)
    lens <- lengths(coords)
    if (length(unique(lens)) != 1) {
      showNotification("❌ All structures must have the same number of atoms.", type = "error")
      return()
    }
    
    withProgress(message = "Clustering in progress...", value = 0, {
      n <- length(coords)
      rmsd_mat <- matrix(0, n, n)
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          rms <- rmsd(coords[[i]], coords[[j]])
          rmsd_mat[i, j] <- rms
          rmsd_mat[j, i] <- rms
        }
      }
      incProgress(0.3)
      
      if (any(is.na(rmsd_mat))) {
        showNotification("❌ RMSD matrix contains NA values.", type = "error")
        return()
      }
      
      mds <- cmdscale(rmsd_mat, k = 3)
      k <- input$k
      km <- kmeans(mds, centers = k, nstart = 20)
      incProgress(0.3)
      
      tmpdir <- tempdir()
      rep_files <- character(k)
      for (cluster_id in 1:k) {
        inds <- which(km$cluster == cluster_id)
        submat <- rmsd_mat[inds, inds]
        medoid_idx <- inds[which.min(colMeans(submat))]
        rep_name <- paste0("cluster_", cluster_id, "_rep.cif")
        rep_path <- file.path(tmpdir, rep_name)
        file.copy(paths[[medoid_idx]], rep_path, overwrite = TRUE)
        rep_files[cluster_id] <- rep_path
      }
      incProgress(0.2)
      
      results$mds_df <- data.frame(X = mds[,1], Y = mds[,2], Z = mds[,3], Cluster = as.factor(km$cluster), File = names(paths))
      results$km <- km
      results$rmsd_mat <- rmsd_mat
      results$paths <- paths
      results$rep_files <- setNames(rep_files, names(rep_files))
      
      rep_files_pdb <- convert_structure_with_biopython(
        setNames(rep_files, basename(rep_files)),
        output_dir = file.path(tmpdir, "pdbs")
      )
      
      results$rep_files_pdb <- rep_files_pdb
      
      stats <- lapply(1:k, function(cluster_id) {
        inds <- which(km$cluster == cluster_id)
        submat <- rmsd_mat[inds, inds]
        upper <- submat[upper.tri(submat)]
        mds_cluster <- mds[inds, , drop = FALSE]
        centroid <- colMeans(mds_cluster)
        dists_to_centroid <- sqrt(rowSums((mds_cluster - matrix(centroid, nrow = nrow(mds_cluster), ncol = 3, byrow = TRUE))^2))
        
        data.frame(
          Cluster = cluster_id,
          Size = length(inds),
          Avg_RMSD = mean(upper),
          Max_RMSD = max(upper),
          SD_RMSD = sd(upper),
          Dist_to_Centroid = mean(dists_to_centroid)
        )
      })
      results$stats <- do.call(rbind, stats)
    })
  })
  
  output$mdsPlot <- renderPlotly({
    req(results$mds_df)
    df <- results$mds_df
    
    marker_style <- list(
      size = 10,                 
      opacity = 0.7
    )
    
    if (input$use3D) {
      plot_ly(
        df,
        x = ~X, y = ~Y, z = ~Z,
        color = ~Cluster,
        type = 'scatter3d',
        mode = 'markers',
        text = ~File,
        marker = marker_style
      ) %>%
        layout(
          legend = list(title = list(text = "Cluster"))
        )
      
    } else {
      plot_ly(
        df,
        x = ~X, y = ~Y,
        color = ~Cluster,
        type = 'scatter',
        mode = 'markers',
        text = ~File,
        marker = marker_style
      ) %>%
        layout(
          legend = list(title = list(text = "Cluster"))
        )
    }
  })
  
  
  output$barPlot <- renderPlotly({
    req(results$km)
    
    df <- data.frame(Cluster = factor(results$km$cluster)) %>%
      count(Cluster)
    
    cluster_colors <- c(
      "#A8D5BA", "#79B7A4", "#5D9B9B", "#4F7E87",
      "#3B6072", "#2F4858", "#A4CBB1", "#8FC1A9", "#C1E1C1", "#6BA292"
    )
    
    plot_ly(
      data = df,
      x = ~Cluster,
      y = ~n,
      type = 'bar',
      textposition = 'auto',
      marker = list(
        color = cluster_colors[seq_along(df$Cluster)]
      )
    ) %>%
      layout(
        title = list(
          text = "Number of Structures per Cluster",
          font = list(size = 20, color = "#2c3e50")
        ),
        xaxis = list(
          title = "Cluster",
          tickfont = list(size = 14, color = "#2c3e50"),
          titlefont = list(size = 16, color = "#2c3e50")
        ),
        yaxis = list(
          title = "Count",
          tickfont = list(size = 14, color = "#2c3e50"),
          titlefont = list(size = 16, color = "#2c3e50")
        ),
        plot_bgcolor = "#ffffff",
        paper_bgcolor = "#ffffff",
        margin = list(l = 60, r = 30, b = 60, t = 70)
      )
  })
  
  
  
  output$clusterTable <- renderTable({
    req(results$stats)
    results$stats
  })
  
  output$downloadData <- downloadHandler(
    filename = function() "cluster_representatives.zip",
    content = function(file) {
      zip(file, files = results$rep_files)
    }
  )
  
  output$downloadTable <- downloadHandler(
    filename = function() "cluster_summary.csv",
    content = function(file) {
      write.csv(results$stats, file, row.names = FALSE)
    }
  )
  
  output$selectClusterUI <- renderUI({
    req(results$rep_files)
    selectInput("selectedCluster", "Select Cluster to View:", 
                choices = paste0("Cluster ", seq_along(results$rep_files)), 
                selected = "Cluster 1")
  })
  
  check_python <- function() {
    result <- suppressWarnings(system("python --version", intern = TRUE, ignore.stderr = TRUE))
    if (length(result) == 0) {
      stop("Python not found. Please install Python and ensure it's added to your system PATH.")
    }
  }
  
  convert_structure_with_biopython <- function(input_files, direction = "cif2pdb", output_dir = "converted", 
                                               python_exec = NULL) {
    if (is.null(python_exec)) {
      python_exec <- Sys.which("python")
      if (python_exec == "") {
        stop("❌ Python not found. Please ensure it's installed and in your PATH.")
      }
    }
    
    script_path <- switch(direction,
                          cif2pdb = normalizePath("cif_to_pdb.py", mustWork = TRUE),
                          pdb2cif = normalizePath("pdb_to_cif.py", mustWork = TRUE),
                          stop("Invalid direction."))
    
    output_dir <- normalizePath(output_dir, mustWork = FALSE)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    out_files <- character(length(input_files))
    
    for (i in seq_along(input_files)) {
      input <- normalizePath(input_files[[i]], mustWork = TRUE)
      original_name <- names(input_files)[i]
      base <- tools::file_path_sans_ext(basename(original_name))
      output_ext <- ifelse(direction == "cif2pdb", ".pdb", ".cif")
      output_path <- file.path(output_dir, paste0(base, output_ext))
      
      cmd <- sprintf('"%s" "%s" "%s" "%s"', python_exec, script_path, input, output_path)
      message("[CMD] ", cmd)
      system(cmd)
      
      out_files[i] <- output_path
    }
    
    return(out_files)
  }
  
  output$structureViewer <- renderR3dmol({
    req(input$selectedCluster, results$rep_files_pdb)
    
    cluster_num <- as.integer(gsub("Cluster ", "", input$selectedCluster))
    file_path <- results$rep_files_pdb[[cluster_num]]
    
    r3dmol() %>%
      m_add_model(data = readChar(file_path, file.info(file_path)$size), format = "pdb") %>%
      m_set_style(style = m_style_cartoon()) %>%
      m_zoom_to()
  })
  
  converted_files <- reactiveVal(NULL)
  
  observeEvent(input$convert_button, {
    req(input$convert_files)
    files <- input$convert_files
    paths <- files$datapath
    names(paths) <- files$name
    
    dir_out <- file.path(tempdir(), "converted_files")
    if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
    
    converted_paths <- switch(
      input$convert_direction,
      cif2pdb = convert_structure_with_biopython(paths, direction = "cif2pdb", output_dir = dir_out),
      pdb2cif = convert_structure_with_biopython(paths, direction = "pdb2cif", output_dir = dir_out)
    )
    converted_files(converted_paths)
  })
  
  output$converted_files_ui <- renderUI({
    req(converted_files())
    tags$ul(lapply(converted_files(), function(path) {
      tags$li(basename(path))
    }))
  })
  
  output$downloadConverted <- downloadHandler(
    filename = function() {
      paste0("converted_files_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(converted_files())
      zip::zip(zipfile = file, files = converted_files(), mode = "cherry-pick")
    }
  )
  
} # close server

shinyApp(ui, server)
