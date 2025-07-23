library(shiny)
library(bio3d)
library(ggplot2)
library(plotly)
library(shinybusy)
library(r3dmol)

ui <- fluidPage(
  add_busy_spinner(spin = "fading-circle"),
  
  titlePanel("Structural k-means clustering (.cif / .pdb files)"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload .cif or .pdb files", multiple = TRUE,
                accept = c(".cif", ".pdb")),
      sliderInput("k", "Number of clusters (k)", min = 2, max = 10, value = 3),
      #checkboxInput("showLabels", "Show labels in MDS plot", FALSE),
      checkboxInput("use3D", "3D MDS Plot", FALSE),
      actionButton("run", "Run Clustering"),
      downloadButton("downloadData", "Download Cluster Representatives (.zip)"),
      downloadButton("downloadTable", "Download Cluster Summary (.csv)"),
      uiOutput("selectClusterUI")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("MDS Plot", plotlyOutput("mdsPlot", height = "600px")),
        tabPanel("Cluster Summary", tableOutput("clusterTable")),
        tabPanel("Cluster Sizes", plotOutput("barPlot")),
        tabPanel("3D Structure Viewer", r3dmolOutput("structureViewer", height = "600px")),
        tabPanel("Cif/Pdb Converter",
                 fileInput("convert_files", "Upload .cif or .pdb files", multiple = TRUE,
                           accept = c(".cif", ".pdb")),
                 radioButtons("convert_direction", "Conversion Direction:",
                              choices = c("cif → pdb" = "cif2pdb", "pdb → cif" = "pdb2cif")),
                 actionButton("convert_button", "Convert"),
                 uiOutput("converted_files_ui"),
                 downloadButton("downloadConverted", "Download Converted Files (.zip)")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  results <- reactiveValues()
  
  observeEvent(input$run, {
    req(input$files)
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
    if (input$use3D) {
      plot_ly(df, x = ~X, y = ~Y, z = ~Z, color = ~Cluster, colors = "Set1",
              type = 'scatter3d', #mode = if (input$showLabels) 'markers+text' else 'markers',
              text = ~File)
    } else {
      plot_ly(df, x = ~X, y = ~Y, color = ~Cluster, colors = "Set1",
              type = 'scatter', #mode = if (input$showLabels) 'markers+text' else 'markers',
              text = ~File)
    }
  })
  
  output$barPlot <- renderPlot({
    req(results$km)
    barplot(
      table(results$km$cluster),
      main = "Number of Structures per Cluster",
      xlab = "Cluster",
      ylab = "Count",
      col = rainbow(input$k)
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
