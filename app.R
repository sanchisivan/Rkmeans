
library(shiny)
library(bio3d)
library(ggplot2)
library(plotly)

ui <- fluidPage(
  titlePanel("Structural Clustering (.cif / .pdb files)"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload .cif or .pdb files", multiple = TRUE,
                accept = c(".cif", ".pdb")),
      numericInput("k", "Number of clusters (k)", value = 3, min = 2),
      actionButton("run", "Run Clustering"),
      downloadButton("downloadData", "Download Cluster Representatives (.cif)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("MDS Plot", plotlyOutput("mdsPlot")),
        tabPanel("Cluster Summary", tableOutput("clusterTable")),
        tabPanel("Cluster Sizes", plotOutput("barPlot"))
      )
    )
  )
)

server <- function(input, output, session) {
  results <- reactiveValues()
  
  observeEvent(input$run, {
    req(input$files)
    files <- input$files
    paths <- files$datapath
    names(paths) <- files$name
    
    # Read structures
    structs <- lapply(paths, function(p) {
      ext <- tools::file_ext(p)
      if (tolower(ext) == "cif") read.cif(p) else read.pdb(p)
    })
    
    coords <- lapply(structs, function(p) p$xyz)
    
    # Check same atom count
    lens <- lengths(coords)
    if (length(unique(lens)) != 1) {
      showNotification("❌ All structures must have the same number of atoms.", type = "error")
      return()
    }
    
    # RMSD matrix
    n <- length(coords)
    rmsd_mat <- matrix(0, n, n)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        rms <- rmsd(coords[[i]], coords[[j]])
        rmsd_mat[i, j] <- rms
        rmsd_mat[j, i] <- rms
      }
    }
    
    if (any(is.na(rmsd_mat))) {
      showNotification("❌ RMSD matrix contains NA values.", type = "error")
      return()
    }
    
    # MDS and Clustering
    mds <- cmdscale(rmsd_mat, k = 3)
    k <- input$k
    km <- kmeans(mds, centers = k, nstart = 20)
    
    # Save medoids
    tmpdir <- tempdir()
    rep_files <- character(k)
    for (cluster_id in 1:k) {
      inds <- which(km$cluster == cluster_id)
      submat <- rmsd_mat[inds, inds]
      medoid_idx <- inds[which.min(colMeans(submat))]
      rep_name <- paste0("cluster_", cluster_id, "_rep.cif")
      rep_path <- file.path(tmpdir, rep_name)
      file.copy(paths[medoid_idx], rep_path, overwrite = TRUE)
      rep_files[cluster_id] <- rep_path
    }
    
    # Save reactive results
    results$mds_df <- data.frame(X = mds[,1], Y = mds[,2], Z = mds[,3], Cluster = as.factor(km$cluster))
    results$km <- km
    results$rmsd_mat <- rmsd_mat
    results$paths <- paths
    results$rep_files <- rep_files
    
    # Cluster stats
    stats <- lapply(1:k, function(cluster_id) {
      inds <- which(km$cluster == cluster_id)
      submat <- rmsd_mat[inds, inds]
      upper <- submat[upper.tri(submat)]
      
      # Distancia al centroide en MDS
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
  
  output$mdsPlot <- renderPlotly({
    req(results$mds_df)
    plot_ly(data = results$mds_df, x = ~X, y = ~Y, color = ~Cluster, colors = "Set1",
            type = 'scatter', mode = 'markers') %>%
      layout(title = "MDS of Structures Colored by Cluster")
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
}

shinyApp(ui, server)
