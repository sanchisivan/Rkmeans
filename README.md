**Rkmeans** is an interactive Shiny application for structural clustering of molecular files (`.cif` or `.pdb`) based on RMSD and k-means clustering.

## 🚀 Features

- Upload multiple `.cif` or `.pdb` files
- Compute pairwise RMSD between structures
- Perform multidimensional scaling (MDS)
- Cluster structures using k-means
- Visualize clusters in interactive 2D plots (plotly)
- Download representative structures (medoids) for each cluster
- View summary statistics (cluster size, RMSD dispersion, centroid distances)

## 🧪 Requirements

- R ≥ 4.1
- Packages: `shiny`, `plotly`, `bio3d`, `ggplot2`

You can install required packages with:

```r
install.packages(c("shiny", "plotly", "ggplot2"))
# For bio3d
if (!requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d", repos = "http://bio3d.colorado.edu/CRAN/")
}
```

## 🖥️ Running the App

```r
shiny::runApp()
```

Or clone the repository and open `app.R` in RStudio.

## 📄 License

MIT License.

---

## 👤 Author

**Iván Sanchis, PhD**  
Laboratorio de Péptidos Bioactivos  
Facultad de Bioquímica y Ciencias Biológicas  
Universidad Nacional del Litoral  
Santa Fe, Argentina  
📧 sanchisivan@gmail.com / sanchisivan@fbcb.unl.edu.ar
