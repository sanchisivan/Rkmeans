**Rkmeans** is an interactive Shiny application for structural clustering of molecular files (`.cif` or `.pdb`) based on RMSD and k-means clustering.

This app was born from the need to cluster AlphaFold 3 output structures. At the time, no simple and direct online tool allowed the structural clustering of `.cif` files. **Rkmeans** provides an intuitive and local interface to process and visualize clusters of predicted structures.

## 🚀 Features

- Upload multiple `.cif` or `.pdb` files
- Compute pairwise RMSD between structures
- Perform multidimensional scaling (MDS)
- Cluster structures using k-means
- Visualize clusters in interactive 2D or 3D MDS plots (Plotly)
- View and select representative structures for each cluster
- Convert `.cif` representative structures to `.pdb` format (for visualization compatibility)
- Visualize selected structures in 3D (r3dmol)
- Download representative structures (medoids) for each cluster
- Download summary statistics (cluster size, RMSD dispersion, centroid distances)

## 🧪 Requirements

- R ≥ 4.1
- R packages: `shiny`, `plotly`, `bio3d`, `ggplot2`, `r3dmol`, `shinybusy`
- Python ≥ 3.8
- Python packages: `biopython`

### Installing R packages

```r
install.packages(c("shiny", "plotly", "ggplot2", "shinybusy"))
# For bio3d
if (!requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d", repos = "http://bio3d.colorado.edu/CRAN/")
}
# For r3dmol
if (!requireNamespace("r3dmol", quietly = TRUE)) {
  install.packages("r3dmol")
}
```

### Installing Python and Biopython (Windows)

1. Download Python from: https://www.python.org/downloads/windows/
2. During installation, check **"Add Python to PATH"**.
3. Open a terminal (cmd) and install Biopython:

```bash
pip install biopython
```

## ⚠️ Important

In the file `app.R`, update the following line with the correct path to your Python executable (only required if not added to PATH):

```r
python_exec <- "C:/path/to/python.exe"
```

If Python is in your system PATH, you can simplify this to just `"python"`.

## 💥 Running the App

From R or RStudio:

```r
shiny::runApp()
```

Or clone the repository and open `app.R`.

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
