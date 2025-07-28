# Rkmeans

**Rkmeans** is a small interactive R/Shiny application for structural clustering of molecular files (`.cif` or `.pdb`) based on RMSD and k-means clustering.

This app was born from the need to cluster AlphaFold 3 output structures, especially for R users. At the time, no simple and direct online tool allowed structural clustering of `.cif` files.  
This tool is best suited for **clustering a small number of predicted structures** (e.g., 5–100 models), adjusting the number of clusters (`k`) accordingly.  
For example:
- If you input **10 structures**, trying `k = 2–4` is sensible.  
- For **30 structures**, you might test `k = 3–6`.  
- For **50+ structures**, `k` = 10 may be appropriate.

---

## Features

- Load multiple `.pdb` or `.cif` files from a folder.
- Automatically convert AlphaFold `.cif` files to `.pdb` (via Python and Biopython).
- Calculate pairwise RMSD using alpha carbons (Cα).
- Run multidimensional scaling (MDS) and k-means clustering.
- Visualize 3D structures with `r3dmol` directly in the app.
- Show cluster details (size, representative structure).
- Save cluster representatives.

---

## Requirements

### R packages

Install these packages if you don’t have them already:

```r
install.packages(c(
  "shiny", "shinyFiles", "bio3d", "r3dmol", "tools", 
  "plotly", "dplyr", "readr", "tibble", "ggplot2"
))
```

### (Optional) Python + Biopython

If you want to use the **3D structure viewer** with `.cif` files (especially those downloaded from AlphaFold), you'll need Python and Biopython installed. This is only required for visualization — **not for clustering**.

```bash
pip install biopython
```

Make sure Python is installed and accessible from R via:

```r
Sys.which("python")
```

If you're only interested in the clustering results and don't need structure visualization, you can skip this step.

---

## How to Use

1. Clone this repository or download the app files.
2. Place your `.cif` or `.pdb` files in a folder (e.g., `structures/`).
3. Launch the app in R:

```r
shiny::runApp("path_to_app_folder")
```

4. In the app:
   - Select a folder containing `.cif` or `.pdb` files.
   - The app will align the structures, calculate the RMSD matrix, and perform **k-means clustering**.
   - You can choose the number of clusters (`k`) manually.
   - Once clustering is complete, you'll see:
     - A **summary table** showing cluster sizes, within-cluster RMSD, and medoid filenames.
     - A **3D scatter plot** of the structures in reduced dimensions (via MDS).
     - Optional **3D viewer** for inspecting any selected structure.

---

## Just a Simple App

This is not a polished piece of software — just a tool I made for myself to cluster structural models. Sharing it here in case it helps others using R for protein analysis.

---

## License

MIT License.

---

## 👤 Author

**Ivan Sanchis, PhD**  
Laboratorio de Péptidos Bioactivos  
Facultad de Bioquímica y Ciencias Biológicas  
Universidad Nacional del Litoral  
Santa Fe, Argentina  
📧 sanchisivan@gmail.com / sanchisivan@fbcb.unl.edu.ar
