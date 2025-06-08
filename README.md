
# Quick Start Guide

> 🎉 **Great News!**  
> This repository comes with all GWAS data pre-loaded! No downloading, no waiting—just clone and run for instant access to all Manhattan plots.

---

## Step 1: Quick Setup (2 minutes)

Get up and running in just a few commands:

```bash
# 1) Install & initialize Git LFS (one-time per machine)
brew install git-lfs        # macOS (or see https://git-lfs.github.com/)
git lfs install

# 2) Clone & fetch LFS-tracked .db files
git clone https://github.com/ciarareeve/STRXplorer.git
cd STRXplorer
git lfs pull

# 3) Install Python dependencies
pip install -r requirements.txt

# 4) Run the application
python STRXplorer.py
````

> **Success!**
> The app will start at [http://localhost:5000](http://localhost:5000) with full functionality.

---

## Step 2: What’s Included

This repository comes with everything you need:

```
STRXplorer/
├── STRXplorer.py             ← Main Flask application
├── locus_data.db             ← STR locus information
├── manhattan_data.db         ← Pre-loaded GWAS data (LFS-stored)
├── manhattan_plot.py         ← Plot generation utilities
├── locus_plots.py            ← Additional plotting functions
├── Procfile                  ← EB deployment instruction
├── requirements.txt          ← Python deps
├── .gitattributes            ← Git LFS tracking rules
└── templates/
    ├── home.html
    ├── error.html
    ├── browse_traits.html
    ├── browse_loci.html
    ├── trait_overview.html
    └── … other templates
```

---

## Step 3: Using the Platform

Once running, you can immediately:

1. **Browse Traits:**
   [http://localhost:5000/browse\_traits](http://localhost:5000/browse_traits)
2. **Browse STR Loci:**
   [http://localhost:5000/browse\_loci](http://localhost:5000/browse_loci)
3. **Generate Manhattan Plots:**
   Click any trait in the list to see its Manhattan plot
4. **Check Status:**
   [http://localhost:5000/database\_status](http://localhost:5000/database_status)

```bash
# Example URLs
http://localhost:5000/
http://localhost:5000/browse_traits
http://localhost:5000/browse_loci
http://localhost:5000/trait_overview/mean_platelet_volume
http://localhost:5000/database_status
```

---

## Step 4: System Requirements

* **Git LFS**

  * Install: `brew install git-lfs` or see [https://git-lfs.github.com/](https://git-lfs.github.com/)
  * Initialize: `git lfs install`
* **Python** ≥ 3.8
* **RAM:** ≥ 2 GB (4 GB recommended)
* **Disk:** ≥ 1 GB free
* **Web Browser:** Chrome, Firefox, Safari, Edge

---

## Step 5: Pre-loaded Data Overview

The **manhattan\_data.db** contains:

* **Multiple traits:** Blood traits, anthropometric measures, and more
* **Millions of variants:** Genome-wide association data
* **Statistical summaries:** P-values, effect sizes, confidence intervals
* **Optimized indexes:** Fast querying for real-time plotting

> **Data Source:**
> GWAS data from Margoliash *et al.* (2023) study on STR associations with complex traits.

---

## Step 6: Troubleshooting

### Application won’t start

* Check Python version:

  ```bash
  python --version
  ```
* Install missing packages:

  ```bash
  pip install -r requirements.txt
  ```
* Ensure you’re in the project root directory

### Database errors

* Verify both `.db` files are present (`locus_data.db`, `manhattan_data.db`)
* Check file permissions (should be world-readable)
* Run from the repo root

### No plots showing

* Visit the [Database Status](/database_status) page
* Refresh your browser or clear cache

### Port already in use

* Kill existing Flask processes:

  ```bash
  pkill -f flask
  ```
* Or run on a different port:

  ```bash
  flask run --port 5001
  ```

---

## Git LFS Notes

We use **Git Large File Storage** to host the two big `.db` files:

1. **Track** them in your local clone:

   ```bash
   git lfs track "*.db"
   git add .gitattributes
   git commit -m "Track .db files with Git LFS"
   ```
2. **Clone & pull** as shown in Step 1 to fetch the actual database blobs

Make sure **`.gitattributes`** contains:

```gitattributes
*.db filter=lfs diff=lfs merge=lfs -text
```

---

That’s it — happy exploring!

