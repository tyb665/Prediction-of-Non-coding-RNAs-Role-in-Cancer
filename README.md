

# **Prediction of Non-coding RNAs Role in Cancer through Multi-omics Data Integration**

An integrative analysis to uncover the regulatory mechanisms of lncRNA **MIR100HG** across five cancer types using multi-omics data.

---

## 📑 Table of Contents

* [Project Overview](#project-overview)
* [Data Sources](#data-sources)
* [Installation](#installation)
* [Usage](#usage)
* [Main Features](#main-features)
* [Results Highlights](#results-highlights)
* [Future Work](#future-work)
* [Contributors](#contributors)
* [License](#license)

---

## 🔍 Project Overview

This project investigates the role of lncRNA **MIR100HG** in cancer by integrating gene expression, DNA methylation, and transcription factor–target data.
Five cancer types — **PRAD**, **PAAD**, **LUAD**, **SKCM**, and **STAD** — are systematically analyzed to identify shared and cancer-specific molecular regulatory mechanisms.

Machine learning models, including **Logistic Regression**, **Random Forest**, and **XGBoost**, are employed for prediction and feature discovery.

---

## 🧬 Data Sources

* **TCGA** — Tumor gene expression and methylation data
* **GTEx** — Normal tissue gene expression data
* **ENCODE** — Transcription factor (TF)–target regulatory network
* **MSigDB** — Molecular signature and pathway reference sets
* **MyGene.info API** — Gene symbol and annotation mapping

---

## ⚙️ Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/MIR100HG-multiomics.git
cd MIR100HG-multiomics

# (Optional) Create a virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install required packages
pip install -r requirements.txt
```

> **Main dependencies**: `pandas`, `numpy`, `scikit-learn`, `xgboost`, `lifelines`, `matplotlib`, `seaborn`

---

## 🚀 Usage

```bash
# Step 1: Preprocess multi-omics data
python scripts/data_preprocessing.py

# Step 2: Perform statistical analyses (DEG, DMA, TF identification)
python scripts/statistical_analysis.py

# Step 3: Train machine learning models
python scripts/model_training.py

# Step 4: Generate visualizations
python scripts/visualization.py
```

All results will be saved under the `/outputs` directory.

---

## 🌟 Main Features

* Multi-omics integration: expression, methylation, TF–target relationships
* Differential expression/methylation analysis between high/low MIR100HG groups
* Cancer vs. normal tissue comparative analysis
* Kaplan–Meier and Cox survival analysis
* Machine learning-based prediction and feature importance analysis
* Pathway enrichment and functional interpretation

---

## 📊 Results Highlights

* **Shared pathways**: Enriched pathways such as *Calcium signaling* and *PI3K–Akt signaling*
* **Key TFs**: Regulators like **POLR2A**, **EZH2**, and **CTCF** frequently identified
* **Methylation insights**: Cancer-specific MIR100HG promoter methylation profiles
* **Prognostic value**: MIR100HG serves as a prognostic biomarker in **STAD**

---

## 🔭 Future Work

* Expand analysis to more cancer types and external datasets (e.g., **GEO**, **ICGC**)
* Integrate more omics layers (e.g., **miRNA regulation**, **chromatin accessibility**)
* Investigate **multi-task learning** and **transfer learning** for model enhancement

---

## 👥 Contributors

* [Yibai Tang](mailto:xv24324@bristol.ac.uk)
* [Shuojingrui He](mailto:eo24594@bristol.ac.uk)
* [Weilin He](mailto:te23071@bristol.ac.uk)
* [Kiran Udayakumar](mailto:zt24069@bristol.ac.uk)
* [Mansi Srivastava](mailto:ah24844@bristol.ac.uk)

> University of Bristol, Department of Engineering

---

## 📄 License

This project is licensed under the **MIT License**. See the [LICENSE](./LICENSE) file for details.

---
