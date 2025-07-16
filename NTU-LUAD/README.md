# Group Repository for the Data Science Mini-Project (EMATM0050)

## Please edit the fields below with your information
Group Number: M-28

Problem Assigned:  Prediction of non-coding RNAs role in cancer through multi-omics data 
integration.

Group Members: Kiran Udayakumar, Shuojingrui He, Yibai Tang, Mansi Srivastava, Weilin He.

# Initial Data Description

This document provides an overview of the datasets available for the project **Prediction of non-coding RNAs role in cancer through multi-omics data integration**. The datasets include information on DNA methylation, RNA expression, and clinical phenotypes, sourced from TCGA, GTEx, and other bioinformatics platforms.

---

## Data Summary

### 1. Gene Attribute Edges
- **Filename**: `gene_attribute_edges.txt.gz`
- **Size**: 11.5 MB
- **Description**: 
  - A compressed text file describing relationships between transcription factors with target genes. Useful for building gene networks or pathway analysis. Each line in the data represents the relationship between a transcription factor (source) and its regulated target gene (target).
  - ***Values***:<br>1. `source`: the name of transcription factors.<br>2. `source_desc`: additional description of transcription factors.(NA means this field may not provide additional information)<br>3. `source_id`: ID of transcription factor.<br>4. `target`: the name of the target gene.<br>5. `target_desc`: additional description of the target gene.(NA means this field may not provide additional information)<br>6. `target_id`: ID of the target gene.<br>7. `weight`: the strength of the transcription factor's association with the target gene.
- ***Application***:<br>These data can be used to integrate the analysis of the relationship between transcription factors and target genes, specifically the differences between high and low MIR100HG expression groups in different cancer types.
- ***What we can do***:<br>1. Select the essencial ranscription factors and gene pairs according to the weight.<br>2. Combine the gene expression data and DNA methylation data to analyze whether these transcription factors play a key role in MIR100HG regulation.<br>3. Investigate the differences in the role of these genes in cancer and normal tissue.


    

---

### 2. Gene Annotation (hg19)
- **Filename**: `geneAnnotation_hg19_basicgenes.txt`
- **Size**: 159.6 MB
- **Description**: 
  - This table describes the information of the gene promoter region(promoter region based on human reference genome version hg19 annotation), which is often used to study the transcriptional regulation mechanism of genes, especially in the study of gene expression regulation. Each row of data provides detailed information about a promoter region, including its location on the chromosome, length, chain orientation, and corresponding genes and transcripts. 
  - ***Values***:<br>1. `chr`: the number of the chromosome where the promoter region is located. For example, chr1 means that the promoter is located on chromosome.<br>2. `start`: starting position of the promoter region. For example, 10874 indicates that the region begins at the 10,874th base of chromosome.<br>3. `end`: the end position of the promoter area. For example, 11873 indicates that the region ends at the 11873rd base.<br>4. `width`: the width or length of the promoter region. For example, 1000 indicates that the length of the promoter region is 1000 base pairs.<br>5. `strand`:the chain in which the promoter is located(where + is a positive chain and - is a negative chain)This information indicates the direction of transcription.<br>6. `id`: the unique identifier of the promoter region.<br>7. `tx_id`: a unique identifier representing the transcript associated with the promoter region.<br>8. `gene_id`: the unique identifier of the gene corresponding to the promoter.<br>9. `symbol`: the name(or symbol)of a gene.<br>10. `type`: the data type, all entries are hg19_genes_promoters, indicating that these promoters are promoter regions based on the human reference genome version hg19 annotation.
- ***What we can do***:<br>1. combination with gene expression data and DNA methylation data, we analyzed whether the methylation level of promoter region was significantly different between high and low MIR100HG expression groups.<br>2. To investigate whether MIR100HG affects gene expression by regulating the activity of certain gene promoters.<br>3. By using the transcription factor binding information of ENCODE database, we studied whether there were important transcription factor binding sites in these promoter regions.<br>4. Compare the differences in the properties of promoter regions (such as methylation levels or transcription factor binding) between cancer and normal tissues to look for specific regulatory patterns.










---

### 3. Illumina Probe Map
- **Filename**: `probeMap_illuminaMethylation450_GPL16304_TCGAlegacy`
- **Size**: 18.2 MB
- **Description**: 
  - This table describes data on gene methylation sites and may be used to analyze the relationship between DNA methylation and gene regulation. These data are important for integrating studies of DNA methylation patterns and gene expression, especially in the context of cancer.
  - ***Values***:<br>1. `id`: A unique identifier representing a methylation site, such as cg13332474. Each id corresponds to a specific methylation site (CpG island).<br>2. `gene`: the genes associated with this methylation site.(Multiple gene names are separated by commas, for example, RSPH14, and GNAZ indicates that the locus may regulate both genes. If there is`.`, the locus is not directly annotated to any gene.)<br>3. `chrom`: the chromosome where the methylation site is located, such as chr7.<br>4. `chromStart`: the starting position of the methylation site on the chromosome, for example 25935146.<br>5. `chromEnd`: the end position of the methylation site on the chromosome, such as 25935148.<br>6. `strand`: the chain information (positive or negative) of the methylation site.(If`.`Is displayed, the information is not marked)
- ***What we can do***:<br>1. Examine differences in methylation patterns between high and low MIR100HG expression groups. Analyze whether methylation specifically affects the expression of certain genes (e.g. RSPH14, GNAZ)<br>2. Correlate methylation sites with transcription factor binding sites (such as ENCODE data) to search for potential regulatory networks.<br>3. Compare the differences in the distribution of methylation sites in cancer and normal tissues to determine whether specific methylation patterns exist.




---

### 4. Survival Data
- **Filename**: `Survival_SupplementalTable_S1_20171025_xena_sp`
- **Size**: 2.4 MB
- **Description**: 
  - The table contains clinical information of cancer patients, mainly involving the basic characteristics of patients, tumor stage, treatment results, survival data, etc.
  - ***Values***:<br>***A. Patient & Demographic Information***<br>
`sample`: Unique sample ID, e.g., TCGA-OR-A5J1-01<br>
`_PATIENT`: Patient ID (without sample suffix)<br>
`cancer type abbreviation`: The type of cancer, in this case, ACC (Adrenocortical Carcinoma).<br>
`age_at_initial_pathologic_diagnosis`: Age when the patient was diagnosed.<br>
`gender`: Patient’s gender (MALE / FEMALE).<br>
`race`: Patient’s race (e.g., WHITE, BLACK OR AFRICAN AMERICAN).<br>
***B. Cancer Diagnosis & Staging***<br>
`ajcc_pathologic_tumor_stage`: The tumor stage according to AJCC (American Joint Committee on Cancer) classification (e.g., Stage II, Stage III, Stage IV).<br>
`clinical_stage`: Additional staging information (empty in this dataset).<br>
`histological_type`: Specific subtype of Adrenocortical Carcinoma (ACC).<br>
***C. Treatment & Disease Progression***<br>
`initial_pathologic_dx_year`: Year of diagnosis.<br>
`menopause_status`: Menopause status of female patients (empty in this dataset).<br>
`birth_days_to`: Number of days from birth to diagnosis (negative values).<br>
`vital_status`: Whether the patient is Alive or Dead.<br>
`tumor_status`: Whether the patient has a tumor at the last follow-up (WITH TUMOR or TUMOR FREE).<br>
`last_contact_days_to`: Days from diagnosis to the last contact.<br>
`death_days_to`: Days from diagnosis to death (if applicable).<br>
`cause_of_death`: If applicable, specifies the cause of death.<br>
`new_tumor_event_type`: If a new tumor event occurred (Distant Metastasis, Locoregional Recurrence).<br>
`new_tumor_event_site`: Location of tumor recurrence/metastasis (Peritoneal Surfaces, Soft Tissue, Lung).<br>
`new_tumor_event_dx_days_to`: Days from diagnosis to new tumor event.<br>
***D. Treatment Outcome & Prognosis***<br>
`treatment_outcome_first_course`: Initial treatment response (Complete Remission/Response, Progressive Disease).<br>
`margin_status`: Presence of cancer cells at the tumor margin after surgery (empty in this dataset).<br>
`residual_tumor`: Residual tumor presence after surgery (empty in this dataset).<br>
***E. Survival & Disease-Free Intervals(survival outcomes are represented in different time metrics)***:<br>
`OS` (Overall Survival): Whether the patient is alive (1 = deceased, 0 = alive).<br>
`OS.time`: Time (in days) from diagnosis to death or last follow-up.<br>
`DSS` (Disease-Specific Survival): Survival only considering death due to disease (1 = deceased, 0 = alive).<br>
`DSS.time`: Time (in days) from diagnosis to death caused by disease.<br>
`DFI` (Disease-Free Interval): Whether the patient remained tumor-free after treatment (1 = recurrence, 0 = no recurrence).<br>
`DFI.time`: Time (in days) before first recurrence.<br>
`PFI` (Progression-Free Interval): Similar to DFI, but includes disease progression cases.<br>
`PFI.time`: Time (in days) before disease progression.<br>
- ***What we can do***:<br>1. Stratifying Patients Based on MIR100HG expression: correlate survival outcomes, recurrence rates, and tumor progression with MIR100HG expression levels.<br>2. Identifying Prognostic Factors: Explore whether MIR100HG expression is linked to tumor stage, metastasis, and survival.<br>3. Analyzing Treatment Responses: Determine whether MIR100HG-high vs. low expression impacts response to treatments.(whether treatment_outcome_first_course (treatment outcome) is related to MIR100HG.)




---

### 5. TCGA Methylation Data (by Cancer Type)
#### Example Files:
- `TCGA.LUAD.sampleMap_HumanMethylation450.gz` (Lung Adenocarcinoma, 449.8 MB)
- `TCGA.PAAD.sampleMap_HumanMethylation450.gz` (Pancreatic Adenocarcinoma, 178.4 MB)
- **Description**:
  - Compressed datasets containing genome-wide methylation profiles from TCGA samples. Different files represent distinct cancer types.

---

### 6. GTEx Phenotype Data
- **Filename**: `TOIL_GTEX_PHENOTYPE_5_CANCERS.csv`<br>`TOIL_GTEX_RSEM_TPM_5_CANCERS.csv`
- **Size**: 172 KB
- **Description**:
  - Contains clinical and phenotypic details for samples from GTEx, focusing on five cancer types.

---

### 7. RNA Expression Data
#### Example Files:
- `TOIL_RSEM_TPM_LUAD.csv` (Lung Adenocarcinoma, 225.5 MB, 513 samples)
- `TOIL_RSEM_TPM_PAAD.csv` (Pancreatic ductual Adenocarcinoma, 74.7 MB, 178 samples)
- `TOIL_RSEM_TPM_PRAD.csv` (Prostate Adenocarcinoma, 206 MB, 495 samples)
- `TOIL_RSEM_TPM_SKCM.csv` (Skin cutaneous Adenocarcinoma, 43.5 MB, 102 samples)
- `TOIL_RSEM_TPM_STAD.csv` (Stomach Adenocarcinoma, 180.4 MB, 414 samples)
- **Description**:
  - The five files represent the 56404 gene expression data for five Cancer types(lung, pancreatic ductual, prostate, skin and stomach) from The Cancer Genome Atlas (TCGA) project, which was used to analyze gene expression patterns in different cancer types, specifically the role of MIR100HG in different cancers. The data is organized in the same way: the first column is the gene ID, the second column is the gene name, and the subsequent columns is the sample ID and its gene expression value. Expression values are usually positive (high expression) or negative (low expression/no expression)
  - ***Values***:<br>1. `Unnamed`: the ID of gene.<br>2. `HGNC_symbol`: the symbolic name of a gene.<br>3. `other columns`: each column corresponds to a patient sample with their corresponds gene expression value.
- ***What we can do***:<br>1. Compare the expression data of MIR100HG in different cancers.<br>2. Divide the high/low gene expression group for Differential Expression Analysis. <br>3. Combine the data of transcription factors - target gene data and DNA methylation data to explore the mechanism of MIR100HG.<br>4. Compare with GTEx normal tissue data, search for cancer-specific gene expression.


 
---

## Research Applications
These datasets enable a range of bioinformatics analyses, including:
1. **Non-coding RNA Prediction**:
   - Identify the role of ncRNAs in cancer through methylation and expression data.
2. **Survival Analysis**:
   - Correlate genomic features with patient survival outcomes.
3. **Multi-omics Integration**:
   - Combine RNA, methylation, and clinical data for holistic cancer analysis.
4. **Cross-Cancer Insights**:
   - Uncover common and unique molecular patterns across cancers.

For additional details on individual files or analysis workflows, please contact the group.

**Group Members**: Kiran Udayakumar, Shuojingrui He, Yibai Tang, Mansi Srivastava, Weilin He.
