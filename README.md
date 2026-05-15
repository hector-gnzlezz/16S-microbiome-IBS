# 16S rRNA Microbiome Analysis: Healthy vs IBS

Reproducible DADA2 pipeline for fecal microbiota profiling comparing a **healthy donor** with **IBS-C** (constipation-predominant) and **IBS-D** (diarrhea-predominant) patients using Illumina MiSeq 2×250 bp amplicon sequencing of the **V4 hypervariable region** of the 16S rRNA gene.

Developed as part of the MSc in Bioinformatics and Data Science (ISCIII · ENS · CNIO), session: *"Microbiome — 16S rRNA gene sequencing"*.

---

## Background

The gut microbiota has been consistently altered in IBS patients compared to healthy controls. Key findings from the literature include increased *Firmicutes* and reduced *Bacteroidetes* and *Bifidobacteria*, with specific changes at the family level (e.g., *Lachnospiraceae*, *Enterobacteriaceae*) depending on the IBS subtype.

This pipeline processes paired-end 16S amplicon reads through the full DADA2 workflow, from raw FASTQ files to taxonomic assignment and ecological diversity analysis with phyloseq.

---

## Pipeline Overview

```
Raw reads (FASTQ)
    │
    ├── Quality inspection        plotQualityProfile()
    │
    ├── Primer trimming           cutadapt (recommended, optional)
    │
    ├── Filter & trim             filterAndTrim(truncLen=c(240,160), maxEE=c(2,2))
    │
    ├── Error learning            learnErrors() — Forward & Reverse
    │
    ├── ASV inference             dada() — sample-level denoising
    │
    ├── Paired-end merging        mergePairs()
    │
    ├── ASV table                 makeSequenceTable()
    │
    ├── Chimera removal           removeBimeraDenovo(method="consensus")
    │
    ├── Taxonomic assignment      assignTaxonomy() + addSpecies()
    │                             ├── RDP Training Set 18 (default)
    │                             └── SILVA v138.1 (recommended, more RAM)
    │
    └── Ecological analysis       phyloseq + ggplot2
                                  ├── Alpha diversity (Shannon, Simpson, Chao1, Fisher)
                                  ├── Relative abundance bar plots (Phylum → Family)
                                  └── Heatmaps (Bray-Curtis NMDS ordination)
```

---

## Repository Structure

```
16S-microbiome-IBS/
├── extractASVs_Healthy_IBS.R    # Main DADA2 analysis script
├── env.yml           # Conda environment (all dependencies pinned)
├── microbiomas/
│   ├── inputFQfiles/            # Raw paired-end FASTQ files (not tracked)
│   │   ├── Healthy_F.fq.gz
│   │   ├── Healthy_R.fq.gz
│   │   ├── IBSC_F.fq.gz
│   │   ├── IBSC_R.fq.gz
│   │   ├── IBSD_F.fq.gz
│   │   └── IBSD_R.fq.gz
│   └── dbs/
│       └── RDP/                 # RDP Training Set 18 (download separately)
├── results/
│   ├── counts_16S_RDP.tsv       # ASV count table with taxonomy
│   ├── counts_s16S_sequences_RDP.tsv
│   └── figures/                 # Quality plots, diversity plots, heatmaps
└── README.md
```

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/hector-gnzlezz/16S-microbiome-IBS.git
cd 16S-microbiome-IBS
```

### 2. Create and activate the conda environment

All software dependencies are pinned in `microbiomaLola.yml` (R 4.3.3, DADA2 1.30, phyloseq 1.46, Python 3.12, Kraken2, BLAST, cutadapt, and more):

```bash
conda env create -f microbiomaLola.yml
conda activate microbiomaLola
```

### 3. Download reference databases

**RDP Training Set 18** (used by default — lower RAM requirements):
```bash
mkdir -p microbiomas/dbs/RDP
# Download rdp_train_set_18.fa.gz and rdp_species_assignment_18.fa.gz
# from: https://zenodo.org/record/4587955
```

**SILVA v138.1** (optional — recommended for species-level resolution):
```bash
mkdir -p microbiomas/dbs/silvaDB
# Download silva_nr99_v138.1_wSpecies_train_set.fa.gz
# from: https://zenodo.org/record/4587955
```

### 4. (Optional) Trim primers with cutadapt

```bash
cutadapt -a ^GTGCCAGCMGCCGCGGTAA \
         -A ^GGACTACHVGGGTWTCTAAT \
         -m 200 --discard-untrimmed \
         -o microbiomas/inputFQfiles/Healthy_F_trimmed.fq \
         -p microbiomas/inputFQfiles/Healthy_R_trimmed.fq \
         microbiomas/inputFQfiles/Healthy_F.fq.gz \
         microbiomas/inputFQfiles/Healthy_R.fq.gz
```

---

## Running the analysis

Open R or RStudio from within the activated conda environment:

```bash
conda activate microbiomaLola
Rscript extractASVs_Healthy_IBS.R
# or interactively: rstudio extractASVs_Healthy_IBS.R
```

Key parameters in the script:

| Parameter | Value | Rationale |
|---|---|---|
| `truncLen` | `c(240, 160)` | Fw: last 10 nt trimmed; Rv: quality crashes at pos 160 |
| `maxEE` | `c(2, 2)` | Max expected errors per read |
| `MAX_CONSIST` | `20` | Extended self-consistency loop for Rv convergence |
| `method` (chimera) | `"consensus"` | Recommended for multi-sample datasets |

---

## Key Results

| Step | Healthy | IBS-C | IBS-D |
|---|---|---|---|
| Raw reads (Fw) | 76,086 | 74,892 | 78,035 |
| After filtering | ~55,626 | ~54,964 | ~57,865 |
| After merging | ~51,851 | ~51,349 | ~53,998 |
| Final ASVs (no chimeras) | 446 total across 3 samples |

**Taxonomic findings** consistent with IBS literature:
- ↑ *Firmicutes* (Lachnospiraceae, Ruminococcaceae) in IBS patients
- ↑ *Proteobacteria* / Enterobacteriaceae in IBS-D
- ↓ *Bacteroidetes* overall in IBS vs Healthy

---

## Dependencies

| Tool | Version | Purpose |
|---|---|---|
| R | 4.3.3 | Analysis environment |
| DADA2 | 1.30.0 | ASV inference |
| phyloseq | 1.46.0 | Ecological analysis |
| ggplot2 | 3.5.2 | Visualization |
| cutadapt | 5.2 | Primer trimming |
| BLAST | 2.17.0 | Sequence alignment |

Full dependency list: see `microbiomaLola.yml`

---

## References

- Callahan BJ et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583.
- McMurdie PJ & Holmes S (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLOS ONE*, 8(4).
- Pittayanon R et al. (2019). Gut Microbiota in Patients With Irritable Bowel Syndrome — A Systematic Review. *Gastroenterology*, 157(1), 97–108.

---

## Author

**Héctor González Pérez** — [hectorgonzalezperez15@gmail.com](mailto:hectorgonzalezperez15@gmail.com)  
MSc Bioinformatics & Data Science · ISCIII · ENS · CNIO · Madrid
