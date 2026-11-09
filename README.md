# DASH4D gut microbiome analysis

Code supporting the DASH4D gut microbiome study. The publication citation and DOI will be added when available.

**Data availability:** Raw metagenomic sequencing reads are deposited under NCBI BioProject [PRJNA1505569](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1505569). Study/clinical metadata and processed analysis data are available from the corresponding author upon reasonable request. Data files and unpublished preprocessing scripts are not distributed in this repository.

## Repository structure

```text
analysis/
  DASH4D_microbiome_analysis.Rmd    # Microbiome, SCFA, and functional-gene analyses
  model_helpers.R                   # Reusable statistical modeling functions
environment/
  R_sessionInfo.txt                 # R and system information
  R_package_versions.csv            # R package versions
  Python_environment.txt            # Python and system information
  Python_requirements.txt           # Python package versions
figures/
resources/
  butyrate_candidates_with_literature.tsv  # Host and environment for Butyrate-producing strains
utils/
  rarefaction_depth_mp.py           # Rarefaction and alpha diversity
visualization/
  gutmicrobiome_visualization.ipynb # Main and supplementary figures
```

## Metagenomic sequencing and functional profiling

Read preprocessing and profiling were conducted as described in the manuscript; upstream processing scripts are not included here.

1. **Quality control and host removal:** Reads were processed with **Sunbeam v4.6.0** using **Cutadapt v2.6** for adapter trimming, **Trimmomatic v0.33** for quality filtering, and **komplexity v0.3.6** to remove reads with complexity scores below 0.55. **BWA-MEM v0.7.17** removed reads mapping to **GRCh38** or **phiX174** (minimum 50% sequence identity and 60% aligned).
2. **Taxonomic profiling:** **Kraken2 v2.1.3** classified reads against **UHGG v2.0.2**, using `--minimum-hit-groups 4` and `--confidence 0.05`. Species abundance was estimated using **Bracken v3.1**; samples with fewer than 1,000,000 classified reads were excluded.
3. **Functional profiling:** **DIAMOND v2.1.12** aligned host-removed reads to butyrate synthesis gene references and **KEGG release 89.0** (`evalue: 1e-5`, `alnLen: 30`, `mismatch: 3`). The butyrate reference is available [here](https://raw.githubusercontent.com/zhaoc1/sunbeam_databases/master/dbs). Butyrate gene analyses used the **human/animal–GI** subset (`Isolate_origin == "human/animal"`, `Isolation_site == "GI"`), annotated in `resources/butyrate_candidates_with_literature.tsv`.

## Alpha diversity rarefaction

`utils/rarefaction_depth_mp.py` repeatedly subsamples integer feature counts without replacement and computes alpha diversity metrics using scikit-bio. The input is a TSV with feature IDs in rows and samples in columns by default.

```bash
python utils/rarefaction_depth_mp.py \
  --input_file /path/to/counts.tsv \
  --output_file /path/to/alpha_diversity_rarefaction.tsv \
  --depths 1000,5000,10000 \
  --metrics observed_otus,shannon,simpson \
  --repeats 10 \
  --seed 42 \
  --processes 4
```

Use `--sample-as-rows` for transposed input. Samples below a requested depth are skipped. Run `python utils/rarefaction_depth_mp.py --help` for additional options.

## Sensitivity analysis

Sensitivity analyses are included in `analysis/DASH4D_microbiome_analysis.Rmd` and apply the primary modeling approaches after excluding the relevant participants or visits, as described in the manuscript. No separate sensitivity-analysis script is required.
