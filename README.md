# Madagascar Influenza

*Last updated: 2026-03-31*

**Report: https://edidpasteur.github.io/madagascar_influenza/**

## Table of Contents

- [Goal](#goal)
- [Resources](#resources)
- [Data](#data)
  - [GISAID data](#gisaid-data)
  - [Unpublished avian data (Norosoa)](#unpublished-avian-data-norosoa)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Progress](#progress)

## Goal

Characterize influenza sequences from Madagascar and Africa available in GISAID:

- Which influenza sequences exist from Madagascar and Africa in GISAID?
- What is their geographic location and sampling time?
- Are the genomes complete?
- Which isolates are suitable for phylogenetic analysis (all 8 segments present and of sufficient length)?

## Resources

- **GISAIDR** — R package for programmatic access to GISAID: https://github.com/Wytamma/GISAIDR
  - Note: GISAIDR does not support EpiFlu. Data was downloaded manually from the GISAID EpiFlu interface.
- **Nextstrain avian-flu**: https://nextstrain.org/avian-flu/h5n1/ha/2y
- **Nextstrain avian-flu pipeline**: https://github.com/nextstrain/avian-flu

## Data

### GISAID data

Downloaded from GISAID EpiFlu (manual download, 16 March 2026). Raw GISAID files are in `data/` and are **not tracked by git** (GISAID redistribution policy). Derived outputs (`combined_metadata.tsv`, `madagascar_missing_metadata.tsv`) **are tracked**.

### Metadata files (XLS)

| File tag | Description | Isolates |
|----------|-------------|----------|
| `a_human_africa` | Influenza A — human host — Africa | 19,902 |
| `b_human_africa` | Influenza B — human host — Africa | 5,081 |
| `c_human_africa` | Influenza C — human host — Africa | 10 |
| `animal_environmental_africa` | All animal/environmental hosts — Africa | 5,160 |
| **Total** | After cross-file deduplication | **30,065** |

### FASTA sequence files

The Influenza A metadata is split into two FASTA downloads by subtype (GISAID export limit).

| File tag | Description | Sequences |
|----------|-------------|-----------|
| `a_excepth3_human_africa` | Influenza A (excl. H3) — human host — Africa | 49,494 |
| `a_h3_human_africa` | Influenza A H3 — human host — Africa | 61,407 |
| `b_human_africa` | Influenza B — human host — Africa | 27,900 |
| `c_human_africa` | Influenza C — human host — Africa | 28 |
| `animal_environmental_africa` | All animal/environmental hosts — Africa | 22,314 |
| **Total** | Unique sequences across all files | **161,143** |

### Summary

**30,065 isolates (after deduplication) — 1,849 from Madagascar (6.1% of all African isolates)**

### Genome completeness (Madagascar)

An isolate is **analysis-ready** when it meets all three criteria:

1. **Complete genome (phylo-ready)** — all 8 segments (PB2, PB1, PA, HA, NP, NA, MP, NS) present in the downloaded FASTA and each above minimum length thresholds (~75–80% of full-length reference; e.g. 1,800 nt for PB2/PB1, 1,300 nt for HA, 700 nt for NS). Thresholds are defined in `MIN_SEQ_LEN` in `scripts/analyse_gisaid.py` and have not yet been validated against a specific downstream pipeline (see Progress).
2. **Known geographic location** — Country field non-empty
3. **Known sampling date** — Collection date present and plausible (1990–2026). All 30,065 African isolates (including all 1,849 from Madagascar) have full year-month-day precision; 0 are year-only.

| Level | Count | % |
|-------|-------|---|
| Total isolates | 1,849 | 100% |
| Complete genome (phylo-ready) | 782 | 42.3% |
| Known location | 1,849 | 100% |
| Known date | 1,849 | 100% |
| **Analysis-ready (all 3 criteria)** | **782** | **42.3%** |

Isolates that are not phylo-ready are listed in `data/madagascar_missing_metadata.tsv`.

### Unpublished avian data (Norosoa)

Unpublished avian influenza sequences from Madagascar collected by Norosoa Raharinosy (Institut Pasteur de Madagascar, 2021–2023). Provided as 7 FASTA files organised by HA/NA subtype, covering **109 unique isolates** (830 segment sequences in total). Raw files are in `data/` and are **not tracked by git**.

| File | Subtype(s) | Isolates | Sequences |
|------|-----------|----------|-----------|
| `norosoa_avian_h1_madagascar.fasta` | H1N2 | 8 | 62 |
| `norosoa_avian_h4_madagascar.fasta` | H4 (partial) | 14 | 109 |
| `norosoa_avian_h6_madagascar.fasta` | H6N2, H6N8 | 13 | 103 |
| `norosoa_avian_h7_madagascar.fasta` | H7N7 | 3 | 24 |
| `norosoa_avian_h9_madagascar.fasta` | H9N2 | 55 | 427 |
| `norosoa_avian_n2_madagascar.fasta` | N2 (partial) | 12 | 77 |
| `norosoa_avian_n6_madagascar.fasta` | N6 (partial) | 4 | 28 |
| **Total** | | **109** | **830** |

**Notes on data quality (confirmed by Norosoa):**

- **3 isolates already in GISAID** (IDs: 78609-23, 78830-23, 79149-23 — H9N2, July 2023). The GISAID entries are partial (5 segments, not phylo-ready); these FASTA files contain the complete genomes sequenced by CDC.
- **20 isolates partially typed** (e.g. "N2", "H9"): incomplete sequencing — full subtype cannot be determined.
- **1 co-infection** (sample H1N2/H6N1): confirmed double infection, not a data entry error.
- **1 imported case** (sample 80824-23, Chicken, H4N6, 22 Aug 2023, region "Brésil"): confirmed import from Brazil — part of punctual testing of imported day-old chicks alongside routine market surveillance.

## Reproducing the analysis

### 1. Create the environment

```bash
mamba env create -f environment.yaml
conda activate madagascar_influenza
```

`seqkit` is included in `environment.yaml` (from the `bioconda` channel) — no separate module load is needed.

### 2. Download data from GISAID

Download the XLS metadata exports and FASTA sequence files from [GISAID EpiFlu](https://gisaid.org) and place them in `data/`, following the naming convention:

```
data/gisaid_epiflu_isolates_<tag>.xls
data/gisaid_epiflu_sequence_<tag>.fasta
```

### 3. Run the pipeline

```bash
make analyse   # → data/combined_metadata.tsv
make report    # → docs/index.html
```

## Progress

- [x] GISAID access obtained
- [x] All Africa datasets downloaded (30,065 isolates, all hosts)
- [x] Conda environment created (`madagascar_influenza`)
- [x] Metadata analysis: geographic distribution, sampling timeline, genome completeness
- [x] Quarto report with Plotly interactive figures (`docs/index.html`)
- [x] Missing metadata fields identified for all 1,849 Madagascar isolates
- [x] `analyse_gisaid.py` hardened: logging, argparse, vectorised ops, provenance JSON
- [x] Report published on GitHub Pages (mirrored from GitLab)
- [x] Unpublished avian sequences from Norosoa (109 isolates, 830 segments) added to `data/`
- [ ] Integrate Norosoa sequences into phylogenetic analysis
- [ ] Phylo-ready threshold validation against Nextstrain avian-flu pipeline
