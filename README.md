# Madagascar Influenza

*Last updated: 2026-03-16*

## Table of Contents

- [Goal](#goal)
- [Resources](#resources)
- [Data](#data)
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

An isolate is **phylo-ready** when all 8 core gene segments (PB2, PB1, PA, HA, NP, NA, MP, NS) are present in the downloaded FASTA *and* each segment meets a minimum length threshold (~75–80% of the full-length reference). This ensures every segment has enough sequence for reliable phylogenetic reconstruction.

| Level | Count | % |
|-------|-------|---|
| Total isolates | 1,849 | 100% |
| Complete metadata (8 segment IDs in GISAID) | 782 | 42.3% |
| Phylo-ready (8 segments + sufficient length) | 782 | 42.3% |

Isolates that are not phylo-ready are listed in `data/madagascar_missing_metadata.tsv`.

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
- [x] Report published on QuartoPub
- [ ] Phylo-ready threshold validation against Nextstrain avian-flu pipeline
- [ ] GitLab Pages deployment
