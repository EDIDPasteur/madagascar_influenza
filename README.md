# Madagascar Influenza

*Last updated: 2026-04-10*

**Report: https://edidpasteur.github.io/madagascar_influenza/**

## Table of Contents

- [Goal](#goal)
- [Data](#data)
  - [GISAID data](#gisaid-data)
  - [Unpublished avian data (Norosoa)](#unpublished-avian-data-norosoa)
- [Sequence alignment pipeline](#sequence-alignment-pipeline)
- [Phylogenetic tree pipeline](#phylogenetic-tree-pipeline)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Progress](#progress)

## Goal

Characterize influenza sequences from Madagascar and Africa available in GISAID:

- Which influenza sequences exist from Madagascar and Africa in GISAID?
- What is their geographic location and sampling time?
- Are the genomes complete?
- Which isolates are suitable for phylogenetic analysis (all 8 segments present and of sufficient length)?
- **Does Madagascar sequencing capture unique African genetic diversity?** Do Malagasy viruses form distinct local lineages (monophyletic clades), or are they sporadic importations from the broader African gene pool? This determines whether sustained in-country influenza surveillance in Madagascar is justified from a genomic-epidemiology perspective.


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

Unpublished avian influenza sequences from Madagascar collected by Norosoa Razanajatovo (Institut Pasteur de Madagascar, 2021–2023). Provided as 7 FASTA files organised by HA/NA subtype, covering **109 unique isolates** (830 segment sequences in total). Raw files are in `data/` and are **not tracked by git**.

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

## Sequence alignment pipeline

### Why per-segment, never concatenated

Influenza has 8 independent genomic segments that can **reassort** — meaning different segments can have different evolutionary histories within the same host.

### What is being aligned

Each alignment is a **single influenza segment from a single subtype**. Mixing subtypes is avoided because highly variable segments like HA and NA diverge too much across subtypes to be meaningfully aligned together. Internal segments (PB2, PB1, PA, NP, MP, NS) are more conserved but are also kept per-subtype at the alignment stage; for tree building they are merged across all subtypes into one tree per segment (see Phylogenetic tree pipeline).

Each segment × subtype combination is aligned in **two scopes**:

- **Africa** — all African sequences from GISAID + Norosoa's Madagascar sequences
- **Madagascar** — Madagascar sequences only (GISAID + Norosoa)

This gives a maximum of 16 alignments per subtype (8 segments × 2 scopes), reduced to those with ≥ 3 sequences.

### Inputs

| Source | Files | Sequences |
|--------|-------|-----------|
| GISAID EpiFlu | `data/gisaid_epiflu_sequence_*.fasta` (5 files) | 161,143 |
| Norosoa (unpublished) | `data/norosoa_avian_*.fasta` (7 files) | 830 |

**Phase 1** (`scripts/align_segments.sh`): the Python block reads all input FASTA files, splits them by segment × subtype × scope, and writes intermediate files to `alignments/split/`. Non-IUPAC characters are replaced with `N` at this stage (228 characters corrected across all GISAID sequences).

### Outputs

**358 aligned FASTA files** in `alignments/aligned/<SEG>_<SUBTYPE>.<scope>.aln.fasta` (425 initial split files, 67 removed after filtering sequences with partial/unknown subtypes — only `HxNy`, `B`, and `C` are kept).

For example:
```
alignments/aligned/HA_H9N2.madagascar.aln.fasta   # HA segment, H9N2, Madagascar only
alignments/aligned/HA_H3N2.africa.aln.fasta        # HA segment, H3N2, all Africa
alignments/aligned/PB2_H5N1.africa.aln.fasta       # PB2 segment, H5N1, all Africa
```

### Scripts

| Script | Purpose |
|--------|---------|
| `scripts/align_segments.sh` | Main pipeline: splits sequences, submits one SLURM job per alignment |
| `scripts/benchmark_alignment.sh` | Submits 4 test jobs (one per resource tier) to measure real RAM and runtime |
| `scripts/parse_benchmark.sh` | Parses benchmark results via `reportseff` and `/usr/bin/time -v` |

### Resource tiers (calibrated from benchmarks on HA segment)

| Tier | Sequences | MAFFT mode | CPUs | Memory | Measured peak RAM | Jobs |
|------|-----------|------------|------|--------|-------------------|------|
| tiny | ≤ 500 | `--localpair --maxiterate 1000` (L-INS-i) | 4 | 2G | 168 MB | ~271 |
| small | 501–2,000 | `--auto` | 4 | 2G | 186 MB | ~11 |
| medium | 2,001–8,000 | `--auto` | 8 | 4G | 468 MB | ~16 |
| large | > 8,000 | `--auto` | 16 | 16G | 5.1 GB | ~3 |

All jobs run on the `seqbio` partition with no time limit. Alignments already present in `alignments/aligned/` are skipped (resumable). Job logs are written to `alignments/logs/`.

### Running the pipeline

```bash
# (Optional) benchmark resource usage first
bash scripts/benchmark_alignment.sh
bash scripts/parse_benchmark.sh   # after jobs finish

# Submit all 358 alignment jobs in batches of 20
bash scripts/align_segments.sh --batch-size 20
# Re-run the same command to submit the next batch (completed alignments are skipped)

# Or submit all at once
bash scripts/align_segments.sh

# Monitor
squeue -u $USER
reportseff $(sacct -u $USER --format=JobID --noheader -S today | tr '\n' ',')
```

## Phylogenetic tree pipeline

### Tree design

**16 ML trees** (IQ-TREE 3.1.0, HKY+G, seed 42) are built to assess whether Madagascar sequences form **monophyletic clades** — groups of ≥ 2 Madagascar sequences sharing a common ancestor exclusive to Madagascar, indicating locally-sustained lineages — or are scattered singletons within an otherwise African tree, indicating sporadic importation. Trees are filtered to datasets where Madagascar has ≥ 10 sequences and Africa has ≥ 20. Two strategies apply by segment:

| Category | Segments | Tree scope | # Trees |
|----------|----------|------------|---------|
| **HA per-subtype** | HA | One tree per HA subtype | 5 |
| **NA per-subtype** | NA | One tree per NA subtype | 5 |
| **Internal merged** | PB2, PB1, PA, NP, MP, NS | All subtypes merged — one tree per segment | 6 |
| **Total** | | | **16** |

**HA and NA** are kept per-subtype because HA/NA subtypes represent distinct antigenic lineages separated by long branches — mixing them would produce an unresolvable star topology. **Internal segments** (PB2, PB1, PA, NP, MP, NS) do not cluster by HxNy subtype (Poon, *Virus Evolution* 2024); merging all subtypes into one tree per segment recovers the full evolutionary history of each internal gene, including reassortment signatures.

**Subsampling** (all trees, seed 42):
- **Madagascar**: all sequences kept
- **Africa (non-Madagascar)**: proportionally sampled by country to **4 × n_Madagascar** — large enough that the tree reflects African diversity, small enough for IQ-TREE to run in hours–days

**Filter**: trees with Mdg < 10 or Africa < 20 are skipped (56 subtypes excluded — mostly avian subtypes with no Madagascar sequences, e.g. H5N8, H7N6).

### Tree inventory

| Tree | DB total | Mdg (all) | Africa (all) | Mdg in tree | Africa in tree | Tree size | % Mdg |
|------|----------|-----------|--------------|-------------|----------------|-----------|-------|
| HA_B | 5,015 | 365 | 4,650 | 365 | 1,492 | 1,857 | 19.7% |
| HA_H1N1 | 8,621 | 234 | 8,387 | 234 | 1,027 | 1,261 | 18.6% |
| HA_H3N2 | 9,587 | 375 | 9,212 | 375 | 1,557 | 1,932 | 19.4% |
| HA_H6N2 | 64 | 27 | 37 | 27 | 37 | 64 | 42.2% |
| HA_H9N2 | 1,438 | 71 | 1,367 | 71 | 303 | 374 | 19.0% |
| NA_B | 4,451 | 342 | 4,109 | 342 | 1,395 | 1,737 | 19.7% |
| NA_H1N1 | 7,959 | 203 | 7,756 | 203 | 894 | 1,097 | 18.5% |
| NA_H3N2 | 8,965 | 369 | 8,596 | 369 | 1,526 | 1,895 | 19.5% |
| NA_H6N2 | 65 | 28 | 37 | 28 | 37 | 65 | 43.1% |
| NA_H9N2 | 1,049 | 71 | 978 | 71 | 299 | 370 | 19.2% |
| MP | 20,510 | 826 | 19,684 | 826 | 3,440 | 4,266 | 19.4% |
| NP | 18,451 | 907 | 17,544 | 907 | 3,733 | 4,640 | 19.5% |
| NS | 19,568 | 907 | 18,661 | 907 | 3,747 | 4,654 | 19.5% |
| PA | 16,970 | 773 | 16,197 | 773 | 3,205 | 3,978 | 19.4% |
| PB1 | 16,148 | 687 | 15,461 | 687 | 2,863 | 3,550 | 19.4% |
| PB2 | 16,710 | 736 | 15,974 | 736 | 3,057 | 3,793 | 19.4% |

> **DB total**: sequences in the full Africa alignment before subsampling.  
> **Africa in tree**: sampled proportionally by country (4 × n_Mdg, seed 42).  
> H6N2 is small enough that all sequences are kept (no subsampling applied).  
> Internal segment trees are built from unaligned `alignments/split/` files merged across all subtypes; MAFFT re-aligns them before IQ-TREE.

### Scripts

| Script | Purpose |
|--------|---------|
| `scripts/run_trees.sh` | Subsamples inputs (Phase 1A: HA/NA from aligned files; Phase 1B: internal segments merged from split files), then submits all IQ-TREE jobs via `parallel + srun` |
| `scripts/run_one_tree.sh` | Runs one job: MAFFT re-alignment (internal segments only) followed by IQ-TREE, via `srun` on the `seqbio` partition |

### Running the pipeline

```bash
# Submit all 16 tree jobs
bash scripts/run_trees.sh

# Or in batches
bash scripts/run_trees.sh --batch-size 8

# Monitor
squeue -u $USER
```

Trees are written to `trees/<NAME>.treefile`. Subsampled FASTAs (and MAFFT re-alignments for internal segments) are in `trees/subsampled/`. Completed trees are skipped on re-run (resumable).

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
- [x] All sequences split by segment × subtype × scope; 67 partial-subtype combos excluded
- [x] **358 per-segment alignments completed** (MAFFT, resource-tiered, via `parallel+srun`)
- [x] **16 ML trees completed** (IQ-TREE 3.1.0, HKY+G — 5 HA + 5 NA per-subtype, 6 internal segments merged; all 16 treefiles verified April 9 2026)
- [x] **Monophyletic clade analysis completed** (`scripts/analyse_trees.py`, ete3) — results in `results/clade_summary.tsv` and `results/clade_details.tsv`; 74–100 % of Madagascar sequences form local clades depending on subtype
- [x] Report extended: Norosoa data section + phylogenetic analysis section
- [ ] Render trees as annotated PNG/HTML for visual inspection (ETE3 or Auspice)
- [ ] Assess reassortment signal across segments for key subtypes (H9N2, H3N2, H1N1, B)
