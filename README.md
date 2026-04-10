# Madagascar Influenza

*Last updated: 2026-04-10*

**Report:** <https://edidpasteur.github.io/madagascar_influenza/>

## Contents

- [Goal](#goal)
- [Data](#data)
- [Alignment pipeline](#alignment-pipeline)
- [Phylogenetic tree pipeline](#phylogenetic-tree-pipeline)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Progress](#progress)

## Goal

Do Malagasy influenza viruses form distinct local lineages (monophyletic clades), or are they sporadic
importations from the broader African gene pool?
This determines whether sustained in-country surveillance in Madagascar is justified from a genomic-epidemiology perspective.

## Data

### GISAID (downloaded 16 March 2026)

**30,065 isolates** (all Africa, all hosts, after deduplication) — **1,849 from Madagascar (6.1%)**.
Raw GISAID files are not tracked by git; derived outputs (`combined_metadata.tsv`) are.

| Tag | Description | Isolates | Sequences |
|-----|-------------|----------|-----------|
| `a_human_africa` | Flu A human Africa | 19,902 | 110,901 |
| `b_human_africa` | Flu B human Africa | 5,081 | 27,900 |
| `animal_environmental_africa` | All animal/env Africa | 5,160 | 22,314 |
| **Total** | after deduplication | **30,065** | **161,143** |

**782 / 1,849 Madagascar isolates (42.3%) are analysis-ready** (complete 8-segment genome + known location + known date).

### Norosoa Razanajatovo (unpublished avian, 2021–2023)

109 avian isolates, 830 segment sequences from Institut Pasteur de Madagascar.
Files in `data/norosoa_avian_*.fasta` — not tracked by git.

| Subtypes | Isolates |
|----------|----------|
| H9N2 | 55 |
| H4 (partial) | 14 |
| H6N2, H6N8 | 13 |
| H1N2 | 8 |
| other | 19 |

Notes: 3 isolates already in GISAID (partial); 1 confirmed co-infection (H1N2/H6N1); 1 confirmed import from Brazil (H4N6).

## Alignment pipeline

**358 MAFFT alignments** — one per segment × subtype × scope (Africa / Madagascar), filtered to ≥ 3 sequences.
Input: 161,143 GISAID + 830 Norosoa sequences. Outputs: `alignments/aligned/`.

```bash
bash scripts/align_segments.sh          # submits SLURM jobs (resumable)
squeue -u $USER                         # monitor
```

| Script | Purpose |
|--------|---------|
| `scripts/align_segments.sh` | Splits sequences, submits one SLURM job per alignment |
| `scripts/benchmark_alignment.sh` | Benchmarks resource usage (4 test jobs) |
| `scripts/parse_benchmark.sh` | Parses benchmark results |

## Phylogenetic tree pipeline

**16 ML trees** (IQ-TREE 3.1.0, HKY+G, seed 42) — 5 HA + 5 NA per-subtype, 6 internal segments merged across subtypes.
Africa subsampled to 4 × n_Madagascar per tree (proportional by country, seed 42).

```bash
bash scripts/run_trees.sh               # submits all 16 jobs (resumable)
```

Trees: `trees/<NAME>.treefile`. Clade analysis results: `results/clade_summary.tsv`.

**Key results:**

| Pattern | Subtypes | Interpretation |
|---------|----------|----------------|
| 1–2 large clades, stem BL ~0.02–0.03 | **H6N2, H9N2** (avian) | Founder event → sustained local circulation; distinct from continental pool |
| 31–43 small clades, stem BL ~0.001–0.002 | **H1N1, H3N2, Flu B** (seasonal) | Annual reintroductions, no lasting local lineage |
| 51–71 fragmented clades | **Internal segments** | Reassortment signature |

Surveillance is justified: avian lineages (H6N2, H9N2) would be invisible to any other African country's programme.

| Script | Purpose |
|--------|---------|
| `scripts/run_trees.sh` | Subsamples inputs, submits IQ-TREE jobs |
| `scripts/run_one_tree.sh` | Single job: re-align (internal segments) + IQ-TREE |
| `scripts/analyse_trees.py` | Detects monophyletic Madagascar clades (ete3) |
| `scripts/render_trees.py` | Renders annotated PNGs (ETE3, headless) |

```bash
# Render tree images
mamba run -p ./env python scripts/render_trees.py --force
```

## Reproducing the analysis

```bash
# 1. Create environment
mamba env create -f environment.yaml

# 2. Download GISAID data → data/gisaid_epiflu_*.{xls,fasta}

# 3. Run
make analyse   # → data/combined_metadata.tsv
make report    # → docs/index.html
```

## Progress

- [x] GISAID data downloaded (30,065 isolates)
- [x] Metadata analysis + Quarto report published
- [x] Norosoa avian data (109 isolates) added
- [x] 358 alignments completed (MAFFT)
- [x] 16 ML trees completed (IQ-TREE 3.1.0, HKY+G)
- [x] Monophyletic clade analysis (`results/clade_summary.tsv`)
- [x] Report: Norosoa section + phylogenetic analysis section
- [ ] Render trees as annotated PNGs (`scripts/render_trees.py`)
- [ ] Reassortment analysis (H9N2, H3N2, H1N1, B)
