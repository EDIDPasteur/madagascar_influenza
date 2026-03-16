# Madagascar Influenza

*Last updated: 2026-03-16*

## Goal

Characterize influenza sequences from Madagascar and Africa available in GISAID:

- Which influenza sequences exist from Madagascar and Africa in GISAID?
- What is their geographic location, sampling time, and whether the genomes are complete?
- Which isolates are suitable for phylogenetic analysis (all 8 segments present and of sufficient length)?

## Resources

- **GISAIDR** — R package for programmatic access to GISAID: https://github.com/Wytamma/GISAIDR
  - Note: GISAIDR does not support EpiFlu. Data was downloaded manually from the GISAID EpiFlu interface.
- **Nextstrain avian-flu**: https://nextstrain.org/avian-flu/h5n1/ha/2y
- **Nextstrain avian-flu pipeline**: https://github.com/nextstrain/avian-flu

## Data

Downloaded from GISAID EpiFlu (manual download, 16 March 2026). Raw files are in `data/` (not tracked by git).

| File tag | Description | Isolates |
|----------|-------------|----------|
| `a_human_africa` | Influenza A — human host — Africa (excl. H3) | 19,902 |
| `a_excepth3_human_africa` + `a_h3_human_africa` | FASTA sequences for the above (split by subtype) | 161,143 seqs |
| `b_human_africa` | Influenza B — human host — Africa | 5,081 |
| `c_human_africa` | Influenza C — human host — Africa | 10 |
| `animal_environmental_africa` | All animal/environmental hosts — Africa | 5,160 |

**Combined (deduplicated): 30,065 isolates — 1,849 from Madagascar**

### Genome completeness (Madagascar)

| Level | Count | % |
|-------|-------|---|
| Total isolates | 1,849 | 100% |
| Complete metadata (8 segment IDs in GISAID) | 782 | 42.3% |
| Phylo-ready (8 segments + sufficient length) | 782 | 42.3% |

Isolates that are not phylo-ready are listed in `data/madagascar_incomplete_genomes.tsv`
with columns indicating which segments are missing from metadata, absent from FASTA, or too short.

## Repository structure

```
├── data/                        # GISAID EpiFlu downloads (not tracked by git)
│   └── madagascar_incomplete_genomes.tsv   # Isolates needing completion
├── scripts/
│   ├── analyse_gisaid.py        # Main pipeline: load XLS + FASTA → combined_metadata.tsv
│   └── incomplete_genomes.py   # Export list of incomplete Madagascar isolates
├── report/
│   └── report.qmd              # Quarto report source
├── docs/                        # Rendered HTML report (GitLab Pages)
├── environment.yaml             # Conda environment (Python 3.11, plotly, quarto)
├── Makefile                     # make analyse | make report
└── README.md
```

## Reproducing the analysis

```bash
conda activate madagascar_influenza
make analyse   # → data/combined_metadata.tsv
make report    # → docs/index.html
```

## Progress

- [x] GISAID access obtained
- [x] All Africa datasets downloaded (30,065 isolates, all hosts)
- [x] Conda environment created (`madagascar_influenza`)
- [x] Metadata analysis: geographic distribution, sampling timeline, genome completeness
- [x] Quarto report with Plotly interactive figures (`docs/index.html`)
- [x] Incomplete Madagascar genomes identified and exported for collaborators
- [ ] Phylo-ready threshold validation against Nextstrain avian-flu pipeline
- [ ] GitLab Pages deployment
