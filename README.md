# Madagascar Influenza

*Last updated: 2026-03-16*

## Goal

Characterize avian influenza sequences from Madagascar and Africa available in GISAID:

- Which avian influenza sequences exist from Madagascar and Africa in GISAID?
- What is their geographic location, sampling time, and whether the genomes are complete or not?

## Resources

- **GISAIDR** — R package for programmatic access to GISAID: https://github.com/Wytamma/GISAIDR
  - Note: GISAIDR does not support EpiFlu (avian influenza). Data was downloaded manually from the GISAID EpiFlu interface.
- **Nextstrain H5N1 HA (2y)**: https://nextstrain.org/avian-flu/h5n1/ha/2y

## Data

Downloaded from GISAID EpiFlu (manual download). Files are in `data/`:

| File | Description | Isolates | Sequences |
|------|-------------|----------|-----------|
| `gisaid_epiflu_isolates.xls` | Metadata — all hosts, Madagascar | 1,849 | 10,805 segments |
| `gisaid_epiflu_sequence.fasta` | Sequences — all hosts, Madagascar | — | 10,805 segments |
| `gisaid_epiflu_isolates_avian_africa.xls` | Metadata — **avian host, all Africa** | 4,923 | 21,455 segments |
| `gisaid_epiflu_sequence_avian_africa.fasta` | Sequences — **avian host, all Africa** | — | 21,455 segments |

> Files containing `avian_africa` in their name correspond to sequences where the **host is avian** and the **location is Africa** (all countries).

## Repository structure

```
├── data/                  # GISAID EpiFlu downloads (not tracked by git)
├── scripts/               # Analysis scripts
│   └── query_gisaid.R     # Attempted programmatic GISAID query (GISAIDR, EpiFlu unsupported)
├── docs/                  # Quarto-rendered HTML report (GitHub Pages)
├── environment.yaml       # Conda environment (Python 3.11 + plotly + quarto)
└── README.md
```

## Progress

- [x] GISAID access obtained
- [x] Madagascar sequences downloaded (1,849 isolates, 10,805 segments)
- [x] Avian Africa sequences downloaded (4,923 isolates, 21,455 segments)
- [x] Conda environment created (`madagascar_influenza`)
- [ ] Metadata analysis: geographic distribution, sampling timeline, genome completeness
- [ ] Quarto report with Plotly interactive figures published to GitHub Pages
