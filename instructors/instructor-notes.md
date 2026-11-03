---
title: 'Instructor Notes'
---

## Data used on Day 2

The files in `episodes/data/` were made from the organiser run of nf-core/taxprofiler 2.0.1 on the **full** five Sun *et al.* (2021) gut samples (about 10 million read pairs each), run on 27 October 2026. Day 1 participants run the same pipeline on a 100,000-pair subset (`Nectar_setup/make_participant_subset.sh`).

| Profiler | Version | Database |
|---|---|---|
| MetaPhlAn | 4.1.1 | `mpa_vJun23_CHOCOPhlAnSGB_202403` |
| mOTUs | 3.1.0 | `db_mOTU_v3.1.0`, read counts (`-c`) |
| Kraken2 / Bracken | 2.1.5 / 3.0.1 | `k2_standard_08_GB_20260626` (Standard-8, capped at 8 GB) |
| sylph / sylph-tax | 0.7.0 / 1.2.0 | GTDB R220, `c200` |
| Taxpasta | 0.7.0 | |

Pre-processing: fastp and FastQC only. Host removal was off (simulated reads).

### Why the lesson does not use the Taxpasta tables directly

In this run Taxpasta was called without a taxonomy, so `taxpasta/*.tsv` contain taxonomy IDs only, with all ranks mixed in one table. Taxpasta 0.7.0 also does not support sylph. The script in `data-prep/` reads each profiler's per-sample reports instead (they contain names) and writes species-level tables with the same layout.

To get Taxpasta tables with names straight from the pipeline, add these to the Day 1 parameters (they need an NCBI taxdump folder with `nodes.dmp`, `names.dmp` and `merged.dmp`):

```
--taxpasta_taxonomy_dir /path/to/ncbi_taxdump
--taxpasta_add_name
--taxpasta_add_rank
--taxpasta_add_lineage
```

### Regenerating the data

`--truth` is the folder with the Sun *et al.* tables `5_gut_taxonomic_abd.txt` and `5_gut_sequence_abd.txt`.

```bash
python3 data-prep/prepare_day2_data.py \
    --results  /path/to/taxprofiler_results \
    --truth    /path/to/truth \
    --cami     /path/to/sun2021-cami-profiles/taxonomic_abd \
    --gtdb-map data-prep/species_map_gtdb_r220.tsv \
    --out      episodes/data
```

The CAMI profiles are from <https://github.com/shenwei356/sun2021-cami-profiles>. They provide the NCBI taxonomy IDs of the truth species (NCBI taxonomy of 2021-12-06). If you regenerate the data with different pipeline results, check the numbers quoted in the text of episodes 02 and 03.

### GTDB names of the truth species

The truth species were first mapped to GTDB R220 by name (the GTDB species of each NCBI type strain): column `gtdb_r220_species_by_name` in `ncbi_to_gtdb_r220.tsv`. For 12 species, GTDB splits the NCBI species, and sylph reported the simulated genome in a sibling GTDB species (e.g. *Fusobacterium nucleatum* → *Fusobacterium animalis*). These 12 were reassigned when sylph's species had the same presence pattern across the five samples, abundances that differed by at most 0.15 percentage points in total, and ANI ≥ 97% (column `gtdb_r220_species`, `evidence` = "genome"). Because this uses sylph's output, mention it when discussing sylph's near-perfect scores (Episode 03, "Is sylph really that good?"). A fully independent check would extract the simulated genome accessions from the read names and look them up in the GTDB R220 metadata.

## Key results to expect in Episode 03

- Bracken reports ~500–660 species per sample (truth: 56–76). Precision 0.11 without a threshold, 0.93 at 0.1%.
- Mean Bray–Curtis to truth: Bracken 0.09 (sequence) vs 0.21 (taxonomic); MetaPhlAn 0.17 vs 0.05; mOTUs 0.16 vs 0.03; sylph 0.02 / 0.01 against its matching truth.
- Missed species: MetaPhlAn misses *Rothia dentocariosa* (reported as the unnamed `Rothia SGB49305`); Bracken misses *Faecalitalea cylindroides* (no *Faecalitalea* in the Standard-8 database).
- PERMANOVA on Bray–Curtis: sample R² ≈ 0.95, method R² ≈ 0.045 (p = 0.001); `betadisper` not significant.
- Shannon index: taxonomic and sequence truth are similar here (paired Wilcoxon p ≈ 0.81), unlike in Sun *et al.*. With n = 5 pairs, the smallest possible exact p-value is 0.0625.

## Common problems

- `left_join()` error *"Input includes an object that is not 'SummarizedExperiment'"*: `mia` was loaded after `dplyr`. Load `mia` first, or use `dplyr::left_join()`.
- File-not-found errors: the RStudio project must be the folder that **contains** `data/`. Episode 3 reads from `/data/workshop_data/day2_ep3` (the workshop RStudio server) by default; remote participants set `data_dir <- "data"`. When the website is built, a hidden chunk falls back to `episodes/data/`.
