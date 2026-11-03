---
title: "Recap from Day 1"
teaching: 10
exercises: 0
---

:::::::::::::::::::::::::::::::::::::: questions

- What did we produce on Day 1, and where does Day 2 pick up?
- Which `nf-core/taxprofiler` and `nf-core/funcprofiler` outputs will we analyse today?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Review the main stages from FASTQ files to nf-core pipeline outputs
- Identify the files from Day 1 that are the inputs for today's analysis

::::::::::::::::::::::::::::::::::::::::::::::::

## From reads to abundance tables

On Day 1 we went from raw shotgun metagenomic reads to standardised abundance tables:

| Stage | What happened | Tool | Output |
|---|---|---|---|
| 1. Data | Downloaded paired-end FASTQ files and sample metadata | Figshare, `wget`, MD5 check | Raw reads, `metadata.tsv` |
| 2. Pipeline setup | Built sample and database sheets and a reproducible launch command | Nextflow, nf-core | `samplesheet.csv`, `databases.csv`, `params.json` |
| 3. Pre-processing | Read QC and adapter trimming. Host-read removal was switched off, because the reads are simulated and contain no human DNA | nf-core/taxprofiler (FastQC, fastp) | MultiQC report |
| 4. Taxonomic profiling | The same reads were classified by four profilers | MetaPhlAn, mOTUs, Kraken2/Bracken, sylph | One profile per tool and sample |
| 5. Standardisation | Profiles were converted to a common format and merged | Taxpasta (MetaPhlAn, mOTUs, Kraken2/Bracken), sylph-tax (sylph) | One merged table per profiler |
| 6. Import into R | Tables and metadata were loaded into a single object | `mia`, `TreeSummarizedExperiment` | `TreeSummarizedExperiment` object |

<!-- FIXME: replace with a workflow figure (e.g. episodes/fig/workflow.png) shared with the Day 1 lesson -->

## The dataset

We are using the five simulated human gut metagenomes from the benchmarking study by Sun *et al.* (2021), [*Challenges in benchmarking metagenomic profilers*](https://doi.org/10.1038/s41592-021-01141-3). The reads were simulated from known genomes in known proportions, so we know the **true** composition of each sample.

| | |
|---|---|
| Samples | 5 (`Gut_sample1` to `Gut_sample5`) |
| Species | 95 in total, 56–76 per sample |
| Reads | about 10 million read pairs per sample, Illumina 2 × 150 bp |
| Ground truth | taxonomic abundance (fraction of cells) and sequence abundance (fraction of reads) for every species |

To avoid favouring one profiler, Sun *et al.* only used species present in the reference databases of MetaPhlAn, mOTUs and Kraken2 at the time. The community is therefore a benchmark rather than a realistic gut: it includes oral, vaginal and skin species, and lacks some typical gut dominants.

On Day 1, you ran the pipeline on a random subset of 100,000 read pairs per sample, so that it finishes within the session. A profile from 1–2% of the reads misses many low-abundance species. Today we analyse the profiles from the **full** samples, produced with the same pipeline (nf-core/taxprofiler 2.0.1), the same profilers and the same databases.

## What we do today

1. **Concepts:** why the four profilers give different answers for the same reads.
2. **Practice:** compare the profilers against the ground truth, and see how the choice of profiler changes alpha and beta diversity.
3. **Differential abundance:** test which taxa differ between groups of samples.

::::::::::::::::::::::::::::::::::::: discussion

### Quick check

On Day 1, every sample was run through four profilers. Before we look at the results, do you expect them to agree? On what: the list of species present, their abundances, or both?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints

- nf-core/taxprofiler ran four profilers on the same reads, and Taxpasta and sylph-tax standardised their outputs into comparable tables.
- Our samples have a known composition (the ground truth), so we can judge how accurate each profiler is.
- Today's analysis uses the profiles of the full samples, not the subset from Day 1.

::::::::::::::::::::::::::::::::::::::::::::::::
