---
title: "Comparing taxonomic profilers: concepts"
teaching: 35
exercises: 10
---

:::::::::::::::::::::::::::::::::::::: questions

- Why do different taxonomic profilers give different answers for the same reads?
- What does "abundance" actually mean in each profiler's output?
- How to choose the right profiler and database for the research question?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Describe the two main reference taxonomies used in metagenomics: NCBI RefSeq and GTDB
- Distinguish between sequence abundance, relative taxonomic abundance, and absolute abundance
- Explain Average Nucleotide Identity (ANI) and containment ANI
- Describe the database and approach of MetaPhlAn, mOTUs, Kraken2/Bracken and sylph
- Weigh the pros and cons of each approach for a given scenario
- Choose an appropriate output for a research question, and use the other outputs as complementary information

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

A taxonomic profiler answers two questions: *who is there?* and *how much of each?* Each profiler's answer depends on three choices:

1. the **reference database** used for profiling;
2. the **method** used to assign reads, or parts of reads, to taxa;
3. the **type of abundance** being reported.

Two profilers that differ in any of these can give different results for the same sample, and none is wrong! Before we compare the outputs in R (next episode), we need to understand these differences.

## Reference databases: the two pillars

Most metagenomic databases are built on one of two taxonomies:

| | NCBI Taxonomy / RefSeq | GTDB |
|---|---|---|
| What it is | NCBI's curated collection of reference sequences, labelled with the NCBI Taxonomy | The Genome Taxonomy Database: a taxonomy of bacteria and archaea built from genome phylogeny |
| How taxa are defined | Historical, literature-based classification. Ranks are not standardised across lineages | Phylogeny of conserved marker proteins, with ranks normalised by relative evolutionary divergence. Species are clusters of genomes at about 95% ANI |
| Coverage | Bacteria, archaea, viruses, fungi, other eukaryotes, plasmids, the human genome | Bacteria and archaea only, including many uncultured genomes (MAGs) |
| Names | Familiar, e.g. *Shigella flexneri*, phylum Firmicutes / Bacillota | Can differ from NCBI. For example, *Shigella* species fall within *Escherichia coli*, and genera are split with suffixes such as `Clostridium_A` |
| Updates | Continuously; each tool uses the snapshot from when its database was built | Numbered releases, about once a year. The sylph database used in this workshop is **GTDB R220** (April 2024) |

In the Sun *et al.* data, *Clostridium perfringens* genome is called *Sarcina perfringens* by GTDB, and *[Ruminococcus] gnavus* is *Ruminococcus_B gnavus*. This disparity even exits at a higher taxonomic level. For `Gut_sample1`, MetaPhlAn (NCBI) puts 66.9% of the community in the phylum Firmicutes, while sylph (GTDB R220) splits the same species across four phyla: *Bacillota* (32.7%), *Bacillota_A* (26.9%), *Bacillota_C* (5.8%) and *Bacillota_I* (2.5%).

::::::::::::::::::::::::::::::::::::: callout

### Names are not interchangeable

A profiler using a GTDB database and one using an NCBI database can detect the same organism under **different names**. Before comparing the profilers' outputs, we have to map them to a common taxonomy. Otherwise a correct detection looks like a false positive plus a false negative. Even within NCBI, names change between versions while the **taxonomy ID** stays the same. In our data, taxonomy ID 33038 is *[Ruminococcus] gnavus* for mOTUs, *Ruminococcus gnavus* for MetaPhlAn and *Mediterraneibacter gnavus* for Kraken2. Therefore, taxonomy matching should be done using taxonomy IDs, not names. 

::::::::::::::::::::::::::::::::::::::::::::::::

## What does "abundance" mean?

Sun *et al.* (2021) showed that much of the disagreement between profilers comes from reporting different kinds of abundance:

- **Sequence abundance**: the fraction of *sequenced reads (or bases)* that come from a taxon.
- **Taxonomic abundance**: the fraction of *cells (genomes)* that belong to a taxon.
- **Absolute abundance**: the actual amount of a taxon, e.g. cells per gram of stool. None of the profilers we use can estimate this from sequencing alone. It requires extra measurements such as spike-ins, qPCR or flow cytometry.

Sequence and taxonomic abundance are related through the genome size. A taxon with a genome twice as large contributes roughly twice as many reads per cell:

$$\text{sequence abundance}_i \propto \text{taxonomic abundance}_i \times \text{genome size}_i$$
Polyploidy (a cell or an organism possessing three or more complete sets of chromosomes) further complicates this relationship.

::::::::::::::::::::::::::::::::::::: challenge

### Challenge 1: Same cells, different numbers

A community has equal numbers of cells of species A (genome size 2 Mb) and species B (genome size 6 Mb). Assuming every read is classified, what is each species' taxonomic abundance and sequence abundance?

:::::::::::::::::::::::: solution

Taxonomic abundance: A = 50%, B = 50%.

Sequence abundance: B contributes three times as many bases per cell as A, so A = 2 / (2 + 6) = 25% and B = 75%.

The same community gives very different numbers depending on what is reported.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

## ANI and containment ANI

**Average nucleotide identity (ANI)** measures how similar two genomes are. Genomes with roughly ≥95% ANI are usually considered the same species. GTDB uses this criteria to (re-)classify taxonomy.  

**Containment ANI** estimates the ANI between a *reference genome* and the *reads in a metagenome*. It is based on the fraction of the reference genome's k-mers found (contained) in the sample. If the reference is present, nearly all its k-mers should be found. A lower containment ANI suggests a related strain or species, instead of a reference match. At low sequencing depth, many k-mers are missed simply by chance. `Sylph` corrects for this and reports an **adjusted ANI**.

In our data, most genomes detected by `Sylph` have an adjusted ANI above 99.5%. A high adjusted contaiment ANI means that either the reads came from that genome or they are almost identical to it. 

- *Gemella morbillorum* (99.98% ANI): The sample’s genome is practically identical to the reference genome. The reads were simulated from the same strain that is in the database.

## The four profilers

| | MetaPhlAn | mOTUs | Kraken2 + Bracken | sylph |
|---|---|---|---|---|
| Approach | Maps reads to **clade-specific marker genes** | Maps reads to **universal single-copy marker genes** (10 genes) | **Exact k-mer matching**, each read assigned to the lowest common ancestor (LCA) of matching taxa. Bracken then re-estimates abundance at a chosen rank | **k-mer sketching**. Estimates containment ANI and coverage for each reference genome |
| Database | ChocoPhlAn marker database, organised into species-level genome bins (SGBs) | mOTUs database: marker-gene clusters from reference genomes *and* metagenome-assembled genomes | User-chosen, often built from RefSeq (e.g. the Kraken2 *Standard* database). Bracken files are built per database and read length | Pre-built genome databases, e.g. GTDB |
| Reports | Taxonomic abundance (relative) | Relative abundance from single-copy markers (≈ taxonomic). Read counts optional | Reads per taxon (sequence abundance) | **Both** sequence and taxonomic abundance, plus ANI and coverage |
| Unknown taxa | Unknown SGBs (uSGBs) from MAGs | Reports an *unassigned* (`-1`) fraction for species not in the database | Reads stay unclassified, or are assigned to a higher rank | Genomes below the ANI threshold are not reported |
| **In this workshop** | MetaPhlAn 4.1.1, <code>mpa_<wbr>vJun23_<wbr>CHOCOPhlAnSGB_<wbr>202403</code> | mOTUs 3.1.0, run with read counts (`-c`) | Kraken2 2.1.5 + Bracken 3.0.1, Standard-8 database <code>k2_<wbr>standard_<wbr>08_<wbr>GB_<wbr>20260626</code> (capped at 8 GB) | sylph 0.7.0 + sylph-tax 1.2.0, GTDB R220 (`c200`) |
| **Values in our tables** | Taxonomic abundance (%) | Marker-gene read counts (<code>insert.<wbr>scaled_<wbr>counts</code>) | Reads per species (<code>new_<wbr>est_<wbr>reads</code>) | Taxonomic and sequence abundance (%) |

All four ran inside nf-core/taxprofiler 2.0.1. Taxpasta 0.7.0 standardised the MetaPhlAn, mOTUs and Kraken2/Bracken outputs. Taxpasta does not support sylph, so the pipeline merges sylph's output with sylph-tax instead.

### Pros and cons

- **MetaPhlAn:** fast, with few false positives thanks to the unique marker genes. Covers bacteria, archaea and some eukaryotes and viruses. Can miss taxa that have few markers or are present at low coverage. Gives no read counts by default.
- **mOTUs:** estimates the proportion of *unknown* species, and works well for poorly characterised environments thanks to its MAG-derived markers. Only covers prokaryotes.
- **Kraken2/Bracken:** classifies every read, can include viruses and eukaryotes, and its read counts suit count-based statistics. Prone to false positives, especially with small or incomplete databases, so it often needs abundance or read-count filtering. Its abundance is sequence abundance. With the capped Standard-8 database used here, about a quarter of the reads stay unclassified.
- **sylph:** very fast, accurate at species level, reports both abundance types, and its ANI output helps tell strains apart. There is an option to use NCBI RefSeq fungal reference database, however, whether or not the algorithm works well with eukaryotic genomes is still unclear.

::::::::::::::::::::::::::::::::::::: challenge

### Challenge 2: Pick a profiler

Which profiler (or output) would you choose for each scenario, and why?

1. You want to know what fraction of the community belongs to species that are missing from reference databases.
2. You plan to use a count-based differential abundance method that needs read counts.
3. Your samples may contain DNA viruses and fungi as well as bacteria.
4. You want to compare cell proportions with flow-cytometry estimates (give cells/gram).

:::::::::::::::::::::::: solution

1. **mOTUs.** Its `-1` fraction estimates the unassigned species. sylph with a MAG-rich database is also an option.
2. **Kraken2/Bracken** read counts, or mOTUs run with read counts enabled.
3. **Kraken2** with a database that includes viral and fungal genomes. MetaPhlAn also covers some eukaryotes and viruses.
4. **Taxonomic abundance:** sylph's `Taxonomic_abundance`, MetaPhlAn or mOTUs. Not Bracken's sequence abundance.

These are discussion points rather than single right answers.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

## Using outputs as complementary information

There is rarely one "best" profiler. A practical strategy:

1. **Choose a primary output** that matches your question: taxonomic abundance for questions about cell proportions, read counts for count-based statistical models.
2. **Use the other profilers as a sanity check.** Taxa detected by several tools, built on different databases and methods, are more trustworthy. Taxa detected by only one tool, at low abundance, deserve scrutiny.
3. **Report what you used:** tool, version, database and abundance type. Results from different profilers, or different abundance types, should not be mixed in one analysis.

::::::::::::::::::::::::::::::::::::: keypoints

- Profilers differ in their database (NCBI/RefSeq vs GTDB), method (marker genes, k-mers, sketching) and type of abundance, and all three affect the results.
- Sequence abundance (fraction of reads) and taxonomic abundance (fraction of cells) differ systematically because of variation in genome sizes and polyploidy.
- Containment ANI estimates how much of a reference genome is "contained" in the sample.
- Choose the profiler and abundance type that match your research question, and use agreement between tools as supporting evidence.

::::::::::::::::::::::::::::::::::::::::::::::::
