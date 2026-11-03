---
title: 'Episode 3: background for instructors'
---

This page explains the reasoning behind each statement in [Community analysis](../episodes/03-community-analysis.Rmd), and gives the numbers behind it. All the numbers were recomputed from the files in `episodes/data/` with the episode's own code (the `mia`/`scater` parts from `episodes/data/ep03_precomputed.RData`). Statements attributed to Sun *et al.* (2021) were checked against the bioRxiv preprint ([10.1101/2020.11.14.382994](https://www.biorxiv.org/content/10.1101/2020.11.14.382994)).

## Before teaching: points to fix or soften in the episode

| Where | Statement | What the data show | Suggestion |
|---|---|---|---|
| Distances, after the table | "With Bray–Curtis, L1, L2 and rJSD, the 'winner' depends on the truth we choose" | True for Bray–Curtis, L1 and L2. **Not for rJSD**: against the sequence truth, Bracken (0.173) is *worse* than MetaPhlAn (0.160) and mOTUs (0.152) | "With Bray–Curtis, L1 and L2 …", and mention that rJSD still moves each tool towards its matching truth |
| Distances, first bullet | MetaPhlAn and mOTUs ≈ 0.16 against the sequence truth | MetaPhlAn 0.169, mOTUs 0.157 | "≈ 0.16–0.17" (the instructor notes already say 0.17) |
| Challenge 2 instructor box | the 12 split species were reassigned with "ANI ≥ 97%" | `data-prep/prepare_day2_data.py` checks the presence pattern and the abundance difference only. It does not read the ANI, although its `evidence` column says "sylph ANI" | Check the ANI in sylph's `.sylphmpa` genome (`t__`) rows, or remove "and ANI ≥ 97%" |
| Challenge 3 solution | removing false positives "restores the dependence on abundance type" | Shown by Sun *et al.*. Here it only partly replicates: clear for mOTUs, weak for MetaPhlAn, absent for Bracken (see below) | Keep it as Sun *et al.*'s finding; don't promise learners they will see it here |
| Challenge 4 solution | "mOTUs and sylph change little" | sylph 66.2 → 65.4, but mOTUs 72.4 → 67.0 | "sylph changes little, and mOTUs drops by about 5" |

## Load packages and data

**Load `mia` before `dplyr`.** `mia` exports its own `left_join()`, `inner_join()`, `full_join()` and `right_join()` for merging `SummarizedExperiment` objects (see its `NAMESPACE`). Whichever package is attached last wins, so with `dplyr` first, `left_join()` on two data frames calls `mia`'s version. It then stops in `mergeSEs()` with *"Input includes an object that is not 'SummarizedExperiment'"*.

**The data table** (versions and databases) comes from the organiser run of nf-core/taxprofiler 2.0.1 (27 October 2026); see `instructors/instructor-notes.md`. How each file was made from the pipeline output is in `data-prep/` and on the learner page [From pipeline output to analysis tables](../learners/pipeline-to-tables.md).

## A first look at the outputs

**"The truth sums to 1, MetaPhlAn and sylph to 100, Bracken and mOTUs report counts."** Column sums per sample: truth 1, MetaPhlAn 100, sylph 100, mOTUs 6,577–7,514, Bracken 7.39–7.88 million. Each tool keeps its native unit: the Sun *et al.* truth tables are fractions, MetaPhlAn and sylph report percentages, and Bracken and mOTUs report counts.

**"Bracken counts reads, mOTUs counts reads mapped to its 10 marker genes."** Bracken redistributes Kraken2's classified reads to species. With paired-end data Kraken2 counts *fragments* (read pairs): 9,991,772 QC pairs × (1 − 0.2347 unclassified) ≈ 7.65 million, which matches Bracken's 7.63 million for Gut_sample1. mOTUs only counts reads that align to 10 universal single-copy marker genes. Those genes are a tiny part of each genome, hence about 7,000 counts per sample.

**"The samples contain 56–76 species. MetaPhlAn and mOTUs report a handful more, sylph almost exactly the right number, and Bracken more than 500."**

| | Sample 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Truth | 67 | 76 | 56 | 67 | 64 |
| MetaPhlAn | 71 | 84 | 64 | 80 | 75 |
| mOTUs | 71 | 89 | 61 | 74 | 72 |
| sylph | 68 | 76 | 56 | 67 | 64 |
| Bracken | 534 | 571 | 508 | 662 | 576 |

MetaPhlAn and mOTUs report 4–13 extra species. sylph's one extra species is *Bifidobacterium* sp959020155 in sample 1, the only false positive it reports.

### Challenge 1

**"About 21–26% of the reads were not classified."** `kraken2_unclassified_pct` is 21.06–25.88%.

**Why, if the genomes are in the databases?** Sun *et al.* picked the simulated genomes from the overlap of the profilers' reference databases, to avoid database bias. The species are therefore there, at least in the database versions of 2020. Kraken2 classifies a read only if its k-mers (minimizers) match the database. The **Standard-8** index is the standard database *capped at 8 GB*, so many k-mers are left out to make it fit, and reads whose k-mers were all dropped stay unclassified. A full-size database would classify more of them.

**"Hundreds of species with tens to a few thousand reads each."** In Gut_sample1, Bracken reports 534 species. Species ranked 60–75 have 1,604–11,446 reads, and the last ones have exactly 10 reads. That is Bracken's default minimum, so the tail is cut at 10. On average, 3.3% of Bracken's relative abundance sits on species that are not in the truth (MetaPhlAn 0.7%, mOTUs 1.8%, sylph 0.15%).

**Why false positives?** A read is assigned to the lowest taxon whose genome contains its k-mers. Short k-mers shared between related genomes, plus Bracken's redistribution of genus-level reads to species, spread reads thinly onto relatives of the true species.

## Same species, different names

**Taxonomy ID 33038** appears as *[Ruminococcus] gnavus* in the truth and mOTUs, *Ruminococcus gnavus* in MetaPhlAn, and *Mediterraneibacter gnavus* in Bracken. Each tool ships with a snapshot of the NCBI taxonomy from when its database was built: mOTUs and the truth with a 2021 snapshot, Kraken2's database with one from 2026, after NCBI moved the species to *Mediterraneibacter*. The square brackets are NCBI's old way of marking a genus name that is known to be wrong. **NCBI keeps the taxonomy ID when a name changes**, which is why the NCBI-based tools are joined by ID.

**sylph uses GTDB R220**, where the species is *Ruminococcus_B gnavus* in phylum *Bacillota_A*. GTDB defines taxa from genome phylogeny and species as ≥ 95% ANI clusters. It splits polyphyletic NCBI genera with suffixes (`_B`) and splits phyla (Firmicutes → Bacillota, Bacillota_A, …). GTDB has no taxonomy IDs, so sylph's names have to be mapped.

**"All but one of sylph's species map to a truth species."** After `add_ncbi_id()`, the only sylph species without an NCBI ID is *Bifidobacterium* sp959020155, at 0.77% in sample 1 and absent from the other samples.

### Challenge 2

**"12 species"** differ between the two GTDB columns: *Clostridium botulinum*, *Desulfovibrio piger*, *Faecalibacterium prausnitzii*, *Finegoldia magna*, *Fusobacterium nucleatum*, *F. periodonticum*, *Gardnerella vaginalis*, *Isoptericola variabilis*, *Oxalobacter formigenes*, *Peptoniphilus harei*, *Rothia mucilaginosa* and *Streptococcus cristatus*. For each, GTDB splits the NCBI species into several ANI clusters. The type strain falls in one of them, which gives the name-based mapping, while the genome that was simulated falls in another (e.g. *F. nucleatum* → *Fusobacterium animalis*).

**"13 species"** fail to map by name: the 12 above plus the real false positive.

**Answer 3** follows from how precision and recall are counted. The detected GTDB name has no NCBI ID, so it counts as a false positive. The truth species then has no detection, so it also counts as a false negative.

**Instructor box: how the 12 were checked.** The summed difference between the truth and sylph's abundance over the five samples is 0.016–0.134 percentage points for the 12 species, so "at most 0.15" holds. The presence pattern matches for all 12, which the script checks. The script's own cut-off is < 0.5 percentage points. The ANI claim is not checked by the script (see the table at the top).

## Combine everything in one `TreeSummarizedExperiment`

**mOTUs' `-1` row** is 0.13–0.32% of the counts per sample ("about 0.1–0.3%").

**The object** has 1,101 taxa and 35 columns: 5 samples × 7 "methods" (two truths and five profiler outputs). Most of the 1,101 rows are Bracken's false positives. `table(tse$method, tse$abundance_type)` shows 5 columns per method, each with a single abundance type.

## Sequence abundance versus taxonomic abundance

**"If all genomes were the same size, the two would be identical."** Reads come from DNA, so a species' share of the reads is its share of the cells weighted by its genome size *G*:

sequence~i~ = taxonomic~i~ × *G*~i~ / Σ~j~ (taxonomic~j~ × *G*~j~)

With equal *G*, the weights cancel.

**"Within a sample, the ratio of the two abundances is proportional to genome size."** From the formula, sequence~i~ / taxonomic~i~ = *G*~i~ / *Ḡ*, where *Ḡ* is the cell-weighted mean genome size of that sample. *Ḡ* is the same for every species in a sample, so the ratio is proportional to *G*~i~. It changes slightly between samples, which is why the episode takes the median over samples.

**The extremes.** *Mycoplasma hominis* has a median ratio of 0.24 (0.24–0.26 across samples), "about a quarter". It has a very small genome (about 0.7 Mb). *Hungatella hathewayi* has a ratio of 2.74 (2.6–2.8) and a genome of about 7 Mb. The two agree: 0.7 / 0.24 ≈ 2.9 Mb and 7 / 2.74 ≈ 2.6 Mb for the average genome.

## How close is each profiler to the truth?

**The five measures.**

- **Bray–Curtis** = Σ|x − y| / Σ(x + y). Large absolute differences, which occur in abundant taxa, dominate it.
- **L1** = Σ|x − y|. When both profiles sum to 1, Bray–Curtis = L1 / 2 exactly (largest difference here 5 × 10⁻¹⁶), so the two always give the **same ranking**. L1 is in the list because Sun *et al.* used it.
- **L2** squares each difference before adding, so a few large errors dominate.
- **rJSD** is the square root of the Jensen–Shannon divergence: the average Kullback–Leibler divergence of each profile from their mean. It is symmetric and bounded (≤ √ln 2 ≈ 0.83 with the natural log used here).
- **rAD** applies a *robust* centred log-ratio transform to each profile: the log of each non-zero value minus the mean log of the non-zero values, with zeros left at 0. It then takes the Euclidean distance. On a log scale, a taxon at 0.001% weighs as much as one at 10%.

**Mean distance to each truth (5 samples):**

| Measure | Method | vs sequence truth | vs taxonomic truth |
|---|---|---|---|
| Bray–Curtis | Bracken | **0.092** | 0.209 |
| | MetaPhlAn | 0.169 | **0.050** |
| | mOTUs | 0.157 | **0.034** |
| | sylph (sequence) | **0.017** | 0.153 |
| | sylph (taxonomic) | 0.155 | **0.011** |
| L2 | Bracken / MetaPhlAn / mOTUs | 0.048 / 0.079 / 0.072 | 0.094 / 0.026 / 0.018 |
| rJSD | Bracken / MetaPhlAn / mOTUs | 0.173 / 0.160 / 0.152 | 0.227 / 0.098 / 0.076 |
| rAD | Bracken / MetaPhlAn / mOTUs | 62.7 / 17.3 / 7.8 | 62.9 / 17.1 / 7.2 |

**Why the ranking flips with the truth.** Bracken estimates the share of *reads*, which is sequence abundance. MetaPhlAn and mOTUs count hits to marker genes that occur once per genome, so they estimate the share of *genomes*, which is taxonomic abundance. Sun *et al.* make the same split: DNA-to-DNA tools (Kraken2, Bracken) give sequence abundance, and DNA-to-marker tools (MetaPhlAn, mOTUs) give taxonomic abundance. sylph reports both, and each output is closest to its matching truth.

### Challenge 3

**rAD** puts Bracken, MetaPhlAn and mOTUs about the same distance from both truths (62.7/62.9, 17.3/17.1, 7.8/7.2), and Bracken far behind. The order follows the number of false positives per sample: Bracken 505, MetaPhlAn 9.8, mOTUs 6.4, sylph 0.2. sylph, which has almost no false positives, still shows the abundance-type effect under rAD (0.30 vs 2.99), which supports the explanation.

**Sun *et al.*** found that rAD favoured MetaPhlAn2 and mOTUs2 whatever the truth, and that after removing the false positives rAD gave the same benchmark result as the other measures. Repeating that here (keeping only truth species):

| Method | vs sequence truth | vs taxonomic truth |
|---|---|---|
| Bracken | 25.6 | 25.8 |
| MetaPhlAn | 7.3 | 6.8 |
| mOTUs | 3.0 | 0.8 |

The dependence returns clearly for mOTUs, weakly for MetaPhlAn, and not for Bracken. Bracken's remaining distance most likely comes from errors in its low-abundance true species, which rAD's log scale also weights heavily.

## Precision and recall

**"Both truth tables contain the same species."** True: the same 95 taxonomy IDs, and the same species present in each sample. Detection therefore doesn't depend on which truth is used.

**Mean over the 5 samples:**

| Method | Threshold | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| Bracken | 0 | 65.2 | 505.0 | 0.8 | 0.11 | 0.99 | 0.21 |
| Bracken | 0.01% | 63.6 | 28.2 | 2.4 | 0.69 | 0.96 | 0.81 |
| Bracken | 0.1% | 60.8 | 4.8 | 5.2 | 0.93 | 0.92 | 0.92 |
| MetaPhlAn | 0 | 65.0 | 9.8 | 1.0 | 0.87 | 0.99 | 0.92 |
| mOTUs | 0 | 66.0 | 6.4 | 0.0 | 0.91 | 1.00 | 0.96 |
| sylph (either) | 0 | 66.0 | 0.2 | 0.0 | 1.00 | 1.00 | 1.00 |

"Only about 1 in 9 of the species it reports is real" is precision 0.11. Bracken's missed species is *Faecalitalea cylindroides*, in the 4 samples that contain it (no *Faecalitalea* in Standard-8). MetaPhlAn misses *Rothia dentocariosa* in 4 samples, which it reports as the unnamed SGB *Rothia* SGB49305 (instructor notes), and *Lachnospira eligens* in sample 5.

**Why a threshold helps Bracken most.** Its false positives are numerous but each is tiny. 505 of them fall to 28 at 0.01% and to 5 at 0.1%. The cost is the true species below the threshold.

**"Is sylph really that good?"** The reasons in the discussion box are all real:

- The simulated genomes, or genomes within GTDB's 95% ANI species clusters, are in sylph's GTDB database.
- The reads are simulated, with no host DNA or contamination.
- The GTDB names of 12 truth species were assigned using sylph's own output.

## Alpha diversity

**The three indices** are Hill numbers of order 0, 1 and 2. Observed richness counts every taxon equally. exp(Shannon) weights taxa by abundance. Inverse Simpson weights the abundant ones even more. The higher the order, the less a rare false positive matters.

**Means by method:**

| Method | Observed | Shannon | Inverse Simpson |
|---|---|---|---|
| Truth (taxonomic) | 66 | 3.66 | 24.4 |
| Truth (sequence) | 66 | 3.66 | 26.9 |
| MetaPhlAn | 74.8 | 3.63 | 23.2 |
| mOTUs | 72.4 | 3.70 | 25.4 |
| Bracken | 570 | 3.74 | 28.0 |
| sylph (either) | 66.2 | 3.68 | 24.8 / 27.3 |

**"Richness is very sensitive; Shannon changes much less."** Bracken has 570 species against 66, but a Shannon of 3.74 against 3.66. Its ~500 extra species carry only about 3% of the abundance, so they add little to Shannon.

**"Sun *et al.* found that the Shannon index is systematically higher for taxonomic abundance."** The preprint reports Shannon and Simpson significantly higher for taxonomic abundance (p < 0.001), across many simulated samples and habitats.

**"Here the two abundance types give similar Shannon values."** Per sample (taxonomic vs sequence): 3.49/3.48, 3.95/4.00, 3.49/3.35, 3.72/3.81, 3.67/3.69. Taxonomic is higher in only 2 of 5 samples. Paired Wilcoxon gives V = 6, p = 0.81.

**"The smallest p-value with 5 pairs is 0.0625."** Under the null hypothesis, each of the 5 differences is equally likely to be positive or negative, so there are 2⁵ = 32 equally likely sign patterns. The most extreme result (all the same sign) has probability 1/32 in each direction, so the two-sided p-value is at least 2/32 = 0.0625. With 5 pairs the test can never reach 0.05.

### Challenge 4

Observed richness before → after the 0.1% threshold: Bracken 570 → 65.6, MetaPhlAn 74.8 → 63.6, mOTUs 72.4 → 67.0, sylph 66.2 → 65.4, truth 66 → 65.0/65.4. The truth itself loses species because 5 truth species-in-sample entries are at or below 0.1% (2, 1, 0, 1 and 1 per sample). That is why "a filtered richness is no longer the richness of the community".

## Beta diversity

**"All 35 columns"** are 5 samples × 7 methods.

**"Five tight groups, one per sample"** in the PCoA. PERMANOVA (marginal R², 999 permutations): sample R² = 0.954 (F = 4,890), method R² = 0.045 (F = 152), both p = 0.001. That is the smallest possible p-value with 999 permutations (1/1000). The design is balanced (every sample × method once), so marginal and sequential R² agree.

**Something to point out:** the `method` factor includes the two truths, which differ only in abundance type. Part of the "4.5% due to method" is therefore the sequence-versus-taxonomic difference, not profiler error.

**`betadisper`:** F = 0.01, p = 1. PERMANOVA can be significant either because group centroids differ or because groups differ in spread. Equal spread means the method effect is a shift.

## Keypoints

Each keypoint summarises one section above: units and mapping (first look, names), joining by ID (Challenge 2), abundance type (sequence vs taxonomic, distances), k-mer false positives and thresholds (Challenge 1, precision and recall), robustness of alpha diversity indices (alpha diversity), and the profiler as a source of variation (PERMANOVA).
