# Prepare and screen candidate datasets for Day 2, episode 4 (differential abundance)
#
# What it does, for each candidate study in curatedMetagenomicData:
#   1. downloads the MetaPhlAn profile (as estimated counts) and the HUMAnN pathway
#      abundances for the stool samples, with the sample metadata;
#   2. saves them as one small .rds file per study in data-prep/day2_ep4/;
#   3. runs a quick ANCOM-BC2 screen (genus level, disease vs control, adjusted for
#      age, sex and BMI) and records run time, number of hits and the top genera;
#   4. lists the pathways that the top case-enriched genera contribute to.
#
# Everything it finds is written to data-prep/day2_ep4/screen_report.txt.
# Please send that file (and the sessionInfo at its end) back to Claude.
#
# How: from the top of the Day 2 lesson folder, run
#   Rscript data-prep/prepare_day2_ep4_data.R
# Downloads are about 1 GB in total the first time (cached by ExperimentHub
# afterwards). The ANCOM-BC2 runs take a few minutes in total.

# ---- packages ------------------------------------------------------------------

needed <- c("curatedMetagenomicData", "mia", "ANCOMBC", "dplyr")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Please install first: BiocManager::install(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(mia)
  library(ANCOMBC)
  library(dplyr)
})

# ---- settings ------------------------------------------------------------------

candidates <- list(
  QinN_2014       = "cirrhosis",
  ThomasAM_2019_c = "CRC",
  WirbelJ_2018    = "CRC"
)
out_dir <- file.path("data-prep", "day2_ep4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
report_file <- file.path(out_dir, "screen_report.txt")

# Write to the console and to the report at the same time
report <- file(report_file, open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  writeLines(txt, report)
}
say_print <- function(x) {
  txt <- capture.output(print(x))
  cat(txt, sep = "\n")
  writeLines(txt, report)
}

say("Day 2 episode 4 dataset screen, ", format(Sys.time()))
say("curatedMetagenomicData ", as.character(packageVersion("curatedMetagenomicData")),
    " | mia ", as.character(packageVersion("mia")),
    " | ANCOMBC ", as.character(packageVersion("ANCOMBC")))

# ANCOMBC renamed some arguments between versions; use whichever this version has
anc_args <- names(formals(ANCOMBC::ancombc2))
assay_arg <- if ("assay.type" %in% anc_args) "assay.type" else "assay_name"

# ---- main loop -----------------------------------------------------------------

for (study in names(candidates)) {
  case <- candidates[[study]]
  say("\n==================== ", study, " (", case, " vs control) ====================")

  meta <- sampleMetadata |>
    filter(study_name == study, body_site == "stool",
           study_condition %in% c(case, "control")) |>
    filter(!is.na(age), !is.na(gender), !is.na(BMI)) |>
    distinct(subject_id, .keep_all = TRUE)       # one sample per subject
  say("Samples with complete age, sex and BMI: ", nrow(meta))
  say_print(table(meta$study_condition))

  # -- download -----------------------------------------------------------------
  tse_tax <- returnSamples(meta, "relative_abundance", counts = TRUE, rownames = "short")
  tse_path <- returnSamples(meta, "pathway_abundance", rownames = "long")

  # use predictable assay names, whatever this package version calls them
  assayNames(tse_tax)[1] <- "counts"
  assayNames(tse_path)[1] <- "pathway_abundance"

  # tidy the sample metadata we need
  colData(tse_tax)$group <- factor(ifelse(colData(tse_tax)$study_condition == case,
                                          case, "control"),
                                   levels = c("control", case))
  colData(tse_path)$group <- colData(tse_tax)$group[match(colnames(tse_path),
                                                          colnames(tse_tax))]

  # drop pathway rows that are zero in every sample to keep the file small
  tse_path <- tse_path[rowSums(assay(tse_path) > 0) > 0, ]

  say("Taxonomic TSE: ", nrow(tse_tax), " taxa x ", ncol(tse_tax), " samples; assays: ",
      paste(assayNames(tse_tax), collapse = ", "), "; ranks: ",
      paste(taxonomyRanks(tse_tax), collapse = ", "))
  say("Pathway TSE: ", nrow(tse_path), " rows (", sum(!grepl("|", rownames(tse_path), fixed = TRUE)),
      " unstratified) x ", ncol(tse_path), " samples; assays: ",
      paste(assayNames(tse_path), collapse = ", "))
  say("Pathway column sums (first 3 samples), to check the units:")
  say_print(head(colSums(assay(tse_path)), 3))
  say("First pathway rows:")
  say_print(head(rownames(tse_path), 4))

  rds_file <- file.path(out_dir, paste0(study, ".rds"))
  saveRDS(list(taxa = tse_tax, pathways = tse_path), rds_file, compress = "xz")
  say("Saved ", rds_file, " (", round(file.size(rds_file) / 1e6, 1), " MB)")

  # -- ANCOM-BC2 screen -----------------------------------------------------------
  tse_genus <- agglomerateByRank(tse_tax, rank = "genus")
  anc_call <- list(
    data         = tse_genus,
    tax_level    = NULL,
    fix_formula  = "group + age + gender + BMI",
    rand_formula = NULL,
    p_adj_method = "holm",
    pseudo_sens  = TRUE,
    prv_cut      = 0.10,
    lib_cut      = 0,
    s0_perc      = 0.05,
    group        = NULL,
    struc_zero   = FALSE,
    neg_lb       = FALSE,
    alpha        = 0.05,
    global       = FALSE,
    n_cl         = 1,
    verbose      = FALSE
  )
  anc_call[[assay_arg]] <- "counts"

  set.seed(123)
  t0 <- Sys.time()
  out <- tryCatch(do.call(ancombc2, anc_call), error = function(e) e)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  if (inherits(out, "error")) {
    say("ANCOM-BC2 FAILED after ", elapsed, " s: ", conditionMessage(out))
    next
  }
  say("ANCOM-BC2 run time: ", elapsed, " s")

  res <- out$res
  say("Result columns: ", paste(colnames(res), collapse = ", "))

  lfc  <- paste0("lfc_group", case)
  se   <- paste0("se_group", case)
  q    <- paste0("q_group", case)
  diff <- paste0("diff_group", case)
  pss  <- paste0("passed_ss_group", case)

  res_g <- res |>
    transmute(genus = taxon, lfc = .data[[lfc]], se = .data[[se]], q = .data[[q]],
              diff = .data[[diff]], passed_ss = .data[[pss]]) |>
    arrange(q)

  say("Genera tested: ", nrow(res_g),
      " | q < 0.05: ", sum(res_g$diff, na.rm = TRUE),
      " | q < 0.05 and passed sensitivity: ", sum(res_g$diff & res_g$passed_ss, na.rm = TRUE),
      " (up in ", case, ": ", sum(res_g$diff & res_g$passed_ss & res_g$lfc > 0, na.rm = TRUE),
      ", down: ", sum(res_g$diff & res_g$passed_ss & res_g$lfc < 0, na.rm = TRUE), ")")
  say("Top 20 genera by q-value:")
  say_print(as.data.frame(mutate(head(res_g, 20), across(c(lfc, se), \(x) round(x, 2)),
                                 q = signif(q, 2))))

  # -- which pathways do the top case-enriched genera contribute to? -------------
  top_up <- res_g |>
    filter(diff, passed_ss, lfc > 0) |>
    arrange(desc(lfc)) |>
    head(5) |>
    pull(genus)
  say("Top case-enriched genera (by LFC, among significant): ",
      paste(top_up, collapse = ", "))

  is_case <- colData(tse_path)$group == case
  for (g in top_up) {
    g_short <- sub("^g__", "", g)
    rows <- grepl(paste0("|g__", g_short, "."), rownames(tse_path), fixed = TRUE)
    if (!any(rows)) {
      say("  ", g_short, ": no stratified pathway rows")
      next
    }
    contrib <- rowMeans(assay(tse_path)[rows, is_case, drop = FALSE])
    contrib <- sort(tapply(contrib, sub("\\|.*$", "", names(contrib)), sum),
                    decreasing = TRUE)
    say("  ", g_short, ": contributes to ", length(contrib),
        " pathways; top 8 by mean contribution in ", case, ":")
    say_print(head(round(contrib, 2), 8))
  }
}

say("\n==================== sessionInfo ====================")
say_print(sessionInfo())
close(report)
cat("\nDone. Please send", report_file, "back.\n")
