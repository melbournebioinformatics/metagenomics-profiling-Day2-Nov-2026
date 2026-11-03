# Make the episode 4 data file from the curatedMetagenomicData download
#
# Input : data-prep/day2_ep4/QinN_2014.rds, made by data-prep/prepare_day2_ep4_data.R
# Output: episodes/data/qin2014_cirrhosis.rds, a list with two TreeSummarizedExperiments
#   $taxa     MetaPhlAn estimated counts (species x samples), taxonomy in rowData
#   $pathways HUMAnN pathway abundances: all unstratified pathways, plus the
#             per-taxon (stratified) rows of the pathways used in the episode
# Both share the same colData: group, age, gender, BMI, number_reads.
#
# How: from the top of the Day 2 lesson folder, run
#   Rscript data-prep/make_day2_ep4_episode_data.R
# Then copy episodes/data/qin2014_cirrhosis.rds to /data/workshop_data/day2_ep4/
# on the workshop RStudio Server.

suppressPackageStartupMessages(library(TreeSummarizedExperiment))

x <- readRDS(file.path("data-prep", "day2_ep4", "QinN_2014.rds"))
taxa <- x$taxa
pathways <- x$pathways

# Assay names differ between curatedMetagenomicData versions; use fixed ones
assayNames(taxa)[1] <- "counts"
assayNames(pathways)[1] <- "abundance"

# Plain species names as row names (curatedMetagenomicData adds "species:")
rownames(taxa) <- sub("^species:", "", rownames(taxa))

# Keep only the sample information used in the episode
keep_cols <- c("group", "age", "gender", "BMI", "number_reads")
col_df <- colData(taxa)[, keep_cols]
col_df$gender <- factor(col_df$gender)
colData(taxa) <- col_df
pathways <- pathways[, colnames(taxa)]
colData(pathways) <- col_df

# Pathways: all unstratified rows plus the stratified rows of the focus pathways
focus <- c("P125-PWY", "HEME-BIOSYNTHESIS-II", "PWY-6737")
stratified <- grepl("|", rownames(pathways), fixed = TRUE)
pathway_id <- sub(":.*$", "", rownames(pathways))
pathways <- pathways[!stratified | pathway_id %in% focus, ]

# Taxa: drop species that are absent from every sample
taxa <- taxa[rowSums(assay(taxa, "counts")) > 0, ]

dir.create(file.path("episodes", "data"), showWarnings = FALSE)
out <- file.path("episodes", "data", "qin2014_cirrhosis.rds")
saveRDS(list(taxa = taxa, pathways = pathways), out, compress = "xz")
message("Saved ", out, " (", round(file.size(out) / 1e6, 2), " MB): ",
        nrow(taxa), " species, ", nrow(pathways), " pathway rows, ",
        ncol(taxa), " samples")
