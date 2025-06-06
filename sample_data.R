# Read the sample data and functions
markers_df_hhc <- readRDS("data/markers_df_hhc.rds")
markers_df_hth <- readRDS("data/markers_df_hth.rds")
source("define_functions.R")

# Main executable lines
markers_df_common <- common_de_marker_selection(
  de_list_1 = markers_df_hhc,
  de_list_2 = markers_df_hth,
  suffix_list = c("_hhc", "_hth"),
  xtitle = "LogFC (HIV vs HC)",
  ytitle = "LogFC (HIV_MTB vs HIV)",
  plot_title = "Common Gene LogFC Distribution among Three Groups"
)
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
converted_common_markers <- common_de_geneid_converter(
  common_de_list = markers_df_common$common_markers, mart = mart
)
common_enrich <- common_markers_enrichment(converted_common_markers)
p1 <- common_plot_kegg_results(
  common_enrich,
  hjust = 1.7,
  title = "KEGG for CD4 Naive in four quadrants"
)
