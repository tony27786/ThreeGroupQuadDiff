# Load the functions
source("define_functions.R")

# Step 1: Select the common features of the two datasets.
markers_df_common <- common_de_marker_selection(
  de_list_1 = markers_df_hhc,
  de_list_2 = markers_df_hth,
  suffix_list = c("_hhc", "_hth"),
  xtitle = "LogFC (HIV vs HC)",
  ytitle = "LogFC (HIV_MTB vs HIV)",
  plot_title = "Common Gene LogFC Distribution among Three Groups"
)
# Step 2: Convert the combined list.
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
converted_common_markers <- common_de_geneid_converter(
  common_de_list = markers_df_common$common_markers, mart = mart
)
# Step 3: Enrich the common list using databases.
# You might consider integrating additional databases at this stage, such as REACTOME and WikiPathways, among others.
common_enrich <- common_markers_enrichment(converted_common_markers)
# Step 4: Plot the KEGG pathway enrichment results.
c1 <- common_plot_kegg_results_v5(cd4naive_common_enrich, hjust = 1.7, title = "KEGG for CD4 Naive in four quadrants")

### Suggest Citations
print("Please consider cite the packages 'clusterProfiler' used in this repository.")
print("Please consider cite the article of this repository: https://doi.org/10.3389/fimmu.2025.1680538")
