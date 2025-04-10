# ThreeGroupQuadDiff
## Main Executable Lines

**Note:** `markers_df_hhc` & `markers_df_hth` are data frames from the `Seurat::FindMarkers` function.

### Step 1: Select the common features of the two datasets.
```r
markers_df_common <- common_de_marker_selection(
  de_list_1 = markers_df_hhc,
  de_list_2 = markers_df_hth,
  suffix_list = c("_hhc", "_hth"),
  xtitle = "LogFC (HIV vs HC)",
  ytitle = "LogFC (HIV+TB vs HIV)"
)
```
### Step 2: Convert the combined list.
```r
converted_common_markers <- common_de_geneid_converter(
  common_de_list = markers_df_common$common_markers
)
```
### Step 3: Enrich the common list using databases.
```r
common_enrich <- common_markers_enrichment(converted_common_markers)
```
**Hint:** You might consider integrating additional databases at this stage, such as `REACTOME` and `WikiPathways`, among others.

### Step 4: Select the common features of the two datasets.
```r
c1 <- common_plot_kegg_results_v5(
  cd4naive_common_enrich,
  hjust = 1.7,
  title = "KEGG for CD4 Naive in four quadrants"
)
```
### Suggest Citations
```r
print("Please consider cite the package 'clusterProfiler' used in this repository.")
```
