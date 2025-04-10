# ThreeGroupQuadDiff
# Main Executable Lines

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
