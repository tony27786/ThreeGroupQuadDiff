# Define the functions
common_de_marker_selection <- function(de_list_1, de_list_2, 
                                        suffix_list = c("_1", "_2"),
                                        plot_title = "Common Gene LogFC Distribution among Three Groups", 
                                        xtitle, 
                                        ytitle) {
  de_list_1 <- de_list_1[de_list_1$p_val_adj <= 0.05, ]
  de_list_2 <- de_list_2[de_list_2$p_val_adj <= 0.05, ]
  de_list_1$GeneSymbol <- rownames(de_list_1)
  de_list_2$GeneSymbol <- rownames(de_list_2)
  
  common_markers <- merge(de_list_1, de_list_2, by = "GeneSymbol", suffixes = suffix_list)
  
  avg_log2FC_col_1 <- paste0("avg_log2FC", suffix_list[1])
  avg_log2FC_col_2 <- paste0("avg_log2FC", suffix_list[2])
  p_val_cols <- paste0("p_val", suffix_list)
  
  common_markers <- common_markers %>%
    dplyr::select(-dplyr::contains(c("pct", p_val_cols))) %>%
    dplyr::mutate(
      class = dplyr::case_when(
        .[[avg_log2FC_col_1]] > 0 & .[[avg_log2FC_col_2]] > 0 ~ "1",
        .[[avg_log2FC_col_1]] < 0 & .[[avg_log2FC_col_2]] > 0 ~ "2",
        .[[avg_log2FC_col_1]] < 0 & .[[avg_log2FC_col_2]] < 0 ~ "3",
        .[[avg_log2FC_col_1]] > 0 & .[[avg_log2FC_col_2]] < 0 ~ "4"
      )
    )
  
  summary <- common_markers %>% dplyr::count(class) %>% print()
  
  library(ggplot2)
  p <- ggplot(common_markers, 
              aes(x = !!sym(avg_log2FC_col_1), 
                  y = !!sym(avg_log2FC_col_2), 
                  color = class)) +
    geom_point(alpha = 0.6, size = 3) +
    ggrepel::geom_text_repel(aes(label = GeneSymbol),
                             size = 3, force = 0.3,
                             max.overlaps = 22,
                             box.padding = 0.2,
                             point.padding = 0.2) +
    labs(title = plot_title,
         x = xtitle,
         y = ytitle) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    annotate("text", x = 2,  y = 2,  label = paste0("Class 1: ", sum(common_markers$class == "1"))) +
    annotate("text", x = -2, y = 2,  label = paste0("Class 2: ", sum(common_markers$class == "2"))) +
    annotate("text", x = -2, y = -2, label = paste0("Class 3: ", sum(common_markers$class == "3"))) +
    annotate("text", x = 2,  y = -2, label = paste0("Class 4: ", sum(common_markers$class == "4")))
  print(p)
  return(list(common_markers = common_markers, plot = p))
}

common_de_geneid_converter <- function(common_de_list) {
  # Load library
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(biomaRt)
  
  if (is.null(mart)) {
    stop("Please download biomaRt reference first! For example:\nmart <- useEnsembl('ensembl', dataset='hsapiens_gene_ensembl')")
  }
  
  common_gene <- common_de_list$GeneSymbol
  de_data <- data.frame(
    OriginalID = common_de_list$GeneSymbol,
    Case = common_de_list$class,
    stringsAsFactors = FALSE
  )
  
  gene_symbols <- de_data$OriginalID[!grepl("^ENSG", de_data$OriginalID)]
  ensembl_ids_with_version <- de_data$OriginalID[grepl("^ENSG", de_data$OriginalID)]
  ensembl_ids_no_version <- sub("\\..*$", "", ensembl_ids_with_version)
  
  combined_data <- data.frame(
    GeneSymbol = character(),
    EnsemblID  = character(),
    EntrezID   = character(),
    Case       = character(),
    stringsAsFactors = FALSE
  )
  
  # (1) Deal with gene symbol part
  if (length(gene_symbols) > 0) {
    gene_symbols_to_entrez <- clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
    
    gene_symbols_to_ensembl <- clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType   = "ENSEMBL",
      OrgDb    = org.Hs.eg.db
    )
    
    data_symbols <- de_data %>%
      dplyr::filter(OriginalID %in% gene_symbols) %>%
      dplyr::left_join(gene_symbols_to_entrez, by = c("OriginalID" = "SYMBOL")) %>%
      dplyr::left_join(gene_symbols_to_ensembl, by = c("OriginalID" = "SYMBOL")) %>%
      dplyr::rename(
        GeneSymbol = OriginalID,
        EntrezID   = ENTREZID,
        EnsemblID  = ENSEMBL
      ) %>%
      dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
      dplyr::select(GeneSymbol, EnsemblID, EntrezID, Case)
    
    combined_data <- dplyr::bind_rows(combined_data, data_symbols)
  } else {
    warning("Warning: Gene symbol dataset is empty")
  }
  
  # (2) Deal with EnsemblID part
  if (length(ensembl_ids_with_version) > 0) {
    ensembl_to_entrez <- getBM(
      attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = unique(ensembl_ids_no_version),
      mart = mart
    )
    
    # Ensure entrezgene_id is character
    ensembl_to_entrez$entrezgene_id <- as.character(ensembl_to_entrez$entrezgene_id)
    
    data_ensembl <- de_data %>%
      dplyr::mutate(EnsemblID_NoVer = sub("\\..*$", "", OriginalID)) %>%
      dplyr::filter(OriginalID %in% ensembl_ids_with_version) %>%
      dplyr::left_join(ensembl_to_entrez, by = c("EnsemblID_NoVer" = "ensembl_gene_id")) %>%
      dplyr::rename(
        EnsemblID  = OriginalID,
        EntrezID   = entrezgene_id,
        GeneSymbol = external_gene_name
      ) %>%
      dplyr::select(GeneSymbol, EnsemblID, EntrezID, Case)
    
    combined_data <- dplyr::bind_rows(combined_data, data_ensembl)
  } else {
    warning("Warning: Ensembl ids dataset is empty")
  }
  
  # Final clean up
  combined_data <- combined_data %>%
    dplyr::filter(!is.na(EntrezID) & EntrezID != "")
  combined_data <- combined_data %>%
    dplyr::filter(!is.na(GeneSymbol) & GeneSymbol != "")
  
  row.names(combined_data) <- NULL
  return(combined_data)
}

common_markers_enrichment <- function(common_entrezid_list) {
  # Load Library
  library(clusterProfiler)
  library(org.Hs.eg.db)
  # Check input data structure
  if (!"Case" %in% colnames(common_entrezid_list) || !"EntrezID" %in% colnames(common_entrezid_list)) {
    stop("The input data must contain 'Case' and 'EntrezID' columns.")
  }
  # Initial result list
  results <- list()
  
  # Loop for cases
  for (case_id in 1:4) {
    # extract case list
    case_genes <- common_entrezid_list %>%
      dplyr::filter(Case == case_id) %>%
      pull(EntrezID)
    # If some case with no gene, skip the case
    if (length(case_genes) == 0) {
      print(paste0("Case ", case_id, " is empty, skipping..."))
      results[[paste0("Case_", case_id)]] <- list(KEGG = NULL, GO = NULL)
      next
    }
    
    # KEGG
    kegg_result <- enrichKEGG(gene = case_genes, organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH")
    # GO
    go_result <- enrichGO(gene = case_genes, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH")
    
    # Save result
    results[[paste0("Class_", case_id)]] <- list(KEGG = kegg_result, GO = go_result)
    print(paste0("Class ", case_id, " KEGG and GO analysis completed."))
  }
  
  return(results)
}

common_plot_kegg_results <- function(results_list, 
                                        cases = c("Class_1", "Class_2", "Class_3", "Class_4"), hjust = 1.8, title = "KEGG Enriched Pathways in four quadrants") {
  library(clusterProfiler)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(RColorBrewer)
  
  # RColorBrewer for four colours
  brewer_cols <- RColorBrewer::brewer.pal(4, "Set1")
  case_colors <- setNames(brewer_cols, cases)
  
  all_dfs <- list()
  for(case_name in cases) {
    kegg_result <- results_list[[case_name]]$KEGG
    if (is.null(kegg_result) || nrow(as.data.frame(kegg_result)) == 0) {
      dummy <- data.frame(
        Description = paste0("No KEGG enrichment for ", case_name),
        GeneRatio = 0.01,
        Count = 0,
        p.adjust = NA,
        Case = case_name,
        stringsAsFactors = FALSE
      )
      all_dfs[[case_name]] <- dummy
    } else {
      df <- as.data.frame(kegg_result)
      df$Case <- case_name
      if ("GeneRatio" %in% colnames(df)) {
        df$GeneRatio <- sapply(df$GeneRatio, function(x) {
          if (is.character(x) && grepl("/", x)) {
            nums <- unlist(strsplit(x, "/"))
            return(as.numeric(nums[1]) / as.numeric(nums[2]))
          } else {
            return(as.numeric(x))
          }
        })
      }
      all_dfs[[case_name]] <- df
    }
  }
  combined_df <- bind_rows(all_dfs)
  
  ## Case_1, Case_2, Case_3, Case_4 from top to bottom
  # Reverse factor levels
  reversed_cases <- rev(cases)
  combined_df$Case <- factor(combined_df$Case, levels = reversed_cases)
  # Number the paths within each case
  combined_df <- combined_df %>%
    dplyr::group_by(Case) %>%
    dplyr::mutate(row_in_case = row_number()) %>%
    dplyr::ungroup()
  
  # Calculate offset in each case
  case_counts <- combined_df %>%
    dplyr::group_by(Case) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::arrange(Case) %>%
    dplyr::mutate(offset = lag(cumsum(n), default = 0))
  # Get y_value based on offset calculated above
  combined_df <- combined_df %>%
    dplyr::left_join(case_counts %>% dplyr::select(Case, offset), by = "Case") %>%
    dplyr::mutate(y_value = offset + row_in_case)
  
  ## Plotting
  # Base plot
  p <- ggplot(combined_df, aes(x = GeneRatio, y = y_value, color = Case, shape = Case)) +
    geom_point(aes(size = Count), na.rm = TRUE) +
    scale_color_manual(values = case_colors) +
    scale_shape_manual(values = c(0, 1, 2, 3)) +
    labs(x = "Gene Ratio", 
         y = "", 
         size = "Gene Count", 
         color = "Case") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  # Add case label in graph
  p <- p + geom_text_repel(aes(label = Description), 
                           na.rm = TRUE, size = 3, show.legend = FALSE)
  
  # Add dotted line
  for (i in 1:(nrow(case_counts) - 1)) {
    y_line <- sum(case_counts$n[1:i]) + 0.5
    p <- p + geom_hline(yintercept = y_line,
                        linetype = "dashed", 
                        color = "grey")
  }
  
  # Calculate case label position on Y-axis
  case_labels <- combined_df %>%
    dplyr::group_by(Case) %>%
    dplyr::summarize(y_mid = mean(y_value, na.rm = TRUE))
  
  # Increase left margin
  p <- p + 
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(t = 4, r = 4, b = 4, l = 80)
    )
  min_x <- min(combined_df$GeneRatio, na.rm = TRUE)
  p <- p + geom_text(
    data = case_labels,
    aes(x = min_x - 0.1, y = y_mid, label = Case),
    color = "black",
    size = 4,
    fontface = "bold",
    hjust = 1.6
  ) +
    ggtitle(title)
  print(p)
  
  return(p)
}
