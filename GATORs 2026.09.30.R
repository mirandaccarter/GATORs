#install and load packages function
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org/")
      library(pkg, character.only = TRUE)
    }
  }
}

# Required packages then load them
required_packages <- c("shiny", "DT", "dplyr", "ggplot2", "purrr", "reactable", 
                       "readr", "ggrepel", "ggtext", "tidyr",
                       "ggiraph","ggVennDiagram", "paletteer", "sf", "RSQLite")
install_and_load(required_packages)



# Ensure BiocManager is loaded to use version()
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org/")
}
if (BiocManager::version() < "3.23") {
  message("Current Bioconductor version is: ", BiocManager::version())
  message("Please update to at least version 3.23 to continue")
}

if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
  library(DESeq2)
}

if (!require("limma", quietly = TRUE)) {
  BiocManager::install("limma")
  library(limma)
}

if (!require("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
  library(enrichplot)
}

if (!require("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
}

# Define UI
ui <- fluidPage(
  titlePanel("GATORs: Gene enrichment Analysis Tool Of RNAseq"),
  fluidRow( #data input sections (directory)
    column(3, 
           wellPanel(
             h4("Data Input"),
             textInput("data_path",
                       "Enter Path to Directory containing CSV files:",
                       value = ""),
             helpText("Each file should contain counts data with locustags, columns designated as controls, and columns designated as experimental"),
             helpText("Treatment replicates should use underscore separator (e.g., 'drugA_1', 'drugA_2'). Columns with the same prefix will be grouped as replicates.")
             
           )),
    column(3, #data input sections (gene info)
           wellPanel(
             h4("Gene Library Input"),
             textInput("gene_library_path",
                       "Enter path to file:",
                       value = "")
           )),
    
    column(3, #data input sections (labels)
           wellPanel(
             h4("Label Input"),
             textInput("labels_path",
                       "Enter path to file:",
                       value = "")
           )),
    column(3, 
           wellPanel(
             h4("Control identifier"),
             textInput("control_pattern",
                       "Enter control identifier:",
                       placeholder = "e.g., 'control', 'WT', 'untreated")
           ))
  ),
  
  fluidRow( #button to load data 
    column(4, offset = 1,
           actionButton("load_data", "Load Data & Run Deseq2 analysis", class = "btn-primary"),),
    column(5, offset = 1,
           conditionalPanel(
             condition = "output.data_loaded",
             textOutput("filtering_summary") 
           ))),
  
  
  br(), br(),
  
  
  fluidRow( #enter row label section
    column(5, offset = 1,
           conditionalPanel(
             condition = "output.data_loaded",
             wellPanel(
               h4("Row Selection"),
               textInput("selected_label", 
                         "Enter Gene Number:", 
                         value = "",
                         placeholder = "Type a label name..."),
               helpText("Please use the SMU_#### format")
             )
           )
    ),
    column(5,
           conditionalPanel(
             condition = "output.data_loaded",
             wellPanel(
               h4("Graph Options"),
               checkboxGroupInput("columns_to_graph",
                                  "Select Conditions to Graph:",
                                  choices = NULL,
                                  selected = NULL)
           )))
  ),
  
  br(),
  
  fluidRow( 
    column(10, offset = 1,    #bar plot spot
           conditionalPanel(
             condition = "output.data_loaded",
             tabsetPanel(
               tabPanel("Regulation Across Experiments",
                        plotOutput("bar_plot", height = "800px"),
                        br(),
                        fluidRow(
                          
                          column(1,
                                 downloadButton("bars_jpg", 
                                                "Download Bar Plot (JPG)", 
                                                class = "btn-secondary")),
                          br(),
                          br(),
                        )),
               
               tabPanel("Additional Information", #table for p values and gene info
                        tableOutput("info_table")
               ))
           )),
    br(),
    br()
  ),
  fluidRow(
    column(10, offset = 1,
           conditionalPanel(
             condition = "output.data_loaded",
             tabsetPanel(
               tabPanel("All Data", #spot for big table
                        h4("All data"),
                        DT::dataTableOutput("data_table"),
                        br(),
                        div(style = "text-align:center; margin-top:20px;",
                            downloadButton("bigexport_csv",
                                           "Export Data to CSV",
                                           class = "btn-success",
                                           icon = icon("download"))),
                        br(),
                        br(),
                        br()),
               tabPanel("Volcano Plots",  #volacno plots
                        fluidRow(
                          column(6,
                                 wellPanel(
                                   h4("Volcano Plot Options"),  #choose volcano plot parameters
                                   selectInput("volcano_dataset", 
                                               "Select Dataset:", 
                                               choices = NULL),
                                   textInput("plot_title", 
                                             "Plot Title:", 
                                             value = "ex: Volcano Plot"),
                                   numericInput("label_top_n", 
                                                "Number of top/bottom genes to label:", 
                                                value = 18, min = 0, max = 50),
                                   checkboxInput("show_labels", "Show top/bottom gene labels", value = TRUE),
                                 )),
                          column(6,
                                 wellPanel(
                                   h4("Threshold Settings"),
                                   numericInput("logfc_threshold", 
                                                "LogFC Threshold:", 
                                                value = 2, min = 0.5, max = 5, step = 0.5),
                                   selectInput(
                                     "volcano_pvalue_type",
                                     "P-value displayed and used for significance:",
                                     choices = c(
                                       "Raw P-value" = "PValue",
                                       "Adjusted P-value (FDR)" = "padj"
                                     ),
                                     selected = "PValue"
                                   ),
                                   
                                   numericInput(
                                     "volcano_p_threshold",
                                     "Significance threshold:",
                                     value = 0.05,
                                     min = 0.0001,
                                     max = 0.1,
                                     step = 0.001
                                   )
                                  
                                 ))
                        ),
                        girafeOutput("volcano_plot", height = "600px", width = "1000px"), #display volcano plot
                        br(),
                        fluidRow(
                          column(6,
                                 wellPanel(
                                   h4("Select a log2 Fold Change cutoff:"),  #set parameters for differential genes table
                                   numericInput("graph_threshold",
                                                "",
                                                value = 2, min = 0.5, max = 15, step = 0.1))),
                          column(6,
                                 wellPanel(
                                   h4("Select a Significance threshold:"),
                                   numericInput("padj_threshold", 
                                                "P-Value Threshold:", 
                                                value = 0.05, min = 0.0001, max = 0.1, step = 0.001)
                                 
                          ))),
                        fluidRow(
                          column (10,
                                  reactable::reactableOutput("diff_genes_table"),
                                  br(),
                                  div(style = "text-align:center; margin-top:20px;",
                                      downloadButton("export_csv",
                                                     "Export Filtered Data to CSV",
                                                     class = "btn-success",
                                                     icon = icon("download"))),
                                  br()))
               ),
               
               tabPanel("Venn Diagrams",  #venn diagrams
                        fluidRow(
                          column(6,
                                 wellPanel(
                                   h4("Venn Diagram Options:"),  #choose venn diagram conditions
                                   checkboxGroupInput("conditions_to_overlap",
                                                      "Select Conditions to Include:",
                                                      choices = NULL,
                                                      selected=NULL))
                                  ),
                          column(6,    #select threshold for venn diagrams
                                 wellPanel(
                                     numericInput("venn_threshold", 
                                                "LogFC Threshold:", 
                                                value = 2, min = 0.5, max = 5, step = 0.5),
                                    actionButton("load_venn",
                                                 "Load Venn Diagrams",
                                                 class = "btn-primary"))  
                                 )),
                        fluidRow(
                          column(6,       #place for venn diagrams to go
                                 h4("Upregulated Genes:"),
                                 girafeOutput("venn_plot_up", 
                                            height = "500px", 
                                            width = "500px")),
                          column(6,
                                 h4("Downregulated Genes:"),
                                 girafeOutput("venn_plot_down", 
                                            height = "500px", 
                                            width = "500px"))
                        ),
                        fluidRow(
                          column(6,
                                 h4("Selected Upregulated Genes:"),    #When a section is clicked on table with genes apppears here
                                 reactable::reactableOutput("venn_up_table")),
                          column(6,
                                 h4("Selected Downregulated Genes:"),
                                 reactable::reactableOutput("venn_down_table"))
                        ),
                        br(),
                        fluidRow(  #download buttons
                          column(2,offset = 1,
                                 downloadButton("selected_venn_up_csv",
                                                "Export Data to CSV",
                                                class = "btn-success",
                                                icon = icon("download"))),
                          column(2,offset = 5,
                                 downloadButton("selected_venn_down_csv",
                                                "Export Data to CSV",
                                                class = "btn-success",
                                                icon = icon("download")))
                        )
                       ),
               tabPanel("Kegg Analysis",   #kegg analysis
                        fluidRow(
                          column(4, 
                                 wellPanel(
                                   h4("KEGG Enrichment Options"),    #give kegg organism code here
                                   textInput("kegg_organism",   
                                             "KEGG Organism Code:",
                                             placeholder = "e.g. smu, eco, bsu"),
                                   helpText("Find your code at genome.jp/kegg/catalog/org_list.html"),
                                   radioButtons("kegg_direction",      #choose to show up or downregulated genes
                                                "Gene regulation:",
                                                choices = c("Upregulated" = "up", 
                                                            "Downregulated" = "down"),
                                                selected = "up"),
                                   numericInput("kegg_lfc_cutoff",     #set a logfc cutoff
                                                "Log2FC Cutoff (upregulated):",
                                                value= 1, min =0, max = 5, step=0.5),
                                   numericInput("kegg_padj_cutoff",       #select a padj cutoff
                                                "Adjusted P-Value Cutoff:",
                                                value = 0.05, min = 0.001, max = 0.1, step = 0.005),
                                   numericInput("kegg_top_n",    #max pathways 
                                                "Max pathways to display(0 = all):",
                                                value = 20, min = 0, max = 100, step = 5),
                                   checkboxGroupInput("kegg_conditions",     #choose conditions to show
                                                      "Select Conditions to Include:",
                                                      choices = NULL,
                                                      selected = NULL),
                                   actionButton("run_kegg", "Run KEGG Analysis", class = "btn-primary"),
                                   br(), br(),
                                   helpText("Note: Requires internet connection to query KEGG database."),
                                   helpText("Only genes with an Adjusted P-value < 0.05 and greater than 1 Log2 FC will be included in the analysis")
                                 )),
                          column(8,
                                 conditionalPanel(    #show bubble plot
                                   condition = "output.kegg_ready",
                                   girafeOutput("kegg_bubble_plot", height = "700px", width = "100%"),
                                   br(),
                                   downloadButton("kegg_plot_download",   #download button
                                                  "Download Plot (JPG)",
                                                  class="btn-secondary"),
                                   br(), br(),
                                   h4("Enrichment Results Table (includes significant and non-significant results)"),
                                   reactable::reactableOutput("kegg_results_table"),   #show results table
                                   br(),
                                   downloadButton("kegg_csv_download",    #download button
                                                  "Export Results to CSV",
                                                  class = "btn-success")
                                   
                                 ),
                                 conditionalPanel(
                                   condition = "!output.kegg_ready",
                                   div(style = "padding: 40px; text-align: center; color:gray;",
                                       h4("Configure options and click 'Run KEGG Analysis'"))
                                 ))
                          
                        ))
             ))),
           conditionalPanel(
             condition = "!output.data_loaded",   #if data hasn't been loaded display this help text
             div(class = "text-left",
                 h4("Instructions:", 
                    style = "font-size: 28px"),
                 tags$ul(
                   tags$li("Place all CSV files in a single directory"),
                   tags$li("Each CSV should have a gene ID column labelled Locustag and count columns"),
                   tags$li("Control columns should contain the pattern you specify (ie. 'control')"),
                   tags$li("Treatment replicates should be named with the same prefix (e.g., 'drugA_1', 'drugA_2', 'drugA_3')"),
                   tags$li("The separator character '_' splits treatment name from replicate number, ensure there are no additional undescores in labels"),
                   tags$li("The app will use DESeq2 to compare each treatment group (all replicates) against the control within the same csv file"),
                   tags$li("DESeq2 will calculate normalized log2 fold changes and adjusted p-values"),
                   tags$li(strong("Gene library file should be tab-separated with UTF-8 encoding")),
                   tags$li(strong("All count values should be raw, non-normalized integer counts")),
                   tags$li(strong("Example column names: control_1, control_2, drugA_1, drugA_2, drugB_1, drugB_2"))
                 ),
                 style = "padding:50px; background-color: #f8f9fa; border-radius: 10px; margin-top: 20px; font-size: 20px")
           )
    )
  )
  


# Define the server for the shiny object
server <- function(input, output, session) {
  
  #assign variables to reactive values
  #these are needed soeverything in the server can access them
  values <- reactiveValues(
    big_table = NULL,
    numeric_cols = NULL,
    info_cols = NULL,
    label_mapping = NULL,
    reverse_label_mapping = NULL,  
    filtering_log = NULL,
    data_loaded = FALSE,
    individual_datasets = NULL,
    diff_genes_table = NULL,
    venn_ready = FALSE,
    venn_up_genes = NULL,
    venn_down_genes=NULL,
    kegg_results_df = NULL,
    all_kegg_results_df = NULL,
    kegg_ready = FALSE,
    gene_info = NULL
  )
  
  
  
  ######functions##########################################################
  
  
  
  #create a function to make labels 
  create_label_mapping <- function(numeric_cols, graph_labels) {
    #as part of this function these give info of whats happening
    cat("Graph labels column names:", paste(colnames(graph_labels), collapse = ", "), "\n")
    cat("Graph labels shortname column:", paste(graph_labels$shortname, collapse = ", "), "\n")
    cat("Graph labels longname column:", paste(graph_labels$longname, collapse = ", "), "\n")
    
    
    label_mapping <- setNames(numeric_cols, numeric_cols)  
    
    #loop to cycle through columns that have been defined as numeric (the logFC columns)
    for (i in seq_along(numeric_cols)) {
      col_name <- numeric_cols[i]
      cat("Processing column:", col_name, "\n")
      
      #if the column name starts with logFC_ find the shortname
      if (startsWith(col_name, "logFC_")) {
        partial_name <- sub("^logFC_", "", col_name)
        cat("Looking for shortname:", partial_name, "\n")
        
        #Assign the shortname in graphlabels to matching row if it equals the partial name
        matching_row <- which(graph_labels$shortname == partial_name)
        cat("Matching rows found:", length(matching_row), "\n")
        
        #if the length of the matching row is greater than 0 assign the longname value to longname_value
        if (length(matching_row) > 0) {
          longname_value <- graph_labels$longname[matching_row[1]]
          cat("Found longname:", longname_value, "\n")
          
          #if long name value is not blank map longname value to the column name in label mapping
          if (!is.na(longname_value) && longname_value != "") {
            label_mapping[col_name] <- longname_value
            cat("Successfully mapped", col_name, "to", longname_value, "\n")
          } else {
            cat("Warning: Empty longname for", col_name, "\n")
          }
        } else {
          cat("Warning: No matching shortname found for", partial_name, "\n")
        }
      }
    }
    return(label_mapping)
  }
  
  #Function to create reverse mapping (longname -> shortname)
  create_reverse_mapping <- function(label_mapping) {
    reverse_mapping <- setNames(names(label_mapping), unname(label_mapping))
    return(reverse_mapping)
  }
  
  
  #function for deseq2 analysis
  
  #function requires the all_counts_list to function
  run_deseq2 <- function(all_counts_list) {
    cat("Running DESeq analysis per dataset...\n")
    results_list <- list()
    
    #for loop to go through the list of counts for every unique file name
    for (file_name in names(all_counts_list)) {
      #assign the file in the list thats being worked with to a variable
      file_data <- all_counts_list[[file_name]]
      cat("Processing experiment file:", file_name, "\n")
      
      #create a matrix with the counts data in the working file
      count_matrix <- round(file_data$counts)
      
      #assign gene ids to the rownames and sample namens to the columns
      rownames(count_matrix) <- file_data$gene_ids
      colnames(count_matrix) <- file_data$sample_names
      
      # Identify genes with zero counts across all samples
      all_genes_in_file <- rownames(count_matrix)
      zero_count_genes <- rownames(count_matrix[rowSums(count_matrix) == 0, , drop = FALSE])
      
      # Log the loss
      cat("Removing", length(zero_count_genes), "genes with zero counts in file:", file_name, "\n")
      values$filtering_log <- c(values$filtering_log, 
                                paste0(file_name, ": Removed ", length(zero_count_genes), " zero-count genes"))
      
      #filter the matrix
      count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]
      
      #create a dataframe containing conditions and samples names only
      colData <- data.frame(
        condition = factor(file_data$conditions),
        row.names = file_data$sample_names
      )
      
      #if there is no identified control condition throw an error
      if (!"control" %in% levels(colData$condition)) {
        warning(paste("No 'control' condition found in file:", file_name))
        next
      }
      
      #relevel the conditions to have control first
      colData$condition <- relevel(colData$condition, ref = "control")
      
      # Build and fit the DESeq2 model 
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = colData,
        design = ~ condition
      )
      
      # Using fitType = "local" and minmu prevents dispersion collapse on small counts
      dds <- DESeq(dds, 
                   fitType = "parametric",
                   quiet = TRUE)
      
      #create a variable using only the non-control variables
      treatment_groups <- unique(file_data$conditions[file_data$conditions != "control"])
      
      #for loop to go through unique treatment types
      for (treatment in treatment_groups) {
        cat("  Extracting results for:", treatment, "\n")
        
        #for a given treatment pull out the deseq2 results
        res <- results(dds, contrast = c("condition", treatment, "control"))
        
        #gets the baseMean cutoff deseq used
        #this is where the padjust values can be set to null if they don't meet the basemena score 
        filter_threshold <- metadata(res)$filterThreshold
        filter_threshold <- if (length(filter_threshold) > 0) unname(filter_threshold) else NA
        
        #puts the deseq results into a dataframe
        res_df <- as.data.frame(res)
        res_df$gene_id <- rownames(res_df)
        
        
        #this is supposed to tell if the glm parameter estimation converged, (i think this was added for troubleshootingsomething)
        beta_converged <- if ("betaConv" %in% colnames(mcols(dds))) {
          mcols(dds)$betaConv
        } else {
          rep(TRUE, nrow(dds))
        }
        names(beta_converged) <- rownames(dds)
        
        
        #finds samples involved in contrast
        samples_in_comparison <- rownames(colData(dds))[colData(dds)$condition %in% c(treatment, "control")]
        #finds cooks distance
        cooks_mat <- assays(dds)[["cooks"]][, samples_in_comparison, drop = FALSE]
        
        max_cooks <- apply(cooks_mat, 1, max, na.rm = TRUE)
        
        #sets outlier boudnary
        n_reps <- length(samples_in_comparison)
        p_params <- ncol(model.matrix(design(dds), colData(dds)))
        cooks_cutoff <- if (n_reps > p_params) qf(0.99, p_params, n_reps - p_params) else Inf
        is_cooks_outlier <- max_cooks > cooks_cutoff
        
        #finds genes with too wide of a dispersion
        disp_outlier <- if ("dispOutlier" %in% colnames(mcols(dds))) {
          mcols(dds)$dispOutlier
        } else {
          rep(FALSE, nrow(dds))
        }
        names(disp_outlier) <- rownames(dds)
        
        
        #these are all the possible reasons why a padjust value might be na and puts them into the results list
        res_df$na_reason <- NA_character_
        is_na_pval <- is.na(res_df$pvalue)
        is_na_padj <- is.na(res_df$padj)
        
        needs_reason <- is_na_pval | is_na_padj
        for (i in which(needs_reason)) {
          gid <- res_df$gene_id[i]
          
          if (res_df$baseMean[i] == 0) {
            res_df$na_reason[i] <- "Not detected - zero counts in this comparison"
          } else if (!isTRUE(beta_converged[gid])) {
            res_df$na_reason[i] <- "Model fit did not converge (unstable count pattern - check raw counts for this gene)"
          } else if (isTRUE(is_cooks_outlier[gid])) {
            res_df$na_reason[i] <- paste0(
              "Outlier - one replicate's count is extreme relative to others (max Cook's distance: ",
              round(max_cooks[gid], 2), ", cutoff: ", round(cooks_cutoff, 2), ")"
            )
          } else if (is_na_padj[i] && !is_na_pval[i]) {
            res_df$na_reason[i] <- paste0(
              "Low expression - filtered out by independent filtering (baseMean = ",
              round(res_df$baseMean[i], 2), ", threshold = ", round(filter_threshold, 2), ")"
            )
          } else {
            res_df$na_reason[i] <- paste0(
              "Could not compute a valid test (baseMean = ", round(res_df$baseMean[i], 2),
              ", dispersion outlier: ", isTRUE(disp_outlier[gid]),
              ", max Cook's distance: ", round(max_cooks[gid], 2), ")"
            )
          }
        }
        
        results_list[[treatment]] <- data.frame(
          gene_id = res_df$gene_id,
          logFC = res_df$log2FoldChange,
          PValue = res_df$pvalue,
          padj = res_df$padj,
          baseMean = res_df$baseMean,
          na_reason = res_df$na_reason
        )
      }
    }
    return(results_list)
  }

  
#function to use the kegg analysis button  
  normalize2kegg<- function(ids, gene_info) {
    
    cat("\n=== NORMALIZE2KEGG ===\n")
    cat("Input IDs (first 5):", paste(head(ids, 5), collapse = ", "), "\n")
    
    
    # Debug gene_info input
    cat("gene_info class:", class(gene_info), "\n")
    cat("gene_info is null:", is.null(gene_info), "\n")
    
    #if the gene info is empty throw this error
    if (is.null(gene_info)) {
      cat("gene_info is NULL — returning IDs unchanged\n")
      return(ids)
    }
    
    #give info on whats happening
    cat("gene_info dims:", nrow(gene_info), "x", ncol(gene_info), "\n")
    cat("gene_info columns:", paste(colnames(gene_info), collapse = ", "), "\n")
    cat("kegg_gene column exists:", "kegg_gene" %in% colnames(gene_info), "\n")
    
    #give info about whats happening if the kegg gene is in th gene info file
    if ("kegg_gene" %in% colnames(gene_info)) {
      cat("kegg_gene sample values:", paste(head(gene_info$kegg_gene, 5), collapse = ", "), "\n")
      cat("kegg_gene non-empty:", sum(!is.na(gene_info$kegg_gene) & gene_info$kegg_gene != ""), "\n")
    }
    
    #explain that the kegg gene is not found in the gene info sheet
    if (!"kegg_gene" %in% colnames(gene_info)) {
      cat("No kegg_gene column found in gene info — returning IDs unchanged\n")
      return(ids)
    }
    
    #pick out specific columns for a lookup dataframe
    lookup <- gene_info %>%
      filter(!is.na(kegg_gene), kegg_gene != "") %>%
      select(Locustag, kegg_gene) %>%
      distinct()
    
    cat("Lookup entries available:", nrow(lookup), "\n")
    
    #look up the locustag in the matching row as the kegg_gene
    normalized <- sapply(ids, function(id) {
      match_row <- lookup$kegg_gene[lookup$Locustag == id]
      if (length(match_row) > 0 && !is.na(match_row[1])) {
        return(match_row[1])
      }
      return(id)  # no mapping found, return original
    })
    
    #find number of genes 
    after_matches <- sum(normalized %in% lookup$kegg_gene)
    cat("Successfully mapped:", after_matches, "/", length(ids), "\n")
    cat("Coverage:", round(after_matches / length(ids) * 100, 1), "%\n")
    
    #genes with no kegg match
    still_missing <- ids[!normalized %in% lookup$kegg_gene]
    if (length(still_missing) > 0) {
      cat("Unmapped (", length(still_missing), "):",
          paste(head(still_missing, 5), collapse = ", "),
          if (length(still_missing) > 5) "..." else "", "\n")
      cat("These likely have no KEGG annotation\n")
    }
    cat("=== END NORMALIZE2KEGG ===\n\n")
    
    return(unname(normalized))
  }
  
  #create the differentially expressed gene table at the bottom of the volcano plot
  diff_genes_reactive <- reactive({
    req(input$volcano_dataset, values$individual_datasets) # Make sure these inputs exist
    req(input$graph_threshold, input$padj_threshold)  # Make sure these inputs exist
    
    # Get the selected dataset
    selected_data <- values$individual_datasets[[input$volcano_dataset]]
    
    if (is.null(selected_data)) {
      return(NULL)
    }
    
    # Create a copy to avoid modifying the original
    diff_genes <- selected_data
    
    # Calculate differential expression based on current thresholds using padj
    diff_genes$diffexpressed <- "NO"
    if ("padj" %in% colnames(diff_genes)) {
      diff_genes$diffexpressed[!is.na(diff_genes$padj) & 
                                 diff_genes$logFC > input$graph_threshold & 
                                 diff_genes$padj < input$padj_threshold] <- "UP"
      diff_genes$diffexpressed[!is.na(diff_genes$padj) & 
                                 diff_genes$logFC < -input$graph_threshold & 
                                 diff_genes$padj < input$padj_threshold] <- "DOWN"
    }
    
    # Filter for only UP and DOWN regulated genes
    diff_genes_filtered <- subset(diff_genes, diffexpressed == "UP" | diffexpressed == "DOWN")
    
    return(diff_genes_filtered)
  })
  
  
  
  # Function to create volcano plot
  create_volcano_plot <- function() {
    req(input$volcano_dataset, values$individual_datasets)   # Make sure these inputs exist
    req(input$volcano_pvalue_type, input$volcano_p_threshold)   # Make sure these inputs exist
    
    #assign selected dataset to a value
    selected_data <- values$individual_datasets[[input$volcano_dataset]]   
    
    #if this is empty give nothing
    if (is.null(selected_data) || nrow(selected_data) == 0) {
      return(NULL)
    }
    
    #choose a p-value type (padj or p value)
    p_column <- input$volcano_pvalue_type
    p_values <- selected_data[[p_column]]
    
    #cleans up data to plot normally to remove NAs and - p values
    selected_data <- selected_data[
      !is.na(selected_data$logFC) &
        is.finite(selected_data$logFC) &
        !is.na(p_values) &
        is.finite(p_values) &
        p_values >= 0,
      ,
      drop = FALSE
    ]
    #give the p value column to a variable
    p_values <- selected_data[[p_column]]
    
    #manages any 0 p values and replaces with a minimum p value
    nonzero_p <- p_values[p_values > 0]
    min_p <- if (length(nonzero_p) > 0) min(nonzero_p) else 1e-50
    p_values_adjusted <- ifelse(p_values == 0, min_p, p_values)
    
 
    
    # Classify genes using whichever statistic the user selected
    selected_data$diffexpressed <- "NO"
    
    is_significant <- p_values < input$volcano_p_threshold
    
    selected_data$diffexpressed[
      is_significant &
        selected_data$logFC >= input$logfc_threshold
    ] <- "UP"
    
    selected_data$diffexpressed[
      is_significant &
        selected_data$logFC > 1 &
        selected_data$logFC < input$logfc_threshold
    ] <- "increase"
    
    selected_data$diffexpressed[
      is_significant &
        selected_data$logFC <= -input$logfc_threshold
    ] <- "DOWN"
    
    selected_data$diffexpressed[
      is_significant &
        selected_data$logFC < -1 &
        selected_data$logFC > -input$logfc_threshold
    ] <- "decrease"
    
    #set up gene labels to be gene name or locustag
    gene_labels <- ifelse(
      !is.na(selected_data$Gene) &
        selected_data$Gene != "" &
        selected_data$Gene != "Unknown",
      selected_data$Gene,
      selected_data$Locustag
    )
    
    # Label only significant up- and downregulated genes
    selected_data$label <- ""
    
    #checks that show labels box is checked
    if (isTRUE(input$show_labels) && input$label_top_n > 0) {
      significant_rows <- which(
        selected_data$diffexpressed %in% c("UP", "DOWN"))
      
      #finds rows where genes pass threshholds
      if (length(significant_rows) > 0) {
        up_rows <- significant_rows[
          selected_data$diffexpressed[significant_rows] == "UP"]
        
        down_rows <- significant_rows[
          selected_data$diffexpressed[significant_rows] == "DOWN"]
        
        #sorts genes in descending order
        top_up <- if (length(up_rows) > 0) {
          ordered_up <- up_rows[
            order(selected_data$logFC[up_rows], decreasing = TRUE)]
          
          ordered_up[
            seq_len(min(input$label_top_n, length(ordered_up)))]
        } else { integer(0)}
        
        #sorts genes in ascending order
        top_down <- if (length(down_rows) > 0) {
          ordered_down <- down_rows[
            order(selected_data$logFC[down_rows], decreasing = FALSE)]
          
          ordered_down[
            seq_len(min(input$label_top_n, length(ordered_down)))]
        } else { integer(0)}
        
        #creates final label rows
        label_rows <- c(top_up, top_down)
        selected_data$label[label_rows] <- gene_labels[label_rows]
      }
    }
    
    #build hover box
    selected_data$tooltip_text <- paste0(
      "<b>", selected_data$Locustag, "</b>\n",
      "Gene: ",    selected_data$Gene, "\n",
      "Log2FC: ",  round(selected_data$logFC, 3), "\n",
      "P-value: ", formatC(selected_data$PValue, format = "e", digits = 3), "<br>",
      "Adjusted P-value: ", formatC(selected_data$padj, format = "e", digits = 3), "<br>",
      selected_data$Description
    )
    
    #set up y axis label to chnage based on pvalue
    y_axis_label <- if (p_column == "PValue") {
      expression(bold("-log"[10] * "(P-value)"))
    } else {
      expression(bold("-log"[10] * "(adjusted P-value)"))
    }
    
    # Create the volcano plot
    volcano_plot <- ggplot(data = selected_data, 
                           aes(x = logFC, 
                               y = -log10(pmax(.data[[p_column]], min_p)),
                               col = diffexpressed)) +
      geom_vline(xintercept = c(-input$logfc_threshold, -1, 1, input$logfc_threshold), 
                 col = "gray", 
                 linetype = 'dashed') +
      geom_hline(yintercept = -log10(input$volcano_p_threshold), 
                 col = "gray", 
                 linetype = 'dashed') + 
      ggtitle(input$plot_title) +
      labs(
        color = "Expression",
        x = expression(bold("log"[2] * "FC")),
        y = if (p_column == "PValue") {
          expression(bold("-log"[10] * "(P-value)"))
        } else {
          expression(bold("-log"[10] * "(adjusted P-value)"))
        }
      )  +
      geom_point_interactive(aes(tooltip=tooltip_text,
                                 data_id = Locustag),
                             size = 2) +
      geom_text_repel(data = selected_data[selected_data$label != "",],
                      aes(x = logFC, 
                          y = -log10(pmax(.data[[p_column]], min_p)),
                          label = label,
                          col = diffexpressed),
                      inherit.aes = FALSE,
                      max.overlaps = 20,
                      size =4)+
      scale_y_continuous(expand = expansion(mult= c(0.05,0.1)))+
      
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5,
                                      face = "bold",
                                      size = 20),
            legend.position = "top",
            plot.clip = "off",
            legend.title = element_blank(),
            legend.text = element_text(size = 10, face= "bold"),
            axis.title = element_text(size = 16, color = "black", face="bold"),
            axis.text = element_text(size = 12, color = "black", face="bold")) +
      scale_color_manual(values = c("DOWN" = "blue", 
                                    "increase" = "palevioletred1",
                                    "NO" = "grey", 
                                    "decrease" = "lightskyblue",
                                    "UP" = "red"),
                         labels = c("DOWN" = "Highly Downregulated", 
                                    "NO" = "Not Significant", 
                                    "UP" = "Highly Upregulated",
                                    "increase" = "Slightly Upregulated",
                                    "decrease" = "Slightly Downregulated"),
                         breaks = c("DOWN", "decrease", "NO", "increase", "UP"))
    
    return(volcano_plot)
  }
  
  
  
  
  # Function to create bar plot
  create_bar_plot <- function() {
    req(input$selected_label, input$columns_to_graph)
    
    #if the selected row name is blank return nothing
    if (!input$selected_label %in% values$big_table$Locustag) {
      return(NULL)
    }
    
    #assign the row in the big table which matches with the locustag in the selected label
    selected_row <- values$big_table[values$big_table$Locustag == input$selected_label, ]
    #variable containg checked boxes from columns to rgaph
    selected_columns <- input$columns_to_graph
    
    # Check if all selected columns exist in the data
    missing_cols <- selected_columns[!selected_columns %in% colnames(values$big_table)]
    if (length(missing_cols) > 0) {
      return(NULL)
    }
    
    # Get data for selected columns
    plot_data <- selected_row[, selected_columns, drop = FALSE]
    
    # Convert to long format for plotting
    plot_data_long <- data.frame(
      Variable = names(plot_data),
      Value = as.numeric(plot_data[1, ])
    )
    
    # Map variable names to longnames for display
    plot_data_long$DisplayName <- sapply(plot_data_long$Variable, function(x) {
      if (x %in% names(values$label_mapping)) {
        return(values$label_mapping[[x]])
      } else {
        return(x)
      }
    })
    
    # Create ggplot bar chart using DisplayName for x-axis
    bar_plot <- ggplot(plot_data_long, aes(x = DisplayName, y = Value)) +
      geom_col(fill = "steelblue", color = "darkblue", width = 0.7) +
      geom_text(aes(y = Value + 0.2 * sign(Value), label = round(Value, 4), fontface = "bold")) +
      labs(
        title = paste("Data for:", input$selected_label),
        x = "Experiment",
        y = "log2FC"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))
    
    return(bar_plot)
  }
  
  #function to create venn diagram
  
  #make the upregulated venn diagram
  create_venn_up<- function(selected_keys, venn_threshold) {
    
    #if user picked less than two conditions throw an error
    if (length(selected_keys)< 2 ){
      showNotification("Please select at least 2 conditions", type = "warning")
      return(NULL)
    }
    selected_data <- values$individual_datasets[selected_keys]
    
    print(names(selected_data))
    
    venn_up <-list()
    
    #for loop to go through all the names in selected data for upregulated genes
    for (name in names(selected_data)) {
      df <- selected_data[[name]]
      
      df$diffexpressed <- "NO"
      df$diffexpressed[!is.na(df$padj) &
                         df$logFC > input$venn_threshold &
                         df$padj < 0.05] <- "UP"
      
      up_genes <- df$Locustag[df$diffexpressed == "UP"]
      
      cat(name, "UP genes:", length(up_genes), "\n")
      
      if (length(up_genes) > 0) {
        venn_up[[name]] <- up_genes
      }
      
    }
    
    cat("Selected conditions:\n")
    print(input$conditions_to_overlap)
    
    cat("Names in selected_data:\n")
    print(names(selected_data))
    
    cat("Venn UP set sizes:\n")
    print(sapply(venn_up, length))
    
    print(input$conditions_to_overlap)
    print(names(values$individual_datasets))
    
    # Check again after filtering empty sets
    if (length(venn_up) < 2) {
      showNotification("Not enough overlapping gene sets (check thresholds)", type = "warning")
      return(NULL)
    }
    
    # Store sets for click observer to use
    values$venn_up_sets <- venn_up
    
    # Build the venn data manually so we can make labels interactive
    venn_data <- Venn(venn_up)
    venn_processed <- process_data(venn_data)
    
    
    # Extract all data as regular dataframes
    regions <- venn_regionedge(venn_processed)
    region_labels <- venn_regionlabel(venn_processed)
    set_labels <- venn_setlabel(venn_processed)
    
    cat("region_labels columns:", paste(colnames(region_labels), collapse=", "), "\n")
    cat("set_labels columns:", paste(colnames(set_labels), collapse=", "), "\n")
    
    regions <- venn_regionedge(venn_processed)
    cat("regionedge class:", class(regions), "\n")
    cat("regionedge columns:", paste(colnames(regions), collapse=", "), "\n")
    cat("regionedge head:\n")
    print(head(regions))
    
    venn_plot_up <- ggplot() +
      # Polygon regions - regionedge has all boundary points already
      geom_polygon(data = venn_regionedge(venn_processed),
                   aes(x = X, y = Y, group = id, fill = count),
                   color = "black", alpha = 0.6) +
      # Set name labels
      geom_text(data = venn_setlabel(venn_processed),
                aes(x = X, y = Y, label = name),
                size = 4, fontface = "bold") +
      # Interactive clickable label boxes
      geom_label_interactive(
        data = venn_regionlabel(venn_processed),
        aes(x = X, y = Y,
            label = count,
            tooltip = paste0(name, "\n", count, " genes\nClick to view"),
            data_id = id),
        size = 6,
        fill = "white",
        label.size = 1.5
      ) +
      scale_fill_paletteer_c("grDevices::Pastel 1", direction = 1) +
      theme_void()
    
    return(venn_plot_up)
    
  } 
  
   #same as above except with downregulated genes
  create_venn_down<- function(selected_keys, venn_threshold) {
    
    
    if (length(selected_keys)< 2 ){
      showNotification("Please select at least 2 conditions", type = "warning")
      return(NULL)
    }
    selected_data <- values$individual_datasets[selected_keys]
    
    print(names(selected_data))
    
    
    venn_down<-list()
    
    for (name in names(selected_data)) {
      df <- selected_data[[name]]
      
      df$diffexpressed <- "NO"
      
      df$diffexpressed[!is.na(df$padj) &
                         df$logFC < -input$venn_threshold &
                         df$padj < 0.05] <- "DOWN"
      
      
      down_genes <- df$Locustag[df$diffexpressed == "DOWN"]
      
      cat(name, "DOWN genes:", length(down_genes), "\n")
      
      
      if (length(down_genes) > 0) {
        venn_down[[name]] <- down_genes
      }
    }
    
    cat("Selected conditions:\n")
    print(input$conditions_to_overlap)
    
    cat("Names in selected_data:\n")
    print(names(selected_data))
    
    
    cat("Venn DOWN set sizes:\n")
    print(sapply(venn_down, length))
    
    print(input$conditions_to_overlap)
    print(names(values$individual_datasets))
    
    # Check again after filtering empty sets
    if (length(venn_down) < 2) {
      showNotification("Not enough overlapping gene sets (check thresholds)", type = "warning")
      return(NULL)
    }
    
    
    # Store sets for click observer to use
    values$venn_down_sets <- venn_down
    
    # Build the venn data manually so we can make labels interactive
    venn_data <- Venn(venn_down)
    venn_processed <- process_data(venn_data)
    
    
    # Extract all data as regular dataframes - no sf
    regions <- venn_regionedge(venn_processed)
    region_labels <- venn_regionlabel(venn_processed)
    set_labels <- venn_setlabel(venn_processed)
    
    cat("region_labels columns:", paste(colnames(region_labels), collapse=", "), "\n")
    cat("set_labels columns:", paste(colnames(set_labels), collapse=", "), "\n")
    
    regions <- venn_regionedge(venn_processed)
    cat("regionedge class:", class(regions), "\n")
    cat("regionedge columns:", paste(colnames(regions), collapse=", "), "\n")
    cat("regionedge head:\n")
    print(head(regions))
    
    
    
    venn_plot_down<- ggplot() +
      # Polygon regions - regionedge has all boundary points already
      geom_polygon(data = venn_regionedge(venn_processed),
                   aes(x = X, y = Y, group = id, fill = count),
                   color = "black", alpha = 0.6) +
      # Set name labels
      geom_text(data = venn_setlabel(venn_processed),
                aes(x = X, y = Y, label = name),
                size = 4, fontface = "bold") +
      # Interactive clickable label boxes
      geom_label_interactive(
        data = venn_regionlabel(venn_processed),
        aes(x = X, y = Y,
            label = count,
            tooltip = paste0(name, "\n", count, " genes\nClick to view"),
            data_id = id),
        size = 6,
        fill = "white",
        label.size = 1.5
      ) +
      scale_fill_paletteer_c("grDevices::Pastel 1", direction = 1) +
      theme_void() 
    
    
    return(venn_plot_down)
    
  } 
  
  
  
  
  #load in the data this is where things really start
  observeEvent(input$load_data, {
    cat("=== STARTING DATA LOAD PROCESS ===\n")
    
    
    #tryCatch helps find errors this one looks for if the inputs are blank
    tryCatch({
      if (input$data_path == "" || input$gene_library_path == "" || input$labels_path == "") {
        cat("ERROR: Missing file paths\n")
        cat("Data path:", input$data_path, "\n")
        cat("Gene library path:", input$gene_library_path, "\n") 
        cat("Labels path:", input$labels_path, "\n")
        showNotification("Please provide all file paths", type = "error")
        return()
      }
      
      cat("All paths provided, proceeding...\n")
      
      
      #read in graph labels with a potential error message if there is an issue
      graph_labels <- tryCatch({
        read.csv(input$labels_path, stringsAsFactors = FALSE, strip.white= TRUE, encoding = "UTF-8")
      }, error = function(e) {
        showNotification(paste("Error reading labels file:", e$message), type = "error")
        return(NULL)
      })
      
      #if graph labels are null return nothing
      if (is.null(graph_labels)) return()
      
      #look at the gene library file, and give a few details about the beginning of the file
      gene_info <- tryCatch({
        cat("Examining gene library file:", input$gene_library_path, "\n")
        first_few_lines <- readLines(input$gene_library_path, n = 5, encoding = "UTF-8")
        cat("First few lines of gene library:\n")
        for(line_num in 1:length(first_few_lines)) {
          cat("Line", line_num, ":", first_few_lines[line_num], "\n")
          comma_count <- lengths(regmatches(first_few_lines[line_num], gregexpr(",", first_few_lines[line_num])))
          tab_count <- lengths(regmatches(first_few_lines[line_num], gregexpr("\t", first_few_lines[line_num])))
          cat(" Comma count:", comma_count, " Tab count:", tab_count, "\n")
        }
        
        file_ext <- tools::file_ext(input$gene_library_path)
        cat("File extension:", file_ext, "\n")
        
        header_line <- first_few_lines[1]
        has_tabs <- grepl("\t", header_line)
        has_commas <- grepl(",", header_line)
        
        
        #if the gene info file is a .tsv file (WHICH IT SHOULD BE) read as a .tsv file
        if (file_ext %in% c("tsv", "txt") || (has_tabs && !has_commas)) {
          cat("Reading as TSV file (tab-separated)\n")
          read.delim(input$gene_library_path, 
                     stringsAsFactors = FALSE, 
                     strip.white = TRUE,
                     fill = TRUE,
                     blank.lines.skip = TRUE,
                     check.names = TRUE,
                     comment.char = "", 
                     encoding = "UTF-8")
          #otherwise read it as a .csv file
        } else {
          cat("Reading as CSV file (comma-separated)\n")
          read.csv(input$gene_library_path, 
                   stringsAsFactors = FALSE, 
                   strip.white = TRUE,
                   fill = TRUE,
                   blank.lines.skip = TRUE,
                   check.names = TRUE,
                   comment.char = "", 
                   encoding = "UTF-8")
        }
      }, error = function(e) {
        #if it has issues then spit out an error
        showNotification(paste("Error reading gene library file:", e$message), type = "error")  
        return(NULL) 
      }) 
      
      
      
      if (is.null(gene_info)) return()
      
      
      
      #make a list of all the csv files in the data_path input directory
      files <- list.files(path = input$data_path, pattern = "*.csv", full.names = TRUE)
      
      
      #find the label filename and remove labels file from data directory if it is in the same file
      labels_filename <- basename(input$labels_path)
      files <- files[!basename(files) %in% labels_filename]
      
      
      #list CSV files to process
      cat("Found", length(files), "CSV files to process (excluding labels file):\n")
      for(i in seq_along(files)) {
        cat(" ", i, ":", basename(files[i]), "\n")
      } 
      if (length(files) == 0) {
        showNotification("No CSV files found in directory", type = "error")
        return()
      }
      
      
      #if the number of labels found is not the same as number of files spit out a warning
      if(nrow(graph_labels) < length(files)) {  
        showNotification(paste("Warning: Only", nrow(graph_labels), "labels found for", length(files), "files"), type = "warning") 
      }
      
      # First pass: collect all data 
      all_counts_list <- list()
      
      for (i in 1:length(files)) {
        file_name <- tools::file_path_sans_ext(basename(files[i]))
        cat("Processing file", i, ":", file_name, "\n")  
        
        counts_data <- read.csv(files[i],  
                                stringsAsFactors = FALSE, 
                                check.names = FALSE,
                                strip.white = TRUE,  
                                fill = TRUE,
                                blank.lines.skip = TRUE, 
                                encoding = "UTF-8")
        
        cat("Successfully read file with", nrow(counts_data), "rows and",  ncol(counts_data), "columns\n")  
        cat("Column names:", paste(colnames(counts_data), collapse = ", "), "\n\n") 
        
        # Check for required columns
        if (!"Locustag" %in% colnames(counts_data)) {  
          cat("ERROR: No 'Locustag' column found in file:", files[i], "\n")  
          stop("Missing required Locustag column")
        }
        # Get count columns
        count_cols <- colnames(counts_data)[sapply(counts_data, is.numeric)]
        control_cols <- grep(input$control_pattern, count_cols, ignore.case = TRUE, value = TRUE)
        treatment_cols <- setdiff(count_cols, control_cols)
        
        cat("Control columns:", paste(control_cols, collapse = ", "), "\n")
        cat("Treatment columns:", paste(treatment_cols, collapse = ", "), "\n")
        
        # Group treatments by prefix
        get_treatment_group <- function(col_name) {
          parts <- strsplit(col_name, "_", fixed = TRUE)[[1]]
          if (length(parts) > 1) {
            return(paste(parts[-length(parts)], collapse = "_"))
          } else {
            return(col_name)
          }
        }
        
        # Store data 
        sample_names <- c(control_cols, treatment_cols)
        conditions <- c(rep("control", length(control_cols)),
                        sapply(treatment_cols, get_treatment_group))
        
        all_counts_list[[file_name]] <- list(
          gene_ids = counts_data$Locustag,
          counts = as.matrix(counts_data[, sample_names]),
          sample_names = sample_names,
          conditions = conditions,
          control_cols = control_cols,
          treatment_cols = treatment_cols
        )
      }
      
      # Run DESeq2 
      cat("\n=== Running DESeq2 analysis ===\n")
      deseq_results <- run_deseq2(all_counts_list)
      
      cat("DESeq2 result names:", paste(names(deseq_results), collapse = ", "), "\n")
      cat("Expected shortnames:", paste(graph_labels$shortname, collapse = ", "), "\n")
      
      # Create individual datasets and combined table
      individual_datasets <- list()
      data_list <- list()
      
      for (treatment_name in names(deseq_results)) {
        cat("Processing results for:", treatment_name, "\n")
        
        deseq_res <- deseq_results[[treatment_name]]
        colnames(deseq_res)[1] <- "Locustag"
        
        # Join with gene info
        dataset <- left_join(deseq_res, gene_info, by = "Locustag")
        
        # Reorder columns
        dataset <- dataset %>%
          select(Locustag, Gene, Description, Functional_Group, 
                 logFC, PValue, padj, na_reason, baseMean, everything())
        
        # Add differential expression categories
        dataset$diffexpressed <- "NO"
        dataset$diffexpressed[dataset$logFC > 2 & dataset$padj < 0.05] <- "UP"
        dataset$diffexpressed[dataset$logFC < 2 & dataset$logFC > 1 & dataset$padj < 0.05] <- "increase"
        dataset$diffexpressed[dataset$logFC < -2 & dataset$padj < 0.05] <- "DOWN"
        dataset$diffexpressed[dataset$logFC > -2 & dataset$logFC < -1 & dataset$padj < 0.05] <- "decrease"
        
        # Fill missing gene info
        dataset <- dataset %>%
          mutate(
            Gene = ifelse(is.na(Gene) | Gene == "", "Unknown", Gene),
            Description = ifelse(is.na(Description) | Description == "", "No description available", Description),
            Functional_Group = ifelse(is.na(Functional_Group) | Functional_Group == "", "Unclassified", Functional_Group)
          )
        
        individual_datasets[[treatment_name]] <- dataset
        
        # Prepare for combined table
        temp_df <- data.frame(Locustag = deseq_res$Locustag)
        temp_df[[paste0("logFC_", treatment_name)]] <- deseq_res$logFC
        temp_df[[paste0("PValue_", treatment_name)]] <- deseq_res$PValue
        temp_df[[paste0("padj_", treatment_name)]] <- deseq_res$padj
        temp_df[[paste0("no_padj_bc_", treatment_name)]] <-deseq_res$na_reason
        
        data_list[[treatment_name]] <- temp_df
      }
      
      # Combine all data
      big_table <- data_list %>%
        purrr::reduce(full_join, by = "Locustag")
      
      # Join with gene info
      big_table <- left_join(big_table, gene_info, by = "Locustag") %>%
        dplyr::select(Locustag, Gene, Description, Functional_Group, everything())
      
      cat("Filling missing gene information...\n")
      
      # Fill missing info
      big_table <- big_table %>%
        mutate(
          Gene = ifelse(is.na(Gene) | Gene == "", "Unknown", Gene),
          Description = ifelse(is.na(Description) | Description == "", "No description available", Description),
          Functional_Group = ifelse(is.na(Functional_Group) | Functional_Group == "", "Unclassified", Functional_Group)
        )
      
      #identify columns for graphing
      cat("Identifying numeric and info columns...\n")
      numeric_cols <- grep("^logFC_", colnames(big_table), value = TRUE)
      padj_cols <- grep("^padj_", colnames(big_table), value = TRUE)
      info_cols <- c("Locustag", "Gene", "Description", "Functional_Group", padj_cols)
      
      cat("Numeric columns found:", paste(numeric_cols, collapse = ", "), "\n")
      cat("Info columns found:", paste(info_cols, collapse = ", "), "\n")
      
      # Create label mapping using graph_labels
      cat("Creating label mapping...\n")
      label_mapping <- NULL
      tryCatch({
        label_mapping <- create_label_mapping(numeric_cols, graph_labels)
        cat("Label mapping created successfully. Length:", length(label_mapping), "\n")
        cat("Label mapping names:", paste(names(label_mapping), collapse = ", "), "\n")
      }, error = function(e) {
        cat("ERROR in label mapping creation:", e$message, "\n")
        cat("Creating fallback label mapping...\n")
        label_mapping <<- setNames(numeric_cols, numeric_cols)
        cat("Fallback label mapping created with length:", length(label_mapping), "\n")
      })
      
      if (is.null(label_mapping)) {
        cat("label_mapping is NULL, creating default mapping...\n")
        label_mapping <- setNames(numeric_cols, numeric_cols)
      }
      
      # Create reverse mapping to put the short names back
      reverse_label_mapping <- create_reverse_mapping(label_mapping)
      cat("Reverse label mapping created with length:", length(reverse_label_mapping), "\n")
      
      cat("Storing values in reactive values...\n")
      values$big_table <- big_table
      values$numeric_cols <- numeric_cols
      values$info_cols <- info_cols
      values$label_mapping <- label_mapping
      values$reverse_label_mapping <- reverse_label_mapping
      values$individual_datasets <- individual_datasets
      values$data_loaded <- TRUE  
      values$gene_info <- gene_info
      
      cat("Updating UI elements...\n")
      # Update UI
      tryCatch({
        #updates names for checkboxes to show datasets in bar graph
        if (length(values$label_mapping) > 0 && length(names(values$label_mapping)) > 0) {
          checkbox_choices <- setNames(names(values$label_mapping), unname(values$label_mapping))
          #Use longname as the display label and column name as the value
          updateCheckboxGroupInput(session, "conditions_to_overlap",
                                   choices = checkbox_choices,
                                   selected = names(values$label_mapping))
          
          updateCheckboxGroupInput(session, "columns_to_graph",
                                   choices = checkbox_choices,
                                   selected = names(values$label_mapping))
          
          updateCheckboxGroupInput(session, "kegg_conditions",
                                   choices = checkbox_choices,
                                   selected = names(values$label_mapping))
          cat("Checkbox group updated successfully\n")
        } else {
          cat("WARNING: Empty label mapping, skipping checkbox update\n")
        }
        
        #Volcano plot selection
        if (length(individual_datasets) > 0) {
          dataset_choices <- setNames(names(individual_datasets), 
                                      sapply(names(individual_datasets), function(x) {
                                        # Try to find in graph_labels
                                        matching_shortname <- NULL
                                        for (i in 1:nrow(graph_labels)) {
                                          if (graph_labels$shortname[i] ==x) {
                                            matching_shortname <- graph_labels$longname[i]
                                            break
                                          }
                                        }
                                        if (!is.null(matching_shortname)) {
                                          return(matching_shortname)
                                        }
                                        return(x)
                                      }))
          updateSelectInput(session, "volcano_dataset",
                            choices = dataset_choices,
                            selected = names(individual_datasets)[1])
          cat("Volcano dataset selector updated successfully\n")
        }
        
       
        
        if (nrow(big_table) > 0) {
          updateTextInput(session, "selected_label", value = big_table$Locustag[1])
          cat("Text input updated successfully\n")
        }
      }, error = function(e) {
        cat("ERROR updating UI elements:", e$message, "\n")
        stop(paste("Failed to update UI:", e$message))
      })
      
      showNotification("Data loaded and DESeq2 analysis completed!", type = "message")
      
    }, error = function(e){
      cat("=== DETAILED ERROR INFORMATION ===\n")
      cat("Error message:", e$message, "\n")
      cat("Error class:", class(e), "\n")
      if(exists("e$call")) cat("Error call:", deparse(e$call), "\n")
      cat("Traceback:\n")
      traceback()
      cat("=====================================\n")
      
      showNotification(paste("Error loading data:", e$message), type = "error")
    })
  })
  #load venn diagram event
  observeEvent(input$load_venn, {
    cat("=== LOAD VENN BUTTON CLICKED ===\n")
    
    selected_keys <- sub("^logFC_", "", input$conditions_to_overlap)
    cat("Selected keys after stripping:", paste(selected_keys, collapse=", "), "\n")
    
    values$venn_selected_keys <- selected_keys
    values$venn_threshold_saved <- input$venn_threshold
    values$venn_ready <- TRUE
    
    cat("venn_ready set to:", values$venn_ready, "\n")
    cat("venn_selected_keys set to:", paste(values$venn_selected_keys, collapse=", "), "\n")
    cat("venn_threshold_saved set to:", values$venn_threshold_saved, "\n")
  })
  
  observeEvent(list(input$conditions_to_overlap, input$venn_threshold), {
    values$venn_ready <- FALSE
  })
  
  #a region of the venn diagram is selected
  observeEvent(input$venn_plot_up_selected, {
    req(values$venn_up_sets)
    
    selected_id <- input$venn_plot_up_selected
    cat("Venn UP region selected, id:", selected_id, "\n")
    
    if (!is.null(selected_id) && length(selected_id) > 0) {
      venn_data <- Venn(values$venn_up_sets)
      venn_processed <- process_data(venn_data)
      regions <- venn_regionlabel(venn_processed)
      
      cat("Available region ids:", paste(regions$id, collapse=", "), "\n")
      
      # Match by id directly
      matched_row <- which(regions$id == selected_id)
      cat("Matched row:", matched_row, "\n")
      
      if (length(matched_row) > 0) {
        clicked_genes <- regions$item[[matched_row]]
        cat("Genes in region:", length(clicked_genes), "\n")
        
        if (length(clicked_genes) > 0) {
          gene_table <- values$big_table %>%
            filter(Locustag %in% clicked_genes) %>%
            select(Locustag, Gene, Description, Functional_Group)
          values$venn_up_genes <- gene_table
        } else {
          values$venn_up_genes <- data.frame(Message = "No genes in this region")
        }
      }
    }
  })
  
  #same as above except for venn down
  observeEvent(input$venn_plot_down_selected, {
    req(values$venn_down_sets)
    
    selected_id <- input$venn_plot_down_selected
    cat("Venn UP region selected, id:", selected_id, "\n")
    
    if (!is.null(selected_id) && length(selected_id) > 0) {
      venn_data <- Venn(values$venn_down_sets)
      venn_processed <- process_data(venn_data)
      regions <- venn_regionlabel(venn_processed)
      
      cat("Available region ids:", paste(regions$id, collapse=", "), "\n")
      
      # Match by id directly
      matched_row <- which(regions$id == selected_id)
      cat("Matched row:", matched_row, "\n")
      
      if (length(matched_row) > 0) {
        clicked_genes <- regions$item[[matched_row]]
        cat("Genes in region:", length(clicked_genes), "\n")
        
        if (length(clicked_genes) > 0) {
          gene_table <- values$big_table %>%
            filter(Locustag %in% clicked_genes) %>%
            select(Locustag, Gene, Description, Functional_Group)
          values$venn_down_genes <- gene_table
        } else {
          values$venn_down_genes <- data.frame(Message = "No genes in this region")
        }
      }
    }
  })
  
  #run kegg analysis
  observeEvent(input$run_kegg, {
    req(values$individual_datasets, input$kegg_organism)
    values$kegg_ready <- FALSE
    
    if (input$kegg_organism == "") {
      showNotification("Please enter a KEGG organism code", type = "error")
      return()
    }
    
    selected_keys <- sub("^logFC_","",input$kegg_conditions)
    
    if (length(selected_keys) ==0) {
      showNotification("Please select at least one conditions", type = "error")
      return()
    }
    
    showNotification("Running KEGG analysis... this may take a moment", type = "message")
    
    tryCatch({
    
      results_per_condition <- map(selected_keys, function(key) {  #iterates over selected_keys and defines each individual one as the key
        cat("Running KEGG ORA for:", key, "\n")
        res_df <- values$individual_datasets[[key]]
        
        if (is.null(res_df)) {
          cat("No data found for:", key, "\n")
          return(NULL)
        }
        
        up_genes <- res_df %>%
          filter(!is.na(padj), padj < 0.05) %>%
          { if (input$kegg_direction == "up")   filter(., logFC > 1)
            else if (input$kegg_direction == "down") filter(., logFC < -1)} %>%
          pull(Locustag) %>%
          as.character()
        
        universe_genes <- res_df %>%
          filter(!is.na(padj)) %>%
          pull(Locustag) %>%
          as.character()
        
        up_genes<- normalize2kegg(up_genes, values$gene_info)
        universe_genes <- normalize2kegg(universe_genes, values$gene_info)
        
          mapped <- sum(up_genes %in% universe_genes)
          cat("  Genes recognised by KEGG:", mapped, "out of", length(up_genes), "\n")
          cat("  Unmapped genes:", paste(head(up_genes[!up_genes %in% universe_genes], 10), collapse = ", "), "\n")
        
        
        cat("  Upregulated genes found:", length(up_genes), "\n")
        cat("  Sample gene IDs being sent to KEGG:", paste(head(up_genes, 5), collapse = ", "), "\n")
        
        if (length(up_genes) <5) {
          showNotification(paste("Too few upregulated genes for condition:",key),
                           type = "warning")
          return(NULL)
        }

       
                
        result <- withCallingHandlers(
          tryCatch({
            cat("  Calling enrichKEGG...\n")
            res <- enrichKEGG(
              gene          = up_genes,
              universe      = universe_genes,
              organism      = input$kegg_organism,
              pAdjustMethod = "BH",
              pvalueCutoff  = 1,
              qvalueCutoff  = 1
            )
            cat("  enrichKEGG completed successfully\n")
            cat("  Result rows:", nrow(as.data.frame(res)), "\n")
            res
            
          }, error = function(e) {
            cat("  ERROR in enrichKEGG:", conditionMessage(e), "\n")
            print(e)
            return(NULL)
          }),
          warning = function(w) {
            cat("  WARNING in enrichKEGG:", conditionMessage(w), "\n")
            invokeRestart("muffleWarning")
          })
        
        if(is.null(result) || nrow(as.data.frame(result)) ==0) {
          cat("  No significant KEGG pathways for:", key, "\n")
          return(NULL)
        }
        
        cat("  Building output dataframe...\n")
        
        
        tryCatch({
          display_name <- key
          if (!is.null(values$label_mapping)) {
            logfc_col <- paste0("logFC_", key)
            cat("  Looking for label:", logfc_col, "\n")
            if (logfc_col %in% names(values$label_mapping)) {
              display_name <- values$label_mapping[[logfc_col]]
              cat("  Display name resolved to:", display_name, "\n")
            }
          }
          
          as.data.frame(result) %>%
            mutate(Condition = display_name, GeneCount = Count)
          
      
          
        }, error = function(e) {
          cat("  ERROR building output dataframe:", conditionMessage(e), "\n")
          print(e)
          return(NULL)
        })
        
      })
      
      
      results_per_condition <- Filter(Negate(is.null), results_per_condition)
      
      if (length(results_per_condition) == 0) {
        showNotification("No significant KEGG pathways found. Try relaxing your thresholds.", 
                         type = "warning")
        return()
      }
      
      combined_df <- bind_rows(results_per_condition)
      
      values$all_kegg_results_df <- combined_df

      
      sig_pathways <- combined_df %>% #only provides pathways significantly upregulated
        filter(p.adjust<0.05) %>%
        pull(Description) %>%
        unique()
      
      
      cat("
          
          printing unique significant pwathways...")
      print(sig_pathways)
      

      
      
      if (length(sig_pathways) == 0) {
        showNotification(
          "No pathways passed p.adjust < 0.05. Try relaxing your thresholds.",
          type = "warning"
        )
        return()
      }
      
    
      plot_df<- combined_df %>%
        filter(p.adjust<0.05, Description %in% sig_pathways) #only take significant ones with unique descriptions
      
      
      cat("
          
          
          printing dataframe for plot")
      print(plot_df)
      
      top_n <- input$kegg_top_n
      
      if(top_n>0) {
        top_ids<- plot_df %>%
          group_by(Description) %>%
          summarise(best_padj = min(p.adjust)) %>%
          slice_min(best_padj, n=top_n) %>%
          pull(Description)
        plot_df<- plot_df %>%filter(Description %in% top_ids)
      }
      
      values$kegg_results_df<- plot_df
      values$kegg_ready<- TRUE
      showNotification("KEGG analysis complete!", type = "message")
      
  }, error = function(e) {
    cat("KEGG ERROR:", e$message, "\n")
    showNotification(paste("KEGG analysis failed:", e$message), type = "error")
    
  }) 
  }) 
      
    
  
  

  
  ######Outputs
  
  
  #give an output for the loaded data
  output$data_loaded <- reactive({
    return(values$data_loaded)
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  
  #output the big data table
  output$data_table <- DT::renderDataTable({
    DT::datatable(
      values$big_table,
      extensions  = 'FixedColumns',
      options = list(
        pageLength = 10,
        lengthMenu = c(5, 10, 25, 50, 100),
        searching = TRUE,
        ordering = TRUE,
        info = TRUE,
        autoWidth = TRUE,
        scrollX = TRUE,
        fixedColumns = list(leftColumns = 2),
        scrollCollapse = TRUE),
      rownames = FALSE,
      class = 'cell-border stripe hover',
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: center; font-size: 1.2em; font-weight: bold;'
      )
    )
  })
  
  
  
  # Output differential genes table (2+/- fold change)
  output$diff_genes_table <- reactable::renderReactable({
    
    diff_genes_data <- diff_genes_reactive()
    
    if (is.null(diff_genes_data) || nrow(diff_genes_data) == 0) {
      return(reactable(data.frame(Message = "No differentially expressed genes found with current thresholds")))
    }
    
    # Store in reactive values for export
    values$diff_genes_table <- diff_genes_data
    
    reactable(diff_genes_data, 
              pagination = TRUE,
              showPageSizeOptions = TRUE,
              searchable = TRUE,
              columns = list(
                Locustag = colDef(name = "Locustag", width = 110),
                Gene = colDef(width = 80),
                Description = colDef(width = 275),
                Functional_Group = colDef(name = "Functional Groups", width = 275),
                logFC = colDef(name = "Log2 Fold Change", format = colFormat(digits = 3)),
                padj = colDef(name = "P-Value", format = colFormat(digits = 6)),
                diffexpressed = colDef(name = "Expression", width = 100)
              ),
              theme = reactableTheme(headerStyle = list(background = "lightgray")))
  })
  
  
  
  # Volcano plot output
  output$volcano_plot <- renderGirafe({
    girafe(
      ggobj = create_volcano_plot(),
      width_svg = 12,
      height_svg = 6,
      options = list(
        opts_hover(css = "fill-opacity:1; stroke:black; stroke-width:1.5px;"),
        opts_tooltip(
          css = "background-color:white; color:black; border:1px solid gray;
               border-radius:5px; padding:8px; font-size:13px;",
          use_fill = FALSE
        ),
        opts_zoom(max = 5)   # bonus: lets users scroll-zoom into dense regions
      )
    )
  })
  
  #barplot graphing
  output$bar_plot <- renderPlot({
    create_bar_plot()
  })
  
  #venn diagram
  output$venn_plot_up <-renderGirafe({
    
    req(values$data_loaded)   
    req(values$venn_ready)
    req(values$venn_selected_keys)
    
    tryCatch({
      cat("=== ATTEMPTING TO CREATE VENN UP ===\n")
      plot <- create_venn_up(values$venn_selected_keys, values$venn_threshold_saved)
      
      girafe(ggobj = plot,
             options = list(
               opts_selection(type = "single", css = girafe_css(css = "fill:white",
                                                                text="fill:black; stroke:black; font-weight:bold")),
               opts_hover(css = "fill-opacity:0.7; cursor:pointer;")
             ))
    }, error = function(e) {
      cat("ERROR in venn_plot_up:", e$message, "\n")
      print(traceback())
    })
  })
  
  output$venn_plot_down <-renderGirafe({
    
    req(values$data_loaded)   
    req(values$venn_ready)
    req(values$venn_selected_keys)
    
    tryCatch({
      cat("=== ATTEMPTING TO CREATE VENN DOWN ===\n")
      plot <- create_venn_down(values$venn_selected_keys, values$venn_threshold_saved)
      
      girafe(ggobj = plot,
             options = list(
               opts_selection(type = "single", css = girafe_css(css = "fill:white",
                                                                text="fill:black; stroke:black; font-weight:bold")),
               opts_hover(css = "fill-opacity:0.7; cursor:pointer;")
             ))
    }, error = function(e) {
      cat("ERROR in venn_plot_down:", e$message, "\n")
      print(traceback())
    })
  })
  
  
  outputOptions(output, "venn_plot_up", suspendWhenHidden = FALSE)
  outputOptions(output, "venn_plot_down", suspendWhenHidden = FALSE)
  
  
  output$venn_up_table <- reactable::renderReactable({
    req(values$venn_up_genes)
    reactable(values$venn_up_genes,
              searchable = TRUE,
              columns = list(
                Locustag = colDef(name = "Locustag", width = 110),
                Gene = colDef(width = 80),
                Description = colDef(width = 275),
                Functional_Group = colDef(name = "Functional Group", width = 200)
              ),
              theme = reactableTheme(headerStyle = list(background = "lightgray")))
  })
  
  output$venn_down_table <- reactable::renderReactable({
    req(values$venn_down_genes)
    reactable(values$venn_down_genes,
              searchable = TRUE,
              columns = list(
                Locustag = colDef(name = "Locustag", width = 110),
                Gene = colDef(width = 80),
                Description = colDef(width = 275),
                Functional_Group = colDef(name = "Functional Group", width = 200)
              ),
              theme = reactableTheme(headerStyle = list(background = "lightgray")))
  })
  
  
  #KEGG outputs
  
  output$kegg_ready <- reactive({ values$kegg_ready})
  outputOptions(output, "kegg_ready", suspendWhenHidden = FALSE)
  
  
  output$kegg_bubble_plot <- renderGirafe({
    req(values$kegg_results_df)
    
    plot_df <- values$kegg_results_df
    
    plot_df <- plot_df %>%
      mutate(
        tooltip_text = paste0(
          "<b>", Description,"</b>\n",
          "Condition: ", Condition, "\n",
          "Genes in pathway: ", GeneCount, "\n",
          "Adj. p-value: ", formatC(p.adjust,format = "e", digits = 3), "\n\n",
          "Genes:\n",
          gsub("/", "\n", geneID)
        )
      )
    
    
    p<- ggplot(plot_df, 
           aes(
             x= Condition,
             y = Description,
             size = GeneCount,
             fill = p.adjust))+
      geom_point_interactive(aes(
        tooltip = tooltip_text,
        data_id = paste(Condition, Description, sep = "_")),
        shape = 21,
                 color = "grey30", 
                 alpha = 0.85)+
      scale_size_continuous(name = "Upregulated\nGene Count",
                            range = c(3,14))+
      scale_fill_gradient(name = "Adjusted\np-value",
                          low = "red",
                          high = "blue",
                          trans = "log10")+
      labs(title = paste("Significantly", input$kegg_direction ,"regulated KEGG Pathways[", toupper(input$kegg_organism), "]"),
           x = "Condition",
           y = "KEGG Pathway",
           caption = "Bubble size = # of upregulated genes; Color = adjusted p-value")+
      theme_classic()+
      theme(
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        panel.grid.major = element_line(color = "grey90"),
        plot.title = element_text(face = "bold",
                                  hjust= 0.5)
      )
    
    girafe(
      ggobj = p,
      width_svg = 10,
      height_svg = 7,
      options = list(
        opts_hover(css = "fill-opacity:1; 
                   stroke:black; 
                   stroke-width:2px; 
                   cursor:pointer;"),
        opts_tooltip(
          css = "background-color:white;
          color:black;
          border:1px solid gray;
          border-radius:5px;
          padding:10px;
          font-size:12px;
          white-space:pre;",
          use_fill = FALSE,
          offx = 10,
          offy = 10
        ),
        opts_zoom(max = 5)
      )
    )
    
  })
  
  output$kegg_results_table  <-reactable::renderReactable({
    req(values$all_kegg_results_df)
    
    cat("Columns in kegg_results_df:", paste(colnames(values$all_kegg_results_df), collapse = ", "), "\n")
    
    desired_cols <- c("Condition", "Description", "GeneCount", "GeneRatio", 
                      "BgRatio", "p.adjust", "geneID")
    
    available_cols <- desired_cols[desired_cols %in% colnames(values$all_kegg_results_df)]
    
    cat("Available cols selected:", paste(available_cols, collapse = ", "), "\n")
    
    display_df <- values$all_kegg_results_df %>% select(all_of(available_cols))
    
    all_col_defs <- list(
      Condition   = colDef(width = 150),
      Description = colDef(width = 250),
      GeneCount   = colDef(name = "Gene Count", width = 100),
      GeneRatio   = colDef(width = 100),
      BgRatio     = colDef(width = 100),
      p.adjust    = colDef(name = "Adj. P-value", format = colFormat(digits = 5)),
      geneID      = colDef(name = "Gene IDs", width = 300)
    )
    
    active_col_defs <- all_col_defs[names(all_col_defs) %in% available_cols]
    
    reactable(
      display_df,
      searchable         = TRUE,
      showPageSizeOptions = TRUE,
      columns            = active_col_defs,
      theme              = reactableTheme(headerStyle = list(background = "lightgray"))
    )
  })
  
  
  output$kegg_plot_download <- downloadHandler(
    filename = function() paste0("keggbubble_", input$kegg_organism, "_", Sys.Date(),".pdf"),
    content = function(file) {
      req(values$kegg_results_df)
      plot_df <- values$kegg_results_df
      p <- ggplot(plot_df, 
                  aes(x = Condition, 
                      y = Description, 
                      size = GeneCount, 
                      fill = p.adjust)) +
        geom_point(shape = 21, 
                   colour = "grey30", 
                   alpha = 0.85) +
        scale_size_continuous(name = "Upregulated\nGene Count", 
                              range = c(3, 14)) +
        scale_fill_gradient(name = "Adjusted\np-value", 
                            low = "red", 
                            high = "blue", 
                            trans = "log10") +
        labs(title = paste("Significantly", input$kegg_direction ,"regulated KEGG Pathways[", toupper(input$kegg_organism), "]"),
             x = "Condition", 
             y = "KEGG Pathway") +
        theme_bw(base_size = 12) +
        theme(axis.text.y = element_text(size = 9),
              axis.text.x = element_text(angle = 30, 
                                         hjust = 1),
              plot.title  = element_text(face = "bold", 
                                         hjust = 0.5))
      
      ggsave(file, plot = p, device = "pdf", width = 10, height = 8)
    }
  )
  
  output$kegg_csv_download <- downloadHandler(
    filename = function() paste0("kegg_results_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(values$kegg_results_df, file, row.names = FALSE)
    }
  )
  
  
  # Display additional information table
  output$info_table <- renderTable({
    
    #input the selected label
    req(input$selected_label)
    
    #if the label is not in the big table return this message
    if (!input$selected_label %in% values$big_table$Locustag) {
      return(data.frame(Message = "Label not found"))
    }
    
    #define the selected row and into data to show
    selected_row <- values$big_table[values$big_table$Locustag == input$selected_label, ]
    info_data <- selected_row[, values$info_cols, drop = FALSE]
    
    #info display
    info_display <- data.frame(
      Attribute = names(info_data),
      Value = as.character(info_data[1, ])
    )
    
    #relabel things to look nicer
    info_display$Attribute <- sapply(info_display$Attribute, function(attr_name) {
      if (grepl("padj_", attr_name)) {
        experiment_shortname <- sub("padj_", "", attr_name)
        logfc_col <- paste0("logFC_", experiment_shortname)
        if (logfc_col %in% names(values$label_mapping)) {
          return(paste("Adjusted P-Value:", values$label_mapping[[logfc_col]]))
        }
      }
      return(attr_name)
    })
    
    
    
    #conditional formatting for significance
    for (i in 1:nrow(info_display)) {
      if (grepl("Adjusted P-Value:", info_display$Attribute[i])) {
        raw_val <- info_display$Value[i]
        pval <- suppressWarnings(as.numeric(raw_val))
        
        if (is.na(pval)) {
          experiment_label <- sub("Adjusted P-Value: ", "", info_display$Attribute[i])
          shortname <-values$reverse_label_mapping[[experiment_label]]
          shortname <- sub("^logFC_", "", shortname)
          reason_col <- paste0("reason_", shortname)
          
          reason_text <- "Inconclusive (reason unavailable)"
          if (reason_col %in% names(selected_row)) {
            r<- selected_row[[reason_col]][1]
            if (!is.na(r)) reason_text <- r
          }
          
          info_display$Value[i] <- paste0(
            '<span style = "background-color: #fff3cd; font-weight: bold; padding: 2px 4px; border-radius:3px;">',
            "Inconclusive: ", reason_text,
            '</span>'
          )
        } else if (pval < 0.05) {
          info_display$Value[i] <- paste0('<span style="background-color: #ffcccc; font-weight: bold; padding: 2px 4px; border-radius: 3px;">', 
                                          info_display$Value[i], '</span>')
        }
      }
    }
    info_display
  }, striped = TRUE, bordered = TRUE, sanitize.text.function = function(x) x)
  
  #filtering summary
  output$filtering_summary <- renderText({
    req(values$filtering_log)
    paste("Analysis Log:", paste(values$filtering_log, collapse = " | "))
  })
  
  ###Download handlers
  
  output$export_csv <- downloadHandler(
    filename = function() {
      dataset_name <- if(!is.null(input$volcano_dataset)) input$volcano_dataset else "unknown"
      paste("differential_expression_", dataset_name, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Get the current differential genes data
      export_data <- diff_genes_reactive()
      write.csv(export_data, file, row.names = FALSE)
    }
  )
  
  output$bigexport_csv <- downloadHandler(
    filename = function() {
      paste("AllData_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(values$big_table, file, row.names= FALSE)
    }
  )
  
  output$selected_venn_up_csv <- downloadHandler(
    filename = function() {
      paste("Selected upreguated genes", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(values$venn_up_genes, file, row.names= FALSE)
    }
  )
  
  output$selected_venn_down_csv <- downloadHandler(
    filename = function() {
      paste("Selected downreguated genes", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(values$venn_down_genes, file, row.names= FALSE)
    }
  )
  
  output$bars_jpg <- downloadHandler(
    filename = function() {
      dataset_name <- if(!is.null(input$selected_label)) input$selected_label else "unknown"
      paste("bar_plot_", dataset_name, "_", Sys.Date(), ".jpg", sep = "")
    },
    content = function(file) {
      plot<- create_bar_plot()
      if (!is.null(plot)) {
        ggsave(file, plot = plot, device = "jpeg", width = 10, height = 8, units = "in", dpi = 300)
      }
    }
  )
  
  
  
}


#the shiny app!
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE, height = 1200, width = 1000))

