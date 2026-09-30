# 1. Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr) # For adding the correlation coefficient (R) and p-value
library(magrittr)
library(ggVennDiagram)
        
        # 2. Prepare the data
        # Assuming 'edger_results' is the dataframe from the company 
        # and 'deseq2_results' is the dataframe from your tool
        # Both should have a column called 'Locustag' and 'logFC'
        
setwd("/home/mirandacarter/HOME/RNAseq tool/validation/")

    deseq2_results <- read.csv("/home/mirandacarter/HOME/RNAseq tool/validation/candida results/deseq candida results.csv")
    edgeR_results <- read.csv("/home/mirandacarter/HOME/RNAseq tool/validation/candida results/candida study results.csv")

        comparison_df <- inner_join(
          # Rename columns to avoid conflict during the join
          edgeR_results %>% dplyr::select(Locustag, edger_lfc = logFC), 
          deseq2_results %>% dplyr::select(Locustag, deseq_lfc = logFC), 
          by = "Locustag"
        )
        


        
        
        # 3. Create the plot
        correlation_plot <- ggplot(comparison_df, aes(x = edger_lfc, y = deseq_lfc)) +
          # Add a density contour to show where most genes cluster (optional but looks great)
          geom_bin2d() + 
          scale_fill_gradient(low = "lightblue", high = "darkblue") +
          # Add a linear regression line (the "Trend" line)
          geom_smooth(method = "lm", color = "red", se = TRUE) +
          # Add the Correlation Coefficient (R) and p-value automatically
          stat_cor(aes(label = paste(after_stat(rr.label),
                                     after_stat(p.label), 
                                     sep = "~`,`~")),
                   label.x = -8, 
                   label.y = 9,
                   size = 6) + 
          labs(
            title = "Comparison of Log2FoldChange:\nPublished vs. GATORs",
            x = "Published Log2FC",
            y = "GATORs Log2FC",
            fill = "Gene Density"
          ) +
          theme_classic() +
          theme(
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.title = element_text(face = "bold",
                                      size = 14),
            axis.text = element_text(face = "bold",
                                     size = 10),
            legend.text = element_text( size = 10,
                                        face = "bold"),
            legend.title = element_text(size = 15,
                                        face = "bold",
                                        hjust = 0.5))
        
        # Display the plot
        print(correlation_plot)
        
        ggsave("correlation_plot.png", plot = correlation_plot, width = 6, height = 6, dpi = 600)
        
       #get different genes 
        my_genes <- deseq2_results %>% 
          filter(Pvalue < 0.05 & abs(logFC) > 1) %>% 
          pull(Locustag)
        
        published_genes <- edgeR_results %>% 
          filter(pvalue < 0.05 & abs(logFC) > 1) %>% 
          pull(Locustag)
        
        # 1. Genes found in BOTH (The Overlap)
        both_genes <- intersect(my_genes, published_genes)
        
        # 2. Genes found ONLY in your tool (Unique to GATORs)
        unique_to_me <- setdiff(my_genes, published_genes)
        
        # 3. Genes found ONLY in the published data
        unique_to_published <- setdiff(published_genes, my_genes)
        
        
        # Combine the gene lists into a named list
        venn_list <- list(
          GATORs = my_genes,
          Published = published_genes
        )
        
        # Create the Venn Diagram
        venn_plot <- ggVennDiagram(venn_list,
                                   label_size = 8,
                                   category.names = c("","")) +
          coord_flip()+
          annotate("text", 
                   label = c("GATORs", "Published"),
                   x = c(4.5, 4.5), 
                   y = c(-0.1,4.1),
                   size = 7)+
          scale_fill_gradient(low = "white", high = "white") +
          labs(title = "Differentially Expressed Genes",
               fill = "Gene Number") +
          theme(legend.position = "inside",
                legend.position.inside = c(0,-4),
                legend.text = element_text(size = 10,
                                           face = "bold"),
                legend.title = element_text(size = 12,
                                            face = "bold"),
                plot.title = element_text(face = "bold",
                                          size = 18,
                                          hjust = 0.5,
                                          margin = margin( t = 20))
                
                )
        
        print(venn_plot)
        
        # Save the plot for your paper
        ggsave("venn_validation.png", plot = venn_plot, width = 6, height = 5, dpi = 300)
        
        
        