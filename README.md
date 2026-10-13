# GATORs

Thank you for choosing to use GATORs for your RNAseq analysis! GATORs is an all in one tool created to generate log2 Fold Change from RNAseq derived gene counts, perform DeSeq2 powered differential expression analysis, and produce publication-quality visuals, including volcano plots, Venn Diagrams, bar graphs, and perform functional analysis through KEGG.

# How to use GATORs
## Prior to download
GATORs is an R-based shiny app. Prior to downloading the GATORs script ensure you have installed both R and Rstudio on your device 
#### R: https://cran.r-project.org/ 
#### Rstudio: https://posit.co/products/open-source/rstudio 

## After download
Once you have downloaded the script (GATORS 2026.09.30.R), open the file in Rstudio. Upon running GATORs for the first time there will be many packages that need to be downloaded, which GATORs will do the first time it runs. Every time afterwards GATORs will continue to check for updates and routinely update where needed. 
#####
To begin first run, click the button in the top right hand corner of the window containing the script labeled "Run App" with a green arrow next to it. If this is your first time using GATORs it will take several minutes to download and install necessary packages. This may take several minutes. If you are asked to update packages during this process indicate yes and that all packages should be updated.

<img width="407" height="289" alt="image" src="https://github.com/user-attachments/assets/16944a6c-8f49-4239-8ae5-4aa1c65a5c3c" />


## File preparation
Three sample datasets have been provided here to demonstrate the correct file set up. 

GATORs has been developed with the intent to compare trends in RNAseq data, and will accept multiple spreadsheets within the data folder. An example of this set up has been provided in the _S. mutans_ zinc transcriptome.

### Main requirements
1. Ensure all column names are formatted without underscores until the separator between name and replicate number (ie drug1.A_1, where "_1" indicates the replicate number)
2. GATORs was developed with the intention to compare multiple datasets to each other, in order to do so each spreadsheet included in the data file should have ONE set of control replicates which all the experimental conditions will be compared to. If you would like to use multiple different datasets with different controls, include multiple spreadsheets within the data folder, but all controls must be labeled the same
3. All files MUST be saved as a .csv, the only exception is the library file which may be saved as a .tsv
4. Library set up must include the following columns: Locustag (containing gene id), Gene (containing gene name), Description (containing a brief description of the gene), Functional_Group (containing relevant functional group info), and kegg_gene (containing the gene locus as found in KEGG).
6. The labels file contains two columns, shortname and longname. Shortname is the column  indicating the label as shown on the data. Long name is the corresponding label to be used for display purposes.

## Get started!
Once GATORs is ready, a page will appear asking for three pathway inputs and a control identifier. Paste in the pathway to the folder containing data spreadsheets, the pathway to the gene library, and pathway to the labels file. 
<img width="936" height="249" alt="image" src="https://github.com/user-attachments/assets/021d315e-6297-4630-9728-47b1f5f5bff2" />
- When providing the path to the files, ensure that all of the backwards slashes are changed to forward slashes
#####
You will then need to provide a control identifier. This can be whatever you would like it to be, but it must be the same across all spreadsheets. Click Run.
#####
After a few seconds GATORs will load a new page.
#### Comparative Analysis
The first section to appear includes two tabs. First a bar graph with the log2 Fold Change of a chosen gene under each experimental condition provided. Change the gene displayed by typing in the locustag for a different gene in the box above. The second tab includes additional information on the gene selected and the significance of the chosen gene's regulation differences, with significant values highlighted. This section will remain visible regardless of which lower tabs are chosen.

### Lower tabs
#### Spreadsheet of all values 
The first tab shows a downloadable spreadsheet of all resulting data. This spreadsheet includes the log2FC, p-value, adjusted p-value, and all information data for each gene. Additionally, for those genes which do not display an adjusted p-value, the reason why is indicated in the following column. This spreadsheet contains information for all of the datasets included in the analysis, and is  downloadable as a .csv file.
#### Volcano Plots
The second tab provides ready-made volcano plots, allowing the user to select the dataset to display, cut-off values for p-value and log2FC, number of top and bottom genes to be labeled, and if you would like the p-value or adjusted p-value to be used. 

Included in this tab is a table of the differentially regulated genes in the dataset chosen for the volcano plot. This is a downloadable table displaying only genes above or below a certain threshold of log2FC and p-value.
#### Venn Diagrams
The third tab provides a tool to build Venn Diagrams of gene overlap between up to six different datasets. After selecting datasets to compare and a log2FC threshold for genes to include, GATORs will build two Venn diagrams. One for comparison of upregulated genes and one for comparison of downregulated genes. Each overlapping section will have the number of overlapping genes, which, if clicked on, will produce a downloadable spreadsheet of all the genes included in that section.
#### KEGG Analysis
The final tab provides a tool analyze KEGG pathways associated with your results. Powered by the R package clusterProfiler, users can select which datasets they would like to include, choose upregulated or downregulated genes, and GATORs will produce a bubble chart to show which KEGG pathways are significantly represented. Additionally a downloadable spreadsheet is provided with specific pathway identifiers, Gene Ratios, and genes which are differentially regulated in each pathway.

