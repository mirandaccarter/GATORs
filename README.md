# GATORs

Thank you for choosing to use GATORs for your RNAseq analysis! GATORs is an all in one tool created to generate log2 Fold Change from RNAseq derived gene counts, perform DeSeq2 powered differential expression analysis, and produce publication-quality visuals, including volcano plots, Venn Diagrams, bar graphs, and perform functional analysis through KEGG.

# How to use GATORs
## Prior to download
GATORs is an R-based shiny app. Prior to downloading the GATORs script ensure you have installed both R and Rstudio on your device 
#### R: https://cran.r-project.org/ 
#### Rstudio: https://posit.co/products/open-source/rstudio 

## After download
Once you have downloaded the script, open the file in Rstudio. Upon running GATORs for the first time there will be many packages that need to be downloaded, which GATORs will do the first time it runs. Every time afterwards GATORs will continue to check for updates and routinely update where needed. 
#####
To begin first run, click the button in the top right hand corner of the window containing the script labeled "Run App" with a green arrow next to it. If this is your first time using GATORs it will take several minutes to download and install necessary packages. This may take several minutes. If you are asked to update packages during this process indicate yes.

## File preparation
Three sample datasets have been provided here to demonstrate the correct file set up. 
### Main requirements
1. Ensure all column names are formatted without underscores until the separator between name and replicate number (ie drug1.A_1, where "_1" indicates the replicate number)
2. GATORs was developed with the intention to compare multiple datasets to each other, in order to do so each spreadsheet included in the data file should have ONE set of control replicates which all the experimental conditions will be compared to. If you would like to use multiple different datasets with different controls, include multiple spreadsheets within the data folder, but all controls must be labeled the same
3. All files MUST be saved as a .csv, the only exception is the library file which may be saved as a .tsv
4. Library set up must include the following columns: Locustag (containing gene id), Gene (containing gene name), Description (containing a brief description of the gene), Functional_Group (containing relevant functional group info), and kegg_gene (containing the gene locus as found in KEGG).
6. The labels file contains two columns, shortname and longname. Shortname is the column  indicating the label as shown on the data. Long name is the corresponding label to be used for display purposes.

## Get started!
Once GATORs is ready, a page will appear asking for three pathway inputs and a control identifier. Paste in the pathway to the folder containing data spreadsheets, the pathway to the gene library, and pathway to the labels file. 
- When providing the path to the files, ensure that all of the backwards slashes are changed to forward slashes
#####
You will then need to provide a control identifier. This can be whatever you would like it to be, but it must be the same across all spreadsheets. Click Run.
#####
After a few seconds GATORs will load a new page with several tabs. 
#### Comparative Analysis
The first section to appear includes two tabs. First a bar graph with the log2 Fold Change of a chosen gene under each experimental condition provided. CHange the gene displayed by typing in the locustag for a different gene in the box above. The second tab includes additional information on the gene selected and the significance of the chosen gene's regulation differences, with significant values highlighted.
#### Spreadsheet of all values 
#### Volcano Plots
#### Venn Diagrams
#### KEGG Analysis


