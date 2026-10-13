# GATORs

Thank you for choosing to use GATORs for your RNAseq analysis! GATORs is an all in one tool created to generate log2 Fold Change from RNAseq derived gene counts, perform DeSeq2 powered differential expression analysis, and produce publication-quality visuals, including volcano plots, Venn Diagrams, bar graphs, and perform functional analysis through KEGG.

# How to use GATORs
## Prior to download
GATORs is an R-based shiny app. Prior to downloading the GATORs script ensure you have installed both R and Rstudio on your device 
#### R: https://cran.r-project.org/ 
#### Rstudio: https://posit.co/products/open-source/rstudio 

## After download
Once you have downloaded the script (GATORS 2026.10.13.R), open the file in Rstudio. Upon running GATORs for the first time there will be many packages that need to be downloaded, which GATORs will do the first time it runs. Every time afterwards GATORs will continue to check for updates and routinely update where needed. 
#####
To begin first run, click the button in the top right hand corner of the window containing the script labeled "Run App" with a green arrow next to it. If this is your first time using GATORs it will take several minutes to download and install necessary packages. This may take several minutes. If you are asked to update packages during this process indicate yes and that all packages should be updated.

<img width="407" height="289" alt="image" src="https://github.com/user-attachments/assets/16944a6c-8f49-4239-8ae5-4aa1c65a5c3c" />


## File preparation
Three sample datasets have been provided here to demonstrate the correct file set up. 

GATORs has been developed with the intent to compare trends in RNAseq data, and will accept multiple spreadsheets within the data folder. An example of this set up has been provided in the _S. mutans_ zinc transcriptome.

For each spreadsheet provided to GATORs in the data folder, there must only be one control to which every distinct experimental condition will be compared to. If you wish to compare the results of different experiments where different controls were used, you must ensure that the control condition is indicated in the same way across all spreadsheets. **Comparisons across experiments should NOT be used to make conclusions without additional experimental validation, and should ONLY be useful for observation of trends.**

In the case of the provided example datasets, the Control condition in each spreadsheet indicates the following: 
- _S. mutans_ 4 mM Zn treatment: Control = untreated sample taken prior to zinc exposure for 15 or 90 minutes
- _S. mutans_ HP29 + Zn: Control = untreated sample taken after the 30 minute exposure period
- _C. albicans_ fluconazole treatment: Nodrug-rna serves as the control and is an untreated sample collected at the same time as flucoonazole treated samples
- Macrophage outer membrane vesicle treatment: Control = Macrophages not treated with either wild type or Δspd OMVs


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

You will then need to provide a control identifier. This can be whatever you would like it to be, but it must be the same across all spreadsheets. Click Run.

After a few seconds GATORs will load a new page. For demonstration purposes, the _S. mutans_ zinc transcriptome will be pictured below. 

Immediately to the right of the Load and Run DESeq2 button a notifier of zero count genes will appear. As part of the innate filtering with DESeq2 genes with a zero-count will be removed. If this is applicable to your dataset you will find the number of excluded genes here.

<img width="850" height="113" alt="image" src="https://github.com/user-attachments/assets/91d16e0c-aa62-4455-a472-6b7f2f962ae9" />

#### Comparative Analysis
The first section to appear includes two tabs. First a bar graph with the log2 Fold Change of a chosen gene under each experimental condition provided. Change the gene displayed by typing in the locustag for a different gene in the box above. The second tab includes additional information on the gene selected and the significance of the chosen gene's regulation differences, with significant values highlighted. This section will remain visible regardless of which lower tabs are chosen.

##### Alterable bar graph displaying gene regulation for each experimental condition
<img width="829" height="424" alt="image" src="https://github.com/user-attachments/assets/0aaeb831-8b72-4775-93ea-3a560ef07828" />

##### Additional information tab displaying significance of gene regulation and relevant genetic info
<img width="349" height="134" alt="image" src="https://github.com/user-attachments/assets/44dca445-8f23-4eaf-8a8f-1d0a5e582ebe" />

### Lower tabs
#### Spreadsheet of all values 
The first tab shows a downloadable spreadsheet of all resulting data. This spreadsheet includes the log2FC, p-value, adjusted p-value, and all information data for each gene. Additionally, for those genes which do not display an adjusted p-value, the reason why is indicated in the following column. This spreadsheet contains information for all of the datasets included in the analysis, and is  downloadable as a .csv file.

##### Searchable Spreadsheet of all data produced after differential expression analysis of all files 
<img width="798" height="420" alt="image" src="https://github.com/user-attachments/assets/47efa1b9-0428-4e8f-b76e-e6d7629528e6" />

#### Volcano Plots
The second tab provides ready-made volcano plots, allowing the user to select the dataset to display, cut-off values for p-value and log2FC, number of top and bottom genes to be labeled, and if you would like the p-value or adjusted p-value to be used. 
##### Adjustable volcano plot of 4 mM Zn treatment after 15 minutes 
<img width="785" height="497" alt="image" src="https://github.com/user-attachments/assets/4018a98c-f63c-4cc8-8742-a257e45271bf" />

Included in this tab is a table of the differentially regulated genes in the dataset chosen for the volcano plot. This is a downloadable table displaying only genes above or below a certain threshold of log2FC and p-value.

##### Searchable table of differentially regulated genes for dataset shown in volcano plot above
<img width="800" height="481" alt="image" src="https://github.com/user-attachments/assets/3d389edd-60ad-47f8-b615-07dcc9f4bff9" />

#### Venn Diagrams
The third tab provides a tool to build Venn Diagrams of gene overlap between up to six different datasets. After selecting datasets to compare and a log2FC threshold for genes to include, GATORs will build two Venn diagrams. One for comparison of upregulated genes and one for comparison of downregulated genes.
If a given dataset does not have any genes that meet the threshold requirements set by the user, no section will be displayed for that dataset (as seen below where HP29 treatment alone does not produce sufficient genes downregulated past 2 log2FC).

##### Venn Diagram Builder
<img width="798" height="386" alt="image" src="https://github.com/user-attachments/assets/ad116c5c-07dc-4cf2-95e2-50280621b1ef" />

Each overlapping section will have the number of overlapping genes, which, if clicked on, will produce a downloadable spreadsheet of all the genes included in that section.

##### Venn Diagram overlap spreadsheet
<img width="401" height="385" alt="image" src="https://github.com/user-attachments/assets/a3fb28c9-996b-405e-bdd2-ee4e19c51046" />

#### KEGG Analysis
The final tab provides a tool analyze KEGG pathways associated with your results. Powered by the R package clusterProfiler, users can select which datasets they would like to include, choose upregulated or downregulated genes, and GATORs will produce a bubble chart to show which KEGG pathways are significantly represented. Additionally a downloadable spreadsheet is provided with specific pathway identifiers, Gene Ratios, and genes which are differentially regulated in each pathway.

##### Bubble chart
<img width="803" height="407" alt="image" src="https://github.com/user-attachments/assets/21b3d571-a41c-4a00-9f5e-cfe8f7bb7c73" />

##### Gene enrichment table
<img width="551" height="418" alt="image" src="https://github.com/user-attachments/assets/cb86746c-cf0c-44f7-b177-3477e97fa065" />

