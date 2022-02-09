############################################################
### README - Datasets and scripts supporting the manuscript
### 
### "Genetic analysis reveals three novel QTLs underpinning 
### a butterfly egg-induced hypersensitive response-like 
### cell death in Brassica rapa" 
### 
###	available at ResearchSquare - DOI: 10.21203/rs.3.rs-1107347/v1
###
### Author:   Niccolo` Bassetti - niccolo.bassetti@protonmail.com
### 
### Co-authors: Lotte Caarls, Gabriella Bukovinszkine’Kiss, Mohamed El-Soda, Jeroen van Veen, 
###				Klaas Bouwmeester, Bas J. Zwaan, M. Eric Schranz, Guusje Bonnema, Nina E. Fatouros
###
### Principal Investigator:   Nina Fatouros - nina.fatouros@wur.nl


1) This .zip archive contains one folder for each figure of the manuscript. 
The folders contain datasets and scripts used for data analysis and/or plotting of data. 
Final versions of figures have been edited in PowerPoint and are not included in this archive.


2) The raw data used for data analysis are attached to the manuscript in as tables within the file "Additional File 3 - Original datasets.xlsx":
	Table S11:  Phenotypic data of the 56 B. rapa accessions tested for the germplasm screening.
	Table S12:	Phenotypic data of B. rapa homozygous lines re-evaluated for egg- and egg wash-induced cell death.
	Table S13:	Phenotypic data of parents and L58 x R-o-18 RIL population across three repeated QTL experiments.
	Table S14:	Marker data used to construct the L58 x R-o-18 RIL population genetic map.
	Table S15:	Phenotypic data of twelve selected RIL lines that were re-evaluated for egg-induced cell death.

NOTE: The raw data are provided also as .txt files within the folders of the figures.


3) These folders contains .R scripts to replicate data analysis and plot.
   Each folder contains .md5 files to check the integrity of datasets after download (MD5 hashes).

	Figure_1 - Phenotypic data of ten selected B. rapa accessions: 
		Data analysis & plot.
		This figures requires raw data from "dataset_Table_S12.txt".
	Figure_2 - Phenotypic data of RIL population L58 x R-o-18:
		Data analysis & plot.
		This figures requires raw data from "dataset_Table_S13.txt".
	Figure_3 - QTL mapping of HR cell death in RIL population: 
		Data analysis & plot.
		This figures requires raw data from "dataset_Table_S13.txt" and "dataset_Table_S14.txt".
	Figure_4 - Validation of QTL effects in twelve select RIL lines: 
		Data analysis and plot.
		This figures requires raw data from "dataset_Table_S15.txt".
	Sup_Figure_S1 - Germplasm screening of 56 B. rapa accessions:
		Data analysis.
		This figures requires raw data from "dataset_Table_S11.txt".
	Sup_Figure_S3 - Image-based phenotyping of cell death
	Sup_Figure_S5 - QTL LOD scores across all 10 B. rapa chromosomes: plot
	Sup_Figure_S6 - Heatmap of LOD scores from two-QTL model analysis: plot
	Sup_Figure_S7 - Phenotypic data of twelve selected RIL lines for validation of QTL effects: plot

4) These folders contain .pdf of figures that were generated elsewhere:
	Sup_Figure_S2 - Image-based phenotyping protocol:
		Image was assembled in PowerPoint.
	Sup_Figure_S4 - Genetic map of L58 x R-o-18:
		Image was generated with JoinMap v4.0
	Sup_Figure_S8 - Syntheny between QTLs Pbc1,2,3 and A. thaliana:
		Image was assembled in PowerPoint.
		Analysis was done in CoGe (https://genomevolution.org/coge/) and can be reproduced with this link:
		https://genomevolution.org/r/1kheq

5) Folder "z_z_Image_based_phenotyping" contains FiJi/Image script to carry out image analysis of HR-like cell death