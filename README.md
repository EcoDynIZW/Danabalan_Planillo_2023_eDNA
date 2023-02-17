# Danabalan_Planillo_2013_eDNA
Comparison of mosquito and fly derived DNA as a tool for sampling vertebrate biodiversity in suburban forests in Berlin, Germany


##Abstract
The use of invertebrate-derived DNA (iDNA) is a promising non-invasive tool to monitor wildlife. While most studies have been carried out in dense tropical and sub-tropical forests and have focused on the use of a single category of invertebrates, this study compares the use of flies and mosquitoes derived DNA to assess vertebrate diversity in semi-urban environments. We conducted our sampling in four different forest plots in Berlin, Germany. Pools of flies and non-bloodfed mosquitoes were metabarcoded using 108 bp vertebrate-specific 12S rRNA (12S-V5) and 94 bp mammal-specific 16S rRNA (16Smam) mitochondrial markers, and individual bloodfed mosquitoes were sequenced using the 340 bp vertebrate-specific 12S rRNA fragment (Mam-12S-340). Most sequencing was only successful for mammal species. From the fly pools, we detected 10 mammal species using 16Smam, and six species using 12S-V5. From the non-bloodfed mosquito pools, we only amplified putative contaminant DNA, indicating that mosquito females without visual signs of a blood meal carry no traces of vertebrate DNA. Finally, in the bloodfed mosquitoes we identified four mammal species. We did not find significant differences in the proportion of mammal species detected regarding the total available number of species between sampling localities. Fly samples were easier to obtain and more abundant over the sampled localities compared to mosquito samples. We conclude that, while there are a few advantages in using mosquito blood meals, the use of flies in the detection of wildlife in a suburban environment is more effective in terms of collection of samples and detection of vertebrates, although this technique is limited to few mammal species in the urban environment. 

--------------------------------------------------------------
## Scripts
This project contains the R scripts to run the statistical analysis and create the figures included in the above paper. The scripts are prepared to be run setting the project folder as working directory. All packages required are included in the file "source_packages". The data required to run the scrips is provided in the "data" folder. The created figures are included in the "output" folder. There are three scripts that can be ran independently:
* _R_script_map_study_area_: creates map of the study area or Figure 1 of the paper
* _R_script_plot_summaries_: creates the figures 2, 3 and S1 
* _R_script_GLM_PCOA_: runs the statistical analyses and creates the associated figures 4 and 5

