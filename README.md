# Pant-diversity-2025
This is a computationally reproducible documents that contain the R code, verbal descriptions, and outputs. With this code all analyses in the paper can be reproduced, and the graphs, except the maps, can be remade. Code to calculate the effect size log ratio response(lnRR), to deal with pseudo-duplication of data and to model phylogenetic correction are provided.

***Please run Data_clean.R file to pre-process data set before you run main analysis
***Please run Supplementary Table 3-10&Supplementary Figs.9-42.R file before you run Fig 2.R, and run Supplementary Table 17 before Supplementary Figs 1-8.

Description of datasets
data of all trophic groups_2025_1001.xlsx: the whole data used for main analyse, there are 15 sheets inside;
0903Paired data for herbivore and crop.xlsx: the data for investigating the interactions between plant and herbivore performance;
0831Paired data for enemies and herbivore and crop.xlsx: the data for investigating the interactions between plant, herbivore, and natural enemy performance;
Fig1a_wordmap.csv: this dataset provides the geographic coordinates (latitude and longitude) of each study location. The sites were categorized according to trophic groups (i.e., predators, parasitoids, natural enemies [NEs], herbivores, and crops) and used to generate Fig. 1a in the main text.
Fig1c_number_of_species.csv: this dataset contains information on the number of added plant species, crop life form (herbaceous or woody), crop production purpose (food or cash), and climate zone for each study, which were used to generate Fig. 1b in the main text.
Figs2a-2f.csv: These six datasets contain summary information on the number of studies, number of effect sizes, mean effect size values, and their 95% confidence intervals across all studies, categorized by crop life form, crop production purpose, herbivore diet breadth, natural enemy host range, diversification strategy, experimental design, and climate zone.

Description of columns in data of all trophic groups_2025_1001.xlsx
Study_ID:
Variable name: Study_ID
Abbreviation (in scripts): S_ID
Description:  the paper serial number
Type: Continuous
Values: according to the extracted paper

CP_Category: 
Variable name: CP_Category
Abbreviation (in scripts): CP_Ct
Description: Purpose of the crop grown (Food crop or Cash crop)
Type: Categorical
Values: Food_crop, Cash_crop

CP_Category3: 
Variable name: CP_Category3
Abbreviation (in scripts): CP_C3
Description: Crop life form (Herbaceous plant or Woody plant)
Type: Categorical
Values: Herbaceous_plant, Woody_plant

Main_herbivore_category:
Variable name: Main_herbivore_category
Abbreviation (in scripts): Mn__
Description: Herbivore diet breath (Generalist or Specialist)
Type: Categorical
Values: Generalist, Specialist

NE_category:
Variable name: NE_category
Abbreviation (in scripts): NE_c
Description: Natural enemy host range (Generalist or Specialist)
Type: Categorical
Values: Generalist, Specialist

NCP_Pattern:
Variable name: NCP_Pattern
Abbreviation (in scripts): NCP_P
Description: Crop diversification strategy (Intercropping, Cover cropping, and Sown field margins)
Type: Categorical
Values: Intercropping, Cover_cropping, Sown_field_margins

Study_design:
Variable name: Study_design
Abbreviation (in scripts): Std_
Description: Experimental study design (Natural, Semi-natural, and Controlled)
Type: Categorical
Values: Natural, Semi_natural, Controlled

Zone:
Variable name: Zone
Abbreviation (in scripts): Zone
Description: Study carried climatic zone (Temperate or Tropical zones)
Type: Categorical
Values: Temperate, Tropic

NCP_Number:
Variable name: NCP_Number
Abbreviation (in scripts): NCP_N
Description: The number of added plant species in the treatment (e.g, 1, 2, 3)
Type: Continuous
Values: 1, 2, 3,...

TREATMENT_REPLICATION_NUMBER:
Variable name: TREATMENT_REPLICATION_NUMBER
Abbreviation (in scripts): TREATMENT_R
Description: The sample size in treatment group.
Type: Continuous
Values: according to the extracted data in the paper.

TREATMENT_MEAN:
Variable name: TREATMENT_MEAN
Abbreviation (in scripts): TREATMENT_M
Description: The indicator value in treatment group
Type: Continuous
Values: according to the extracted data in the paper

TREATMENT_SD：
Variable name: TREATMENT_MEAN
Abbreviation (in scripts): TREATMENT_M
Description: The indicator value in treatment group
Type: Continuous
Values: according to the extracted data in the paper

CONTROL_REPLICATION_NUMBER：
Variable name: CONTROL_REPLICATION_NUMBER
Abbreviation (in scripts): CONTROL_R
Description: The sample size in control group.
Type: Continuous
Values: according to the extracted data in the paper.

CONTROL_MEAN:
Variable name: CONTROL_MEAN
Abbreviation (in scripts): CONTROL_M
Description: The indicator value in control group
Type: Continuous
Values: according to the extracted data in the paper

CONTROL_SD:
Variable name: CONTROL_MEAN
Abbreviation (in scripts): CONTROL_M
Description: The indicator value in control group
Type: Continuous
Values: according to the extracted data in the paper

CP_Latin_name_after_check：
Variable name: CP_Latin_name_after_check
Abbreviation (in scripts): CP_Ltn_nm_f_
Description: The Latin name of main crop and checked in www.theplantlist.com.
Type: Taxonomic
Values: according to the extracted data in the paper, e.g., Zea_mays


LnR: the indicator value in treatment groups divided by the indicator value in the control group.

Indicator: the concrete indicator used to measure invertebrate predator reproduction, predator diversity, predator predation, invertebrate parasitoid growth, parasitoid reproduction, parasitoid diversity, parasitoid parasitism, invertebrate natural enemies (NEs) reproduction, NEs diversity, invertebrate herbivore growth, herbivore reproduction, herbivore damage, crop growth, crop reproduction, and crop quality.

Effect size identity
To handle non-independence in effect sizes, we have now added a random effect as a unique identifier for each effect size. Our use of random effects now allowed to account for different variances between the studies. Above R code was used for generating a sequence for effect sizes.
# generate a sequence for effect sizes for subsequent calculation of variance-covariance matrixData <- Data[order(Data$S_ID),]Cod<-table(Data$S_ID)esid<-NULLfor (i in 1:length(Cod)){  esid <- c(esid,seq(1,Cod[i],1))}Data$EsID <- esid

Generating new columns for subsequent analysis
Trophic_group: including Predator_performance, Parasitoid_performance, NEs_performance, Herbivore_performance, Crop_performance.

Trophic_group_response: including Predator_reproduction_response, Predator_diversity_response, Predator_predation_response, Parasitoid_growth_response, Parasitoid_reproduction_response, Parasitoid_diversity_response, Parasitoid_parasitism_response, NEs_reproduction_response, NEs_diversity_response, Herbivore_growth_response, Herbivore_reproduction_response, Herbivore_damage_response, Crop_growth_response, Crop_reproduction_response, Crop_quality_response.

added_plant_species: the logarithm of the number of added plant species in the plant diversity treatment over the control.
