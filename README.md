# Goodson Shelter Archaeofaunal Data and Analyses

This repository provides the raw data and R code for data cleaning, organization, and statistical analyses of the archaeofauna recovered from Goodson Shelter.

<i><b>01_data_fauna.csv</b></i>: Catalogued faunal specimens for the site.
<i><b>01_data_provenience.csv</b></i>: Vertical and horizontal provenience information for the faunal specimens.
<i><b>03_DeerSEA.csv</b></i>: Bone density and dietary utility index values for deer elements.


<i><b>01_Goodson_data_cleanup.R</b></i>: This script requires and reads the first two csv files. It joins the data tables and produces cleaned datasets summarized by taxon and provenience. It is necessary to run this script prior to running the other two R scripts. After it completes, it outputs <i>Goodson_cleaned.RData</i> to the working directory.
<i><b>02_Goodson_faunal_analyses</b></i>: This script requires and reads <i>Goodson_cleaned.RData</i>. It analyzes spatial patterns in burning and fragmentation, taxonomic richness, and an artiodactyl index (Artiodactyla_NISP/[Artiodactyla_NISP + Lagomorpha_NISP]). It outputs figures summarizing results to the subdirectory <i>working directory/Figures</i>.
<i><b>03_Goodson_deer_element_analyses</b></i>: This script requires and reads <i>Goodson_cleaned.RData</i> and <i>03_DeerSEA.csv</i>. It examines the relationships between deer element representation and several predictors (bone density and dietary utility indices). It outputs figures summarizing results to the subdirectory <i>working directory/Figures</i>.
