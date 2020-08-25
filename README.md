# Goodson Shelter Archaeofaunal Data and Analyses

This repository provides the raw data and R code for data cleaning, organization, and statistical analyses of the archaeofauna recovered from Goodson Shelter.

<i>01_data_fauna.csv</i>: Catalogued faunal specimens for the site.

<i>01_data_provenience.csv</i>: Vertical and horizontal provenience information for the faunal specimens.


<i>01_Goodson_data_cleanup.R</i>: This script requires and reads both csv files. It joins the data tables and produces cleaned datasets summarized by taxon and provenience. It is necessary to run this script prior to running the other two R scripts. After it completes, it outputs <i>Goodson_cleaned.RData</i> to the working directory.

<i>02_Goodson_faunal_analyses</i>: This script requires and reads <i>Goodson_cleaned.RData</i>. It analyzes spatial patterns in burning and fragmentation, taxonomic richness, and an artiodactyl index (Artiodactyla_NISP/[Artiodactyla_NISP + Lagomorpha_NISP]). It outputs figures summarizing results to the subdirectory <i>working directory/Figures</i>.

<i>03_Goodson_deer_element_analyses</i>: This script requires and reads <i>Goodson_cleaned.RData</i>. It examines the relationships between deer element representation and several predictors (bone density and dietary utility indices). It outputs figures summarizing results to the subdirectory <i>working directory/Figures</i>.
