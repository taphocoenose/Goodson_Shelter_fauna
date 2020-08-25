# This script cleans and organizes faunal data for Goodson
# Shelter.
# 
# Last updated on a Windows 10 machine, August 19, 2020.
# Questions? rbreslawski@smu.edu

# Load libraries
library(parallel)
library(data.table)
library(reshape2)
library(stringr)

# Import data
fd <- read.csv("01_data_fauna.csv", header=TRUE,
               stringsAsFactors=FALSE)
pd <- read.csv("01_data_provenience.csv", header=TRUE,
               stringsAsFactors=FALSE)

# Replace non-numeric values with NA
pd$Level <- suppressWarnings(as.numeric(pd$Level))
pd$EDM <- suppressWarnings(as.numeric(pd$EDM))
pd$N <- suppressWarnings(as.numeric(pd$N))
pd$E <- suppressWarnings(as.numeric(pd$E))
pd$Z <- suppressWarnings(as.numeric(pd$Z))

# Create data frame of levels and elevations
Levels <- data.frame(Level=seq(from=30, to=179, by=1),
                     Starting=seq(from=102.5, to=95.05, by=-0.05))
Levels$Ending <- Levels$Starting-0.05

# Remove Zs from specimens without EDM points
pd$Z <- sapply(1:nrow(pd), function(x){
  ifelse(is.na(pd$Z[x]), return(NA), 
         ifelse((pd$Z[x]*100 %% 5 == 0) & is.na(pd$EDM[x]),
                return(NA), return(pd$Z[x])))
})

# For specimens with Zs, replace the reported level with the
# level from the Level's data frame.
pd$Level <- sapply(1:nrow(pd), function(x){
  ifelse(!is.na(pd$Z[x]), 
         return(max(Levels$Level[which(Levels$Starting > pd$Z[x])])),
         return(pd$Level[x]))
})

# Create data frame of E and N by block number
BlockLetters <- paste(sort(rep(LETTERS[1:26], 30)), rep(seq(1:30), 26), sep="")
BlockLetters <- lapply(BlockLetters, function(x){
  sapply(1:25, function(z) paste(x,"-",z,sep=""))
})
BlockLetters <- do.call("c", BlockLetters)
GridDF <- data.frame(Block=BlockLetters, stringsAsFactors=FALSE)
GridDF$E <- rep(seq(from=0, to=4, by=1), 5*26*30)
GridDF$E <- GridDF$E + sort(rep(seq(from=0, to=25*5, by=5), 30*25))
GridDF$N <- rep(sort(rep(seq(from=0, to=4, by=1), 5)), 30*26)
GridDF$N <- GridDF$N + rep(sort(rep(seq(from=0, to=29*5, by=5), 25)), 26)
# Adjust for Goodson grid
GridDF$E <- GridDF$E + 940
GridDF$N <- GridDF$N + 930

# Fill in missing Ns by block
pd$N <- sapply(1:nrow(pd), function(x){
  if(!is.na(pd$N[x])){return(pd$N[x])}
  else{
    if(str_count(pd$Bag..[x], "-")==2){
      EndChar <- gregexpr("-", pd$Bag..[x])[[1]][2] - 1
      Block <- substr(pd$Bag..[x], 1, EndChar)
      if(Block %in% GridDF$Block){
        return(GridDF$N[which(GridDF$Block==Block)])
      }else{return(NA)}
    }else{return(NA)}
  }
})

# Fill in missing Es by block
pd$E <- sapply(1:nrow(pd), function(x){
  if(!is.na(pd$E[x])){return(pd$E[x])}
  else{
    if(str_count(pd$Bag..[x], "-")==2){
      EndChar <- gregexpr("-", pd$Bag..[x])[[1]][2] - 1
      Block <- substr(pd$Bag..[x], 1, EndChar)
      if(Block %in% GridDF$Block){
        return(GridDF$E[which(GridDF$Block==Block)])
      }else{return(NA)}
    }else{return(NA)}
  }
})

##################################################
############--- MISSING INFO ---##################

# Missing bags
uniqueBags <- unique(fd$Bag.)
missingBags <- uniqueBags[!(uniqueBags %in% pd$Bag..)]
missingIDs <- sapply(missingBags, function(x){
  nums <- fd$aID[which(fd$Bag.==x)]
  return(paste(nums[1],"-",nums[length(nums)],sep=""))
})


# Vector of elements not comprised of bone
nonbone <- c("SHELL","CTMR","CTMX","CT","CUN","CUUN","dPM2MR","dPM4MR",
             "dPM3MR","I1MX","I1U","I1UN","I2MR","IMR","IMX","IUN","IUS",
             "IUUN","M1MX","M1MR","M2MX","M2MR","M3MX","M3MR","MMR","MMX",
             "MUMR","MUMX","PM2MR","PM2MX","PM3MR","PM3MX","PM4MR","PM4MX",
             "PMMX","PMUMX", "PMUMR", "TFR", "IUMX")

#############################################################
###############------ APPEND PROV INFO TO FAUNAL DATA ----###

# Make cluster
cl <- makeCluster(detectCores() - 1)
# Export p for parallel processing
clusterExport(cl=cl, varlist="pd")

# Append prov info from pd to fd using 5 sapply functions
fd$Level <- as.numeric((parSapply(cl, fd$Bag., function(x){
  if(x %in% pd$Bag..){pd$Level[which(pd$Bag..==x)][1]}
  else{-1}})))
fd$N <- parSapply(cl, fd$Bag., function(x){
  if(x %in% pd$Bag..){pd$N[which(pd$Bag..==x)][1]}
  else{-1}})
fd$E <- parSapply(cl, fd$Bag., function(x){
  if(x %in% pd$Bag..){pd$E[which(pd$Bag..==x)][1]}
  else{-1}})
fd$Z <- parSapply(cl, fd$Bag., function(x){
  if(x %in% pd$Bag..){pd$Z[which(pd$Bag..==x)][1]}
  else{-1}})
fd$EDM <- parSapply(cl, fd$Bag., function(x){
  if(x %in% pd$Bag..){pd$EDM[which(pd$Bag..==x)][1]}
  else{-1}})

stopCluster(cl)

# Calulate fresh fracture index
fd_FFI <- fd[which(fd$FFI.A != "-" & fd$FFI.O != "-" & fd$FFI.T != "-" ),]
fd_FFI$FFI.A <- as.numeric(fd_FFI$FFI.A)
fd_FFI$FFI.T <- as.numeric(fd_FFI$FFI.T)
fd_FFI$FFI.O <- as.numeric(fd_FFI$FFI.O)
fd_FFI$FFI.A <- sapply(fd_FFI$FFI.A, function(x) ifelse(x>2, 2, x))
fd_FFI$FFI.T <- sapply(fd_FFI$FFI.T, function(x) ifelse(x>2, 2, x))
fd_FFI$FFI.O <- sapply(fd_FFI$FFI.O, function(x) ifelse(x>2, 2, x))
fd_FFI$FFI <- sapply(1:nrow(fd_FFI), function(z) mean(c(fd_FFI$FFI.A[z], 
                                                        fd_FFI$FFI.O[z],
                                                        fd_FFI$FFI.T[z])))

##############################################################
#####---- CHECK FOR MISSING TAXA AND TALLY TAXA BY LEVEL --###

# Create frequency tables to check for misspelled taxa and tally NISP
TaxonSize <- data.frame(table(fd$Size))
colnames(TaxonSize) <- c("Size.Class","NISP")
TaxonSize$NISP <- sapply(TaxonSize$Size.Class, function(x){
  sum(fd$count[which(fd$Size==x)])})
TaxonClass <- data.frame(table(fd$Class))
colnames(TaxonClass) <- c("Class","NISP")
TaxonClass$NISP <- sapply(TaxonClass$Class, function(x){
  sum(fd$count[which(fd$Class==x)])})
TaxonOrder <- data.frame(table(fd$Order))
colnames(TaxonOrder) <- c("Order","NISP")
TaxonOrder$NISP <- sapply(TaxonOrder$Order, function(x){
  sum(fd$count[which(fd$Order==x)])})
TaxonFamily <- data.frame(table(fd$Family))
colnames(TaxonFamily) <- c("Family","NISP")
TaxonFamily$NISP <- sapply(TaxonFamily$Family, function(x){
  sum(fd$count[which(fd$Family==x)])})
TaxonGenus <- data.frame(table(fd$Genus))
colnames(TaxonGenus) <- c("Genus","NISP")
TaxonGenus$NISP <- sapply(TaxonGenus$Genus, function(x){
  sum(fd$count[which(fd$Genus==x)])})
TaxonSpecies <- data.frame(table(fd$Species))
colnames(TaxonSpecies) <- c("Species","NISP")
TaxonSpecies$NISP <- sapply(TaxonSpecies$Species, function(x){
  sum(fd$count[which(fd$Species==x)])})

# Store taxa in a list
taxa <- list(TaxonSize,TaxonClass,TaxonOrder,TaxonFamily,
             TaxonGenus,TaxonSpecies)

# Remove taxonomic data frames
rm(TaxonSize, TaxonClass, TaxonOrder, TaxonFamily, 
   TaxonGenus, TaxonSpecies)

# Find total number of taxa
ntaxa <- 0
for(y in 1:length(taxa)){ntaxa <- ntaxa + nrow(taxa[[y]])}

# Remove possible NA counts
fd <- fd[!is.na(fd$count),]

# NISP with missing prov info
NISPnoProv <- data.frame(Reason=c("Missing bag info", "No prov on bag"),
                         Numbags=c(length(missingBags)-1,
                          length(unique(fd$Bag.[is.na(fd$Level)]))), 
                         NISP=c(sum(fd$count[fd$Bag. %in% missingBags]),
                          sum(fd$count[is.na(fd$Level)])), stringsAsFactors=FALSE)

# Replace level NAs with -2
fd$Level[is.na(fd$Level)] <- -2

# Create df to store data by level, allocate columns for each taxon
ldat <- data.frame(level=sort(unique(as.numeric(fd$Level))))

# Create list of taxonomic IDs and create column names for ldat
taxaindex <- data.frame(rbindlist(taxa, idcol="tLevel"), stringsAsFactors=FALSE)
colnames(taxaindex)[2] <- "taxon"
taxaindex$taxon <- as.character(taxaindex$taxon)
taxaindex <-  taxaindex[!is.na(taxaindex$NISP),]
ldat <- data.frame(cbind(ldat, matrix(0,nrow(ldat),nrow(taxaindex))))
colnames(ldat)[2:ncol(ldat)] <- taxaindex$taxon

# NISP by taxon by level
# 
# # Make cluster
cl <- makeCluster(detectCores() - 1)
# Export p for parallel processing
clusterExport(cl=cl, varlist="pd")

for(t in 1:nrow(taxaindex)){
  t_lev <- fd[,taxaindex$tLevel[t]+2]
  clusterExport(cl=cl, varlist=c("t_lev", "fd", "taxaindex", "t"))
  ldat[,t+1] <- parSapply(cl, ldat[,1], function(x){
    sum(fd$count[which((fd$Level==x) & (t_lev==taxaindex$taxon[t]))])
  })
}

# Stop cluster
stopCluster(cl)

# Save cleaned data
save(fd, fd_FFI, ldat, Levels, NISPnoProv, pd,
     taxa, taxaindex, nonbone, missingIDs,
     uniqueBags, file="Goodson_cleaned.RData")
