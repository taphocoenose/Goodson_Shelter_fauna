# This script analyzes deer element frequencies at Goodson Shelter.
# 
# Last updated on a Windows 10 machine, August 20, 2020.
# Questions? rbreslawski@smu.edu

# Load libraries
library(ggplot2)
library(rethinking)
library(patchwork)

# Load deer and caribou skeletal part density values
DeerSEA <- read.csv("03_DeerSEA.csv", header=TRUE,
                    stringsAsFactors=FALSE, 
                    na.strings=c("", "NA"))
# Load data
load("Goodson_cleaned.RData")

# Vector of observed scan deer scan sites from Goodson
DeerScans <- c(fd$Scan1[which(fd$Size=="m" & fd$Order=="Artiodactyla")],
               fd$Scan2[which(fd$Size=="m" & fd$Order=="Artiodactyla")],
               fd$Scan3[which(fd$Size=="m" & fd$Order=="Artiodactyla")],
               fd$Scan4[which(fd$Size=="m" & fd$Order=="Artiodactyla")])
DeerScans <- data.frame(table(DeerScans[which(DeerScans!="-")]))
DeerScans$Var1 <- as.character(DeerScans$Var1)

# Construct element data frame based on unique elements defined in Madrigal and
# Holt's 2002 paper on deer element return rates
DeerElements <- lapply(unique(DeerSEA$UtilElement[!is.na(DeerSEA$UtilElement)]),
                              function(x){
                                
                              # Subset scan site freq table by element x
                              ScanSub <- DeerScans[DeerScans$Var1 %in% 
                                         DeerSEA$ScanSite[which(DeerSEA$UtilElement==x)],]
                              
                              # Get scan site for element x with highest frequency
                              SiteX <- ScanSub$Var1[which.max(ScanSub$Freq)]
                              
                              # Return 1-case data frame corresponding to SiteX
                              return(DeerSEA[which(DeerSEA$ScanSite==SiteX),])
                              
                              })

# Combine list of element into a single data frame for modelling
DeerElements <- do.call("rbind", DeerElements)
# Generate scan site based MNE for each element
DeerElements$MNE <- sapply(DeerElements$ScanSite, function(x){
  DeerScans$Freq[which(DeerScans$Var1==x)]})
# Generate MAU values
DeerElements$MAU <- DeerElements$MNE/DeerElements$AnatomFreq
# Create log anatomical frequency for Poisson GLM offset
DeerElements$LogAnFreq <- log(DeerElements$AnatomFreq)
DeerElements$EIndex <- seq(1:nrow(DeerElements))

# Cycle through predictors for simple correlations
for(z in c(5,8,9,10,11)){
  
  # Creata date frame for specified variables
  maucor <- data.frame(MAU=DeerElements$MAU,
                       util=DeerElements[,z])
  
  # Store correlation as variable
  corresults <- cor.test(y=maucor$MAU, 
                         x=maucor$util,
                         method="kendall", exact=FALSE)
  
  # print correlation
  cat(paste(colnames(DeerElements)[z], "x MAU: tau =", 
                     round(corresults$estimate, 3),
                     ", p =", round(corresults$p.value, 3),"\n"))
  
  ifelse(z==5, xlab <- "Volume density", 
         ifelse(z==8, xlab <- "Meat kcal/hr",
                ifelse(z==9, xlab <- "Marrow kcal/hr",
                       ifelse(z==10, xlab <- "Meat kcal",
                              xlab <- "Marrow kcal"))))
  
  tempplot <- ggplot(maucor, aes(x=util, y=MAU))+
    geom_smooth(method=lm, se=FALSE, color="grey")+
    geom_point(color="black")+
    annotate("text", label=paste("Kendall's tau =",
                                 round(corresults$estimate, 3)),
             hjust=1, x=0.99*max(maucor$util),
             y=3.8)+
    annotate("text", label=paste("p =",
                                 round(corresults$p.value, 3)),
             hjust=1, x=0.99*max(maucor$util),
             y=3.5)+
    labs(x=xlab, y="MAU")+
    theme(panel.background=element_rect(color="grey", 
                                        fill="white"),
          panel.grid=element_blank())
  
  assign(paste0("tp",z), tempplot)
}
cor_plot <- tp5+tp8+tp9+tp10+tp11+plot_layout(nrow=2)

if(!dir.exists("Figures")) dir.create("Figures")
ggsave("Figures/Deer_Element_Correlations.pdf", plot=cor_plot,
       device="pdf", width=10, height=6, units="in")

# Standardize predictors
DeerElements$CaribouDensity <- (DeerElements$CaribouDensity-
                                  mean(DeerElements$CaribouDensity))/
                                  sd(DeerElements$CaribouDensity)
DeerElements$MeatRR <- (DeerElements$MeatRR-
                                  mean(DeerElements$MeatRR))/
                                  sd(DeerElements$MeatRR)
DeerElements$MarRR <- (DeerElements$MarRR-
                          mean(DeerElements$MarRR))/
                          sd(DeerElements$MarRR)
DeerElements$MeatKcal <- (DeerElements$MeatKcal-
                          mean(DeerElements$MeatKcal))/
                          sd(DeerElements$MeatKcal)
DeerElements$MarKcal <- (DeerElements$MarKcal-
                         mean(DeerElements$MarKcal))/
                         sd(DeerElements$MarKcal)

### Fit GLMs ###
################

# Predictors: Bone density
ModelDeer_D <- map2stan(alist(
  MNE ~ dpois(lambda),
  log(lambda) <- a + a_unique[EIndex] + LogAnFreq + 
    Bd*CaribouDensity,
  a ~ dnorm(0, 5),
  Bd ~ dnorm(0, 5),
  a_unique[EIndex] ~ dnorm(0, sigma_unique),
  sigma_unique ~ dexp(2)),
  data=DeerElements, chains=4, iter=4e3, cores=4,
  start=list(a=0, Bd=0),
  control=list(adapt_delta=0.99))

# Predictors: Bone density + MeatRR
ModelDeer_D_MeRR <- map2stan(alist(
  MNE ~ dpois(lambda),
  log(lambda) <- a + a_unique[EIndex] + LogAnFreq + 
    Bd*CaribouDensity + Bmeat*MeatRR,
  a ~ dnorm(0, 5),
  c(Bd, Bmeat) ~ dnorm(0, 5),
  a_unique[EIndex] ~ dnorm(0, sigma_unique),
  sigma_unique ~ dexp(2)),
  data=DeerElements, chains=4, iter=8e3, cores=4,
  start=list(a=0, Bd=0, Bmeat=0),
  control=list(adapt_delta=0.9999))

# Predictors: Bone density + MarRR
ModelDeer_D_MaRR <- map2stan(alist(
  MNE ~ dpois(lambda),
  log(lambda) <- a + a_unique[EIndex] + LogAnFreq + 
    Bd*CaribouDensity + Bmar*MarRR,
  a ~ dnorm(0, 5),
  c(Bd, Bmar) ~ dnorm(0, 5),
  a_unique[EIndex] ~ dnorm(0, sigma_unique),
  sigma_unique ~ dexp(2)),
  data=DeerElements, chains=4, iter=8e3, cores=4,
  start=list(a=0, Bd=0, Bmar=0),
  control=list(adapt_delta=0.99))

# Predictors: Bone density + MeatKcal
ModelDeer_D_MeKc <- map2stan(alist(
  MNE ~ dpois(lambda),
  log(lambda) <- a + a_unique[EIndex] + LogAnFreq + 
    Bd*CaribouDensity + Bmeat*MeatKcal,
  a ~ dnorm(0, 5),
  c(Bd, Bmeat) ~ dnorm(0, 5),
  a_unique[EIndex] ~ dnorm(0, sigma_unique),
  sigma_unique ~ dexp(2)),
  data=DeerElements, chains=4, iter=4e3, cores=4,
  start=list(a=0, Bd=0, Bmeat=0),
  control=list(adapt_delta=0.99))

# Predictors: Bone density + MarKcal
ModelDeer_D_MaKc <- map2stan(alist(
  MNE ~ dpois(lambda),
  log(lambda) <- a + a_unique[EIndex] + LogAnFreq + 
    Bd*CaribouDensity + Bmar*MarKcal,
  a ~ dnorm(0, 5),
  c(Bd, Bmar) ~ dnorm(0, 5),
  a_unique[EIndex] ~ dnorm(0, sigma_unique),
  sigma_unique ~ dexp(2)),
  data=DeerElements, chains=4, iter=4e3, cores=4,
  start=list(a=0, Bd=0, Bmar=0),
  control=list(adapt_delta=0.99))

# Predictors: Bone density + MeatRR + MarRR
ModelDeer_D_MeRR_MaRR <- map2stan(alist(
  MNE ~ dpois(lambda),
  log(lambda) <- a + a_unique[EIndex] + LogAnFreq + 
    Bd*CaribouDensity + Bmeat*MeatRR + Bmar*MarRR,
  a ~ dnorm(0, 5),
  c(Bd, Bmeat, Bmar) ~ dnorm(0, 5),
  a_unique[EIndex] ~ dnorm(0, sigma_unique),
  sigma_unique ~ dexp(2)),
  data=DeerElements, chains=4, iter=8e3, cores=4,
  start=list(a=0, Bd=0, Bmeat=0, Bmar=0),
  control=list(adapt_delta=0.9999))

# Predictors: Bone density + MeatKcal + MarKcal
ModelDeer_D_MeKc_MaKc <- map2stan(alist(
  MNE ~ dpois(lambda),
  log(lambda) <- a + a_unique[EIndex] + LogAnFreq + 
    Bd*CaribouDensity + Bmeat*MeatKcal + Bmar*MarKcal,
  a ~ dnorm(0, 10),
  c(Bd, Bmeat, Bmar) ~ dnorm(0, 10),
  a_unique[EIndex] ~ dnorm(0, sigma_unique),
  sigma_unique ~ dexp(2)),
  data=DeerElements, chains=4, iter=4e3, cores=4,
  start=list(a=0, Bd=0, Bmeat=0, Bmar=0),
  control=list(adapt_delta=0.99))

# Display fitted model parameters and HMC diagnostics
print("Predictors: Density")
print(precis(ModelDeer_D, depth=2, prob=0.95))
print("Predictors: Density + Marrow RR")
print(precis(ModelDeer_D_MaRR, depth=2, prob=0.95))
print("Predictors: Density + Meat RR")
print(precis(ModelDeer_D_MeRR, depth=2, prob=0.95))
print("Predictors: Density + Marrow Kcal")
print(precis(ModelDeer_D_MaKc, depth=2, prob=0.95))
print("Predictors: Density + Meat Kcal")
print(precis(ModelDeer_D_MeKc, depth=2, prob=0.95))
print("Predictors: Density + Meat RR + Marrow RR")
print(precis(ModelDeer_D_MeRR_MaRR, depth=2, prob=0.95))
print("Predictors: Density + Meat Kcal + Marrow Kcal")
print(precis(ModelDeer_D_MeKc_MaKc, depth=2, prob=0.95))

# WAIC model comparison
print(compare(ModelDeer_D, ModelDeer_D_MaRR, ModelDeer_D_MeRR,
        ModelDeer_D_MaKc, ModelDeer_D_MeKc, 
        ModelDeer_D_MeRR_MaRR, ModelDeer_D_MeKc_MaKc))

comparison <- compare(ModelDeer_D, ModelDeer_D_MaRR, ModelDeer_D_MeRR,
              ModelDeer_D_MaKc, ModelDeer_D_MeKc, 
              ModelDeer_D_MeRR_MaRR, ModelDeer_D_MeKc_MaKc)


# Create data frame of posterior parameter estimates
postpars <- data.frame(Model=sort(rep(1:7, 6), decreasing=TRUE), 
                       ParIndex=rep(1:6, 7),
                       PointEst=c(-0.59, 0.80, NA, NA, NA, NA,
                                  -0.56, 0.93, 0.61, NA, NA, NA,
                                  -0.53, 0.99, NA, 0.87, NA, NA,
                                  -0.61, 0.55, NA, NA, -0.53, NA,
                                  -0.60, 0.79, NA, NA, NA, -0.04,
                                  -0.60, 0.77, 0.52, NA, -0.33, NA,
                                  -0.53, 1.03, NA, 0.88, NA, 0.11),
                       Lower=c(-1.31, 0.10, NA, NA, NA, NA,
                               -1.24, 0.26, -0.13, NA, NA, NA,
                               -1.13, 0.41, NA, 0.25, NA, NA,
                               -1.33, -0.23, NA, NA, -1.35, NA,
                               -1.37, -0.03, NA, NA, NA, -0.89,
                               -1.30, -0.07, -0.30, NA, -1.20, NA,
                               -1.16, 0.38, NA, 0.25, NA, -0.55),
                       Upper=c(0.04, 1.53, NA, NA, NA, NA,
                               0.06, 1.64, 1.31, NA, NA, NA,
                               -0.03, 1.61, NA, 1.45, NA, NA,
                               0.06, 1.34, NA, NA, 0.25, NA,
                               0.09, 1.62, NA, NA, NA, 0.73,
                               0.04, 1.62, 1.29, NA, 0.50, NA,
                               0.00, 1.72, NA, 1.49, NA, 0.72))

# Cycle through each column of predictors and create a plot
for(u in 1:6){
  
  ptitle <- c("alpha", "beta[Density]", "beta[Marrow~Kcal]",
             "beta[Marrow~Kcal/hr]", "beta[Meat~Kcal]",
             "beta[Meat~Kcal/hr]")[u]
  
  # Subset data frame for plotting
  plotdf <- postpars[which(postpars$ParIndex==u),]
  
  tempparplot <- ggplot(data=plotdf, aes(y=Model))+
    annotate("segment", y=-Inf, yend=7.5, x=0, xend=0,
             color="white", linetype="longdash", size=1)+
    annotate("segment", y=-Inf, yend=7.5, x=c(-1,1), xend=c(-1,1),
             color="white", linetype="longdash", size=0.5)+
    geom_segment(aes(x=Lower, xend=Upper, yend=Model))+
    annotate("text", label=ptitle, parse=TRUE, x=0, y=7.5,
             size=4.5, vjust=0)+
    geom_point(aes(x=PointEst))
  
  if(u==1){
    
    tempparplot <- tempparplot +
      scale_y_continuous(breaks=seq(7,1,-1),
                         limits=c(0,8),
                         labels=c("7"="A", "6"="B", "5"="C",
                                  "4"="D", "3"="E", "2"="F",
                                  "1"="G"),
                         expand=c(0,0))+
      scale_x_continuous(limits=c(-2,2), 
                         breaks=seq(-2, 2, 1), expand=c(0,0))+
      theme(panel.grid=element_blank(),
            axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title=element_blank(),
            axis.line.x=element_line(color="grey"))
    
    parplot <- tempparplot
    
  }else{
    
    tempparplot <- tempparplot +
      scale_x_continuous(limits=c(-2,2), 
                         breaks=seq(-2, 2, 1),
                         expand=c(0,0))+
      scale_y_continuous(breaks=seq(7,1,-1), limits=c(0,8),
                         expand=c(0,0))+
      theme(panel.grid=element_blank(),
            axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title=element_blank(),
            axis.line.x=element_line(color="grey"))
    
    parplot <- parplot | tempparplot
    
  }
  
}

if(!dir.exists("Figures")) dir.create("Figures")
ggsave("Figures/Deer_Element_Models.pdf", plot=parplot,
       device="pdf", width=10, height=4, units="in")