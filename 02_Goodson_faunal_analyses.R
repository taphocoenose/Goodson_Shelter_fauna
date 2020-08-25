# This script performs analyses and creates plots based
# on the cleaned and organized Goodson faunal data.
# 
# Last updated on a Windows 10 machine, August 22, 2020.
# Questions? rbreslawski@smu.edu

# Load libraries
library(ggplot2)
library(rethinking)
library(parallel)
library(reshape2)
library(patchwork)
library(stringr)
library(data.table)

# Load cleaned data
load("Goodson_cleaned.RData")

# Fix randomization
set.seed(67896)

### ARTIODACTYL INDEX ANALYSIS ###
##################################

# Artiodactyl index data frame
AIdf <- data.frame(level=ldat$level, 
                   a=ldat$Artiodactyla, 
                   l=ldat$Lagomorpha)
AIdf$nisp <- AIdf$a + AIdf$l
# Subset levels for analysis (some do not have adequate
# NISP values for analysis, some lack prov. infomations)
AIdf <- AIdf[which(AIdf$level>85 & AIdf$nisp >= 1 & 
                     AIdf$level <= 130),]

# Calculate AI and confidence intervals by level
AIdf$AI_mu <- AIdf$a/AIdf$nisp
AIdf$AI_lower <- sapply(1:nrow(AIdf), function(x){
  binom.test(AIdf$a[x], AIdf$nisp[x], conf.level=0.95)$conf.int[1]})
AIdf$AI_upper <- sapply(1:nrow(AIdf), function(x){
  binom.test(AIdf$a[x], AIdf$nisp[x], conf.level=0.95)$conf.int[2]})

# Fit Bayesian binomial GLM to artiodactyl index
AImodel <- map2stan(alist(
  # Likelihood
  a ~ dbinom(nisp, p),
  # Linear model
  logit(p) <- al + Bl*level,
  # priors
  al ~ dnorm(0,10),
  Bl ~ dnorm(0,10)),
  # Stan args
  data=AIdf, iter=1e4, warmup=2e3, chains=4, cores=4)

# Get HPDI for relationship
LevelSeq1 <- seq(from=min(AIdf$level), to=max(AIdf$level), length.out=100)
LinkedMean1 <- link(AImodel, data=data.frame(level=LevelSeq1))
LinkedMean.mean <- apply(LinkedMean1, 2, mean)
LevelSeq2 <- seq(from=min(AIdf$level), to=max(AIdf$level), length.out=100)
# The following line takes 1 million samples for each LevelSeq2 value,
# which may take some time.
LinkedMean2 <- link(AImodel, data=data.frame(level=LevelSeq2), n=1e6)
LinkedMean.HPDI <- apply(LinkedMean2, 2, HPDI, prob=0.95)

# Plot artiodactyl index
AIplot <- ggplot(data=AIdf)+
  annotate("line", x=LevelSeq2, y=LinkedMean.HPDI[1,], color="blue", 
           size=1, linetype="dotted")+
  annotate("line", x=LevelSeq2, y=LinkedMean.HPDI[2,], color="blue", 
           size=1, linetype="dotted")+
  annotate("line", x=LevelSeq1, y=LinkedMean.mean, color="blue", 
           size=0.5)+
  geom_segment(aes(y=AI_lower, yend=AI_upper, x=level, xend=level),
               size=0.5)+
  geom_point(aes(y=AI_mu, x=level), shape=16, size=3)+
  labs(y="Artiodactyl index", x="Levels")+
  scale_x_reverse(breaks=seq(from=85, to=125, by=5))+
  coord_flip()+
  theme(panel.border=element_rect(fill=NA,size=0.5,color="grey"),
        panel.background=element_blank())

if(!dir.exists("Figures")) dir.create("Figures")
ggsave("Figures/Artiodactyl_Index.pdf", plot=AIplot,
       device="pdf", width=7, height=7, units="in")

############################################
###### CREATE ARTIODACTYL DATA FRAMES ######
artio_medium <- fd[which(fd$Order=="Artiodactyla" & fd$Size=="m"),]
artio_medium$elementside <- paste(artio_medium$Element, artio_medium$Side)
artio_m <- rep(artio_medium$elementside, artio_medium$count)
artio_m <- data.frame(table(artio_m))

####------ GASTRODPOD ANALYSIS ---------######
##############################################

# Create data frame for gastropod shell size analysis.
gastro_df <- suppressWarnings(fd[which((fd$Class=="Gastropoda") & 
                                         (fd$Portion=="CO")),])

# Test for dimension lists with bad formatting
gastro_bool <- sapply(gastro_df$MaxD, function(x){
  m <- unlist(strsplit(x,","))
  metrics <- suppressWarnings(as.numeric(m))
  if(anyNA(metrics)){
    return(1)
  }else{
    return(0)
  }
})

if(1 %in% gastro_bool){
  print("ROWS WITH MAX DIMENSION FORMATTING ERRORS:")
  print(gastro_df[which(gastro_bool==1),1:2])
  # If all lists are properly formatted, create list of data frames
  # with dimensions for each aID
}else{
  rm(gastro_bool)
  gastro_list <- lapply(1:nrow(gastro_df), function(x){
    fvec <- as.numeric(unlist(strsplit(gastro_df$MaxD[x],",")))
    data.frame(aID=rep(gastro_df$aID[x],length(fvec)),
               Bag=rep(gastro_df$Bag.[x],length(fvec)),
               Element=rep(gastro_df$Element[x],length(fvec)),
               Burn=rep(gastro_df$Burn[x],length(fvec)),
               Level=rep(gastro_df$Level[x],length(fvec)),
               N=rep(gastro_df$N[x],length(fvec)),
               E=rep(gastro_df$E[x],length(fvec)),
               Z=rep(gastro_df$Z[x],length(fvec)),
               MaxD=fvec)})
  
  gastropods <- data.frame(rbindlist(gastro_list,idcol="listID"))
  rm(gastro_list)
}

# Remove specimens with unusual levels (probably entry errors)
gastropods <- gastropods[which(gastropods$Level>10),]


###-------- FRAGMENTATION ANALYSIS --------###
##############################################

# Create data frame for fragment size analysis, 
# remove shell, teeth, and complete specimens
d_frag <- suppressWarnings(fd[which(!is.na(fd$MaxD) & !is.na(fd$Level) &
                                      !(fd$Element %in% nonbone) & 
                                      !(fd$Portion=="CO")), ])

# Test for dimension lists with bad formatting
frag_bool <- sapply(d_frag$MaxD, function(x){
  m <- unlist(strsplit(x,","))
  metrics <- suppressWarnings(as.numeric(m))
  if(anyNA(metrics)){
    return(1)
  }else{
    return(0)
  }
})

if(1 %in% frag_bool){
  print("ROWS WITH MAX DIMENSION FORMATTING ERRORS:")
  print(d_frag[which(frag_bool==1),1:2])
  # If all lists are properly formatted, create list of data frames
  # with dimensions for each aID
}else{
  rm(frag_bool)
  frags_list <- lapply(1:nrow(d_frag), function(x){
    fvec <- as.numeric(unlist(strsplit(d_frag$MaxD[x],",")))
    data.frame(aID=rep(d_frag$aID[x],length(fvec)),
               Bag=rep(d_frag$Bag.[x],length(fvec)),
               Element=rep(d_frag$Element[x],length(fvec)),
               Burn=rep(d_frag$Burn[x],length(fvec)),
               Level=rep(d_frag$Level[x],length(fvec)),
               N=rep(d_frag$N[x],length(fvec)),
               E=rep(d_frag$E[x],length(fvec)),
               Z=rep(d_frag$Z[x],length(fvec)),
               MaxD=fvec)})
  
  frags_df <- data.frame(rbindlist(frags_list,idcol="listID"))
  rm(frags_list)
}


# Test for dimension lists with bad formatting for FFI
frag_boolffi <- sapply(fd_FFI$MaxD, function(x){
  m <- unlist(strsplit(x,","))
  metrics <- suppressWarnings(as.numeric(m))
  if(anyNA(metrics)){
    return(1)
  }else{
    return(0)
  }
})

if(1 %in% frag_boolffi){
  print("ROWS WITH MAX DIMENSION FORMATTING ERRORS:")
  print(fd_FFI[which(frag_boolffi==1),1:2])
  # If all lists are properly formatted, create list of data frames
  # with dimensions for each aID
}else{
  rm(frag_boolffi)
  frags_listffi <- lapply(1:nrow(fd_FFI), function(x){
    fvec <- rep(fd_FFI$FFI[x], fd_FFI$count[x])
    data.frame(aID=rep(fd_FFI$aID[x],length(fvec)),
               Bag=rep(fd_FFI$Bag.[x],length(fvec)),
               Element=rep(fd_FFI$Element[x],length(fvec)),
               Burn=rep(fd_FFI$Burn[x],length(fvec)),
               Level=rep(fd_FFI$Level[x],length(fvec)),
               N=rep(fd_FFI$N[x],length(fvec)),
               E=rep(fd_FFI$E[x],length(fvec)),
               Z=rep(fd_FFI$Z[x],length(fvec)),
               FFI=fvec)})
  
  FFI_df <- data.frame(rbindlist(frags_listffi,idcol="listID"))
  rm(frags_listffi)
}

# Subset data for potential erroneous level entries, then jitter
# the x and y coordinates for plotting.
FFI_df <- FFI_df[which(FFI_df$Level>50),]
FFI_df$ljit <- sapply(FFI_df$Level, function(x) x + runif(1, -0.3,0.3))
FFI_df$xjit <- sapply(FFI_df$FFI, function(x) x + runif(1, -0.02,0.02))

# Remove specimens under 5 mm and unusual levels (probably entry errors)
frags_df <- frags_df[which(frags_df$MaxD>=5 & 
                             frags_df$Level <= 140 &
                             frags_df$Level >= 80),]
# Convert partially carbonized to general carbonized category
frags_df$Burn2 <- as.character(sapply(frags_df$Burn, 
                                      function(x) ifelse(x==2,1,x)))
# jitter level designation
frags_df$ljit <- sapply(frags_df$Level, function(x) x + runif(1, -0.3,0.3))

# Create data frame of fragmentation data by level
frLVL <- data.frame(Level=sort(unique(frags_df$Level)))
frLVL$NISP <- sapply(frLVL$Level, 
                     function(x) nrow(frags_df[which(frags_df$Level==x),]))
frLVL$propUnburned <- sapply(frLVL$Level, 
                             function(x) nrow(frags_df[which(frags_df$Burn2=="0" &
                                                               frags_df$Level==x),])/frLVL$NISP[which(frLVL$Level==x)])
frLVL$propCarbonized <- sapply(frLVL$Level, 
                               function(x) nrow(frags_df[which(frags_df$Burn2=="1" &
                                                                 frags_df$Level==x),])/frLVL$NISP[which(frLVL$Level==x)])
frLVL$propCalcined <- sapply(frLVL$Level, 
                             function(x) nrow(frags_df[which(frags_df$Burn2=="3" &
                                                               frags_df$Level==x),])/frLVL$NISP[which(frLVL$Level==x)])
frLVL$muUnburned <- sapply(frLVL$Level, function(x){
  mean(frags_df$MaxD[which(frags_df$Burn2=="0" & frags_df$Level==x)])})
frLVL$seUnburned <- sapply(frLVL$Level, function(x){
  sd(frags_df$MaxD[which(frags_df$Level==x & frags_df$Burn2=="0")])/
    sqrt(nrow(frags_df[which(frags_df$Level==x & frags_df$Burn2=="0"),])
    )})
frLVL$upperUnburned <- frLVL$muUnburned + 1.96*frLVL$seUnburned
frLVL$lowerUnburned <- frLVL$muUnburned - 1.96*frLVL$seUnburned
frLVL$muCarbonized<- sapply(frLVL$Level, function(x){
  mean(frags_df$MaxD[which(frags_df$Burn2=="1" & frags_df$Level==x)])})
frLVL$seCarbonized <- sapply(frLVL$Level, function(x){
  sd(frags_df$MaxD[which(frags_df$Level==x & frags_df$Burn2=="1")])/
    sqrt(nrow(frags_df[which(frags_df$Level==x & frags_df$Burn2=="1"),])
    )})
frLVL$upperCarbonized <- frLVL$muCarbonized + 1.96*frLVL$seCarbonized
frLVL$lowerCarbonized <- frLVL$muCarbonized - 1.96*frLVL$seCarbonized
frLVL$muCalcined <- sapply(frLVL$Level, function(x){
  mean(frags_df$MaxD[which(frags_df$Burn2=="3" & frags_df$Level==x)])})
frLVL$seCalcined <- sapply(frLVL$Level, function(x){
  sd(frags_df$MaxD[which(frags_df$Level==x & frags_df$Burn2=="3")])/
    sqrt(nrow(frags_df[which(frags_df$Level==x & frags_df$Burn2=="3"),])
    )})
frLVL$upperCalcined <- frLVL$muCalcined + 1.96*frLVL$seCalcined
frLVL$lowerCalcined <- frLVL$muCalcined - 1.96*frLVL$seCalcined
frLVL$FFImean <- sapply(frLVL$Level, function(z){
  if(z %in% FFI_df$Level){
    return(mean(FFI_df$FFI[which(FFI_df$Level==z)]))
  }else{return(NA)}
})

frLVL$FFImean <- sapply(frLVL$Level, function(z){
  if(z %in% FFI_df$Level){
    return(mean(FFI_df$FFI[which(FFI_df$Level==z)]))
  }else{return(NA)}
})

frLVL$FFIlower <- sapply(frLVL$Level, function(z){
  if(length(FFI_df$Level[which(FFI_df$Level==z)])>10){
    
    l <- FFI_df$FFI[which(FFI_df$Level==z)]
    
    lmeans <- sapply(1:1e3, function(v){
      return(mean(sample(l, size=length(l), replace=TRUE)))
    }) 
    
    return(quantile(lmeans, probs=0.025))
    
  }else{return(0)}
})


frLVL$FFIupper <- sapply(frLVL$Level, function(z){
  if(length(FFI_df$Level[which(FFI_df$Level==z)])>10){
    
    l <- FFI_df$FFI[which(FFI_df$Level==z)]
    
    lmeans <- sapply(1:1e3, function(v){
      return(mean(sample(l, size=length(l), replace=TRUE)))
    }) 
    
    return(quantile(lmeans, probs=0.975))
    
  }else{return(2)}
})


# Remove any levels in frLVL that contain NAs
frLVL <- frLVL[2:42,]

# Melt data frame for plotting polygons
frMELT <- melt(frLVL[,c(1,8,9,12,13,16,17)], id.vars="Level")
frMELT$bgroup <- substr(frMELT$variable, 6, length(frMELT$variable))
frMELT$variable <- substr(frMELT$variable, 1, 5)
frMELT <- frMELT[complete.cases(frMELT),]
# Order frMELT for plotting
for(y in c("Unburned","Carbonized","Calcined")){
  burntemp <- frMELT[which(frMELT$bgroup==y & frMELT$variable=="lower"),]
  frMELT[which(frMELT$bgroup==y & frMELT$variable=="lower"),] <- burntemp[rev(rownames(burntemp)),]
}

# Levels for plotting
ys <- unique(frags_df$Level)

#### PLOTS ####
###############

# Plot fragmentation
fragplot <- ggplot(data=frags_df, aes(x=MaxD, y=ljit))+
  geom_point(aes(color=factor(Burn2)), alpha=0.7, shape=16, size=0.85)+
  annotate("segment", x=rep(-Inf,44), xend=rep(Inf,44), 
           y=seq(from=84.5, to=127.5, by=1), yend=seq(from=84.5, to=127.5, by=1),
           color=rep("white",44), size=rep(0.5,44))+
  annotate("text",label="a", size=6, x=90.4, y=85.2, vjust=0)+
  labs(x="Max dimension (mm)", y="Level")+
  scale_x_continuous(limits=c(4,100), expand=c(0,0), breaks=seq(from=10,to=90,by=10))+
  scale_y_reverse(expand=c(0,0), limits=c(max(ys)+0.5, min(ys)-0.5), breaks=seq(130, 80, -2))+
  scale_colour_manual(values = c(c("0"="tan4","1"="black","3"="cadetblue4")))+
  theme(panel.grid=element_blank(), legend.position="none", 
        plot.background=element_blank(),
        panel.border=element_rect(fill=NA, color="grey", size=1),
        plot.margin=unit(c(0.1,-0.1,0.1,0.1), "cm"))

mufragplot <- ggplot(data=frLVL, aes(y=Level))+
  annotate("segment", x=rep(-Inf,44), xend=rep(Inf,44), 
           y=seq(from=84.5, to=127.5, by=1), 
           yend=seq(from=84.5, to=127.5, by=1),
           color=rep("white",44), size=rep(0.5,44))+
  geom_polygon(data=frMELT, aes(x=value, y=Level, fill=bgroup), alpha=0.2)+
  geom_path(aes(x=muUnburned), color="tan4", alpha=0.6, size=1)+
  geom_path(aes(x=muCarbonized), color="black", alpha=0.6, size=1)+
  geom_path(aes(x=muCalcined), color="cadetblue4", alpha=0.6, size=1)+
  labs(x="Mean max dimension (mm)")+
  annotate("text",label="b", size=6,x=29.1, y=85.2, vjust=0)+
  scale_x_continuous(expand=c(0,0), limits=c(3,32), 
                     breaks=seq(from=5, to=30, by=5))+
  scale_y_reverse(expand=c(0,0), limits=c(max(ys)+0.5, min(ys)-0.5))+
  scale_fill_manual(values = c("Unburned"="tan4","Carbonized"="black",
                               "Calcined"="cadetblue4"))+
  theme(panel.grid=element_blank(), legend.position="none", 
        plot.background=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border=element_rect(fill=NA, color="grey", size=1),
        plot.margin=unit(c(0.1,0.1,0.1,-0.1), "cm"))

propfragplot <- ggplot(data=frLVL, aes(ymin=Level-0.5, ymax=Level+0.5))+
  geom_rect(aes(xmin=rep(0, 41), xmax=propUnburned),fill="tan4",alpha=0.8)+
  geom_rect(aes(xmin=propUnburned, xmax=propUnburned+propCarbonized),
            fill="black",alpha=0.8)+
  geom_rect(aes(xmin=1-propCalcined, xmax=rep(1,41)),
            fill="cadetblue4",alpha=0.8)+
  annotate("segment", x=rep(-Inf,44), xend=rep(Inf,44), 
           y=seq(from=84.5, to=127.5, by=1), yend=seq(from=84.5, to=127.5, by=1),
           color=rep("white",44), size=rep(0.5,44))+
  annotate("text",label="c",size=6, x=0.9, y=85.2, vjust=0)+
  labs(x="Proportion burn stage")+
  scale_x_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(from=0.2, to=0.8, by=0.2))+
  scale_y_reverse(expand=c(0,0), limits=c(max(ys)+0.5, min(ys)-0.5))+
  theme(panel.grid=element_blank(), legend.position="none", 
        plot.background=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border=element_rect(fill=NA, color="grey", size=1),
        plot.margin=unit(c(0.1,0.1,0.1,-0.2), "cm"))

ffiplot <- ggplot(data=frLVL)+
  annotate("segment", x=rep(-Inf,44), xend=rep(Inf,44), 
           y=seq(from=84.5, to=127.5, by=1), yend=seq(from=84.5, to=127.5, by=1),
           color=rep("white",44), size=rep(0.5,44))+
  geom_rect(aes(xmin=FFIlower, xmax=FFIupper, 
                ymin=Level-0.5, ymax=Level+0.5),
            fill="red", alpha=0.2)+
  geom_segment(aes(y=Level-0.5, yend=Level+0.5, 
                   x=FFImean, xend=FFImean), color="red", size=0.6)+
  geom_point(data=FFI_df, aes(x=xjit, y=ljit), alpha=0.5, size=0.6, shape=16)+
  annotate("text",label="d",size=6, x=1.7, y=85.2, vjust=0)+
  labs(x="FFI")+
  scale_x_continuous(expand=c(0,0), limits=c(-0.1, 2.1), 
                     breaks=seq(from=0, to=2, by=0.5))+
  scale_y_reverse(expand=c(0,0), limits=c(max(ys)+0.5, min(ys)-0.5))+
  theme(panel.grid=element_blank(), legend.position="none", 
        plot.background=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border=element_rect(fill=NA, color="grey", size=1),
        plot.margin=unit(c(0.1,0.5,0.1,-0.2), "cm"))

frag_plots <- fragplot + mufragplot + propfragplot + ffiplot + 
  plot_layout(nrow=1, widths=c(5,2,2,2))

if(!dir.exists("Figures")) dir.create("Figures")
ggsave("Figures/Vertical_Patterns.pdf", plot=frag_plots,
       device="pdf", width=15, height=7, units="in")


######################################################################
########------- Horizontal Fragmentaion Model -------#################
######################################################################

# Create unit strings
frags_df$Unit <- paste(floor(frags_df$N), floor(frags_df$E))
frags_df$burn_ind <- coerce_index(factor(frags_df$Burn2))
frags_df$unit_ind <- coerce_index(factor(frags_df$Unit))
frags_df$llMaxD <- log(log(frags_df$MaxD))

# Create data frame with one row per level/unit
unitlevellist <- lapply(unique(frags_df$Level), function(v){
  
  levelstats <- t(sapply(unique(frags_df$Unit), function(z){
    
    if(v %in% frags_df$Level[which(frags_df$Unit==z)]){
      
      temp <- frags_df[which(frags_df$Unit==z & frags_df$Level==v),]
      
      if(nrow(temp[which(temp$burn_ind==1),])>5){
        
        u <- temp$MaxD[which(temp$burn_ind==1)]
        UnburnedMean <- mean(u)
        UnburnedSE <- sd(u)/sqrt(length(u))
        UnburnedLogmean <- mean(log(u))
        UnburnedLogSE <- sd(log(u))/sqrt(length(u))
        
      }else{
        
        UnburnedMean <- UnburnedSE <- 
          UnburnedLogmean <- UnburnedLogSE <- NA
        
      }
      
      if(nrow(temp[which(temp$burn_ind==2),])>5){
        
        b <- temp$MaxD[which(temp$burn_ind==2)]
        BurnedMean <- mean(b)
        BurnedSE <- sd(b)/sqrt(length(b))
        BurnedLogmean <- mean(log(b))
        BurnedLogSE <- sd(log(b))/sqrt(length(b))
        
      }else{
        
        BurnedMean <- BurnedSE <- 
          BurnedLogmean <- BurnedLogSE <- NA
        
      }
      
      if(nrow(temp[which(temp$burn_ind==3),])>5){
        
        c <- temp$MaxD[which(temp$burn_ind==3)]
        CalcinedMean <- mean(c)
        CalcinedSE <- sd(c)/sqrt(length(c))
        CalcinedLogmean <- mean(log(c))
        CalcinedLogSE <- sd(log(c))/sqrt(length(c))
        
      }else{
        
        CalcinedMean <- CalcinedSE <- 
          CalcinedLogmean <- CalcinedLogSE <- NA
        
      }
      
      return(c(UnburnedMean, UnburnedSE, UnburnedLogmean, 
               UnburnedLogSE, BurnedMean, BurnedSE, BurnedLogmean,
               BurnedLogSE, CalcinedMean, CalcinedSE, CalcinedLogmean,
               CalcinedLogSE, z))
    }else{return(c(rep(NA,12), z))}
    
    
  }))
  
  return(cbind(levelstats, rep(v, nrow(levelstats))))
  
})


modeldata <- data.frame(do.call("rbind", unitlevellist), stringsAsFactors=FALSE)
colnames(modeldata) <- c("UnburnedMean", "UnburnedSE", "UnburnedLogmean", 
                         "UnburnedLogSE", "BurnedMean", "BurnedSE", "BurnedLogmean",
                         "BurnedLogSE", "CalcinedMean", "CalcinedSE", "CalcinedLogmean",
                         "CalcinedLogSE", "Unit", "Level")
for(i in c(1:12, 14)){
  modeldata[,i] <- as.numeric(modeldata[,i])
}
modeldata$North <- sapply(modeldata$Unit, function(y) as.numeric(strsplit(y, " ")[[1]][1]))
modeldata$East <- sapply(modeldata$Unit, function(y) as.numeric(strsplit(y, " ")[[1]][2]))

u_df <- modeldata[!is.na(modeldata$UnburnedMean),]
b_df <- modeldata[!is.na(modeldata$BurnedMean),]
c_df <- modeldata[!is.na(modeldata$CalcinedMean),]

# Create datasets of logged and standardized dimension data
# and standardized provenience data
frags_df2 <- frags_df
frags_df2$MaxD <- log(frags_df2$MaxD)
frags_df2$Level <- (frags_df2$Level - mean(frags_df2$Level))/
  sd(frags_df2$Level)
frags_df2$N <- (frags_df2$N - mean(frags_df2$N))/
  sd(frags_df2$N)
frags_df2$E <- (frags_df2$E - mean(frags_df2$E))/
  sd(frags_df2$E)
frags_df2$MaxD <- (frags_df2$MaxD - mean(frags_df2$MaxD))/
  sd(frags_df2$MaxD)

frags_unburned <- list(MaxD=frags_df2$MaxD[which(frags_df2$Burn==0)],
                       Level=frags_df2$Level[which(frags_df2$Burn==0)],
                       No=frags_df2$N[which(frags_df2$Burn==0)],
                       Ea=frags_df2$E[which(frags_df2$Burn==0)])
frags_carbonized <- list(MaxD=frags_df2$MaxD[which(frags_df2$Burn==1)],
                       Level=frags_df2$Level[which(frags_df2$Burn==1)],
                       No=frags_df2$N[which(frags_df2$Burn==1)],
                       Ea=frags_df2$E[which(frags_df2$Burn==1)])
frags_calcined <- list(MaxD=frags_df2$MaxD[which(frags_df2$Burn==3)],
                         Level=frags_df2$Level[which(frags_df2$Burn==3)],
                         No=frags_df2$N[which(frags_df2$Burn==3)],
                         Ea=frags_df2$E[which(frags_df2$Burn==3)])

# Fit three spatial models for fragmentation (max dimension),
# one for each burning stage.
Model_unburned <- map2stan(alist(
  MaxD ~ dnorm(mu, sigma),
  mu <- a + bl*Level + bn*No + be*Ea,
  c(a, bl, bn, be) ~ dnorm(0, 1),
  sigma ~ dexp(1)),
  data=frags_unburned,
  chains=4, cores=4, iter=4e3, warmup=2e3)
Model_carbonized <- map2stan(alist(
  MaxD ~ dnorm(mu, sigma),
  mu <- a + bl*Level + bn*No + be*Ea,
  c(a, bl, bn, be) ~ dnorm(0, 1),
  sigma ~ dexp(1)),
  data=frags_carbonized,
  chains=4, cores=4, iter=4e3, warmup=2e3)
Model_calcined <- map2stan(alist(
  MaxD ~ dnorm(mu, sigma),
  mu <- a + bl*Level + bn*No + be*Ea,
  c(a, bl, bn, be) ~ dnorm(0, 1),
  sigma ~ dexp(1)),
  data=frags_calcined,
  chains=4, cores=4, iter=4e3, warmup=2e3)



######################################################################
############### PLOT FRAG AND BURNING ACROSS UNITS ###################
######################################################################

# Create data frame with one row per level/unit

fd$Unit <- paste(floor(fd$N), floor(fd$E))

unitstats <- t(sapply(unique(frags_df$Unit), function(z){
  
  temp <- frags_df[which(frags_df$Unit==z),]
  
  if(nrow(temp[which(temp$burn_ind==1),])>5){
    
    u <- temp$MaxD[which(temp$burn_ind==1)]
    UnburnedMean <- mean(u)
    UnburnedSE <- sd(u)/sqrt(length(u))
    UnburnedLogmean <- mean(log(u))
    UnburnedLogSE <- sd(log(u))/sqrt(length(u))
    
  }else{
    
    UnburnedMean <- UnburnedSE <- 
      UnburnedLogmean <- UnburnedLogSE <- NA
    
  }
  
  if(nrow(temp[which(temp$burn_ind==2),])>5){
    
    b <- temp$MaxD[which(temp$burn_ind==2)]
    BurnedMean <- mean(b)
    BurnedSE <- sd(b)/sqrt(length(b))
    BurnedLogmean <- mean(log(b))
    BurnedLogSE <- sd(log(b))/sqrt(length(b))
    
  }else{
    
    BurnedMean <- BurnedSE <- 
      BurnedLogmean <- BurnedLogSE <- NA
    
  }
  
  if(nrow(temp[which(temp$burn_ind==3),])>5){
    
    c <- temp$MaxD[which(temp$burn_ind==3)]
    CalcinedMean <- mean(c)
    CalcinedSE <- sd(c)/sqrt(length(c))
    CalcinedLogmean <- mean(log(c))
    CalcinedLogSE <- sd(log(c))/sqrt(length(c))
    
  }else{
    
    CalcinedMean <- CalcinedSE <- 
      CalcinedLogmean <- CalcinedLogSE <- NA
    
  }
  
  
  
  burnstats <- c(UnburnedMean, UnburnedSE, UnburnedLogmean, 
                 UnburnedLogSE, BurnedMean, BurnedSE, BurnedLogmean,
                 BurnedLogSE, CalcinedMean, CalcinedSE, CalcinedLogmean,
                 CalcinedLogSE)
  
  if(z %in% paste(floor(FFI_df$N), floor(FFI_df$E))){
    ffi <- mean(FFI_df$FFI[which(paste(floor(FFI_df$N), 
                                       floor(FFI_df$E))==z)])
  }else{ffi <- NA}
  
  nisp <- sum(fd$count[which(fd$Unit==z)])
  
  ubcount <- sum(fd$count[which(!(fd$Element %in% nonbone) & 
                                  fd$Unit==z &
                                  fd$Burn==0)])
  
  bcount <- sum(fd$count[which(!(fd$Element %in% nonbone) & 
                                 fd$Unit==z &
                                 fd$Burn!=0)])
  
  exvol <- length(unique(fd$Level[which(fd$Unit==z &
                                          fd$Level>50)]))/20
  
  return(c(burnstats, ffi, nisp, ubcount, bcount, exvol, z))
  
}))

unitstats <- rbind(unitstats, c(rep(NA, 17), "1009 998"))
unitstats <- as.data.frame(unitstats, stringsAsFactors = FALSE)
for(j in 1:(ncol(unitstats)-1)) unitstats[,j] <- as.numeric(unitstats[,j])
colnames(unitstats) <- c("UnburnedMean", "UnburnedSE", "UnburnedLogmean", 
                         "UnburnedLogSE", "BurnedMean", "BurnedSE", "BurnedLogmean",
                         "BurnedLogSE", "CalcinedMean", "CalcinedSE", "CalcinedLogmean",
                         "CalcinedLogSE", "FFImean", "NISP", "NonBurnedNISP", "BurnedNISP",
                         "ExVol", "Unit")

unitstats$North <- sapply(unitstats$Unit, function(y) as.numeric(strsplit(y, " ")[[1]][1]))
unitstats$East <- sapply(unitstats$Unit, function(y) as.numeric(strsplit(y, " ")[[1]][2]))
unitstats$PropBurned <- unitstats$BurnedNISP/(unitstats$BurnedNISP + unitstats$NonBurnedNISP)
unitstats$MaterialDensity <- unitstats$NISP/unitstats$ExVol


Plot_Unit_Ufrag <- ggplot(data=unitstats, aes(x=East+0.5, y=North+0.5))+
  geom_tile(aes(fill=UnburnedMean), color="black")+
  geom_text(aes(label=round(UnburnedMean,2)))+
  scale_fill_gradient(low="white", high="blue", limits=c(8,22))+
  annotate("text", label="Not burned", y=1009.5, x=1000.5, vjust=0, hjust=0)+
  annotate("text", label="a", x=995.5, y=1009.5, size=5)+
  scale_y_continuous(breaks=seq(1007, 1010, 1))+
  labs(y="North")+
  theme(axis.line=element_blank(), axis.text.x=element_blank(),
        panel.background=element_blank(), panel.grid=element_blank(),
        axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="none")

Plot_Unit_Bfrag <- ggplot(data=unitstats, aes(x=East+0.5, y=North+0.5))+
  geom_tile(aes(fill=BurnedMean), color="black")+
  geom_text(aes(label=round(BurnedMean,2)))+
  scale_fill_gradient(low="white", high="blue", limits=c(8,22), name="mm")+
  annotate("text", label="Carbonized", y=1009.5, x=1000.5, vjust=0, hjust=0)+
  annotate("text", label="b", x=995.5, y=1009.5, size=5)+
  scale_y_continuous(breaks=seq(1007, 1010, 1))+
  labs(y="North")+
  theme(axis.line=element_blank(), axis.text.x=element_blank(),
        panel.background=element_blank(), panel.grid=element_blank(),
        axis.title.x=element_blank(), axis.ticks.x=element_blank())

Plot_Unit_Cfrag <- ggplot(data=unitstats, aes(x=East+0.5, y=North+0.5))+
  geom_tile(aes(fill=CalcinedMean), color="black")+
  geom_text(aes(label=round(CalcinedMean,2)))+
  scale_fill_gradient(low="white", high="blue", limits=c(8,22))+
  annotate("text", label="Calcined", y=1009.5, x=1000.5, vjust=0, hjust=0)+
  annotate("text", label="c", x=995.5, y=1009.5, size=5)+
  scale_y_continuous(breaks=seq(1007, 1010, 1))+
  scale_x_continuous(breaks=seq(995, 1004, 1))+
  labs(y="North", x="East")+
  theme(axis.line=element_blank(), panel.background=element_blank(), 
        panel.grid=element_blank(),legend.position="none")

Plot_Unit_Propburn <- ggplot(data=unitstats, aes(x=East+0.5, y=North+0.5))+
  geom_tile(aes(fill=PropBurned), color="black")+
  geom_text(aes(label=round(PropBurned,2)))+
  annotate("text", label="d", x=995.5, y=1009.5, size=5)+
  scale_fill_gradient(low="white", high="orange", limits=c(0,1),
                      name="Prop. Burned")+
  theme(axis.line=element_blank(), axis.text=element_blank(),
        panel.background=element_blank(), panel.grid=element_blank(),
        axis.title=element_blank(), axis.ticks=element_blank())

Plot_Unit_FFI <- ggplot(data=unitstats, aes(x=East+0.5, y=North+0.5))+
  geom_tile(aes(fill=FFImean), color="black")+
  geom_text(aes(label=round(FFImean,2)))+
  annotate("text", label="e", x=995.5, y=1009.5, size=5)+
  scale_fill_gradient(low="white", high="forestgreen", limits=c(0,2),
                      name="FFI")+
  theme(axis.line=element_blank(), axis.text=element_blank(),
        panel.background=element_blank(), panel.grid=element_blank(),
        axis.title=element_blank(), axis.ticks=element_blank())

Plot_Unit_dens <- ggplot(data=unitstats, aes(x=East+0.5, y=North+0.5))+
  geom_tile(aes(fill=MaterialDensity), color="black")+
  geom_text(aes(label=round(MaterialDensity,0)))+
  annotate("text", label="f", x=995.5, y=1009.5, size=5)+
  scale_fill_gradient(low="white", high="purple", limits=c(30,4800),
                      name="NISP/m")+
  scale_x_continuous(breaks=seq(995, 1004, 1))+
  labs(x="East")+
  theme(axis.line=element_blank(), axis.text.y=element_blank(),
        panel.background=element_blank(), panel.grid=element_blank(),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())


horiz_frag_plot <- (Plot_Unit_Ufrag / Plot_Unit_Bfrag / Plot_Unit_Cfrag) | 
                   (Plot_Unit_Propburn / Plot_Unit_FFI / Plot_Unit_dens)

if(!dir.exists("Figures")) dir.create("Figures")
ggsave("Figures/Horizontal_Patterns.pdf", device="pdf",
       width=13, height=6, units="in")

######################################################################
###########--- MNIs for lagomorphs and artiodctyls ---################
######################################################################

ArtioDF <- fd[which(fd$Order=="Artiodactyla" & fd$count==1 & fd$Size=="m"),]
LagoDF <- fd[which(fd$Order=="Lagomorpha" & fd$count==1),]

# Collapse scan site columns for each taxon
ArtioDF <- lapply(seq(from=15, to=18, by=1), function(x){
  df <- ArtioDF[c(seq(from=1, to=14, by=1), x)]
  colnames(df)[ncol(df)] <- "Scan"
  return(df)
})
ArtioDF <- do.call("rbind", ArtioDF)
ArtioDF$SidedScans <- paste(ArtioDF$Scan, ArtioDF$Element, ArtioDF$Side, sep="")
ArtioTable <- data.frame(table(ArtioDF$SidedScans), stringsAsFactors=FALSE)

LagoDF <- lapply(seq(from=15, to=18, by=1), function(x){
  df <- LagoDF[c(seq(from=1, to=14, by=1), x)]
  colnames(df)[ncol(df)] <- "Scan"
  return(df)
})
LagoDF <- do.call("rbind", LagoDF)
LagoDF$SidedScans <- paste(LagoDF$Scan, LagoDF$Element, LagoDF$Side, sep="")
LagoTable <- data.frame(table(LagoDF$SidedScans), stringsAsFactors=FALSE)


##################################################################
#########------- Taxonomic Richness -------#######################
##################################################################

# FUNCTION: Chao1 estimator
C1 <- function(m){
  
  if((1 %in% m[,1]) & (nrow(m) > 2)){
    
    Sobs <- nrow(m)
    n <- sum(m[,2]*m[,1])
    k <- 1-1/n
    f1 <- m[which(m[,1]==1),2]
    
    if(2 %in% m[,1]){
      
      f2 <- m[which(m[,1]==2),2]
      Shat <- Sobs + ((n-1)/n)*((f1^2)/(2*f2))
      VarS <- f2*((k/2*(f1/f2)^2) + ((k^2)*(f1/f2)^3) +
                    (1/4)*(k^2)*(f1/f2)^4)
      
    }else{
      
      Shat <- Sobs + ((n-1)/n)*f1*(f1-1)/2
      VarS <- ((k*f1*(f1-1))/2) + ((k^2*f1*(2*f1-1)^2)/4) -
        ((k^2*f1^4)/(4*Shat))
      
    }
    
    C <- exp(1.96*(log(1+VarS/(Shat-Sobs)^2))^(1/2))
    LB <- Sobs+(Shat-Sobs)/C
    UB <- Sobs+(Shat-Sobs)*C
    
    if(is.nan(UB)) Shat <- LB <- UB <- NA
    
    return(c(Sobs, Shat, LB, UB, n))
  } else {
    return(c(nrow(m), rep(NA,3), sum(m[,2]*m[,1])))
  }
}

# Obtain family level taxonomic richness values by level--convert
# to table for Chao estimator function in SPECIES package.

# Select only those cases in fd with family level IDs
fdFAMILY <- fd[which(fd$Family != "-"),]

# Levels for following apply function
LVLs <- data.frame(u=frLVL$Level[frLVL$Level%%2==1])
LVLs$l <- LVLs$u+1

# Return matrix of Chao CIs, one for each level
ChaoCIs <- data.frame(t(sapply(1:nrow(LVLs), function(x){
  
  # Find unique families for level x
  families <- unique(fdFAMILY$Family[fdFAMILY$Level %in% LVLs[x,]])
  fdFAMILYsub <- fdFAMILY[which(fdFAMILY$Level %in% LVLs[x,]),]
  
  # Sum NISP values for each family z in level x
  CountsFamilies <- sapply(families, function(z){
    sum(fdFAMILYsub$count[which(fdFAMILYsub$Family==z)])
  })
  
  if(length(CountsFamilies)>1){
    
    # Create 'frequency of frequencies' matrix
    FoF <- as.data.frame(table(CountsFamilies), stringsAsFactors=FALSE)
    FoF <- as.matrix(cbind(as.numeric(FoF[,1]), FoF[,2]))
    
    # If counts are appropriate, return Chao estimator CI
    return(C1(FoF))
    
  }else{if(length(families)==1){
    return(c(1,rep(NA,3),sum(fdFAMILYsub$count)))}else(return(rep(NA,5)))}
  
})))
colnames(ChaoCIs) <- c("ObsRichness", "Chao1", "Chao1LB", "Chao1UB", "NISP")
ChaoCIs$Level <- paste(LVLs[,1], "-", LVLs[,2])
ChaoCIs$Index <- seq(from=nrow(ChaoCIs), to=1, by=-1)

# Plot Richness measures
Plot_R_obs <- ggplot(data=ChaoCIs, aes(y=ObsRichness, x=Index))+
  geom_line(color="blue", size=0.5)+
  geom_point(color="blue", size=1)+
  scale_x_continuous(breaks=ChaoCIs$Index, labels=ChaoCIs$Level,
                     limits=c(min(ChaoCIs$Index)-0.5, 
                              max(ChaoCIs$Index)+0.5),
                     expand=c(0,0))+
  annotate("text", label="a", x=1, y=4.4, size=5, vjust=0)+
  labs(y="Observed richness", x="Level")+
  coord_flip()+
  theme(panel.background=element_blank(), 
        axis.line.x=element_line(color="grey"),
        axis.line.y=element_line(color="grey"),
        axis.ticks.y=element_line(color="grey"),
        axis.ticks.x=element_line(color="grey"))

Plot_R_chao <- ggplot(data=ChaoCIs, aes(y=Chao1, x=Index))+
  geom_point(color="red", size=1)+
  geom_segment(aes(y=Chao1LB, yend=Chao1UB, xend=Index), 
               size=0.75, color="red", alpha=0.5)+
  scale_x_continuous(breaks=ChaoCIs$Index, labels=ChaoCIs$Level,
                     limits=c(min(ChaoCIs$Index)-0.5, 
                              max(ChaoCIs$Index)+0.5),
                     expand=c(0,0))+
  annotate("text", label="c", x=1, y=25, size=5, vjust=0)+
  labs(y="Chao1 estimator")+
  coord_flip()+
  theme(panel.background=element_blank(), 
        axis.line.x=element_line(color="grey"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(color="grey"))

Plot_R_nisp <- ggplot(data=ChaoCIs, aes(y=NISP, x=Index))+
  geom_bar(fill="purple", stat="identity")+
  scale_x_continuous(breaks=ChaoCIs$Index, labels=ChaoCIs$Level,
                     limits=c(min(ChaoCIs$Index)-0.5, 
                              max(ChaoCIs$Index)+0.5),
                     expand=c(0,0))+
  annotate("text", label="b", x=1, y=120, size=5, vjust=0)+
  labs(y="NISP")+
  coord_flip()+
  theme(panel.background=element_blank(), 
        axis.line.x=element_line(color="grey"),        
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(color="grey"))

Richness_plot <- Plot_R_obs + Plot_R_nisp + Plot_R_chao + 
                 plot_layout(ncol=3, widths=c(3/8, 2/8, 3/8))

if(!dir.exists("Figures")) dir.create("Figures")
ggsave("Figures/Richness_Patterns.pdf", plot=Richness_plot,
       device="pdf", width=8, height=5, units="in")

#########################################
########### Burning by taxon ############
Burn_df <- fd[!(fd$Element %in% nonbone),]

Burn_order <- unique(Burn_df$Order)

Burn_taxa <- t(sapply(Burn_order, function(z){
  
  nb <- sum(Burn_df$count[which(Burn_df$Order==z & Burn_df$Burn==0)])
  b <- sum(Burn_df$count[which(Burn_df$Order==z & (Burn_df$Burn==1 | Burn_df$Burn==2))])
  c <- sum(Burn_df$count[which(Burn_df$Order==z & Burn_df$Burn==3)])
  
  pnb <- nb/(nb + b + c)
  pb <- b/(nb + b + c)
  pc <- c/(nb + b + c)
  
  return(c(nb, pnb, b, pb, c, pc))
}))

Burn_taxadf <- data.frame(Taxon=Burn_order, 
                          nb_NISP=Burn_taxa[,1], 
                          nb_per=Burn_taxa[,2],
                          b_NISP=Burn_taxa[,3],
                          b_per=Burn_taxa[,4],
                          c_NISP=Burn_taxa[,5],
                          c_per=Burn_taxa[,6])
