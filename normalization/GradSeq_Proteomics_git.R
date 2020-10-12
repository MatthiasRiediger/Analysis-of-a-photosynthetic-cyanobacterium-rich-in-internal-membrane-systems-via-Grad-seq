setwd("C:/Users/matth/OneDrive/Dokumente/HessLab_Data_PhD_Backup_130320/projects/3_GradSeq_MS/Results/Proteomics/Analysis/Proteomics_summary/Raw_data_processing")

###### Quality control ######
##############################################################################################

require(reshape2)
require(ggplot2)
require(gridExtra)
require(lme4)


### Input = total MS intensities ###
input <- read.csv("GradSeqInput_proteomics.csv", sep=";", header=TRUE) # DF: ID, Repl, Fractionwise RC



  
x <- input
x <- cbind(x[2], x[,7:70])
names(x)[1] <- "ID"
x<- melt(x, id=c("ID"))
names(x)[2:3]<-c("Replicate_Fraction", "Intensity") #Minimum Readcount / ID and Fraction from all Replicates
Replicate_Fraction <- colsplit(as.factor(x$Replicate_Fraction), "_", c("Replicate", "Fraction"))
x <- cbind(x[1], Replicate_Fraction[1:2], x[3])
x <- droplevels(subset(x, Replicate != "Reference"))
  
  
### Pre-format unnormalized data: Filter features with <1 Intensities in all fractions ###
  
Min_Filter <- function(xPre){
  
  xmax<- aggregate(xPre$Intensity, by = list(xPre$ID, xPre$Replicate), max)
    names(xmax)[1:3]<-c("ID", "Replicate", "maxInt") #Maximum...
  xPre <- dcast(xPre, ID + Replicate ~ Fraction)   
  xPre <- merge(xPre, xmax, by=c("ID", "Replicate"))
  xPre <- subset(xPre, maxInt >=1)
  xPre <- xPre[,1:18]
  xPre <- melt(xPre, id=c("ID", "Replicate"))
    names(xPre)[3:4]<-c("Fraction", "Intensity")
  xPre <- dcast(xPre, ID + Fraction ~ Replicate)   
  xPre <- xPre[complete.cases(xPre), ]
  xPre <- melt(xPre, id=c("ID", "Fraction"))
    names(xPre)[3:4]<-c("Replicate", "Intensity")
xPre <<- xPre

}
Min_Filter(x) # Output:  xPre(Formatted x)
  
  
### Normalization to Median Intensity -> Not applied!!! ###

Normalization <- function(xNorm){ 
  
  xNorm$Intensity <- log2(xNorm$Intensity+1)
  xMed <- droplevels(subset(xNorm, Intensity >= .1))
  xMed<- aggregate(xMed$Intensity, by = list(xMed$Fraction, xMed$Replicate), median)
    names(xMed) <- c("Fraction", "Replicate", "MedInt")
  xMedMin <- aggregate(xMed$MedInt, by = list(xMed$Fraction), min)
    names(xMedMin) <- c("Fraction", "MedIntMin")
  xMed <- merge(xMed, xMedMin, by=c("Fraction"))
  xMed$MedInt <- xMed$MedIntMin/xMed$MedInt
    
  xNorm <- merge(xNorm, xMed[1:3], by=c("Fraction", "Replicate"))
  xNorm$Intensity <- xNorm$Intensity*xNorm$MedInt
  xNorm <- xNorm[,1:4]
  xNorm <- dcast(xNorm, ID + Fraction ~ Replicate)   
xNorm <<- xNorm

}#Output: xNorm
  Normalization(xPre)
  
  xNorm <- xPre
  xNorm$Intensity <- log2(xNorm$Intensity+1)
  xNorm <- dcast(xNorm, ID + Fraction ~ Replicate)   
  

### Compute lin Regr Coeff ### 
  
LinRegrCoeff <- function(x) {require(lme4)
    Repl12 <- coef(lmList(formula = Repl2 ~ Repl1 | Fraction, data= xNorm))
    Repl13 <- coef(lmList(formula = Repl3 ~ Repl1 | Fraction, data= xNorm))
    Repl23 <- coef(lmList(formula = Repl3 ~ Repl2 | Fraction, data= xNorm))
    names(Repl12)[1:2]<-c("Repl12_intercept", "Repl12")
    names(Repl13)[1:2]<-c("Repl13_intercept", "Repl13")
    names(Repl23)[1:2]<-c("Repl23_intercept", "Repl23")
    Repl12 <- data.frame(Subject=rownames(Repl12),Repl12,check.names=FALSE)
    Repl13 <- data.frame(Subject=rownames(Repl13),Repl13,check.names=FALSE)
    Repl23 <- data.frame(Subject=rownames(Repl23),Repl23,check.names=FALSE)
    rownames(Repl12) <- NULL 
    rownames(Repl13) <- NULL 
    rownames(Repl23) <- NULL 
    Rep_norm <- merge(Repl12, Repl13, by=c("Subject"))
    Rep_norm <- merge(Rep_norm, Repl23, by=c("Subject"))
    names(Rep_norm)[1]<-c("Fraction")
    Coeff <<- Rep_norm
  } #Output: Coeff
  LinRegrCoeff(xNorm)
  
  
### Regression and Boxplot of RC distribution between Replicate Fractions ###
 
  R2Rep12 <- function(xNorm){
    Repl12 = lm(Repl2 ~ Repl1, xNorm);
    eqRepl12 <- substitute(~~R^2~"="~r2, 
                           list(r2 = format(summary(Repl12)$r.squared, digits = 2)))
    as.character(as.expression(eqRepl12))
  }
  R2Rep13 <- function(xNorm){
    Repl13 = lm(Repl3 ~ Repl1, xNorm);
    eqRepl13 <- substitute(~~R^2~"="~r2, 
                           list(r2 = format(summary(Repl13)$r.squared, digits = 2)))
    as.character(as.expression(eqRepl13))
  }
  R2Rep23 <- function(xNorm){
    Repl23 = lm(Repl3 ~ Repl2, xNorm);
    eqRepl23 <- substitute(~~R^2~"="~r2, 
                           list(r2 = format(summary(Repl23)$r.squared, digits = 2)))
    as.character(as.expression(eqRepl23))
  }
  R2Rep12 <- by(xNorm, xNorm$Fraction, R2Rep12)
  R2Rep13 <- by(xNorm, xNorm$Fraction, R2Rep13)
  R2Rep23 <- by(xNorm, xNorm$Fraction, R2Rep23)
  R2Rep12 <- data.frame(eq = unclass(R2Rep12), Fraction = names(R2Rep12))  
  R2Rep13 <- data.frame(eq = unclass(R2Rep13), Fraction = names(R2Rep13))  
  R2Rep23 <- data.frame(eq = unclass(R2Rep23), Fraction = names(R2Rep23))  

    corrplot <- function (dat, xvar, yvar, R2df, x_name, y_name) {
    ggplot(dat,aes(x= xvar, y=yvar) , group=Fraction)+
      geom_point(alpha=0.25, size=0.1)+
      #geom_smooth(method='lm',formula=y~x)+
      geom_segment(x=15, y=15, xend=38, yend=38, linetype="dashed", size=0.2, colour="grey50")+  
      facet_grid(. ~ Fraction)+
      geom_text(data = R2df, aes(x = 22, y = 35, label = eq), 
                color = 'black',  parse = TRUE) +
      scale_x_continuous(name = paste("log2(norm. intensities [Replicate ", x_name, "])"), limits=c(12,38),  expand=c(0,0)) +
      scale_y_continuous(name = paste("log2(norm. intensities [Replicate ", y_name, "])"), limits=c(12,38), expand=c(0,0)) +
      theme_bw()
  } #InputDF: ID, Fraction, L2Repl1, L2Repl2, L2Repl3
    Rep12 <- corrplot(xNorm, xNorm[,3], xNorm[,4], R2Rep12, names(xNorm[3]), names(xNorm[4]))
    Rep13 <- corrplot(xNorm, xNorm[,3], xNorm[,5], R2Rep13, names(xNorm[3]), names(xNorm[5]))
    Rep23 <- corrplot(xNorm, xNorm[,4], xNorm[,5], R2Rep23, names(xNorm[4]), names(xNorm[5]))
 
      #pdf("Normalized_Fractionwise_linear_regression_transposed.pdf", width=20, height=6)
      grid.arrange(Rep12, Rep13, Rep23, nrow=3)
      #dev.off()
  
    xNorm <- melt(xNorm, id=c("ID", "Fraction"))
    names(xNorm)[3:4]<-c("Replicate", "Intensity") 
    xMed <- droplevels(subset(xNorm, Intensity > 0))
    boxP <- function (dat, xvar, yvar, xgroup) {
    ggplot(dat,aes(x= xgroup, y=yvar, fill=xvar))+
      geom_boxplot()+
      #facet_grid(xgroup ~ .)+
      #geom_segment(x=0, y=0, xend=25, yend=25)+  
      #scale_x_continuous(name = "",  expand=c(0,0)) +
      scale_y_continuous(name = paste("log2(norm. reads)"), expand=c(0,0)) +
      theme_bw()
  } #InputDF: ID, Fraction, Replicate, Log2RC
  
      #pdf("Normalized_fractionwise_intensity_distribution.pdf", width=12, height=6)
      boxP(xMed, xMed$Replicate, xMed$Intensity, xMed$Fraction)
      #dev.off()
  
  
### Correlation between Features/IDs across Replicates ###
  SpCorrelation <- function(xNorm){      
    require(reshape2)
    xNorm<- dcast(xNorm, ID + Fraction ~ Replicate)
    require(plyr)
    SpCor12 <- ddply(xNorm, ~ID, summarise, "corr" = cor(Repl1, Repl2, method = "spearman"))
    SpCor13 <- ddply(xNorm, ~ID, summarise, "corr" = cor(Repl1, Repl3, method = "spearman"))
    SpCor23 <- ddply(xNorm, ~ID, summarise, "corr" = cor(Repl2, Repl3, method = "spearman"))
    SpCor123 <- merge( merge( SpCor12, SpCor13, by=c("ID") ), SpCor23, by=c("ID") )
    names(SpCor123)[2:4]<-c("SpCor12", "SpCor13", "SpCor23")
    SpCor <- melt(SpCor123, id=c("ID"))
    names(SpCor)[2:3]<-c("Replicate", "Corr_coeff")
    SpCorAv <- aggregate(as.numeric(SpCor$Corr_coeff), by = list(SpCor$ID), mean)
    names(SpCorAv)[1:2]<-c("ID", "SpCorAv")
    SpCor <- merge( merge( merge(SpCor12, SpCor13, by=c("ID")), 
                           SpCor23, by=c("ID")), SpCorAv, by=c("ID"))
    names(SpCor)[2:4]<-c("SpCor12", "SpCor13", "SpCor23")
    SpCor <- melt(SpCor, id=c("ID"))
    names(SpCor)[2:3]<-c("Replicate", "Corr_coeff")
    SpCor <<- SpCor
  } #Output: SpCor
  SpCorrelation(xNorm)
  require(ggplot2)
  SpCorHis <- function(x, Corr_coeff, Replicate) {ggplot(x,aes(x= Corr_coeff, fill=Replicate))+
      geom_histogram(binwidth=0.01, position="stack", color="black") +
      facet_grid(Replicate ~ .)+
      scale_x_continuous(name = "Correlation Coefficient", limits=c(0,1), expand=c(0,0)) +
      scale_y_continuous(name = "Number of Genes", expand=c(0,0)) +
      theme_bw()+
      theme(panel.grid.major = element_blank())
  }
  boxP <- function (dat, xvar, yvar) {
    ggplot(dat,aes(x= xvar, y=yvar, fill=xvar))+
      geom_boxplot()+
      #facet_grid(xgroup ~ .)+
      #geom_segment(x=0, y=0, xend=25, yend=25)+  
      #scale_x_continuous(name = "",  expand=c(0,0)) +
      scale_y_continuous(name = paste("log2(norm. Intensity)"), expand=c(0,0)) +
      theme_bw()
  } #InputDF: ID, Fraction, Replicate, Log2RC
  #pdf("SpearmanCorrelationHistogram.pdf", width=6, height=6)
  SpCorHis(SpCor, SpCor$Corr_coeff, SpCor$Replicate)
  #dev.off()
  #pdf("SpearmanCorrelationBoxplot.pdf", width=12, height=6)
  boxP(SpCor, SpCor$Replicate, SpCor$Corr_coeff)
  #dev.off()
  
  
  
  ### Prepare Output ###
  SpCor<- dcast(SpCor, ID  ~ Replicate) 
  xNormmax <- dcast(xNorm, ID + Fraction ~ Replicate)
  xNormmax <- transform(xNormmax, AvInt = rowMeans(xNormmax[,3:5]))
  xNormmax <- aggregate(xNormmax$AvInt, by = list(xNormmax$ID), max)
  names(xNormmax)[1:2]<-c("ID", "maxAvInt")
  Output <<- merge(xNormmax, SpCor, by=c("ID"))
  
  corrplot <- function (dat, xvar, yvar) {ggplot(dat,aes(x= xvar, y=yvar) , group=Fraction)+
      geom_point(alpha=0.25, size=0.1)+
      geom_smooth()+
      scale_x_continuous(name = paste("Average Spearman Correlation Coefficient"),  expand=c(0,0)) +
      scale_y_continuous(name = paste("Average maximum Log2(Normalized Intensity)"), expand=c(0,0)) +
      theme_bw()
  } #InputDF: ID, maxAvInt, SpCorAv
  
  #pdf("Regression_AvSpCor2AvInt.pdf", width=6, height=6)
  corrplot(Output, Output$SpCorAv, Output$maxAvInt)
  #dev.off()
  


xNorm$Intensity <- (2^xNorm$Intensity)-1
x <- xNorm
#write.csv(x, "xNorm.csv")

xNorm <- aggregate(xNorm$Intensity, by = list(xNorm$ID, xNorm$Fraction), mean)  
names(xNorm)[1:3]<-c("ID","Fraction", "Intensity") 
xNorm <-  dcast(xNorm, ID  ~ Fraction)
xNorm <- merge(xNorm, Output, by="ID")

#xNorm[2:17] <- xNorm[2:17]/rowSums(xNorm[2:17]) # Relative Readcounts

write.csv(xNorm, "Summary_Output.csv")

x <- read.csv("output/xNorm.csv")

x <- dcast(x, ID ~ Replicate + Fraction)
write.csv(x, "WGCNA_input_total_proteomics.csv")
x[2:17] <- x[2:17]/rowSums(x[2:17]) # Relative Intensities Repl1
x[18:34] <- x[18:34]/rowSums(x[18:34]) # Relative Intensities Repl2
x[35:49] <- x[35:49]/rowSums(x[35:49]) # Relative Intensities Repl3
write.csv(x, "WGCNA_input_proteomics.csv")









