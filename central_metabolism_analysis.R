library(dplyr)
library(tidyr)
library(readxl)
library(randomForest)

cm_data <- read_excel("~/Documents/Misc Projects/Albert_Pollen/Data/Updated_Ideas/PollenBiomarker_DATASET_ForBryan.xlsx",sheet="CentralMetabolism")

#Get the data ready for analysis (transpose, shorten names, make factors xnumeric)
cnames_orig <- cm_data$`row identity (all IDs)`
cnames <- unlist(lapply(strsplit(cnames_orig,"_"),function(x)return(x[1])))
cm_data <- data.frame(t(cm_data[,-1]))
colnames(cm_data) <- cnames
rnames <- rownames(cm_data)
rownames(cm_data) <- NULL
rnames <- gsub(".raw baseline-corrected Peak area","",rnames)
rnames <- gsub("LCMS_","",rnames)
cm_data$Sample <- rnames
for(i in 4:(ncol(cm_data)-1)){
  cm_data[,i] <- as.numeric(as.character(cm_data[,i]))
}

#Monosaccharides column cant be duplicated
colnames(cm_data)[which(colnames(cm_data)=="Monosaccharides")] <- paste("Monosaccharides",unlist(lapply(strsplit(cnames_orig[which(colnames(cm_data)=="Monosaccharides")],"_"),function(x)return(x[2]))),sep="_")

##-----------------------##
# Get the data objects ready for filtering, modelling

pollen_dat <- cm_data%>%filter(SAMPLE.TYPE=="Pollen")
filter_dat <- cm_data%>%filter(SAMPLE.TYPE!="Pollen")
xmat_pollen <- pollen_dat%>%dplyr::select(`2-Deoxy-D Ribose`:Tryptophan)%>%data.matrix()
rownames(xmat_pollen) <- pollen_dat$ID2
xmat_filter <- filter_dat%>%dplyr::select(`2-Deoxy-D Ribose`:Tryptophan)%>%data.matrix()

yvec_pollen <- as.character(pollen_dat$CELL.FACTOR)
yvec_filter <- gsub("Peak","",filter_dat$CELL.FACTOR)


##----------------------##
#Scale the data
maxes <- matrix(apply(xmat_pollen,2,max,na.rm=T),nrow=nrow(xmat_pollen),ncol=ncol(xmat_pollen),byrow=TRUE)
mins <- matrix(apply(xmat_pollen,2,min,na.rm=T),nrow=nrow(xmat_pollen),ncol=ncol(xmat_pollen),byrow=TRUE)
xmat_pollen_scaled <- (xmat_pollen-mins)/(maxes-mins)

maxes <- matrix(apply(xmat_filter,2,max,na.rm=T),nrow=nrow(xmat_filter),ncol=ncol(xmat_filter),byrow=TRUE)
mins <- matrix(apply(xmat_filter,2,min,na.rm=T),nrow=nrow(xmat_filter),ncol=ncol(xmat_filter),byrow=TRUE)
xmat_filter_scaled <- (xmat_filter-mins)/(maxes-mins)


#Replace missing values with zeros
xmat_pollen_scaled[is.na(xmat_pollen_scaled)] <- 0
xmat_filter_scaled[is.na(xmat_filter_scaled)] <- 0

##----------------------##
# Modeling options
set.seed(30)
#Random forest
rf_fit <- randomForest(x=xmat_pollen_scaled,y=as.factor(yvec_pollen),xtest=xmat_filter_scaled,ytest=as.factor(yvec_filter),importance=TRUE,ntree=1000,mtry=10,keep.forest = TRUE)
rf_fit

#Make plot for ambient samples (Figure 3)
amb_pred <- predict(rf_fit, newdata=xmat_filter_scaled,type='prob')
amb_p_orig <- data.frame(amb_pred,Truth=yvec_filter,Sample=as.character(1:13))
amb_p <- amb_p_orig%>%gather(Prediction,Prob,-Truth,-Sample)
amb_p$Sample <- factor(amb_p$Sample,levels=as.character(c(1:13)))
qplot(Sample,Prob,colour=Prediction,data=amb_p)+geom_line(aes(group=Prediction))+theme_bw()+#theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab("Sample Number")+ylab("Probability of Each Species")+geom_vline(aes(xintercept=3.5),colour='gray50')+geom_vline(aes(xintercept=10.5),colour='gray50')


