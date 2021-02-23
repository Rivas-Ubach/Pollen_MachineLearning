library(dplyr)
library(tidyr)
library(readxl)
library(mixOmics)

pas_data <- read_excel("~/Documents/Misc Projects/Albert_Pollen/Data/Updated_Ideas/PollenBiomarker_DATASET_ForBryan.xlsx",sheet="NEG Fingerprint")
#Get the data ready for analysis (transpose, shorten names, make factors xnumeric)
cnames_orig <- pas_data$`row identity (all IDs)`

cnames <- cnames_orig#unlist(lapply(strsplit(cnames_orig,"_"),function(x)return(x[1])))
pas_data <- data.frame(t(pas_data[,-1]))
colnames(pas_data) <- cnames
rnames <- rownames(pas_data)
rownames(pas_data) <- NULL
rnames <- gsub(".raw baseline-corrected Peak area","",rnames)
rnames <- gsub("LCMS_","",rnames)
pas_data$Sample <- rnames
for(i in 4:(ncol(pas_data)-1)){
  pas_data[,i] <- as.numeric(as.character(pas_data[,i]))
}

#Fix duplicated column names
dups <- names(which(table(colnames(pas_data))>1))
for(i in 1:length(dups)){
  colnames(pas_data)[which(colnames(pas_data)==dups[i])] <- paste0(colnames(pas_data)[which(colnames(pas_data)==dups[i])],1:2)
}

##-----------------------##
# Get the data objects ready for filtering, modelling

pollen_dat <- pas_data%>%filter(SAMPLE.TYPE=="Pollen")
filter_dat <- pas_data%>%filter(SAMPLE.TYPE!="Pollen")
xmat_pollen <- pollen_dat%>%dplyr::select(-SAMPLE.TYPE,-ID2,-CELL.FACTOR,-Sample)%>%data.matrix()
rownames(xmat_pollen) <- pollen_dat$ID2
xmat_filter <- filter_dat%>%dplyr::select(-SAMPLE.TYPE,-ID2,-CELL.FACTOR,-Sample)%>%data.matrix()

yvec_pollen <- as.character(pollen_dat$CELL.FACTOR)
yvec_filter <- gsub("Peak","",filter_dat$CELL.FACTOR)


##----------------------##
#Try different types of scalings/imputation

maxes <- matrix(apply(xmat_pollen,2,max,na.rm=T),nrow=nrow(xmat_pollen),ncol=ncol(xmat_pollen),byrow=TRUE)
mins <- matrix(apply(xmat_pollen,2,min,na.rm=T),nrow=nrow(xmat_pollen),ncol=ncol(xmat_pollen),byrow=TRUE)
xmat_pollen_scaled <- (xmat_pollen-mins)/(maxes-mins)
    
maxes <- matrix(apply(xmat_filter,2,max,na.rm=T),nrow=nrow(xmat_filter),ncol=ncol(xmat_filter),byrow=TRUE)
mins <- matrix(apply(xmat_filter,2,min,na.rm=T),nrow=nrow(xmat_filter),ncol=ncol(xmat_filter),byrow=TRUE)
xmat_filter_scaled <- (xmat_filter-mins)/(maxes-mins)



##----------------------##
# Modeling options (SPLS-DA since there are so many)

#SPLS DA
ncs <- 2
spls_model <- splsda(X=xmat_pollen_scaled, Y=yvec_pollen,ncomp=ncs, keepX=rep(30,ncs))
filter_hat <- predict(spls_model,newdata = xmat_filter_scaled, dist="max.dist")
pred_hat <- filter_hat$class$max.dist[,ncol(filter_hat$class$max.dist)]
con_mat <- get.confusion_matrix(truth=yvec_filter, predicted = pred_hat)
rowSums(abs(con_mat-diag(diag(con_mat)))) #Miss classified counts by species
con_mat

#Plot the proabilities for each class
cps <- data.frame(filter_hat$predict[,,dim(filter_hat$predict)[3]])
cps$Sample <- rownames(cps)
cps$Truth <- yvec_filter
class_probs <- cps%>%gather(Species,Value,-Sample,-Truth)%>%arrange(Truth,desc(Value))
class_probs$Sample <- factor(class_probs$Sample,levels=as.character(c(1:13)))
qplot(Sample,Value,colour=Species,data=class_probs)+geom_line(aes(group=Species))+theme_bw()+#theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab("Sample Number")+ylab("Probability of Each Species")+geom_vline(aes(xintercept=3.5),colour='gray50')+geom_vline(aes(xintercept=10.5),colour='gray50')




