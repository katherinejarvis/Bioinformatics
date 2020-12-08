# load the affy library
library(affy) #to read in CEL files
library(readxl) #to read biomart excel file
library(dplyr) # for functions
library(affyPLM) #for MA plot
library(limma)

#Reading in the spreadsheet from Biomart listing the gene names and their cooresponding probesets
mart_export <- read_excel('C:/Users/jarvi/Box/2020 spring/BMB502/Project/mart_export.xls')

# Read in the CEL files in the directory, then normalize the data
data <- ReadAffy()

#Convert to Expressionset
eset.data <- rma(data)

#Name of the samples
Samples<-c("Control_1","Control_2","Control_3","FrozenTissue_1", "FrozenTissue_2","FrozenTissue_3",
           "BL23_1","BL23_2","BL23_3","LP299v_1", "LP299v_2", "LP299v_3",
           "LP299v(A-)_1","LP299v(A-)_2","LP299v(A-)_3","Inflammed_1","Inflammed_2","Inflammed_3")

#Pairs figure to compare all samples
pairs(exprs(eset.data),labels=Samples,pch = '.', xlim=c(0,16), ylim=c(0,16))

#Name of the groups of samples
groups<-c("Control", "DirectFrozen","BL23","LP299v","LP299vA","Inflammed")

#Assigning the sample names to the expressionset
eset.data$sample<-Samples

#Removing duplicates of probe sets to have a array of the probesets of interest=geneinfo
genesinfo<-mart_export[!duplicated(cbind(mart_export$`Gene name`, mart_export$`AFFY HG U133 Plus 2 probe`)),]

#Making a expressionset of just the probeset of interest and adding Gene Names
probesetofinterest<-eset.data[featureNames(eset.data) %in% genesinfo$`AFFY HG U133 Plus 2 probe`,]

#Adding Gene Names to ExpressionSet
ID <- featureNames(probesetofinterest)
Symbol <- genesinfo$`Gene name`[genesinfo$`AFFY HG U133 Plus 2 probe`%in% ID ]
fData(probesetofinterest) <- data.frame(Symbol=Symbol)

#Convert Expressionset to Dataframe
write.exprs(probesetofinterest,file="probesofinterestdata.txt")
series_matrix_data <- read.table("probesofinterestdata.txt",sep="\t",nrows=9,quote="\"",header=TRUE)
series_dataframe<-data.frame(series_matrix_data[,2:19])
rownames(series_dataframe)<-series_matrix_data[,1]
colnames(series_dataframe)<-Samples

#Extract the Expression Values of the genes of interest
design <- model.matrix(~ 0+factor(c(rep(1:6,each=3))))
colnames(design) <- groups
fit <- lmFit(probesetofinterest, design)
fit <- eBayes(fit)
individualgenesstats<-topTableF(fit, adjust="BH")
write.csv(individualgenesstats,"genestats.csv", row.names =TRUE)

#Use the values found in 'fit' to find the fold change
contrast.matrix <- makeContrasts(Control- BL23, Control-LP299v, Control- LP299vA, Control-Inflammed, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
pvaluescontrast<-cbind(fit2$p.value,data.frame(fit2$F.p.value))
pvalues<-cbind(fit2$p.value,data.frame(fit$F.p.value))
contraststatsAll<-topTable(fit2)
write.csv(contraststatsAll,"foldchangestats.csv", row.names =TRUE)

#Assign palette for the heatmap
my_palette <- colorRampPalette(c("orange", "white", "blue"))(n = 32)
library(ComplexHeatmap)
Heatmap(as.matrix(individualgenesstats[2:7]), name ="Expression Values", row_labels =rownames(individualgenesstats),
        col=my_palette,row_split = individualgenesstats$Symbol)


# Barplot of Fold Change
library(ggplot2)
specie<- rep(genesinfo$`Gene name`,4)
condition <- rep(groups[3:6],each=9)
value <- c(contraststatsAll[,2],contraststatsAll[,3],contraststatsAll[,4],contraststatsAll[,5])
data <- data.frame(specie,condition,value)

ggplot(data, aes(x=value,y=condition, fill=condition)) + 
      geom_boxplot() + coord_flip() + 
      facet_wrap(~specie) + 
       geom_jitter(shape=16, position=position_jitter(0.1))+
       theme(legend.position="none") +
      ggtitle("Fold Change of Gene Expression Relative to Control")+
      ylab("Condition")+xlab("Expression Fold Change (log2)")

