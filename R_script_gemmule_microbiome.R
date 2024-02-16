
####  Project name : S.lacustris_gemmules_microbiome ####

####  Description of the project ####
#Sponges and their microbiome are playing important role in their ecosystems, such as the cycling of the organic matter. 
#However, the mechanisms associated with the transmission of the bacterial community from one generation to another is still unexplored, especially for freshwater sponges. 
#Spongilla lacustris is a model widely distributed in European rivers and lakes, producing before winter dormant cysts (named gemmules) for their asexual reproduction. 
#Through an in vitro experiment, this study aims to describe the dynamics of the bacterial community and its transmission modes following the hatching of these gemmules.

####  Objective of the script ####
# To analyse the 16S rRNA gene metabarcoding reads obtained from the DADA2 pipeline. 
# The main package used for this analysis will be phyloseq (for more details see tutorials such as https://joey711.github.io/phyloseq/index.html)

#### Author of the script : Dr. Benoît PAIX

#### Citation : not yet available (paper submitted)


####  1. Setting the working directory ####
setwd("C:/Users/benoit.paix/Mon Drive/2021_2023 Postdoc/8. Gemmules growth project/4. Lab notebook/3. Metabarcoding analyses/4. Phyloseq analysis v2 with 2")


####  2. Installation and loading of the packages ####


# example of few package installations : 
#install.packages("randomcoloR")
#install.packages("ggplot2")
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
#install.packages("stringi")
#install.packages("Rcpp")
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install(version='devel')
#BiocManager::install("microbiome")
#install.packages("hrbrthemes")
#install.packages("ggthemes")
#install.packages("RColorBrewer")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")



library(stringi)
library(vegan)
library(Rcpp)
library(ggplot2)
library(randomcoloR)
library(rlang)
library(phyloseq)
library(ape)
library(dplyr)
library(agricolae) #for statistical tests
library(RVAideMemoire) #for statistical tests
library(microbiome)
library(hrbrthemes)
library(cowplot) #for figures
#library(ggthemes)
library(RColorBrewer) #for figures
library(decontam) #for decontamination of the dataset
library(ggrepel) #for figures
library(microbiomeMarker) #for lefse analysis
library(multcompView)
library(rcompanion)
library(pairwiseAdonis) 
library(ggpubr)
library(ggh4x)


####  3. Creation of the output folders


#dir.create("1. Data_prep results")
#dir.create("2. Alpha_div results")
#dir.create("3. Beta_div results")
#dir.create("4. Compositional results")


####  4. Reading of the csv files (ASV, Taxonomy and Metadata tables)   ####

ASV_table = read.csv(file = "ASV_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(ASV_table)

TAX_table = read.csv(file = "Taxonomy_table.csv" , sep = ";" , header = T , row.names = 1)
dim(TAX_table)

META_table = read.csv(file = "Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(META_table)


TAX_table = as.data.frame(TAX_table)
TAX_table

TAX_table.v2 <- TAX_table

TAX_table.v2$Kingdom <- TAX_table.v2$Kingdom
TAX_table.v2$Phylum <- paste(TAX_table.v2$Kingdom,"|" ,TAX_table$Phylum)
TAX_table.v2$Class <- paste(TAX_table.v2$Phylum,"|" ,TAX_table$Class)
TAX_table.v2$Order <- paste(TAX_table.v2$Class ,"|" ,TAX_table$Order)
TAX_table.v2$Family <- paste(TAX_table.v2$Order,"|" ,TAX_table$Family)
TAX_table.v2$Genus <- paste(TAX_table.v2$Family,"|" ,TAX_table$Genus)

head(TAX_table.v2)

####  5. Creation of the main phyloseq object  ####

ASV_phylo = otu_table(ASV_table, taxa_are_rows = TRUE)
dim(ASV_phylo)


TAX_table = as.matrix(TAX_table.v2)

TAX_phylo = tax_table(TAX_table)
dim(TAX_phylo)

META_phylo = sample_data(META_table)
dim(META_phylo)

physeq_raw = phyloseq(ASV_phylo, TAX_phylo, META_phylo)
physeq_raw


####  6. Decontamination of the dataset #####


#inspection of the library sizes
df <- as.data.frame(sample_data(physeq_raw)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq_raw)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

plot_lib_size = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control))
plot_lib_size = plot_lib_size + geom_point()
plot_lib_size

ggsave(filename = "Plot_lib_size.pdf", 
       plot = plot_lib_size, 
       device = "pdf" , 
       width = 15 , height = 10, units = "cm", 
       path = "./1. Data_prep results")


#identification of the contaminants (prevalence method)

sample_data(physeq_raw)$is.neg <- sample_data(physeq_raw)$Sample_or_Control == "Control"

contamdf.prev05 <- isContaminant(physeq_raw, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)


contamdf.prev05 = cbind(as.data.frame(tax_table(physeq_raw)) , contamdf.prev05)
contamdf.prev05

write.csv(contamdf.prev05, file.path("./1. Data_prep results" , "Contamination_table_prev05.csv"))


# Make phyloseq object of presence-absence in negative controls and true samples

ps.pa <- transform_sample_counts(physeq_raw, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)

plot_prevalence = ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) 
plot_prevalence = plot_prevalence + geom_point() 
plot_prevalence = plot_prevalence + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
plot_prevalence


ggsave(filename = "Plot_prevalence.pdf", 
       plot = plot_prevalence, 
       device = "pdf" , 
       width = 15 , height = 10, units = "cm", 
       path = "./1. Data_prep results")




#create a phyloseq_decontam object without the contaminants
physeq_decontam = prune_taxa(!contamdf.prev05$contaminant, physeq_raw)
physeq_decontam


####  7. Removing unwanted taxa (EuK, Chloroplast, Mitochondria) ####



physeq_filtered = subset_taxa(physeq_decontam, 
                              (Kingdom != "Eukaryota" | is.na(Kingdom)) &
                              (Order!= "Bacteria | Cyanobacteria | Cyanobacteriia | Chloroplast"| is.na(Order)) &
                              (Family != "Bacteria | Proteobacteria | Alphaproteobacteria | Rickettsiales | Mitochondria"| is.na(Family)))

physeq_filtered


####  8.  Verification of the rarefaction curves ######

as.data.frame(t(otu_table(physeq_filtered)))

rarecurve(as.data.frame(t(otu_table(physeq_filtered))), step = 20, cex = 0.5)
#save manually the plot to ./1. Data_prep results

#set a minimum reads threshold based on the rarefaction curves results



####  9. Removing unwanted samples for the analysis (control samples and those below the threshold)  ######

sample_sums(physeq_filtered)

physeq_subsampled = prune_samples(sample_sums(physeq_filtered)>=10000, physeq_filtered)
physeq_subsampled


saveRDS(physeq_subsampled, "Physeq_subsampled.RDS")


####  10. Creating phylose_subsampled objects ######

physeq_subsampled = readRDS("Physeq_subsampled.RDS") #######
physeq_subsampled

sample_sums(physeq_subsampled)

physeq_subsampled_t0 = subset_samples(physeq_subsampled,  Sampling_time == "04. T0")
physeq_subsampled_t1 = subset_samples(physeq_subsampled,  Sampling_time == "05. T1" & Sample_type == "Cultivated sponge")
physeq_subsampled_t2 = subset_samples(physeq_subsampled,  Sampling_time == "06. T2" & Sample_type == "Cultivated sponge")
physeq_subsampled_t3 = subset_samples(physeq_subsampled,  Sampling_time == "07. T3" & Sample_type == "Cultivated sponge")
physeq_subsampled_t4 = subset_samples(physeq_subsampled,  Sampling_time == "08. T4" & Sample_type == "Cultivated sponge")
physeq_subsampled_t5 = subset_samples(physeq_subsampled,  Sampling_time == "09. T5" & Sample_type == "Cultivated sponge")
physeq_subsampled_t6 = subset_samples(physeq_subsampled,  Sampling_time == "10. T6" & Sample_type == "Cultivated sponge")

  


#### 11. ALPHA DIVERSITY ANALYSES ##################



# Create a phyloseq rarefied objects for alpha div analyses

physeq_rarefied = rarefy_even_depth(physeq_subsampled)
physeq_rarefied

sample_sums(physeq_rarefied)


# ALPHA DIVERSITY INDEX


data_alpha = estimate_richness(physeq_rarefied , measures = c("Observed", "Shannon", "Chao1"))
data_alpha


Pielou = data_alpha$Shannon / log(data_alpha$Observed)
Pielou

data_alpha = cbind(sample_data(physeq_rarefied), data_alpha , Pielou)
data_alpha = data_alpha[, c("Sample_type","Sampling_time" , "Cross_conditions",	"Washing"	, "Medium", "Shannon" ,"Chao1","Pielou")]
data_alpha

data_alpha$Sample_type[data_alpha$Sample_type == "Gemmule"] <- "2. Gemmules"
data_alpha$Sample_type[data_alpha$Sample_type == "Cultivated sponge"] <- "3. Cultivated sponges"
data_alpha$Sample_type[data_alpha$Sample_type == "In situ maternal sponge"] <- "1. In situ maternal sponge"
data_alpha$Sample_type[data_alpha$Sample_type == "Filtered freshwater"] <- "4. Filtered freshwater"


data_alpha

#Q1 : differences per samples types 

color_vector_1 = c("G+EM" = "lightcoral",
                   "G+EM_F-PB" = "plum3",
                   "G+EM_F+PB" = "sandybrown",
                   "G-EM" = "lightgreen",
                   "G-EM_F-PB" = "mediumaquamarine",
                   "G-EM_F+PB" = "yellow1")

plot_chao1_Q1 = ggplot(data_alpha, aes(Cross_conditions ,Chao1, fill = Cross_conditions))
plot_chao1_Q1 = plot_chao1_Q1 + geom_boxplot(alpha = 0.8, size = 1 ) + facet_grid(  ~ Sample_type, scales="free", space = "free") 
plot_chao1_Q1 = plot_chao1_Q1 + theme_bw(base_size = 15) 
plot_chao1_Q1 = plot_chao1_Q1 + theme(legend.position="none") +  theme(axis.text.x=element_blank())
plot_chao1_Q1 = plot_chao1_Q1 + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_chao1_Q1 = plot_chao1_Q1 + scale_fill_manual(values = c(color_vector_1), na.value = c("pink3","slateblue"))
plot_chao1_Q1 = plot_chao1_Q1 + geom_point(fill = "white" , size = 3, alpha = 0.5, pch = 21, stroke = 1)
plot_chao1_Q1

plot_shannon_Q1 = ggplot(data_alpha, aes(Cross_conditions ,Shannon, fill = Cross_conditions))
plot_shannon_Q1 = plot_shannon_Q1 + geom_boxplot(alpha = 0.8, size = 1 ) + facet_grid(  ~ Sample_type, scales="free", space = "free") 
plot_shannon_Q1 = plot_shannon_Q1 + theme_bw(base_size = 15) 
plot_shannon_Q1 = plot_shannon_Q1 + theme(legend.position="none") +  theme(axis.text.x=element_blank())
plot_shannon_Q1 = plot_shannon_Q1 + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_shannon_Q1 = plot_shannon_Q1 + scale_fill_manual(values = c(color_vector_1), na.value = c("pink3","slateblue"))
plot_shannon_Q1 = plot_shannon_Q1 + geom_point(fill = "white" , size = 3, alpha = 0.5, pch = 21, stroke = 1)
plot_shannon_Q1

plot_pielou_Q1 = ggplot(data_alpha, aes(Cross_conditions ,Pielou, fill = Cross_conditions))
plot_pielou_Q1 = plot_pielou_Q1 + geom_boxplot(alpha = 0.8, size = 1 ) + facet_grid(  ~ Sample_type, scales="free", space = "free") 
plot_pielou_Q1 = plot_pielou_Q1 + theme_bw(base_size = 15) 
plot_pielou_Q1 = plot_pielou_Q1 + theme(legend.position="none") +  theme(axis.text.x=element_text(angle = 45, hjust = 1))
plot_pielou_Q1 = plot_pielou_Q1 + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_pielou_Q1 = plot_pielou_Q1 + scale_fill_manual(values = c(color_vector_1), na.value = c("pink3","slateblue"))
plot_pielou_Q1 = plot_pielou_Q1 + geom_point(fill = "white" , size = 3, alpha = 0.5, pch = 21, stroke = 1)
plot_pielou_Q1


plot_alpha_Q1 =  plot_grid( plot_chao1_Q1, plot_shannon_Q1 ,plot_pielou_Q1 , labels=c("A", "B", "C"),rel_heights = c(1, 1, 1.6) , ncol = 1, nrow = 3 ,label_size = 20)
plot_alpha_Q1


ggsave(filename = "Plot_alpha_Q1.pdf", 
       plot = plot_alpha_Q1, 
       device = "pdf" , 
       width = 30 , height = 25, units = "cm", 
       path = "./2. Alpha_div results")


#Q2 : within gemmules & cultivated gemmules : per time ###########


data_alpha_Q2 = subset(data_alpha, Sample_type == "3. Cultivated sponges" | Sample_type == "2. Gemmules"  )
data_alpha_Q2 


plot_chao1_Q2 = ggplot(data_alpha_Q2, aes(Cross_conditions ,Chao1, fill = Cross_conditions))
plot_chao1_Q2 = plot_chao1_Q2 + geom_boxplot(alpha = 0.8, size = 1 ) + facet_grid(  ~ Sampling_time, scales="free", space = "free") 
plot_chao1_Q2 = plot_chao1_Q2 + theme_bw(base_size = 15) 
plot_chao1_Q2 = plot_chao1_Q2 + theme(legend.position="none") +  theme(axis.text.x=element_blank())
plot_chao1_Q2 = plot_chao1_Q2 + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_chao1_Q2 = plot_chao1_Q2 + scale_fill_manual(values = c(color_vector_1), na.value = c("pink3","slateblue"))
plot_chao1_Q2 = plot_chao1_Q2 + geom_point(fill = "white" , size = 3, alpha = 0.5, pch = 21, stroke = 1)
plot_chao1_Q2

plot_shannon_Q2 = ggplot(data_alpha_Q2, aes(Cross_conditions ,Shannon, fill = Cross_conditions))
plot_shannon_Q2 = plot_shannon_Q2 + geom_boxplot(alpha = 0.8, size = 1 ) + facet_grid(  ~ Sampling_time, scales="free", space = "free") 
plot_shannon_Q2 = plot_shannon_Q2 + theme_bw(base_size = 15) 
plot_shannon_Q2 = plot_shannon_Q2 + theme(legend.position="none") +  theme(axis.text.x=element_text(angle = 90, hjust = 1))
plot_shannon_Q2 = plot_shannon_Q2 + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_shannon_Q2 = plot_shannon_Q2 + scale_fill_manual(values = c(color_vector_1), na.value = c("pink3","slateblue"))
plot_shannon_Q2 = plot_shannon_Q2 + geom_point(fill = "white" , size = 3, alpha = 0.5, pch = 21, stroke = 1)
plot_shannon_Q2

plot_pielou_Q2 = ggplot(data_alpha_Q2, aes(Cross_conditions ,Pielou, fill = Cross_conditions))
plot_pielou_Q2 = plot_pielou_Q2 + geom_boxplot(alpha = 0.8, size = 1 ) + facet_grid(  ~ Sampling_time, scales="free", space = "free") 
plot_pielou_Q2 = plot_pielou_Q2 + theme_bw(base_size = 15) 
plot_pielou_Q2 = plot_pielou_Q2 + theme(legend.position="none") +  theme(axis.text.x=element_text(angle = 90, hjust = 1))
plot_pielou_Q2 = plot_pielou_Q2 + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_pielou_Q2 = plot_pielou_Q2 + scale_fill_manual(values = c(color_vector_1), na.value = c("pink3","slateblue"))
plot_pielou_Q2 = plot_pielou_Q2 + geom_point(fill = "white" , size = 3, alpha = 0.5, pch = 21, stroke = 1)
plot_pielou_Q2


plot_shannon_Q2

ggsave(filename = "Plot_shannon_Q2.pdf", 
       plot = plot_shannon_Q2, 
       device = "pdf" , 
       width = 25 , height = 15, units = "cm", 
       path = "./2. Alpha_div results")


plot_chao1_pielou_Q2 =  plot_grid( plot_chao1_Q2 ,plot_pielou_Q2 , labels=c("A", "B"),rel_heights = c(1,  1.6) , ncol = 1, nrow = 2 ,label_size = 20)
plot_chao1_pielou_Q2


ggsave(filename = "Plot_chao1_pielou_Q2.pdf", 
       plot = plot_chao1_pielou_Q2, 
       device = "pdf" , 
       width = 30 , height = 20, units = "cm", 
       path = "./2. Alpha_div results")




### Univariate statistics for alpha div(anova) ####


#Q1 : differences between sample types 



#Shapiro tests from the checking of the normal distribution

shapiro_data_shannon_Q1 = shapiro.test(data_alpha$Shannon)
shapiro_data_shannon_Q1$p.value
#From the output, the p-value < 0.05 implying that the distribution of the data is significantly different from normal distribution. In other words, we cannot assume the normality.

shapiro_data_chao1_Q1 = shapiro.test(data_alpha$Chao1)
shapiro_data_chao1_Q1
#From the output, the p-value < 0.05 implying that the distribution of the data is significantly different from normal distribution. In other words, we cannot assume the normality.

shapiro_data_pielou_Q1 = shapiro.test(data_alpha$Pielou)
shapiro_data_pielou_Q1
#From the output, the p-value < 0.05 implying that the distribution of the data is significantly different from normal distribution. In other words, we cannot assume the normality.



#writing a table with the results from the shapiro tests with the 3 indexes
shapiro_data_alpha_Q1 <- matrix(nrow = 3 ,  ncol=2, byrow=TRUE)
colnames(shapiro_data_alpha_Q1) = c("W","p-value")
rownames(shapiro_data_alpha_Q1) = c("Shannon","Chao1","Pielou")

shapiro_data_alpha_Q1[,1] <- c(shapiro_data_shannon_Q1$statistic,
                              shapiro_data_chao1_Q1$statistic,
                              shapiro_data_pielou_Q1$statistic)

shapiro_data_alpha_Q1[,2]<- c(shapiro_data_shannon_Q1$p.value,
                              shapiro_data_chao1_Q1$p.value,
                              shapiro_data_pielou_Q1$p.value)

shapiro_data_alpha_Q1


write.csv(shapiro_data_alpha_Q1, file.path("./2. Alpha_div results" , "Shapiro_data_alpha_Q1.csv"))


# Kruskal tests (non parametric anova) with Sample type factor : Q1

krustal_data_shannon_Q1 = kruskal.test(Shannon ~ Sample_type, data_alpha)
krustal_data_shannon_Q1

krustal_data_chao1_Q1 = kruskal.test(Chao1 ~ Sample_type, data_alpha)
krustal_data_chao1_Q1

krustal_data_pielou_Q1 = kruskal.test(Pielou ~ Sample_type, data_alpha)
krustal_data_pielou_Q1


kruskal_data_alpha_Q1  <- matrix(nrow = 3 ,  ncol=3, byrow=TRUE)
colnames(kruskal_data_alpha_Q1) = c("Chi-square","Df","p-value")
rownames(kruskal_data_alpha_Q1) = c("Shannon","Chao1","Pielou")

kruskal_data_alpha_Q1

kruskal_data_alpha_Q1[,1] <- c(krustal_data_shannon_Q1$statistic,
                               krustal_data_chao1_Q1$statistic,
                               krustal_data_pielou_Q1$statistic)

kruskal_data_alpha_Q1[,2]<- c(krustal_data_shannon_Q1$parameter,
                              krustal_data_chao1_Q1$parameter,
                              krustal_data_pielou_Q1$parameter)

kruskal_data_alpha_Q1[,3]<- c(krustal_data_shannon_Q1$p.value,
                              krustal_data_chao1_Q1$p.value,
                              krustal_data_pielou_Q1$p.value)

kruskal_data_alpha_Q1


write.csv(kruskal_data_alpha_Q1, file.path("./2. Alpha_div results" , "Kruskal_data_alpha_Q1.csv"))


# Wilcox pairwise comparison test for anova with p < 0.05 : here for Pielou

data_wilcox_pielou_Q1 = pairwise.wilcox.test(data_alpha$Pielou, data_alpha$Sample_type, p.adjust.method = "bonf")
data_wilcox_pielou_Q1
data_wilcox_pielou_Q1$p.value

data_wilcox_pielou_Q1_full = fullPTable(data_wilcox_pielou_Q1$p.value)
data_wilcox_pielou_Q1_full

indices_wilcox_pielou_Q1 = multcompLetters(data_wilcox_pielou_Q1_full,
                            compare="<",
                            threshold=0.05,
                            Letters=letters,
                            reversed = FALSE)
indices_wilcox_pielou_Q1$Letters



data_wilcox_pielou_Q1_full2 = cbind(data_wilcox_pielou_Q1_full, indices_wilcox_pielou_Q1$Letters )
data_wilcox_pielou_Q1_full2

write.csv(data_wilcox_pielou_Q1_full2, file.path("./2. Alpha_div results" , "Wilcox_data_pielou_Q1.csv"))



#Q2 : differences between conditions over time


#Shapiro tests 

shapiro_data_shannon_Q2 = shapiro.test(data_alpha_Q2$Shannon)
shapiro_data_shannon_Q2
#From the output, the p-value < 0.05 implying that the distribution of the data is significantly different from normal distribution. In other words, we cannot assume the normality.

shapiro_data_chao1_Q2 = shapiro.test(data_alpha_Q2$Chao1)
shapiro_data_chao1_Q2
#From the output, the p-value < 0.05 implying that the distribution of the data is significantly different from normal distribution. In other words, we cannot assume the normality.

shapiro_data_pielou_Q2 = shapiro.test(data_alpha_Q2$Pielou)
shapiro_data_pielou_Q2

#writing a table with the results from the shapiro tests with the 3 indexes
shapiro_data_alpha_Q2 <- matrix(nrow = 3 ,  ncol=2, byrow=TRUE)
colnames(shapiro_data_alpha_Q2) = c("W","p-value")
rownames(shapiro_data_alpha_Q2) = c("Shannon","Chao1","Pielou")

shapiro_data_alpha_Q2[,1] <- c(shapiro_data_shannon_Q2$statistic,
                               shapiro_data_chao1_Q2$statistic,
                               shapiro_data_pielou_Q2$statistic)

shapiro_data_alpha_Q2[,2]<- c(shapiro_data_shannon_Q2$p.value,
                              shapiro_data_chao1_Q2$p.value,
                              shapiro_data_pielou_Q2$p.value)

shapiro_data_alpha_Q2


write.csv(shapiro_data_alpha_Q2, file.path("./2. Alpha_div results" , "Shapiro_data_alpha_Q2.csv"))

#From the output, the p-value < 0.05 implying that the distribution of the data is significantly different from normal distribution. In other words, we cannot assume the normality.
#for Q2 I need to create a new factor column in the dataset : Sampling_time_cross_condition. This factor will be used for the Kruskal-wallis tests
# --> solution used : or (to simplify the analysis, but harder to justify in m&m ) a factor Sampling_time_washing (since the medium type seems to have no effect)

St_Cc = paste(data_alpha_Q2$Sampling_time , 
              data_alpha_Q2$Washing)

St_Cc

data_alpha_Q2.2 = cbind(data_alpha_Q2, St_Cc)
head(data_alpha_Q2.2)

# Kruskal tests (non parametric anova) with Sample_time and Cross_conditions factors : Q2


krustal_data_chao1_Q2 = kruskal.test(Chao1 ~ St_Cc  , data_alpha_Q2.2)
krustal_data_chao1_Q2

krustal_data_shannon_Q2 = kruskal.test(Shannon ~ St_Cc  , data_alpha_Q2.2)
krustal_data_shannon_Q2

krustal_data_pielou_Q2 = kruskal.test(Pielou ~ St_Cc  , data_alpha_Q2.2)
krustal_data_pielou_Q2


kruskal_data_alpha_Q2  <- matrix(nrow = 3 ,  ncol=3, byrow=TRUE)
colnames(kruskal_data_alpha_Q2) = c("Chi-square","Df","p-value")
rownames(kruskal_data_alpha_Q2) = c("Shannon","Chao1","Pielou")

kruskal_data_alpha_Q2

kruskal_data_alpha_Q2[,1] <- c(krustal_data_shannon_Q2$statistic,
                               krustal_data_chao1_Q2$statistic,
                               krustal_data_pielou_Q2$statistic)

kruskal_data_alpha_Q2[,2]<- c(krustal_data_shannon_Q2$parameter,
                              krustal_data_chao1_Q2$parameter,
                              krustal_data_pielou_Q2$parameter)

kruskal_data_alpha_Q2[,3]<- c(krustal_data_shannon_Q2$p.value,
                              krustal_data_chao1_Q2$p.value,
                              krustal_data_pielou_Q2$p.value)

kruskal_data_alpha_Q2


write.csv(kruskal_data_alpha_Q2, file.path("./2. Alpha_div results" , "Kruskal_data_alpha_Q2.csv"))


data_wilcox_shannon_Q2 = pairwise.wilcox.test(data_alpha_Q2.2$Shannon, data_alpha_Q2.2$St_Cc, p.adjust.method = "bonferroni")
data_wilcox_shannon_Q2
data_wilcox_shannon_Q2$p.value

data_wilcox_chao1_Q2 = pairwise.wilcox.test(data_alpha_Q2.2$Chao1, data_alpha_Q2.2$St_Cc, p.adjust.method = "bonferroni")
data_wilcox_chao1_Q2
data_wilcox_chao1_Q2$p.value
  
data_wilcox_pielou_Q2 = pairwise.wilcox.test(data_alpha_Q2.2$Pielou, data_alpha_Q2.2$St_Cc, p.adjust.method = "bonferroni")
data_wilcox_pielou_Q2
data_wilcox_pielou_Q2$p.value


data_wilcox_shannon_Q2_full = fullPTable(data_wilcox_shannon_Q2$p.value)
data_wilcox_shannon_Q2_full

indices_wilcox_shannon_Q2 = multcompLetters(data_wilcox_shannon_Q2_full,
                                           compare="<",
                                           threshold=0.05,
                                           Letters=letters,
                                           reversed = FALSE)
indices_wilcox_shannon_Q2$Letters



data_wilcox_shannon_Q2_full2 = cbind(data_wilcox_shannon_Q2_full, indices_wilcox_shannon_Q2$Letters )
data_wilcox_shannon_Q2_full2


write.csv(data_wilcox_shannon_Q2_full2, file.path("./2. Alpha_div results" , "Wilcox_data_shannon_Q2.csv"))



data_wilcox_chao1_Q2_full = fullPTable(data_wilcox_chao1_Q2$p.value)
data_wilcox_chao1_Q2_full

indices_wilcox_chao1_Q2 = multcompLetters(data_wilcox_chao1_Q2_full,
                                            compare="<",
                                            threshold=0.05,
                                            Letters=letters,
                                            reversed = FALSE)
indices_wilcox_chao1_Q2$Letters



data_wilcox_chao1_Q2_full2 = cbind(data_wilcox_chao1_Q2_full, indices_wilcox_chao1_Q2$Letters )
data_wilcox_chao1_Q2_full2


write.csv(data_wilcox_chao1_Q2_full2, file.path("./2. Alpha_div results" , "Wilcox_data_chao1_Q2.csv"))


data_wilcox_pielou_Q2_full = fullPTable(data_wilcox_pielou_Q2$p.value)
data_wilcox_pielou_Q2_full

indices_wilcox_pielou_Q2 = multcompLetters(data_wilcox_pielou_Q2_full,
                                          compare="<",
                                          threshold=0.05,
                                          Letters=letters,
                                          reversed = FALSE)
indices_wilcox_pielou_Q2$Letters



data_wilcox_pielou_Q2_full2 = cbind(data_wilcox_pielou_Q2_full, indices_wilcox_pielou_Q2$Letters )
data_wilcox_pielou_Q2_full2


write.csv(data_wilcox_pielou_Q2_full2, file.path("./2. Alpha_div results" , "Wilcox_data_pielou_Q2.csv"))



#### 12. BETA DIVERSITY : NMDS ############

physeq_compo = transform(physeq_subsampled, "compositional")
physeq_compo



nmds_Q1 = ordinate(physeq_compo, "NMDS", "bray")
nmds_Q1$points
nmds_Q1$stress
#write.csv(nmds$points , file = "nmds.points.t3.csv")




### Sample plot with phyloseq with all samples

palette_nmdsplot = c("pink",
                     "pink3",
                     "tan4",
                     "orchid", 
                     "cornflowerblue" , 
                     "aquamarine3" , 
                     "green" , 
                     "yellow" , 
                     "darkorange" , 
                     "brown1")

palette_nmdsplot

data_nmds = plot_ordination(physeq_compo, nmds_Q1, type="samples", justDF = TRUE)
data_nmds

plot_nmds_Q1_samples = ggplot(data_nmds, aes(NMDS1, NMDS2))
plot_nmds_Q1_samples = plot_nmds_Q1_samples + geom_point(aes(fill = Sampling_time, pch = Sample_type) , stroke = 1.5 , size = 7, alpha = 0.7, color = "black" , show.legend = T)
plot_nmds_Q1_samples = plot_nmds_Q1_samples + scale_shape_manual(values = c(21, 25 ,22, 23))
plot_nmds_Q1_samples = plot_nmds_Q1_samples + theme_bw(base_size = 20) 
plot_nmds_Q1_samples = plot_nmds_Q1_samples + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_nmds_Q1_samples = plot_nmds_Q1_samples + scale_fill_manual(values = palette_nmdsplot)
plot_nmds_Q1_samples = plot_nmds_Q1_samples + annotate("text", label="2D stress = 0.15", x=-3.4, y=1.6, size = 5)
plot_nmds_Q1_time_samples = plot_nmds_Q1_samples
plot_nmds_Q1_time_samples

ggsave(filename = "Plot_NMDS_Q1_time_samples.pdf", 
       plot = plot_nmds_Q1_time_samples, 
       device = "pdf" , 
       width = 35 , height = 20, units = "cm", 
       path = "./3. Beta_div_results")


#same plot but to visualize the different cross-conditions
color_vector_1

plot_nmds_Q1_samples = ggplot(data_nmds, aes(NMDS1, NMDS2))
plot_nmds_Q1_samples = plot_nmds_Q1_samples + geom_point(aes(fill = Cross_conditions, pch = Sample_type) , stroke = 1.5 , size = 7, alpha = 0.7, color = "black" , show.legend = T)
plot_nmds_Q1_samples = plot_nmds_Q1_samples + scale_shape_manual(values = c(21, 25 ,22, 23))
plot_nmds_Q1_samples = plot_nmds_Q1_samples + theme_bw(base_size = 20) 
plot_nmds_Q1_samples = plot_nmds_Q1_samples + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_nmds_Q1_samples = plot_nmds_Q1_samples + scale_fill_manual(values = color_vector_1)
plot_nmds_Q1_samples = plot_nmds_Q1_samples + annotate("text", label="2D stress = 0.16", x=-3.4, y=1.6, size = 5)
plot_nmds_Q1_crossco_samples = plot_nmds_Q1_samples
plot_nmds_Q1_crossco_samples

ggsave(filename = "Plot_NMDS_Q1_crossco_samples.pdf", 
       plot = plot_nmds_Q1_crossco_samples, 
       device = "pdf" , 
       width = 35 , height = 20, units = "cm", 
       path = "./3. Beta_div_results")

# ASVs plots all samples (Q1) : is it necessary? 


sum_taxa = taxa_sums(physeq_compo)
sum_taxa
select.taxa = plot_ordination(physeq_compo, nmds, type="taxa", justDF = TRUE)
select.taxa
select.taxa.2 = cbind(select.taxa, sum_taxa)
select.taxa.2

select.taxa.4 = subset(select.taxa.2, sum_taxa > 0.1 )
select.taxa.4

check_asv_names <- select.taxa.4[order(select.taxa.4$sum_taxa,decreasing=TRUE),]
head(check_asv_names)

head(select.taxa.4)
select.taxa.4$Phylum
levels(factor(select.taxa.4$Phylum))



colorpal_phylum = c( "Bacteria | Acidobacteriota"   ="darkgoldenrod1",
                    "Bacteria | Actinobacteriota"  ="darkslateblue",
                    "Bacteria | Armatimonadota" = "orchid4",
                    "Bacteria | Bacteroidota" = "yellow1",
                    "Bacteria | Bdellovibrionota"     =  "olivedrab1",
                    "Bacteria | Cyanobacteria"     = "plum1",
                    "Bacteria | Gemmatimonadota" ="orangered1",
                    "Bacteria | Myxococcota"       ="slategray2",
                    "Bacteria | Patescibacteria"     = "darkblue",
                    "Bacteria | Planctomycetota"    ="aquamarine4",
                    "Bacteria | Proteobacteria"    ="dodgerblue",
                    "Bacteria | Spirochaetota"     ="peachpuff",
                    "Bacteria | Verrucomicrobiota" = "brown")


plot_nmds_Q1_ASV = ggplot(select.taxa.4, aes(NMDS1, NMDS2))
plot_nmds_Q1_ASV = plot_nmds_Q1_ASV +  geom_point(aes(fill = Phylum, size = sum_taxa), stroke = 1 , color = "black", alpha = 0.6, pch = 21) #+ facet_wrap(vars(Family))
plot_nmds_Q1_ASV = plot_nmds_Q1_ASV + scale_size_continuous(range = c(1, 15)) 
plot_nmds_Q1_ASV = plot_nmds_Q1_ASV  + scale_fill_manual(values = colorpal_phylum)
plot_nmds_Q1_ASV = plot_nmds_Q1_ASV + theme_bw(base_size = 17)  + theme(legend.position="right") 
plot_nmds_Q1_ASV = plot_nmds_Q1_ASV + guides(fill = guide_legend(override.aes = list(size = 5)))
plot_nmds_Q1_ASV = plot_nmds_Q1_ASV + geom_point(data = data_nmds, aes(data_nmds$NMDS1, data_nmds$NMDS2), pch = 22 , size = 7, alpha = 0)
plot_nmds_Q1_ASV


ggsave(filename = "Plot_NMDS_Q1_ASV.pdf", 
       plot = plot_nmds_Q1_ASV, 
       device = "pdf" , 
       width = 35 , height = 20, units = "cm", 
       path = "./3. Beta_div_results")


##### TEST STAT PERMANOVA et pairwise for Q1.1
metadata <- as(sample_data(physeq_compo), "data.frame")
metadata

data_distbeta = as.matrix(distance(physeq_compo, method="bray"))
data_distbeta

metadata$Sample_type

data_permnova_Q1.1 = adonis2(data_distbeta ~ Sample_type, data = metadata)
data_permnova_Q1.1
write.csv(as.data.frame(data_permnova_Q1.1), 
          file.path("./3. Beta_div_results" , "Permanova_Q1.1_sampletype.csv"))
 
data_pairwiseadonis_Q1.1 = pairwise.adonis(data_distbeta, metadata$Sample_type)
data_pairwiseadonis_Q1.1

write.csv(data_pairwiseadonis_Q1.1, 
          file.path("./3. Beta_div_results" , "Pairwise_adonis_Q1.1_sampletype.csv"))


# dataset for Q1.2 : cross condition * time comparison, so only cultivated juveniles sponges

metadata
physeq_subsampled
physeq_subsampled_cultivated = subset_samples(physeq_subsampled,  Sample_type  == "Cultivated sponge")
physeq_subsampled_cultivated

physeq_subsampled_cultivated_compo = transform(physeq_subsampled_cultivated, "compositional")
physeq_subsampled_cultivated_compo

metadata <- as(sample_data(physeq_subsampled_cultivated_compo), "data.frame")
metadata
head(metadata)

data_distbeta = as.matrix(distance(physeq_subsampled_cultivated_compo, method="bray"))
data_distbeta


data_permnova_Q1.2 = adonis2(data_distbeta ~ Sampling_time*Washing*Medium , data = metadata)
data_permnova_Q1.2
write.csv(as.data.frame(data_permnova_Q1.2), 
          file.path("./3. Beta_div_results" , "Permanova_Q1.2_sampletime_crossco.csv"))

data_pairwiseadonis_Q1.2 = pairwise.adonis(data_distbeta, metadata$Sampling_time)
data_pairwiseadonis_Q1.2

write.csv(data_pairwiseadonis_Q1.2, 
          file.path("./3. Beta_div_results" , "Pairwise_adonis_Q1.2_sampletime.csv"))

data_pairwiseadonis_Q1.2 = pairwise.adonis(data_distbeta, metadata$Cross_conditions)
data_pairwiseadonis_Q1.2

write.csv(data_pairwiseadonis_Q1.2, 
          file.path("./3. Beta_div_results" , "Pairwise_adonis_Q1.2_crossco.csv"))



### NMDS for Q2 : Sample plot with phyloseq by time

color_vector_1 = c("G+EM" = "lightcoral",
                   "G+EM_F-PB" = "plum3",
                   "G+EM_F+PB" = "sandybrown",
                   "G-EM" = "lightgreen",
                   "G-EM_F-PB" = "mediumaquamarine",
                   "G-EM_F+PB" = "yellow1")
# phyloseq compo object for all sampling times 

physeq_compo_t0 = transform(physeq_subsampled_t0, "compositional")
physeq_compo_t1 = transform(physeq_subsampled_t1, "compositional")
physeq_compo_t2 = transform(physeq_subsampled_t2, "compositional")
physeq_compo_t3 = transform(physeq_subsampled_t3, "compositional")
physeq_compo_t4 = transform(physeq_subsampled_t4, "compositional")
physeq_compo_t5 = transform(physeq_subsampled_t5, "compositional")
physeq_compo_t6 = transform(physeq_subsampled_t6, "compositional")

#just for a test
physeq_compo_maternal_and_t0 = transform(physeq_subsampled_maternal_and_t0, "compositional")
physeq_compo_maternal_and_t0

nmds_t1 = ordinate(physeq_compo_t1, "NMDS", "bray")
nmds_t1$points
nmds_t1$stress

data_nmds_t1 = plot_ordination(physeq_compo_t1, nmds_t1, type="samples", justDF = TRUE)
data_nmds_t1


plot.nmds = ggplot(data_nmds_t1, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 21, stroke = 1.5 , alpha = 0.8, size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
#plot.nmds = plot.nmds + geom_text_repel(aes(label = Samples_name))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t1 = plot.nmds



nmds_t2 = ordinate(physeq_compo_t2, "NMDS", "bray")
nmds_t2$points
nmds_t2$stress

data_nmds_t2 = plot_ordination(physeq_compo_t2, nmds_t2, type="samples", justDF = TRUE)
data_nmds_t2


plot.nmds = ggplot(data_nmds_t2, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 21, stroke = 1.5 , alpha = 0.8, size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
#plot.nmds = plot.nmds + geom_text_repel(aes(label = Samples_name))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t2 = plot.nmds



nmds_t3 = ordinate(physeq_compo_t3, "NMDS", "bray")
nmds_t3$points
nmds_t3$stress

data_nmds_t3 = plot_ordination(physeq_compo_t3, nmds_t3, type="samples", justDF = TRUE)
data_nmds_t3


plot.nmds = ggplot(data_nmds_t3, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 21, stroke = 1.5 , alpha = 0.8, size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
#plot.nmds = plot.nmds + geom_text_repel(aes(label = Samples_name))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t3 = plot.nmds




nmds_t4 = ordinate(physeq_compo_t4, "NMDS", "bray")
nmds_t4$points
nmds_t4$stress

data_nmds_t4 = plot_ordination(physeq_compo_t4, nmds_t4, type="samples", justDF = TRUE)
data_nmds_t4


plot.nmds = ggplot(data_nmds_t4, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 21, stroke = 1.5 , alpha = 0.8, size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
#plot.nmds = plot.nmds + geom_text_repel(aes(label = Samples_name))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t4 = plot.nmds




nmds_t5 = ordinate(physeq_compo_t5, "NMDS", "bray")
nmds_t5$points
nmds_t5$stress

data_nmds_t5 = plot_ordination(physeq_compo_t5, nmds_t5, type="samples", justDF = TRUE)
data_nmds_t5


plot.nmds = ggplot(data_nmds_t5, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 21, stroke = 1.5 , alpha = 0.8, size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
#plot.nmds = plot.nmds + geom_text_repel(aes(label = Samples_name))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t5 = plot.nmds




nmds_t6 = ordinate(physeq_compo_t6, "NMDS", "bray")
nmds_t6$points
nmds_t6$stress

data_nmds_t6 = plot_ordination(physeq_compo_t6, nmds_t6, type="samples", justDF = TRUE)
data_nmds_t6


plot.nmds = ggplot(data_nmds_t6, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 21, stroke = 1.5 , alpha = 0.8, size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
#plot.nmds = plot.nmds + geom_text_repel(aes(label = Samples_name))
legend_t1.6 = as_ggplot(get_legend(plot.nmds+ theme_bw(base_size = 11)))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t6 = plot.nmds
legend_t1.6



nmds_t0 = ordinate(physeq_compo_t0, "NMDS", "bray")
nmds_t0$points
nmds_t0$stress

data_nmds_t0 = plot_ordination(physeq_compo_t0, nmds_t0, type="samples", justDF = TRUE)
data_nmds_t0

plot.nmds = ggplot(data_nmds_t0, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 22, stroke = 1.5 , alpha = 0.8,size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
#plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
#plot.nmds = plot.nmds + geom_text_repel(aes(label = Samples_name))
legend_t0 = as_ggplot(get_legend(plot.nmds+ theme_bw(base_size = 11)))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t0 = plot.nmds

legend_t0


#just for a test
physeq_compo_maternal_and_t0
nmds_t0maternal_test = ordinate(physeq_compo_maternal_and_t0, "NMDS", "bray")
nmds_t0maternal_test$points
nmds_t0maternal_test$stress

data_nmds_t0maternal_test = plot_ordination(physeq_compo_maternal_and_t0, nmds_t0maternal_test, type="samples", justDF = TRUE)
data_nmds_t0maternal_test

plot.nmds = ggplot(data_nmds_t0maternal_test, aes(NMDS1, NMDS2))
plot.nmds = plot.nmds + geom_point(aes(fill = Cross_conditions) , pch = 22, stroke = 1.5 , alpha = 0.8,size = 5, color = "black" , show.legend = T)
plot.nmds = plot.nmds + theme_bw(base_size = 12)
plot.nmds = plot.nmds + scale_fill_manual(values = color_vector_1)
#plot.nmds = plot.nmds + guides(fill = guide_legend(override.aes = list(shape = 21)))
#plot.nmds = plot.nmds + scale_fill_manual(values = palette.boxplot.2)
#plot.nmds = plot.nmds + annotate("text", label="2D stress = 0.129 ", x=-0.4, y=0.4, size = 6)
plot.nmds = plot.nmds + geom_text_repel(aes(label = Sampling_time))
legend_t0 = as_ggplot(get_legend(plot.nmds+ theme_bw(base_size = 11)))
plot.nmds = plot.nmds + theme(legend.position="none",
                              axis.title=element_blank())
plot.nmds
plot.nmds.t0maternal_test = plot.nmds
### end of test : distance t-1 to +EM lower than t-1 to -EM  ! ok



plot_all_time = plot_grid(plot.nmds.t0, 
                          legend_t0, 
                          legend_t1.6,
                          plot.nmds.t1, 
                          plot.nmds.t2, 
                          plot.nmds.t3, 
                          plot.nmds.t4, 
                          plot.nmds.t5, 
                          plot.nmds.t6,
                          labels = c('T0', '', '', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6'), label_size = 12, ncol = 3, nrow = 3)
plot_all_time
plot_NMDS_Q2_all_bytime = plot_all_time
plot_NMDS_Q2_all_bytime

ggsave(filename = "Plot_NMDS_Q2_all_bytime.pdf", 
       plot = plot_NMDS_Q2_all_bytime, 
       device = "pdf" , 
       width = 25 , height = 20, units = "cm", 
       path = "./3. Beta_div_results")


nmds_t0$stress
nmds_t1$stress
nmds_t2$stress
nmds_t3$stress
nmds_t4$stress
nmds_t5$stress
nmds_t6$stress


### permanova analyses for each times


#t1
metadata_t1 <- as(sample_data(physeq_compo_t1), "data.frame")
metadata_t1

data_distbeta_t1 = as.matrix(distance(physeq_compo_t1, method="bray"))
data_distbeta_t1

data_permnova_Q2.t1 = adonis2(data_distbeta_t1 ~ Washing*Medium , data = metadata_t1)
data_permnova_Q2.t1

data_pairwiseadonis_Q2_t1 = pairwise.adonis(data_distbeta_t1, metadata_t1$Cross_conditions)
data_pairwiseadonis_Q2_t1

#t2
metadata_t2 <- as(sample_data(physeq_compo_t2), "data.frame")
metadata_t2

data_distbeta_t2 = as.matrix(distance(physeq_compo_t2, method="bray"))
data_distbeta_t2

data_permnova_Q2.t2 = adonis2(data_distbeta_t2 ~ Washing*Medium , data = metadata_t2)
data_permnova_Q2.t2

data_pairwiseadonis_Q2_t2 = pairwise.adonis(data_distbeta_t2, metadata_t2$Cross_conditions)
data_pairwiseadonis_Q2_t2

#t3
metadata_t3 <- as(sample_data(physeq_compo_t3), "data.frame")
metadata_t3

data_distbeta_t3 = as.matrix(distance(physeq_compo_t3, method="bray"))
data_distbeta_t3

data_permnova_Q2.t3 = adonis2(data_distbeta_t3 ~ Washing*Medium , data = metadata_t3)
data_permnova_Q2.t3

data_pairwiseadonis_Q2_t3 = pairwise.adonis(data_distbeta_t3, metadata_t3$Cross_conditions)
data_pairwiseadonis_Q2_t3

#t4
metadata_t4 <- as(sample_data(physeq_compo_t4), "data.frame")
metadata_t4

data_distbeta_t4 = as.matrix(distance(physeq_compo_t4, method="bray"))
data_distbeta_t4

data_permnova_Q2.t4 = adonis2(data_distbeta_t4 ~ Washing*Medium , data = metadata_t4)
data_permnova_Q2.t4

data_pairwiseadonis_Q2_t4 = pairwise.adonis(data_distbeta_t4, metadata_t4$Cross_conditions)
data_pairwiseadonis_Q2_t4

#t5
metadata_t5 <- as(sample_data(physeq_compo_t5), "data.frame")
metadata_t5

data_distbeta_t5 = as.matrix(distance(physeq_compo_t5, method="bray"))
data_distbeta_t5

data_permnova_Q2.t5 = adonis2(data_distbeta_t5 ~ Washing*Medium , data = metadata_t5)
data_permnova_Q2.t5

data_pairwiseadonis_Q2_t5 = pairwise.adonis(data_distbeta_t5, metadata_t5$Cross_conditions)
data_pairwiseadonis_Q2_t5

#t6
metadata_t6 <- as(sample_data(physeq_compo_t6), "data.frame")
metadata_t6

data_distbeta_t6 = as.matrix(distance(physeq_compo_t6, method="bray"))
data_distbeta_t6

data_permnova_Q2.t6 = adonis2(data_distbeta_t6 ~ Washing*Medium , data = metadata_t6)
data_permnova_Q2.t6

data_pairwiseadonis_Q2_t6 = pairwise.adonis(data_distbeta_t6, metadata_t6$Cross_conditions)
data_pairwiseadonis_Q2_t6

#t0
metadata_t0 <- as(sample_data(physeq_compo_t0), "data.frame")
metadata_t0

data_distbeta_t0 = as.matrix(distance(physeq_compo_t0, method="bray"))
data_distbeta_t0

data_permnova_Q2.t0 = adonis2(data_distbeta_t0 ~ Washing , data = metadata_t0)
data_permnova_Q2.t0

data_pairwiseadonis_Q2_t0 = pairwise.adonis(data_distbeta_t0, metadata_t0$Washing)
data_pairwiseadonis_Q2_t0

permanova_Q2_alltime =    rbind(c("result t0", "", "", "", ""),
                          as.data.frame(data_permnova_Q2.t0),
                          c("result t1", "", "", "", ""),
                          as.data.frame(data_permnova_Q2.t1),
                          c("result t2", "", "", "", ""),
                          as.data.frame(data_permnova_Q2.t2),
                          c("result t3", "", "", "", ""),
                          as.data.frame(data_permnova_Q2.t3),
                          c("result t4", "", "", "", ""),
                          as.data.frame(data_permnova_Q2.t4),
                          c("result t5", "", "", "", ""),
                          as.data.frame(data_permnova_Q2.t5),
                          c("result t6", "", "", "", ""),
                          as.data.frame(data_permnova_Q2.t6))
permanova_Q2_alltime


write.csv(permanova_Q2_alltime, 
          file.path("./3. Beta_div_results" , "Permanova_Q2_alltime.csv"))

data_pairwiseadonis_Q2_t0
data_pairwiseadonis_Q2_t1
data_pairwiseadonis_Q2_t2
data_pairwiseadonis_Q2_t3
data_pairwiseadonis_Q2_t4
data_pairwiseadonis_Q2_t5
data_pairwiseadonis_Q2_t6


pairwise_adonis_Q2_alltime = rbind(c("result t0", "", "", "", "", "", "", ""),
                                   data_pairwiseadonis_Q2_t0,
                                 c("result t1", "", "", "", "", "", "", ""),
                                 data_pairwiseadonis_Q2_t1,
                                 c("result t2", "", "", "", "", "", "", ""),
                                 data_pairwiseadonis_Q2_t2,
                                 c("result t3", "", "", "", "", "", "", ""),
                                 data_pairwiseadonis_Q2_t3,
                                 c("result t4", "", "", "", "", "", "", ""),
                                 data_pairwiseadonis_Q2_t4,
                                 c("result t5", "", "", "", "", "", "", ""),
                                 data_pairwiseadonis_Q2_t5,
                                 c("result t6", "", "", "", "", "", "", ""),
                                 data_pairwiseadonis_Q2_t6)
pairwise_adonis_Q2_alltime

write.csv(pairwise_adonis_Q2_alltime, 
          file.path("./3. Beta_div_results" , "Pairwise_adonis_Q2_alltime.csv"))


#### 13. BETA DIVERSITY : BETA-DISPER analyses per times (Q2) ########



#t1
metadata_t1 <- as(sample_data(physeq_compo_t1), "data.frame")
metadata_t1

data_distbeta_t1 = distance(physeq_compo_t1, method="bray")
data_distbeta_t1

permdisp_t1 <- betadisper(data_distbeta_t1, data_nmds_t1$Cross_conditions)
dist_centroid = permdisp_t1$distances
dist_centroid

data_nmds_t1_dist_centroid = cbind(data_nmds_t1, dist_centroid)
data_nmds_t1_dist_centroid

#t2
metadata_t2 <- as(sample_data(physeq_compo_t2), "data.frame")
metadata_t2

data_distbeta_t2 = distance(physeq_compo_t2, method="bray")
data_distbeta_t2

permdisp_t2 <- betadisper(data_distbeta_t2, data_nmds_t2$Cross_conditions)
dist_centroid = permdisp_t2$distances
dist_centroid

data_nmds_t2_dist_centroid = cbind(data_nmds_t2, dist_centroid)
data_nmds_t2_dist_centroid

#t3
metadata_t3 <- as(sample_data(physeq_compo_t3), "data.frame")
metadata_t3

data_distbeta_t3 = distance(physeq_compo_t3, method="bray")
data_distbeta_t3

permdisp_t3 <- betadisper(data_distbeta_t3, data_nmds_t3$Cross_conditions)
dist_centroid = permdisp_t3$distances
dist_centroid

data_nmds_t3_dist_centroid = cbind(data_nmds_t3, dist_centroid)
data_nmds_t3_dist_centroid

#t4
metadata_t4 <- as(sample_data(physeq_compo_t4), "data.frame")
metadata_t4

data_distbeta_t4 = distance(physeq_compo_t4, method="bray")
data_distbeta_t4

permdisp_t4 <- betadisper(data_distbeta_t4, data_nmds_t4$Cross_conditions)
dist_centroid = permdisp_t4$distances
dist_centroid

data_nmds_t4_dist_centroid = cbind(data_nmds_t4, dist_centroid)
data_nmds_t4_dist_centroid

#t5
metadata_t5 <- as(sample_data(physeq_compo_t5), "data.frame")
metadata_t5

data_distbeta_t5 = distance(physeq_compo_t5, method="bray")
data_distbeta_t5

permdisp_t5 <- betadisper(data_distbeta_t5, data_nmds_t5$Cross_conditions)
dist_centroid = permdisp_t5$distances
dist_centroid

data_nmds_t5_dist_centroid = cbind(data_nmds_t5, dist_centroid)
data_nmds_t5_dist_centroid

#t6
metadata_t6 <- as(sample_data(physeq_compo_t6), "data.frame")
metadata_t6

data_distbeta_t6 = distance(physeq_compo_t6, method="bray")
data_distbeta_t6

permdisp_t6 <- betadisper(data_distbeta_t6, data_nmds_t6$Cross_conditions)
dist_centroid = permdisp_t6$distances
dist_centroid

data_nmds_t6_dist_centroid = cbind(data_nmds_t6, dist_centroid)
data_nmds_t6_dist_centroid

#t0
metadata_t0 <- as(sample_data(physeq_compo_t0), "data.frame")
metadata_t0

data_distbeta_t0 = distance(physeq_compo_t0, method="bray")
data_distbeta_t0

permdisp_t0 <- betadisper(data_distbeta_t0, data_nmds_t0$Cross_conditions)
dist_centroid = permdisp_t0$distances
dist_centroid

data_nmds_t0_dist_centroid = cbind(data_nmds_t0, dist_centroid)
data_nmds_t0_dist_centroid


### all data_nmds_dist_centroid in one table

data_nmds_dist_centroid_alltime = rbind(data_nmds_t0_dist_centroid,
                                        data_nmds_t1_dist_centroid,
                                        data_nmds_t2_dist_centroid,
                                        data_nmds_t3_dist_centroid,
                                        data_nmds_t4_dist_centroid,
                                        data_nmds_t5_dist_centroid,
                                        data_nmds_t6_dist_centroid)

data_nmds_dist_centroid_alltime


plot_distcentr = ggplot(data_nmds_dist_centroid_alltime, aes(Cross_conditions ,dist_centroid, fill = Cross_conditions))
plot_distcentr = plot_distcentr + geom_boxplot(alpha = 0.8, size = 1)   + facet_grid(  ~ Sampling_time, scales="free", space = "free")
plot_distcentr = plot_distcentr + geom_point(size = 2, alpha = 0.8, pch = 21, stroke = 1)
plot_distcentr = plot_distcentr + theme_bw(base_size = 15) 
plot_distcentr = plot_distcentr + theme(legend.position="left") +  theme(axis.text.x=element_blank())
plot_distcentr = plot_distcentr + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_distcentr = plot_distcentr + scale_fill_manual(values = color_vector_1)
plot_distcentr

ggsave(filename = "Plot_distcentr_Q2_all_bytime.pdf", 
       plot = plot_distcentr, 
       device = "pdf" , 
       width = 35 , height = 15, units = "cm", 
       path = "./3. Beta_div_results")



## Shapiro test for the betadisper analyse (for the dist to centroids) at each time


shapiro_data_distcentr_Q2_t1 = shapiro.test(data_nmds_t1_dist_centroid$dist_centroid)
shapiro_data_distcentr_Q2_t1
#From the output, the p-value > 0.05 implying that the distribution of the data is not significantly different from normal distribution. 
#In other words, we can assume the normality. so ANOVA test + Tukey's HSD test

shapiro_data_distcentr_Q2_t2 = shapiro.test(data_nmds_t2_dist_centroid$dist_centroid)
shapiro_data_distcentr_Q2_t2
#From the output, the p-value > 0.05 implying that the distribution of the data is not significantly different from normal distribution. 
#In other words, we can assume the normality. so ANOVA test + Tukey's HSD test

shapiro_data_distcentr_Q2_t3 = shapiro.test(data_nmds_t3_dist_centroid$dist_centroid)
shapiro_data_distcentr_Q2_t3
#From the output, the p-value < 0.05 implying that the distribution of the data is  significantly different from normal distribution. 
#In other words, we cannot assume the normality. so Kruskal-Wallis test + Wilcoxon test

shapiro_data_distcentr_Q2_t4 = shapiro.test(data_nmds_t4_dist_centroid$dist_centroid)
shapiro_data_distcentr_Q2_t4
#From the output, the p-value < 0.05 implying that the distribution of the data is  significantly different from normal distribution. 
#In other words, we cannot assume the normality. so Kruskal-Wallis test + Wilcoxon test

shapiro_data_distcentr_Q2_t5 = shapiro.test(data_nmds_t5_dist_centroid$dist_centroid)
shapiro_data_distcentr_Q2_t5
#From the output, the p-value > 0.05 implying that the distribution of the data is not significantly different from normal distribution. 
#In other words, we can assume the normality. so ANOVA test + Tukey's HSD test

shapiro_data_distcentr_Q2_t6 = shapiro.test(data_nmds_t6_dist_centroid$dist_centroid)
shapiro_data_distcentr_Q2_t6
#From the output, the p-value > 0.05 implying that the distribution of the data is not significantly different from normal distribution. 
#In other words, we can assume the normality. so ANOVA test + Tukey's HSD test

shapiro_data_distcentr_Q2_t0 = shapiro.test(data_nmds_t0_dist_centroid$dist_centroid)
shapiro_data_distcentr_Q2_t0
#From the output, the p-value > 0.05 implying that the distribution of the data is not significantly different from normal distribution. 
#In other words, we can assume the normality. so ANOVA test + Tukey's HSD test


#writing a table with the results from the shapiro tests with the 3 indexes
shapiro_data_distcentr_Q2_alltime <- matrix(nrow = 7 ,  ncol=2, byrow=TRUE)
colnames(shapiro_data_distcentr_Q2_alltime) = c("W","p-value")
rownames(shapiro_data_distcentr_Q2_alltime) = c("t0","t1","t2", "t3", "t4", "t5", "t6")

shapiro_data_distcentr_Q2_alltime[,1] <- c(shapiro_data_distcentr_Q2_t0$statistic,
                                           shapiro_data_distcentr_Q2_t1$statistic,
                                           shapiro_data_distcentr_Q2_t2$statistic,
                                           shapiro_data_distcentr_Q2_t3$statistic,
                                           shapiro_data_distcentr_Q2_t4$statistic,
                                           shapiro_data_distcentr_Q2_t5$statistic,
                                           shapiro_data_distcentr_Q2_t6$statistic)

shapiro_data_distcentr_Q2_alltime[,2]<- c(shapiro_data_distcentr_Q2_t0$p.value,
                                          shapiro_data_distcentr_Q2_t1$p.value,
                                          shapiro_data_distcentr_Q2_t2$p.value,
                                          shapiro_data_distcentr_Q2_t3$p.value,
                                          shapiro_data_distcentr_Q2_t4$p.value,
                                          shapiro_data_distcentr_Q2_t5$p.value,
                                          shapiro_data_distcentr_Q2_t6$p.value)

shapiro_data_distcentr_Q2_alltime


write.csv(shapiro_data_distcentr_Q2_alltime, 
          file.path("./3. Beta_div_results" , "Shapiro_data_distcentr_Q2_alltime.csv"))


#### ANOVA tests for dist to centroids from betadisper, for each time (Q2)


#t1
anovadata_t1_dist_centroid <- aov(dist_centroid ~  Washing*Medium, data_nmds_t1_dist_centroid)
anovadata_t1_dist_centroid = as.data.frame(summary(anovadata_t1_dist_centroid)[[1]])
anovadata_t1_dist_centroid

data_hsd_t1_dist_centroid = HSD.test(aov(dist_centroid ~ Cross_conditions, data_nmds_t1_dist_centroid), "Cross_conditions", group=T)
data_hsd_t1_dist_centroid

#t2
anovadata_t2_dist_centroid <- aov(dist_centroid ~  Washing*Medium, data_nmds_t2_dist_centroid)
anovadata_t2_dist_centroid = as.data.frame(summary(anovadata_t2_dist_centroid)[[1]])
anovadata_t2_dist_centroid

data_hsd_t2_dist_centroid = HSD.test(aov(dist_centroid ~ Cross_conditions, data_nmds_t2_dist_centroid), "Cross_conditions", group=T)
data_hsd_t2_dist_centroid

#t3 ###### no anova but kruskal
krustaldata_t3_dist_centroid = kruskal.test(dist_centroid ~ Cross_conditions, data_nmds_t3_dist_centroid)
krustaldata_t3_dist_centroid


data_wilcox_t3_dist_centroid = pairwise.wilcox.test(data_nmds_t3_dist_centroid$dist_centroid, data_nmds_t3_dist_centroid$Cross_conditions, p.adjust.method = "bonferroni")
data_wilcox_t3_dist_centroid
data_wilcox_t3_dist_centroid$p.value


#t4  ###### no anova but kruskal
krustaldata_t4_dist_centroid = kruskal.test(dist_centroid ~ Cross_conditions, data_nmds_t4_dist_centroid)
krustaldata_t4_dist_centroid

data_wilcox_t3_dist_centroid = pairwise.wilcox.test(data_nmds_t4_dist_centroid$dist_centroid, data_nmds_t4_dist_centroid$Cross_conditions, p.adjust.method = "bonferroni")
data_wilcox_t3_dist_centroid
data_wilcox_t3_dist_centroid$p.value


#t5
anovadata_t5_dist_centroid <- aov(dist_centroid ~  Washing*Medium, data_nmds_t5_dist_centroid)
anovadata_t5_dist_centroid = as.data.frame(summary(anovadata_t5_dist_centroid)[[1]])
anovadata_t5_dist_centroid

data_hsd_t5_dist_centroid = HSD.test(aov(dist_centroid ~ Cross_conditions, data_nmds_t5_dist_centroid), "Cross_conditions", group=T)
data_hsd_t5_dist_centroid

#t6
anovadata_t6_dist_centroid <- aov(dist_centroid ~  Washing*Medium, data_nmds_t6_dist_centroid)
anovadata_t6_dist_centroid = as.data.frame(summary(anovadata_t6_dist_centroid)[[1]])
anovadata_t6_dist_centroid

data_hsd_t6_dist_centroid = HSD.test(aov(dist_centroid ~ Cross_conditions, data_nmds_t6_dist_centroid), "Cross_conditions", group=T)
data_hsd_t6_dist_centroid


#t0
anovadata_t0_dist_centroid <- aov(dist_centroid ~  Washing, data_nmds_t0_dist_centroid)
anovadata_t0_dist_centroid = as.data.frame(summary(anovadata_t0_dist_centroid)[[1]])
anovadata_t0_dist_centroid

data_hsd_t0_dist_centroid = HSD.test(aov(dist_centroid ~ Cross_conditions, data_nmds_t0_dist_centroid), "Cross_conditions", group=T)
data_hsd_t0_dist_centroid


#### table with all the anova results from the dist centroid analyse (betadisper for Q2) per times
anovadata_alltime_dist_centroid = rbind(c("result t0", "", "", "", ""),
                                          anovadata_t0_dist_centroid,
                                        c("result t1", "", "", "", ""),
                                          anovadata_t1_dist_centroid,
                                        c("result t2", "", "", "", ""),
                                          anovadata_t2_dist_centroid,
                                        c("result t5", "", "", "", ""),
                                          anovadata_t5_dist_centroid,
                                        c("result t6", "", "", "", ""),
                                          anovadata_t6_dist_centroid)
anovadata_alltime_dist_centroid

write.csv(anovadata_alltime_dist_centroid, 
          file.path("./3. Beta_div_results" , "Anova_data_distcentr_Q2_alltime.csv"))


HSD_data_alltime_dist_centroid = rbind(c("result t0", ""),
                                       data_hsd_t0_dist_centroid$groups,
                                       c("result t1", ""),
                                       data_hsd_t1_dist_centroid$groups,
                                       c("result t2", ""),
                                       data_hsd_t2_dist_centroid$groups,
                                       c("result t5", ""),
                                       data_hsd_t5_dist_centroid$groups,
                                       c("result t6", ""),
                                       data_hsd_t6_dist_centroid$groups)
HSD_data_alltime_dist_centroid

write.csv(HSD_data_alltime_dist_centroid, 
          file.path("./3. Beta_div_results" , "HSD_data_distcentr_Q2_alltime.csv"))


#### 13. COMPOSITIONAL ANALYSES  BAR PLOT analysis  at the family level ######

physeq_compo_family = aggregate_rare(physeq_compo, level = "Family", detection = 10/100, prevalence = 0/100)
physeq_compo_family

taxtable_fam = as.data.frame(tax_table(physeq_compo_family))
print(as.data.frame(taxtable_fam$Family), row.names = FALSE)

ntaxa(physeq_compo_family)

  
colorpal_family = c(  "Bacteria | Actinobacteriota | Actinobacteria | Frankiales | Sporichthyaceae"=   "violet",
                      "Bacteria | Actinobacteriota | Actinobacteria | Micrococcales | Microbacteriaceae"=  "plum4",
                      "Bacteria | Actinobacteriota | Actinobacteria | Micrococcales | Micrococcaceae" =  "violetred1" , 
                      "Bacteria | Bacteroidota | Bacteroidia | Chitinophagales | Chitinophagaceae" =   "sandybrown",
                      "Bacteria | Bacteroidota | Bacteroidia | Chitinophagales | Saprospiraceae"=  "orange" , 
                      "Bacteria | Bacteroidota | Bacteroidia | Cytophagales | Hymenobacteraceae"=  "sienna2" , 
                      "Bacteria | Bacteroidota | Bacteroidia | Cytophagales | Microscillaceae"=   "peru", 
                      "Bacteria | Bacteroidota | Bacteroidia | Cytophagales | Spirosomaceae"=   "gold1", 
                      "Bacteria | Bacteroidota | Bacteroidia | Flavobacteriales | Flavobacteriaceae"=   "tan1", 
                      "Bacteria | Cyanobacteria | Cyanobacteriia | Synechococcales | Cyanobiaceae"=   "purple", 
                      "Bacteria | Planctomycetota | Planctomycetes | Pirellulales | Pirellulaceae"=   "red", 
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Caulobacterales | Caulobacteraceae"= "palegreen4"  , 
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Elsterales | Elsteraceae"=   "greenyellow",
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Rhizobiales | Devosiaceae"=  "palegreen" , 
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Rhizobiales | Rhizobiaceae"= "darkgreen"  , 
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Rhodobacterales | Rhodobacteraceae"=  "olivedrab" , 
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Rhodospirillales | Terasakiellaceae"=  "springgreen3" , 
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Sphingomonadales | Sphingomonadaceae"=  "yellowgreen" , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Burkholderiaceae"=  "turquoise1" , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Comamonadaceae"= "lightblue2"  , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Methylophilaceae"=   "turquoise3",
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Oxalobacteraceae"= "cornflowerblue"  , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Rhodocyclaceae"= "darkblue"  , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Enterobacterales | Alteromonadaceae"=  "deepskyblue" , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | NA | NA" = "slategrey" ,
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | Moraxellaceae"=  "turquoise4" , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | Pseudomonadaceae"= "royalblue4"  , 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Salinisphaerales | Solimonadaceae"=   "blue", 
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Xanthomonadales | Xanthomonadaceae"= "darkslateblue"  , 
                      "Bacteria | Verrucomicrobiota | Chlamydiae | Chlamydiales | Parachlamydiaceae"= "yellow"  , 
                      "Other"= "grey"   )


# for a barplot without border the phyloseq function plot_bar need to be changed. 
#Use the plot_bar_2 function as follows

plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


physeq_compo_family_cultivated = subset_samples(physeq_compo_family,  Sample_type  == "Cultivated sponge")
physeq_compo_family_cultivated


physeq_compo_family_nocultivated = subset_samples(physeq_compo_family,  Sample_type  != "Cultivated sponge")
physeq_compo_family_nocultivated

sample_data(physeq_compo_family_nocultivated)$Sample_type[sample_data(physeq_compo_family_nocultivated)$Sample_type== "Gemmule"] <- "3. Gemmules" 
sample_data(physeq_compo_family_nocultivated)$Sample_type[sample_data(physeq_compo_family_nocultivated)$Sample_type== "In situ maternal sponge"] <- "2. In situ maternal sponge" 
sample_data(physeq_compo_family_nocultivated)$Sample_type[sample_data(physeq_compo_family_nocultivated)$Sample_type== "Filtered freshwater"] <- "1. Filtered freshwater" 


barplot = plot_bar_2(physeq_compo_family_cultivated ,"Samples_name" ,fill = "Family") + theme_bw()
barplot = barplot + geom_bar(stat = "identity"  , position="stack", width=1)
barplot = barplot + scale_fill_manual(values = colorpal_family)
barplot = barplot + theme(legend.text = element_text(size=11),
                          axis.ticks.y = element_blank(),
                          legend.position="bottom",
                          axis.text.x = element_blank(),
                          legend.title = element_blank(),
                          axis.ticks.x=element_blank())
barplot = barplot + guides(fill = guide_legend(ncol = 2)) 
barplot = barplot +  scale_y_percent() #+ coord_flip()
barplot = barplot + theme(axis.title = element_blank())
barplot = barplot + facet_nested(~ Sampling_time+Washing+Medium , scales = "free", space = "free")
barplot = barplot + theme(panel.spacing = unit(0, "cm", data = NULL), 
                          panel.border = element_rect(color = "black", fill = NA, size = 1))
barplot_family_cultivated = barplot
barplot_family_cultivated

ggsave(filename = "Barplot_family_cultivated.pdf", 
       plot = barplot_family_cultivated, 
       device = "pdf" , 
       width = 40 , height = 30, units = "cm", 
       path = "./4. Compositional results")



barplot = plot_bar_2(physeq_compo_family_nocultivated ,"Samples_name" ,fill = "Family") + theme_bw()
barplot = barplot + geom_bar(stat = "identity"  , position="stack", width=1)
barplot = barplot + scale_fill_manual(values = colorpal_family)
barplot = barplot + theme(legend.text = element_text(size=10),
                          axis.ticks.y = element_blank(),
                          legend.position="right",
                          axis.text.x = element_blank(),
                          legend.title = element_blank(),
                          axis.ticks.x=element_blank())
barplot = barplot + guides(fill = guide_legend(ncol = 1)) 
barplot = barplot +  scale_y_percent() #+ coord_flip()
barplot = barplot + theme(axis.title = element_blank())
barplot = barplot + facet_nested(~ Sample_type+Sampling_time+Washing, scales = "free", space = "free")
barplot = barplot + theme(panel.spacing = unit(0, "cm", data = NULL), 
                          panel.border = element_rect(color = "black", fill = NA, size = 1))
barplot_family_nocultivated = barplot
barplot_family_nocultivated



ggsave(filename = "Barplot_family_nocultivated.pdf", 
       plot = barplot_family_nocultivated, 
       device = "pdf" , 
       width = 40 , height = 22, units = "cm", 
       path = "./4. Compositional results")


