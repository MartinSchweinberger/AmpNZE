##################################################################
# Titel:      The Amplifier System of New Zealand English - Part 3
# R version:  3.4.1 (2017-06-30) -- "Single Candle"
# Autor:      Martin Schweinberger
# Date:       2018-11-06
# Contact:    martin.schweinberger.hh@gmail.com
##################################################################
# Disclaimer: If you have questions,suggestions or you found errors
#             or in case you would to provide feedback, questions
#             write an email to martin.schweinberger.hh@gmail.com.
# Citation:   If you use this script or results thereof, please cite it as:
#             Schweinberger, Martin. 2018. "The Amplifier System of New Zealand English, Part 3",
#             unpublished R script, The University of Queensland.
###############################################################
#                   START
###############################################################
# remove all lists from the current workspace
rm(list=ls(all=T))
# set wd
setwd("D:\\Uni\\Projekte\\02-Intensification\\AmpNZE")
# load libraries
library(ape)
library(ca)
library(car)
library(cfa)
library(cluster)
library(dplyr)
library(ggnetwork)
library(ggplot2)
library(gsubfn)
library(Hmisc)
library(MASS)
library(languageR)
library(lme4)
library(plyr)
library(psych)
library(pvclust)
library(QuantPsyc)
library(readr)
library(reshape)
#install.packAges("C:/R/Rling_1.0.tar.gz", repos = NULL, type = "source")
library(Rling)
library(rms)
library(stringr)
library(tm)
library(tnet)
library(network) # keep after tnet
library(zoo)
# load self written Function
source("D:\\R/multiplot_ggplot2.R") # for multiple ggplot2 plots in one window
###############################################################
# Setting options
options(stringsAsFactors = F)
options(scipen = 999)
options(max.prAmplified=10000)
# define imAge dnzectors
imAgeDirectory<-"images"
###############################################################
# load data
ampnze <- read.table("ampnze04_clean.txt", sep = "\t", header = T)
###############################################################
# inspect data
str(ampnze)

###############################################################
# recode text type
ampnze$Genre <- factor(ampnze$Genre, levels = c("PrivateDialogue", "PublicDialogue",
                                                "UnscriptedMonologue", "ScriptedMonologue"))
###############################################################
#              SEMANTIC VECTOR SPACE MODEL 1
# tabulate data
t1 <- tapply(ampnze$Amplified, list(ampnze$Adjective, ampnze$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3nze <- t3
t3nze <- t3nze[, 2: ncol(t3nze)]
# remove Adjectives that were not Amplifiedensified
t3nze <- t3nze[rowSums(t3nze) > 0, ]
# save row and column names
colnamesnze <- colnames(t3nze)
rownamesnze <- rownames(t3nze)
# turn dataframe Amplifiedo matrix
svsmnze <- as.matrix(t3nze)
# convert token frequency to type frequency
#svsmnze <- apply(svsmnze, 1, Function(x) { x <- ifelse(x > 1, 1, x) } )
svsmnze <- t(svsmnze)
#svsmnze <- svsmnze[, colSums(svsmnze) >= 2]
#svsmnze <- svsmnze[rowSums(svsmnze) >= 2, ]
svsmnze

# compute expected values
svsmnze.exp <- chisq.test(svsmnze)$expected
# calculate PMI and PPMI
svsmnze.PMI <- log2(svsmnze/svsmnze.exp)
svsmnze.PPMI <- ifelse(svsmnze.PMI < 0, 0, svsmnze.PMI)
# calculate cosine similarity
svsmnze.tmp1 <- svsmnze.PPMI
svsmnze.cos <- cossim(svsmnze.tmp1)
#round(svsmnze.cos, 2)
###############################################################
#               CLUSTER SEMANTIC VECTORS
# find max value that is not 1
svsmnze.cos.test <- apply(svsmnze.cos, 1, function(x){
  x <- ifelse(x == 1, 0, x) } )
maxval <- max(svsmnze.cos.test)
# create distance matrix
svsmnze.dist <- 1 - (svsmnze.cos/maxval)
clustd <- as.dist(svsmnze.dist)
# create distance matrix
clustd <- dist(svsmnze.cos, method = "manhattan") 
# alternative methods
# eucledian - not good when dealing with many dimensions
# manhattan - most popular choice
# method - here the difference between poAmplifieds dominates
# canberra - for count data
# binary - for binary data only!
# minkowski - is not a true distance measure

# find optimal number of clusters
asw <- as.vector(unlist(sapply(2:nrow(svsmnze)-1, function(x) pam(clustd, k = x)$silinfo$avg.width)))
# determine the optimal number of clusters (max width is optimal)
optclust <- which(asw == max(asw))+1 # optimal number of clusters

# inspect clustering with optimal number of clusters
svsmnze.clust <- pam(clustd, optclust)
svsmnze.clust$clustering

# create cluster object
# alternative methods: "single", "ward.D2", "averAge", "mcquitty", "median", "centroid"
ampnzehclust <- hclust(clustd, method="ward.D")    
# plot cluster solution
png("images/Clustnze.png",  width = 480, height = 480) # save plot
plot(ampnzehclust, main = "", xlab = "", ylab = "")
rect.hclust(ampnzehclust, k = optclust)
dev.off()
# load libraries for nicer dendrograms
library(factoextra)
library(dendextend)
# plot with colored clusters
png("images/Clustnze.png",  width = 800, height = 500) # save plot
fviz_dend(ampnzehclust, k = optclust, cex = 1, horiz = F,  
          k_colors = c("grey30"), 
          rect_border = c("grey30"), 
          rect_fill = F, main = "", labels_track_height=2, rect = T)
dev.off()
# plot as unrooted tree
png("images/PhyClustAmpnze.png",  width = 680, height = 480) 
fviz_dend(ampnzehclust, k = optclust, color_labels_by_k = T, type = "phylogenic", repel = TRUE, cex = .9,
          k_colors = c("grey70", "grey50", "grey30"))
dev.off()
###############################################################
# Unrooted clustering
# library ape
library(ape)
# convert 'hclust' to 'phylo' object
phylo_tree = as.phylo(ampnzehclust)
# get edges
graph_edges = phylo_tree$edge
# library igraph
library(igraph)
# get graph from edge list
graph_net = graph.edgelist(graph_edges)
# extract layout (x-y coords)
graph_layout = layout.auto(graph_net)
# number of observations
nobs = nrow(svsmnze.cos)
# save plot
png("images/UClustAmpnze.png",  width = 680, height = 480) 
# start plot
plot(graph_layout[,1], graph_layout[,2], type = "n", axes = FALSE,
     xlab = "", ylab = "")
# draw tree branches
segments(
  x0 = graph_layout[graph_edges[,1],1], 
  y0 = graph_layout[graph_edges[,1],2],
  x1 = graph_layout[graph_edges[,2],1],
  y1 = graph_layout[graph_edges[,2],2],
  col = "gray90", lwd = 2
)
# add labels
text(graph_layout[1:nobs,1], graph_layout[1:nobs,2],
     phylo_tree$tip.label, cex = .9, xpd = TRUE, font = 1)
dev.off()
###############################################################
#                 WARNING
#             DATA REDUCTION
# exclude amplifiers that are very dissimilar to main group of amplifiers
rmvamp <- c("especially", "exceptionally", "excruciatingly", "profoundly", "mighty", "decidedly", "much")
nrow(ampnze)

ampnze <- ampnze[!ampnze$Variant %in% rmvamp, ]
nrow(ampnze)

###############################################################
# check if variable levels need to be collapsed (only PrivateDialogue data)
tstdt <- ampnze[ampnze$Genre == "PrivateDialogue",]
FileSpeakerAgetfb <- ftable(tstdt$FileSpeaker, tstdt$Age)
FileSpeakerAgetfb <- apply(FileSpeakerAgetfb, 1, function(x) ifelse(x > 1, 1, x)) 
rowSums(FileSpeakerAgetfb)

###############################################################
# prepare data for plotting
# create data frame with relevant variables
pd <- data.frame(ampnze$Age, ampnze$Genre, ampnze$Function, ampnze$Amplified, ampnze$Variant)
# clean col names
colnames(pd) <- gsub("ampnze.", "", colnames(pd))
colnames(pd)[5] <- "Variant"
# convert Age column
Agelbs <- names(table(pd$Age))
pd$Age <- ifelse(pd$Age == "16-24", 4,
                 ifelse(pd$Age == "25-39", 3, 
                        ifelse(pd$Age == "40-49", 2, 
                               ifelse(pd$Age == "50+", 1, pd$Age))))
# multiply Amplified * 100 to get percent for Variant
pd$Amplified <- ifelse(pd$Amplified == 1, 100, 0)
# convert Age and Amplified Amplifiedo a numeric variables
clnm <- c("Age", "Amplified")
pd[clnm] <- lapply(pd[clnm], as.numeric)
famps <- names(table(pd$Variant))[which(table(pd$Variant) > 20)]
# reclassify Adjectives - infreq. Adjectives are collapsed Amplifiedo category other
pd$Variant <- ifelse(pd$Variant  %in% famps, pd$Variant , "other")
# create variables 
pd$other <- ifelse(pd$Variant == "other", 100, 0)
pd$pretty <- ifelse(pd$Variant == "pretty", 100, 0) 
pd$really <- ifelse(pd$Variant == "really", 100, 0) 
pd$so <- ifelse(pd$Variant == "so", 100, 0) 
pd$very <- ifelse(pd$Variant == "very", 100, 0)
pd$zero <- ifelse(pd$Variant == "0", 100, 0)
###############################################################
# p1
p1d <- pd
# increase n
p1d <- rbind(p1d, p1d, p1d)
# start plot: Amplified
p1 <- ggplot(p1d, aes(x = jitter(Age), y = Amplified)) +
  geom_smooth(aes(y = Amplified), size=.5, col = "gray30", lty = "longdash") +
  facet_grid(vars(Function), vars(Genre)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  theme_light(base_size = 10) +
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(Agelbs))
ggsave(file = paste(imAgeDirectory,"Amplified_Genre_Function.png",sep="/"))
p1

###############################################################
# p2
p2d <- pd
# remove non-amplified instances
p2d <- p2d[p2d$Amplified != 0,]
# increase n
p2d <- rbind(p2d, p2d, p2d, p2d, p2d, p2d, p2d, p2d, p2d, p2d)
# start plot: all
p2 <- ggplot(p2d, aes(x = jitter(Age), y = very)) +
  facet_grid(vars(Function), vars(Genre)) +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), size=.25) +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), size=.25) +
  geom_smooth(aes(y = so, color = "so", linetype = "so"), size=.25) +
  geom_smooth(aes(y = pretty, color = "pretty", linetype = "pretty"), size=.25) +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), size=.25) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed","solid", "dashed","solid", "dotdash"),
                        name="Variant",
                        breaks = c("other", "pretty", "really", "so", "very"), 
                        labels = c("other", "pretty", "really", "so", "very")) +
  scale_colour_manual(values=c("indianred4","goldenrod2", "grey30", "indianred4", "goldenrod2"),
                      name="Variant", 
                      breaks=c("other", "pretty", "really", "so", "very"), 
                      labels = c("other", "pretty", "really", "so", "very")) +
  theme_set(theme_light(base_size = 10)) +
  theme(legend.position="top") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE)+
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(Agelbs))
ggsave(file = paste(imAgeDirectory,"Variant_Genre_Function.png",sep="/"))
p2

###############################################################
# p3
p3d <- pd
# remove non-amplified instances
p3d <- p3d[p3d$Amplified != 0,]
# increase n
p3d <- rbind(p3d, p3d, p3d, p3d, p3d, p3d, p3d, p3d, p3d, p3d)
# start plot: all with zero
p3 <- ggplot(p3d, aes(x = jitter(Age), y = very)) +
  facet_grid(vars(Function)) +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), size=.25) +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), size=.25) +
  geom_smooth(aes(y = so, color = "so", linetype = "so"), size=.25) +
  geom_smooth(aes(y = pretty, color = "pretty", linetype = "pretty"), size=.25) +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), size=.25) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed","solid", "dashed","solid", "dotdash"),
                        name="Variant",
                        breaks = c("other", "pretty", "really", "so", "very"), 
                        labels = c("other", "pretty", "really", "so", "very")) +
  scale_colour_manual(values=c("indianred4","goldenrod2", "grey30", "indianred4", "goldenrod2"),
                      name="Variant", 
                      breaks=c("other", "pretty", "really", "so", "very"), 
                      labels = c("other", "pretty", "really", "so", "very")) +
  theme_set(theme_light(base_size = 10)) +
  theme(legend.position="top") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE)+
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(Agelbs))
ggsave(file = paste(imAgeDirectory,"Variant_pd_Function.png",sep="/"), width = 10, height = 10, units = c("cm"),  dpi = 320)
p3

###############################################################
# p4
p4d <- pd
p4d <- subset(p4d, Genre == "PrivateDialogue")
# start plot: all
p4 <- ggplot(p4d, aes(x = jitter(Age), y = zero)) +
  facet_grid(vars(Function)) +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), size=.25) +
  geom_smooth(aes(y = pretty, color = "pretty", linetype = "pretty"), size=.25) +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), size=.25) +
  geom_smooth(aes(y = so, color = "so", linetype = "so"), size=.25) +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), size=.25) +
  geom_smooth(aes(y = zero, color = "zero", linetype = "zero"), size=.25) +  
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed", "dotdash","solid", "dashed", "dotdash","solid"),
                        name="Variant",
                        breaks = c("other", "pretty", "really", "so", "very", "zero"), 
                        labels = c("other", "pretty", "really", "so", "very", "zero")) +
  scale_colour_manual(values=c("indianred4","goldenrod2", "indianred4", "grey30", "goldenrod2", "grey30"),
                      name="Variant", 
                      breaks=c("other", "pretty", "really", "so", "very", "zero"), 
                      labels = c("other", "pretty", "really", "so", "very", "zero")) +
  theme_set(theme_light(base_size = 10)) +
  theme(legend.position="top") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE)+
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(Agelbs))
ggsave(file = paste(imAgeDirectory,"Variantwz_Genre.png",sep="/"), width = 10, height = 10, units = c("cm"),  dpi = 320)
p4

###############################################################
#            WARNING: DATA REDUCTION
###############################################################
#ampnze <- ampnze[ampnze$Function == "Attributive",]
ampnze <- ampnze[ampnze$Genre == "PrivateDialogue",]
# recode adjectives
ntfrqadj <- names(table(ampnze$Adjective))[which(table(ampnze$Adjective) <= 10)]
ampnze$Adjective <- ifelse(ampnze$Adjective %in% ntfrqadj, "other", ampnze$Adjective)
###############################################################
###             TABULARIZATION
###############################################################
# tb 1
Varianttbnze <- table(ampnze$Variant)
Varianttbnze <- Varianttbnze[order(table(ampnze$Variant), decreasing = T)]
Variantnames <- as.vector(names(Varianttbnze))
Variantn <- as.vector(Varianttbnze)
Variantprcnt <- round(Variantn/sum(Variantn)*100, 2)
Variantprcnt2 <-  c(0, round(Variantn[2:length(Variantn)]/sum(Variantn[2:length(Variantn)])*100, 2))
Varianttbnze <- data.frame(Variantnames, Variantn, Variantprcnt, Variantprcnt2)
colnames(Varianttbnze) <- c("Amplifiedensifier", "TokenFrequency", "PercentAgeSlots", "PercentAgeAmplifiedensifiers")
Varianttbnze <- rbind(Varianttbnze, c("Total", sum(as.vector(Varianttbnze$TokenFrequency)), "", ""))
rownames(Varianttbnze) <- NULL
# inspect data
head(Varianttbnze)

# save data to disc
write.table(Varianttbnze, "Varianttbnze.txt", sep = "\t", row.names = F)
###############################################################
# tb 2
ampnzeVariantdf <- ampnze[ampnze$Variant != "0",]
tbAdjective <- table(ampnzeVariantdf$Adjective)
tbAdjective <- tbAdjective[order(tbAdjective, decreasing = T)]
tbAdjective <- tbAdjective[which(tbAdjective >= 5)]
freqAdjective <- names(tbAdjective)
dffAdjective <- ampnzeVariantdf[ampnzeVariantdf$Adjective %in% freqAdjective,]
dffAdjective <- dffAdjective[, c(1, 14, 2)] # Age Variant Adjective
tst5 <- table(ampnze$Variant)
infreqAmplified <- names(tst5[which(tst5 < 15)])
dffAdjective$Variant <- ifelse(dffAdjective$Variant %in% infreqAmplified, "other", dffAdjective$Variant)
dffAdjective <- dffAdjective[dffAdjective$Variant != "0",]
# tabulate by Age
ampnzeVariantdf <- ftable(dffAdjective$Adjective, dffAdjective$Variant, dffAdjective$Age)
ampnzeVariantdf

# save data to disc
write.table(ampnzeVariantdf, "AmplifiedAdjectivedateabs.txt", sep = "\t", row.names = T)
###############################################################
# tb 3
# reorder data
relfreqtb_abs <- ftable(dffAdjective$Adjective, dffAdjective$Age, dffAdjective$Variant)
relfreqtb_rel <- round(prop.table(relfreqtb_abs, 1)*100, 1) # row percentAges
# replave NA with 0
relfreqtb_relwona <- apply(relfreqtb_rel, 1, function(x) ifelse(is.na(x) == T, 0, x))
# convert Amplifiedo matrix
relfreqmx <- matrix(relfreqtb_relwona, ncol = 4, byrow = F)
# convert Amplifiedo data frame
relfreqdf <- as.data.frame(relfreqmx)
# create a vector with Adjectives
Adjective <- rep(unlist(attr(relfreqtb_rel, "row.vars")[1]), each =  
                   (nrow(relfreqmx)/length(unlist(attr(relfreqtb_rel, "row.vars")[1])))) # Age groups
# create a vector with amplifiers
amp <- rep(unlist(attr(relfreqtb_rel, "col.vars")[1]),  
           (nrow(relfreqmx)/length(unlist(attr(relfreqtb_rel, "col.vars")[1])))) # Age groups
# create a data frame
relfreqdf <- data.frame(Adjective, amp, relfreqdf)
# add colnames
colnames(relfreqdf) <- c("Adjective", "amp", unlist(attr(relfreqtb_rel, "row.vars")[2]))
# inspect data
str(relfreqdf); head(relfreqdf)

# save data to disc
write.table(relfreqdf, "iampAdjectiveampAgepcnt.txt",  sep = "\t", row.names = T)
###############################################################
#              SEMANTIC VECTOR SPACE MODEL 2
# evaluation how strongly really and very correlate
# tabulate data
t1 <- tapply(ampnze$Amplified, list(ampnze$Adjective, ampnze$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3nze <- t3
#t3nze <- t3nze[, 2: ncol(t3nze)]
# remove Adjectives that were not Amplifiedensified
t3nze <- t3nze[rowSums(t3nze) > 0, ]
# save row and column names
colnamesnze <- colnames(t3nze)
rownamesnze <- rownames(t3nze)
# turn dataframe Amplifiedo matrix
svsmnze <- as.matrix(t3nze)
# convert token frequency to type frequency
#svsmnze <- apply(svsmnze, 1, Function(x) { x <- ifelse(x > 1, 1, x) } )
svsmnze <- t(svsmnze)
#svsmnze <- svsmnze[, colSums(svsmnze) >= 2]
#svsmnze <- svsmnze[rowSums(svsmnze) >= 2, ]
svsmnze

# determine overall n in data
n_nze <- sum(svsmnze)
n_nze

# correlate amplifiers based on collocation
r_nze <- cor(t(svsmnze))
r_nze

# extract correlation coefficient r for really and very
r_reallyverynze <- r_nze[which(attr(r_nze, "dimnames")[[1]] == "very"), which(attr(r_nze, "dimnames")[[2]] == "really")]
r_reallyverynze

z_reallyverynze <- fisherz(r_reallyverynze)
z_reallyverynze

# the z value can be tested for significance using the r-test from the psych library
#r.test(n=100,r12=.5,r34=.4, n2=80) 
###############################################################
#               PLOTTING LEXICAL DIVESRITY
# Function for extracting lexdiv values
lexdiv <- function(x){
  Varianttokentbnze <- table(x$Variant, x$Adjective)
  Varianttokentbnze <- Varianttokentbnze[2:nrow(Varianttokentbnze), ]
  Varianttokentbnze <- Varianttokentbnze[rowSums(Varianttokentbnze) > 1, ]
  # extract typefrequency of tokenectives
  Varianttokentbnzetyp <- t(apply(Varianttokentbnze, 1, function(x) ifelse(x > 1, 1, x)  ))
  # claculate lexical diversity measure
  lexdivnze <- rowSums(Varianttokentbnzetyp)/rowSums(Varianttokentbnze)
  lexdivnze <- lexdivnze[order(lexdivnze)]
  return(lexdivnze)
}
# apply Function to data
lexdivnze <- lexdiv(ampnze)

# ggplot2 p5
lexdivdf <- data.frame(1:length(lexdivnze), names(lexdivnze), round(lexdivnze, 2))
colnames(lexdivdf) <- c("id", "amp", "lexdiv")

# start plot: Amplified
p5 <- ggplot(lexdivdf, aes(x = jitter(id), y = lexdiv, label = lexdivdf$lexdiv), size = 8) +
  geom_line(aes(y = lexdiv), col = "indianred4", lwd = .5) + 
  geom_smooth(aes(y = lexdiv), size=.5, col = "gray30", lty = "dotted") +
  geom_text(label = lexdivdf$lexdiv, hjust = 0.1, nudge_y = -0.1, size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Amplifier Type", y = "Lexical Diversity") +
  theme_light(base_size = 8) +
  scale_x_continuous(name = "Amplifier Type",
                     breaks = c(1:length(lexdivdf$amp)),
                     labels=lexdivdf$amp)
ggsave(file = paste(imAgeDirectory,"lexdiv_thememinimal.png", sep="/"), width = 15,  height = 7.5, units = c("cm"),  dpi = 320)
p5

###############################################################
#            PLOT LEXICAL DIVERSITY OVER TIME
# extract tokenfrequency of Amplifiedensifiers
a1624 <- subset(ampnze, Age == "16-24")
a2539 <- subset(ampnze, Age == "25-39")
a4049 <- subset(ampnze, Age == "40-49")
a5080<- subset(ampnze, Age == "50+")
# apply Function to data sets
lexdiv1624 <- lexdiv(a1624)
lexdiv2539 <- lexdiv(a2539)
lexdiv4049 <- lexdiv(a4049)
lexdiv5080 <- lexdiv(a5080)
# find common items
cmnamps <- Reduce(intersect, list(names(lexdiv1624),names(lexdiv2539),names(lexdiv4049), names(lexdiv5080)))
# extract lex div values for amps which occur in all Age groups
lexdivvls <- data.frame(lexdiv1624[which(names(lexdiv1624) %in% cmnamps)][order(names(lexdiv1624[which(names(lexdiv1624) %in% cmnamps)]))], 
                        lexdiv2539[which(names(lexdiv2539) %in% cmnamps)][order(names(lexdiv2539[which(names(lexdiv2539) %in% cmnamps)]))], 
                        lexdiv4049[which(names(lexdiv4049) %in% cmnamps)][order(names(lexdiv4049[which(names(lexdiv4049) %in% cmnamps)]))], 
                        lexdiv5080[which(names(lexdiv5080) %in% cmnamps)][order(names(lexdiv5080[which(names(lexdiv5080) %in% cmnamps)]))])# transpose data
lexdivvlst <- t(lexdivvls)
# combine lexdiv tables
p6d <- data.frame(1:length(cmnamps), names(table(ampnze$Age)), lexdivvlst)
colnames(p6d)[1:2] <- c("id", "Age")
rownames(p6d) <- 1:nrow(p6d)
Agelbs <- names(table(p6d$Age))
p6d$Age <- ifelse(p6d$Age == "16-24", 4,
                  ifelse(p6d$Age == "25-39", 3, 
                         ifelse(p6d$Age == "40-49", 2, 
                                ifelse(p6d$Age == "50+", 1, pd$Age))))

# start plot: Amplified
p6 <- ggplot(p6d, aes(x = jitter(Age), y = really, label = Age), size = 8) +
  geom_smooth(aes(y = really, color = "really", lty = "really"), size=.5) +
  geom_smooth(aes(y = very, color = "very", lty = "very"), size=.5) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed","solid"),
                        name="",
                        breaks = c("really",  "very"), 
                        labels = c("really",  "very")) +
  scale_colour_manual(values=c("indianred4", "grey30"),
                      name="", 
                      breaks=c("really",  "very"), 
                      labels = c("really",  "very")) +
  theme(legend.position="top") +
  theme_light(base_size = 15) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Amplifier Type", y = "Lexical Diversity") +
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(Agelbs))
ggsave(file = paste(imAgeDirectory,"lexdivAge_thememinimal.png", sep="/"), width = 15,  height = 7.5, units = c("cm"),  dpi = 320)
p6

###############################################################
#       CORRESPONDENCE ANALYSIS
# load packAges for CA
library("FactoMineR")
library("gplots")
# convert the data as a table
catbnze <- t3nze[rowSums(t3nze) >= 5, ]
catbnze <- catbnze[, colSums(catbnze) >= 2]
# only variable items
catbnze2 <- t(apply(catbnze, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
catbnze <- catbnze[, colSums(catbnze2) >= 2]
# convert data Amplifiedo table
dt <- as.table(as.matrix(catbnze))
# check wehther rows and columns correlate significantly
chisq <- chisq.test(catbnze)
# use only if X2 is significant
chisq$p.value

# create correspondence object
res.ca <- CA(catbnze, graph = FALSE)
# extract eigencvalues
eig.val <- get_eigenvalue(res.ca)
# show percentAge of explained variances per dimension
#fviz_screeplot(res.ca, addlabels = TRUE, ylim = c(0, 60))

# show percentAge of explained variances per dimension with threshold
#fviz_screeplot(res.ca) +  
#geom_hline(yAmplifiedercept=mean(eig.val), linetype=2, color="red")

# Contributions of variables to PC1
#fviz_contrib(res.ca, choice = "row", axes = 1, top = 10)

# Contributions of variables to PC2
#fviz_contrib(res.ca, choice = "col", axes = 2, top = 10)

# plot correspondence analysis results
png("images/CAFactorMapFvisAmpnze.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray")
dev.off()
# alternative plot
png("images/CAFactorMapAmpnze.png",  width = 680, height = 480) 
plot(res.ca, shadow = T, cex = 1, selectRow = "cos2 0.1", selectCol = "cos2 0.1", col.row = "gray50", title = "")
dev.off()
# add dentrogram of results
res <- hcut(catbnze, k = 6, stand = TRUE)
# Visualize
#fviz_dend(res, rect = TRUE, cex = 0.5, k_colors = c("indianred4","goldenrod2", "grey30", "indianred4", "goldenrod2", "grey30"))

# Optimal number of clusters for k-means
my_data <- scale(catbnze)
fviz_nbclust(my_data, kmeans, method = "gap_stat")
###############################################################
#           COVARYING COLLEXEME ANALYSIS
# collex Function
# collex Function
collex <- function(data = data, cv1 = cv1){
  # set up rslttb
  rslttb <- matrix(c("amp", "Adjective", "namp", "nAdjective", "obs", "exp", "prob", "cs", "or", "p"), ncol = 10)
  colnames(rslttb) <- c("Amp", "Adjective", "N(Amp)", "N(Adjective)", "OBS", "EXP", 
                        "Probability", "CollStrength", "OddsRatio", "p")
  rvs <- 1:nrow(t3)
  # define column values
  cv0 <- 1
  # set up table
  sapply(rvs, function(x){
    # extract values
    obs <- t3[x,cv1] # freq Adjective with amp
    fAdjective <- sum(t3[x,]) # freq Adjective
    n <- sum(t3[,cv1]) # freq amp
    fall <- sum(t3) # freq amps and Adjectives
    # calculate exp
    exp <- (fAdjective*n)/fall
    prob <- exp/n
    # create table to extract odds ratios
    m <- matrix(c(obs, (n-obs), (fAdjective-obs), (fall-fAdjective)), nrow = 2, byrow = T)
    o <- fisher.test(m)
    # perform binomial test
    rslt <- binom.test(obs, n, prob)
    # set up table with results
    rslttb <- list(c(colnames(data)[cv1], 
                     rownames(data)[x], 
                     n, 
                     fAdjective,
                     obs,
                     round(exp, 1),
                     round(prob, 2),
                     round(abs(log(as.vector(unlist(rslt$p.value, 10)), 10)), 2), 
                     round(as.vector(unlist(o$estimate)), 2),
                     round(as.vector(unlist(rslt$p.value)), 6)
    ))
    # return results
    return(rslttb)
  } )
}
###############################################################
#                 CCLA ON ALL Age GROUPS COMBINED
# rename data
cclad <- ampnze
# define infrequent adjectives
ntfrqadj <- names(table(cclad$Adjective))[which(table(cclad$Adjective) <= 10)]
# recode adjectives
cclad$Adjective <- ifelse(cclad$Adjective %in% ntfrqadj, "other", cclad$Adjective)

# create table
t1 <- tapply(cclad$Amplified, list(cclad$Adjective, cclad$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3 <- t3[, 2:ncol(t3)]
t3 <- t3[rowSums(t3) > 0, ]
### WARNING!
# apply collex Function (Amplifieds >= 250: 
colnames(t3)[which(colSums(t3) >= 10)]

which(colSums(t3) >= 10)

Amplifiers <- names(which(colSums(t3) >= 10))

pretty  <- collex(data = t3, cv1 = which(colnames(t3) == "pretty"))
really  <- collex(data = t3, cv1 = which(colnames(t3) == "really"))
so  <- collex(data = t3, cv1 = which(colnames(t3) == "so"))
very  <- collex(data = t3, cv1 = which(colnames(t3) == "very"))
# extract informaltion
pretty <- matrix(unlist(pretty),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
collextab <- rbind(pretty,  really, so, very)
# write Function to process collextab (input = collextab)
collextbedit <- function(collextab){
  # convert Amplifiedo data frame
  collexdf <- as.data.frame(collextab)
  # add colnames
  colnames(collexdf) <- c("Amp", "Adjective", "N(Amp)", "N(Adjective)", "OBS", "EXP", 
                          "Probability", "CollStrength", "OddsRatio", "p")
  # add attraction column
  collexdf$attr <- ifelse(as.numeric(collexdf$OBS) > as.numeric(collexdf$EXP), "attr", "repel")
  # modify CollStrength column 
  collexdf$CollStrength <- ifelse(collexdf$attr == "repel", 
                                  paste("-", collexdf$CollStrength, sep =""), collexdf$CollStrength)
  # perform bonferroni correction
  corr05 <- 0.05/nrow(collexdf)
  collexdf$corr05 <- rep(corr05, nrow(collexdf))
  corr01 <- 0.01/nrow(collexdf)
  collexdf$corr01 <- rep(corr01, nrow(collexdf))
  corr001 <- 0.001/nrow(collexdf)
  collexdf$corr001 <- rep(corr001, nrow(collexdf))
  # calculate corrected significance status
  collexdf$sig <- as.vector(unlist(sapply(collexdf$p, function(x){
    x <- ifelse(x <= corr001, "p<.001",
                ifelse(x <= corr01, "p<.01",
                       ifelse(x <= corr001, "p<.001", "n.s."))) } )))
  return(collexdf)
}
# apply collextbedit Function
collex_all <- collextbedit(collextab)
# inspect results
subset(collex_all, sig != "n.s.")

###############################################################
#                  CVCLA : 16-24
# collocation analysis Age 16-24
df1624 <- cclad[cclad$Age == "16-24",]
t1 <- tapply(df1624$Amplified, list(df1624$Adjective, df1624$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3 <- t3[, 2:ncol(t3)]
t3 <- t3[rowSums(t3) > 0, ]
### WARNING!
# apply collex Function (colnames(t3)[which(colSums(t3) >= 10)]
which(colnames(t3) %in% Amplifiers)
pretty  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[1])
really  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[2])
so  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[3])
very  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[4])
# extract informaltion
pretty <- matrix(unlist(pretty),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
collextab <- rbind(pretty, really, so, very)
# apply collextbedit Function
collex1624 <- collextbedit(collextab)
# inspect results
subset(collex1624, sig != "n.s.")

###############################################################
#                  CVCLA: 25-39
df2539 <- cclad[cclad$Age == "25-39",]
t1 <- tapply(df2539$Amplified, list(df2539$Adjective, df2539$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3 <- t3[, 2:ncol(t3)]
t3 <- t3[rowSums(t3) > 0, ]
### WARNING!
# apply collex Function (colnames(t3)[which(colSums(t3) >= 3)]
which(colnames(t3) %in% Amplifiers)
pretty  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[1])
really  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[2])
so  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[3])
very  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[4])
# extract informaltion
pretty <- matrix(unlist(pretty),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
collextab <- rbind(pretty, really, so, very)
# apply collextbedit Function
collex2539 <- collextbedit(collextab)
# inspect results
subset(collex2539, sig != "n.s.")

###############################################################
#                  CVCLA: 40-49
df4049 <- cclad[cclad$Age == "40-49",]
t1 <- tapply(df4049$Amplified, list(df4049$Adjective, df4049$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3 <- t3[, 2:ncol(t3)]
t3 <- t3[rowSums(t3) > 0, ]
### WARNING!
# apply collex Function (colnames(t3)[which(colSums(t3) >= 3)]
which(colnames(t3) %in% Amplifiers)
pretty  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[1])
really  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[2])
so  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[3])
very  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[4])
# extract informaltion
pretty <- matrix(unlist(pretty),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
collextab <- rbind(pretty, really, so, very)
# apply collextbedit Function
collex4049 <- collextbedit(collextab)
# inspect results
subset(collex4049, sig != "n.s.")

###############################################################
#                  CVCLA: 50-80
df5080 <- cclad[cclad$Age == "50+",]
t1 <- tapply(df5080$Amplified, list(df5080$Adjective, df5080$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3 <- t3[, 2:ncol(t3)]
t3 <- t3[rowSums(t3) > 0, ]
### WARNING!
# apply collex Function (colnames(t3)[which(colSums(t3) >= 3)]
which(colnames(t3) %in% Amplifiers)
pretty  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[1])
really  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[2])
so  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[3])
very  <- collex(data = t3, cv1 = which(colnames(t3) %in% Amplifiers)[4])
# extract informaltion
pretty <- matrix(unlist(pretty),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
collextab <- rbind(pretty, really, so, very)
# apply collextbedit Function
collex5080 <- collextbedit(collextab)
# inspect results
subset(collex5080, sig != "n.s.")

###########################################################################
# combine covar collex data frames
collexdfAge <- rbind(collex1624[,c(1:11,15)], collex2539[,c(1:11,15)], 
                     collex4049[,c(1:11,15)], collex5080[,c(1:11,15)])
Age <- c(rep("1624", nrow(collex1624)), rep("2539", nrow(collex2539)), 
         rep("4049", nrow(collex4049)), rep("5080", nrow(collex5080)))
# create data frame
covarcoldf <- data.frame(Age, collexdfAge)
#convert Amplifiedo numeric
covarcoldf$Probability <- as.numeric(covarcoldf$Probability)
covarcoldf$CollStrength <- as.numeric(covarcoldf$CollStrength)
covarcoldf$OddsRatio <- as.numeric(covarcoldf$OddsRatio)
# inspect data
str(covarcoldf); head(covarcoldf); summary(covarcoldf$CollStrength)

###########################################################################
# rename data
p10d <- covarcoldf[, c(1:3, 9)] # Age, Amp, Adjective, CollStrength
#extract Adjectives present in all ages
allagesAdjective <- which(rowSums(ftable(p10d$Amp, p10d$Adjective, p10d$Age)) == 4)
ampAdjectiveageftb <- ftable(p10d$Amp, p10d$Adjective, p10d$Age)
amp <- unlist(attr(ampAdjectiveageftb, "row.vars")[1])
Adjective <- unlist(attr(ampAdjectiveageftb, "row.vars")[2])
Age <- unlist(attr(ampAdjectiveageftb, "col.vars")[1])
Adjectiver <- rep(Adjective, length(amp))
ampr <- rep(amp, each = length(Adjective))
freqAdjectives1 <- unique(Adjectiver[allagesAdjective])
freqAdjectives2 <- names(table(cclad$Adjective)[order(table(cclad$Adjective), decreasing = T)])[1:20]
freqAdjectives <- intersect(freqAdjectives1, freqAdjectives2)
# use only data with Adjectives that are present among all amps and all age groups
p10d <- subset(p10d, Adjective %in% freqAdjectives)
# create new id variable
p10d$AgeAdjective <- paste(p10d$Age, "_", p10d$Adjective, sep = "")
# reorder data frame
p10tb <- reshape(p10d, idvar = "AgeAdjective", timevar = "Amp",direction = "wide")
# select relevant column
# age, Adjective, collstrength:quite, collstrength:so, collstrength:really, collstrength:very
p10tb <- p10tb[, c(2:4, 7, 10, 13)] 
colnames(p10tb) <- c("Age", "Adjective", "pretty", "really", "so", "very")
p10tb$Adjective <- as.factor(p10tb$Adjective)
# recode Age
p10tb$Age <- ifelse(p10tb$Age == "1624", 4, p10tb$Age)
p10tb$Age <- ifelse(p10tb$Age == "2539", 3, p10tb$Age)
p10tb$Age <- ifelse(p10tb$Age == "4049", 2, p10tb$Age)
p10tb$Age <- ifelse(p10tb$Age == "5080", 1, p10tb$Age)
p10tb$Age <- as.numeric(p10tb$Age)
# create vector with Age groups
Agelbs <- c("16-24", "25-39", "40-49", "50-80")
head(p10tb)

# start plot: all
p10 <- ggplot(p10tb, aes(x = jitter(Age), y = pretty)) +
  facet_grid(vars(Adjective)) +
  geom_smooth(aes(y = pretty, color = "pretty", linetype = "pretty"), size=.5, se = F) +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), size=.5, se = F) +
  geom_smooth(aes(y = so, color = "so", linetype = "so"), size=.5, se = F) +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), size=.5, se = F) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed","twodash", "longdash", "solid"),
                        name="Variant",
                        breaks = c("pretty", "really", "so", "very"), 
                        labels = c("pretty", "really", "so", "very")) +
  scale_colour_manual(values=c("grey70","goldenrod2", "indianred4", "grey30"),
                      name="Variant", 
                      breaks=c("pretty", "really", "so", "very"), 
                      labels = c("pretty", "really", "so", "very")) +
  theme_set(theme_light(base_size = 8)) +
  theme(legend.position="top") +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  labs(x = "Year", y = "Collocation Strength (LOG(p), 10)") +
  guides(size = FALSE)+
  guides(alpha = FALSE)+
  #  guides(linetype = FALSE)+
  #theme(legend.title=element_blank())+
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(Agelbs))
ggsave(file = paste(imAgeDirectory,"covarcoll_coha_freqAdjective1.png",sep="/"), width = 7.5, height = 20, units = c("cm"),  dpi = 320)
p10

###########################################################################
#                  CHANGES IN Adjective FREQ
# tabulate data
ftbadjnze <- ftable(cclad$Adjective, cclad$Age)
rwnms <- as.vector(unlist(attr(ftbadjnze, "row.vars")))
ftbadjnze <- ftbadjnze[2:nrow(ftbadjnze),]
rownames(ftbadjnze) <- rwnms[2:length(rwnms)]
svrwnms <- as.vector(unlist(attr(ftbadjnze, "dimnames")))[which(rowSums(ftbadjnze) >= 35)]
ftbadjnze <- ftbadjnze[which(rowSums(ftbadjnze) >= 35),]
rownames(ftbadjnze) <- svrwnms
colnames(ftbadjnze) <- names(table(cclad$Age))
ptbadjnze <- prop.table(ftbadjnze, margin=2)*100
#ptbampnze <- ptbampnze[rowSums(ptbampnze) > 1, ]
ptbadjnze

# save data to disc
write.table(ptbadjnze, "ptbadjnze.txt", sep = "\t", row.names = F)
###########################################################################
p11d <- cclad
famp <- names(table(p11d$Adjective))[which(table(p11d$Adjective) > 35)]
p11d$Adjective <- ifelse(p11d$Adjective %in% famp, p11d$Adjective, "other")
#table(p11d$Adjective)[order(table(p11d$Adjective), decreasing = T)]

# create vars for Variant
p11d$other <- ifelse(p11d$Adjective == "other", 100, 0)
p11d$big <- ifelse(p11d$Adjective == "big", 100, 0)
p11d$good <- ifelse(p11d$Adjective == "good", 100, 0)
p11d$new <- ifelse(p11d$Adjective == "new", 100, 0)
p11d$nice <- ifelse(p11d$Adjective == "nice", 100, 0)
p11d$old <- ifelse(p11d$Adjective == "old", 100, 0)
p11d$right <- ifelse(p11d$Adjective == "right", 100, 0)
# recode Age
p11d$Age <- ifelse(p11d$Age == "16-24", 4, p11d$Age)
p11d$Age <- ifelse(p11d$Age == "25-39", 3, p11d$Age)
p11d$Age <- ifelse(p11d$Age == "40-49", 2, p11d$Age)
p11d$Age <- ifelse(p11d$Age == "50+", 1, p11d$Age)
p11d$Age <- as.numeric(p11d$Age)
# create vector with Age groups
Agelbs <- c("16-24", "25-39", "40-49", "50-80")
head(p11d)

table(p11d$Adjective)[order(table(p11d$Adjective), decreasing = T)]

# start plot: Adjective
p11 <- ggplot(p11d, aes(x = jitter(Age), y = other, label = Age), size = 8) +
  geom_smooth(aes(y = other, color = "other", lty = "other"), size=.5) +
  geom_smooth(aes(y = good, color = "good", lty = "good"), size=.5) +
  geom_smooth(aes(y = nice, color = "nice", lty = "nice"), size=.5) +
  geom_smooth(aes(y = new, color = "new", lty = "new"), size=.5) +
  geom_smooth(aes(y = big, color = "big", lty = "big"), size=.5) +
  geom_smooth(aes(y = old, color = "old", lty = "old"), size=.5) +
  geom_smooth(aes(y = right, color = "right", lty = "right"), size=.5) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed", "longdash", "solid","dotted", "dashed","longdash", "solid"),
                        name="",
                        breaks = c("other", "good", "nice", "new", "big", "old", "right"), 
                        labels = c("other", "good", "nice", "new", "big", "old", "right")) +
  scale_colour_manual(values=c("grey30", "grey60", "goldenrod2",  "indianred4", "grey30", "goldenrod2", "grey60"),
                      name="", 
                      breaks=c("other", "good", "nice", "new", "big", "old", "right"), 
                      labels = c("other", "good", "nice", "new", "big", "old", "right")) +
  theme(legend.position="top") +
  theme_light(base_size = 8) +
  coord_cartesian(ylim = c(0, 70)) +
  labs(x = "Decade", y = "Percent of Adjectives") +
  scale_x_continuous(name = "Age",
                     breaks = 4:1,
                     labels=Agelbs)
ggsave(file = paste(imAgeDirectory,"AdjectivefreqAge.png", sep="/"), width = 15,  height = 7.5, units = c("cm"),  dpi = 320)
p11

###########################################################################
AdjectiveAge <- 1:4
Adjectivelm <- ptbadjnze
head(Adjectivelm)

str(Adjectivelm)

nrow(Adjectivelm)
sigAdjective <- apply(Adjectivelm, 1, function(x){
  x <- lm(x ~ AdjectiveAge)
  x <- summary(x)[4][[1]][[8]]})

sigAdjectives <- which(sigAdjective < .05)
sigAdjectives

###########################################################################
#                  REGRESSION DATA SET
###########################################################################
ampnze <- ampnze[ampnze$Amplified == 1,]
# remove superfluous columns
ampnze$ID <- NULL
ampnze$SpeechUnit <- NULL
ampnze$CleanSpeechUnit <- NULL
ampnze$PosTaggedSpeechUnit <- NULL
ampnze$OriginalString <- NULL
ampnze$PreContext <- NULL
ampnze$PostContext <- NULL
ampnze$speech.unit.count <- NULL
ampnze$word.count <- NULL
ampnze$txtspk <- NULL
ampnze$PreContextLong <- NULL
ampnze$Genre <- NULL
# inspect data
nrow(ampnze); str(ampnze)

###############################################################
write.table(ampnze, "icenzeamp03_regdat.txt", row.names= F, sep = "\t")
###############################################################
#                   END PART 3
###############################################################
