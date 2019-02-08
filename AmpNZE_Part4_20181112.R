##################################################################
# Titel:      The Amplifier System of New Zealand English - Part 4
# R version:  3.5.1 (2018-07-02) -- "Feather Spray"
# Autor:      Martin Schweinberger
# Date:       2018-11-06
# Contact:    martin.schweinberger.hh@gmail.com
##################################################################
# Disclaimer: If you have questions,suggestions or you found errors
#             or in case you would to provide feedback, questions
#             write an email to martin.schweinberger.hh@gmail.com.
# Citation:   If you use this script or results thereof, please cite it as:
#             Schweinberger, Martin. 2018. "The Amplifier System of New Zealand English, Part 4",
#             unpublished R script, The University of Queensland.
###############################################################
#                   START
###############################################################
# remove all lists from the current workspace
rm(list=ls(all=T))
# set wd
setwd("D:\\Uni\\Projekte\\02-Intensification\\AmpNZE")
# load libraries
library(Boruta)
library(effects)
library(dplyr)
library(ggplot2)
library(gsubfn)
library(Hmisc)
library(MASS)
library(mlogit)
#library(nlme)
library(languageR)
library(lme4)
library(plyr)
library(QuantPsyc)
library(reshape)
#install.packages("D:/R/Rling_1.0.tar.gz", repos = NULL, type = "source")
library(Rling)
library(RLRsim)
library(rms)
library(sjPlot)
library(visreg)
# load self written function
source("D:\\R/multiplot_ggplot2.R") # for multiple ggplot2 plots in one window
source("D:\\R/PseudoR2lmerBinomial.R")
source("D:\\R/mlr.summary.R")
source("D:\\R/blr.summary.R")
source("D:\\R/meblr.summary.R")
source("D:\\R/ModelFittingSummarySWSU.R") # for Mixed Effects Model fitting (step-wise step-up): Binary Logistic Mixed Effects Models
source("D:\\R/ModelFittingSummarySWSD.R") # for Mixed Effects Model fitting (step-wise step-down): Binary Logistic Mixed Effects Models
source("D:\\R/ModelFittingSummarySWSULogReg.R") # for Fixed Effects Model fitting: Binary Logistic Models
###############################################################
# Setting options
options(stringsAsFactors = F)
options(scipen = 999)
options(max.print=10000)
# define image dausctors
imageDirectory<-"images"
###############################################################
###                   START
###############################################################
# load data
reallynze <- read.table("icenzeamp03_regdat.txt", sep = "\t", header = T)
###############################################################
# inspect data
str(reallynze); nrow(reallynze)

# recode variables 
reallynze$AudienceSize <- ifelse(reallynze$AudienceSize == "2" | reallynze$AudienceSize == "3", "2-3", reallynze$AudienceSize)
reallynze$AudienceSize <- ifelse(reallynze$AudienceSize == "4" | reallynze$AudienceSize == "5" | reallynze$AudienceSize == "6", "4+", reallynze$AudienceSize)
# factorize variables
clfct <- c("FileSpeaker", "Function", "Gender", "Occupation", "Emotionality", "Gradabilty", "Ethnicity", 
           "SemanticCategory", "Priming", "ConversationType", "AudienceSize", "Age")
reallynze[clfct] <- lapply(reallynze[clfct], factor)
# inspect tabulated data
lapply(reallynze[clfct], table)

tbampnze <- lapply(reallynze[clfct], table)
lapply(tbampnze, sum)

# convert dep var into numeric
reallynze$really <- as.numeric(reallynze$really)
# test which variables can be excluded
sum(table(reallynze$age)); nrow(reallynze)

# inspect data
str(reallynze); head(reallynze)

# create frequent adjectives column for random intercepts
tbadj <- table(reallynze$Adjective)
fradj <- names(tbadj)[tbadj >= 10] 
reallynze$Adjective <- as.vector(unlist(sapply(reallynze$Adjective, function(x){
  x <- ifelse(x %in% fradj, x, "other")})))
reallynze$Adjective <- factor(reallynze$Adjective)
# inspect adj
table(reallynze$Adjective)

# provide an overview of the data
head(reallynze); str(reallynze)

# remove superfluous variables
reallynze$Variant <- NULL
reallynze$File <- NULL
reallynze$Subfile <- NULL
reallynze$Speaker <- NULL
reallynze$very <- NULL
reallynze$very <- NULL
reallynze$so <- NULL
reallynze$pretty <- NULL
reallynze$Amplified <- NULL
reallynze$AgeOriginalClassification <- NULL

################################################################
#               CONDITION INFERENCE TREES
library(partykit)
# create data
citd <- reallynze
# set.seed
set.seed(201811081) 
# apply bonferroni correction (1 minus alpha multiplied by n of predictors)
control = ctree_control(mincriterion = 1-(.05*14))
# create initial conditional inference tree model
citd.ctree <- ctree(really ~ Age + Adjective + FileSpeaker + Function + Priming + Gender + 
                      Ethnicity + Occupation + ConversationType + AudienceSize + Freq + 
                      Gradabilty + SemanticCategory + Emotionality,
                    data = citd)
# plot final ctree
png("images/final_ctree.png",  width = 680, height = 480) 
plot(citd.ctree, gp = gpar(fontsize = 8))
dev.off()
# test prediction accuracy
ptb <- table(predict(citd.ctree), citd$really)
(((ptb[1]+ptb[4])+(ptb[2]+ptb[3]))/sum(table(predict(citd.ctree), citd$really)))*100
##100

# determine baseline
(table(citd$really)[[2]]/sum(table(citd$really)))*100
## 45.53

###############################################################
#                   RANDOM FOREST I
# prepare data
rfd <- reallynze
# remove FileSpeaker from data (too many levels)
rfd$FileSpeaker <- NULL
# convert really into a factor
rfd$really <- as.factor(rfd$really)
# start with random forest
# set seed
set.seed(111)
# partition data for evaluating rf 
id <- sample(2, nrow(rfd), replace = T, prob = c(.7, .3))
train <- rfd[id == 1, ]
test <- rfd[id == 2,]
# load library
library(randomForest)
# create initial model
reallynze_rf1 <- randomForest(really~., data = train)
# inspect model
print(reallynze_rf1)

# inspect attibutes
attributes(reallynze_rf1)

# start model evaluation
# install package
#source("https://bioconductor.org/biocLite.R"); biocLite(); library(Biobase)
#install.packages("Biobase", repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com", 
#                                      "http://cran.rstudio.com/", dependencies=TRUE))
#install.packages("dimRed", dependencies = TRUE)
#install.packages('caret', dependencies = TRUE)

# load caret library
library(caret) # because initially caret did not work, the libraries above had to be installed
# extract prediction for training data
ptrain1 <- predict(reallynze_rf1, train)
# inspect predictions
head(ptrain1); head(train$really)

# create confusionMatrix
confusionMatrix(ptrain1, train$really)

# extract prediction for test data
ptest1 <- predict(reallynze_rf1, test)
# create confusionMatrix
confusionMatrix(ptest1, test$really)

# determine errorrate of random forest model
plot(reallynze_rf1, main = "")

# tune model
reallynze_rf2 <- tuneRF(train[, !colnames(train)== "really"], train[, colnames(train)== "really"], 
       stepFactor = .1, # for most values 3 appears to be optimal
       plot = T,
       ntreeTry = 1000,
       trace = T,
       improve = .05
       )
# create improved model
# create initial model
reallynze_rf2 <- randomForest(really~., data = train, 
                              ntree = 1000,
                              ntry = 3,
                              importance= T,
                              proximity = T)
# inspect model
print(reallynze_rf2)

# predict based on improved model
ptrain2 <- predict(reallynze_rf2, train)
# create confusionMatrix
confusionMatrix(ptrain2, train$really)

# extract prediction for test data
ptest2 <- predict(reallynze_rf2, test)
# create confusionMatrix
confusionMatrix(ptest2, test$really)

# inspect number of nodes for trees
hist(treesize(reallynze_rf2), main = "", col = "lightgray")

# check variable importance
varImpPlot(reallynze_rf2, main = "") 

# left plot (Accuracy): how much accuracy decreases if factor is left out
# left plot (Gini/Pureness): how much more unpure (ambigious) the distributions become if fector is left out
# extract variable importance values
importance(reallynze_rf2)

#which variables have been used in the trees
varUsed(reallynze_rf2)

# partial dependence plot
partialPlot(reallynze_rf2, train, Freq, 1)

partialPlot(reallynze_rf2, train, ConversationType, 1)
            
partialPlot(reallynze_rf2, train, Age, 1)

partialPlot(reallynze_rf2, train, Gender, 1)


# extract tree
getTree(reallynze_rf2, 1, labelVar = T)

# mds plot
MDSplot(reallynze_rf2, test$really)

###############################################################
#                   RANDOM FOREST I
# detach partykit
detach("package:partykit", unload=TRUE)
# load package party
library(party)
# prepare data
rfd <- reallynze
# set seed
set.seed(222)

# create initial model
reallynze.rf <- cforest(really ~ really ~ Age + Adjective + FileSpeaker + Function + Priming + Gender + 
                          Ethnicity + Occupation + ConversationType + AudienceSize + Freq + 
                          Gradabilty + SemanticCategory + Emotionality,
                        data = rfd, controls = cforest_unbiased(ntree = 50, mtry = 3))
# determine importance of factors
reallynze.varimp <- varimp(reallynze.rf, conditional = T)
round(reallynze.varimp, 3)

# plot result
png("images/randomforest_factorimportance.png",  width = 680, height = 480) 
dotchart(sort(reallynze.varimp), pch = 20, main = "Conditional importance of variables")
dev.off()

library(Hmisc)
# evaluate random forst
reallynze.rf.pred <- unlist(treeresponse(reallynze.rf))[c(FALSE,TRUE)]
somers2(reallynze.rf.pred, as.numeric(rfd$really) - 1)
##           C          Dxy            n      Missing 
##      0.7216286    0.4432572 2264.0000000    0.0000000 
###############################################################
#                     RANDOM FOREST 2
# load library
library(party)
# create data
randomforestdata <- reallynze

cf1 <- cforest(really ~ . , data= randomforestdata, control=cforest_unbiased(mtry=2,ntree=100)) # fit the random forest
varimp(cf1) # get variable importance, based on mean decrease in accuracy

varimp(cf1, conditional=TRUE) # conditional=True, adjusts for correlations between predict

varimpAUC(cf1)  # more robust towards class imbalance.

###############################################################
#                  BORUTA
# create dada for boruta
borutadata <- reallynze
# run 1
boruta.ampnze <- Boruta(really~.,data=borutadata)
print(boruta.ampnze)

getConfirmedFormula(boruta.ampnze)

png("images/BorutaAmpnze.png",  width = 1500, height = 300)
plot(boruta.ampnze, cex = .75)
dev.off()
plot(boruta.ampnze)

png("images/BorutaAmpnzeHistory.png",  width = 680, height = 480)
plotImpHistory(boruta.ampnze)
dev.off()
plotImpHistory(boruta.ampnze)

# remove superfluous variables
borutadata$Emotionality <- NULL
borutadata$Priming <- NULL
borutadata$Gradabilty <- NULL
borutadata$Gradabilty <- NULL
borutadata$Ethnicity <- NULL
borutadata$Occupation <- NULL
# run2
boruta.ampnze <- Boruta(really~.,data=borutadata)
print(boruta.ampnze)

getConfirmedFormula(boruta.ampnze)

png("images/BorutaAmpnze.png",  width = 1500, height = 300)
plot(boruta.ampnze, cex = .75)
dev.off()
plot(boruta.ampnze)

png("images/BorutaAmpnzeHistory.png",  width = 680, height = 480)
plotImpHistory(boruta.ampnze)
dev.off()
plotImpHistory(boruta.ampnze)

# remove superfluous variables
borutadata$Function <- NULL
borutadata$SemanticCategory <- NULL
# run3
boruta.ampnze <- Boruta(really~.,data=borutadata)
print(boruta.ampnze)

getConfirmedFormula(boruta.ampnze)

png("images/BorutaAmpnze.png",  width = 1500, height = 750)
plot(boruta.ampnze, cex.axis=.7, las=2, xlab="", col = c("grey30", "grey30", "grey30",  "grey30", "grey30", "grey30", "grey30","grey90","grey90","grey90"))
dev.off()
plot(boruta.ampnze)

png("images/BorutaAmpnzeHistory.png",  width = 680, height = 480)
plotImpHistory(boruta.ampnze)
dev.off()
plotImpHistory(boruta.ampnze)






###############################################################
# prepare plot data (pd) 
pd <- reallynze
pd$really <- pd$really *100
agelbs <- names(table(pd$age))

# start plotting
p10 <- ggplot(pd, aes(gender, really, color = gender)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Gender", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
p10

# activate (remove #) to show
p11 <- ggplot(pd, aes(gender, really, color = gender)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Gender", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyGender.png",sep="/"))
p11

p12 <- ggplot(pd, aes(age, really, color = age)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Age", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyAge.png",sep="/"))
p12

# prepare data for p13
p13d <- pd
p13d$age <- as.numeric(ifelse(p13d$age == "17-25", 4,
                              ifelse(p13d$age == "26-30", 3, 
                                     ifelse(p13d$age == "31-40", 2, 
                                            ifelse(p13d$age == "41-80", 1, p13d$age)))))
p13 <- ggplot(p13d, aes(x = jitter(age), y = really)) +
  geom_smooth(aes(y = really), size=1, col = "gray40", lty = "longdash") +
  geom_smooth(method='lm', se = FALSE, size=1, col = "blue", lty = "dotted") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  theme_light(base_size = 20) +
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(agelbs))
ggsave(file = paste(imageDirectory,"ReallyAgeSmooth.png",sep="/"))
p13

p14 <- ggplot(pd, aes(priming, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Priming", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyPriming.png",sep="/"))
p14

p15 <- ggplot(pd, aes(fun, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Syntactic function", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyFun.png",sep="/"))
p15

p16 <- ggplot(pd, aes(emo, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Emotionality of Adj.", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyEmo.png",sep="/"))
p16

p17 <- ggplot(pd, aes(grad, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Gradability of Adj.", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyGrad.png",sep="/"))
p17

p18 <- ggplot(pd, aes(sem, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Semantic class of Adj.", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50", "grey50", "grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallySem.png",sep="/"), width = 14)
p18

p19 <- ggplot(pd, aes(smsx, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Conversation type.", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50", "grey50", "grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyConvTyp.png",sep="/"), width = 14)
p19

# start plot: audience size
p20 <- ggplot(pd, aes(audiencen, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_light(base_size = 20)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Conversation type.", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyAudienceSize.png",sep="/"), width = 14)
p20

p21 <- ggplot(pd, aes(fradj, really)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  theme_set(theme_bw(base_size = 21)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position="none") +
  labs(x = "Semantic class of Adj.", y = "Percent (REALLY in Ampl. Adj. Slots)") +
  scale_color_manual(values = c("grey50", "grey50", "grey50", "grey50", "grey50"))
ggsave(file = paste(imageDirectory,"ReallyFradj.png",sep="/"), width = 14)
p21

# interaction plots
# p100
pd100 <- data.frame(pd$age, pd$gender, pd$really)
colnames(pd100) <- gsub("pd.", "", colnames(pd100))
pd100 <- na.omit(pd100)
p100 <- ggplot(pd100, aes(age, really, colour = gender)) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 1.25) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_set(theme_light(base_size = 20)) +
  #  theme(legend.position="none") +
  labs(x = "Gender", y = "Percent (REALLY in Ampl. Adj. Slots)", colour = "gender") +
  scale_color_manual(values = c("grey40", "grey60", "grey80"))
ggsave(file = paste(imageDirectory,"ReallyAgeGender.png",sep="/"))
p100

# prepare data p101
p101d <- tapply(pd$really, list(pd$age, pd$gender), mean)
p101d <- data.frame(rownames(p101d), p101d)
colnames(p101d) <- c("age", "Women", "Men")
p101d$age <- as.numeric(ifelse(p101d$age == "17-25", 4,
                               ifelse(p101d$age == "26-30", 3, 
                                      ifelse(p101d$age == "31-40", 2, 
                                             ifelse(p101d$age == "41-80", 1, p101d$age)))))
# increase n
p101d <- p101d[rep(1:nrow(p101d), 3), ]
# p101
p101 <- ggplot(p101d, aes(jitter(age), jitter(Women))) +
  geom_smooth(aes(y = jitter(Women), color = "Women", linetype = "Women"), size=1) +
  geom_smooth(aes(y = jitter(Men), color = "Men", linetype = "Men"), size=1) +
  geom_smooth(aes(y = jitter(Women), color = "Women (Reg. Line)", linetype = "Women (Reg. Line)"), method='lm', se = FALSE, size=1) +
  geom_smooth(aes(y = jitter(Men), color = "Men (Reg. Line)", linetype = "Men (Reg. Line)"), method='lm', se = FALSE, size=1) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed","dotted", "solid", "dotted"),
                        name="Gender",
                        breaks = c("Women", "Men", "Women (Reg. Line)", "Men (Reg. Line)"), 
                        labels = c("Women", "Men", "Women (Reg. Line)", "Men (Reg. Line)")) +
  scale_colour_manual(values=c("grey40", "blue", "grey40", "red"),
                      name="Gender", 
                      breaks=c("Women", "Men", "Women (Reg. Line)", "Men (Reg. Line)"), 
                      labels = c("Women", "Men", "Women (Reg. Line)", "Men (Reg. Line)")) +
  theme_set(theme_light(base_size = 10)) +
  theme(legend.position="top") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE)+
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4),
                     labels=rev(agelbs))
ggsave(file = paste(imageDirectory,"ReallyAgeGenderSmooth.png",sep="/"))
p101

# p102
pd102 <- data.frame(pd$sem, pd$gender, pd$really)
colnames(pd102) <- gsub("pd.", "", colnames(pd102))
pd102 <- na.omit(pd102)
p102 <- ggplot(pd102, aes(sem, really, colour = gender)) +
  stat_summary(fun.y = mean, geom = "point", size = 1) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = .5) +
  coord_cartesian(ylim = c(0, 102)) +
  theme_set(theme_light(base_size = 10)) +
  #  theme(legend.position="none") +
  labs(x = "Semantic Category", y = "Percent (REALLY in Ampl. Adj. Slots)", colour = "gender") +
  scale_color_manual(values = c("grey40", "grey80"))
ggsave(file = paste(imageDirectory,"ReallyGenderSem.png",sep="/"), width = 15, height = 7.5, units = c("cm"))
p102

# p103
pd103 <- data.frame(pd$grad, pd$audiencen, pd$really)
colnames(pd103) <- gsub("pd.", "", colnames(pd103))
pd103 <- na.omit(pd103)
p103 <- ggplot(pd103, aes(grad, really, colour = audiencen)) +
  stat_summary(fun.y = mean, geom = "point", size = 1) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = .5) +
  coord_cartesian(ylim = c(0, 103)) +
  theme_set(theme_light(base_size = 10)) +
  #  theme(legend.position="none") +
  labs(x = "Gradability", y = "Percent (REALLY in Ampl. Adj. Slots)", colour = "Audience Size") +
  scale_color_manual(values = c("grey40", "grey80"))
ggsave(file = paste(imageDirectory,"ReallyGradAudienceSize.png",sep="/"), width = 15, height = 7.5, units = c("cm"))
p103

# p104
pd104 <- data.frame(pd$smsx, pd$audiencen, pd$really)
colnames(pd104) <- gsub("pd.", "", colnames(pd104))
pd104 <- na.omit(pd104)
p104 <- ggplot(pd104, aes(smsx, really, colour = audiencen)) +
  stat_summary(fun.y = mean, geom = "point", size = 1) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = .5) +
  coord_cartesian(ylim = c(0, 104)) +
  theme_set(theme_light(base_size = 10)) +
  #  theme(legend.position="none") +
  labs(x = "Conversation Type", y = "Percent (REALLY in Ampl. Adj. Slots)", colour = "Audience Size") +
  scale_color_manual(values = c("grey40", "grey80"))
ggsave(file = paste(imageDirectory,"ReallySmsxAudienceSize.png",sep="/"), width = 15, height = 7.5, units = c("cm"))
p104

###############################################################
#        Fixed Effects Binomial Logistic Regression
###############################################################
# set options
options(contrasts  =c("contr.treatment", "contr.poly"))
reallynze.dist <- datadist(reallynze)
options(datadist = "reallynze.dist")
# a few words on glm vs lrm: Baayen (2008:196-197) states that lrm should be
# the function of choice in cases where each row contains
# exactly 1 success OR failure (1 or 0) while glm is preferrable if there are two
# columns holding the number of successes and the number of failures
# respectively. i have tried it both ways and both functions work fine if
# each row contains exactly 1 success OR failure but only glm can handle the
# latter case.
# generate initial saturated regression model including
# all variables and their interactions
m0.glm = glm(really ~ 1, family = binomial, data = reallynze) # baseline model glm
m0.lrm = lrm(really ~ 1, data = reallynze, x = T, y = T) # baseline model lrm
# inspect results
summary(m0.glm)

m0.lrm

###############################################################
#        Mixed Effects Binomial Logistic Regression
###############################################################
# create model with a random intercept for token
m0.lmer <- lmer(really ~ (1|fradj), data = reallynze, family = binomial)
# Baayen (2008:278-284) uses the call above but the this call is now longer
# up-to-date because the "family" parameter is deprecated
# we switch to glmer (suggested by R) instead but we will also
# create a lmer object of the final minimal adequate model as some functions
# will not (yet) work on glmer
m0.glmer = glmer(really ~ (1|fradj), data = reallynze, family = binomial)

# results of the lmer object
print(m0.lmer, corr = F)

# check if including the random effect is permitted by comparing the aic from the glm to aic from the glmer model
aic.glmer <- AIC(logLik(m0.glmer))
aic.glm <- AIC(logLik(m0.glm))
aic.glmer; aic.glm

# the aic of the glmer object is smaller which shows that including the random
# intercepts is justified

# test random effects
null.id = -2 * logLik(m0.glm) + 2 * logLik(m0.glmer)
pchisq(as.numeric(null.id), df=1, lower.tail=F) # sig m0.glmer better than m0.glm

# inspect results
summary(m0.glm)

summary(m0.glmer)

###########################################################################
# model fitting
# fit the model to find the "best" model, i.e. the minimal adequate model
# we will use a step-wise step up procedure
# we need to add "control = glmerControl(optimizer = "bobyqa")" 
# because otherwise R fails to converge
#	manual modelfitting
m0.glmer <- glmer(really ~ 1+ (1|fradj), family = binomial, data = reallynze, 
                  control=glmerControl(optimizer="bobyqa"))
# add age
m1.glmer <- update(m0.glmer, .~.+age)
m1.glm <- update(m0.glm, .~.+age)
vif(m1.glm) # VIFs ok
anova(m0.glmer, m1.glmer, test = "Chi") # SIG (p<.05)

# add gender
ftable(reallynze$grad, reallynze$really)
m2.glmer <- update(m1.glmer, .~.+gender)
m2.glm <- update(m1.glm, .~.+gender)
vif(m2.glm) # VIFs ok
anova(m1.glmer, m2.glmer, test = "Chi") # only mar sig (p=.05) # check AIC(-2), BIC(+2), logLik(-2): no reason to include gender

# add fun
#ftable(reallynze$fun, reallynze$really)
m3.glmer <- update(m1.glmer, .~.+fun)
#m3.glm <- update(m1.glm, .~.+fun)
#vif(m3.glm) # VIFs ok
#anova(m1.glmer, m3.glmer, test = "Chi") # not sig (p=.81)

# add emo
#ftable(reallynze$emo, reallynze$really)
m4.glmer <- update(m1.glmer, .~.+emo)
#m4.glm <- update(m1.glm, .~.+emo)
#vif(m4.glm) # VIFs ok
#anova(m1.glmer, m4.glmer, test = "Chi") # not sig (p=.72)

# add grad
#ftable(reallynze$sem, reallynze$really)
m5.glmer <- update(m1.glmer, .~.+grad)
#m5.glm <- update(m1.glm, .~.+grad)
#vif(m5.glm) # VIFs ok
#anova(m1.glmer, m5.glmer, test = "Chi") # only mar sig. (p<.059)

# add sem
#ftable(reallynze$sem, reallynze$really)
m6.glmer <- update(m1.glmer, .~.+sem)
#m6.glm <- update(m1.glm, .~.+sem)
#vif(m6.glm) # VIFs ok
#anova(m1.glmer, m6.glmer, test = "Chi") # not sig (p=.74)

# add smsx
#ftable(reallynze$smsx, reallynze$really)
m7.glmer <- update(m1.glmer, .~.+smsx)
#m7.glm <- update(m1.glm, .~.+smsx)
#vif(m7.glm) # VIFs ok
#anova(m1.glmer, m7.glmer, test = "Chi") # not sig (p=.39)

# add audiencesize
#ftable(reallynze$audiencen, reallynze$really)
m8.glmer <- update(m1.glmer, .~.+audiencen)
#m8.glm <- update(m1.glm, .~.+audiencen)
#vif(m8.glm) # VIFs ok
#anova(m1.glmer, m8.glmer, test = "Chi") # not sig (p=.24)


# add priming
ftable(reallynze$priming, reallynze$really)
m9.glmer <- update(m1.glmer, .~.+priming)
m9.glm <- update(m1.glm, .~.+priming)
vif(m9.glm) # VIFs ok
anova(m1.glmer, m9.glmer, test = "Chi") # SIG (p=.00088)

meblrm_m9 <- meblrm.summary(m0.glm, m9.glm, m0.glmer, m9.glmer, reallynze$really) #
meblrm_m9

###########################################################################
# find all 2-way interactions
#install.packages("utils")
library(utils)
vars <- c("age", "gender", "fun", "emo", "grad", "sem", "smsx", "audiencen", "priming")
intac <- t(combn(vars, 2))
intac

# add interactions
# add age*gender
ftable(reallynze$age, reallynze$gender)		
m10.glmer <- update(m9.glmer, .~.+age*gender)
m10.glm <- update(m9.glm, .~.+age*gender)
vif(m10.glm) # VIFs ok: max 3.8 # Zuur et al (2009: 386): cut off point either 3 or 5 (here 5)
anova(m10.glmer, m9.glmer, test = "Chi") # SIG (p=.027)

meblrm_m10 <- meblrm.summary(m0.glm, m10.glm, m0.glmer, m10.glmer, reallynze$really) #
meblrm_m10

# add age*fun
#ftable(reallynze$age, reallynze$fun)			
m11.glmer <- update(m10.glmer, .~.+age*fun)
#m11.glm <- update(m10.glm, .~.+age*fun)
#vif(m11.glm) # VIFs too high (max = 6.5)
#anova(m11.glmer, m10.glmer, test = "Chi") # not sig (p=.19)

# add age*emo
#ftable(reallynze$age, reallynze$emo)			
m12.glmer <- update(m10.glmer, .~.+age*emo)
#m12.glm <- update(m10.glm, .~.+age*emo)
#vif(m12.glm) # VIFs too high (max = )
#anova(m12.glmer, m10.glmer, test = "Chi") # not sig (p=.92)

# add age*grad
#ftable(reallynze$age, reallynze$grad)			
m13.glmer <- update(m10.glmer, .~.+age*grad)
#m13.glm <- update(m10.glm, .~.+age*grad)
#vif(m13.glm) # VIFs too high (max = )
#anova(m13.glmer, m10.glmer, test = "Chi") # not sig (p=.82)

# add age*sem
#ftable(reallynze$age, reallynze$sem)			
m14.glmer <- update(m10.glmer, .~.+age*sem) # WARNING: Model is nearly unidentifiable: large eigenvalue ratio
#m14.glm <- update(m10.glm, .~.+age*sem)
#vif(m14.glm) # VIFs too high (max = )
#anova(m14.glmer, m10.glmer, test = "Chi") # not sig (p=.51)

# add age*smsx
#ftable(reallynze$age, reallynze$smsx)			
m15.glmer <- update(m10.glmer, .~.+age*smsx)
#m15.glm <- update(m10.glm, .~.+age*smsx)
#vif(m15.glm) # VIFs too high (max = )
#anova(m15.glmer, m10.glmer, test = "Chi") # not sig (p=.57)

# add age*audiencen
#ftable(reallynze$age, reallynze$audiencen)		
m16.glmer <- update(m10.glmer, .~.+age*audiencen)
#m16.glm <- update(m10.glm, .~.+age*audiencen)
#vif(m16.glm) # VIFs too high (max = )
#anova(m16.glmer, m10.glmer, test = "Chi") # not sig (p=.17)

# add age*priming
#ftable(reallynze$age, reallynze$priming)		
m17.glmer <- update(m10.glmer, .~.+age*priming)
#m17.glm <- update(m10.glm, .~.+age*priming)
#vif(m17.glm) # VIFs too high (max = )
#anova(m17.glmer, m10.glmer, test = "Chi") # not sig (p=.21)

# add gender*fun
#ftable(reallynze$gender, reallynze$fun)		
m18.glmer <- update(m10.glmer, .~.+gender*fun)
#m18.glm <- update(m10.glm, .~.+gender*fun)
#vif(m18.glm) # VIFs too high (max = )
#anova(m18.glmer, m10.glmer, test = "Chi") # not sig (p=.94)

# add gender*emo
#ftable(reallynze$gender, reallynze$emo)		
m19.glmer <- update(m10.glmer, .~.+gender*emo)
#m19.glm <- update(m10.glm, .~.+gender*emo)
#vif(m19.glm) # VIFs too high (max = )
#anova(m19.glmer, m10.glmer, test = "Chi") # not sig (p=.86)

# add gender*grad
#ftable(reallynze$gender, reallynze$grad)		
m20.glmer <- update(m10.glmer, .~.+gender*grad)
#m20.glm <- update(m10.glm, .~.+gender*grad)
#vif(m20.glm) # VIFs too high (max = )
#anova(m20.glmer, m10.glmer, test = "Chi") # not sig (p=.74)

# add gender*sem
#ftable(reallynze$gender, reallynze$sem)		
m21.glmer <- update(m10.glmer, .~.+gender*sem)
#m21.glm <- update(m10.glm, .~.+gender*sem)
#vif(m21.glm) # VIFs too high (max = )
#anova(m21.glmer, m10.glmer, test = "Chi") # not sig (p=.066) = marg sig (AIC+1, BIC+32, logLik -8)

# add gender*smsx
#ftable(reallynze$gender, reallynze$smsx)		
m22.glmer <- update(m10.glmer, .~.+gender*smsx)
#m22.glm <- update(m10.glm, .~.+gender*smsx)
#vif(m22.glm) # VIFs too high (max = )
#anova(m22.glmer, m10.glmer, test = "Chi") # not sig (p=.33)

# add gender*audiencen
#ftable(reallynze$gender, reallynze$audiencen)		
m23.glmer <- update(m10.glmer, .~.+gender*audiencen)
#m23.glm <- update(m10.glm, .~.+gender*audiencen)
#vif(m23.glm) # VIFs too high (max = )
#anova(m23.glmer, m10.glmer, test = "Chi") # not sig (p=.1)

# add gender*priming
#ftable(reallynze$gender, reallynze$priming)		
m24.glmer <- update(m10.glmer, .~.+gender*priming)
#m24.glm <- update(m10.glm, .~.+gender*priming)
#vif(m24.glm) # VIFs too high (max = )
#anova(m24.glmer, m10.glmer, test = "Chi") # not sig (p=.12)

# add fun*emo
#ftable(reallynze$fun, reallynze$emo)			
m25.glmer <- update(m10.glmer, .~.+fun*emo)
#m25.glm <- update(m10.glm, .~.+fun*emo)
#vif(m25.glm) # VIFs too high (max = )
#anova(m25.glmer, m10.glmer, test = "Chi") # not sig (p=.83)

# add fun*grad
#ftable(reallynze$fun, reallynze$grad)			
m26.glmer <- update(m10.glmer, .~.+fun*grad)
#m26.glm <- update(m10.glm, .~.+fun*grad)
#vif(m26.glm) # VIFs too high (max = )
#anova(m26.glmer, m10.glmer, test = "Chi") # not sig (p=.12)

# add fun*sem
#ftable(reallynze$fun, reallynze$sem)			
m27.glmer <- update(m10.glmer, .~.+fun*sem)
#m27.glm <- update(m10.glm, .~.+fun*sem)
#vif(m27.glm) # VIFs too high (max = )
#anova(m27.glmer, m10.glmer, test = "Chi") # not sig (p=.79)

# add fun*smsx
#ftable(reallynze$fun, reallynze$smsx)			 
m28.glmer <- update(m10.glmer, .~.+fun*smsx)
#m28.glm <- update(m10.glm, .~.+fun*smsx)
#vif(m28.glm) # VIFs too high (max = )
#anova(m28.glmer, m10.glmer, test = "Chi") # not sig (p=.22)

# add fun*audiencen
#ftable(reallynze$fun, reallynze$audiencen)
m29.glmer <- update(m10.glmer, .~.+fun*audiencen)
#m29.glm <- update(m10.glm, .~.+fun*audiencen)
#vif(m29.glm) # VIFs too high (max = )
#anova(m29.glmer, m10.glmer, test = "Chi") # not sig (p=.091) = mar sig  (AIC-1, BIC+11, logLik -7)

# add emo*grad
#ftable(reallynze$emo, reallynze$grad)
m30.glmer <- update(m10.glmer, .~.+emo*grad)
#m30.glm <- update(m10.glm, .~.+emo*grad)
#vif(m30.glm) # VIFs too high (max = )
#anova(m30.glmer, m10.glmer, test = "Chi") # not sig (p=.9)

# add emo*sem
#ftable(reallynze$emo, reallynze$sem)
m31.glmer <- update(m10.glmer, .~.+emo*sem)
#m31.glm <- update(m10.glm, .~.+emo*sem)
#vif(m31.glm) # VIFs too high (max = )
#anova(m31.glmer, m10.glmer, test = "Chi") # not sig (p=.42)

# add emo*smsx
#ftable(reallynze$emo, reallynze$smsx)
m32.glmer <- update(m10.glmer, .~.+emo*smsx)
#m32.glm <- update(m10.glm, .~.+emo*smsx)
#vif(m32.glm) # VIFs too high (max = )
#anova(m32.glmer, m10.glmer, test = "Chi") # not sig (p=.75)

# add emo*audiencen
#ftable(reallynze$emo, reallynze$audiencen)
m33.glmer <- update(m10.glmer, .~.+emo*audiencen)
#m33.glm <- update(m10.glm, .~.+emo*audiencen)
#vif(m33.glm) # VIFs too high (max = )
#anova(m33.glmer, m10.glmer, test = "Chi") # not sig (p=.15)

# add emo*priming
#ftable(reallynze$emo, reallynze$priming)
m34.glmer <- update(m10.glmer, .~.+emo*priming)
#m34.glm <- update(m10.glm, .~.+emo*priming)
#vif(m34.glm) # VIFs too high (max = )
#anova(m34.glmer, m10.glmer, test = "Chi") # not sig (p=.025)

# add grad*sem: WARNING: not üpossible!
#ftable(reallynze$grad, reallynze$sem) 
#m35.glmer <- update(m10.glmer, .~.+grad*sem)
#m35.glm <- update(m10.glm, .~.+grad*sem)
#vif(m35.glm) # VIFs too high (max = )
#anova(m35.glmer, m10.glmer, test = "Chi")

# add grad*smsx
#ftable(reallynze$grad, reallynze$smsx)
m36.glmer <- update(m10.glmer, .~.+grad*smsx)
#m36.glm <- update(m10.glm, .~.+grad*smsx)
#vif(m36.glm) # VIFs too high (max = )
#anova(m36.glmer, m10.glmer, test = "Chi") # not sig (p=.56)

# add grad*audiencen
#ftable(reallynze$grad, reallynze$audiencen)
m37.glmer <- update(m10.glmer, .~.+grad*audiencen)
#m37.glm <- update(m10.glm, .~.+grad*audiencen)
#vif(m37.glm) # VIFs too high (max = )
#anova(m37.glmer, m10.glmer, test = "Chi") # not sig (p=.063) = mar sig  (AIC-1, BIC+19, logLik -7)

# add grad*priming
#ftable(reallynze$grad, reallynze$priming)
m38.glmer <- update(m10.glmer, .~.+grad*priming)
#m38.glm <- update(m10.glm, .~.+grad*priming)
#vif(m38.glm) # VIFs too high (max = )
#anova(m38.glmer, m10.glmer, test = "Chi") # not sig (p=.11)

# add sem*smsx
#ftable(reallynze$sem, reallynze$smsx)
m39.glmer <- update(m10.glmer, .~.+sem*smsx)
#m39.glm <- update(m10.glm, .~.+sem*smsx)
#vif(m39.glm) # VIFs too high (max = )
#anova(m39.glmer, m10.glmer, test = "Chi") # not sig (p=.76)

# add sem*audiencen
#ftable(reallynze$sem, reallynze$audiencen)
m40.glmer <- update(m10.glmer, .~.+sem*audiencen)
#m40.glm <- update(m10.glm, .~.+sem*audiencen)
#vif(m40.glm) # VIFs too high (max = )
#anova(m40.glmer, m10.glmer, test = "Chi") # not sig (p=.57)

# add sem*priming
#ftable(reallynze$sem, reallynze$priming)
m41.glmer <- update(m10.glmer, .~.+sem*priming)
#m41.glm <- update(m10.glm, .~.+sem*priming)
#vif(m41.glm) # VIFs too high (max = )
#anova(m41.glmer, m10.glmer, test = "Chi") # not sig (p=.96)

# add smsx*audiencen
ftable(reallynze$smsx, reallynze$audiencen)
m42.glmer <- update(m10.glmer, .~.+smsx*audiencen)
m42.glm <- update(m10.glm, .~.+smsx*audiencen)
vif(m42.glm) # VIFs okay (max = 4.022 < 5, cf. Zuur et al. 2009)
anova(m42.glmer, m10.glmer, test = "Chi") # SIG (p=.022) BUT AIC-4, BIC +8(!), logLik -5

meblrm_m42 <- meblrm.summary(m0.glm, m42.glm, m0.glmer, m42.glmer, reallynze$really) #
meblrm_m42 # does not indicate sig! 

library(car)
meblrm_m42_Anova <- Anova(m42.glmer, type = "III", test = "Chi")
meblrm_m42_Anova # does show that smsx*audiencen is NOT SIGNIFICANT

# add smsx*priming
#ftable(reallynze$smsx, reallynze$priming)
m43.glmer <- update(m10.glmer, .~.+smsx*priming)
#m43.glm <- update(m10.glm, .~.+smsx*priming)
#vif(m43.glm) # VIFs ok.
#anova(m43.glmer, m10.glmer, test = "Chi") # not sig (p=.41)

###########################################################################
###              3-way interactions
###########################################################################
# find all 3-way interactions
#install.packages("utils")
library(utils)
vars <- c("age", "gender", "fun", "emo", "grad", "sem", "smsx", "audiencen", "priming")
intac <- t(combn(vars, 3))
intac

# check if interaction is possible
agegenderfun <- min(ftable(reallynze$age, reallynze$gender, reallynze$fun, reallynze$really)) 
agegenderemo <- min(ftable(reallynze$age, reallynze$gender, reallynze$emo, reallynze$really)) 
agegendergrad <- min(ftable(reallynze$age, reallynze$gender, reallynze$grad, reallynze$really)) 
agegendersem <- min(ftable(reallynze$age, reallynze$gender, reallynze$sem, reallynze$really)) 
agegendersmsx <- min(ftable(reallynze$age, reallynze$gender, reallynze$smsx, reallynze$really)) 
agegenderaudiencen <- min(ftable(reallynze$age, reallynze$gender, reallynze$audiencen, reallynze$really)) 
agegenderpriming <- min(ftable(reallynze$age, reallynze$gender, reallynze$priming, reallynze$really)) 
agefunemo <- min(ftable(reallynze$age, reallynze$fun, reallynze$emo, reallynze$really)) 
agefungrad <- min(ftable(reallynze$age, reallynze$fun, reallynze$grad, reallynze$really)) 
agefunsem <- min(ftable(reallynze$age, reallynze$fun, reallynze$sem, reallynze$really)) 
agefunsmsx <- min(ftable(reallynze$age, reallynze$fun, reallynze$smsx, reallynze$really)) 
agefunaudiencen <- min(ftable(reallynze$age, reallynze$fun, reallynze$audiencen, reallynze$really)) 
agefunpriming <- min(ftable(reallynze$age, reallynze$fun, reallynze$priming, reallynze$really)) 
ageemograd <- min(ftable(reallynze$age, reallynze$emo, reallynze$grad, reallynze$really)) 
ageemosem <- min(ftable(reallynze$age, reallynze$emo, reallynze$sem, reallynze$really)) 
ageemosmsx <- min(ftable(reallynze$age, reallynze$emo, reallynze$smsx, reallynze$really)) 
ageemoaudiencen <- min(ftable(reallynze$age, reallynze$emo, reallynze$audiencen, reallynze$really)) 
ageemopriming <- min(ftable(reallynze$age, reallynze$emo, reallynze$priming, reallynze$really)) 
agegradsem <- min(ftable(reallynze$age, reallynze$grad, reallynze$sem, reallynze$really)) 
agegradsmsx <- min(ftable(reallynze$age, reallynze$grad, reallynze$smsx, reallynze$really)) 
agegradaudiencen <- min(ftable(reallynze$age, reallynze$grad, reallynze$audiencen, reallynze$really)) 
agegradpriming <- min(ftable(reallynze$age, reallynze$grad, reallynze$priming, reallynze$really)) 
agesemsmsx <- min(ftable(reallynze$age, reallynze$sem, reallynze$smsx, reallynze$really)) 
agesemaudiencen <- min(ftable(reallynze$age, reallynze$sem, reallynze$audiencen, reallynze$really)) 
agesempriming <- min(ftable(reallynze$age, reallynze$sem, reallynze$priming, reallynze$really)) 
agesmsxaudiencen <- min(ftable(reallynze$age, reallynze$smsx, reallynze$audiencen, reallynze$really)) 
agesmsxpriming <- min(ftable(reallynze$age, reallynze$smsx, reallynze$priming, reallynze$really)) 
ageaudiencenpriming <- min(ftable(reallynze$age, reallynze$audiencen, reallynze$priming, reallynze$really)) 
genderfunemo <- min(ftable(reallynze$gender, reallynze$fun, reallynze$emo, reallynze$really)) 
genderfungrad <- min(ftable(reallynze$gender, reallynze$fun, reallynze$grad, reallynze$really)) 
genderfunsem <- min(ftable(reallynze$gender, reallynze$fun, reallynze$sem, reallynze$really)) 
genderfunsmsx <- min(ftable(reallynze$gender, reallynze$fun, reallynze$smsx, reallynze$really)) 
genderfunaudiencen <- min(ftable(reallynze$gender, reallynze$fun, reallynze$audiencen, reallynze$really)) 
genderfunpriming <- min(ftable(reallynze$gender, reallynze$fun, reallynze$priming, reallynze$really)) 
genderemograd <- min(ftable(reallynze$gender, reallynze$emo, reallynze$grad, reallynze$really)) 
genderemosem <- min(ftable(reallynze$gender, reallynze$emo, reallynze$sem, reallynze$really)) 
genderemosmsx <- min(ftable(reallynze$gender, reallynze$emo, reallynze$smsx, reallynze$really)) 
genderemoaudiencen <- min(ftable(reallynze$gender, reallynze$emo, reallynze$audiencen, reallynze$really)) 
genderemopriming <- min(ftable(reallynze$gender, reallynze$emo, reallynze$priming, reallynze$really)) 
gendergradsem <- min(ftable(reallynze$gender, reallynze$grad, reallynze$sem, reallynze$really)) 
gendergradsmsx <- min(ftable(reallynze$gender, reallynze$grad, reallynze$smsx, reallynze$really)) 
gendergradaudiencen <- min(ftable(reallynze$gender, reallynze$grad, reallynze$audiencen, reallynze$really)) 
gendergradpriming <- min(ftable(reallynze$gender, reallynze$grad, reallynze$priming, reallynze$really)) 
gendersemsmsx <- min(ftable(reallynze$gender, reallynze$sem, reallynze$smsx, reallynze$really)) 
gendersemaudiencen <- min(ftable(reallynze$gender, reallynze$sem, reallynze$audiencen, reallynze$really)) 
gendersempriming <- min(ftable(reallynze$gender, reallynze$sem, reallynze$priming, reallynze$really)) 
gendersmsxaudiencen <- min(ftable(reallynze$gender, reallynze$smsx, reallynze$audiencen, reallynze$really)) 
gendersmsxpriming <- min(ftable(reallynze$gender, reallynze$smsx, reallynze$priming, reallynze$really)) 
genderaudiencenpriming <- min(ftable(reallynze$gender, reallynze$audiencen, reallynze$priming, reallynze$really)) 
funemograd <- min(ftable(reallynze$fun, reallynze$emo, reallynze$grad, reallynze$really)) 
funemosem <- min(ftable(reallynze$fun, reallynze$emo, reallynze$sem, reallynze$really)) 
funemosmsx <- min(ftable(reallynze$fun, reallynze$emo, reallynze$smsx, reallynze$really)) 
funemoaudiencen <- min(ftable(reallynze$fun, reallynze$emo, reallynze$audiencen, reallynze$really)) 
funemopriming <- min(ftable(reallynze$fun, reallynze$emo, reallynze$priming, reallynze$really)) 
fungradsem <- min(ftable(reallynze$fun, reallynze$grad, reallynze$sem, reallynze$really)) 
fungradsmsx <- min(ftable(reallynze$fun, reallynze$grad, reallynze$smsx, reallynze$really)) 
fungradaudiencen <- min(ftable(reallynze$fun, reallynze$grad, reallynze$audiencen, reallynze$really)) 
fungradpriming <- min(ftable(reallynze$fun, reallynze$grad, reallynze$priming, reallynze$really)) 
funsemsmsx <- min(ftable(reallynze$fun, reallynze$sem, reallynze$smsx, reallynze$really)) 
funsemaudiencen <- min(ftable(reallynze$fun, reallynze$sem, reallynze$audiencen, reallynze$really)) 
funsempriming <- min(ftable(reallynze$fun, reallynze$sem, reallynze$priming, reallynze$really)) 
funsmsxaudiencen <- min(ftable(reallynze$fun, reallynze$smsx, reallynze$audiencen, reallynze$really)) 
funsmsxpriming <- min(ftable(reallynze$fun, reallynze$smsx, reallynze$priming, reallynze$really)) 
funaudiencenpriming <- min(ftable(reallynze$fun, reallynze$audiencen, reallynze$priming, reallynze$really)) 
emogradsem <- min(ftable(reallynze$emo, reallynze$grad, reallynze$sem, reallynze$really)) 
emogradsmsx <- min(ftable(reallynze$emo, reallynze$grad, reallynze$smsx, reallynze$really)) 
emogradaudiencen <- min(ftable(reallynze$emo, reallynze$grad, reallynze$audiencen, reallynze$really)) 
emogradpriming <- min(ftable(reallynze$emo, reallynze$grad, reallynze$priming, reallynze$really)) 
emosemsmsx <- min(ftable(reallynze$emo, reallynze$sem, reallynze$smsx, reallynze$really)) 
emosemaudiencen <- min(ftable(reallynze$emo, reallynze$sem, reallynze$audiencen, reallynze$really)) 
emosempriming <- min(ftable(reallynze$emo, reallynze$sem, reallynze$priming, reallynze$really)) 
emosmsxaudiencen <- min(ftable(reallynze$emo, reallynze$smsx, reallynze$audiencen, reallynze$really)) 
emosmsxpriming <- min(ftable(reallynze$emo, reallynze$smsx, reallynze$priming, reallynze$really)) 
emoaudiencenpriming <- min(ftable(reallynze$emo, reallynze$audiencen, reallynze$priming, reallynze$really)) 
gradsemsmsx <- min(ftable(reallynze$grad, reallynze$sem, reallynze$smsx, reallynze$really)) 
gradsemaudiencen <- min(ftable(reallynze$grad, reallynze$sem, reallynze$audiencen, reallynze$really)) 
gradsempriming <- min(ftable(reallynze$grad, reallynze$sem, reallynze$priming, reallynze$really)) 
gradsmsxaudiencen <- min(ftable(reallynze$grad, reallynze$smsx, reallynze$audiencen, reallynze$really)) 
gradsmsxpriming <- min(ftable(reallynze$grad, reallynze$smsx, reallynze$priming, reallynze$really)) 
gradaudiencenpriming <- min(ftable(reallynze$grad, reallynze$audiencen, reallynze$priming, reallynze$really)) 
semsmsxaudiencen <- min(ftable(reallynze$sem, reallynze$smsx, reallynze$audiencen, reallynze$really)) 
semsmsxpriming <- min(ftable(reallynze$sem, reallynze$smsx, reallynze$priming, reallynze$really)) 
semaudiencenpriming <- min(ftable(reallynze$sem, reallynze$audiencen, reallynze$priming, reallynze$really)) 
smsxaudiencenpriming <- min(ftable(reallynze$smsx, reallynze$audiencen, reallynze$priming, reallynze$really)) 
# test which interactions are possible
testpos3intact <- c(agegenderfun, agegenderemo, agegendergrad, agegendersem, agegendersmsx, 
                    agegenderaudiencen, agegenderpriming, agefunemo, agefungrad, agefunsem, 
                    agefunsmsx, agefunaudiencen, agefunpriming, ageemograd, ageemosem, ageemosmsx, 
                    ageemoaudiencen, ageemopriming, agegradsem, agegradsmsx, agegradaudiencen, 
                    agegradpriming, agesemsmsx, agesemaudiencen, agesempriming, agesmsxaudiencen, 
                    agesmsxpriming, ageaudiencenpriming, genderfunemo, genderfungrad, genderfunsem, 
                    genderfunsmsx, genderfunaudiencen, genderfunpriming, genderemograd, genderemosem, 
                    genderemosmsx, genderemoaudiencen, genderemopriming, gendergradsem, gendergradsmsx, 
                    gendergradaudiencen, gendergradpriming, gendersemsmsx, gendersemaudiencen, 
                    gendersempriming, gendersmsxaudiencen, gendersmsxpriming, genderaudiencenpriming, 
                    funemograd, funemosem, funemosmsx, funemoaudiencen, funemopriming, fungradsem, 
                    fungradsmsx, fungradaudiencen, fungradpriming, funsemsmsx, funsemaudiencen, 
                    funsempriming, funsmsxaudiencen, funsmsxpriming, funaudiencenpriming, emogradsem, 
                    emogradsmsx, emogradaudiencen, emogradpriming, emosemsmsx, emosemaudiencen, 
                    emosempriming, emosmsxaudiencen, emosmsxpriming, emoaudiencenpriming, gradsemsmsx, 
                    gradsemaudiencen, gradsempriming, gradsmsxaudiencen, gradsmsxpriming, 
                    gradaudiencenpriming, semsmsxaudiencen, semsmsxpriming, semaudiencenpriming)
names(testpos3intact) <- c("agegenderfun", "agegenderemo", "agegendergrad", "agegendersem", "agegendersmsx", 
                           "agegenderaudiencen", "agegenderpriming", "agefunemo", "agefungrad", "agefunsem", 
                           "agefunsmsx", "agefunaudiencen", "agefunpriming", "ageemograd", "ageemosem", "ageemosmsx", 
                           "ageemoaudiencen", "ageemopriming", "agegradsem", "agegradsmsx", "agegradaudiencen", 
                           "agegradpriming", "agesemsmsx", "agesemaudiencen", "agesempriming", "agesmsxaudiencen", 
                           "agesmsxpriming", "ageaudiencenpriming", "genderfunemo", "genderfungrad", "genderfunsem", 
                           "genderfunsmsx", "genderfunaudiencen", "genderfunpriming", "genderemograd", "genderemosem", 
                           "genderemosmsx", "genderemoaudiencen", "genderemopriming", "gendergradsem", "gendergradsmsx", 
                           "gendergradaudiencen", "gendergradpriming", "gendersemsmsx", "gendersemaudiencen", 
                           "gendersempriming", "gendersmsxaudiencen", "gendersmsxpriming", "genderaudiencenpriming", 
                           "funemograd", "funemosem", "funemosmsx", "funemoaudiencen", "funemopriming", "fungradsem", 
                           "fungradsmsx", "fungradaudiencen", "fungradpriming", "funsemsmsx", "funsemaudiencen", 
                           "funsempriming", "funsmsxaudiencen", "funsmsxpriming", "funaudiencenpriming", "emogradsem", 
                           "emogradsmsx", "emogradaudiencen", "emogradpriming", "emosemsmsx", "emosemaudiencen", 
                           "emosempriming", "emosmsxaudiencen", "emosmsxpriming", "emoaudiencenpriming", "gradsemsmsx", 
                           "gradsemaudiencen", "gradsempriming", "gradsmsxaudiencen", "gradsmsxpriming", 
                           "gradaudiencenpriming", "semsmsxaudiencen", "semsmsxpriming", "semaudiencenpriming")
tstintact3 <- names(testpos3intact)[which(testpos3intact >= 1)]
tstintact3; length(tstintact3)

# add age*fun*emo
#ftable(reallynze$age, reallynze$fun, reallynze$emo)
m100.glmer <- update(m10.glmer, .~.+age*fun*emo)
#m100.glm <- update(m10.glm, .~.+age*fun*emo)
#vif(m100.glm) # vifs too high max = 
#anova(m100.glmer, m10.glmer, test = "Chi") # not sig (p=.68)

# add age*fun*smsx
#ftable(reallynze$age, reallynze$fun, reallynze$smsx)
m101.glmer <- update(m10.glmer, .~.+age*fun*smsx)
#m101.glm <- update(m10.glm, .~.+age*fun*smsx)
#vif(m101.glm) # vifs too high max = 
#anova(m101.glmer, m10.glmer, test = "Chi") # not sig (p=.34)

# add gender*fun*emo
#ftable(reallynze$gender, reallynze$fun, reallynze$emo)
m102.glmer <- update(m10.glmer, .~.+gender*fun*emo)
#m102.glm <- update(m10.glm, .~.+gender*fun*emo)
#vif(m102.glm) # vifs too high max = 
#anova(m102.glmer, m10.glmer, test = "Chi") # not sig (p=.62)

# add gender*fun*grad
#ftable(reallynze$gender, reallynze$fun, reallynze$grad)
m103.glmer <- update(m10.glmer, .~.+gender*fun*grad)
#m103.glm <- update(m10.glm, .~.+gender*fun*grad)
#vif(m103.glm) # vifs too high max = 
#anova(m103.glmer, m10.glmer, test = "Chi") # not sig (p=.41)

# add gender*fun*smsx
#ftable(reallynze$gender, reallynze$fun, reallynze$smsx)
m104.glmer <- update(m10.glmer, .~.+gender*fun*smsx)
#m104.glm <- update(m10.glm, .~.+gender*fun*smsx)
#vif(m104.glm) # vifs too high max = 
#anova(m104.glmer, m10.glmer, test = "Chi") # not sig (p=.55)

# add gender*fun*priming
#ftable(reallynze$gender, reallynze$fun, reallynze$priming)
m105.glmer <- update(m10.glmer, .~.+gender*fun*priming)
#m105.glm <- update(m10.glm, .~.+gender*fun*priming)
#vif(m105.glm) # vifs too high max = 
#anova(m105.glmer, m10.glmer, test = "Chi") # not sig (p=.66)

# add gender*emo*smsx
#ftable(reallynze$gender, reallynze$emo, reallynze$smsx)
m106.glmer <- update(m10.glmer, .~.+gender*emo*smsx)
#m106.glm <- update(m10.glm, .~.+gender*emo*smsx)
#vif(m106.glm) # vifs too high max = 
#anova(m106.glmer, m10.glmer, test = "Chi") # not sig (p=.47)

# add gender*emo*priming
#ftable(reallynze$gender, reallynze$emo, reallynze$priming)
m107.glmer <- update(m10.glmer, .~.+gender*emo*priming)
#m107.glm <- update(m10.glm, .~.+gender*emo*priming)
#vif(m107.glm) # vifs too high max = 
#anova(m107.glmer, m10.glmer, test = "Chi") # not sig (p=.21)

# add gender*grad*smsx
#ftable(reallynze$gender, reallynze$grad, reallynze$smsx)
m108.glmer <- update(m10.glmer, .~.+gender*grad*smsx)
#m108.glm <- update(m10.glm, .~.+gender*grad*smsx)
#vif(m108.glm) # vifs too high max = 
#anova(m108.glmer, m10.glmer, test = "Chi") # not sig (p=.25)

# add gender*grad*priming
#ftable(reallynze$gender, reallynze$grad, reallynze$priming)
m109.glmer <- update(m10.glmer, .~.+gender*grad*priming)
#m109.glm <- update(m10.glm, .~.+gender*grad*priming)
#vif(m109.glm) # vifs too high max = 
#anova(m109.glmer, m10.glmer, test = "Chi") # not sig (p=.22)

# add gender*sem*smsx
#ftable(reallynze$gender, reallynze$sem, reallynze$smsx)
m110.glmer <- update(m10.glmer, .~.+gender*sem*smsx)
#m110.glm <- update(m10.glm, .~.+gender*sem*smsx)
#vif(m110.glm) # vifs too high max = 
#anova(m110.glmer, m10.glmer, test = "Chi") # not sig (p=.12)

# add gender*smsx*priming
#ftable(reallynze$gender, reallynze$smsx, reallynze$priming)
m111.glmer <- update(m10.glmer, .~.+gender*smsx*priming)
#m111.glm <- update(m10.glm, .~.+gender*smsx*priming)
#vif(m111.glm) # vifs ok
#anova(m111.glmer, m10.glmer, test = "Chi") # not sig (p=.11)

# add fun*emo*smsx
#ftable(reallynze$fun, reallynze$emo, reallynze$smsx)
m112.glmer <- update(m10.glmer, .~.+fun*emo*smsx)
#m112.glm <- update(m10.glm, .~.+fun*emo*smsx)
#vif(m112.glm) # vifs too high max = 
#anova(m112.glmer, m10.glmer, test = "Chi") # not sig (p=.69)

# add fun*emo*audiencen
#ftable(reallynze$fun, reallynze$emo, reallynze$audiencen)
m113.glmer <- update(m10.glmer, .~.+fun*emo*audiencen)
#m113.glm <- update(m10.glm, .~.+fun*emo*audiencen)
#vif(m113.glm) # vifs too high max = 
#anova(m113.glmer, m10.glmer, test = "Chi") # not sig (p=.23)

# add fun*emo*priming
#ftable(reallynze$fun, reallynze$emo, reallynze$priming)
m114.glmer <- update(m10.glmer, .~.+fun*emo*priming)
#m114.glm <- update(m10.glm, .~.+fun*emo*priming)
#vif(m114.glm) # vifs too high max = 
#anova(m114.glmer, m10.glmer, test = "Chi") # not sig (p=.6)

# add fun*grad*priming
#ftable(reallynze$fun, reallynze$grad, reallynze$priming)
m115.glmer <- update(m10.glmer, .~.+fun*grad*priming)
#m115.glm <- update(m10.glm, .~.+fun*grad*priming)
#vif(m115.glm) # vifs too high max = 
#anova(m115.glmer, m10.glmer, test = "Chi") # not sig (p=.1) = mar sig  (AIC+4, BIC+42, logLik -8)

# add fun*sem*priming
#ftable(reallynze$fun, reallynze$sem, reallynze$priming)
m116.glmer <- update(m10.glmer, .~.+fun*sem*priming)
#m116.glm <- update(m10.glm, .~.+fun*sem*priming)
#vif(m116.glm) # vifs too high max = 
#anova(m116.glmer, m10.glmer, test = "Chi") # not sig (p=.99)

# add fun*smsx*priming
#ftable(reallynze$fun, reallynze$smsx, reallynze$priming)
m117.glmer <- update(m10.glmer, .~.+fun*smsx*priming)
#m117.glm <- update(m10.glm, .~.+fun*smsx*priming)
#vif(m117.glm) # vifs too high max = 
#anova(m117.glmer, m10.glmer, test = "Chi") # not sig (p=.56)

# add fun*audiencen*priming
#ftable(reallynze$fun, reallynze$audiencen, reallynze$priming)
m118.glmer <- update(m10.glmer, .~.+fun*audiencen*priming)
#m118.glm <- update(m10.glm, .~.+fun*audiencen*priming)
#vif(m118.glm) # vifs too high max = 
#anova(m118.glmer, m10.glmer, test = "Chi") # not sig (p=.12)

# add emo*sem*smsx
#ftable(reallynze$emo, reallynze$sem, reallynze$smsx)
m119.glmer <- update(m10.glmer, .~.+emo*sem*smsx)
#m119.glm <- update(m10.glm, .~.+emo*sem*smsx)
#vif(m119.glm) # vifs too high max = 
#anova(m119.glmer, m10.glmer, test = "Chi") # not sig (p=.89)

# add emo*smsx*priming
#ftable(reallynze$emo, reallynze$smsx, reallynze$priming)
m120.glmer <- update(m10.glmer, .~.+emo*smsx*priming)
#m120.glm <- update(m10.glm, .~.+emo*smsx*priming)
#vif(m120.glm) # vifs too high max = 
#anova(m120.glmer, m10.glmer, test = "Chi") # not sig (p=.47)

# add emo*audiencen*priming
#ftable(reallynze$emo, reallynze$audiencen, reallynze$priming)
m121.glmer <- update(m10.glmer, .~.+emo*audiencen*priming)
#m121.glm <- update(m10.glm, .~.+emo*audiencen*priming)
#vif(m121.glm) # vifs too high max = 
#anova(m121.glmer, m10.glmer, test = "Chi") # not sig (p=.29)

# add grad*smsx*priming
#ftable(reallynze$grad, reallynze$smsx, reallynze$priming)
m122.glmer <- update(m10.glmer, .~.+grad*smsx*priming)
#m122.glm <- update(m10.glm, .~.+grad*smsx*priming)
#vif(m122.glm) # vifs too high max = 
#anova(m122.glmer, m10.glmer, test = "Chi") # not sig (p=.43)

#########################################
# set up summary table
meblrm_mapnze <- meblrm.summary(m0.glm, m10.glm, m0.glmer, m10.glmer, reallynze$really) #
meblrm_mapnze

# save results to disc
write.table(meblrm_mapnze, "meblrm_mapnze.txt", sep="\t")

library(car)
meblrm_mapnze_Anova <- Anova(m10.glmer, type = "III", test = "Chi")
meblrm_mapnze_Anova

# save results to disc
write.table(meblrm_mapnze_Anova, "meblrm_mapnze_Anova.txt", sep="\t")

effectage <- anova(m0.glmer, m1.glmer, test = "Chi")

effectpriming <- anova(m1.glmer, m9.glmer, test = "Chi")

effectagegender <- anova(m9.glmer, m10.glmer, test = "Chi")

# use customized model comparireallyn function
# create comparireallyns
m1.m0 <- anova(m1.glmer, m0.glmer, test = "Chi")
m2.m1 <- anova(m2.glmer, m1.glmer, test = "Chi")
m3.m1 <- anova(m3.glmer, m1.glmer, test = "Chi")
m4.m1 <- anova(m4.glmer, m1.glmer, test = "Chi")
m5.m1 <- anova(m5.glmer, m1.glmer, test = "Chi")
m6.m1 <- anova(m6.glmer, m1.glmer, test = "Chi")
m7.m1 <- anova(m7.glmer, m1.glmer, test = "Chi")
m8.m1 <- anova(m8.glmer, m1.glmer, test = "Chi")
m9.m1 <- anova(m9.glmer, m1.glmer, test = "Chi")
m10.m9 <- anova(m10.glmer, m9.glmer, test = "Chi")
m11.m10 <- anova(m11.glmer, m10.glmer, test = "Chi")
m12.m10 <- anova(m12.glmer, m10.glmer, test = "Chi")
m13.m10 <- anova(m13.glmer, m10.glmer, test = "Chi")
m14.m10 <- anova(m14.glmer, m10.glmer, test = "Chi")
m15.m10 <- anova(m15.glmer, m10.glmer, test = "Chi")
m16.m10 <- anova(m16.glmer, m10.glmer, test = "Chi")
m17.m10 <- anova(m17.glmer, m10.glmer, test = "Chi")
m18.m10 <- anova(m18.glmer, m10.glmer, test = "Chi")
m19.m10 <- anova(m19.glmer, m10.glmer, test = "Chi")
m20.m10 <- anova(m20.glmer, m10.glmer, test = "Chi")
m21.m10 <- anova(m21.glmer, m10.glmer, test = "Chi")
m22.m10 <- anova(m22.glmer, m10.glmer, test = "Chi")
m23.m10 <- anova(m23.glmer, m10.glmer, test = "Chi")
m24.m10 <- anova(m24.glmer, m10.glmer, test = "Chi")
m25.m10 <- anova(m25.glmer, m10.glmer, test = "Chi")
m26.m10 <- anova(m26.glmer, m10.glmer, test = "Chi")
m27.m10 <- anova(m27.glmer, m10.glmer, test = "Chi")
m28.m10 <- anova(m28.glmer, m10.glmer, test = "Chi")
m29.m10 <- anova(m29.glmer, m10.glmer, test = "Chi")
m30.m10 <- anova(m30.glmer, m10.glmer, test = "Chi")
m31.m10 <- anova(m31.glmer, m10.glmer, test = "Chi")
m32.m10 <- anova(m32.glmer, m10.glmer, test = "Chi")
m33.m10 <- anova(m33.glmer, m10.glmer, test = "Chi")
m34.m10 <- anova(m30.glmer, m10.glmer, test = "Chi")
m36.m10 <- anova(m36.glmer, m10.glmer, test = "Chi")
m37.m10 <- anova(m37.glmer, m10.glmer, test = "Chi")
m38.m10 <- anova(m38.glmer, m10.glmer, test = "Chi")
m39.m10 <- anova(m39.glmer, m10.glmer, test = "Chi")
m40.m10 <- anova(m40.glmer, m10.glmer, test = "Chi")
m41.m10 <- anova(m41.glmer, m10.glmer, test = "Chi")
m10.m10 <- anova(m10.glmer, m10.glmer, test = "Chi")
m43.m10 <- anova(m43.glmer, m10.glmer, test = "Chi")
m101.m10 <- anova(m101.glmer, m10.glmer, test = "Chi")
m102.m10 <- anova(m102.glmer, m10.glmer, test = "Chi")
m103.m10 <- anova(m103.glmer, m10.glmer, test = "Chi")
m104.m10 <- anova(m104.glmer, m10.glmer, test = "Chi")
m105.m10 <- anova(m105.glmer, m10.glmer, test = "Chi")
m106.m10 <- anova(m106.glmer, m10.glmer, test = "Chi")
m107.m10 <- anova(m107.glmer, m10.glmer, test = "Chi")
m108.m10 <- anova(m108.glmer, m10.glmer, test = "Chi")
m109.m10 <- anova(m109.glmer, m10.glmer, test = "Chi")
m110.m10 <- anova(m110.glmer, m10.glmer, test = "Chi")
m111.m10 <- anova(m111.glmer, m10.glmer, test = "Chi")
m112.m10 <- anova(m112.glmer, m10.glmer, test = "Chi")
m113.m10 <- anova(m113.glmer, m10.glmer, test = "Chi")
m114.m10 <- anova(m114.glmer, m10.glmer, test = "Chi")
m115.m10 <- anova(m115.glmer, m10.glmer, test = "Chi")
m116.m10 <- anova(m116.glmer, m10.glmer, test = "Chi")
m117.m10 <- anova(m117.glmer, m10.glmer, test = "Chi")
m118.m10 <- anova(m118.glmer, m10.glmer, test = "Chi")
m119.m10 <- anova(m119.glmer, m10.glmer, test = "Chi")
m120.m10 <- anova(m120.glmer, m10.glmer, test = "Chi")
m121.m10 <- anova(m121.glmer, m10.glmer, test = "Chi")
m122.m10 <- anova(m122.glmer, m10.glmer, test = "Chi")
# create a list of the model comparireallyns
mdlcmp <- list(m1.m0, m2.m1, m3.m1, m4.m1, m5.m1, m6.m1, m7.m1, m8.m1, m9.m1, m10.m9, 
               m11.m10, m12.m10, m13.m10, m14.m10, m15.m10, m16.m10, m17.m10, m18.m10, 
               m19.m10, m20.m10, m21.m10, m22.m10, m23.m10, m24.m10, m25.m10, m26.m10, 
               m27.m10, m28.m10, m29.m10, m30.m10, m31.m10, m32.m10, m33.m10, m34.m10, 
               m36.m10, m37.m10, m38.m10, m39.m10, m40.m10, m41.m10, m10.m10, m43.m10, 
               m101.m10, m102.m10, m103.m10, m104.m10, m105.m10, m106.m10, m107.m10, 
               m108.m10, m109.m10, m110.m10, m111.m10, m112.m10, m113.m10, m114.m10, 
               m115.m10, m116.m10, m117.m10, m118.m10, m119.m10, m120.m10, m121.m10, 
               m122.m10)
# apply function
mdl.cmp.glmersc.swsu.dm <- mdl.fttng.swsu(mdlcmp)
# inspect output
mdl.cmp.glmersc.swsu.dm

write.table(mdl.cmp.glmersc.swsu.dm, "mdl_cmp_glmersc_swsu_dm_20181025.txt", sep="\t")
###########################################################
# Post-hoc analysis
library (multcomp)
summary(glht(m42.glmer, mcp(grad="Tukey")))

################################################################
#                 IMPORTANT OBJECTS
################################################################
# inspect very important objects
head(reallynze)

# glmer
effectage

effectpriming

effectagegender

meblrm_mapnze

meblrm_mapnze_Anova

###############################################################
###              END PART 4
###############################################################


