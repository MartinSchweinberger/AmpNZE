##################################################################
# Titel:      The Amplifier System of New Zealand English - Part 1
# R version:  3.4.1 (2017-06-30) -- "Single Candle"
# Autor:      Martin Schweinberger
# Date:       2018-11-06
# Contact:    martin.schweinberger.hh@gmail.com
##################################################################
# Disclaimer: If you have questions,suggestions or you found errors
#             or in case you would to provide feedback, questions
#             write an email to martin.schweinberger.hh@gmail.com.
# Citation:   If you use this script or results thereof, please cite it as:
#             Schweinberger, Martin. 2018. "The Amplifier System of New Zealand English, Part 1",
#             unpublished R script, The University of Queensland.
###############################################################
#                   START
###############################################################
# remove all lists from the current workspace
rm(list=ls(all=T))
# set wd
setwd("D:\\Uni\\Projekte\\02-Intensification\\AmpNZE")
# load packages
library(car)
library(gsubfn)
library(Hmisc)
library(plyr)
library(reshape)
library(rms)
library(stringr)
library(tm)
source("D:\\R/POStagObject.R") # for pos-tagging objects in R
###############################################################
# set options
options(stringsAsFactors = F)
options(scipen = 999)
# define image directory
imageDirectory<-"images"
# specify path to corpra
corpusnzepath <- "D:\\Uni\\Korpora\\Original\\ICE New Zealand\\Spoken"
bionzepath <- "D:\\Uni\\Korpora\\Metadata/BiodataIceNewZealand.txt"
# define corpus files
corpus.files = list.files(path = corpusnzepath, pattern = NULL, all.files = T,
                          full.names = T, recursive = T, ignore.case = T, include.dirs = T)
###############################################################
# load and start processing corpus
corpusnze <- sapply(corpus.files, function(x) {
  x <- scan(x, what = "char", sep = "", quote = "", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  x <- str_replace_all(x, fixed("\n"), " ")
  x <- paste(x, sep = " ", collapse = " ")
  x <- strsplit(gsub("(<I>)", "~\\1", x), "~" )
  x <- unlist(x)
  x <- x[2:length(x)]
  x <- gsub("<I> ", "", x)
  x <- as.vector(unlist(x))
} )
sf <- as.vector(unlist(sapply(corpusnze, function(x){
  x <- length(x)
  x <- str_replace_all(x, "2","1 2")
  x <- str_replace_all(x, "3","1 2 3")
  x <- str_replace_all(x, "4","1 2 3 4")
  x <- strsplit(x, " ")
})))
corpusnze01 <- unlist(corpusnze)
fl <- as.vector(unlist(sapply(corpusnze01, function(x){
  x <- gsub("\\#.*", "", x)
  x <- gsub(".*:", "", x)
})))
suspk <- lapply(corpusnze01, function(x){
  x <- strsplit(gsub("(<ICE-NZ:[A-Z][0-9][A-Z])", "~\\1", x), "~" )
  x <- unlist(x)
  x <- x[2:length(x)]
} )
l <- as.vector(unlist(sapply(suspk, function(x){
  x <- length(x)
})))
fln <- rep(fl, l)
sfn <- rep(sf, l)
suspkn <- as.vector(unlist(suspk))
spk <- as.vector(unlist(sapply(suspkn, function(x){
  x <- gsub(">.*", "", x)
  x <- gsub(".*\\:", "", x)
})))
# create a df out of the vectors
df <- data.frame(fln, sfn, spk, suspkn)
# create a full-id vector
flid <- as.vector(unlist(apply(df, 1, FUN=function(x){
  x <- paste("<", x[1], ":", x[2], "$", x[3], ">", sep = "", collapse = "")
} )))
suspknn <- sub(" ", " <#>", suspkn)
suspknn <- gsub("<#><#>", "<#>", suspknn)
# speech unit
su <- sapply(suspknn, function(x){
  x <- strsplit(gsub("(<#>)", "~\\1", x), "~" )
})
sun <- sapply(su, function(x){
  x <- x[2:length(x)]
})
sunn <- as.vector(unlist(sun))
# extract the number of speech units per turn
n <- as.vector(unlist(sapply(sun, function(x){
  x <- length(x)
})))
# create a clean vector of the speech units
sucl <- sunn
sucl <- gsub("<#>", " ", sucl)
sucl <- gsub("</{0,1}\\?>", " ", sucl)
sucl <- gsub("<&> {0,1}.* {0,1}</&>", " ", sucl)
sucl <- gsub("<O> {0,1}.* {0,1}</O>", " ", sucl)
sucl <- gsub("<unclear> {0,1}.* {0,1}</unclear>", " ", sucl)
sucl <- gsub("<.{0,1}/{0,1}.{0,1}[A-Z]{0,1}[a-z]{1,}>", " ", sucl)
sucl <- gsub("<,{1,3}>", " ", sucl)
sucl <- gsub("<\\.>[a-z]{1,}</\\.>", " ", sucl)
sucl <- gsub("</{0,1}\\{[0-9]{0,2}>", " ", sucl)
sucl <- gsub("</{0,1}[0-9]{0,1}\\[{0,1}\\{{0,1}[0-9]{0,2}>", " ", sucl)
sucl <- gsub("<\\[/[0-9]{0,2}>", " ", sucl)
sucl <- gsub("</{0,1}\\[[0-9]{0,2}>", " ", sucl)
sucl <- gsub("</{0,1}\\}[0-9]{0,2}>", " ", sucl)
sucl <- gsub("</{0,1}\\][0-9]{0,2}>", " ", sucl)
sucl <- gsub("</{0,1}\\.[0-9]{0,2}>", " ", sucl)
sucl <- gsub("</{0,1}[A-Z]{0,2}[0-9]{0,2}>", " ", sucl)
sucl <- gsub("</{0,1}[a-z]{1,}=[A-Z]{0,1}[a-z]{1,}>", " ", sucl)
sucl <- gsub(" {2,}", " ", sucl)
sucl <- str_trim(sucl, side = "both")
# test if non-words are still present
#t1 <- sucl[grep("\\[", sucl)]
# additional cleaning
sucl <- gsub("<&> laughter <&>", " ", sucl)
sucl <- gsub("<", " ", sucl)
sucl <- gsub("&", " ", sucl)
sucl <- gsub(">", " ", sucl)
sucl <- gsub("/", "", sucl)
sucl <- gsub("[", " ", sucl, fisucled = T)
sucl <- gsub("{", " ", sucl, fisucled = T)
sucl <- gsub("?", "", sucl, fisucled = T)
sucl <- gsub(" {2,}", " ", sucl)
sucl <- str_trim(sucl, side = "both")
# rep vectors number of times of the speech units per turn
flnn <- rep(fln, n)
sfnn <- rep(sfn, n)
spkn <- rep(spk, n)
flidn <- rep(flid, n)
corpusnzedf <- data.frame(flidn, flnn, sfnn, spkn, sunn, sucl)
#t1 <- corpusnzedf$sucl[grep("\\[", corpusnzedf$sucl)] #test
###############################################################
# save raw data to disc
write.table(corpusnzedf, "corpusnzeraw.txt", sep = "\t", row.names = F, col.names = T)
corpusnzedf <- read.table("corpusnzeraw.txt", sep = "\t", header = T)
###############################################################
# remove empty speech units
corpusnzedf <- corpusnzedf[corpusnzedf$sucl != "",]
# split data into smaller chunks
pos01 <- corpusnzedf$sucl[1:5000]
pos02 <- corpusnzedf$sucl[5001:10000]
pos03 <- corpusnzedf$sucl[10001:15000]
pos04 <- corpusnzedf$sucl[15001:20000]
pos05 <- corpusnzedf$sucl[20001:25000]
pos06 <- corpusnzedf$sucl[25001:30000]
pos07 <- corpusnzedf$sucl[30001:35000]
pos08 <- corpusnzedf$sucl[35001:40000]
pos09 <- corpusnzedf$sucl[40001:45000]
pos10 <- corpusnzedf$sucl[45001:50000]
pos11 <- corpusnzedf$sucl[50001:55000]
pos12 <- corpusnzedf$sucl[55001:60000]
pos13 <- corpusnzedf$sucl[60001:65000]
pos14 <- corpusnzedf$sucl[65001:nrow(corpusnzedf)]
# reload libraries
source("D:\\R/POStagObject.R") # for pos-tagging objects in R
library(NLP)
library(openNLP)
library(openNLPmodels.en)
# pos tagging data
nzepos01 <- POStag(object = pos01)
nzepos01 <- as.vector(unlist(nzepos01))
writeLines(nzepos01, con = "nzepos01.txt", sep = "\n", useBytes = FALSE)
# chunk 2
nzepos02 <- POStag(object = pos02)
nzepos02 <- as.vector(unlist(nzepos02))
writeLines(nzepos02, con = "nzepos02.txt", sep = "\n", useBytes = FALSE)
# chunk 03
nzepos03 <- POStag(object = pos03)
nzepos03 <- as.vector(unlist(nzepos03))
writeLines(nzepos03, con = "nzepos03.txt", sep = "\n", useBytes = FALSE)
# chunk 04
nzepos04 <- POStag(object = pos04)
nzepos04 <- as.vector(unlist(nzepos04))
writeLines(nzepos04, con = "nzepos04.txt", sep = "\n", useBytes = FALSE)
# chunk 05
nzepos05 <- POStag(object = pos05)
nzepos05 <- as.vector(unlist(nzepos05))
writeLines(nzepos05, con = "nzepos05.txt", sep = "\n", useBytes = FALSE)
# chunk 06
nzepos06 <- POStag(object = pos06)
nzepos06 <- as.vector(unlist(nzepos06))
writeLines(nzepos06, con = "nzepos06.txt", sep = "\n", useBytes = FALSE)
# chunk 07
nzepos07 <- POStag(object = pos07)
nzepos07 <- as.vector(unlist(nzepos07))
writeLines(nzepos07, con = "nzepos07.txt", sep = "\n", useBytes = FALSE)
# chunk 08
nzepos08 <- POStag(object = pos08)
nzepos08 <- as.vector(unlist(nzepos08))
writeLines(nzepos08, con = "nzepos08.txt", sep = "\n", useBytes = FALSE)
# chunk 09
nzepos09 <- POStag(object = pos09)
nzepos09 <- as.vector(unlist(nzepos09))
writeLines(nzepos09, con = "nzepos09.txt", sep = "\n", useBytes = FALSE)
# chunk 10
nzepos10 <- POStag(object = pos10)
nzepos10 <- as.vector(unlist(nzepos10))
writeLines(nzepos10, con = "nzepos10.txt", sep = "\n", useBytes = FALSE)
# chunk 11
nzepos11 <- POStag(object = pos11)
nzepos11 <- as.vector(unlist(nzepos11))
writeLines(nzepos11, con = "nzepos11.txt", sep = "\n", useBytes = FALSE)
# chunk 12
nzepos12 <- POStag(object = pos12)
nzepos12 <- as.vector(unlist(nzepos12))
writeLines(nzepos11, con = "nzepos12.txt", sep = "\n", useBytes = FALSE)
# chunk 13
nzepos13 <- POStag(object = pos13)
nzepos13 <- as.vector(unlist(nzepos13))
writeLines(nzepos13, con = "nzepos13.txt", sep = "\n", useBytes = FALSE)
# chunk 14
nzepos14 <- POStag(object = pos14)
nzepos14 <- as.vector(unlist(nzepos14))
writeLines(nzepos14, con = "nzepos14.txt", sep = "\n", useBytes = FALSE)
# list pos tagged elements
postag.files = c("nzepos01.txt", "nzepos02.txt", "nzepos03.txt",  "nzepos04.txt", "nzepos05.txt",
                 "nzepos06.txt",  "nzepos07.txt", "nzepos08.txt", "nzepos09.txt",  "nzepos10.txt",
                 "nzepos11.txt", "nzepos12.txt", "nzepos13.txt", "nzepos14.txt")
# load pos tagged elements
nzepos <- sapply(postag.files, function(x) {
  x <- scan(x, what = "char", sep = "\n", quote = "", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  x <- str_replace_all(x, fixed("\n"), " ")
})
# unlist pos tagged elements
corpusnzedf$nzepos <- unlist(nzepos)
###############################################################
# extract number of adjs per line
pstggd <- corpusnzedf$nzepos
lpstggd <- strsplit(pstggd, " ")
nlpstggd <- as.vector(unlist(sapply(lpstggd, function(x){
  x <- x[grep("JJ", x)]
  x <- length(x) } )))
rp <- nlpstggd
rp <- ifelse(rp == 0, 1, rp)
# load function for concordancing
source("D:\\R/ConcR_2.3_loadedfiles.R")
# set parameters for concordancing
pattern <- "[A-Z]{0,1}[a-z]{1,}\\/JJ[A-Z]{0,1}"
context <- 50
# extract all adjectives (concordance)
concjjnze <- ConcR(corpusnzedf$nzepos, pattern, context, all.pre = FALSE)
# repeat rows in data farem as often as there are adjectives in it (if 0 adj, repeat once)
corpusnzeadjdf <- corpusnzedf[rep(seq(nrow(corpusnzedf)), rp),]
# combine data sets
corpusnzeadj <- data.frame(1:nrow(corpusnzeadjdf), corpusnzeadjdf, concjjnze)
# remove rows without Tokens
ampnze <- corpusnzeadj[is.na(corpusnzeadj$Token) == F,]
# add clean column names
colnames(ampnze) <- c("ID", "FileSpeaker", "File", "Subfile", "Speaker", "SpeechUnit", "CleanSpeechUnit",
                      "PosTaggedSpeechUnit", "OriginalString", "PreContext", "Adjective", "PostContext")
# clean adjectives
ampnze$Adjective <- str_replace_all(ampnze$Adjective, fixed("/JJ"), "")
# add Vraiant column
ampnze$Variant <- gsub(".* ", "", str_trim(ampnze$PreContext, side = "both")) 
# inspect data
#nrow(ampnze); head(ampnze)

###############################################################
# define amplifiers
amplifiers <- c("absolutely", "actually", "aggressively", 
                "amazingly", "appallingly", "awful", "awfully", 
                "badly", "bloody", "certainly", "clearly",
                "complete", "dead", "completely", "considerably", 
                "crazy", "decidedly", "definitely",  "distinctly", 
                "dreadfully", "enormously", "entirely", "especially", 
                "exactly", "exceedingly", "exceptionally", 
                "excruciatingly", "extraordinarily", "extremely",
                "fiercely", "firmly", "frightfully", "fucking", 
                "fully", "genuinely", "greatly",
                "grossly", "heavily", "highly", "hopelessly", 
                "horrendously", "hugely",
                "immediately", "immensely", "incredibly", 
                "infinitely", "intensely", "irrevocably",
                "mad", "mega", "mighty", "most", "much", 
                "obviously", "openly", "overwhelmingly", "particularly", 
                "perfectly", "plenty", "positively", "precisely", 
                "pretty", "profoundly", "purely", 
                #"quite", 
                "real", "really", "remarkably", "seriously", 
                "shocking",   "significant", "significantly", "so", 
                "specially", "specifically", "strikingly",
                "strongly", "substantially", "super", "surely", 
                "terribly", "terrifically", 
                #"too",
                "total", "totally", "traditionally", "true", 
                "truly", "ultra", "utterly", "very",
                "viciously", 
                #"well", 
                "wholly", "wicked", "wildly")
# clean ice nze data
ampnze$Function <- str_trim(ampnze$PostContext, side = "both")
ampnze$Function <- tolower(ampnze$Function)
ampnze$Function <- gsub(" {2,}", " ", ampnze$Function)
ampnze$Function <- gsub(".* ", "", ampnze$Function)
ampnze$Function <- gsub(".*/n.*", "Attributive", ampnze$Function)
ampnze$Function <- ifelse(ampnze$Function == "", "Predicative", ampnze$Function)
ampnze$Function <- ifelse(ampnze$Function == "Attributive" | ampnze$Function == "Predicative", ampnze$Function, "remove")
# remove items for which Function could not be clearly determined
ampnze <- ampnze[ampnze$Function != "remove",]
# register
ampnze$Genre <- gsub("_.*", "", ampnze$File)
ampnze$Genre <- gsub("-.*", "", ampnze$Genre)
ampnze$Genre <- ifelse(ampnze$Genre == "S1A", "PrivateDialogue", ampnze$Genre)
ampnze$Genre <- ifelse(ampnze$Genre == "S1B", "PublicDialogue", ampnze$Genre)
ampnze$Genre <- ifelse(ampnze$Genre == "S2A", "UnscriptedMonologue", ampnze$Genre)
ampnze$Genre <- ifelse(ampnze$Genre == "S2B", "ScriptedMonologue", ampnze$Genre)
ampnze$Genre <- ifelse(ampnze$Genre == "PrivateDialogue" | ampnze$Genre == "PublicDialogue" |  
                         ampnze$Genre == "UnscriptedMonologue" | ampnze$Genre == "ScriptedMonologue", ampnze$Genre, "remove")
ampnze <- ampnze[ampnze$Genre != "remove",]
# shorten post Context
ampnze$PostContext <- substr(ampnze$PostContext, 1, ifelse((nchar(ampnze$PostContext)+25) <25, maampnze(nchar(ampnze$PostContext)), 25))
# pre Context
ampnze$PreContext <- str_trim(ampnze$PreContext, side = "both")
ampnze$PreContextLong <- ampnze$PreContext
ampnze$PreContextLong <- substr(ampnze$PreContextLong, ifelse(nchar(ampnze$PreContextLong)-25 <=0, 1, 
                                                              nchar(ampnze$PreContextLong)-25), nchar(ampnze$PreContextLong))
ampnze$PreContext <- gsub(".* ", "", ampnze$PreContext)
# amplifier variant
ampnze$PreContext <- gsub("\\/.*", "", ampnze$PreContext)
ampnze$Variant <- ifelse(ampnze$PreContext %in% amplifiers, ampnze$PreContext, "0")
# amplified y/n
ampnze$Amplified <- ifelse(ampnze$Variant == "0", 0, 1) 
# adjective
ampnze$Adjective <- tolower(ampnze$Adjective)
# inspect data
nrow(ampnze); head(ampnze); table(ampnze$Variant)

# define forms that require removal
sups <- c(".*most.*", ".*more.*") 
negs <- c(".*not.*", ".*never.*", ".*n't.*")
downtoners <- c(".*sort/.*", ".*kind/.*", ".* bit/.*", ".*somewhat.*", ".*fairly.*", 
                ".*rather.*", ".*reasonably.*", ".*slightly.*", ".*comparatively.*", ".*semi.*", 
                ".*relatively.*", ".*little.*", ".*somehow.*", ".*almost.*", ".*partly.*", 
                ".*hardly.*", ".* less.*", ".*barely.*", ".* just/.*")
specialforms <- c(".* too.*", ".*quite.*")
PostContextdowntoners <- c(".*enough.*")
nonpropadj <- c("only", "much", "many", "cheaper", "cheaperr", "bests", "larger")
# check length of dataset
str(ampnze); head(ampnze); nrow(ampnze)#; table(ampnze$pint); head(ampnze$PreContextLong); head(ampnze$PreContextLong)

# find items to be removed
supsidx <- unique(grep(paste(sups,collapse="|"), ampnze$PreContextLong, value=F))
negsidx <- unique(grep(paste(negs,collapse="|"), ampnze$PreContextLong, value=F))
downtonersidx <- unique(grep(paste(downtoners,collapse="|"), ampnze$PreContextLong, value=F))
specialformsidx <- unique(grep(paste(specialforms,collapse="|"), ampnze$PreContextLong, value=F))
PostContextdowntonersidx <- unique(grep(paste(PostContextdowntoners,collapse="|"), ampnze$PostContext, value=F))
nonpropadjidx <- unique(grep(paste(nonpropadj,collapse="|"), ampnze$Adjective, value=F))
# combine indices
idxs <- unique(c(supsidx, negsidx, downtonersidx, specialformsidx, PostContextdowntonersidx, nonpropadjidx))
# remove forms that require removal
ampnze <- ampnze[-idxs,]
# remove empty values
ampnze <- ampnze[!ampnze$Variant == "", ]
###############################################################
# save raw data to disc
write.table(ampnze, "ampnze02_wo_neg.txt", sep = "\t", row.names = F)
###############################################################
# code priming
prim1 <- c(rep(0, 1), ampnze$Variant[1:length(ampnze$Variant)-1])
prim2 <- c(rep(0, 2), ampnze$Variant[1:(length(ampnze$Variant)-2)])
prim3 <- c(rep(0, 3), ampnze$Variant[1:(length(ampnze$Variant)-3)])
primtb <- cbind(ampnze$Variant, prim1, prim2, prim3)

ampnze$Priming <- as.vector(unlist(apply(primtb, 1, function(x){
  x <- ifelse(x[1]== "0" , "noprime",
              ifelse(x[1] == x[2] | x[1] == x[3] | x[1] == x[4], "prime", "noprime"))
})))
# remove items that were not intensified by a minimum of 2 intensifier variants
nrow(ampnze)

pintadjtb <- table(ampnze$Adjective, ampnze$Variant)
#pintadjtb <- pintadjtb[2:nrow(pintadjtb),]
pintadjtb <- pintadjtb[,2:ncol(pintadjtb)]
pintadjtb2 <- apply(pintadjtb, 1, function(x){
  x <- ifelse(x > 1, 1, x)})
pintadjtb3 <- colSums(pintadjtb2)
pintadjschildes <- names(pintadjtb3)[which(pintadjtb3 >=2 )]
ampnze <- ampnze[ampnze$Adjective %in% pintadjschildes, ]
nrow(ampnze)

# inspect adjectives
names(table(ampnze$Adjective))

# clean adjectives
ampnze$Adjective <- ifelse(ampnze$Adjective == "cheaperr", "cheaper",
                           ifelse(ampnze$Adjective == "largerr", "larger", ampnze$Adjective))
# create vector with false adjectives
rmvadj <- c("er", "okay")
ampnze$remove <- ifelse(ampnze$Adjective %in% rmvadj, "remove", ampnze$Adjective)
ampnze <- ampnze[ampnze$remove != "remove",]
ampnze$remove <- NULL
# inspecta data
nrow(ampnze); length(table(ampnze$Adjective)); head(ampnze)

###############################################################
# save raw data to disc
write.table(ampnze, "ampnze03_semiclean.txt", sep = "\t", row.names = F)
###############################################################
#                        END PART 1
###############################################################
