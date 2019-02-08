##################################################################
# Titel:      The Amplifier System of New Zealand English - Part 2
# R version:  3.4.1 (2017-06-30) -- "Single Candle"
# Autor:      Martin Schweinberger
# Date:       2018-11-06
# Contact:    martin.schweinberger.hh@gmail.com
##################################################################
# Disclaimer: If you have questions,suggestions or you found errors
#             or in case you would to provide feedback, questions
#             write an email to martin.schweinberger.hh@gmail.com.
# Citation:   If you use this script or results thereof, please cite it as:
#             Schweinberger, Martin. 2018. "The Amplifier System of New Zealand English, Part 2",
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
library(syuzhet)
library(tm)
###############################################################
# Setting options
options(stringsAsFactors = F)
options(scipen = 999)
# define image dausctors
imageDnzectory<-"images"
# specify path to biodata
bionzepath <- "D:\\Uni\\Korpora\\Metadata/BiodataIceNewZealand.txt"
# load corpus data
ampnze <- read.table("ampnze03_semiclean.txt", sep = "\t", header = T)
###############################################################
# add biodata
###############################################################
# read in data
bio <- read.table(bionzepath, sep = "\t", header=TRUE)
# clean file for joining
bio$File <- gsub("([A|B])", "\\1-", bio$text.id)
bio$Subfile <- bio$subfile.id
bio$Speaker <- bio$spk.ref
################################################################
# add conversation type (same vs mixed) and audience size (number of interlocutors)
# add info on samesex vs diffsex files
# create a "file plus speaker" variable
bio$txtspk <- paste(bio$text.id, "$", bio$spk.ref, sep = "")
sstb <- table(bio$txtspk, bio$sex)
fl <- gsub("\\$.*", "", rownames(sstb))
spk <- gsub(".*\\$([A-Z]{0,1}\\?{0,1}).*", "\\1", rownames(sstb))
fml <- as.vector(unlist(sapply(sstb[,1], function(x){
  x <- ifelse(x >0, 1, 0) } )))
ml <- as.vector(unlist(sapply(sstb[,2], function(x){
  x <- ifelse(x >0, 1, 0) } )))
# set up data frame
ssdf <- data.frame(rownames(sstb), fl, spk, fml, ml)
colnames(ssdf) <- c("full", "file", "spk", "female", "male")
# clean speakers
ssdf$spk <- gsub("?", "", ssdf$spk, fixed = T)
ssdf <- ssdf[ssdf$spk != "", ]
# inspect data
#head(ssdf)

ssdf <- ssdf[rowSums(ssdf[, c(4:5)]) > 0, ]
ssdf$sex <- as.vector(unlist(sapply(ssdf[, 5], function(x){
  x <- ifelse(x == 1, "male", "female") } )))
# inspect data
#head(ssdf)

tst <- table(ssdf$file, ssdf$sex)
intloc <- apply(tst, 1, function(x){
  x <- rowSums(tst) } )
intloc <- intloc[,1]
tst2 <- t(apply(tst, 1, function(x){
  x <- ifelse(x > 0, 1, 0) } ))
tst2 <- as.data.frame(tst2)
rwsm <- as.vector(unlist(rowSums(tst2)))
tst2$ConversationType <- as.vector(unlist(sapply(rwsm, function(x){
  x <- ifelse(x == 1, "samesex", "mixedsex") } )))
tst2$text.id <- rownames(tst2)
tst2$ints <- intloc
tst3 <- data.frame(tst2$ConversationType, tst2$text.id, tst2$ints)
colnames(tst3) <- gsub("tst2.", "", colnames(tst3))
# chnage colnames so that they match the bio colnames
colnames(tst3) <- c("ConversationType", "text.id", "AudienceSize") 
# inspect data
#head(tst3)

# combine number of interlocutors and emoire data frame
bio <- join(bio, tst3, by = c("text.id"), type = "left")
# inspect data
#str(bio); head(bio)

################################################################
# homogenize column names
bio$Subfile <- as.factor(bio$Subfile)
# join data
ampnze <- join(ampnze, bio, by = c("File", "Subfile", "Speaker"), type = "left")
# remove superfluous columns
ampnze$id <- NULL
ampnze$id.orig <- NULL
ampnze$file.speaker.id <- NULL
ampnze$text.id <- NULL
ampnze$subfile.id <- NULL
ampnze$spk.ref <- NULL
# rename columns
colnames(ampnze) <- ifelse(colnames(ampnze) == "age", "Age", 
                           ifelse(colnames(ampnze) == "sex", "Gender",
                                  ifelse(colnames(ampnze) == "ethnicity", "Ethnicity",
                                         ifelse(colnames(ampnze) == "occupation", "Occupation", colnames(ampnze)))))
# inspect data
head(ampnze); str(ampnze); colnames(ampnze)

################################################################
# remove speakers below the age of 18
ampnze <- ampnze[ampnze$Age != "0-18", ]
################################################################
# recode ethnicity
Maori <- c("Cook Island Maori",  "Maori",  "Maori/French",  "Maori/Nuiean/Samoan",  "Maori/Samoan",  "Samoan",  "Tokelauan")
NAN <- c("Chinese",  "European Asian",  "Lebanese",  "Maori/Negro/Pakeha",
         "Maori/Pakeha",  "Other",  "Pakeha/Maori",  "Pakeha/Samoan",  "Samoan/Dutch",
         "Samoan/Pakeha",  "Scots/Maori",  "Semitic")
Pakeha <- c("Cook Island/Pakeha",  "Dutch",  "Greek",  "Irish",  "Jewish",
            "NZ Greek",  "Pakeha",  "Pakeha/Arabic",  "Pakeha/Asian",  "Pakeha/Danish",
            "Pakeha/Dutch",  "Pakeha/German",  "Spanish/German")
ampnze$Ethnicity <- ifelse(ampnze$Ethnicity %in% Maori, "Maori",
                           ifelse(ampnze$Ethnicity %in% NAN, NA,
                                  ifelse(ampnze$Ethnicity %in% Pakeha, "Pakeha", ampnze$Ethnicity)))
################################################################
# recode occupation 1
ampnze$Occupation <- sapply(ampnze$Occupation, function(x) {
  x <- gsub(".*[T|t]eache.*", "acmp", x)     # teachers
  x <- gsub(".*[D|d]nzect.*", "acmp", x)     # dnzectors
  x <- gsub(".*[L|l]awye.*", "acmp", x)      # lawyers
  x <- gsub(".*[B|b]arrist.*", "acmp", x)    # barristers
  x <- gsub(".*[M|m]anage.*", "acmp", x)     # managers
  x <- gsub(".*[C|c]ler[i]{0,1}[c|k].*", "acmp", x)  # clerks
  x <- gsub(".*[A|a]ccounta.*", "acmp", x)   # accountants
  x <- gsub(".*[A|a]naly.*", "acmp", x)      # analysts
  x <- gsub(".*[M|m]iniste.*", "acmp", x)    # ministers
  x <- gsub(".*[A|a]dministr.*", "acmp", x)  # administartors
  x <- gsub(".*[J|j]udge.*", "acmp", x)      # judges
  x <- gsub(".*[A|a]pprent.*", "sml", x)     # apprentices
  x <- gsub(".*[J|j]ournal.*", "acmp", x)    # journalists
  x <- gsub(".*[W|w]rite.*", "acmp", x)      # writers
  x <- gsub(".*[S|s]olicit.*", "acmp", x)    # solicitors
  x <- gsub(".*[P|p]arliam.*", "acmp", x)    # pariamentary speakers/members of parliament
  x <- gsub(".*[P|p]rincip.*", "acmp", x)    # principals
  x <- gsub(".*[R|r]esearch.*", "acmp", x)   # researchers/research assistants
  x <- gsub(".*[C|c]omment.*", "acmp", x)    # political commentators
  x <- gsub(".*MP.*", "acmp", x)             # members of parliament
  x <- gsub(".*[L|l]ectur.*", "acmp", x)     # lecturers
  x <- gsub(".*[P|p]rogramm.*", "acmp", x)   # programmers & programm producers/dnzectors
  x <- gsub(".*[E|e]dit.*", "acmp", x)       # editors
  x <- gsub(".*[B|b]ank.*", "acmp", x)       # bankers
  x <- gsub(".*[P|p]rof.*", "acmp", x)       # professors
}  )
acmp <- c("Actor", "Anglican Youth Worker", "Auctioneer", "Book Shop Assistant",
          "Broadcaster", "Choreographer", "Company Chairman", "Composer/Music Education Advisor",
          "Constable", "Consultant", "Cultural Consultant", "Dbase Admin/Checkout Operator",
          "Delivery Contractor","Director Education Centre",  "Documentary Producer", "Economist", "Education Officer",
          "Entertainer", "ESL Tutor", "Executive Radio Producer", "Faculty", "Knowledge Engineer",
          "Leader of the Opposition", "Librarian", "Library Assistant", "Marketing Rep",
          "Musician", "Occupational Therapist", "P/T Counsellor", "P/T Japanese Tutor",
          "P/T Telemarketing Trainer/Student", "P/T Tutor", "Police Constable", "PR Office Assistant",
          "Priest", "Privacy Commissioner", "Processing Officer", "Radio Broadcaster",
          "Reservations Consultant", "Restauranteur", "Self-Employed Consultant", "Self Employed Musician",
          "Senior Policy Advisor", "Shop Owner", "Student", "Student Support Person",
          "Student/Checkout Operator", "Student/Tutor", "Talkback Host", "Teaching Assistant",
          "Tutor", "Tutor/Fulltime Student", "Tutor/Student", "TV Quiz Show Host",
          "Uni Student", "Advisor to Race Relations", "Board Secretary", "Broadcaster/Financial Consultant",
          "Cartage Contractor", "Catholic Priest", "Company Executive", "Computer Consultant",
          "Dean of Commerce", "Deputy Academic Registrar", "Deputy Secretary to Treasury",
          "Distributor/Student", "Education Offiver/Tupperware", "Environmental Geochemist",
          "Freelance Broadcast", "FREELANCE Broadcaster", "Governor General of NZ",
          "Home Economist", "Launchmaster", "Liaison Officer", "Media", "Musician/Broadcaster",
          "News Reader", "Newscaster", "Newsreader", "NZ High Commission London",
          "NZRFU Resource Coach", "Oral Historian", "P/T Cinema Attendant/Student",
          "Patent Attorney", "Post Doctoral Fellow", "Publisher", "Radio Announcer",
          "Radio Announcer/News Reader", "Radio News Reader", "Radio Newsreader",
          "Radio Producer/News Reader", "Radio Producer/Newsreader", "Radio/Television Presenter",
          "Radio/TV Presenter", "Receptionist/Telephonist", "Regional Councillor", "Reporter",
          "Retnzed Public Servant/Academic", "Scientist", "Scientist Enviroment Officer",
          "Secretary", "Seismologist", "Sports Broadcaster", "Sports Broadcaster RNZ",
          "Sports Education Consultant", "Student/Teaching Fellow", "Television News Presenter",
          "Television Newscaster", "Television Presenter", "Television Presenter/Radio Announcer",
          "Television Presenter/Reporter", "TV/Film Producer", "University Reader")
NAN <- c("Baker/Student", "Barmaid/Student", "Barman/Student", "Building Supervisor",
         "Chef/Student", "Cleaner/Student", "Clerk/Cleaner/Student", "Detective NZ Police",
         "ECE Worker/Student", "Factory Labourer/Student", "Gymnastics Coach/Student",
         "Jewellery Salesperson", "Kitchenhand/Student", "Machinist in Proof Centre/Student",
         "P/T Barman/Student", "P/T Cleaner/Student", "P/T Pool Attendant/Student",
         "P/T Tutor/KFC/Student", "P/T Tutor/Retail Assistant", "P/T Vet's Nurse/Student",
         "P/T Waitress/Student", "Police Officer", "Registered Nurse/Student",
         "Retail/Student", "Self-Employed", "Self Employed", "Shop Assistant/Student",
         "Steward/Student", "Student/P/T Waitress", "Telephonist/Student", "Temp",
         "Tutor/Shop Assistant", "Retnzed", "Unemployed", "Waiter/PR Officer", "Waitress/Student",
         "(no BI)")
sml <- c("Ambulance Officer", "Baker/Patissier", "Builder", "Car Salesman", "Caregiver",
         "Catering", "Cleaner", "Cook", "Courier", "Doing Catering Course", "Firefighter",
         "Florist", "Food Delivery", "Joiner", "Nanny", "P/T Cafe Worker", "P/T Cleaner",
         "P/T Film Technician", "P/T House Cleaner", "P/T Shop Assistant", "Primary Health Social Worker",
         "Property Developer", "Quantity Surv & Carpentry", "Receptionist",
         "Receptionist/Secretary", "Service Station Attendant", "Shop Assistant",
         "Sports Shop Assistant", "Video Hnze", "Waitress/Bartender", "Farm Tour Guide",
         "Fisheries Technician", "Former Airline Serviceman", "Kentucky Fried Chicken",
         "Technician")
# recode occupation 2
ampnze$Occupation <- ifelse(ampnze$Occupation %in% acmp == T, "acmp",
                            ifelse(ampnze$Occupation %in% NAN == T, NA,
                                   ifelse(ampnze$Occupation %in% sml == T, "sml", ampnze$Occupation))) 
# create columns with intensifier freqs
ampnze$very <- ifelse(ampnze$Variant == "very", 1, 0) 
ampnze$really <- ifelse(ampnze$Variant == "really", 1, 0) 
ampnze$so <- ifelse(ampnze$Variant == "so", 1, 0) 
ampnze$pretty <- ifelse(ampnze$Variant == "pretty", 1, 0) 
# clean data
ampnze <- ampnze[complete.cases(ampnze), ]
# inspect results
head(ampnze)

###############################################################
# recode age
ampnze$AgeOriginalClassification <- ampnze$Age
ampnze$Age <- ifelse(ampnze$Age == "16-19", "16-24",
                     ifelse(ampnze$Age == "20-24", "16-24",
                            ifelse(ampnze$Age == "25-29", "25-39",
                                   ifelse(ampnze$Age == "30-34", "25-39",
                                          ifelse(ampnze$Age == "35-39", "25-39",
                                                 ifelse(ampnze$Age == "40-44", "40-49",
                                                        ifelse(ampnze$Age == "45-49", "40-49",
                                                               ifelse(ampnze$Age == "50-54", "50+",
                                                                      ifelse(ampnze$Age == "55-59", "50+",
                                                                             ifelse(ampnze$Age == "60-64", "50+",
                                                                                    ifelse(ampnze$Age == "65-69", "50+",
                                                                                           ifelse(ampnze$Age == "70-74", "50+", ampnze$Age)))))))))))) 
###############################################################
# code freq of adj type by age group
frqadjtb <- table(ampnze$Age, ampnze$Adjective)
relfreqadjtb <- round(prop.table(frqadjtb, margin = 1)*100, 5)
relfreqadjdf <- as.data.frame(relfreqadjtb)
colnames(relfreqadjdf)[1:2] <- c("Age", "Adjective")
# add freq by date to data
ampnze <- merge(ampnze, relfreqadjdf, by=c("Age", "Adjective"))
# reorder data
ampnze <- ampnze[order(ampnze$ID),]
# inspect data
head(ampnze)

###############################################################
# code gradability
# gradability - manual classification
# done by trained grad students
# add gradability
ngrd_manual <- c("abject", "able", "abrasive", "abstract", "absurd", "abundant", "abusive",
                 "accurate", "acrimonious", "active", "advanced", "adverse", "affectionate", "afraid",
                 "aged", "aggressive", "agile", "agitated", "aimless", "airy", "alert", "alleged",
                 "allusive", "amazing", "ambitious", "amused", "amusing", "ancient", "angry",
                 "annoying", "anxious", "appalling", "apparent", "appealing", "applicable", "applied",
                 "appreciative", "apprehensive", "approachable", "appropriate", "approving",
                 "arduous", "arrogant", "ashamed", "associated", "astute", "athletic", "atrocious",
                 "attitudinal", "attractive", "authentic", "authoritarian   ", "authoritative",
                 "available", "aware", "awesome", "awful", "awkward", "awry", "bad", "bare", "base",
                 "battered", "beautiful", "beloved", "benevolent", "benign", "besetting", "bitter",
                 "bizarre", "bleak", "bleary", "bloody", "blotchy", "bold", "boppy", "bored",
                 "boring", "bossy", "brave", "brief", "bright", "brilliant", "broad", "browsing",
                 "brutal", "bubbly", "burly", "buzzy", "callous", "calm", "campy", "candid",
                 "capable", "careful", "careless", "casual", "cautious", "ceremonial", "challenging",
                 "changed", "charismatic", "charming", "cheap", "circumspect", "civic", "civil",
                 "civilised", "classy", "clever", "cocky", "cold", "collective", "colossal", "colourful",
                 "comfortable", "commandeered", "committed", "compatible", "compelling", "competent",
                 "competitive", "complex", "complicated", "conceivable", "concentrated", "concerned",
                 "confident", "confidential", "confused", "confusing", "considerable", "constructive",
                 "consultative", "contrived", "controversial", "convenient", "conventional", "converted",
                 "convinced", "cool", "corrupt", "cosy", "coy", "cramped", "crass", "crazy", "creative",
                 "criminal", "crippling", "critical", "cross", "crowded", "crucial", "cruel", "cumbersome",
                 "curious", "cushy", "cute", "cynical", "damaged", "damaging", "damp", "dangerous",
                 "daring", "darkened", "darn", "daunting", "dear", "debatable", "decent", "dedicated",
                 "deep", "defective", "defensive", "delicate", "delicious", "delighted", "delightful",
                 "dense", "dependent", "depressed", "desirable", "despairing", "desperate", "despicable",
                 "despondent", "destructive", "detailed", "detrimental", "devilish", "difficult",
                 "dirty", "disabled", "disadvantaged", "disappointed", "disappointing", "disastrous",
                 "disenchanted", "disgraceful", "disgusting", "dishonest", "disparaging", "distant",
                 "distinguished   ", "distorted", "distressed", "disturbed", "disturbing", "dizzy",
                 "dodgy", "dominant", "dotty", "double", "doubtful", "downhill", "dramatic", "dreadful",
                 "driving", "drunk", "drunken", "ductile", "dull", "dumb", "dusty", "dylan", "dynamic",
                 "dynamical", "eager", "early", "earnest", "earthy", "easterly", "eastern", "easy",
                 "eccentric", "economic", "edible", "effective", "efficient", "elderly", "elegant",
                 "eligible", "elitist", "elusive", "embarrassed", "embarrassing", "emergent", "eminent",
                 "emotional", "emotive", "encouraging", "energetic", "enlightening", "enormous",
                 "entertaining", "enthusiastic", "epic", "erudite", "estimated", "estranged", "everyday",
                 "evil", "exact", "exceptional", "excessive", "excited", "exciting", "expensive",
                 "experienced", "expert", "explicit", "express", "expressive", "extended", "extensive",
                 "extraordinary", "extravagant", "extroverted", "fabulous", "facile", "factual", "faint",
                 "familiar", "famous", "fanatic", "fancy", "fantastic", "fascinating", "fast", "fastidious",
                 "fat", "favourable", "favoured", "fearful", "feisty", "fergal", "ferocious", "fertile",
                 "fierce", "fiery", "filthy", "fine", "finished", "finite", "firm", "fitting", "fizzy",
                 "flexible", "fluffy", "fluttering", "foggy", "foolish", "forceful", "formalised",
                 "formidable", "fortunate", "frank", "frantic", "fraudulent", "fraught", "frenzied",
                 "frequent", "friendly", "frightening", "frightful", "frustrated", "frustrating",
                 "fulsome", "fun", "funny", "furious", "generous", "gentle", "giant", "gifted",
                 "gigantic", "glad", "glib", "glorious", "glossy", "good", "goodhearted", "gorgeous",
                 "gracious", "gradual", "grand", "grandiose", "grateful", "grave", "greasy", "great",
                 "grim", "groggy", "groovy", "gross", "grubby", "guilty", "gutless", "habitual",
                 "handsome", "handy", "hapless", "happy", "hard", "hardy", "harmful", "harmless",
                 "harmonic", "harsh", "hazardous", "hazy", "heavy", "hectic", "helpful", "hideous",
                 "high", "hilarious", "holy", "honest", "honorable", "honorary", "honourable",
                 "hooked", "hopeful", "hopeless", "horrendous", "horrible", "horrific", "hostile",
                 "hot", "huge", "humble", "humorous", "hungry", "hurt", "hysterical", "idealistic",
                 "igneous", "ignorant", "imaginative", "immature", "immediate", "immense", "imperative",
                 "important", "impotent", "impressive", "inane", "incompetent", "inconsistent",
                 "incorporate", "incorporated", "increased", "incredible", "incredulous", "indecent",
                 "independent", "individual", "individualistic ", "ineffective", "ineffectual",
                 "inept", "inevitable", "inexorable", "inexpensive", "inexperienced", "infamous",
                 "infertile", "informal", "infuriating", "injured", "innovative", "insatiable",
                 "insecure", "insidious", "inspirational   ", "inspired", "instructive", "insuperable",
                 "integrated", "intellectual", "intelligent", "intense", "intensive", "intimate",
                 "intolerant", "invaluable", "inventive", "ironic", "irresponsible", "irritable",
                 "irritating", "itchy", "jealous", "joyful", "justified", "justifying", "keen",
                 "labour", "ladylike", "lame", "large", "late", "layered", "lazy", "lean", "legitimate",
                 "leisurely", "less", "liberal", "liberating", "light", "likely", "limp", "little",
                 "loath", "locating", "lone", "lonely", "long", "loony", "loud", "lousy", "lovely",
                 "low", "loyal", "lucky", "lumbering", "luminous", "lumpy", "lunatic", "lush",
                 "mad", "magic", "magnificent", "major", "mandatory", "manipulated", "marginal",
                 "marvellous", "massive", "matrimonial", "mean", "meaningful", "measurable", "medical",
                 "medicinal", "mediocre", "mere", "mighty", "mild", "minatory", "minded", "minor",
                 "minted", "miraculous", "miscellaneous", "misleading", "mixed", "mock", "modal",
                 "modern", "modest", "modesty", "momentous", "monetary", "monstrous", "moral", "motivating",
                 "muddy", "muggy", "multiple", "mutual", "mystical", "mythical", "naive", "narrow", "nasty",
                 "naughty", "near", "nearby", "neat", "necessary", "neglected", "negligent", "nervous",
                 "net", "new", "nice", "noble", "noisy", "normal", "northern", "nostalgic", "notable",
                 "noted", "noteworthy", "noxious", "numerous", "objective", "obnoxious", "obscure",
                 "observant", "odd", "off", "oily", "okay", "old", "oldfashioned", "operatic", "optimistic",
                 "orderly", "ordinary", "orientated", "oriented", "other", "outdated", "outrageous",
                 "outstanding", "over", "overhanging", "overwhelming", "painful", "parky", "parlous",
                 "passionate", "pathetic", "patronising", "patterned", "peaked", "peculiar", "perforated",
                 "perishable", "pernicious", "perplexed", "perplexing", "persistent", "personal",
                 "persuasive", "perverted", "pessimistic", "petite", "petty", "phenomenal", "picturesque",
                 "pinkish", "plain", "pleasant", "pleased", "pleasing", "pleasurable", "plenty",
                 "poetic", "polite", "poor", "popular", "possessive", "potent", "potential", "powerful",
                 "practical", "pragmatic", "preachy", "precarious", "precious", "precise", "predatory",
                 "predictable", "prepared", "prescriptive", "pressing", "prestigious", "presumptuous",
                 "pretentious", "pretty", "prevalent", "primitive", "privileged", "prodigious", "productive",
                 "professional", "profitable", "profligate", "progressive", "prominent", "promotional",
                 "prone", "proper", "proportionate", "prospective", "prosperous", "protective", "proud",
                 "provocative", "prudential", "psycho", "psychotic", "public", "puerile", "purposeful",
                 "quaint", "qualitative", "queer", "quick", "quiet", "racist", "radical", "rainy",
                 "rampant", "rank", "rapid", "rapt", "rare", "rational", "rattled", "raw", "reactionary",
                 "reactive", "ready", "realistic", "reasonable", "recognisable", "recognised",
                 "recreational", "reddish", "reduced", "refreshing", "regretful", "regular", "relaxed",
                 "relaxing", "relentless", "relevant", "reliable", "reluctant", "remote", "required",
                 "resourceful", "respected", "responsible", "restless", "revealing", "rich", "ridiculous",
                 "risky", "robust", "rocky", "romantic", "rotten", "rough", "rowdy", "rude", "rumbling",
                 "rusty", "sacred", "sad", "safe", "sandy", "sarcastic", "satisfied", "satisfying",
                 "savage", "scarce", "scared", "sceptical", "scientific", "scrappy", "scratchy",
                 "scruffy", "scurrilous", "secret", "secular", "secure", "sedate", "seduced", "seedy",
                 "seismic", "selfconfessed", "selfish", "selfreliant", "senile", "sensational",
                 "sensible", "sensitive", "sentimental", "serious", "severe", "sexist", "sexual",
                 "sexy", "shadowy", "shaky", "shaped", "sharp", "shiny", "shitty", "shocking",
                 "short", "sick", "sickly", "silly", "simple", "sizeable", "skilful", "skilled",
                 "sleepy", "slight", "slim", "slippery", "sloppy", "slow", "small", "smart", "snoopy",
                 "snotty", "sociable", "social", "soft", "soggy", "solid", "sophisticated", "sore",
                 "sorry", "sour", "south", "spare", "sparkling", "spectacular", "spectral", "spiritual",
                 "spiteful", "splendid", "sporting", "starkly", "startling", "staunch", "steady",
                 "steamy", "steep", "stellar", "sticky", "stiff", "stimulating", "stoical", "stormy",
                 "strange", "strategic", "stressful", "stretched", "strict", "striking", "strong",
                 "structured", "stubborn", "stunning", "stupid", "subject", "subtle", "successful",
                 "suffering", "suitable", "sunny", "super", "superficial", "superior", "supernatural",
                 "supportive", "suppressed", "sure", "surplus", "surprised", "surprising", "susceptible",
                 "suspicious", "sustainable", "sweaty", "sweet", "swift", "sympathetic", "tacky",
                 "tactic", "talented", "tall", "tantalising", "tasteful", "tedious", "teensy", "temperate",
                 "tempting", "tended", "tense", "tentative", "terrible", "terrific", "theatrical",
                 "theoretical", "thermal", "thick", "thickened", "thin", "thirsty", "thoughtful",
                 "threatening", "thriving", "tight", "tiny", "tired", "titanic", "tony", "top", "topical",
                 "torrential", "tortious", "tortured", "torturous", "tough", "touring", "tragic",
                 "transcendental", "transferable", "traumatic", "treacherous", "tremendous", "trendy",
                 "tricky", "trim", "triumphal", "trivial", "troubled", "twee", "twisted", "typical",
                 "ugly", "ulterior", "unable", "unattractive", "unaware", "unbeknown", "unbelievable",
                 "uncaring", "uncertain", "unclear", "unctuous", "undecided", "undeniable", "undifferentiated",
                 "undignified", "uneven", "unexpected", "unfair", "unfamiliar", "unfavourable",
                 "unfit", "unflattering", "unforced", "unfortunate", "ungrateful", "unhappy", "unholy",
                 "unified", "unknown", "unlikely", "unlucky", "unorthodox", "unpleasant", "unreal",
                 "unseemly", "unsmiling", "unsocial", "unsound", "unstable", "unusual", "upset",
                 "uptight", "urban", "urbanised", "urgent", "useful", "vague", "vain", "valiant",
                 "valuable", "variable", "varied", "vast", "venerated", "vengeful", "versatile",
                 "vested", "veteran", "viable", "vigorous", "vile", "violent", "virtual", "visionary",
                 "visual", "vital", "vivid", "vocal", "volatile", "vulnerable", "wakeful", "warm",
                 "wayward", "weak", "weakly", "wealthy", "weary", "wee", "weird", "wet", "wicked",
                 "wide", "widespread", "wild", "willing", "wise", "wishful", "witty", "wobbly",
                 "wonderful", "wondrous", "worried", "worthwhile", "worthy", "wounded", "young",
                 "yukky", "yummy")
grd_manual <- c("delusive", "abdominal", "aboriginal", "absent", "absolute", "academic",
                "accented", "acceptable", "accessible", "accomplished", "accountable", "acoustical",
                "acrylic", "actual", "additional", "adequate", "adjacent", "administrative",
                "adolescent", "advantageous", "aerial", "affected", "affirmative", "affordable",
                "african", "aggregate", "agricultural", "albanian", "alive", "allergic", "alternative",
                "ambiguous", "american", "analogous", "analytical", "ancestral", "anecdotal",
                "angled", "anglican", "announced", "annual", "anonymous", "antarctic", "apocryphal",
                "aqueous", "arbitrary", "archaeological", "archaic", "arctic", "armed", "armoured",
                "artificial", "artistic", "asian", "asthmatic", "atmospheric", "atomic", "aussie",
                "australian", "austrian", "authorised", "automatic", "autonomous", "average", "awake",
                "back", "backward", "baked", "balanced", "bald", "bankrupt", "basic", "bearded",
                "beneficial", "best", "biblical", "bibliographic", "binding", "biodegradeable",
                "biographical", "biological", "black", "blank", "blatant", "blind", "blonde", "blue",
                "bodily", "booed", "botanical", "bottom", "british", "broke", "broken", "brown",
                "bucketful", "budgetary", "bureaucratic", "burnt", "businesslike", "busy", "californian",
                "canonical", "capitalistic", "captive", "cardiac", "catholic", "cellular", "central",
                "centralised", "centred", "certain", "characteristic  ", "chartered", "cheated",
                "chemical", "chilean", "chinese", "chivalrous", "christian", "chromatic", "chronological",
                "churchy", "classic", "classical", "clean", "clear", "close", "closed", "coarse",
                "coated", "coherent", "cohesive", "coincidental", "colloquial", "coloured", "coming",
                "commercial", "common", "compact", "comparable", "complete", "compound", "comprehensive",
                "compulsory", "computerised", "conceptual", "concrete", "confessional", "confirmed",
                "conscious", "conservative", "consistent", "constant", "constituent", "contemporary",
                "contestable", "continual", "contraceptive", "contrary", "cooked", "cooking",
                "corporate", "correct", "cracked", "crushed", "cubic", "cultural", "curly", "current",
                "customary", "cut", "daily", "dark", "dead", "deadly", "deaf", "decisive", "definite",
                "definitive", "deliberate", "democratic", "demographic", "determined", "diagnostic",
                "diagonal", "dietetic", "different", "digestive", "digital", "diplomatic", "direct",
                "discursive", "displaced", "disqualified", "distinct", "distinctive", "diverse", "divine",
                "domestic", "down", "downward", "dry", "dual", "dubious", "dummy", "dutch", "east",
                "educational", "effluent", "egalitarian", "electable", "electric", "electrical", "electronic",
                "elemental", "empty", "endemic", "endless", "english", "enough", "enrolled", "entailed",
                "entire", "equal", "equatorial", "equestrian", "equitable", "equivalent", "eritrean",
                "essential", "estonian", "ethic", "ethiopian", "ethnic", "european", "ewen", "exalted",
                "excellent", "executive", "exiguous", "existent", "existing", "exotic", "expected",
                "experimental", "explosive", "exponential", "external", "extinct", "extra", "extreme",
                "fair", "fake", "false", "far", "fatal", "favourite", "federal", "federated",
                "fellow", "female", "feminist", "feudal", "few", "fictional", "final", "financial",
                "first", "fixed", "flagged", "flannelled", "flat", "fleet", "flowing", "fluent",
                "fluid", "focused", "folded", "folding", "following", "foreign", "foremost", "formal",
                "forthcoming", "forward", "fossil", "foster", "founding", "fragile", "free", "french",
                "fresh", "frisian", "front", "frontal", "frosted", "fucking", "full", "fundamental",
                "funded", "further", "future", "gaelic", "gay", "general", "generational", "generic",
                "genuine", "geographical", "geological", "geotechnical", "german", "germanic", "glandular",
                "global", "gold", "golden", "governmental", "granulitic", "graphical", "gray", "greek",
                "green", "grey", "guaranteed", "half", "halved", "halving", "healthy", "hereditary",
                "heterogeneous   ", "heterogenious", "hidden", "historic", "historical", "holistic",
                "homosexual", "hooped", "horizontal", "hourly", "human", "humanitarian", "humiliating",
                "hungary", "hydroplaning", "hypocritical", "hypothetical", "iambic", "ideal", "identical",
                "ideological", "idle", "ill", "illegal", "imaginable", "immune", "imperial", "implicit",
                "implied", "impossible", "improved", "inaccessible", "inaccurate", "inadequate", "inclusive",
                "incoming", "incorrect", "incumbent", "indian", "indifferent", "indigenous", "indispensable",
                "indisputable", "industrial", "inefficient", "inescapable", "inexplicable", "infallible",
                "inflatable", "inflated", "informed", "infrequent", "inherent", "initial", "innate",
                "inner", "innocent", "innumerable", "inorganic", "inside", "insignificant", "instant",
                "instrumental", "insufficient", "intact", "integral", "intentional", "interactive",
                "intercultural", "interested", "interesting", "internal", "international", "interrupted",
                "intervening", "intriguing", "intrinsic", "inverted", "iraq", "irish", "irrelevant",
                "irrespective", "islamic", "italian", "japanese", "jewish", "joint", "journalistic",
                "judicial", "judicious", "junior", "just", "last", "latter", "leading", "learned",
                "learnt", "left", "lefthand", "legal", "legged", "legislative", "lesbian", "liable",
                "lime", "limited", "linear", "linguistic", "liquid", "literary", "liturgical", "live",
                "loaded", "local", "logarithmic", "logical", "logistic", "lost", "macrocyclic",
                "magnetic", "main", "male", "marine", "marked", "married", "masqueraded", "masterly",
                "materialistic   ", "maternal", "mathematical", "mature", "maximum", "mechanistic",
                "medieval", "mega", "melodic", "mental", "messy", "metamorphic", "meterological",
                "metrical", "metropolitan", "mexican", "micro", "microeconomic", "mid", "middle",
                "militaristic", "military", "milky", "minimal", "minimalist", "minimum", "ministerial",
                "missionary", "mobile", "moderate", "molecular", "molten", "monotonous", "mundane",
                "muscovite", "musical", "mutant", "naked", "narrative", "nasal", "natal", "national",
                "nationwide", "native", "natural", "nautical", "naval", "nazi", "needy", "negative",
                "neurotic", "next", "nitric", "north", "noticeable", "now", "nuclear", "obligatory",
                "obvious", "occasional", "octave", "official", "olympic", "ongoing", "only", "onward",
                "open", "operational", "opposed", "opposite", "optical", "optimum", "optional", "oral",
                "orange", "orchestral", "orchestrated", "organic", "original", "outside", "overlapping",
                "pacific", "painless", "pakistani", "parallel", "paramount", "parental", "parliamentary",
                "partial", "particular", "partisan", "passive", "past", "pastoral", "paternal",
                "paternalistic", "patriarchal", "patriotic", "perfect", "peripheral", "permanent",
                "permissive", "pertinent", "peruvian", "philosophical", "phonetic", "physical", "pink",
                "plastic", "pluralistic", "polar", "political", "politicised", "polynesian", "pornographic",
                "portable", "positive", "possible", "practicable", "preconceived", "preferential",
                "preferred", "pregnant", "preliminary", "presbyterian", "present", "presidential",
                "previous", "prewarned", "priceless", "primary", "prime", "principal", "prior", "pristine",
                "private", "privatised", "probable", "procedural", "programmed", "prolonged", "pronged",
                "proportional", "provincial", "psychiatric", "psychic", "pure", "purple", "quantifiable",
                "quantitative", "racial", "radioactive", "random", "readable", "real", "rear", "recent",
                "recycled", "red", "redundant", "reformed", "regional", "registered", "regulated",
                "regulatory", "reissued", "related", "relational", "relative", "remarkable", "remedial",
                "reportable", "reported", "residential", "respective", "resulting", "retrospective",
                "reusable", "reverse", "revolutionary", "ridged", "right", "rightful", "righthand",
                "rigid", "rigorous", "romanian", "rotary", "round", "royal", "ruined", "rural",
                "russian", "same", "samoan", "sane", "saturated", "scandinavian", "scholastic",
                "scottish", "scriptural", "seasonal", "secondary", "securing", "selected", "selective",
                "selfstyled", "senior", "senseless", "separate", "separated", "serial", "sheer",
                "siberian", "significant", "silver", "similar", "simultaneous", "sincere", "singaporean",
                "single", "skinned", "sleeveless", "sliced", "smokefree", "smooth", "sober", "socialist",
                "sociodemographic", "socioeconomic", "sole", "solitary", "soluble", "southern",
                "southwest", "sovereign", "soviet", "spanish", "special", "specialised", "specific",
                "spinal", "spontaneous", "spurious", "square", "stable", "stagnant", "standard",
                "stated", "stationary", "statistical", "statutory", "steely", "stereo", "stolen",
                "straight", "stratospheric", "striped", "structural", "subconscious", "subordinate",
                "subset", "substantial", "substantive", "suburban", "sudden", "sufficient", "suggestive",
                "sundry", "superheated", "supplementary", "supreme", "surgical", "sustained", "swedish",
                "swiss", "swollen", "symbolic", "synthetic", "technical", "technological", "temporary",
                "terminal", "territorial", "textual", "textural", "thematic", "thorough", "thoroughgoing",
                "timely", "total", "totalitarian", "toxic", "traditional", "transmitted", "traversable",
                "true", "twin", "ultimate", "unacceptable", "unaffected", "unallocated", "unannounced",
                "unanswered", "unbeaten", "unbiased", "unblemished", "unchanged", "uncoordinated",
                "under", "undisclosed", "undone", "unemployed", "unequal", "unexpired", "unfilled",
                "unfurnished", "unique", "universal", "unlimited", "unnatural", "unnecessary", "unoccupied",
                "unofficial", "unplayable", "unpopular", "unprecedented", "unprejudiced", "unpretentious",
                "unpromising", "unreceptive", "unregulated", "unrelated", "unresolved", "unrhymed",
                "unseeded", "unseen", "unselective", "unselfish", "unspecified", "unspoilt", "unstressed",
                "unsubsidised", "untold", "untrue", "unvarnished", "unwanted", "unwarranted", "unwilling",
                "unwrinkled", "upward", "usable", "useless", "usual", "utter", "vacant", "valid",
                "various", "veiled", "venetian", "verbal", "verbatim", "verifiable", "vertical",
                "volcanic", "voluntary", "weekly", "west", "western", "white", "whole", "wilful",
                "wooden", "woollen", "written", "wrong", "yellow", "youthful", "religious")
# gradability - data driven classification
# determine which adjectives are gradable
alladj <- names(table(ampnze$Adjective)) 
grd1 <- names(table(ampnze$Adjective[ampnze$Variant == "very"])) # adjs with very are gradable
grd2 <- names(table(ampnze$Adjective[ampnze$Variant == "extremely"])) # adjs with extremely are gradable
ngrd1 <- names(table(ampnze$Adjective[ampnze$Variant == "completely"])) # adjs with completely are gradable
ngrd2 <- names(table(ampnze$Adjective[ampnze$Variant == "total"])) # adjs with total are gradable
ngrd3 <- names(table(ampnze$Adjective[ampnze$Variant == "totally"])) # adjs with totally are gradable
ngrd4 <- names(table(ampnze$Adjective[ampnze$Variant == "utterly"])) # adjs with utterly are gradable
ngrd5 <- names(table(ampnze$Adjective[ampnze$Variant == "absolutely"])) # adjs with absolutely are gradable
# create vector with gradable adjectives
grdadj <- intersect(grd1, grd2)
# create vector with non-gradable adjectives
ngrdadj <- c(ngrd1, ngrd2, ngrd3, ngrd4, ngrd5)
ngrdadj <- names(table(ngrdadj))
# find elements that occur in both groups
bthgrdngrd1 <- grdadj[grdadj %in% ngrdadj]
bthgrdngrd2 <- ngrdadj[ngrdadj %in% grdadj]
bthgrdngrd <- names(table(c(bthgrdngrd1, bthgrdngrd2)))
# extract adjs that are clearly gradable 
grd_datadriven <- grdadj[!grdadj %in% bthgrdngrd]
# extract adjs that are clearly not gradable 
ngrd_datadriven <- ngrdadj[!ngrdadj %in% bthgrdngrd]
# find elements that are neither in gradable nor in nongradable
gradnongrad <- names(table(c(grd_datadriven, ngrd_datadriven)))
nagrad <- alladj[!alladj %in% gradnongrad]
naadjs <- names(table(c(nagrad, bthgrdngrd)))
###############################################################
# combine data driven and manual classification
# adj unclassified by datadriven now assigned value based on manual coding
grd_add <- naadjs[naadjs %in% grd_manual]
ngrd_add <- naadjs[naadjs %in% ngrd_manual]
# combine data driven and manual coding
grdnze <- names(table(c(grd_datadriven, grd_add)))
ngrdnze <- names(table(c(ngrd_datadriven, ngrd_add)))
# check which adjs are still unclassified
nagradnze1 <- alladj[!alladj %in% grdnze]
nagradnze <- nagradnze1[!nagradnze1 %in% ngrdnze]
# inspect unclassified adj
#nagradnze

# inspect length of vectors
length(alladj); length(grdnze); length(ngrdnze); length(nagradnze)

# add gradability coding
ampnze$Gradabilty <- ifelse(ampnze$Adjective %in% grdnze, "Gradable", ampnze$Adjective)
ampnze$Gradabilty <- ifelse(ampnze$Gradabilty %in% ngrdnze, "NotGradable", ampnze$Gradabilty)
ampnze$Gradabilty <- ifelse(ampnze$Gradabilty  == "Gradable" | ampnze$Gradabilty  == "NotGradable", 
                            ampnze$Gradabilty, "GradabilityUndetermined")
# inspect data
head(ampnze); table(ampnze$Gradabilty)

# save table of gradable and non-gradable adj
gradtb <- table(ampnze$Adjective, ampnze$Gradabilty)
gradtb <- data.frame(gradtb)
colnames(gradtb) <- c("adj", "grad", "freq")
head(gradtb)

# dave table to file
write.table(gradtb, "gradtb.txt", sep = "\t", row.names = T)
###############################################################
# add semantic types (tagliamonte 2008, based on dixon 1977)
# dimension = semdim (e.g. big, large, little, small, long, short, wide, narrow, thick)
# difficulty = semdif (e.g. difficult, simple)
# physical property = (e.g. hard, soft, heavy, light, rough, smooth, hot, sweet)
# color = semcol (e.g. black, white, red)
# human propensity: semhup (e.g. jealous, happy, kind, clever, generous, gay, rude)
# age = semage (e.g. new, young, old) 
# value (e.g. good, bad, proper, perfect, excellent, delicious, poor), 
# speed Speed (fast, quick, slow)
# position (e.g. right, left, near, far)
# other

# age
semage <- c("actual", "adolescent", "aged", "ancestral", "ancient", "annual", "archaeological", "archaic", 
            "biographical", "contemporary", "elderly", "foster", "generational", "historic", "historical", 
            "immature", "junior", "late", "mature", "medieval", "modern", "old", "oldfashioned", "outdated", 
            "past", "preliminary", "present", "primary", "prime", "prior", "puerile", "recent", "seasonal", 
            "senile", "senior", "temporary", "topical", "veteran", "young", "youthful")
# color
semcol <- c("colourful", "darkened", "pinkish", "reddish", "black", "blue", "brown",
            "coloured", "dark", "gold", "golden", "gray", "green", "grey", "lime", "marine",
            "orange", "pink", "purple", "red", "silver", "white", "yellow")
semdif <- c("basic", "complicated", "difficult", "easy", "elusive", "facile", "precarious",
            "risky", "simple", "stressful", "tricky", "twisted", "unpromising")
# dimension
semdim <- c("adjacent", "angled", "arctic", "back", "backward", "big", "bottom", "brief", "bright", 
            "broad", "central", "centralised", "centred", "close", "compact", "deep", "diagonal", 
            "direct", "distant", "distorted", "down", "downward", "early", "east", "easterly", "eastern", 
            "endemic", "endless", "equatorial", "european", "ewen", "far", "few", "first", "flat", 
            "foreign", "foremost", "forthcoming", "forward", "free", "front", "frontal", "further", 
            "geographical", "giant", "gigantic", "global", "grand", "half", "halved", "halving", 
            "high", "horizontal", "huge", "inner", "inside", "internal", "international", "large", 
            "last", "latter", "left", "linear", "little", "local", "locating", "long", "low", "massive", 
            "micro", "mid", "middle", "minimal", "minimalist", "minimum", "minor", "misleading", 
            "narrow", "national", "nationwide", "native", "near", "nearby", "next", "north", "northern", 
            "off", "onward", "orientated", "outside", "over", "overhanging", "overlapping", "pacific", 
            "parallel", "paramount", "peripheral", "petite", "polar", "proportional", "provincial", 
            "public", "rear", "regional", "remote", "reverse", "round", "rural", "separate", 
            "separated", "short", "sizeable", "slight", "small", "south", "southern", "southwest", 
            "spinal", "square", "steep", "stratospheric", "suburban", "super", "tall", "teensy", 
            "terminal", "territorial", "thick", "thickened", "thin", "tight", "tiny", "titanic", 
            "top", "torrential", "touring", "tremendous", "under", "universal", "unseeded", "upward", 
            "urban", "urbanised", "vast", "vertical", "warped", "wee", "west", "western", "wide", "widespread")
semhup <- c("able", "abrasive", "abusive", "academic", "accomplished", "advanced", "adverse", "afraid", 
            "aggressive", "aimless", "amused", "amusing", "analytical", "angry", "anxious", "appreciative", 
            "apprehensive", "ashamed", "astute", "aware", "benevolent", "besetting", "bold", "bossy", 
            "brave", "brutal", "busy", "callous", "candid", "capable", "careful", "challenging", 
            "charismatic", "cheated", "clever", "cocky", "compelling", "competent", "competitive", 
            "concerned", "confident", "consultative", "convinced", "creative", "cross", "cruel", "cute", 
            "cynical", "delighted", "depressed", "despairing", "desperate", "despondent", "disappointed", 
            "dodgy", "dotty", "dubious", "dull", "dumb", "eager", "elitist", "embarrassed", "encouraging", 
            "entertaining", "enthusiastic", "erudite", "evil", "excited", "fanatic", "fearful", "ferocious", 
            "fierce", "foolish", "forceful", "fortunate", "fraudulent", "friendly", "frustrated", "fun", 
            "funny", "furious", "generous", "gifted", "glad", "goodhearted", "gracious", "grateful", 
            "grim", "gross", "gutless", "hapless", "happy", "hopeful", "hopeless", "horrible", "hostile", 
            "hysterical", "ignorant", "ill", "imperative", "incompetent", "inexorable", "inexperienced", 
            "infallible", "informed", "insatiable", "insidious", "intellectual", "intelligent", "intriguing", 
            "inventive", "jealous", "joyful", "keen", "lazy", "learned", "learnt", "loath", "lone", "lonely", 
            "lucky", "lunatic", "mad", "mean", "minded", "motivating", "nasty", "nervous", "nice", 
            "optimistic", "passionate", "patronising", "pessimistic", "pleased", "polite", "poor", 
            "preachy", "prepared", "presumptuous", "primitive", "procedural", "professional", "promotional", 
            "proud", "prudential", "psycho", "puzzled", "rapt", "rational", "regretful", "relentless", 
            "resourceful", "respected", "rich", "romantic", "rowdy", "rude", "sad", "sane", "sarcastic", 
            "satisfied", "satisfying", "scared", "sceptical", "selective", "selfish", "sensitive", 
            "sentimental", "sick", "silly", "skilful", "skilled", "smart", "snotty", "sociable", "sophisticated", 
            "sorry", "sovereign", "spiteful", "staunch", "strategic", "strict", "stubborn", "stupid", 
            "suffering", "superior", "supportive", "suspicious", "tactic", "talented", "technical", "treacherous", 
            "troubled", "unable", "unanswered", "unaware", "uncaring", "ungrateful", "unhappy", "unsmiling", 
            "unsocial", "upset", "valiant", "valid", "vengeful", "vile", "wicked", "willing", "wise", "witty", 
            "worried")
# physical property
semphy <- c("cheap", "clear", "cold", "comfortable", "cool", "dark", "different", "dry", "flexible", "hard", 
            "heavy", "hot", "light", "neat", "obvious", "quick", "quiet", "real", "same", "scarce", "similar", 
            "slow", "strong", "sweet", "tidy", "warm")
# Value
semval <- c("amazing", "appropriate", "awful", "bad", "beautiful", "bizarre", "boring", "brilliant", 
            "competitive", "counterproductive", "easy", "effective", "efficient", "essential", "excellent", 
            "exciting", "expensive", "fantastic", "fat", "good", "great", "important", "interesting", 
            "new", "original", "painful", "pathetic", "popular", "relevant", "ridiculous", "right", 
            "scary", "serious", "simple", "special", "strange", "sure", "tough", "trendy", "true", 
            "unrealistic", "unusual", "useful", "useless", "weird", "worser", "worthwhile", "wrong", "yummy")
# add semantic type classification
ampnze$SemanticCategory <- ifelse(ampnze$Adjective %in% semage, "Age", ampnze$Adjective)
ampnze$SemanticCategory <- ifelse(ampnze$SemanticCategory %in% semcol, "Color", ampnze$SemanticCategory)
ampnze$SemanticCategory <- ifelse(ampnze$SemanticCategory %in% semdif, "Difficulty", ampnze$SemanticCategory)
ampnze$SemanticCategory <- ifelse(ampnze$SemanticCategory %in% semdim, "Dimension", ampnze$SemanticCategory)
ampnze$SemanticCategory <- ifelse(ampnze$SemanticCategory %in% semhup, "HumanPropensity", ampnze$SemanticCategory)
ampnze$SemanticCategory <- ifelse(ampnze$SemanticCategory %in% semphy, "PhysicalProperty", ampnze$SemanticCategory)
ampnze$SemanticCategory <- ifelse(ampnze$SemanticCategory %in% semval, "Value", ampnze$SemanticCategory)
ampnze$SemanticCategory <- ifelse(ampnze$SemanticCategory == "Age" | ampnze$SemanticCategory == "Color" | ampnze$SemanticCategory == "Difficulty" | 
                                    ampnze$SemanticCategory == "Dimension" | ampnze$SemanticCategory == "HumanPropensity" |
                                    ampnze$SemanticCategory == "PhysicalProperty" | ampnze$SemanticCategory == "Value",  ampnze$SemanticCategory, "NoSemType")
# table sem class of tokens
table(ampnze$SemanticCategory)

# check classification
names(table(ampnze$Adjective[which(ampnze$SemanticCategory == "NoSemType")]))

# inspect data
head(ampnze); nrow(ampnze); length(table(ampnze$Adjective))

###############################################################
# code emotion
class_emo <- get_nrc_sentiment(ampnze$Adjective)
# process sentiment
ampnze$Emotionality <- as.vector(unlist(apply(class_emo, 1, function(x){
  x <- ifelse(x[9] == 1, "emotional",
              ifelse(x[10] == 1, "emotional", "nonemotional")) } )))
# revert order of factor Emotionality
ampnze$Emotionality <- factor(ampnze$Emotionality, levels = c("nonemotional", "emotional" ))
# inspect data
head(ampnze)

###############################################################
# save raw data to disc
write.table(ampnze, "ampnze04_clean.txt", sep = "\t", row.names = F)
###############################################################
#                        END PART 2
###############################################################
