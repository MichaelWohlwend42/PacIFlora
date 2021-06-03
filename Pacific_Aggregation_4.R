###########################################
###    Pacific Data Aggregation Script  ###
### by M.R.Wohlwend v.4.0.5. 03.06.2021 ###
###########################################

setwd("C:/Users/mwohl/Documents/PacificInvaders")

## Packages
library(devtools)
library(lcvplants)
library(data.table)
library(LCVP)
library(dplyr)
library(tidyr)
# devtools::install_github("idiv-biodiversity/lcvplants")
# devtools::install_github("idiv-biodiversity/LCVP",force=T)



## Operators
'%!in%' <- function(x,y)!('%in%'(x,y))

## Functions
# source("pacific_island_unifier.txt") #run manualy
# source("pacific_species_unifier.txt")


#### Chapter 1: Harmonizing Rawdata - Pier ####

pierref<-read.csv2("PierRef2.csv",header=F)
pierref$refnr<-(sub(".*PIER reference ", "", pierref[,1]))
pierref$REF<-(sub("PIER reference .*", "", pierref[,1]))

pierref2<-read.csv2("wohlwend_ref_not_found_pier_1.csv",header=T)
pierref2$REF<-paste(pierref2$Author.s.,pierref2$Year,pierref2$Remainder.of.citation,sep=",")
pierref3<-read.csv2("wohlwend_ref_not_found_pier_2.csv",header=T)
pierref3$REF<-paste(pierref3$Author.s.,pierref3$Year,pierref3$Remainder.of.citation,sep=",")

pier2<-fread("Old_Data/PIER_SppPresences.csv") 


for(i in which(pier2$PierBibliogReference%in%pierref$refnr)){
  pier2$PierBibliogReference[i]<-pierref$REF[which(pierref$refnr==pier2$PierBibliogReference[i])]
} #gives the child a name

for(i in which(pier2$PierBibliogReference%in%as.character(pierref3$PIER.reference.number))){
  pier2$PierBibliogReference[i]<-pierref3$REF[which(pierref3$PIER.reference.number==pier2$PierBibliogReference[i])]
} #gives the child a name 2


## Cleaning & Selections

pier2<-pier2[which(!is.na(pier2$ProduceSpList4ThisRegion.)),]         # Removes corrupted/non referenced entries (n=41)
pier2<-pier2[which(pier2$PIER_region=="Pacific"),]                                # Removes non-Pacific entries
pier2<-pier2[which(pier2$Ref.Non.NativeStatus%!in%c("Eradicated","Extirpated")),] # Removes expired records

glonselect<-read.csv2("glonselect.csv")# a lookup table for ids
dat_list$tdwg4_name<-"NA"
for(i in (which(dat_list$tdwg4_id%in%glonselect$tdwg4_id))){
  dat_list$tdwg4_name[i]<-unique(glonselect$tdwg4_name[which(glonselect$tdwg4_id==dat_list$tdwg4_id[i])])
}# asigning groups from old file. Ignore warnings.

## Apply functions to standardize species & islands

pier_is<-island.unifier(pier2$FinestResolutionGeogLoc,pier2$IslandGroup) # standardizes island and group name. Ignore warnings
pier2$island<-pier_is[,1]
pier2$islandgroup<-pier_is[,2]
pier2<-pier2[which(pier2$islandgroup%!in%c("Cook")),]
pier2$islandgroup[which(pier2$islandgroup%in%c("Kiribati"))]<-"Line" #bug
pier_species<-species.unifier(pier2$TaxonCode) # standardizes species

pier2$species<-"NA"
pier2$family<-"NA"
pier2$order<-"NA"

for(i in 1:nrow(pier2)){
  pier2$species[i]<-unique(pier_species$Assigned_Name[which(pier_species$Submitted_Name==pier2$TaxonCode[i])])
  pier2$family[i]<-unique(pier_species$Family[which(pier_species$Submitted_Name==pier2$TaxonCode[i])])
  pier2$order[i]<-unique(pier_species$Order[which(pier_species$Submitted_Name==pier2$TaxonCode[i])])
} # applies species

### applying cultivation & native status

#for cultivation
pier2$cult_stat<-0
pier2$cult_stat[which(pier2$RefCitesAsCultivated.=="T")]<-1 
pier2$cult_stat[which(pier2$RefCitesAsCultivated.=="")]<-0.5 

#for native status
nat_rep<-matrix(0,11,2)
nat_rep[,1]<-c("1","Native","Probably native","Possibly native","Uncertain if native","NULL","?",
               "Possibly introduced","Probably introduced","Introduced","0")
nat_rep[,2]<-c(1,1,0.8,0.6,0.5,0.5,0.5,0.4,0.2,0,0)
colnames(nat_rep)<-c("Name","Score")
nat_rep<-as.data.frame(nat_rep)
nat_rep$Score<-as.numeric(as.character(nat_rep$Score))
nat_rep$Name<-as.character(nat_rep$Name)

pier2$nat_score<-0
for (i in 1:nrow(nat_rep)){
  pier2$nat_score[which(pier2$Ref.Non.NativeStatus==nat_rep$Name[i])]<-nat_rep$Score[i]
}
#this two stati have to be adjusted later as at this point in rare cases unique species/island combinations can have two different stati



#### Chapter 2: Rawdata Glonaf ####


dat_list2<-dat_list[which(dat_list$region_id%in%glonselect$region_id),]#which islands do you want
dat_flora_orig2<-dat_flora_orig[which(dat_flora_orig$list_id%in%dat_list2$list_id),]# which species lists do you want
glon_lit<-read.csv2("glonafdb_2021-05-31_recnum.csv",header=T)



dat_taxon2<-dat_taxon[which(dat_taxon$taxon_orig_id%in%dat_flora_orig2$taxon_orig_id),]# which species are in that list
hyb<-str_count(dat_taxon2$taxon_orig, pattern = " x ")
dat_taxon2$Usage[which(hyb==0)]<-sub("^(\\S*\\s+\\S+).*", "\\1", dat_taxon2$taxon_orig[which(hyb==0)]) #shortens output to binomial style
dat_taxon2$Usage[which(hyb==1)]<-word(dat_taxon2$taxon_orig[which(hyb==1)], 1, 3)
dat_taxon2$Usage[which(hyb==1)][3]<-"Kyllinga nemoralis" # one author had a " x " in the name

glon_3000<-species.unifier(dat_taxon2$Usage) # so it can get lcvp-bashed
for(i in 1:nrow(dat_taxon2)){
  dat_taxon2$Assigned_Name[i]<-unique(glon_3000$Assigned_Name[which(glon_3000$Submitted_Name==dat_taxon2$Usage[i])])
} #assigns the new name

for(i in 1:nrow(dat_flora_orig2)){
  dat_flora_orig2$Original_name[i]<-dat_taxon2$taxon_orig[which(dat_taxon2$taxon_orig_id==dat_flora_orig2$taxon_orig_id[i])]
  dat_flora_orig2$Assigned_Name[i]<-dat_taxon2$Assigned_Name[which(dat_taxon2$taxon_orig_id==dat_flora_orig2$taxon_orig_id[i])]
} # add that to the island + list file

Is_gr_glon<-island.unifier(dat_list2$region_name,dat_list2$tdwg4_name)
dat_list2$tdwg4_name2<-Is_gr_glon[,2]
dat_list2$region_name2<-Is_gr_glon[,1] # extracting, harmonizing and assigning island and group name


for(i in 1:nrow(dat_flora_orig2)){
  dat_flora_orig2$island[i]<-dat_list2$region_name2[which(dat_list2$list_id==dat_flora_orig2$list_id[i])]
  dat_flora_orig2$island_group[i]<-dat_list2$tdwg4_name2[which(dat_list2$list_id==dat_flora_orig2$list_id[i])]
} # and merging this to the complete file

### Asigning a standardized value for cultivated, native, invasive
dat_flora_orig2$invasive<-dat_flora_orig2$status
dat_flora_orig2$cultivated<-dat_flora_orig2$status
dat_flora_orig2$cultivated[which(dat_flora_orig2$cultivated%in%c("alien","casual"))]<-0.5
dat_flora_orig2$cultivated[which(dat_flora_orig2$cultivated%in%c("cultivated"))]<-1
dat_flora_orig2$cultivated[which(dat_flora_orig2$cultivated%!in%c("1","0.5"))]<-0 
dat_flora_orig2$invasive[which(dat_flora_orig2$invasive=="invasive")]<-1
dat_flora_orig2$invasive[which(dat_flora_orig2$invasive!="1")]<-0 # splitting status into invasive and cultivated information
dat_flora_orig2$presence_uncertain[which(is.na(dat_flora_orig2$presence_uncertain))]<-0.5 # rare, as communicated with anke stein
dat_flora_orig2<-dat_flora_orig2[which(dat_flora_orig2$presence_uncertain!=1),] # removing uncertainties

table(dat_list2$endnote_rec_num%in%glon_lit$Nr)

for (i in 1:nrow(dat_list2)){
  dat_list2$source[i]<-glon_lit$Source[which(glon_lit$Nr==dat_list2$endnote_rec_num[i])]
}
unique(dat_list2$source)





for(i in 1:nrow(dat_flora_orig2)){
dat_flora_orig2$reference[i]<-dat_list2$source[which(dat_list2$list_id==dat_flora_orig2$list_id[i])]
}
glon2<-dat_flora_orig2[,c("Assigned_Name","island","island_group","invasive","cultivated","questionably_native","Original_name","reference")]
glon2$questionably_native[which(is.na(glon2$questionably_native))]<-0 # as communicated with anke stein
glon2$sp_is<-paste(glon2$Assigned_Name,glon2$island) #simplifies identifying unique entries later
glon3<-glon2

for(i in which(glon3$Assigned_Name%in%glon_3000$Assigned_Name)){
  glon3$family[i]<-unique(glon_3000$Family[which(glon_3000$Assigned_Name==glon3$Assigned_Name[i])])
  glon3$order[i]<-unique(glon_3000$Order[which(glon_3000$Assigned_Name==glon3$Assigned_Name[i])])
}

#### - Your additional Data - not coded for efficiency/large datasets####
#your_names <- species.unifier(yourdata$species)
#yourdata$species<-your_names[,1]
#yourdata$family<-your_names[,2]
#yourdata$order<-your_names[,3]
#your_islands<-island.unifier(your_data$island,yourdata$islandgroup) # if you don't have island group, enter island twice, it should work
#yourdata$island<-your_islands[,1]
#yourdata$islandgroup<-your_islands[,2]


#### Chapter 3: Combination ####

glon3$dataset<-"glon"
pier2$dataset<-"pier"

glon3<-glon3[,c("Assigned_Name","island","island_group","questionably_native","cultivated","family","order","dataset","invasive","Original_name","reference")]
pier3<-pier2[,c("species","island","islandgroup","nat_score","cult_stat","family","order","dataset","RefCitesAsInvasive.","TaxonCode","PierBibliogReference")]
colnames(pier3)[9]<-"inva_stat"
colnames(pier3)[10]<-"orginal_name"
colnames(pier3)[11]<-"reference"

colnames(glon3)<-colnames(pier3)

pacific_is4<-rbind(glon3,pier3) # I recommend to bind in your data at this point if you have some information of cultivation nativeness etc. and later when you have non.
pacific_is4$islandgroup[which(pacific_is4$islandgroup=="Midway")]<-"NW Hawaiian"#bug
apply_this_4<-island.unifier(pacific_is4$island,pacific_is4$islandgroup) #security measures
pacific_is4$island<-apply_this_4[,1]
pacific_is4$islandgroup<-apply_this_4[,2]

pacific_is4$inva_stat[which(pacific_is4$inva_stat=="T")]<-1
pacific_is4$inva_stat[which(pacific_is4$inva_stat=="F")]<-0
pacific_is4$inva_stat[which(pacific_is4$inva_stat=="")]<-0.5
pacific_is4$inva_stat<-as.numeric(pacific_is4$inva_stat)

pacific_is5<-pacific_is4[which(!is.na(pacific_is4$islandgroup)),]
#islands can be NA, but if also no group is present we have no idea where this location is situated
#there should never be the case where an island has no group. if so it needs to be included in the .csv file for cleaning manually, or this step must be adjusted accordingly

pacific_is5$sp_ar<-paste(pacific_is5$species,pacific_is5$islandgroup,sep = "_")
pacific_is5$sp_is<-paste(pacific_is5$species,pacific_is5$island,sep = "_")

# apply standardized dataset to be able to safely sort out duplicates without loosing information later
pacific_is5$dataset_ar<-pacific_is5$dataset
for( i in 1:nrow(pacific_is5)){
  ifelse(length(unique(pacific_is5$dataset[which(pacific_is5$sp_ar==pacific_is5$sp_ar[i])]))>1,pacific_is5$dataset_ar[i]<-"glonpier",NA)
}

for( i in 1:nrow(pacific_is5)){
  ifelse(length(unique(pacific_is5$dataset[which(pacific_is5$sp_ar==pacific_is5$sp_ar[i])]))>1,pacific_is5$dataset[i]<-"glonpier",NA)
}

## Form lowest degree of cultivation & mean degree of "nativeness" for each island/islandgroup

doubles<-pacific_is5$species[which(duplicated(pacific_is5$sp_is))]
doubles<-doubles[which(!is.na(doubles))]
doubles<-doubles[which(doubles!="NA")]

cult_sorter<-table((pacific_is5[which(pacific_is5$species%in%doubles),c("sp_is","cult_stat")]))
naturalized<-names(cult_sorter[which(cult_sorter[,1]!=0),1])
pacific_is5$cult_stat[which(pacific_is5$sp_is%in%naturalized)]<-0 
# so what this does is reassign a "certainly naturalized status to duplicated records who may/may not have this information
# if you look at the steps in between you will see that it is rarely the case that a species is recored as maybe cultivated & not cultivated on an island group 
# this had to be done because otherwise in the next step, a "random" naturalization status would have been selected

pacific_is5$nat_score2<-pacific_is5$nat_score
pacific_is5$nat_score<-pacific_is5$nat_score2

for(i in 1:nrow(pacific_is5)){
  pacific_is5$nat_score[i]<-mean(pacific_is5$nat_score2[which(pacific_is5$sp_is==pacific_is5$sp_is[i])])
} # reach consensus among authors on native status. For simplicity, everything over 0.5 will be treated as native. 

natives<-pacific_is5$species[which(pacific_is5$nat_score>0.51)]
exotics<-pacific_is5$species[which(pacific_is5$nat_score<0.51)]# this ensures that uncertain cases are included 

true_natives<-unique(natives[which(natives%!in%exotics)])
pacific_is5<-pacific_is5[which(pacific_is5$species%!in%true_natives),]#eliminates species with only native records

cultivars<-pacific_is5$species[which(pacific_is5$cult_stat>0.51)]
wilds<-pacific_is5$species[which(pacific_is5$cult_stat<0.51)]# this ensures that uncertain cases are included 

true_cults<-unique(cultivars[which(cultivars%!in%wilds)])
pacific_is5<-pacific_is5[which(pacific_is5$species%!in%true_cults),]# eliminates all species with only cultivated records

pacific_is5$inva_stat2<-pacific_is5$inva_stat
pacific_is5$inva_stat<-pacific_is5$inva_stat2
for(i in 1:nrow(pacific_is5)){
  pacific_is5$inva_stat[i]<-max(pacific_is5$inva_stat2[which(pacific_is5$sp_is==pacific_is5$sp_is[i])])
} # reach consensus among authors on invasive status. 

pacific_is5<-pacific_is5[which(pacific_is5$islandgroup!=""),]

pacific_is5_rebind<-pacific_is5[which((pacific_is5$species=="NA")),]# takes NAs out to rebind later
mynewdoubles<-pacific_is5$sp_is[which(duplicated(pacific_is5$sp_is))]
b_pacific_is5<-pacific_is5[which(pacific_is5$sp_is%in%mynewdoubles), ]
b_pacific_is5$rspis<-paste(b_pacific_is5$sp_is,b_pacific_is5$reference)
c_pacific_is5<-b_pacific_is5[which(!duplicated(b_pacific_is5$rspis)),]
c_pacific_is5<-c_pacific_is5[which(c_pacific_is5$species!="NA"),]
c_pacific_is5<-c_pacific_is5[which(!is.na(c_pacific_is5$species)),]

duplicated_dataset<-c_pacific_is5[,c(1:11)]

ref_table<-matrix(0,22,3)
ref_table[,1]<-names(head(sort(table(pacific_is5$reference),decreasing = T),22))
ref_table[,2]<-as.numeric(head(sort(table(pacific_is5$reference),decreasing = T),22))
for (i in 1:nrow(ref_table)){
  ref_table[i,3]<-paste(names(sort(table(pacific_is5$islandgroup[which(pacific_is5$reference==ref_table[i,1])]),decreasing = T)),collapse = ",")
} # what this does: Sort out the island groups the authors we use the most deal with

write.csv2(ref_table,"Reference_Overview.csv",row.names = F)
myref_2<-matrix(0,length(unique(pacific_is5$reference)),2)
myref_2[,1]<-unique(pacific_is5$reference)
myref_2[,2]<-c(1:nrow(myref_2))



for (i in 1:nrow(pacific_is5)){
  pacific_is5$refid[i]<-paste(c(myref_2[which(myref_2[,1]==pacific_is5$reference[i]),2]),collapse = ",")
}

testref<-matrix(0,nrow(pacific_is5),4)
testref[,1]<-pacific_is5$sp_is
testref[,2]<-pacific_is5$refid
testref[,3]<-pacific_is5$reference
for(i in 1:nrow(pacific_is5)){
pacific_is5$refid[i]<-paste(c(unique(testref[which(testref[,1]==pacific_is5$sp_is[i]),2])),collapse = ",")
}

for (i in 1:nrow(testref)){
  testref[i]<-length((myref_2[which(myref_2[,1]==pacific_is5$reference[i]),2]))
} 


write.csv2(testref[which(!duplicated(testref[,2]))],"Reference_table_0106.csv")
write.table(testref[which(!duplicated(testref[,2]))],"Reference_table_0306.txt",row.names=F)
#okay so what this all was for to asign each entry all references that where found for it, but without pasting the text & profiding an extra table to look up where else they where found


pacific_is5_b<-pacific_is5[which(!duplicated(pacific_is5$sp_is)),]

pacific_is5_rebind<-pacific_is5_rebind[which(rownames(pacific_is5_rebind)%!in%rownames(pacific_is5)),]

pacific_is6<-rbind(pacific_is5_b,pacific_is5_rebind)# positive sideeffect: Now all NA species are at the bottom of the list
##### Here, you could bind your data when you have no status information ####
## pacific_is7<-rbind(pacific_is6,yourdata)

#now every unique species island combination has one value


## apply for islandgroups
pacific_ar5<-pacific_is6 # from now on you have to work with two different datasets.

pacific_ar6<-pacific_ar5[,c("sp_ar","dataset")] #dummy

for(i in which(duplicated(pacific_ar6$sp_ar))){
  ifelse(length(unique(pacific_ar6$dataset[which(pacific_ar6$sp_ar==pacific_ar6$sp_ar[i])]))>1,pacific_ar5$dataset[i]<-"glonpier",NA)
} # assigns double dataset status for island

doubles<-pacific_ar5$species[which(duplicated(pacific_ar5$sp_ar))]
doubles<-doubles[which(doubles!="NA")]
doubles<-doubles[which(!is.na(doubles))]# to make sure species that could not be assigned an accepted name are not excluded.
cult_sorter<-table((pacific_ar5[which(pacific_ar5$species%in%doubles),c("sp_ar","cult_stat")]))
naturalized<-names(cult_sorter[which(cult_sorter[,1]!=0),1])
pacific_ar5$cult_stat[which(pacific_ar5$sp_ar%in%naturalized)]<-0 
# so what this does is reassign a "certainly naturalized status to duplicated records who may/may not have this information
# if you look at the steps in between you will see that it is rarely the case that a species is recorded as maybe cultivated & not cultivated on an island group 
# this had to be done because otherwise in the next step, a "random" naturalization status would have been selected

pacific_ar5$nat_score2<-pacific_ar5$nat_score
pacific_ar5$nat_score<-pacific_ar5$nat_score2 # safety

for(i in 1:nrow(pacific_ar5)){
  pacific_ar5$nat_score[i]<-mean(pacific_ar5$nat_score2[which(pacific_ar5$sp_ar==pacific_ar5$sp_ar[i])])
} # reach consensus among authors on native status. For simplicity, everything over 0.5 will be treated as native. 

natives<-pacific_ar5$species[which(pacific_ar5$nat_score>0.51)]
exotics<-pacific_ar5$species[which(pacific_ar5$nat_score<0.51)]# this ensures that uncertain cases (0.5) are included 

true_natives<-unique(natives[which(natives%!in%exotics)])
pacific_ar5<-pacific_ar5[which(pacific_ar5$species%!in%true_natives),]

cultivars<-pacific_ar5$species[which(pacific_ar5$cult_stat>0.51)]
cultivars<-cultivars[which(cultivars!="NA")]
wilds<-pacific_ar5$species[which(pacific_ar5$cult_stat<0.51)]# this ensures that uncertain cases (0.5) which included 

true_cults<-unique(cultivars[which(cultivars%!in%wilds)])
pacific_ar5<-pacific_ar5[which(pacific_ar5$species%!in%true_cults),]


pacific_ar5_rebind<-pacific_ar5[which(is.na(pacific_ar5$species)),]
pacific_ar5<-pacific_ar5[which(!duplicated(pacific_ar5$sp_ar)),]
pacific_ar5_rebind<-pacific_ar5_rebind[which(rownames(pacific_ar5_rebind)%!in%rownames(pacific_ar5)),]
pacific_ar6<-rbind(pacific_ar5,pacific_ar5_rebind)
pacific_ar6<-pacific_ar6[which(pacific_ar6$islandgroup!=""),]

pacific_is6$nat_score<-round(pacific_is6$nat_score,3)#due to technical reviewer request
pacific_is6$species[which(pacific_is6$species=="unresolved")]<-NA #due to technical reviewer request
pacific_is6$species[which(pacific_is6$species=="Agropogon_x lutosus")]<-"Agropogon lutosus_x" # altered  technical reviewer request

aberhallo<-pacific_is6$orginal_name
pacific_is6$orginal_name<-gsub("\r?\n|\r", " ", pacific_is6$orginal_name)# thanks @ technical review


pacific_is7<-pacific_is6[which(!is.na(pacific_is6$island)),]

## With that, pacific_ar6 & pacific_is7 are complete lists. pacific_is6 contains all information 


#### Chapter 4: Matrix aggregation ####
library(fuzzySim) # makes life easier
pacgr<-splist2presabs(data=pacific_ar6,sites.col = "islandgroup",sp.col = "species")
rownames(pacgr)<-pacgr[,1]
pacgr<-pacgr[,c(2:ncol(pacgr))]


pacis<-splist2presabs(data=pacific_is7,sites.col = "island",sp.col = "species")
rownames(pacis)<-pacis[,1]
pacis<-pacis[,c(2:ncol(pacis))]

# ### smaller lists without native/cultivated records - to be provided as code only. Maybe I should write a function
# 
# pacific_ar_s<-pacific_ar[which(pacific_ar$nat_score<0.5),]
# 
# pacgr_s<-splist2presabs(data=pacific_ar_s,sites.col = "islandgroup",sp.col = "species")
# rownames(pacgr_s)<-pacgr_s[,1]
# pacgr_s<-pacgr_s[(1:nrow(pacgr_s)),c(2:ncol(pacgr_s))]


#### Chapter 5: coordinates ####

pacific_df <- read.csv("Pacific_Df.csv",sep=";",dec=",")

gift_holder<-island.unifier(pacific_df$geo_entity,pacific_df$archipelago)
pacific_df$island<-gift_holder[,1]
pacific_df$islandgroup<-gift_holder[,2]

pacific_df<-pacific_df[which(!is.na(pacific_df$island)),]


myapp<-island.unifier(pacific_is7$island,pacific_is7$islandgroup)

pacific_is7$island<-myapp[,1]


pacific_is7$latitude<-"NA"
pacific_is7$longitude<-"NA"
for (i in which(pacific_is7$island%in%pacific_df$island)){
  pacific_is7$longitude[i]<-pacific_df$point_x[which(pacific_df$island==pacific_is7$island[i])]
  pacific_is7$latitude[i]<-pacific_df$point_y[which(pacific_df$island==pacific_is7$island[i])]
}

google_coords<-read.csv2("to_coord_table_2.csv")
for( i in which(pacific_is7$island%in%google_coords$island)){
  pacific_is7$longitude[i]<-google_coords$latitude[which(google_coords$island==pacific_is7$island[i])]
  pacific_is7$latitude[i]<-google_coords$longitude[which(google_coords$island==pacific_is7$island[i])]
}



myapp<-island.unifier(pacific_is6$island,pacific_is6$islandgroup)

pacific_is6$island<-myapp[,1]


pacific_is6$latitude<-"NA"
pacific_is6$longitude<-"NA"
for (i in which(pacific_is6$island%in%pacific_df$island)){
  pacific_is6$longitude[i]<-pacific_df$point_x[which(pacific_df$island==pacific_is6$island[i])]
  pacific_is6$latitude[i]<-pacific_df$point_y[which(pacific_df$island==pacific_is6$island[i])]
}

google_coords<-read.csv2("to_coord_table_2.csv")
for( i in which(pacific_is6$island%in%google_coords$island)){
  pacific_is6$longitude[i]<-google_coords$latitude[which(google_coords$island==pacific_is6$island[i])]
  pacific_is6$latitude[i]<-google_coords$longitude[which(google_coords$island==pacific_is6$island[i])]
}


###
myapp<-island.unifier(pacific_is6$island,pacific_is6$islandgroup)

pacific_is6$island<-myapp[,1]
pacific_is6$islandgroup[which(pacific_is6$island=="Temotu")]

pacific_is6$latitude<-"NA"
pacific_is6$longitude<-"NA"


for (i in which(pacific_is6$island%in%pacific_df$island)){
  pacific_is6$longitude[i]<-pacific_df$point_x[which(pacific_df$island==pacific_is6$island[i])]
  pacific_is6$latitude[i]<-pacific_df$point_y[which(pacific_df$island==pacific_is6$island[i])]
}

google_coords<-read.csv2("to_coord_table_2.csv")
for( i in which(pacific_is6$island%in%google_coords$island)){
  pacific_is6$longitude[i]<-google_coords$longitude[which(google_coords$island==pacific_is6$island[i])]
  pacific_is6$latitude[i]<-google_coords$latitude[which(google_coords$island==pacific_is6$island[i])]
}

pacific_is6<-pacific_is6[,which(colnames(pacific_is6)%!in%c("nat_score2","inva_stat2","sp_is","sp_ar","dataset_ar"))]
colnames(pacific_is6)[4:5]<-c("Native","Cultivated")
colnames(pacific_is6)[9]<-"Invasive"
colnames(pacific_is6)[10]<-"Original_name"

unused_ref<-unique(pier3$reference)[which((unique(pier3$reference)%!in%unique(pacific_is6$reference)))] # all glonaf references are included. 

write.table(unused_ref,"Additional_references.txt")
write.csv2(pacific_ar6,"Species_list_group_2905.csv")
write.csv2(pacific_is6,"Species_list_full_2905.csv",row.names = F)


colnames(pacific_is6)<-c("Species","Island","IslandGroup","Native","Cultivated","Family","Order","Database","Invasive","Original_name","Source","Source_ID","Latitude","Longitude")

write.table(pacific_is6,"PacIFlora.txt",row.names = F)
write.csv2(pacific_is7,"Species_list_island_2905.csv")
#write.csv2(pacgr_s,"Species_mat_pacific_ars_2905.csv")
write.csv2(pacgr,"Species_mat_group_2905.csv")
write.csv2(pacis,"Species_mat_island_2905.csv")


nrow(pacific_is6[which(pacific_is6$Invasive>0.51),])
length(unique(pacific_is6$Species[which(pacific_is6$Cultivated>0.51)]))
length(unique(pacific_is6$Species[which(pacific_is6$Native>0.51)]))
length(unique(pacific_is6$Species[which(pacific_is6$Invasive>0.51)]))
length(unique(pacific_is6$Species))
nrow(pacific_is6)
nrow(pacific_is6[which(is.na(pacific_is6$Island)),])
nrow(pacific_is6[which(is.na(pacific_is6$Species)),])
nrow(pacific_is6[which((pacific_is6$Species=="NA")),])
length(unique(pacific_is6$Island))
hist(pacific_is6$Native)
