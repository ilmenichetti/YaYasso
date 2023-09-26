
##################################################################################################################
# Reconstruct figures of Tupek et al. 2023: 
# Modeling boreal forest’s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier
#
###################################################################################################################

# Boris Tupek
# Natural Resources Institute Finland (LUKE)
# June 2023
# boris.tupek@luke.fi

################################################################################################################
#code below produces LITTER INPUT based on the tree stand inventory (data to plot Fig S2 and Fig S3)############
#########################
#alternatively load files:
#load modeled biomass and litter by method of finnish GHG described in Tupek et al 2019
#load(file=paste(paste0(path.data,"biomass.litter_trees.understory.RData"),sep=""))
#includes "ul" understory,"tl"- trees, "tl.msd" agregated for species, "tul.sum" agre. for sites,
#"tul.bp" for plotting (NOTE tree finerots excluded from the litte input as the tree roots were cut)

#load litter input from trenching for 2004 JULY Time series input!!!
#load(paste(paste0(path.data,"litter.trench.cut_coarse.fineroots.RData"),sep=""))
#includes "tree.w.coarseroot.cut04.awen.sum", "treeunder.nw.fineroot.cut04.awen"

#load lists of simulated AWEN data woody and nowoody litter 
#load(file=paste(paste0(path.data,"litter.eco_equilibrium_months.RData"),sep=""))
#includes lists "litter.eq.awen", "litter.natural.month.awen","litter.fc.month.awen"

## Tree stand measurements ##########################################################################################################
stands <- read.csv(paste0(path.data,'Stands.ecotone.csv'), header=T, sep = ",",stringsAsFactors = FALSE, na.strings=c("NA",""))

# CN ratio
c.cont.s <- c(43.17,24.22,49.63,47.09,45.36,48.68,50.30,45.76,48.20)
n.cont.s <- c(1.02,0.61,1.18,1.59,2.19,1.47,1.12,1.29,0.96)
cn.o <-  c.cont.s/n.cont.s 
cn.o.df <-data.frame(plot=1:9,cn.o=cn.o )
stands <-merge(stands,cn.o.df,by="plot")
 #View(stands)

#reorganize variables of the stand table ##########
stands4 <- stands[,c(c("plot","species","n.400","n.ha","cn.o","dmed_cm","hmed_m","hcb_m"), names(stands)[grep("4_",names(stands))])]
stands5 <- stands[,c(c("plot","species","n.400","n.ha","cn.o","dmed_cm","hmed_m","hcb_m"), names(stands)[grep("5",names(stands))])]
stands6 <- stands[,c(c("plot","species","n.400","n.ha","cn.o","dmed_cm","hmed_m","hcb_m"), names(stands)[grep("6",names(stands))])]

stands4$year<-2004
stands5$year<-2005
stands6$year<-2006

names(stands4)<-sub("4_","",names(stands4))
names(stands5)<-sub("5_","",names(stands4))
names(stands6)<-sub("6_","",names(stands4))

stands3<-stands4
stands3$year<-2003
names(stands3)
stands3$d13.cm <- stands3$dmed_cm
stands3$hm <- stands3$hmed_m

stands.y <- rbind(stands3,stands4,stands5,stands6)
#add basal area m2/ha
stands.y$basal.m2ha <- round(stands.y$n.ha*pi*(stands.y$d13.cm/200)^2,2)

names(stands.y)
stands.y<-stands.y[,c("plot","species","year","cn.o","n.ha","basal.m2ha","d13.cm","hm","hcb_m")]
stands.y<-stands.y[order(stands.y$plot),]
#View(stands.y)

#for tree biomass models separate stand data by species 
#1 pine, 2 spruce, 3 deciduous
stands.pi <- subset(stands.y,species==1)
stands.su <- subset(stands.y,species==2)
stands.de <- subset(stands.y,species==3)

#for understory biomass models plot level means
library(doBy)
#stands weighted age and height
stands.wmean <- stands[,c("plot","species","age","n.400","n.ha","dmed_cm","hmed_m")]
#View(stands.wmean)

n.sum <-summaryBy(n.ha ~ plot, stands.wmean, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
stands.wmean <- merge(stands.wmean, n.sum, by="plot")
#weighted means
stands.wmean$age.w0 <- with(stands.wmean,age*n.ha/n.ha.sum)
stands.wmean$d.w0 <- with(stands.wmean,dmed_cm*n.ha/n.ha.sum)
stands.wmean$h.w0 <- with(stands.wmean,hmed_m*n.ha/n.ha.sum)

## TREE BIOMASS REPOLA ###########################################################################################
#
#Functions for biomass  calculation
source(paste0(path.scripts,"biomass.core.r")) #biomass.core.r is needed for BIOMASS_functions.r
source(paste0(path.scripts,"BIOMASS_functions.r"))

#foliage #kg/ha
repola.pine.needles
repola.spruce.needles
stands.pi$needles <- mapply(repola.pine.needles, stands.pi$d13.cm,stands.pi$hm,stands.pi$hcb_m)*stands.pi$n.ha
stands.su$needles <- mapply(repola.spruce.needles, stands.su$d13.cm,stands.su$hm,stands.su$hcb_m)*stands.su$n.ha
stands.de$needles <- mapply(repola.birch.foliage, stands.de$d13.cm,stands.de$hm,stands.de$hcb_m)*stands.de$n.ha

#living branch
repola.pine.livbranch
repola.spruce.livbranch
stands.pi$branch <- mapply(repola.pine.livbranch, stands.pi$d13.cm,stands.pi$hm,stands.pi$hcb_m)*stands.pi$n.ha
stands.su$branch <- mapply(repola.spruce.livbranch, stands.su$d13.cm,stands.su$hm,stands.su$hcb_m)*stands.su$n.ha
stands.de$branch <- mapply(repola.birch.livbranch, stands.de$d13.cm,stands.de$hm,stands.de$hcb_m)*stands.de$n.ha

#stem
repola.pine.stembark.simp
repola.spruce.stembark.simp
stands.pi$stembark <- mapply(repola.pine.stembark.simp, stands.pi$d13.cm,stands.pi$hm)*stands.pi$n.ha
stands.su$stembark <- mapply(repola.spruce.stembark.simp, stands.su$d13.cm,stands.su$hm)*stands.su$n.ha
stands.de$stembark <- mapply(repola.birch.stembark.simp, stands.de$d13.cm,stands.de$hm)*stands.de$n.ha

#stump
repola.pine.stump #(d13rm)
repola.spruce.stump 
stands.pi$stump <- mapply(repola.pine.stump, stands.pi$d13.cm)*stands.pi$n.ha
stands.su$stump <- mapply(repola.spruce.stump, stands.su$d13.cm)*stands.su$n.ha
stands.de$stump <- mapply(repola.birch.stump, stands.de$d13.cm)*stands.de$n.ha

#coarse.roots
repola.pine.roots #(d13rm)
repola.spruce.roots
stands.pi$coarseroots <- mapply(repola.pine.roots, stands.pi$d13.cm)*stands.pi$n.ha
stands.su$coarseroots <- mapply(repola.spruce.roots, stands.su$d13.cm)*stands.su$n.ha
stands.de$coarseroots <- mapply(repola.birch.roots, stands.de$d13.cm,stands.de$hm)*stands.de$n.ha

#fine.roots
#Lethtonen et al. 2015 FORECO 
fineroots # (basal,decid,cno) decid = 0/1

stands.pi$fineroots <- mapply(fineroots, stands.pi$basal.m2ha,0,stands.pi$cn.o)
stands.su$fineroots <- mapply(fineroots, stands.su$basal.m2ha,0,stands.su$cn.o)
stands.de$fineroots <- mapply(fineroots, stands.de$basal.m2ha,1,stands.de$cn.o)

#merge data
stands.bio<- rbind(stands.pi,stands.su,stands.de)

#sum biomass components 
stands.bio$biomass <- rowSums(stands.bio[,c("needles","branch","stembark","stump","coarseroots","fineroots")])

### CONVERT TO CARBON ##########################

carbon <- function(M) { ###  function for carbon!!!
  return(M*0.5)
}
#convert all to carbon!!!
names(stands.bio)
stands.bio[,c("needles","branch","stembark","stump","coarseroots","fineroots","biomass")] <- 
  carbon(stands.bio[,c("needles","branch","stembark","stump","coarseroots","fineroots","biomass")])

biomass.c <- stands.bio[,c("plot","species", "year","needles","branch","stembark","stump","coarseroots","fineroots","biomass")]
#View(stands.bio)

#stope

################################################################################################
## Functions for LITTER calculation ##############################################################
#load litter functions
source(paste0(path.scripts,"LITTER_functions.r"))
source(paste0(path.scripts,"litter_updated.r"))

lsf.str()[grep("litter",lsf.str())]

#add litter data to biomass table
biolit<-biomass.c

#update foliage species specific litter/biomass turnover rates from Tupek et al. 2015
#for southern finland
#needle.turnover.spruce <- 0.139
#needle.turnover.pine <- 0.278
foliage.litter <- function(Mb,spec){
  if (spec==1) {#pine
    return(Mb*0.278)}
  if (spec==2) {#spruce
    return(Mb*0.139)}
  if (spec==3) {
    return(Mb*0.79)}
}
biolit$lit.needle <- mapply(foliage.litter,biolit$needles,biolit$species)

#branch litter 
biolit$lit.branch <- mapply(branch.litter, biolit$branch,biolit$species)

#stem 
#"stem.pine","stem.spruce"
biolit$lit.stem <- mapply(rep.stem.litter, biolit$stembark,biolit$species)

#stump
# "stump.pine","stump.spruce"
biolit$lit.stump <-mapply(stump.litter, biolit$stump, biolit$species)

#root 
#"root.pine","root.spruce"
#sum coarse.roots-repola 2009 estimate and coroot.t - measured 2-5mm
biolit$lit.coarseroot <- mapply(root.litter, biolit$coarseroots, biolit$species)

#fine.root
#proportion of maximum temperature(1=max,0=min) needed for fine.root litter
#temp.ratio <- (c$Tair.mean+abs(min(c$Tair.mean)))/diff(range(c$Tair.mean))
#use function with NO correction for temperature 
fineroot.litter.old <- function(Mf,spec) {
  if (spec==1) {
    return(Mf*0.868)} #uncorrected Mf*0.868
  if (spec==2) {
    return(Mf*0.811)} #uncorrected Mf*0.811
  if (spec==3) {
    return(Mf*1)} #uncorrected Mf*1
}
biolit$lit.fineroot <- mapply(fineroot.litter.old,biolit$fineroots,biolit$species)

#sum litter components and convert to carbon
biolit$litter <- rowSums(biolit[,c("lit.needle","lit.branch","lit.stem","lit.stump","lit.coarseroot","lit.fineroot")])

#total biomass+litter 
biolit$tot.biolit <- biolit$biomass + biolit$litter

#NPP from annual differences
biolit$c.npp.kgha <- c(NA,diff(biolit$biomass))+biolit$litter
summary(biolit$c.npp.kgha/10000*1000)
#441 g/m2/yr #NPP Boreal Forest: Kuusamo, Finland, 1967-1972, R1 #https://daac.ornl.gov/cgi-bin/dataset_lister.pl?p=13
d.plot <- c(0,diff(biolit$plot))
biolit$c.npp.kgha[which(d.plot !=0)] <-NA
biolit$c.npp.kgha[which(biolit$c.npp.kgha < 0)] <-NA

stands.biolit <- merge(stands.wmean[,c("plot","species","age")],biolit,by=c("plot","species"))
names(stands.biolit)
#equilibrium stem/biomass litter (as mean of annual increment)
stands.biolit$needles.eq.lit <- stands.biolit$needles/stands.biolit$age
stands.biolit$branch.eq.lit <- stands.biolit$branch/stands.biolit$age
stands.biolit$stembark.eq.lit <- stands.biolit$stembark/stands.biolit$age
stands.biolit$stump.eq.lit <- stands.biolit$stump/stands.biolit$age
stands.biolit$coarseroots.eq.lit <- stands.biolit$coarseroots/stands.biolit$age
stands.biolit$fineroots.eq.lit <- stands.biolit$fineroots/stands.biolit$age

d.plot <- c(0,diff(biolit$plot))
#biolit$eq.litter[which(d.plot !=0)] <-NA
#biolit$eq.litter[which(biolit$eq.litter< 0)] <-NA

#View(biolit)
#View(stands.wmean)
#View(stands.biolit)

#natural forest stem litter input (equilibrium forest stem litter = mean increment)
#equilibrium stem/biomass litter (as mean of annual increment)
biolit$eq.stem.litter <- c(NA,diff(biolit$stembark))
d.plot <- c(0,diff(biolit$plot))
biolit$eq.stem.litter[which(d.plot !=0)] <-NA
biolit$eq.stem.litter[which(biolit$eq.stem.litter< 0)] <-NA

biolit$eq.litter <- c(NA,diff(biolit$biomass))
d.plot <- c(0,diff(biolit$plot))
biolit$eq.litter[which(d.plot !=0)] <-NA
biolit$eq.litter[which(biolit$eq.litter< 0)] <-NA
cbind(biolit$stembark,biolit$eq.stem.litter,biolit$eq.litter)

#stope

## BIOMASS LITTER TREESTANDS - OUTPUT1 #UNITS Kg C ha #############################################################################
#save(biolit, file=paste0(path.data,"ECO_biomass.litter.c.RData"))
###################################
tl <-biolit #RENAME TO TREE LITTER#
###################################


#add understory to stands with weighted age and height
stands.wmean$basal <- round(stands.wmean$n.ha*pi*(stands.wmean$dmed_cm/200)^2,2)
stands.wmean.specsum <- summaryBy(basal+age.w0 ~ plot, stands.wmean, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
stands.underst <- stands.wmean.specsum  #rename

## add understory measurements from Tupek 2008 gm2
stands.underst$underst.eco.above.gm2 <- c(181,130,165,52,135,190,162,360,550)
#convert to carbon kg/ha!
stands.underst$underst.eco.above.kgha <- carbon(stands.underst$underst.eco.above.gm2)/1000*10000
stands.underst$underst.eco.below.kgha <- stands.underst$underst.eco.above.kgha*c(0.7,0.7,0.7,0.7,0.5,0.5,0.3,0.3,0.3)
#20% as fine woody
#80% as nonwoody
names(stands.underst)
#[1] "plot"                   "basal.sum"              "age.w0.sum"            
#[4] "underst.tot"            "underst.below"          "underst.above"         
#[7] "underst.eco.above.gm2"  "underst.eco.above.kgha" "underst.eco.below.kgha"
#
understory.biolit <- stands.underst[,c("plot","basal.sum","age.w0.sum","underst.eco.above.gm2","underst.eco.above.kgha", "underst.eco.below.kgha")]
names(understory.biolit) <- c("plot","basal0","age0","ubio.eco.above.gm2","ubio.above.kgha", "ubio.below.kgha")
# add understory litter
lsf.str()[grep("litter",lsf.str())]
# "bryof.litter"         "dwarfshrub.litter"   
#"fineroot.litter"      "fineroot.litter.old"  "fineroot.litter.reg"  "fineroot.litter.tsum"
#  "grass.litter"         "herb.litter"         
# "lichen.litter"  
dwarfshrub.litter <-function(mass, comp) {
  if (comp=='abv'){
    return(0.37*mass) }
  if (comp=='bel'){
    return( 0.08116103*mass) }
  
}
herb.litter <- function(mass,comp){
  if (comp=='abv'){
    return(1*mass) }
  if (comp=='bel'){
    return((1/1.7)*mass) }
}
bryof.litter<-function(mass) {
  return(0.417*mass)
}
#distribute understory biomass
#Forest 0.6 of dwarfshrabs in 0.1 herbs and 0.3 mosses
#Transition 0.4 of dwarfshrabs in 0.1 herbs and 0.5 mosses
#Mire 0.3 of dwarfshrabs in 0.3 herbs/grass and 0.4 mosses
#dwarfshrubs
understory.biolit$dwarfs.abv <- understory.biolit$ubio.above.kgha*c(rep(0.6,4),rep(0.4,3),rep(0.2,2))
understory.biolit$dwarfs.bel <- understory.biolit$ubio.below.kgha*c(rep(0.6,4),rep(0.4,3),rep(0.2,2))
#herbs/grass
understory.biolit$herbs.abv <- understory.biolit$ubio.above.kgha*c(rep(0.1,4),rep(0.1,3),rep(0.3,2))
understory.biolit$herbs.bel <- understory.biolit$ubio.below.kgha*c(rep(0.1,4),rep(0.1,3),rep(0.3,2))
#mosses
understory.biolit$moss <- understory.biolit$ubio.above.kgha*c(rep(0.3,4),rep(0.5,3),rep(0.5,2))

###litter understory 
#dwarfshrubs
understory.biolit$litt.dwarfs.abv <- dwarfshrub.litter(understory.biolit$dwarfs.abv,"abv")
understory.biolit$litt.dwarfs.bel <- dwarfshrub.litter(understory.biolit$dwarfs.abv,"bel")
#herbs/grass
understory.biolit$litt.herbs.abo <- herb.litter(understory.biolit$herbs.abv,"abv")
understory.biolit$litt.herbs.bel <- herb.litter(understory.biolit$herbs.bel,"bel")
#mosses
understory.biolit$litt.moss <- bryof.litter(understory.biolit$moss)
#View(understory.biolit)

## BIOMASS LITTER UNDERSTORY OUTPUT2  #UNITS Kg C ha #############################################################################
#save understory
#save(understory.biolit, file=paste0(path.data,"ECO_understory.biomass.litter.c.kgha.RData"))
###################################################
ul <-understory.biolit #RENAME to understory litter
###################################################

## LITTER AWEN FOR Equilibrium ##############################################################
#load litter functions including AWENS
source(paste0(path.scripts,"litter_updated.r"))
#list AWEN functions
lsf.str()[grep("AWEN",lsf.str())]

names(tl)
names(ul)

View(tl)
View(ul)


#AGGREGATE LITTER site mean LEVEL #include finerrots for equilibrium forest
tl$totlit.rh.eq <- rowSums(tl[,c("lit.needle","lit.branch","lit.stem","lit.stump","lit.coarseroot","lit.fineroot")])
#AGGREGATE LITTER site mean LEVEL #exculde fineroots as they should be cut
tl$totlit.rh <- rowSums(tl[,c("lit.needle","lit.branch","lit.stem","lit.stump","lit.coarseroot")])
tl$lit.stemstump <- tl$lit.stem +  tl$lit.stump

ul$totlit.u <- rowSums(ul[,c("litt.dwarfs.abv","litt.dwarfs.bel","litt.herbs.abo","litt.herbs.bel","litt.moss")])


library(doBy)
#
tl.msd <-summaryBy(lit.needle+lit.branch+lit.stemstump+lit.coarseroot+lit.fineroot+totlit.rh+totlit.rh.eq ~ plot+species, 
                   tl, FUN = function(x) { c(m=mean(x, na.rm=T),
                                             sd = sd(x, na.rm=T)) }) # , se= sd(x,na.rm=T)/sqrt(length(x)))})
tl.m <-summaryBy(lit.needle+lit.branch+lit.stemstump+lit.coarseroot+lit.fineroot+totlit.rh+totlit.rh.eq ~ plot+species, 
                 tl, FUN = function(x) { c(m=mean(x, na.rm=T)) }) 

tl.sum <-summaryBy(lit.needle.m+lit.branch.m+lit.stemstump.m+lit.coarseroot.m+lit.fineroot.m+totlit.rh.m+totlit.rh.sd+totlit.rh.eq.m+totlit.rh.eq.sd ~ plot, 
                   tl.msd, FUN = function(x) { c(sum = sum(x, na.rm=T)) }) 
tul.sum <- tl.sum
tul.sum$totlit.u<- ul$totlit.u 
names(tul.sum)
tul.sum$totlit.lu <- tul.sum$totlit.rh.m.sum+tul.sum$totlit.u
tul.sum$totlit.eq.lu <- tul.sum$totlit.rh.eq.m.sum+tul.sum$totlit.u

#data for boxplots bp
tul.bp <- tul.sum[,c("totlit.u","lit.needle.m.sum","lit.branch.m.sum","lit.stemstump.m.sum",
                     "lit.coarseroot.m.sum")]
#litter for eqilibrium
tul.bp.eq <- tul.sum[,c("totlit.u","lit.needle.m.sum","lit.branch.m.sum","lit.stemstump.m.sum",
                        "lit.coarseroot.m.sum","lit.fineroot.m.sum")]

###############################################################################################################################
#save biomass and litter
save(list = c("ul","tl","tl.msd","tul.sum", "tul.bp", "tul.bp.eq"), 
     file=paste(path.data,"biomass.litter_trees.understory_years.RData",sep=""))
################################################################################################################################

## separate total litter to AWEN fractions 
#
#TREES
tree.woody.awen <- data.frame(branches.AWEN(tl$lit.branch+
                                              tl$lit.coarseroot+
                                              tl$lit.stump))
names(tree.woody.awen)<-c("At","Wt","Et","Nt")
tree.woody.awen$plot <- tl$plot
tree.woody.awen$species <- tl$species
tree.woody.awen$year <- tl$year

#during filed campaing exclude coarse roots which were cut off during collar installation 
#and entered soil as litter-biomass (see below litter cut04)
tree.woody.fc.awen <- data.frame(branches.AWEN(tl$lit.branch+
                                                 tl$lit.stump))
names(tree.woody.fc.awen)<-c("At","Wt","Et","Nt")
tree.woody.fc.awen$plot <- tl$plot
tree.woody.fc.awen$species <- tl$species
tree.woody.fc.awen$year <- tl$year

tree.woody.eq.awen <- data.frame(branches.AWEN(tl$lit.branch+
                                                 tl$lit.coarseroot+
                                                 tl$lit.stump+
                                                 tl$eq.stem.litter))
names(tree.woody.eq.awen)<-c("At","Wt","Et","Nt")
tree.woody.eq.awen$plot <- tl$plot
tree.woody.eq.awen$year <- tl$year

#index 1pine 2spruce 3birch
itl1<-which(tl$species==1)
itl2<-which(tl$species==2)
itl3<-which(tl$species==3)
tree.nonwoody.awen <- data.frame(matrix(NA,dim(tree.woody.awen)[1],4))
names(tree.nonwoody.awen)<-c("At","Wt","Et","Nt")
for(i in 1:3){
  # i = 1
  tree.nonwoody.awen[get(c("itl1","itl2","itl3")[i]),]<-
    foliage.AWEN(tl$lit.needle[get(c("itl1","itl2","itl3")[i])], i) + 
    fineroot.AWEN(tl$lit.fineroot[get(c("itl1","itl2","itl3")[i])], i)
}

tree.nonwoody.awen$plot <- tl$plot
tree.nonwoody.awen$year <- tl$year
str(tree.nonwoody.awen)

## trenching cut of coarseroots and fine roots of trees and understory in JULY 2004
tree.w.coarseroot.cut.awen <- data.frame(branches.AWEN(tl$coarseroots))
names(tree.w.coarseroot.cut.awen)<-c("At","Wt","Et","Nt")
tree.w.coarseroot.cut.awen$plot <- tl$plot
tree.w.coarseroot.cut.awen$species <- tl$species
tree.w.coarseroot.cut.awen$year <- tl$year

tree.nw.fineroot.cut.awen <- data.frame(matrix(NA,dim(tree.woody.awen)[1],4))
names(tree.nw.fineroot.cut.awen)<-c("At","Wt","Et","Nt")
for(i in 1:3){
  # i = 1
  tree.nw.fineroot.cut.awen[get(c("itl1","itl2","itl3")[i]),]<- 
    fineroot.AWEN(tl$fineroots[get(c("itl1","itl2","itl3")[i])], i)
}
names(tree.nw.fineroot.cut.awen)<-c("At","Wt","Et","Nt")
tree.nw.fineroot.cut.awen$plot <- tl$plot
tree.nw.fineroot.cut.awen$species <- tl$species
tree.nw.fineroot.cut.awen$year <- tl$year
# trenched coarseroots and fineroots biomass=litter
tree.w.coarseroot.cut04.awen <- subset(tree.w.coarseroot.cut.awen, year== 2004)
tree.nw.fineroot.cut04.awen <- subset(tree.nw.fineroot.cut.awen, year== 2004)

#sum species awen to plot level per year
library(doBy)

tree.w.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year, tree.woody.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
tree.w.fc.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year, tree.woody.fc.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
tree.nw.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year, tree.nonwoody.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})

#equilibrium #woody awen
tree.w.eq.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year, tree.woody.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
tree.w.eq.awen <- summaryBy(At.sum+Wt.sum+Et.sum+Nt.sum ~ plot, tree.w.eq.awen.sum, FUN = function(x) { c(mean = mean(x, na.rm=TRUE))})
#equilibrium #non-woody awen
tree.nw.eq.awen <-summaryBy(At.sum+Wt.sum+Et.sum+Nt.sum ~ plot, tree.nw.awen.sum, FUN = function(x) { c(mean = mean(x, na.rm=TRUE))})

#increase eq.awen of 3  and 5 site due to logs on ground 3(1dead) 5(2 dead)
tree.w.eq.awen[c(3,5),c(2:5)] <- tree.w.eq.awen[c(3,5),c(2:5)]*c(1.3,1.5)
tree.nw.eq.awen[c(3,5),c(2:5)] <- tree.nw.eq.awen[c(3,5),c(2:5)]*c(1.3,1.5)

#remove sum from names
#names(tree.nw.awen.sum ) <-sub(".sum","",names(tree.nw.awen.sum ))

#list AWEN functions
lsf.str()[grep("AWEN",lsf.str())]


#UNDERSTORY
#ul
understory.nw.awen <- data.frame(grass.AWEN(ul$litt.dwarfs.abv+ul$litt.herbs.abo,1)+
                                   grass.AWEN(ul$litt.dwarfs.bel+ul$litt.herbs.bel,0)+ 
                                   moss.AWEN(ul$litt.moss))
names(understory.nw.awen)<-c("Au","Wu","Eu","Nu")
understory.nw.awen$plot <- ul$plot


#trenched understory cut
nw.understory.cut04.awen <- data.frame(grass.AWEN(ul$dwarfs.bel+ul$herbs.bel,0))
names(nw.understory.cut04.awen)<-c("Au","Wu","Eu","Nu")

#LITTER AWEN TREES + UNDERSTORY 
#equilibrium forest
#annual level no seasonal distribution!
treeunder.w.eq.awen <- tree.w.eq.awen
treeunder.nw.eq.awen <- tree.nw.eq.awen[,2:5]+understory.nw.awen[,1:4]
treeunder.nw.eq.awen$plot <- ul$plot
names(treeunder.nw.eq.awen)
treeunder.nw.eq.awen <- treeunder.nw.eq.awen[,c( "plot","At.sum.mean","Wt.sum.mean","Et.sum.mean","Nt.sum.mean")]

#LITTER AWEN TREES + UNDERSTORY BIOMASS input after TRENCHING for time series
# woody tree coarse roots

tree.w.coarseroot.cut04.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year, tree.w.coarseroot.cut04.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
tree.nw.fineroot.cut04.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year, tree.nw.fineroot.cut04.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})

#add fineroots of trees and understory
treeunder.nw.fineroot.cut04.awen <- tree.nw.fineroot.cut04.awen.sum
treeunder.nw.fineroot.cut04.awen[,3:6] <- tree.nw.fineroot.cut04.awen.sum[,3:6]+nw.understory.cut04.awen[,1:4]

#saved litter from trenching
litter.trenchcut04 <- list()
litter.trenchcut04[[1]] <- tree.w.coarseroot.cut04.awen.sum
litter.trenchcut04[[2]] <- treeunder.nw.fineroot.cut04.awen
names(litter.trenchcut04) <- c("tree.w.coarseroot.cut04.awen.sum", "treeunder.nw.fineroot.cut04.awen")
#########################################################################################################################
save(litter.trenchcut04, 
     file=paste0(path.data,"litter.trench.cut_coarse.fineroots.RData"))
#########################################################################################################################

## LITTER AWEN FOR TIME SERIES ##############################################################
#LITTER AWEN TREES + UNDERSTORY 
#time series 2004,2005,2006

#seasonally distributed fineroots,  foliage, understory
#foliage RATIO (applied in Tupek EJSS 2019) saved from icp2.litter.meteo.model.input_07.02.20.r
fol.ratio <- read.csv(paste(path.data,"foliar.monthly.ratios_spruce.pine.csv" , sep=""),
                      header=T, sep = ",",stringsAsFactors = FALSE, na.strings=c("NA",""))
fol.ratio$birch <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.05, 0.10, 0.16, 0.35, 0.29, 0.05, 0.00)
sum(fol.ratio$birch)
par(mfrow=c(1,1))
matplot(fol.ratio, type="l", xlab="months")
legend("topleft", c("spruce","pine","birch"),
       lty=c(1,2,3),lwd=1,col=c(1,2,3),bty="n")

#index 1pine 2spruce 3birch
itl1<-which(tl$species==1)
itl2<-which(tl$species==2)
itl3<-which(tl$species==3)
tree.nonwoody.fol.awen <- data.frame(matrix(NA,dim(tree.woody.awen)[1],4))
names(tree.nonwoody.fol.awen)<-c("At","Wt","Et","Nt")
for(i in 1:3){ # 1pine 2spruce 3birch
  # i = 1
  tree.nonwoody.fol.awen[get(c("itl1","itl2","itl3")[i]),]<-
    foliage.AWEN(tl$lit.needle[get(c("itl1","itl2","itl3")[i])], i)
}
tree.nonwoody.fol.awen$plot <- tl$plot
tree.nonwoody.fol.awen$species <- tl$species
tree.nonwoody.fol.awen$year <- tl$year

#OPTIONAL/Don't use in case of Rh (roots were trenched, Zhiyansku_2014 90% of fineroots in organic and OM layer
#fineroots in first round (fineroots were cut in July 2004) test only against 2005, 2006 data?
tree.nonwoody.fir.awen <- data.frame(matrix(NA,dim(tree.woody.awen)[1],4))
names(tree.nonwoody.fir.awen)<-c("At","Wt","Et","Nt")
for(i in 1:3){
  # i = 1
  tree.nonwoody.fir.awen[get(c("itl1","itl2","itl3")[i]),]<-
    fineroot.AWEN(tl$lit.fineroot[get(c("itl1","itl2","itl3")[i])], i)
}
tree.nonwoody.fir.awen$plot <- tl$plot
tree.nonwoody.fir.awen$species <- tl$species
tree.nonwoody.fir.awen$year <- tl$year
#fineroot RATIO #Zhiyansku_2014 Seasonal-dynamics-of-fine-root-biomass-in-selected-forest-stands.pdf

#apply foliar Modified birch #fol.ratio$birch #
understory.seas.ratio <- c(0, 0,0, 0.04, 0.10, 0.15, 0.21, 0.20, 0.15, 0.10, 0.05, 0.00)
fol.ratio$birch
sum(understory.seas.ratio)
plot(understory.seas.ratio, type = "l")

understory.nw.awen

#replicate UNDERSTORY litter awens to monthly time step
library(dplyr)
#replicate
understory.nw.awen$species <- 0
understory.nw.awen$year <- NA
for(j in 1:4){ #years
  for(i in 1:9){
    #i = 1
    rep.understory.nw.awen.i <- bind_rows(replicate(12, understory.nw.awen[i,], simplify = FALSE))
    
    rep.understory.nw.awen.i[,1:4]<-rep.understory.nw.awen.i[,1:4]*understory.seas.ratio
    rep.understory.nw.awen.i$year <- c(2003:2006)[j] #include 2003 to match data in the other tables
    rep.understory.nw.awen.i$month<- 1:12
    
    
    if(i == 1 & j == 1){
      understory.month.awen <- rep.understory.nw.awen.i
    }else{
      understory.month.awen <- rbind(understory.month.awen,rep.understory.nw.awen.i)
    }
  }
}

#replicate TREES litter awens to monthly time step
for(i in 1:dim(tree.nonwoody.fol.awen)[1]){
  rep.woody.awen.i <- bind_rows(replicate(12, tree.woody.awen[i,], simplify = FALSE))
  rep.woody.fc.awen.i <- bind_rows(replicate(12, tree.woody.fc.awen[i,], simplify = FALSE)) #fc field campaign cutting roots
  rep.fol.awen.i <- bind_rows(replicate(12, tree.nonwoody.fol.awen[i,], simplify = FALSE))
  rep.fir.awen.i <- bind_rows(replicate(12, tree.nonwoody.fir.awen[i,], simplify = FALSE))
  
  #foliage distribute litter to months by species specific fractions
  spec.ni <- unique(rep.fol.awen.i$species)
  rep.fol.awen.i[,1:4]<-rep.fol.awen.i[,1:4]*(as.numeric(fol.ratio[,spec.ni]))
  rep.fol.awen.i$month<-1:12
  #fineroot distribute litter to months evenly by 1/12 
  #(should be also seasonal Zhiyansku_2014 but fineroots were cut in July 2004) 
  # Calibrate only against 2005, 2006 data?
  rep.woody.awen.i[,1:4]<-rep.woody.awen.i[,1:4]*1/12
  rep.woody.awen.i$month<-1:12
  rep.woody.fc.awen.i[,1:4]<-rep.woody.fc.awen.i[,1:4]*1/12 #fc field campaign cutting roots
  rep.woody.fc.awen.i$month<-1:12
  rep.fir.awen.i[,1:4]<-rep.fir.awen.i[,1:4]*1/12
  rep.fir.awen.i$month<-1:12
  
  #str(rep.fol.awen.i)
  if(i == 1){
    tree.woody.month.awen <- rep.woody.awen.i
    tree.woody.fc.month.awen <- rep.woody.fc.awen.i
    tree.nonwoody.fol.month.awen <- rep.fol.awen.i
    tree.nonwoody.fir.month.awen <- rep.fir.awen.i
  }else{
    tree.woody.month.awen <- rbind(tree.woody.month.awen,rep.woody.awen.i)
    tree.woody.fc.month.awen <- rbind(tree.woody.fc.month.awen,rep.woody.fc.awen.i)
    tree.nonwoody.fol.month.awen <- rbind(tree.nonwoody.fol.month.awen,rep.fol.awen.i)
    tree.nonwoody.fir.month.awen <- rbind(tree.nonwoody.fir.month.awen,rep.fir.awen.i)
  }
}

str(tree.woody.month.awen)
str(tree.woody.fc.month.awen)
str(tree.nonwoody.fol.month.awen)
str(tree.nonwoody.fir.month.awen)

str(understory.month.awen)

#sum species to site and year - month level
tree.woody.month.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year+month, tree.woody.month.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
tree.woody.fc.month.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year+month, tree.woody.fc.month.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
tree.nonwoody.fol.month.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year+month, tree.nonwoody.fol.month.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
tree.nonwoody.fir.month.awen.sum <-summaryBy(At+Wt+Et+Nt ~ plot+year+month, tree.nonwoody.fir.month.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})
# understory matching other tables
understory.month.awen.sum <-summaryBy(Au+Wu+Eu+Nu ~ plot+year+month, understory.month.awen, FUN = function(x) { c(sum = sum(x, na.rm=TRUE))})

##save awen for eqilibrium and monthly simulations
#equilibrium
#ls()[grep("eq.awen",ls())]
litter.eq.awen <- list()
litter.eq.awen[[1]] <- treeunder.nw.eq.awen
litter.eq.awen[[2]] <- treeunder.w.eq.awen
names(litter.eq.awen) <- c("treeunder.nw.eq.awen", "treeunder.w.eq.awen")

#View(treeunder.nw.eq.awen)

#monthly
ls()[grep("month.awen",ls())]
#[1] "tree.nonwoody.fir.month.awen" "tree.nonwoody.fol.month.awen" "tree.woody.month.awen"   
litter.natural.month.awen <- list()
litter.natural.month.awen[[1]] <- tree.woody.month.awen.sum
litter.natural.month.awen[[2]] <- tree.nonwoody.fol.month.awen.sum
litter.natural.month.awen[[3]] <- tree.nonwoody.fir.month.awen.sum
litter.natural.month.awen[[4]] <- understory.month.awen.sum 
names(litter.natural.month.awen) <- c("tree.woody.month.awen",
                                      "tree.nonwoody.fol.month.awen",
                                      "tree.nonwoody.fir.month.awen",
                                      "understory.month.awen.sum")


#account for trenching, harvesting - disturbed litter input
#index period before ternching
before.fc.2003.ix <- which(tree.woody.month.awen.sum$year < 2004)
before.fc.2004.ix <- which(tree.woody.month.awen.sum$year > 2003 & tree.woody.month.awen.sum$year < 2005 &
                             tree.woody.month.awen.sum$month < 7)
before.fc.ix <- c(before.fc.2003.ix, before.fc.2004.ix)
#index.july.2004
july.2004.ix <- which(tree.woody.month.awen.sum$year == 2004 &
                        tree.woody.month.awen.sum$month == 7)
#index days after by excluding before and july
after.fc.ix <- -before.fc.ix
j_after.fc.ix <- -c(before.fc.ix ,july.2004.ix) #including july

tree.woody.month.awen.sum[j_after.fc.ix,]

#1) natural monhtly litter to period before JULY 2004 ()
#replace by natural litter
#woody
tree.woody.fc0.month.awen.sum <- tree.woody.fc.month.awen.sum #reduced woody litter branches,setem+stump during all period
tree.woody.fc.month.awen.sum[before.fc.ix,] <- tree.woody.month.awen.sum[before.fc.ix,]
#View(tree.woody.fc.month.awen.sum)
#non-woody
tree.nonwoody.fc.month.fir.awen.sum <- tree.nonwoody.fir.month.awen.sum
tree.nonwoody.fc.month.fir.awen.sum[after.fc.ix,
                                    -match(c("plot","year","month"),
                                           names(tree.nonwoody.fc.month.fir.awen.sum))] <- 0

understory.fc.month.awen.sum <- understory.month.awen.sum
understory.fc.month.awen.sum[after.fc.ix,
                             -match(c("plot","year","month"),names(understory.fc.month.awen.sum))] <- 0

#2) add biomass/cut litter in JULY 2004 to woody litter in only branches+stemstump litter

#coarse roots to woody
tree.woody.fc.month.awen.sum[july.2004.ix,c( "At.sum","Wt.sum","Et.sum","Nt.sum")] <-
  tree.w.coarseroot.cut04.awen.sum[,c( "At.sum","Wt.sum","Et.sum","Nt.sum")] 

#fine-roots to non-woody
tree.nonwoody.fc.month.fir.awen.sum[july.2004.ix,c( "At.sum","Wt.sum","Et.sum","Nt.sum")] <-
  treeunder.nw.fineroot.cut04.awen[,c( "At.sum","Wt.sum","Et.sum","Nt.sum")] 

#3 foliage+branches+stemstump litter after trenching #stays same
#tree.woody.fc.month.awen.sum #branches+stemstump
#tree.nonwoody.fol.month.awen.sum #foliage

#monthly disturbed by filed campaign trenching
ls()[grep("fc",ls())]
litter.fc.month.awen <- list()
litter.fc.month.awen[[1]] <- tree.woody.fc.month.awen.sum
litter.fc.month.awen[[2]] <- tree.nonwoody.fol.month.awen.sum #stays same
litter.fc.month.awen[[3]] <- tree.nonwoody.fc.month.fir.awen.sum
litter.fc.month.awen[[4]] <- understory.fc.month.awen.sum
names(litter.fc.month.awen) <- c("tree.woody.fc.month.awen.sum",
                                 "tree.nonwoody.fol.month.awen.sum",
                                 "tree.nonwoody.fc.month.fir.awen.sum",
                                 "understory.fc.month.awen.sum")
#########################################################################################################################
#save lists of simulated monthly litter awens
#units Kg C ha!!!
save(list= c("litter.eq.awen", "litter.natural.month.awen","litter.fc.month.awen"), 
     file=paste0(path.data,"litter.eco_equilibrium_months.RData"))
#########################################################################################################################

