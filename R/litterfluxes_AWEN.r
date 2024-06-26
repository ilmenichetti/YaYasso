
#  This piece of code defines litter production rates and the AWEN conversion factors to feed the Yasso decomposition model (5 pools)

#  29.01.2009 by Aleksi Lehtonen

#  update functions 14.10.2019 by Boris Tupek

# new update by Lorenzo Menichetti, Sept-Octo 2023

# updated according to Lehtonen et al 2016 (Tupek et al  2015)

# most functions should be in Silva fennica 2008 annd 2009, Repola

#' Calculates foliage litter
#'
#' @param Mf mass of living foliage (unitless, usually in Mg ha$^{-1}$)
#' @param spec tree species, 1=pine 2=spruce, 3=residues (broadleaves)
#' @param reg Region (1= Northern Finland or 2= Southern Finland)
#' @param min mineral soil, TRUE/FALSE
#' @references Liski, J. et al. Carbon accumulation in Finland’s forests 1922–2004 – an estimate obtained by combination of forest inventory data with modelling of biomass, litter and soil. Ann. For. Sci. 63, 687–697 (2006).
#' @return litter inputs (mass) from a mass of leaves, to be passed to the AWEN partitioning function
#' @seealso \code{\link{foliage.AWEN}}
foliage.litter <- function(Mf, spec, reg, min) {
  if (spec==1 & reg==1 & min==T) {
    return(Mf*0.278)}
  if (spec==1 & reg==2 & min==T) {
    return(Mf*0.213) }
  if (spec==1 & min==F) {
    return(Mf*0.33) }

  if (spec==2) {
    if (reg==1) {
          return(Mf*0.139)}
    if (reg==2) {
          return(Mf*0.1)}
  }
    if (spec==3) {
           return(Mf*0.79)}
}

#' Calculates foliage litter based on the older approach, with one less option (min). Deprecated.
#'
#' @param Mf mass of living foliage (usually in Mg ha$^{-1}$)
#' @param spec tree species, 1=pine 2=spruce, 3=residues (broadleaves)
#' @param reg Region (1= Northern Finland or 2= Southern Finland)
#' @references Liski, J. et al. Carbon accumulation in Finland’s forests 1922–2004 – an estimate obtained by combination of forest inventory data with modelling of biomass, litter and soil. Ann. For. Sci. 63, 687–697 (2006).
#' @return litter inputs (mass) from a mass of leaves, to be passed to the AWEN partitioning function
#' @seealso \code{\link{foliage.AWEN}}
foliage.litter.Liski2006 <- function(Mf, spec, reg) {
  if (spec==1) {
    if (reg==1) {
      return(Mf*0.22)}
    if (reg==2) {
      return(Mf*0.1)}
  }

  if (spec==2) {
    if (reg==1) {
      return(Mf*0.1)}
    if (reg==2) {
      return(Mf*0.05)}
  }
  if (spec==3) {
    return(Mf*0.70)}
}



#' Calculates branch litter
#'
#' @param Mb Mass living branches (usually in Mg ha$^{-1}$)
#' @inheritParams foliage.litter
#' @return litter inputs (mass) from a mass of branches, to be passed to the AWEN partitioning function
#' @seealso \code{\link{branches.AWEN}}
branches.litter <- function(Mb, spec) {
  if (spec==1) {
          return(Mb*0.02)}
  if (spec==2) {
          return(Mb*0.0125)}
  if (spec==3) {
          return(Mb*0.0135)}
}


# the bark of the stump
#' Calculates stump bark litter
#'
#' @param Mst Mass stumps (usually in Mg ha$^{-1}$)
#' @inheritParams foliage.litter
#' @return litter inputs (mass) from a mass of stumps, to be passed to the AWEN partitioning function
#' @seealso \code{\link{stem.AWEN}}
stump.litter <- function(Mst, spec) {
  if (spec==1) {
          return(Mst*0.0029)}
  if (spec==2) {
          return(0)}
  if (spec==3) {
          return(Mst*0.0001)}
}

#' Calculates root litter
#'
#' @param Mr Mass roots (usually in Mg ha$^{-1}$)
#' @inheritParams foliage.litter
#' @return litter inputs (mass) from a mass of roots, to be passed to the AWEN partitioning function
#' @seealso \code{\link{branches.AWEN}}
root.litter <- function(Mr, spec) {
  if (spec==1) {
          return(Mr*0.0184)}
  if (spec==2) {
          return(Mr*0.0125)}
  if (spec==3) {
          return(Mr*0.0135)}
}

#' Calculates bark litter
#'
#' @param Ms Mass bark (usually in Mg ha$^{-1}$)
#' @inheritParams foliage.litter
#' @return litter inputs (mass) from a mass of bark, to be passed to the AWEN partitioning function
#' @seealso \code{\link{stem.AWEN}}
bark.litter <- function(Ms, spec) {
  if (spec==1) {
          return(Ms*0.0052)}
  if (spec==2) {
          return(Ms*0.0027)}
  if (spec==3) {
          return(Ms*0.0029)}
}




#' Calculates fine root litter
#'
#' @param Mfr mass of fine roots
#' @inheritParams foliage.litter
#' @return litter inputs (mass) from a mass of fine roots, to be passed to the AWEN partitioning function
#' @seealso \code{\link{fineroot.AWEN}}
fineroot.litter.reg <- function(Mfr,reg) {
  if (reg==1) {
          return(Mfr*0.85)}
  if (reg==2) {
          return(Mfr*0.5)}
}
#
# Based on assumptions
#

#' Calculates fine root litter with three different options based on tsum
#'
#' @inheritParams fineroot.litter.reg
#' @param tsum temperature sum
#' @return litter inputs (mass) from a mass of fine roots, to be passed to the AWEN partitioning function
#' @seealso \code{\link{fineroot.AWEN}}
fineroot.litter.tsum <- function(Mfr, tsum) {
  if (tsum>=1200) {
          return(Mfr*0.85)}
  if (tsum<=700) {
          return(Mfr*0.5)}
  if (tsum<1200 && tsum>700){
    return(Mfr*(0.0007*tsum+0.001))}
}



#########################
#' Converts biomass in carbon
#'
#' @description
#' The function assumes 50\% C content of organic matter
#'
#' @param mass mass of organic matter
#' @returns a scalar, mass of C
carbon <- function(mass) {
  return(mass*0.5)
}





###############################################3
####################### HERE is a piece to calculate AWEN (acid, water, ethanol & non-solubles) for Yasso input
###########################

# Note for Birch Betula pubenscens and brown leaves is used
#########################
#' Calculates the AWEN proportins of foliage biomass
#'
#' @param Lf leaves biomass
#' @param spec tree species
#' @returns the mass of the four AWEN components
#' @references Yasso07 user-interface manual. J Liski, M Tuomi, J Rasinmäki - , Helsinki, 2009
#' @seealso \code{\link{foliage.litter}}
#' @author Boris Tupek
#'
foliage.AWEN <- function(Lf, spec) {

  fol.AWEN <- matrix(0,nrow=length(Lf), ncol=4)
  ma <- (1:length(Lf))[spec==1]
  ku <- (1:length(Lf))[spec==2]
  ko <- (1:length(Lf))[spec==3]

fol.AWEN[,1][ma] <- 0.518*Lf[ma]
fol.AWEN[,1][ku] <- 0.4826*Lf[ku]
fol.AWEN[,1][ko] <- 0.4079*Lf[ko]

fol.AWEN[,2][ma] <- 0.1773*Lf[ma]
fol.AWEN[,2][ku] <- 0.1317*Lf[ku]
fol.AWEN[,2][ko] <- 0.198*Lf[ko]

fol.AWEN[,3][ma] <- 0.0887*Lf[ma]
fol.AWEN[,3][ku] <- 0.0658*Lf[ku]
fol.AWEN[,3][ko] <- 0.099*Lf[ko]

fol.AWEN[,4][ma] <- 0.216*Lf[ma]
fol.AWEN[,4][ku] <- 0.3199*Lf[ku]
fol.AWEN[,4][ko] <- 0.2951*Lf[ko]

  return(fol.AWEN)
}


#' Calculates the AWEN proportins of fine roots biomass
#'
#' @param Lfr fine roots biomass
#' @param spec tree species
#' @inherit foliage.AWEN
#' @seealso \code{\link{fineroot.litter.reg}} \code{\link{fineroot.litter.tsum}}
#' @author Boris Tupek
fineroot.AWEN <- function(Lfr, spec) {

  fr.AWEN <- matrix(0,nrow=length(Lfr), ncol=4)
  ma <- (1:length(Lfr))[spec==1]
  ku <- (1:length(Lfr))[spec==2]
  ko <- (1:length(Lfr))[spec==3]

fr.AWEN[,1][ma] <- 0.5791*Lfr[ma]
fr.AWEN[,1][ku] <- 0.5508*Lfr[ku]
fr.AWEN[,1][ko] <- 0.4079*Lfr[ko]

fr.AWEN[,2][ma] <- 0.1286*Lfr[ma]
fr.AWEN[,2][ku] <- 0.1331*Lfr[ku]
fr.AWEN[,2][ko] <- 0.198*Lfr[ko]

fr.AWEN[,3][ma] <- 0.0643*Lfr[ma]
fr.AWEN[,3][ku] <- 0.0665*Lfr[ku]
fr.AWEN[,3][ko] <- 0.099*Lfr[ko]

fr.AWEN[,4][ma] <- 0.228*Lfr[ma]
fr.AWEN[,4][ku] <- 0.2496*Lfr[ku]
fr.AWEN[,4][ko] <- 0.2951*Lfr[ko]

  return(fr.AWEN)
}


## Branches are here
# It seems that there is only values for pine (these are applied for others as well)
#' Calculates the AWEN proportions of branches biomass
#'
#' @param Lb branches biomass
#' @inherit foliage.AWEN
#' @seealso \code{\link{branches.litter}}
#' @author Boris Tupek
branches.AWEN <- function(Lb) {
   fb.AWEN <- matrix(0,nrow=length(Lb), ncol=4)

a <- c(0.4763,0.4933,0.4289,0.5068,0.4607,0.5047,0.4642,0.5307,0.5256,0.4661,0.5060,
0.4941,0.4848,0.4158,0.4605,0.4423,0.4811,0.4434,0.5141,0.4312,0.4867,0.3997,0.4758,0.4741,0.4996)

w <- c(0.0196,0.0105,0.0197,0.0120,0.0107,0.0106,0.0130,0.0126,0.0116,0.0195,0.0180,
0.0257,0.0219,0.0295,0.0242,0.0198,0.0242,0.0263,0.0188,0.0218,0.0207,0.0234,0.0176,0.0248,0.0188)

e <- c(0.0870,0.0659,0.1309,0.0506,0.0874,0.0519,0.0840,0.0382,0.0394,0.0996,0.0647,
0.0905,0.0633,0.1131,0.0874,0.1101,0.0681,0.1108,0.0561,0.1128,0.0452,0.1161,0.0678,0.0698,0.0470)

n <- c(0.4170,0.4303,0.4205,0.4306,0.4412,0.4328,0.4388,0.4186,0.4234,0.4148,0.4112,
0.4456,0.4300,0.4416,0.4279,0.4278,0.4266,0.4195,0.4110,0.4341,0.4474,0.4608,0.4388,0.4313,0.4346)

fb.AWEN[,1] <- mean(a)*Lb
fb.AWEN[,2] <- mean(w)*Lb
fb.AWEN[,3] <- mean(e)*Lb
fb.AWEN[,4] <- mean(n)*Lb

return(fb.AWEN)
 }

#' Calculates the AWEN proportions of stems biomass
#'
#' @param Lst stem biomass
#' @param spec tree species
#' @return AWEN proportions of stem biomass
#' @inherit foliage.AWEN
stem.AWEN <- function(Lst, spec) {

  st.AWEN <- matrix(0,nrow=length(Lst), ncol=4)
  ma <- (1:length(Lst))[spec==1]
  ku <- (1:length(Lst))[spec==2]
  ko <- (1:length(Lst))[spec==3]

st.AWEN[,1][ma] <- 0.5*(0.66+0.68)*Lst[ma]
st.AWEN[,1][ku] <- 0.5*(0.63+0.7)*Lst[ku]
st.AWEN[,1][ko] <- 0.5*(0.65+0.78)*Lst[ko]

st.AWEN[,2][ma] <- 0.5*(0.03+0.015)*Lst[ma]
st.AWEN[,2][ku] <- 0.5*(0.03+0.005)*Lst[ku]
st.AWEN[,2][ko] <- 0.5*(0.03+0)*Lst[ko]

st.AWEN[,3][ma] <- 0.5*(0+0.015)*Lst[ma]
st.AWEN[,3][ku] <- 0.5*(0+0.005)*Lst[ku]
st.AWEN[,3][ko] <- 0

st.AWEN[,4][ma] <- 0.5*(0.28+0.29)*Lst[ma]
st.AWEN[,4][ku] <- 0.5*(0.33+0.28)*Lst[ku]
st.AWEN[,4][ko] <- 0.5*(0.32+0.22)*Lst[ko]

  return(st.AWEN)
}

#' Calculates the AWEN proportions of grass litter
#'
#' @param Lg grass biomass
#' @param comp compartment, abv = above, bel = below
#' @inherit foliage.AWEN
#' @inherit understorey.dwarfshrub
#' @seealso \code{\link{understorey.grass.litter}}
#' @author Boris Tupek
understorey.grass.AWEN <- function(Lg, comp) {
  grass.AWEN <- matrix(0, nrow = length(Lg), ncol = 4)

  a <- NULL
  b <- NULL

  if (comp == "abv") {
    a <- 1:length(Lg)
    b <- numeric(0)
  } else if (comp == "bel") {
    b <- 1:length(Lg)
    a <- numeric(0)
  } else {
    stop("Invalid 'above' value. It should be either 0 or 1.")
  }

  grass.AWEN[, 1][a] <- 0.273 * Lg[a]
  grass.AWEN[, 2][a] <- 0.427518 * Lg[a]
  grass.AWEN[, 3][a] <- 0.274482 * Lg[a]
  grass.AWEN[, 4][a] <- 0.025 * Lg[a]

  grass.AWEN[, 1][b] <- 0.273 * Lg[b]
  grass.AWEN[, 2][b] <- 0.506844 * Lg[b]
  grass.AWEN[, 3][b] <- 0.195156 * Lg[b]
  grass.AWEN[, 4][b] <- 0.025 * Lg[b]

  return(grass.AWEN)
}



### Note this is lichen (jäkälä)
#' Calculates the AWEN proportions of lichens biomass
#'
#' @param Ll lichens litter biomass
#' @seealso \code{\link{understorey.lichen.litter}}
understorey.lichen.AWEN <- function(Ll) {
 lichen.AWEN <- matrix(0,nrow=length(Ll), ncol=4)

lichen.AWEN[,1] <- 0.836*Ll
lichen.AWEN[,2] <- 0.080864*Ll
lichen.AWEN[,3] <- 0.031136*Ll
lichen.AWEN[,4] <- 0.052*Ll
return(lichen.AWEN)
}

### Note this is moss (sammal)
#' Calculates the AWEN proportions of bryophytes (was: mosses) biomass
#'
#' @param Lm moss litter biomass
#' @seealso \code{\link{understorey.bryoph.litter}}
#' @author Boris Tupek
understorey.bryoph.AWEN <- function(Lm) {
 moss.AWEN <- matrix(0,nrow=length(Lm), ncol=4)

moss.AWEN[,1] <- 0.736*Lm
moss.AWEN[,2] <- 0.096026*Lm
moss.AWEN[,3] <- 0.036974*Lm
moss.AWEN[,4] <- 0.131*Lm
return(moss.AWEN)
}

#' Calculates the AWEN proportions of wheat biomass
#'
#' @param Lwheat wheat litter biomass
#' @author Boris Tupek
wheat.AWEN <- function(Lwheat) {
  wh.AWEN <- matrix(0,nrow=length(Lwheat), ncol=4)

wh.AWEN[,1] <- 0.712956*Lwheat
wh.AWEN[,2] <- 0.102822*Lwheat
wh.AWEN[,3] <- 0.052463*Lwheat
wh.AWEN[,4] <- 0.131759*Lwheat

  return(wh.AWEN)
}

#' Calculates the AWEN proportions of barley biomass
#'
#' @param Lbarley barley litter biomass
#' @author Boris Tupek
barley.AWEN <- function(Lbarley) {
  bar.AWEN <- matrix(0,nrow=length(Lbarley), ncol=4)

bar.AWEN[,1] <- 0.728503*Lbarley
bar.AWEN[,2] <- 0.133199*Lbarley
bar.AWEN[,3] <- 0.062327*Lbarley
bar.AWEN[,4] <- 0.075972*Lbarley

  return(bar.AWEN)
}

#' Calculates the AWEN proportions of cow dung biomass
#'
#' @param Lshit cow dung biomass
#' @author Boris Tupek
shit.AWEN <- function(Lshit) {
  shi.AWEN <- matrix(0,nrow=length(Lshit), ncol=4)

shi.AWEN[,1] <- 0.701348*Lshit
shi.AWEN[,2] <- 0.102378*Lshit
shi.AWEN[,3] <- 0.067075*Lshit
shi.AWEN[,4] <- 0.129199*Lshit

  return(shi.AWEN)
}


##############################################
# Here adding understorey litter ratios into file
##############################################
#' Calculates the biomass of the understorey dwarf shrubs
#'
#' @param mass mass of shrubs
#' @param comp compartment, abv = above, bel = below
#' @return biomass of dwarf shrubs
#' @seealso \code{\link{understorey.dwarfshrub}}
understorey.dwarfshrub.litter <- function(mass, comp) {
     if (comp=='abv'){
    return(0.37*mass) }
     if (comp=='bel'){
    return( 0.08116103*mass) }

 }

# here calculating bel litter
# Varpujen hienojuurien osuus maanalaisesta kokonaisbiomassasta on keskimäärin 7\%.
# 1/1.7 = 0.5882353
# maavarret 20v 0.05
# tr (93*0.05*0.975+7*0.87*0.5882353) / 100 = 0.08116103
# Hei,
# tässä ehdotus, en tiedä onko tarpeellinen.

# Tutkin Heljä-Siskon kairanäytteitä ICP level2 koealoilta. Totaali maanalainen (jonka mukaan cover biomassa mallit tehtiin) sisältää 7 \% hienojuuria ja 93\% maavarsia. Hienojuurista 13\% on kuolleita ja 87\% eläviä. Maavarsista 2.5 \% on kuolleita. Mitä mieltä olet, pitääkö tämä kuolleiden osuus poistaa aineistosta.
# Maija

#' Calculates the litter input of the understorey bryophytes
#'
#' @param mass mass of bryophytes
#' @return litter inputs (mass) from a mass of bryophytes, to be passed to the AWEN partitioning function
#' @seealso \code{\link{understorey.bryoph}}  \code{\link{understorey.bryoph.AWEN}}
understorey.bryoph.litter <- function(mass) {
    return(0.417*mass)
}

#' Calculates the litter input of the understorey lichens
#'
#' @param mass mass of lichens
#' @return litter inputs (mass) from a mass of lichens, to be passed to the AWEN partitioning function
#' @seealso \code{\link{understorey.lichen}}  \code{\link{understorey.lichen.AWEN}}
understorey.lichen.litter <- function(mass) {
    return(0.1*mass)
}

#' Calculates the litter input of the understorey grasses
#'
#' @param mass mass of grasses
#' @param comp compartment, abv = above, bel = below
#' @return litter inputs (mass) from a mass of grasses, to be passed to the AWEN partitioning function
#' @seealso \code{\link{understorey.grass}}  \code{\link{understorey.grass.AWEN}}
understorey.grass.litter <- function(mass,comp) {
     if (comp=='abv'){
    return(0.33*mass) }
     if (comp=='bel'){
    return((1/1.7)*mass) }
}

#' Calculates the litter input of the understorey herbs
#'
#' @param mass mass of herbs
#' @param comp compartment, abv = above, bel = below
#' @return litter inputs (mass) from a mass of herbs, to be passed to the AWEN partitioning function
#' @seealso \code{\link{understorey.herb}} \code{\link{understorey.grass.AWEN}}
understorey.herb.litter <- function(mass,comp){
     if (comp=='abv'){
    return(1*mass) }
     if (comp=='bel'){
    return((1/1.7)*mass) }
 }



#' Calculates the biomass of the understorey based on \% cover, uplands soils.
#' @description Calculates the biomass of the understorey based on \% cover, uplands soils.
#'
#' @param p_lich coverage of lichens (ratio)
#' @param p_bryoph coverage of bryphites (ratio)
#' @param p_shrubs coverage of shrubs (ratio)
#' @param p_grasses coverage of grasses (ratio)
#' @param spec tree species, 1 = pine 2 = spruce, 3 = residues (broadleaves). For broadleaves the functions are missing, it relies on the same functions than in pine forests
#' @return a vector with biomass of briophytes, lichens, shrubs and grasses, values in g m^-2
#' @references reference Muukkonen, P. et al. Relationship between biomass and percentage cover in understorey vegetation of boreal coniferous forests. Silva Fenn. 40, (2006).
#'
#' @author author Boris Tupek, Lorenzo Menichetti
understorey.biomass <- function(p_lich, p_bryoph, p_shrubs, p_grasses, spec){

  #function introduced by Lorenzo Menichetti, 2023
  eq1=function(p, b_0, b_1) {p^2/(b_0+b_1*p)^2}
  eq2=function(p, b_0) {b_0*p}

      #pine forest (and birch, since it is missing)
      if(spec==1 | spec==3){
        #bryophytes
        b_bri=eq1(p_lich, b_0=1.1833, b_1=0.0334)
        #lychens
        b_lich=eq1(p_bryoph, b_0=4.3369, b_1=0.0128)

        #shrubs
        b_shrubs=eq2(p_shrubs, b_0=2.1262)
        #grasses and herbs
        b_grass=eq2(p_grasses, b_0=0.8416)
      } else if (spec==2){
      #spruce fores
        #bryophytes
        b_bri=eq1(p_lich, b_0=1.8304, b_1=0.0482)
        #lychens
        b_lich=eq1(p_bryoph, b_0=4.3369, b_1=0.0128) #unchanged from pines since it is missing in the paper

        #shrubs
        b_shrubs=eq2(p_shrubs, b_0=1.3169)
        #grasses and herbs
        b_grass=eq2(p_grasses, b_0=0.6552)
      }


return(c(brophytes_biomass=b_bri, lichens_biomass=b_lich, shrubs_biomass=b_shrubs, grasses_herbs_biomass=b_grass))
}


