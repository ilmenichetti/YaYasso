
# This piece of code defines biomass- and stem volume functions #  16.12.2008 by Aleksi Lehtonen (originally by Mikko + Anett) # ################################################################################
# code was documented and collected in a package by Lorenzo Menichetti, 2023, and modified so that all allometric functions are self-contained and not relying any longer on an external file (called biomass_core.r)

## defining functions

## stem volume [m3] on tree level  =============================== FUNCTIONS
#' Function for calculating the stem volume, pine
#' @param d13rm tree diameter at breast heigh, cm
#' @param ht  tree height, m
#' @return the volme of stem (\eqn{m^3})
#' @references Laasasenaho, J. Taper curve and volume functions for pine, spruce and birch. at (1982).
stem.vol.pine <- function(d13rm, ht) {                                 #[dm3 =l]
  b0 <- 0.036089
  b1 <- 2.01395
  b2 <- 2.07025
  b3 <- -1.07209
  b4 <- 0.99676
  return(b0* (d13rm^b1)* (b4^d13rm)* (ht^b2)* ((ht-1.3)^b3)) }

#' Function for calculating the stem volume, spruce
#' @inherit stem.vol.pine
stem.vol.spruce <- function(d13rm, ht) {
  b0 <- 0.022927
  b1 <- 1.91505
  b2 <- 2.82541
  b3 <- -1.53547
  b4 <- 0.99146
  return(b0* (d13rm^b1)* (b4^d13rm)* (ht^b2)* ((ht-1.3)^b3)) }

#' Function for calculating the stem volume, birch
#' @inherit stem.vol.pine
stem.vol.birch <- function(d13rm, ht) {
  b0 <- 0.011197
  b1 <- 2.10253
  b2 <- 3.98519
  b3 <- -2.65900
  b4 <- 0.98600
  return(b0* (d13rm^b1)* (b4^d13rm)* (ht^b2)* ((ht-1.3)^b3)) } #################here functions for stem volume, including dbh, d6 & h

#' Function for calculating the stem volume, pine
#' @inherit stem.vol.pine
#' @param d6  diameter at 6 m height (cm)
stem.vol.pine.comp <- function(d13rm, ht, d6) {              #[dm3 =l]
  b0 <- 0.268621
  b1 <- -0.0145543
  b2 <- -0.0000478628
  b3 <- 0.000334101
  b4 <- 0.0973148
  b5 <- 0.0440716

  return(b0*d13rm^2+b1*d13rm^2*ht+b2*d13rm^3*ht+b3*d13rm^2*ht^2+b4*(d13rm^2+d13rm*d6+d6^2)+b5*d6^2*(ht-6))

}

#' Function for calculating the stem volume, spruce
#' @inherit stem.vol.pine.comp
stem.vol.spruce.comp <- function(d13rm, ht, d6) {           #[dm3 =l]
  b0 <-  0.208043000
  b1 <- -0.014956700
  b2 <- -0.000114406
  b3 <-  0.000436781
  b4 <- 0.133947000
  b5 <- 0.037459900

  return(b0*d13rm^2+b1*d13rm^2*ht+b2*d13rm^3*ht+b3*d13rm^2*ht^2+b4*(d13rm^2+d13rm*d6+d6^2)+b5*d6^2*(ht-6))

}

#' Function for calculating the stem volume, birch
#' @inherit stem.vol.pine.comp
stem.vol.birch.comp <- function(d13rm, ht, d6) {            #[dm3 =l]
  b0 <- 0.226547000
  b1 <- -0.010469100
  b2 <- -0.000122258
  b3 <- 0.000438033
  b4 <- 0.099162000
  b5 <- 0.033483600

  return(b0*d13rm^2+b1*d13rm^2*ht+b2*d13rm^3*ht+b3*d13rm^2*ht^2+b4*(d13rm^2+d13rm*d6+d6^2)+b5*d6^2*(ht-6))
}



## stem wood density [kg*m-3] on tree level      ================= FUNCTIONS
#' Repola function for wood density, pine
#' @param d13rm tree diameter at breast heigh, (cm). The function converts it to dk=2+1.25*d13rm.
#' @param t13 tree age at breast height (years)
#' @param tsum average annual effective temperature sum (>5 °C)
#' @return bulk density of wood (dimensionless, ratio)
#' @references  Repola, J. & Kukkola, R. O. and M. Biomass functions for Scots pine, Norway spruce and birch in Finland. Working Papers of the Finnish Forest Research Institute (2007).
repola.pine.dens <- function(d13rm,t13,tsum) {                          # Pine
  b0 <- 378.39
  b1 <- -78.829
  b2 <- 0.039
  return(b0 + b1*d13rm/t13 + b2*tsum) }

#' Repola function for wood density, spruce
#' @inherit repola.pine.dens
repola.spruce.dens <- function(d13rm,t13) {                             # Spruce
  b0 <- 442.03
  b1 <- -0.904
  b2 <- -82.695
  dk <- 2+1.25*d13rm
  return(b0 + b1*dk + b2*d13rm/t13)
}

#' Repola function for wood density, spruce
#' @inherit repola.pine.dens
repola.birch.dens <- function(d13rm,t13) {                               # Birch
  b0 <- 431.43
  b1 <- 28.054
  b2 <- -52.203
  dk <- 2+1.25*d13rm
  return(b0 + b1*log(dk) + b2*d13rm/t13) }



## stem biomass [kg] on tree level      =================  SIMPLE FUNCTIONS
# Functions based only on dbh & h
#' Repola function for pine, stem wood, model 1
#' @param d13rm tree diameter at breast heigh, cm (the function converts it to dk=2+1.25*d13rm)
#' @param ht  tree height, m
#' @return The biomass component defined by the function (kg)
#' @references  Repola, J. Biomass equations for birch in Finland. Silva Fenn. 42, (2008) and Repola, J. Biomass equations for Scots pine and Norway spruce in Finland. Silva Fenn. 43, (2009).
repola.pine.stem.simp <- function(d13rm, ht) {
  ma.ru.m1 = function(dk, ht) exp(-3.721+8.103*dk/(dk+14)+5.066*ht/(ht+12)+(0.002+0.009)/2)  #Model 1
  dk <- 2+1.25*d13rm
  return(ma.ru.m1(dk,ht))
}

#' Repola function for spruce, stem wood, model 1
#' @inherit repola.pine.stem.simp
repola.spruce.stem.simp <- function(d13rm, ht) {
  ku.ru.m1 = function(dk, ht) exp(-3.555+8.042*dk/(dk+14)+0.869*log(ht)+0.015*ht+(0.009+0.009)/2)   #Model 1
  dk <- 2+1.25*d13rm
  return(ku.ru.m1(dk,ht))
}
#' Repola function for birch, stem wood, model 1
#' @inherit repola.pine.stem.simp
repola.birch.stem.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ko.ru.m1 = function(dk,ht) exp(-4.879+9.651*dk/(dk+12)+1.012*log(ht)+(0.00263+0.00544)/2)#Model 1
  return(ko.ru.m1(dk,ht))

}


## stem bark biomass [kg] on tree level      =================  SIMPLE FUNCTIONS
# Functions based only on dbh & h
#' Repola function for pine, stem bark, model 2
#' @inherit repola.pine.stem.simp
repola.pine.stembark.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ma.kuo.m2 = function(dk, ht) exp(-4.695+8.727*dk/(dk+12)+0.228*log(ht)+(0.014+0.057)/2)             #Model 2
  return(ma.kuo.m2(dk,ht))

}

#' Repola function for spruce, stem bark, model 2
#' @inherit repola.pine.stem.simp
repola.spruce.stembark.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ku.kuo.m2 = function(dk, h) exp(-4.437+10.071*dk/(dk+18)+0.261*log(h)+(0.019+0.039)/2)#Model 2
  return(ku.kuo.m2(dk,ht))
}

#' Repola function for birch, stem bark, model 2
#' @inherit repola.pine.stem.simp
repola.birch.stembark.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ko.kuo.m2 = function(dk, h) exp(-5.433+10.121*dk/(dk+12)+2.647*h/(h+20)+(0.01059+0.04419)/2)#Model 2
  return(ko.kuo.m2(dk,ht))
}


##  weight of foliage [kg] on tree level    ======================= FUNCTIONS
#' Repola function for pine, needles  model 2
#' @inherit repola.pine.stem.simp
#' @param hlc  Height of living crown (m)
repola.pine.needles <- function(d13rm, ht, hlc) {
  cl <- ht- hlc
  dk <- 2+1.25*d13rm
  ma.neul.m2 = function(dk, h, cl) exp(-1.748+14.824*dk/(dk+4)-12.684*h/(h+1)+1.209*log(cl)+(0.032+0.093)/2)                #Model 2
  return(ma.neul.m2(dk,ht,cl))

}
#' Repola function for spruce, needles  model 2
#' @inherit repola.pine.needles
repola.spruce.needles <- function(d13rm, ht, hlc) {
  cl <- ht- hlc
  dk <- 2+1.25*d13rm
  ku.neul.m2 = function(dk, h, cl) exp(-0.085+15.222*dk/(dk+4)-14.446*h/(h+1)+1.273*log(cl)+(0.028+0.087)/2)    #Model 2
  return(ku.neul.m2(dk,ht,cl))
}

#' Repola function for birch, leaves  model 2
#' @inherit repola.pine.needles
repola.birch.foliage <- function(d13rm, ht, hlc) {
  cl <- ht- hlc
  cr <- cl/ht
  dk <- 2+1.25*d13rm
  ko.lehdet.m2= function(dk, cr) exp(-20.856+22.320*dk/(dk+2)+2.819*cr+(0.01082+0.04355)/2)#Model 2
  return(ko.lehdet.m2(dk,cr))
}


######Marklundß
#' Marklund function for pine, needles
#' @inherit repola.pine.needles
#' @references Marklund L.G., Biomassafunktioner för tall, gran och björk i Sve-rige, Sveriges lantbruksuniversitet, Rapporter-Skog 45 (1988) 1–73
marklund.pine.needles <-  function(d13rm, ht) {
  return(exp(-3.4781 + 12.1095*d13rm/(d13rm+7) -1.565*log(ht) +0.0413*ht ) ) }
marklund.spruce.needles <-  function(d13rm, ht) {
  return(exp(-1.8551 + 9.7809*d13rm/(d13rm+12)-0.4873*log(ht)))
}

#' Marklund function for spruce, needles
#' @inherit marklund.pine.needles
#' @inherit repola.pine.needles
marklund.spruce.needles.comp <-  function(d13rm, ht, hlc) {
  cl <- ht- hlc
  return(exp(-1.5732 + 8.4127*d13rm/(d13rm+12) -1.5628*log(ht) + 1.4032*log(cl))) }
marklund.birch.livbranch <- function(d13rm) {
  return(exp(-3.3633+ 10.2806*d13rm/(d13rm+ 10) )) }




# Functions based only on dbh & hß
#' Repola function for pine, needles  model 1
#' @inherit repola.pine.stem.simp
repola.pine.needles.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ma.neul.m1 = function(dk, h) exp(-6.303+14.472*dk/(dk+6)-3.976*h/(h+1)+(0.109+0.118)/2)             #Model 1
  return(ma.neul.m1(dk,ht))
}

#' Repola function for spruce, needles  model 1ß
#' @inherit repola.pine.stem.simp
repola.spruce.needles.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ku.neul.m1 = function(dk, h) exp(-2.994+12.251*dk/(dk+10)-3.415*h/(h+1)+(0.107+0.089/2))                #Model 1
  return(ku.neul.m1(dk,ht))
}

#' Repola function for birch, leaves  model 1
#' @inherit repola.pine.stem.simp
repola.birch.foliage.simp <- function(d13rm) {
  dk <- 2+1.25*d13rm
  ko.lehdet.m1= function(dk) exp(-29.556+33.372*dk/(dk+2)+(0.077)/2) #Model 1
  return(ko.lehdet.m1(dk))
}


## weight of living branches [kg] on tree level    =============== FUNCTIONS
#' Repola function for pine, living branches,  model 2
#' @inherit repola.pine.needles
repola.pine.livbranch <- function(d13rm, ht, hlc) {
  cl <- ht- hlc
  cr <- cl/ht
  dk <- 2+1.25*d13rm
  ma.ran.m2 = function(dk, h, cl) exp(-5.166+13.085*dk/(dk+12)-5.189*h/(h+8)+1.110*log(cl)+(0.020+0.063)/2)#Model 2
  return(ma.ran.m2(dk,ht,cl))
}

#' Repola function for spruce, living branches,  model 2
#' @inherit repola.pine.needles
repola.spruce.livbranch <- function(d13rm, ht, hlc) {
  cl <- ht- hlc
  cr <- cl/ht
  dk <- 2+1.25*d13rm
  ku.ran.m2 = function(dk, h, cl) exp(-3.023+12.017*dk/(dk+14)-5.722*h/(h+5)+1.033*log(cl)+(0.017+0.068)/2)  #Model 2
  return(ku.ran.m2(dk,ht,cl))
}

#' Repola function for birch, living branches,  model 2
#' @inherit repola.pine.needles
repola.birch.livbranch <- function(d13rm, ht, hlc) {
  cl <- ht- hlc
  cr <- cl/ht
  dk <- 2+1.25*d13rm
  ko.ran.m2= function(dk, h, cl) exp(-5.067+14.614*dk/(dk+12)-5.074*h/(h+12)+0.092*cl+(0.01508+0.05663)/2)      #Model 2
  return(ko.ran.m2(dk,ht,cl))
}



# Functions based only on dbh & h
#' Repola function for pine, living branches,  model 1
#' @inherit repola.pine.needles
repola.pine.livbranch.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  return(ma.ran.m1(dk,ht))
}

#' Repola function for spruce, living branches,  model 1
#' @inherit repola.pine.needles
repola.spruce.livbranch.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  return(ku.ran.m1(dk,ht))
}

#' Repola function for birch, living branches,  model 1
#' @inherit repola.pine.needles
repola.birch.livbranch.simp <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  return(ko.ran.m1(dk,ht))
}



##  aboveground biomass [kg] on tree level   ====================== FUNCTIONS
#' Repola function for pine, aboveground biomass,  model 1
#' @inherit repola.pine.needles
repola.pine.above <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ma.tot.m1 = function(dk, h) exp(-3.198+9.547*dk/(dk+12)+3.241*h/(h+20)+(0.009+0.010)/2)            #Model 1
  return(ma.tot.m1(dk,ht))
}

#' Repola function for spruce, aboveground biomass,  model 1
#' @inherit repola.pine.needles
repola.spruce.above <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ku.tot.m1 = function(dk, h) exp(-1.808+9.482*dk/(dk+20)+0.469*log(h)+(0.006+0.013)/2)           #Model 1
  return(ku.tot.m1(dk,ht))
}

#' Repola function for birch, aboveground biomass, model 1
#' @inherit repola.pine.needles
repola.birch.above <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ko.tot.m1 = function(dk, h) exp(-3.654+10.582*dk/(dk+12)+3.018*h/(h+22)+(0.00068+0.00727)/2) #Model 1
  return(ko.tot.m1(dk,ht))
}




##  belowground biomass [kg] on tree level    ===================== FUNCTIONS
# 1) REPOLA
#' Repola function for pine, stump biomass, model 1
#' @inherit repola.pine.needles
repola.pine.stump <- function(d13rm) {                                 # Pine
  dk <- 2+1.25*d13rm
  ma.kanto.m1 = function(dk) exp(-6.753+12.681*dk/(dk+12)+(0.010+0.044)/2)              #Model 1
  return(ma.kanto.m1(dk))
}

#' Repola function for pine, roots biomass, model 1
#' @inherit repola.pine.needles
repola.pine.roots <- function(d13rm) {
  dk <- 2+1.25*d13rm
  ma.juuret.m1 = function(dk) exp(-5.550+13.408*dk/(dk+15)+0.079/2)   #Model 1
  return(ma.juuret.m1(dk))
}

#' Repola function for spruce, stump biomass, model 1
#' @inherit repola.pine.needles
repola.spruce.stump <- function(d13rm) {                                # Spruce
  dk <- 2+1.25*d13rm
  ku.kanto.m1 = function(dk) exp(-3.964+11.730*dk/(dk+26)+(0.065+0.058)/2)        #Model 1
  return(ku.kanto.m1(dk))
}

#' Repola function for spruce, roots biomass, model 1
#' @inherit repola.pine.needles
repola.spruce.roots <- function(d13rm) {
  dk <- 2+1.25*d13rm
  ku.juuret.m1 = function(dk) exp(-2.294+10.646*dk/(dk+24)+(0.105+0.114)/2)       #Model 1
  #return(ma.juuret.m1(dk)) korjattu 2.10.2012 Aleksi Lehtone
  return(ku.juuret.m1(dk))
}

#' Repola function for birch, stump biomass, model 1
#' @inherit repola.pine.needles
repola.birch.stump <- function(d13rm) {                                 # Birch
  dk <- 2+1.25*d13rm
  ko.kanto.m1 = function(dk) exp(-3.574+11.304*dk/(dk+26)+(0.02154+0.04542)/2)               #Model 1
  return(ko.kanto.m1(dk))
}

#' Repola function for birch, roots biomass, model 1
#' @inherit repola.pine.needles
repola.birch.roots <- function(d13rm, ht) {
  dk <- 2+1.25*d13rm
  ko.juuret.m1 = function(dk, h) exp(-3.223+6.497*dk/(dk+22)+1.033*log(h)+(0.048+0.02677)/2)    #Model 1
  return(ko.juuret.m1(dk, ht))
  }




# 2) PETERSSON

#' Petersson function for pine
#' @description
#' The function is an alternative approach to Repola, and it calculated the below-ground biomass to roots >2 mm. Please note that the unit is different from repola (g instead of kg).
#'
#' @inherit repola.pine.needles
#' @return thbe belowground biomass (>2mm) in g
#' @references Petersson, H. & Ståhl, G. Functions for below-ground biomass of Pinus sylvestris, Picea abies, Betula pendula and Betula pubescens in Sweden. Scand. J. For. Res. 21, 84–93 (2006).
petersson.pine.below <-  function(d13rm) {
  return(exp(3.39014 + 11.06822*(d13rm*10)/((d13rm*10)+113) + ((0.35702)^2)/2)) }
#' Petersson function for spruce
#' @inherit petersson.pine.below
petersson.spruce.below <-  function(d13rm) {
    return(exp(4.52965 + 10.57571*(d13rm*10)/((d13rm*10)+142) + ((0.31487)^2)/2)) }
#' Petersson function for birch
#' @inherit petersson.pine.below
petersson.birch.below <-  function(d13rm) {
      return(exp(4.90864 + 9.91194*(d13rm*10)/((d13rm*10)+138) + ((0.44180)^2)/2)) }




##  weight of fine roots [kg] on tree level    ======================= FUNCTIONS

#' Fine root based on root mass
#' @param Mf mass of foliage (any mass unit)
#' @param spec tree species, 1 = Pine, 2 = Spruce, 3 = Birch
#' @references ?????????? not sure at all!!!  Helmisaari, H.-S., Derome, J., Nöjd, P. & Kukkola, M. Fine root biomass in relation to site and stand characteristics in Norway spruce and Scots pine stands. Tree Physiol. 27, 1493–1504 (2007).
fineroot.linear <- function(Mf,spec) {
  if (spec==1) {
    return(Mf*0.676)}
  if (spec==2) {
    return(Mf*0.256)}
  if (spec==3) {
    return(Mf*0.5)}
}

#' Fine root based on root mass and fertility
#' @inherit fineroot.linear
#' @param type Forest type, it is a proxy for soil fertility. 3 = Myrtillus, 4 = Vaccinium, 5 = Calluna
fineroot.linear.fertility <- function(Mf,spec,type) {
  if (spec==1 & type<3) {
    return(Mf*0.2)}
  if (spec==2 & type<3) {
    return(Mf*0.18)}
  if (spec==3 & type<3) {
    return(Mf*1)}

  if (spec==1 & type==3) {
    return(Mf*0.36)}
  if (spec==2 & type==3) {
    return(Mf*0.3)}
  if (spec==3 & type==3) {
    return(Mf*1.5)}

  if (spec==1 & type==4) {
    return(Mf*0.51)}
  if (spec==2 & type==4) {
    return(Mf*0.42)}
  if (spec==3 & type==4) {
    return(Mf*2)}

  if (spec==1 & type>4) {
    return(Mf*0.7)}
  if (spec==2 & type>4) {
    return(Mf*0.54)}
  if (spec==3 & type>4) {
    return(Mf*2.5)}
}


# includes all fineroots, also understorey # tsum from 1961 - 1990 # Helmisaari et al. 2007 # tonnia biomassaa
#' Fine root of the whole stand (including understorey), based on tsum
#' @inherit fineroot.linear
#' @inheritParams repola.pine.dens
fineroot.total.tsum <- function(tsum) {
  return((-0.396*tsum+771.4)/100)
}



#Lethtonen et al. 2015 FORECO
#fine.root model? email from Aleksi 27.2.2015 Subject: fineroot model
#> summary(f3)

#Call:
#  lm(formula = ln_fr ~ ln_ppa + koivu + log(CN), data = fineroot)

#Residuals:
#  Min       1Q   Median       3Q      Max
#-1.09852 -0.32911  0.02857  0.39351  1.01399

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)  3.10397    0.76350   4.065 0.000103 ***
#  ln_ppa       0.73447    0.08698   8.444  5.3e-13 ***
#  koivu        0.66066    0.21455   3.079 0.002760 **
#  log(CN)      0.72086    0.19287   3.738 0.000328 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 0.4488 on 89 degrees of freedom
#(2 observations deleted due to missingness)
#Multiple R-squared:  0.467,  Adjusted R-squared:  0.449
#F-statistic: 25.99 on 3 and 89 DF,  p-value: 3.64e-12




#' Fine root biomass function (basal area)
#' @description
#' Part of a series of functions predicting the fine root biomass of a stand
#'
#' @param basal stand basal area (\eqn{m^2 ha^{-1}})
#' @param decid dominance of birch, factor (0= dominance of conifer, 1 = dominance of birch????)
#' @param CNo carbon:nitrogen ratio of organic layer or upper 0–20 cm peat layer
#' @return fine root biomass (<2mm) in g
#' @references Lehtonen, A. et al. Modelling fine root biomass of boreal tree stands using site and stand variables. For. Ecol. Manag. 359, 361–369 (2016).
fineroot.stand <- function(basal,decid,CNo) {
  return(exp(3.10397 +0.73447*log(basal) + 0.66066*decid + 0.72086*log(CNo) + (0.4488^2)/2))
}


# fineroote model based on stem volume updated 27.2.2015
#' Fine root biomass function (volume)
#' @description
#' Part of a series of functions predicting the fine root biomass of a stand
#'
#' @param vol stand volume (\eqn{m^3 ha^{-1}})
#' @return fine root biomass (<2mm) in g
#' @references Lehtonen, A. et al. Modelling fine root biomass of boreal tree stands using site and stand variables. For. Ecol. Manag. 359, 361–369 (2016).

fineroot.stand.volume <- function(vol) {
  return(exp(6.203+0.320*log(vol)+(0.530^2)/2))
}

fineroot.stand.volume.OLD <- function(vol) {
  return(exp(6.44043+0.27071*log(vol)+(0.542^2)/2))

}
### above coding is from Sanna Harkonen biomassa.r





##########################################
### Understorey vegetation, based data from M.Salemaa ### predictions are kg per ha ##########################################

###is it this ????
# Salemaa, M., Hamberg, L., Kalinauskaite, N., Korpela, L., Lindroos, A-J., Nöjd, P., & Tonteri, T. (2013). Understorey vegetation on level II plots during 1998-2009. In L. Merilä, & S. Jortikka (Eds.), Forest Condition Monitoring in Finland - National Report. The Finnish Forest Research Institute The Finnish Forest Research Institute. http://urn.fi/URN:NBN:fi:metla-201305087568

#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @param cov coverage, % ???? it's probably not % but 0 to 1 ratio!!!!
#' @param comp compartment, abv = above, bel = below ????
#' @param reg region (1 = south 2 = north Finland ????)
#' @references ?????? probably missing, data source could be Salemaa, M., Hamberg, L., Kalinauskaite, N., Korpela, L., Lindroos, A-J., Nöjd, P., & Tonteri, T. (2013). Understorey vegetation on level II plots during 1998-2009. In L. Merilä, & S. Jortikka (Eds.), Forest Condition Monitoring in Finland - National Report. The Finnish Forest Research Institute The Finnish Forest Research Institute. http://urn.fi/URN:NBN:fi:metla-201305087568
understorey.dwarfshrub <- function(cov, comp, reg) {
  if (comp=='abv' & reg==1){
    return(exp(2.555+1.087*log(cov+1)+0.179))}
  if (comp=='abv' & reg==2){
    return(exp(3.494+0.997*log(cov+1)+0.122))}
  if (comp=='bel' & reg==1){
    return(exp(3.377+0.893*log(cov+1)+0.488))}
  if (comp=='bel' & reg==2){
    return(exp(6.177+0.479*log(cov+1)+0.158))}
}

#' Understorey vegetation functions, probably developed by Lehtonen, origin unknown
#' @inheritParams understorey.dwarfshrub
understorey.dwarfshrub.OLD <- function(cov, comp, reg) {
  if (comp=='abv' & reg==1){
    return(cov*19.195)}
  if (comp=='abv' & reg==2){
    return(cov*30.232)}
  if (comp=='bel' & reg==1){
    return(54.531*cov^0.7327)}
  if (comp=='bel' & reg==2){
    return(324.08*cov^0.06074)}
}

# Huom! etelan kerroin molemmille
#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inheritParams understorey.dwarfshrub
understorey.herb <- function(cov, comp, reg) {
  if (comp=='abv' & reg==1){
    return(exp(1.116+1.08*log(cov+1)+0.445)) }
  if (comp=='abv' & reg==2){
    return(exp(1.116+1.08*log(cov+1)+0.445)) }
  if (comp=='bel' & reg==1){
    return(exp(1.87+0.906*log(cov+1)+0.709))}
  if (comp=='bel' & reg==2){
    return(exp(1.87+0.906*log(cov+1)+0.709))}

}

#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inherit understorey.dwarfshrub.OLD
understorey.herb.OLD <- function(cov, comp, reg) {
  if (comp=='abv' & reg==1){
    return(cov*2.5845)}
  if (comp=='abv' & reg==2){
    return(cov*2.5845)}
  if (comp=='bel' & reg==1){
    return(6.9353*cov)}
  if (comp=='bel' & reg==2){
    return(6.9353*cov)}
}

#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inheritParams understorey.dwarfshrub
understorey.grass <- function(cov, comp, reg) {
  if (comp=='abv' & reg==1){
    return(exp(2.75+0.918*log(cov+1)+0.292))}
  if (comp=='abv' & reg==2){
    return(exp(2.75+0.918*log(cov+1)+0.292))}
  if (comp=='bel' & reg==1){
    return(exp(3.368+0.589*log(cov+1)+0.639))}
  if (comp=='bel' & reg==2){
    return(exp(3.368+0.589*log(cov+1)+0.639))}
}

#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inherit understorey.dwarfshrub.OLD
understorey.grass.OLD <- function(cov, comp, reg) {
  if (comp=='abv' & reg==1){
    return(cov*14.128)}
  if (comp=='abv' & reg==2){
    return(cov*13.787)}
  if (comp=='bel' & reg==1){
    return(7.7022*cov)}
  if (comp=='bel' & reg==2){
    return(21.452*cov)}
}

#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inheritParams understorey.dwarfshrub
understorey.bryoph <- function(cov, reg) {
  if (cov < 80) {
    if (reg==1){
      return(exp(3.13+0.795*log(cov+1)+0.281))}
    if (reg==2){
      return(exp(3.13+0.795*log(cov+1)+0.281))}
  }
  if (cov > 79.99) {
    if (reg==1){
      return(1128.517)
    }
    if (reg==2){
      return(1060.17)
    }
  }
}


#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inherit understorey.dwarfshrub.OLD
understorey.bryoph.OLD <- function(cov, reg) {
  if (reg==1){
    return(126.78*exp(cov*0.0229))}
  if (reg==2){
    return(80.079*exp(cov*0.0272))}
}

#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inheritParams understorey.dwarfshrub
understorey.lichen <- function(cov) {
    return(exp(3.69+0.894*log(cov+1)+0.109))
}

#' Understorey vegetation functions, developed by Lehtonen based on data from Salemaa, M., et. al 2013
#' @inherit understorey.dwarfshrub.OLD
understorey.lichen.OLD <- function(cov) {
  return(60.087*cov^0.7015)}

