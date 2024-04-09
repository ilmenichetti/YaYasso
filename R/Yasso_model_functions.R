## Corrected SoilR Yasso07Model #######################################################################
#install.packages("SoilR", dependencies = T)
#library(SoilR)
#Yasso07Model #original version

#the original SoilR version IS MISSINNG FRACTIONATION TO AWEN!!!!!
#as a result C input is not distributed to A W E N pools but all enters pool A

# FIXED SoilR Yasso07model function includes:
# 1) AWEN fractionation
# 2) decomopsition dependency on size of litter
# 3) original environmental functions used for model calibration

# by Boris Tupek boris.tupek@luke.fi
# Dec 2019


#' Yasso07 model as in Tuomi et al. 2011 ecological applications.
#'
#' @description
#' The function relies on the soilR implementation, and it generates a model object that needs then to be solved by soilR functions (see examples)
#'
#' @param t time (years)
#' @param ksY a vector with the decomposition constants of the five Yasso decomposing pools, default c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, kH = 0.0017)
#' @param pY a vector with the transfers (including feedback) between the five
#'   Yasso pools, default  pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99, p5
#'   = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0, p11 = 0.01, p12 = 0.92,
#'   pH = 0.0045),
#' @param beta1  temperature dependence parameter
#' @param beta2 temperature dependence parameter
#' @param gamma  precipitation dependence parameter
#' @param delta1  Woody litter size dependence parameters (cm^-1)
#' @param delta2  Woody litter size dependence parameters (cm^-1)
#' @param r  Woody litter size dependence parameters (adimensional)
#' @param C0  initial C , 5 element vector, e.g. c(0,0,0,0,0)
#' @param In  litter C input, same length as years, if AWEN fractionation not provided it has to be a data.frame (years, 5 AWEN pools), as long as the time steps
#' @param AWEN 5 element vector, fractionation of plant C litter input to yasso AWEN pools
#' @param xi  climate scaling, if xi = FALSE then the model will scale decomposition based on temperature and precipitation data, otherwise it will ignore it (xi = TRUE or xi = 1 is equivalent to setting the xi scaling factor to 1)
#' @param xi_modifier  linear modifier of the climatic effect on decomposition, if = 1 effect is null. Used mostly for model calibration.
#' @param MT Annual mean temperature (vector of length t) (cat('\eqn{^\circ}0') C)
#' @param TA  Temperature Amplitude = (mothly temp. range)/2  (cat('\eqn{^\circ}0') C) (vector of length t)
#' @param PR_mm, Annual precipitation_mm (vector of length t)
#' @param WS woody size, for example 0 no effect  for nonwoody, 2 finewoody, 20 coarse woody (cm)
#' @param solver the solver used by deSolve to solve the linear ODE system, called by SoiLR, default is = deSolve.lsoda.wrapper,
#' @param pass = FALSE
#' @references Ťupek, B. et al. Modeling boreal forest&rsquo;s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier. EGUsphere 2023, 1–34 (2023).
#' @author Boris Tupek, boris.tupek@@luke.fi
#' @return
#'
#' @examples
#' years=seq(1:10)
#' Litter= cbind(years, matrix(rep(c(2,2,2,2, 0), each = 10), nrow = 10))
#' PR_mm=rep(200, 10)
#' MT= rep(10,10)
#' TA= rep(2,10)
#'
#' yassofix <- Yasso07.Modelfi(years,
#'                            C0=rep(0,5), #initial carbon
#'                            AWEN = c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
#'                            In=Litter,#litter C input (same length as years)
#'                            xi = 1, # only xi = 1  will replace climate data no climate effect,
#'                            MT=MT,#MeanTemperature
#'                            TA=TA, #TemperatureAmplitude
#'                            PR_mm=PR_mm,#Precipitation_mm)
#'                            WS=2)
#'
#' Ct=getC(yassofix)
#' Rt=getReleaseFlux(yassofix) #respiration
#'
#' @seealso \link{Yasso07.Modelfi.month}
Yasso07.Modelfi <- function(t, #years
                            #decomposition rates
                            ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031,
                                    kH = 0.0017),
                            #transfers and feedbacks
                            pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99,
                                   p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0,
                                   p11 = 0.01, p12 = 0.92, pH = 0.0045),
                            # climate dependence parameters
                            beta1 = 0.096,
                            beta2 = -0.0014,
                            gamma = -1.21,
                            # Woody litter size dependence parameters
                            delta1 = -1.7,
                            delta2 = 0.86,
                            r =  -0.306,
                            C0, #initial C , 5 element vector, e.g. c(0,0,0,0,0)
                            In, # litter C input, data.frame(years, litter), same length as years, if AWEN faractionatoin not provided it has to be in 5 element form for each time step
                            AWEN, #5 element vector, fractionation of plant C litter input to yasso AWEN pools
                            xi = T , # x1 = T will use climate data
                            xi_modifier = 1, # linear modifier of the climatic effect on decomposition
                            MT = NULL ,# MeanTemperature
                            TA = NULL , # TemperatureAmplitude = (mothly temp. range)/2
                            PR_mm = NULL , # Precipitation_mm
                            WS, # woody size, 0 no effect  for nonwoody, 2 finewoody, 20 coarse woody
                            solver = deSolve.lsoda.wrapper,
                            pass = FALSE){
  # AWEN fractionation # see https://en.ilmatieteenlaitos.fi/yasso-download-and-support,
  # Yasso07 user-interface manual (pdf 450 kB) p.13-14
  t_start = min(t)
  t_end = max(t)

  #structural matrix
  Ap = diag(-1, 5, 5)
  Ap[1, 2] = pY[1]
  Ap[1, 3] = pY[2]
  Ap[1, 4] = pY[3]
  Ap[2, 1] = pY[4]
  Ap[2, 3] = pY[5]
  Ap[2, 4] = pY[6]
  Ap[3, 1] = pY[7]
  Ap[3, 2] = pY[8]
  Ap[3, 4] = pY[9]
  Ap[4, 1] = pY[10]
  Ap[4, 2] = pY[11]
  Ap[4, 3] = pY[12]
  Ap[5, 1:4] = pY[13]

  # add Woody litter size dependence to structural matrix
  # WS in cm, e.g. 0 for nonwoody, 2 for finewoody, 20 - for coarsewoody litter
  # the effect of wl size as in Eq. 3.1 in model description
  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),ksY[5])
    #print(ksY.ws)
    A1.ws = abs(diag(ksY.ws))
    AYS.ws = Ap %*% A1.ws
    #print(AYS.ws)
    return(AYS.ws)
  }
  AYS.ws <- AYS.wsfun(WS)

  #LitterInput
  if (length(In[1,]) == 6) {
    LI = as.matrix(In[,2:6]) #first column for years, 2:6 for AWEN
  } else if (length(AWEN) != 5) {
    stop("the AWEN litter fractionation 5 element vector must be provided if C litter input is not fractionated to AWEN")
  } else {
    #fractionate liter input to yasso07 A W E N pools
    LA =  matrix(AWEN , nrow = length(t),
                 ncol = 5, byrow = TRUE)
    LI = LA*as.vector(In[,2])  #first column for years
  }

  inputFluxes = function(t) {
    matrix(nrow = 5, ncol = 1, LI[t,] )
  }
  inputFluxes_tm = SoilR::BoundInFluxes(
    inputFluxes,
    t_start,
    t_end)

  #environmental effect as function for matrix multiplication
  #this can be replaced by any soil TEMP-SWC function
  y07.Efun <- function(MT,TA,PR_mm){
    #MT = Temp
    #TA =TempAmplit
    PR = PR_mm/1000  # conversion from mm to meters

    #seasonal approximation of temperature (if annual mean)
    T1 = MT + 4*TA/pi*(1/sqrt(2) - 1)          # Eq. 2.4 in model description
    T2 = MT - 4*TA/(sqrt(2)*pi)              # Eq. 2.5 in model description
    T3 = MT + 4*TA/pi*(1 - 1/sqrt(2))          # Eq. 2.6 in model description
    T4 = MT + 4*TA/(sqrt(2)*pi)
    TS =  (exp(beta1*cbind(T1,T2,T3,T4) + beta2*(cbind(T1,T2,T3,T4)^2))*(1 - exp(gamma*PR))) * xi_modifier
    TS
    apply(TS,1,mean)
  }
  xiE = function(t){
    y07.Efun(MT[t], TA[t], PR_mm[t])
  }
  #if xi = 1 replace climate by no effect
  if (xi == T) {
    #woody size but NO environtal effect
    #if woody size 0 than no woody size too
    AYS.Ews_t = SoilR::BoundLinDecompOp(
      function(t){AYS.ws},
      t_start,
      t_end
    )
  } else{
    #woody size and environtal effect
    AYS.Ews_t = SoilR::BoundLinDecompOp(
      function(t){xiE(t)*AYS.ws},
      t_start,
      t_end
    )
  }

  #Yasso07 on Carlos Sierra's general model
  Mod = SoilR::GeneralModel(t, A = AYS.Ews_t, ivList = C0, inputFluxes = inputFluxes_tm,
                   solver, pass)
  return(Mod)
}



#' Yasso07 model implemented as a matrix (to solve its steady state)
#'
#' @description
#' The function generates the structural matrix of the model
#'
#' @param clim a vector containing the mean values of MT, TA, PR_mm as in \link{Yasso07.Modelfi}. If clim = 1 ignored. Since we are considering steady states, climate must be constant
#' @param wetlands if wetlands = TRUE, apply 35% reduction of decomposition se Kleinen et al. (2021). default is FALSE
#' @param A.print if A.print = "y", prints the structural matrix, "n" ignored
#'
#' @inherit Yasso07.Modelfi
#'
#' @references
#' Tuomi, M., Rasinmäki, J., Repo, A., Vanhala, P., & Liski, J. (2011). soil carbon model yasso07 graphical user interface.. https://doi.org/10.48550/arxiv.1105.4961
#'
#' Ťupek, B. et al. Modeling boreal forest&rsquo;s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier. EGUsphere 2023, 1–34 (2023). #'
#'
#' Kleinen, T., Gromov, S., Steil, B., & Brovkin, V. (2021). atmospheric methane underestimated in future climate projections. Environmental Research Letters, 16(9), 094006. https://doi.org/10.1088/1748-9326/ac1814
#'
#' @author Boris Tupek, boris.tupek@@luke.fi
#' @return
#'
#' @example
#' ### steady state calculation
#' AYS.ws <- yasso.matrix.fun(WS = 2, clim = c(5, 500, 7), A.print = "n")
#' xss=-1*solve(AYS.ws)%*%u.Li #inverse of matrix solve(B)

yasso.matrix.fun <- function(WS,# 0,2,20 cm, when 0 ignored
                             clim, #MT, TA, PR_mm, if clim = 1 ignored
                             wetlands = FALSE, # if wetlands = "y", apply 35% reduction of decomposition se Kleinen et al. (2021)
                             A.print, #if A.print = "y", prints the structural matrix, "n" ignored
                             #decomposition rates
                             ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031,
                                     kH = 0.0017),
                             #transfers and feedbacks
                             pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99,
                                    p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0,
                                    p11 = 0.01, p12 = 0.92, pH = 0.0045),
                             # climate dependence parameters
                             beta1 = 0.096,
                             beta2 = -0.0014,
                             gamma = -1.21,
                             # Woody litter size dependence parameters
                             delta1 = -1.7,
                             delta2 = 0.86,
                             r =  -0.306
){
  # for below to work
  # read the structure functions to the environment from the model!


  #structural matrix
  Ap = diag(-1, 5, 5)
  Ap[1, 2] = pY[1]
  Ap[1, 3] = pY[2]
  Ap[1, 4] = pY[3]
  Ap[2, 1] = pY[4]
  Ap[2, 3] = pY[5]
  Ap[2, 4] = pY[6]
  Ap[3, 1] = pY[7]
  Ap[3, 2] = pY[8]
  Ap[3, 4] = pY[9]
  Ap[4, 1] = pY[10]
  Ap[4, 2] = pY[11]
  Ap[4, 3] = pY[12]
  Ap[5, 1:4] = pY[13]

  AYS = Ap %*% abs(diag(ksY))

  if(wetlands==TRUE){
    AYS = Ap %*% abs(diag(0.35*ksY))
    return(AYS)
  }

  if(A.print=="y"){
    print("default structural matrix")
    print(AYS)
  }

  # add Woody litter size dependence to structural matrix
  # WS in cm, e.g. 0 for nonwoody, 2 for finewoody, 20 - for coarsewoody litter
  # the effect of wl size as in Eq. 3.1 in model description

  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),ksY[5])
    A1.ws = abs(diag(ksY.ws))
    AYS.ws = Ap %*% A1.ws
    if(A.print=="y"){
      print("structural matrix modified by woody size")
      print(AYS.ws)
    }
    return(AYS.ws)
  }
  AYS.ws <- AYS.wsfun(WS)

  #yasso environmental function
  if(length(clim) ==3){

    MT = clim[1]
    TA = clim[3]
    PR_mm = clim[2]

    y07.Efun <- function(MT,TA,PR_mm){
      #MT = Temp
      #TA =TempAmplit
      PR = PR_mm/1000  # conversion from mm to meters

      #seasonal approximation of temperature (if annual mean)
      T1=MT+4*TA/pi*(1/sqrt(2)-1)          # Eq. 2.4 in model description
      T2=MT-4*TA/(sqrt(2)*pi)              # Eq. 2.5 in model description
      T3=MT+4*TA/pi*(1-1/sqrt(2))          # Eq. 2.6 in model description
      T4=MT+4*TA/(sqrt(2)*pi)
      TS =  exp(beta1*cbind(T1,T2,T3,T4)+beta2*(cbind(T1,T2,T3,T4)^2))*(1-exp(gamma*PR))
      TS
      apply(TS,1,mean)
    }
    xiE <- y07.Efun(MT, TA, PR_mm)
    AYS.wsE <- xiE*AYS.ws
  } else{
    AYS.wsE <- 1*AYS.ws
  }
  if(A.print=="y"){
    print("returns structural matrix modified by woody size and climate")
    print(AYS.wsE)
  }
  return(AYS.wsE)
}



#' Yasso07 steady state calculation
#'
#' @description
#' The function is just a wrapper for the example in \link{yasso.matrix.fun}, to calculate the steady state by solving the matrix model definition
#'
#' @inherit yasso.matrix.fun
#' @inherit Yasso07.Modelfi
#' @param In_ave average AWEN inputs
#' @author Lorenzo Menichetti, lorenzo.menichetti@@luke.fi
#' @return
#'
#' @example
#' yasso07.SS(WS = 2, clim = c(5, 500, 7), In_ave = c(2, 2, 2, 2, 0))
yasso07.SS <- function(WS,# 0,2,20 cm, when 0 ignored
                       clim, #MT,  PR_mm, TA, if clim = 1 ignored
                       wetlands = FALSE, # if wetlands = T, apply 35% reduction of decomposition se Kleinen et al. (2021)
                       In_ave, #average inputs
                       #decomposition rates
                       ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031,
                               kH = 0.0017),
                       #transfers and feedbacks
                       pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99,
                              p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0,
                              p11 = 0.01, p12 = 0.92, pH = 0.0045),
                       # climate dependence parameters
                       beta1 = 0.096,
                       beta2 = -0.0014,
                       gamma = -1.21,
                       # Woody litter size dependence parameters
                       delta1 = -1.7,
                       delta2 = 0.86,
                       r =  -0.306){

  AYS.ws <- yasso.matrix.fun(WS = WS, clim = clim, A.print = "n",
                             ksY = ksY,
                             #transfers and feedbacks
                             pY = pY ,
                             # climate dependence parameters
                             beta1 = beta1,
                             beta2 = beta2,
                             gamma = gamma,
                             # Woody litter size dependence parameters
                             delta1 = delta1,
                             delta2 = delta2,
                             r =  r)


  yasso_ss=-1*solve(AYS.ws)%*%In_ave #inverse of matrix solve(B)

  rownames(yasso_ss) <- c("A", "W", "E", "N", "O")
  colnames(yasso_ss) <- c("Steady State")

  return(yasso_ss)

}




#' Yasso07 model as in Tuomi et al. 2011 ecological applications, modified to run in monthly steps
#'
#' @inherit Yasso07.Modelfi
#' @param t months (1/12 of year)
#' @return
#' @examples
#' @seealso \link{Yasso07.Modelfi}
#' @references
#' Tuomi, M., Rasinmäki, J., Repo, A., Vanhala, P., & Liski, J. (2011). soil carbon model yasso07 graphical user interface.. https://doi.org/10.48550/arxiv.1105.4961
#'
#' Kleinen, T., Gromov, S., Steil, B., & Brovkin, V. (2021). atmospheric methane underestimated in future climate projections. Environmental Research Letters, 16(9), 094006. https://doi.org/10.1088/1748-9326/ac1814
#Yasso07 model as in Tuomi et al. 2011 ecological applications
Yasso07.Modelfi.month <- function(t, #months 1 month 1/12 of year
                                  #decomposition rates
                                  ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031,
                                          kH = 0.0017)/12, #decrease decomposition rates by 12!!!
                                  wetlands = FALSE, # if wetlands = T decrease kSY down to 35% (Goll et al. 2015, Kleinen et al. 2021)
                                  #transfers and feedbacks
                                  pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99,
                                         p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0,
                                         p11 = 0.01, p12 = 0.92, pH = 0.0045),
                                  # climate dependence parameters
                                  beta1 = 0.096,
                                  beta2 = -0.0014,
                                  gamma = -1.21,
                                  # Woody litter size dependence parameters
                                  delta1 = -1.7,
                                  delta2 = 0.86,
                                  r =  -0.306,
                                  C0, #initial C , 5 element vector, e.g. c(0,0,0,0,0)
                                  In, # litter C input, data.frame(years, litter), same length as years, if AWEN faractionatoin not provided it has to be in 5 element form for each time step
                                  AWEN, #5 element vector, fractionation of plant C litter input to yasso AWEN pools
                                  xi = 1, # x1 = T will use climate data, otherwise will ignore climate data, no climate effect,
                                  MT = NULL,# MeanTemperature
                                  TA = NULL, # TemperatureAmplitude = (mothly temp. range)/2
                                  PR_mm = NULL, # Precipitation_mm
                                  WS, # woody size, 0 no effect  for nonwoody, 2 finewoody, 20 coarse woody
                                  solver = deSolve.lsoda.wrapper,
                                  pass = FALSE){
  # AWEN fractionation # see https://en.ilmatieteenlaitos.fi/yasso-download-and-support,
  # Yasso07 user-interface manual (pdf 450 kB) p.13-14
  t_start = min(t)
  t_end = max(t)


  if (wetlands == TRUE) {
    ksY = 0.35*ksY
    #return(ksY)
  }

  #structural matrix
  Ap = diag(-1, 5, 5)
  Ap[1, 2] = pY[1]
  Ap[1, 3] = pY[2]
  Ap[1, 4] = pY[3]
  Ap[2, 1] = pY[4]
  Ap[2, 3] = pY[5]
  Ap[2, 4] = pY[6]
  Ap[3, 1] = pY[7]
  Ap[3, 2] = pY[8]
  Ap[3, 4] = pY[9]
  Ap[4, 1] = pY[10]
  Ap[4, 2] = pY[11]
  Ap[4, 3] = pY[12]
  Ap[5, 1:4] = pY[13]

  # add Woody litter size dependence to structural matrix
  # WS in cm, e.g. 0 for nonwoody, 2 for finewoody, 20 - for coarsewoody litter
  # the effect of wl size as in Eq. 3.1 in model description
  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1 + delta1*WS + delta2*(WS^2))^(r),ksY[5])
    #print(ksY.ws)
    A1.ws = abs(diag(ksY.ws))
    AYS.ws = Ap %*% A1.ws
    #print(AYS.ws)
    return(AYS.ws)
  }
  AYS.ws <- AYS.wsfun(WS)

  #LitterInput
  if (length(In[1,]) == 6) {
    LI = as.matrix(In[,2:6]) #first column for years, 2:6 for AWEN
  } else if (length(AWEN) != 5) {
    stop("the AWEN litter fractionation 5 element vector must be provided if C litter input is not fractionated to AWEN")
  } else {
    #fractionate liter input to yasso07 A W E N pools
    LA =  matrix(AWEN , nrow = length(t),
                 ncol = 5, byrow = TRUE)
    LI = LA*as.vector(In[,2])  #first column for years
  }

  inputFluxes = function(t) {
    matrix(nrow = 5, ncol = 1, LI[t,] )
  }

  inputFluxes_tm = BoundInFluxes(
    inputFluxes,
    t_start,
    t_end)

  #environmental effect as function for matrix multiplication
  #this can be replaced by any soil TEMP-SWC function
  y07.Efun.month <- function(MT, PR_mm){
    #MT = Temp
    PR = 12*PR_mm/1000  # *12 convert monthly precipitation to annual sum # convert  from mm to meters

    #monthly temperature and upscaled monthly precipitation to annual level
    TS =  exp(beta1*MT + beta2*MT^2)*(1 - exp(gamma*PR))
    TS
  }

  xiE = function(t){
    y07.Efun.month(MT[t], PR_mm[t])
  }

  #if xi = 1 replace climate by no effect
  if (xi == T) {
    #woody size but NO environmental effect
    #if woody size 0 than no woody size too
    AYS.Ews_t = BoundLinDecompOp(
      function(t){AYS.ws},
      t_start,
      t_end
    )
  } else{
    #woody size and environtal effect
    AYS.Ews_t = BoundLinDecompOp(
      function(t){xiE(t)*AYS.ws},
      t_start,
      t_end
    )
  }

  #Yasso07 on Carlos Sierra's general model
  Mod = GeneralModel(t, A=AYS.Ews_t, ivList = C0, inputFluxes = inputFluxes_tm,
                   solver, pass)
  return(Mod)
}
