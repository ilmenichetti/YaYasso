
timesim=seq(1:28)
PR_mm=rep(600, length(timesim))
MT= rep(10, length(timesim))
TA=rep(2, length(timesim))

Litter= cbind(years=timesim, matrix(rep(c(2,2,2,2, 0), each = 28), nrow = 28))
colnames(Litter)=c("years", "A", "W", "E", "N", "O")
str(Litter)

### steady state calculation
AYS.ws <- yasso.matrix.fun(WS = 2, clim = c(5, 500, 7), A.print = "n")
xss=-1*solve(AYS.ws)%*%u.Li #inverse of matrix solve(B)


yassofix <- Yasso07Modelfi(t = years,
                           C0=rep(0,5), #initial carbon
                           AWEN = NULL, #to separate litter to yasso AWEN pools, this depends on plant organ and species
                           In = Litter,#litter C input (same length as years)
                           xi = F, # only xi = 1  will replace climate data no climate effect,
                           MT = MT,#MeanTemperature
                           TA = TA, #TemperatureAmplitude
                           PR_mm = PR_mm,
                           WS=2)


Ct=getC(yassofix)
Rt=getReleaseFlux(yassofix) #respiration





yasso07.SS(WS = 2, clim = c(5, 500, 7), In_ave = c(2, 2, 2, 2, 0))


AYS.ws <- yasso.matrix.fun(WS = 2, clim = c(5, 500, 7), A.print = "n")
xss=-1*solve(AYS.ws)%*%u.Li #inverse of matrix solve(B)
xss
