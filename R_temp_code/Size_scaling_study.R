

ksY=c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, kH = 0.0017)
delta1 = -1.7
delta2 = 0.86
r =  -0.306
WS=seq(from=0, to=20, by=0.1)
kinetic_scalinng = c((1+delta1*WS+delta2*(WS^2))^(r))

plot(WS, kinetic_scalinng, type="l", xlab="Diameter of the residual (cm)", ylab="Scaling of the AWEN kinetic", main="size dependence of kinetic in Yasso07")
abline(v=WS[which.max(kinetic_scalinng)])


plot(kinetic_scalinng, WS, type="l")
