% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Yasso_model_functions.R
\name{Yasso07Modelfi}
\alias{Yasso07Modelfi}
\title{Yasso07 model as in Tuomi et al. 2011 ecological applications}
\usage{
Yasso07Modelfi(
  t,
  ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, kH = 0.0017),
  pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99, p5 = 0, p6 = 0.01, p7 = 0, p8 = 0,
    p9 = 0.03, p10 = 0, p11 = 0.01, p12 = 0.92, pH = 0.0045),
  beta1 = 0.096,
  beta2 = -0.0014,
  gamma = -1.21,
  delta1 = -1.7,
  delta2 = 0.86,
  r = -0.306,
  C0,
  In,
  AWEN,
  xi = T,
  MT,
  TA,
  PR_mm,
  WS,
  solver = deSolve.lsoda.wrapper,
  pass = FALSE
)
}
\arguments{
\item{t}{time (years)}

\item{ksY}{a vector with the decomposition constants of the five Yasso decomposing pools, default c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, kH = 0.0017)}

\item{pY}{a vector with the transfers (including feedbacks) between the five
Yasso pools, defailt  pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99, p5
= 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0, p11 = 0.01, p12 = 0.92,
pH = 0.0045),}

\item{beta1}{climate dependence parameter}

\item{beta2}{climate dependence parameter}

\item{delta1}{Woody litter size dependence parameters}

\item{delta2}{Woody litter size dependence parameters}

\item{r}{TO BE DONE}

\item{C0}{initial C , 5 element vector, e.g. c(0,0,0,0,0)}

\item{In}{litter C input, data.frame(years, litter), same length as years, if AWEN faractionatoin not provided it has to be in 5 element form for each time step}

\item{AWEN}{5 element vector, fractionation of plant C litter input to yasso AWEN pools}

\item{xi}{climate scaling, if x=TRUE then the model will scale decompositon based on temperature and precipitation data, otherwise it will ignore it (equivalent to setting the xi scaling factor to 1)}

\item{MT}{Annual mean temperature}

\item{TA}{TemperatureAmplitude = (mothly temp. range)/2}

\item{PR_mm, }{Annual precipitation_mm}

\item{WS}{woody size, 0 no effect  for nonwoody, 2 finewoody, 20 coarse woody}

\item{solver}{the solver used by deSolve to solve the linear ODE system, called by SoiLR, default is = deSolve.lsoda.wrapper,}

\item{pass}{= FALSE}

\item{gammma}{climate dependence parameter}
}
\description{
Yasso07 model as in Tuomi et al. 2011 ecological applications
}
\references{
Ťupek, B. et al. Modeling boreal forest&rsquo;s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier. EGUsphere 2023, 1–34 (2023).
}
\author{
Boris Tupek, boris.tupek@luke.fi
}
