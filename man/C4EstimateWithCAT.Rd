% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/C4EstimateWithCAT.R
\name{C4EstimateWithCAT}
\alias{C4EstimateWithCAT}
\title{Estimate photosynthesis parameters for C4 species using Sharkey's fitting procedure}
\usage{
C4EstimateWithCAT(ACi, Tleaf, Patm, alpha1, x, CaBreakL, CaBreakH, startp,
  TresponseKc, TresponseKo, Tresponsegammastar, TresponseKp, Tresponsegbs,
  TresponseVcmax, TresponseJ, TresponseRd, Tresponsegm, TresponseVpmax,
  TresponseKp2, Tresponsekr, Tresponsekf, TresponseKH, TresponsekCA)
}
\arguments{
\item{ACi}{Gas exchange measurement from Li6400 or other equipment. It is a
dataframe iput. Ci with the unit of ppm. You can prepare the data in Excel file
like the given example and save it as "tab delimited text". Then import data by
ACi <- read.table(file = "/Users/haoranzhou/Desktop/Learn R/ACi curve.txt",
header = TRUE)}

\item{Tleaf}{Leaf temperature when A_Ci curve is measured.}

\item{Patm}{Atmosphere pressure when A_Ci curve is measured.}

\item{alpha1}{The fraction of O2 evolution occurring in the bundle sheath.
Unless you have enough information, input it as the 0.15.}

\item{x}{the fraction of total electron transport that are confined to be
used for the PEP regeneration out of J, which is the total electron transport.}

\item{CaBreakL}{Higher bound of Ci below which A is thought to be controled by
Rubisco Carboxylation (Ac). Start with 10.}

\item{CaBreakH}{Lower bound of Ci above which A is thought to be controled by
RuBP regeneration (Aj). Start with 50. If the estimation results showed
"inadmissible fits", change the CaBreakL and CaBreakH until "inadmissible fits"
disappear.}

\item{startp}{A vector that gives the start points for the estimation
(c(Vcmax,J,Rd,gm,Vpmax,KCA))}

\item{TresponseKc}{A vector that gives the temperature response parameters for the
Kc (c(Kc25,deltaHaKc))}

\item{TresponseKo}{A vector that gives the temperature response parameters for the
Ko (c(Ko25,deltaHaKo))}

\item{Tresponsegammastar}{A vector that gives the temperature response parameters for the
gammastar (c(gammastar25,deltaHagammastar))}

\item{TresponseKp}{A vector that gives the temperature response parameters for the
Kp for CO2 (c(Kp25,deltaHaKp))}

\item{Tresponsegbs}{A vector that gives the temperature response parameters for the
gbs (c(gbs25,deltaHagbs,deltaHdgbs,deltaSgbs))}

\item{TresponseVcmax}{A vector that gives the temperature response parameters for the
Vcmax (c(deltaHaVcmax))}

\item{TresponseJ}{A vector that gives the temperature response parameters for the
J (c(deltaHaJ,deltaHdJ,deltaSJ))}

\item{TresponseRd}{A vector that gives the temperature response parameters for the
Rd (c(deltaHaRd))}

\item{Tresponsegm}{A vector that gives the temperature response parameters for the
gm (c(deltaHagm,deltaHdgm,deltaSgm))}

\item{TresponseVpmax}{A vector that gives the temperature response parameters for the
Vpmax (c(deltaHaVpmax,deltaHdVpmax,deltaSVpmax))}

\item{TresponseKp2}{A vector that gives the temperature response parameters for the
Kp for HCO3- (c(Kp225,deltaHaKp2))}

\item{Tresponsekr}{A vector that gives the temperature response parameters for the
kr (c(kr25,deltaHakr))}

\item{Tresponsekf}{A vector that gives the temperature response parameters for the
kf (c(kf25,deltaHakf))}

\item{TresponseKH}{A vector that gives the temperature response parameters for the
KH (c(KH25))}

\item{TresponsekCA}{A vector that gives the temperature response parameters for the
kCA (c(deltaHakCA,deltaHdkCA,deltaSkCA))}
}
\value{
This package will return a dataframe that contains the following values
(c(Vcmax,J,Rd,gm,Vpmax and KCA)). You can try with c(30, 150, 3, 10, 50, 80).

Parameter at leaf temperature:      A vector (c(Vcmax,J,Rd,gm,Vpmax and
KCA)) returns the estimation parameters at leaf temperature.

Parameter at 25°C:                  A vector (c(Vcmax,J,Rd,gm,Vpmax and
KCA)) returns the estimation parameters at leaf temperature.

Objective:                          The final objective value based on
the estimation results.

Convergence:                        An integer code. 0 indicates
successful convergence.

Message:	                          A character string giving any
additional information returned by the optimizer, or NULL. For details, see PORT
documentation.

Iterations:	                        Number of iterations performed.

Evaluations:	                      Number of objective function and
gradient function evaluations.
}
\description{
Using the gas exchange measurement (A_Ci curve), C4 photosynthesis model with
carbonic anhydrase reaction and Sharkey et al. (2007) fitting processure to do
nonlinear curve fitting (using nlminb package) for estimating photosynthesis
parameters (Vcmax,J,Rd,gm and Vpmax) for C4 species. The difference
between this method with C4EstimateWithCA is that temperature response
parameters need to be provided by the users. Thus, this method provides the option
to alter temperature response parameters. If only planing to alter several
parameters, not all of them, one can use the other parameters provided by
Table S1 in Zhou et al. (2019) ("Deriving C4 photosynthesis parameters by fitting
intensive A/Ci curves"). Make sure to load the "stats"
package before intstalling and using the "C4Estimation" package.
}
