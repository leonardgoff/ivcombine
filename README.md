## Stata package ivcombine for IV estimation with multiple instruments

This is a public repository for the package ```ivcombine``` for Stata, which implements the estimation procedure for combining multiple instrumental variables described in the paper [A Vector Monotonicity Assumption for Multiple Instruments](https://arxiv.org/abs/2009.00553 "Paper").

This is a preliminary version of the code and is offered without warranty. I appreciate any feedback or issues noted as I continue to develop it. A Stata help file is forthcoming; for now, please see below as a guide to usage.

## Installation

The file [ivcombine.ado](ivcombine.ado) is the main code file. The package can be installed directly from within Stata by running
```
net from https://raw.githubusercontent.com/leonardgoff/ivcombine/master
net describe ivcombine
net install ivcombine
```

Alternatively, you can install the package by downloading ```ivcombine.ado``` to your ado directory for Stata, e.g. C:\ado\personal. 

## Usage

The function delivers a point estimate and standard error for a conditional average treatment effect parameter of the form $\Delta_c = \mathbb{E}[Y_i(1)-Y_i(0)|c(G_i,Z_i)=1]$ where $G_i$ indicates the selection group of unit $i$ and $Z_i$ is their realization of the instruments. The user must choose a function $c$ that indicates a group of generalized ``compliers'' for which $\Delta_c$ is identified. By default, ```ivcombine``` estimates the "all-compliers LATE" (ACLATE), which is the largest group of such compliers. Other choices are described in the paper and can be implemented as below.

Here's some test code to get you going, once ```ivcombine``` is installed. This uses a simulated dataset ```sampledata.dta```, which can also be downloaded from this repository.

```
use "sampledata.dta", clear

*To estimate e.g. the all compliers late on the treated:
ivcombine Y D Z1 Z2, vary(1 2) treated

*To estimate e.g. the ``Set LATE'' for instrument 1:
ivcombine Y D Z1 Z2, vary(1)
```

## Syntax

The basic syntax of the command is as follows:

``` ivcombine outcomevar treatmentvar instrument1 instrument2 ..., [vary(numlist) treated untreated covs(varlist) vce(string)]```

where:

- ```outcomevar``` is an outcome variable (Y, in the notation of the paper)
- ```treatmentvar``` is a binary endogenous treatment variable (D, in the notation of the paper)
- ```instrument1``` is a binary instrumental variable for ```treatmentvar``` (denoted Z<sub>1</sub> in the paper), ```instrument2``` is another binary instrumental variable for `treatmentvar` (denoted Z<sub>2</sub> in the paper),  etc., where all of the instruments are assumed to be jointly valid and satisfy vector monotonicity

Optional arguments to ``ivcombine'' include:

- ```vary``` helps specify the parameter of interest to estimate, by giving a list of instrument indices that are to be varied. If one has three instruments and wants the ACLATE, for example, one would specify ```vary(1 2 3)```. If not set by the user, ```vary``` will default to the set of all instruments.
- ```treated``` helps specify the parameter of interest to estimate, by allowing the user to condition on units that receive treatment. For example, if ```treated``` is set along with ```vary(1)```, the ```ivcombine``` function will return the parameter $SLATT_{\\{1\\}}$, the set LATE among the treated for instrument one.
- ```untreated``` helps specify the parameter of interest to estimate, by allowing the user to condition on units that do not receive treatment. For example, if ```treated``` is set along with ```vary(1)```, the ```ivcombine``` function will return the parameter $SLATU_{\\{1\\}}$, the set LATE among the untreated for instrument one.
- ```covs``` allows a set of control variables (denoted as $X$ in the paper). ```ivcombine``` simply uses these controls as additional regressors when estimating linear projections. This allows the user to recover the desired parameter of interest, assuming the instruments are valid conditional on $X$, under the approximation that the conditional expectation functions $E[Y_i|Z_{i},X_i]$ and $E[D_i|Z_{i},X_i]$ are both additively separable between the vector of instruments $Z_i$ and $X_i$, and linear in the components of $X_i$.
- ```vce``` allows the user to specify a desired variance estimator. Options include ```robust``` and ```cluster```. ```ivcombine``` defaults to heteroskedasticity-robust standard errors, while specifying ```cluster(clustervar)``` would allow cluster-robust standard errors by ```clustervar```.


## Output

In addition to outputting the numerical results to the screen, The ```ivcombine``` function sets the following scalars for use in your code:
- ```param_ptest``` returns the point estimate of the parameter of interest $\Delta_c$
- ```param_se``` returns a standard error for the above estimate of the parameter of interest $\Delta_c$
- ```fstage_ptest``` returns a point estimate of $P(c(G_i,Z_i)=1)$ for the parameter of interest (intuitively, this is the size of the group for whom average treatment effects are being estimated, and serves as the ``denominator'' in $\Delta_c$)
- ```fstage_se``` returns a standard error for the above estimate of $P(c(G_i,Z_i)=1)$
- ```ivcombine``` also stores ```param_ptest``` and the square of ```param_se``` as the Stata matrices ```e(b)``` and ```e(V)``` (each is a $1 \times 1$ matrix). This can be useful when reporting results using e.g. the ```esttab``` command

## Current version of ```ivcombine``` code: 1.1
