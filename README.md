## Stata package ivcombine for IV estimation with multiple instruments

This is a public repository for the package ```ivcombine``` for Stata, which implements the estimation procedure for combining multiple instrumental variables in the paper [A Vector Monotonicity Assumption for Multiple Instruments](https://arxiv.org/abs/2009.00553 "Paper").

This is a preliminary version of the code and is offered without warranty. I appreciate any feedback or issues noted as I continue to develop it. A Stata help file is forthcoming; for now, please see below as a guide to usage.

## Installation

The file [ivcombine.ado](ivcombine.ado) is the main code file. The package can be installed directly from within Stata by running
```
net from https://raw.githubusercontent.com/leonardgoff/ivcombine/
net describe ivcombine
net install ivcombine
```

Alternatively, you can install the package by downloading ```ivcombine.ado``` to your ado directory for Stata, e.g. C:\ado\personal. 

Here's some test code to get you going, once ```ivcombine``` is installed. This uses a simulated dataset ```sampledata.dta```, which can also be downloaded from this repository.

```
use "sampledata.dta", clear

## To estimate e.g. the all compliers late (ACL):
ivcombine Y D Z1 Z2, covs(`controls') vary(1 2)

## To estimate e.g. the ``Set LATE'' for instrument 1, among the treated:
ivcombine Y D Z1 Z2, covs(`controls') vary(1) treated
```
