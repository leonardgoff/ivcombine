/*
********************************************************************************
VERSION 1.00 (AUGUST 29, 2022)
DISCLAIMER: THIS CODE SHOULD BE CONSIDERED PRELIMINARY AND IS OFFERED WITHOUT WARRANTY.
I APPRECIATE YOUR COMMENTS AND ANY ISSUES NOTICED, AT LEONARD.GOFF AT UCALGARY DOT CA.
Copyright (C) 2022 Leonard Goff

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
********************************************************************************
*/

*ssc install desmat
*ssc install coefplot

capture program drop ivcombine
program define ivcombine, eclass

timer clear 51
timer on 51

set more off

qui set type float, permanently

*Preparation and checking input
********************************************************************************

	*Check syntax
	syntax varlist [if] [in],[vary(numlist) treated untreated covs(varlist) vce(string)] 
	

	*Limit sample if requested
		marksample touse
		preserve
		tempfile tmpfilesample tmpfileinputs tmpfileoutcomes_b tmpfileTEs tmpfile tmpfile2
		qui keep if `touse'
		
		gen temporder=_n
	
	*Define main variables
		tokenize "`varlist'"
		local w : word count `varlist'
		if `w' < 3 {
			di as error  "{err}{cmd:varlist()} requires an outcome variable name followed by treatment followed by at least one binary instrument"  
			exit 102
		}
		if `w' >= 3 {
			gen outcomevar =`1'
			gen treatmentvar =`2'
			foreach n of numlist 3/`w' {
				local j = `n'-2
				qui gen instrument`j'=``n''
			}
			local numZ = `w'-2
		}
		
		gen Zprod0 = 1
		gen Zprod1 = 1
		foreach j of numlist 2/`numZ' {
			foreach v of varlist Zprod* {
				gen `v'0 = 1
				gen `v'1 = 1
				drop `v'
			}
		}
		
		foreach v of varlist Zprod* {	
			foreach j of numlist 1/`numZ' {
				if(substr("`v'",5+`j',1)=="1"){
					qui replace `v'=`v'*instrument`j'
				}
			}
		}
	
	*Load standard-errors style
		local w : word count `vce'
		if `w' == 0 {
			local clusterstring = "robust"
		}
		if `w' > 0 {
			local clusterstring = "`vce'"
		}
	
	*Read in covariates
		local w : word count `covs'
		if `w' == 0 {
			local gotcovs = 0
			local ncovlevels=0
		}
		else if `w' >= 1{
			local gotcovs = 1
			local controls = "`covs'"
		}
		
	*Define parameter of interest
	
		*Check if user wants to condition on treated or untreated
		local paramtreated = 0
		if "`treated'"=="treated" {
			local paramtreated = 1
		}
		if "`untreated'"=="untreated" {
			local paramtreated = -1
		}
		if ("`treated'"=="treated") & ("`untreated'"=="untreated") {
			di as error "Cannot set both treated and untreated"
			exit 102
		}
		
		*Load in script J set
		if "`vary'"=="" {
			matrix scriptJ = J(`numZ', 1, 1)
			if(`paramtreated'==0){
				di as result "Parameter of interest: vary() is not set, estimating ACL"
			}
			if(`paramtreated'==1){
				di as result "Parameter of interest: vary() is not set, estimating ACL on the treated"
			}
			if(`paramtreated'==-1){
				di as result "Parameter of interest: vary() is not set, estimating ACL on the untreated"
			}
		}
		else{
			matrix scriptJ = J(`numZ', 1, 0)
			
			foreach j of numlist `vary' {
				matrix scriptJ[`j',1] = 1
			}
			
			if(`paramtreated'==0){
				di as result "Parameter of interest: SLATE on set {`vary'}"
			}
			if(`paramtreated'==1){
				di as result "Parameter of interest: SLATT on set {`vary'}"
			}
			if(`paramtreated'==-1){
				di as result "Parameter of interest: SLATU on set {`vary'}"
			}
		}

	estimates clear
	qui reg outcomevar Zprod* `controls', nocons
	estimates store yvar
	qui reg treatmentvar Zprod* `controls', nocons
	estimates store dvar
	
	cap gen ones=1
	foreach v of varlist Zprod* {	
		qui reg `v' ones, nocons
		estimates store avg_`v'
	}
			
	local firstterm = 1
	
	foreach v of varlist Zprod* {
		local allzeroes = 1
		foreach j of numlist 1/`numZ' {
			if(substr("`v'",5+`j',1)=="1"){
				local allzeroes = 0
			}
		}
		
		*Drop avg_Zprod0000 from memory so it's not included in avg_* in suest
		if `allzeroes'==1{
			estimates drop avg_`v'
		}
				
		*First: check for overlap between scriptJ and `v'
			local overlap = 0
			foreach j of numlist 1/`numZ' {
				if(substr("`v'",5+`j',1)=="1"){
					if(scriptJ[`j',1]==1){
						local overlap = 1
					}
				}
			}
			
			if (`overlap'==1) {
				local newV = "Zprod"
				local newVallzeroes = 1
				
				if(`paramtreated'==0){
					foreach j of numlist 1/`numZ' {
						if(scriptJ[`j',1]==1){
							local newV = "`newV'"+"0"
						}
						else{
							local newV = "`newV'"+substr("`v'",5+`j',1)
							if(substr("`v'",5+`j',1)=="1"){
								local newVallzeroes = 0
							}
						}
					}
					
					if `newVallzeroes'==1{
						local cval = "1"
					}
					else {
						local cval = "_b[avg_`newV'_mean:ones]"
					}
				}
				if(`paramtreated'==1){
					local newV = "`v'"
					local newVallzeroes = `allzeroes'
					
					if `newVallzeroes'==1{
						local cval = "1"
					}
					else {
						local cval = "_b[avg_`newV'_mean:ones]"
					}
				}
				if(`paramtreated'==-1){
					foreach j of numlist 1/`numZ' {
						if(scriptJ[`j',1]==1){
							local newV = "`newV'"+"0"
						}
						else{
							local newV = "`newV'"+substr("`v'",5+`j',1)
							if(substr("`v'",5+`j',1)=="1"){
								local newVallzeroes = 0
							}
						}
					}
															
					if `newVallzeroes'==1{
						local cval = "(1-_b[avg_`v'_mean:ones])"
					}
					else {
						local cval = "(_b[avg_`newV'_mean:ones]-_b[avg_`v'_mean:ones])"
					}
				}
			}
			else {
				local cval = "0"
			}
		
		if `firstterm'== 1 {
			local numerator = "`cval'*_b[yvar_mean:`v']"
			local denominator = "`cval'*_b[dvar_mean:`v']"
			local firstterm = 0
		}
		else{
			local numerator = "`numerator'" + "+ `cval'*_b[yvar_mean:`v']"
			local denominator = "`denominator'" + "+ `cval'*_b[dvar_mean:`v']"
		}
	}
	local nlcomstring = "(`numerator')/(`denominator')"
	
	qui suest yvar dvar avg_*, vce(`clusterstring')
		qui nlcom "`nlcomstring'", post
		matrix temp = e(V)
		scalar param_se = sqrt(temp[1,1])
		matrix temp = e(b)
		scalar param_ptest = temp[1,1]
		
	
	qui suest yvar dvar avg_*, vce(`clusterstring')
		qui nlcom "`denominator'", post
		matrix temp = e(V)
		scalar fstage_se = sqrt(temp[1,1])
		matrix temp = e(b)
		scalar fstage_ptest = temp[1,1]
		
	local param_se: di %4.2f scalar(param_se)
	local param_ptest: di %4.2f scalar(param_ptest)
	local fstage_se: di %4.2f scalar(fstage_se)
	local fstage_ptest: di %4.1f scalar(fstage_ptest)*100 "%"	
	
	di as result "...Estimated proportion of population that are compliers for this parameter: `fstage_ptest', standard error `fstage_se'"
	di as result "...Estimated value of target parameter: `param_ptest', standard error `param_se'"
	
	drop outcomevar treatmentvar instrument* Zprod* ones
	
	estimates clear
	ereturn scalar param_ptest = scalar(param_ptest)
	ereturn scalar param_se = scalar(param_se)
	ereturn scalar fstage_ptest = scalar(fstage_ptest)
	ereturn scalar fstage_se = scalar(fstage_se)
	
	timer off 51
	qui timer list 51
	ereturn scalar seconds_taken = r(t51)
end

