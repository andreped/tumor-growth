// clear existing data and results (if any)
cls
clear

// load data from CSV
. import delimited "/Users/andreped/workspace/tumor-growth/data/merged_longitudinal_data_090123.csv"

// set variable of interest
//stset volume_change_relative

// replace -999 with .a for relevant vectors
replace t2 = . if t2 == -999
replace oedema = . if oedema == -999

//replace t2 = 4 if t2 == .
//replace oedema = 4 if oedema == .

// convert string variable genders to categorical variable
//encode gender, generate(gender_n)
encode patient, generate(patient_n)

// explore data
//graph twoway (line volume follow_up_months, connect(ascending)), by(gender_n) xtitle(Follow-up [Months]) ytitle(Volume)


** Linear regression **
// - using random intercept and mixed-effect modelling
meglm volume follow_up_months initial_volume current_age_years ///
	i.gender_bin i.multifocality i.t2 i.oedema ///
	|| patient:, family(gaussian) vce(robust)
	
estat ic
estat icc





** Quadratic regression **
gen follow_up_months2 = follow_up_months^2
meglm volume follow_up_months follow_up_months2 initial_volume ///
	current_age_years i.gender_bin i.multifocality i.multifocality ///
	i.t2 i.oedema || patient:, family(gaussian) vce(robust)
	
estat ic





** Exponential regression **
gen log_volume = log(volume)
//xtmixed log_volume follow_up_months current_age_years initial_volume || patient:, mle
meglm log_volume follow_up_months initial_volume current_age_years ///
	i.gender_bin i.multifocality i.t2 i.oedema ///
	|| patient:, family(gaussian) vce(robust)

// NOTE: Exp only better when T2 and Oedema is added because missing values in
//  t2 and oedema were removed!!! Hence, why "more" variance was explained!

estat ic
estat icc




** Linear radial growth regression **
// @TODO: Are we introducing the multi-level information correct? That is
//   is it sufficient to write {U[patient_n]} to introduce it as a random
//   intercept to introduce a two-level regression model? Not sure...
gen const_ = log(4 * c(pi) / 3)
menl log_volume = {U[patient_n]} + log(4 * c(pi) / 3) + /// 
	3 * log({r0} + {a0} * follow_up_months + {b1} * initial_volume + ///
	{b2} * current_age_years + {b3} * gender_bin + {b4} * multifocality + ///
	{b5} * t2 + {b6} * oedema)
	
estat ic


// TEST factors with menl
// DOES NOT WORK -> menl y = (({b0}+{U0[Id]})*conc)/(({b1}+ i.Group+ {U1[Id]})+conc)
// WORKS ->         menl y = (({b0}+{U0[Id]})*conc)/( {xb: i.Group U1[Id]}+conc)

** Gompertzian regression **
// @TODO: Adding factors outside of function (at y level), degrades log-likelihood!
//   what's the correct place to put them? It fails to start if I put it next to
//   the dependent variables inside the exp/log stuff
 
menl log_volume = {U0[patient_n]} + {k0} + log({v0} / {k0}) * ///
	exp(-{xb:}), ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years ///
	) initial(k0 4.7993 v0 0.0509) mle

// NOTE: We wish to maximize log-likelihood - hence -1600 is poorer than -750

estat ic
estat icc


// initial values gompertz, found in Python through rpy2: Am/Rd/LCP: 5.2992923110723495 -0.15090776979923248 4.450200080871582


** Plot Gompertzian for different genders **

twoway function 





