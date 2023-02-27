// clear existing data and results (if any)
cls
clear

// load data from CSV
. import delimited "/Users/andreped/workspace/tumor-growth/data/merged_longitudinal_data_090123.csv"

// set variable of interest
//stset volume_change_relative

// replace -999 with . for relevant vectors
replace t2 = . if t2 == -999
replace oedema = . if oedema == -999

//replace t2 = 4 if t2 == .
//replace oedema = 4 if oedema == .

// convert string variable genders to categorical variable
//encode gender, generate(gender_n)
encode patient, generate(patient_n)

gen log_volume = log(volume)

// explore data
//graph twoway (line volume follow_up_months, connect(ascending)), by(gender_n) xtitle(Follow-up [Months]) ytitle(Volume)


** Linear regression **
// - using random intercept and mixed-effect modelling
menl volume = {U0[patient_n]} + {xb:}, ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)
	
estat ic




** Quadratic regression **
//gen follow_up_months2 = follow_up_months^2
//meglm volume follow_up_months follow_up_months2 initial_volume ///
//	current_age_years i.gender_bin i.multifocality i.multifocality ///
//	i.t2 i.oedema || patient:, family(gaussian) vce(robust)
	
//estat ic




** Exponential regression **
//xtmixed log_volume follow_up_months current_age_years initial_volume || patient:, mle
menl log_volume = {U0[patient_n]} + {xb:}, ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)

estat ic

// NOTE: Exp only better when T2 and Oedema is added because missing values in
//  t2 and oedema were removed!!! Hence, why "more" variance was explained!




** Linear radial growth regression **
// @TODO: Are we introducing the multi-level information correct? That is
//   is it sufficient to write {U[patient_n]} to introduce it as a random
//   intercept to introduce a two-level regression model? Not sure...

menl log_volume = {U[patient_n]} + log(4 * c(pi) / 3) + 3 * log({xb:}), ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)
	
estat ic




// TEST factors with menl
// DOES NOT WORK -> menl y = (({b0}+{U0[Id]})*conc)/(({b1}+ i.Group+ {U1[Id]})+conc)
// WORKS ->         menl y = (({b0}+{U0[Id]})*conc)/( {xb: i.Group U1[Id]}+conc)

** Gompertzian regression **
// @TODO: Adding factors outside of function (at y level), degrades log-likelihood!
//   what's the correct place to put them? It fails to start if I put it next to
//   the dependent variables inside the exp/log stuff

// Should group variable patient_n be outside function, as a bias term?

/* 
menl log_volume = {k0} + log({v0} / {k0}) * exp(-{xb:}), ///
	define(xb: U0[patient_n] i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years) ///
	initial(k0 4.1993 v0 0.0569)
 */
menl log_volume = {U0[patient_n]} + log({a0}) - exp(-{xb:}), ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)

// NOTE: We wish to maximize log-likelihood - hence -1600 is poorer than -750

estat ic
//estat icc


// initial values gompertz, found in Python through rpy2: Am/Rd/LCP: 5.2992923110723495 -0.15090776979923248 4.450200080871582


** Plot Gompertzian for different genders **

//twoway function 





