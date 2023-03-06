// clear existing data and results (if any)
cls
clear

// load data from CSV
. import delimited "/Users/andreped/workspace/tumor-growth/data/merged_longitudinal_data_090123.csv"

// replace -999 with . for relevant vectors
replace t2 = . if t2 == -999
replace oedema = . if oedema == -999

// convert string variable genders to categorical variable
encode patient, generate(patient_n)

gen log_volume = log(volume)


** Linear regression **
menl volume = {U0[patient_n]} + {xb:}, ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)
	
estat ic


** Exponential regression **
menl log_volume = {U0[patient_n]} + {xb:}, ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)

estat ic

// NOTE: Exp only better when T2 and Oedema is added because missing values in
//  t2 and oedema were removed!!! Hence, why "more" variance was explained!


** Linear radial growth regression **
menl log_volume = {U[patient_n]} + log(4 * c(pi) / 3) + 3 * log({xb:}), ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)
	
estat ic


** Gompertzian regression **
menl log_volume = {U0[patient_n]} + log({a0}) - exp(-{xb:}), ///
	define(xb: i.gender_bin i.multifocality i.oedema i.t2 ///
	follow_up_months initial_volume current_age_years)

estat ic

// NOTE: We wish to maximize log-likelihood - hence -1600 is poorer than -750
// initial values gompertz, found in Python through rpy2: Am/Rd/LCP: 5.2992923110723495 -0.15090776979923248 4.450200080871582
