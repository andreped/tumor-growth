// clear existing data and results (if any)
cls
clear

// load data from CSV
. import delimited "/Users/andreped/workspace/tumor-growth/data/merged_processed_regression_data_080123.csv"

// set variable of interest
//stset volume_change_relative

// convert string variable genders to categorical variable
encode genders, generate(genders_n)

// perform linear regression
nl (volume_change_relative = {b0} + {b1} * total_follow_up_months + ///
	{b2} * init_volume + {b3} * age_at_t1 + {b4} * genders_n)
	//{b6} * age_at_t1 + {b7} * spacing3 + {b8} * multifocality)
	//{b3} * group(genders))
. estat ic

exit  // used to stop DO file from computing

// perform exponential regression
nl (volume_change_relative = {V0} * exp({alpha} * total_follow_up_months))
. estat ic

// perform radial growth regression
nl (log(volume_change_relative) = log(4 * c(pi) / 3) + 3 * log({r0} + {alpha} * total_follow_up_months))
. estat ic

// perform gompertz growth regression
nl (log(volume_change_relative) = log({K}) + log({V0}/{K}) * exp(-{alpha} * total_follow_up_months))
. estat ic
