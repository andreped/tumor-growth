// clear existing data (if any)
clear

// load data from CSV
. import delimited "/Users/andreped/workspace/tumor-growth/data/merged_processed_regression_data_080123.csv"

// set independent variable
stset volume_change_relative

// perform exponential regression
streg volume_change_relative total_follow_up_months age_at_t1 init_volume, d(exponential)

// perform radial growth regression
//streg volume_change_relative total_follow_up_months age_at_t1 init_volume, d(exponential)

// perform gompertz growth regression
streg volume_change_relative total_follow_up_months age_at_t1 init_volume, d(gompertz)
