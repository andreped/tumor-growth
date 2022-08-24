import pandas as pd
import os
import numpy as np
from utils import sort_timestamps, remove_surgery_patients, plot_graphs, str2datetime, get_earliest_timestamp
import seaborn as sns
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
from lmfit import Minimizer, Parameters, report_fit


# R-related imports
stats = importr('stats')
base = importr('base')
nlme = importr('nlme')
car = importr('car')
outliers = importr('outliers')
pandas2ri.activate()
R = ro.r


def preprocess(data_path):
    cohort_personal_info = pd.read_csv(os.path.join(data_path, "cohort_personal_info.csv"))
    cohort_volumes_quality = pd.read_csv(os.path.join(data_path, "cohort_volumes_quality-removed_P09016-T6.csv"))
    # interrater_variability_study = pd.read_csv(os.path.join(data_path, "interrater_variability_study.csv"))  # not considered for now
    volumes = pd.read_csv(os.path.join(data_path, "Volumes.csv"))

    # get unique patients
    patients = cohort_personal_info["Patient"]

    # remove all patients that have had surgery
    patients_no_surgery, patient_filter = remove_surgery_patients(patients)

    # filter other datasets by selected patients
    cohort_personal_info_filtered = cohort_personal_info[patient_filter]

    filter_ = [x in patients_no_surgery for x in cohort_volumes_quality["Patient"]]
    filtered_cohort_volumes_quality = cohort_volumes_quality[filter_]

    filter_ = [x in patients_no_surgery for x in volumes["OP_ID"]]
    filtered_volumes = volumes[filter_]

    print("----\nFiltered data:")
    print(cohort_personal_info_filtered)
    print(filtered_cohort_volumes_quality)
    print(filtered_volumes)

    # 1) First assumption (lets sum all fragmented tumors together into one - total tumor volume in patient,
    #  for each time point). Use cohort volumes quality to catch all patients and time points
    data = []
    unique_patients = np.unique(list(filtered_volumes["OP_ID"]))
    # np.random.shuffle(unique_patients)
    iter = 0
    for pat in unique_patients:
        # get all data for current patient
        curr_data = filtered_volumes[filtered_volumes["OP_ID"] == pat]

        # get unique timestamps
        curr_timestamps = curr_data["Timestamp"]
        curr_timestamps = list(curr_timestamps)

        unique_timestamps = np.unique(list(curr_timestamps))
        unique_timestamps = sort_timestamps(unique_timestamps)

        # get earliest timestamp with non-NaN volume
        for t in unique_timestamps:
            tmp = curr_data[curr_data["Timestamp"] == t]
            curr_v = np.array(tmp["Volume"]).astype("float32")
            curr_v = sum(curr_v)
            if not pd.isnull(curr_v):
                break
        init_timestamp = t

        # get date of first timestamp of patient (T1 might not be the earliest!! If T1 is NaN, then T2 might be, etc.
        first_timestamp_date = curr_data[curr_data["Timestamp"] == init_timestamp]["Date"]
        first_timestamp_date = str2datetime(np.asarray(first_timestamp_date)[0])

        # get initial volume size at first scan
        initial_volume = curr_data[curr_data["Timestamp"] == init_timestamp]["Volume"]
        initial_volume = np.asarray(initial_volume)[0]

        # for each time stamp, all clusters and sum these into one value (total tumor amount in ml)
        for timestamp in unique_timestamps:
            # get volume for current timestamp - if mulitple, sum these (total tumor volume)
            times = curr_data[curr_data["Timestamp"] == timestamp]
            tmp = np.nan_to_num(times["Volume"])  # convert 0 to NaN for summation
            curr_volume = sum(list(tmp))

            # get current date for timestamp - if multiple, select the first (should have save date for the same tumors in the same timestamp)
            curr_dates = times["Date"]
            curr_dates = np.asarray(curr_dates)
            curr_date = curr_dates[0]

            # translate current date to datetime format
            curr_date = str2datetime(curr_date)

            data.append([pat, timestamp, initial_volume, first_timestamp_date, curr_date, curr_volume])
            iter += 1
    data = np.array(data)

    # merge this with the cohort volumes quality stuff
    full_data = filtered_cohort_volumes_quality.copy()

    full_data["Volumes"] = data[:, -1].astype("float32")
    full_data["Date"] = data[:, -2]
    full_data["First_Timestamp_Date"] = data[:, -3]
    full_data["Initial_Volume"] = data[:, -4].astype("float32")

    # need to filter NaN volumes on merged data frame
    full_data = full_data[full_data.Volumes != 0]

    # reset indices in full_data to go 0:1:N
    full_data.index = list(range(len(full_data)))

    # add patient characteristics to full data frame
    full_data["Birth_Year"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    full_data["Gender"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    for pat, gender, byear in zip(cohort_personal_info_filtered["Patient"], cohort_personal_info_filtered["Gender"], cohort_personal_info_filtered["Birth_Year"]):
        row_ids = np.where(full_data["Patient"] == pat)[0]
        for r in row_ids:
            byear_new_format = str(byear) + "-07-01"
            byear_new_format = str2datetime(byear_new_format)

            full_data.loc[r, 'Birth_Year'] = byear_new_format
            full_data.loc[r, 'Gender'] = gender

    # remove all occurences where Volumes=0 (not possible -> tumor was not annotated)
    filter_zero_volumes = full_data["Volumes"] != str(0.0)
    full_data_nonzero = full_data[filter_zero_volumes]

    # get current age at scan and add to data frame
    curr_age_at_scan = full_data_nonzero["Date"] - full_data_nonzero["Birth_Year"]
    curr_age_at_scan = curr_age_at_scan.dt.days
    full_data_nonzero["Current_Age"] = curr_age_at_scan.astype(float)

    # get relative difference in days between scans
    relative_difference_between_scans = full_data_nonzero["Date"] - full_data_nonzero["First_Timestamp_Date"]
    relative_difference_between_scans = relative_difference_between_scans.dt.days
    full_data_nonzero["Relative_Days_Difference"] = relative_difference_between_scans.astype("float32")
    full_data_nonzero["Follow_Up_Months"] = relative_difference_between_scans.astype("float32") / 30

    # get relative volume ratios between scans
    relative_volume_ratio = full_data_nonzero["Volumes"] / full_data_nonzero["Initial_Volume"]
    full_data_nonzero["Relative_Volume_Ratio"] = relative_volume_ratio.astype("float32")

    # @TODO: Should normalize the variables in some way to avoid exploding stuff issues
    #   - perhaps age in years (instead of in days) make more sense?

    print(full_data_nonzero)

    # capture outliers and remove them - assuming normality using quantiles
    lower = R.quantile(full_data_nonzero["Relative_Volume_Ratio"], 0.025)[0]
    higher = R.quantile(full_data_nonzero["Relative_Volume_Ratio"], 0.975)[0]

    x = np.array(full_data_nonzero["Relative_Volume_Ratio"])
    filter_ = (x > lower) & (x < higher)
    full_data_nonzero = full_data_nonzero[filter_]

    tmp_df = pd.DataFrame({
        "Relative_Volume_Ratio": full_data_nonzero["Relative_Volume_Ratio"],
        "Follow_Up_Months": full_data_nonzero["Follow_Up_Months"],
        "Gender": full_data_nonzero["Gender"],
    })

    # test linear regression
    model = R.lm('Relative_Volume_Ratio ~ Follow_Up_Months + factor(Gender)', data=tmp_df)
    summary_model = R.summary(model)
    print(summary_model)  # .rx2('coefficients'))
    print("ANOVA:")
    print(R.anova(model))

    #exit()

    # get (x, y) variables for both genders
    gender_study_df_man = tmp_df[tmp_df["Gender"] == "man"]
    gender_study_df_woman = tmp_df[tmp_df["Gender"] == "woman"]

    # make regression curves for each gender
    model_man = R.lm("Relative_Volume_Ratio ~ Follow_Up_Months", data=gender_study_df_man)
    model_woman = R.lm("Relative_Volume_Ratio ~ Follow_Up_Months", data=gender_study_df_woman)
    # summary_man = R.summary(model_man)
    #R.abline(model_man)
    # print(summary_man)
    #exit()


    # get goodness of fit measures (AIC, Mallows Cp, R^2)
    aic_ = stats.AIC(model)[0]
    bic_ = stats.BIC(model)[0]

    print("R^2 | AIC | BIC:", summary_model.rx2("adj.r.squared")[0], aic_, bic_)

    # sns.scatterplot(x=x_filtered, y=y_filtered, data=full_data_nonzero)
    sns.scatterplot(x="Follow_Up_Months", y="Relative_Volume_Ratio", hue="Gender",
                    data=full_data_nonzero, legend=True)
    sns.lineplot(x=[min(full_data_nonzero["Follow_Up_Months"]), max(full_data_nonzero["Follow_Up_Months"])], y=[1, 1],
                 palette="g")
    plt.grid("on")
    plt.tight_layout()

    print(len(model.rx2('fitted.values')))
    print(len(full_data_nonzero["Relative_Days_Difference"]))

    # draw regression curves
    plt.plot(full_data_nonzero["Follow_Up_Months"], model.rx2('fitted.values'), 'k')

    # plot regression curves for both genders
    plt.plot(gender_study_df_man["Follow_Up_Months"], model_man.rx2('fitted.values'), 'b')
    plt.plot(gender_study_df_woman["Follow_Up_Months"], model_woman.rx2('fitted.values'), 'r')

    # show figure
    plt.show()

    exit()

    # apply log to y
    y_log_filtered = np.log(y_filtered)

    # get mean exponential growth for each groups
    y_log_filtered_man = y_log_filtered[gender_filtered == "man"]
    y_log_filtered_woman = y_log_filtered[gender_filtered == "woman"]

    # get coefficients
    mu_y_log_man = np.mean(y_log_filtered_man)
    mu_y_log_woman = np.mean(y_log_filtered_woman)

    print("\nmeans (man/woman):")
    print(mu_y_log_man)
    print(mu_y_log_woman)

    # print(gender_study_df.head(60))
    print(full_data)
    exit()
    exit()

    #sns.scatterplot(x=x_filtered, y=y_log_filtered, hue=gender_filtered, legend=True)
    #plt.show()
    # -> looks like white noise (but appears to be some trends). Should apply

    #sns.lineplot(x=[min(x_filtered), max(x_filtered)], y=[mu_y_log_man, mu_y_log_man])
    #sns.lineplot(x=[min(x_filtered), max(x_filtered)], y=[mu_y_log_woman, mu_y_log_woman])

    # @FIXME: Probably interaction between "Relative_Days_Difference" and "age"

    # perform linear regression to get coefficients
    # @TODO: lm or glm?
    model = R.lm('log(Relative_Volume_Ratio) ~ Relative_Days_Difference * age + factor(gender)', data=gender_study_df)
    summary = R.summary(model)
    coeffs = summary.rx2("coefficients")
    print(summary)  # .rx2('coefficients'))

    # Build a DataFrame from the coefficients tables
    #print(coeffs.names)
    #result = pd.DataFrame(pandas2ri.py2ri(coeffs), index=coeffs.names[0], columns=coeffs.names[1])
    #print(result)

    # Non-linear regression using lmfit (Python package)


    # Gompertzian growth curve

    sns.scatterplot(x=gender_study_df["Relative_Days_Difference"], y=gender_study_df["Relative_Volume_Ratio"],
                    hue=gender_filtered, legend=True)
    plt.show()

    # @TODO: Challenging to set good initial values - optimization just fails...
    gender_study_df['a'] = 1.0
    gender_study_df['b'] = 13.0
    gender_study_df['c'] = 0.09
    startvector = ro.ListVector({'a': 1.0, 'b': 13.0, 'c': 0.09})
    gomp_formula = stats.as_formula('Volumes ~ a * exp(b * exp(c * Relative_Days_Difference))')
    # gomp_model = nlme.nls('Volumes ~ SSgompertz(Relative_Days_Difference,)')
    gomp_model = stats.nls(gomp_formula, data=gender_study_df, start=startvector)

    print(gomp_model)


if __name__ == "__main__":
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../data/")
    preprocess(data_path)
