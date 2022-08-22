import pandas as pd
import os
import numpy as np
from utils import sort_timestamps, remove_surgery_patients, plot_graphs
import seaborn as sns
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
import matplotlib.pyplot as plt
from datetime import datetime


pandas2ri.activate()
R = ro.r


def translate_date_to_years(data):
    return


def str2datetime(str_):
    return datetime.strptime(str_, '%Y-%m-%d')


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

        # get date of first timestamp of patient
        first_timestamp_date = curr_data[curr_data["Timestamp"] == "T1"]["Date"]
        first_timestamp_date = str2datetime(np.asarray(first_timestamp_date)[0])

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

            # translate current date to year
            #curr_date = int(curr_date.split("-")[0])

            data.append([pat, timestamp, first_timestamp_date, curr_date, curr_volume])
            iter += 1
    data = np.array(data)

    # merge this with the cohort volumes quality stuff
    #full_data = np.concatenate([filtered_cohort_volumes_quality, data[:, 2:3]], axis=1)
    full_data = filtered_cohort_volumes_quality.copy()
    full_data["Volumes"] = data[:, -1].astype(float)
    full_data["Date"] = data[:, -2]
    full_data["First_Timestamp_Date"] = data[:, -3]

    # indeces = full_data.index

    # add patient characteristics to full data frame
    full_data["Birth_Year"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    full_data["Gender"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    for pat, gender, byear in zip(cohort_personal_info_filtered["Patient"], cohort_personal_info_filtered["Gender"], cohort_personal_info_filtered["Birth_Year"]):
        row_ids = np.where(full_data["Patient"] == pat)[0]
        #if type(row_ids) is not list():
        #    row_ids = [row_ids]
        #print(full_data["Patient"] == pat)
        #print("row ids:", row_ids)
        for r in row_ids:
            byear_new_format = str(byear) + "-07-01"
            byear_new_format = str2datetime(byear_new_format)

            #print("r:", r)
            full_data.loc[r + 60, 'Birth_Year'] = byear_new_format  # @FIXME: Need better solution, quick-fix: 60 off in row ID...
            full_data.loc[r + 60, 'Gender'] = gender

    print(full_data.shape)
    print(full_data.head(30))

    # remove all occurences where Volumes=0 (not possible!)
    filter_zero_volumes = full_data["Volumes"] != str(0.0)
    full_data_nonzero = full_data[filter_zero_volumes]

    print("\nold vs new shape after filtering:", full_data.shape, full_data_nonzero.shape)
    del full_data  # sanity check, to avoid using the wrong volume

    print("\nRemoved all volumes=0 incidents:")
    print(full_data_nonzero.head(30))

    # get current age at scan and add to data frame
    print(full_data_nonzero["Date"])
    print(full_data_nonzero["Birth_Year"])

    curr_age_at_scan = full_data_nonzero["Date"] - full_data_nonzero["Birth_Year"]
    curr_age_at_scan = curr_age_at_scan.dt.days
    print(curr_age_at_scan)

    # curr_age_at_scan = full_data_nonzero["Date"].astype(int) - full_data_nonzero["Birth_Year"].astype(int)
    print(curr_age_at_scan)
    full_data_nonzero["Current_Age"] = curr_age_at_scan.astype(float)

    print(full_data_nonzero.head(30))

    # get relative difference in days between scans
    print(full_data_nonzero["First_Timestamp_Date"])
    relative_difference_between_scans = full_data_nonzero["Date"] - full_data_nonzero["First_Timestamp_Date"]
    relative_difference_between_scans = relative_difference_between_scans.dt.days
    print(relative_difference_between_scans)

    full_data_nonzero["Relative_Days_Difference"] = relative_difference_between_scans

    print(full_data_nonzero.head(30))

    ### Make scatter plots and fit curve (regression)
    # - 1) test linear regression
    #model = sm.GLM(np.array(full_data.Volumes).astype(float), pd.get_dummies(full_data.Gender))  # , family=sm.families.Gamma())
    #results = model.fit()
    #print(results.summary())

    # - 1) Regression modelling using rpy2 to run R code from Python
    #model = R.lm('Volumes ~ Current_Age', data=full_data_nonzero)
    #print(R.summary(model))  # .rx2('coefficients'))

    # - Plot scatter plot

    # sns.scatterplot(x="Relative_Days_Difference", y="Volumes", data=full_data_nonzero, hue="Gender", legend=True)
    #current_age = np.asarray(full_data_nonzero["Current_Age"]).astype(float)
    #vols = np.asarray(full_data_nonzero["Volumes"]).astype(float)
    #sns.scatterplot(x=list(range(len(vols))), y=vols)
    #x = sns.scatterplot(x=current_age, y=vols)
    #x.set_xlabel("he")

    x = np.asarray(full_data_nonzero["Relative_Days_Difference"])
    y = np.asarray(full_data_nonzero["Volumes"])
    gender = np.asarray(full_data_nonzero["Gender"])
    age = np.asarray(full_data_nonzero["Current_Age"])

    # @TODO: Should normalize the variables in some way to avoid exploding stuff issues
    #   - perhaps age in years (instead of in days) make more sense?

    filter_ = y > 0
    x_filtered = x[filter_]
    y_filtered = y[filter_]
    gender_filtered = gender[filter_]
    age_filtered = age[filter_]

    # new simplified, filtered data frame
    gender_study_df = pd.DataFrame({
        "gender": gender_filtered,
        "Volumes": y_filtered,
        "Relative_Days_Difference": x_filtered,
        "age": age_filtered
    })

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

    sns.scatterplot(
        x=x_filtered,
        y=y_log_filtered,
        hue=gender_filtered,
        legend=True
    )
    # -> looks like white noise (but appears to be some trends). Should apply

    #sns.lineplot(x=[min(x_filtered), max(x_filtered)], y=[mu_y_log_man, mu_y_log_man])
    #sns.lineplot(x=[min(x_filtered), max(x_filtered)], y=[mu_y_log_woman, mu_y_log_woman])

    # @FIXME: Probably interaction between "Relative_Days_Difference" and "age"

    # perform linear regression to get coefficients
    # @TODO: lm or glm?
    model = R.lm('log(Volumes) ~ Relative_Days_Difference * age + factor(gender)', data=gender_study_df)
    print(R.summary(model))  # .rx2('coefficients'))

    plt.show()


if __name__ == "__main__":
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../data/")
    preprocess(data_path)
