import pandas as pd
import os
import numpy as np
import statsmodels.api as sm
from utils import sort_timestamps, remove_surgery_patients, plot_graphs


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
    np.random.shuffle(unique_patients)
    iter = 0
    for pat in unique_patients:
        # get all data for current patient
        curr_data = filtered_volumes[filtered_volumes["OP_ID"] == pat]

        # get unique timestamps
        curr_timestamps = curr_data["Timestamp"]
        curr_timestamps = list(curr_timestamps)
        unique_timestamps = np.unique(list(curr_timestamps))
        unique_timestamps = sort_timestamps(unique_timestamps)

        # for each time stamp, all clusters and sum these into one value (total tumor amount in ml)
        for timestamp in unique_timestamps:
            tmp = curr_data[curr_data["Timestamp"] == timestamp]
            tmp = np.nan_to_num(tmp["Volume"])  # convert 0 to NaN for summation
            curr_volume = sum(list(tmp))
            data.append([pat, timestamp, curr_volume])
            iter += 1
    data = np.array(data)

    # merge this with the cohort volumes quality stuff
    #full_data = np.concatenate([filtered_cohort_volumes_quality, data[:, 2:3]], axis=1)
    full_data = filtered_cohort_volumes_quality.copy()
    full_data["Volumes"] = data[:, -1]

    print(full_data.shape)
    print(full_data)

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
            #print("r:", r)
            full_data.loc[r + 60, 'Birth_Year'] = byear  # @FIXME: Need better solution, quick-fix: 60 off in row ID...
            full_data.loc[r + 60, 'Gender'] = gender

    print(full_data)
    print(full_data.head(30))

    # remove all occurences where Volumes=0 (not possible!)
    filter_ = full_data["Volumes"] != str(0.0)
    full_data_nonzero = full_data[filter_]

    print("\nRemoved all volumes=0 incidents:")
    print(full_data_nonzero.head(30))
    print("\nold vs new shape after filtering:", full_data.shape, full_data_nonzero.shape)

    ### Make scatter plots and fit curve (regression)
    # - 1) test linear regression
    #model = sm.GLM(np.array(full_data.Volumes).astype(float), pd.get_dummies(full_data.Gender))  # , family=sm.families.Gamma())
    #results = model.fit()
    #print(results.summary())


if __name__ == "__main__":
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../data/")
    preprocess(data_path)
