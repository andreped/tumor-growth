import pandas as pd
import os
import numpy as np
from utils import sort_timestamps, remove_surgery_patients, str2datetime, get_earliest_timestamp,\
    get_last_timestamp
from statistics import kruskal_wallis_test_prompt, test_univariate_normality, wilcox_test_custom
import seaborn as sns
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import scipy
from scipy.stats import spearmanr


# R-related imports
stats = importr('stats')
base = importr('base')
nlme = importr('nlme')
car = importr('car')
outliers = importr('outliers')
minpack_lm = importr('minpack.lm')
# nlsfit = importr('nlsfit')  # not compatible with R-4.2.2
pandas2ri.activate()
R = ro.r


def study_inter_rater_variability(data_path):
    interrater_variability_study = pd.read_csv(os.path.join(data_path, "interrater_variability_study.csv"))

    print(interrater_variability_study)

    volume = interrater_variability_study["Raw volume"]
    dice = interrater_variability_study.Dice

    print(volume)
    print(dice)

    x = interrater_variability_study["Raw volume"].astype("float32")
    y = interrater_variability_study["Per volume"].astype("float32")

    mu = np.mean([x, y], axis=0)
    diff = x - y
    md = np.mean(diff)
    sd = np.std(diff, axis=0)

    print(md, sd)
    print(mu - 1.96 * sd, mu + 1.96 * sd)

    # outlier detection
    outliers = (diff > md + 1.96 * sd) | (diff < md - 1.96 * sd)
    outliers = np.array(outliers)

    print(diff > md + 1.96 * sd)
    print(outliers)
    print(sum(outliers))

    x = np.array(x)
    y = np.array(y)

    x_filtered = x[np.logical_not(outliers)]
    y_filtered = y[np.logical_not(outliers)]

    # redo bland-altman plot
    mu = np.mean([x_filtered, y_filtered], axis=0)
    diff = x_filtered - y_filtered
    md = np.mean(diff)
    sd = np.std(diff, axis=0)

    plt.scatter(mu, diff)
    plt.axhline(md, color='gray', linestyle='--')
    plt.axhline(md + 1.96 * sd, color='gray', linestyle='--')
    plt.axhline(md - 1.96 * sd, color='gray', linestyle='--')

    # print((diff < md - 1.96 * sd) or (diff > md + 1.96 * sd))

    plt.title('Bland-Altman Plot')
    plt.show()

    #sns.scatterplot(x, y)
    #plt.show()


def get_patient_characteristics(data_path):
    cohort_personal_info = pd.read_csv(os.path.join(data_path, "cohort_personal_info.csv"))
    print(cohort_personal_info)


def preprocess(data_path):
    cohort_personal_info = pd.read_csv(os.path.join(data_path, "cohort_personal_info.csv"))
    cohort_volumes_quality = pd.read_csv(os.path.join(data_path, "cohort_volumes_quality-removed_P09016-T6.csv"))
    # interrater_variability_study = pd.read_csv(os.path.join(data_path, "interrater_variability_study.csv"))  # not considered for now
    volumes = pd.read_csv(os.path.join(data_path, "Volumes.csv"))
    t2_oedema = pd.read_csv(os.path.join(data_path, "T2_and_peritumorial_oedema.csv"), sep=";")
    scanner_info = pd.read_csv(os.path.join(data_path, "scanners_info.csv"), sep=",")

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

    filter_ = [x in patients_no_surgery for x in t2_oedema["Patient"]]
    filtered_t2_oedema = t2_oedema[filter_]

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

        # get earliest timestamp with non-NaN or non-zero volume
        for t in unique_timestamps:
            tmp = curr_data[curr_data["Timestamp"] == t]
            curr_v = np.array(tmp["Volume"]).astype("float32")
            curr_v = sum(curr_v)
            if (curr_v != 0) and not pd.isnull(curr_v):
                break
        init_timestamp = t
        # print("final timestamp (and size):", init_timestamp, curr_v)

        # get last timestamp with non-NaN volume
        for t in unique_timestamps[::-1]:  # reversed ordered timestamp list
            tmp = curr_data[curr_data["Timestamp"] == t]
            curr_v = np.array(tmp["Volume"]).astype("float32")
            curr_v = sum(curr_v)
            #if not pd.isnull(curr_v):
            if (curr_v != 0) and not pd.isnull(curr_v):
                break
        last_timestamp = t

        # get date of first timestamp of patient (T1 might not be the earliest!! If T1 is NaN or 0, then T2 might be...
        first_timestamp_date = curr_data[curr_data["Timestamp"] == init_timestamp]["Date"]
        first_timestamp_date = str2datetime(np.asarray(first_timestamp_date)[0])

        # get last timestamp
        last_timestamp_date = curr_data[curr_data["Timestamp"] == last_timestamp]["Date"]
        last_timestamp_date = str2datetime(np.asarray(last_timestamp_date)[0])

        # get initial volume size at first scan
        initial_volume = curr_data[curr_data["Timestamp"] == init_timestamp]["Volume"]
        initial_volume = np.asarray(initial_volume)[0]

        # get final volume size at last scan
        last_volume = curr_data[curr_data["Timestamp"] == last_timestamp]["Volume"]
        last_volume = np.asarray(last_volume)[0]

        # get relative volume change
        relative_volume_change = (last_volume - initial_volume) / initial_volume

        # get cluster numbers for current patient (if above 1, multifocal "by definition")
        multifocality = int(np.any(curr_data["Clusters total"] > 1))
        clusters_total = max(curr_data["Clusters total"])

        # counter number of timestamps
        # nb_timestamps = len(curr_data["Timestamp"])
        nb_timestamps = len(unique_timestamps)

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

            # check if earliest timestamp, store categorical value
            earliest_timestamp = int(timestamp == init_timestamp)
            #if earliest_timestamp == 1:
            #    print(timestamp)
            #    print(init_timestamp)
            #    print(earliest_timestamp)
            #    print(curr_volume)
            #    print()

            # translate current date to datetime format
            curr_date = str2datetime(curr_date)

            data.append([pat, timestamp, relative_volume_change, clusters_total, multifocality, earliest_timestamp,
                         nb_timestamps, initial_volume, last_volume, first_timestamp_date,
                         last_timestamp_date, curr_date, curr_volume])
            iter += 1
    data = np.array(data)

    # merge this with the cohort volumes quality stuff
    full_data = filtered_cohort_volumes_quality.copy()

    full_data["Volume"] = data[:, -1].astype("float32")
    full_data["Date"] = data[:, -2]
    full_data["Last_Timestamp_Date"] = data[:, -3]
    full_data["First_Timestamp_Date"] = data[:, -4]
    full_data["Final_Volume"] = data[:, -5].astype("float32")
    full_data["Initial_Volume"] = data[:, -6].astype("float32")
    full_data["Number_Of_Timestamps"] = data[:, -7].astype("float32")
    full_data["Earliest_Timestamp"] = data[:, -8]
    full_data["Multifocality"] = data[:, -9]
    full_data["Clusters_total"] = data[:, -10]
    full_data["Relative_Volume_Change"] = data[:, -11].astype("float32")

    # need to filter NaN volumes on merged data frame
    full_data = full_data[full_data.Volume != 0]

    # reset indices in full_data to go 0:1:N
    full_data.index = list(range(len(full_data)))

    # add patient characteristics to full data frame
    full_data["Birth_Year"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    full_data["Gender"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    for pat, gender, byear in zip(cohort_personal_info_filtered["Patient"], cohort_personal_info_filtered["Gender"],
                                  cohort_personal_info_filtered["Birth_Year"]):
        row_ids = np.where(full_data["Patient"] == pat)[0]
        for r in row_ids:
            byear_new_format = str(byear) + "-07-01"
            byear_new_format = str2datetime(byear_new_format)

            full_data.loc[r, 'Birth_Year'] = byear_new_format
            full_data.loc[r, 'Gender'] = gender

    # add T2 and oedema information to full data frame
    full_data["T2"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    full_data["Oedema"] = (-999 * np.ones(full_data.shape[0])).astype("str")
    for patient in filtered_t2_oedema["Patient"]:
        row_ids = np.where(full_data["Patient"] == patient)[0]
        curr = filtered_t2_oedema[filtered_t2_oedema["Patient"] == patient]
        for r in row_ids:
            full_data.loc[r, "T2"] = np.array(curr["T2"])
            full_data.loc[r, "Oedema"] = np.array(curr["peritumorial_oedema"])

    # add scanner info to the full data frame
    full_data["Manufacturer"] = (np.nan * np.ones(full_data.shape[0])).astype("str")
    full_data["Model_Name"] = (np.nan * np.ones(full_data.shape[0])).astype("str")
    full_data["Tesla"] = (np.nan * np.ones(full_data.shape[0]))
    for i in range(len(scanner_info)):
        patient, timestamp, manufacturer, model_name, tesla = scanner_info.loc[i]
        row_id = np.where((full_data["Patient"] == patient) & (full_data["Timestamp"] == timestamp))[0]

        if len(row_id) == 0:
            continue

        full_data.loc[row_id[0], "Manufacturer"] = manufacturer
        full_data.loc[row_id[0], "Model_Name"] = model_name
        full_data.loc[row_id[0], "Tesla"] = float(tesla)

    # remove all occurences where Volumes=0 (not possible -> tumor was not annotated)
    filter_zero_volumes = full_data["Volume"] != str(0.0)
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
    relative_volume_ratio = full_data_nonzero["Volume"] / full_data_nonzero["Initial_Volume"]
    full_data_nonzero["Relative_Volume_Ratio"] = relative_volume_ratio.astype("float32")

    # @TODO: Should normalize the variables in some way to avoid exploding stuff issues
    #   - perhaps age in years (instead of in days) make more sense?

    print(full_data_nonzero.head(40))

    # capture outliers and remove them - assuming normality using quantiles
    #lower = R.quantile(full_data_nonzero["Relative_Volume_Ratio"], 0.025)[0]
    #higher = R.quantile(full_data_nonzero["Relative_Volume_Ratio"], 0.975)[0]

    #x = np.array(full_data_nonzero["Relative_Volume_Ratio"])
    #filter_ = (x > lower) & (x < higher)
    #full_data_nonzero = full_data_nonzero[filter_]

    # filter patients that show no growth? - how to determine if tumor has grown?
    # Look at first and last timestep volume size?
    volume_change = full_data_nonzero["Final_Volume"] - full_data_nonzero["Initial_Volume"]
    print(volume_change)

    # remove patients with slice thickness higher or equal than X
    # slice_thickness_filter = np.array(full_data_nonzero["Spacing3"]) < 2
    # full_data_nonzero = full_data_nonzero[slice_thickness_filter]

    # create summary statistics for study - Table 1
    full_data_nonzero["Current_Age_Years"] = np.array(full_data_nonzero["Current_Age"]).astype("float32") / 365.25

    # patient_filter_ = full_data_nonzero["Timestamp"] == "T1"
    # @TODO: After removing volumes with 0 size, some T1 points are now missing (FIXED BELOW)
    patient_filter_ = np.array(full_data_nonzero["Earliest_Timestamp"]) == 1

    # multifocality
    multifocality = np.array(full_data_nonzero["Multifocality"][patient_filter_])
    Clusters_total = np.array(full_data_nonzero["Clusters_total"][patient_filter_])

    # age_at_T1 = np.array(full_data_nonzero["Current_Age"][patient_filter_]).astype("float32") / 365.25
    genders = np.array(full_data_nonzero["Gender"][patient_filter_])

    # init_volume_size = np.array(full_data_nonzero["Volume"][patient_filter_])
    number_of_mri_scans = np.array(full_data_nonzero["Number_Of_Timestamps"][patient_filter_])
    slice_thickness = np.array(full_data_nonzero["Spacing3"][patient_filter_])

    t2_hyperintense_orig = np.array(full_data_nonzero["T2"][patient_filter_])
    t2_hyperintense_orig[t2_hyperintense_orig == '-999.0'] = np.nan
    #t2_hyperintense = t2_hyperintense_orig[t2_hyperintense_orig != '-999.0']
    t2_hyperintense = t2_hyperintense_orig[~pd.isnull(t2_hyperintense_orig)]
    t2_hyperintense = np.array([int(x) for x in t2_hyperintense])

    oedema_orig = np.array(full_data_nonzero["Oedema"][patient_filter_])
    oedema_orig[oedema_orig == "-999.0"] = np.nan
    #oedema = oedema_orig[oedema_orig != '-999.0']
    oedema = oedema_orig[~pd.isnull(oedema_orig)]
    oedema = np.array([int(x) for x in oedema])

    earliest_timestamp_filter = np.array(full_data_nonzero["Earliest_Timestamp"] == 1)
    init_volume_size = np.array(full_data_nonzero["Volume"][earliest_timestamp_filter])
    age_at_T1 = np.array(full_data_nonzero["Current_Age"][earliest_timestamp_filter]).astype("float32") / 365.25

    patients = np.unique(full_data_nonzero["Patient"])
    total_follow_up_days = []
    for patient in patients:
        curr = full_data_nonzero[full_data_nonzero["Patient"] == patient]["Relative_Days_Difference"]
        total_follow_up_days.append(max(curr))
    total_follow_up_months = np.array(total_follow_up_days) / 30

    # @TODO: Which threshold to use? Base it on measurements error (largest error, quantiles?) in
    #   inter-rater variability study? 15 % makes sense as it corresponds to the largest error in the inter-rater study
    relative_growth_threshold = 0.15

    full_data_nonzero_grew_only = full_data_nonzero.copy()
    volume_change = []
    volume_change_relative = []
    volume_grew = []
    volume_shrank = []
    volume_no_change = []
    volume_change_categorical = []
    grow_patients = []
    for patient in patients:
        curr = full_data_nonzero[full_data_nonzero["Patient"] == patient]
        first_timestamp = get_earliest_timestamp(curr["Timestamp"])
        # first_timestamp = np.array(curr["Timestamp"][curr["Earliest_Timestamp"] == 1])
        last_timestamp = get_last_timestamp(curr["Timestamp"])
        initial_size = np.array(curr[curr["Timestamp"] == first_timestamp]["Volume"])[0]
        final_size = np.array(curr[curr["Timestamp"] == last_timestamp]["Volume"])[0]

        print("--")
        print(initial_size, final_size)
        print(curr["Initial_Volume"])
        print(curr["Final_Volume"])

        relative_change = (final_size - initial_size) / initial_size
        initial_size = float(initial_size)
        final_size = float(final_size)
        volume_change.append(final_size - initial_size)
        volume_change_relative.append(relative_change)

        if relative_change > relative_growth_threshold:
            volume_change_categorical.append(1)
            volume_grew.append([patient, initial_size, relative_change])
        elif relative_change < - relative_growth_threshold:
            volume_change_categorical.append(-1)
            volume_shrank.append([patient, initial_size, relative_change])
        else:
            volume_change_categorical.append(0)
            volume_no_change.append([patient, initial_size, relative_change])

    volume_change = np.array(volume_change)
    volume_change_relative = np.array(volume_change_relative)

    # get yearly growth
    yearly_growth = np.array(volume_change_relative) / (np.array(total_follow_up_days) / 356.25)

    print(full_data_nonzero)
    N = len(age_at_T1)
    # print(age_at_T1)
    print("age: median/IQR/min/max:", np.round(np.median(age_at_T1), 1), np.round(scipy.stats.iqr(age_at_T1), 1),
          np.round(np.min(age_at_T1), 1), np.round(np.max(age_at_T1), 1))
    print("gender (count/%):", sum(genders == "woman"), np.mean(genders == "woman"))
    print("initial volume size at T1: (median/IQR/min/max):", np.round(np.median(init_volume_size), 1),
          np.round(scipy.stats.iqr(init_volume_size), 1), np.round(np.min(init_volume_size), 1),
          np.round(np.max(init_volume_size), 1))
    print("number of MRI scans per patient: (median/IQR/min/max):", np.round(np.median(number_of_mri_scans), 1),
          np.round(scipy.stats.iqr(number_of_mri_scans), 1), np.round(np.min(number_of_mri_scans), 1),
          np.round(np.max(number_of_mri_scans), 1))
    print("total follow up in months: (median/IQR/min/max):", np.round(np.median(total_follow_up_months), 1),
          np.round(scipy.stats.iqr(total_follow_up_months), 1), np.round(np.min(total_follow_up_months), 1),
          np.round(np.max(total_follow_up_months), 1))
    print("volume change (T1 to T-last): (median/IQR/min/max):", np.round(np.median(volume_change), 1),
          np.round(scipy.stats.iqr(volume_change), 1), np.round(np.min(volume_change), 1),
          np.round(np.max(volume_change), 1))
    print("relative volume change (T1 to T-last): (median/IQR/min/max):", np.round(np.median(volume_change_relative), 3),
          np.round(scipy.stats.iqr(volume_change_relative), 3), np.round(np.min(volume_change_relative), 3),
          np.round(np.max(volume_change_relative), 3))
    print("yearly relative growth (T1 to T-last): (median/IQR/min/max):",
          np.round(np.median(yearly_growth), 3),
          np.round(scipy.stats.iqr(yearly_growth), 3), np.round(np.min(yearly_growth), 3),
          np.round(np.max(yearly_growth), 3))
    print("number of patients with tumors that grew/no change/shrank:",
          len(volume_grew), len(volume_no_change), len(volume_shrank), "| % |",
          np.round(len(volume_grew) / N, 3),
          np.round(len(volume_no_change) / N, 3),
          np.round(len(volume_shrank) / N), 3)
    print("slice thickness: (median/IQR/min/max):",
          np.round(np.median(slice_thickness), 3),
          np.round(scipy.stats.iqr(slice_thickness), 3), np.round(np.min(slice_thickness), 3),
          np.round(np.max(slice_thickness), 3))
    print("multifocality (count + %):", sum(multifocality), sum(multifocality) / len(multifocality))
    print("clusters total: (median/IQR/min/max):",
          np.round(np.median(Clusters_total), 3),
          np.round(scipy.stats.iqr(Clusters_total), 3), np.round(np.min(Clusters_total), 3),
          np.round(np.max(Clusters_total), 3))
    print("Manufacturer (count + %):", np.unique(full_data_nonzero["Manufacturer"], return_counts=True),
          np.unique(full_data_nonzero["Manufacturer"], return_counts=True)[1] / len(full_data_nonzero))
    print("Model_Name (count + %):", np.unique(full_data_nonzero["Model_Name"], return_counts=True),
          np.unique(full_data_nonzero["Model_Name"], return_counts=True)[1] / len(full_data_nonzero))
    print("Tesla (count + %):", np.unique(full_data_nonzero["Tesla"], return_counts=True),
          np.unique(full_data_nonzero["Tesla"], return_counts=True)[1] / len(full_data_nonzero))

    print(wilcox_test_custom(yearly_growth))

    # T2 hyperintense signal - calculate summary statistics
    tmp = np.unique(t2_hyperintense, return_counts=True)
    print("T2 hyperintense signal (counts for 1/2/3 categories):",
          tmp, tmp[1] / len(t2_hyperintense))
    tmp = np.unique(oedema, return_counts=True)
    print("Oedema (counts for 0/1 categories):",
          tmp, tmp[1] / len(oedema))

    # perform ANOVA to assess the association between tumor growth and different variables
    # 1) first test for normality
    print("Univariate normality tests:")
    print("p-value tumor growth:", test_univariate_normality(volume_change, full_data_nonzero))
    print("p-value Age:", test_univariate_normality(volume_change_relative, full_data_nonzero))
    print("p-value Slice thickness:", test_univariate_normality(full_data_nonzero["Spacing3"], full_data_nonzero))
    # print("p-value T2:", test_univariate_normality(full_data_nonzero["T2"]))
    # print("p-value Oedema:", test_univariate_normality(full_data_nonzero["Oedema"]))
    # -> Data is NOT normal! Cannot perform ordinary ANOVA analysis

    # create temporary dataframe to store data relevant for statistical analysis
    df_association = pd.DataFrame({
        "volume_change": volume_change,
        "volume_change_relative": volume_change_relative,
        "init_volume_size": init_volume_size,
        "age_at_T1": age_at_T1,
        "total_follow_up_months": total_follow_up_months,
        "T2": t2_hyperintense_orig.astype("float32"),
        "oedema": oedema_orig.astype("float32"),
        "genders": genders,
        "Spacing3": slice_thickness,
        "Multifocality": multifocality,
        "yearly_growth": yearly_growth,
        "volume_change_categorical": volume_change_categorical,
    })

    df_association_dropped_na = df_association.copy()
    df_association_dropped_na["T2"] = t2_hyperintense_orig
    df_association_dropped_na["oedema"] = oedema_orig
    df_association_dropped_na = df_association_dropped_na.dropna()
    df_association_dropped_na["T2"] = df_association_dropped_na["T2"].astype(int)
    df_association_dropped_na["oedema"] = df_association_dropped_na["oedema"].astype(int)
    df_association_dropped_na["yearly_growth"] = df_association["yearly_growth"]

    for var in ["age_at_T1", "total_follow_up_months", "init_volume_size", "Spacing3"]:
        corr_ = spearmanr(df_association["yearly_growth"], df_association[var])
        print("spearman correlation (" + var + "):", corr_, "\n")

    for var in ["genders", "T2", "oedema", "Multifocality"]:
        kruskal_wallis_test_prompt(var, data=df_association_dropped_na)

    print("\nHandle growth as a categorical variable:\n")
    for var in ["age_at_T1", "total_follow_up_months", "init_volume_size", "Spacing3"] + \
               ["genders", "T2", "oedema", "Multifocality"]:
        print("var:", var)
        if var in ["age_at_T1", "total_follow_up_months", "init_volume_size", "Spacing3"]:
            curr_formula = stats.as_formula('factor(' + var + ') ~ factor(volume_change_categorical)')
            result = stats.kruskal_test(curr_formula, data=df_association, **{'na.action': stats.na_omit})
            print("Kruskal wallis test, volume_change_categorical vs " + var + ". p-value:", result.rx2('p.value')[0])
        else:
            result = stats.chisq_test(df_association["yearly_growth"], df_association[var])
            print("Chisq test, volume_change_categorical vs " + var + ". p-value:", result.rx2('p.value')[0])


    print("\n######\n", "Now, lets study the growth pattern only for the tumours that grew")

    # df_association_grew_only = df_association[df_association["volume_change_categorical"] == 1]

    full_data_nonzero_grew_only = full_data_nonzero[full_data_nonzero["Relative_Volume_Change"] > 0.15]

    print(list(full_data_nonzero_grew_only.keys()))
    print(full_data_nonzero_grew_only)

    df_regression_association = pd.DataFrame({
        "Relative_Volume_Change": full_data_nonzero_grew_only["Relative_Volume_Change"],
        "Follow_Up_Months": full_data_nonzero_grew_only["Follow_Up_Months"],
    })

    # linear growth
    print("\nLINEAR:")
    model = R.lm('Relative_Volume_Change ~ Follow_Up_Months', data=df_regression_association)
    summary_model = R.summary(model)
    print(summary_model)  # .rx2('coefficients'))
    # print("AIC:", stats.extractAIC(model))

    # exponential growth
    print("\nEXPONENTIAL:")
    model = R.lm('log(Relative_Volume_Change) ~ Follow_Up_Months', data=df_regression_association)
    summary_model = R.summary(model)
    print(summary_model)  # .rx2('coefficients'))
    # print("AIC:", stats.extractAIC(model))

    # linear radial growth
    print("\nLinear radial growth")
    df_regression_association['Am'] = max(df_regression_association["Relative_Volume_Change"]) - \
                                      min(df_regression_association["Relative_Volume_Change"])
    df_regression_association['Rd'] = - min(df_regression_association["Relative_Volume_Change"])
    df_regression_association['LCP'] = max(df_regression_association["Relative_Volume_Change"]) - 1
    startvector = ro.ListVector({
        'Am': max(df_regression_association["Relative_Volume_Change"]) - min(
            df_regression_association["Relative_Volume_Change"]),
        'Rd': - min(df_regression_association["Relative_Volume_Change"]),
        'LCP': max(df_regression_association["Relative_Volume_Change"]) - 1
    })
    gomp_model = minpack_lm.nlsLM("Relative_Volume_Change ~ Am * (1 - ((1 - Rd/Am))^(1 - (Follow_Up_Months / LCP)))",
                                  start=startvector, data=df_regression_association)

    summary_model = R.summary(gomp_model)
    print(summary_model)
    # print("AIC:", stats.extractAIC(gomp_model))

    # gompertzian growth
    df_regression_association['Am'] = max(df_regression_association["Relative_Volume_Change"]) - \
                                      min(df_regression_association["Relative_Volume_Change"])
    df_regression_association['Rd'] = - min(df_regression_association["Relative_Volume_Change"])
    df_regression_association['LCP'] = max(df_regression_association["Relative_Volume_Change"]) - 1

    startvector = ro.ListVector({
        'Am': max(df_regression_association["Relative_Volume_Change"]) - min(df_regression_association["Relative_Volume_Change"]),
        'Rd': - min(df_regression_association["Relative_Volume_Change"]),
        'LCP': max(df_regression_association["Relative_Volume_Change"]) - 1
    })

    print(np.unique(df_regression_association["Relative_Volume_Change"]))
    print(np.unique(df_regression_association["Follow_Up_Months"]))
    #gomp_model = stats.nls("log(Relative_Volume_Change) ~ SSgompertz(log(Follow_Up_Months+1), Asym, b2, b3)",
    #                       data=df_regression_association)
    #data = dat,
    #       start = list(y0=0.15, ymax=7, k=0.5, lag=3))
    # gomp_model = stats.nls(gomp_formula, data=df_regression_association, start=startvector)

    gomp_model = minpack_lm.nlsLM("Relative_Volume_Change ~ Am * (1 - ((1 - Rd/Am))^(1 - (Follow_Up_Months / LCP)))",
                                  start=startvector, data=df_regression_association)

    summary_model = R.summary(gomp_model)
    print(summary_model)
    # print("AIC:", stats.extractAIC(gomp_model))

    # save processed data frame as CSV on disk
    df_association.to_csv("./data/merged_processed_regression_data_080123.csv")


if __name__ == "__main__":
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../data/")
    preprocess(data_path)
