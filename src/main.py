import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def sort_timestamps(input_):
    tmp = np.array([int(x[1:]) for x in input_])
    return input_[np.argsort(tmp)]


def remove_surgery_patients(patients):
    filter_ = []
    for pat in patients:
        filter_.append(not pat.startswith("O"))
    filter_ = np.array(filter_)
    patients_no_surgery = list(patients[filter_])
    N = len(patients)
    print("Number of patients with surgery vs total:", N - len(patients_no_surgery), "out of", N)
    return patients_no_surgery, filter_


def plot_graphs(data):
    sns.set_style('darkgrid')
    sns.set(rc={'figure.figsize': (14, 8)})

    ax = sns.lineplot(data=data, x='OP_ID', y='Volume',
                      #hue='District',
                      palette='viridis',
                      legend='full', lw=3)

    #ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
    #plt.legend(bbox_to_anchor=(1, 1))
    #plt.ylabel('PM2.5 (µg/m3)')
    #plt.xlabel('Year-Month')
    plt.show()


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
    #print(len(unique_patients))
    #print(len(np.unique(filtered_cohort_volumes_quality["Patient"])))
    #exit()
    allowed = np.array(["T" + str(i) for i in range(9000)])
    iter = 0
    for pat in unique_patients:
        # get all data for current patient
        curr_data = filtered_volumes[filtered_volumes["OP_ID"] == pat]

        # get unique timestamps
        curr_timestamps = curr_data["Timestamp"]
        #print(curr_timestamps)
        curr_timestamps = list(curr_timestamps)
        for t in curr_timestamps:
            #print(t, allowed)
            #print()
            if t not in allowed:
                raise ValueError("jail time")
        #exit()
        unique_timestamps = np.unique(list(curr_timestamps))
        unique_timestamps = sort_timestamps(unique_timestamps)

        #print(unique_timestamps)

        #if len(unique_timestamps) != len(curr_timestamps):
        #    raise ValueError("you fucked up")
        print(len(unique_timestamps), len(curr_timestamps))

        # for each time stamp, all clusters and sum these into one value (total tumor amount in ml)
        #vols = []
        for timestamp in unique_timestamps:
            tmp = curr_data[curr_data["Timestamp"] == timestamp]
            curr_volume = sum(list(tmp["Volume"]))

            data.append([pat, timestamp, curr_volume])
            iter += 1
        print()

    print("iter:", iter)

    print("\n\n\n---")
    data = np.array(data)

    # verify if patient - timepoint pair exist in cohort volumes quality
    for i in range(len(filtered_cohort_volumes_quality["Patient"])):
        pat = list(filtered_cohort_volumes_quality["Patient"])[i]
        ts = list(filtered_cohort_volumes_quality["Timestamp"])[i]
        print(i, pat, ts)
        #if ()
        #print(data)
        success = False
        for j in range(data.shape[0]):
            curr = data[j]
            #print(curr)
            #print(curr[0], pat)
            #print(curr[1], ts)
            #print("-")
            if (curr[0] == pat) and (curr[1] == ts):
                success = True

            #if curr[0] == pat:
            #pr
        print(success)
        if not success:
            raise ValueError("here there is something wrong...")
        #exit()
        #print("####")

    exit()
    print(data)

    print(data.shape)
    print(filtered_cohort_volumes_quality.shape)

    print(filtered_cohort_volumes_quality)

    # compare patients and timestamps to see if there is a patient there is something wrong with
    # for pat in np.unique(data[])

    # merge this with the cohort volumes quality stuff
    # full_data = np.stack([filtered_cohort_volumes_quality, data[..., 1:]], axis=1)

    full_data = []
    counter = 0
    patters = np.unique(list(data[:, 0]))
    for pat in np.unique(list(filtered_cohort_volumes_quality["Patient"])):
        #tmp = filtered_cohort_volumes_quality[pat]
        #full_data.append(tmp)
        print(pat)
        ts = filtered_cohort_volumes_quality["Timestamp"][filtered_cohort_volumes_quality["Patient"] == pat]
        uniques = np.unique(list(ts))
        print(filtered_cohort_volumes_quality["Timestamp"])
        print(ts)
        print(len(uniques))
        print(len(ts))
        if len(uniques) != len(ts):
            raise ValueError("you fckd up")
        #print(pat, patters)
        if pat not in patters:
            raise ValueError("one patient is sus")
        counter += len(uniques)
        #print(ts)
        print()

    print("patient counts:")
    print()

    print("counts:")
    print(filtered_cohort_volumes_quality.shape)
    print(counter)
    print(data.shape)

    # start merging one by one and see where it fails
    new = filtered_cohort_volumes_quality.copy()
    for i in range(data.shape[0]):
        print(i, data[i])
    vols = -999 * np.ones(filtered_cohort_volumes_quality.shape[0])
    new['Volumes'] = vols


    #new.loc[(new['Patient'] > ), 'Match'] = df['Percent'].shift(-1)

    print(new)
    exit()

    print(full_data)





if __name__ == "__main__":
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../data/")
    preprocess(data_path)
