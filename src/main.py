import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def preprocess(data_path):
    cohort_personal_info = pd.read_csv(os.path.join(data_path, "cohort_personal_info.csv"))
    cohort_volumes_quality = pd.read_csv(os.path.join(data_path, "cohort_volumes_quality.csv"))
    interrater_variability_study = pd.read_csv(os.path.join(data_path, "interrater_variability_study.csv"))
    volumes = pd.read_csv(os.path.join(data_path, "Volumes.csv"))

    print(cohort_personal_info.head(40))
    print(cohort_volumes_quality.head(40))
    print(interrater_variability_study.head(40))
    print(volumes.head(40))

    # get unique patients
    patients = cohort_personal_info["Patient"]
    print(patients)


    exit()

    plot_graphs(volumes)


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


if __name__ == "__main__":
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../data/")
    preprocess(data_path)

