import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime


def sort_timestamps(input_):
    tmp = np.array([int(x[1:]) for x in input_])
    return input_[np.argsort(tmp)]


def get_earliest_timestamp(input_):
    value = min([int(x[1:]) for x in input_])
    return "T" + str(value)


def get_last_timestamp(input_):
    value = max([int(x[1:]) for x in input_])
    return "T" + str(value)


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


def str2datetime(str_):
    return datetime.strptime(str_, '%Y-%m-%d')


def calculate_volumetric_diameter(volume):
    return (6 * volume / np.pi) ** (1 / 3)
