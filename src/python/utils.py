import numpy as np
from datetime import datetime
from scipy.stats import norm


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


def str2datetime(str_):
    return datetime.strptime(str_, '%Y-%m-%d')


def calculate_volumetric_diameter(volume):
    return (6 * volume / np.pi) ** (1 / 3)


def BCa_interval(X, func, B=1000, q=0.975):
    theta_hat = func(X)

    N = len(X)
    order = np.array(range(N))
    order_boot = np.random.choice(order, size=(B, N), replace=True)
    X_boot = X[order_boot]
        
    # bootstrap
    theta_hat_boot = np.array([func(X_boot[i]) for i in range(X_boot.shape[0])])

    # 1) find jackknife estimates
    tmp = np.transpose(np.reshape(np.repeat(order, repeats=len(order)), (len(order), len(order))))  # make NxN matrix
    tmp_mat = tmp[~np.eye(tmp.shape[0], dtype=bool)].reshape(tmp.shape[0], -1)
    X_tmp_mat = X[tmp_mat]

    jk_theta = np.array([func(X_tmp_mat[i]) for i in range(X_tmp_mat.shape[0])])
    phi_jk = np.mean(jk_theta) - jk_theta  # jackknife estimates

    # 2) Find a
    a = 1 / 6 * np.sum(phi_jk ** 3) / np.sum(phi_jk ** 2) ** (3 / 2)

    # 3) Find b
    b = norm.ppf(np.sum(theta_hat_boot < theta_hat) / B)  # inverse standard normal

    # 4) Find gamma values -> limits
    def gamma1_func(a, b, q):
        return norm.cdf(b + (b + norm.ppf(1 - q)) / (1 - a * (b + norm.ppf(1 - q))))

    def gamma2_func(a, b, q):
        return norm.cdf(b + (b + norm.ppf(q)) / (1 - a * (b + norm.ppf(q))))

    # 5) get confidence interval of BCa
    CI_BCa = np.percentile(theta_hat_boot, [100 * gamma1_func(a, b, q), 100 * gamma2_func(a, b, q)])

    return CI_BCa, theta_hat_boot
