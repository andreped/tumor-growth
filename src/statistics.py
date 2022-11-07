from rpy2.robjects.packages import importr
from rpy2 import robjects as ro
import numpy as np


# R-related imports
stats = importr('stats')
R = ro.r


def kruskal_wallis_test_prompt(dependent_variable, data):
    if dependent_variable in ["genders", "T2", "oedema", "Multifocality"]:
        curr_formula = stats.as_formula('yearly_growth ~ factor(' + dependent_variable + ')')
    else:
        curr_formula = stats.as_formula('yearly_growth ~ ' + dependent_variable)
    result = stats.kruskal_test(curr_formula, data=data, **{'na.action': stats.na_omit})
    print("kruskal wallis test, yearly_growth vs " + dependent_variable + ". p-value:", result.rx2('p.value')[0])


def test_univariate_normality(x, full_data_nonzero):
    """
    Shapiro-Wilk suitable for datasets with N < 50, whereas Kolmogorov-Smirnov tests more suitable for N > 50
    """
    x = full_data_nonzero["Volume"]
    result = stats.ks_test(x, "pnorm", mean=np.mean(x), sd=np.std(x))
    # result = stats.shapiro_test(x)
    return result.rx2('p.value')[0]


def wilcox_test_custom(x):
    R.assign("variable", x)
    R("res <- wilcox.test(variable, exact=TRUE, conf.int=TRUE, conf.level=0.95)")
    r_result = R("res")
    return r_result.rx2('p.value')[0], r_result.rx2('conf.int').tolist()
