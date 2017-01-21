#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.signal import argrelmax
from sklearn.neighbors.kde import KernelDensity
from sklearn import cluster


tp = [210.0, 353.0, 403.0, 26.0, 238.0, 409.0, 75.0, 317.0, 279.0, 410.0, 96.0, 311.0, 108.0, 235.0, 375.0, 283.0, 264.0, 193.0, 258.0, 308.0, 355.0, 83.0, 415.0, 22.0, 278.0, 276.0, 442.0, 362.0, 342.0, 450.0, 404.0, 196.0, 372.0, 283.0, 291.0, 194.0, 560.0, 174.0, 251.0, 210.0, 481.0, 272.0, 416.0, 242.0, 313.0, 199.0, 411.0, 322.0, 300.0, 443.0, 237.0, 310.0, 405.0, 223.0, 223.0, 425.0, 227.0, 301.0, 317.0, 166.0, 287.0, 261.0, 253.0, 296.0, 175.0, 311.0, 483.0, 483.0, 141.0, 270.0, 471.0, 390.0, 285.0, 98.0, 279.0, 242.0, 420.0, 119.0, 310.0, 285.0, 412.0, 212.0, 194.0, 279.0, 103.0, 403.0, 228.0, 270.0, 298.0, 283.0, 220.0, 75.0, 411.0, 204.0, 363.0, 283.0, 248.0, 332.0, 326.0, 323.0, 307.0, 95.0, 110.0, 51.0, 336.0, 60.0, 354.0, 248.0, 456.0]
tn = [257.0, 282.0, 230.0, 147.0, 248.0, 189.0, 248.0, 301.0, 106.0, 271.0, 246.0, 296.0, 138.0, 100.0, 262.0, 166.0]


# Fit KDE
def kde_sklearn(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return kde_skl, np.exp(log_pdf)


# Find intersection variant one
def get_intersection(fun1, fun2, xo):
    return fsolve(lambda x : fun1(x) - fun2(x), xo, xtol=1e-01, maxfev=500)


def find_density_drop_idx(pdf_tn, x_axis):
    read_flag = 0
    for density_x in pdf_tn:
        if read_flag == 0 and density_x > 0.0001:
            read_flag = 1
        if read_flag == 1 and density_x < 0.0001:
            density_drop_on_x_axis = [x for x, y in enumerate(list(pdf_tn)) if y == density_x][0]
            return x_axis[density_drop_on_x_axis]
    return None


def best_intersection_if_curves_overlap(intersec_points, lower, upper):
    try:
        testlist = [w for w in intersec_points if lower < w < upper]
        print(testlist)
        for point in sorted(intersec_points):
            if lower < point < upper:
                return point
    except (TypeError, IndexError):
        return None


def find_density_intersection(tp, tn):
    x_axis = np.linspace(0, max(tp), 4000)
    tp_array = np.array(tp)
    tn_array = np.array(tn)
    bandwidth = cluster.estimate_bandwidth(np.array([[x] for x in x_axis]), quantile=0.1)
    kde_tp, pdf_tp = kde_sklearn(tp_array, x_axis, bandwidth=bandwidth)
    kde_tn, pdf_tn = kde_sklearn(tn_array, x_axis, bandwidth=bandwidth)
    highest_maxima_tn = [[pdf_tn[first_tn], x_axis[first_tn]] for first_tn in argrelmax(pdf_tn)[0]]
    maxima_x_axis = int(round(sorted(highest_maxima_tn, key=lambda x: x[0], reverse=True)[0][1], -1))
    funcA = lambda x: np.exp(kde_tp.score_samples([[x]][0]))
    funcB = lambda x: np.exp(kde_tn.score_samples([[x]][0]))
    intersection_interval = range(maxima_x_axis, int(max(tn)) + 10, 10)
    point_guesses = [round(float(get_intersection(funcA, funcB, guess))) for guess in intersection_interval]
    intersec_points = sorted(set(point_guesses), reverse=True)
    if int(min(tp)) > int(max(tn)):
        return round((max(tn) + min(tp)) / 2)
    else:
        tn_curve_end = find_density_drop_idx(pdf_tn, x_axis)
        best_intersec = best_intersection_if_curves_overlap(intersec_points, maxima_x_axis, tn_curve_end) #float(min(tn))
        #plt.plot(x_axis, pdf_tp, color='green')
        #plt.plot(x_axis, pdf_tn, color='blue')
        #plt.axvline(best_intersec, color='red')
        #plt.show()
        return best_intersection_if_curves_overlap(intersec_points, maxima_x_axis, tn_curve_end) #float(min(tn))

#print(find_density_intersection(tp, tn))
