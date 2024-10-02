# ########
# Credits:
# The fit with uncertainties was developed following the book by Peter Young called
# "Everything you wanted to know about data analysis and fitting but were afraid
# to ask" (ISBN 978-3-319-19050-1 978-3-319-19051-8), and more specifically of
# the section 3.1 "Fitting to a straight line".


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


g = 9.8
f = lambda x, a, b: b + a * x
ff = lambda x, a: a * x


df = pd.read_csv("data_buckling.csv")
data = df[df['m'] != 0]

# L = np.array(data['L'])
# t = np.array(data['t'])
# w = np.array(data['w'])
# m = np.array(data['m'])
# err_m = np.array(data['err_m'])
# F_app = m * g
# err_F_app = err_m * g

L_vals = [25, 30, 35, 40]
t_vals = [1, 2, 3]
w_vals = [1, 2]

def length_dependence(log=True):
    # To study the influence on length, we have to obtain all data sets where t and w are held constant and only L varies:
    i = 0
    for tt in t_vals:
        for ww in w_vals:
            mini_df = data[(data['t'] == tt)*(data['w'] == ww)]
            print(mini_df)
            print()
            F = mini_df['m'] * g / 1000
            err_F = mini_df['err_m'] * g / 1000
            L = mini_df['L']
            L_m = L / 100

            if not log:
                if len(mini_df) > 1:
                    i += 1
                    plt.errorbar(L_m, F, yerr=err_F, capsize=3, label=r"$F_{\mathrm{applied}} = f(L)$")
                    plt.xticks([x for x in [0.25, 0.30, 0.35, 0.40] if x*100 in L.values])
                    plt.ylabel(r"$F_{\mathrm{applied}}$ (N)")
                    plt.xlabel(r"length $L$ (m)")
                    plt.legend()
                    plt.title(f"$w = {ww}w_0$, $t = {tt}\\, \\mathrm{{ mm}}$")
                    plt.savefig(f"Fapp_length_nonlog_{i}.pdf")
                    plt.show()
            else:
                if len(mini_df) > 1:
                    i += 1
                    x = np.log(L_m)
                    y = np.log(F)
                    yerr = err_F / F

                    # %% /!\ DO NOT EDIT THIS SECTION
                    n = len(x)  # n is the number of data points
                    p, covm = curve_fit(f, x, y, sigma=yerr)  # p contains the parameters while covm will be used to derive the uncertainties
                    NDF = n - len(p)  # number of degrees of freedom
                    chisq = np.sum(((f(x, *p) - y) / yerr) ** 2) / NDF  # chi-squared of the fit
                    errp = np.sqrt(np.diag(covm) / chisq)  # uncertainty on the parameters of the fit

                    # print the results
                    print('********************')
                    print("slope = %10.4f +/- %7.4f " % (p[0], errp[0]))
                    if len(p) > 1:  # i.e. if the fit is affine and has an intercept.
                        print("intercept = %10.4f +/- %7.4f" % (p[1], errp[1]))
                    print('********************')

                    beta = "%10.2f" % (2 - p[0])
                    err_slope = "%10.2f" % errp[0]

                    err_slope = "%10.2f" % errp[0]
                    if np.isnan(errp[0]):
                        err_slope = ""
                    else:
                        err_slope = "\\pm" + err_slope

                    # plot the results
                    plt.plot(x, f(x, *p), '--', label=f'linear fit: $y= c+(2-\\beta)x$ with $\\beta={beta}{err_slope}$')
                    plt.errorbar(x, y, yerr=yerr, marker='o', capsize=3, linestyle='',
                                 label=r"$\log F_{\mathrm{applied}} = f(\log L)$")
                    plt.ylabel(r"$\log\left(F_{\mathrm{applied}}/1 \mathrm{ N}\right)$")
                    plt.xlabel(r"$\log\left(L/1 \mathrm{ m}\right)$")
                    plt.legend()
                    plt.title(f"$w = {ww}w_0$, $t = {tt}\\, \\mathrm{{ mm}}$")
                    plt.savefig(f"Fapp_length_log_{i}.pdf")
                    plt.show()


def thickness_dependence(log=True):
    # To study the influence on thickness, we have to obtain all data sets where L and w are held constant and only t varies:
    i = 0
    for LL in L_vals:
        for ww in w_vals:
            mini_df = data[(data['L'] == LL)*(data['w'] == ww)]
            print(mini_df)
            print()
            F = mini_df['m'] * g / 1000
            err_F = mini_df['err_m'] * g / 1000
            t = mini_df['t']
            t_m = t / 1000

            if not log:
                if len(mini_df) > 1:
                    i += 1
                    plt.errorbar(t_m, F, yerr=err_F, capsize=3, label=r"$F_{\mathrm{applied}} = f(t)$")
                    plt.xticks([x for x in [0.001, 0.002, 0.003] if x*1000 in t.values])
                    plt.ylabel(r"$F_{\mathrm{applied}}$ (N)")
                    plt.xlabel(r"thickness $t$ (m)")
                    plt.legend()
                    plt.title(f"$w = {ww}w_0$, $L = {LL}\\, \\mathrm{{ cm}}$")
                    plt.savefig(f"Fapp_thickness_nonlog_{i}.pdf")
                    plt.show()
            else:
                if len(mini_df) > 1:
                    i += 1
                    x = np.log(t_m)
                    y = np.log(F)
                    yerr = err_F / F

                    # %% /!\ DO NOT EDIT THIS SECTION
                    n = len(x)  # n is the number of data points
                    p, covm = curve_fit(f, x, y, sigma=yerr)  # p contains the parameters while covm will be used to derive the uncertainties
                    NDF = n - len(p)  # number of degrees of freedom
                    chisq = np.sum(((f(x, *p) - y) / yerr) ** 2) / NDF  # chi-squared of the fit
                    errp = np.sqrt(np.diag(covm) / chisq)  # uncertainty on the parameters of the fit

                    # print the results
                    print('********************')
                    print("slope = %10.4f +/- %7.4f " % (p[0], errp[0]))
                    if len(p) > 1:  # i.e. if the fit is affine and has an intercept.
                        print("intercept = %10.4f +/- %7.4f" % (p[1], errp[1]))
                    print('********************')

                    slope = "%10.2f" % p[0]
                    err_slope = "%10.2f" % errp[0]
                    if np.isnan(errp[0]):
                        err_slope = ""
                    else:
                        err_slope = "\\pm" + err_slope

                    # plot the results
                    plt.plot(x, f(x, *p), '--', label=f'linear fit: $y= c+(\\beta-\\alpha)x$ with $\\beta-\\alpha={slope}{err_slope}$')
                    plt.errorbar(x, y, yerr=yerr, marker='o', capsize=3, linestyle='',
                                 label=r"$\log F_{\mathrm{applied}} = f(\log t)$")
                    plt.ylabel(r"$\log\left(F_{\mathrm{applied}}/1 \mathrm{ N}\right)$")
                    plt.xlabel(r"$\log\left(t/1 \mathrm{ m}\right)$")
                    plt.legend()
                    plt.title(f"$w = {ww}w_0$, $L = {LL}\\, \\mathrm{{ cm}}$")
                    plt.savefig(f"Fapp_thickness_log_{i}.pdf")
                    plt.show()

def width_dependence(log=True):
    # To study the influence on width, we have to obtain all data sets where L and t are held constant and only w varies:
    i = 0
    for LL in L_vals:
        for tt in t_vals:
            mini_df = data[(data['L'] == LL)*(data['t'] == tt)]
            print(mini_df)
            print()
            F = mini_df['m'] * g / 1000
            err_F = mini_df['err_m'] * g / 1000
            w = mini_df['w']
            w_m = w / 100

            if not log:
                if len(mini_df) > 1:
                    i += 1
                    plt.errorbar(w_m, F, yerr=err_F, capsize=3, label=r"$F_{\mathrm{applied}} = f(w)$")
                    plt.xticks([x for x in [0.01, 0.02] if x * 100 in w.values])
                    plt.ylabel(r"$F_{\mathrm{applied}}$ (N)")
                    plt.xlabel(r"width $w$ (units of $w_0$)")
                    plt.legend()
                    plt.title(f"$t = {tt}\\, \\mathrm{{ mm}}$, $L = {LL}\\, \\mathrm{{ cm}}$")
                    plt.savefig(f"Fapp_width_nonlog_{i}.pdf")
                    plt.show()
            else:
                if len(mini_df) > 1:
                    i += 1
                    x = np.log(w_m)
                    y = np.log(F)
                    yerr = err_F / F

                    # %% /!\ DO NOT EDIT THIS SECTION
                    n = len(x)  # n is the number of data points
                    p, covm = curve_fit(f, x, y, sigma=yerr)  # p contains the parameters while covm will be used to derive the uncertainties
                    NDF = n - len(p)  # number of degrees of freedom
                    chisq = np.sum(((f(x, *p) - y) / yerr) ** 2) / NDF  # chi-squared of the fit
                    errp = np.sqrt(np.diag(covm) / chisq)  # uncertainty on the parameters of the fit

                    # print the results
                    print('********************')
                    print("slope = %10.4f +/- %7.4f " % (p[0], errp[0]))
                    if len(p) > 1:  # i.e. if the fit is affine and has an intercept.
                        print("intercept = %10.4f +/- %7.4f" % (p[1], errp[1]))
                    print('********************')

                    slope = "%10.2f" % p[0]
                    err_slope = "%10.2f" % errp[0]
                    if np.isnan(errp[0]):
                        err_slope = ""
                    else:
                        err_slope = "\\pm" + err_slope

                    # plot the results
                    plt.plot(x, f(x, *p), '--', label=f'linear fit: $y= c+\\alpha x$ with $\\alpha={slope}{err_slope}$')
                    plt.errorbar(x, y, yerr=yerr, marker='o', capsize=3, linestyle='',
                                 label=r"$\log F_{\mathrm{applied}} = f(\log w)$")
                    plt.ylabel(r"$\log\left(F_{\mathrm{applied}}/1 \mathrm{ N}\right)$")
                    plt.xlabel(r"$\log\left(w/w_0\right)$")
                    plt.legend()
                    plt.title(f"$t = {tt}\\, \\mathrm{{ mm}}$, $L = {LL}\\, \\mathrm{{ cm}}$")
                    plt.savefig(f"Fapp_width_log_{i}.pdf")
                    plt.show()

def param_dependence():
    L = np.array(data['L']) / 100
    t = np.array(data['t']) / 1000
    w = np.array(data['w']) * 1.5 / 100
    m = np.array(data['m'])
    err_m = np.array(data['err_m'])
    F_app = m * g / 1000
    err_F_app = err_m * g / 1000

    x = w * t**3 / L**2 * 1e6
    y = F_app
    yerr = err_F_app

    # %% /!\ DO NOT EDIT THIS SECTION
    n = len(x)  # n is the number of data points
    p, covm = curve_fit(ff, x, y, sigma=yerr)  # p contains the parameters while covm will be used to derive the uncertainties
    NDF = n - len(p)  # number of degrees of freedom
    chisq = np.sum(((ff(x, *p) - y) / yerr) ** 2) / NDF  # chi-squared of the fit
    errp = np.sqrt(np.diag(covm) / chisq)  # uncertainty on the parameters of the fit

    # print the results
    print('********************')
    print("slope = %10.4f +/- %7.4f " % (p[0], errp[0]))
    if len(p) > 1:  # i.e. if the fit is affine and has an intercept.
        print("intercept = %10.4f +/- %7.4f" % (p[1], errp[1]))
    print('********************')

    slope = "%10.2f" % round(p[0]/1000, 2)
    err_slope = "%10.2f" % round(errp[0]/1000, 2)
    if np.isnan(errp[0]):
        err_slope = ""
    else:
        err_slope = "\\pm" + err_slope

    # plot the results
    # plt.plot(x, f(x, *p), label=f'...')
    plt.plot(x, ff(x, *p), label=r'fit: $F= \alpha\left(\dfrac{wt^3}{L^2}\right)$ with '+f"$\\alpha=({slope}{err_slope})$ GPa")
    plt.errorbar(x, y, yerr=yerr, marker='o', capsize=3, linestyle='',
                 label=r"$F=F\left(\dfrac{wt^3}{L^2}\right)$")
    plt.ylabel(r"$F$ (N)")
    plt.xlabel(r"$wt^3/L^2$ (mm$^2$)")
    plt.legend()
    plt.title(r"Dependence of the force on the parameter $\frac{wt^3}{L^2}$")
    plt.savefig(f"F_app_linear_parameter.pdf")
    plt.show()

length_dependence()
length_dependence(log=False)
thickness_dependence()
thickness_dependence(log=False)
width_dependence()
width_dependence(log=False)

param_dependence()
