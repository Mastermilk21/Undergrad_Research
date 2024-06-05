import numpy as np
import matplotlib.pyplot as plt

from boundry_analysis import BoundryLightCurve

filename = r"c:\Users\eitan\code repos\data\extracted tarballs\sector54\lc_2022pix_cleaned.txt"

x_min = 2781
x_max = 2782
cadence = 0.5
analyze = BoundryLightCurve(filename, x_min=x_min, x_max=x_max)
analyze.plot_lightcurve()

analyze.plot_lombscargle(cadence)
analyze.plot_combined()
analyze.phasefold()