import numpy as np
from test_object import Plotter
import matplotlib.pyplot as plt
from analysis import Lightcurve

cadence = 0.5
file_path = r"c:\Users\eitan\code repos\data\extracted tarballs\sector26\cam4_ccd1\lc_discovery\lc_2020nww_cleaned.txt"
analyze = Lightcurve(file_path)
#analyze.plot_lightcurve()
analyze.lombscargle(0.5)
#analyze.plot_lombscargle()
analyze.plot_combined()