import numpy as np
import matplotlib.pyplot as plt

from boundry_analysis import BoundryLightCurve

filename = r"c:\Users\eitan\code repos\data\extracted tarballs\sector28\lc_2020rtt_cleaned.txt"

x_min = 2075    
x_max = 2078
height = 0.015
analyze = BoundryLightCurve(filename, x_min=x_min, x_max=x_max, height=height)
analyze.plot_lightcurve()
analyze.plot_combined()
#analyze.phasefold()
analyze.plotposter()
#analyze.bounded_lightcurve()
