import numpy as np
import matplotlib.pyplot as plt

from boundry_analysis import BoundryLightCurve

sector = "65"
object_name = "2023iwl"

filename = f"c:/Users/eitan/code repos/data/extracted tarballs/sector{sector}/lc_{object_name}_cleaned.txt"


x_min = 3069.6
x_max = 3074
height = 0.015
analyze = BoundryLightCurve(filename, x_min=x_min, x_max=x_max, height=height)
analyze.plot_lightcurve()
analyze.newfastplot()
#analyze.plot_combined()
#analyze.phasefold()
#analyze.plotposter()
#analyze.bounded_lightcurve()
