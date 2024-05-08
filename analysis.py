import os
import numpy as np
import matplotlib.pyplot as plt
#from astropy.timeseries import LombScargle
#from scipy.signal import find_peaks



class Lightcurve:
    def __init__(self, file_path):
        self.file_path = file_path

    def plot_lightcurve(self):
        data = np.loadtxt(self.file_path)
        self.time = data[:, 1]  
        self.magnitude = data[:, 4] 
        plt.plot(self.time, self.magnitude, color='blue')
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        plt.title('Light Curve')
        plt.grid(True)
        plt.show()

    def lombscargle(self, )

    