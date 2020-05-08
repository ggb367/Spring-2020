import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

calibration_df = pd.read_excel('Data/Lab7_LoadCell_cal.xlsx')

load = np.array(calibration_df['Load, g'])
voltage = np.array(calibration_df['Voltage, V'])

slope, intercept, r_value, p_value, std_err = linregress(voltage, load)
print(slope, intercept, r_value)
