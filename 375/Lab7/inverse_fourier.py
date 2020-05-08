import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import fft

data = pd.read_csv('Data/Lab7_start0.csv')
freq = np.array(data['freq'])
accel = np.array(data['Accel'])
force = np.array(data['Force'])

accel_td = fft.irfft(accel)
force_td = fft.irfft(force)

time_domain = {'Accel': accel_td, 'Force': force_td}

time_df = pd.DataFrame(data=time_domain)

time_df.to_csv('Data/lab7_TimeDomain.csv')

print('complete')
