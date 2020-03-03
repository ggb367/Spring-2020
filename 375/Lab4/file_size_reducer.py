import time

import numpy as np
import pandas as pd

start = time.time()
file = pd.read_csv('data/part1.csv')
end = time.time()
print("Data Read! It took:", end-start, "seconds")
new_df = []
# count = 0
for index, row in file.iterrows():
    if np.mod(index, 5) == 0:
        new_df.append(row)
print("reduction complete, writing to csv....")
new_df = pd.DataFrame(new_df)
new_df.to_csv('data/part2_reduced.csv')
print("Operation Complete!")
