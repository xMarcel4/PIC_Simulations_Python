# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 23:49:15 2024

@author: marce
"""

import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Set the path to your project directory
project_dir = r'C:\Users\Marcel\OneDrive\Desktop\Uni\Semester 8\Studienprojekt\GeerdeteBox_Debugging\\'

# Define the path to the CSV file
csv_file_path = os.path.join(project_dir, 'runtime_diags.csv')

# Load the CSV data into a DataFrame
df = pd.read_csv(csv_file_path, index_col=1)
print(df)


# Plot the data
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(df.loc[:, 'total_E'], label='E_total')
ax1.plot(df.loc[:, 'KE.e-'], label='KE.e-')
ax1.plot(df.loc[:, 'KE.O+'], label='KE.O+')
ax1.plot(df.loc[:, 'PE'], label='PE')

# Label the axes
ax1.set_xlabel('Zeit [s]', fontsize=15)
ax1.set_ylabel('Energie [J]', fontsize=15)

# Set axis limits and scale
#ax1.set_xlim(0, 2e-6)
#ax1.set_ylim(0, 8e-11)

# Format the y-axis ticks
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2e'))

# Add a legend
plt.legend(loc=1)

# Show the plot
plt.show()
