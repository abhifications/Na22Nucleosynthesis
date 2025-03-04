import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Replace 'your_data_file.csv' with the path to your data file
file_path = '20221030_000000.calo'


df = pd.read_csv(file_path, sep='\s+', decimal=',', header=0, names=['Timestamp', 'RTD_Selection', 'T_Setpoint', 'RTD0', 'RTD1', 'RTD2', 'RTD3', 'RTD4', 'Unknown', 'POWER', 'V', 'I'])

# Adjusting for the possibility of mixed space separation in the timestamp
# Since 'Timestamp' now may not accurately represent separate date and time due to space splitting,
# it's considered as a single initial column for alignment purposes

# Print the first value in the 'POWER' column to check if it's correct
print("First value in 'POWER' column:", df['POWER'].iloc[0])

# Generating a simple sequence of integers based on the DataFrame's length to represent time in seconds
time_in_seconds = range(len(df))
time_labels = [str(pd.Timedelta(seconds=s)) for s in time_in_seconds]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_labels, df['POWER'], marker='o', linestyle='-')
plt.title('Power vs Time in 1-second increment')
plt.xlabel('Time (seconds)')
plt.ylabel('Power')
plt.xticks(rotation=45, ha="right")

plt.gca().set_xticks(np.linspace(0, len(time_labels) - 1, num=10, dtype=int))  # Adjust the number of ticks as necessary
plt.tight_layout()  # Adjust layout to not cut off labels
plt.show()

avg_pwr = sum(df['POWER']/len(df))
print(avg_pwr)