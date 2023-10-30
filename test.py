import numpy as np
import matplotlib.pyplot as plt
from math import *
import csv
import pandas as pd
from datetime import datetime

def load_column_data_with_filter(csv_file_path, column_name, filter_value):
    column_data = []
    with open(csv_file_path, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['option_type'] == filter_value:
                column_data.append(float(row[column_name]))
    return column_data

file_path = 'data_SPY.csv'

filter_value = 'call'


column_name_to_extract = 'theo'
theo_price = load_column_data_with_filter(file_path, column_name_to_extract, filter_value)

column_name_to_extract = 'underlying_price'
stock_price = load_column_data_with_filter(file_path, column_name_to_extract, filter_value)

column_name_to_extract = 'strike'
strike_price = load_column_data_with_filter(file_path, column_name_to_extract, filter_value)

r = 0.05

column_name_to_extract = 'implied_volatility'
sigma = load_column_data_with_filter(file_path, column_name_to_extract, filter_value)

def calculate_date_differences_with_filter(csv_file_path, date_column1, date_column2, filter_value):
    date_format1 = '%Y-%m-%d'
    date_format2 = '%m/%d/%y'

    date_differences = []
    with open(csv_file_path, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['option_type'] == filter_value:
                date1 = datetime.strptime(row[date_column1].split("T")[0], date_format1)
                date2 = datetime.strptime(row[date_column2], date_format2)
                difference = (date2 - date1).days
                date_differences.append(float(difference - 1))
    return date_differences

# Replace 'date_column1' and 'date_column2' with the actual column names in your CSV file
date_column1 = 'executed_at'
date_column2 = 'expiry'

date_differences_list = calculate_date_differences_with_filter(file_path, date_column1, date_column2, filter_value)

T = []

for i in range(len(date_differences_list)):
    T.append((sigma[i]**2)*(date_differences_list[i]/365)/2)

calc_price = []

for n in range(len(theo_price)):

    print(n)

    k=2*r/sigma[n]**2

    # Parameters
    L = 2 # [-L,L]
    alpha = 1  # Thermal diffusivity
    N = 1000#1000  # Number of grid points
    M = 10000#10000  # Number of time steps

    # Discretization
    dx = 2*L / N  # Grid spacing
    dt = T[n] / M  # Time increment

    # Initialize the temperature matrix
    u = np.zeros((N+1, M+1))

    # Set initial conditions
    for i in range(N+1):
        u[i,0]=max(e**(0.5*(k+1)*(-L+i*dx))-e**(0.5*(k-1)*(-L+i*dx)),0)

    # Set Dirichlet boundary conditions
    u[0, :] = 0  # Left boundary
    i = 0
    for m in np.arange(0, T[n] + dt, dt): ### ugly here
        if not i > 10000:
            u[N,i]=(e**L-e**(-m*k))*e**(0.5*(k-1)*L+0.25*((k+1)**2)*m)
            i += 1

    # Finite difference method
    for j in range(M):
        for i in range(1, N):
            u[i, j+1] = u[i, j] + alpha * dt / dx**2 * (u[i+1, j] - 2*u[i, j] + u[i-1, j])

    # find the index
    min_distance=100
    for i in range(N+1):
        if abs(-L+i*dx-log(stock_price[n]/strike_price[n]))<min_distance:
            min_distance=abs(-L+i*dx-log(stock_price[n]/strike_price[n]))
            index=i

    C = strike_price[n]*u[index,M]*e**((1-k)*log(stock_price[n]/strike_price[n])/2-(k+1)**2*T[n]/4)
    calc_price.append(C)

print("Theoretical price:", theo_price)
print("Calculated price:", calc_price)

x_axis = list(range(0, 26))  # x-axis from 1 to 26

plt.plot(x_axis, theo_price, 'b-o', label='Theoretical Prices')  # 'b-o' means blue line with circle markers
plt.plot(x_axis, calc_price, 'r-o', label='Calculated Prices')   # 'r-o' means red line with circle markers

plt.xlabel('Points')  # Add a label for the x-axis
plt.ylabel('Prices')  # Add a label for the y-axis

plt.title('Theoretical vs Calculated Prices')

plt.legend()

plt.grid(True)  # Show grid lines

plt.show()  # Display the plot.
