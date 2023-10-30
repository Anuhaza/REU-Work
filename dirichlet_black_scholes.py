import numpy as np
import matplotlib.pyplot as plt
from math import *

### the folllwing is another set of parameters for verification
# S = 8.5
# K = 8
# r = 0.02
# sigma = 0.2
# T = (sigma**2)*(1)/2#(sigma**2)*(1/365)/2  # Final time
### C=1.0336,0.504 P=0.3751


S = 90
K = 90
r = 0.3
sigma = 0.5
T = (sigma**2)*(1/365)/2  # Final time
#save C and P here: C=0.977 P=0.98 use the put-call parity

k=2*r/sigma**2

# Parameters
L = 2 # [-L,L]
alpha = 1  # Thermal diffusivity
N = 1000#1000  # Number of grid points
M = 10000#10000  # Number of time steps

# Discretization
dx = 2*L / N  # Grid spacing
dt = T / M  # Time increment

# Initialize the temperature matrix
u = np.zeros((N+1, M+1))

# Set initial conditions
for i in range(N+1):
    u[i,0]=max(e**(0.5*(k+1)*(-L+i*dx))-e**(0.5*(k-1)*(-L+i*dx)),0)

# Set Dirichlet boundary conditions
u[0, :] = 0  # Left boundary
i = 0;
for m in np.arange(0, T + dt, dt): ### ugly here
   u[N,i]=(e**L-e**(-m*k))*e**(0.5*(k-1)*L+0.25*((k+1)**2)*m)
   i += 1

# Finite difference method
for j in range(M):
    for i in range(1, N):
        u[i, j+1] = u[i, j] + alpha * dt / dx**2 * (u[i+1, j] - 2*u[i, j] + u[i-1, j])

# find the index
min_distance=100
for i in range(N+1):
    if abs(-L+i*dx-log(S/K))<min_distance:
        min_distance=abs(-L+i*dx-log(S/K))
        index=i

# print(u[index,M])
C = K*u[index,M]*e**((1-k)*log(S/K)/2-(k+1)**2*T/4)
print(C)

# Plotting
x = np.linspace(0, L, N+1)
t = np.linspace(0, T, M+1)
X, T = np.meshgrid(x, t)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, T, u.T)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('Temperature')
plt.show()
