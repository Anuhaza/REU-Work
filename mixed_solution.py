import numpy as np
import matplotlib.pyplot as plt
import math

# Parameters
L = 1.0  # Length of the domain
T = 0.1  # Final time
alpha = 0.01 * 300  # Thermal diffusivity
N = 100  # Number of grid points
M = 10000  # Number of time steps

# Discretization
dx = L / N  # Grid spacing
dt = T / M  # Time increment

# Initialize the temperature matrix
u = np.zeros((N+1, M+1))

# Set initial conditions
for i in range(N+1):
    u[i, 0] = math.sin(3 * i * dx)  # Initial temperature

# Set boundary conditions
A = math.sin(3 * 0.0)  # Temperature at left boundary
B = math.cos(3 * 1.0)*3 # Heat flux at right boundary

# Apply boundary conditions
u[0, :] = A  # Dirichlet boundary condition
u[N, :] = u[N-1, :] + B * dx  # Mixed boundary condition

# Finite difference method
for j in range(M):
    for i in range(1, N):
        u[i, j+1] = u[i, j] + alpha * dt / dx**2 * (u[i+1, j] - 2*u[i, j] + u[i-1, j])

# Plotting
x = np.linspace(0, L, N+1)
t = np.linspace(0, T, M+1)

X, T = np.meshgrid(x, t)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, T, u.T, cmap='coolwarm')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('Temperature')
plt.show()
