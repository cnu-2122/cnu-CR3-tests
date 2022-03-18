import numpy as np
import matplotlib.pyplot as plt

def freefall(u0, v0, k, nmax, dt):
    '''
    Simulates the position and velocity of an object in freefall over time,
    with initial displacement u0 and velocity v0, a drag parameter k (>0),
    using a finite difference method with step size dt,
    for nmax time steps.
    '''
    # Set gravitational constant
    g = 9.81
    
    # Initialise empty vectors to store the result
    U = np.zeros(nmax)
    V = np.zeros(nmax)
    
    # Impose the initial conditions
    U[0] = u0
    V[0] = v0

    # Compute the solution for the remaining nmax-1 time steps
    for n in range(1, nmax):

        # Stop the simulation if the object falls to the ground
        if U[n-1] <= 0:
            break

        U[n] = dt * V[n-1] + U[n-1]
        V[n] = dt * (-k*V[n-1] - g) + V[n-1]
    
    # Return the (possibly truncated) vectors
    return U[:n], V[:n]


# TESTING

def plot_solution(U, V, title):
    '''
    Convenience function for plotting computed
    position and velocity of the object.
    '''
    n = len(U)
    time = np.linspace(0, (n - 1) * dt, n)

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(time, U, 'r.', markersize=2)
    ax[0].set(ylabel='Position (m)')
    ax[1].plot(time, V, 'b.', markersize=2)
    ax[1].set(xlabel='Time (s)', ylabel='Velocity (m/s)')

    fig.suptitle(title)

    # Return the figure and axes objects for further plotting
    return fig, ax


# Suggested test
u0, v0 = 100, 10
k = 0.5
nmax = 500
dt = 0.02

# Compute and plot the solution
U, V = freefall(u0, v0, k, nmax, dt)
plot_solution(U, V, 'Suggested test')
plt.show()


# Further tests -- not exhaustive!

# Set the initial position to zero
U, V = freefall(0, v0, k, nmax, dt)
plot_solution(U, V, 'Initial position set to zero')
plt.show()

# Try different values of initial velocity
v0_values = [-100, 0, 100]
for v0 in v0_values:
    U, V = freefall(u0, v0, k, nmax, dt)
    fig, ax = plot_solution(U, V, f'Different values of v0: v0 = {v0}')
    plt.show()
    

# Set the drag to zero
v0 = 10
U, V = freefall(u0, v0, 0, nmax, dt)
fig, ax = plot_solution(U, V, 'Drag set to zero')

# Compute and plot the exact solution
def exact_solution(u0, v0, time):
    g = 9.81
    U_exact = -0.5*g*time**2 + v0*time + u0
    V_exact = -g*time + v0
    return U_exact, V_exact

n = len(U)
time = np.linspace(0, (n - 1) * dt, n)
U_exact, V_exact = exact_solution(u0, v0, time)
ax[0].plot(time, U_exact, 'k-', label='Exact solution')
ax[1].plot(time, V_exact, 'k-', label='Exact solution')
ax[0].legend()
ax[1].legend()
plt.show()


# Check that increasing dt decreases accuracy
dt_values = [0.1, 0.5, 1, 2]
tmax = 10

# Plot the exact solution (using the previous time vector --
# dt doesn't matter here as long as we have enough points for plotting)
U_exact, V_exact = exact_solution(u0, v0, time)
fig, ax = plt.subplots()
ax.plot(time, U_exact, 'k-', label='Exact solution')

for dt in dt_values:
    # Calculate the max number of steps to cover a simulation time of tmax
    nmax = int(tmax / dt)

    # Compute the approximate solution
    U, V = freefall(u0, v0, 0, nmax, dt)

    # Plot the computed solution
    n = len(U)
    time = np.linspace(0, (n-1)*dt, n)
    ax.plot(time, U, 'x', label=f'dt = {dt}')

ax.set(xlabel='Time (s)', ylabel='Position (m)', title='Different values of dt')
ax.legend()
plt.show()
