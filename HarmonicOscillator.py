import numpy as np
import matplotlib.pyplot as plt

# This code solves the ODE mx'' + bx' + kx = F0*cos(Wd*t)
# m is the mass of the object in kg, b is the damping constant in Ns/m
# k is the spring constant in N/m, F0 is the driving force in N,
# Wd is the frequency of the driving force and x is the position

# Setting up

timeFinal= 300.0   # This is how far the graph will go in seconds
steps = 10000     # Number of steps
dT = timeFinal/steps      # Step length
time = np.linspace(0, timeFinal, steps+1)
# Creates an array with steps+1 values from 0 to timeFinal

# Allocating arrays for velocity and position
vel = np.zeros(steps+1)
pos = np.zeros(steps+1)

# Setting constants and initial values for vel. and pos.
k = 0.1
m = 0.01
vel0 = 0.05
pos0 = 0.0
freqNatural = 10.0**0.5
b = 0.0
F0 = 0.0
Wd = 0.0
vel[0] = vel0    #Sets the initial velocity
pos[0] = pos0    #Sets the initial position



# Numerical solution using Euler's
# Splitting the ODE into two first order ones
# v'(t) = -(k/m)*x(t) - (b/m)*v(t) + (F0/m)*cos(Wd*t)
# x'(t) = v(t)
# Using the definition of the derivative we get
# (v(t+dT) - v(t))/dT on the left side of the first equation
# (x(t+dT) - x(t))/dT on the left side of the second
# In the for loop t and dT will be replaced by i and 1

for i in range(0, steps):
    vel[i+1] = (-k/m)*dT*pos[i] + vel[i]*(1-dT*b/m) + (dT/m)*F0*np.cos(Wd*i)
    pos[i+1] = dT*vel[i] + pos[i]

# Ploting
#----------------
# With no damping
plt.plot(time, pos, 'g-', label='Undampened')

# Damping set to 10% of critical damping
#b = (freqNatural/50)*0.1
b = 0.3

# Using Euler's again to compute new values for new damping
for i in range(0, steps):
    vel[i+1] = (-k/m)*dT*pos[i] + vel[i]*(1-(dT*(b/m))) + (F0*dT/m)*np.cos(Wd*i*dT)
    pos[i+1] = dT*vel[i] + pos[i]

plt.plot(time, pos, 'b-', label = '10% of crit. damping')
plt.plot(time, 0*time, 'k-')      # This plots the x-axis
plt.legend(loc = 'upper right')

#---------------
plt.show()