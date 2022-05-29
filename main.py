import numpy as np
import matplotlib.pyplot as plt
from Body import *
from functions import *
from scipy.integrate import odeint
from solver_rk4method import rk4method

o1 = np.array([4, 10, -3])
mom1 = np.array([0, 0, 0])
ang_mom1 = np.array([0, 0, 0])
R1 = np.eye(3)
# quad = Quad(o1, inertia, mom1, mass1, d, R1, ang_mom1)

mom2 = np.array([0, 0, 0])
R2 = np.eye(3)
ang_mom2 = np.array([0, 0, 0])
# pendulum = Pendulum(o2 ,inertia, v_mom2, mass2, length, R2, ang_mom2)

ref1 = np.array([[1, 1, 1], [1, 1, 0]])
ref2 = np.array([[1, 1, 0], [2, 1, 0]])
k33 = 5


def system_ode(t, x):
    quad = Quad(x[0], x[1], x[2], x[3])
    pendulum = Pendulum(x[4], x[5], x[6], quad)
    x_dot = dynamics(quad, pendulum, ref1, ref2, k33)
    return x_dot


d = np.array([0, 0, -0.5])
length = 1
c3 = np.array([0, 0, 1])
o2 = o1 + R1.dot(d) - (R2.dot(c3) * (length/2))
initial_cond = [o1, R1, mom1, ang_mom1, R2, mom2, ang_mom2]
initial_x = np.empty(7, dtype='object')
for i in range(len(initial_cond)):
    initial_x[i] = initial_cond[i]
time = np.linspace(1, 10, 100)
y = rk4method(system_ode, initial_x, time, 7)
# for i in y:
#     R1 = i[1]
#     print(np.linalg.det(R1))
