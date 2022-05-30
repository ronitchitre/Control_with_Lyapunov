import numpy as np
import matplotlib.pyplot as plt
from Body import *
from functions import *
from scipy.integrate import solve_ivp
from solver_rk4method import rk4method

o1 = np.array([0, 0, 0])
mom1 = np.array([0, 0, 0])
ang_mom1 = np.array([0, 0, 0])
Q1 = np.array([1, 0, 0, 0])
# R1 = quaternion_rotation_matrix(Q1)
R1 = np.eye(3)
quad = Quad(o1, R1, mom1, ang_mom1)

mom2 = np.array([0, 0, 0])
Q2 = np.array([1, 0, 0, 0])
# R2 = quaternion_rotation_matrix(Q2)
R2 = np.eye(3)
ang_mom2 = np.array([0, 0, 0])
pendulum = Pendulum(R2, ang_mom2, quad)

ref1 = np.array([[1, 1, 1], [1, 1, 0]])
ref2 = np.array([[1, 1, 0], [1, 1, 0]])
k33 = 10


def system_ode(t, x):
    # R1 = quaternion_rotation_matrix(x[1])
    quad = Quad(x[0], x[1], x[2], x[3])
    # R2 = quaternion_rotation_matrix(x[4])
    pendulum = Pendulum(x[4], x[5], quad)
    x_dot = dynamics(quad, pendulum, ref1, ref2, k33)
    return x_dot


d = np.array([0, 0, -0.5])
length = 1
c3 = np.array([0, 0, 1])
o2 = o1 + R1.dot(d) - (R2.dot(c3) * (length/2))
initial_cond = [o1, R1, mom1, ang_mom1, R2, ang_mom2]
initial_x = np.empty(6, dtype='object')
for i in range(len(initial_cond)):
    initial_x[i] = initial_cond[i]
time = np.linspace(1, 5.01, 100)
y = rk4method(system_ode, initial_x, time, 6)
iter = 0
for i in y:
    # print(iter)
    quad = Quad(i[0], i[1], i[2], i[3])
    pendulum = Pendulum(i[4], i[5], quad)
    # print(W_dot(quad, pendulum, ref1, ref2, k33))
    print(W(quad, pendulum))
    iter +=1
# s = y[11]
# quad_fin = Quad(s[0], s[1], s[2], s[3])
# pend_fin = Pendulum(s[4], s[5], quad_fin)
# print(control(quad, pend_fin, ref1, ref2, k33)[2])
