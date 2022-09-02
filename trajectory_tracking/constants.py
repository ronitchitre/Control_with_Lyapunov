import numpy as np

mass1 = 6
mass2 = 1
mew = (mass1 + mass2) / (mass1 * mass2)
length_0 = 1
d = np.array([0, 0, 0])
k_spring = 10
c_damper = 10
g = np.array([0, 0, 9.8])
inertia = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
Kp_o = 1
Kd_o = 1
Kd_ang1 = 1
