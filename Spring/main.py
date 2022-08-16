import numpy as np
import matplotlib.pyplot as plt
from Body import *
from functions import *
from solver_rk4method import rk4method
import constants
from scipy.spatial.transform import Rotation as R

print('done')

o1 = np.array([0, 0, 0])
mom1 = np.array([0, 0, 0])
ang_mom1 = np.array([0., 0., 0.])
# R1 = quaternion_rotation_matrix(Q1)
# r1 = R.from_euler('xyz', [20, 10, 0], degrees=True)
# R1 = r1.as_matrix()
R1 = np.eye(3)
quad = Quad(o1, R1, mom1, ang_mom1)

mom2 = np.array([0, 0, 0])
r2 = R.from_euler('xyz', [40, 20, 30], degrees=True)
R2 = r2.as_matrix()
# R2 = np.eye(3)
c3 = np.array([0, 0, 1])
d = np.array([0, 0, -0.5])
o2 = o1 + R1.dot(d) - (R2.dot(c3) * constants.length / 2)
ang_mom2 = np.array([0., 0., 0.])
pendulum = Pendulum(o2, R2, mom2, ang_mom2)

ref1 = np.array([[3, 2, 10], [0, 0, 0]])
refo2 = ref1[0] + d - (c3 * constants.length / 2) - (constants.mass2 * 9.8 / constants.k_spring) * c3
print(refo2)
refp2 = np.array([0, 0, 0])
ref2 = np.array([refo2, refp2])
k33 = 0.1
k11 = 0.1
tf = 20


def system_ode(t, x):
    quad = Quad(x[0], x[1], x[2], x[3])
    pendulum = Pendulum(x[4], x[5], x[6], x[7])
    x_dot = dynamics(quad, pendulum, ref1, ref2, k33, k11)
    return x_dot


initial_cond = [o1, R1, mom1, ang_mom1, o2, R2, mom2, ang_mom2]
initial_x = np.empty(8, dtype='object')
for i in range(len(initial_cond)):
    initial_x[i] = initial_cond[i]
time = np.linspace(0, tf, 10000)
y = rk4method(system_ode, initial_x, time, 8)

# for i in y:
#     # print(iter)
#     quad = Quad(i[0], i[1], i[2], i[3])
#     pendulum = Pendulum(i[4], i[5], i[6], i[7])
#     # inertia_2_inv = np.linalg.norm(pendulum.inertia2)
#     # inertia_2_inv_spatial = pendulum.R2.dot(inertia_2_inv).dot(pendulum.R2.T)
#     R1 = quad.R1
#     R2 = pendulum.R2
#     c3 = np.array([0, 0, 1])
#     print(pendulum.position2 + R2.dot(c3) - quad.position1 - R1.dot(quad.d))
#     # print(o_2e, W(quad, pendulum, ref1, ref2, k11, k33))


# Plots
def plots(y, time):
    N = y.shape[0]
    o1 = np.zeros([N, 3])
    p1 = np.zeros([N, 3])
    pi1 = np.zeros([N, 3])
    o2 = np.zeros([N, 3])
    p2 = np.zeros([N, 3])
    pi2 = np.zeros([N, 3])
    a = np.zeros([N, 3])
    W_val = np.zeros([N, 1])
    W_dot_val = np.zeros([N, 1])
    f_u_1 = np.zeros([N, 3])
    f_u_2 = np.zeros([N, 3])
    torq_u_1 = np.zeros([N, 3])
    torq_u_2 = np.zeros([N, 3])
    for i in range(N):
        quad = Quad(y[i][0], y[i][1], y[i][2], y[i][3])
        pendulum = Pendulum(y[i][4], y[i][5], y[i][6], y[i][7])
        o1[i] = y[i][0]
        p1[i] = y[i][2]
        pi1[i] = y[i][3]
        o2[i] = y[i][4]
        p2[i] = y[i][6]
        pi2[i] = y[i][7]
        a[i] = pendulum.position2 - quad.position1 - quad.R1.dot(quad.d)
        W_val[i] = W(quad, pendulum, ref1, ref2, k33, k11)
        W_dot_val[i] = W_dot(quad, pendulum, ref1, ref2, k33, k11)
        control_app = control(quad, pendulum, ref1, ref2, k33, k11)
        f_u_1[i] = control_app[0]
        f_u_2[i] = control_app[3]
        torq_u_1[i] = control_app[1]
        torq_u_2[i] = control_app[2]
    plt.plot(time, o1[:, 0])
    plt.plot(time, o1[:, 1])
    plt.plot(time, o1[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('position of quad')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.savefig('o1.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, p1[:, 0])
    plt.plot(time, p1[:, 1])
    plt.plot(time, p1[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('linear momentum of quad')
    plt.xlabel('time')
    plt.ylabel('momentum')
    plt.savefig('p1.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, pi1[:, 0])
    plt.plot(time, pi1[:, 1])
    plt.plot(time, pi1[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('ang momentum of quad')
    plt.xlabel('time')
    plt.ylabel('ang momentum')
    plt.savefig('pi1.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, a[:, 0])
    plt.plot(time, a[:, 1])
    plt.plot(time, a[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('o2 - o1 - R1d')
    plt.xlabel('time')
    plt.ylabel('o2 - o1 - R1d')
    plt.savefig('ori.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, o2[:, 0])
    plt.plot(time, o2[:, 1])
    plt.plot(time, o2[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('position of pendulum')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.savefig('o2.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, p2[:, 0])
    plt.plot(time, p2[:, 1])
    plt.plot(time, p2[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('linear momentum of pendulum')
    plt.xlabel('time')
    plt.ylabel('momentum')
    plt.savefig('p2.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, pi2[:, 0])
    plt.plot(time, pi2[:, 1])
    plt.plot(time, pi2[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('ang momentum of pendulum')
    plt.xlabel('time')
    plt.ylabel('ang momentum')
    plt.savefig('pi2.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, f_u_1[:, 0])
    plt.plot(time, f_u_1[:, 1])
    plt.plot(time, f_u_1[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('control force on quad')
    plt.xlabel('time')
    plt.ylabel('force')
    plt.savefig('force_quad.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, f_u_2[:, 0])
    plt.plot(time, f_u_2[:, 1])
    plt.plot(time, f_u_2[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('control force on pendulum')
    plt.xlabel('time')
    plt.ylabel('force')
    plt.savefig('force_pend.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, torq_u_1[:, 0])
    plt.plot(time, torq_u_1[:, 1])
    plt.plot(time, torq_u_1[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('control torque on quad')
    plt.xlabel('time')
    plt.ylabel('torque')
    plt.savefig('torq_quad.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, torq_u_2[:, 0])
    plt.plot(time, torq_u_2[:, 1])
    plt.plot(time, torq_u_2[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('control torque on pendulum')
    plt.xlabel('time')
    plt.ylabel('torque')
    plt.savefig('torq_pend.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, W_val)
    plt.title('W')
    plt.xlabel('time')
    plt.ylabel('W')
    plt.savefig('W.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, W_dot_val)
    plt.title('W_dot')
    plt.xlabel('time')
    plt.ylabel('W_dot')
    plt.savefig('W_dot.jpeg', dpi=1200)
    plt.show()


plots(y, time)
