import numpy as np
import matplotlib.pyplot as plt
from Body import *
from functions import *
from solver_rk4method import rk4method
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
Q2 = np.array([1, 0, 0, 0])
r2 = R.from_euler('xyz', [10, 20, 30], degrees=True)
R2 = r2.as_matrix()
# R2 = np.eye(3)
ang_mom2 = np.array([0., 0., 0.])
pendulum = Pendulum(R2, ang_mom2, quad)

ref1 = np.array([[1, 1, 0], [0, 0, 0]])
ref2 = np.array([[1, 1, 0], [0, 0, 0]])
k33 = 0.001
tf = 20


def system_ode(t, x):
    quad = Quad(x[0], x[1], x[2], x[3])
    pendulum = Pendulum(x[4], x[5], quad)
    x_dot = dynamics(quad, pendulum, ref1, ref2, k33)
    # print(np.linalg.det(x[4]))
    return x_dot


d = np.array([0, 0, -0.5])
length = 1
c3 = np.array([0, 0, 1])
o2 = o1 + R1.dot(d) - (R2.dot(c3) * (length / 2))
initial_cond = [o1, R1, mom1, ang_mom1, R2, ang_mom2]
initial_x = np.empty(6, dtype='object')
for i in range(len(initial_cond)):
    initial_x[i] = initial_cond[i]
time = np.linspace(1, tf, 10000)
y = rk4method(system_ode, initial_x, time, 6)

gamma = np.array([0, 0, 1])
iter = 0



# for i in y:
#     # print(iter)
#     quad = Quad(i[0], i[1], i[2], i[3])
#     pendulum = Pendulum(i[4], i[5], quad)
#     R2 = i[4]
#     R1 = quad.m_R1
#     # omega2 = R2.dot(np.linalg.inv(pendulum.m_inertia2)).dot(R2.T).dot(pendulum.v_ang_mom2)
#     # gamma = np.array([0, 0, 1])
#     # W2 = -1 * np.dot(gamma, hat(omega2).dot(R2).dot(gamma))
#     a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
#     H_p_2 = pendulum.v_ang_mom2 + np.cross(a, pendulum.v_mom2)
#     print(pendulum.v_ang_mom2, H_p_2)
# s = y[1]
# quad_fin = Quad(s[0], s[1], s[2], s[3])
# pend_fin = Pendulum(s[4], s[5], quad_fin)
# print(control(quad_fin, pend_fin, ref1, ref2, k33)[0])
# print(control(quad, pend_fin, ref1, ref2, k33)[1])

# Plots
def plot_quantity(y, k, tf, title):
    N = y.shape[0]
    result = np.zeros([N, 3])
    for i in range(N):
        result[i] = y[i][k]
    result1 = result[:, 0]
    result2 = result[:, 1]
    result3 = result[:, 2]
    time = np.linspace(0, tf, N)
    plt.plot(time, result1)
    plt.plot(time, result2)
    plt.plot(time, result3)
    plt.legend(['x', 'y', 'z'])
    plt.title(title)
    plt.ylabel(title)
    plt.xlabel('time')
    # plt.savefig(f'{title}.jpeg')
    plt.show()


def plot_euler(y, k, tf, title):
    N = y.shape[0]
    result = np.zeros([N, 3])
    for i in range(N):
        rot_mat = y[i][k]
        r = R.from_matrix(rot_mat)
        euler_ang = r.as_euler('zxz', degrees=True)
        result[i] = np.array(euler_ang)
    result1 = result[:, 0]
    result2 = result[:, 1]
    result3 = result[:, 2]
    time = np.linspace(0, tf, N)
    plt.plot(time, result1)
    plt.plot(time, result2)
    plt.plot(time, result3)
    plt.legend(['x', 'y', 'z'])
    plt.title(title)
    plt.ylabel('degree')
    plt.xlabel('time')
    plt.show()


def plot_pend(y, tf, ref1, ref2, k33):
    N = y.shape[0]
    resulta = np.zeros([N, 3])
    resulto2 = np.zeros([N, 3])
    resultp2 = np.zeros([N, 3])
    resultW = np.zeros([N, 1])
    resultW_dot = np.zeros([N, 1])
    resultfu = np.zeros([N, 3])
    resulttorqu1 = np.zeros([N, 3])
    resulttorqu2 = np.zeros([N, 3])
    # resultdet = np.zeros([N, 3])
    resultpi2 = np.zeros([N, 3])
    # resulto1 = np.zeros([N, 3])
    for i in range(N):
        quad = Quad(y[i][0], y[i][1], y[i][2], y[i][3])
        pendulum = Pendulum(y[i][4], y[i][5], quad)
        R1 = quad.m_R1
        R2 = pendulum.m_R2
        a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
        control_app = control(quad, pendulum, ref1, ref2, k33)
        resulta[i] = a
        resulto2[i] = pendulum.v_position2
        # resulto1[i] = quad.v_position1
        resultp2[i] = pendulum.v_mom2
        resultW[i] = W(quad, pendulum, ref1, ref2, k33)
        resultW_dot[i] = W_dot(quad, pendulum, ref1, ref2, k33)
        resultfu[i] = control_app[0]
        resulttorqu1[i] = control_app[1]
        resulttorqu2[i] = control_app[2]
        # resultdet[i] = np.linalg.det(R2)
        resultpi2[i] = pendulum.v_ang_mom2

    # fig = plt.figure(figsize=(10, 7))
    # ax = plt.axes(projection="3d")
    # ax.scatter3D(resulto2[:, 0], resulto2[:, 1], resulto2[:, 2], color='green')
    # ax.scatter3D(resulto1[:, 0], resulto1[:, 1], resulto1[:, 2], color='red')
    # plt.show()

    resulto21 = resulto2[:, 0]
    resulto22 = resulto2[:, 1]
    resulto23 = resulto2[:, 2]
    time = np.linspace(0, tf, N)
    plt.plot(time, resulto21)
    plt.plot(time, resulto22)
    plt.plot(time, resulto23)
    plt.legend(['x', 'y', 'z'])
    plt.title('pos_pendulum')
    plt.ylabel('position')
    plt.xlabel('time')
    # plt.savefig('o2.jpeg')
    plt.show()

    resultp21 = resultp2[:, 0]
    resultp22 = resultp2[:, 1]
    resultp23 = resultp2[:, 2]
    time = np.linspace(0, tf, N)
    plt.plot(time, resultp21)
    plt.plot(time, resultp22)
    plt.plot(time, resultp23)
    plt.legend(['x', 'y', 'z'])
    plt.title('mom_pendulum')
    plt.ylabel('mom_pendulum')
    plt.xlabel('time')
    # plt.savefig('p2.jpeg')
    plt.show()

    resultpi21 = resultpi2[:, 0]
    resultpi22 = resultpi2[:, 1]
    resultpi23 = resultpi2[:, 2]
    time = np.linspace(0, tf, N)
    plt.plot(time, resultpi21)
    plt.plot(time, resultpi22)
    plt.plot(time, resultpi23)
    plt.legend(['x', 'y', 'z'])
    plt.title('ang_mom_pendulum')
    plt.ylabel('ang_mom_pendulum')
    plt.xlabel('time')
    # plt.savefig('ang_mom_pend2.jpeg')
    plt.show()

    resulta1 = resulta[:, 0]
    resulta2 = resulta[:, 1]
    resulta3 = resulta[:, 2]
    time = np.linspace(0, tf, N)
    plt.plot(time, resulta1)
    plt.plot(time, resulta2)
    plt.plot(time, resulta3)
    plt.legend(['x', 'y', 'z'])
    plt.title('o2 - o1 - R1d')
    plt.ylabel('position')
    plt.xlabel('time')
    # plt.savefig('ori.jpeg')
    plt.show()

    plt.plot(time, resultfu[:, 0])
    plt.plot(time, resultfu[:, 1])
    plt.plot(time, resultfu[:, 2])
    plt.title('Control force')
    plt.legend(['x', 'y', 'z'])
    plt.ylabel('force')
    plt.xlabel('time')
    # plt.savefig('fu.jpeg')
    plt.show()

    plt.plot(time, resulttorqu1[:, 0])
    plt.plot(time, resulttorqu1[:, 1])
    plt.plot(time, resulttorqu1[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('Control torq on quad')
    plt.ylabel('torq')
    plt.xlabel('time')
    # plt.savefig('torq1.jpeg')
    plt.show()
    #
    plt.plot(time, resulttorqu2[:, 0])
    plt.plot(time, resulttorqu2[:, 1])
    plt.plot(time, resulttorqu2[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('Control torq on pend')
    plt.ylabel('torq')
    plt.xlabel('time')
    # plt.savefig('torq2.jpeg')
    plt.show()

    plt.plot(time, resultW)
    plt.title('W')
    plt.xlabel('time')
    # plt.savefig('W.jpeg')
    plt.show()
    #
    plt.plot(time, resultW_dot)
    plt.title('W_dot')
    plt.xlabel('time')
    # plt.savefig('W_dot.jpeg')
    plt.show()

    # plt.plot(time, resultdet)
    # plt.title('Det of R2')
    # plt.savefig('det_R2.jpeg')
    # plt.show()


plot_quantity(y, 0, tf, 'position_quad')
plot_quantity(y, 2, tf, 'mom_quad')
plot_quantity(y, 3, tf, 'ang_mom_quad')
plot_pend(y, tf, ref1, ref2, k33)
# plot_euler(y, 1, 20, 'euler_quad')
# plot_euler(y, 4, 20, 'euler_pend')
