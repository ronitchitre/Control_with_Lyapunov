
import matplotlib.pyplot as plt
from functions import *
from solver_rk4method import rk4method
import constants
from scipy.spatial.transform import Rotation as R

print('done')

o1 = np.array([0, 0, 0])
vel1 = np.array([0, 0, 0])
ang_vel1 = np.array([0., 0., 0.])
# R1 = quaternion_rotation_matrix(Q1)
# r1 = R.from_euler('xyz', [20, 10, 0], degrees=True)
# R1 = r1.as_matrix()
R1 = np.eye(3)
quad = Quad(o1, vel1, R1, ang_vel1)

vel2 = np.array([0, 0, 0])
ang_mom2 = np.array([0., 0., 0.])
l = 1
l_dot = 0
q = unit_vec(np.array([0, 0, -1]))
ang_vel2 = np.array([0, 0, 0])
pendulum = Pendulum(quad, l, l_dot, q, ang_vel2)

ref1 = np.array([[3, 2, 10], [0, 0, 0]])
ref2 = np.array([2, 0])
tf = 20


def system_ode(t, x):
    quad = Quad(x[0], x[1], x[2], x[3])
    pendulum = Pendulum(quad, x[4], x[5], x[6], x[7])
    x_dot = dynamics(quad, pendulum, ref1, ref2)
    return x_dot


initial_cond = [o1, vel1, R1, ang_vel1, l, l_dot, q, ang_vel2]
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
    o1_list = np.zeros([N, 3])
    o2_list = np.zeros([N, 3])
    v1_list = np.zeros([N, 3])
    ang_vel1_list = np.zeros([N, 3])
    l_list = np.zeros([N, 1])
    l_dot_list = np.zeros([N, 1])
    q_list = np.zeros([N, 3])
    ang_vel2_list = np.zeros([N, 3])
    f_u_1_list = np.zeros([N, 3])
    torq_u_1_list = np.zeros([N, 3])
    torq_u_2_list = np.zeros([N, 3])
    for i in range(N):
        quad = Quad(y[i][0], y[i][1], y[i][2], y[i][3])
        pendulum = Pendulum(quad, y[i][4], y[i][5], y[i][6], y[i][7])
        o1_list[i] = y[i][0]
        v1_list[i] = y[i][1]
        ang_vel1_list[i] = y[i][3]
        l_list[i] = y[i][4]
        l_dot_list[i] = y[i][5]
        q_list[i] = y[i][6]
        o2_list[i] = o1_list[i] + l_list[i] * q_list[i]
        ang_vel2_list[i] = y[i][7]
        control_app = control(quad, pendulum, ref1, ref2)
        f_u_1_list[i] = control_app[0]
        torq_u_1_list[i] = control_app[1]
        torq_u_2_list[i] = control_app[2]
    plt.plot(time, o1_list[:, 0])
    plt.plot(time, o1_list[:, 1])
    plt.plot(time, o1_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('position of quad')
    plt.xlabel('time')
    plt.ylabel('position')
    # plt.savefig('o1.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, o2_list[:, 0])
    plt.plot(time, o2_list[:, 1])
    plt.plot(time, o2_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('position of pendulum')
    plt.xlabel('time')
    plt.ylabel('position')
    # plt.savefig('o1.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, v1_list[:, 0])
    plt.plot(time, v1_list[:, 1])
    plt.plot(time, v1_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('linear momentum of quad')
    plt.xlabel('time')
    plt.ylabel('momentum')
    # plt.savefig('p1.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, ang_vel1_list[:, 0])
    plt.plot(time, ang_vel1_list[:, 1])
    plt.plot(time, ang_vel1_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('ang momentum of quad')
    plt.xlabel('time')
    plt.ylabel('ang momentum')
    # plt.savefig('pi1.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, l_list)
    plt.title('length')
    plt.xlabel('time')
    plt.ylabel('l')
    # plt.savefig('W.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, l_dot_list)
    plt.title('l_dot')
    plt.xlabel('time')
    plt.ylabel('l_dot')
    # plt.savefig('W.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, q_list[:, 0])
    plt.plot(time, q_list[:, 1])
    plt.plot(time, q_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('q')
    plt.xlabel('time')
    plt.ylabel('q')
    # plt.savefig('p2.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, ang_vel2_list[:, 0])
    plt.plot(time, ang_vel2_list[:, 1])
    plt.plot(time, ang_vel2_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('ang velocity of pendulum')
    plt.xlabel('time')
    plt.ylabel('ang velocity')
    # plt.savefig('pi2.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, f_u_1_list[:, 0])
    plt.plot(time, f_u_1_list[:, 1])
    plt.plot(time, f_u_1_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('control force on quad')
    plt.xlabel('time')
    plt.ylabel('force')
    plt.savefig('force_quad.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, torq_u_1_list[:, 0])
    plt.plot(time, torq_u_1_list[:, 1])
    plt.plot(time, torq_u_1_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('control torque on quad')
    plt.xlabel('time')
    plt.ylabel('torque')
    plt.savefig('torq_quad.jpeg', dpi=1200)
    plt.show()

    plt.plot(time, torq_u_2_list[:, 0])
    plt.plot(time, torq_u_2_list[:, 1])
    plt.plot(time, torq_u_2_list[:, 2])
    plt.legend(['x', 'y', 'z'])
    plt.title('control torque on pendulum')
    plt.xlabel('time')
    plt.ylabel('torque')
    plt.savefig('torq_pend.jpeg', dpi=1200)
    plt.show()


plots(y, time)
