import matplotlib.pyplot as plt
from functions import *
import Body
from solver_rk4method import rk4method
import constants
from scipy.spatial.transform import Rotation as R

print('done')


class Main:
    def __init__(self, quad, pendulum, h, pend_des_traj):
        self.quad = quad
        self.pendulum = pendulum
        self.h = h
        self.quad_traj = [quad]
        self.pend_traj = [pendulum]
        self.pend_des_traj = pend_des_traj

    def estimate_o1r(self, o_2r, o_2r_dot, o_2r_ddot):
        o_1r = o_2r - (self.pendulum.l * self.pendulum.q)
        o_1r_dot = o_2r_dot - (self.pendulum.l_dot * self.pendulum.q)
        o_1r_dot += -1 * self.pendulum.l * np.cross(self.pendulum.ang_vel2, self.pendulum.q)
        l_ddot_est = (self.pend_traj[-1].l_dot - self.pend_traj[-2].l_dot) / self.h
        omega_dot_est = (self.pend_traj[-1].ang_vel2 - self.pend_traj[-2].ang_vel2) / self.h
        o_1r_ddot = o_2r_ddot - (l_ddot_est * self.pendulum.q)
        o_1r_ddot += -2 * self.pendulum.l_dot * np.cross(self.pendulum.ang_vel2, self.pendulum.q)
        o_1r_ddot += -1 * self.pendulum.l * np.cross(omega_dot_est, self.pendulum.q)
        o_1r_ddot += self.pendulum.l * (np.linalg.norm(self.pendulum.ang_vel2) ** 2) * self.pendulum.q
        return o_1r, o_1r_dot, o_1r_ddot

    def controller(self, t):
        o_2r = self.pend_des_traj[0](t)
        o_2r_dot = self.pend_des_traj[1](t)
        o_2r_ddot = self.pend_des_traj[2](t)
        o1_r, o1r_dot, o1r_ddot = self.estimate_o1r(o_2r, o_2r_dot, o_2r_ddot)
        o1_e = o1_r - self.quad.position1
        o1_dot_e = o1r_dot - self.quad.vel1
        f_u_1 = self.quad.mass1 * o1r_ddot + self.quad.mass1 * constants.g
        f_u_1 += -1 * self.pendulum.spring_force() * self.pendulum.q
        f_u_1 += -1 * self.pendulum.damper() * self.pendulum.q
        f_u_1 += (self.quad.mass1 * constants.Kp_o * o1_e) + (self.quad.mass1 * constants.Kd_o * o1_dot_e)
        f_u_1 += self.quad.mass1 * o1r_ddot

        torq_u_1 = -1 * np.cross(self.quad.inertia1.dot(self.quad.ang_vel1), self.quad.ang_vel1)
        torq_u_1 += constants.Kd_ang1 * self.quad.ang_vel1

        return [f_u_1, torq_u_1]

    def dynamics(self, t):
        x_dot = np.empty(8, dtype='object')
        R1 = self.quad.R1
        control_app = control(self.quad, self.pendulum, t)
        f_s = spring_force(self.pendulum)
        o1_dot = self.quad.vel1
        x_dot[0] = o1_dot
        v1_dot = -1 * constants.g + (self.pendulum.spring_force() / self.quad.mass1)
        v1_dot += (self.pendulum.damper() / self.quad.mass1) + (control_app[0] / self.quad.mass1)
        x_dot[1] = v1_dot
        R1_dot = R1.dot(hat(self.quad.ang_vel1))
        x_dot[2] = R1_dot
        ang_vel1_dot = np.cross(self.quad.inertia1.dot(self.quad.ang_vel1), quad.ang_vel1) + control_app[1]
        ang_vel1_dot = np.linalg.inv(self.quad.inertia1).dot(ang_vel1_dot)
        x_dot[3] = ang_vel1_dot

        l_dot = self.pendulum.l_dot
        x_dot[4] = l_dot
        l_dot_dot = self.pendulum.l * (np.linalg.norm(self.pendulum.ang_vel2) ** 2)
        l_dot_dot += (-1 * constants.mew * np.dot(self.pendulum.q, self.pendulum.spring_force()))
        l_dot_dot += (-1 * constants.mew * np.dot(self.pendulum.q, self.pendulum.damper()))
        l_dot_dot += (-1 * control_app[0] / self.quad.mass1) * self.pendulum.q
        l_dot_dot += -1 * (np.dot(control_app[0], self.pendulum.q)) / self.quad.mass1
        x_dot[5] = l_dot_dot
        q_dot = np.cross(self.pendulum.ang_vel2, self.pendulum.q)
        x_dot[6] = q_dot
        ang_vel2_dot = (-2 * self.pendulum.l_dot * self.pendulum.ang_vel2 / self.pendulum.l)
        ang_vel2_dot += -1 * np.cross(self.pendulum.q, control_app[0]) / (self.quad.mass1 * self.pendulum.l)
        x_dot[7] = ang_vel2_dot
        return x_dot


def o_2r(t):
    return np.array([np.sin(t), np.cos(t), 0])


def o_2r_dot(t):
    return np.array([np.cos(t), -1 * np.sin(t), 0])


def o_2r_ddot(t):
    return np.array([-1 * np.sin(t), -1 * np.cos(t), 0])


def func(t, y):
    quad = Body.Quad(y[0], y[1], y[2], y[3])
    pend = Body.Pendulum(y[4], y[5], y[6], y[7])
    system = Main(quad, pendulum, h)


h = 0.001
desired_traj = [o_2r, o_2r_dot, o_2r_ddot]

o1_initial = np.array([0, 0, 1])
v1_initial = np.array([0, 0, 0])
R1_initial = np.eye(3)
ang_vel1_initial = np.array([0, 0, 0])
quad_initial = Body.Quad(o1_initial, v1_initial, R1_initial, ang_vel1_initial)

l_initial = 1
l_dot_initial = 0
q_initial = np.array([0, 0, -1])
omega_initial = np.array([0, 0, 0])
pend_initial = Body.Pendulum(l_initial, l_dot_initial, q_initial, omega_initial)

system = Main(quad_initial, pend_initial, h, desired_traj)


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
