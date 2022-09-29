import matplotlib.pyplot as plt
import numpy as np
import Body
from solver_rk4method import rk4method_update
import constants
from scipy.spatial.transform import Rotation as R

print('done')


def hat(v_u):
    return np.array([[0, -1 * v_u[2], v_u[1]],
                     [v_u[2], 0, -1 * v_u[0]],
                     [-1 * v_u[1], v_u[0], 0]])


def unit_vec(v):
    if np.linalg.norm(v) != 0:
        return v / np.linalg.norm(v)
    else:
        return np.array([0, 0, 0])


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
        if len(self.pend_traj) == 1:
            l_ddot_est = 0
            omega_dot_est = np.array([0, 0, 0])
        else:
            l_ddot_est = (self.pend_traj[-1].l_dot - self.pend_traj[-2].l_dot) / self.h
            omega_dot_est = (self.pend_traj[-1].ang_vel2 - self.pend_traj[-2].ang_vel2) / self.h
        o_1r_ddot = o_2r_ddot - (l_ddot_est * self.pendulum.q)
        o_1r_ddot += -2 * self.pendulum.l_dot * np.cross(self.pendulum.ang_vel2, self.pendulum.q)
        o_1r_ddot += -1 * self.pendulum.l * np.cross(omega_dot_est, self.pendulum.q)
        o_1r_ddot += self.pendulum.l * (np.linalg.norm(self.pendulum.ang_vel2) ** 2) * self.pendulum.q
        return o_1r, o_1r_dot, o_2r_ddot

    def controller(self, t):
        o_2r = self.pend_des_traj[0](t)
        o_2r_dot = self.pend_des_traj[1](t)
        o_2r_ddot = self.pend_des_traj[2](t)
        o_1r, o_1r_dot, o_1r_ddot = self.estimate_o1r(o_2r, o_2r_dot, o_2r_ddot)
        o1_e = o_1r - self.quad.position1
        o1_dot_e = o_1r_dot - self.quad.vel1
        f_u_1 = (self.quad.mass1 * o_1r_ddot) + (self.quad.mass1 * constants.g)
        f_u_1 += -1 * self.pendulum.spring_force()
        f_u_1 += -1 * self.pendulum.damper()
        f_u_1 += (self.quad.mass1 * constants.Kp_o * o1_e) + (self.quad.mass1 * constants.Kd_o * o1_dot_e)

        torq_u_1 = -1 * np.cross(self.quad.inertia1.dot(self.quad.ang_vel1), self.quad.ang_vel1)
        torq_u_1 += constants.Kd_ang1 * self.quad.ang_vel1

        q_e = 1 - np.dot(self.pendulum.q, np.array([0, 0, -1]))
        q_e_cross = np.cross(self.pendulum.q, np.array([0, 0, -1]))
        torq_u_2 = 2 * self.quad.mass1 * self.pendulum.l * self.pendulum.l_dot * self.pendulum.ang_vel2
        torq_u_2 += self.pendulum.l * np.cross(self.pendulum.q, f_u_1)
        torq_u_2 += self.quad.mass1 * (self.pendulum.l ** 2) * constants.Kp_ang2 * q_e * q_e_cross
        torq_u_2 += -1 * self.quad.mass1 * (self.pendulum.l ** 2) * constants.Kd_ang2 * self.pendulum.ang_vel2
        return [f_u_1, torq_u_1, torq_u_2]

    def dynamics(self, t):
        x_dot = np.empty(8, dtype='object')
        R1 = self.quad.R1
        f_u_1, torq_u_1, torq_u_2 = self.controller(t)
        f_s = self.pendulum.spring_force()
        o1_dot = self.quad.vel1
        x_dot[0] = o1_dot
        v1_dot = -1 * constants.g + (self.pendulum.spring_force() / self.quad.mass1)
        v1_dot += (self.pendulum.damper() / self.quad.mass1) + (f_u_1 / self.quad.mass1)
        x_dot[1] = v1_dot
        R1_dot = R1.dot(hat(self.quad.ang_vel1))
        x_dot[2] = R1_dot
        ang_vel1_dot = np.cross(self.quad.inertia1.dot(self.quad.ang_vel1), self.quad.ang_vel1) + torq_u_1
        ang_vel1_dot = np.linalg.inv(self.quad.inertia1).dot(ang_vel1_dot)
        x_dot[3] = ang_vel1_dot

        l_dot = self.pendulum.l_dot
        x_dot[4] = l_dot
        l_dot_dot = self.pendulum.l * (np.linalg.norm(self.pendulum.ang_vel2) ** 2)
        l_dot_dot += (-1 * constants.mew * np.dot(self.pendulum.q, self.pendulum.spring_force()))
        l_dot_dot += (-1 * constants.mew * np.dot(self.pendulum.q, self.pendulum.damper()))
        l_dot_dot += (-1 * np.dot(f_u_1, self.pendulum.q) / self.quad.mass1)
        x_dot[5] = l_dot_dot
        q_dot = np.cross(self.pendulum.ang_vel2, self.pendulum.q)
        x_dot[6] = q_dot
        ang_vel2_dot = (-2 * self.pendulum.l_dot * self.pendulum.ang_vel2 / self.pendulum.l)
        ang_vel2_dot += -1 * np.cross(self.pendulum.q, f_u_1) / (self.quad.mass1 * self.pendulum.l)
        ang_vel2_dot += torq_u_2 / (self.quad.mass1 * self.pendulum.l * self.pendulum.l)
        x_dot[7] = ang_vel2_dot
        return x_dot


def o_2r(t):
    return np.array([np.sin(t), np.cos(t), 0])


def o_2r_dot(t):
    return np.array([np.cos(t), -1 * np.sin(t), 0])


def o_2r_ddot(t):
    return np.array([-1 * np.sin(t), -1 * np.cos(t), 0])


h = 0.01
desired_traj = [o_2r, o_2r_dot, o_2r_ddot]

o1_initial = np.array([0, 2, constants.length_0 + (constants.mass2 * 9.8 / constants.k_spring)])
v1_initial = np.array([0, 0, 0])
R1_initial = np.eye(3)
ang_vel1_initial = np.array([0, 0, 0])
quad_initial = Body.Quad(o1_initial, v1_initial, R1_initial, ang_vel1_initial)

l_initial = constants.length_0 + (constants.mass2 * 9.8 / constants.k_spring) + 5
l_dot_initial = 0.0
q_initial = unit_vec(np.array([0.5, 0.5, -0.9]))
omega_initial = np.array([0, 0, 0])
pend_initial = Body.Pendulum(l_initial, l_dot_initial, q_initial, omega_initial)

system = Main(quad_initial, pend_initial, h, desired_traj)
t_mesh = np.linspace(0, 20, int(20 / h))


def func(t, y):
    return system.dynamics(t)


for i in range(len(t_mesh) - 1):
    x_quad = system.quad.get_state()
    x_pend = system.pendulum.get_state()
    current_state = np.concatenate([x_quad, x_pend])
    new_state = rk4method_update(func, current_state, t_mesh[i], h)
    quad_new = Body.Quad(new_state[0], new_state[1], new_state[2], new_state[3])
    pend_new = Body.Pendulum(new_state[4], new_state[5], new_state[6], new_state[7])
    system.quad = quad_new
    system.pendulum = pend_new
    system.quad_traj.append(quad_new)
    system.pend_traj.append(pend_new)

o1_list = np.array([i.position1 for i in system.quad_traj])
l_list = np.array([i.l for i in system.pend_traj])
q_list = np.array([i.q for i in system.pend_traj])

o2_list = o1_list + (q_list * l_list[:, np.newaxis])


o_2r_list = np.array([o_2r(i) for i in t_mesh])

plt.plot(t_mesh, o2_list[:, 0])
plt.plot(t_mesh, o2_list[:, 1])
plt.plot(t_mesh, o2_list[:, 2])
# plt.plot(t_mesh, o_2r_list[:, 0], '-')
# plt.plot(t_mesh, o_2r_list[:, 1], '-')
# plt.plot(t_mesh, o_2r_list[:, 2], '-')
plt.title('position of pend')
plt.legend(['x', 'y', 'z'])
plt.xlabel('time (s)')
plt.ylabel('position')
plt.show()

plt.plot(o2_list[:, 0], o2_list[:, 1])
plt.plot(o_2r_list[:, 0], o_2r_list[:, 1], 'r-')
plt.title('locus')
plt.legend(['real', 'desired'])
plt.xlabel('time (s)')
plt.ylabel('position')
plt.show()

plt.plot(t_mesh, o2_list[:, 0])
plt.plot(t_mesh, o_2r_list[:, 0])
plt.title('compare trajectory x')
plt.legend(['actual', 'desired'])
plt.xlabel('time (s)')
plt.ylabel('position')
plt.show()

plt.plot(t_mesh, o2_list[:, 1])
plt.plot(t_mesh, o_2r_list[:, 1])
plt.title('compare trajectory y')
plt.legend(['actual', 'desired'])
plt.xlabel('time (s)')
plt.ylabel('position')
plt.show()

plt.plot(t_mesh, o2_list[:, 2])
plt.plot(t_mesh, o_2r_list[:, 2])
plt.title('compare trajectory z')
plt.legend(['actual', 'desired'])
plt.xlabel('time (s)')
plt.ylabel('position')
plt.show()

plt.plot(t_mesh, o1_list[:, 0])
plt.plot(t_mesh, o1_list[:, 1])
plt.plot(t_mesh, o1_list[:, 2])
plt.title('position of quad')
plt.legend(['x', 'y', 'z'])
plt.xlabel('time (s)')
plt.ylabel('position')
plt.show()

plt.plot(t_mesh, l_list)
plt.title('length')
plt.xlabel('time (s)')
plt.ylabel('length')
plt.show()

plt.plot(t_mesh, q_list[:, 0])
plt.plot(t_mesh, q_list[:, 1])
plt.plot(t_mesh, q_list[:, 2])
plt.title('q')
plt.legend(['x', 'y', 'z'])
plt.show()

# propogate system


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
