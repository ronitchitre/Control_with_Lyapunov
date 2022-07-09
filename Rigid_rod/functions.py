import numpy as np
from Body import *


def hat(v_u):
    return np.array([[0, -1 * v_u[2], v_u[1]],
                     [v_u[2], 0, -1 * v_u[0]],
                     [-1 * v_u[1], v_u[0], 0]])


def unit_vec(v):
    if np.linalg.norm(v) != 0:
        return v / np.linalg.norm(v)
    else:
        return np.array([0, 0, 0])


def vec_to_mat(v_u):
    v_u = np.array([v_u])
    return v_u.T.dot(v_u)


def beta(quad, pendulum, k33):
    gamma3 = np.array([0, 0, 1])
    R1 = quad.m_R1
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    R2 = pendulum.m_R2
    inertia2 = pendulum.m_inertia2
    inertia2_inv = np.linalg.inv(inertia2)
    inertia2_inv_spatial = R2.dot(inertia2_inv).dot(R2.T)
    b1 = inertia2_inv_spatial.dot(pendulum.v_ang_mom2)
    b2 = np.cross(R2.dot(gamma3), gamma3)
    b3 = pendulum.v_ang_mom2 + np.cross(a, pendulum.v_mom2)
    if abs(b3[2]) >= 10 ** -3:
        return -1 * np.dot(b1, b2) / (k33 * b3[2]) * np.array([0, 0, 1])
    elif abs(b3[1]) >= 10 ** -3:
        return -1 * np.dot(b1, b2) / (k33 * b3[1]) * np.array([0, 1, 0])
    elif abs(b3[0]) >= 10 ** -3:
        return -1 * np.dot(b1, b2) / (k33 * b3[0]) * np.array([1, 0, 0])
    else:
        return np.array([0, 0, 0])


def beta_force(quad, pendulum, ref1, ref2, k33):
    p_1r = ref1[1]
    p_2r = ref2[1]
    p_1e = p_1r - quad.v_mom1
    p_2e = p_2r - pendulum.v_mom2
    if p_1e[2] + p_2e[2] == 0:
        return np.array([0, 0, 0])
    gamma3 = np.array([0, 0, 1])
    R1 = quad.m_R1
    # a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    R2 = pendulum.m_R2
    inertia2 = pendulum.m_inertia2
    inertia2_inv = np.linalg.inv(inertia2)
    inertia2_inv_spatial = R2.dot(inertia2_inv).dot(R2.T)
    omega2 = inertia2_inv_spatial.dot(pendulum.v_ang_mom2)
    inertia1 = quad.m_inertia1
    inertia1_inv = np.linalg.inv(inertia1)
    inertia1_inv_spatial = R1.dot(inertia1_inv).dot(R1.T)
    omega1 = inertia1_inv_spatial.dot(quad.v_ang_mom1)
    # c1 = np.cross(R1.dot(gamma3), gamma3)
    # c2 = -1 * np.dot(omega1, c1)
    b1 = np.cross(R2.dot(gamma3), gamma3)
    b2 = -1 * np.dot(omega2, b1)
    scal = b2 / (k33 * (p_1e[2] + p_2e[2]))
    # scal += c2 / (k33 * (p_1e[2] + p_2e[2]))
    return scal * np.array([0, 0, 1])


# def beta_new(quad, k33):
#     gamma3 = np.array([0, 0, 1])
#     R1 = quad.m_R1
#     inertia1 = quad.m_inertia1
#     inertia1_inv = np.linalg.inv(inertia1)
#     inertia1_inv_spatial = R1.dot(inertia1_inv).dot(R1.T)
#     omega1 = inertia1_inv_spatial.dot(quad.v_ang_mom1)
#     Hp1 = quad.v_ang_mom1 - np.cross(R1.dot(quad.v_d), quad.v_mom1)
#     temp = np.cross(R1.dot(gamma3), gamma3)
#     if abs(Hp1[2]) >= 10 ** -3:
#         return (np.dot(omega1, temp) / (k33 * Hp1[2])) * np.array([0, 0, 1])
#     elif abs(Hp1[1]) >= 10 ** -3:
#         return (np.dot(omega1, temp) / (k33 * Hp1[1])) * np.array([0, 1, 0])
#     elif abs(Hp1[0]) >= 10 ** -3:
#         return (np.dot(omega1, temp) / (k33 * Hp1[0])) * np.array([1, 0, 0])
#     else:
#         return np.array([0, 0, 0])


def constraint_force(quad, pendulum, control):
    R1 = quad.m_R1
    R2 = pendulum.m_R2
    f_u_1 = control[0]
    torq_u_1 = control[1]
    torq_u_2 = control[2]
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    inertia1_spatial = R1.dot(quad.m_inertia1).dot(R1.T)
    inertia2_spatial = R2.dot(pendulum.m_inertia2).dot(R2.T)
    inertia1_spatial_inv = np.linalg.inv(inertia1_spatial)
    inertia2_spatial_inv = np.linalg.inv(inertia2_spatial)
    v_omega1 = inertia1_spatial_inv.dot(quad.v_ang_mom1)
    v_omega2 = inertia2_spatial_inv.dot(pendulum.v_ang_mom2)
    v_omega1_hat = hat(v_omega1)
    v_omega2_hat = hat(v_omega2)
    b1 = -1 * (f_u_1 + quad.f_e_1) / quad.mass1
    b2 = inertia1_spatial_inv.dot(hat(torq_u_1 - torq_u_2 + np.cross(R1.dot(quad.pos_of_control), f_u_1))).dot(
        R1.dot(quad.v_d))
    b2 = -1 * b2
    b3 = pendulum.f_e_2 / pendulum.mass2
    b4 = -1 * v_omega1_hat.dot(v_omega1_hat).dot(R1).dot(quad.v_d)
    b5 = -1 * inertia2_spatial_inv.dot(hat(torq_u_2)).dot(a)
    b6 = (quad.v_mom1 / quad.mass1) + (v_omega1_hat.dot(R1).dot(quad.v_d)) - (pendulum.v_mom2 / pendulum.mass2)
    b6 = v_omega2_hat.dot(b6)
    b = b1 + b2 + b3 + b4 + b5 + b6

    A1 = np.eye(3) / quad.mass1 + np.eye(3) / pendulum.mass2
    A2 = inertia1_spatial_inv * (np.linalg.norm(quad.v_d)) ** 2
    A3 = inertia2_spatial_inv * (pendulum.length / 2) ** 2
    A4 = -1 * inertia1_spatial_inv.dot(vec_to_mat(R1.dot(quad.v_d)))
    A5 = -1 * inertia2_spatial_inv.dot(vec_to_mat(a))
    A = A1 + A2 + A3 + A4 + A5

    A_inv = np.linalg.inv(A)
    f_c = A_inv.dot(b)
    return f_c


def control(quad, pendulum, ref1, ref2, k33):
    gamma = np.array([0, 0, 1])
    o_1r = ref1[0]
    p_1r = ref1[1]
    o_2r = ref2[0]
    p_2r = ref2[1]
    o_1e = o_1r - quad.v_position1
    o_2e = o_2r - pendulum.v_position2
    p_1e = p_1r - quad.v_mom1
    p_2e = p_2r - pendulum.v_mom2
    R1 = quad.m_R1
    R2 = pendulum.m_R2
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    e3 = np.array([0, 0, 1])
    inertia1_spatial = R1.dot(quad.m_inertia1).dot(R1.T)
    inertia1_spatial_inv = np.linalg.inv(inertia1_spatial)
    omega1 = inertia1_spatial_inv.dot(quad.v_ang_mom1)
    # if p_1e[2] + p_2e[2] != 0:
    #     f1 = np.dot(o_1e, e3) * p_1e[2] / ((p_1e[2] + p_2e[2]) * quad.mass1 * k33) * np.array([0, 0, 1])
    # if p_1e[1] + p_2e[1] != 0:
    #     f1 = np.dot(o_1e, e3) * p_1e[2] / ((p_1e[1] + p_2e[1]) * quad.mass1 * k33) * np.array([0, 1, 0])
    # if p_1e[0] + p_2e[0] != 0:
    #     f1 = np.dot(o_1e, e3) * p_1e[2] / ((p_1e[0] + p_2e[0]) * quad.mass1 * k33) * np.array([1, 0, 0])
    # else:
    #     f1 = np.array([0, 0, 0])
    # f_u_1 = (-1 * quad.f_e_1) + (-1 * pendulum.f_e_2) + f1 + p_1e + p_2e
    H_p_1 = quad.v_ang_mom1 + np.cross(-1 * R1.dot(quad.v_d), quad.v_mom1)
    H_p_2 = pendulum.v_ang_mom2 + np.cross(a, pendulum.v_mom2)
    if p_1e[2] + p_1e[2] != 0:
        f1 = o_1e[2] * p_1e[2] / (k33 * (p_1e[2] + p_2e[2]) * quad.mass1)
    else:
        f1 = 0

    f_u_1 = (-1 * quad.f_e_1[2]) + (-1 * pendulum.f_e_2[2]) + f1
    f_u_1 = np.array([0, 0, f_u_1]) + beta_force(quad, pendulum, ref1, ref2, k33) + p_1e + p_2e

    # f_u_1 = -1 * quad.f_e_1 - pendulum.f_e_2

    # f_u_1 = (-1 * quad.f_e_1) - pendulum.f_e_2 + (100 * abs(o_1e[2]) * (p_1e[2] + p_2e[2])) * np.array([0, 0, 1])
    # f_u_1 += ((quad.mass1 * o_1e[2] + pendulum.mass2 * o_2e[2]) / k33) * np.array([0, 0, 1])

    # H_p_2_unit = unit_vec(H_p_2)
    # H_p_2_scal = 100 * np.linalg.norm(R2.dot(gamma) - gamma) * H_p_2_unit
    #
    torq_u_2 = -1 * H_p_2 - np.cross(a, pendulum.f_e_2)
    c = quad.v_mom1 / quad.mass1 + hat(omega1).dot(R1).dot(quad.v_d) - (pendulum.v_mom2 / pendulum.mass2)
    torq_u_2 += np.cross(c, pendulum.v_mom2)
    # torq_u_2 = -1 * pendulum.v_ang_mom2

    r = R1.dot(quad.pos_of_control - quad.v_d)
    # H_p_1_unit = unit_vec(H_p_1)
    # H_p_1_scal = 100 * np.linalg.norm(R2.dot(gamma) - gamma) * H_p_1_unit
    #
    torq_u_1 = torq_u_2 - np.cross(r, f_u_1) - H_p_1 + np.cross(R1.dot(quad.v_d), quad.f_e_1)
    torq_u_1 += np.cross(hat(omega1).dot(R1).dot(quad.v_d), quad.v_mom1)

    # torq_u_1 = torq_u_2 - np.cross(R1.dot(quad.pos_of_control), f_u_1) - quad.v_ang_mom1
    # torq_u_1 += -1 * np.cross(R1.dot(quad.v_d), pendulum.f_e_2)

    # gamma = np.array([0, 0, 1])
    # torq_u_1 = torq_u_2 - np.cross(quad.pos_of_control, f_u_1) - quad.v_ang_mom1 - np.cross(R1.dot(gamma), gamma)

    # torq_u_1 = torq_u_2 - np.cross(r, f_u_1)
    return [f_u_1, torq_u_1, torq_u_2]


def dynamics(quad, pendulum, ref1, ref2, k33):
    x_dot = np.empty(6, dtype='object')
    R1 = quad.m_R1
    R2 = pendulum.m_R2
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    control_app = control(quad, pendulum, ref1, ref2, k33)
    f_c = constraint_force(quad, pendulum, control_app)
    torq_c_1 = np.cross(R1.dot(quad.v_d), f_c)
    torq_c_2 = np.cross(a, f_c)
    torq_fu = np.cross(R1.dot(quad.pos_of_control), control_app[0])
    o1_dot = quad.v_mom1 / quad.mass1
    x_dot[0] = o1_dot
    omega1 = R1.dot(np.linalg.inv(quad.m_inertia1)).dot(R1.T).dot(quad.v_ang_mom1)
    omega2 = R2.dot(np.linalg.inv(pendulum.m_inertia2)).dot(R2.T).dot(pendulum.v_ang_mom2)
    R1_dot = hat(omega1).dot(R1)
    x_dot[1] = R1_dot
    p1_dot = quad.f_e_1 + control_app[0] + f_c
    x_dot[2] = p1_dot
    ang_mom1_dot = control_app[1] - control_app[2] + torq_c_1 + torq_fu
    x_dot[3] = ang_mom1_dot

    # o2_dot = pendulum.v_mom2 / pendulum.mass2
    # x_dot[4] = o2_dot
    R2_dot = hat(omega2).dot(R2)
    x_dot[4] = R2_dot
    p2_dot = pendulum.f_e_2 - f_c
    # x_dot[5] = p2_dot
    ang_mom2_dot = control_app[2] + torq_c_2
    # a_dot = pendulum.v_mom2 / pendulum.mass2 - quad.v_mom1 / quad.mass1 - hat(omega1).dot(R1).dot(quad.v_d)
    # ang_mom2_dot = H_P2_dot - np.cross(a_dot, pendulum.v_mom2) - np.cross(a, p2_dot)
    x_dot[5] = ang_mom2_dot
    return x_dot


def W_dot(quad, pendulum, ref1, ref2, k33):
    R2 = pendulum.m_R2
    R1 = quad.m_R1
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    control_app = control(quad, pendulum, ref1, ref2, k33)
    o_1r = ref1[0]
    p_1r = ref1[1]
    o_2r = ref2[0]
    p_2r = ref2[1]
    o_1e = o_1r - quad.v_position1
    o_2e = o_2r - pendulum.v_position2
    p_1e = p_1r - quad.v_mom1
    p_2e = p_2r - pendulum.v_mom2
    inertia1_spatial = R1.dot(quad.m_inertia1).dot(R1.T)
    inertia1_spatial_inv = np.linalg.inv(inertia1_spatial)
    omega1 = inertia1_spatial_inv.dot(quad.v_ang_mom1)
    omega2 = R2.dot(np.linalg.inv(pendulum.m_inertia2)).dot(R2.T).dot(pendulum.v_ang_mom2)
    W1 = o_1e[2] * p_1e[2] / quad.mass1
    # W1 = (quad.mass1 * o_1e[2] + pendulum.mass2 * o_2e[2]) * (p_1e[2] + p_2e[2])
    gamma = np.array([0, 0, 1])
    W2 = -1 * np.dot(gamma, hat(omega2).dot(R2).dot(gamma))
    W3 = k33 * np.dot(p_1e + p_2e, -1 * quad.f_e_1 - pendulum.f_e_2 - control_app[0])
    # W3 = k33 * ((p_1e[2] + p_2e[2]) * (-1 * quad.f_e_1[2] - pendulum.f_e_2[2] - control_app[0][2]))
    H_p_1 = quad.v_ang_mom1 + np.cross(-1 * R1.dot(quad.v_d), quad.v_mom1)
    H_p_2 = pendulum.v_ang_mom2 + np.cross(a, pendulum.v_mom2)
    r = quad.m_R1.dot(quad.pos_of_control - quad.v_d)
    W4 = k33 * np.dot(H_p_1, control_app[1] - control_app[2] + np.cross(r, control_app[0]) - np.cross(R1.dot(quad.v_d),
                                                                                                      quad.f_e_1) - np.cross(
        hat(omega1).dot(R1).dot(quad.v_d), quad.v_mom1))
    c = quad.v_mom1 / quad.mass1 + hat(omega1).dot(R1).dot(quad.v_d) - (pendulum.v_mom2 / pendulum.mass2)
    W5 = k33 * np.dot(H_p_2, control_app[2] + np.cross(a, pendulum.f_e_2) - np.cross(c, pendulum.v_mom2))
    ref1 = k33 * (p_1e[2] + p_2e[2]) ** 2
    ref2 = k33 * np.linalg.norm(H_p_1) ** 2
    ref3 = k33 * np.linalg.norm(H_p_2) ** 2
    return W1 + W2 + W3 + W4 + W5


def W(quad, pendulum, ref1, ref2, k33):
    R1 = quad.m_R1
    R2 = pendulum.m_R2
    o_1r = ref1[0]
    o_2r = ref2[0]
    o_1e = o_1r - quad.v_position1
    o_2e = o_2r - pendulum.v_position2
    e3 = np.array([0, 0, 1])
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    W1 = (np.dot(o_1e, e3) ** 2) / 2
    # W1 = (np.dot(quad.mass1 * o_1e + pendulum.mass2 * o_2e, e3) ** 2) / 2
    W2 = (np.linalg.norm(e3 - R2.dot(e3)) ** 2) / 2
    p_1e = ref1[1] - quad.v_mom1
    p_2e = ref2[1] - pendulum.v_mom2
    H_p_1 = quad.v_ang_mom1 + np.cross(-1 * R1.dot(quad.v_d), quad.v_mom1)
    H_p_2 = pendulum.v_ang_mom2 + np.cross(a, pendulum.v_mom2)
    k33 = k33 / 2
    W3 = (k33 * np.linalg.norm(H_p_1) ** 2) + (k33 * np.linalg.norm(H_p_2) ** 2)
    W4 = k33 * np.linalg.norm(p_1e + p_2e) ** 2
    return W1 + W2 + W3 + W4


if __name__ == "__main__":
    o1 = np.array([4, 10, -3])
    d = np.array([0, 0, -0.5])
    inertia = np.eye(3)
    mass1 = 6
    mom1 = np.array([1, -5, 3])
    ang_mom1 = np.array([-10, 1, 3])
    R1 = np.eye(3)
    quad = Quad(o1, R1, mom1, ang_mom1)

    v_mom2 = np.array([5, 9, -3])
    mass2 = 1
    length = 1
    R2 = np.eye(3)
    ang_mom2 = np.array([4, 2, 1])
    c3 = np.array([0, 0, 1])
    o2 = o1 + R1.dot(d) - (R2.dot(c3) * (length / 2))
    pendulum = Pendulum(R2, ang_mom2, quad)

    f_u_1 = np.array([8, -1, 2])
    torq_u_1 = np.array([1, 5, 1])
    torq_u_2 = np.array([-4, -9, -1])

    ref1 = np.array([[1, 1, 1], [1, 1, 0]])
    ref2 = np.array([[1, 1, 0], [2, 1, 0]])
