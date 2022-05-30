import numpy as np
from Body import *


def hat(v_u):
    return np.array([[0, -1 * v_u[2], v_u[1]],
                     [v_u[2], 0, -1 * v_u[0]],
                     [-1 * v_u[1], v_u[0], 0]])


def vec_to_mat(v_u):
    v_u = np.array([v_u])
    return v_u.T.dot(v_u)


def rotation_angles(matrix, order):
    """
    input
        matrix = 3x3 rotation matrix (numpy array)
        oreder(str) = rotation order of x, y, z : e.g, rotation XZY -- 'xzy'
    output
        theta1, theta2, theta3 = rotation angles in rotation order
    """
    r11, r12, r13 = matrix[0]
    r21, r22, r23 = matrix[1]
    r31, r32, r33 = matrix[2]

    if order == 'xzx':
        theta1 = np.arctan(r31 / r21)
        theta2 = np.arctan(r21 / (r11 * np.cos(theta1)))
        theta3 = np.arctan(-r13 / r12)

    elif order == 'xyx':
        theta1 = np.arctan(-r21 / r31)
        theta2 = np.arctan(-r31 / (r11 * np.cos(theta1)))
        theta3 = np.arctan(r12 / r13)

    elif order == 'yxy':
        theta1 = np.arctan(r12 / r32)
        theta2 = np.arctan(r32 / (r22 * np.cos(theta1)))
        theta3 = np.arctan(-r21 / r23)

    elif order == 'yzy':
        theta1 = np.arctan(-r32 / r12)
        theta2 = np.arctan(-r12 / (r22 * np.cos(theta1)))
        theta3 = np.arctan(r23 / r21)

    elif order == 'zyz':
        theta1 = np.arctan(r23 / r13)
        theta2 = np.arctan(r13 / (r33 * np.cos(theta1)))
        theta3 = np.arctan(-r32 / r31)

    elif order == 'zxz':
        theta1 = np.arctan(-r13 / r23)
        theta2 = np.arctan(-r23 / (r33 * np.cos(theta1)))
        theta3 = np.arctan(r31 / r32)

    elif order == 'xzy':
        theta1 = np.arctan(r32 / r22)
        theta2 = np.arctan(-r12 * np.cos(theta1) / r22)
        theta3 = np.arctan(r13 / r11)

    elif order == 'xyz':
        theta1 = np.arctan(-r23 / r33)
        theta2 = np.arctan(r13 * np.cos(theta1) / r33)
        theta3 = np.arctan(-r12 / r11)

    elif order == 'yxz':
        theta1 = np.arctan(r13 / r33)
        theta2 = np.arctan(-r23 * np.cos(theta1) / r33)
        theta3 = np.arctan(r21 / r22)

    elif order == 'yzx':
        theta1 = np.arctan(-r31 / r11)
        theta2 = np.arctan(r21 * np.cos(theta1) / r11)
        theta3 = np.arctan(-r23 / r22)

    elif order == 'zyx':
        theta1 = np.arctan(r21 / r11)
        theta2 = np.arctan(-r31 * np.cos(theta1) / r11)
        theta3 = np.arctan(r32 / r33)

    elif order == 'zxy':
        theta1 = np.arctan(-r12 / r22)
        theta2 = np.arctan(r32 * np.cos(theta1) / r22)
        theta3 = np.arctan(-r31 / r33)

    theta1 = theta1 * 180 / np.pi
    theta2 = theta2 * 180 / np.pi
    theta3 = theta3 * 180 / np.pi

    return (theta1, theta2, theta3)


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
    if abs(b3[2]) >= 10 ** -5:
        return -1 * np.dot(b1, b2) / (k33 * b3[2]) * np.array([0, 0, 1])
    elif abs(b3[1]) >= 10 ** -5:
        return -1 * np.dot(b1, b2) / (k33 * b3[1]) * np.array([0, 1, 0])
    if abs(b3[1]) >= 10 ** -5:
        return -1 * np.dot(b1, b2) / (k33 * b3[0]) * np.array([1, 0, 0])
    else:
        return 0


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
    b2 = inertia1_spatial_inv.dot(hat(torq_u_1 - torq_u_2 + np.cross(quad.pos_of_control, f_u_1))).dot(R1.dot(quad.v_d))
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
    o_1r = ref1[0]
    p_1r = ref1[1]
    o_2r = ref2[0]
    p_2r = ref2[1]
    o_1e = o_1r - quad.v_position1
    p_1e = p_1r - quad.v_mom1
    p_2e = p_2r - pendulum.v_mom2
    R1 = quad.m_R1
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    e3 = np.array([0, 0, 1])
    if abs(p_1e[2] + p_2e[2]) >= 10 ** -3:
        f1 = np.dot(o_1e, e3) * p_1e[2] / ((p_1e[2] + p_2e[2]) * quad.mass1) * np.array([0, 0, 1])
        f1 = np.dot(o_1e, e3) * p_1e[2] / ((p_1e[2] + p_2e[2]) * quad.mass1) * np.array([0, 0, 1])
    if abs(p_1e[1] + p_2e[1]) >= 10 ** -3:
        f1 = np.dot(o_1e, e3) * p_1e[2] / ((p_1e[1] + p_2e[1]) * quad.mass1) * np.array([0, 1, 0])
    if abs(p_1e[0] + p_2e[0]) >= 10 ** -3:
        f1 = np.dot(o_1e, e3) * p_1e[2] / ((p_1e[0] + p_2e[0]) * quad.mass1) * np.array([1, 0, 0])
    else:
        f1 = np.array([0, 0, 0])
    f_u_1 = (-1 * quad.f_e_1) + (-1 * pendulum.f_e_2) + f1 + p_1e + p_2e
    H_p_2 = pendulum.v_ang_mom2 + np.cross(a, pendulum.v_mom2)
    torq_u_2 = beta(quad, pendulum, k33) - H_p_2
    H_p_1 = quad.v_ang_mom1 + np.cross(-1 * R1.dot(quad.v_d), quad.v_mom1)
    r = quad.m_R1.dot(quad.pos_of_control - quad.v_d)
    torq_u_1 = torq_u_2 - np.cross(r, f_u_1) - H_p_1
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
    # p2_dot = pendulum.f_e_2 - f_c
    # x_dot[5] = p2_dot
    ang_mom2_dot = control_app[2] + torq_c_2
    x_dot[5] = ang_mom2_dot
    return x_dot


def W(quad, pendulum):
    R1 = quad.m_R1
    R2 = pendulum.m_R2
    o_1r = np.array([1, 1, 1])
    o_1e = o_1r - quad.v_position1
    e3 = np.array([0, 0, 1])
    a = pendulum.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
    W1 = (np.dot(o_1e, e3) ** 2) / 2
    W2 = (np.linalg.norm(e3 + R2.dot(e3)) ** 2) / 2
    p_1e = np.array([1, 1, 0]) - quad.v_mom1
    p_2e = np.array([1, 1, 0]) - pendulum.v_mom2
    H_p_1 = quad.v_ang_mom1 + np.cross(-1 * quad.v_d, quad.v_mom1)
    H_p_2 = pendulum.v_ang_mom2 + np.cross(a, pendulum.v_mom2)
    W3 = 5 * np.linalg.norm(p_1e + p_2e) ** 2 + 5 * np.linalg.norm(H_p_1) ** 2 + 5 * np.linalg.norm(H_p_2) ** 2
    return W1 + W2 + W3


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
    pendulum = Pendulum(o2, R2, v_mom2, ang_mom2)

    f_u_1 = np.array([8, -1, 2])
    torq_u_1 = np.array([1, 5, 1])
    torq_u_2 = np.array([-4, -9, -1])

    ref1 = np.array([[1, 1, 1], [1, 1, 0]])
    ref2 = np.array([[1, 1, 0], [2, 1, 0]])
