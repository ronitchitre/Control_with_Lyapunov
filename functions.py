import numpy as np
from Body import Quad, Pendulum


def hat(v_u):
    return np.array([[0, -1 * v_u[2], v_u[1]],
                     [v_u[2], 0, -1 * v_u[0]],
                     [-1 * v_u[1], v_u[0], 0]])


def vec_to_mat(v_u):
    v_u = np.array([v_u])
    return v_u.T.dot(v_u)


def constraint_force(quad, pendulum, f_u_1, torq_u_1, torq_u_2):
    R1 = quad.m_R1
    R2 = pendulum.m_R2
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
    A2 = inertia1_spatial_inv * (np.linalg.norm(R1.dot(quad.v_d))) ** 2
    A3 = inertia2_spatial_inv * (np.linalg.norm(a)) ** 2
    A4 = -1 * inertia1_spatial_inv.dot(vec_to_mat(R1.dot(quad.v_d)))
    A5 = -1 * inertia2_spatial_inv.dot(vec_to_mat(a))
    A = A1 + A2 + A3 + A4 + A5

    A_inv = np.linalg.inv(A)
    f_c = A_inv.dot(b)
    return f_c


if __name__ == "__main__":
    o1 = np.array([0, 0, 0])
    d = np.array([0, 0, -0.5])
    inertia = np.eye(3)
    mass1 = 6
    mom1 = np.array([1, 0, 0])
    ang_mom1 = np.array([0, 1, 0])
    R1 = np.eye(3)
    quad = Quad(o1, inertia, mom1, mass1, d, R1, ang_mom1)

    inertia = np.array([[1, 0, 0], [0, 1 / 12, 0], [0, 0, 1 / 12]])
    v_mom2 = np.array([0.05, 0, 0])
    mass2 = 1
    length = 1
    R2 = np.eye(3)
    ang_mom2 = np.array([1, 0, 0])
    pendulum = Pendulum(inertia, v_mom2, mass2, length, R2, ang_mom2, quad)

    f_u_1 = np.array([0, 1, 1])
    torq_u_1 = np.array([1, 0, 0])
    torq_u_2 = np.array([0, 1, 0])

    print(constraint_force(quad, pendulum, f_u_1, torq_u_1, torq_u_2))
