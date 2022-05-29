import numpy as np


def hat(v_u):
    return np.array([[0, -1 * v_u[2], v_u[1]],
                     [v_u[2], 0, -1 * v_u[0]],
                     [-1 * v_u[1], v_u[0], 0]])


class Quad:
    def __init__(self, v_position1, m_R1, v_mom1, v_ang_mom1):
        self.v_position1 = v_position1
        self.m_inertia1 = np.eye(3)
        self.v_mom1 = v_mom1
        self.mass1 = 6
        self.v_d = np.array([0, 0, -0.5])
        self.m_R1 = m_R1
        self.v_ang_mom1 = v_ang_mom1
        gravity = np.array([0, 0, -9.8])
        self.f_e_1 = self.m_R1.dot(self.mass1 * gravity)
        self.pos_of_control = [0, 0.5, 0]


class Pendulum:
    def __init__(self, m_R2, v_ang_mom2, quad):
        self.m_inertia2 = np.array([[1 / 12, 0, 0], [0, 1 / 12, 0], [0, 0, 1 / 120]])
        self.mass2 = 1
        self.length = 1
        self.m_R2 = m_R2
        self.v_ang_mom2 = v_ang_mom2
        gravity = np.array([0, 0, -9.8])
        self.f_e_2 = self.m_R2.dot(self.mass2 * gravity)
        c3 = np.array([0, 0, 1])
        self.v_position2 = quad.v_position1 + quad.m_R1.dot(quad.v_d) - (self.m_R2.dot(c3) * (self.length / 2))
        R1 = quad.m_R1
        a = self.v_position2 - quad.v_position1 - R1.dot(quad.v_d)
        inertia1_spatial = R1.dot(quad.m_inertia1).dot(R1.T)
        inertia1_spatial_inv = np.linalg.inv(inertia1_spatial)
        inertia2_spatial = m_R2.dot(self.m_inertia2).dot(m_R2.T)
        inertia2_spatial_inv = np.linalg.inv(inertia2_spatial)
        omega1 = inertia1_spatial_inv.dot(quad.v_ang_mom1)
        omega2 = inertia2_spatial_inv.dot(self.v_ang_mom2)
        self.v_mom2 = quad.v_mom1 / quad.mass1 + hat(omega1).dot(R1).dot(quad.v_d) + hat(omega2).dot(a)
