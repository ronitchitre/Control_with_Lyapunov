import numpy as np

class Quad:
    def __init__(self, v_position1, m_inertia1, v_mom1, mass1, v_d, m_R1, v_ang_mom1):
        self.v_position1 = v_position1
        self.m_inertia1 = m_inertia1
        self.v_mom1 = v_mom1
        self.mass1 = mass1
        self.v_d = v_d
        self.m_R1 = m_R1
        self.v_ang_mom1 = v_ang_mom1
        gravity = np.array([0, 0, -9.8])
        self.f_e_1 = self.m_R1.dot(self.mass1*gravity)
        self.pos_of_control = [0, 0.5, 0]


class Pendulum:
    def __init__(self, m_inertia2, v_mom2, mass2, length, m_R2, v_ang_mom2, quad):
        self.m_inertia2 = m_inertia2
        self.v_mom2 = v_mom2
        self.mass2 = mass2
        self.length = length
        self.m_R2 = m_R2
        self.v_ang_mom2 = v_ang_mom2
        gravity = np.array([0, 0, -9.8])
        self.f_e_2 = self.m_R2.dot(self.mass2 * gravity)
        c3 = np.array([0, 0, 1])
        self.v_position2 = quad.v_position1 + quad.m_R1.dot(quad.v_d) - (length / 2) * m_R2.dot(c3)
