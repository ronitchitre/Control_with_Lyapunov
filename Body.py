import numpy as np

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
        self.f_e_1 = self.m_R1.dot(self.mass1*gravity)
        self.pos_of_control = [0, 0.5, 0]


class Pendulum:
    def __init__(self, v_position2, m_R2, v_mom2, v_ang_mom2):
        self.v_position2 = v_position2
        self.m_inertia2 = np.array([[1, 0, 0], [0, 1 / 12, 0], [0, 0, 1 / 12]])
        self.v_mom2 = v_mom2
        self.mass2 = 1
        self.length = 1
        self.m_R2 = m_R2
        self.v_ang_mom2 = v_ang_mom2
        gravity = np.array([0, 0, -9.8])
        self.f_e_2 = self.m_R2.dot(self.mass2 * gravity)
