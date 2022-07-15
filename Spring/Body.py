import numpy as np
import constants


def hat(v_u):
    return np.array([[0, -1 * v_u[2], v_u[1]],
                     [v_u[2], 0, -1 * v_u[0]],
                     [-1 * v_u[1], v_u[0], 0]])


class Quad:
    def __init__(self, position1, R1, mom1, ang_mom1):
        self.position1 = position1
        self.inertia1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.mom1 = mom1
        self.mass1 = constants.mass1
        self.d = constants.d
        self.R1 = R1
        self.ang_mom1 = ang_mom1
        gravity = np.array([0, 0, -9.8])
        self.f_e_1 = self.mass1 * gravity
        self.pos_of_control = [0, 0, 0]


class Pendulum:
    def __init__(self, position2, R2, mom2, ang_mom2):
        self.inertia2 = np.array([[1 / 12, 0, 0], [0, 1 / 12, 0], [0, 0, 1 / 120]])
        self.mass2 = constants.mass2
        self.length = constants.length
        self.ks = constants.k_spring
        self.position2 = position2
        self.mom2 = mom2
        self.R2 = R2
        self.ang_mom2 = ang_mom2
        gravity = np.array([0, 0, -9.8])
        self.f_e_2 = self.mass2 * gravity
