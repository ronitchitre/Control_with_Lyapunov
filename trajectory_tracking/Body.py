import numpy as np
import constants


def hat(v_u):
    return np.array([[0, -1 * v_u[2], v_u[1]],
                     [v_u[2], 0, -1 * v_u[0]],
                     [-1 * v_u[1], v_u[0], 0]])


class Quad:
    def __init__(self, position1, vel1, R1, ang_vel1):
        self.position1 = position1
        self.inertia1 = constants.inertia
        self.vel1 = vel1
        self.mass1 = constants.mass1
        self.d = constants.d
        self.R1 = R1
        self.ang_vel1 = ang_vel1


class Pendulum:
    def __init__(self, l, l_dot, q, ang_vel2):
        self.mass2 = constants.mass2
        self.length_0 = constants.length_0
        self.l = l
        self.q = q
        self.l_dot = l_dot
        self.ang_vel2 = ang_vel2

    def spring_force(self):
        return constants.k_spring * (self.l - constants.length_0) * self.q

    def damper(self):
        return constants.c_damper * self.l_dot * self.q
