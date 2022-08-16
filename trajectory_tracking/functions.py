
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


def spring_force(pendulum):
    return constants.k_spring * (pendulum.l - constants.length_0)


def control(quad, pendulum, ref1, ref2):
    o_1r = ref1[0]
    v_1r = ref1[1]
    o_1e = o_1r - quad.position1
    v_1e = v_1r - quad.vel1
    q_r = np.array([0, 0, -1])
    q_e = np.dot(pendulum.q, q_r) - 1
    f_s = spring_force(pendulum)
    u_perp = quad.mass1 * pendulum.l * (np.linalg.norm(pendulum.ang_vel2) ** 2)
    u_perp += -1 * constants.mew * f_s * quad.mass1
    u_perp += quad.mass1 * (pendulum.l - ref2[0])
    u_perp += quad.mass1 * constants.k_l * (pendulum.l_dot - ref2[1])
    u_par = (quad.mass1 * constants.g) + (quad.mass1 * o_1e) + (quad.mass1 * constants.k_o * v_1e)
    f_u_1 = (u_perp * pendulum.q) - np.cross(pendulum.q, np.cross(pendulum.q, u_par))

    torq_u_2 = pendulum.l * np.cross(pendulum.q, f_u_1)
    torq_u_2 += (2 * quad.mass1 * pendulum.l * pendulum.l_dot * pendulum.ang_vel2)
    torq_u_2 += -1 * quad.mass1 * pendulum.l * pendulum.l * q_e * np.cross(pendulum.q, q_r)
    torq_u_2 += -1 * quad.mass1 * pendulum.l * pendulum.l * pendulum.ang_vel2 * constants.k_q

    torq_u_1 = -1 * np.cross(quad.inertia1.dot(quad.ang_vel1), quad.ang_vel1) - (constants.k_omega * quad.ang_vel1)
    return [f_u_1, torq_u_1, torq_u_2]


def dynamics(quad, pendulum, ref1, ref2):
    x_dot = np.empty(8, dtype='object')
    R1 = quad.R1
    control_app = control(quad, pendulum, ref1, ref2)
    f_s = spring_force(pendulum)
    o1_dot = quad.vel1
    x_dot[0] = o1_dot
    v1_dot = -1 * constants.g + constants.k_spring * (pendulum.l - constants.length_0) * pendulum.q
    v1_dot += (control_app[0] / quad.mass1)
    x_dot[1] = v1_dot
    R1_dot = R1.dot(hat(quad.ang_vel1))
    x_dot[2] = R1_dot
    ang_vel1_dot = np.cross(quad.inertia1.dot(quad.ang_vel1), quad.ang_vel1) + control_app[1]
    ang_vel1_dot = np.linalg.inv(quad.inertia1).dot(ang_vel1_dot)
    x_dot[3] = ang_vel1_dot

    l_dot = pendulum.l_dot
    x_dot[4] = l_dot
    l_dot_dot = pendulum.l * (np.linalg.norm(pendulum.ang_vel2) ** 2) - constants.mew * f_s
    l_dot_dot += -1 * (np.dot(control_app[0], pendulum.q)) / quad.mass1
    x_dot[5] = l_dot_dot
    q_dot = np.cross(pendulum.ang_vel2, pendulum.q)
    x_dot[6] = q_dot
    ang_vel2_dot = (-2 * pendulum.l_dot * pendulum.ang_vel2 / pendulum.l)
    ang_vel2_dot += -1 * np.cross(pendulum.q, control_app[0]) / (quad.mass1 * pendulum.l)
    ang_vel2_dot += control_app[2] / (quad.mass1 * pendulum.l * pendulum.l)
    x_dot[7] = ang_vel2_dot
    return x_dot

# def W_dot(quad, pendulum, ref1, ref2, k11, k33):
#     R2 = pendulum.R2
#     R1 = quad.R1
#     control_app = control(quad, pendulum, ref1, ref2, k11, k33)
#     o_1r = ref1[0]
#     p_1r = ref1[1]
#     o_2r = ref2[0]
#     p_2r = ref2[1]
#     o_1e = o_1r - quad.position1
#     o_2e = o_2r - pendulum.position2
#     p_1e = p_1r - quad.mom1
#     p_2e = p_2r - pendulum.mom2
#     # inertia1_spatial = R1.dot(quad.inertia1).dot(R1.T)
#     # inertia1_spatial_inv = np.linalg.inv(inertia1_spatial)
#     # omega1 = inertia1_spatial_inv.dot(quad.ang_mom1)
#     omega2 = R2.dot(np.linalg.inv(pendulum.inertia2)).dot(R2.T).dot(pendulum.ang_mom2)
#     W1 = np.dot(o_1e, p_1e) / quad.mass1
#     W2 = np.dot(o_2e, p_2e) / pendulum.mass2
#     gamma = np.array([0, 0, 1])
#     W3 = -1 * np.dot(omega2, np.cross(R2.dot(gamma), gamma))
#     f_s = spring_force(quad, pendulum, ref1, ref2)
#     torq_s_1 = np.cross(R1.dot(quad.d), f_s)
#     torq_fu = np.cross(R1.dot(quad.pos_of_control), control_app[0])
#     W4 = k11 * np.dot(p_1e, -1 * quad.f_e_1 - f_s - control_app[0])
#     W5 = k11 * np.dot(p_2e, -1 * pendulum.f_e_2 + f_s - control_app[3])
#     W6 = k33 * np.dot(quad.ang_mom1, control_app[1] - control_app[2] + torq_s_1 + torq_fu)
#     W7 = k33 * np.dot(pendulum.ang_mom2, control_app[2])
#     dec1 = k11 * np.linalg.norm(p_1e) ** 2
#     dec2 = k33 * np.linalg.norm(quad.ang_mom1) ** 2
#     dec3 = k33 * np.linalg.norm(pendulum.ang_mom2) ** 2
#     return W1 + W3 + W2 + W4 + W5 + W6 + W7
#
#
# def W(quad, pendulum, ref1, ref2, k11, k33):
#     # R1 = quad.R1
#     R2 = pendulum.R2
#     o_1r = ref1[0]
#     o_2r = ref2[0]
#     p_2r = ref2[1]
#     p_1r = ref1[1]
#     o_1e = o_1r - quad.position1
#     o_2e = o_2r - pendulum.position2
#     p_2e = p_2r - pendulum.mom2
#     p_1e = p_1r - quad.mom1
#     e3 = np.array([0, 0, 1])
#     W1 = (np.linalg.norm(o_1e) ** 2) / 2
#     W2 = (np.linalg.norm(o_2e) ** 2) / 2
#     W3 = (np.linalg.norm(R2.dot(e3) - e3) ** 2) / 2
#     W4 = ((k11 * np.linalg.norm(p_1e) ** 2) / 2) + ((k11 * np.linalg.norm(p_2e) ** 2) / 2)
#     W5 = ((k33 * np.linalg.norm(quad.ang_mom1) ** 2) / 2) + ((k33 * np.linalg.norm(pendulum.ang_mom2) ** 2) / 2)
#     return W1 + W2 + W3 + W4 + W5


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
