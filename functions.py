import numpy as np

def hat(v_u):
    return np.array([[0, -1*v_u[2], v_u[1]],
                     [v_u[2], 0, -1*v_u[0]],
                     [-1*v_u[1], v_u[0], 0]])

def constraint_force(quad, pendulum, f_u_1, torq_u_1, torq_u_2):
    R1 = quad.m_R1
    R2 = quad.m_R2
    inertia1_spatial = R1.dot(quad.m_inertia1).dot(R1.T)
    inertia2_spatial = R2.dot(quad.m_inertia2).dot(R2.T)
    inertia1_spatial_inv = np.linalg.inv(inertia1_spatial)
    inertia2_spatial_inv = np.linalg.inv(inertia2_spatial)
    b1 = -1 * (f_u_1 + quad.f_e_1) / quad.mass1
    b2 = inertia1_spatial_inv(hat(torq_u_1 - torq_u_2 + np.cross(quad.pos_of_control, f_u_1))).dot(R1.dot(quad.v_d))
    b2 = -1 * b2

