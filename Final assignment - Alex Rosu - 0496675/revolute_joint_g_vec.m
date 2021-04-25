% Alex Rosu
% 25.04.2021

% Supplementary function for assembling the part of g-vector related to
% revolute joints

function g_r = revolute_joint_g_vec(i, j, s_i, s_j, q, qp)

idx_i = body_idx(i);
phi_i = q(idx_i(3));
phi_i_dot = qp(idx_i(3));
idx_j = body_idx(j);
phi_j = q(idx_j(3));
phi_j_dot = qp(idx_j(3));

g_r = rot(phi_i) * s_i * phi_i_dot^2 - rot(phi_j) * s_j * phi_j_dot^2;
