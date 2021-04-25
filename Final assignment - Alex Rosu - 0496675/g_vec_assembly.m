% Alex Rosu
% 25.04.2021

% A function which assembles the needed g-vector to solve accelerations

function G = g_vec_assembly(revolute, simple, driving, t, q, qp)

r_len = length(revolute);
s_len = length(simple);
d_len = length(driving);

n_constr = 2 * r_len + s_len + d_len;

G = zeros(n_constr, 1);

g_idx = 0;

for r = revolute
    G(g_idx + (1:2)) = revolute_joint_g_vec(r.i, r.j, r.s_i, r.s_j, q, qp);
    g_idx = g_idx + 2;
end

for s = simple
    g_idx = g_idx + 1;
    G(g_idx) = 0;
end

for d = driving
    g_idx = g_idx + 1;
    G(g_idx) = d.d_k_tt(t);
end

end