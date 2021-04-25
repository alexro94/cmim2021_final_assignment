% Alex Rosu
% 25.04.2021

% A function which solves the kinematics of the system by utilizing
% Newton-Rhapson method and the related theory of multibody systems.

function [T, Q, QP, QPP] = pos_vel_acc_NR(C_fun, Cq_fun, Ct_fun, g_vec, tend, q_0, dt)

    % Initializing
    N_t = floor(round(tend/dt));
    T = linspace(0, N_t*dt, N_t+1)';
    Q = zeros(N_t+1, length(q_0));
    QP = zeros(N_t+1, length(q_0));
    QPP = zeros(N_t+1, length(q_0));

    % [x, iteration_counter] = NR_method(F, J, x, eps)

    % Solving initial position, velocity and acceleration of the system
    [qi, icnt] = NR_method(@(q) C_fun(0, q), @(q) Cq_fun(0, q), q_0, 1e-8);
    Cqi = Cq_fun(0, qi);
    qip = -Cqi\Ct_fun(0, qi);
    g_vec_i = g_vec(0, qi, qip);
    qipp = Cqi\g_vec_i;
    Q(1, :) = qi';
    QP(1, :) = qip';
    QPP(1, :) = qipp';

    % Step equations forward in time
    for n = 1 : N_t
        [qi, icnt] = NR_method(@(q) C_fun(T(n + 1), q), ...
            @(q) Cq_fun(T(n + 1), q), ...
            qi, 1e-8);
        Cqi = Cq_fun(T(n + 1), qi);
        qip = -Cqi\Ct_fun(T(n + 1), qi);
        g_vec_i = g_vec(T(n + 1), qi, qip);
        qipp = Cqi\g_vec_i;
        disp(icnt)
        Q(n + 1, :) = qi';
        QP(n + 1, :) = qip';
        QPP(n + 1, :) = qipp';
    end

end