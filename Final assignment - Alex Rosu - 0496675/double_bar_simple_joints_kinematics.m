%% Double-bar mechanism with two simple joint - kinematic analysis
% Alex Rosu 
% 25.04.2021

%% Coordinates
% ground
q1 = [0; 0; 0];

% slider 1
q2 = [0
    1 * sind(30) + 1 * sind(60)
    0];

% bar 1
q3 = [0.5 * cosd(30)
    (1 * sind(30) + 1 * sind(60)) - sind(30)*0.5
    -deg2rad(30)];

% bar 2
q4 = [1 * cosd(30) + 0.5 * cosd(60)
    0.5 * sind(60)
    -deg2rad(60)];

% slider 2
q5 = [1 * cosd(30) + 1 * cosd(60)
    0
    0];

q_0 = [q1; q2; q3; q4; q5]; % initial coordinates

%% We need two constraint types (geometric ones)
% - revolute
% - simple constraints

%% Revolute joints
% 1 connects slider 1 and bar 1
revolute(1).i = 2;
revolute(1).j = 3;
revolute(1).s_i = [0; 0];
revolute(1).s_j = [-0.5; 0];

% 2 connects bar 1 and bar 2
revolute(2).i = 3;
revolute(2).j = 4;
revolute(2).s_i = [0.5; 0];
revolute(2).s_j = [-0.5; 0];

% 3 connects bar 2 and slider 2
revolute(3).i = 4;
revolute(3).j = 5;
revolute(3).s_i = [0.5; 0];
revolute(3).s_j = [0; 0];

% % Check revolute joint constraints
% r = revolute(3);
% C_r_i = revolute_joint(r.i, r.j, r.s_i, r.s_j, q_0)

%% Simple constraints

% Three simple joints to fix the ground origin
simple(1).i = 1;
simple(1).k = 1;
simple(1).c_k = 0;

simple(2).i = 1;
simple(2).k = 2;
simple(2).c_k = 0;

simple(3).i = 1;
simple(3).k = 3;
simple(3).c_k = 0;

% slider - use simple joints instead of translational
simple(4).i = 2;
simple(4).k = 1;
simple(4).c_k = 0;

simple(5).i = 2;
simple(5).k = 3;
simple(5).c_k = 0;

simple(6).i = 5;
simple(6).k = 2;
simple(6).c_k = 0;

simple(7).i = 5;
simple(7).k = 3;
simple(7).c_k = 0;

% % check simple constraints
% for s = simple
%     C_s_i = simple_joint(s.i, s.k, s.c_k, q_0)
% end

%% Add some driving constraints
driving(1).i = 2;
driving(1).k = 2;
driving(1).d_k = @(t) sin(t)/4 + 1.2;
driving(1).d_k_t = @(t) cos(t)/4;
driving(1).d_k_tt = @(t) -sin(t)/4;

driving(2).i = 5;
driving(2).k = 1;
driving(2).d_k = @(t) cos(t)/4 + 1.2;
driving(2).d_k_t = @(t) -sin(t)/4;
driving(2).d_k_tt = @(t) -cos(t)/4;

% % Verify
% d = driving(1);
% C_d_i = driving_joint(d.i, d.k, d.d_k, 0, q_0)

% %% Verify constraint function
% clc
% C = constraint(revolute, simple, driving, 0, q_0)

%% Solve constraint equation using fsolve
C_fun = @(t, q) constraint(revolute, simple, driving, t, q);
tic
[T, Q] = position_fsolve(C_fun, 2, q_0, 0.1);
toc
%% Some verification plots
plot(Q(:, 4), Q(:, 5), ...
    Q(:, 7), Q(:, 8), ...
    Q(:, 10), Q(:, 11), ...
    Q(:, 13), Q(:, 14), ...
    0, 0, '*', 'LineWidth', 2);
axis equal

%% Jacobian of our constraints
Cq = constraint_dq(revolute, simple, driving, 0, q_0);

%% Solve constraint equation using NR
C_fun = @(t, q) constraint(revolute, simple, driving, t, q);
Cq_fun = @(t, q) constraint_dq(revolute, simple, driving, t, q);
tic
[T, Q] = position_NR(C_fun, Cq_fun, 2, q_0, 0.1);
toc

%% Some verification plots
plot(Q(:, 4), Q(:, 5), ...
    Q(:, 7), Q(:, 8), ...
    Q(:, 10), Q(:, 11), ...
    Q(:, 13), Q(:, 14), ...
    0, 0, '*', 'LineWidth', 2);
axis equal

%% Verify Ct
Ct = constraint_dt(revolute, simple, driving, 0, q_0);

%% Solve constraint equation using NR for position and velocity
C_fun = @(t, q) constraint(revolute, simple, driving, t, q);
Cq_fun = @(t, q) constraint_dq(revolute, simple, driving, t, q);
Ct_fun = @(t, q) constraint_dt(revolute, simple, driving, t, q);
[T, Q, QP] = pos_vel_NR(C_fun, Cq_fun, Ct_fun, 1, q_0, 0.1);

%% Some verification plots
plot(Q(:, 4), Q(:, 5), ...
    Q(:, 7), Q(:, 8), ...
    Q(:, 10), Q(:, 11), ...
    Q(:, 13), Q(:, 14), ...
    0, 0, '*', 'LineWidth', 2);
grid on;
axis equal;
xlabel('X [m]');
ylabel('Y [m]');
set(gca,'FontSize',14);

%% Some verification plots
plot(QP(:, 4), QP(:, 5), ...
    QP(:, 7), QP(:, 8), ...
    QP(:, 10), QP(:, 11), ...
    QP(:, 13), QP(:, 14), ...
    0, 0, '*', 'LineWidth', 2);
grid on;
axis equal;
xlabel('v_x [m/s]');
ylabel('v_y [m/s]');
set(gca,'FontSize',14)

%% Solve constraint equation using NR for position, velocity and acceleration
C_fun = @(t, q) constraint(revolute, simple, driving, t, q);
Cq_fun = @(t, q) constraint_dq(revolute, simple, driving, t, q);
Ct_fun = @(t, q) constraint_dt(revolute, simple, driving, t, q);
g_vec = @(t, q, qp) g_vec_assembly(revolute, simple, driving, t, q, qp);
[T, Q, QP, QPP] = pos_vel_acc_NR(C_fun, Cq_fun, Ct_fun, g_vec, 1, q_0, 0.1);

%% Some verification plots
plot(QPP(:, 4), QPP(:, 5), ...
    QPP(:, 7), QPP(:, 8), ...
    QPP(:, 10), QPP(:, 11), ...
    QPP(:, 13), QPP(:, 14), ...
    0, 0, '*', 'LineWidth', 2);
grid on;
axis equal;
xlabel('a_x [m/s]');
ylabel('a_y [m/s]');
set(gca,'FontSize',14)