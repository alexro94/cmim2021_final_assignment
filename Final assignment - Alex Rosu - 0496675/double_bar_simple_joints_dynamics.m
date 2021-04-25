%% Double-bar mechanism with two simple joints - dynamic analysis
% Alex Rosu 
% 25.04.2021

clc; clear;

%% Defining mass properties and coordinates of the system:
% Mass properties of the system:
body(1).m = 1;          % Mass of ground body
body(2).m = 1e-8;       % Mass of slider 1 (simplified)
body(3).m = 1;          % Mass of bar 1
body(4).m = 1;          % Mass of bar 1
body(5).m = 1e-8;       % Mass of slider 2 (simplified)

l = 1;                  % The length of all bodies is assumed to be 1 m

% Mass moment of inertia properties of the system:
body(1).Ic = body(1).m * l^2 / 12;  % Moment of inertia of ground body
body(2).Ic = body(2).m * l^2 / 12;  % Moment of inertia of slider 1
body(3).Ic = body(3).m * l^2 / 12;  % Moment of inertia of bar 1
body(4).Ic = body(4).m * l^2 / 12;  % Moment of inertia of bar 2
body(5).Ic = body(5).m * l^2 / 12;  % Moment of inertia of slider 2

% Generalized coordinates of each bodies in the system
body(1).q = [0; 0; 0];
body(2).q = [0; 1 * sind(30) + 1 * sind(60); 0];
body(3).q = [0.5 * cosd(30)
            (1 * sind(30) + 1 * sind(60)) - sind(30)*0.5
            -deg2rad(30)];
body(4).q = [1 * cosd(30) + 0.5 * cosd(60)
            0.5 * sind(60)
            -deg2rad(60)];
body(5).q = [1 * cosd(30) + 1 * cosd(60)
            0
            0];

%% Get mass matrix
M = mass_matrix(body);

%% Adding a force to effect bar 1 and bar 2 (the force is distributed)
sforce(1).f = [0; -1];
sforce(1).i = 3;
sforce(1).u_i = [0.5; 0];

sforce(2).f = [0; -1];
sforce(2).i = 4;
sforce(2).u_i = [-0.5; 0];

grav = [0; -9.81]; % gravitational acceleration

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

%% Dynamic solution

q0 = system_coordinates(body);
q0p = zeros(size(q0));

y0 = [q0; q0p];

tspan = 0 : 0.01 : 0.5;

[T, Y] = ode45(@(t,y) odefun(t,y,body,revolute,simple,sforce,grav), tspan, y0);

% Plotting solved positions
plot(Y(:, 4), Y(:, 5), ...
    Y(:, 7), Y(:, 8), ...
    Y(:, 10), Y(:, 11), ...
    Y(:, 13), Y(:, 14), ...
    0, 0, '*', 'LineWidth', 2);
grid on;
axis equal;
xlabel('X [m]');
ylabel('Y [m]');
set(gca,'FontSize',14);

% Plotting solved velocities
figure();
plot(Y(:, 19), Y(:, 20), ...
    Y(:, 22), Y(:, 23), ...
    Y(:, 25), Y(:, 26), ...
    Y(:, 28), Y(:, 29), ...
    0, 0, '*', 'LineWidth', 2);
grid on;
axis equal;
xlabel('v_x [m/s]');
ylabel('v_y [m/s]');
set(gca,'FontSize',14)

% A local function which returns an ode-function which is used for solving
% positions and velocities of the studied constrained system
function dydt = odefun(t, y, bodies, revolute, simple, sforce, grav)

    nq = length(y) / 2;

    q = y(1 : nq);
    qp = y(nq + 1 : end);
    M = mass_matrix(bodies);
    Cq = constraint_dq(revolute, simple, [], t, q);

    C = constraint(revolute, simple, [], t, q);

    lhs = [M, Cq'
        Cq, zeros(length(C))];

    Cdot = Cq * qp;

    f = force_vector(grav, sforce, bodies, q);

    g = g_vec_assembly(revolute, simple, [], t, q, qp);

    alpha = 20; beta = 20;
    
    g_mod = g - 2 * alpha * Cdot - beta ^ 2 * C;

    rhs = [f
        g_mod];

    x = lhs\rhs;

    qpp = x(1 : nq);

    dydt = [qp; qpp];

end

