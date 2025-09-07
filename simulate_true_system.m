function [x_true, z_meas, F_true] = simulate_true_system(params, x0, u, dt, T, sigma_meas)
%% SIMULATE_TRUE_SYSTEM Simulates the true nonlinear mass-spring-damper system
%
% Inputs:
%   params     - Structure with m, b, k1, k3, q, tau (true parameters)
%   x0         - Initial state [position; velocity; force]
%   u          - Input force vector
%   dt         - Sample time
%   T          - Total simulation time
%   sigma_meas - Measurement noise standard deviation
%
% Outputs:
%   x_true - True state trajectory (3 x N)
%   z_meas - Noisy position measurements (1 x N)
%   F_true - True stochastic force (1 x N)

% Extract parameters
m = params.m;
b = params.b;
k1 = params.k1;
k3 = params.k3;
q = params.q;
tau = params.tau;

% Time vector
t = 0:dt:T;
N = length(t);

% Initialize state [position; velocity; force]
x_true = zeros(3, N);
x_true(:, 1) = x0;
F_true = zeros(1, N);

% Process noise (for stochastic force)
dbeta = randn(N, 1) * sqrt(q * dt);

% Measurement noise
v = sigma_meas * randn(1, N);

% Simulate nonlinear system with Euler integration
for k = 1:N-1
    x = x_true(1, k);  % position
    xdot = x_true(2, k);  % velocity
    F = x_true(3, k);  % stochastic force
    
    % Nonlinear dynamics
    xddot = (-k1*x - k3*x^3 - b*xdot + F + u(k)) / m;
    
    % Update states
    x_true(1, k+1) = x + xdot * dt;
    x_true(2, k+1) = xdot + xddot * dt;
    x_true(3, k+1) = F - F/tau * dt + dbeta(k);  % Ornstein-Uhlenbeck process
    
    F_true(k) = F;
end
F_true(N) = x_true(3, N);

% Generate measurements (position only)
z_meas = x_true(1, :) + v;

end