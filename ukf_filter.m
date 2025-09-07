function [x_est, P_est] = ukf_filter(z, u, params, x0, P0, R, dt, meas_schedule)
%% Unscented Kalman Filter with sparse measurements
%
% Handles sparse/intermittent measurements for nonlinear MSD system
%
% Inputs:
%   z      - Measurements (1 x N) with NaN for no measurement
%   u      - Input force (N x 1)
%   params - Filter parameters (m, b, k1, k3, tau, q)
%   x0     - Initial state estimate [3x1]
%   P0     - Initial covariance [3x3]
%   R      - Measurement noise variance
%   dt     - Sample time
%   meas_schedule - Boolean array indicating when measurements available
%
% Outputs:
%   x_est - State estimates (3 x N)
%   P_est - Covariance history (3 x 3 x N)

% Extract parameters
m = params.m;
b = params.b;
k1 = params.k1;
k3 = params.k3;
tau = params.tau;
q = params.q;

% Number of samples
N = length(z);

% State dimension
n = 3;

% UKF parameters - tuned for strong nonlinearity
alpha = 0.01;  % Increased for better nonlinearity capture
beta = 2;
kappa = 3 - n;
lambda = alpha^2 * (n + kappa) - n;

% Weights for mean and covariance
Wm = zeros(2*n+1, 1);
Wc = zeros(2*n+1, 1);
Wm(1) = lambda / (n + lambda);
Wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);
for i = 2:2*n+1
    Wm(i) = 1 / (2*(n + lambda));
    Wc(i) = 1 / (2*(n + lambda));
end

% Initialize
x_est = zeros(3, N);
P_est = zeros(3, 3, N);
x_est(:, 1) = x0;
P_est(:, :, 1) = P0;

% Process noise covariance
G = [0; 0; 1];
Q = G * q * G' * dt;

% Measurement matrix
H = [1, 0, 0];

% UKF loop
for k = 1:N-1
    % Current state estimate
    x_k = x_est(:, k);
    P_k = P_est(:, :, k);
    
    % === PREDICTION STEP ===
    
    % Generate sigma points
    try
        sqrt_P = chol((n + lambda) * P_k)';
    catch
        % If not positive definite, use SVD
        [U,S,~] = svd(P_k);
        sqrt_P = U * sqrt(S) * sqrt(n + lambda);
    end
    
    chi = zeros(n, 2*n+1);
    chi(:, 1) = x_k;
    for i = 1:n
        chi(:, i+1) = x_k + sqrt_P(:, i);
        chi(:, i+n+1) = x_k - sqrt_P(:, i);
    end
    
    % Propagate sigma points through nonlinear dynamics
    chi_pred = zeros(n, 2*n+1);
    for i = 1:2*n+1
        pos = chi(1, i);
        vel = chi(2, i);
        force = chi(3, i);
        
        % Nonlinear state equations
        pos_dot = vel;
        vel_dot = (-k1*pos - k3*pos^3 - b*vel + force + u(k))/m;
        force_dot = -force/tau;
        
        % Euler integration
        chi_pred(1, i) = pos + pos_dot * dt;
        chi_pred(2, i) = vel + vel_dot * dt;
        chi_pred(3, i) = force + force_dot * dt;
    end
    
    % Calculate predicted mean
    x_pred = zeros(n, 1);
    for i = 1:2*n+1
        x_pred = x_pred + Wm(i) * chi_pred(:, i);
    end
    
    % Calculate predicted covariance
    P_pred = Q;
    for i = 1:2*n+1
        dx = chi_pred(:, i) - x_pred;
        P_pred = P_pred + Wc(i) * (dx * dx');
    end
    
    % Ensure symmetry
    P_pred = (P_pred + P_pred') / 2;
    
    % === UPDATE STEP ===
    if k < N && meas_schedule(k+1) && ~isnan(z(k+1))
        % Measurement is available
        
        % Propagate sigma points through measurement model
        z_sigma = zeros(1, 2*n+1);
        for i = 1:2*n+1
            z_sigma(i) = H * chi_pred(:, i);
        end
        
        % Predicted measurement
        z_pred = 0;
        for i = 1:2*n+1
            z_pred = z_pred + Wm(i) * z_sigma(i);
        end
        
        % Innovation
        y = z(k+1) - z_pred;
        
        % Calculate innovation covariance and cross-covariance
        Pyy = R;
        Pxy = zeros(n, 1);
        for i = 1:2*n+1
            dz = z_sigma(i) - z_pred;
            dx = chi_pred(:, i) - x_pred;
            Pyy = Pyy + Wc(i) * (dz * dz');
            Pxy = Pxy + Wc(i) * (dx * dz');
        end
        
        % Kalman gain
        K = Pxy / Pyy;
        
        % State update
        x_est(:, k+1) = x_pred + K * y;
        
        % Covariance update
        P_est(:, :, k+1) = P_pred - K * Pyy * K';
    else
        % No measurement - prediction only
        x_est(:, k+1) = x_pred;
        P_est(:, :, k+1) = P_pred;
    end
    
    % Ensure positive definiteness
    P_est(:, :, k+1) = (P_est(:, :, k+1) + P_est(:, :, k+1)') / 2;
    eigvals = eig(P_est(:, :, k+1));
    if min(eigvals) < 1e-10
        P_est(:, :, k+1) = P_est(:, :, k+1) + eye(3) * 1e-9;
    end
end

end