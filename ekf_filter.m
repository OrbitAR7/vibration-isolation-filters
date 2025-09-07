function [x_est, P_est] = ekf_filter(z, u, params, x0, P0, R, dt, meas_schedule)
%% Extended Kalman Filter with sparse measurements
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

% EKF loop
for k = 1:N-1
    % Current state estimate
    x_k = x_est(:, k);
    P_k = P_est(:, :, k);
    
    % Extract states
    pos = x_k(1);
    vel = x_k(2);
    force = x_k(3);
    
    % === PREDICTION STEP ===
    % Nonlinear state transition
    f = [vel; 
         (-k1*pos - k3*pos^3 - b*vel + force + u(k))/m;
         -force/tau];
    
    x_pred = x_k + f * dt;
    
    % Jacobian of state transition (linearization at current estimate)
    F = [0, 1, 0;
         (-k1 - 3*k3*pos^2)/m, -b/m, 1/m;
         0, 0, -1/tau];
    
    % Discrete-time state transition matrix
    Phi = eye(3) + F * dt;
    
    % Predicted covariance
    P_pred = Phi * P_k * Phi' + Q;
    
    % === UPDATE STEP ===
    if k < N && meas_schedule(k+1) && ~isnan(z(k+1))
        % Measurement is available
        
        % Innovation
        y = z(k+1) - H * x_pred;
        
        % Innovation covariance
        S = H * P_pred * H' + R;
        
        % Kalman gain
        K = P_pred * H' / S;
        
        % Updated estimate
        x_est(:, k+1) = x_pred + K * y;
        
        % Updated covariance (Joseph form for numerical stability)
        I_KH = eye(3) - K * H;
        P_est(:, :, k+1) = I_KH * P_pred * I_KH' + K * R * K';
    else
        % No measurement - prediction only
        x_est(:, k+1) = x_pred;
        P_est(:, :, k+1) = P_pred;
    end
    
    % Ensure covariance remains positive definite
    P_est(:, :, k+1) = (P_est(:, :, k+1) + P_est(:, :, k+1)') / 2;
    eigvals = eig(P_est(:, :, k+1));
    if min(eigvals) < 1e-10
        P_est(:, :, k+1) = P_est(:, :, k+1) + eye(3) * 1e-9;
    end
end

end