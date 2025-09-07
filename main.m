%% Vibration Isolation State Estimation for Sensor Pod
% Main script for comparing EKF, UKF, and SLF performance
% Application: Inspection robot sensor pod vibration isolation with nonlinear stiffness

clear; clc; close all;

%% System Parameters (Truth)
params_true.m = 1.0;          % Mass [kg]
params_true.b = 5;            % Damping [N·s/m] - light for sustained oscillation
params_true.k1 = 400;         % Linear stiffness [N/m] - softened spring
params_true.k3 = 4e5;         % Cubic stiffness [N/m³] - strong nonlinearity
params_true.tau = 0.4;        % Force correlation time [s] - fairly fast OU force
params_true.sigma = 2.5;      % Noise strength [N] - energetic disturbance
params_true.q = 2*params_true.sigma^2/params_true.tau;  % Process noise intensity

%% Filter Parameters (with significant mismatch)
mismatch = 0.3;  % 30% parameter uncertainty

params_filter.m   = params_true.m   * (1 - mismatch*0.3);   % slightly underestimated mass
params_filter.b   = params_true.b   * (1 + mismatch);       % damping overestimated
params_filter.k1  = params_true.k1  * (1 - mismatch*0.5);   % linear stiffness underestimated
params_filter.k3  = params_true.k3  * (1 + mismatch*1.5);   % big error in nonlinear term!
params_filter.tau = params_true.tau * (1 + mismatch*0.5);   % thinks force is more correlated
params_filter.q   = params_true.q   * 0.5;                  % underestimates process noise


%% Simulation Settings
dt = 0.001;         % Sample time [s]
dt_meas = 0.01;     % Measurement sample time [s] - sparse measurements
T = 5;              % Total simulation time [s]
t = 0:dt:T;         % Time vector
N = length(t);      % Number of samples

%% Measurement Settings
R = 1e-5;           % Position measurement variance [m²] - noisy measurements
sigma_meas = sqrt(R);

% Create measurement schedule (only every dt_meas)
meas_schedule = false(1, N);
meas_schedule(1:round(dt_meas/dt):end) = true;

%% Initial Conditions
x0_true = [0.002; 0.1; 0.5];     % True initial state with velocity and force
x0_est = [0; 0; 0];              % Filter starts from zero - poor initialization
P0 = diag([1e-3, 1, 1]);         % Large initial uncertainty

%% Generate Input Signal - Controlled excitation
u = zeros(N, 1);

% Sinusoidal base vibration
freqs = [1, 3, 7];  % Hz
amps = [0.002, 0.001, 0.0005];  % Small amplitudes
for i = 1:length(freqs)
    u = u + amps(i) * sin(2*pi*freqs(i)*t' + pi/4*i);
end

% Series of moderate pulses
pulse_times = [0.5, 1.5, 2.5, 3.5];
pulse_amps = [0.5, -0.7, 0.6, -0.4];  % Reduced amplitudes
for i = 1:length(pulse_times)
    pulse_idx = find(t >= pulse_times(i) & t < pulse_times(i) + 0.05);
    u(pulse_idx) = u(pulse_idx) + pulse_amps(i);
end

%% Simulate True System
fprintf('Simulating strongly nonlinear system...\n');
[x_true, z_meas_all, F_true] = simulate_true_system(params_true, x0_true, u, dt, T, sigma_meas);

% Apply measurement schedule (sparse measurements)
z_meas = z_meas_all;
z_meas(~meas_schedule) = NaN;  % No measurement available

%% Run Filters
fprintf('Running EKF with %.0f%% parameter mismatch...\n', mismatch*100);
tic;
[x_ekf, P_ekf] = ekf_filter(z_meas, u, params_filter, x0_est, P0, R, dt, meas_schedule);
t_ekf = toc;

fprintf('Running UKF with %.0f%% parameter mismatch...\n', mismatch*100);
tic;
[x_ukf, P_ukf] = ukf_filter(z_meas, u, params_filter, x0_est, P0, R, dt, meas_schedule);
t_ukf = toc;

fprintf('Running SLF with %.0f%% parameter mismatch...\n', mismatch*100);
tic;
[x_slf, P_slf] = slf_filter(z_meas, u, params_filter, x0_est, P0, R, dt, meas_schedule);
t_slf = toc;

fprintf('\nComputation times: EKF=%.3fs, UKF=%.3fs, SLF=%.3fs\n', t_ekf, t_ukf, t_slf);

%% Calculate RMSE (exclude initial transient)
start_idx = round(1.0/dt);  % Start after 1 second
rmse_ekf_pos = sqrt(mean((x_true(1,start_idx:end) - x_ekf(1,start_idx:end)).^2));
rmse_ekf_vel = sqrt(mean((x_true(2,start_idx:end) - x_ekf(2,start_idx:end)).^2));
rmse_ekf_force = sqrt(mean((x_true(3,start_idx:end) - x_ekf(3,start_idx:end)).^2));

rmse_ukf_pos = sqrt(mean((x_true(1,start_idx:end) - x_ukf(1,start_idx:end)).^2));
rmse_ukf_vel = sqrt(mean((x_true(2,start_idx:end) - x_ukf(2,start_idx:end)).^2));
rmse_ukf_force = sqrt(mean((x_true(3,start_idx:end) - x_ukf(3,start_idx:end)).^2));

rmse_slf_pos = sqrt(mean((x_true(1,start_idx:end) - x_slf(1,start_idx:end)).^2));
rmse_slf_vel = sqrt(mean((x_true(2,start_idx:end) - x_slf(2,start_idx:end)).^2));
rmse_slf_force = sqrt(mean((x_true(3,start_idx:end) - x_slf(3,start_idx:end)).^2));

%% Calculate peak errors
peak_ekf_pos = max(abs(x_true(1,:) - x_ekf(1,:)));
peak_ukf_pos = max(abs(x_true(1,:) - x_ukf(1,:)));
peak_slf_pos = max(abs(x_true(1,:) - x_slf(1,:)));

%% Display Results
fprintf('\n=== RMSE Results (after 1s transient) ===\n');
fprintf('Position RMSE [μm]:\n');
fprintf('  EKF: %.2f\n', rmse_ekf_pos*1e6);
fprintf('  UKF: %.2f\n', rmse_ukf_pos*1e6);
fprintf('  SLF: %.2f\n', rmse_slf_pos*1e6);
fprintf('Velocity RMSE [mm/s]:\n');
fprintf('  EKF: %.3f\n', rmse_ekf_vel*1000);
fprintf('  UKF: %.3f\n', rmse_ukf_vel*1000);
fprintf('  SLF: %.3f\n', rmse_slf_vel*1000);
fprintf('Force RMSE [N]:\n');
fprintf('  EKF: %.3f\n', rmse_ekf_force);
fprintf('  UKF: %.3f\n', rmse_ukf_force);
fprintf('  SLF: %.3f\n', rmse_slf_force);

fprintf('\n=== Peak Position Errors [μm] ===\n');
fprintf('  EKF: %.2f\n', peak_ekf_pos*1e6);
fprintf('  UKF: %.2f\n', peak_ukf_pos*1e6);
fprintf('  SLF: %.2f\n', peak_slf_pos*1e6);

%% Calculate improvements
improv_ukf_pos = (rmse_ekf_pos - rmse_ukf_pos) / rmse_ekf_pos * 100;
improv_slf_pos = (rmse_ekf_pos - rmse_slf_pos) / rmse_ekf_pos * 100;
improv_ukf_vel = (rmse_ekf_vel - rmse_ukf_vel) / rmse_ekf_vel * 100;
improv_slf_vel = (rmse_ekf_vel - rmse_slf_vel) / rmse_ekf_vel * 100;

fprintf('\n=== Improvements over EKF ===\n');
fprintf('  UKF: %.1f%% (pos), %.1f%% (vel)\n', improv_ukf_pos, improv_ukf_vel);
fprintf('  SLF: %.1f%% (pos), %.1f%% (vel)\n', improv_slf_pos, improv_slf_vel);

%% Downsample for plotting
downsample_factor = 10;
t_plot = t(1:downsample_factor:end);
x_true_plot = x_true(:, 1:downsample_factor:end);
x_ekf_plot = x_ekf(:, 1:downsample_factor:end);
x_ukf_plot = x_ukf(:, 1:downsample_factor:end);
x_slf_plot = x_slf(:, 1:downsample_factor:end);
z_plot = z_meas_all(1:downsample_factor:end);
meas_plot = meas_schedule(1:downsample_factor:end);

%% Plot Results with Error Analysis
plot_results(t_plot, x_true_plot, x_ekf_plot, x_ukf_plot, x_slf_plot, ...
    z_plot, meas_plot, rmse_ekf_pos, rmse_ekf_vel, rmse_ukf_pos, rmse_ukf_vel, ...
    rmse_slf_pos, rmse_slf_vel);

fprintf('\nSimulation complete. Figure saved as comparison_results.png\n');