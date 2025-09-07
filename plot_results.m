function plot_results(t, x_true, x_ekf, x_ukf, x_slf, z_meas, meas_schedule, ...
    rmse_ekf_pos, rmse_ekf_vel, rmse_ukf_pos, rmse_ukf_vel, ...
    rmse_slf_pos, rmse_slf_vel)
% Shows both state estimates and estimation errors to highlight differences

% Inputs:
%   t                 [1xN] time
%   x_true, x_ekf, x_ukf, x_slf  [3xN] state trajectories
%   z_meas            [1xN] position measurements (NaN where none)
%   meas_schedule     [1xN] logical indices of available measurements
%   rmse_*_pos/vel    scalars for summary panel
%
% Force RMSE is computed directly from trajectories.


% ---- Colors / styles ----
col_true = [0 0 0];        % black
col_ekf  = [0 0.4 0.8];    % blue
col_ukf  = [0.8 0.2 0.2];  % red
col_slf  = [0.2 0.7 0.2];  % green
col_meas = [0.6 0.6 0.6];  % light gray

ls_ekf = '-';  ls_ukf = '--';  ls_slf = ':';

% ---- Prep ----
if nargin < 13
    error('plot_results_with_errors:NotEnoughInputs', 'Function signature changed? Provide 13 inputs.');
end

meas_idx = find(meas_schedule);
if isempty(meas_idx), meas_idx = []; end

% Errors
pos_err_ekf = (x_true(1,:) - x_ekf(1,:)) * 1e6;    % [um]
pos_err_ukf = (x_true(1,:) - x_ukf(1,:)) * 1e6;
pos_err_slf = (x_true(1,:) - x_slf(1,:)) * 1e6;

vel_err_ekf = (x_true(2,:) - x_ekf(2,:)) * 1000;   % [mm/s]
vel_err_ukf = (x_true(2,:) - x_ukf(2,:)) * 1000;
vel_err_slf = (x_true(2,:) - x_slf(2,:)) * 1000;

F_err_ekf   = (x_true(3,:) - x_ekf(3,:));          % [N]
F_err_ukf   = (x_true(3,:) - x_ukf(3,:));
F_err_slf   = (x_true(3,:) - x_slf(3,:));

% Compute force RMSE for summary panel
rmse_ekf_F = rms(F_err_ekf);
rmse_ukf_F = rms(F_err_ukf);
rmse_slf_F = rms(F_err_slf);

% ---- Figure ----
figure('Position', [50, 50, 1600, 1000], 'Color', 'w');

%% TOP ROW: State Estimates
% Position
subplot(3,3,1); hold on; grid on; box on;
if ~isempty(meas_idx)
    h_meas = plot(t(meas_idx), z_meas(meas_idx)*1e6, '.', 'Color', col_meas, 'MarkerSize', 4);
    set(h_meas, 'HandleVisibility','off'); % hide from legend
end
h1 = plot(t, x_true(1,:)*1e6, '-',  'Color', col_true, 'LineWidth', 2,   'DisplayName','Truth');
h2 = plot(t, x_ekf(1,:)*1e6,  ls_ekf,'Color', col_ekf,  'LineWidth', 1.5, 'DisplayName','EKF');
h3 = plot(t, x_ukf(1,:)*1e6,  ls_ukf,'Color', col_ukf,  'LineWidth', 1.5, 'DisplayName','UKF');
h4 = plot(t, x_slf(1,:)*1e6,  ls_slf,'Color', col_slf,  'LineWidth', 2.0, 'DisplayName','SLF');
ylabel('Position [\mum]'); title('Position Estimates');
xlim([t(1) t(end)]); legend([h1 h2 h3 h4], 'Location','northeast');

% Velocity
subplot(3,3,2); hold on; grid on; box on;
plot(t, x_true(2,:)*1000, '-',   'Color', col_true, 'LineWidth', 2,   'DisplayName','Truth');
plot(t, x_ekf(2,:)*1000,  ls_ekf,'Color', col_ekf,  'LineWidth', 1.5, 'DisplayName','EKF');
plot(t, x_ukf(2,:)*1000,  ls_ukf,'Color', col_ukf,  'LineWidth', 1.5, 'DisplayName','UKF');
plot(t, x_slf(2,:)*1000,  ls_slf,'Color', col_slf,  'LineWidth', 2.0, 'DisplayName','SLF');
ylabel('Velocity [mm/s]'); title('Velocity Estimates');
xlim([t(1) t(end)]); legend('Location','northeast');

% Force
subplot(3,3,3); hold on; grid on; box on;
plot(t, x_true(3,:), '-',    'Color', col_true, 'LineWidth', 2,   'DisplayName','Truth');
plot(t, x_ekf(3,:),  ls_ekf, 'Color', col_ekf,  'LineWidth', 1.5, 'DisplayName','EKF');
plot(t, x_ukf(3,:),  ls_ukf, 'Color', col_ukf,  'LineWidth', 1.5, 'DisplayName','UKF');
plot(t, x_slf(3,:),  ls_slf, 'Color', col_slf,  'LineWidth', 2.0, 'DisplayName','SLF');
ylabel('Force [N]'); title('Stochastic Force Estimates');
xlim([t(1) t(end)]); legend('Location','northeast');

%% MIDDLE ROW: Estimation Errors
% Position error
subplot(3,3,4); hold on; grid on; box on;
plot(t, pos_err_ekf, ls_ekf, 'Color', col_ekf, 'LineWidth', 1.5, 'DisplayName','EKF');
plot(t, pos_err_ukf, ls_ukf, 'Color', col_ukf, 'LineWidth', 1.5, 'DisplayName','UKF');
plot(t, pos_err_slf, ls_slf, 'Color', col_slf, 'LineWidth', 2.0, 'DisplayName','SLF');
h0 = yline(0,'k-','LineWidth',0.5); set(h0,'HandleVisibility','off');
if ~isempty(meas_idx)
    use_idx = meas_idx(1:max(1,ceil(numel(meas_idx)/20)):end);
    for ii = 1:numel(use_idx)
        hx = xline(t(use_idx(ii)), 'k:', 'LineWidth', 0.5); set(hx,'HandleVisibility','off');
    end
end
ylabel('Position Error [\mum]'); xlabel('Time [s]'); title('Position Error');
xlim([t(1) t(end)]); legend('Location','northeast');

% Velocity error
subplot(3,3,5); hold on; grid on; box on;
plot(t, vel_err_ekf, ls_ekf, 'Color', col_ekf, 'LineWidth', 1.5, 'DisplayName','EKF');
plot(t, vel_err_ukf, ls_ukf, 'Color', col_ukf, 'LineWidth', 1.5, 'DisplayName','UKF');
plot(t, vel_err_slf, ls_slf, 'Color', col_slf, 'LineWidth', 2.0, 'DisplayName','SLF');
h0 = yline(0,'k-','LineWidth',0.5); set(h0,'HandleVisibility','off');
ylabel('Velocity Error [mm/s]'); xlabel('Time [s]'); title('Velocity Error');
xlim([t(1) t(end)]); legend('Location','northeast');

% Force error
subplot(3,3,6); hold on; grid on; box on;
plot(t, F_err_ekf, ls_ekf, 'Color', col_ekf, 'LineWidth', 1.5, 'DisplayName','EKF');
plot(t, F_err_ukf, ls_ukf, 'Color', col_ukf, 'LineWidth', 1.5, 'DisplayName','UKF');
plot(t, F_err_slf, ls_slf, 'Color', col_slf, 'LineWidth', 2.0, 'DisplayName','SLF');
h0 = yline(0,'k-','LineWidth',0.5); set(h0,'HandleVisibility','off');
ylabel('Force Error [N]'); xlabel('Time [s]'); title('Force Error');
xlim([t(1) t(end)]); legend('Location','northeast');

%% BOTTOM ROW: Error Statistics
% Running RMSE (position) with simple moving window
subplot(3,3,7); hold on; grid on; box on;
window_size = max(50, round(numel(t)/200));
rms_ekf = simple_mov_rms((x_true(1,:) - x_ekf(1,:))*1e6, window_size);
rms_ukf = simple_mov_rms((x_true(1,:) - x_ukf(1,:))*1e6, window_size);
rms_slf = simple_mov_rms((x_true(1,:) - x_slf(1,:))*1e6, window_size);
plot(t, rms_ekf, ls_ekf, 'Color', col_ekf, 'LineWidth', 1.5, 'DisplayName','EKF');
plot(t, rms_ukf, ls_ukf, 'Color', col_ukf, 'LineWidth', 1.5, 'DisplayName','UKF');
plot(t, rms_slf, ls_slf, 'Color', col_slf, 'LineWidth', 2.0, 'DisplayName','SLF');
ylabel('Running RMSE [\mum]'); xlabel('Time [s]');
title(sprintf('Position RMSE (window=%d)', window_size));
xlim([t(1) t(end)]); legend('Location','northeast');

% Error histogram (position) — compact, non-overlapping
subplot(3,3,8); hold on; grid on; box on;
% Common edges based on pooled std
sig_pool = max([std(pos_err_ekf), std(pos_err_ukf), std(pos_err_slf), eps]);
edges = linspace(-3, 3, 30) * sig_pool;
[N1, ~] = histcounts(pos_err_ekf, edges);
[N2, ~] = histcounts(pos_err_ukf, edges);
[N3, ~] = histcounts(pos_err_slf, edges);
centers = (edges(1:end-1) + edges(2:end))/2;
w = centers(2) - centers(1);
bar(centers - w*0.33, N1, 0.3, 'FaceColor', col_ekf, 'EdgeColor','none', 'FaceAlpha',0.6, 'DisplayName','EKF');
bar(centers,           N2, 0.3, 'FaceColor', col_ukf, 'EdgeColor','none', 'FaceAlpha',0.6, 'DisplayName','UKF');
bar(centers + w*0.33, N3, 0.3, 'FaceColor', col_slf, 'EdgeColor','none', 'FaceAlpha',0.6, 'DisplayName','SLF');
xlabel('Position Error [\mum]'); ylabel('Count'); title('Error Distribution');
legend('Location','northeast');

% RMSE comparison bar chart (pos/vel/force) with percentage improvements
subplot(3,3,9); hold on; grid on; box on;
rmse_pos = [rmse_ekf_pos, rmse_ukf_pos, rmse_slf_pos]*1e6;   % [um]
rmse_vel = [rmse_ekf_vel, rmse_ukf_vel, rmse_slf_vel]*1000;  % [mm/s]
rmse_for = [rmse_ekf_F, rmse_ukf_F, rmse_slf_F];             % [N]
vals = [rmse_pos; rmse_vel; rmse_for];  % 3x3
bh = bar(vals.', 'LineWidth', 0.5);     % groups: EKF, UKF, SLF
bh(1).FaceColor = [0.85 0.85 0.85];     % Position
bh(2).FaceColor = [0.65 0.65 0.65];     % Velocity
bh(3).FaceColor = [0.45 0.45 0.45];     % Force
xticklabels({'EKF','UKF','SLF'});
ylabel('RMSE (units vary)');
title('RMSE Summary: Pos [\mum], Vel [mm/s], Force [N]');
legend({'Pos','Vel','Force'}, 'Location','northoutside','Orientation','horizontal','Box','off');

% annotate improvements vs EKF over each metric
for alg = 2:3  % UKF, SLF
    for metric = 1:3
        ekf_val = vals(metric,1);
        val     = vals(metric,alg);
        if ekf_val > 0
            imp = (ekf_val - val)/ekf_val * 100;
        else
            imp = 0;
        end
        % approximate label position above the bar
        x = alg + (metric-2)*0.25;
        y = val;
        txt = sprintf('%+.1f%%', imp);
        text(x, y*1.02, txt, 'HorizontalAlignment','center', 'FontSize',8);
    end
end

% Overall title
sgtitle('Nonlinear Vibration Isolation — Filter Comparison', 'FontSize', 14, 'FontWeight', 'bold');

% ---- helper: simple moving RMS (centered window) ----
function y = simple_mov_rms(x, win)
    x = x(:).';
    n = numel(x);
    y = zeros(1,n);
    half = floor(win/2);
    for i = 1:n
        i1 = max(1, i-half);
        i2 = min(n, i+half);
        seg = x(i1:i2);
        y(i) = sqrt(mean(seg.^2));
    end
end

end