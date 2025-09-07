function [x_est, P_est] = slf_filter(z, u, params, x0, P0, R, dt, meas_schedule)
%%   Statistically Linearized Filter (sparse measurements, OU force)

%
% Inputs
%   z      [1xN]   position measurements (NaN for missing)
%   u      [Nx1]   control force sequence (same dt)
%   params struct  with fields: m, b, k1, k3, tau, q
%   x0     [3x1]   initial state estimate
%   P0     [3x3]   initial covariance
%   R      scalar  measurement variance
%   dt     scalar  sample time
%   meas_schedule [1xN] logical, true where a measurement is available
%
% Outputs
%   x_est  [3xN]   state estimates over time
%   P_est  [3x3xN] covariances over time


% Unpack parameters
if nargin < 9 || isempty(opts), opts = struct; end
if ~isfield(opts,'use_rk4'),   opts.use_rk4 = true; end
if ~isfield(opts,'Qxx'),       opts.Qxx = 0; end
if ~isfield(opts,'Qvv'),       opts.Qvv = 0; end
if ~isfield(opts,'gate_chi2'), opts.gate_chi2 = 25; end

% Unpack parameters
m   = params.m;
b   = params.b;
k1  = params.k1;
k3  = params.k3;
tau = params.tau;
q   = params.q;   % OU noise intensity

% Sizes
N = numel(z);
x_est = zeros(3, N);
P_est = zeros(3, 3, N);

% Initialize
x_k = x0(:);
P_k = P0;

% OU discretization
phiF = exp(-dt/tau);
Q33  = q * (tau/2) * (1 - exp(-2*dt/tau));

% Measurement model
H = [1 0 0];

% Covariance floor
eps_floor = 1e-12;

for k = 1:N
    %====================
    % Prediction (SLF)
    %====================
    x  = x_k(1);  v = x_k(2);  F = x_k(3);
    uk = u(k);

    % Statistical linearization pieces
    mu_x    = x;
    sig2_x  = max(P_k(1,1), 0);
    Ex3     = mu_x^3 + 3*mu_x*sig2_x;
    dEx3_dm = 3*mu_x^2 + 3*sig2_x;

    % Define accel function using Ex3 frozen over the step (common SLF choice)
    % If opts.use_rk4, we update Ex3 on intermediate stages with a local sigma proxy.
    acc_fun = @(x_local, v_local, F_local, Ex3_local) ...
        (1/m) * (uk - b*v_local - k1*x_local - k3*Ex3_local + F_local);

    % State prediction
    x_pred = zeros(3,1);

    if opts.use_rk4
        % Local proxy for sigma at stages: use same P_k(1,1)
        f1_x = v;
        f1_v = acc_fun(x, v, F, Ex3);

        x2   = x + 0.5*dt*f1_x;
        v2   = v + 0.5*dt*f1_v;
        Ex3_2 = x2^3 + 3*x2*sig2_x;

        f2_x = v2;
        f2_v = acc_fun(x2, v2, F, Ex3_2);

        x3   = x + 0.5*dt*f2_x;
        v3   = v + 0.5*dt*f2_v;
        Ex3_3 = x3^3 + 3*x3*sig2_x;

        f3_x = v3;
        f3_v = acc_fun(x3, v3, F, Ex3_3);

        x4   = x + dt*f3_x;
        v4   = v + dt*f3_v;
        Ex3_4 = x4^3 + 3*x4*sig2_x;

        f4_x = v4;
        f4_v = acc_fun(x4, v4, F, Ex3_4);

        x_pred(1) = x + (dt/6)*(f1_x + 2*f2_x + 2*f3_x + f4_x);
        x_pred(2) = v + (dt/6)*(f1_v + 2*f2_v + 2*f3_v + f4_v);
    else
        acc = acc_fun(x, v, F, Ex3);
        x_pred(1) = x + dt * v;
        x_pred(2) = v + dt * acc;
    end

    % Exact OU mean for force
    x_pred(3) = phiF * F;

    % Linearized A matrix about (x,v,F) with Ex3 derivative
    A = [ 1,             dt,                 0;
         -dt*(k1 + k3*dEx3_dm)/m, 1 - dt*(b/m),  dt*(1/m);
          0,              0,                 phiF ];

    % Process noise
    Q = diag([opts.Qxx, opts.Qvv, Q33]);

    % Covariance prediction
    P_pred = A * P_k * A.' + Q;
    P_pred = 0.5*(P_pred + P_pred.');

    %====================
    % Update (if meas)
    %====================
    if meas_schedule(k) && ~isnan(z(k))
        y = z(k) - H*x_pred;
        S = H*P_pred*H.' + R;
        if (y^2)/S <= opts.gate_chi2
            K = (P_pred*H.') / S;
            x_upd = x_pred + K*y;
            P_upd = (eye(3)-K*H)*P_pred;
        else
            x_upd = x_pred;
            P_upd = P_pred;
        end
    else
        x_upd = x_pred;
        P_upd = P_pred;
    end

    % PSD floor
    [V,D] = eig(0.5*(P_upd+P_upd.'));
    D = diag(D);
    D = max(D, eps_floor);
    P_upd = V*diag(D)*V.';

    % Store
    x_est(:,k)   = x_upd;
    P_est(:,:,k) = P_upd;
    x_k = x_upd;
    P_k = P_upd;
end
end
