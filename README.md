
# Vibration Isolation State Estimation

Real-time state estimation for a sensor-pod isolation system using EKF, UKF, and SLF on a nonlinear mass–spring–damper with OU forcing.

## Quick Start

```matlab
>> main
```

Runs the simulation and outputs comparison figures + RMSE.

## Application

Inspection-robot sensor pod isolating base vibrations/external disturbances.

## Repo

* `main.m` – entry point
* `simulate_true_system.m` – truth MSD with OU force
* `ekf_filter.m`, `ukf_filter.m`, `slf_filter.m` – filters
* `plot_results.m`  – plots & metrics

## Features

* Position-only measurements (sparse supported)
* Parameter mismatch handling (≈±20–30%)
* Estimates position, velocity, and stochastic force

## Results

![output](https://github.com/user-attachments/assets/0a99c6fa-a127-4f9d-972f-68bc9bc6e79f)


RMSE after a 1 s transient, three random seeds (MATLAB R2025b):

| | Position [µm] | Velocity [mm/s] | Force [N] |
|---|---|---|---|
| EKF | 1939–2083 | 43.8–56.7 | 1.27–1.49 |
| UKF | within 0.2% of EKF | within 0.2% | within 0.3% |
| SLF | within 0.4% of EKF | within 0.4% | within 0.3% |

For this system the three filters give essentially the same answer. The nonlinearity is mild at these displacement amplitudes, so first-order linearization, the unscented transform and statistical linearization all see nearly the same model. A case with larger excursions into the nonlinear spring region would be needed to separate them.
