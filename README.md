
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


UKF often yields the lowest peaks and \~10–15% lower velocity RMSE vs EKF under mismatch; SLF sits between, with robust tuning narrowing the gap.
