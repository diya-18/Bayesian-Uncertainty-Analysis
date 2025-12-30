# Bayesian Uncertainty Analysis for Scientific Measurements

This project demonstrates how uncertainty is handled in real scientific and engineering measurements using statistical and Bayesian methods.

## What this project is about
In real-world experiments, sensor measurements are noisy and the true value of a physical quantity is unknown. This project shows how to:
- Simulate noisy measurements
- Estimate true values
- Quantify uncertainty
- Compare frequentist and Bayesian approaches

## Project Structure

### 1. Measurement Simulation
Simulated repeated measurements of a physical quantity with sensor noise to represent real instrumentation data.

### 2. Statistical Estimation
Computed mean, standard deviation, standard error, and confidence intervals to estimate the measured quantity using classical statistics.

### 3. Monte Carlo Uncertainty Propagation
Used Monte Carlo simulation to propagate uncertainty through a derived quantity (power = voltage × current).

### 4. Bayesian Estimation (Known Noise)
Applied Bayesian inference to combine prior knowledge with measured data and compute a posterior probability distribution.

### 5. Frequentist vs Bayesian Comparison
Compared point estimates from frequentist statistics with full probability distributions from Bayesian inference.

### 6. Bayesian Estimation with Unknown Noise
Extended the Bayesian model by treating measurement noise as an unknown parameter and estimating it directly from data.

### 7. Model Comparison
Compared Bayesian posteriors assuming known and unknown noise levels.

## Why this project matters
- Demonstrates uncertainty-aware thinking essential in physics and engineering
- Uses Monte Carlo methods widely applied in scientific research
- Shows Bayesian reasoning used in detector calibration and experimental analysis
- Reflects real measurement challenges in large-scale research facilities

## Tools Used
- R
- Statistical modeling
- Monte Carlo simulation
- Bayesian inference

This project focuses on scientific rigor rather than black-box machine learning models.
