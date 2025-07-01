# Predicting Hospital Encounters with Infectious Diseases (COVID-19, Influenza, RSV) in Regions and Counties of South Carolina Utilizing State Predictions

## Overview

This repository provides a flexible and scalable framework to distribute infectious disease–related hospitalizations (COVID-19, influenza, RSV) from statewide projections to smaller geographic units: counties and regions in South Carolina.

The approach integrates:
- Electronic health records (EHRs) from MUSC and Prisma Health,
- Hospital admission data from the SC Revenue and Fiscal Affairs (RFA) office,
- CDC reported statewidhospitalization (2024–2025),
- Machine learning (XGBoost) models optimized using Bayesian hyperparameter tuning.

---

## Repository Structure


├── scripts/               # R scripts for preprocessing, modeling, visualization
├── figures/               # Generated graphs 
├── output/                # Prediction results and performance metrics
├── docs/                  # Reports and supplementary documentation
├── README.md              # Project overview and usage guide
├── LICENSE                # License file
├── .gitignore             # Files to ignore in git

---

## Methods Overview

- Data Integration: EHR + RFA data (2021–2023)
- Smoothing: Peak-preserving moving average
- Modeling: XGBoost regression
- Optimization: Bayesian tuning of hyperparameters
- Evaluation: Percent Agreement (PA) between distributed and observed counts
  



