# Drain-Based Irrigation Scheduling with Offline Reinforcement Learning

This repository contains the data pipeline, training and evaluation code for drain-based irrigation
scheduling in protected strawberry cultivation using offline reinforcement learning. A policy decides
when the next irrigation should take place, based on the current greenhouse state and the measured
drainage. The work is fueled in the paper *Drain-Based Irrigation Scheduling with Offline Reinforcement
Learning in Protected Strawberry Cultivation* (Heilbronn University of Applied Sciences).

## Overview

The goal is to keep the daily drainage close to a target value (about 20 percent of the supplied water,
roughly 1 liter) without continuous human intervention. Instead of fixed irrigation timers, an offline
RL policy is trained on episodes that describe complete irrigation cycles or cultivation days. The
training data is generated from a FAO-56 evapotranspiration-drainage model that is calibrated on weather
station and sensor data, so policies can be trained without online interaction with the real system.

The offline RL algorithms are implemented with [d3rlpy](https://github.com/takuseno/d3rlpy). The project
benchmarks a wide range of algorithms and uses BCQ, PLAS and IQL for the final models.

## Problem formulation

- Observation (7 values): `doy` (day of year), `abstime` (absolute time of day), `time` (time since the
  last irrigation), `temp` (temperature), `humi` (humidity), `light` (radiation), `drain` (drainage).
- Action: `timeToIrri`, a continuous scalar giving the predicted minutes until the next irrigation.
- Reward: derived from the drainage target, highest near the target band and lower outside of it.
- Sampling and trigger: states are sampled every 10 minutes; an irrigation is triggered when the policy
  predicts a value of 10 minutes or less.

## Episode strategies

Two ways of structuring the episodes are compared:

- Irrigation-interval-based (referred to as *single episode*): each episode covers one irrigation cycle
  and the cumulative state variables are reset at every irrigation. This setup reaches the highest
  scheduling accuracy.
- Full-day-based (referred to as *day* or *multi episode*): each episode covers a full cultivation day
  with several irrigations and no reset between them.

## Repository structure

```
TEST_WATERS/
  data/            Datasets and data generation
    CSV AgrarMeteo/  Raw weather station data
    episFormula/     Generated episode datasets (CSV per episode)
    script/          Data preprocessing and episode generation
  docs/            Paper, presentations and sources
  script/          Offline RL training and evaluation notebooks
  requirements.txt Python dependencies
  README.md
```

Training and evaluation notebooks write their outputs to a local `logs/` directory (models and metrics),
which is created on first run.

## Data pipeline (`data/script`)

The model is calibrated first and then used to generate the episode datasets. The simulated episodes
(from weather station data) and the measured episodes (from InfluxDB sensor data) are produced in two
parallel paths:

| File | Purpose |
| --- | --- |
| `getFormula.ipynb` | Calibrate the FAO-56 evapotranspiration-drainage model against the measured InfluxDB data (radiation from lux, `Kr`, clear-sky radiation, `Kc`) and validate it against the measured drainage. |
| `formulaEvo.py` | The calibrated FAO-56 model used for episode generation (clear-sky radiation, crop coefficient spline, evapotranspiration). |
| `creatEpis.ipynb` | Generate the simulated episodes from the AgrarMeteo weather station CSVs using `formulaEvo`. |
| `preprocessInfluxDB.ipynb` | Query the raw InfluxDB sensor data, downsample to 10-minute intervals and build the measured episodes (the measured counterpart to `creatEpis`). |
| `plot_epi_data.ipynb` | Load and visualize a single episode CSV. |

Each episode is stored as a CSV file (named by weather station, month, day and episode index). The
columns are `doy`, `abstime`, `time`, `temp`, `humi`, `light`, `drain`, `timeToIrri`, `reward`,
`terminal`, `idealIrri`, `evoCumsum`; the first seven columns form the observation.

Each notebook starts with a short description of its purpose and steps.

## Offline RL workflow (`script`)

| Notebook | Purpose |
| --- | --- |
| `test d3rlpy.ipynb` | Verify the GPU stack and the d3rlpy installation. |
| `sendItEveryAlgo.ipynb` | Benchmark all d3rlpy continuous-control algorithms to select the best ones. |
| `train Single Episode.ipynb` | Train BCQ, PLAS and IQL on the irrigation-interval-based dataset. |
| `train Multi Episode.ipynb` | Train policies on the full-day-based dataset. |
| `testModelsMultiEpisode.ipynb` | Evaluate the interval-based models on simulated and measured data. |
| `testModelsDayEpisode.ipynb` | Roll out the full-day models over complete days and plot the schedule. |
| `testModelsDayEpisodeEvaluation.ipynb` | Quantitative evaluation of the full-day models (offset metrics and tables). |
| `analyseModels.ipynb` | Visualize how predictions evolve over training and compare the `diff_eval` metric. |

A typical run goes through verification, algorithm benchmarking, training and finally evaluation. Each
notebook starts with a short description of its purpose and steps.

## Setup

The notebooks require Python with a CUDA-capable GPU (the configurations use `device='cuda:0'`). Install
the dependencies with:

```
pip install -r requirements.txt
```

The main dependencies are `d3rlpy` (offline RL), `torch`, `numpy`, `pandas`, `matplotlib` and `scipy`.