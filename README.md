# QBHM-SPBO

MATLAB implementation of **QBHM-SPBO**, the quantum-inspired black hole model based
Student Psychology-Based Optimization algorithm, from:

> K. D. Bodha, V. Arun, A. Awasthi, B. Mahato, G. Fotis,
> *Quantum Black Hole-Driven Enhancement of Student Psychology-Based Optimization Algorithm
> for Benchmark and Power System Dispatch Solutions*, SN Computer Science (under review).

## How the algorithm works

SPBO explores well but exploits weakly. QBHM-SPBO adds a quantum-inspired black-hole operator:

* The **best student is the quantum black hole**.
* Each student's **normalised fitness** is `s_i = (f_i - f_best) / (f_worst - f_best + eps)`, which lies in [0, 1]
  and is unchanged by shifting or scaling the objective, so it works for positive, zero and negative values.
* **Event horizon** (grows linearly): `R_t = R_min + (R_max - R_min) * t/T`.
* **Absorption** (exploitation): a student with `s_i <= R_t` is resampled in a quantum delta-potential well
  between itself and the black hole,
  `x = p + sign * beta_t * |x_mean - x_i| * ln(1/u)`,  with `p = phi*x_best + (1-phi)*x_i`,
  and contraction-expansion coefficient `beta_t = beta_max - (beta_max - beta_min) * t/T`.
* Students outside the horizon follow the **original SPBO rules** (exploration).
* **Hawking re-emission**: a student that has not improved for `L` iterations is re-emitted around the black hole.
* Greedy selection; one function evaluation per student per iteration (same budget as SPBO).

| Parameter | Value used in the paper |
|---|---|
| `Rmin`, `Rmax` | 0.10, 0.50 |
| `bmax`, `bmin` (beta) | 1.00, 0.50 |
| `L` (stagnation limit) | max(20, round(T/50)) |

## Usage

```matlab
[x_best, f_best, curve, info] = QBHM_SPBO(pop, T, UB, LB, D, fobj, maxFE, opts);
```

| Input | Meaning |
|---|---|
| `pop` | number of students (population size) |
| `T` | maximum number of iterations |
| `UB`, `LB` | upper / lower bounds (scalar or 1-by-D vector) |
| `D` | dimension |
| `fobj` | handle to the objective, called with a 1-by-D row vector |
| `maxFE` | optional evaluation budget (default `pop*T + pop`) |
| `opts` | optional struct: `absorb`, `hawking` (true/false), `Rmin`, `Rmax`, `bmax`, `bmin`, `L` |

`info` contains the number of evaluations, the run time, the population diversity per iteration and the
number of absorbed / re-emitted students per iteration.

Run `demo.m` for a complete example. No toolbox is required.


## Reproducing the experiments of the paper

| Folder | Content |
|---|---|
| `QBHM_SPBO.m`, `demo.m` | the algorithm and a minimal example |
| `benchmarks_CECfn/` | the 15 CEC-inspired test functions CECfn1–CECfn15 (D = 30) |
| `dispatch/cl_dispatch.m` | all dispatch test systems (3-, 5- and 26-unit ED; 6- and 10-unit CEED) with their data, slack-unit handling and penalty |
| `experiments/` | `cl_job.m` (one run under an exact budget of function evaluations, counted by a wrapper around the objective), `run_experiments.m` (all experiments), SPBO (`cl_SPBO.m`), the single-jump variant used in the ablation study (`cl_origQBHM.m`) and the eight competitors (`bl_*.m`, `EO.m`) with the parameter settings of Table 3 |
| `results/` | raw results of every run: `results_CECfn.csv`, `results_CEC22.csv`, `results_DISP.csv` (best value of each run; for dispatch also the full dispatch, fuel cost, emission, losses, power-balance mismatch and limit violation) |

Algorithm labels in the result files: `QBHM` = QBHM-SPBO, `QBHM_orig` = single-jump variant QBHM-SPBO-SJ, `QBHM_noAbs` / `QBHM_noHawk` = without absorption / without re-emission. Run *r* uses the random seed 2026 + *r* for every algorithm.

The competitor implementations follow the original publications cited in the paper; `EO.m` is adapted from the public code of Faramarzi et al. (2020) with the memory-saving step removed, as described in the paper. The official CEC 2022 code (`cec22_test_func` and its `input_data`) is not redistributed here; it is available from https://github.com/P-N-Suganthan/2022-SO-BO.

## Citation

If you use this code, please cite the paper above.

## License

MIT (see `LICENSE`).
