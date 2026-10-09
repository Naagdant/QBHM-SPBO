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

## Citation

If you use this code, please cite the paper above.

## License

MIT (see `LICENSE`).
