//! Shared forced pseudo-arclength corrector in physical or normalized units.

use super::scaling::{AdmittedScales, length};
use super::{
    AftWorkspace, ContinuableProblem, ContinuationBudget, ContinuationPoint, HbAnchor, HbOrbit,
    HbScales, MAX_AFT_REFINEMENTS, OrbitError, Packing, ScaledHbOrbit, check_sampling, det,
    hb_residual, lu, norm, pack_orbit, solve_hb, solve_hb_scaled_seeded,
};

/// Physical conditioning of a forced branch, including its swept parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct ContinuationScales {
    /// Original state/equation units and angular-frequency conditioning.
    pub hb: HbScales,
    /// Positive finite magnitude in the continuation parameter's own units.
    /// Pressure, tension and angular frequency need distinct physical scales.
    pub parameter: f64,
}

/// One physical branch point with independently checked normalized balance.
#[derive(Debug, Clone)]
pub struct ScaledContinuationPoint {
    /// Original parameter and Fourier coefficients. Here `step` is normalized
    /// arclength and orbit residuals contain physical HB rows only.
    pub point: ContinuationPoint,
    /// Verified Parseval RMS norm of HB defects divided by equation scales.
    pub scaled_residual: f64,
    /// Same norm of the retained balance's independent AFT-grid change.
    pub scaled_sampling_defect: f64,
    /// Absolute dimensionless arclength-equation defect; zero at the seed.
    pub arclength_residual: f64,
}

fn invalid(what: &'static str) -> OrbitError {
    OrbitError::BadParameter { what }
}

/// Continue a selected full forced waveform through parameter folds.
///
/// The canonical physical seed contains every DC and complex harmonic
/// coefficient. The initial point is solved again with the original source;
/// previous convergence metadata is never trusted. Parameter changes remain
/// owned by `set_parameter`, and `omega_of` supplies the physical clock.
///
/// The same corrector as [`super::continue_branch`] uses Parseval-normalized
/// state coordinates and a separately scaled parameter for its predictor,
/// secant and arclength row. Newton columns, equation rows and the independent
/// finer AFT check use the same units. Both augmented normalized balance and
/// sampling defect must be below `budget.hb.tolerance`, irrespective of DC
/// preload. `initial_step`/`min_step` are dimensionless arclength magnitudes.
/// Returned coefficients, parameter, clock and physical work remain unchanged.
///
/// This fixed-order numerical branch is not a branch-identity, stability,
/// fold-location or infinite-harmonic certificate. Newton can find another
/// branch; source parameter state on refusal remains caller-owned.
///
/// # Errors
/// Malformed seeds, scales and budgets refuse before source callbacks.
/// Numerical and sampling failures retain the original [`OrbitError`].
#[allow(clippy::too_many_arguments)]
pub fn continue_branch_scaled_seeded<P: ContinuableProblem>(
    problem: &mut P,
    omega_of: impl Fn(f64) -> f64,
    lambda0: f64,
    lambda1: f64,
    seed: Vec<f64>,
    budget: &ContinuationBudget,
    scales: &ContinuationScales,
) -> Result<Vec<ScaledContinuationPoint>, OrbitError> {
    let positive = |v: f64| v.is_finite() && v > 0.0;
    if !lambda0.is_finite()
        || !lambda1.is_finite()
        || lambda0 == lambda1
        || budget.max_steps == 0
        || !positive(budget.initial_step)
        || !positive(budget.min_step)
        || budget.min_step > budget.initial_step
        || !positive(scales.parameter)
    {
        return Err(invalid(
            "finite distinct branch endpoints and positive continuation scales/budgets required",
        ));
    }
    let anchor = HbAnchor::Forced {
        omega: omega_of(lambda0),
    };
    let admitted = AdmittedScales::new(problem.dim(), anchor, &seed, &budget.hb, &scales.hb)?;
    admitted
        .columns
        .len()
        .checked_add(1)
        .and_then(|n| n.checked_mul(n))
        .ok_or_else(|| invalid("bordered continuation dimensions overflow"))?;
    problem.set_parameter(lambda0);
    let first = solve_hb_scaled_seeded(problem, anchor, seed, &budget.hb, &scales.hb)?;
    continue_impl(
        problem,
        omega_of,
        lambda0,
        lambda1,
        first,
        budget,
        Metric {
            hb: Some(admitted),
            parameter: scales.parameter,
        },
    )
}

pub(super) fn legacy<P: ContinuableProblem>(
    problem: &mut P,
    omega_of: impl Fn(f64) -> f64,
    lambda0: f64,
    lambda1: f64,
    guess_amplitude: f64,
    budget: &ContinuationBudget,
) -> Result<Vec<ContinuationPoint>, OrbitError> {
    if !(budget.initial_step > 0.0 && budget.initial_step.is_finite()) {
        return Err(invalid("initial step must be positive and finite"));
    }
    problem.set_parameter(lambda0);
    let orbit = solve_hb(
        problem,
        HbAnchor::Forced {
            omega: omega_of(lambda0),
        },
        guess_amplitude,
        &budget.hb,
    )?;
    let first = ScaledHbOrbit {
        scaled_residual: orbit.residual,
        scaled_sampling_defect: orbit.sampling_defect,
        orbit,
    };
    Ok(continue_impl(
        problem,
        omega_of,
        lambda0,
        lambda1,
        first,
        budget,
        Metric {
            hb: None,
            parameter: 1.0,
        },
    )?
    .into_iter()
    .map(|p| p.point)
    .collect())
}

struct Metric {
    hb: Option<AdmittedScales>,
    parameter: f64,
}

impl Metric {
    fn column(&self, col: usize, nc: usize) -> f64 {
        if col == nc {
            self.parameter
        } else {
            self.hb.as_ref().map_or(1.0, |s| s.columns[col])
        }
    }
    fn row(&self, row: usize) -> f64 {
        self.hb.as_ref().map_or(1.0, |s| s.rows[row])
    }
    fn balance(&self, residual: &[f64]) -> f64 {
        self.hb
            .as_ref()
            .map_or_else(|| norm(residual), |s| s.row_norm(residual))
    }
    fn augmented(&self, residual: &[f64], arc: f64) -> f64 {
        if self.hb.is_some() {
            det::hypot(self.balance(residual), arc)
        } else {
            det::sqrt(residual.iter().map(|v| v * v).sum::<f64>() + arc * arc)
        }
    }
    fn raw(&self, residual: &[f64]) -> f64 {
        if self.hb.is_some() {
            length(residual.iter().copied())
        } else {
            norm(residual)
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn continue_impl<P: ContinuableProblem>(
    problem: &mut P,
    omega_of: impl Fn(f64) -> f64,
    lambda0: f64,
    lambda1: f64,
    first: ScaledHbOrbit,
    budget: &ContinuationBudget,
    metric: Metric,
) -> Result<Vec<ScaledContinuationPoint>, OrbitError> {
    let pack = Packing {
        d: problem.dim(),
        n_harm: budget.hb.harmonics,
        omega_slot: None,
    };
    let nc = pack.coeff_len();
    let nz = nc + 1;
    let mut z = pack_orbit(&first.orbit, &HbAnchor::Forced { omega: 1.0 });
    z.push(lambda0);
    let mut points = vec![ScaledContinuationPoint {
        point: ContinuationPoint {
            lambda: lambda0,
            orbit: first.orbit,
            step: 0.0,
        },
        scaled_residual: first.scaled_residual,
        scaled_sampling_defect: first.scaled_sampling_defect,
        arclength_residual: 0.0,
    }];
    let dir = (lambda1 - lambda0).signum();
    let mut tangent = vec![0.0; nz];
    tangent[nc] = dir;
    let mut step = budget.initial_step;
    let mut residual = vec![0.0; nc];
    let mut aft = AftWorkspace::new(&pack)?;
    let state_jacobian = problem.island_state_jacobian();
    let mut tangent_workspace = state_jacobian
        .map(|_| super::tangent::Workspace::new(&pack))
        .transpose()?;
    let eval = |problem: &mut P, z: &[f64], aft: &mut AftWorkspace, out: &mut [f64]| {
        let omega = omega_of(z[nc]);
        if metric.hb.is_some()
            && (z.iter().any(|v| !v.is_finite())
                || !omega.is_finite()
                || omega <= 0.0
                || !(super::TAU / omega).is_finite())
        {
            return false;
        }
        problem.set_parameter(z[nc]);
        hb_residual(problem, &pack, &z[..nc], omega, aft, out);
        metric.hb.is_none() || out.iter().all(|v| v.is_finite())
    };
    for _ in 0..budget.max_steps {
        let lam = z[nc];
        if (dir > 0.0 && lam >= lambda1) || (dir < 0.0 && lam <= lambda1) {
            return Ok(points);
        }
        let mut zp: Vec<f64> = z
            .iter()
            .zip(&tangent)
            .enumerate()
            .map(|(i, (a, t))| a + step * t * metric.column(i, nc))
            .collect();
        let mut converged = None;
        let mut residual_trace = Vec::new();
        for _ in 0..budget.hb.max_newton {
            if !eval(problem, &zp, &mut aft, &mut residual) {
                break;
            }
            let arc = zp
                .iter()
                .zip(&z)
                .zip(&tangent)
                .enumerate()
                .map(|(i, ((a, b), t))| ((a - b) / metric.column(i, nc)) * t)
                .sum::<f64>()
                - step;
            let rn = metric.augmented(&residual, arc);
            residual_trace.push(if metric.hb.is_some() {
                metric.raw(&residual)
            } else {
                rn
            });
            let tolerance = budget.hb.tolerance
                * if metric.hb.is_some() {
                    1.0
                } else {
                    norm(&zp[..nc]).max(1.0)
                };
            if rn < tolerance {
                let (finer, checked, raw_defect) =
                    check_sampling(problem, &pack, &zp[..nc], omega_of(zp[nc]), &aft, &residual)?;
                let checked_norm = metric.augmented(&checked, arc);
                let defect = metric
                    .hb
                    .as_ref()
                    .map_or(raw_defect, |s| s.difference_norm(&checked, &residual));
                if checked_norm < tolerance && defect < tolerance {
                    let raw = if metric.hb.is_some() {
                        metric.raw(&checked)
                    } else {
                        checked_norm
                    };
                    *residual_trace.last_mut().expect("current corrector") = raw;
                    let sampling_defect = if metric.hb.is_some() {
                        length(checked.iter().zip(&residual).map(|(a, b)| a - b))
                    } else {
                        raw_defect
                    };
                    let (coeffs, _) = pack.unpack(&zp[..nc]);
                    converged = Some(ScaledContinuationPoint {
                        point: ContinuationPoint {
                            lambda: zp[nc],
                            step,
                            orbit: HbOrbit {
                                omega: omega_of(zp[nc]),
                                coeffs,
                                residual_trace: std::mem::take(&mut residual_trace),
                                residual: raw,
                                aft_samples: finer.samples,
                                sampling_defect,
                            },
                        },
                        scaled_residual: metric.balance(&checked),
                        scaled_sampling_defect: defect,
                        arclength_residual: arc.abs(),
                    });
                    break;
                }
                if finer.refinement >= MAX_AFT_REFINEMENTS {
                    let physical_residual = if metric.hb.is_some() {
                        metric.raw(&checked)
                    } else {
                        checked_norm
                    };
                    let defect = if metric.hb.is_some() {
                        length(checked.iter().zip(&residual).map(|(a, b)| a - b))
                    } else {
                        raw_defect
                    };
                    return Err(OrbitError::SamplingUnresolved {
                        samples: finer.samples,
                        residual: physical_residual,
                        defect,
                    });
                }
                aft = finer;
                continue;
            }
            let mut jac = vec![0.0; nz * nz];
            let mut r2 = vec![0.0; nc];
            let mut valid = true;
            if let (Some(callback), Some(workspace)) = (state_jacobian, tangent_workspace.as_mut())
            {
                workspace.prepare(problem, callback, &pack, omega_of(zp[nc]), &aft)?;
            }
            for col in 0..nz {
                let unit = metric.column(col, nc);
                if let Some(workspace) = tangent_workspace.as_ref().filter(|_| col < nc) {
                    for row in 0..nc {
                        jac[row * nz + col] =
                            (workspace.coefficient_entry(row, col, &pack) / metric.row(row)) * unit;
                    }
                    jac[nc * nz + col] = tangent[col];
                    continue;
                }
                let h = 1.0e-7 * zp[col].abs().max(unit);
                let saved = zp[col];
                zp[col] = saved + h;
                let difference = if metric.hb.is_some() {
                    zp[col] - saved
                } else {
                    h
                };
                let normalized_difference = difference / unit;
                if metric.hb.is_some()
                    && (!normalized_difference.is_finite() || normalized_difference == 0.0)
                {
                    return Err(invalid(
                        "continuation finite-difference column is not representable",
                    ));
                }
                valid = eval(problem, &zp, &mut aft, &mut r2);
                zp[col] = saved;
                if !valid {
                    break;
                }
                for row in 0..nc {
                    jac[row * nz + col] = if metric.hb.is_some() {
                        (r2[row] / metric.row(row) - residual[row] / metric.row(row))
                            / normalized_difference
                    } else {
                        (r2[row] - residual[row]) / normalized_difference
                    };
                }
                jac[nc * nz + col] = tangent[col];
            }
            eval(problem, &zp, &mut aft, &mut residual);
            if !valid {
                break;
            }
            let Ok(factorized) = lu(&jac, nz) else {
                return Err(OrbitError::SingularJacobian);
            };
            let mut delta: Vec<f64> = residual
                .iter()
                .enumerate()
                .map(|(i, r)| r / metric.row(i))
                .collect();
            delta.push(arc);
            factorized.solve(&mut delta);
            for col in 0..nz {
                zp[col] -= delta[col] * metric.column(col, nc);
            }
        }
        if let Some(point) = converged {
            let diff: Vec<f64> = zp
                .iter()
                .zip(&z)
                .enumerate()
                .map(|(i, (a, b))| (a - b) / metric.column(i, nc))
                .collect();
            let nrm = if metric.hb.is_some() {
                length(diff.iter().copied())
            } else {
                norm(&diff).max(1.0e-300)
            };
            if metric.hb.is_some() && (!nrm.is_finite() || nrm == 0.0) {
                return Err(invalid("continuation secant is not representable"));
            }
            tangent = diff.iter().map(|v| v / nrm).collect();
            z = zp;
            points.push(point);
            step = (step * 1.3).min(budget.initial_step * 4.0);
        } else if step > budget.min_step {
            step *= 0.5;
        } else {
            return Err(OrbitError::ContinuationExhausted {
                steps: points.len(),
            });
        }
    }
    Err(OrbitError::ContinuationExhausted {
        steps: points.len(),
    })
}
