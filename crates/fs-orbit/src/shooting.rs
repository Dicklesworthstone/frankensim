//! Shooting with the exact derivative of its discrete RK4 period map.
//!
//! State and tangent use the same four source stages. The autonomous period
//! column differentiates the step size at every stage, rather than replacing
//! the discrete derivative with the continuous endpoint vector field. A fixed
//! full-state phase plane closes the autonomous equations without discarding
//! a physical row or assuming that the first coordinate moves.

use crate::{OrbitError, OrbitProblem, ShootBudget, ShootOrbit, TimeDomainRhs};
use fs_la::{eigen_complex::eig, factor::lu};
use fs_math::{c64::C64, dd::Dd, det};

mod flow;

/// Explicit physical scales and bounds for a complete variational query.
#[derive(Clone, Debug)]
pub struct VariationalShootControls {
    /// Positive original-state magnitudes. Residual acceptance is absolute in
    /// this chart; a large static preload never relaxes another state's limit.
    pub state_scales: Vec<f64>,
    /// Positive period magnitude [s] for the autonomous Newton column.
    pub period_scale: f64,
    /// Minimum observed scaled excursion of an autonomous cycle, strictly
    /// larger than the normalized closure tolerance. Forced zero and static
    /// periodic responses remain permitted.
    pub minimum_autonomous_excursion: f64,
    /// Complete original state count, additionally bounded by 256.
    pub maximum_dimension: usize,
    /// Conservative live numerical bytes, excluding caller source and allocator.
    pub maximum_payload_bytes: usize,
    /// Admission score, not measured FLOPs: maps*steps*(16d^3+128d^2+128d)
    /// + (128*max_newton+8192)*d^3. Shared dense LU/eigen stages are atomic.
    pub maximum_work: u64,
    /// Counts every complete period map, including rejected line-search trials.
    pub maximum_period_maps: usize,
}

/// Complete numerical periodic response with a consistent discrete tangent.
#[derive(Clone, Debug)]
pub struct VariationalShootOrbit {
    /// Physical period and original state. Its trace retains physical residual
    /// norms; acceptance uses the separate normalized residual below.
    pub orbit: ShootOrbit,
    /// S^-1 D(Phi_T) S, where S=diag(state_scales). This disclosed similarity
    /// chart preserves multipliers without forcing mixed-unit entries to overflow.
    pub monodromy: Vec<f64>,
    pub state_scales: Vec<f64>,
    /// Original d(Phi_T)/dT [state/s], for the autonomous period column only.
    pub period_derivative: Option<Vec<f64>>,
    pub normalized_residual: f64,
    pub phase_residual: f64,
    pub maximum_scaled_excursion: f64,
    pub period_maps: usize,
    pub rhs_evaluations: usize,
    pub jacobian_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum VariationalShootError {
    Invalid(&'static str),
    Cancelled,
    Budget {
        quantity: &'static str,
        requested: u128,
        cap: u128,
    },
    TimeDomainUnavailable,
    TangentUnavailable,
    /// A callback refused, left an entry unwritten, or supplied nonfinite data.
    SourceRefused {
        time: f64,
        tangent: bool,
    },
    Numerical(OrbitError),
}
impl core::fmt::Display for VariationalShootError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "variational shooting refused: {self:?}")
    }
}
impl std::error::Error for VariationalShootError {}
impl From<OrbitError> for VariationalShootError {
    fn from(error: OrbitError) -> Self {
        Self::Numerical(error)
    }
}
type Result<T> = core::result::Result<T, VariationalShootError>;
fn invalid(what: &'static str) -> VariationalShootError {
    VariationalShootError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(VariationalShootError::Cancelled)
    } else {
        Ok(())
    }
}
fn cap(quantity: &'static str, requested: u128, limit: u128) -> Result<()> {
    if requested > limit {
        Err(VariationalShootError::Budget {
            quantity,
            requested,
            cap: limit,
        })
    } else {
        Ok(())
    }
}
fn length(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |norm, &v| det::hypot(norm, v))
}
fn phase(normal: &[f64], value: &[f64], origin: &[f64]) -> f64 {
    normal
        .iter()
        .zip(value)
        .zip(origin)
        .fold(Dd::ZERO, |sum, ((&n, &x), &r)| {
            sum + Dd::from_f64(n) * (Dd::from_f64(x) - Dd::from_f64(r))
        })
        .to_f64()
}
fn defect(end: &[f64], start: &[f64]) -> f64 {
    end.iter()
        .zip(start)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

struct Counts {
    maps: usize,
    maximum: usize,
    rhs: usize,
    jacobian: usize,
}
impl Counts {
    fn begin(&mut self) -> Result<()> {
        cap("period maps", self.maps as u128 + 1, self.maximum as u128)?;
        self.maps += 1;
        Ok(())
    }
}

/// Resolve a smooth forced or autonomous periodic response and its Floquet map.
/// Requires the COMPLETE explicit RHS and its original exact state Jacobian;
/// no frequency port is inverted or assumed to have an available realization.
///
/// The fixed-step RK4 tangent is obtained by differentiating the very same
/// stages as the primal map. State and unknown-period columns share these
/// evaluations. Newton uses the existing dense LU, and multipliers use the
/// existing complex eigensolver. There are no trajectory perturbations, extra
/// model copies, state reduction or artificial regularization.
///
/// Autonomous phase is a fixed plane through the seed, normal to its scaled
/// vector field, using ALL original states. The source must actually be
/// autonomous: time derivatives in the period column are then zero. A declared
/// excursion floor excludes equilibrium. Forced sources use their prescribed
/// physical period and can have exact-zero response or a turning point at t=0.
///
/// Results concern this discrete smooth map. Temporal refinement, primitive
/// period, source equivalence, continuous-orbit existence, Floquet conditioning
/// and physical validity require their own evidence. Multipliers are reported
/// without inferring torus existence, global attraction or hybrid saltations.
/// Cancellation is polled at source stages and matrix-product boundaries; the
/// separately bounded shared dense LU/eigen calls remain atomic.
#[allow(clippy::too_many_arguments)]
pub fn solve_shooting_variational<P: OrbitProblem>(
    problem: &P,
    seed: &[f64],
    period_guess: f64,
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit> {
    poll(&mut cancelled)?;
    let n = problem.dim();
    if n == 0
        || n > 256
        || n > controls.maximum_dimension
        || seed.len() != n
        || controls.state_scales.len() != n
        || seed.iter().any(|v| !v.is_finite())
        || controls
            .state_scales
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        || !period_guess.is_finite()
        || period_guess <= 0.0
        || !controls.period_scale.is_finite()
        || controls.period_scale <= 0.0
        || !controls.minimum_autonomous_excursion.is_finite()
        || controls.minimum_autonomous_excursion <= 0.0
        || budget.steps_per_period == 0
        || budget.max_newton == 0
        || budget.max_newton > 4096
        || !budget.tolerance.is_finite()
        || budget.tolerance <= 0.0
        || budget.tolerance >= 1.0
        || controls.maximum_period_maps == 0
    {
        return Err(invalid(
            "complete finite states, positive physical scales and bounded shooting controls required",
        ));
    }
    let d = n as u128;
    let bytes = 384 * d * d + 1024 * d + 64 * (budget.max_newton as u128 + 1);
    let work = (controls.maximum_period_maps as u128)
        .checked_mul(budget.steps_per_period as u128)
        .and_then(|v| v.checked_mul(16 * d * d * d + 128 * d * d + 128 * d))
        .and_then(|v| v.checked_add((128 * budget.max_newton as u128 + 8192) * d * d * d))
        .ok_or_else(|| invalid("variational work envelope overflows"))?;
    cap(
        "numerical payload bytes",
        bytes,
        controls.maximum_payload_bytes as u128,
    )?;
    cap("variational work", work, controls.maximum_work as u128)?;
    // Counters must remain representable on 32-bit hosts as well. The checked
    // work envelope above already bounds this smaller product in u128.
    cap(
        "source evaluations",
        4 * controls.maximum_period_maps as u128 * budget.steps_per_period as u128 + 1,
        usize::MAX as u128,
    )?;
    let auto = problem.autonomous();
    if auto && controls.minimum_autonomous_excursion <= budget.tolerance {
        return Err(invalid(
            "autonomous excursion floor must exceed the normalized closure tolerance",
        ));
    }
    let rhs = problem
        .time_domain_rhs()
        .ok_or(VariationalShootError::TimeDomainUnavailable)?;
    let jacobian = problem
        .time_domain_state_jacobian()
        .ok_or(VariationalShootError::TangentUnavailable)?;
    let columns = n + usize::from(auto);
    let units = &controls.state_scales;
    let origin: Vec<_> = seed.iter().zip(units).map(|(x, s)| x / s).collect();
    let mut z = origin.clone();
    let mut period = period_guess;
    if z.iter().any(|v| !v.is_finite()) || !(period / controls.period_scale).is_finite() {
        return Err(invalid(
            "normalized physical seed and period must be representable",
        ));
    }
    let mut normal = vec![0.0; n];
    let mut counts = Counts {
        maps: 0,
        maximum: controls.maximum_period_maps,
        rhs: 0,
        jacobian: 0,
    };
    if auto {
        normal.fill(f64::NAN);
        rhs(problem, 0.0, seed, &mut normal);
        counts.rhs += 1;
        poll(&mut cancelled)?;
        if normal.iter().any(|v| !v.is_finite()) {
            return Err(VariationalShootError::SourceRefused {
                time: 0.0,
                tangent: false,
            });
        }
        for (v, &s) in normal.iter_mut().zip(units) {
            *v = (*v / s) * controls.period_scale;
        }
        let norm = length(&normal);
        if norm == 0.0 {
            return Err(OrbitError::TrivialCollapse.into());
        }
        if !norm.is_finite() {
            return Err(invalid("scaled phase direction is not representable"));
        }
        for v in &mut normal {
            *v /= norm;
        }
    }
    let mut workspace = flow::Workspace::new(n, columns);
    let mut trace = Vec::with_capacity(budget.max_newton + 1);
    let mut last = f64::INFINITY;
    for iteration in 0..=budget.max_newton {
        let image = workspace.integrate(
            problem,
            rhs,
            Some(jacobian),
            &z,
            period,
            budget.steps_per_period,
            units,
            controls.period_scale,
            &mut counts,
            &mut cancelled,
        )?;
        let residual: Vec<_> = image.end.iter().zip(&z).map(|(a, b)| a - b).collect();
        let physical: Vec<_> = residual.iter().zip(units).map(|(r, s)| r * s).collect();
        if physical.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "original periodicity residual is not representable",
            ));
        }
        trace.push(length(&physical));
        let phase_value = if auto {
            phase(&normal, &z, &origin)
        } else {
            0.0
        };
        let normalized_residual = defect(&image.end, &z);
        last = normalized_residual.max(phase_value.abs());
        if last <= budget.tolerance {
            if auto && image.excursion < controls.minimum_autonomous_excursion {
                return Err(OrbitError::TrivialCollapse.into());
            }
            let monodromy: Vec<_> = (0..n)
                .flat_map(|i| image.tangent[i * columns..i * columns + n].iter().copied())
                .collect();
            let complex: Vec<_> = monodromy.iter().map(|&v| C64::from_re(v)).collect();
            poll(&mut cancelled)?;
            let multipliers = eig(&complex, n).map_err(OrbitError::Eigen)?;
            poll(&mut cancelled)?;
            let x0: Vec<_> = z.iter().zip(units).map(|(x, s)| x * s).collect();
            let period_derivative = auto.then(|| {
                (0..n)
                    .map(|i| (image.tangent[i * columns + n] / controls.period_scale) * units[i])
                    .collect::<Vec<_>>()
            });
            if x0
                .iter()
                .chain(period_derivative.iter().flatten())
                .any(|v| !v.is_finite())
                || multipliers
                    .iter()
                    .any(|v| !v.re.is_finite() || !v.im.is_finite())
            {
                return Err(invalid(
                    "original orbit or Floquet result is not representable",
                ));
            }
            return Ok(VariationalShootOrbit {
                orbit: ShootOrbit {
                    period,
                    x0,
                    multipliers: multipliers.iter().map(|v| (v.re, v.im)).collect(),
                    residual_trace: trace,
                },
                monodromy,
                state_scales: units.clone(),
                period_derivative,
                normalized_residual,
                phase_residual: phase_value.abs(),
                maximum_scaled_excursion: image.excursion,
                period_maps: counts.maps,
                rhs_evaluations: counts.rhs,
                jacobian_evaluations: counts.jacobian,
            });
        }
        if iteration == budget.max_newton {
            break;
        }
        let mut matrix = vec![0.0; columns * columns];
        let mut correction = vec![0.0; columns];
        for i in 0..n {
            poll(&mut cancelled)?;
            correction[i] = residual[i];
            for j in 0..columns {
                matrix[i * columns + j] = image.tangent[i * columns + j];
            }
            matrix[i * columns + i] -= 1.0;
        }
        if auto {
            matrix[n * columns..n * columns + n].copy_from_slice(&normal);
            correction[n] = phase_value;
        }
        poll(&mut cancelled)?;
        let factor = lu(&matrix, columns).map_err(|_| OrbitError::SingularJacobian)?;
        factor.solve(&mut correction);
        poll(&mut cancelled)?;
        if correction.iter().any(|v| !v.is_finite()) {
            return Err(invalid("shooting correction is not representable"));
        }
        let mut accepted = false;
        for trial in 0..8 {
            let alpha = 2.0_f64.powi(-trial);
            let next: Vec<_> = z
                .iter()
                .zip(&correction)
                .map(|(x, dx)| x - alpha * dx)
                .collect();
            let next_period = if auto {
                period - alpha * correction[n] * controls.period_scale
            } else {
                period
            };
            if next.iter().any(|v| !v.is_finite()) || !next_period.is_finite() || next_period <= 0.0
            {
                continue;
            }
            let candidate = workspace.integrate(
                problem,
                rhs,
                None,
                &next,
                next_period,
                budget.steps_per_period,
                units,
                controls.period_scale,
                &mut counts,
                &mut cancelled,
            )?;
            let merit = defect(&candidate.end, &next).max(if auto {
                phase(&normal, &next, &origin).abs()
            } else {
                0.0
            });
            if merit < last {
                z = next;
                period = next_period;
                accepted = true;
                break;
            }
        }
        if !accepted {
            return Err(OrbitError::NewtonStalled {
                residual: last,
                iterations: iteration + 1,
                trace,
            }
            .into());
        }
    }
    Err(OrbitError::NewtonStalled {
        residual: last,
        iterations: budget.max_newton,
        trace,
    }
    .into())
}

// The legacy shooter uses exactly the same primal RK4 stages, without tangent
// storage or a change to its original phase/finite-difference policy.
pub(super) fn legacy_flow<P: OrbitProblem>(
    problem: &P,
    rhs: TimeDomainRhs<P>,
    initial: &[f64],
    period: f64,
    steps: usize,
) -> core::result::Result<Vec<f64>, OrbitError> {
    let n = problem.dim();
    if n == 0 || initial.len() != n || steps == 0 {
        return Err(OrbitError::BadParameter {
            what: "complete nonempty time-domain state and positive step count required",
        });
    }
    let mut workspace = flow::Workspace::new(n, 0);
    let mut counts = Counts {
        maps: 0,
        maximum: usize::MAX,
        rhs: 0,
        jacobian: 0,
    };
    workspace
        .integrate(
            problem,
            rhs,
            None,
            initial,
            period,
            steps,
            &vec![1.0; n],
            1.0,
            &mut counts,
            &mut || false,
        )
        .map(|image| image.end)
        .map_err(|_| OrbitError::BadParameter {
            what: "time-domain period map must be finite and representable",
        })
}
