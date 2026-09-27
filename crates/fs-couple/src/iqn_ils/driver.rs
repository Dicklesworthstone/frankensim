//! Bounded, transactional partitioned physical steps over an actual domain map.
//!
//! Every trial starts from the same committed physical state. Acceptance checks
//! the unrelaxed `G(x) - x` coordinate by coordinate, and each declared balance
//! separately. Neither tiny relaxation nor cancelling signed residuals can
//! establish convergence. IQN-ILS uses fixed, caller-declared coordinate scales.
//!
//! The producer must return a disposable, owned trial state without mutating
//! committed state through interior mutability or performing irreversible writes.
//! It owns domain admissibility and balance evaluation. The cancellation callback
//! should poll the owning execution context, e.g. `|| cx.checkpoint().is_err()`;
//! producers must also poll inside their kernels. This driver polls around each
//! producer call and acceleration, not inside either operation.
//!
//! Passing these numerical gates is not a physical-validity, passivity, stability,
//! discretization-error, or port-schema certificate. Empty balance controls are
//! an explicit decision to perform no conservation gate.

/// Bounded physical-time continuation with committed-step checkpoints.
pub mod march;

use core::fmt;
use std::collections::BTreeSet;

use super::{IqnIls, IqnIlsConfig, IqnIlsError};

/// Explicit physical interval, validated before any producer runs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepInterval {
    /// Start time in seconds.
    pub start_s: f64,
    /// End time in seconds, strictly later than the start.
    pub end_s: f64,
}

impl StepInterval {
    /// Physical duration in seconds.
    #[must_use]
    pub fn duration_s(self) -> f64 {
        self.end_s - self.start_s
    }
}

/// Fixed units and tolerance for one interface coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceControl {
    /// Finite positive scale in the coordinate's native units.
    pub scale: f64,
    /// Finite nonnegative absolute tolerance in native units.
    pub absolute_tolerance: f64,
    /// Finite nonnegative tolerance relative to `scale`, not a growing iterate.
    /// The combined native-unit tolerance must be representable as a finite f64.
    pub relative_tolerance: f64,
}

/// One separately checked accounting equation.
#[derive(Debug, Clone, PartialEq)]
pub struct BalanceControl {
    /// Unique nonempty name; callers should include the physical quantity/units.
    pub name: String,
    /// Finite nonnegative tolerance in this balance's native units.
    pub absolute_tolerance: f64,
}

/// Numerical update policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CouplingMethod {
    /// Fixed under-relaxation, without secant history.
    RelaxedPicard,
    /// Bounded vector IQN-ILS with rank-filtered history.
    IqnIls(IqnIlsConfig),
}

/// Explicit work budget and acceptance policy.
#[derive(Debug, Clone, PartialEq)]
pub struct CouplingControls {
    /// Maximum actual producer calls, including the initial trial; nonzero.
    pub max_evaluations: usize,
    /// Picard/startup relaxation, finite and in `(0, 1]`.
    pub relaxation: f64,
    /// Update policy.
    pub method: CouplingMethod,
    /// One rule per interface coordinate, in producer order.
    pub interfaces: Vec<InterfaceControl>,
    /// Required balances, in producer order; empty means no balance gate.
    pub balances: Vec<BalanceControl>,
}

/// A tentative physical state and the map/balances evaluated with it.
#[derive(Debug, Clone, PartialEq)]
pub struct CouplingTrial<S> {
    /// Owned state, published only after all gates pass.
    pub state: S,
    /// Actual `G(x)` in native coordinate units.
    pub image: Vec<f64>,
    /// Signed residuals, one per declared balance, in declared order.
    pub balance_residuals: Vec<f64>,
}

/// Diagnostics for one completely measured trial.
#[derive(Debug, Clone, PartialEq)]
pub struct IterationResidual {
    /// One-based producer call number.
    pub evaluation: usize,
    /// Rounded infinity norm of the residual divided by fixed scales.
    /// The acceptance gate compares native-unit residuals, so scaling underflow
    /// here cannot promote a nonzero residual to exact convergence.
    pub maximum_scaled_residual: f64,
    /// Number of coordinates outside their tolerances.
    pub unconverged_coordinates: usize,
    /// Number of independently failed balance equations.
    pub failed_balances: usize,
    /// Independent secants used for the subsequent proposal, or zero.
    pub acceleration_columns: usize,
}

/// Numerical evidence only, not a physical certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct CouplingReport {
    /// Attempted physical interval.
    pub interval: StepInterval,
    /// Actual producer calls, including a call that refused or was cancelled.
    pub evaluations: usize,
    /// Completely measured trials, bounded by `max_evaluations`.
    pub iterations: Vec<IterationResidual>,
    /// Last completely measured `G(x) - x`, in native units.
    pub interface_residuals: Vec<f64>,
    /// Last completely measured balance residuals.
    pub balance_residuals: Vec<f64>,
}

/// Invalid controls, malformed producer output or unrepresentable arithmetic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CouplingInputError {
    /// A control violated its documented domain.
    Control {
        /// Responsible field.
        field: &'static str,
        /// Coordinate/rule index, or zero for a scalar control.
        index: usize,
    },
    /// A vector had the wrong shape.
    Shape {
        /// Responsible field.
        field: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        found: usize,
    },
    /// A nonfinite input or arithmetic result.
    NonFinite {
        /// Responsible stage.
        field: &'static str,
        /// Coordinate index.
        index: usize,
        /// IEEE-754 representation of the rejected value.
        bits: u64,
    },
}

/// Why a physical step was not committed.
#[derive(Debug, Clone, PartialEq)]
pub enum CouplingFailure<E> {
    /// Invalid controls, producer shape or arithmetic.
    Input(CouplingInputError),
    /// Original domain-producer error.
    Operator(E),
    /// Original IQN-ILS error.
    Accelerator(IqnIlsError),
    /// Cancellation was observed at an explicit checkpoint.
    Cancelled,
    /// The work budget ended without satisfying every gate.
    NotConverged,
}

/// A refusal and its measured prefix; never a partially published state.
#[derive(Debug, Clone, PartialEq)]
pub struct CouplingError<E> {
    /// Original typed failure.
    pub reason: CouplingFailure<E>,
    /// Diagnostics from the attempted step.
    pub report: CouplingReport,
}

impl<E: fmt::Debug> fmt::Display for CouplingError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "coupled step refused after {} evaluations: {:?}", self.report.evaluations, self.reason)
    }
}

impl<E: std::error::Error + 'static> std::error::Error for CouplingError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.reason {
            CouplingFailure::Operator(error) => Some(error),
            CouplingFailure::Accelerator(error) => Some(error),
            _ => None,
        }
    }
}

impl CouplingReport {
    fn error<E>(&self, reason: CouplingFailure<E>) -> CouplingError<E> {
        CouplingError { reason, report: self.clone() }
    }
}

fn control(valid: bool, field: &'static str, index: usize) -> Result<(), CouplingInputError> {
    if valid { Ok(()) } else { Err(CouplingInputError::Control { field, index }) }
}

fn finite(value: f64, field: &'static str, index: usize) -> Result<f64, CouplingInputError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CouplingInputError::NonFinite { field, index, bits: value.to_bits() })
    }
}

fn shape(field: &'static str, expected: usize, found: usize) -> Result<(), CouplingInputError> {
    if expected == found {
        Ok(())
    } else {
        Err(CouplingInputError::Shape { field, expected, found })
    }
}

fn validate(
    interval: StepInterval,
    interface: &[f64],
    controls: &CouplingControls,
) -> Result<Vec<f64>, CouplingInputError> {
    finite(interval.start_s, "start time", 0)?;
    finite(interval.end_s, "end time", 0)?;
    finite(interval.duration_s(), "duration", 0)?;
    control(interval.end_s > interval.start_s, "duration", 0)?;
    control(controls.max_evaluations > 0, "max evaluations", 0)?;
    control(
        controls.relaxation.is_finite() && controls.relaxation > 0.0 && controls.relaxation <= 1.0,
        "relaxation", 0,
    )?;
    control(!interface.is_empty(), "interface", 0)?;
    shape("interface controls", interface.len(), controls.interfaces.len())?;
    let mut thresholds = Vec::with_capacity(interface.len());
    for (index, (&value, rule)) in interface.iter().zip(&controls.interfaces).enumerate() {
        finite(value, "initial interface", index)?;
        control(rule.scale.is_finite() && rule.scale > 0.0, "scale", index)?;
        control(rule.absolute_tolerance.is_finite() && rule.absolute_tolerance >= 0.0, "absolute tolerance", index)?;
        control(rule.relative_tolerance.is_finite() && rule.relative_tolerance >= 0.0, "relative tolerance", index)?;
        // Gate in native units: dividing a tiny nonzero residual by a large
        // scale can underflow to zero, which must not satisfy a zero tolerance.
        thresholds.push(finite(
            rule.absolute_tolerance + rule.relative_tolerance * rule.scale,
            "combined tolerance", index,
        )?);
    }
    let mut names = BTreeSet::new();
    for (index, rule) in controls.balances.iter().enumerate() {
        control(!rule.name.trim().is_empty() && names.insert(rule.name.as_str()), "balance name", index)?;
        control(rule.absolute_tolerance.is_finite() && rule.absolute_tolerance >= 0.0, "balance tolerance", index)?;
    }
    Ok(thresholds)
}

/// Execute one physical step without publishing any unsuccessful trial.
///
/// `evaluate` receives the same committed state at every nonlinear iteration.
/// It must evaluate the supplied interface, not advance from a rejected trial.
/// `cancelled` returns true when the owning execution context requests a stop.
///
/// On success, `state` becomes the accepted trial state and `interface` becomes
/// the INPUT used to compute it, not an unevaluated proposal or the slightly
/// different map image. The report retains that remaining defect. On every
/// returned error both caller-owned values remain unchanged, provided callbacks
/// obey the no-interior-mutation/no-irreversible-side-effects contract.
///
/// # Errors
/// Invalid controls/output, nonfinite arithmetic, producer or IQN-ILS refusal,
/// observed cancellation, or exhausted producer-call budget.
pub fn coupled_step<S, E, F, C>(
    state: &mut S,
    interface: &mut [f64],
    interval: StepInterval,
    controls: &CouplingControls,
    evaluate: &mut F,
    cancelled: &mut C,
) -> Result<CouplingReport, CouplingError<E>>
where
    F: FnMut(&S, StepInterval, &[f64]) -> Result<CouplingTrial<S>, E>,
    C: FnMut() -> bool,
{
    let mut report = CouplingReport {
        interval, evaluations: 0, iterations: Vec::new(),
        interface_residuals: Vec::new(), balance_residuals: Vec::new(),
    };
    let thresholds = validate(interval, interface, controls)
        .map_err(|error| report.error(CouplingFailure::Input(error)))?;
    let mut accelerator = match controls.method {
        CouplingMethod::RelaxedPicard => None,
        CouplingMethod::IqnIls(config) => Some(IqnIls::new(interface.len(), config)
            .map_err(|error| report.error(CouplingFailure::Accelerator(error)))?),
    };
    let origin = interface.to_vec();
    let mut current = origin.clone();
    let mut normalized_current = vec![0.0; current.len()];
    let mut normalized_image = vec![0.0; current.len()];
    for iteration in 0..controls.max_evaluations {
        if cancelled() { return Err(report.error(CouplingFailure::Cancelled)); }
        report.evaluations += 1;
        let trial = evaluate(state, interval, &current)
            .map_err(|error| report.error(CouplingFailure::Operator(error)))?;
        if cancelled() { return Err(report.error(CouplingFailure::Cancelled)); }

        // Malformed trials must not replace diagnostics from a valid prefix.
        let measure = || -> Result<(Vec<f64>, f64, usize, usize), CouplingInputError> {
            shape("map image", current.len(), trial.image.len())?;
            shape("balances", controls.balances.len(), trial.balance_residuals.len())?;
            let mut residuals = Vec::with_capacity(current.len());
            let mut maximum = 0.0_f64;
            let mut missed = 0;
            for (index, (&mapped, &x)) in trial.image.iter().zip(&current).enumerate() {
                finite(mapped, "map image", index)?;
                let residual = finite(mapped - x, "interface residual", index)?;
                let scaled = finite(residual / controls.interfaces[index].scale, "scaled residual", index)?;
                maximum = maximum.max(scaled.abs());
                missed += usize::from(residual.abs() > thresholds[index]);
                residuals.push(residual);
            }
            let mut failed = 0;
            for (index, (&value, rule)) in trial.balance_residuals.iter().zip(&controls.balances).enumerate() {
                finite(value, "balance residual", index)?;
                failed += usize::from(value.abs() > rule.absolute_tolerance);
            }
            Ok((residuals, maximum, missed, failed))
        };
        let (residuals, maximum, missed, failed_balances) = measure()
            .map_err(|error| report.error(CouplingFailure::Input(error)))?;
        report.interface_residuals = residuals;
        report.balance_residuals = trial.balance_residuals.clone();
        report.iterations.push(IterationResidual {
            evaluation: iteration + 1, maximum_scaled_residual: maximum,
            unconverged_coordinates: missed, failed_balances, acceleration_columns: 0,
        });
        if cancelled() { return Err(report.error(CouplingFailure::Cancelled)); }
        if missed == 0 && failed_balances == 0 {
            *state = trial.state;
            interface.copy_from_slice(&current);
            return Ok(report);
        }
        if iteration + 1 == controls.max_evaluations { break; }
        for index in 0..current.len() {
            normalized_current[index] = finite(
                (current[index] - origin[index]) / controls.interfaces[index].scale,
                "normalized iterate", index,
            ).map_err(|error| report.error(CouplingFailure::Input(error)))?;
            normalized_image[index] = finite(
                (trial.image[index] - origin[index]) / controls.interfaces[index].scale,
                "normalized image", index,
            ).map_err(|error| report.error(CouplingFailure::Input(error)))?;
        }
        let (proposal, columns) = if let Some(accelerator) = &mut accelerator {
            let step = accelerator.step(&normalized_current, &normalized_image, controls.relaxation)
                .map_err(|error| report.error(CouplingFailure::Accelerator(error)))?;
            (step.values, step.used_columns)
        } else {
            let values = normalized_current.iter().zip(&normalized_image).enumerate()
                .map(|(index, (&x, &g))| finite(x + controls.relaxation * (g - x), "Picard proposal", index))
                .collect::<Result<Vec<_>, _>>()
                .map_err(|error| report.error(CouplingFailure::Input(error)))?;
            (values, 0)
        };
        if cancelled() { return Err(report.error(CouplingFailure::Cancelled)); }
        for (index, value) in proposal.into_iter().enumerate() {
            current[index] = finite(origin[index] + value * controls.interfaces[index].scale, "physical proposal", index)
                .map_err(|error| report.error(CouplingFailure::Input(error)))?;
        }
        if let Some(last) = report.iterations.last_mut() { last.acceleration_columns = columns; }
    }
    Err(report.error(CouplingFailure::NotConverged))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    pub(super) fn controls(n: usize, evaluations: usize) -> CouplingControls {
        CouplingControls {
            max_evaluations: evaluations, relaxation: 0.5,
            method: CouplingMethod::IqnIls(IqnIlsConfig::default()),
            interfaces: vec![InterfaceControl {
                scale: 1.0, absolute_tolerance: 1.0e-12, relative_tolerance: 0.0,
            }; n],
            balances: Vec::new(),
        }
    }

    fn interval() -> StepInterval { StepInterval { start_s: 0.0, end_s: 0.5 } }

    #[test]
    fn thermal_exchange_matches_monolithic_backward_euler() {
        let initial = [400.0, 300.0];
        let mut state = initial;
        let mut interface = initial;
        let mut policy = controls(2, 12);
        for rule in &mut policy.interfaces { rule.scale = 100.0; rule.absolute_tolerance = 1.0e-10; }
        policy.balances.push(BalanceControl { name: "thermal-energy-j".into(), absolute_tolerance: 1.0e-8 });
        let mut calls = 0;
        let report = coupled_step(&mut state, &mut interface, interval(), &policy,
            &mut |old: &[f64; 2], step, x: &[f64]| {
                assert_eq!(*old, initial, "rejected trials must not advance physical time");
                calls += 1;
                let dt = step.duration_s();
                let next = [
                    (2.0 / dt * old[0] + 10.0 * x[1]) / (2.0 / dt + 10.0),
                    (3.0 / dt * old[1] + 10.0 * x[0]) / (3.0 / dt + 10.0),
                ];
                Ok::<_, &'static str>(CouplingTrial { state: next, image: next.to_vec(),
                    balance_residuals: vec![2.0 * (next[0] - old[0]) + 3.0 * (next[1] - old[1])] })
            }, &mut || false,
        ).unwrap();
        let delta = 100.0 / (1.0 + 0.5 * 10.0 * (0.5 + 1.0 / 3.0));
        assert!((state[0] - (340.0 + 0.6 * delta)).abs() < 1.0e-9);
        assert!((state[1] - (340.0 - 0.4 * delta)).abs() < 1.0e-9);
        assert_eq!(calls, report.evaluations);
        assert!(report.iterations.iter().any(|row| row.acceleration_columns > 0));
        assert!(report.balance_residuals[0].abs() <= 1.0e-8);
    }

    #[test]
    fn opposing_modes_cannot_cancel_and_budget_failure_rolls_back() {
        let mut state = 42;
        let mut interface = [0.0, 0.0];
        let error = coupled_step(&mut state, &mut interface, interval(), &controls(2, 1),
            &mut |_, _, _| Ok::<_, &'static str>(CouplingTrial {
                state: 43, image: vec![1.0, -1.0], balance_residuals: vec![],
            }), &mut || false,
        ).unwrap_err();
        assert!(matches!(error.reason, CouplingFailure::NotConverged));
        assert_eq!(error.report.iterations[0].unconverged_coordinates, 2);
        assert_eq!((state, interface), (42, [0.0, 0.0]));
    }

    #[test]
    fn tiny_relaxation_cannot_hide_a_large_fixed_point_defect() {
        let mut policy = controls(1, 3);
        policy.relaxation = 1.0e-20;
        policy.method = CouplingMethod::RelaxedPicard;
        let mut state = 0;
        let mut interface = [0.0];
        let error = coupled_step(&mut state, &mut interface, interval(), &policy,
            &mut |_, _, _| Ok::<_, &'static str>(CouplingTrial {
                state: 1, image: vec![1.0], balance_residuals: vec![],
            }), &mut || false,
        ).unwrap_err();
        assert!(matches!(error.reason, CouplingFailure::NotConverged));
        assert_eq!(error.report.evaluations, 3);
        assert_eq!((state, interface), (0, [0.0]));
    }

    #[test]
    fn independently_bad_balances_cannot_cancel_at_a_fixed_point() {
        let mut policy = controls(1, 2);
        policy.balances = ["left-j", "right-j"].into_iter().map(|name| BalanceControl {
            name: name.into(), absolute_tolerance: 0.0,
        }).collect();
        let mut state = 0;
        let mut interface = [0.0];
        let error = coupled_step(&mut state, &mut interface, interval(), &policy,
            &mut |_, _, x| Ok::<_, &'static str>(CouplingTrial {
                state: 1, image: x.to_vec(), balance_residuals: vec![1.0, -1.0],
            }), &mut || false,
        ).unwrap_err();
        assert!(matches!(error.reason, CouplingFailure::NotConverged));
        assert!(error.report.iterations.iter().all(|row| row.unconverged_coordinates == 0 && row.failed_balances == 2));
        assert_eq!(state, 0);
    }

    #[test]
    fn producer_failure_preserves_original_error_and_committed_values() {
        let mut state = 7;
        let mut interface = [0.0];
        let mut calls = 0;
        let error = coupled_step(&mut state, &mut interface, interval(), &controls(1, 4),
            &mut |old, _, _| {
                assert_eq!(*old, 7);
                calls += 1;
                if calls == 2 { Err("domain refusal") } else {
                    Ok(CouplingTrial { state: 8, image: vec![1.0], balance_residuals: vec![] })
                }
            }, &mut || false,
        ).unwrap_err();
        assert!(matches!(error.reason, CouplingFailure::Operator("domain refusal")));
        assert_eq!(error.report.evaluations, 2);
        assert_eq!(error.report.iterations.len(), 1);
        assert_eq!((state, interface), (7, [0.0]));
    }

    #[test]
    fn malformed_or_nonfinite_results_never_publish_state() {
        for image in [vec![], vec![f64::NAN], vec![f64::INFINITY]] {
            let mut state = 0;
            let mut interface = [0.0];
            let error = coupled_step(&mut state, &mut interface, interval(), &controls(1, 2),
                &mut |_, _, _| Ok::<_, &'static str>(CouplingTrial {
                    state: 1, image: image.clone(), balance_residuals: vec![],
                }), &mut || false,
            ).unwrap_err();
            assert!(matches!(error.reason, CouplingFailure::Input(_)));
            assert!(error.report.iterations.is_empty());
            assert_eq!((state, interface), (0, [0.0]));
        }
        for balances in [vec![], vec![f64::NAN], vec![f64::INFINITY]] {
            let mut policy = controls(1, 1);
            policy.balances.push(BalanceControl { name: "energy-j".into(), absolute_tolerance: 0.0 });
            let mut state = 0;
            let mut interface = [0.0];
            let error = coupled_step(&mut state, &mut interface, interval(), &policy,
                &mut |_, _, x| Ok::<_, &'static str>(CouplingTrial {
                    state: 1, image: x.to_vec(), balance_residuals: balances.clone(),
                }), &mut || false,
            ).unwrap_err();
            assert!(matches!(error.reason, CouplingFailure::Input(_)));
            assert_eq!(state, 0);
        }
    }

    #[test]
    fn accepted_interface_is_the_one_that_produced_the_trial() {
        let mut policy = controls(1, 1);
        policy.interfaces[0].absolute_tolerance = 0.1;
        let mut state = 0.0;
        let mut interface = [1.0];
        let report = coupled_step(&mut state, &mut interface, interval(), &policy,
            &mut |_, _, x| Ok::<_, &'static str>(CouplingTrial {
                state: x[0], image: vec![1.05], balance_residuals: vec![],
            }), &mut || false,
        ).unwrap();
        assert_eq!((state, interface), (1.0, [1.0]));
        assert!(report.interface_residuals[0] > 0.0);
    }

    #[test]
    fn invalid_controls_are_refused_before_producer_evaluation() {
        for invalid in 0..7 {
            let mut policy = controls(1, 2);
            match invalid {
                0 => policy.max_evaluations = 0,
                1 => policy.relaxation = 0.0,
                2 => policy.interfaces[0].scale = 0.0,
                3 => policy.interfaces[0].absolute_tolerance = f64::NAN,
                4 => policy.interfaces.clear(),
                5 => policy.interfaces[0].relative_tolerance = -1.0,
                _ => policy.balances = vec![BalanceControl { name: "same".into(), absolute_tolerance: 0.0 }; 2],
            }
            let mut state = 0;
            let mut interface = [0.0];
            let error = coupled_step(&mut state, &mut interface, interval(), &policy,
                &mut |_, _, _| -> Result<CouplingTrial<i32>, &'static str> { panic!("invalid controls reached physics") },
                &mut || false,
            ).unwrap_err();
            assert!(matches!(error.reason, CouplingFailure::Input(_)));
            assert_eq!(error.report.evaluations, 0);
        }
    }

    #[test]
    fn iqn_respects_disparate_declared_coordinate_scales() {
        let scales = [1.0e9, 1.0e-9];
        let mut policy = controls(2, 12);
        for (rule, scale) in policy.interfaces.iter_mut().zip(scales) {
            rule.scale = scale; rule.absolute_tolerance = 0.0; rule.relative_tolerance = 1.0e-11;
        }
        let mut state = [0.0, 0.0];
        let mut interface = state;
        coupled_step(&mut state, &mut interface, interval(), &policy,
            &mut |_, _, x| {
                let next = [(0.7 * x[0] / scales[0] + 2.0) * scales[0], (0.2 * x[1] / scales[1] - 3.0) * scales[1]];
                Ok::<_, &'static str>(CouplingTrial { state: next, image: next.to_vec(), balance_residuals: vec![] })
            }, &mut || false,
        ).unwrap();
        assert!((state[0] / scales[0] - 2.0 / 0.3).abs() < 1.0e-9);
        assert!((state[1] / scales[1] + 3.0 / 0.8).abs() < 1.0e-9);
    }

    #[test]
    fn cancellation_before_or_during_producer_never_commits_a_trial() {
        for initially_cancelled in [true, false] {
            let requested = Cell::new(initially_cancelled);
            let calls = Cell::new(0);
            let mut state = 9;
            let mut interface = [0.0];
            let error = coupled_step(&mut state, &mut interface, interval(), &controls(1, 2),
                &mut |_, _, x| {
                    calls.set(calls.get() + 1);
                    requested.set(true);
                    Ok::<_, &'static str>(CouplingTrial { state: 10, image: x.to_vec(), balance_residuals: vec![] })
                }, &mut || requested.get(),
            ).unwrap_err();
            assert!(matches!(error.reason, CouplingFailure::Cancelled));
            assert_eq!(calls.get(), usize::from(!initially_cancelled));
            assert_eq!(error.report.evaluations, calls.get());
            assert_eq!((state, interface), (9, [0.0]));
        }
    }

    #[test]
    fn scaling_underflow_cannot_satisfy_an_exact_zero_tolerance() {
        let mut policy = controls(1, 1);
        policy.interfaces[0] = InterfaceControl { scale: 1.0e300, absolute_tolerance: 0.0, relative_tolerance: 0.0 };
        let mut state = 7;
        let mut interface = [0.0];
        let error = coupled_step(&mut state, &mut interface, interval(), &policy,
            &mut |_, _, _| Ok::<_, &'static str>(CouplingTrial {
                state: 8, image: vec![1.0e-300], balance_residuals: vec![],
            }), &mut || false,
        ).unwrap_err();
        assert!(matches!(error.reason, CouplingFailure::NotConverged));
        assert_eq!(error.report.iterations[0].unconverged_coordinates, 1);
        assert_eq!((state, interface), (7, [0.0]));
    }
}
