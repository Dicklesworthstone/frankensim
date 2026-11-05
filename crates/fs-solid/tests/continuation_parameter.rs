//! Current parameter derivatives for follower loads and nonlinear load amplitudes.

use std::cell::Cell;

use fs_solid::SolidError;
use fs_solid::continuation::{ArcSettings, PathEvent, PathResidual, PathState, advance};
use fs_sparse::{Coo, Csr};

fn controls() -> ArcSettings {
    ArcSettings {
        ds: 0.08,
        ds_min: 0.08,
        ds_max: 0.08,
        max_corrector_iters: 4,
        tol: 2e-12,
    }
}

/// R = u - (lambda + alpha lambda^2)(1 + u^2).
/// Its exact equilibrium curve has a simple fold at u = 1.
struct Follower {
    alpha: f64,
}

impl Follower {
    fn amplitude(&self, lambda: f64) -> f64 {
        lambda + self.alpha * lambda * lambda
    }

    fn exact_load(&self, u: f64) -> f64 {
        let h = u / (1.0 + u * u);
        // Stable positive root, including alpha = 0.
        2.0 * h / (1.0 + (1.0 + 4.0 * self.alpha * h).sqrt())
    }
}

impl PathResidual for Follower {
    fn ndof(&self) -> usize {
        1
    }
    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        Ok(vec![u[0] - self.amplitude(lambda) * (1.0 + u[0] * u[0])])
    }
    fn tangent(&self, u: &[f64], lambda: f64) -> Result<Csr, SolidError> {
        let mut coo = Coo::new(1, 1);
        coo.push(0, 0, 1.0 - 2.0 * u[0] * self.amplitude(lambda));
        Ok(coo.assemble())
    }
    fn load_vector(&self) -> Vec<f64> {
        vec![1.0]
    }
    fn load_vector_at(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        Ok(vec![
            (1.0 + 2.0 * self.alpha * lambda) * (1.0 + u[0] * u[0]),
        ])
    }
}

#[test]
fn current_load_derivative_crosses_the_exact_follower_fold() {
    for alpha in [0.0, 0.75] {
        let problem = Follower { alpha };
        let settings = controls();
        let mut path = PathState::start(1, &settings);
        advance(&problem, &mut path, &settings, 44).expect("current bordered Newton derivative");
        assert!(path.u[0] > 2.0 && path.tangent.1 < 0.0);
        assert!(
            path.events
                .iter()
                .any(|e| matches!(e, PathEvent::LimitPoint { .. }))
        );
        assert!(
            !path
                .events
                .iter()
                .any(|e| matches!(e, PathEvent::BranchPoint { .. }))
        );
        for &(lambda, u) in &path.trace {
            assert!((lambda - problem.exact_load(u)).abs() < 2e-11);
        }
        let peak = path.trace.iter().map(|row| row.0).fold(0.0, f64::max);
        assert!((peak - problem.exact_load(1.0)).abs() < 2e-3);
        assert!(path.lambda < peak - 0.02);
    }
}

#[test]
fn follower_checkpoint_resume_matches_a_straight_run_bit_for_bit() {
    let problem = Follower { alpha: 0.75 };
    let settings = controls();
    let mut straight = PathState::start(1, &settings);
    let mut split = straight.clone();
    advance(&problem, &mut straight, &settings, 44).unwrap();
    advance(&problem, &mut split, &settings, 17).unwrap();
    advance(&problem, &mut split, &settings, 27).unwrap();
    assert_eq!(format!("{straight:?}"), format!("{split:?}"));
}

#[derive(Clone, Copy, Debug)]
enum Fault {
    Empty,
    Nonfinite,
    Refused,
    Cancelled,
    CancelledAfterReturn,
}

struct RefusingDerivative {
    fault: Fault,
    only_trial: bool,
    calls: Cell<usize>,
    cancelled: Cell<bool>,
}

impl PathResidual for RefusingDerivative {
    fn checkpoint(&self) -> Result<(), SolidError> {
        if self.cancelled.get() {
            Err(SolidError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn ndof(&self) -> usize {
        1
    }
    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        Follower { alpha: 0.0 }.residual(u, lambda)
    }
    fn tangent(&self, u: &[f64], lambda: f64) -> Result<Csr, SolidError> {
        Follower { alpha: 0.0 }.tangent(u, lambda)
    }
    fn load_vector(&self) -> Vec<f64> {
        vec![1.0]
    }
    fn load_vector_at(&self, u: &[f64], _lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.calls.set(self.calls.get() + 1);
        if self.only_trial && u[0] == 0.0 {
            return Ok(vec![1.0]);
        }
        match self.fault {
            Fault::Empty => Ok(Vec::new()),
            Fault::Nonfinite => Ok(vec![f64::NAN]),
            Fault::Refused => Err(SolidError::MaterialRefused {
                what: "parameter derivative refused the trial".into(),
            }),
            Fault::Cancelled => Err(SolidError::Cancelled),
            Fault::CancelledAfterReturn => {
                self.cancelled.set(true);
                Ok(vec![1.0])
            }
        }
    }
}

#[test]
fn parameter_derivative_validation_includes_zero_step_calls() {
    for fault in [Fault::Empty, Fault::Nonfinite] {
        for steps in [0, 1] {
            let problem = RefusingDerivative {
                fault,
                only_trial: false,
                calls: Cell::new(0),
                cancelled: Cell::new(false),
            };
            let settings = controls();
            let mut path = PathState::start(1, &settings);
            let before = format!("{path:?}");
            assert!(matches!(
                advance(&problem, &mut path, &settings, steps),
                Err(SolidError::InvalidInput { .. })
            ));
            assert_eq!(problem.calls.get(), 1);
            assert_eq!(format!("{path:?}"), before);
        }
    }
}

#[test]
fn trial_parameter_derivative_refusals_preserve_the_complete_checkpoint() {
    for fault in [
        Fault::Empty,
        Fault::Nonfinite,
        Fault::Refused,
        Fault::Cancelled,
        Fault::CancelledAfterReturn,
    ] {
        let problem = RefusingDerivative {
            fault,
            only_trial: true,
            calls: Cell::new(0),
            cancelled: Cell::new(false),
        };
        let mut settings = controls();
        if matches!(fault, Fault::Cancelled | Fault::CancelledAfterReturn) {
            // Cancellation must return immediately, even when smaller steps exist.
            settings.ds_min = 1e-5;
        }
        let mut path = PathState::start(1, &settings);
        let before = format!("{path:?}");
        let error = advance(&problem, &mut path, &settings, 1).unwrap_err();
        match fault {
            Fault::Empty | Fault::Nonfinite => {
                assert!(matches!(error, SolidError::InvalidInput { .. }));
            }
            Fault::Refused => assert_eq!(
                error,
                SolidError::MaterialRefused {
                    what: "parameter derivative refused the trial".into(),
                }
            ),
            Fault::Cancelled | Fault::CancelledAfterReturn => {
                assert_eq!(error, SolidError::Cancelled);
            }
        }
        assert_eq!(problem.calls.get(), 2, "{fault:?}");
        assert_eq!(format!("{path:?}"), before);
    }
}

/// A regular nonsymmetric Jacobian with both eigenvalues exactly one:
/// R = [u0 - lambda, u1 + u0^2]. Cholesky of its lower triangle alone
/// fails once 2*u0 > 1, despite the absence of any bifurcation.
struct RegularNonsymmetric;

impl PathResidual for RegularNonsymmetric {
    fn ndof(&self) -> usize {
        2
    }
    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        Ok(vec![u[0] - lambda, u[1] + u[0] * u[0]])
    }
    fn tangent(&self, u: &[f64], _lambda: f64) -> Result<Csr, SolidError> {
        let mut coo = Coo::new(2, 2);
        coo.push(0, 0, 1.0);
        coo.push(1, 0, 2.0 * u[0]);
        coo.push(1, 1, 1.0);
        Ok(coo.assemble())
    }
    fn load_vector(&self) -> Vec<f64> {
        vec![1.0, 0.0]
    }
}

#[test]
fn a_regular_nonsymmetric_tangent_does_not_report_a_branch_point() {
    let settings = controls();
    let mut path = PathState::start(2, &settings);
    advance(&RegularNonsymmetric, &mut path, &settings, 44).unwrap();
    assert!(path.u[0] > 1.0 && path.tangent.1 > 0.0);
    assert!((path.u[0] - path.lambda).abs() < settings.tol);
    assert!((path.u[1] + path.u[0] * path.u[0]).abs() < settings.tol);
    assert!(path.events.is_empty());
}
