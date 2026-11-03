use super::*;
use fs_solver::{CheckedCgConfig, LinearOp};
use fs_sparse::precond::IdentityPrecond;
use std::cell::Cell;

/// An explicit SPD test mass with a nonzero off-diagonal coupling. Its 2x2
/// inverse and endpoint kinetic energies are independent analytical oracles.
struct Mass {
    entries: [[f64; 2]; 2],
    applications: Cell<usize>,
}

impl Mass {
    fn coupled() -> Self {
        Self {
            entries: [[2.0, 0.5], [0.5, 1.0]],
            applications: Cell::new(0),
        }
    }

    fn identity() -> Self {
        Self {
            entries: [[1.0, 0.0], [0.0, 1.0]],
            applications: Cell::new(0),
        }
    }
}

impl LinearOp for Mass {
    fn n(&self) -> usize {
        2
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.applications.set(self.applications.get() + 1);
        y[0] = self.entries[0][0] * x[0] + self.entries[0][1] * x[1];
        y[1] = self.entries[1][0] * x[0] + self.entries[1][1] * x[1];
    }
}

fn controls() -> ImpulseControls {
    ImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 16,
            max_corrections: 2,
        },
        max_workspace_entries: 512,
        absolute_momentum_tolerance_ns: 1e-13,
        relative_momentum_tolerance: 1e-13,
        absolute_work_tolerance_j: 1e-13,
        relative_work_tolerance: 1e-13,
    }
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 3e-14, "{actual} != {expected}");
}

#[test]
fn coupled_original_mass_preserves_signed_impulse_and_reversible_work() {
    let mass = Mass::coupled();
    let before = [1.0, -2.0];
    let impulse = [-0.3, 0.8];
    let (after, report) = velocity_candidate(
        &mass,
        &before,
        &impulse,
        &IdentityPrecond,
        10.0,
        controls(),
        &mut || false,
    )
    .unwrap();

    // det(M)=7/4, so M^-1 J=(-2/5,1), independently of checked CG.
    close(after[0], 0.6);
    close(after[1], -1.0);
    close(
        2.0 * (after[0] - before[0]) + 0.5 * (after[1] - before[1]),
        impulse[0],
    );
    close(
        0.5 * (after[0] - before[0]) + (after[1] - before[1]),
        impulse[1],
    );
    assert!(report.solve.converged_euclidean());
    close(report.kinetic_energy_before_j, 2.0);
    close(report.kinetic_energy_after_j, 0.56);
    close(report.kinetic_energy_change_j, -1.44);
    close(report.applied_impulse_work_j, -1.44);
    close(report.work_residual_j, 0.0);
    assert!(report.max_momentum_residual_ns < 1e-14);

    let (restored, reverse) = velocity_candidate(
        &mass,
        &after,
        &[0.3, -0.8],
        &IdentityPrecond,
        10.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(restored[0], before[0]);
    close(restored[1], before[1]);
    close(reverse.applied_impulse_work_j, 1.44);
    close(
        reverse.kinetic_energy_change_j + report.kinetic_energy_change_j,
        0.0,
    );

    let (unchanged, zero) = velocity_candidate(
        &mass,
        &before,
        &[0.0, 0.0],
        &IdentityPrecond,
        10.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(unchanged, before);
    assert_eq!(zero.solve.iters, 0);
    assert_eq!(zero.kinetic_energy_change_j, 0.0);
    assert_eq!(zero.applied_impulse_work_j, 0.0);
}

#[test]
fn momentum_gate_checks_the_representable_velocity_after_addition() {
    let mass = Mass::identity();
    let before = [1e20, 0.0];
    let error = velocity_candidate(
        &mass,
        &before,
        &[1.0, 0.0],
        &IdentityPrecond,
        1e21,
        controls(),
        &mut || false,
    )
    .unwrap_err();
    // CG solves delta_v=1 exactly, but adding it to 1e20 loses the increment.
    assert!(matches!(
        error,
        ImpulseError::Momentum {
            component: 0,
            residual_ns: -1.0,
            ..
        }
    ));
    assert_eq!(before, [1e20, 0.0]);
}

#[test]
fn work_gate_uses_event_scale_even_with_large_background_kinetic_energy() {
    let mass = Mass::identity();
    let mut loose_momentum = controls();
    loose_momentum.absolute_momentum_tolerance_ns = 2.0;
    loose_momentum.absolute_work_tolerance_j = 0.5;
    let error = velocity_candidate(
        &mass,
        &[1e16, 0.0],
        &[1.0, 0.0],
        &IdentityPrecond,
        1e17,
        loose_momentum,
        &mut || false,
    )
    .unwrap_err();
    // Deliberately allow the lost unit impulse through its momentum gate.
    // Its 1e16 J work defect must not disappear against the 5e31 J background.
    match error {
        ImpulseError::Work {
            residual_j,
            tolerance_j,
        } => {
            assert_eq!(residual_j, -1e16);
            assert!(tolerance_j < 2e3);
        }
        other => panic!("expected independent work refusal, received {other:?}"),
    }
}

#[test]
fn budget_cancellation_solver_and_speed_refusals_publish_no_velocity() {
    let mass = Mass::coupled();
    let before = [1.0, -2.0];
    let impulse = [-0.3, 0.8];
    let mut bounded = controls();
    bounded.max_workspace_entries = 0;
    assert!(matches!(
        velocity_candidate(
            &mass,
            &before,
            &impulse,
            &IdentityPrecond,
            10.0,
            bounded,
            &mut || false,
        ),
        Err(ImpulseError::Budget(_))
    ));
    assert_eq!(mass.applications.get(), 0);

    bounded = controls();
    bounded.solve.max_iterations = usize::MAX;
    assert!(matches!(
        velocity_candidate(
            &mass,
            &before,
            &impulse,
            &IdentityPrecond,
            10.0,
            bounded,
            &mut || false,
        ),
        Err(ImpulseError::Budget(_))
    ));
    assert_eq!(mass.applications.get(), 0);

    assert!(matches!(
        velocity_candidate(
            &mass,
            &before,
            &impulse,
            &IdentityPrecond,
            10.0,
            controls(),
            &mut || true,
        ),
        Err(ImpulseError::Cancelled)
    ));
    assert_eq!(mass.applications.get(), 0);
    assert!(matches!(
        velocity_candidate(
            &mass,
            &before,
            &impulse,
            &IdentityPrecond,
            10.0,
            controls(),
            &mut || mass.applications.get() != 0,
        ),
        Err(ImpulseError::Cancelled)
    ));
    assert!(mass.applications.get() > 0);

    let mut too_few = controls();
    too_few.solve.max_iterations = 1;
    match velocity_candidate(
        &mass,
        &before,
        &impulse,
        &IdentityPrecond,
        10.0,
        too_few,
        &mut || false,
    ) {
        Err(ImpulseError::Solve(report)) => {
            assert_eq!(report.iters, 1);
            assert!(report.euclidean_rel_residual().is_some());
            assert!(!report.converged_euclidean());
        }
        other => panic!("expected solver budget refusal, received {other:?}"),
    }
    assert!(matches!(
        velocity_candidate(
            &mass,
            &before,
            &[30.0, 0.0],
            &IdentityPrecond,
            10.0,
            controls(),
            &mut || false,
        ),
        Err(ImpulseError::Velocity { .. })
    ));
    assert_eq!(before, [1.0, -2.0]);
}
