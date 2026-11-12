//! Original cohesive mass-owner admission and compensated exchange regressions.

use super::kinetic_change;
use fs_solver::LinearOp;

#[path = "kinetic_tests.rs"]
mod tests;

struct UncalledMass;
impl LinearOp for UncalledMass {
    fn n(&self) -> usize {
        2
    }

    fn apply(&self, _: &[f64], _: &mut [f64]) {
        panic!("refused state reached the mass action");
    }

    fn apply_transpose(&self, _: &[f64], _: &mut [f64]) {
        panic!("refused state reached the mass action");
    }
}

#[test]
fn incomplete_original_coordinates_refuse_before_mass_action() {
    assert!(matches!(
        kinetic_change(&UncalledMass, &[1.0], &[2.0], &mut || false),
        Err(super::Error::Invalid(_))
    ));
}

#[test]
fn native_cancellation_refuses_before_mass_action() {
    assert!(matches!(
        kinetic_change(&UncalledMass, &[1.0, 2.0], &[2.0, 3.0], &mut || true),
        Err(super::Error::Cancelled)
    ));
}
