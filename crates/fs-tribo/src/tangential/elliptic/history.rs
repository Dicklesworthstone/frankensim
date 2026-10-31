//! Native history acceptance after an independently admitted mechanical path.
//!
//! This operation records already resolved slip in the checkpoint's current
//! tangent frame. It performs no contact, friction, compliance or time solve.
//! No normal load, traction, impulse, pressure, work or heat is inferred from
//! the displacement. Mechanical admissibility and dissipation remain with the
//! owner that resolved the motion. Objective transport is a separate operation.

use super::{EllipticCoulomb, EllipticCoulombState, TriboError, finite, invalid};
use crate::{nonnegative_finite, positive_finite};

/// A private candidate history and the actual representable vector increment.
/// There is deliberately no mechanical force or dissipation receipt here.
#[derive(Clone, Debug, PartialEq)]
pub struct EllipticSlipHistoryReport {
    pub state: EllipticCoulombState,
    pub slip_increment_m: [f64; 2],
    pub represented_slip_increment_m: [f64; 2],
    /// Represented minus supplied increment, in the current tangent coordinates.
    pub slip_residual_m: [f64; 2],
}

impl EllipticCoulomb {
    /// Record actual resolved slip and an explicitly supplied thermal boundary.
    ///
    /// The mechanical owner must have admitted the complete displacement before
    /// publishing this candidate. Coefficients and texture do not change. This
    /// zero-storage law has no extra slip-dependent internal force to evolve;
    /// it does not authorize the same operation for another constitutive model.
    /// The state carries no clock; its time interval belongs to the caller.
    ///
    /// `absolute_tolerance_m` independently bounds each component's addition
    /// roundoff. A nonzero component that disappears entirely refuses even if
    /// the requested tolerance would otherwise allow it. All failures preserve
    /// the accepted checkpoint.
    pub fn accept_resolved_slip(
        &self,
        state: &EllipticCoulombState,
        slip_increment_m: [f64; 2],
        contact_temperature_k: Option<f64>,
        absolute_tolerance_m: f64,
    ) -> Result<EllipticSlipHistoryReport, TriboError> {
        self.check(state)?;
        nonnegative_finite(absolute_tolerance_m, "elliptic_slip_history_tolerance_m")?;
        if let Some(temperature) = contact_temperature_k {
            positive_finite(temperature, "contact_temperature_k")?;
        }
        let mut next = state.clone();
        let mut represented_slip_increment_m = [0.0; 2];
        let mut slip_residual_m = [0.0; 2];
        for axis in 0..2 {
            let increment = slip_increment_m[axis];
            finite(increment, "elliptic_resolved_slip_increment_m")?;
            let value = state.slip_m[axis] + increment;
            finite(value, "elliptic_accumulated_slip_m")?;
            let represented = value - state.slip_m[axis];
            let residual = represented - increment;
            finite(represented, "elliptic_represented_slip_increment_m")?;
            finite(residual, "elliptic_slip_history_residual_m")?;
            if (increment != 0.0 && represented == 0.0) || residual.abs() > absolute_tolerance_m {
                return Err(invalid("elliptic_resolved_slip_history_resolution"));
            }
            next.slip_m[axis] = value;
            represented_slip_increment_m[axis] = represented;
            slip_residual_m[axis] = residual;
        }
        next.contact_temperature_k = contact_temperature_k;
        Ok(EllipticSlipHistoryReport {
            state: next,
            slip_increment_m,
            represented_slip_increment_m,
            slip_residual_m,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TangentialLaw;
    use crate::partial_slip::TangentFrame;
    use crate::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};

    fn law() -> EllipticCoulomb {
        EllipticCoulomb::new([0.8, 0.4], [0.6, 0.2], [1.0, 0.0]).unwrap()
    }

    #[test]
    fn accepted_vector_history_reverses_without_a_second_mechanical_solve() {
        let law = law();
        let initial = law.initial_state();
        let first = law
            .accept_resolved_slip(&initial, [0.25, -0.5], Some(300.0), 0.0)
            .unwrap();
        let reverse = law
            .accept_resolved_slip(&first.state, [-0.25, 0.5], Some(310.0), 0.0)
            .unwrap();
        assert_eq!(first.state.slip_m(), [0.25, -0.5]);
        assert_eq!(first.represented_slip_increment_m, first.slip_increment_m);
        assert_eq!(first.slip_residual_m, [0.0; 2]);
        assert_eq!(reverse.state.slip_m(), [0.0; 2]);
        assert_eq!(reverse.state.principal_axis(), initial.principal_axis());
        assert_eq!(reverse.state.contact_temperature_k(), Some(310.0));
        assert_eq!(law.stored_energy_j(&reverse.state), 0.0);
        assert_eq!(initial, law.initial_state());
        assert_eq!(
            law.accept_resolved_slip(&initial, [0.0; 2], None, 0.0)
                .unwrap()
                .state,
            initial,
        );
    }

    #[test]
    fn accepted_slip_uses_the_existing_objective_history_isometry() {
        let law = law();
        let initial = law.initial_state();
        let old_frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
        let new_frame = TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap();
        let rotation = TangentialHistoryRotation::between(old_frame, new_frame, 0.0).unwrap();
        let slip = [0.25, -0.5];
        let accepted = law
            .accept_resolved_slip(&initial, slip, Some(300.0), 0.0)
            .unwrap();
        let transported = law.transport_state(&accepted.state, rotation).unwrap();
        let start = law.transport_state(&initial, rotation).unwrap();
        let in_new_frame = law
            .accept_resolved_slip(&start, rotation.apply(slip).unwrap(), Some(300.0), 0.0)
            .unwrap();
        assert_eq!(transported, in_new_frame.state);
        assert_eq!(transported.slip_m(), [-0.5, -0.25]);
        assert_eq!(transported.principal_axis(), [0.0, -1.0]);
    }

    #[test]
    fn foreign_nonfinite_and_unrepresentable_updates_preserve_accepted_history() {
        let law = law();
        let initial = law.initial_state();
        let foreign = EllipticCoulomb::new([0.8, 0.4], [0.5, 0.2], [1.0, 0.0]).unwrap();
        assert!(
            foreign
                .accept_resolved_slip(&initial, [0.0; 2], None, 0.0)
                .is_err()
        );
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(
                law.accept_resolved_slip(&initial, [bad, 0.0], None, 0.0)
                    .is_err()
            );
            assert!(
                law.accept_resolved_slip(&initial, [0.0; 2], Some(bad), 0.0)
                    .is_err()
            );
            assert!(
                law.accept_resolved_slip(&initial, [0.0; 2], None, bad)
                    .is_err()
            );
        }
        assert!(
            law.accept_resolved_slip(&initial, [0.0; 2], Some(0.0), 0.0)
                .is_err()
        );
        assert!(
            law.accept_resolved_slip(&initial, [0.0; 2], None, -1.0)
                .is_err()
        );
        let large = law
            .accept_resolved_slip(&initial, [1e16, 0.0], None, 0.0)
            .unwrap()
            .state;
        let saved = large.clone();
        assert!(
            law.accept_resolved_slip(&large, [0.5, 0.0], None, 100.0)
                .is_err()
        );
        assert!(
            law.accept_resolved_slip(&large, [1.5, 0.0], None, 0.25)
                .is_err()
        );
        assert_eq!(large, saved);
        assert_eq!(initial, law.initial_state());
    }
}
