//! Fixed native face-mass storage and the independent backward-Euler time loss.
use fs_math::dd::Dd;

/// Retain the actual represented endpoint difference for both physical storage
/// and time-discretization loss. Neither quantity is inferred from heat/work.
pub(super) struct FaceKineticIncrement {
    mass: Dd,
    difference: Dd,
    mean: Dd,
}

impl FaceKineticIncrement {
    pub(super) fn new(mass: f64, before: f64, after: f64) -> Self {
        Self {
            mass: Dd::from_f64(mass),
            difference: Dd::from_f64(after) - Dd::from_f64(before),
            mean: Dd::from_f64(0.5 * before) + Dd::from_f64(0.5 * after),
        }
    }

    pub(super) fn stored_j(&self) -> Dd {
        self.mass * self.difference * self.mean
    }

    pub(super) fn backward_euler_loss_j(&self) -> Dd {
        Dd::from_f64(0.5) * self.mass * self.difference * self.difference
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn endpoint_identity_retains_a_starting_speed_lost_by_native_subtraction() {
        let speed = 2.0_f64.powi(54);
        assert_eq!(speed - 1.0, speed);
        let change = FaceKineticIncrement::new(1.0, 1.0, speed);
        // Independent BE identity: u1*(u1-u0) = 2^108 - 2^54.
        let work = Dd::from_f64(2.0_f64.powi(108)) - Dd::from_f64(speed);
        let defect = change.stored_j() + change.backward_euler_loss_j() - work;
        assert!(defect.to_f64().abs() <= 1.0);
        let reverse = FaceKineticIncrement::new(1.0, speed, 1.0);
        assert_eq!(reverse.stored_j(), -change.stored_j());
        assert_eq!(
            reverse.backward_euler_loss_j(),
            change.backward_euler_loss_j()
        );
    }

    #[test]
    fn opposite_large_face_transfers_cancel_before_reduction_rounding() {
        let speed = 2.0_f64.powi(54);
        let positive = FaceKineticIncrement::new(3.0, speed, speed + 4.0);
        let negative = FaceKineticIncrement::new(3.0, speed, speed - 4.0);
        assert_eq!((positive.stored_j() + negative.stored_j()).to_f64(), 48.0);
        assert_eq!(
            (positive.backward_euler_loss_j() + negative.backward_euler_loss_j()).to_f64(),
            48.0,
        );
    }

    #[test]
    fn stationary_background_cannot_erase_a_separate_captured_face() {
        let background = FaceKineticIncrement::new(1.0, 2.0_f64.powi(54), 2.0_f64.powi(54));
        let captured = FaceKineticIncrement::new(3.0, 0.5, 0.0);
        assert_eq!(
            (background.stored_j() + captured.stored_j()).to_f64(),
            -3.0 / 8.0
        );
        assert_eq!(background.backward_euler_loss_j(), Dd::ZERO);
    }
}
