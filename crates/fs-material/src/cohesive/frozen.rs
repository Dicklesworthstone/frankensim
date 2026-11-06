//! Recoverable cohesive stiffness at an immutable accepted history.
//!
//! Differentiating stored energy while HOLDING damage and the onset mode mix
//! fixed gives the unloading stiffness. It omits the derivative of damage on
//! the loading envelope, and therefore is not that envelope's algorithmic
//! tangent. This is the stiffness required by a declared frozen-history modal
//! approximation; it neither advances nor suppresses fracture in a native
//! evolution step. A finite oscillation still needs to remain within the same
//! opening/compression branch and avoid additional damage.
//!
//! All coefficients and remaining integrity come from the original cohesive
//! owner. Complete failure has exact zero opening/shear transfer; no residual
//! stiffness is added. Normal compression retains its original stiffness.
//! Physical bond dashpots, if present, are not part of a storage Hessian.

use super::{CohesiveError, CohesivePoint, finite};

impl CohesivePoint<'_> {
    /// Hessian of recoverable energy with irreversible history held fixed,
    /// in the original local `[normal, shear 1, shear 2]` frame, Pa/m.
    ///
    /// At zero normal jump, unequal opening/compression stiffnesses produce a
    /// typed refusal instead of choosing a branch silently. Equal stiffnesses
    /// have a unique tangent, including a virgin brittle band's zero state.
    /// A state on its damage envelope still returns the FROZEN unloading
    /// tangent, not a claim of a classical derivative of irreversible loading.
    /// This fixed-size observation allocates nothing and changes no history.
    pub fn frozen_tangent_pa_m(&self) -> Result<[[f64; 3]; 3], CohesiveError> {
        self.require_jump(self.jump, self.law.maximum_separation)?;
        // Read native integrity directly: reconstructing it as 1 - damage
        // can erase a small but representable surviving ligament near failure.
        let integrity = if self.law.mixed_mode.is_some() {
            self.mixed_integrity()
        } else {
            self.law.intact(self.maximum)
        };
        let (normal, shear) = self
            .law
            .mixed_mode
            .map_or((self.law.stiffness, self.law.stiffness), |p| {
                (p.normal_stiffness_pa_m, p.tangential_stiffness_pa_m)
            });
        let opening = finite(integrity * normal)?;
        let shear = finite(integrity * shear)?;
        let compression = self.law.compression_stiffness;
        if self.jump[0] == 0.0 && opening != compression {
            return Err(CohesiveError::NormalClosureSwitch {
                opening_stiffness_pa_m: opening,
                compression_stiffness_pa_m: compression,
            });
        }
        Ok([
            [
                if self.jump[0] < 0.0 {
                    compression
                } else {
                    opening
                },
                0.0,
                0.0,
            ],
            [0.0, shear, 0.0],
            [0.0, 0.0, shear],
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohesive::BilinearCohesive;
    use crate::cohesive::mixed_mode::MixedModeParameters;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g1_frozen_unloading_stiffness_differs_from_the_negative_loading_tangent() {
        let law = BilinearCohesive::new(1e9, 1e6, 2000.0, 2e9, 0.02).unwrap();
        let loading = law.bonded_point().trial([0.002, 0.0, 0.0], 1e-9).unwrap();
        assert!(loading.branch_tangent_pa_m[0][0] < 0.0);
        let before = loading.point;
        let frozen = loading.point.frozen_tangent_pa_m().unwrap();
        near(frozen[0][0], 1e9 / 3.0);
        assert_eq!(loading.point, before);
        let unloaded = loading.point.trial([0.0015, 0.0, 0.0], 1e-9).unwrap();
        near(unloaded.traction_pa[0], frozen[0][0] * 0.0015);
        near(
            unloaded.point.stored_energy_j_m2(),
            0.5 * frozen[0][0] * 0.0015_f64.powi(2),
        );
        assert_eq!(unloaded.point.damage(), before.damage());
        assert_eq!(unloaded.point.dissipated_j_m2(), before.dissipated_j_m2());
    }

    #[test]
    fn g0_closure_is_ambiguous_only_when_the_two_normal_stiffnesses_differ() {
        let equal = BilinearCohesive::new(1e9, 1e6, 2000.0, 1e9, 0.02).unwrap();
        assert_eq!(
            equal.bonded_point().frozen_tangent_pa_m().unwrap()[0][0],
            1e9
        );
        let law = BilinearCohesive::new(1e9, 1e6, 2000.0, 2e9, 0.02).unwrap();
        assert!(matches!(
            law.bonded_point().frozen_tangent_pa_m(),
            Err(CohesiveError::NormalClosureSwitch { .. })
        ));
        let failed = law
            .bonded_point()
            .trial([0.006, 0.0, 0.0], 1e-9)
            .unwrap()
            .point;
        assert_eq!(failed.frozen_tangent_pa_m().unwrap(), [[0.0; 3]; 3]);
        let closed = failed.trial([-0.0001, 0.0, 0.0], 1e-9).unwrap().point;
        assert_eq!(
            closed.frozen_tangent_pa_m().unwrap(),
            [[2e9, 0.0, 0.0], [0.0; 3], [0.0; 3]]
        );
        let switch = closed.trial([0.0; 3], 1e-9).unwrap().point;
        assert!(matches!(
            switch.frozen_tangent_pa_m(),
            Err(CohesiveError::NormalClosureSwitch { .. })
        ));
        assert!(switch.failed());
        assert_eq!(switch.dissipated_j_m2(), failed.dissipated_j_m2());
    }

    #[test]
    fn g0_representable_remaining_integrity_is_not_rounded_away_by_damage() {
        let law = BilinearCohesive::new(1.0, 1e-8, 5e-9, 1.0, 2.0).unwrap();
        let jump = law.failure_separation_m().next_down();
        let point = law
            .bonded_point()
            .trial([jump, 0.0, 0.0], 1e-20)
            .unwrap()
            .point;
        assert!(!point.failed());
        assert_eq!(
            point.damage(),
            1.0,
            "1 - integrity rounds to one near failure"
        );
        let stiffness = point.frozen_tangent_pa_m().unwrap();
        assert!(stiffness[0][0] > 0.0);
        assert_eq!(stiffness[0][0] * jump, point.traction_pa()[0]);
    }

    #[test]
    fn g1_mixed_mode_stiffness_retains_native_anisotropy_and_quadratic_storage() {
        let parameters = MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 4e9,
            normal_strength_pa: 1e6,
            tangential_strength_pa: 2e6,
            mode_i_fracture_energy_j_m2: 2000.0,
            mode_ii_fracture_energy_j_m2: 6000.0,
            bk_exponent: 1.6,
        };
        let law = BilinearCohesive::new_mixed_mode(parameters, 1e9, 0.02).unwrap();
        let point = law
            .bonded_point()
            .trial([0.002, 0.001, 0.0003], 1e-8)
            .unwrap()
            .point;
        assert!(point.damage() > 0.0 && !point.failed());
        let stiffness = point.frozen_tangent_pa_m().unwrap();
        near(stiffness[1][1], 4.0 * stiffness[0][0]);
        let jump = point.jump_m().map(|x| 0.6 * x);
        let unloaded = point.trial(jump, 1e-8).unwrap().point;
        let energy: f64 = (0..3)
            .map(|i| 0.5 * stiffness[i][i] * jump[i] * jump[i])
            .sum();
        near(unloaded.stored_energy_j_m2(), energy);
        for (i, traction) in unloaded.traction_pa().into_iter().enumerate() {
            near(traction, stiffness[i][i] * jump[i]);
        }
        assert_eq!(unloaded.mixed_mode_history(), point.mixed_mode_history());
        assert_eq!(unloaded.dissipated_j_m2(), point.dissipated_j_m2());
    }
}
