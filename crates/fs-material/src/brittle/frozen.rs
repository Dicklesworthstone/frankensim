//! Recoverable brittle stiffness at the original accepted damage state.
//!
//! The native cohesive point supplies its frozen unloading/compression
//! stiffness. Projection uses the same fixed modes and parallel elastic
//! remainder as every existing brittle trial. There is no second damage law,
//! isotropic damage multiplier or numerical residual stiffness.
//!
//! A frozen-history spectrum is conditional on fixed damage, coefficients and
//! the current opening/compression branch. It does not establish stability of
//! evolving fracture, a finite vibration amplitude, or contact between released
//! geometric surfaces. In particular a state on its loading envelope may grow
//! damage under an arbitrarily small outward increment.

use super::{BrittleCrackBand, BrittleError, BrittleState};
use crate::Tangent6;

impl BrittleCrackBand {
    /// Recoverable `d(stress)/d(strain)` in tensor-shear Voigt coordinates,
    /// holding the original accepted damage and onset mix fixed.
    ///
    /// This is the native storage Hessian. It is positive semidefinite under
    /// tensor-work weights, including exact zero transverse transfer after
    /// complete failure. In-plane elasticity and normal compression survive.
    /// The matrix is not necessarily symmetric in unweighted tensor-shear
    /// coordinates; `diag(1,1,1,2,2,2) * C` is symmetric.
    ///
    /// Owner/strain refusals remain identical to native evolution. At damaged
    /// zero opening the original cohesive owner returns
    /// `CohesiveError::NormalClosureSwitch`, since no unique two-sided tangent
    /// exists. Nothing in this operation advances history, work or time.
    pub fn frozen_tangent(&self, state: &BrittleState<'_>) -> Result<Tangent6, BrittleError> {
        self.jump_for_trial(state.strain(), state)?;
        self.tangent_from_modes(state.cohesive_point().frozen_tangent_pa_m()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brittle::{BrittleCrackBandParameters, BrittleCrackBandSettings};
    use crate::cohesive::CohesiveError;
    use crate::elastic::IsotropicElastic;
    use crate::tensor::{contract, rotate, rotation};
    use crate::{SmallStrainLaw, Voigt};

    fn law(normal: [f64; 3]) -> BrittleCrackBand {
        BrittleCrackBand::new(
            IsotropicElastic::new(3e9, 0.25, 0.05).unwrap(),
            BrittleCrackBandParameters {
                tensile_strength_pa: 3e6,
                shear_strength_pa: 6e6,
                mode_i_fracture_energy_j_m2: 120.0,
                mode_ii_fracture_energy_j_m2: 360.0,
                bk_exponent: 1.6,
            },
            BrittleCrackBandSettings {
                band_width_m: 0.004,
                normal,
                maximum_strain: 0.05,
                energy_tolerance_j_m3: 1e-6,
            },
        )
        .unwrap()
    }

    fn apply(c: &Tangent6, e: &Voigt) -> Voigt {
        c.map(|row| row.iter().zip(e).map(|(a, b)| a * b).sum())
    }

    fn near(a: f64, b: f64, floor: f64) {
        assert!(
            (a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(floor),
            "{a:e} != {b:e}"
        );
    }

    #[test]
    fn g1_frozen_band_is_native_unloading_storage_and_preserves_irreversible_history() {
        let law = law([1.0, 0.0, 0.0]);
        let loading = law
            .trial(&[0.012, 0.0, 0.0, 0.0, 0.0, 0.0], &law.initial_state())
            .unwrap();
        assert!(loading.branch_tangent_pa[0][0] < 0.0);
        let old = loading.state;
        let c = law.frozen_tangent(&old).unwrap();
        assert!(c[0][0] > 0.0);
        let target = [0.004, -0.0001, 0.0002, 0.0003, 0.0001, -0.0002];
        let unloaded = law.trial(&target, &old).unwrap();
        for (actual, expected) in unloaded.stress_pa.into_iter().zip(apply(&c, &target)) {
            near(actual, expected, 1e4);
        }
        near(
            unloaded.state.stored_energy_j_m3(),
            0.5 * contract(&target, &apply(&c, &target)),
            1.0,
        );
        assert_eq!(unloaded.state.damage(), old.damage());
        assert_eq!(unloaded.state.dissipated_j_m3(), old.dissipated_j_m3());
        assert_eq!(loading.state, old);
    }

    #[test]
    fn g0_failure_keeps_exact_transverse_null_modes_and_undamaged_compression() {
        let law = law([1.0, 0.0, 0.0]);
        let failed = law
            .trial(&[0.022, 0.0, 0.0, 0.0, 0.0, 0.0], &law.initial_state())
            .unwrap()
            .state;
        assert!(failed.failed());
        let c = law.frozen_tangent(&failed).unwrap();
        for index in [0, 3, 5] {
            let mut e = [0.0; 6];
            e[index] = 1.0;
            assert_eq!(apply(&c, &e), [0.0; 6]);
        }
        assert!(c[1][1] > 0.0 && c[4][4] > 0.0);
        let closed = law
            .trial(&[-0.001, 0.0, 0.0, 0.0, 0.0, 0.0], &failed)
            .unwrap()
            .state;
        let compression = law.frozen_tangent(&closed).unwrap();
        near(compression[0][0], law.elastic().stiffness()[0][0], 1.0);
        assert_eq!(compression[3][3], 0.0);
        assert_eq!(closed.dissipated_j_m3(), failed.dissipated_j_m3());
    }

    #[test]
    fn g0_native_owner_and_closure_refusals_do_not_reset_damage() {
        let law = law([1.0, 0.0, 0.0]);
        let c = law.frozen_tangent(&law.initial_state()).unwrap();
        for (row, expected) in c.iter().zip(law.elastic().stiffness()) {
            for (&a, b) in row.iter().zip(expected) {
                near(a, b, 3e9);
            }
        }
        let damaged = law
            .trial(&[0.012, 0.0, 0.0, 0.0, 0.0, 0.0], &law.initial_state())
            .unwrap()
            .state;
        let zero = law.trial(&[0.0; 6], &damaged).unwrap().state;
        assert!(matches!(
            law.frozen_tangent(&zero),
            Err(BrittleError::Cohesive(
                CohesiveError::NormalClosureSwitch { .. }
            ))
        ));
        let identical = law.clone();
        assert_eq!(
            identical.frozen_tangent(&damaged),
            Err(BrittleError::StateMismatch)
        );
        assert_eq!(zero.damage(), damaged.damage());
        assert_eq!(zero.dissipated_j_m3(), damaged.dissipated_j_m3());
    }

    #[test]
    fn g3_rotated_frozen_stiffness_preserves_tensor_work_symmetry_and_psd() {
        let original = law([1.0, 0.0, 0.0]);
        let q = rotation([0.3, -0.5, 0.8], 0.73);
        let rotated = law([q[0][0], q[1][0], q[2][0]]);
        let e = [0.012, 0.0, 0.0, 0.0, 0.0, 0.0];
        let a = original.trial(&e, &original.initial_state()).unwrap().state;
        let b = rotated
            .trial(&rotate(&e, &q), &rotated.initial_state())
            .unwrap()
            .state;
        let ca = original.frozen_tangent(&a).unwrap();
        let cb = rotated.frozen_tangent(&b).unwrap();
        let w = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];
        for i in 0..6 {
            for j in 0..6 {
                near(w[i] * cb[i][j], w[j] * cb[j][i], 3e9);
            }
            let mut direction = [0.0; 6];
            direction[i] = 0.001;
            let rotated_direction = rotate(&direction, &q);
            let expected = rotate(&apply(&ca, &direction), &q);
            let actual = apply(&cb, &rotated_direction);
            for k in 0..6 {
                near(actual[k], expected[k], 1e6);
            }
            assert!(contract(&rotated_direction, &actual) >= 0.0);
        }
    }
}
