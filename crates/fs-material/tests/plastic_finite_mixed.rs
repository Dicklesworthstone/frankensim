//! G0/G1/G3: the original multiplicative plastic return in a mixed potential.
//! Synthetic coefficients verify the numerical law, not a calibrated specimen.

use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const F: [f64; 9] = [1.08, 0.04, 0.01, -0.02, 0.96, 0.03, 0.01, 0.0, 1.02];

fn law(kappa: f64, hardening: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: kappa,
            yield_stress_pa: Y,
            hardening_modulus_pa: hardening,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.05,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 0.8,
        },
    )
    .unwrap()
}

fn laws() -> Vec<FiniteJ2> {
    vec![
        law(K, 4e6),
        law(K, 4e6).with_voce_hardening(6e6, 12.0).unwrap(),
        law(K, 0.0)
            .with_tabulated_hardening(
                TabulatedHardening::new(vec![
                    HardeningPoint {
                        plastic_strain: 0.0,
                        yield_stress_pa: Y,
                    },
                    HardeningPoint {
                        plastic_strain: 0.2,
                        yield_stress_pa: 5e6,
                    },
                    HardeningPoint {
                        plastic_strain: 0.8,
                        yield_stress_pa: 9e6,
                    },
                ])
                .unwrap(),
            )
            .unwrap(),
    ]
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= scale * tolerance,
        "{a:e} versus {b:e}, scale {scale:e}, tolerance {tolerance:e}"
    );
}

fn coaxial(axial: f64, dilation: f64) -> [f64; 9] {
    [
        (dilation + axial).exp(),
        0.0,
        0.0,
        0.0,
        (dilation - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (dilation - 0.5 * axial).exp(),
    ]
}

#[test]
fn g1_original_pressure_recovers_native_plastic_stress_energy_and_condensed_tangent() {
    for material in laws() {
        let old = material.initial_state();
        let native = material.response(&F, &old).unwrap();
        let remainder = material.remainder_response(&F, &old).unwrap();
        let j = determinant(&F);
        let pressure = K * j.ln() / j;
        let mixed = material
            .freeze_remainder(&old)
            .unwrap()
            .mixed_response(&F, pressure, VolumetricTreatment::OriginalMaterial)
            .unwrap();
        assert!(native.plastic_increment > 0.0);
        assert_eq!(remainder.state, native.state);
        assert_eq!(remainder.plastic_increment, native.plastic_increment);
        assert_eq!(remainder.flow_potential_j_m3, native.flow_potential_j_m3);
        assert_eq!(
            remainder.plastic_dissipation_j_m3,
            native.plastic_dissipation_j_m3
        );
        near(mixed.state.pressure_residual, 0.0, 1.0, 2e-14);
        near(
            mixed.state.potential_pa,
            native.incremental_potential_j_m3,
            MU,
            2e-12,
        );
        near(
            remainder.stored_energy_j_m3 + mixed.state.bulk_stored_pa,
            native.stored_energy_j_m3,
            MU,
            2e-12,
        );
        for row in 0..9 {
            near(
                mixed.state.first_piola_pa[row],
                native.first_piola_pa[row],
                MU,
                2e-12,
            );
            for col in 0..9 {
                // Independently eliminate the pressure unknown from the full
                // 2x2 block Hessian. The native penalty tangent must emerge.
                let condensed = mixed.tangent_pa[row][col]
                    - mixed.pressure_coupling[row] * mixed.pressure_coupling[col]
                        / mixed.pressure_tangent_per_pa;
                near(condensed, native.tangent_pa[row][col], K, 3e-11);
            }
        }
    }
}

#[test]
fn g0_mixed_energy_and_all_tangent_blocks_match_independent_differences() {
    for material in laws() {
        let old = material
            .response(&coaxial(0.04, 0.0), &material.initial_state())
            .unwrap()
            .state;
        let frozen = material.freeze_remainder(&old).unwrap();
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            let pressure = 3e6;
            let response = frozen.mixed_response(&F, pressure, treatment).unwrap();
            for col in 0..9 {
                let h = 2e-6;
                let (mut plus, mut minus) = (F, F);
                plus[col] += h;
                minus[col] -= h;
                let plus = frozen.mixed_state(&plus, pressure, treatment).unwrap();
                let minus = frozen.mixed_state(&minus, pressure, treatment).unwrap();
                near(
                    (plus.potential_pa - minus.potential_pa) / (2.0 * h),
                    response.state.first_piola_pa[col],
                    MU,
                    2e-8,
                );
                near(
                    (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
                    response.pressure_coupling[col],
                    1.0,
                    2e-9,
                );
                for row in 0..9 {
                    near(
                        (plus.first_piola_pa[row] - minus.first_piola_pa[row]) / (2.0 * h),
                        response.tangent_pa[row][col],
                        MU,
                        2e-7,
                    );
                    near(
                        response.tangent_pa[row][col],
                        response.tangent_pa[col][row],
                        MU,
                        2e-10,
                    );
                }
            }
            let h = 10.0;
            let plus = frozen.mixed_state(&F, pressure + h, treatment).unwrap();
            let minus = frozen.mixed_state(&F, pressure - h, treatment).unwrap();
            near(
                (plus.potential_pa - minus.potential_pa) / (2.0 * h),
                response.state.pressure_residual,
                1.0,
                2e-9,
            );
            near(
                (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
                response.pressure_tangent_per_pa,
                1.0 / K,
                2e-8,
            );
            for row in 0..9 {
                near(
                    (plus.first_piola_pa[row] - minus.first_piola_pa[row]) / (2.0 * h),
                    response.pressure_coupling[row],
                    1.0,
                    2e-9,
                );
            }
        }
    }
}

#[test]
fn g3_loading_unloading_and_reverse_flow_preserve_original_history_and_yield_heat() {
    for material in laws() {
        let mut old = material.initial_state();
        let mut accumulated = 0.0;
        for (axial, dilation) in [
            (0.03, 0.01),
            (0.08, 0.015),
            (0.075, 0.0),
            (-0.05, -0.01),
            (0.03, 0.0),
        ] {
            let f = coaxial(axial, dilation);
            let before = old.clone();
            let full = material.response(&f, &old).unwrap();
            let rem = material.remainder_response(&f, &old).unwrap();
            assert_eq!(old, before);
            assert_eq!(rem.state, full.state);
            assert_eq!(rem.plastic_dissipation_j_m3, full.plastic_dissipation_j_m3);
            near(
                rem.plastic_dissipation_j_m3,
                Y * (rem.state.alpha - old.alpha),
                MU,
                1e-14,
            );
            material
                .validate_remainder_committed(&f, &rem.state, 1e-4)
                .unwrap();
            near(
                material.stored_remainder_energy(&f, &rem.state).unwrap(),
                rem.stored_energy_j_m3,
                MU,
                2e-12,
            );
            assert!(rem.state.alpha >= accumulated);
            accumulated = rem.state.alpha;
            old = rem.state;
        }
        assert!(accumulated > 0.1);
    }
}

#[test]
fn g3_remainder_is_bulk_independent_and_exact_constraint_never_evaluates_penalty() {
    let moderate = law(K, 4e6);
    let enormous = law(1e308, 4e6);
    let old = moderate.initial_state();
    for f in [F, [0.7, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.7]] {
        let expected = moderate.remainder_response(&f, &old).unwrap();
        let actual = enormous.remainder_response(&f, &old).unwrap();
        assert_eq!(actual, expected);
        let a = moderate
            .freeze_remainder(&old)
            .unwrap()
            .mixed_response(&f, 3e6, VolumetricTreatment::ExactIncompressible)
            .unwrap();
        let b = enormous
            .freeze_remainder(&old)
            .unwrap()
            .mixed_response(&f, 3e6, VolumetricTreatment::ExactIncompressible)
            .unwrap();
        assert_eq!(a.state.first_piola_pa, b.state.first_piola_pa);
        assert_eq!(a.tangent_pa, b.tangent_pa);
        assert_eq!(b.pressure_tangent_per_pa, 0.0);
        assert_eq!(b.state.bulk_stored_pa, 0.0);
        enormous
            .validate_remainder_committed(&f, &actual.state, 1e-4)
            .unwrap();
    }
    // The full penalty Hessian here is unrepresentable. The exact mixed
    // formulation is still usable because it never computes that Hessian.
    assert!(
        enormous
            .response(&[0.7, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.7], &old)
            .is_err()
    );
}

#[test]
fn g0_native_domain_and_original_pressure_branch_refuse_without_history_changes() {
    let material = law(K, 4e6);
    let old = material.initial_state();
    let frozen = material.freeze_remainder(&old).unwrap();
    assert!(
        frozen
            .mixed_response(&I, K, VolumetricTreatment::OriginalMaterial)
            .is_err()
    );
    assert!(
        frozen
            .mixed_response(&I, -100.0 * K, VolumetricTreatment::OriginalMaterial)
            .is_err()
    );
    assert!(
        frozen
            .mixed_response(&I, f64::NAN, VolumetricTreatment::ExactIncompressible)
            .is_err()
    );
    let mut inverted = I;
    inverted[0] = -1.0;
    assert!(material.remainder_response(&inverted, &old).is_err());
    let mut invalid = old.clone();
    invalid.inverse_plastic_cauchy_green[0] = 2.0;
    assert!(material.freeze_remainder(&invalid).is_err());
    assert_eq!(old, material.initial_state());
    assert!(material.remainder_response(&F, &old).is_ok());
}
