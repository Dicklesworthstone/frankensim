//! Original-energy pressure duals and finite-strain mixed constitutive blocks.

use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::mixed::{
    MixedFiniteStrainPotential, VolumetricTreatment, pressure_potential,
};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const F: [f64; 9] = [1.16, 0.07, 0.02, -0.01, 0.97, 0.03, 0.01, 0.0, 0.93];
const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} vs {b:e}; scale {scale:e}, tolerance {tolerance:e}"
    );
}
fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}
fn models() -> Vec<Hyperelastic> {
    vec![
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa: 8e5,
            },
            4.0,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![
                    OgdenTerm {
                        mu_pa: 2e5,
                        alpha: 1.3,
                    },
                    OgdenTerm {
                        mu_pa: -1e4,
                        alpha: -2.0,
                    },
                ],
                volumetric: Volumetric::LogSquared { kappa: 8e5 },
            },
            4.0,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 5e4,
                k1: 1.2e5,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0.0, 0.6], [0.0, 0.6, 0.8]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 8e5 },
            },
            4.0,
        )
        .unwrap(),
    ]
}
fn visco() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        models().remove(0),
        vec![
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3e4 },
                tau_s: 0.2,
            },
            FiniteStrainBranch {
                energy: BranchEnergy::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 2e4,
                        alpha: 1.7,
                    }],
                },
                tau_s: 0.07,
            },
        ],
    )
    .unwrap()
}
fn pressure(law: &Hyperelastic, j: f64) -> f64 {
    match &law.model {
        HyperelasticModel::MooneyRivlin { kappa, .. } => kappa * (j - 1.0),
        HyperelasticModel::Ogden { volumetric, .. }
        | HyperelasticModel::HolzapfelGasserOgden { volumetric, .. } => match volumetric {
            Volumetric::Quadratic { kappa } => kappa * (j - 1.0),
            Volumetric::LogSquared { kappa } => kappa * j.ln() / j,
        },
        HyperelasticModel::NeoHookean { .. } => panic!("coupled fixture not admitted"),
    }
}

#[test]
fn g0_original_pressure_duals_recover_volume_storage_and_independent_derivatives() {
    let kappa = 8e5;
    for volume in [
        Volumetric::Quadratic { kappa },
        Volumetric::LogSquared { kappa },
    ] {
        for j in [1e-6_f64, 0.2, 0.8, 1.0, 1.000000001, 1.2, 2.5] {
            let (p, u, curvature) = match volume {
                Volumetric::Quadratic { .. } => {
                    (kappa * (j - 1.0), 0.5 * kappa * (j - 1.0).powi(2), kappa)
                }
                Volumetric::LogSquared { .. } => (
                    kappa * j.ln() / j,
                    0.5 * kappa * j.ln().powi(2),
                    kappa * (1.0 - j.ln()) / (j * j),
                ),
            };
            let result = pressure_potential(volume, p).unwrap();
            near(result.volume_ratio, j, j, 2e-10);
            near(result.bulk_stored_pa, u, kappa.max(u.abs()), 3e-12);
            near(
                result.complementary_pa,
                p * (j - 1.0) - u,
                kappa.max(p.abs()),
                3e-12,
            );
            near(
                result.compliance_per_pa,
                1.0 / curvature,
                1.0 / curvature,
                2e-10,
            );
            let h = 1e-7 * kappa.max(p.abs());
            let plus = pressure_potential(volume, p + h).unwrap();
            let minus = pressure_potential(volume, p - h).unwrap();
            near(
                (plus.complementary_pa - minus.complementary_pa) / (2.0 * h),
                j - 1.0,
                j.max(1.0),
                2e-7,
            );
            near(
                (plus.volume_ratio - minus.volume_ratio) / (2.0 * h),
                result.compliance_per_pa,
                result.compliance_per_pa,
                2e-5,
            );
        }
    }
    // The log inverse remains on the lower, stable branch close to its fold.
    let j = std::f64::consts::E * (1.0 - 1e-5);
    let result = pressure_potential(Volumetric::LogSquared { kappa: 1.0 }, j.ln() / j).unwrap();
    near(result.volume_ratio, j, j, 5e-10);
    assert!(result.volume_ratio < std::f64::consts::E);
}

fn derivative_gate(law: &impl MixedFiniteStrainPotential, treatment: VolumetricTreatment) {
    law.validate_mixed().unwrap();
    let p = 1.1e5;
    let full = law.mixed_response(&F, p, treatment).unwrap();
    let state = law.mixed_state(&F, p, treatment).unwrap();
    near(state.potential_pa, full.state.potential_pa, 1e6, 1e-14);
    for i in 0..9 {
        near(
            state.first_piola_pa[i],
            full.state.first_piola_pa[i],
            1e6,
            1e-13,
        );
    }
    for col in 0..9 {
        let h = 2e-6;
        let mut plus = F;
        let mut minus = F;
        plus[col] += h;
        minus[col] -= h;
        let plus = law.mixed_state(&plus, p, treatment).unwrap();
        let minus = law.mixed_state(&minus, p, treatment).unwrap();
        near(
            (plus.potential_pa - minus.potential_pa) / (2.0 * h),
            state.first_piola_pa[col],
            1e6,
            2e-7,
        );
        near(
            (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
            full.pressure_coupling[col],
            1.0,
            2e-8,
        );
        for row in 0..9 {
            near(
                (plus.first_piola_pa[row] - minus.first_piola_pa[row]) / (2.0 * h),
                full.tangent_pa[row][col],
                1e6,
                2e-6,
            );
            near(
                full.tangent_pa[row][col],
                full.tangent_pa[col][row],
                1e6,
                2e-12,
            );
        }
    }
    let h = 0.2;
    let plus = law.mixed_state(&F, p + h, treatment).unwrap();
    let minus = law.mixed_state(&F, p - h, treatment).unwrap();
    near(
        (plus.potential_pa - minus.potential_pa) / (2.0 * h),
        state.pressure_residual,
        1.0,
        2e-8,
    );
    near(
        (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
        full.pressure_tangent_per_pa,
        1e-6,
        2e-8,
    );
    for i in 0..9 {
        near(
            (plus.first_piola_pa[i] - minus.first_piola_pa[i]) / (2.0 * h),
            full.pressure_coupling[i],
            1.0,
            2e-8,
        );
    }
}

#[test]
fn g3_elastic_and_frozen_history_mixed_blocks_are_the_full_saddle_hessian() {
    for law in models() {
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            derivative_gate(&law, treatment);
        }
    }
    let law = visco();
    let previous = [1.2, 0.1, 0.0, 0.02, 0.94, 0.0, 0.0, 0.0, 0.91];
    let (accepted, _) = law.step(&law.initial_state(), &previous, 0.04).unwrap();
    let snapshot = accepted.clone();
    for dt in [0.0, 0.025, 0.3] {
        let frozen = law.freeze(&accepted, dt).unwrap();
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            derivative_gate(&frozen, treatment);
        }
    }
    assert_eq!(accepted, snapshot);
}

#[test]
fn g0_pressure_elimination_recovers_each_original_energy_stress_and_tangent() {
    for law in models() {
        let p = pressure(&law, determinant(&F));
        let result = law
            .mixed_response(&F, p, VolumetricTreatment::OriginalMaterial)
            .unwrap();
        near(result.state.pressure_residual, 0.0, 1.0, 2e-14);
        near(result.state.potential_pa, law.energy(&F), 1e6, 2e-13);
        near(
            result.state.bulk_stored_pa,
            law.volumetric_energy(determinant(&F)),
            1e6,
            2e-13,
        );
        let original_p = law.piola(&F).unwrap();
        let original_a = law.tangent(&F).unwrap();
        for i in 0..9 {
            near(result.state.first_piola_pa[i], original_p[i], 1e6, 3e-13);
            for j in 0..9 {
                let condensed = result.tangent_pa[i][j]
                    - result.pressure_coupling[i] * result.pressure_coupling[j]
                        / result.pressure_tangent_per_pa;
                near(condensed, original_a[i][j], 1e6, 3e-12);
            }
        }
    }
    let law = visco();
    let (accepted, _) = law.step(&law.initial_state(), &F, 0.03).unwrap();
    let frozen = law.freeze(&accepted, 0.12).unwrap();
    let p = pressure(&law.equilibrium, determinant(&F));
    let mixed = frozen
        .mixed_response(&F, p, VolumetricTreatment::OriginalMaterial)
        .unwrap();
    let (original_p, original_value) = frozen.stress_and_potential(&F).unwrap();
    let original_a = frozen.consistent_tangent(&F).unwrap();
    near(mixed.state.potential_pa, original_value, 1e6, 1e-12);
    for i in 0..9 {
        near(mixed.state.first_piola_pa[i], original_p[i], 1e6, 1e-12);
        for j in 0..9 {
            let condensed = mixed.tangent_pa[i][j]
                - mixed.pressure_coupling[i] * mixed.pressure_coupling[j]
                    / mixed.pressure_tangent_per_pa;
            near(condensed, original_a[i][j], 1e6, 5e-12);
        }
    }
}

#[test]
fn g0_exact_constraint_bulk_independence_and_pressure_recovered_storage_are_distinct() {
    let a = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 6e4,
            c01: 4e4,
            kappa: 8e5,
        },
        4.0,
    )
    .unwrap();
    let b = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 6e4,
            c01: 4e4,
            kappa: 1e300,
        },
        4.0,
    )
    .unwrap();
    let x = a
        .mixed_response(&F, 1.1e5, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    let y = b
        .mixed_response(&F, 1.1e5, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    assert_eq!(
        x.state.potential_pa.to_bits(),
        y.state.potential_pa.to_bits()
    );
    assert_eq!(x.state.first_piola_pa, y.state.first_piola_pa);
    assert_eq!(x.tangent_pa, y.tangent_pa);
    assert_eq!(x.pressure_tangent_per_pa, 0.0);
    assert_eq!(x.state.bulk_stored_pa, 0.0);
    assert_eq!(x.state.bulk_volume_ratio, 1.0);
    let p = 1.1e5;
    let state = a
        .mixed_state(&F, p, VolumetricTreatment::OriginalMaterial)
        .unwrap();
    near(state.bulk_stored_pa, p * p / (2.0 * 8e5), 1e6, 1e-14);
    assert!((state.bulk_stored_pa - a.volumetric_energy(determinant(&F))).abs() > 100.0);
    near(
        state.potential_pa - (state.remainder_potential_pa + state.bulk_stored_pa),
        p * state.pressure_residual,
        1e6,
        1e-14,
    );
}

#[test]
fn g0_original_branches_remain_isochoric_and_physical_storage_is_not_incremental() {
    let law = visco();
    let jump = law.step_response(&law.initial_state(), &F, 0.0).unwrap();
    let before = jump.state.clone();
    let frozen = law.freeze(&before, 0.12).unwrap();
    let reference = frozen.remainder_response(&F).unwrap();
    let equilibrium = law.equilibrium.remainder_response(&F).unwrap();
    let scale = 1.3;
    let dilated = F.map(|v| scale * v);
    let response = frozen.remainder_response(&dilated).unwrap();
    let equilibrium_dilated = law.equilibrium.remainder_response(&dilated).unwrap();
    near(
        response.potential_pa - equilibrium_dilated.potential_pa,
        reference.potential_pa - equilibrium.potential_pa,
        1e5,
        2e-12,
    );
    for i in 0..9 {
        near(
            scale * (response.first_piola_pa[i] - equilibrium_dilated.first_piola_pa[i]),
            reference.first_piola_pa[i] - equilibrium.first_piola_pa[i],
            1e5,
            3e-12,
        );
    }
    let held = law.step_response(&before, &F, 0.12).unwrap();
    assert!(held.dissipated_pa > 0.0);
    let physical = law.mixed_remainder_free_energy(&held.state, &F).unwrap();
    near(
        physical,
        held.state.stored - law.equilibrium.volumetric_energy(determinant(&F)),
        1e6,
        2e-12,
    );
    assert!(reference.potential_pa > physical);
    assert_eq!(before, jump.state);
    for cv in &held.state.inverse_viscous_cauchy_green {
        let expanded = [
            cv[0], cv[3], cv[4], cv[3], cv[1], cv[5], cv[4], cv[5], cv[2],
        ];
        near(determinant(&expanded), 1.0, 1.0, 2e-12);
    }
    // A large bulk penalty cannot erase the independently evaluated remainder.
    let mut huge = law.clone();
    if let HyperelasticModel::MooneyRivlin { kappa, .. } = &mut huge.equilibrium.model {
        *kappa = 1e300;
    }
    let retained = huge.mixed_remainder_free_energy(&held.state, &F).unwrap();
    assert_eq!(retained.to_bits(), physical.to_bits());
}

#[test]
fn g3_mixed_rotation_and_pressure_unit_scaling_preserve_original_work() {
    let rotate = |f: [f64; 9]| {
        core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            match i {
                0 => 0.8 * f[j] - 0.6 * f[3 + j],
                1 => 0.6 * f[j] + 0.8 * f[3 + j],
                _ => f[6 + j],
            }
        })
    };
    for law in models() {
        let a = law
            .mixed_state(&F, 1.1e5, VolumetricTreatment::OriginalMaterial)
            .unwrap();
        let b = law
            .mixed_state(&rotate(F), 1.1e5, VolumetricTreatment::OriginalMaterial)
            .unwrap();
        near(a.potential_pa, b.potential_pa, 1e6, 2e-13);
        near(a.pressure_residual, b.pressure_residual, 1.0, 2e-14);
        for (actual, expected) in b.first_piola_pa.into_iter().zip(rotate(a.first_piola_pa)) {
            near(actual, expected, 1e6, 2e-12);
        }
    }
    let reference = pressure_potential(Volumetric::LogSquared { kappa: 8e5 }, 1.1e5).unwrap();
    for scale in [1e-12, 1e12] {
        let value =
            pressure_potential(Volumetric::LogSquared { kappa: 8e5 * scale }, 1.1e5 * scale)
                .unwrap();
        near(value.volume_ratio, reference.volume_ratio, 1.0, 1e-13);
        near(
            value.complementary_pa / scale,
            reference.complementary_pa,
            1e6,
            1e-13,
        );
        near(
            value.compliance_per_pa * scale,
            reference.compliance_per_pa,
            1e-5,
            1e-13,
        );
    }
}

#[test]
fn g0_coupled_laws_invalid_pressure_branches_and_stretch_domains_refuse() {
    let coupled = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2e5,
            lambda: 8e5,
        },
        4.0,
    )
    .unwrap();
    assert!(matches!(
        coupled.validate_mixed(),
        Err(MaterialError::Parameters { .. })
    ));
    for treatment in [
        VolumetricTreatment::OriginalMaterial,
        VolumetricTreatment::ExactIncompressible,
    ] {
        assert!(matches!(
            coupled.mixed_state(&I, 0.0, treatment),
            Err(MaterialError::Parameters { .. })
        ));
    }
    assert!(pressure_potential(Volumetric::Quadratic { kappa: 1.0 }, -1.0).is_err());
    assert!(pressure_potential(Volumetric::Quadratic { kappa: 0.0 }, 0.0).is_err());
    assert!(
        pressure_potential(
            Volumetric::LogSquared { kappa: 1.0 },
            fs_math::det::exp(-1.0)
        )
        .is_err()
    );
    assert!(pressure_potential(Volumetric::LogSquared { kappa: 1.0 }, 1.0).is_err());
    let log = models().remove(1);
    let expansion = [2.8, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert!(
        log.mixed_state(&expansion, 0.0, VolumetricTreatment::OriginalMaterial)
            .is_err()
    );
    assert!(
        log.mixed_state(&expansion, 0.0, VolumetricTreatment::ExactIncompressible)
            .is_ok()
    );
    for f in [
        [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        [4.0 * (1.0 + 1e-12), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    ] {
        assert!(
            log.mixed_state(&f, 0.0, VolumetricTreatment::ExactIncompressible)
                .is_err()
        );
    }
    assert!(
        log.mixed_state(&I, f64::NAN, VolumetricTreatment::OriginalMaterial)
            .is_err()
    );
}
