//! G0/G3: objective finite-step work, true endpoint derivatives and roundoff admission.

use fs_material::hyper::discrete_gradient::{
    DiscreteGradientPoint, DiscreteGradientSettings, EnergyPart,
};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};

const ID: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
const F0: [f64; 9] = [1.17, 0.08, -0.02, 0.03, 0.91, 0.04, -0.01, 0.02, 1.06];
const F1: [f64; 9] = [1.36, 0.17, 0.05, -0.04, 0.83, 0.09, 0.03, -0.02, 1.14];

fn settings() -> DiscreteGradientSettings {
    DiscreteGradientSettings {
        minimum_relative_strain_increment: 1e-7,
        absolute_work_tolerance_pa: 1e-8,
        relative_work_tolerance: 1e-12,
    }
}

fn laws() -> Vec<Hyperelastic> {
    vec![
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2e5,
                lambda: 8e5,
            },
            4.,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa: 9e5,
            },
            4.,
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
            4.,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 5e4,
                k1: 1.2e5,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0., 0.6], [0., 0.6, 0.8]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 8e5 },
            },
            4.,
        )
        .unwrap(),
    ]
}

fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:.17e} vs {b:.17e}, scale={scale:e}, tolerance={tolerance:e}"
    );
}
fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}
fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}
fn symmetric_spatial_stress(p: &[f64; 9], f0: &[f64; 9], f1: &[f64; 9]) {
    let fm: [f64; 9] = core::array::from_fn(|i| 0.5 * f0[i] + 0.5 * f1[i]);
    let stress: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..3).map(|a| p[3 * i + a] * fm[3 * j + a]).sum())
    });
    let scale = stress
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(1.0_f64, f64::max);
    for i in 0..3 {
        for j in 0..3 {
            close(stress[i][j], stress[j][i], scale, 2e-13);
        }
    }
}

#[test]
fn g0_original_energies_have_exact_algorithmic_work_and_symmetric_material_moment() {
    let mut correction_is_live = false;
    for (index, law) in laws().iter().enumerate() {
        for part in [EnergyPart::Full, EnergyPart::Remainder] {
            if index == 0 && part == EnergyPart::Remainder {
                continue;
            }
            let point = DiscreteGradientPoint::try_new(law, part, F0, settings()).unwrap();
            let response = point.response(&F1).unwrap();
            let trial = point.state(&F1).unwrap();
            assert!(!response.state.used_midpoint_limit);
            let original = match part {
                EnergyPart::Full => law.energy(&F1),
                EnergyPart::Remainder => law.isochoric_energy(&F1),
            };
            close(
                response.state.energy_after_pa,
                original,
                original.abs().max(1.),
                2e-13,
            );
            close(
                response.state.work_pa,
                response.state.energy_after_pa - response.state.energy_before_pa,
                original.abs().max(1.),
                2e-13,
            );
            for i in 0..9 {
                close(
                    trial.first_piola_pa[i],
                    response.state.first_piola_pa[i],
                    1e6,
                    2e-13,
                );
            }
            symmetric_spatial_stress(&response.state.first_piola_pa, &F0, &F1);
            let fm = core::array::from_fn(|i| 0.5 * F0[i] + 0.5 * F1[i]);
            if part == EnergyPart::Full {
                let ordinary = law.piola(&fm).unwrap();
                correction_is_live |=
                    (0..9).any(|i| (ordinary[i] - response.state.first_piola_pa[i]).abs() > 1.0);
            }
            assert_eq!(point.accepted_deformation(), &F0);
        }
    }
    assert!(
        correction_is_live,
        "finite-step algorithm must not collapse to F-midpoint elasticity"
    );
}

#[test]
fn g0_consistent_endpoint_jacobians_match_finite_differences_without_symmetrizing() {
    let h = 2e-6;
    let mut nonsymmetric = false;
    let mut worst = 0.0_f64;
    for law in laws() {
        let point = DiscreteGradientPoint::try_new(&law, EnergyPart::Full, F0, settings()).unwrap();
        let response = point.response(&F1).unwrap();
        let volume = point.volume_response(&F1, 1e-12).unwrap();
        symmetric_spatial_stress(&volume.state.discrete_cofactor, &F0, &F1);
        close(
            volume.state.volume_ratio_before,
            determinant(&F0),
            1.,
            2e-14,
        );
        close(volume.state.volume_ratio_after, determinant(&F1), 1., 2e-14);
        close(
            volume.state.volume_work,
            determinant(&F1) - determinant(&F0),
            1.,
            2e-14,
        );
        for col in 0..9 {
            let mut plus = F1;
            let mut minus = F1;
            plus[col] += h;
            minus[col] -= h;
            let pp = point.state(&plus).unwrap().first_piola_pa;
            let pm = point.state(&minus).unwrap().first_piola_pa;
            let jp = point.volume_state(&plus, 1e-12).unwrap().discrete_cofactor;
            let jm = point.volume_state(&minus, 1e-12).unwrap().discrete_cofactor;
            close(
                volume.state.endpoint_cofactor[col],
                (determinant(&plus) - determinant(&minus)) / (2. * h),
                1.,
                2e-9,
            );
            for row in 0..9 {
                let fd = (pp[row] - pm[row]) / (2. * h);
                worst = worst.max((fd - response.endpoint_tangent_pa[row][col]).abs() / 1e6);
                close(response.endpoint_tangent_pa[row][col], fd, 1e6, 2e-7);
                close(
                    volume.endpoint_jacobian[row][col],
                    (jp[row] - jm[row]) / (2. * h),
                    1.,
                    2e-8,
                );
                nonsymmetric |= (response.endpoint_tangent_pa[row][col]
                    - response.endpoint_tangent_pa[col][row])
                    .abs()
                    > 1.0;
            }
        }
    }
    assert!(
        nonsymmetric,
        "the true algorithmic tangent is not an elastic Hessian"
    );
    println!("objective discrete-gradient tangent worst normalized difference={worst:e}");
}

#[test]
fn g3_frame_covariance_and_equal_strains_under_finite_rotation_keep_true_c_chain() {
    let q = [0.8, -0.6, 0., 0.6, 0.8, 0., 0., 0., 1.];
    for law in laws() {
        let before = DiscreteGradientPoint::try_new(&law, EnergyPart::Full, F0, settings())
            .unwrap()
            .state(&F1)
            .unwrap();
        let f0 = multiply(&q, &F0);
        let f1 = multiply(&q, &F1);
        let rotated = DiscreteGradientPoint::try_new(&law, EnergyPart::Full, f0, settings())
            .unwrap()
            .state(&f1)
            .unwrap();
        let expected = multiply(&q, &before.first_piola_pa);
        for i in 0..9 {
            close(rotated.first_piola_pa[i], expected[i], 1e6, 3e-12);
        }
        close(rotated.energy_after_pa, before.energy_after_pa, 1e6, 3e-13);

        let coincident = DiscreteGradientPoint::try_new(&law, EnergyPart::Full, F0, settings())
            .unwrap()
            .response(&F0)
            .unwrap();
        let original = law.piola(&F0).unwrap();
        let tangent = law.tangent(&F0).unwrap();
        assert!(coincident.state.used_midpoint_limit);
        for i in 0..9 {
            close(coincident.state.first_piola_pa[i], original[i], 1e6, 2e-13);
            for j in 0..9 {
                close(
                    coincident.endpoint_tangent_pa[i][j],
                    0.5 * tangent[i][j],
                    1e6,
                    2e-12,
                );
            }
        }
    }

    // A half-turn leaves C unchanged and makes F_mid singular. The derivative
    // still uses F_mid and dC1/dF1; half the native endpoint Hessian is WRONG.
    let law = &laws()[0];
    let f0 = [1.2, 0., 0., 0., 0.8, 0., 0., 0., 1.05];
    let f1 = [-1.2, 0., 0., 0., -0.8, 0., 0., 0., 1.05];
    let settings = DiscreteGradientSettings {
        minimum_relative_strain_increment: 1e-4,
        ..settings()
    };
    let point = DiscreteGradientPoint::try_new(law, EnergyPart::Full, f0, settings).unwrap();
    let response = point.response(&f1).unwrap();
    assert!(response.state.used_midpoint_limit);
    close(response.state.work_pa, 0.0, 1., 1e-12);
    symmetric_spatial_stress(&response.state.first_piola_pa, &f0, &f1);
    let h = 1e-7;
    for col in 0..9 {
        let mut fp = f1;
        fp[col] += h;
        let mut fm = f1;
        fm[col] -= h;
        let p = point.state(&fp).unwrap().first_piola_pa;
        let m = point.state(&fm).unwrap().first_piola_pa;
        for row in 0..9 {
            close(
                response.endpoint_tangent_pa[row][col],
                (p[row] - m[row]) / (2. * h),
                1e6,
                2e-8,
            );
        }
    }
}

#[test]
fn g0_small_increment_work_is_reported_and_refused_when_its_budget_is_insufficient() {
    let law = &laws()[0];
    let declared = DiscreteGradientSettings {
        minimum_relative_strain_increment: 1e-3,
        absolute_work_tolerance_pa: 1e-3,
        relative_work_tolerance: 0.0,
    };
    let mut tiny = F0;
    tiny[0] += 1e-4;
    let point = DiscreteGradientPoint::try_new(law, EnergyPart::Full, F0, declared).unwrap();
    let limit = point.state(&tiny).unwrap();
    assert!(limit.used_midpoint_limit);
    assert!(limit.work_defect_pa.abs() > 1e-10 && limit.work_defect_pa.abs() < 1e-3);
    let strict = DiscreteGradientSettings {
        absolute_work_tolerance_pa: 0.0,
        ..declared
    };
    assert!(
        DiscreteGradientPoint::try_new(law, EnergyPart::Full, F0, strict)
            .unwrap()
            .state(&tiny)
            .is_err()
    );
    tiny = F0;
    tiny[0] += 1e-10;
    let small = point.state(&tiny).unwrap();
    let native = law.piola(&F0).unwrap();
    for i in 0..9 {
        close(small.first_piola_pa[i], native[i], 1e6, 2e-9);
    }
    assert_eq!(point.accepted_deformation(), &F0);
    for invalid in [
        DiscreteGradientSettings {
            minimum_relative_strain_increment: f64::NAN,
            ..declared
        },
        DiscreteGradientSettings {
            absolute_work_tolerance_pa: -1.,
            ..declared
        },
        DiscreteGradientSettings {
            relative_work_tolerance: f64::INFINITY,
            ..declared
        },
    ] {
        assert!(DiscreteGradientPoint::try_new(law, EnergyPart::Full, F0, invalid).is_err());
    }
    for f in [
        [0.; 9],
        [-1., 0., 0., 0., 1., 0., 0., 0., 1.],
        [5., 0., 0., 0., 1., 0., 0., 0., 1.],
        [f64::NAN; 9],
    ] {
        assert!(point.state(&f).is_err());
    }
    assert!(point.volume_state(&F1, -1.).is_err());
    assert!(DiscreteGradientPoint::try_new(law, EnergyPart::Remainder, F0, settings()).is_err());
}

#[test]
fn g0_remainder_never_evaluates_unused_bulk_and_full_invariant_fibres_retain_volume() {
    let mr = |kappa| {
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 1e3,
                c01: 0.,
                kappa,
            },
            4.,
        )
        .unwrap()
    };
    let (small, huge) = (mr(1e3), mr(1e308));
    let f = [3., 0., 0., 0., 3., 0., 0., 0., 3.];
    assert!(
        DiscreteGradientPoint::try_new(&huge, EnergyPart::Full, ID, settings())
            .unwrap()
            .response(&f)
            .is_err()
    );
    let a = DiscreteGradientPoint::try_new(&small, EnergyPart::Remainder, ID, settings())
        .unwrap()
        .response(&f)
        .unwrap();
    let b = DiscreteGradientPoint::try_new(&huge, EnergyPart::Remainder, ID, settings())
        .unwrap()
        .response(&f)
        .unwrap();
    assert_eq!(a, b);

    let make = |invariants| {
        Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 1e3,
                k1: 2e3,
                k2: 1.,
                dispersion: 0.,
                fibres: vec![[1., 0., 0.]],
                invariants,
                volumetric: Volumetric::Quadratic { kappa: 1e300 },
            },
            4.,
        )
        .unwrap()
    };
    let (full, iso) = (
        make(FibreInvariants::Full),
        make(FibreInvariants::Isochoric),
    );
    let f = [1.2, 0., 0., 0., 1.2, 0., 0., 0., 1.2];
    let full = DiscreteGradientPoint::try_new(&full, EnergyPart::Remainder, ID, settings())
        .unwrap()
        .state(&f)
        .unwrap();
    let iso = DiscreteGradientPoint::try_new(&iso, EnergyPart::Remainder, ID, settings())
        .unwrap()
        .state(&f)
        .unwrap();
    assert!(full.energy_after_pa > 100.);
    close(iso.energy_after_pa, 0., 1., 1e-10);
    assert!(full.work_pa > 100.);
}
