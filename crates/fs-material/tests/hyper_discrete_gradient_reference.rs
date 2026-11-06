//! Held-reference conservative strain work, angular momentum and exact derivatives.

use fs_material::hyper::discrete_gradient::{
    DiscreteGradientPoint, DiscreteGradientSettings, EnergyPart,
};
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};

const F0E: [f64; 9] = [1.17, 0.08, -0.02, 0.03, 0.91, 0.04, -0.01, 0.02, 1.06];
const F1E: [f64; 9] = [1.36, 0.17, 0.05, -0.04, 0.83, 0.09, 0.03, -0.02, 1.14];

fn settings() -> DiscreteGradientSettings {
    DiscreteGradientSettings {
        minimum_relative_strain_increment: 1e-7,
        absolute_work_tolerance_pa: 1e-8,
        relative_work_tolerance: 1e-12,
    }
}

fn reference() -> ReferenceDeformation {
    ReferenceDeformation::try_new([1.7, 0.24, -0.17, 0.08, 1.3, 0.19, 0.14, -0.11, 1.5]).unwrap()
}

fn laws() -> Vec<Hyperelastic> {
    vec![
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa: 9e5,
            },
            4.0,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: 2e5,
                    alpha: 1.3,
                }],
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

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn close(a: f64, b: f64, allowance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= allowance,
        "{a:.17e} versus {b:.17e}, allowance={allowance:e}"
    );
}

fn symmetric_moment(p: &[f64; 9], f0: &[f64; 9], f1: &[f64; 9]) {
    let fm: [f64; 9] = core::array::from_fn(|i| 0.5 * f0[i] + 0.5 * f1[i]);
    let moment: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..3).map(|a| p[3 * i + a] * fm[3 * j + a]).sum())
    });
    let scale = moment
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(1.0_f64, f64::max);
    for i in 0..3 {
        for j in 0..3 {
            close(moment[i][j], moment[j][i], 2e-13 * scale);
        }
    }
}

#[test]
fn accepted_reference_change_retains_original_small_work_and_source_part() {
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let law = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 900.0,
            c01: 600.0,
            kappa: 1e24,
        },
        4.0,
    )
    .unwrap();
    let same_law = law.clone();
    let before = ReferenceDeformation::try_new(identity.map(|x| 2.0 * x)).unwrap();
    let mut gradient = *before.deformation_gradient();
    gradient[1] = 4e-7;
    let after = ReferenceDeformation::try_new(gradient).unwrap();
    let total = identity.map(|x| 4.0 * x);
    let expected = 8.0 * 1500.0 * (2e-7_f64).powi(2);
    for part in [EnergyPart::Full, EnergyPart::Remainder] {
        let p0 =
            DiscreteGradientPoint::try_new_with_reference(&law, part, &before, total, settings())
                .unwrap();
        let p1 = DiscreteGradientPoint::try_new_with_reference(
            &same_law,
            part,
            &after,
            total,
            settings(),
        )
        .unwrap();
        if part == EnergyPart::Full {
            assert_eq!(p0.energy_before_pa(), p1.energy_before_pa());
        }
        close(
            p0.energy_increment_to(&p1).unwrap(),
            expected,
            3e-12 * expected,
        );
        close(
            p1.energy_increment_to(&p0).unwrap(),
            -expected,
            3e-12 * expected,
        );
        assert_eq!(p0.energy_increment_to(&p0).unwrap(), 0.0);
        assert_eq!(p0.accepted_deformation(), &total);
    }
    let full = DiscreteGradientPoint::try_new(&law, EnergyPart::Full, total, settings()).unwrap();
    let remainder =
        DiscreteGradientPoint::try_new(&law, EnergyPart::Remainder, total, settings()).unwrap();
    assert!(full.energy_increment_to(&remainder).is_err());
}

#[test]
fn accepted_source_changes_retain_both_endpoint_potentials() {
    let laws = laws();
    let g = reference();
    let before_f = multiply(&F0E, g.deformation_gradient());
    for part in [EnergyPart::Full, EnergyPart::Remainder] {
        let before =
            DiscreteGradientPoint::try_new_with_reference(&laws[0], part, &g, before_f, settings())
                .unwrap();
        let after = DiscreteGradientPoint::try_new(&laws[1], part, F1E, settings()).unwrap();
        let expected = after.energy_before_pa() - before.energy_before_pa();
        close(before.energy_increment_to(&after).unwrap(), expected, 2e-8);
        close(after.energy_increment_to(&before).unwrap(), -expected, 2e-8);
    }
    let bare0 =
        DiscreteGradientPoint::try_new(&laws[0], EnergyPart::Full, F0E, settings()).unwrap();
    let bare1 =
        DiscreteGradientPoint::try_new(&laws[0], EnergyPart::Full, F1E, settings()).unwrap();
    assert_eq!(
        bare0.energy_increment_to(&bare1).unwrap(),
        laws[0].energy_increment(&F0E, &F1E)
    );
}

#[test]
fn simultaneous_reference_and_stiffness_changes_retain_weak_signed_source_work() {
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let material = |factor| {
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: factor * 900.0,
                c01: factor * 600.0,
                kappa: 1e24,
            },
            4.0,
        )
        .unwrap()
    };
    let original = material(1.0);
    let changed = material(2.0);
    let g0 = ReferenceDeformation::try_new(identity.map(|x| 2.0 * x)).unwrap();
    let mut gradient = *g0.deformation_gradient();
    gradient[1] = 4e-7;
    let g1 = ReferenceDeformation::try_new(gradient).unwrap();
    let total = identity.map(|x| 4.0 * x);
    let expected = 8.0 * 3000.0 * (2e-7_f64).powi(2);
    for part in [EnergyPart::Full, EnergyPart::Remainder] {
        let before =
            DiscreteGradientPoint::try_new_with_reference(&original, part, &g0, total, settings())
                .unwrap();
        let after =
            DiscreteGradientPoint::try_new_with_reference(&changed, part, &g1, total, settings())
                .unwrap();
        if part == EnergyPart::Full {
            assert_eq!(before.energy_before_pa(), after.energy_before_pa());
        }
        close(
            before.energy_increment_to(&after).unwrap(),
            expected,
            3e-12 * expected,
        );
        close(
            after.energy_increment_to(&before).unwrap(),
            -expected,
            3e-12 * expected,
        );
    }
}

#[test]
fn reference_work_torque_and_nonsymmetric_endpoint_tangent_retain_native_laws() {
    let g = reference();
    let f0 = multiply(&F0E, g.deformation_gradient());
    let f1 = multiply(&F1E, g.deformation_gradient());
    let jg = g.volume_ratio();
    let mut nonsymmetric = false;
    for law in laws() {
        for part in [EnergyPart::Full, EnergyPart::Remainder] {
            let point =
                DiscreteGradientPoint::try_new_with_reference(&law, part, &g, f0, settings())
                    .unwrap();
            let response = point.response(&f1).unwrap();
            let energy = |f: &[f64; 9]| {
                jg * match part {
                    EnergyPart::Full => law.energy(&g.elastic_deformation(f)),
                    EnergyPart::Remainder => law.isochoric_energy(&g.elastic_deformation(f)),
                }
            };
            close(response.state.energy_before_pa, energy(&f0), 2e-7);
            close(response.state.energy_after_pa, energy(&f1), 2e-7);
            let work: f64 = (0..9)
                .map(|i| response.state.first_piola_pa[i] * (f1[i] - f0[i]))
                .sum();
            close(work, energy(&f1) - energy(&f0), 2e-7);
            close(work, response.state.work_pa, 1e-9);
            let allowance = settings().absolute_work_tolerance_pa
                + settings().relative_work_tolerance
                    * response
                        .state
                        .energy_before_pa
                        .abs()
                        .max(response.state.energy_after_pa.abs())
                        .max(response.state.work_pa.abs());
            assert!(response.state.work_defect_pa.abs() <= allowance);
            symmetric_moment(&response.state.first_piola_pa, &f0, &f1);
            for col in 0..9 {
                let h = 2e-6;
                let mut plus = f1;
                let mut minus = f1;
                plus[col] += h;
                minus[col] -= h;
                let pp = point.state(&plus).unwrap().first_piola_pa;
                let pm = point.state(&minus).unwrap().first_piola_pa;
                for row in 0..9 {
                    close(
                        response.endpoint_tangent_pa[row][col],
                        (pp[row] - pm[row]) / (2.0 * h),
                        0.8,
                    );
                    nonsymmetric |= (response.endpoint_tangent_pa[row][col]
                        - response.endpoint_tangent_pa[col][row])
                        .abs()
                        > 1.0;
                }
            }
            assert_eq!(point.accepted_deformation(), &f0);
        }
    }
    assert!(
        nonsymmetric,
        "the algorithmic endpoint derivative must not become a symmetric elastic Hessian"
    );
}

#[test]
fn reference_volume_work_and_derivative_use_total_volume_and_total_cofactor() {
    let g = reference();
    let law = laws().remove(0);
    let f0 = multiply(&F0E, g.deformation_gradient());
    let f1 = multiply(&F1E, g.deformation_gradient());
    let point = DiscreteGradientPoint::try_new_with_reference(
        &law,
        EnergyPart::Remainder,
        &g,
        f0,
        settings(),
    )
    .unwrap();
    let response = point.volume_response(&f1, 1e-12).unwrap();
    close(response.state.volume_ratio_before, determinant(&f0), 3e-14);
    close(response.state.volume_ratio_after, determinant(&f1), 3e-14);
    close(
        response.state.volume_work,
        determinant(&f1) - determinant(&f0),
        3e-14,
    );
    symmetric_moment(&response.state.discrete_cofactor, &f0, &f1);
    let h = 2e-6;
    for col in 0..9 {
        let mut plus = f1;
        let mut minus = f1;
        plus[col] += h;
        minus[col] -= h;
        let pp = point.volume_state(&plus, 1e-12).unwrap().discrete_cofactor;
        let pm = point.volume_state(&minus, 1e-12).unwrap().discrete_cofactor;
        close(
            response.state.endpoint_cofactor[col],
            (determinant(&plus) - determinant(&minus)) / (2.0 * h),
            2e-8,
        );
        for row in 0..9 {
            close(
                response.endpoint_jacobian[row][col],
                (pp[row] - pm[row]) / (2.0 * h),
                2e-8,
            );
        }
    }
    let natural = DiscreteGradientPoint::try_new_with_reference(
        &law,
        EnergyPart::Remainder,
        &g,
        *g.deformation_gradient(),
        settings(),
    )
    .unwrap()
    .volume_state(g.deformation_gradient(), 1e-12)
    .unwrap();
    assert_eq!(natural.volume_ratio_after, g.volume_ratio());
    assert_eq!(natural.volume_work, 0.0);
}

#[test]
fn reference_spatial_rotation_and_singular_midpoint_preserve_objective_work() {
    let g = reference();
    let q = [0.8, -0.6, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 1.0];
    let half_turn = [-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0];
    let f0 = multiply(&F0E, g.deformation_gradient());
    let f1 = multiply(&F1E, g.deformation_gradient());
    for law in laws() {
        let point = DiscreteGradientPoint::try_new_with_reference(
            &law,
            EnergyPart::Full,
            &g,
            f0,
            settings(),
        )
        .unwrap();
        let state = point.state(&f1).unwrap();
        let rotated = DiscreteGradientPoint::try_new_with_reference(
            &law,
            EnergyPart::Full,
            &g,
            multiply(&q, &f0),
            settings(),
        )
        .unwrap()
        .state(&multiply(&q, &f1))
        .unwrap();
        let expected = multiply(&q, &state.first_piola_pa);
        for i in 0..9 {
            close(rotated.first_piola_pa[i], expected[i], 3e-6);
        }
        close(rotated.energy_before_pa, state.energy_before_pa, 3e-7);
        close(rotated.energy_after_pa, state.energy_after_pa, 3e-7);
        close(rotated.work_pa, state.work_pa, 3e-7);
        let same_strain = point.response(&multiply(&half_turn, &f0)).unwrap();
        assert!(same_strain.state.used_midpoint_limit);
        close(same_strain.state.work_pa, 0.0, 1e-7);
        close(
            same_strain.state.energy_after_pa,
            same_strain.state.energy_before_pa,
            1e-7,
        );
        assert!(
            same_strain
                .endpoint_tangent_pa
                .iter()
                .flatten()
                .all(|v| v.is_finite())
        );
    }
}

#[test]
fn identity_is_exact_and_admission_uses_elastic_strain_without_losing_original_work_units() {
    let identity = ReferenceDeformation::identity();
    for law in laws() {
        for part in [EnergyPart::Full, EnergyPart::Remainder] {
            let bare = DiscreteGradientPoint::try_new(&law, part, F0E, settings()).unwrap();
            let referenced = DiscreteGradientPoint::try_new_with_reference(
                &law,
                part,
                &identity,
                F0E,
                settings(),
            )
            .unwrap();
            assert_eq!(
                bare.response(&F1E).unwrap(),
                referenced.response(&F1E).unwrap()
            );
            assert_eq!(
                bare.volume_response(&F1E, 1e-12).unwrap(),
                referenced.volume_response(&F1E, 1e-12).unwrap()
            );
        }
    }
    let law = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 6e4,
            c01: 4e4,
            kappa: 9e5,
        },
        1.15,
    )
    .unwrap();
    let g = ReferenceDeformation::try_new([5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0]).unwrap();
    assert!(
        DiscreteGradientPoint::try_new(
            &law,
            EnergyPart::Full,
            *g.deformation_gradient(),
            settings(),
        )
        .is_err()
    );
    let point = DiscreteGradientPoint::try_new_with_reference(
        &law,
        EnergyPart::Full,
        &g,
        *g.deformation_gradient(),
        settings(),
    )
    .unwrap();
    let relaxed = point.state(g.deformation_gradient()).unwrap();
    assert_eq!(relaxed.energy_after_pa, 0.0);
    close(
        relaxed.first_piola_pa.iter().map(|v| v.abs()).sum(),
        0.0,
        1e-7,
    );
    let invalid = [6.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0];
    assert!(point.response(&invalid).is_err());
    assert!(
        point
            .volume_response(g.deformation_gradient(), f64::NAN)
            .is_err()
    );
    assert_eq!(point.accepted_deformation(), g.deformation_gradient());
    let huge =
        ReferenceDeformation::try_new([1e100, 0.0, 0.0, 0.0, 1e100, 0.0, 0.0, 0.0, 1e100]).unwrap();
    let unrepresentable_allowance = DiscreteGradientSettings {
        absolute_work_tolerance_pa: f64::MIN_POSITIVE,
        ..settings()
    };
    assert!(
        DiscreteGradientPoint::try_new_with_reference(
            &law,
            EnergyPart::Full,
            &huge,
            *huge.deformation_gradient(),
            unrepresentable_allowance,
        )
        .is_err()
    );
    let coupled = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2e5,
            lambda: 8e5,
        },
        4.0,
    )
    .unwrap();
    assert!(
        DiscreteGradientPoint::try_new_with_reference(
            &coupled,
            EnergyPart::Remainder,
            &g,
            *g.deformation_gradient(),
            settings(),
        )
        .is_err()
    );
}
