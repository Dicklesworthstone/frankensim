//! Original full-tensor harmonic sections in an explicitly fixed shear profile.

use fs_material::elastic::OrthotropicElastic;
use fs_material::laminate::harmonic::{
    HarmonicLaminate, HarmonicLaminateBudget, HarmonicLaminateError,
};
use fs_material::laminate::{
    LaminaLayer, LaminateSection, LaminateShearRecovery, LayerLaw, ReferenceSurface,
};
use fs_material::visco::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use std::sync::Arc;

fn budget() -> HarmonicLaminateBudget {
    HarmonicLaminateBudget {
        maximum_layers: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_work: 20_000_000,
    }
}
fn elastic(scale: f64, limit: f64) -> OrthotropicElastic {
    OrthotropicElastic::new(
        [14e9 * scale, 7e9 * scale, 10e9 * scale],
        [0.21, 0.27, 0.18],
        [3e9 * scale, 2e9 * scale, 4e9 * scale],
        limit,
    ).unwrap()
}
fn matrix(diagonal: [f64; 6], vector: [f64; 6], scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| {
        scale * (vector[i] * vector[j] + if i == j { diagonal[i] } else { 0.0 })
    }))
}
fn prony(scale: f64) -> TensorProny {
    let equilibrium = matrix(
        [17.0, 13.0, 11.0, 7.0, 5.0, 9.0],
        [1.3, -0.7, 0.9, 0.4, -0.5, 0.6], 1e9 * scale,
    );
    let branch = matrix(
        [1.1, 0.9, 0.7, 0.5, 0.4, 0.3],
        [0.7, 0.4, -0.8, 0.9, -0.3, 0.5], 1e9 * scale,
    );
    let kelvin = matrix(
        [0.7, 0.8, 0.9, 0.5, 0.3, 0.4],
        [0.3, 0.9, -0.4, 0.7, 0.8, -0.6], 1e5 * scale,
    );
    TensorProny::new_mandel(equilibrium, vec![PronyBranch::new(branch, 0.014).unwrap()], 0.03)
        .unwrap().with_kelvin_viscosity_mandel(kelvin).unwrap()
        .with_band((0.0, 100.0)).unwrap()
}
fn layer(law: LayerLaw, thickness_m: f64, density_kg_m3: f64, orientation_rad: f64) -> LaminaLayer {
    LaminaLayer { law, thickness_m, density_kg_m3, orientation_rad, material: None }
}
fn reference(angle: f64) -> Arc<LaminateSection> {
    Arc::new(LaminateSection::try_new_equilibrated_shear(
        vec![
            layer(LayerLaw::Elastic(elastic(1.0, 0.03)), 0.0026, 450.0, 0.31 + angle),
            layer(LayerLaw::Elastic(elastic(0.6, 0.03)), 0.0002, 1100.0, -0.47 + angle),
        ],
        ReferenceSurface::Offset(0.0008),
        0.27 + angle,
    ).unwrap())
}
fn actual(angle: f64) -> Arc<LaminateSection> {
    Arc::new(LaminateSection::try_new(
        vec![
            layer(LayerLaw::Prony(prony(1.0)), 0.0026, 490.0, 0.31 + angle),
            layer(LayerLaw::Prony(prony(0.6)), 0.0002, 1300.0, -0.47 + angle),
        ],
        ReferenceSurface::Offset(0.0008),
    ).unwrap())
}
fn source(angle: f64) -> HarmonicLaminate {
    let shear = Arc::new(LaminateShearRecovery::try_new(reference(angle)).unwrap());
    HarmonicLaminate::try_new(actual(angle), shear, budget(), || false).unwrap()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual {actual:?}, expected {expected:?}");
}
fn motion() -> [C64; 8] {
    [
        C64::new(2e-6, -3e-7), C64::new(-1e-6, 4e-7), C64::new(7e-7, 9e-7),
        C64::new(8e-4, -2e-4), C64::new(-3e-4, 1e-4), C64::new(5e-4, 4e-4),
        C64::new(3e-6, -1e-6), C64::new(-2e-6, 4e-7),
    ]
}

#[test]
fn original_elastic_section_recovers_full_abd_h_and_actual_density_moments() {
    let section = reference(0.0);
    let shear = Arc::new(LaminateShearRecovery::try_new(Arc::clone(&section)).unwrap());
    let model = HarmonicLaminate::try_new(Arc::clone(&section), Arc::clone(&shear), budget(), || false).unwrap();
    assert!(core::ptr::eq(model.section(), section.as_ref()));
    assert!(core::ptr::eq(model.shear_trial(), shear.as_ref()));
    assert_eq!(model.inertia_moments(), [section.i0_kg_m2(), section.i1_kg_m(), section.i2_kg()]);
    let prepared = model.prepare_harmonic(37.0, budget(), || false).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            for (row, column, expected) in [
                (i, j, section.a_pa_m()[3 * i + j]),
                (i, j + 3, section.b_pa_m2()[3 * i + j]),
                (i + 3, j, section.b_pa_m2()[3 * i + j]),
                (i + 3, j + 3, section.d_pa_m3()[3 * i + j]),
            ] {
                near(C64::from_re(prepared.storage()[8 * row + column]), C64::from_re(expected), 2e-11);
            }
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            near(C64::from_re(prepared.storage()[8 * (i + 6) + j + 6]),
                C64::from_re(section.shear_pa_m()[2 * i + j]), 2e-11);
        }
    }
    assert_eq!(prepared.loss(), &[0.0; 64]);
    for i in 0..6 {
        for j in 6..8 { assert_eq!(prepared.storage()[8 * i + j], 0.0); }
    }
    for layer in 0..2 {
        for fraction in [0.0, 0.23, 0.5, 1.0] {
            let gamma = [1e-4, -6e-5];
            let original = shear.from_section_strain(layer, fraction, gamma).unwrap();
            let mut generalized = [C64::ZERO; 8];
            generalized[6] = C64::from_re(gamma[0]);
            generalized[7] = C64::from_re(gamma[1]);
            let strain = prepared.recover_strain(layer, fraction, &generalized).unwrap();
            near(strain[5], C64::from_re(original.material_engineering_strain[0] / 2.0_f64.sqrt()), 3e-15);
            near(strain[4], C64::from_re(original.material_engineering_strain[1] / 2.0_f64.sqrt()), 3e-15);
        }
    }
}

// Independent reference integration: use the original physical traction/compliance
// observer at FIVE Gauss points, original full tensor stress and one scalar
// normal-stress constraint. No new section map, Bernstein controls or reduced
// material matrix is used to construct this oracle.
fn direct_operator(model: &HarmonicLaminate, omega: f64) -> [[C64; 8]; 8] {
    let nodes = [-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640];
    let weights = [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891];
    let mut result = [[C64::ZERO; 8]; 8];
    for (layer, actual) in model.section().layers().iter().enumerate() {
        let LayerLaw::Prony(law) = &actual.law else { panic!("oracle fixture requires original Prony") };
        let (storage, loss) = law.complex_modulus(omega).unwrap();
        let c: [[C64; 6]; 6] = core::array::from_fn(|i| core::array::from_fn(|j| C64::new(storage[i][j], -loss[i][j])));
        let transform = model.section().principal_strain_transform(layer).unwrap();
        let geometry = &model.section().layer_integrals()[layer];
        let thickness = geometry.z_top_m - geometry.z_bottom_m;
        for point in 0..5 {
            let fraction = 0.5 * (1.0 + nodes[point]);
            let z = geometry.z_bottom_m + fraction * thickness;
            let mut strains = [[C64::ZERO; 6]; 8];
            for column in 0..8 {
                let mut generalized = [0.0; 8];
                generalized[column] = 1.0;
                let membrane: [f64; 3] = core::array::from_fn(|i| generalized[i] + z * generalized[i + 3]);
                let principal: [f64; 3] = core::array::from_fn(|i|
                    (0..3).map(|j| transform[i][j] * membrane[j]).sum()
                );
                let shear = model.shear_trial().from_section_strain(layer, fraction, [generalized[6], generalized[7]])
                    .unwrap().material_engineering_strain;
                strains[column] = [
                    C64::from_re(principal[0]), C64::from_re(principal[1]), C64::ZERO,
                    C64::from_re(principal[2] / 2.0_f64.sqrt()),
                    C64::from_re(shear[1] / 2.0_f64.sqrt()),
                    C64::from_re(shear[0] / 2.0_f64.sqrt()),
                ];
                strains[column][2] = -(0..6).filter(|&i| i != 2)
                    .fold(C64::ZERO, |sum, i| sum + c[2][i] * strains[column][i]) / c[2][2];
            }
            let stress: [[C64; 6]; 8] = strains.map(|strain| core::array::from_fn(|i|
                (0..6).fold(C64::ZERO, |sum, j| sum + c[i][j] * strain[j])
            ));
            for i in 0..8 {
                for j in 0..8 {
                    let work = (0..6).fold(C64::ZERO, |sum, k| sum + strains[i][k].conj() * stress[j][k]);
                    result[i][j] = result[i][j] + work.scale(0.5 * thickness * weights[point]);
                }
            }
        }
    }
    result
}

#[test]
fn full_anisotropy_keeps_shear_cross_work_and_matches_original_physical_stress_integrals() {
    let model = source(0.0);
    let prepared = model.prepare_harmonic(37.0, budget(), || false).unwrap();
    let oracle = direct_operator(&model, 37.0);
    for i in 0..8 {
        for j in 0..8 {
            near(C64::new(prepared.storage()[8 * i + j], -prepared.loss()[8 * i + j]), oracle[i][j], 3e-11);
            assert_eq!(prepared.storage()[8 * i + j], prepared.storage()[8 * j + i]);
            assert_eq!(prepared.loss()[8 * i + j], prepared.loss()[8 * j + i]);
        }
    }
    assert!(prepared.storage()[6].abs() > 1e3);
    assert!(prepared.storage()[8 * 3 + 6].abs() > 1.0);
    assert!(prepared.storage()[3].abs() > 1e3);
    assert_eq!(model.inertia_moments(), [model.section().i0_kg_m2(), model.section().i1_kg_m(), model.section().i2_kg()]);
    assert_ne!(model.inertia_moments()[0], model.shear_trial().section().i0_kg_m2());
    let generalized = motion();
    let resultant = prepared.resultants(&generalized, || false).unwrap();
    let mut work = C64::ZERO;
    for i in 0..8 {
        let expected = (0..8).fold(C64::ZERO, |sum, j| sum + oracle[i][j] * generalized[j]);
        near(resultant[i], expected, 3e-11);
        work = work + generalized[i].conj() * resultant[i];
    }
    assert!(work.re > 0.0 && work.im < 0.0);
    for layer in 0..2 {
        for fraction in [0.0, 0.37, 1.0] {
            let strain = prepared.recover_strain(layer, fraction, &generalized).unwrap();
            let material = prepared.layer_material(layer).unwrap();
            let stress = material.source().prepare_harmonic_material(37.0).unwrap().evaluate(&strain).unwrap().stress_mandel_pa;
            let scale = stress.iter().map(|z| z.abs()).fold(1.0_f64, f64::max);
            assert!(stress[2].abs() < 1e-12 * scale);
        }
    }
    // General anisotropic coupling remains at the physical exterior: a fixed
    // reference profile is not silently relabelled a new traction-free solve.
    let strain = prepared.recover_strain(0, 0.0, &generalized).unwrap();
    let stress = prepared.layer_material(0).unwrap().source()
        .prepare_harmonic_material(37.0).unwrap().evaluate(&strain).unwrap().stress_mandel_pa;
    assert!(stress[4].abs() + stress[5].abs() > 1.0);
}

#[test]
fn exact_section_frequency_tangent_holds_the_declared_profile_and_geometry_fixed() {
    let model = source(0.0);
    let omega = 37.0;
    let h = 1e-3;
    let center = model.prepare_harmonic(omega, budget(), || false).unwrap();
    let (storage_rate, loss_rate) = center.frequency_derivative(budget(), || false).unwrap();
    let left = model.prepare_harmonic(omega - h, budget(), || false).unwrap();
    let right = model.prepare_harmonic(omega + h, budget(), || false).unwrap();
    for i in 0..64 {
        let difference = C64::new(
            (right.storage()[i] - left.storage()[i]) / (2.0 * h),
            (right.loss()[i] - left.loss()[i]) / (2.0 * h),
        );
        near(C64::new(storage_rate[i], loss_rate[i]), difference, 2e-7);
    }
    assert!(model.prepare_harmonic(101.0, budget(), || false).is_err());
    assert!(model.prepare_harmonic(-1.0, budget(), || false).is_err());
    let dc = model.prepare_harmonic(0.0, budget(), || false).unwrap();
    assert_eq!(dc.loss(), &[0.0; 64]);
    assert!(dc.frequency_derivative(budget(), || false).is_ok());
    let mut imaginary = [C64::ZERO; 8];
    imaginary[6].im = 1e-6;
    assert!(dc.resultants(&imaginary, || false).is_err());
}

#[test]
fn whole_thickness_admission_detects_interior_shear_hidden_by_traction_free_faces() {
    let section = Arc::new(LaminateSection::try_new_equilibrated_shear(
        vec![layer(LayerLaw::Elastic(elastic(1.0, 0.01)), 0.003, 450.0, 0.0)],
        ReferenceSurface::Midsurface, 0.0,
    ).unwrap());
    let shear = Arc::new(LaminateShearRecovery::try_new(Arc::clone(&section)).unwrap());
    let model = HarmonicLaminate::try_new(section, shear, budget(), || false).unwrap();
    let prepared = model.prepare_harmonic(10.0, budget(), || false).unwrap();
    let mut generalized = [C64::ZERO; 8];
    generalized[6] = C64::from_re(0.016);
    assert!(prepared.recover_strain(0, 0.0, &generalized).is_ok());
    assert!(prepared.recover_strain(0, 1.0, &generalized).is_ok());
    assert!(prepared.recover_strain(0, 0.5, &generalized).is_err());
    assert!(prepared.admit_generalized_strain(&generalized, || false).is_err());
    assert!(prepared.resultants(&generalized, || false).is_err());
    generalized[6] = C64::new(0.002, -0.001);
    assert!(prepared.admit_generalized_strain(&generalized, || false).is_ok());
    assert!(prepared.recover_strain(1, 0.5, &generalized).is_err());
    assert!(prepared.recover_strain(0, 1.01, &generalized).is_err());
}

fn rotated_motion(original: [C64; 8], angle: f64) -> [C64; 8] {
    let (s, c) = angle.sin_cos();
    let rotate = |v: [C64; 3]| [
        v[0].scale(c * c) + v[1].scale(s * s) - v[2].scale(c * s),
        v[0].scale(s * s) + v[1].scale(c * c) + v[2].scale(c * s),
        (v[0] - v[1]).scale(2.0 * c * s) + v[2].scale(c * c - s * s),
    ];
    let a = rotate([original[0], original[1], original[2]]);
    let b = rotate([original[3], original[4], original[5]]);
    [
        a[0], a[1], a[2], b[0], b[1], b[2],
        original[6].scale(c) - original[7].scale(s),
        original[6].scale(s) + original[7].scale(c),
    ]
}

#[test]
fn corotated_source_axes_probe_direction_and_physical_strain_preserve_work() {
    let original = source(0.0);
    let angle = 0.63;
    let rotated = source(angle);
    let a = original.prepare_harmonic(37.0, budget(), || false).unwrap();
    let b = rotated.prepare_harmonic(37.0, budget(), || false).unwrap();
    let qa = motion();
    let qb = rotated_motion(qa, angle);
    let fa = a.resultants(&qa, || false).unwrap();
    let fb = b.resultants(&qb, || false).unwrap();
    let work = |q: [C64; 8], f: [C64; 8]| (0..8).fold(C64::ZERO, |sum, i| sum + q[i].conj() * f[i]);
    near(work(qa, fa), work(qb, fb), 3e-11);
    for layer in 0..2 {
        for fraction in [0.0, 0.29, 0.81, 1.0] {
            let ea = a.recover_strain(layer, fraction, &qa).unwrap();
            let eb = b.recover_strain(layer, fraction, &qb).unwrap();
            for i in 0..6 { near(ea[i], eb[i], 3e-14); }
        }
    }
}

#[test]
fn incompatible_sources_budgets_and_cancellation_refuse_without_losing_reuse() {
    let actual = actual(0.0);
    let reference = reference(0.0);
    let shear = Arc::new(LaminateShearRecovery::try_new(Arc::clone(&reference)).unwrap());
    for b in [
        HarmonicLaminateBudget { maximum_layers: 1, ..budget() },
        HarmonicLaminateBudget { maximum_payload_bytes: 0, ..budget() },
        HarmonicLaminateBudget { maximum_work: 0, ..budget() },
    ] {
        assert!(HarmonicLaminate::try_new(Arc::clone(&actual), Arc::clone(&shear), b, || false).is_err());
    }
    assert!(matches!(
        HarmonicLaminate::try_new(Arc::clone(&actual), Arc::clone(&shear), budget(), || true),
        Err(HarmonicLaminateError::Cancelled)
    ));
    let mut changed = actual.layers().to_vec();
    changed[0].orientation_rad += 0.1;
    let mismatch = Arc::new(LaminateSection::try_new(changed, actual.reference()).unwrap());
    assert!(HarmonicLaminate::try_new(mismatch, Arc::clone(&shear), budget(), || false).is_err());
    let mut reduced = reference.layers().to_vec();
    reduced[0].law = LayerLaw::PlaneStress {
        stiffness_pa: [14e9, 1e9, 0.0, 1e9, 7e9, 0.0, 0.0, 0.0, 3e9],
        transverse_shear_pa: [4e9, 2e9],
    };
    let reduced = Arc::new(LaminateSection::try_new(reduced, reference.reference()).unwrap());
    assert!(HarmonicLaminate::try_new(reduced, Arc::clone(&shear), budget(), || false).is_err());
    let model = HarmonicLaminate::try_new(actual, shear, budget(), || false).unwrap();
    assert!(model.prepare_harmonic(37.0, HarmonicLaminateBudget { maximum_work: 0, ..budget() }, || false).is_err());
    let prepared = model.prepare_harmonic(37.0, budget(), || false).unwrap();
    let safe = prepared.resultants(&motion(), || false).unwrap();
    assert!(matches!(prepared.resultants(&motion(), || true), Err(HarmonicLaminateError::Cancelled)));
    assert!(matches!(prepared.frequency_derivative(budget(), || true), Err(HarmonicLaminateError::Cancelled)));
    assert!(prepared.frequency_derivative(HarmonicLaminateBudget { maximum_work: 0, ..budget() }, || false).is_err());
    let mut bad = motion();
    bad[0].re = f64::NAN;
    assert!(prepared.resultants(&bad, || false).is_err());
    assert_eq!(safe, prepared.resultants(&motion(), || false).unwrap());
}
