//! Independent finite-deformation checks for a prescribed material reference.

use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};

const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
// Nonsymmetric prescribed shear catches inverse/cofactor transpose mistakes
// that a symmetric thermal reference cannot expose.
const G: [f64; 9] = [1.12, 0.17, -0.04, 0.0, 0.93, 0.08, 0.0, 0.0, 1.07];
const FE: [f64; 9] = [1.06, 0.09, -0.03, 0.02, 0.97, 0.06, -0.01, 0.03, 1.04];

#[test]
fn native_source_amplitudes_retain_stiffness_work_beneath_unchanged_bulk() {
    let delta = 2.0_f64.powi(-20);
    for gamma in [1e-5_f64, 0.4] {
        let mut f = I.map(|x| 2.0 * x);
        f[1] = 2.0 * gamma;
        let volumetric = Volumetric::Quadratic { kappa: 1e24 };
        let models = [
            (
                HyperelasticModel::NeoHookean {
                    mu: 1024.0,
                    lambda: 1e24,
                },
                delta * (4.5 + 2.0 * gamma * gamma - 8.0_f64.ln()),
            ),
            (
                HyperelasticModel::MooneyRivlin {
                    c10: 1024.0,
                    c01: 512.0,
                    kappa: 1e24,
                },
                3.0 * delta * gamma * gamma,
            ),
            (
                HyperelasticModel::Gent {
                    mu: 1024.0,
                    locking_invariant: 10.0,
                    volumetric,
                },
                -5.0 * delta * (-gamma * gamma / 10.0).ln_1p(),
            ),
            (
                HyperelasticModel::Ogden {
                    terms: vec![
                        OgdenTerm {
                            mu_pa: 1024.0,
                            alpha: 1.3,
                        },
                        OgdenTerm {
                            mu_pa: -512.0,
                            alpha: -2.0,
                        },
                    ],
                    volumetric,
                },
                // Isochoric simple-shear stretches are exp(+-asinh(gamma/2)), 1.
                4.0 * delta / 1.3 * (0.65 * (0.5 * gamma).asinh()).sinh().powi(2)
                    + delta * gamma * gamma,
            ),
            (
                HyperelasticModel::HolzapfelGasserOgden {
                    c10: 1024.0,
                    k1: 512.0,
                    k2: 100.0,
                    dispersion: 0.0,
                    fibres: vec![[0.0, 1.0, 0.0]],
                    invariants: FibreInvariants::Isochoric,
                    volumetric,
                },
                delta * gamma * gamma + delta / 100.0 * (100.0 * gamma.powi(4)).exp_m1(),
            ),
        ];
        for (model, expected) in models {
            let original = Hyperelastic::new(model, 4.0).unwrap();
            let mut changed = original.clone();
            match &mut changed.model {
                HyperelasticModel::NeoHookean { mu, .. } | HyperelasticModel::Gent { mu, .. } => {
                    *mu += delta;
                }
                HyperelasticModel::MooneyRivlin { c10, c01, .. } => {
                    *c10 += delta;
                    *c01 += 2.0 * delta;
                }
                HyperelasticModel::Ogden { terms, .. } => {
                    terms[0].mu_pa += delta;
                    terms[1].mu_pa -= 2.0 * delta;
                }
                HyperelasticModel::HolzapfelGasserOgden { c10, k1, .. } => {
                    *c10 += delta;
                    *k1 += 2.0 * delta;
                }
            }
            assert_eq!(changed.energy(&f), original.energy(&f));
            let change = original.energy_increment_to(&f, &changed, &f);
            near(change, expected, expected.abs(), 3e-11);
            near(
                changed.energy_increment_to(&f, &original, &f),
                -change,
                expected.abs(),
                3e-14,
            );
            assert_eq!(original.energy_increment_to(&f, &original, &f), 0.0);
        }
    }
}

#[test]
fn changed_nonlinear_shapes_keep_their_remainders_separate_from_common_bulk() {
    let gamma = 0.4_f64;
    let mut f = I.map(|x| 2.0 * x);
    f[1] = 2.0 * gamma;
    let volumetric = Volumetric::Quadratic { kappa: 1e24 };
    let gent = |locking_invariant| {
        Hyperelastic::new(
            HyperelasticModel::Gent {
                mu: 1024.0,
                locking_invariant,
                volumetric,
            },
            4.0,
        )
        .unwrap()
    };
    let ogden = |alpha| {
        Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: 1024.0,
                    alpha,
                }],
                volumetric,
            },
            4.0,
        )
        .unwrap()
    };
    let hgo = |direction| {
        Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 1024.0,
                k1: 512.0,
                k2: 1.7,
                dispersion: 0.0,
                fibres: vec![direction],
                invariants: FibreInvariants::Isochoric,
                volumetric,
            },
            4.0,
        )
        .unwrap()
    };
    let gent_energy = |limit: f64| -512.0 * limit * (-gamma * gamma / limit).ln_1p();
    let ogden_energy =
        |alpha: f64| 4096.0 / alpha * (0.5 * alpha * (0.5 * gamma).asinh()).sinh().powi(2);
    let cases = [
        (gent(1.0), gent(10.0), gent_energy(10.0) - gent_energy(1.0)),
        (
            ogden(1.3),
            ogden(2.4),
            ogden_energy(2.4) - ogden_energy(1.3),
        ),
        (
            hgo([0.0, 1.0, 0.0]),
            hgo([1.0, 0.0, 0.0]),
            -512.0 / 3.4 * (1.7 * gamma.powi(4)).exp_m1(),
        ),
        (
            gent(10.0),
            ogden(2.4),
            ogden_energy(2.4) - gent_energy(10.0),
        ),
    ];
    for (original, changed, expected) in cases {
        assert_eq!(original.energy(&f), changed.energy(&f));
        let change = original.energy_increment_to(&f, &changed, &f);
        near(change, expected, expected.abs(), 5e-12);
        near(
            changed.energy_increment_to(&f, &original, &f),
            -change,
            expected.abs(),
            5e-14,
        );
    }
}

#[test]
fn native_source_change_admits_only_each_sources_actual_endpoint() {
    let volumetric = Volumetric::Quadratic { kappa: 1e24 };
    let before = Hyperelastic::new(
        HyperelasticModel::Gent {
            mu: 1024.0,
            locking_invariant: 0.01,
            volumetric,
        },
        1.1,
    )
    .unwrap();
    let after = Hyperelastic::new(
        HyperelasticModel::Gent {
            mu: 2048.0,
            locking_invariant: 10.0,
            volumetric,
        },
        4.0,
    )
    .unwrap();
    let mut f = I;
    f[1] = 0.4;
    assert!(before.check_state(&f).is_err());
    let expected = -10240.0 * (-0.016_f64).ln_1p();
    near(
        before.energy_increment_to(&I, &after, &f),
        expected,
        expected,
        3e-13,
    );
    near(
        after.energy_increment_to(&f, &before, &I),
        -expected,
        expected,
        3e-13,
    );
    assert!(before.energy_increment_to(&f, &after, &f).is_nan());
    assert!(after.energy_increment_to(&f, &before, &f).is_nan());

    // A switched-off family must never evaluate its unused overflowing
    // exponential at the new source's otherwise admissible deformation.
    let fibre = |k1| {
        Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 1024.0,
                k1,
                k2: 1000.0,
                dispersion: 0.0,
                fibres: vec![[0.0, 1.0, 0.0]],
                invariants: FibreInvariants::Isochoric,
                volumetric,
            },
            4.0,
        )
        .unwrap()
    };
    let active = fibre(512.0);
    let inactive = fibre(0.0);
    f[1] = 3.0;
    assert!(!active.energy(&f).is_finite());
    let expected = 1024.0 * 9.0;
    near(
        active.energy_increment_to(&I, &inactive, &f),
        expected,
        expected,
        5e-14,
    );
    near(
        inactive.energy_increment_to(&f, &active, &I),
        -expected,
        expected,
        5e-14,
    );
}

#[test]
fn opposing_stiffness_and_strain_scales_keep_finite_native_source_work() {
    let law = |c10| {
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10,
                c01: 0.0,
                kappa: 60.0,
            },
            1e110,
        )
        .unwrap()
    };
    let before = law(1e200);
    let after = law(1e-140);
    let f = [1e70, 0.0, 0.0, 0.0, 1e-35, 0.0, 0.0, 0.0, 1e-35];
    before.check_state(&I).unwrap();
    after.check_state(&f).unwrap();
    assert_eq!(before.energy(&I), 0.0);
    near(after.energy(&f), 1.0, 1.0, 3e-14);
    near(before.energy_increment_to(&I, &after, &f), 1.0, 1.0, 3e-14);
    near(after.energy_increment_to(&f, &before, &I), -1.0, 1.0, 3e-14);
}

#[test]
fn changing_reference_retains_weak_shear_work_beneath_bulk_storage() {
    let law = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 900.0,
            c01: 600.0,
            kappa: 1e24,
        },
        4.0,
    )
    .unwrap();
    let before = ReferenceDeformation::try_new(I.map(|x| 2.0 * x)).unwrap();
    let mut gradient = *before.deformation_gradient();
    gradient[1] = 4e-7;
    let after = ReferenceDeformation::try_new(gradient).unwrap();
    let total = I.map(|x| 4.0 * x);
    let original = ReferencedHyperelastic::new(&law, before);
    let changed = ReferencedHyperelastic::new(&law, after);
    assert_eq!(
        changed.energy_density(&total) - original.energy_density(&total),
        0.0
    );
    // Fe_after = 2I - 4e-7 e_x (x) e_y, Jg=8, Je=8. The
    // original isochoric Mooney-Rivlin energy is (c10+c01)*(2e-7)^2.
    let expected = 8.0 * 1500.0 * (2e-7_f64).powi(2);
    near(
        original.energy_density_increment_to_reference(&total, &total, &after),
        expected,
        expected,
        3e-12,
    );
    near(
        changed.energy_density_increment_to_reference(&total, &total, &before),
        -expected,
        expected,
        3e-12,
    );
    let invalid = [0.0; 9];
    assert!(
        original
            .energy_density_increment_to_reference(&total, &invalid, &after)
            .is_nan()
    );
    assert_eq!(
        original.energy_density_increment_to_reference(&total, &total, &before),
        original.energy_density_increment(&total, &total)
    );
}

#[test]
fn reference_volume_change_retains_original_lengths_and_density_source() {
    let before = ReferenceDeformation::identity();
    let epsilon = 2.0_f64.powi(-30);
    let after =
        ReferenceDeformation::try_new([1.0, epsilon, 0.0, epsilon, 1.0, 0.0, 0.0, 0.0, 1.0])
            .unwrap();
    assert_eq!(after.volume_ratio(), before.volume_ratio());
    assert_eq!(before.volume_ratio_increment_to(&after), -epsilon * epsilon);
    assert_eq!(after.volume_ratio_increment_to(&before), epsilon * epsilon);
    let lengths = [1.0 + epsilon, 1.0 - epsilon, 1.0];
    let axes = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
    let principal = ReferenceDeformation::from_principal_stretches(lengths, axes).unwrap();
    let unrotated = ReferenceDeformation::from_principal_stretches(
        lengths,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    assert_eq!(principal.volume_ratio(), 1.0);
    assert_eq!(
        before.volume_ratio_increment_to(&principal),
        -epsilon * epsilon
    );
    assert_eq!(principal.volume_ratio_increment_to(&unrotated), 0.0);

    // The additive density belongs to this material, not its bare elastic
    // coefficients. It does work when the prescribed relaxed volume changes.
    struct OffsetEnergy<'a>(&'a Hyperelastic);
    impl HyperelasticMaterial for OffsetEnergy<'_> {
        fn hyperelastic(&self) -> &Hyperelastic {
            self.0
        }
        fn energy_density(&self, f: &[f64; 9]) -> f64 {
            2.0_f64.powi(60) + self.0.energy(f)
        }
        fn energy_density_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
            self.0.energy_increment(before, after)
        }
    }
    let law = law();
    let source = OffsetEnergy(&law);
    let material = ReferencedHyperelastic::new(&source, before);
    let changed = ReferencedHyperelastic::new(&source, after);
    assert_eq!(material.energy_density(&I), changed.energy_density(&I));
    let work = material.energy_density_increment_to_reference(&I, &I, &after);
    near(work, -1.0, 1.0, 2e-12);
    near(
        changed.energy_density_increment_to_reference(&I, &I, &before),
        -work,
        1.0,
        2e-15,
    );
}

#[test]
fn changed_reference_potential_matches_finite_original_volume_work() {
    let law = law();
    let before = ReferenceDeformation::try_new(G).unwrap();
    let after =
        ReferenceDeformation::try_new([0.98, -0.11, 0.03, 0.0, 1.21, 0.14, 0.02, 0.0, 1.06])
            .unwrap();
    let original = ReferencedHyperelastic::new(&law, before);
    let changed = ReferencedHyperelastic::new(&law, after);
    let f0 = multiply(&FE, &G);
    let f1 = multiply(&FE, after.deformation_gradient());
    for end in [&f0, &f1] {
        let expected = changed.energy_density(end) - original.energy_density(&f0);
        let change = original.energy_density_increment_to_reference(&f0, end, &after);
        near(change, expected, 1e6, 3e-14);
        near(
            changed.energy_density_increment_to_reference(end, &f0, &before),
            -change,
            1e6,
            3e-14,
        );
    }
}

#[test]
fn changing_reference_avoids_overflow_between_finite_original_energies() {
    struct Prestressed<'a>(&'a Hyperelastic);
    impl HyperelasticMaterial for Prestressed<'_> {
        fn hyperelastic(&self) -> &Hyperelastic {
            self.0
        }
        fn energy_density(&self, f: &[f64; 9]) -> f64 {
            self.0.energy(f) + 1e300 * f[1]
        }
        fn energy_density_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
            self.0.energy_increment(before, after) + 1e300 * (after[1] - before[1])
        }
        fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], fs_material::MaterialError> {
            let mut p = self.0.piola(f)?;
            p[1] += 1e300;
            Ok(p)
        }
    }
    let law = law();
    let source = Prestressed(&law);
    let before = ReferenceDeformation::try_new(I.map(|x| 1e100 * x)).unwrap();
    let after = ReferenceDeformation::try_new(I.map(|x| 1e-100 * x)).unwrap();
    let material = ReferencedHyperelastic::new(&source, before);
    let changed = ReferencedHyperelastic::new(&source, after);
    let f0 = *before.deformation_gradient();
    let mut f1 = *after.deformation_gradient();
    f1[1] = 1e-100;
    assert_eq!(material.energy_density(&f0), 0.0);
    let expected = changed.energy_density(&f1);
    near(expected, 1.0, 1.0, 2e-15);
    near(
        material.energy_density_increment_to_reference(&f0, &f1, &after),
        expected,
        1.0,
        2e-15,
    );
    near(
        changed.energy_density_increment_to_reference(&f1, &f0, &before),
        -expected,
        1.0,
        2e-15,
    );
}

#[test]
fn reference_energy_changes_preserve_original_volume_and_weak_work() {
    let law = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 900.0,
            c01: 600.0,
            kappa: 1e24,
        },
        4.0,
    )
    .unwrap();
    let reference = ReferenceDeformation::try_new(I.map(|x| 2.0 * x)).unwrap();
    let material = ReferencedHyperelastic::new(&law, reference);
    let before = I.map(|x| 4.0 * x);
    let mut after = before;
    after[1] = 8e-7;
    assert_eq!(
        material.energy_density(&after) - material.energy_density(&before),
        0.0
    );
    let increment = material.energy_density_increment(&before, &after);
    let expected = 8.0 * 1500.0 * (2e-7_f64).powi(2);
    near(increment, expected, expected, 3e-12);
    near(
        material.energy_density_increment(&after, &before),
        -expected,
        expected,
        3e-12,
    );
    let unchanged = ReferencedHyperelastic::new(&law, ReferenceDeformation::identity());
    assert_eq!(
        unchanged.energy_density_increment(&I, &FE),
        law.energy_increment(&I, &FE)
    );
    let invalid = I.map(|x| 0.1 * x);
    assert!(
        material
            .energy_density_increment(&invalid, &invalid)
            .is_nan()
    );
}

#[test]
fn default_energy_change_and_tangent_action_keep_an_adapters_own_potential() {
    struct AffineEnergy<'a>(&'a Hyperelastic);
    impl HyperelasticMaterial for AffineEnergy<'_> {
        fn hyperelastic(&self) -> &Hyperelastic {
            self.0
        }
        fn energy_density(&self, f: &[f64; 9]) -> f64 {
            11.0 + 2.0 * self.0.energy(f) + 7.0 * f[1]
        }
        fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], fs_material::MaterialError> {
            let mut p = self.0.piola(f)?.map(|x| 2.0 * x);
            p[1] += 7.0;
            Ok(p)
        }
        fn consistent_tangent(
            &self,
            f: &[f64; 9],
        ) -> Result<[[f64; 9]; 9], fs_material::MaterialError> {
            Ok(self.0.tangent(f)?.map(|row| row.map(|x| 2.0 * x)))
        }
    }
    let law = law();
    let material = AffineEnergy(&law);
    let expected = material.energy_density(&FE) - material.energy_density(&I);
    assert_eq!(material.energy_density_increment(&I, &FE), expected);
    assert_eq!(
        fs_material::hyper::cards::FiniteStrainPotential::energy_density_increment(
            &material, &I, &FE
        ),
        expected
    );
    assert!((expected - law.energy_increment(&I, &FE)).abs() > 1.0);
    let direction = [0.3, -0.1, 0.2, 0.04, 0.17, -0.12, 0.07, -0.08, 0.11];
    let expected_action = material
        .consistent_tangent(&FE)
        .unwrap()
        .map(|row| (0..9).fold(0.0, |sum, j| row[j].mul_add(direction[j], sum)));
    assert_eq!(
        material.consistent_tangent_action(&FE, &direction).unwrap(),
        expected_action
    );
    assert_eq!(
        fs_material::hyper::cards::FiniteStrainPotential::consistent_tangent_action(
            &material, &FE, &direction
        )
        .unwrap(),
        expected_action
    );
    let reference = ReferenceDeformation::try_new(G).unwrap();
    let held = ReferencedHyperelastic::new(&material, reference);
    let before = multiply(&I, &G);
    let after = multiply(&FE, &G);
    let original = held.energy_density(&after) - held.energy_density(&before);
    near(
        held.energy_density_increment(&before, &after),
        original,
        original.abs(),
        2e-13,
    );
    let expected_action = held
        .consistent_tangent(&after)
        .unwrap()
        .map(|row| (0..9).fold(0.0, |sum, j| row[j].mul_add(direction[j], sum)));
    let action = held.consistent_tangent_action(&after, &direction).unwrap();
    for i in 0..9 {
        near(action[i], expected_action[i], 1e6, 3e-13);
    }
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

fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}

fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 2.3e5,
            c01: 0.7e5,
            kappa: 3.1e6,
        },
        3.0,
    )
    .unwrap()
}

#[test]
fn reference_energy_gradient_and_tangent_obey_original_volume_chain_rule() {
    let law = law();
    let reference = ReferenceDeformation::try_new(G).unwrap();
    let material = ReferencedHyperelastic::new(&law, reference);
    let f = multiply(&FE, &G);
    let energy = material.energy_density(&f);
    near(energy, determinant(&G) * law.energy(&FE), 1.0, 1e-13);
    let stress = material.first_piola(&f).unwrap();
    let tangent = material.consistent_tangent(&f).unwrap();
    let h = 2e-6;
    for j in 0..9 {
        let (mut plus, mut minus) = (f, f);
        plus[j] += h;
        minus[j] -= h;
        near(
            stress[j],
            (material.energy_density(&plus) - material.energy_density(&minus)) / (2.0 * h),
            1e6,
            2e-8,
        );
        let pp = material.first_piola(&plus).unwrap();
        let pm = material.first_piola(&minus).unwrap();
        for i in 0..9 {
            near(tangent[i][j], (pp[i] - pm[i]) / (2.0 * h), 1e6, 2e-8);
            near(tangent[i][j], tangent[j][i], 1e6, 2e-14);
        }
    }
}

#[test]
fn reference_directional_tangent_preserves_shear_and_sub_ulp_increments() {
    let law = law();
    let direction = [0.3, -0.1, 0.2, 0.04, 0.17, -0.12, 0.07, -0.08, 0.11];
    for g in [I, G, I.map(|value| 1.25 * value)] {
        let reference = ReferenceDeformation::try_new(g).unwrap();
        let material = ReferencedHyperelastic::new(&law, reference);
        let f = multiply(&FE, &g);
        let tangent = material.consistent_tangent(&f).unwrap();
        for factor in [1e-24, 1.0, -4.0] {
            let increment = direction.map(|value| factor * value);
            let elastic_increment = reference.elastic_deformation_increment(&increment);
            for (actual, expected) in multiply(&elastic_increment, &g).into_iter().zip(increment) {
                near(actual, expected, factor.abs(), 1e-15);
            }
            let action = material.consistent_tangent_action(&f, &increment).unwrap();
            let expected: [f64; 9] =
                tangent.map(|row| (0..9).fold(0.0, |sum, j| row[j].mul_add(increment[j], sum)));
            let scale = expected.iter().map(|x| x.abs()).fold(0.0, f64::max);
            for i in 0..9 {
                near(action[i], expected[i], scale, 3e-13);
            }
        }
        let h = 2e-6;
        let plus = core::array::from_fn(|i| f[i] + h * direction[i]);
        let minus = core::array::from_fn(|i| f[i] - h * direction[i]);
        let pp = material.first_piola(&plus).unwrap();
        let pm = material.first_piola(&minus).unwrap();
        let action = material.consistent_tangent_action(&f, &direction).unwrap();
        for i in 0..9 {
            near(action[i], (pp[i] - pm[i]) / (2.0 * h), 1e6, 2e-8);
        }
        assert!(
            material
                .consistent_tangent_action(&f, &[f64::NAN; 9])
                .is_err()
        );
        assert!(
            material
                .consistent_tangent_action(&[0.0; 9], &[0.0; 9])
                .is_err()
        );
    }
}

#[test]
fn reference_relaxes_prescribed_shear_and_retains_resolved_weak_strain() {
    let reference = ReferenceDeformation::try_new(G).unwrap();
    assert_eq!(reference.elastic_deformation(&G), I);
    let law = law();
    let material = ReferencedHyperelastic::new(&law, reference);
    assert_eq!(material.energy_density(&G), 0.0);
    assert!(
        material
            .first_piola(&G)
            .unwrap()
            .iter()
            .all(|x| x.abs() < 1e-9)
    );
    // Perturb an originally zero component so the input increment is exactly
    // represented despite the finite, sheared reference configuration.
    let mut f = G;
    f[6] = 1e-12;
    let fe = reference.elastic_deformation(&f);
    near(fe[6], 1e-12 / G[0], 1e-12, 5e-15);
    let stress = material.first_piola(&f).unwrap();
    let tangent = material.consistent_tangent(&G).unwrap();
    // A resolved shear component must not be swallowed by roundoff in G G^-1.
    near(stress[6], tangent[6][6] * 1e-12, 1e-6, 5e-4);
    assert!(stress[6].abs() > 1e-7);
}

#[test]
fn reference_preserves_spatial_objectivity_for_sheared_material_axes() {
    let law = law();
    let material = ReferencedHyperelastic::new(&law, ReferenceDeformation::try_new(G).unwrap());
    let rotation = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let f = multiply(&FE, &G);
    let qf = multiply(&rotation, &f);
    near(
        material.energy_density(&f),
        material.energy_density(&qf),
        1e6,
        1e-14,
    );
    let p = material.first_piola(&f).unwrap();
    let qp = material.first_piola(&qf).unwrap();
    for (a, b) in qp.into_iter().zip(multiply(&rotation, &p)) {
        near(a, b, 1e6, 1e-14);
    }
    let direction = [0.3, -0.1, 0.2, 0.04, 0.17, -0.12, 0.07, -0.08, 0.11];
    let tangent = material.consistent_tangent(&f).unwrap();
    let qtangent = material.consistent_tangent(&qf).unwrap();
    let qdirection = multiply(&rotation, &direction);
    let action = core::array::from_fn(|i| (0..9).map(|j| tangent[i][j] * direction[j]).sum());
    let qaction: [f64; 9] =
        core::array::from_fn(|i| (0..9).map(|j| qtangent[i][j] * qdirection[j]).sum());
    for (a, b) in qaction.into_iter().zip(multiply(&rotation, &action)) {
        near(a, b, 1e6, 1e-13);
    }
}

#[test]
fn reference_mixed_blocks_enforce_elastic_volume_and_condense_to_original_law() {
    let law = law();
    let reference = ReferenceDeformation::try_new(G).unwrap();
    let material = ReferencedHyperelastic::new(&law, reference);
    let f = multiply(&FE, &G);
    let p = 3.1e6 * (determinant(&FE) - 1.0);
    let response = material
        .mixed_response(&f, p, VolumetricTreatment::OriginalMaterial)
        .unwrap();
    near(response.state.pressure_residual, 0.0, 1.0, 1e-14);
    near(
        response.state.potential_pa,
        material.energy_density(&f),
        1e6,
        1e-14,
    );
    let full_stress = material.first_piola(&f).unwrap();
    let full_tangent = material.consistent_tangent(&f).unwrap();
    for i in 0..9 {
        near(response.state.first_piola_pa[i], full_stress[i], 1e6, 1e-14);
        for j in 0..9 {
            let condensed = response.tangent_pa[i][j]
                - response.pressure_coupling[i] * response.pressure_coupling[j]
                    / response.pressure_tangent_per_pa;
            near(condensed, full_tangent[i][j], 1e6, 2e-14);
        }
    }
    let relaxed = material
        .mixed_state(&G, 0.0, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    near(relaxed.pressure_residual, 0.0, 1.0, 1e-15);
    let constrained = material
        .mixed_state(&I, 0.0, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    near(
        constrained.pressure_residual,
        1.0 - determinant(&G),
        1.0,
        1e-15,
    );
    for treatment in [
        VolumetricTreatment::OriginalMaterial,
        VolumetricTreatment::ExactIncompressible,
    ] {
        let r = material.mixed_response(&f, 8.0e4, treatment).unwrap();
        let h = 2e-6;
        for j in 0..9 {
            let (mut plus, mut minus) = (f, f);
            plus[j] += h;
            minus[j] -= h;
            let a = material.mixed_state(&plus, 8.0e4, treatment).unwrap();
            let b = material.mixed_state(&minus, 8.0e4, treatment).unwrap();
            near(
                (a.potential_pa - b.potential_pa) / (2.0 * h),
                r.state.first_piola_pa[j],
                1e6,
                2e-8,
            );
            near(
                (a.pressure_residual - b.pressure_residual) / (2.0 * h),
                r.pressure_coupling[j],
                1.0,
                2e-9,
            );
            for i in 0..9 {
                near(
                    (a.first_piola_pa[i] - b.first_piola_pa[i]) / (2.0 * h),
                    r.tangent_pa[i][j],
                    1e6,
                    2e-8,
                );
            }
        }
        let a = material.mixed_state(&f, 8.0e4 + 1.0, treatment).unwrap();
        let b = material.mixed_state(&f, 8.0e4 - 1.0, treatment).unwrap();
        near(
            (a.potential_pa - b.potential_pa) / 2.0,
            r.state.pressure_residual,
            1.0,
            1e-10,
        );
        near(
            (a.pressure_residual - b.pressure_residual) / 2.0,
            r.pressure_tangent_per_pa,
            1e-6,
            1e-10,
        );
        for i in 0..9 {
            near(
                (a.first_piola_pa[i] - b.first_piola_pa[i]) / 2.0,
                r.pressure_coupling[i],
                1.0,
                1e-10,
            );
        }
    }
}

#[test]
fn reference_constitutive_domains_follow_elastic_stretches_and_volume() {
    let law = Hyperelastic::new(
        HyperelasticModel::Gent {
            mu: 2e5,
            locking_invariant: 4.0,
            volumetric: Volumetric::LogSquared { kappa: 8e5 },
        },
        1.5,
    )
    .unwrap();
    let reference =
        ReferenceDeformation::try_new([2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0]).unwrap();
    let material = ReferencedHyperelastic::new(&law, reference);
    let g = reference.deformation_gradient();
    assert!(law.piola(g).is_err());
    assert!(material.first_piola(g).is_ok());
    assert!(
        material
            .mixed_state(g, 0.0, VolumetricTreatment::OriginalMaterial)
            .is_ok()
    );
    let mut invalid = *g;
    invalid[0] = 3.01;
    assert!(material.first_piola(&invalid).is_err());
    assert!(
        material
            .mixed_state(&invalid, 0.0, VolumetricTreatment::OriginalMaterial)
            .is_err()
    );
    let zero = [0.0; 9];
    assert!(material.first_piola(&zero).is_err());
}

#[test]
fn oriented_principal_reference_retains_lengths_and_exact_relaxed_configuration() {
    let q = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
    let s = [1.17, 0.91, 1.06];
    let reference = ReferenceDeformation::from_principal_stretches(s, q).unwrap();
    let g = reference.deformation_gradient();
    for a in 0..3 {
        for i in 0..3 {
            let actual: f64 = (0..3).map(|j| g[3 * i + j] * q[j][a]).sum();
            near(actual, s[a] * q[i][a], 1.0, 1e-15);
        }
    }
    near(reference.volume_ratio(), s[0] * s[1] * s[2], 1.0, 1e-15);
    assert_eq!(reference.elastic_deformation(g), I);
    let isotropic = ReferenceDeformation::from_principal_stretches([1.2; 3], q).unwrap();
    assert_eq!(
        isotropic.deformation_gradient(),
        &[1.2, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, 1.2]
    );
    assert_eq!(
        isotropic.elastic_deformation(isotropic.deformation_gradient()),
        I
    );
    assert!(
        ReferenceDeformation::from_principal_stretches(
            s,
            [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        )
        .is_err()
    );
    assert!(
        ReferenceDeformation::from_principal_stretches(
            s,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]
        )
        .is_err()
    );
}

#[test]
fn identity_reference_is_bitwise_compatible_and_invalid_maps_refuse() {
    let law = law();
    let reference = ReferenceDeformation::identity();
    let material = ReferencedHyperelastic::new(&law, reference);
    assert_eq!(
        material.energy_density(&FE).to_bits(),
        law.energy(&FE).to_bits()
    );
    assert_eq!(material.first_piola(&FE).unwrap(), law.piola(&FE).unwrap());
    assert_eq!(
        material.consistent_tangent(&FE).unwrap(),
        law.tangent(&FE).unwrap()
    );
    let original = law
        .mixed_response(&FE, 8e4, VolumetricTreatment::OriginalMaterial)
        .unwrap();
    let mapped = material
        .mixed_response(&FE, 8e4, VolumetricTreatment::OriginalMaterial)
        .unwrap();
    assert_eq!(original.tangent_pa, mapped.tangent_pa);
    assert_eq!(
        original.state.potential_pa.to_bits(),
        mapped.state.potential_pa.to_bits()
    );
    assert_eq!(original.state.first_piola_pa, mapped.state.first_piola_pa);
    for bad in [
        [0.0; 9],
        [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        [f64::NAN; 9],
        [f64::INFINITY; 9],
    ] {
        assert!(ReferenceDeformation::try_new(bad).is_err());
    }
}
