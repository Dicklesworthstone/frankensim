//! Fixed natural configurations in the original finite-strain Maxwell owner.
//! Independent matrix contractions and finite differences use total F, while
//! an unmapped history supplies the same elastic loading/relaxation trajectory.

use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
// A nonsymmetric map discriminates the inverse from its transpose.
const G: [f64; 9] = [1.25, 0.18, 0.0, 0.03, 0.85, 0.12, 0.0, -0.04, 1.1];
const PREVIOUS: [f64; 9] = [1.18, 0.1, 0.0, 0.02, 0.94, 0.0, 0.0, 0.0, 0.93];
const TRIAL: [f64; 9] = [1.12, 0.06, 0.02, -0.02, 0.98, 0.03, 0.01, 0.0, 0.95];

fn determinant(a: &[f64; 9]) -> f64 {
    a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6])
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn inverse(a: &[f64; 9]) -> [f64; 9] {
    let j = determinant(a);
    [
        a[4] * a[8] - a[5] * a[7],
        a[2] * a[7] - a[1] * a[8],
        a[1] * a[5] - a[2] * a[4],
        a[5] * a[6] - a[3] * a[8],
        a[0] * a[8] - a[2] * a[6],
        a[2] * a[3] - a[0] * a[5],
        a[3] * a[7] - a[4] * a[6],
        a[1] * a[6] - a[0] * a[7],
        a[0] * a[4] - a[1] * a[3],
    ]
    .map(|x| x / j)
}

fn pullback_stress(p: &[f64; 9]) -> [f64; 9] {
    let h = inverse(&G);
    let j = determinant(&G);
    core::array::from_fn(|ij| {
        j * (0..3)
            .map(|a| p[3 * (ij / 3) + a] * h[3 * (ij % 3) + a])
            .sum::<f64>()
    })
}

fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite(), "nonfinite {a} or {b}");
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} vs {b:e}, scale {scale:e}, tolerance {tolerance:e}"
    );
}

fn law(log_bulk: bool) -> ViscoHyperelastic {
    let model = if log_bulk {
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 1.1e5,
                alpha: 1.7,
            }],
            volumetric: Volumetric::LogSquared { kappa: 8e5 },
        }
    } else {
        HyperelasticModel::MooneyRivlin {
            c10: 3e4,
            c01: 2e4,
            kappa: 8e5,
        }
    };
    ViscoHyperelastic::new(
        Hyperelastic::new(model, 2.0).unwrap(),
        vec![
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3e4 },
                tau_s: 0.2,
            },
            FiniteStrainBranch {
                energy: BranchEnergy::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 2e4,
                        alpha: 1.6,
                    }],
                },
                tau_s: 0.07,
            },
        ],
    )
    .unwrap()
}

#[test]
fn natural_configuration_is_relaxed_and_limits_apply_to_elastic_stretch() {
    let natural = [4.0, 0.5, 0.0, 0.0, 0.75, 0.25, 0.0, 0.0, 1.25];
    let base = law(false);
    assert!(base.equilibrium.check_state(&natural).is_err());
    let mapped = base.with_reference_deformation(ReferenceDeformation::try_new(natural).unwrap());
    let initial = mapped.initial_state();
    let remainder = mapped.initial_remainder_state();
    assert_eq!(initial.deformation_gradient, natural);
    assert_eq!(remainder.deformation_gradient, natural);
    for dt in [0.0, 0.2] {
        let response = mapped.step_response(&initial, &natural, dt).unwrap();
        let rem = mapped
            .step_remainder_response(&remainder, &natural, dt)
            .unwrap();
        for value in response.first_piola_pa.iter().chain(&rem.first_piola_pa) {
            near(*value, 0.0, 1e6, 1e-13);
        }
        assert_eq!(
            response.state.inverse_viscous_cauchy_green,
            initial.inverse_viscous_cauchy_green
        );
        assert_eq!(
            rem.state.inverse_viscous_cauchy_green,
            remainder.inverse_viscous_cauchy_green
        );
        assert_eq!(response.state.stored, 0.0);
        assert_eq!(response.state.work, 0.0);
        assert_eq!(response.dissipated_pa, 0.0);
        assert_eq!(rem.state.stored, 0.0);
        assert_eq!(rem.dissipated_pa, 0.0);
    }
    let frozen = mapped.freeze_remainder(&remainder, 0.1).unwrap();
    let exact = frozen
        .mixed_response(&natural, 0.0, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    near(exact.state.pressure_residual, 0.0, 1.0, 2e-14);
    near(
        exact.state.bulk_volume_ratio,
        determinant(&natural),
        10.0,
        2e-14,
    );
    assert_eq!(exact.pressure_tangent_per_pa, 0.0);
    assert_eq!(exact.state.bulk_stored_pa, 0.0);

    let rotation = [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let rotated = multiply(&rotation, &natural);
    let response = mapped.step_response(&initial, &rotated, 0.2).unwrap();
    for p in response.first_piola_pa {
        near(p, 0.0, 1e6, 1e-12);
    }
    near(response.state.stored, 0.0, 1e6, 1e-12);
    near(response.dissipated_pa, 0.0, 1e6, 1e-12);

    let inadmissible = multiply(&[2.1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], &natural);
    assert!(mapped.step_response(&initial, &inadmissible, 0.1).is_err());
    assert!(
        mapped
            .step_remainder_response(&remainder, &inadmissible, 0.1)
            .is_err()
    );
    assert!(mapped.free_energy(&initial, &inadmissible).is_err());
}

#[test]
fn noncommuting_reference_preserves_relaxation_and_original_volume_work() {
    let base = law(false);
    let mapped = base
        .clone()
        .with_reference_deformation(ReferenceDeformation::try_new(G).unwrap());
    let jg = determinant(&G);
    let (mut plain, mut pulled) = (base.initial_state(), mapped.initial_state());
    let (mut plain_rem, mut pulled_rem) = (
        base.initial_remainder_state(),
        mapped.initial_remainder_state(),
    );
    for (elastic, dt) in [
        (PREVIOUS, 0.0),
        (PREVIOUS, 0.07),
        (TRIAL, 0.08),
        (TRIAL, 0.12),
        (I, 0.0),
    ] {
        let total = multiply(&elastic, &G);
        let a = base.step_response(&plain, &elastic, dt).unwrap();
        let b = mapped.step_response(&pulled, &total, dt).unwrap();
        let ar = base
            .step_remainder_response(&plain_rem, &elastic, dt)
            .unwrap();
        let br = mapped
            .step_remainder_response(&pulled_rem, &total, dt)
            .unwrap();
        assert_eq!(b.state.deformation_gradient, total);
        assert_eq!(br.state.deformation_gradient, total);
        for (actual, expected) in b
            .first_piola_pa
            .iter()
            .zip(pullback_stress(&a.first_piola_pa))
        {
            near(*actual, expected, 1e6, 2e-11);
        }
        for (actual, expected) in br
            .first_piola_pa
            .iter()
            .zip(pullback_stress(&ar.first_piola_pa))
        {
            near(*actual, expected, 1e6, 2e-11);
        }
        for (actual, expected) in [
            (b.state.work, a.state.work),
            (b.state.stored, a.state.stored),
            (b.state.dissipated, a.state.dissipated),
            (b.loading_work_pa, a.loading_work_pa),
            (b.incremental_potential_pa, a.incremental_potential_pa),
            (b.dissipated_pa, a.dissipated_pa),
            (br.state.work, ar.state.work),
            (br.state.stored, ar.state.stored),
            (br.state.dissipated, ar.state.dissipated),
            (br.loading_work_pa, ar.loading_work_pa),
            (br.incremental_potential_pa, ar.incremental_potential_pa),
            (br.dissipated_pa, ar.dissipated_pa),
        ] {
            near(actual, jg * expected, 1e6, 2e-11);
        }
        for (actual, expected) in b
            .state
            .inverse_viscous_cauchy_green
            .iter()
            .flatten()
            .zip(a.state.inverse_viscous_cauchy_green.iter().flatten())
        {
            near(*actual, *expected, 1.0, 2e-12);
        }
        for (actual, expected) in br
            .state
            .inverse_viscous_cauchy_green
            .iter()
            .flatten()
            .zip(b.state.inverse_viscous_cauchy_green.iter().flatten())
        {
            near(*actual, *expected, 1.0, 2e-13);
        }
        near(
            mapped.free_energy(&b.state, &total).unwrap(),
            b.state.stored,
            1e6,
            2e-11,
        );
        near(
            mapped.remainder_free_energy(&br.state, &total).unwrap(),
            br.state.stored,
            1e6,
            2e-11,
        );
        near(
            mapped
                .mixed_remainder_free_energy(&b.state, &total)
                .unwrap(),
            br.state.stored,
            1e6,
            2e-11,
        );
        near(
            b.state.work - b.state.stored - b.state.dissipated,
            0.0,
            1e6,
            2e-11,
        );
        near(
            br.state.work - br.state.stored - br.state.dissipated,
            0.0,
            1e6,
            2e-11,
        );
        assert!(b.dissipated_pa >= 0.0 && br.dissipated_pa >= 0.0);
        if dt == 0.0 {
            assert_eq!(b.dissipated_pa, 0.0);
        }
        plain = a.state;
        pulled = b.state;
        plain_rem = ar.state;
        pulled_rem = br.state;
    }
    assert!(pulled.dissipated > 0.0);
}

#[test]
fn mapped_incremental_gradient_and_hessian_are_conjugate_to_total_deformation() {
    let mapped = law(false).with_reference_deformation(ReferenceDeformation::try_new(G).unwrap());
    let accepted = mapped
        .step_response(&mapped.initial_state(), &multiply(&PREVIOUS, &G), 0.03)
        .unwrap()
        .state;
    let snapshot = accepted.clone();
    let frozen = mapped.freeze(&accepted, 0.08).unwrap();
    let f = multiply(&TRIAL, &G);
    let (p, potential) = frozen.stress_and_potential(&f).unwrap();
    let a = frozen.consistent_tangent(&f).unwrap();
    near(
        potential,
        mapped
            .step_response(&accepted, &f, 0.08)
            .unwrap()
            .incremental_potential_pa,
        1e6,
        1e-14,
    );
    for col in 0..9 {
        let h = 2e-6;
        let (mut plus, mut minus) = (f, f);
        plus[col] += h;
        minus[col] -= h;
        let (pp, wp) = frozen.stress_and_potential(&plus).unwrap();
        let (pm, wm) = frozen.stress_and_potential(&minus).unwrap();
        near((wp - wm) / (2.0 * h), p[col], 1e6, 2e-7);
        for row in 0..9 {
            near((pp[row] - pm[row]) / (2.0 * h), a[row][col], 1e6, 2e-7);
            near(a[row][col], a[col][row], 1e6, 2e-12);
        }
    }
    assert_eq!(accepted, snapshot);
}

fn mixed_derivatives(
    law: &impl MixedFiniteStrainPotential,
    f: &[f64; 9],
    treatment: VolumetricTreatment,
) {
    let p = 35000.0;
    let response = law.mixed_response(f, p, treatment).unwrap();
    for col in 0..9 {
        let h = 2e-6;
        let (mut fp, mut fm) = (*f, *f);
        fp[col] += h;
        fm[col] -= h;
        let plus = law.mixed_state(&fp, p, treatment).unwrap();
        let minus = law.mixed_state(&fm, p, treatment).unwrap();
        near(
            (plus.potential_pa - minus.potential_pa) / (2.0 * h),
            response.state.first_piola_pa[col],
            1e6,
            2e-7,
        );
        near(
            (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
            response.pressure_coupling[col],
            1.0,
            2e-8,
        );
        for row in 0..9 {
            near(
                (plus.first_piola_pa[row] - minus.first_piola_pa[row]) / (2.0 * h),
                response.tangent_pa[row][col],
                1e6,
                2e-7,
            );
            near(
                response.tangent_pa[row][col],
                response.tangent_pa[col][row],
                1e6,
                2e-12,
            );
        }
    }
    let h = 0.2;
    let plus = law.mixed_state(f, p + h, treatment).unwrap();
    let minus = law.mixed_state(f, p - h, treatment).unwrap();
    near(
        (plus.potential_pa - minus.potential_pa) / (2.0 * h),
        response.state.pressure_residual,
        1.0,
        2e-8,
    );
    near(
        (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
        response.pressure_tangent_per_pa,
        1e-6,
        2e-8,
    );
    for i in 0..9 {
        near(
            (plus.first_piola_pa[i] - minus.first_piola_pa[i]) / (2.0 * h),
            response.pressure_coupling[i],
            1.0,
            2e-8,
        );
    }
}

#[test]
fn mapped_mixed_constraint_and_pressure_elimination_use_elastic_volume() {
    for log_bulk in [false, true] {
        let mapped =
            law(log_bulk).with_reference_deformation(ReferenceDeformation::try_new(G).unwrap());
        let previous = multiply(&PREVIOUS, &G);
        let accepted = mapped
            .step_response(&mapped.initial_state(), &previous, 0.03)
            .unwrap()
            .state;
        let accepted_rem = mapped
            .step_remainder_response(&mapped.initial_remainder_state(), &previous, 0.03)
            .unwrap()
            .state;
        let full = mapped.freeze(&accepted, 0.08).unwrap();
        let remainder = mapped.freeze_remainder(&accepted_rem, 0.08).unwrap();
        let f = multiply(&TRIAL, &G);
        let je = determinant(&TRIAL);
        let jg = determinant(&G);
        let pressure = if log_bulk {
            8e5 * je.ln() / je
        } else {
            8e5 * (je - 1.0)
        };
        let mixed = full
            .mixed_response(&f, pressure, VolumetricTreatment::OriginalMaterial)
            .unwrap();
        let rem = remainder
            .mixed_response(&f, pressure, VolumetricTreatment::OriginalMaterial)
            .unwrap();
        let (p, value) = full.stress_and_potential(&f).unwrap();
        let tangent = full.consistent_tangent(&f).unwrap();
        near(mixed.state.pressure_residual, 0.0, 1.0, 2e-13);
        near(mixed.state.bulk_volume_ratio, jg * je, 1.0, 2e-13);
        near(mixed.state.potential_pa, value, 1e6, 2e-12);
        near(rem.state.potential_pa, value, 1e6, 2e-12);
        for i in 0..9 {
            near(mixed.state.first_piola_pa[i], p[i], 1e6, 2e-12);
            near(rem.state.first_piola_pa[i], p[i], 1e6, 2e-12);
            for j in 0..9 {
                let condensed = mixed.tangent_pa[i][j]
                    - mixed.pressure_coupling[i] * mixed.pressure_coupling[j]
                        / mixed.pressure_tangent_per_pa;
                near(condensed, tangent[i][j], 1e6, 3e-12);
            }
        }
        let exact = remainder
            .mixed_response(&f, 35000.0, VolumetricTreatment::ExactIncompressible)
            .unwrap();
        near(exact.state.pressure_residual, jg * (je - 1.0), 1.0, 2e-13);
        near(exact.state.bulk_volume_ratio, jg, 1.0, 2e-14);
        assert_eq!(exact.state.bulk_stored_pa, 0.0);
        assert_eq!(exact.pressure_tangent_per_pa, 0.0);
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            mixed_derivatives(&full, &f, treatment);
            mixed_derivatives(&remainder, &f, treatment);
        }
    }
}

#[test]
fn a_fixed_reference_cannot_be_changed_under_retained_history() {
    let mapped = law(false).with_reference_deformation(ReferenceDeformation::try_new(G).unwrap());
    let total = multiply(&PREVIOUS, &G);
    let state = mapped
        .step_response(&mapped.initial_state(), &total, 0.03)
        .unwrap()
        .state;
    let remainder = mapped
        .step_remainder_response(&mapped.initial_remainder_state(), &total, 0.03)
        .unwrap()
        .state;
    let snapshot = (state.clone(), remainder.clone());
    let other = mapped
        .clone()
        .with_reference_deformation(ReferenceDeformation::identity());
    assert!(matches!(
        other.step_response(&state, &total, 0.04),
        Err(MaterialError::State { .. })
    ));
    assert!(matches!(
        other.step_remainder_response(&remainder, &total, 0.04),
        Err(MaterialError::State { .. })
    ));
    assert!(other.freeze(&state, 0.04).is_err());
    assert!(other.freeze_remainder(&remainder, 0.04).is_err());
    assert!(other.free_energy(&state, &total).is_err());
    assert!(other.mixed_remainder_free_energy(&state, &total).is_err());
    assert!(other.remainder_free_energy(&remainder, &total).is_err());
    assert_eq!((state, remainder), snapshot);
}
