//! Remainder histories advance the original metrics without any bulk evaluation.

use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const F: [f64; 9] = [1.2, 0.1, 0.0, 0.02, 0.94, 0.0, 0.0, 0.0, 0.91];
const G: [f64; 9] = [1.16, 0.07, 0.02, -0.01, 0.97, 0.03, 0.01, 0.0, 0.93];

fn law(kappa: f64) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6000.0,
                c01: 4000.0,
                kappa,
            },
            4.0,
        )
        .unwrap(),
        vec![
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 30_000.0 },
                tau_s: 0.2,
            },
            FiniteStrainBranch {
                energy: BranchEnergy::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 18_000.0,
                        alpha: 1.7,
                    }],
                },
                tau_s: 0.07,
            },
        ],
    )
    .unwrap()
}
fn det(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}
fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= scale * tolerance,
        "{a:e} vs {b:e}, scale {scale:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g0_remainder_history_matches_original_branch_evolution_and_keeps_its_own_ledger() {
    let law = law(8e5);
    let mut original = law.initial_state();
    let mut remainder = law.initial_remainder_state();
    for (f, dt) in [(F, 0.0), (F, 0.04), (G, 0.03), (G, 0.2)] {
        let old = remainder.clone();
        let full = law.step_response(&original, &f, dt).unwrap();
        let reduced = law.step_remainder_response(&remainder, &f, dt).unwrap();
        assert_eq!(
            full.state.inverse_viscous_cauchy_green,
            reduced.state.inverse_viscous_cauchy_green
        );
        assert_eq!(
            full.state.deformation_gradient,
            reduced.state.deformation_gradient
        );
        assert_eq!(
            full.dissipated_pa.to_bits(),
            reduced.dissipated_pa.to_bits()
        );
        assert_eq!(
            full.state.dissipated.to_bits(),
            reduced.state.dissipated.to_bits()
        );
        let bulk = law.equilibrium.volumetric_energy(det(&f));
        near(full.state.stored - reduced.state.stored, bulk, 1e5, 2e-12);
        near(full.state.work - reduced.state.work, bulk, 1e5, 2e-12);
        near(
            reduced.loading_work_pa - (reduced.state.stored - old.stored) - reduced.dissipated_pa,
            0.0,
            1e5,
            2e-11,
        );
        near(
            law.remainder_free_energy(&reduced.state, &f).unwrap(),
            reduced.state.stored,
            1e5,
            2e-11,
        );
        assert_eq!(remainder, old);
        let frozen_full = law.freeze(&original, dt).unwrap();
        let frozen_remainder = law.freeze_remainder(&remainder, dt).unwrap();
        frozen_remainder.validate_mixed().unwrap();
        let a = frozen_full
            .mixed_response(&f, 3500.0, VolumetricTreatment::ExactIncompressible)
            .unwrap();
        let b = frozen_remainder
            .mixed_response(&f, 3500.0, VolumetricTreatment::ExactIncompressible)
            .unwrap();
        assert_eq!(a.state.first_piola_pa, b.state.first_piola_pa);
        assert_eq!(
            a.state.potential_pa.to_bits(),
            b.state.potential_pa.to_bits()
        );
        assert_eq!(a.tangent_pa, b.tangent_pa);
        original = full.state;
        remainder = reduced.state;
    }
    assert!(remainder.dissipated > 0.0);
}

#[test]
fn g0_remainder_updates_are_bulk_independent_even_when_original_bulk_energy_overflows() {
    let reference = law(1e3);
    for kappa in [1e300, 1e308] {
        let extreme = law(kappa);
        let mut a = reference.initial_remainder_state();
        let mut b = extreme.initial_remainder_state();
        for (f, dt) in [
            (F, 0.0),
            (G, 0.08),
            ([3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0], 0.1),
        ] {
            let ra = reference.step_remainder_response(&a, &f, dt).unwrap();
            let rb = extreme.step_remainder_response(&b, &f, dt).unwrap();
            assert_eq!(ra, rb);
            for (x, y) in [
                ra.state.work,
                ra.state.stored,
                ra.state.dissipated,
                ra.incremental_potential_pa,
                ra.loading_work_pa,
                ra.dissipated_pa,
            ]
            .into_iter()
            .zip([
                rb.state.work,
                rb.state.stored,
                rb.state.dissipated,
                rb.incremental_potential_pa,
                rb.loading_work_pa,
                rb.dissipated_pa,
            ]) {
                assert_eq!(x.to_bits(), y.to_bits());
            }
            for (x, y) in ra
                .state
                .inverse_viscous_cauchy_green
                .iter()
                .flatten()
                .zip(rb.state.inverse_viscous_cauchy_green.iter().flatten())
            {
                assert_eq!(x.to_bits(), y.to_bits());
            }
            let frozen = extreme.freeze_remainder(&b, dt).unwrap();
            let response = frozen
                .mixed_response(&f, 1800.0, VolumetricTreatment::ExactIncompressible)
                .unwrap();
            assert!(response.state.potential_pa.is_finite());
            assert_eq!(response.state.bulk_stored_pa, 0.0);
            a = ra.state;
            b = rb.state;
        }
        if kappa == 1e308 {
            let expansion = [3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0];
            assert!(
                !extreme
                    .equilibrium
                    .volumetric_energy(det(&expansion))
                    .is_finite()
            );
            assert!(
                extreme
                    .step_response(&extreme.initial_state(), &expansion, 0.1)
                    .is_err()
            );
            // The finite remainder result above therefore cannot have depended
            // on evaluating and subtracting the original full bulk energy.
        }
    }
}

#[test]
fn g3_new_frozen_history_retains_incremental_work_derivatives_and_accepted_state() {
    let law = law(1e300);
    let accepted = law
        .step_remainder_response(&law.initial_remainder_state(), &F, 0.03)
        .unwrap()
        .state;
    let snapshot = accepted.clone();
    let frozen = law.freeze_remainder(&accepted, 0.08).unwrap();
    let direction = [0.3, -0.2, 0.1, 0.07, -0.13, 0.17, -0.08, 0.11, 0.09];
    let h = 2e-6;
    let plus = core::array::from_fn(|i| G[i] + h * direction[i]);
    let minus = core::array::from_fn(|i| G[i] - h * direction[i]);
    let response = frozen.remainder_response(&G).unwrap();
    let (pp, wp) = frozen.remainder_stress_and_potential(&plus).unwrap();
    let (pm, wm) = frozen.remainder_stress_and_potential(&minus).unwrap();
    let work: f64 = response
        .first_piola_pa
        .iter()
        .zip(direction)
        .map(|(p, d)| p * d)
        .sum();
    near((wp - wm) / (2.0 * h), work, 1e5, 2e-8);
    for row in 0..9 {
        let derivative: f64 = response.tangent_pa[row]
            .iter()
            .zip(direction)
            .map(|(k, d)| k * d)
            .sum();
        near((pp[row] - pm[row]) / (2.0 * h), derivative, 1e5, 2e-7);
    }
    let endpoint = law.step_remainder_response(&accepted, &G, 0.08).unwrap();
    near(
        response.potential_pa,
        endpoint.incremental_potential_pa,
        1e5,
        1e-14,
    );
    assert!(endpoint.incremental_potential_pa > endpoint.state.stored);
    assert_eq!(accepted, snapshot);
}

#[test]
fn g0_full_invariant_fibres_see_the_actual_volume_without_a_deformation_projection() {
    let (c10, k1, k2, mu) = (6000.0, 40_000.0, 1.2, 30_000.0);
    let law = ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10,
                k1,
                k2,
                dispersion: 0.0,
                fibres: vec![[1.0, 0.0, 0.0]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 1e300 },
            },
            4.0,
        )
        .unwrap(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: mu },
            tau_s: 0.2,
        }],
    )
    .unwrap();
    let f = [1.2, 0.0, 0.0, 0.0, 1.1, 0.0, 0.0, 0.0, 0.95];
    let result = law
        .step_remainder_response(&law.initial_remainder_state(), &f, 0.0)
        .unwrap();
    let i1 = f.iter().map(|x| x * x).sum::<f64>() * det(&f).powf(-2.0 / 3.0);
    let fibre_strain = f[0] * f[0] - 1.0;
    let expected = (c10 + 0.5 * mu) * (i1 - 3.0)
        + k1 / (2.0 * k2) * ((k2 * fibre_strain * fibre_strain).exp() - 1.0);
    near(result.state.stored, expected, 1e5, 2e-13);
    assert_eq!(result.state.deformation_gradient, f);
    assert!(det(&result.state.deformation_gradient) > 1.2);
    assert_eq!(result.dissipated_pa, 0.0);
    assert_eq!(
        result.state.inverse_viscous_cauchy_green[0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
}

#[test]
fn g0_invalid_remainder_ledgers_histories_durations_and_domains_refuse() {
    let law = law(1e300);
    let initial = law.initial_remainder_state();
    assert!(law.freeze_remainder(&initial, -0.1).is_err());
    let mut invalid = initial.clone();
    invalid.stored = f64::NAN;
    assert!(law.freeze_remainder(&invalid, 0.1).is_err());
    assert!(law.step_remainder_response(&invalid, &F, 0.1).is_err());
    invalid = initial.clone();
    invalid.dissipated = -1.0;
    assert!(law.freeze_remainder(&invalid, 0.1).is_err());
    invalid = initial.clone();
    invalid.inverse_viscous_cauchy_green[0][0] = -1.0;
    assert!(law.step_remainder_response(&invalid, &F, 0.1).is_err());
    assert!(
        law.step_remainder_response(
            &initial,
            &[4.0 * (1.0 + 1e-12), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            0.1
        )
        .is_err()
    );
    let mut coupled = law;
    coupled.equilibrium = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 6000.0,
            lambda: 24_000.0,
        },
        4.0,
    )
    .unwrap();
    assert!(coupled.freeze_remainder(&initial, 0.1).is_err());
    assert!(coupled.step_remainder_response(&initial, &F, 0.1).is_err());
}
