//! The shared temperature calibration drives native finite-strain updates.
//! All energies, viscous metrics and tangents come from the original owner.

use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm};
use fs_material::visco::temperature::{
    HyperelasticTemperatureError, TemperatureShift, TemperatureShiftModel,
};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use std::error::Error;

const F: [f64; 9] = [1.2, 0.1, 0.0, 0.02, 0.94, 0.0, 0.0, 0.0, 0.91];
const G: [f64; 9] = [1.16, 0.07, 0.02, -0.01, 0.97, 0.03, 0.01, 0.0, 0.93];

fn reference() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6000.0,
                c01: 4000.0,
                kappa: 8e5,
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
    .with_local_newton(1e-13, 96)
    .unwrap()
}

fn shifts() -> [TemperatureShift; 2] {
    [
        TemperatureShift::new(
            TemperatureShiftModel::WilliamsLandelFerry {
                c1: 2.0,
                c2_k: 70.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 35_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
    ]
}

fn near(a: f64, b: f64, floor: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= 2e-10 * a.abs().max(b.abs()).max(floor),
        "{a:e} versus {b:e}"
    );
}

#[test]
fn calibrated_kinetics_reach_native_full_and_mixed_tangents() {
    let original = reference();
    let untouched = original.clone();
    let accepted = original
        .step_response(&original.initial_state(), &F, 0.01)
        .unwrap()
        .state;
    let accepted_remainder = original
        .step_remainder_response(&original.initial_remainder_state(), &F, 0.01)
        .unwrap()
        .state;
    for shift in shifts() {
        assert_eq!(
            shift.shift_hyperelastic(&original, 300.0).unwrap().law(),
            &original
        );
        for temperature in [280.0, 325.0, 360.0] {
            let shifted = shift.shift_hyperelastic(&original, temperature).unwrap();
            let mut explicit = original.clone();
            for branch in &mut explicit.branches {
                branch.tau_s *= shifted.state().factor();
            }
            assert_eq!(shifted.law(), &explicit);
            assert_eq!(shifted.state().temperature_k(), temperature);
            assert_eq!(shifted.law().local_newton_max_iterations(), 96);
            assert_eq!(
                shifted.law().step_response(&accepted, &G, 0.03).unwrap(),
                explicit.step_response(&accepted, &G, 0.03).unwrap()
            );
            let actual = shifted.law().freeze(&accepted, 0.03).unwrap();
            let expected = explicit.freeze(&accepted, 0.03).unwrap();
            assert_eq!(
                actual.consistent_tangent(&G).unwrap(),
                expected.consistent_tangent(&G).unwrap()
            );
            let mixed = shifted
                .law()
                .freeze_remainder(&accepted_remainder, 0.03)
                .unwrap();
            let mixed_reference = explicit
                .freeze_remainder(&accepted_remainder, 0.03)
                .unwrap();
            let a = mixed
                .mixed_response(&G, 3500.0, VolumetricTreatment::ExactIncompressible)
                .unwrap();
            let b = mixed_reference
                .mixed_response(&G, 3500.0, VolumetricTreatment::ExactIncompressible)
                .unwrap();
            assert_eq!(a.state.first_piola_pa, b.state.first_piola_pa);
            assert_eq!(a.tangent_pa, b.tangent_pa);
            assert_eq!(
                shifted
                    .law()
                    .step_remainder_response(&accepted_remainder, &G, 0.03)
                    .unwrap(),
                explicit
                    .step_remainder_response(&accepted_remainder, &G, 0.03)
                    .unwrap()
            );
        }
    }
    assert_eq!(original, untouched);
}

#[test]
fn temperature_history_carries_metrics_and_obeys_reduced_time() {
    let original = reference();
    for shift in shifts() {
        let loaded = original
            .step_response(&original.initial_state(), &F, 0.015)
            .unwrap()
            .state;
        let mut physical = loaded.clone();
        let mut reduced = loaded.clone();
        // Each interval freezes its explicitly supplied temperature. Both
        // trajectories retain all accepted metrics and accumulated ledgers.
        for (temperature, deformation, dt) in [
            (340.0, G, 0.007),
            (285.0, F, 0.009),
            (320.0, G, 0.011),
            (300.0, F, 0.006),
        ] {
            let shifted = shift.shift_hyperelastic(&original, temperature).unwrap();
            let old = physical.clone();
            // A temperature jump by itself neither resets memory nor adds heat.
            let jump = shifted
                .law()
                .step_response(&physical, &physical.deformation_gradient, 0.0)
                .unwrap();
            let same = original
                .step_response(&physical, &physical.deformation_gradient, 0.0)
                .unwrap();
            assert_eq!(jump, same);
            assert_eq!(jump.dissipated_pa, 0.0);

            let a = shifted
                .law()
                .step_response(&physical, &deformation, dt)
                .unwrap();
            let b = original
                .step_response(&reduced, &deformation, dt / shifted.state().factor())
                .unwrap();
            for (x, y) in a.first_piola_pa.iter().zip(b.first_piola_pa) {
                near(*x, y, 1e5);
            }
            for (x, y) in a
                .state
                .inverse_viscous_cauchy_green
                .iter()
                .flatten()
                .zip(b.state.inverse_viscous_cauchy_green.iter().flatten())
            {
                near(*x, *y, 1.0);
            }
            for (x, y) in [
                (a.state.work, b.state.work),
                (a.state.stored, b.state.stored),
                (a.state.dissipated, b.state.dissipated),
                (a.dissipated_pa, b.dissipated_pa),
            ] {
                near(x, y, 1e5);
            }
            assert!(a.dissipated_pa >= 0.0);
            assert!(a.state.dissipated >= old.dissipated);
            near(
                a.loading_work_pa - (a.state.stored - old.stored) - a.dissipated_pa,
                0.0,
                1e5,
            );
            assert_eq!(physical, old);
            physical = a.state;
            reduced = b.state;
        }
        assert!(physical.dissipated > loaded.dissipated);
    }
}

#[test]
fn heating_accelerates_relaxation_without_changing_elastic_storage() {
    let original = reference();
    let loaded = original
        .step_response(&original.initial_state(), &F, 0.0)
        .unwrap()
        .state;
    for shift in shifts() {
        let cold = shift.shift_hyperelastic(&original, 280.0).unwrap();
        let hot = shift.shift_hyperelastic(&original, 360.0).unwrap();
        assert!(hot.state().factor() < cold.state().factor());
        assert_eq!(
            hot.law().free_energy(&loaded, &F).unwrap(),
            cold.law().free_energy(&loaded, &F).unwrap()
        );
        let slow = cold.law().step_response(&loaded, &F, 0.02).unwrap();
        let fast = hot.law().step_response(&loaded, &F, 0.02).unwrap();
        assert!(fast.state.stored < slow.state.stored);
        assert!(fast.dissipated_pa > slow.dissipated_pa);
        near(fast.loading_work_pa, 0.0, 1e5);
        near(fast.state.stored + fast.dissipated_pa, loaded.stored, 1e5);
        assert_eq!(loaded.dissipated, 0.0);
    }
}

#[test]
fn invalid_reference_or_shift_refuses_before_any_history_change() {
    let shift = shifts()[0];
    let original = reference();
    assert!(matches!(
        shift.shift_hyperelastic(&original, 279.0),
        Err(HyperelasticTemperatureError::Shift(_))
    ));
    for tau in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = original.clone();
        invalid.branches[0].tau_s = tau;
        let error = shift.shift_hyperelastic(&invalid, 300.0).unwrap_err();
        assert!(
            error
                .source()
                .unwrap()
                .downcast_ref::<MaterialError>()
                .is_some()
        );
        assert!(matches!(error, HyperelasticTemperatureError::Material(_)));
    }
    let mut invalid_equilibrium = original.clone();
    invalid_equilibrium.equilibrium.stretch_limit = -1.0;
    assert!(matches!(
        shift.shift_hyperelastic(&invalid_equilibrium, 300.0),
        Err(HyperelasticTemperatureError::Material(_))
    ));
    for (tau, temperature) in [(f64::MAX, 280.0), (f64::from_bits(1), 360.0)] {
        let mut unrepresentable = original.clone();
        unrepresentable.branches[0].tau_s = tau;
        let before = unrepresentable.clone();
        assert!(matches!(
            shift.shift_hyperelastic(&unrepresentable, temperature),
            Err(HyperelasticTemperatureError::Shift(_))
        ));
        assert_eq!(unrepresentable, before);
    }
}
