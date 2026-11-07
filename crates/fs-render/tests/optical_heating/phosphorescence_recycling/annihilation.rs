//! Nonlinear triplet encounters on the original optical/thermal history.
//! All encounter coefficients and active volumes below are manufactured inputs.

use super::*;
use fs_render::volumes::heating::luminescence::phosphorescence::PhosphorescenceAnnihilation;
use fs_thermal::luminescence::transfer::annihilation::TripletAnnihilation;
use fs_thermal::luminescence::transfer::recycling::annihilation::AnnihilationLimits;

const PAIR_RATE_M3_S: f64 = 1e-20;
const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;

fn controls() -> AnnihilationLimits {
    AnnihilationLimits {
        absolute_energy_tolerance_j: 1e-6,
        relative_tolerance: 1e-5,
        // Includes the graded endpoint segments as well as the base grid.
        maximum_substeps: 1024,
    }
}

fn annihilator(
    grid_cell: usize,
    active_volume_m3: f64,
    singlet_yield: f64,
) -> PhosphorescenceAnnihilation {
    PhosphorescenceAnnihilation {
        grid_cell,
        law: TripletAnnihilation::try_new(PAIR_RATE_M3_S, active_volume_m3, singlet_yield).unwrap(),
    }
}

fn within(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}, relative={relative:e}, absolute={absolute:e}"
    );
}

#[test]
fn g1_zero_encounter_rate_preserves_linear_results_and_no_return_paths_agree() {
    with_cx(|cx| {
        let (thermal, owner) = setup(cx, [500.0, 600.0], &[vec![[0.0; 2]]], &[[0.0; 2]], 3.0, 0.5);
        let zero = PhosphorescenceAnnihilation {
            grid_cell: 0,
            law: TripletAnnihilation::try_new(0.0, 1.0, 0.6).unwrap(),
        };
        let unchanged = owner
            .clone()
            .with_triplet_annihilation(&[zero], controls())
            .unwrap();
        let (_, linear, _) = thermal
            .advance_observed(&pulse(owner.clone(), 1.0), 3.0, config(1.0), || Ok(()))
            .unwrap();
        let (_, zero_rate, _) = thermal
            .advance_observed(&pulse(unchanged, 1.0), 3.0, config(1.0), || Ok(()))
            .unwrap();
        // A zero authored rate takes the original linear fast path exactly.
        assert_eq!(linear.energy(), zero_rate.energy());
        assert_eq!(linear.cell_emission_w(), zero_rate.cell_emission_w());
        assert_eq!(
            linear.state_stored_energy_j(0),
            zero_rate.state_stored_energy_j(0)
        );
        assert_eq!(zero_rate.energy().annihilation_heat_j, 0.0);

        let law = annihilator(0, 1.0, 0.6);
        let independent = owner
            .clone()
            .with_triplet_annihilation(&[law], controls())
            .unwrap();
        let network = owner
            .with_singlet_triplet_recycling(&[], &thermal, || Ok(()))
            .unwrap()
            .with_triplet_annihilation(&[law], controls())
            .unwrap();
        let (independent_thermal, independent, _) = thermal
            .advance_observed(&pulse(independent, 1.0), 3.0, config(1.0), || Ok(()))
            .unwrap();
        let (network_thermal, network, _) = thermal
            .advance_observed(&pulse(network, 1.0), 3.0, config(1.0), || Ok(()))
            .unwrap();
        // No-return transport must not change the independently composed law.
        // Compare at the explicitly authored nonlinear integration allowance.
        for state in 0..2 {
            within(
                independent.state_stored_energy_j(0).unwrap()[state],
                network.state_stored_energy_j(0).unwrap()[state],
                1e-5,
                3e-6,
            );
            within(
                independent.energy().emitted_j[state],
                network.energy().emitted_j[state],
                1e-5,
                3e-6,
            );
        }
        within(
            independent.energy().annihilation_heat_j,
            network.energy().annihilation_heat_j,
            1e-5,
            3e-6,
        );
        within(
            independent_thermal.temperature_k(CellId(1)).unwrap(),
            network_thermal.temperature_k(CellId(1)).unwrap(),
            0.0,
            1e-5,
        );
        assert!(independent.energy().annihilation_heat_j > 0.0);
        partition(&independent);
        partition(&network);
    });
}

#[test]
fn g1_dark_pairs_reach_the_independent_density_limited_steady_state_and_heat_the_cell() {
    with_cx(|cx| {
        let (thermal, owner) = setup(cx, [500.0, 600.0], &[vec![[0.0; 2]]], &[[0.0; 2]], 3.0, 0.5);
        let volume = 1.0;
        let owner = owner
            .with_triplet_annihilation(&[annihilator(0, volume, 0.0)], controls())
            .unwrap();
        let (heated, accepted, _) = thermal
            .advance_observed(&owner, 20.0, config(1.0), || Ok(()))
            .unwrap();

        // With no productive singlet branch the S equation is independent:
        // Us*=8 W/(1+3)/s=2 J, so the triplet production is G=6/Es per second.
        // The positive root of G - beta*Nt - 2*(k/V)*Nt^2 = 0 is evaluated
        // without subtracting two nearly equal roots.
        let es = HC_J_NM / 500.0;
        let et = HC_J_NM / 600.0;
        let production = 6.0 / es;
        let beta = 0.2;
        let coefficient = PAIR_RATE_M3_S / volume;
        let nt = 2.0 * production / (beta + (beta * beta + 8.0 * coefficient * production).sqrt());
        let stores = accepted.state_stored_energy_j(0).unwrap();
        // Twenty seconds suppress the physical transient below this explicit
        // 1e-4 comparison; the composition uses the tighter controls() above.
        within(stores[0], 2.0, 1e-4, 2e-5);
        within(stores[1], nt * et, 1e-4, 2e-5);
        assert!(stores[1] < 0.5 * (production / beta) * et);
        assert!(accepted.energy().annihilation_heat_j > 1.0);
        partition(&accepted);

        // The native sensible law has 4 J/K for this one-cell fixture. The
        // newly reported TTA heat must also enter the actual enthalpy timeline.
        close(
            heated.temperature_k(CellId(1)).unwrap(),
            300.0 + accepted.pump_energy().thermal_input_j / 4.0,
        );
        let e = accepted.energy();
        let heat_without_pairs =
            e.prompt_heat_j + e.nonradiative_heat_j + e.relaxation_heat_j + e.crossing_heat_j
                - e.reverse_crossing_heat_absorbed_j
                + e.secondary_absorption_heat_j;
        close(
            accepted.pump_energy().thermal_input_j - heat_without_pairs,
            e.annihilation_heat_j,
        );
    });
}

#[test]
fn g1_smaller_active_volume_increases_encounters_and_productive_dark_singlet_light() {
    with_cx(|cx| {
        let (thermal, owner) = setup(cx, [500.0, 600.0], &[vec![[0.0; 2]]], &[[0.0; 2]], 3.0, 0.5);
        let dilute = owner
            .clone()
            .with_triplet_annihilation(&[annihilator(0, 1.0, 0.8)], controls())
            .unwrap();
        let dense = owner
            .with_triplet_annihilation(&[annihilator(0, 0.125, 0.8)], controls())
            .unwrap();
        assert_ne!(
            dilute.emission_source_identity(),
            dense.emission_source_identity()
        );

        let mut afterglows = Vec::new();
        for source in [dilute, dense] {
            let (charged_thermal, charged, _) = thermal
                .advance_observed(&pulse(source, 1.0), 1.0, config(1.0), || Ok(()))
                .unwrap();
            let (_, dark, _) = charged_thermal
                .advance_observed(&charged, 5.0, config(1.0), || Ok(()))
                .unwrap();
            let dark_singlet_light = dark.energy().emitted_j[0] - charged.energy().emitted_j[0];
            assert!(dark_singlet_light > 0.0);
            partition(&charged);
            partition(&dark);
            afterglows.push((dark_singlet_light, dark));
        }
        // Independent manufactured reference, integrated with SciPy solve_ivp
        // (rtol=1e-11, atol=1e-13), not generated by the Rust composition:
        // Us' = P - 4 Us + phi*1.2*q*Ut^2,
        // Ut' = 2.5 Us - 0.2 Ut - 2*q*Ut^2,
        // singlet-light' = Us, TTA-heat' = (2 - 1.2*phi)*q*Ut^2,
        // q=k/(V*Et), phi=0.8, P=8 W on [0,1] and zero until t=5.
        // Entries: Us, Ut, total S light, TTA heat, pump-off S light [J].
        let reference = [
            [
                0.016_102_563_309_489_72,
                1.390_392_688_213_618_6,
                2.219_187_444_535_174,
                0.967_256_703_237_710_1,
                0.698_847_332_966_817_4,
            ],
            [
                0.013_395_616_786_840_786,
                0.436_610_103_681_364_16,
                2.570_325_699_730_817_7,
                2.485_923_283_685_961_5,
                0.992_042_041_289_961_5,
            ],
        ];
        for ((dark_light, afterglow), expected) in afterglows.iter().zip(reference) {
            let stores = afterglow.state_stored_energy_j(0).unwrap();
            for (actual, expected) in [
                stores[0],
                stores[1],
                afterglow.energy().emitted_j[0],
                afterglow.energy().annihilation_heat_j,
                *dark_light,
            ]
            .into_iter()
            .zip(expected)
            {
                // This includes the authored per-interval composition error
                // and its accumulation over the five physical intervals.
                within(actual, expected, 2e-4, 2e-4);
            }
        }
        let (dilute_light, dilute) = &afterglows[0];
        let (dense_light, dense) = &afterglows[1];
        assert!(dense_light > dilute_light);
        assert!(
            dense.state_stored_energy_j(0).unwrap()[1]
                < dilute.state_stored_energy_j(0).unwrap()[1]
        );
        assert!(dense.energy().annihilation_heat_j > dilute.energy().annihilation_heat_j);
        assert!(dense.energy().emitted_j[0] > dilute.energy().emitted_j[0]);
    });
}

#[test]
fn g1_productive_encounters_release_closed_triplet_return_and_heat_an_unpumped_receiver() {
    with_cx(|cx| {
        let rows = [vec![[0.0; 2], [0.0, 1.0]], vec![[0.0, 1.0], [0.0; 2]]];
        let (thermal, owner) = setup(cx, [500.0, 600.0], &rows, &[[0.0; 2]; 2], 3.0, 1.0);
        let owner = owner
            .with_singlet_triplet_recycling(
                &[receiver(0, 1, 1.0), receiver(1, 1, 1.0)],
                &thermal,
                || Ok(()),
            )
            .unwrap();
        let nonlinear = owner
            .clone()
            .with_triplet_annihilation(
                &[annihilator(0, 1.0, 0.8), annihilator(1, 1.0, 0.8)],
                controls(),
            )
            .unwrap();
        let (linear_thermal, linear, _) = thermal
            .advance_observed(&pulse(owner, 1.0), 20.0, config(1.0), || Ok(()))
            .unwrap();
        let (heated, accepted, _) = thermal
            .advance_observed(&pulse(nonlinear, 1.0), 20.0, config(1.0), || Ok(()))
            .unwrap();
        close(linear.energy().stored_excitation_j, 5.0);
        close(linear_thermal.temperature_k(CellId(2)).unwrap(), 300.0);
        assert!(accepted.energy().stored_excitation_j > 0.0);
        assert!(accepted.energy().stored_excitation_j < linear.energy().stored_excitation_j);
        assert!(accepted.energy().emitted_j[0] > linear.energy().emitted_j[0]);
        assert!(accepted.energy().annihilation_heat_j > 0.0);
        close(accepted.energy().escaped_j[1], 0.0);
        close(
            accepted.energy().recycled_excitation_j,
            accepted.energy().emitted_j[1],
        );
        // Cell 2 receives no beam and its resonant photons produce no passive
        // heat. Its new temperature rise requires the populated TTA/S pathways.
        assert!(heated.temperature_k(CellId(2)).unwrap() > 300.0);
        partition(&linear);
        partition(&accepted);
    });
}

#[test]
fn g4_admission_retry_cancellation_and_refinement_limits_preserve_accepted_history() {
    with_cx(|cx| {
        let (thermal, owner) = setup(cx, [500.0, 600.0], &[vec![[0.0; 2]]], &[[0.0; 2]], 3.0, 0.5);
        let law = annihilator(0, 1.0, 0.8);
        assert!(
            owner
                .clone()
                .with_triplet_annihilation(&[annihilator(9, 1.0, 0.8)], controls())
                .is_err()
        );
        assert!(
            owner
                .clone()
                .with_triplet_annihilation(&[law, law], controls())
                .is_err()
        );
        let invalid = AnnihilationLimits {
            absolute_energy_tolerance_j: 0.0,
            relative_tolerance: 0.0,
            maximum_substeps: 512,
        };
        assert!(
            owner
                .clone()
                .with_triplet_annihilation(&[law], invalid)
                .is_err()
        );
        let (_, already_advanced, _) = thermal
            .advance_observed(&owner, 1.0, config(1.0), || Ok(()))
            .unwrap();
        assert!(
            already_advanced
                .with_triplet_annihilation(&[law], controls())
                .is_err()
        );

        let limited = owner
            .clone()
            .with_triplet_annihilation(
                &[law],
                AnnihilationLimits {
                    maximum_substeps: 1,
                    ..controls()
                },
            )
            .unwrap();
        let mut refused = limited.clone();
        assert!(
            refused
                .step_sources(&thermal, 0.0, 1.0, 1.0, &mut || Ok(()))
                .is_err()
        );
        assert_eq!(
            refused.state_stored_energy_j(0),
            limited.state_stored_energy_j(0)
        );
        assert_eq!(refused.energy(), limited.energy());
        assert_eq!(refused.cell_emission_w(), limited.cell_emission_w());

        let bound = owner.with_triplet_annihilation(&[law], controls()).unwrap();
        let mut candidate = bound.clone();
        candidate
            .step_sources(&thermal, 0.0, 1.0, 1.0, &mut || Ok(()))
            .unwrap();
        let retried = candidate
            .step_sources(&thermal, 0.0, 0.5, 0.5, &mut || Ok(()))
            .unwrap();
        let fresh = bound
            .clone()
            .step_sources(&thermal, 0.0, 0.5, 0.5, &mut || Ok(()))
            .unwrap();
        assert_eq!(retried, fresh);
        assert_eq!(
            candidate.state_stored_energy_j(0),
            bound.state_stored_energy_j(0)
        );
        assert_eq!(candidate.energy(), bound.energy());

        let mut cancelled = bound.clone();
        let mut calls = 0;
        let result = cancelled.step_sources(&thermal, 0.0, 1.0, 1.0, &mut || {
            calls += 1;
            if calls >= 24 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(calls >= 24 && result.is_err());
        assert_eq!(
            cancelled.state_stored_energy_j(0),
            bound.state_stored_energy_j(0)
        );
        assert_eq!(cancelled.energy(), bound.energy());
        assert_eq!(cancelled.cell_emission_w(), bound.cell_emission_w());
    });
}
