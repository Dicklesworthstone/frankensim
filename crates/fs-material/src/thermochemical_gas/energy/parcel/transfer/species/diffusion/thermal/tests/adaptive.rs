//! G1/G3/G4 actual finite-gas trajectories through the existing thermal sources.
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::adaptive::{
    DiffusionNetworkAdvanceControls, DiffusionNetworkAdvanceError,
};
use network::{
    GasThermodiffusionLink, GasThermodiffusionNetworkAdvanceControls,
    GasThermodiffusionNetworkError,
};

fn accuracy(temperature: f64, heat: f64) -> GasThermodiffusionNetworkAdvanceControls {
    GasThermodiffusionNetworkAdvanceControls {
        temporal: DiffusionNetworkAdvanceControls {
            maximum_step_s: 1.0,
            minimum_step_s: 1e-9,
            maximum_evaluations: 100_000,
            absolute_species_tolerance_kg: 1e-10,
            relative_species_tolerance: 1e-5,
            temperature_tolerance_k: temperature,
        },
        absolute_species_transfer_tolerance_kg: 1e-10,
        absolute_additional_heat_tolerance_j: heat,
    }
}

fn pair_endpoints() -> BTreeMap<u64, GasEnergySolve> {
    BTreeMap::from([(1, energy_solve()), (2, energy_solve())])
}

#[test]
fn g1_pure_species_heat_relaxation_refines_against_independent_finite_capacity_solution() {
    let bank = model(vec![0.5, 0.5], 0.001);
    let data = cards([0.0; 2], 0.02);
    let context = point(400.0, 100000.0);
    let law = prepared(&bank, &data, &context);
    let original = BTreeMap::from([
        (1, gas(&bank, 0.0, 300.0, 100000.0)),
        (2, gas(&bank, 0.0, 500.0, 100000.0)),
    ]);
    let saved = original.clone();
    let links = [GasThermodiffusionLink {
        id: 7,
        parcels: [1, 2],
        geometry: geometry(),
        model: &law,
        solve: solve(),
    }];
    // Pure B has M=0.004 kg/mol and constant Cv=1.5R; no species can move
    // in the mass-average frame. Both finite caloric capacities remain fixed.
    let capacities = [1, 2].map(|id| original[&id].mass_kg() / 0.004 * 1.5 * R);
    let equilibrium =
        (capacities[0] * 300.0 + capacities[1] * 500.0) / (capacities[0] + capacities[1]);
    let conductance = 0.02 * geometry().area_m2 / geometry().distance_m;
    let duration = 2.0;
    let decay = (-conductance * (1.0 / capacities[0] + 1.0 / capacities[1]) * duration).exp();
    let expected = [
        equilibrium + (300.0 - equilibrium) * decay,
        equilibrium + (500.0 - equilibrium) * decay,
    ];
    let mut previous = f64::INFINITY;
    let mut previous_work = 0;
    for tolerance in [1e-2, 2.5e-3, 6.25e-4] {
        let candidate = network::advance_candidate(
            &original,
            &links,
            &mut [law.workspace().unwrap()],
            GasTransportPaths::default(),
            duration,
            &pair_endpoints(),
            network_solve(),
            accuracy(tolerance, 1e6),
            || false,
        )
        .unwrap();
        let endpoint = [1, 2].map(|id| candidate.parcels()[&id].state().state().temperature);
        let error = (endpoint[0] - expected[0])
            .abs()
            .max((endpoint[1] - expected[1]).abs());
        assert!(error < previous);
        assert!(candidate.report().evaluations > previous_work);
        previous = error;
        previous_work = candidate.report().evaluations;
        assert!(candidate.report().transfers.is_empty());
        assert!(candidate.report().maximum_temperature_error_ratio <= 1.0);
        assert!(candidate.report().maximum_transfer_error_ratio.unwrap() <= 1.0);
        assert!(
            candidate.report().endpoint_energy_change_j.abs()
                <= network_solve().exchange.energy_tolerance_j
        );
        assert!(candidate.report().entropy_change_j_k > 0.0);
        let heat = candidate
            .report()
            .heat
            .iter()
            .find(|row| row.parcels == [1, 2])
            .unwrap();
        near(
            -heat.heat_left_to_right_j,
            capacities[0] * (endpoint[0] - 300.0),
            1e-8,
        );
        near(
            heat.heat_left_to_right_j,
            capacities[1] * (endpoint[1] - 500.0),
            1e-8,
        );
        // Source samples describe the actual START of the final accepted half
        // step, rather than pretending to be endpoint face evaluations.
        let last = &candidate.report().last_substep;
        assert_eq!(last.links[0].id, 7);
        assert!(last.links[0].flux.face_temperature_k != 400.0);
        assert!(
            (candidate.report().minimum_accepted_substep_s
                ..=candidate.report().maximum_accepted_substep_s)
                .contains(&last.transport.duration_s)
        );
    }
    assert!(previous < 0.1);
    assert_eq!(original, saved);
}

#[test]
fn g1_additional_heat_integrals_have_their_own_accuracy_even_with_loose_endpoint_control() {
    let bank = model(vec![0.5, 0.5], 0.001);
    let data = cards([0.0; 2], 0.02);
    let context = point(400.0, 100000.0);
    let law = prepared(&bank, &data, &context);
    let original = BTreeMap::from([
        (1, gas(&bank, 0.0, 300.0, 100000.0)),
        (2, gas(&bank, 0.0, 500.0, 100000.0)),
    ]);
    // Distinct parallel exposures have opposite AUTHORED orientation. Each
    // original link is compared independently before the native heat merge.
    let links = [
        GasThermodiffusionLink {
            id: 9,
            parcels: [1, 2],
            geometry: geometry(),
            model: &law,
            solve: solve(),
        },
        GasThermodiffusionLink {
            id: 2,
            parcels: [2, 1],
            geometry: geometry(),
            model: &law,
            solve: solve(),
        },
    ];
    let run = |tolerance| {
        network::advance_candidate(
            &original,
            &links,
            &mut [law.workspace().unwrap(), law.workspace().unwrap()],
            GasTransportPaths::default(),
            1.0,
            &pair_endpoints(),
            network_solve(),
            accuracy(1e6, tolerance),
            || false,
        )
        .unwrap()
    };
    let loose = run(1e6);
    let tight = run(1e-5);
    assert_eq!(loose.report().evaluations, 3);
    assert!(tight.report().evaluations > loose.report().evaluations);
    assert!(tight.report().maximum_transfer_error_ratio.unwrap() <= 1.0);
    assert_eq!(tight.report().last_substep.links[0].id, 2);
    assert_eq!(tight.report().last_substep.links[1].id, 9);
    assert!(
        tight.report().last_substep.links[0]
            .flux
            .additional_heat_rate_left_to_right_w
            > 0.0
    );
    assert!(
        tight.report().last_substep.links[1]
            .flux
            .additional_heat_rate_left_to_right_w
            < 0.0
    );
    assert_eq!(tight.report().heat.len(), 1);
    assert!(tight.report().heat[0].heat_left_to_right_j < 0.0);
}

#[test]
fn g3_thermal_trajectory_retains_soret_and_caloric_reference_invariance() {
    let data = cards([100000.0, -30000.0], 0.02);
    let context = point(400.0, 100000.0);
    let mut results = Vec::new();
    for bank in [shifted_bank([0.0, 0.0]), shifted_bank([1e24, -3e24])] {
        let law = prepared(&bank, &data, &context);
        let original = BTreeMap::from([
            (1, gas(&bank, 0.4, 300.0, 100000.0)),
            (2, gas(&bank, 0.4, 500.0, 100000.0)),
        ]);
        let links = [GasThermodiffusionLink {
            id: 1,
            parcels: [1, 2],
            geometry: geometry(),
            model: &law,
            solve: solve(),
        }];
        let result = network::advance_candidate(
            &original,
            &links,
            &mut [law.workspace().unwrap()],
            GasTransportPaths::default(),
            0.1,
            &pair_endpoints(),
            network_solve(),
            accuracy(1e-3, 1e-5),
            || false,
        )
        .unwrap();
        assert!(result.report().accepted_substeps > 2);
        assert!(result.report().entropy_change_j_k > 0.0);
        assert!(!result.report().transfers.is_empty());
        assert_ne!(
            result.parcels()[&1]
                .model()
                .thermodynamics()
                .mole_composition(),
            original[&1].model().thermodynamics().mole_composition(),
        );
        assert!(result.report().maximum_species_error_ratio <= 1.0);
        assert!(result.report().maximum_temperature_error_ratio <= 1.0);
        assert!(result.report().maximum_transfer_error_ratio.unwrap() <= 1.0);
        assert!(
            result.report().maximum_species_mass_residual_kg
                <= network_solve().exchange.mass_tolerance_kg
        );
        assert!(
            result.report().endpoint_energy_change_j.abs()
                <= network_solve().exchange.energy_tolerance_j
        );
        results.push(result);
    }
    for id in [1, 2] {
        assert_eq!(
            results[0].parcels()[&id].state().state().temperature,
            results[1].parcels()[&id].state().state().temperature
        );
        assert_eq!(
            results[0].parcels()[&id].mass_kg(),
            results[1].parcels()[&id].mass_kg()
        );
        assert_eq!(
            results[0].parcels()[&id]
                .model()
                .thermodynamics()
                .mole_composition(),
            results[1].parcels()[&id]
                .model()
                .thermodynamics()
                .mole_composition()
        );
    }
    assert_eq!(
        results[0].report().evaluations,
        results[1].report().evaluations
    );
    assert_eq!(results[0].report().heat, results[1].report().heat);
    assert_eq!(
        results[0].report().maximum_transfer_error_ratio,
        results[1].report().maximum_transfer_error_ratio
    );
    for (a, b) in results[0]
        .report()
        .transfers
        .iter()
        .zip(&results[1].report().transfers)
    {
        assert_eq!(
            (a.donor, a.receiver, &a.species, a.mass_kg),
            (b.donor, b.receiver, &b.species, b.mass_kg)
        );
        // The gross original partial enthalpy legitimately carries its datum.
        assert_ne!(a.carried_enthalpy_j, b.carried_enthalpy_j);
    }
}

#[test]
fn g0_zero_coupling_retains_the_original_equilibrium_and_empty_transport() {
    let bank = model(vec![0.4, 0.6], 0.001);
    let data = cards([0.0; 2], 0.0);
    let context = point(400.0, 100000.0);
    let law = prepared(&bank, &data, &context);
    let original = BTreeMap::from([
        (1, gas(&bank, 0.4, 400.0, 100000.0)),
        (2, gas(&bank, 0.4, 400.0, 100000.0)),
    ]);
    let links = [GasThermodiffusionLink {
        id: 1,
        parcels: [1, 2],
        geometry: geometry(),
        model: &law,
        solve: solve(),
    }];
    let candidate = network::advance_candidate(
        &original,
        &links,
        &mut [law.workspace().unwrap()],
        GasTransportPaths::default(),
        1.0,
        &pair_endpoints(),
        network_solve(),
        accuracy(1e-3, 1e-8),
        || false,
    )
    .unwrap();
    assert!(candidate.report().transfers.is_empty());
    assert!(candidate.report().heat.is_empty());
    assert_eq!(candidate.report().maximum_transfer_error_ratio, Some(0.0));
    assert_eq!(candidate.report().evaluations, 3);
    for id in [1, 2] {
        near(
            candidate.parcels()[&id].state().state().temperature,
            400.0,
            1e-10,
        );
        assert_eq!(candidate.parcels()[&id].mass_kg(), original[&id].mass_kg());
        assert_eq!(
            candidate.parcels()[&id]
                .model()
                .thermodynamics()
                .mole_composition(),
            original[&id].model().thermodynamics().mole_composition()
        );
    }
}

#[test]
fn g4_work_source_and_late_cancellation_refusals_publish_no_partial_thermal_network() {
    let bank = model(vec![0.5, 0.5], 0.001);
    let data = cards([0.0; 2], 0.02);
    let context = point(400.0, 100000.0);
    let law = prepared(&bank, &data, &context);
    let original = BTreeMap::from([
        (1, gas(&bank, 0.0, 300.0, 100000.0)),
        (2, gas(&bank, 0.0, 500.0, 100000.0)),
    ]);
    let saved = original.clone();
    let links = [GasThermodiffusionLink {
        id: 1,
        parcels: [1, 2],
        geometry: geometry(),
        model: &law,
        solve: solve(),
    }];
    let mut controls = accuracy(1e6, 1e6);
    controls.temporal.maximum_evaluations = 3;
    assert!(matches!(
        network::advance_candidate(
            &original,
            &links,
            &mut [law.workspace().unwrap()],
            GasTransportPaths::default(),
            2.0,
            &pair_endpoints(),
            network_solve(),
            controls,
            || false,
        ),
        Err(DiffusionNetworkAdvanceError::WorkLimit {
            elapsed_s: 1.0,
            accepted_intervals: 1,
            evaluations: 3
        })
    ));
    let mut calls = 0;
    assert_eq!(
        network::advance_candidate(
            &original,
            &links,
            &mut [law.workspace().unwrap()],
            GasTransportPaths::default(),
            2.0,
            &pair_endpoints(),
            network_solve(),
            accuracy(1e-3, 1e-6),
            || {
                calls += 1;
                calls == 70
            },
        ),
        Err(DiffusionNetworkAdvanceError::Cancelled)
    );
    assert_eq!(calls, 70);
    let mut unsupported = cards([0.0; 2], 0.02);
    unsupported[1] = scalar(
        THERMAL_CONDUCTIVITY_PROPERTY,
        THERMAL_CONDUCTIVITY_DIMS,
        0.02,
        450.0,
    );
    let restricted = prepared(&bank, &unsupported, &context);
    let invalid_links = [GasThermodiffusionLink {
        id: 11,
        parcels: [1, 2],
        geometry: geometry(),
        model: &restricted,
        solve: solve(),
    }];
    assert!(matches!(
        network::advance_candidate(
            &original,
            &invalid_links,
            &mut [restricted.workspace().unwrap()],
            GasTransportPaths::default(),
            2.0,
            &pair_endpoints(),
            network_solve(),
            accuracy(1e-3, 1e-6),
            || false,
        ),
        Err(DiffusionNetworkAdvanceError::ThermalStep {
            elapsed_s: 0.0,
            step_s: 1.0,
            error: GasThermodiffusionNetworkError::Face { id: 11, .. },
        })
    ));
    assert_eq!(original, saved);
}
