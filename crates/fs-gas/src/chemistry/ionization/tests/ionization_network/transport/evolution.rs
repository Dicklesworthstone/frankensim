use super::*;
use crate::chemistry::collisional_radiative::CollisionalRadiativeThermal;
use crate::chemistry::radiation_transport::{
    RadiationTransportCell, RadiationTransportSolve, advance_radiation_transport,
};
use fs_qty::Power;

fn compare_states(actual: &FrozenIonizedGasStateV1, expected: &FrozenIonizedGasStateV1) {
    assert_eq!(actual.thermodynamics(), expected.thermodynamics());
    assert_eq!(actual.volume(), expected.volume());
    for (a, b) in actual.amounts().iter().zip(expected.amounts()) {
        assert!((a.value() - b.value()).abs() < 2e-9, "{a:?} versus {b:?}");
    }
    assert!(
        (actual.electron_temperature().value() - expected.electron_temperature().value()).abs()
            < 2e-3
    );
    assert_eq!(actual.heavy_temperature(), expected.heavy_temperature());
}

#[test]
fn sequential_helium_stages_share_finite_photons_and_preserve_charge_and_energy() {
    let (original, sources) = fixture(true, [0.0; 3]);
    // Seed only H ionization. Both higher helium stages start exactly absent;
    // the second helium edge must become active as the first creates He+.
    let state = FrozenIonizedGasStateV1::new(
        original.thermodynamics(),
        &[0.9, 0.1, 0.7, 0.0, 0.0, 0.1].map(Amount::new),
        original.volume(),
        original.heavy_temperature(),
        Temperature::new(6000.0),
        1e-12,
        || false,
    )
    .unwrap();
    let beams = [beam(60.0, 2e23, 0.65), beam(26.0, 5e22, -0.4)];
    let duration = 1e-6;
    let area = 100.0;
    let mut network = controls(1e-7);
    network.local.initial_step_s = duration / 4.0;
    network.local.maximum_step_s = duration / 4.0;
    let solve = NetworkRadiationSolve::new(network);
    let run = |channels: &[AtomicRadiativeRecombinationV1]| {
        advance_network_radiation_transport(
            &[NetworkRadiationCell::new(&state, channels)],
            Area::new(area),
            &beams,
            Time::new(duration),
            solve,
            || false,
        )
        .unwrap()
    };
    let result = run(&sources);
    let mut reversed = sources.clone();
    reversed.reverse();
    assert_eq!(result, run(&reversed));
    let cell = &result.cells()[0];
    let amounts = cell.state().amounts();
    assert_eq!(state.amounts()[3].value(), 0.0);
    assert_eq!(state.amounts()[4].value(), 0.0);
    assert!(amounts[3].value() > 0.0 && amounts[4].value() > 0.0);
    assert!(cell.channel_balances()[2].photoionization_extent_mol > 0.0);
    assert!(amounts.iter().all(|n| n.value() >= 0.0));
    relative(amounts[0].value() + amounts[1].value(), 1.0, 5e-13);
    relative(
        amounts[2].value() + amounts[3].value() + amounts[4].value(),
        0.7,
        5e-13,
    );
    relative(
        amounts[1].value() + amounts[3].value() + 2.0 * amounts[4].value(),
        amounts[5].value(),
        5e-13,
    );
    assert_eq!(cell.state().thermodynamics(), state.thermodynamics());
    assert_eq!(cell.state().volume(), state.volume());

    let mut absorbed_photons = 0.0;
    let mut absorbed_energy = 0.0;
    for (index, (beam, band)) in beams.iter().zip(result.beam_balances()).enumerate() {
        let incoming =
            beam.photons.number_flux_per_m2_s * area * beam.direction_cosine.abs() * duration;
        relative(band.incident_photons, incoming, 2e-12);
        relative(
            band.absorbed_photons + band.transmitted_photons,
            incoming,
            2e-12,
        );
        relative(
            band.absorbed_energy_j,
            band.absorbed_photons * beam.photons.photon_energy.value(),
            2e-12,
        );
        assert!(band.absorbed_photons <= incoming * (1.0 + 2e-15));
        assert!(band.photon_residual.abs() <= 2e-12 * incoming);
        assert_eq!(cell.beam_balances()[index], *band);
        absorbed_photons += band.absorbed_photons;
        absorbed_energy += band.absorbed_energy_j;
    }
    let balance = cell.balance();
    let channel_photons = cell
        .channel_balances()
        .iter()
        .map(|c| c.photoionization_extent_mol * NA)
        .sum::<f64>();
    relative(channel_photons, absorbed_photons, 2e-12);
    relative(balance.absorbed_photons, absorbed_photons, 2e-12);
    relative(balance.absorbed_radiation_energy_j, absorbed_energy, 2e-12);
    relative(
        balance.net_radiation_energy_j,
        absorbed_energy - balance.emitted_radiation_energy_j,
        2e-11,
    );
    let gas_change = balance.chemical_storage_change_j
        + balance.electron_sensible_energy_change_j
        + balance.heavy_sensible_energy_change_j;
    assert!(
        (gas_change - balance.net_radiation_energy_j).abs() <= network.local.energy_tolerance_j
    );
    assert!(balance.energy_residual_j.abs() <= network.local.energy_tolerance_j);
    assert!(balance.relative_charge_residual < 1e-12);
    assert!(balance.emitted_photons > 0.0);
}

#[test]
fn single_source_and_dark_limits_match_existing_chemistry_owners() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let beams = [beam(14.0, 2e24, 1.0)];
    let duration = Time::new(1e-8);
    let area = Area::new(100.0);
    let mut network = controls(1e-7);
    network.local.initial_step_s = 1e-9;
    network.local.maximum_step_s = 1e-9;
    let solve = NetworkRadiationSolve::new(network);
    let transported = advance_network_radiation_transport(
        &[NetworkRadiationCell::new(&state, &sources[..1])],
        area,
        &beams,
        duration,
        solve,
        || false,
    )
    .unwrap();
    let existing = advance_radiation_transport(
        &[RadiationTransportCell::new(&state, &sources[0])],
        &[beams[0].photons],
        area,
        duration,
        RadiationTransportSolve::new(network.local),
        || false,
    )
    .unwrap();
    compare_states(transported.cells()[0].state(), existing.cells()[0].state());
    relative(
        transported.beam_balances()[0].absorbed_photons,
        existing.group_balances()[0].absorbed_photons,
        2e-6,
    );
    relative(
        transported.beam_balances()[0].transmitted_photons,
        existing.group_balances()[0].transmitted_photons,
        2e-6,
    );

    let dark = advance_network_radiation_transport(
        &[NetworkRadiationCell::new(&state, &sources)],
        area,
        &[],
        duration,
        solve,
        || false,
    )
    .unwrap();
    let local =
        advance_ionization_network(&state, &sources, &[], duration, network, || false).unwrap();
    compare_states(dark.cells()[0].state(), local.state());
    assert!(dark.beam_balances().is_empty());
    let balance = dark.cells()[0].balance();
    assert_eq!(balance.absorbed_photons, 0.0);
    assert!(balance.emitted_photons > 0.0);
    relative(
        balance.net_radiation_energy_j,
        -balance.emitted_radiation_energy_j,
        2e-12,
    );
    relative(
        balance.emitted_photons,
        local.balance().emitted_photons,
        2e-6,
    );
}

#[test]
fn electron_birth_uses_the_combined_photoelectron_energy_without_a_seed() {
    let (state, sources) = fixture(true, [0.0; 3]);
    // Low-excess H photons dominate the birth pool; the separate He-ionizing
    // band also heats H. The combined birth temperature lies inside support,
    // even though H absorbing only 26 eV photons would exceed that support.
    let beams = [beam(14.0, 2e24, 1.0), beam(26.0, 2e23, 1.0)];
    let area = 100.0;
    let mut birth_mol_per_s = 0.0;
    let mut heat_w = 0.0;
    for beam in beams {
        let energy_ev = beam.photons.photon_energy.value() / EV;
        let depths = channel_depths(&state, area, energy_ev);
        let depth = depths.iter().sum::<f64>();
        let photons = area * beam.photons.number_flux_per_m2_s * -(-depth).exp_m1();
        for edge in 0..3 {
            let absorbed = photons * depths[edge] / depth;
            birth_mol_per_s += absorbed / NA;
            heat_w += absorbed * (energy_ev - GAPS_EV[edge]).max(0.0) * EV;
        }
    }
    let expected_temperature = heat_w / (CV * birth_mol_per_s);
    assert!((1000.0..50000.0).contains(&expected_temperature));
    let solve = NetworkRadiationSolve::new(controls(1e-7));
    let run = |channels: &[AtomicRadiativeRecombinationV1]| {
        advance_network_radiation_transport(
            &[NetworkRadiationCell::new(&state, channels)],
            Area::new(area),
            &beams,
            Time::new(1e-10),
            solve,
            || false,
        )
        .unwrap()
    };
    let result = run(&sources);
    let mut reordered = sources.clone();
    reordered.rotate_left(1);
    assert_eq!(result, run(&reordered));
    assert_eq!(state.amounts()[5].value(), 0.0);
    assert!(result.cells()[0].state().amounts()[5].value() > 0.0);
    assert!(result.cells()[0].state().amounts()[3].value() > 0.0);
    assert!(
        (result.cells()[0].state().electron_temperature().value() - expected_temperature).abs()
            < 1e-3
    );
    let balance = result.cells()[0].balance();
    assert!(balance.energy_residual_j.abs() <= solve.network.local.energy_tolerance_j);
    assert!(balance.relative_charge_residual < 1e-12);
}

#[test]
fn cancellation_work_and_source_refusals_leave_inputs_replayable() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let saved_state = state.clone();
    let saved_sources = sources.clone();
    let beams = [beam(60.0, 2e23, 1.0)];
    let cells = [NetworkRadiationCell::new(&state, &sources)];
    let solve = NetworkRadiationSolve::new(controls(1e-6));
    let duration = Time::new(1e-8);
    assert!(matches!(
        advance_network_radiation_transport(
            &cells,
            Area::new(100.0),
            &beams,
            duration,
            solve,
            || true
        ),
        Err(NetworkRadiationError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        advance_network_radiation_transport(
            &cells,
            Area::new(100.0),
            &beams,
            duration,
            solve,
            || {
                polls += 1;
                polls >= 8
            }
        ),
        Err(NetworkRadiationError::Cancelled)
    ));
    assert!(polls >= 8);

    let mut bounded = solve;
    bounded.maximum_radiation_work = 1;
    assert!(matches!(
        advance_network_radiation_transport(
            &cells,
            Area::new(100.0),
            &beams,
            duration,
            bounded,
            || false
        ),
        Err(NetworkRadiationError::WorkBudget)
    ));
    bounded = solve;
    bounded.maximum_radiation_evaluations = 1;
    assert!(matches!(
        advance_network_radiation_transport(
            &cells,
            Area::new(100.0),
            &beams,
            duration,
            bounded,
            || false
        ),
        Err(NetworkRadiationError::WorkBudget)
    ));
    assert!(
        advance_network_radiation_transport(
            &cells,
            Area::new(0.0),
            &beams,
            duration,
            solve,
            || false
        )
        .is_err()
    );
    assert!(
        advance_network_radiation_transport(
            &cells,
            Area::new(100.0),
            &[beam(60.0, 2e23, 0.0)],
            duration,
            solve,
            || false
        )
        .is_err()
    );
    // Unsupported source energies remain errors even for a dark beam: absent
    // or unilluminated channels must not silently acquire fabricated support.
    assert!(matches!(
        advance_network_radiation_transport(
            &cells,
            Area::new(100.0),
            &[beam(2000.0, 0.0, 1.0)],
            duration,
            solve,
            || false
        ),
        Err(NetworkRadiationError::Source { .. })
    ));
    assert_eq!(state, saved_state);
    assert_eq!(sources, saved_sources);
    let run = || {
        advance_network_radiation_transport(
            &cells,
            Area::new(100.0),
            &beams,
            duration,
            solve,
            || false,
        )
        .unwrap()
    };
    assert_eq!(run(), run());
}

#[test]
fn empty_columns_transmit_beams_and_zero_duration_preserves_the_original_state() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let beams = [beam(60.0, 2e23, 0.65), beam(26.0, 5e22, -0.4)];
    let area = Area::new(100.0);
    let duration = Time::new(2e-8);
    let solve = NetworkRadiationSolve::new(controls(1e-6));
    let empty =
        advance_network_radiation_transport(&[], area, &beams, duration, solve, || false).unwrap();
    assert!(empty.cells().is_empty());
    assert_eq!(empty.elapsed_time(), duration);
    for (beam, balance) in beams.iter().zip(empty.beam_balances()) {
        let expected = area.value()
            * beam.direction_cosine.abs()
            * beam.photons.number_flux_per_m2_s
            * duration.value();
        relative(balance.incident_photons, expected, 2e-12);
        relative(balance.transmitted_photons, expected, 2e-12);
        assert_eq!(balance.absorbed_photons, 0.0);
        assert_eq!(balance.absorbed_energy_j, 0.0);
        assert_eq!(balance.photon_residual, 0.0);
    }
    assert_eq!(empty.balance().gas_energy_change_j, 0.0);
    assert_eq!(empty.balance().emitted_photons, 0.0);
    assert_eq!(empty.balance().energy_residual_j, 0.0);

    let zero = advance_network_radiation_transport(
        &[NetworkRadiationCell::new(&state, &sources)],
        area,
        &beams,
        Time::new(0.0),
        solve,
        || false,
    )
    .unwrap();
    let cell = &zero.cells()[0];
    assert_eq!(cell.initial_state(), &state);
    assert_eq!(cell.state(), &state);
    assert_eq!(cell.sources(), sources);
    assert!(
        cell.channel_balances()
            .iter()
            .all(|b| *b == IonizationChannelBalance::default())
    );
    for balance in zero.beam_balances() {
        assert_eq!(balance.incident_photons, 0.0);
        assert_eq!(balance.transmitted_photons, 0.0);
        assert_eq!(balance.absorbed_photons, 0.0);
        assert_eq!(balance.absorbed_energy_j, 0.0);
    }
    assert_eq!(cell.beam_balances(), zero.beam_balances());
    assert_eq!(zero.balance().gas_energy_change_j, 0.0);
    assert_eq!(zero.balance().energy_residual_j, 0.0);
    assert_eq!(zero.balance().accepted_steps, 0);
    assert_eq!(zero.balance().rejected_steps, 0);
}

#[test]
fn held_pool_powers_act_once_during_mixed_species_beam_evolution() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let mut thermal =
        CollisionalRadiativeThermal::new(Temperature::new(1000.0), Temperature::new(50000.0));
    thermal.collisions = None;
    thermal.electron_power = Power::new(2e8);
    thermal.heavy_power = Power::new(3e7);
    let duration = 1e-7;
    let beams = [beam(14.0, 2e22, 1.0), beam(26.0, 5e21, -0.75)];
    let mut network = controls(1e-7);
    network.local.initial_step_s = 1e-8;
    network.local.maximum_step_s = 1e-8;
    let result = advance_network_radiation_transport(
        &[NetworkRadiationCell::new(&state, &sources).with_thermal(&thermal)],
        Area::new(100.0),
        &beams,
        Time::new(duration),
        NetworkRadiationSolve::new(network),
        || false,
    )
    .unwrap();
    let cell = &result.cells()[0];
    let balance = cell.balance();
    let whole = result.balance();
    let tolerance = network.local.energy_tolerance_j;
    let electron_work = thermal.electron_power.value() * duration;
    let heavy_work = thermal.heavy_power.value() * duration;
    assert_eq!(cell.sources(), sources);
    assert_eq!(cell.thermal_inputs(), Some(&thermal));
    assert_eq!(cell.state().thermodynamics(), state.thermodynamics());
    assert_eq!(cell.state().volume(), state.volume());
    assert_eq!(balance.electron_deposited_energy_j, electron_work);
    assert_eq!(balance.heavy_deposited_energy_j, heavy_work);
    assert_eq!(whole.external_work_j, electron_work + heavy_work);
    assert_eq!(balance.electron_to_heavy_energy_j, 0.0);
    assert!(balance.absorbed_photons > 0.0);
    assert!(balance.emitted_photons > 0.0);

    // Every atomic stage has the same monatomic Cv and ionization preserves
    // the total heavy-particle inventory. With no elastic exchange, this is
    // an independent exact heavy-temperature prediction, irrespective of the
    // three reaction edges, changing electron count, and photon absorption.
    let heavy_mol = state.amounts()[..5].iter().map(|n| n.value()).sum::<f64>();
    let capacity = CV * heavy_mol;
    let expected_heavy_temperature = state.heavy_temperature().value() + heavy_work / capacity;
    let actual_heavy_temperature = cell.state().heavy_temperature().value();
    assert!((actual_heavy_temperature - expected_heavy_temperature).abs() <= tolerance / capacity);
    assert!(
        (capacity * (actual_heavy_temperature - state.heavy_temperature().value()) - heavy_work)
            .abs()
            <= tolerance
    );
    assert!((balance.heavy_sensible_energy_change_j - heavy_work).abs() <= tolerance);

    let gas_change = balance.chemical_storage_change_j
        + balance.electron_sensible_energy_change_j
        + balance.heavy_sensible_energy_change_j;
    let cell_residual = gas_change + balance.emitted_radiation_energy_j
        - balance.absorbed_radiation_energy_j
        - electron_work
        - heavy_work;
    let column_residual = whole.gas_energy_change_j + whole.emitted_radiation_energy_j
        - whole.absorbed_radiation_energy_j
        - electron_work
        - heavy_work;
    for residual in [
        cell_residual,
        column_residual,
        balance.energy_residual_j,
        balance.electron_energy_residual_j,
        balance.heavy_energy_residual_j,
        whole.energy_residual_j,
    ] {
        assert!(
            residual.abs() <= tolerance,
            "energy residual {residual:e} J"
        );
    }
    let amounts = cell.state().amounts();
    assert!(amounts.iter().all(|n| n.value() >= 0.0));
    relative(amounts[0].value() + amounts[1].value(), 1.0, 5e-13);
    relative(
        amounts[2].value() + amounts[3].value() + amounts[4].value(),
        0.7,
        5e-13,
    );
    relative(
        amounts[1].value() + amounts[3].value() + 2.0 * amounts[4].value(),
        amounts[5].value(),
        5e-13,
    );
    assert!(balance.relative_charge_residual < 1e-12);
}
