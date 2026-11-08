use super::*;
use crate::chemistry::diffuse_recombination::DiffuseRecombinationError;
use crate::chemistry::diffuse_recombination::evolution::{
    DiffuseRecombinationEvolutionCell, DiffuseRecombinationEvolutionSolve,
    advance_diffuse_recombination,
};

fn live_controls(duration: f64) -> NetworkDiffuseSolve {
    let mut solve = diffuse_controls(1e-8);
    solve.network.network.local.initial_step_s = duration / 4.0;
    solve.network.network.local.maximum_step_s = duration / 4.0;
    solve
}

fn audit(result: &NetworkRadiationResult, solve: NetworkDiffuseSolve) {
    let tolerance = solve.network.network.local.energy_tolerance_j;
    let summary = result.diffuse_radiation().unwrap();
    let mut gas = 0.0;
    for cell in result.cells() {
        let initial = cell.initial_state();
        let state = cell.state();
        let before = initial.amounts();
        let after = state.amounts();
        assert!(after.iter().all(|n| n.value() >= 0.0));
        relative(after[0].value() + after[1].value(), 1.0, 5e-13);
        relative(
            after[2].value() + after[3].value() + after[4].value(),
            0.7,
            5e-13,
        );
        relative(
            after[1].value() + after[3].value() + 2.0 * after[4].value(),
            after[5].value(),
            5e-13,
        );
        assert_eq!(state.thermodynamics(), initial.thermodynamics());
        assert_eq!(state.volume(), initial.volume());
        relative(state.mass().value(), initial.mass().value(), 5e-13);
        let chemical = NA
            * EV
            * (GAPS_EV[0] * (after[1].value() - before[1].value())
                + GAPS_EV[1]
                    * (after[3].value() - before[3].value() + after[4].value()
                        - before[4].value())
                + GAPS_EV[2] * (after[4].value() - before[4].value()));
        let electron = CV
            * (before[5].value()
                * (state.electron_temperature().value() - initial.electron_temperature().value())
                + (after[5].value() - before[5].value()) * state.electron_temperature().value());
        let heavy = CV
            * before[..5].iter().map(|n| n.value()).sum::<f64>()
            * (state.heavy_temperature().value() - initial.heavy_temperature().value());
        let actual = chemical + electron + heavy;
        let exposure = cell.diffuse_radiation().unwrap();
        let balance = cell.balance();
        assert!((actual - exposure.net_gas_energy_j).abs() <= tolerance);
        assert!((balance.net_radiation_energy_j - exposure.net_gas_energy_j).abs() <= tolerance);
        assert!(
            (exposure.diffuse_absorbed_energy_j - balance.absorbed_radiation_energy_j).abs()
                <= tolerance
        );
        assert!(
            (exposure.emitted_energy_j - balance.emitted_radiation_energy_j).abs() <= tolerance
        );
        relative(
            exposure.diffuse_absorbed_photons(),
            balance.absorbed_photons,
            3e-12,
        );
        relative(exposure.emitted_photons, balance.emitted_photons, 3e-12);
        assert!(balance.energy_residual_j.abs() <= tolerance);
        assert!(balance.relative_charge_residual < 1e-12);
        gas += actual;
    }
    let balance = summary.balance();
    let boundary_tolerance = tolerance + balance.boundary_energy_roundoff_j;
    assert_eq!(balance.incident_photons, 0.0);
    assert_eq!(balance.incident_energy_j, 0.0);
    assert!((gas + balance.escaped_energy_j).abs() <= boundary_tolerance);
    assert!(
        (balance.net_radiation_energy_j + balance.escaped_energy_j).abs() <= boundary_tolerance
    );
    assert!(balance.radiation_energy_residual_j.abs() <= boundary_tolerance);
    assert!(balance.energy_residual_j.abs() <= tolerance);
    assert!(balance.photon_residual.abs() <= 3e-12 * balance.emitted_photons);
    relative(
        -balance.net_ionization_extent_mol * NA,
        balance.escaped_photons,
        3e-11,
    );
    let interfaces = summary.interfaces();
    assert_eq!(interfaces[0].right_photons, 0.0);
    assert_eq!(interfaces.last().unwrap().left_photons, 0.0);
    relative(
        interfaces[0].left_photons + interfaces.last().unwrap().right_photons,
        balance.escaped_photons,
        3e-12,
    );
    relative(
        interfaces[0].left_energy_j + interfaces.last().unwrap().right_energy_j,
        balance.escaped_energy_j,
        3e-12,
    );
}

#[test]
fn one_channel_live_recycling_matches_the_existing_diffuse_owner() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let duration = 1e-8;
    let area = Area::new(100.0);
    let solve = live_controls(duration);
    let result = advance_network_diffuse_recombination(
        &[NetworkRadiationCell::new(&state, &sources[..1])],
        area,
        &[],
        Time::new(duration),
        solve,
        || false,
    )
    .unwrap();
    let mut original_controls =
        DiffuseRecombinationEvolutionSolve::new(solve.network.network.local);
    original_controls.radiation = solve.radiation;
    let original = advance_diffuse_recombination(
        &[DiffuseRecombinationEvolutionCell::new(&state, &sources[0])],
        area,
        Time::new(duration),
        original_controls,
        || false,
    )
    .unwrap();
    for (actual, expected) in result.cells()[0]
        .state()
        .amounts()
        .iter()
        .zip(original.cells()[0].state().amounts())
    {
        assert!((actual.value() - expected.value()).abs() < 3e-9);
    }
    assert!(
        (result.cells()[0].state().electron_temperature().value()
            - original.cells()[0].state().electron_temperature().value())
        .abs()
            < 2e-3
    );
    relative(
        result
            .diffuse_radiation()
            .unwrap()
            .balance()
            .escaped_photons,
        original.balance().escaped_photons,
        3e-6,
    );
    relative(
        result
            .diffuse_radiation()
            .unwrap()
            .balance()
            .escaped_energy_j,
        original.balance().escaped_energy_j,
        3e-6,
    );
    relative(
        result.cells()[0].channel_balances()[0].radiative_recombination_extent_mol,
        original.cells()[0]
            .balance()
            .radiative_recombination_extent_mol,
        3e-6,
    );
    assert!(
        result.cells()[0]
            .diffuse_radiation()
            .unwrap()
            .self_absorbed_photons
            > 0.0
    );
    audit(&result, solve);
}

#[test]
fn nonlocal_recycling_creates_receiver_electrons_and_refreshes_multiple_charge_stages() {
    let (neutral, sources) = fixture(true, [0.0; 3]);
    // H dominates emission while both upper helium stages are present. This
    // keeps physical photoelectron birth inside the original 50k K support;
    // it does not introduce any seed electrons into the neutral receiver.
    let donor = mixture_state(
        &neutral,
        [0.8, 0.2, 0.699, 0.000999, 0.000001, 0.201001],
        20000.0,
    );
    let duration = 2e-8;
    let area = Area::new(100.0);
    let solve = live_controls(duration);
    let mut endpoints = Vec::new();
    for (stored_temperature, reverse) in [(1000.0, false), (45000.0, true)] {
        let receiver = mixture_state(&neutral, [1.0, 0.0, 0.7, 0.0, 0.0, 0.0], stored_temperature);
        let mut ordered = sources.clone();
        if reverse {
            ordered.reverse();
        }
        let cells = [
            NetworkRadiationCell::new(&donor, &ordered),
            NetworkRadiationCell::new(&receiver, &ordered),
        ];
        let initial_field =
            evaluate_network_diffuse_recombination(&cells, area, &[], solve, || false).unwrap();
        let result = advance_network_diffuse_recombination(
            &cells,
            area,
            &[],
            Time::new(duration),
            solve,
            || false,
        )
        .unwrap();
        let emitter = &result.cells()[0];
        let born = &result.cells()[1];
        assert_eq!(receiver.amounts()[5].value(), 0.0);
        assert!(born.state().amounts()[1].value() > 0.0);
        assert!(born.state().amounts()[3].value() > 0.0);
        assert!(born.state().amounts()[5].value() > 0.0);
        assert!(born.diffuse_radiation().unwrap().nonlocal_absorbed_photons > 0.0);
        assert!(born.diffuse_radiation().unwrap().net_electron_energy_j > 0.0);
        assert_ne!(emitter.state().amounts()[3], donor.amounts()[3]);
        assert_ne!(emitter.state().amounts()[4], donor.amounts()[4]);
        assert_ne!(
            emitter.state().electron_temperature(),
            donor.electron_temperature()
        );
        for cell in result.cells() {
            assert_eq!(cell.sources(), sources);
        }
        let end_cells = [
            NetworkRadiationCell::new(emitter.state(), &sources),
            NetworkRadiationCell::new(born.state(), &sources),
        ];
        let final_field =
            evaluate_network_diffuse_recombination(&end_cells, area, &[], solve, || false).unwrap();
        assert_ne!(
            initial_field.cells()[0].channel_rates(),
            final_field.cells()[0].channel_rates()
        );
        assert!(
            final_field.cells()[1]
                .channel_rates()
                .iter()
                .any(|r| r.emitted_photons_per_s > 0.0)
        );
        audit(&result, solve);
        endpoints.push([emitter.state().clone(), born.state().clone()]);
    }
    // Neither source permutation nor the meaningless temperature stored in an
    // empty electron pool may change the physical reacting endpoint.
    for (base, replay) in endpoints[0].iter().zip(&endpoints[1]) {
        for (actual, expected) in base.amounts().iter().zip(replay.amounts()) {
            assert!((actual.value() - expected.value()).abs() < 3e-10);
        }
        assert!(
            (base.electron_temperature().value() - replay.electron_temperature().value()).abs()
                < 2e-3
        );
    }
}

#[test]
fn diffuse_source_work_and_cancellation_refusals_do_not_publish_partial_states() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let saved_state = state.clone();
    let saved_sources = sources.clone();
    let duration = 1e-9;
    let area = Area::new(100.0);
    let solve = live_controls(duration);
    let cells = [NetworkRadiationCell::new(&state, &sources)];
    assert!(matches!(
        advance_network_diffuse_recombination(
            &cells,
            area,
            &[],
            Time::new(duration),
            solve,
            || true
        ),
        Err(NetworkRadiationError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        advance_network_diffuse_recombination(
            &cells,
            area,
            &[],
            Time::new(duration),
            solve,
            || {
                polls += 1;
                polls >= 12
            }
        ),
        Err(NetworkRadiationError::Cancelled)
    ));
    assert!(polls >= 12);
    let mut bounded = solve;
    bounded.radiation.maximum_work = 1;
    assert!(matches!(
        advance_network_diffuse_recombination(
            &cells,
            area,
            &[],
            Time::new(duration),
            bounded,
            || false
        ),
        Err(NetworkRadiationError::WorkBudget)
            | Err(NetworkRadiationError::Diffuse(
                DiffuseRecombinationError::WorkBudget
            ))
    ));
    bounded = solve;
    bounded.network.maximum_radiation_evaluations = 1;
    assert!(matches!(
        advance_network_diffuse_recombination(
            &cells,
            area,
            &[],
            Time::new(duration),
            bounded,
            || false
        ),
        Err(NetworkRadiationError::WorkBudget)
    ));

    let gap = sources[0].source().channel().ionization_energy_j_mol() / NA;
    let provenance = Provenance {
        source: "manufactured truncated bound-free support for domain refusal".into(),
        license: "analytic test fixture".into(),
        artifact: None,
    };
    let short = AtomicPhotoionizationV1::new(
        sources[0].source().channel().clone(),
        provenance.clone(),
        [0.0, 0.5]
            .map(|excess| PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap + excess * EV),
                cross_section: Area::new(SIGMA[0]),
            })
            .to_vec(),
    )
    .unwrap();
    let short =
        AtomicRadiativeRecombinationV1::new(short, Area::new(SIGMA[0]), provenance).unwrap();
    let mut unsupported = sources.clone();
    unsupported[0] = short;
    assert!(
        advance_network_diffuse_recombination(
            &[NetworkRadiationCell::new(&state, &unsupported)],
            area,
            &[],
            Time::new(duration),
            solve,
            || false
        )
        .is_err()
    );
    assert!(
        advance_network_diffuse_recombination(
            &cells,
            area,
            &[beam(2000.0, 0.0, 1.0)],
            Time::new(duration),
            solve,
            || false
        )
        .is_err()
    );
    assert_eq!(state, saved_state);
    assert_eq!(sources, saved_sources);
    let run = || {
        advance_network_diffuse_recombination(&cells, area, &[], Time::new(duration), solve, || {
            false
        })
        .unwrap()
    };
    assert_eq!(run(), run());
}
