use super::*;
use network::adaptive::{
    DiffusionNetworkAdvanceControls, DiffusionNetworkAdvanceError as AdvanceError,
};
use network::binary::{
    BinaryDiffusionNetworkAdvanceCandidate, BinaryDiffusionNetworkAdvanceControls,
    advance_adaptive_candidate,
};

fn accuracy(relative: f64, maximum_step_s: f64) -> BinaryDiffusionNetworkAdvanceControls {
    BinaryDiffusionNetworkAdvanceControls {
        temporal: DiffusionNetworkAdvanceControls {
            maximum_step_s,
            minimum_step_s: 1e-6,
            maximum_evaluations: 10_000,
            absolute_species_tolerance_kg: 1e-14,
            relative_species_tolerance: relative,
            temperature_tolerance_k: 1e-6,
        },
        maximum_pair_evaluations: 100_000,
    }
}

fn balances(
    before: &BTreeMap<u64, ClosedGasParcel>,
    after: &BinaryDiffusionNetworkAdvanceCandidate,
    physical: DiffusionNetworkSolve,
) {
    let report = after.report();
    assert!(report.total_mass_residual_kg.abs() <= physical.exchange.mass_tolerance_kg);
    assert!(report.maximum_species_mass_residual_kg <= physical.exchange.mass_tolerance_kg);
    assert!(report.total_energy_residual_j.abs() <= physical.exchange.energy_tolerance_j);
    assert!(report.endpoint_energy_change_j.abs() <= physical.exchange.energy_tolerance_j);
    assert!(report.entropy_change_j_k >= -physical.entropy_tolerance_j_k);
    near(
        after.parcels().values().map(internal_energy).sum(),
        before.values().map(internal_energy).sum(),
        physical.exchange.energy_tolerance_j,
    );
    for name in ["A", "B"] {
        let initial: f64 = before.values().map(|p| species_mass(p, name)).sum();
        let final_mass: f64 = after
            .parcels()
            .values()
            .map(|p| species_mass(p, name))
            .sum();
        near(initial, final_mass, physical.exchange.mass_tolerance_kg);
        for &id in before.keys() {
            let net: f64 = report
                .transfers
                .iter()
                .filter(|row| row.species.as_str() == name)
                .map(|row| {
                    if row.receiver == id {
                        row.mass_kg
                    } else if row.donor == id {
                        -row.mass_kg
                    } else {
                        0.0
                    }
                })
                .sum();
            near(
                species_mass(&after.parcels()[&id], name),
                species_mass(&before[&id], name) + net,
                physical.exchange.mass_tolerance_kg,
            );
        }
    }
}

#[test]
fn g1_stiff_binary_window_uses_three_pair_calls_and_only_fine_transfers() {
    let d = card(2e-5);
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let right = ClosedGasParcel::new(gas(0.25).model().clone(), context.clone(), 0.004).unwrap();
    let initial = BTreeMap::from([(0, gas(0.75)), (1, right)]);
    let saved = initial.clone();
    let bounds = endpoints(&initial);
    let links = [link(7, [1, 0], &pairs)];
    let mut physical = solve();
    physical.exchange.maximum_species_fraction = 0.9;
    assert!(
        network::step_candidate(&initial, &links, 100.0, &bounds, physical, || false,).is_err()
    );
    let next = advance_adaptive_candidate(
        &initial,
        &links,
        100.0,
        &bounds,
        physical,
        accuracy(1e-8, 100.0),
        || false,
    )
    .unwrap();
    let decay = (-2.5_f64).exp();
    near(composition(&next.parcels()[&0]), 0.35 + 0.4 * decay, 3e-12);
    near(composition(&next.parcels()[&1]), 0.35 - 0.1 * decay, 3e-12);
    let report = next.report();
    assert_eq!(report.duration_s, 100.0);
    assert_eq!(report.accepted_intervals, 1);
    assert_eq!(report.accepted_substeps, 2);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.evaluations, 3);
    assert_eq!(report.pair_evaluations, 3);
    assert_eq!(report.last_substep.len(), 1);
    assert_eq!(report.last_substep[0].relaxation.duration_s, 50.0);
    // Counting the coarse proposal would double this actual inventory change.
    balances(&initial, &next, physical);
    assert_eq!(initial, saved);
}

#[test]
fn g1_tighter_binary_network_accuracy_refines_the_independent_chain_mode() {
    let d = card(2e-5);
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let links = [link(9, [1, 2], &pairs), link(2, [0, 1], &pairs)];
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.5)), (2, gas(0.25))]);
    let bounds = endpoints(&initial);
    let mut physical = solve();
    physical.exchange.maximum_species_fraction = 0.9;
    let exact = [
        0.5 + 0.25 * (-1.0_f64).exp(),
        0.5,
        0.5 - 0.25 * (-1.0_f64).exp(),
    ];
    let mut errors = Vec::new();
    for relative in [1e-3, 1e-5] {
        let next = advance_adaptive_candidate(
            &initial,
            &links,
            50.0,
            &bounds,
            physical,
            accuracy(relative, 50.0),
            || false,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(report.duration_s, 50.0);
        assert!(report.accepted_intervals > 1 && report.rejected_intervals > 0);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert_eq!(report.pair_evaluations, 3 * report.evaluations);
        assert_eq!(report.last_substep.len(), 3);
        assert!(report.maximum_species_error_ratio <= 1.0);
        assert!(report.maximum_temperature_error_ratio <= 1.0);
        errors.push(
            (0..3)
                .map(|i| (composition(&next.parcels()[&(i as u64)]) - exact[i]).abs())
                .fold(0.0_f64, f64::max),
        );
        balances(&initial, &next, physical);
    }
    assert!(
        errors[0] < 3e-4 && errors[1] < 0.1 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn g1_nonisothermal_source_refresh_refines_to_an_independent_network_ode() {
    let d = temperature_card();
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let initial: BTreeMap<_, _> = [(0, 0.75, 390.0), (1, 0.5, 400.0), (2, 0.25, 410.0)]
        .into_iter()
        .map(|(id, x, t)| {
            (
                id,
                ClosedGasParcel::new(gas(x).model().clone(), point(t, 100000.0), 0.001).unwrap(),
            )
        })
        .collect();
    let make = |id, nodes| {
        let mut edge = link(id, nodes, &pairs);
        edge.solve.maximum_temperature_difference_k = 50.0;
        edge.solve.maximum_pressure_difference_pa = 50_000.0;
        edge
    };
    let links = [make(3, [0, 1]), make(8, [1, 2]), make(2, [0, 2])];
    let bounds = endpoints(&initial);
    // Independent RK4 integration of simultaneous species/partial-enthalpy
    // rates, with analytic quadratic NASA inversion (8192 steps over 4 s).
    // These are continuum-in-time references, not another pair-sweep oracle.
    let exact = [
        (0.6979513894072313, 390.8109800553009),
        (0.5001289547192785, 399.88329994384014),
        (0.3045855915230852, 408.25146127694865),
    ];
    let mut errors = Vec::new();
    for (relative, temperature) in [(1e-3, 0.002), (2.5e-4, 0.000125)] {
        let mut requested = accuracy(relative, 4.0);
        requested.temporal.temperature_tolerance_k = temperature;
        let next =
            advance_adaptive_candidate(&initial, &links, 4.0, &bounds, solve(), requested, || {
                false
            })
            .unwrap();
        let report = next.report();
        assert!(report.rejected_intervals > 0);
        assert!(report.maximum_species_error_ratio <= 1.0);
        assert!(report.maximum_temperature_error_ratio <= 1.0);
        assert_eq!(report.pair_evaluations, 5 * report.evaluations);
        let error = (0..3)
            .map(|i| (next.parcels()[&(i as u64)].state().state().temperature - exact[i].1).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        for i in 0..3 {
            near(composition(&next.parcels()[&(i as u64)]), exact[i].0, 3e-6);
        }
        for stage in &report.last_substep {
            let flux = &stage.initial_flux;
            near(
                flux.pairs[0].diffusivity_m2_s,
                5e-8 * flux.face_temperature_k,
                1e-17,
            );
            assert!(
                flux.pairs[0]
                    .material
                    .query_point()
                    .contains(&("temperature".into(), flux.face_temperature_k))
            );
        }
        balances(&initial, &next, solve());
    }
    assert!(
        errors[0] < 0.006 && errors[1] < 0.35 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn g0_withdrawal_cutbacks_charge_actual_partial_sweeps_and_preserve_budgets() {
    let d = card(2e-5);
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let initial = BTreeMap::from([(0, gas(0.9)), (1, gas(0.1)), (2, gas(0.1))]);
    let bounds = endpoints(&initial);
    let links = [link(1, [0, 1], &pairs), link(2, [0, 2], &pairs)];
    let next = advance_adaptive_candidate(
        &initial,
        &links,
        25.0,
        &bounds,
        solve(),
        accuracy(0.1, 25.0),
        || false,
    )
    .unwrap();
    let report = next.report();
    // The initial coarse sweep fails at its second edge: two calls, not three.
    assert_eq!(report.rejected_intervals, 1);
    assert_eq!(report.accepted_intervals, 2);
    assert_eq!(report.evaluations, 7);
    assert_eq!(report.pair_evaluations, 20);
    balances(&initial, &next, solve());
    let mut limited = accuracy(0.1, 25.0);
    limited.maximum_pair_evaluations = 5;
    assert!(matches!(advance_adaptive_candidate(
        &initial, &links, 25.0, &bounds, solve(), limited, || false,
    ), Err(AdvanceError::PairWorkLimit {
        elapsed_s, accepted_intervals: 0,
        pair_evaluations: 5, maximum_pair_evaluations: 5,
    }) if elapsed_s == 0.0));
}

#[test]
fn g4_missing_sources_and_late_cancellation_leave_no_partial_network() {
    let d = card(2e-5);
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.5)), (2, gas(0.25))]);
    let saved = initial.clone();
    let bounds = endpoints(&initial);
    let requested = accuracy(1e-3, 0.1);
    let missing = [link(3, [0, 1], &pairs), link(8, [1, 2], &[])];
    assert!(matches!(advance_adaptive_candidate(
        &initial, &missing, 0.1, &bounds, solve(), requested, || false,
    ), Err(AdvanceError::BinaryStep {
        elapsed_s, step_s,
        error: BinaryDiffusionNetworkError::Face { id: 8, .. },
    }) if elapsed_s == 0.0 && step_s == 0.1));
    let links = [link(3, [0, 1], &pairs), link(8, [1, 2], &pairs)];
    let mut count = 0;
    advance_adaptive_candidate(&initial, &links, 0.1, &bounds, solve(), requested, || {
        count += 1;
        false
    })
    .unwrap();
    let mut calls = 0;
    assert_eq!(
        advance_adaptive_candidate(&initial, &links, 0.1, &bounds, solve(), requested, || {
            calls += 1;
            calls == count
        },),
        Err(AdvanceError::Cancelled)
    );
    assert_eq!(initial, saved);
}

#[test]
fn g0_empty_graphs_still_respect_the_complete_trial_budget() {
    let initial = BTreeMap::from([(0, gas(0.5))]);
    let saved = initial.clone();
    let mut requested = accuracy(1e-3, 1.0);
    requested.temporal.maximum_evaluations = 3;
    let one = advance_adaptive_candidate(
        &initial,
        &[],
        1.0,
        &endpoints(&initial),
        solve(),
        requested,
        || false,
    )
    .unwrap();
    assert_eq!(one.parcels(), &initial);
    assert_eq!(one.report().pair_evaluations, 0);
    assert_eq!(one.report().evaluations, 3);
    assert!(matches!(advance_adaptive_candidate(
        &initial, &[], 2.0, &endpoints(&initial), solve(), requested, || false,
    ), Err(AdvanceError::WorkLimit {
        elapsed_s, accepted_intervals: 1, evaluations: 3,
    }) if elapsed_s == 1.0));
    assert_eq!(initial, saved);
}
