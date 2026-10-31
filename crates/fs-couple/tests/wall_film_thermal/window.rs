//! Complete physical windows: heat totals, retries, wet/dry events and time refinement.
use super::*;
use fs_couple::finite_gas_thermal::wall_film::window::{
    WallFilmWindowError, WallFilmWindowReport, WallFilmWindowSolve,
};
use fs_material::phase::{UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

fn window(step: f64) -> WallFilmWindowSolve {
    WallFilmWindowSolve {
        maximum_step_s: step,
        minimum_step_s: 1e-7,
        maximum_attempts: 2000,
        energy_tolerance_j: 1e-6,
        mass_tolerance_kg: 1e-13,
    }
}
fn window_balance(
    input: UniformEnthalpyStepInput<'_>,
    gas: &ClosedGasParcel,
    film: GasWallFilmInventory,
    binding: &fs_material::evaporation::gas::GasLiquidBinding,
    gas_heat: f64,
    next: &UniformEnthalpyStep<WallFilmWindowReport>,
) {
    let report = &next.report;
    let wall_change = input.mass_kg
        * (next.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
    let gas_change = report.next_gas.mass_kg()
        * report.next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let film_change = film_energy(report.next_film, binding) - film_energy(film, binding);
    near(
        wall_change + gas_change + film_change,
        input.internal_heat_j + gas_heat + report.radiation_heat_j,
        1e-6,
    );
    near(
        wall_change,
        input.internal_heat_j + next.external_heat_j,
        1e-6,
    );
    near(
        report.next_film.mass_kg() + report.next_gas.species_mass_kg(&water()).unwrap(),
        film.mass_kg() + gas.species_mass_kg(&water()).unwrap(),
        1e-13,
    );
    near(
        report.next_gas.species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-17,
    );
    near(
        report.next_film.mass_kg() + report.next_gas.mass_kg(),
        film.mass_kg() + gas.mass_kg(),
        1e-13,
    );
    assert!(report.sum_absolute_energy_residual_j <= 1e-6);
    assert_eq!(report.samples.first().unwrap().time_s, 0.0);
    assert_eq!(report.samples.last().unwrap().time_s, input.duration_s);
    assert_eq!(report.samples.len(), report.accepted_intervals + 1);
    assert_eq!(
        report.attempted_intervals,
        report.accepted_intervals + report.rejected_intervals
    );
    assert!(
        report
            .samples
            .windows(2)
            .all(|samples| samples[0].time_s < samples[1].time_s)
    );
}

#[test]
fn uniform_sources_are_total_window_energies_and_final_partial_interval_keeps_its_clock() {
    let source = source();
    let gas = gas(&source, 300.0, 0.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let mut input = input(&curve, 300.0, 0.55);
    input.internal_heat_j = 4.0;
    let mut environment = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    environment.gas_heat_j = 1.0;
    let next = scope(|cx, _| environment.advance_window(cx, input, window(0.1)).unwrap());
    assert_eq!(next.report.accepted_intervals, 6);
    assert_eq!(next.report.rejected_intervals, 0);
    assert_eq!(next.report.next_film, GasWallFilmInventory::Dry);
    assert_eq!(next.report.condensed_mass_kg, 0.0);
    assert_eq!(next.report.evaporated_mass_kg, 0.0);
    assert_eq!(next.report.internal_body_heat_j, 4.0);
    assert_eq!(next.report.external_gas_heat_j, 1.0);
    near(
        next.state.temperature_k(),
        300.0 + 4.0 / (input.mass_kg * CP_WALL),
        1e-7,
    );
    let gas_capacity = 2.5 * 100000.0 * GAS_VOLUME / 300.0;
    near(
        next.report.next_gas.state().state().temperature,
        300.0 + 1.0 / gas_capacity,
        1e-7,
    );
    near(next.report.body_heat_apportionment_residual_j, 0.0, 1e-14);
    near(next.report.gas_heat_apportionment_residual_j, 0.0, 1e-14);
    window_balance(input, &gas, GasWallFilmInventory::Dry, &binding, 1.0, &next);
    assert_eq!(gas, original_gas);
}

#[test]
fn donor_retries_keep_original_heat_totals_and_one_whole_window_balance_allowance() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let mut input = input(&curve, 300.0, 0.2);
    input.internal_heat_j = 0.5;
    let mut environment = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    environment.solve.phase.maximum_inventory_fraction = 0.01;
    environment.gas_heat_j = 0.2;
    let next = scope(|cx, _| environment.advance_window(cx, input, window(0.2)).unwrap());
    assert!(next.report.rejected_intervals > 0);
    assert!(next.report.accepted_intervals > 1);
    assert!(next.report.next_film.mass_kg() > 0.0);
    assert_eq!(next.report.internal_body_heat_j, 0.5);
    assert_eq!(next.report.external_gas_heat_j, 0.2);
    window_balance(input, &gas, GasWallFilmInventory::Dry, &binding, 0.2, &next);
    // The same repeated solves cannot renew the global energy allowance.
    assert!(next.report.sum_absolute_energy_residual_j > 0.0);
    let mut tight = window(0.2);
    tight.energy_tolerance_j = 0.5 * next.report.sum_absolute_energy_residual_j;
    assert!(matches!(
        scope(|cx, _| environment.advance_window(cx, input, tight)),
        Err(WallFilmWindowError::EnergyBalance { .. })
    ));
    let mut exhausted = window(0.2);
    exhausted.maximum_attempts = 1;
    assert!(matches!(
        scope(|cx, _| environment.advance_window(cx, input, exhausted)),
        Err(WallFilmWindowError::WorkLimit { accepted_time_s, .. }) if accepted_time_s == 0.0
    ));
    assert_eq!(gas, original_gas);
    assert_eq!(environment.film, GasWallFilmInventory::Dry);
}

#[test]
fn the_same_wall_and_water_survive_onset_complete_dryout_and_recondensation() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let cold_input = input(&curve, 300.0, 0.1);
    let (cold, hot, recooled) = scope(|cx, _| {
        let cold = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body)
            .advance_window(cx, cold_input, window(0.025))
            .unwrap();
        assert!(cold.report.next_film.mass_kg() > 0.0);
        assert!(cold.report.samples.iter().any(|sample| sample.film_formed));
        window_balance(
            cold_input,
            &gas,
            GasWallFilmInventory::Dry,
            &binding,
            0.0,
            &cold,
        );

        let warm_input = UniformEnthalpyStepInput {
            initial: cold.state,
            duration_s: 2.0,
            internal_heat_j: 40.0,
            ..cold_input
        };
        let hot = environment(
            &cold.report.next_gas,
            &source,
            cold.report.next_film,
            &binding,
            &body,
        )
        .advance_window(cx, warm_input, window(0.025))
        .unwrap();
        assert_eq!(hot.report.next_film, GasWallFilmInventory::Dry);
        assert!(hot.report.evaporated_mass_kg >= cold.report.next_film.mass_kg());
        assert!(
            hot.report
                .samples
                .iter()
                .any(|sample| sample.dry_out_time_s.is_some())
        );
        let event_index = hot
            .report
            .samples
            .iter()
            .position(|sample| sample.dry_out_time_s.is_some())
            .unwrap();
        let event = hot.report.samples[event_index].dry_out_time_s.unwrap();
        assert!(event > 0.0 && event <= hot.report.samples[event_index].time_s);
        assert!(
            hot.report.samples[event_index..]
                .iter()
                .all(|sample| sample.film_mass_kg == 0.0)
        );
        window_balance(
            warm_input,
            &cold.report.next_gas,
            cold.report.next_film,
            &binding,
            0.0,
            &hot,
        );

        let cooling_input = UniformEnthalpyStepInput {
            initial: hot.state,
            duration_s: 2.0,
            internal_heat_j: -40.0,
            ..cold_input
        };
        let recooled = environment(
            &hot.report.next_gas,
            &source,
            hot.report.next_film,
            &binding,
            &body,
        )
        .advance_window(cx, cooling_input, window(0.025))
        .unwrap();
        assert!(recooled.report.next_film.mass_kg() > 0.0);
        assert!(
            recooled
                .report
                .samples
                .iter()
                .any(|sample| sample.film_formed)
        );
        assert!(recooled.state.temperature_k() < hot.state.temperature_k());
        window_balance(
            cooling_input,
            &hot.report.next_gas,
            hot.report.next_film,
            &binding,
            0.0,
            &recooled,
        );
        (cold, hot, recooled)
    });
    for next in [&cold, &hot, &recooled] {
        near(
            next.report.next_gas.species_mass_kg(&water()).unwrap()
                + next.report.next_film.mass_kg(),
            gas.species_mass_kg(&water()).unwrap(),
            1e-13,
        );
    }
    assert_eq!(gas, original_gas);
}

#[test]
fn complete_heat_phase_and_wall_clock_has_a_first_order_refinement_ladder() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let curve = curve();
    let body = body(20.0);
    let input = input(&curve, 300.0, 0.2);
    let environment = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    let endpoints: Vec<_> = [0.04, 0.02, 0.01, 0.005]
        .into_iter()
        .map(|dt| {
            let next = scope(|cx, _| environment.advance_window(cx, input, window(dt)).unwrap());
            assert!(next.report.next_film.mass_kg() > 0.0);
            assert!(next.report.solid_to_gas_heat_j < 0.0);
            assert_eq!(
                next.report.rejected_intervals, 0,
                "refinement rung must stay live"
            );
            assert_eq!(
                next.report.accepted_intervals,
                (input.duration_s / dt).round() as usize
            );
            window_balance(input, &gas, GasWallFilmInventory::Dry, &binding, 0.0, &next);
            [
                next.state.temperature_k(),
                next.report.next_gas.state().state().temperature,
                next.report.next_film.mass_kg() * LATENT / (input.mass_kg * CP_WALL),
            ]
        })
        .collect();
    let differences: Vec<_> = endpoints
        .windows(2)
        .map(|ends| {
            ends[0]
                .iter()
                .zip(ends[1])
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
        })
        .collect();
    assert!(differences.iter().all(|value| *value > 0.0));
    for pair in differences.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (1.4..2.8).contains(&ratio),
            "full-clock differences={differences:?}, ratio={ratio}"
        );
    }
}

#[test]
fn late_work_refusal_and_cancellation_leave_the_entire_input_window_unpublished() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let input = input(&curve, 300.0, 0.1);
    let original_wall = input.initial;
    let environment = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    let mut short = window(0.025);
    short.maximum_attempts = 1;
    assert!(matches!(
        scope(|cx, _| environment.advance_window(cx, input, short)),
        Err(WallFilmWindowError::WorkLimit { attempted_intervals: 1, accepted_time_s }) if accepted_time_s == 0.025
    ));
    scope(|cx, gate| {
        gate.request();
        assert!(matches!(
            environment.advance_window(cx, input, window(0.025)),
            Err(WallFilmWindowError::Cancelled)
        ));
    });
    assert_eq!(gas, original_gas);
    assert_eq!(input.initial, original_wall);
    assert_eq!(environment.film, GasWallFilmInventory::Dry);
}
