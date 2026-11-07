//! G1/G3/G4: thermal reacting cycles, original latent inventory and global refusal.
use super::*;
use crate::material::compression::history::{
    CompressionKnot, MaterialCompressionHistoryError as HistoryError,
};
use crate::material::compression::history::thermal::MaterialThermalCompressionHistorySolve as HistorySolve;
use crate::material::compression::thermal::adaptive::{
    AdaptiveMaterialThermalCompressionError as AdaptiveError,
    AdaptiveMaterialThermalCompressionSolve,
};

fn knot(time: f64, volume: f64) -> CompressionKnot {
    CompressionKnot { elapsed_time: Time::new(time), volume: Volume::new(volume) }
}
fn fixed(substeps: usize, maximum_substeps: usize) -> HistorySolve {
    HistorySolve::Fixed { controls: path(substeps), maximum_substeps }
}
fn adaptive(tolerance: f64) -> AdaptiveMaterialThermalCompressionSolve {
    AdaptiveMaterialThermalCompressionSolve {
        volume: path(1).volume,
        absolute_temperature_tolerance_k: tolerance,
        absolute_amount_tolerance_mol: 1e-10,
        absolute_work_tolerance_j: 0.25 * tolerance,
        absolute_heat_tolerance_j: 0.25 * tolerance,
        absolute_body_enthalpy_tolerance_j: 0.25 * tolerance,
        relative_tolerance: 0.0,
        initial_step_s: 0.05,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.05,
        maximum_intervals: 4096,
        maximum_evaluations: 12288,
    }
}

// Independent exact first law on every log-volume leg:
// T'=(a-b)T+b*Tb+c*exp(-4t), a=ln(Vold/Vnew)/(3*dt), b=G/(n*3R).
// Chemical time is global; turns/dwells never reset the reactant inventory.
fn bath_history_reference(
    knots: &[CompressionKnot], bath: f64, conductance: f64, release: f64,
) -> [f64; 3] {
    let capacity = 0.01 * 3.0 * R;
    let b = conductance / capacity;
    let (mut temperature, mut work, mut heat) = (400.0, 0.0, 0.0);
    for pair in knots.windows(2) {
        let start = pair[0].elapsed_time.value();
        let duration = pair[1].elapsed_time.value() - start;
        let a = (pair[0].volume.value() / pair[1].volume.value()).ln() / (3.0 * duration);
        let lambda = a - b;
        let c = 4.0 * release / (3.0 * R) * (-4.0 * start).exp();
        assert!(lambda != 0.0 && lambda + 4.0 != 0.0);
        let exponential = (lambda * duration).exp();
        let chemical = (-4.0 * duration).exp();
        let integral_exponential = (lambda * duration).exp_m1() / lambda;
        let integral_temperature = temperature * integral_exponential
            + b * bath * (integral_exponential - duration) / lambda
            + c * (integral_exponential - (1.0 - chemical) / 4.0) / (lambda + 4.0);
        temperature = temperature * exponential
            + b * bath * integral_exponential
            + c * (exponential - chemical) / (lambda + 4.0);
        work += capacity * a * integral_temperature;
        heat += conductance * (bath * duration - integral_temperature);
    }
    [temperature, work, heat]
}

#[test]
fn held_bath_closed_reacting_cycle_preserves_turns_dwell_and_work_heat() {
    let knots = [
        knot(0.0, 0.01), knot(0.0625, 0.005),
        knot(0.09375, 0.005), knot(0.15625, 0.01),
    ];
    let expected = bath_history_reference(&knots, 300.0, 0.2, 12_000.0);
    let rigid = bath_reference(0.01, 0.15625, 300.0, 0.2, 12_000.0);
    let (bank, mechanism, initial) = fixture(-12_000.0, 2000.0);
    let mut reactor = MaterialThermalCompressionReactor::new(
        bank, mechanism, kinetics(), || false,
    ).unwrap();
    let workspace = reactor.workspace_bytes();
    let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(300.0)).unwrap();
    let mut errors = Vec::new();
    for substeps in [4, 8, 16] {
        let mut parcel = initial.clone();
        let report = reactor.advance_history_with_heat_bath(
            &mut parcel, &knots, bath, fixed(substeps, 2 * substeps + 1), || false,
        ).unwrap();
        let physical = report.physical;
        errors.push((parcel.state().state().temperature - expected[0]).abs());
        let stored = parcel.mass_kg() * (
            parcel.state().specific_internal_energy_j_kg()
                - initial.state().specific_internal_energy_j_kg()
        );
        close(stored, physical.internal_energy_change_j, 5e-6);
        close(stored, physical.boundary_work_on_gas_j + physical.heat_input_j
            + physical.energy_residual_j, 5e-6);
        close(physical.energy_residual_j, 0.0, 2e-6);
        close(physical.thermal_partner_entropy_change_j_k,
            -physical.heat_input_j / 300.0, 2e-10);
        assert!(physical.entropy_production_j_k >= -1e-7);
        assert!(physical.boundary_work_on_gas_j < -5.0 && physical.heat_input_j < -8.0);
        assert!((parcel.state().state().temperature - rigid[0]).abs() > 25.0);
        close(parcel.model().thermodynamics().mole_composition().fractions()[1],
            1.0 - (-4.0_f64 * 0.15625).exp(), 5e-7);
        assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
        assert_eq!(parcel.volume_m3().to_bits(), initial.volume_m3().to_bits());
        assert_eq!(parcel.point().axes()["sample_id"], 7.0);
        assert_eq!(parcel.model().transport_materials(), reactor.bank().transport_materials());
        assert_eq!(report.completed_segments, 3);
        assert_eq!(report.evaluations, 2 * substeps + 1);
        assert_eq!(physical.substeps, report.evaluations);
        assert_eq!(report.total_chemical_accepted_steps, physical.chemical_accepted_steps);
        assert_eq!(report.maximum_accepted_error_ratios, None);
        assert_eq!(report.suggested_step_s, None);
        assert_eq!(report.elapsed_time_error_s, 0.0);
        assert_eq!(physical.elapsed_time, Time::new(0.15625));
        parcel.state().state().validate_coherence().unwrap();
        if substeps == 16 {
            close(parcel.state().state().temperature, expected[0], 0.003);
            close(physical.boundary_work_on_gas_j, expected[1], 0.0015);
            close(physical.heat_input_j, expected[2], 0.0007);
        }
    }
    assert!(errors.windows(2).all(|pair| pair[0] / pair[1] > 3.7), "{errors:?}");
    assert_eq!(reactor.workspace_bytes(), workspace);
}

// Independent variable gas Cv and finite body Cp balances, integrated without
// any native caloric inversion, chemical stage, volume map or history method.
fn finite_body_history_reference(knots: &[CompressionKnot]) -> [f64; 4] {
    let mut state = [400.0, 800.0, 0.0, 0.0];
    for pair in knots.windows(2) {
        let start = pair[0].elapsed_time.value();
        let duration = pair[1].elapsed_time.value() - start;
        let compression_rate = (pair[0].volume.value() / pair[1].volume.value()).ln() / duration;
        let rate = |time: f64, value: [f64; 4]| {
            let heat = 2.0 * (value[1] - value[0]);
            let work = 0.01 * R * value[0] * compression_rate;
            [
                (work + 0.01 * 24_000.0 * 4.0 * (-4.0 * time).exp() + heat)
                    / (0.01 * R * (3.0 + 0.004 * value[0])),
                -heat / (0.001 * (60.0 + 0.2 * value[1])),
                work,
                heat,
            ]
        };
        let dt = duration / 4096.0;
        for index in 0..4096 {
            let time = start + index as f64 * dt;
            let shift = |rate: [f64; 4], fraction: f64| {
                std::array::from_fn(|i| state[i] + fraction * dt * rate[i])
            };
            let a = rate(time, state);
            let b = rate(time + 0.5 * dt, shift(a, 0.5));
            let c = rate(time + 0.5 * dt, shift(b, 0.5));
            let d = rate(time + dt, shift(c, 1.0));
            for i in 0..4 {
                state[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
            }
        }
    }
    state
}

#[test]
fn finite_source_body_cycles_retain_variable_capacities_density_and_two_work_accounts() {
    let knots = [
        knot(0.0, 0.01), knot(0.05, 0.005),
        knot(0.075, 0.005), knot(0.125, 0.01),
    ];
    let expected = finite_body_history_reference(&knots);
    let (bank, mechanism, mut parcel) = variable_fixture(24_000.0);
    let initial = parcel.clone();
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 800.0, 2.0);
    let mut reactor = MaterialThermalCompressionReactor::new(
        bank, mechanism, kinetics(), || false,
    ).unwrap();
    let report = reactor.advance_history_with_finite_body(
        &mut parcel, &knots, body, fixed(64, 129), || false,
    ).unwrap();
    let physical = report.physical;
    close(parcel.state().state().temperature, expected[0], 0.012);
    close(physical.body.state().temperature_k(), expected[1], 0.012);
    close(physical.gas.boundary_work_on_gas_j, expected[2], 0.006);
    close(physical.gas.heat_input_j, expected[3], 0.006);
    let gas_u = parcel.mass_kg() * (
        parcel.state().specific_internal_energy_j_kg()
            - initial.state().specific_internal_energy_j_kg()
    );
    close(gas_u, physical.gas.internal_energy_change_j, 2e-5);
    close(physical.body_enthalpy_change_j, -physical.gas.heat_input_j, 2e-5);
    close(gas_u + physical.body_enthalpy_change_j,
        physical.gas.boundary_work_on_gas_j, 2e-5);
    close(gas_u + physical.body_internal_energy_change_j + physical.body_boundary_work_j,
        physical.gas.boundary_work_on_gas_j, 2e-5);
    close(physical.combined_heat_residual_j, 0.0, 2e-5);
    close(physical.combined_energy_residual_j, 0.0, 2e-5);
    // Authored rho(T) = 1040 - 0.2*T on the same variable-Cp body source.
    let expected_body_work = 2e5 * 0.001
        * ((1040.0 - 0.2 * expected[1]).recip() - (1040.0_f64 - 0.2 * 800.0).recip());
    close(physical.body_boundary_work_j, expected_body_work, 2e-6);
    close(physical.body_entropy_change_j_k,
        body.mass().value() * curve.specific_entropy_change_j_kg_k(
            body.state().specific_enthalpy_j_kg(),
            physical.body.state().specific_enthalpy_j_kg(),
        ).unwrap(), 3e-8);
    assert_eq!(physical.body.mass(), body.mass());
    assert_eq!(physical.body.pressure(), body.pressure());
    assert_eq!(physical.body.conductance(), body.conductance());
    assert!(std::ptr::eq(physical.body.curve(), &curve));
    assert_eq!(parcel.volume_m3(), initial.volume_m3());
    assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
    assert_eq!(report.completed_segments, 3);
    assert_eq!(report.evaluations, 129);
    parcel.state().state().validate_coherence().unwrap();
}

fn latent_curve(bank: &ThermochemicalGasModel) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, t, fraction, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: fraction,
        bulk_density_kg_m3: rho,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        bank.thermodynamics().components()[0].regions()[0].card().content_hash(),
        vec![
            knot(0.0, 200.0, 0.0, 1100.0),
            knot(20_000.0, 400.0, 0.0, 1000.0),
            knot(120_000.0, 400.0, 1.0, 800.0),
            knot(200_000.0, 1200.0, 1.0, 700.0),
        ],
    ).unwrap()
}

#[test]
fn latent_body_inventory_survives_turns_and_reverses_heat_at_fixed_temperature() {
    // End before ln(2)/4: the manufactured irreversible equal-caloric source
    // remains on the physically increasing side of its mixing entropy.
    let knots = [
        knot(0.0, 0.01), knot(0.05, 0.005), knot(0.075, 0.005),
        knot(0.125, 0.01), knot(0.15625, 0.01),
    ];
    let expected = bath_history_reference(&knots, 400.0, 1.0, 0.0);
    let (bank, mechanism, initial) = fixture(0.0, 2000.0);
    let curve = latent_curve(&bank);
    let body = FiniteHeatBody::new(
        &curve, Mass::new(0.001), Pressure::new(2e5), 70_000.0,
        ThermalConductance::new(1.0),
    ).unwrap();
    let mut reactor = MaterialThermalCompressionReactor::new(
        bank, mechanism, kinetics(), || false,
    ).unwrap();
    let mut prefix = initial.clone();
    let prefix_report = reactor.advance_history_with_finite_body(
        &mut prefix, &knots[..4], body, fixed(64, 130), || false,
    ).unwrap();
    let mut refined = adaptive(1e6);
    refined.absolute_amount_tolerance_mol = 1e6;
    refined.absolute_body_enthalpy_tolerance_j = 1e-5;
    for solve in [fixed(64, 130), HistorySolve::Adaptive(refined)] {
        let mut parcel = initial.clone();
        let report = reactor.advance_history_with_finite_body(
            &mut parcel, &knots, body, solve, || false,
        ).unwrap();
        let physical = report.physical;
        close(parcel.state().state().temperature, expected[0], 0.005);
        close(physical.gas.boundary_work_on_gas_j, expected[1], 0.002);
        close(physical.gas.heat_input_j, expected[2], 0.002);
        assert_eq!(physical.body.state().temperature_k(), 400.0);
        let fraction = physical.body.state().liquid_mass_fraction();
        close(fraction, 0.5 - expected[2] / 100.0, 2e-5);
        assert!(fraction > 0.5 && fraction < 1.0);
        // The last dwell warms the expanded gas and partly solidifies the SAME
        // finite body. Comparing temperatures alone would lose this history.
        assert!(fraction < prefix_report.physical.body.state().liquid_mass_fraction() - 0.004);
        assert!(physical.body.state().specific_enthalpy_j_kg()
            < prefix_report.physical.body.state().specific_enthalpy_j_kg());
        close(physical.body_enthalpy_change_j, -physical.gas.heat_input_j, 2e-5);
        close(physical.body_internal_energy_change_j + physical.body_boundary_work_j,
            physical.body_enthalpy_change_j, 2e-5);
        close(physical.body_entropy_change_j_k,
            physical.body_enthalpy_change_j / 400.0, 3e-8);
        close(physical.combined_heat_residual_j, 0.0, 2e-5);
        close(physical.combined_energy_residual_j, 0.0, 2e-5);
        assert_eq!(physical.body.mass(), body.mass());
        assert!(std::ptr::eq(physical.body.curve(), &curve));
        assert_eq!(parcel.volume_m3().to_bits(), initial.volume_m3().to_bits());
        assert_eq!(report.completed_segments, 4);
        if let Some(errors) = report.maximum_accepted_error_ratios {
            assert!(errors.into_iter().all(|error| error <= 1.0));
            assert!(errors[..4].iter().all(|&error| error < 1e-3));
            assert!(report.rejected_adaptive_intervals > 0);
            assert_eq!(report.evaluations,
                3 * (report.accepted_adaptive_intervals + report.rejected_adaptive_intervals) - 4);
            assert_eq!(physical.gas.substeps, 2 * report.accepted_adaptive_intervals - 2);
            assert!(report.total_chemical_accepted_steps > physical.gas.chemical_accepted_steps);
        } else {
            assert_eq!(report.evaluations, 130);
        }
    }
}

#[test]
fn global_thermal_budgets_include_dwells_and_roll_back_both_endpoints_then_replay() {
    let knots = [
        knot(0.0, 0.01), knot(0.03125, 0.005),
        knot(0.0625, 0.01), knot(0.09375, 0.01),
    ];
    let (bank, mechanism, initial) = fixture(-12_000.0, 2000.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 800.0, 2.0);
    let original_body = body;
    let mut reactor = MaterialThermalCompressionReactor::new(
        bank, mechanism, kinetics(), || false,
    ).unwrap();
    let mut parcel = initial.clone();
    assert!(matches!(reactor.advance_history_with_finite_body(
        &mut parcel, &knots, body, fixed(4, 8), || false,
    ), Err(HistoryError::WorkLimit {
        completed_segments: 2, elapsed_time_s, evaluations: 8,
    }) if elapsed_time_s == 0.0625));
    assert_eq!(parcel, initial);
    assert_eq!(body, original_body);

    let mut accepted = initial.clone();
    let mut polls = 0;
    let expected = reactor.advance_history_with_finite_body(
        &mut accepted, &knots, body, fixed(4, 9), || {
            polls += 1;
            false
        },
    ).unwrap();
    assert_eq!(expected.evaluations, 9);
    for stop in [polls / 2, polls] {
        let mut parcel = initial.clone();
        let mut calls = 0;
        assert_eq!(reactor.advance_history_with_finite_body(
            &mut parcel, &knots, body, fixed(4, 9), || {
                calls += 1;
                calls >= stop
            },
        ).unwrap_err(), HistoryError::Cancelled);
        assert_eq!(parcel, initial);
        assert_eq!(body, original_body);
        let retry = reactor.advance_history_with_finite_body(
            &mut parcel, &knots, body, fixed(4, 9), || false,
        ).unwrap();
        assert_eq!(parcel, accepted);
        assert_eq!(retry, expected);
    }

    let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(700.0)).unwrap();
    let mut solve = adaptive(1e6);
    solve.absolute_amount_tolerance_mol = 1e6;
    solve.initial_step_s = 0.03125;
    solve.maximum_step_s = 0.03125;
    solve.maximum_intervals = 3;
    solve.maximum_evaluations = 7;
    let mut parcel = initial.clone();
    let report = reactor.advance_history_with_heat_bath(
        &mut parcel, &knots, bath, HistorySolve::Adaptive(solve), || false,
    ).unwrap();
    assert_eq!(report.evaluations, 7); // Three + three + final one-evaluation dwell.
    assert_eq!(report.accepted_adaptive_intervals, 3);
    assert_eq!(report.rejected_adaptive_intervals, 0);
    assert_eq!(report.physical.substeps, 5);
    assert_eq!(report.suggested_step_s, Some(0.03125));
    for (evaluations, intervals) in [(6, 3), (7, 2)] {
        let mut parcel = initial.clone();
        let mut refused = solve;
        refused.maximum_evaluations = evaluations;
        refused.maximum_intervals = intervals;
        assert!(matches!(reactor.advance_history_with_heat_bath(
            &mut parcel, &knots, bath, HistorySolve::Adaptive(refused), || false,
        ), Err(HistoryError::WorkLimit {
            completed_segments: 2, elapsed_time_s, evaluations: 6,
        }) if elapsed_time_s == 0.0625));
        assert_eq!(parcel, initial);
    }
}

#[test]
fn adaptive_thermal_errors_report_global_elapsed_time_and_spent_prefix_work() {
    let knots = [knot(0.0, 0.01), knot(0.03125, 0.01), knot(0.0625, 0.005)];
    let (bank, mechanism, initial) = fixture(-12_000.0, 2000.0);
    let mut reactor = MaterialThermalCompressionReactor::new(
        bank, mechanism, kinetics(), || false,
    ).unwrap();
    let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(300.0)).unwrap();
    let mut solve = adaptive(1e-10);
    solve.initial_step_s = 0.03125;
    solve.maximum_step_s = 0.03125;
    solve.maximum_evaluations = 4;
    let mut parcel = initial.clone();
    assert!(matches!(reactor.advance_history_with_heat_bath(
        &mut parcel, &knots, bath, HistorySolve::Adaptive(solve), || false,
    ), Err(HistoryError::ThermalAdaptive {
        segment: 1, start_time_s, end_time_s,
        source: AdaptiveError::WorkLimit {
            elapsed_time_s, accepted_intervals: 1, rejected_intervals: 1, evaluations: 4,
        },
    }) if start_time_s == 0.03125 && end_time_s == 0.0625 && elapsed_time_s == 0.03125));
    assert_eq!(parcel, initial);
    solve.minimum_step_s = 0.03125;
    solve.maximum_evaluations = 100;
    assert!(matches!(reactor.advance_history_with_heat_bath(
        &mut parcel, &knots, bath, HistorySolve::Adaptive(solve), || false,
    ), Err(HistoryError::ThermalAdaptive {
        segment: 1, start_time_s, end_time_s,
        source: AdaptiveError::Accuracy { elapsed_time_s, step_s, .. },
    }) if start_time_s == 0.03125 && end_time_s == 0.0625
        && elapsed_time_s == 0.03125 && step_s == 0.03125));
    assert_eq!(parcel, initial);
}

#[test]
fn thermal_history_source_refusal_after_a_feasible_prefix_and_invalid_knots_are_atomic() {
    let (bank, mechanism, initial) = fixture(0.0, 450.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 300.0, 0.05);
    let mut reactor = MaterialThermalCompressionReactor::new(
        bank, mechanism, kinetics(), || false,
    ).unwrap();
    let knots = [
        knot(0.0, 0.01), knot(0.02, 0.01),
        knot(0.04, 0.009), knot(0.1, 0.005),
    ];
    let mut prefix = initial.clone();
    let accepted = reactor.advance_history_with_finite_body(
        &mut prefix, &knots[..3], body, fixed(8, 17), || false,
    ).unwrap();
    assert!(prefix.state().state().temperature > 400.0);
    assert!(accepted.physical.body.state().temperature_k() > body.state().temperature_k());
    let mut parcel = initial.clone();
    assert!(matches!(reactor.advance_history_with_finite_body(
        &mut parcel, &knots, body, fixed(8, 17), || false,
    ), Err(HistoryError::Segment { segment: 2, .. })));
    assert_eq!(parcel, initial);
    let retry = reactor.advance_history_with_finite_body(
        &mut parcel, &knots[..3], body, fixed(8, 17), || false,
    ).unwrap();
    assert_eq!(parcel, prefix);
    assert_eq!(retry, accepted);
    let mut parcel = initial.clone();
    assert!(matches!(reactor.advance_history_with_finite_body(
        &mut parcel, &[knot(0.0, 0.02), knot(0.1, 0.005)],
        body, fixed(8, 17), || false,
    ), Err(HistoryError::Knot { index: 0, .. })));
    assert_eq!(parcel, initial);
}

#[test]
fn thermal_dwells_report_the_exact_independent_knot_roundoff() {
    let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
    let mut reactor = MaterialThermalCompressionReactor::new(
        bank, mechanism, kinetics(), || false,
    ).unwrap();
    let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(300.0)).unwrap();
    let knots = [
        knot(0.0, 0.01), knot(0.025, 0.01), knot(0.1, 0.01), knot(0.125, 0.01),
    ];
    let report = reactor.advance_history_with_heat_bath(
        &mut parcel, &knots, bath, fixed(4096, 3), || false,
    ).unwrap();
    // Exact real sum of the three represented f64 subtractions is 1/8 + 2^-57.
    assert_eq!(report.elapsed_time_error_s, 2.0_f64.powi(-57));
    assert_eq!(report.physical.elapsed_time, Time::new(0.125));
    assert_eq!(report.evaluations, 3);
    assert_eq!(report.physical.boundary_work_on_gas_j, 0.0);
    assert!(report.physical.heat_input_j < 0.0);
    close(parcel.model().thermodynamics().mole_composition().fractions()[1],
        1.0 - (-0.5_f64).exp(), 5e-7);
}
