//! Coupled temporal accuracy, birth exposure and physical inventory scales.
use super::super::adaptive::{GasCoatingAdaptiveError, advance_coatings_adaptive_with_gas};
use super::*;
use crate::sublimation::gas::timeline::adaptive::GasTimelineAdaptiveSolve;

fn accuracy(tolerance_k: f64) -> GasTimelineAdaptiveSolve {
    GasTimelineAdaptiveSolve {
        absolute_mass_tolerance_kg: tolerance_k * 1e-6,
        absolute_temperature_tolerance_k: tolerance_k,
        absolute_specific_enthalpy_tolerance_j_kg: 800.0 * tolerance_k,
        relative_tolerance: 0.0,
        minimum_dt_s: 1e-10,
        maximum_trials: 20000,
    }
}

// Independent simultaneous ODE: two supported films, a thermal contact and an
// ideal carrier/vapor gas. No production phase step, heat solver or adaptive
// controller is called here. States are T1,T2,Tg,m1,m2,m_vapor.
fn reference(gas: &ClosedGasParcel, duration: f64) -> [f64; 6] {
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let ra = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let carrier = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap();
    let vapor = gas.species_mass_kg(&fixture::species()).unwrap();
    let rates = |y: [f64; 6]| {
        let pressure = y[5] * r * y[2] / 0.01;
        let incoming = 0.1 * 9e-5 * pressure / (2.0 * core::f64::consts::PI * r * y[2]).sqrt();
        let mut derivative = [0.0; 6];
        let mut gas_heat = 0.0;
        for i in 0..2 {
            let saturation = 1000.0 * (200000.0 / r * (1.0 / 400.0 - 1.0 / y[i])).exp();
            let outgoing =
                0.1 * 9e-5 * saturation / (2.0 * core::f64::consts::PI * r * y[i]).sqrt();
            let into_gas = outgoing - incoming;
            let stream = outgoing * 3.5 * r * y[i] - incoming * 3.5 * r * y[2];
            let heat = if i == 0 {
                10.0 - 2.0 * (y[0] - y[1])
            } else {
                2.0 * (y[0] - y[1])
            };
            let deposit_h = 4.0 * r * y[i] - 200000.0;
            derivative[i] = (heat - stream + deposit_h * into_gas) / (8.0 + y[3 + i] * 4.0 * r);
            derivative[3 + i] = -into_gas;
            derivative[5] += into_gas;
            gas_heat += stream;
        }
        derivative[2] =
            (gas_heat - 3.0 * r * y[2] * derivative[5]) / (carrier * 2.5 * ra + y[5] * 3.0 * r);
        derivative
    };
    let mut y = [300.0, 315.0, 350.0, 1e-4, 1e-4, vapor];
    let dt = duration / 16384.0;
    for _ in 0..16384 {
        let a = rates(y);
        let b = rates(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = rates(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = rates(std::array::from_fn(|i| y[i] + dt * c[i]));
        y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn complete_coating_heat_and_finite_gas_refine_against_simultaneous_rates() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate { datum: 0.0 };
    let law = fixture::phase();
    let mut second = cell(2, substrate, 0.0);
    second.specific_h_j_kg = 800.0 * 15.0;
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &bindings,
        vec![
            birth(0.0, cell(1, substrate, 10.0)),
            birth(0.0, second),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1e-5,
                    conductance_w_k: 2.0,
                }),
            },
        ],
        [1e-4; 2],
        100000,
    );
    let selection = selected(&bindings);
    let expected = reference(&gas, 1.0);
    let run = |tolerance| {
        advance_coatings_adaptive_with_gas(
            &original,
            &selection,
            None,
            &gas,
            &catalog,
            1.0,
            heat(1.0),
            controls(),
            accuracy(tolerance),
            || Ok(()),
        )
        .unwrap()
    };
    let loose = run(0.02);
    let tight = run(1e-5);
    let error = |candidate: &GasCoatingTimelineCandidate<'_, '_>| {
        let actual = [
            candidate.coating().temperature_k(CellId(1)).unwrap(),
            candidate.coating().temperature_k(CellId(2)).unwrap(),
            candidate.gas().state().state().temperature,
            candidate.coating().coatings()[&CellId(1)].mass_kg,
            candidate.coating().coatings()[&CellId(2)].mass_kg,
            candidate
                .gas()
                .species_mass_kg(&fixture::species())
                .unwrap(),
        ];
        actual
            .into_iter()
            .zip(expected)
            .enumerate()
            .map(|(i, (a, b))| (a - b).abs() * if i < 3 { 1.0 } else { 1e6 })
            .fold(0.0, f64::max)
    };
    let errors = [error(loose.endpoint()), error(tight.endpoint())];
    assert!(
        errors[1] < 0.15 * errors[0] && errors[1] < 0.05,
        "{errors:?}"
    );
    assert!(tight.audit().accepted_intervals > loose.audit().accepted_intervals);
    assert!(tight.audit().rejected_trials > 0);
    assert!(tight.audit().maximum_accepted_error_ratio <= 1.0);
    assert!(tight.audit().thermal_substeps_evaluated > tight.endpoint().audit().thermal.substeps);
    close(tight.endpoint().audit().thermal.external_heat_j, 10.0, 1e-9);
    close(tight.audit().maximum_species_mass_residual_kg, 0.0, 1e-14);
    close(tight.audit().solid_stream_energy_residual_j, 0.0, 1e-7);
    close(tight.audit().gas_stream_energy_residual_j, 0.0, 1e-7);
    check(tight.endpoint());
    assert!(original.cells().is_empty());
}

#[test]
fn bare_coatings_preserve_birth_exposure_and_canonical_replay() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate { datum: -100000.0 };
    let law = fixture::phase();
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &bindings,
        vec![
            birth(0.01, cell(1, substrate, 2.0)),
            birth(0.02, cell(2, substrate, 0.0)),
            TimedEvent {
                time_s: 0.02,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1e-4,
                    conductance_w_k: 0.0,
                }),
            },
        ],
        [0.0; 2],
        20000,
    );
    let selection = selected(&bindings);
    let run = |selection: &[GasCoatingTimelineSurface<'_>]| {
        advance_coatings_adaptive_with_gas(
            &original,
            selection,
            None,
            &gas,
            &catalog,
            0.03,
            heat(0.02),
            controls(),
            accuracy(0.001),
            || Ok(()),
        )
        .unwrap()
    };
    let accepted = run(&selection);
    let replay = run(&[selection[1], selection[0]]);
    assert_eq!(
        accepted.endpoint().coating().cells(),
        replay.endpoint().coating().cells()
    );
    assert_eq!(
        accepted.endpoint().coating().coatings(),
        replay.endpoint().coating().coatings()
    );
    assert_eq!(accepted.endpoint().gas(), replay.endpoint().gas());
    assert_eq!(accepted.endpoint().audit(), replay.endpoint().audit());
    assert_eq!(accepted.audit(), replay.audit());
    assert_eq!(accepted.audit().dormant_advances, 1);
    let endpoint = accepted.endpoint();
    assert!(endpoint.coating().coatings()[&CellId(1)].mass_kg > 0.0);
    assert_eq!(endpoint.coating().coatings()[&CellId(2)].mass_kg, 0.0);
    assert_eq!(endpoint.audit().cells[&CellId(2)].incoming_mass_kg, 0.0);
    assert_eq!(endpoint.audit().thermal.events, 3);
    close(endpoint.audit().thermal.incoming_mass_kg, 0.02, 1e-14);
    close(endpoint.audit().thermal.external_heat_j, 0.04, 1e-12);
    check(endpoint);
    assert!(original.cells().is_empty());
}

#[test]
fn coarse_work_uses_the_original_budget_and_late_cancellation_preserves_all_state() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let mut bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    for binding in &mut bindings {
        binding.accommodation_coefficient = 0.0;
    }
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate { datum: 0.0 };
    let law = fixture::phase();
    let events = vec![
        birth(0.0, cell(1, substrate, 1.0)),
        birth(0.0, cell(2, substrate, 0.0)),
    ];
    let enough = initial(
        &substrate,
        &law,
        &deposit_h,
        &bindings,
        events.clone(),
        [0.0; 2],
        3,
    );
    let short = initial(&substrate, &law, &deposit_h, &bindings, events, [0.0; 2], 2);
    let selection = selected(&bindings);
    let mut polls = 0;
    let success = advance_coatings_adaptive_with_gas(
        &enough,
        &selection,
        None,
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        accuracy(0.001),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(success.audit().accepted_intervals, 1);
    assert_eq!(success.audit().thermal_substeps_evaluated, 3);
    assert_eq!(success.endpoint().audit().thermal.substeps, 2);
    close(
        success.endpoint().audit().thermal.external_heat_j,
        0.001,
        1e-14,
    );
    assert!(matches!(
        advance_coatings_adaptive_with_gas(
            &short,
            &selection,
            None,
            &gas,
            &catalog,
            0.001,
            heat(0.001),
            controls(),
            accuracy(0.001),
            || Ok(())
        ),
        Err(GasCoatingAdaptiveError::Timeline(
            GasCoatingTimelineError::Thermal(ThermalSublimationError::Thermal(
                DepositionError::Budget(_)
            ))
        ))
    ));
    for stop in [1, polls] {
        let mut current = 0;
        assert!(matches!(
            advance_coatings_adaptive_with_gas(
                &enough,
                &selection,
                None,
                &gas,
                &catalog,
                0.001,
                heat(0.001),
                controls(),
                accuracy(0.001),
                || {
                    current += 1;
                    if current == stop {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                }
            ),
            Err(GasCoatingAdaptiveError::Timeline(
                GasCoatingTimelineError::Thermal(ThermalSublimationError::Thermal(
                    DepositionError::Cancelled
                ))
            ))
        ));
    }
    assert!(enough.cells().is_empty() && enough.coatings().is_empty());
    assert_eq!(enough.pending_events().len(), 2);
}

#[test]
fn deposit_mass_error_is_not_hidden_by_a_heavy_substrate_or_large_datum() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let mut results = Vec::new();
    // The substrate is 100,000 times heavier than the initial coating. Keep
    // datum * mass representable within the original energy audit tolerance.
    for datum in [0.0, 1e5] {
        let substrate = Substrate { datum };
        let events = (1..=2)
            .map(|id| {
                let mut cell = cell(id, substrate, 0.0);
                cell.mass_kg = 10.0;
                cell.geometry.volume_m3 = cell.mass_kg / 2700.0;
                birth(0.0, cell)
            })
            .collect();
        let original = initial(
            &substrate, &law, &deposit_h, &bindings, events, [1e-4; 2], 20000,
        );
        let mut tolerance = accuracy(1.0);
        tolerance.absolute_specific_enthalpy_tolerance_j_kg = 1000.0;
        tolerance.absolute_mass_tolerance_kg = 1e-10;
        let accepted = advance_coatings_adaptive_with_gas(
            &original,
            &selected(&bindings),
            None,
            &gas,
            &catalog,
            0.1,
            heat(0.1),
            controls(),
            tolerance,
            || Ok(()),
        )
        .unwrap();
        assert!(accepted.audit().accepted_intervals > 1);
        assert!(accepted.audit().rejected_trials > 0);
        results.push(accepted.endpoint().coating().coatings()[&CellId(1)].mass_kg);
        close(accepted.endpoint().audit().mass_residual_kg, 0.0, 1e-14);
        close(accepted.endpoint().audit().energy_residual_j, 0.0, 1e-7);
    }
    close(results[0], results[1], 1e-12);
}
