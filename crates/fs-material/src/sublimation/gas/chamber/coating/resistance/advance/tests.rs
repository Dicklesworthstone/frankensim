//! G1 independent heat/ODE limits; G3 datum/order; G4 whole-call refusal.
use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, chamber::GasChamberSolve, tests as fixture};
use fs_qty::chemistry::SpeciesId;
use fs_thermal::sublimation::{
    SolidSurfaceState, SolverConfig, SublimationLaw, coating::resistance::ResistiveSurfaceState,
};

const SUPPORT_MASS: f64 = 0.02;
const SUPPORT_CP: f64 = 500.0;
const FILM_MASS: f64 = 0.002;
const AREA: f64 = 1e-4;

fn controls() -> GasResistiveAdvanceSolve {
    GasResistiveAdvanceSolve {
        step: GasSupportedChamberSolve {
            chamber: GasChamberSolve {
                surface: fixture::controls(),
                maximum_surfaces: 8,
            },
            storage: SolverConfig {
                relative_tolerance: 1e-12,
                absolute_energy_tolerance_j: 1e-10,
                max_iterations: 128,
            },
        },
        absolute_mass_tolerance_kg: 1e-10,
        absolute_temperature_tolerance_k: 1e-3,
        absolute_specific_enthalpy_tolerance_j_kg: 0.25,
        absolute_heat_tolerance_j: 1e-4,
        relative_tolerance: 1e-6,
        minimum_dt_s: 1e-10,
        maximum_dt_s: 0.5,
        maximum_trials: 20000,
    }
}
fn state(ts: f64, td: f64, chart: (f64, f64), support_datum: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: SUPPORT_MASS,
        substrate_specific_enthalpy_j_kg: support_datum + SUPPORT_CP * (ts - 400.0),
        deposit: SolidSurfaceState {
            mass_kg: FILM_MASS,
            area_m2: AREA,
            specific_enthalpy_j_kg: chart.0 + chart.1 * (td - 400.0),
        },
    }
}
fn support(h: f64) -> Result<f64, SublimationError> {
    Ok(400.0 + h / SUPPORT_CP)
}
fn row<'a>(
    id: u64,
    state: ResistiveSurfaceState,
    conductance: f64,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    substrate: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    deposit: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasResistiveSurface<'a> {
    GasResistiveSurface {
        id,
        state,
        conductance_w_k: conductance,
        law,
        binding,
        substrate_temperature: substrate,
        deposit_temperature: deposit,
    }
}
fn loose() -> GasResistiveAdvanceSolve {
    GasResistiveAdvanceSolve {
        absolute_mass_tolerance_kg: 1e-7,
        absolute_temperature_tolerance_k: 0.05,
        absolute_specific_enthalpy_tolerance_j_kg: 100.0,
        absolute_heat_tolerance_j: 0.2,
        relative_tolerance: 1e-4,
        ..controls()
    }
}
fn assert_balanced(candidate: &GasResistiveAdvanceCandidate<'_>) {
    let audit = candidate.audit();
    for residual in [
        audit.chamber.energy_residual_j,
        audit.chamber.gas_stream_energy_residual_j,
        audit.chamber.solid_stream_energy_residual_j,
    ] {
        assert!(residual.abs() <= controls().step.chamber.surface.gas.energy_tolerance_j);
    }
    assert!(
        audit.chamber.mass_residual_kg.abs()
            <= controls().step.chamber.surface.gas.mass_tolerance_kg
    );
    assert!(
        audit.chamber.maximum_species_mass_residual_kg
            <= controls().step.chamber.surface.gas.mass_tolerance_kg
    );
    for row in &audit.thermal {
        assert!(
            row.substrate_energy_residual_j.abs()
                <= controls().step.chamber.surface.gas.energy_tolerance_j
        );
        assert!(
            row.deposit_energy_residual_j.abs()
                <= controls().step.chamber.surface.gas.energy_tolerance_j
        );
    }
}

#[test]
fn tighter_contact_heat_control_approaches_the_independent_two_body_exponential() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let deposit = |h| fixture::temperature(h, chart);
    let rows = [row(
        7,
        state(300.0, 350.0, chart, 0.0),
        2.0,
        &law,
        &binding,
        &support,
        &deposit,
    )];
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let films = BTreeMap::from([(7, 0.0)]);
    let duration = 0.5;
    let cs = SUPPORT_MASS * SUPPORT_CP;
    let cd = FILM_MASS * chart.1;
    let common = (cs * 300.0 + cd * 350.0) / (cs + cd);
    let difference = (300.0 - 350.0) * (-2.0 * (1.0 / cs + 1.0 / cd) * duration).exp();
    let expected = (
        common + cd / (cs + cd) * difference,
        common - cs / (cs + cd) * difference,
    );
    let mut errors = Vec::new();
    let mut intervals = Vec::new();
    for solve in [loose(), controls()] {
        let accepted = advance_resistive_chamber_adaptive_with_films(
            &rows,
            None,
            &films,
            &gas,
            &catalog,
            duration,
            solve,
            || false,
        )
        .unwrap();
        let surface = &accepted.endpoint().surfaces()[0];
        let history = &accepted.audit().thermal[0];
        errors.push(
            (surface.substrate_temperature_k[1] - expected.0).abs()
                + (surface.exchange.surface_temperature_k[1] - expected.1).abs(),
        );
        intervals.push(accepted.audit().chamber.accepted_intervals);
        fixture::near(
            history.interface_heat_j,
            -cs * (surface.substrate_temperature_k[1] - 300.0),
            1e-7,
        );
        fixture::near(
            history.substrate_energy_change_j + history.deposit_energy_change_j,
            0.0,
            1e-7,
        );
        assert_eq!(surface.state.deposit.mass_kg, FILM_MASS);
        assert_eq!(accepted.gas(), &gas);
        assert_eq!(accepted.audit().chamber.elapsed_s, duration);
        assert!(accepted.audit().chamber.maximum_accepted_error_ratio <= 1.0);
        assert!(accepted.audit().chamber.rejected_trials > 0);
        assert_balanced(&accepted);
    }
    assert!(errors[1] < errors[0] * 0.2, "{errors:?}");
    assert!(intervals[1] > intervals[0]);
}

// Independent continuum equations for this manufactured constant-Cp mixture.
// State = [deposit mass, substrate T, deposit T, vapor mass, gas T]. No production
// sublimation/heat endpoint or flux function enters this reference integration.
fn derivative(y: [f64; 5], carrier_capacity: f64, conductance: f64) -> [f64; 5] {
    let [mass, ts, td, vapor, tg] = y;
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let saturation = 1000.0 * (200000.0 / r * (1.0 / 400.0 - 1.0 / td)).exp();
    let partial = vapor * r * tg / 0.01;
    let outgoing = AREA * 0.1 * saturation / (2.0 * core::f64::consts::PI * r * td).sqrt();
    let incoming = AREA * 0.1 * partial / (2.0 * core::f64::consts::PI * r * tg).sqrt();
    let dm = outgoing - incoming;
    let stream = outgoing * 3.5 * r * td - incoming * 3.5 * r * tg;
    let contact = conductance * (ts - td);
    [
        -dm,
        -contact / (SUPPORT_MASS * SUPPORT_CP),
        (-stream + dm * (4.0 * r * td - 200000.0) + contact) / (mass * 4.0 * r),
        dm,
        (stream - dm * 3.0 * r * tg) / (carrier_capacity + vapor * 3.0 * r),
    ]
}
fn rk4(mut y: [f64; 5], duration: f64, carrier_capacity: f64, conductance: f64) -> [f64; 5] {
    let dt = duration / 8192.0;
    let add = |a: [f64; 5], b: [f64; 5], scale| std::array::from_fn(|i| a[i] + scale * b[i]);
    for _ in 0..8192 {
        let a = derivative(y, carrier_capacity, conductance);
        let b = derivative(add(y, a, 0.5 * dt), carrier_capacity, conductance);
        let c = derivative(add(y, b, 0.5 * dt), carrier_capacity, conductance);
        let d = derivative(add(y, c, dt), carrier_capacity, conductance);
        y = std::array::from_fn(|i| y[i] + dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    y
}

#[test]
fn coupled_deposition_and_two_temperatures_converge_to_independent_finite_gas_ode() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let deposit = |h| fixture::temperature(h, chart);
    let rows = [row(
        3,
        state(300.0, 350.0, chart, 0.0),
        2.0,
        &law,
        &binding,
        &support,
        &deposit,
    )];
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let vapor = gas.species_mass_kg(&fixture::species()).unwrap();
    let carrier = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap();
    let carrier_capacity =
        carrier * 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let reference = rk4(
        [FILM_MASS, 300.0, 350.0, vapor, 400.0],
        1.0,
        carrier_capacity,
        2.0,
    );
    let mut errors = Vec::new();
    for solve in [loose(), controls()] {
        let accepted =
            advance_resistive_chamber_adaptive(&rows, None, &gas, &catalog, 1.0, solve, || false)
                .unwrap();
        let end = &accepted.endpoint().surfaces()[0];
        let actual = [
            end.state.deposit.mass_kg,
            end.substrate_temperature_k[1],
            end.exchange.surface_temperature_k[1],
            accepted.gas().species_mass_kg(&fixture::species()).unwrap(),
            accepted.gas().state().state().temperature,
        ];
        let error: f64 = actual
            .iter()
            .zip(reference)
            .zip([1e5, 1.0, 1.0, 1e5, 1.0])
            .map(|((&a, b), scale)| scale * (a - b).abs())
            .sum();
        errors.push(error);
        assert!(actual[0] > FILM_MASS && actual[3] < vapor);
        fixture::near(actual[0] + actual[3], FILM_MASS + vapor, 1e-14);
        fixture::near(
            accepted
                .gas()
                .species_mass_kg(&SpeciesId::new("A").unwrap())
                .unwrap(),
            carrier,
            1e-14,
        );
        fixture::near(
            accepted.audit().chamber.surfaces[0].transferred_mass_kg,
            FILM_MASS - actual[0],
            1e-14,
        );
        assert_balanced(&accepted);
    }
    assert!(errors[1] < errors[0] * 0.2, "{errors:?}");
    assert!(errors[1] < 0.1, "{errors:?}");
}

#[test]
fn accepted_heat_is_exactly_the_fine_history_and_input_order_does_not_change_it() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let deposit = |h| fixture::temperature(h, chart);
    let rows = [
        row(
            42,
            state(300.0, 350.0, chart, 0.0),
            0.5,
            &law,
            &binding,
            &support,
            &deposit,
        ),
        row(
            7,
            state(330.0, 375.0, chart, 0.0),
            3.0,
            &law,
            &binding,
            &support,
            &deposit,
        ),
    ];
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let accommodation = ThermalAccommodation::new(
        &catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.25),
            (fixture::species(), 0.5),
        ],
        [200.0, 900.0],
    )
    .unwrap();
    let reflected = BTreeMap::from([(7, accommodation.clone()), (42, accommodation)]);
    let solve = GasResistiveAdvanceSolve {
        absolute_mass_tolerance_kg: 1.0,
        absolute_temperature_tolerance_k: 1e6,
        absolute_specific_enthalpy_tolerance_j_kg: 1e6,
        absolute_heat_tolerance_j: 1e6,
        relative_tolerance: 0.0,
        ..controls()
    };
    let duration = 0.001;
    let accepted = advance_resistive_chamber_adaptive(
        &rows,
        Some(&reflected),
        &gas,
        &catalog,
        duration,
        solve,
        || false,
    )
    .unwrap();
    let replay = advance_resistive_chamber_adaptive(
        &[rows[1], rows[0]],
        Some(&reflected),
        &gas,
        &catalog,
        duration,
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(accepted, replay);
    assert_eq!(accepted.audit().chamber.accepted_intervals, 1);
    assert_eq!(accepted.audit().chamber.kernel_evaluations, 3);
    let first = super::super::advance_resistive_chamber(
        &rows,
        Some(&reflected),
        &gas,
        &catalog,
        0.5 * duration,
        solve.step,
        || false,
    )
    .unwrap();
    let middle = [
        GasResistiveSurface {
            state: first.surfaces()[0].state,
            ..rows[1]
        },
        GasResistiveSurface {
            state: first.surfaces()[1].state,
            ..rows[0]
        },
    ];
    let second = super::super::advance_resistive_chamber(
        &middle,
        Some(&reflected),
        first.gas(),
        &catalog,
        0.5 * duration,
        solve.step,
        || false,
    )
    .unwrap();
    assert_eq!(accepted.endpoint(), &second);
    let mut expected_reflection = 0.0;
    for (((a, b), mass), heat) in first
        .surfaces()
        .iter()
        .zip(second.surfaces())
        .zip(&accepted.audit().chamber.surfaces)
        .zip(&accepted.audit().thermal)
    {
        fixture::near(
            mass.energy_into_gas_j,
            a.exchange.energy_into_gas_j + b.exchange.energy_into_gas_j,
            1e-12,
        );
        fixture::near(
            heat.interface_heat_j,
            a.interface_heat_j + b.interface_heat_j,
            1e-12,
        );
        fixture::near(
            mass.reflected_heat_into_gas_j,
            a.reflected_heat_into_gas_j + b.reflected_heat_into_gas_j,
            1e-12,
        );
        expected_reflection += mass.reflected_heat_into_gas_j;
    }
    fixture::near(
        accepted.audit().chamber.reflected_heat_into_gas_j,
        expected_reflection,
        1e-12,
    );
    assert!(expected_reflection < 0.0);
    assert_balanced(&accepted);
}

#[test]
fn independent_condensed_and_gas_datums_do_not_relax_the_physical_error_norm() {
    let reference_catalog = fixture::catalog(0.0);
    let chart = fixture::chart(&reference_catalog, &fixture::binding(0.0));
    let law = fixture::phase();
    let deposit = |h| fixture::temperature(h, chart);
    let mut reference: Option<(f64, f64, f64, f64, f64)> = None;
    for (gas_datum, support_datum) in [(0.0, 0.0), (1e13, 1e6)] {
        let catalog = fixture::catalog(gas_datum);
        let binding = fixture::binding(gas_datum);
        let substrate = |h| support(h - support_datum);
        let rows = [row(
            9,
            state(300.0, 350.0, chart, support_datum),
            2.0,
            &law,
            &binding,
            &substrate,
            &deposit,
        )];
        let gas = fixture::gas(&catalog, 400.0, 3000.0);
        let mut solve = controls();
        solve
            .step
            .chamber
            .surface
            .latent_compatibility_tolerance_j_kg = 0.01;
        let accepted =
            advance_resistive_chamber_adaptive(&rows, None, &gas, &catalog, 0.1, solve, || false)
                .unwrap();
        let end = &accepted.endpoint().surfaces()[0];
        let physical = (
            end.state.deposit.mass_kg,
            end.substrate_temperature_k[1],
            end.exchange.surface_temperature_k[1],
            accepted.gas().state().state().temperature,
            accepted.audit().thermal[0].interface_heat_j,
        );
        if let Some(expected) = reference {
            fixture::near(physical.0, expected.0, 1e-11);
            fixture::near(physical.1, expected.1, 2e-6);
            fixture::near(physical.2, expected.2, 2e-6);
            fixture::near(physical.3, expected.3, 2e-6);
            fixture::near(physical.4, expected.4, 1e-5);
            assert!(accepted.audit().chamber.energy_into_gas_j.abs() > 1e6);
        } else {
            reference = Some(physical);
        }
        assert!(accepted.audit().chamber.accepted_intervals > 1);
        assert_balanced(&accepted);
    }
}

#[test]
fn inventory_retry_zero_time_budgets_and_late_cancellation_keep_the_original_inputs() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let deposit = |h| fixture::temperature(h, chart);
    let initial = state(350.0, 400.0, chart, 0.0);
    let rows = [row(5, initial, 2.0, &law, &binding, &support, &deposit)];
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let original_gas = gas.clone();
    let zero =
        advance_resistive_chamber_adaptive(&rows, None, &gas, &catalog, 0.0, controls(), || false)
            .unwrap();
    assert_eq!(zero.endpoint().surfaces()[0].state, initial);
    assert_eq!(zero.gas(), &gas);
    assert_eq!(zero.audit().chamber.trials, 0);
    assert_eq!(zero.audit().chamber.kernel_evaluations, 1);
    assert_eq!(zero.audit().thermal[0].interface_heat_j, 0.0);
    let mut solve = controls();
    solve.step.chamber.surface.maximum_gross_inventory_fraction = 1e-5;
    let mut polls = 0;
    let result =
        advance_resistive_chamber_adaptive(&rows, None, &gas, &catalog, 0.001, solve, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(result.audit().chamber.rejected_trials > 0);
    assert!(
        result.audit().chamber.kernel_evaluations > 2 * result.audit().chamber.accepted_intervals
    );
    assert_balanced(&result);
    for stop in [1, polls] {
        let mut seen = 0;
        assert!(matches!(
            advance_resistive_chamber_adaptive(&rows, None, &gas, &catalog, 0.001, solve, || {
                seen += 1;
                seen == stop
            }),
            Err(GasResistiveAdvanceError::Surface(GasSurfaceError::Surface(
                SublimationError::Cancelled
            )))
        ));
    }
    solve.maximum_trials = 1;
    assert!(matches!(
        advance_resistive_chamber_adaptive(&rows, None, &gas, &catalog, 0.001, solve, || false),
        Err(GasResistiveAdvanceError::TrialBudget { maximum_trials: 1 })
    ));
    solve.minimum_dt_s = 0.0005;
    assert!(matches!(
        advance_resistive_chamber_adaptive(&rows, None, &gas, &catalog, 0.001, solve, || false),
        Err(GasResistiveAdvanceError::MinimumStep { .. })
    ));
    let bare = [GasResistiveSurface {
        state: ResistiveSurfaceState {
            deposit: SolidSurfaceState {
                mass_kg: 0.0,
                ..initial.deposit
            },
            ..initial
        },
        ..rows[0]
    }];
    assert!(matches!(
        advance_resistive_chamber_adaptive(&bare, None, &gas, &catalog, 0.001, controls(), || {
            false
        }),
        Err(GasResistiveAdvanceError::Invalid(_))
    ));
    assert_eq!(rows[0].state, initial);
    assert_eq!(gas, original_gas);
}
