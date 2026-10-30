//! Independent circuit and caloric balances for actual field-dependent sources.
use super::*;
use fs_couple::circuit_enthalpy::conduction::circuit::{
    ConductionCircuitTerminal, response::ConductionResponseCircuit,
};
use fs_phs::circuit::{Branch, CircuitGraph};

const C: f64 = 0.5;
const L: f64 = 0.05;

fn circuit(
    cards: &[MaterialCard; 2],
    coordinate: ConductorMagnitude,
    charge: f64,
    flux: Option<f64>,
) -> ConductionResponseCircuit {
    let graph = CircuitGraph {
        node_count: if flux.is_some() { 3 } else { 2 },
        branches: if flux.is_some() {
            vec![
                (1, 2, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::Inductor { henries: L }),
                (1, 0, Branch::Capacitor { farads: C }),
            ]
        } else {
            vec![
                (1, 0, Branch::Resistor { ohms: 999.0 }),
                (1, 0, Branch::Capacitor { farads: C }),
            ]
        },
        transformers: vec![],
    };
    ConductionResponseCircuit::new(
        problem(cards, coordinate),
        graph,
        ConductionCircuitTerminal {
            branch: 0,
            positive: vec![2, 5],
            negative: vec![0, 3],
        },
        0.0,
        flux.into_iter().collect(),
        vec![charge],
        &[],
        limits(),
        field_controls(),
        controls(),
        || false,
    )
    .unwrap()
}

// Independent elimination for two one-metre layers, each of unit area/mass.
// For a series RLC, V_res(I,T) + (2L/dt+dt/(2C))*I = q0/C + 2*lambda0/dt.
// Thermal backward Euler gives the independent two-cell solve already derived
// in endpoint_at_current. No field solver or constitutive response API is used.
fn rlc_reference(
    dt: f64,
    q0: f64,
    flux0: f64,
    old: [f64; 2],
    r: [f64; 2],
    beta: [f64; 2],
) -> (f64, [f64; 2], [f64; 2]) {
    let rhs = q0 / C + 2.0 * flux0 / dt;
    let mut lo = 0.0;
    let mut hi = rhs / (r[0] + r[1]);
    for _ in 0..100 {
        let j = 0.5 * (lo + hi);
        let (_, e) = endpoint_at_current(j, dt, old, r, beta, true);
        if e[0] + e[1] + (2.0 * L / dt + dt / (2.0 * C)) * j < rhs {
            lo = j;
        } else {
            hi = j;
        }
    }
    let j = 0.5 * (lo + hi);
    let (theta, e) = endpoint_at_current(j, dt, old, r, beta, true);
    (j, theta, e)
}

#[test]
fn material_data_changes_real_rlc_current_local_heat_and_continuing_charge_and_flux() {
    let coordinate = ConductorMagnitude::ElectricField;
    let beta = [0.05, 0.08];
    let mut currents = Vec::new();
    for first in [1.0, 2.0] {
        let r = [first, 2.0];
        let (a, ca) = material(r[0], CP[0], ALPHA[0], beta[0], coordinate, false);
        let (b, cb) = material(r[1], CP[1], ALPHA[1], beta[1], coordinate, false);
        let (cards, curves) = ([a, b], [ca, cb]);
        let context = point(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let layout = network(&bindings);
        for sign in [1.0, -1.0] {
            let mut thermal = initial(&curves, 999.0, None);
            let mut actual = circuit(&cards, coordinate, sign * 3.0, Some(0.0));
            let (mut q, mut flux, mut old) = (3.0, 0.0, [0.0; 2]);
            let total_initial = 3.0_f64.powi(2) / (2.0 * C);
            for step in 1..=3 {
                let (j, theta, e) = rlc_reference(0.1, q, flux, old, r, beta);
                let next = layout
                    .step_response_circuit_candidate(
                        &thermal,
                        step as f64 * 0.1,
                        &actual,
                        &[],
                        &[],
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                let electrical = &next.report().electrical;
                near(electrical.circuit.spatial_current_a, sign * j);
                near(electrical.circuit.spatial_voltage_v, sign * (e[0] + e[1]));
                assert!(electrical.field_solves > 1);
                assert!(next.report().circuit_trials >= electrical.field_solves);
                near(next.report().deposited_joule_j, electrical.circuit.heat_j);
                near(electrical.circuit.energy_residual_j, 0.0);
                near(electrical.circuit.native.supplied, 0.0);
                q -= 0.1 * j;
                flux = 2.0 * L * j - flux;
                near(next.circuit().charges_c()[0], sign * q);
                near(next.circuit().fluxes_wb()[0], sign * flux);
                near(next.circuit().time_s(), step as f64 * 0.1);
                for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                    near(
                        next.thermal().temperature_k(cell).unwrap(),
                        300.0 + theta[i],
                    );
                    for t in 2 * i..2 * i + 2 {
                        near(
                            electrical.circuit.spatial.current_density_a_m2[t][0],
                            -sign * j,
                        );
                        near(electrical.circuit.spatial.cell_power_w[t], 0.5 * e[i] * j);
                        let source = &electrical.evaluated[t].conductor;
                        near(source.material().query().axes()["T"], 300.0 + theta[i]);
                        near(source.material().query().axes()["operating"], e[i]);
                        assert_eq!(source.source_key().name(), "original_transport");
                    }
                }
                near(
                    next.circuit().charges_c()[0].powi(2) / (2.0 * C)
                        + next.circuit().fluxes_wb()[0].powi(2) / (2.0 * L)
                        + CP[0] * (next.thermal().temperature_k(CellId(30)).unwrap() - 300.0)
                        + CP[1] * (next.thermal().temperature_k(CellId(10)).unwrap() - 300.0),
                    total_initial,
                );
                if step == 1 && sign == 1.0 {
                    currents.push(electrical.circuit.spatial_current_a);
                }
                old = theta;
                (thermal, actual, _) = next.into_parts();
            }
        }
    }
    assert!(currents[0] > currents[1]);
}

#[test]
fn fixed_context_rc_uses_actual_midpoint_current_not_initial_or_unit_field() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let (a, _) = material(1.0, CP[0], 0.0, 0.1, coordinate, false);
    let (b, _) = material(2.0, CP[1], 0.0, 0.2, coordinate, false);
    let cards = [a, b];
    let mut actual = circuit(&cards, coordinate, 3.0, None);
    let mut q = 3.0;
    for step in 1..=3 {
        // (3 + dt/(2C))*I + 0.3 I² = q0/C.
        let a = 3.0 + 0.1 / (2.0 * C);
        let j = 2.0 * q / C / (a + (a * a + 1.2 * q / C).sqrt());
        let next = actual
            .step_candidate(
                step as f64 * 0.1,
                &[],
                limits(),
                field_controls(),
                controls(),
                config().circuit,
                || false,
            )
            .unwrap();
        let report = next.step();
        near(report.circuit.spatial_current_a, j);
        near(report.circuit.spatial_voltage_v, (3.0 + 0.3 * j) * j);
        near(report.circuit.heat_j, 0.1 * (3.0 + 0.3 * j) * j * j);
        for source in &report.evaluated {
            near(source.conductor.material().query().axes()["operating"], j);
        }
        q -= 0.1 * j;
        near(next.circuit().charges_c()[0], q);
        assert!(report.field_solves > 1);
        (actual, _) = next.into_parts();
    }
}

#[test]
fn response_trials_share_the_heat_budget_and_refusals_leave_all_histories_unchanged() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], 0.05, coordinate, false);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], 0.08, coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let thermal = initial(&curves, 999.0, None);
    let actual = circuit(&cards, coordinate, 3.0, Some(0.0));
    let context = point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let mut budget = config();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_response_circuit_candidate(
            &thermal,
            0.1,
            &actual,
            &[],
            &[],
            limits(),
            field_controls(),
            controls(),
            budget,
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(matches!(
        layout.step_response_circuit_candidate(
            &thermal,
            0.1,
            &actual,
            &[],
            &[],
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || true
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        actual
            .step_candidate(
                0.1,
                &[],
                limits(),
                field_controls(),
                ConductionResponseControls {
                    maximum_iterations: 1,
                    ..controls()
                },
                config().circuit,
                || false
            )
            .is_err()
    );
    assert_eq!(actual.time_s().to_bits(), 0.0_f64.to_bits());
    assert_eq!(actual.charges_c(), &[3.0]);
    assert_eq!(actual.fluxes_wb(), &[0.0]);
    assert_eq!(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
    let run = || {
        layout
            .step_response_circuit_candidate(
                &thermal,
                0.1,
                &actual,
                &[],
                &[],
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap()
    };
    let a = run();
    let b = run();
    assert_eq!(a.circuit().charges_c(), b.circuit().charges_c());
    assert_eq!(a.circuit().fluxes_wb(), b.circuit().fluxes_wb());
    assert_eq!(
        a.report().electrical.circuit.spatial.cell_power_w,
        b.report().electrical.circuit.spatial.cell_power_w
    );
}
