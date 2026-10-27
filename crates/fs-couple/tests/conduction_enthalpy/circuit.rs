//! Spatial field, native circuit storage and finite heat share one endpoint.
use super::*;
use fs_couple::circuit_enthalpy::{
    ResistorCell,
    conduction::circuit::{ConductionCircuit, ConductionCircuitTerminal},
};
use fs_phs::circuit::{Branch, CircuitGraph};

fn terminal() -> ConductionCircuitTerminal {
    ConductionCircuitTerminal {
        branch: 0,
        positive: vec![2, 5],
        negative: vec![0, 3],
    }
}
fn rc(capacitance: f64) -> CircuitGraph {
    CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Resistor { ohms: 999.0 }),
            (
                1,
                0,
                Branch::Capacitor {
                    farads: capacitance,
                },
            ),
        ],
        transformers: vec![],
    }
}

#[test]
fn g1_capacitor_discharges_through_actual_layered_field_and_local_finite_heat() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material_with_electrical_basis(
        ElectricalPropertyBasis::Conductivity,
        1.0 / 3.0,
        2.0,
        0.0,
        1.0,
    );
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 999.0, None);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let circuit = ConductionCircuit::new(
        problem(
            &cards,
            [
                ElectricalPropertyBasis::Resistivity,
                ElectricalPropertyBasis::Conductivity,
            ],
            300.0,
            1.0,
        ),
        rc(0.5),
        terminal(),
        0.0,
        vec![],
        vec![4.0],
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    let next = layout
        .step_circuit_candidate(
            &thermal,
            0.5,
            &circuit,
            &[],
            &[],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    // R=1+3=4 ohm, C=.5 F. Native midpoint gives
    // q1/q0=(2RC-dt)/(2RC+dt)=7/9 and I=(q0-q1)/dt=16/9 A.
    let q1 = 28.0 / 9.0;
    let current = 16.0 / 9.0;
    let r = next.report();
    let e = &r.electrical;
    near(next.circuit().charges_c()[0], q1);
    near(e.spatial_current_a, current);
    near(e.spatial_voltage_v, 4.0 * current);
    near(e.stored_before_j, 16.0);
    near(e.stored_after_j, q1 * q1);
    near(e.native.supplied, 0.0);
    near(e.heat_j, 16.0 - q1 * q1);
    near(r.non_electrical_heat_j, 0.0);
    near(r.deposited_joule_j, e.heat_j);
    for (i, &p) in e.spatial.cell_power_w.iter().enumerate() {
        near(p, 0.5 * current * current * if i < 2 { 1.0 } else { 3.0 });
        near(e.spatial.current_density_a_m2[i][0], -current);
    }
    // Independent two-store BE conduction solve with heat [I²/2, 3I²/2].
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        300.0 + 2.0 * current * current / 3.5,
    );
    near(
        next.thermal().temperature_k(CellId(10)).unwrap(),
        300.0 + 2.5 * current * current / 3.5,
    );
    near(
        next.thermal().energy_j() + e.stored_after_j,
        thermal.energy_j() + 16.0,
    );
    assert!(r.total_energy_residual_j.abs() < config().total_energy_tolerance_j);
    assert_eq!(next.circuit().time_s(), next.thermal().time_s());
    assert!(e.graph_resistors.is_empty());
    assert_eq!(circuit.charges_c(), &[4.0]);
    assert_eq!(thermal.time_s(), 0.0);
    let continued = layout
        .step_circuit_candidate(
            next.thermal(),
            1.0,
            next.circuit(),
            &[],
            &[],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    near(continued.circuit().charges_c()[0], q1 * 7.0 / 9.0);
    near(
        continued.thermal().energy_j() + continued.report().electrical.stored_after_j,
        thermal.energy_j() + 16.0,
    );
}

#[test]
fn g1_inductor_flux_releases_its_actual_energy_into_the_spatial_conductor() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 0.0, None);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Resistor { ohms: 123.0 }),
            (0, 1, Branch::Inductor { henries: 2.0 }),
        ],
        transformers: vec![],
    };
    let circuit = ConductionCircuit::new(
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        ),
        graph,
        terminal(),
        0.0,
        vec![4.0],
        vec![],
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    let next = network(&bindings)
        .step_circuit_candidate(
            &thermal,
            0.25,
            &circuit,
            &[],
            &[],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    // phi1=phi0(2L-Rdt)/(2L+Rdt), with L=2 H, R=4 ohm and dt=.25 s.
    near(next.circuit().fluxes_wb()[0], 2.4);
    near(next.report().electrical.spatial_current_a, 1.6);
    near(next.report().electrical.stored_before_j, 4.0);
    near(next.report().electrical.stored_after_j, 1.44);
    near(next.report().electrical.heat_j, 2.56);
    near(next.thermal().energy_j() - thermal.energy_j(), 2.56);
    near(next.report().electrical.native.supplied, 0.0);
}

// Independent scalar balance: two identical unit-capacity slabs and C=1 F.
// Heat conservation gives T_end-T_start=(q_start²-q_end²)/4, while midpoint
// RC stepping uses R_end=2*rho0*(1+alpha*(T_end-300)). No production solver,
// material interpolation, field extraction or thermal residual computes this.
fn rc_reference(q0: f64, t0: f64, rho0: f64, alpha: f64, dt: f64) -> (f64, f64) {
    let mut low = 0.0;
    let mut high = q0.abs();
    for _ in 0..90 {
        let q = 0.5 * (low + high);
        let t = t0 + (q0 * q0 - q * q) / 4.0;
        let r = 2.0 * rho0 * (1.0 + alpha * (t - 300.0));
        let rhs = q0.abs() * (2.0 * r - dt) / (2.0 * r + dt);
        if q < rhs {
            low = q;
        } else {
            high = q;
        }
    }
    let q = 0.5 * (low + high);
    (q.copysign(q0), t0 + (q0 * q0 - q * q) / 4.0)
}

#[test]
fn g3_temperature_and_material_data_change_charge_decay_and_continued_fields() {
    let alpha = 0.01;
    let mut positive = Vec::new();
    for rho in [1.0, 2.0] {
        for sign in [1.0, -1.0] {
            let (a, ca) = material(rho, 1.0, alpha);
            let (b, cb) = material(rho, 1.0, alpha);
            let cards = [a, b];
            let curves = [ca, cb];
            let thermal = initial(&curves, 0.0, None);
            let context = point(300.0);
            let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let bindings = cells(&cards, &curves, &context, &key);
            let layout = network(&bindings);
            let circuit = ConductionCircuit::new(
                problem(
                    &cards,
                    [ElectricalPropertyBasis::Resistivity; 2],
                    300.0,
                    1.0,
                ),
                rc(1.0),
                terminal(),
                0.0,
                vec![],
                vec![sign * 4.0],
                &[],
                field_controls(),
                || false,
            )
            .unwrap();
            let first = layout
                .step_circuit_candidate(
                    &thermal,
                    0.5,
                    &circuit,
                    &[],
                    &[],
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let (q1, t1) = rc_reference(sign * 4.0, 300.0, rho, alpha, 0.5);
            near(first.circuit().charges_c()[0], q1);
            near(first.thermal().temperature_k(CellId(10)).unwrap(), t1);
            near(first.thermal().temperature_k(CellId(30)).unwrap(), t1);
            near(
                first.circuit().problem().conductors()[0].resistivity_ohm_m(),
                rho * (1.0 + alpha * (t1 - 300.0)),
            );
            assert_eq!(first.report().electrical.spatial_current_a.signum(), sign);
            let second = layout
                .step_circuit_candidate(
                    first.thermal(),
                    1.0,
                    first.circuit(),
                    &[],
                    &[],
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let (q2, t2) = rc_reference(q1, t1, rho, alpha, 0.5);
            near(second.circuit().charges_c()[0], q2);
            near(second.thermal().temperature_k(CellId(30)).unwrap(), t2);
            near(
                second.thermal().energy_j() + second.report().electrical.stored_after_j,
                thermal.energy_j() + 8.0,
            );
            if sign > 0.0 {
                positive.push((q2, t2));
            }
        }
    }
    assert!(positive[1].0 > positive[0].0);
    assert!(positive[1].1 < positive[0].1);
}

#[test]
fn g1_external_circuit_work_and_all_graph_resistor_heat_close_once() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 999.0, None);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::Resistor { ohms: 999.0 }),
            (1, 2, Branch::Resistor { ohms: 2.0 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let circuit = ConductionCircuit::new(
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        ),
        graph,
        terminal(),
        0.0,
        vec![],
        vec![],
        &[12.0],
        field_controls(),
        || false,
    )
    .unwrap();
    let route = ResistorCell {
        branch: 1,
        cell: CellId(30),
    };
    let layout = network(&bindings);
    let next = layout
        .step_circuit_candidate(
            &thermal,
            0.5,
            &circuit,
            &[12.0],
            &[route],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let e = &next.report().electrical;
    near(e.spatial_current_a, 2.0);
    near(e.spatial.heat_j, 8.0);
    near(e.graph_resistors[0].heat_j, 4.0);
    near(e.heat_j, 12.0);
    near(e.native.supplied, 12.0);
    near(next.thermal().energy_j() - thermal.energy_j(), 12.0);
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        300.0 + 18.0 / 3.5,
    );
    near(
        next.thermal().temperature_k(CellId(10)).unwrap(),
        300.0 + 12.0 / 3.5,
    );
    for bad in [
        vec![],
        vec![route, route],
        vec![ResistorCell { branch: 0, ..route }],
        vec![ResistorCell {
            cell: CellId(999),
            ..route
        }],
    ] {
        assert!(
            layout
                .step_circuit_candidate(
                    &thermal,
                    0.5,
                    &circuit,
                    &[12.0],
                    &bad,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.time_s(), 0.0);
}

#[test]
fn g4_late_refusal_preserves_charge_heat_sources_and_retry() {
    let (a, ca) = material(1.0, 1.0, 0.01);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 0.0, None);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let circuit = ConductionCircuit::new(
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        ),
        rc(0.5),
        terminal(),
        0.0,
        vec![],
        vec![4.0],
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    let saved_cells = thermal.cells().clone();
    let saved_sources = circuit.problem().conductors().to_vec();
    let saved_graph = circuit.graph().clone();
    let mut count = 0;
    let accepted = layout
        .step_circuit_candidate(
            &thermal,
            0.5,
            &circuit,
            &[],
            &[],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        layout.step_circuit_candidate(
            &thermal,
            0.5,
            &circuit,
            &[],
            &[],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                polls == count
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        layout
            .step_circuit_candidate(
                &thermal,
                0.5,
                &circuit,
                &[1.0],
                &[],
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    let mut exhausted = config();
    exhausted.maximum_circuit_trials = 1;
    assert!(
        layout
            .step_circuit_candidate(
                &thermal,
                0.5,
                &circuit,
                &[],
                &[],
                limits(),
                field_controls(),
                exhausted,
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(circuit.problem().conductors(), &saved_sources);
    assert_eq!(circuit.graph(), &saved_graph);
    assert_eq!(circuit.charges_c(), &[4.0]);
    let retry = layout
        .step_circuit_candidate(
            &thermal,
            0.5,
            &circuit,
            &[],
            &[],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(retry.thermal().cells(), accepted.thermal().cells());
    assert_eq!(retry.circuit().charges_c(), accepted.circuit().charges_c());
}
