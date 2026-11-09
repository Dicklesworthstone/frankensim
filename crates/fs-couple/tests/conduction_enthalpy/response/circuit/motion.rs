//! G1/G3/G4: native RLC storage, nonlinear material-frame current and motion.
use super::*;
use fs_em::conduction::motion::ConductionMotion;

fn motion(speed: f64) -> ([[f64; 2]; 4], [f64; 4]) {
    ([[0.0, speed]; 4], [1.0; 4])
}

fn sources(coordinate: ConductorMagnitude) -> [MaterialCard; 2] {
    [
        material(1.0, CP[0], 0.0, 0.05, coordinate, false).0,
        material(2.0, CP[1], 0.0, 0.08, coordinate, false).0,
    ]
}

fn terminal() -> ConductionCircuitTerminal {
    ConductionCircuitTerminal {
        branch: 0,
        positive: vec![2, 5],
        negative: vec![0, 3],
    }
}

// Independent series reduction. This is not a call to either field/response owner.
fn fields(j: f64, coordinate: ConductorMagnitude) -> [f64; 2] {
    [0, 1].map(|i| {
        let r = [1.0, 2.0][i];
        let beta = [0.05, 0.08][i];
        match coordinate {
            ConductorMagnitude::ElectricField => r * j / (1.0 - beta * j.abs()),
            ConductorMagnitude::CurrentDensity => (r + beta * j.abs()) * j,
        }
    })
}

fn current(rhs: f64, extra: f64, coordinate: ConductorMagnitude) -> f64 {
    let (mut lo, mut hi) = (-5.0, 5.0);
    for _ in 0..120 {
        let j = 0.5 * (lo + hi);
        let e = fields(j, coordinate);
        if e[0] + e[1] + extra * j < rhs {
            lo = j;
        } else {
            hi = j;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn moving_rc_retains_actual_nonlinear_sources_and_signed_work_beside_dormant_storage() {
    let dt = 0.1;
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let cards = sources(coordinate);
        for sign in [-1.0, 1.0] {
            for speed in [0.0, 0.5, 3.0, 5.0] {
                let speed = sign * speed;
                let (velocity, flux_density) = motion(speed);
                let drive = ConductionMotion {
                    velocity_m_s: &velocity,
                    normal_flux_density_t: &flux_density,
                };
                let mut circuit = ConductionResponseCircuit::new_with_motion(
                    problem(&cards, coordinate).with_relaxation(0.5).unwrap(),
                    CircuitGraph {
                        node_count: 3,
                        branches: vec![
                            (1, 0, Branch::Resistor { ohms: 999.0 }),
                            (1, 0, Branch::Capacitor { farads: C }),
                            (2, 0, Branch::Capacitor { farads: 1.0 }),
                        ],
                        transformers: vec![],
                    },
                    terminal(),
                    0.0,
                    vec![],
                    vec![sign * 3.0, 1e12],
                    &[],
                    drive,
                    limits(),
                    field_controls(),
                    controls(),
                    300,
                    || false,
                )
                .unwrap();
                let mut q = sign * 3.0;
                for interval in 1..=2 {
                    let j = current(q / C - 2.0 * speed, dt / (2.0 * C), coordinate);
                    let ep = fields(j, coordinate);
                    let next_q = q - dt * j;
                    let voltage = ep[0] + ep[1] + 2.0 * speed;
                    let heat = dt * j * (ep[0] + ep[1]);
                    let mechanical = dt * 2.0 * speed * j;
                    let next = circuit
                        .step_with_motion_candidate(
                            interval as f64 * dt,
                            &[],
                            drive,
                            limits(),
                            field_controls(),
                            controls(),
                            config().circuit,
                            300,
                            || false,
                        )
                        .unwrap();
                    let report = next.step();
                    let native = &report.circuit;
                    let moving = report.motion.as_ref().unwrap();
                    near(next.circuit().charges_c()[0], next_q);
                    assert_eq!(next.circuit().charges_c()[1], 1e12);
                    assert_eq!(native.stored_after_j - native.stored_before_j, 0.0);
                    near(
                        native.stored_change_j,
                        (next_q - q) * (next_q + q) / (2.0 * C),
                    );
                    near(native.spatial_voltage_v, voltage);
                    near(native.spatial_current_a, j);
                    near(native.supplied_work_j, 0.0);
                    near(native.heat_j, heat);
                    near(report.mechanical_work_j, mechanical);
                    near(native.spatial.supplied_work_j, dt * voltage * j);
                    near(native.stored_change_j + heat + mechanical, 0.0);
                    near(native.energy_residual_j, 0.0);
                    assert_eq!(report.field_solves % 3, 0);
                    for triangle in 0..4 {
                        let layer = triangle / 2;
                        near(native.spatial.current_density_a_m2[triangle][0], -j);
                        near(moving.conductor_field_v_m[triangle][0], -ep[layer]);
                        near(native.spatial.cell_power_w[triangle], 0.5 * j * ep[layer]);
                        let expected = match coordinate {
                            ConductorMagnitude::ElectricField => ep[layer].abs(),
                            ConductorMagnitude::CurrentDensity => j.abs(),
                        };
                        near(
                            report.evaluated[triangle]
                                .conductor
                                .material()
                                .query()
                                .axes()["operating"],
                            expected,
                        );
                    }
                    if speed.abs() == 5.0 {
                        assert!(native.stored_change_j > 0.0);
                        assert!(report.mechanical_work_j < 0.0);
                    }
                    if speed.abs() == 3.0 {
                        near(native.heat_j, 0.0);
                        near(native.spatial_current_a, 0.0);
                    }
                    assert_eq!(circuit.time_s(), (interval - 1) as f64 * dt);
                    near(circuit.charges_c()[0], q);
                    (circuit, _) = next.into_parts();
                    q = next_q;
                }
            }
        }
    }
}

#[test]
fn moving_rlc_keeps_authored_ports_graph_heat_and_both_original_storage_coordinates() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let cards = sources(coordinate);
    let dt = 0.1;
    let source_current = 0.4;
    let graph_ohms = 7.0;
    for speed in [0.5, 5.0] {
        let (velocity, flux_density) = motion(speed);
        let drive = ConductionMotion {
            velocity_m_s: &velocity,
            normal_flux_density_t: &flux_density,
        };
        let mut circuit = ConductionResponseCircuit::new_with_motion(
            problem(&cards, coordinate).with_relaxation(0.5).unwrap(),
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (1, 2, Branch::Resistor { ohms: 999.0 }),
                    (2, 0, Branch::Inductor { henries: L }),
                    (1, 0, Branch::Capacitor { farads: C }),
                    (1, 0, Branch::CurrentSource { port: 0 }),
                    (1, 0, Branch::Resistor { ohms: graph_ohms }),
                ],
                transformers: vec![],
            },
            terminal(),
            0.0,
            vec![0.08],
            vec![3.0],
            &[source_current],
            drive,
            limits(),
            field_controls(),
            controls(),
            300,
            || false,
        )
        .unwrap();
        let (mut q, mut flux) = (3.0, 0.08);
        for interval in 1..=2 {
            let a = dt / (2.0 * C);
            let shunt = 1.0 + a / graph_ohms;
            let extra = 2.0 * L / dt + a / shunt;
            let rhs = (q / C + a * source_current) / shunt + 2.0 * flux / dt - 2.0 * speed;
            let j = current(rhs, extra, coordinate);
            let ep = fields(j, coordinate);
            let capacitor_voltage = (q / C + a * (source_current - j)) / shunt;
            let next_q = q + dt * (source_current - j - capacitor_voltage / graph_ohms);
            let next_flux = 2.0 * L * j - flux;
            let stored = (next_q - q) * (next_q + q) / (2.0 * C)
                + (next_flux - flux) * (next_flux + flux) / (2.0 * L);
            let graph_heat = dt * capacitor_voltage * capacitor_voltage / graph_ohms;
            let spatial_heat = dt * j * (ep[0] + ep[1]);
            let mechanical = dt * 2.0 * speed * j;
            let supply = dt * source_current * capacitor_voltage;
            let next = circuit
                .step_with_motion_candidate(
                    interval as f64 * dt,
                    &[source_current],
                    drive,
                    limits(),
                    field_controls(),
                    controls(),
                    config().circuit,
                    300,
                    || false,
                )
                .unwrap();
            let report = next.step();
            let native = &report.circuit;
            near(next.circuit().charges_c()[0], next_q);
            near(next.circuit().fluxes_wb()[0], next_flux);
            near(native.spatial_current_a, j);
            near(native.spatial_voltage_v, ep[0] + ep[1] + 2.0 * speed);
            near(native.stored_change_j, stored);
            near(native.supplied_work_j, supply);
            near(native.graph_resistors[0].heat_j, graph_heat);
            near(native.heat_j, graph_heat + spatial_heat);
            near(report.mechanical_work_j, mechanical);
            near(
                stored + spatial_heat + graph_heat + mechanical - supply,
                0.0,
            );
            near(native.energy_residual_j, 0.0);
            assert_eq!(next.circuit().graph().branches.len(), 5);
            assert_eq!(native.native.y.len(), 2);
            (circuit, _) = next.into_parts();
            q = next_q;
            flux = next_flux;
        }
    }
}

#[test]
fn moving_circuit_work_bounds_and_cancellation_leave_charge_flux_sources_and_time_reusable() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let cards = sources(coordinate);
    let (velocity, flux_density) = motion(5.0);
    let drive = ConductionMotion {
        velocity_m_s: &velocity,
        normal_flux_density_t: &flux_density,
    };
    let build = |maximum| {
        ConductionResponseCircuit::new_with_motion(
            problem(&cards, coordinate).with_relaxation(0.5).unwrap(),
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Resistor { ohms: 999.0 }),
                    (1, 0, Branch::Capacitor { farads: C }),
                ],
                transformers: vec![],
            },
            terminal(),
            0.0,
            vec![],
            vec![3.0],
            &[],
            drive,
            limits(),
            field_controls(),
            controls(),
            maximum,
            || false,
        )
    };
    assert!(matches!(build(2), Err(CircuitEnthalpyError::CircuitBudget)));
    let circuit = build(300).unwrap();
    let original_sources = circuit.conductors().to_vec();
    let mut calls = 0;
    let baseline = circuit
        .step_with_motion_candidate(
            0.1,
            &[],
            drive,
            limits(),
            field_controls(),
            controls(),
            config().circuit,
            300,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(matches!(
        circuit.step_with_motion_candidate(
            0.1,
            &[],
            drive,
            limits(),
            field_controls(),
            controls(),
            config().circuit,
            baseline.step().field_solves - 1,
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    for cutoff in [1, calls / 2, calls] {
        let mut polls = 0;
        assert!(matches!(
            circuit.step_with_motion_candidate(
                0.1,
                &[],
                drive,
                limits(),
                field_controls(),
                controls(),
                config().circuit,
                300,
                || {
                    polls += 1;
                    polls == cutoff
                },
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
    }
    let retry = circuit
        .step_with_motion_candidate(
            0.1,
            &[],
            drive,
            limits(),
            field_controls(),
            controls(),
            config().circuit,
            300,
            || false,
        )
        .unwrap();
    assert_eq!(retry.circuit().charges_c(), baseline.circuit().charges_c());
    assert_eq!(circuit.charges_c(), &[3.0]);
    assert_eq!(circuit.time_s(), 0.0);
    assert_eq!(circuit.conductors(), original_sources);
}
