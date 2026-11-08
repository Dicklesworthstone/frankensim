//! Independently excited induction supplies real passive electrode loads.
use super::*;

struct CubeOracle {
    potential: C64,
    voltage: C64,
    impedance: C64,
    heat: f64,
    impressed_supply: C64,
    terminal_supply: C64,
}

fn cube_oracle(omega: f64, forcing: C64, current: C64) -> CubeOracle {
    // Direct tetrahedral integration: int W=(1,1,1)/6, int |W|^2=1/5,
    // K=2 and sigma=5. A uniform impressed [j,j,j] supplies F=j/2.
    // Current enters the x=1 electrode; phi=V*x. Eliminating V gives
    // D*A=F-I/6, V=I/5+i*omega*A/6, including the open-circuit EMF.
    let d = C64::new(2.0, -omega * 5.0 * (1.0 / 5.0 - 1.0 / 36.0));
    let potential = (forcing - current.scale(1.0 / 6.0)) / d;
    let induced = C64::new(0.0, omega) * potential;
    let voltage = current.scale(0.2) + induced.scale(1.0 / 6.0);
    CubeOracle {
        potential,
        voltage,
        impedance: C64::from_re(0.2) + C64::new(0.0, -omega / 36.0) / d,
        heat: 2.5
            * (voltage.norm_sq() - (voltage * induced.conj()).re / 3.0 + induced.norm_sq() / 5.0),
        impressed_supply: (-induced * forcing.conj()).scale(0.5),
        terminal_supply: (voltage * current.conj()).scale(0.5),
    }
}

#[test]
fn open_electrodes_resolve_independent_faraday_voltage_and_reverse_its_orientation() {
    let (problem, sources) = cube(1, None);
    let j = C64::new(1.6, -0.4);
    let impressed = [[j; 3]; 6];
    let velocity = [[0.0; 3]; 8];
    for reverse in [false, true] {
        let mut faces = electrodes(&problem);
        if reverse {
            faces.reverse();
        }
        let source = current_source(&problem, sources.clone(), &faces);
        for omega in [0.3, 3.0, 40.0] {
            let frequency = SpatialHarmonicTerminalFrequency::with_impressed_current(
                &source,
                &impressed,
                omega,
                &velocity,
                terminal_controls(),
                &mut || false,
            )
            .unwrap();
            let expected = cube_oracle(omega, j.scale(0.5), C64::ZERO);
            let voltage = expected.voltage.scale(if reverse { -1.0 } else { 1.0 });
            complex_near(frequency.open_circuit_voltages_v()[0], voltage, 3e-9);
            complex_near(frequency.impedance_ohm()[0], expected.impedance, 3e-9);
            assert_eq!(frequency.impressed_current_a_m2(), &impressed);
            assert!(frequency.impressed_report().is_some());
            let response = frequency
                .current_response(&[C64::ZERO], &mut || false)
                .unwrap();
            complex_near(response.currents_a[0], C64::ZERO, 3e-9);
            complex_near(response.terminal_voltages_v[0], voltage, 3e-9);
            complex_near(
                response.field.potential_wb()[diagonal(&problem)],
                expected.potential,
                3e-9,
            );
            let balance = response.field.balance();
            near(balance.mean_ohmic_power_w, expected.heat, 3e-9);
            near(
                balance.mean_impressed_power_w,
                expected.impressed_supply.re,
                3e-9,
            );
            near(
                balance.reactive_impressed_power_var,
                expected.impressed_supply.im,
                3e-9,
            );
            near(balance.mean_terminal_power_w, 0.0, 3e-9);
            near(
                balance.mean_magnetic_energy_j,
                0.5 * expected.potential.norm_sq(),
                3e-9,
            );
        }
        assert_eq!(source.time_s(), 0.0);
        assert!(source.potential_wb().iter().all(|&a| a == 0.0));
    }
}

#[test]
fn impressed_induction_powers_passive_resistor_and_rlc_loads_with_separate_source_accounts() {
    let (problem, sources) = cube(1, None);
    let source = current_source(&problem, sources, &electrodes(&problem));
    let j = C64::new(1.6, -0.4);
    let velocity = [[0.0; 3]; 8];
    for omega in [3.0, 40.0] {
        let frequency = SpatialHarmonicTerminalFrequency::with_impressed_current(
            &source,
            &[[j; 3]; 6],
            omega,
            &velocity,
            terminal_controls(),
            &mut || false,
        )
        .unwrap();
        let open = cube_oracle(omega, j.scale(0.5), C64::ZERO);
        for reactive in [false, true] {
            let (graph, branch, load) = if reactive {
                (
                    CircuitGraph {
                        node_count: 4,
                        branches: vec![
                            (1, 2, Branch::Resistor { ohms: 0.4 }),
                            (2, 3, Branch::Inductor { henries: 0.03 }),
                            (3, 0, Branch::Capacitor { farads: 0.7 }),
                            (1, 0, Branch::Inductor { henries: f64::NAN }),
                        ],
                        transformers: vec![],
                    },
                    3,
                    C64::new(0.4, 1.0 / (omega * 0.7) - omega * 0.03),
                )
            } else {
                (
                    CircuitGraph {
                        node_count: 2,
                        branches: vec![
                            (1, 0, Branch::Resistor { ohms: 0.4 }),
                            (1, 0, Branch::Inductor { henries: f64::NAN }),
                        ],
                        transformers: vec![],
                    },
                    1,
                    C64::from_re(0.4),
                )
            };
            // Both branch orientations leave the same graph node. KCL gives
            // I_field=-I_load, hence I_field=-V_open/(Z_field+Z_load).
            let current = -open.voltage / (open.impedance + load);
            let expected = cube_oracle(omega, j.scale(0.5), current);
            for scale in [0.25, 4.0] {
                let response = frequency
                    .circuit_response(
                        &graph,
                        &[branch],
                        &[],
                        &[],
                        HarmonicCircuitControls {
                            reference_resistance_scale: scale,
                            ..circuit_controls()
                        },
                        &mut || false,
                    )
                    .unwrap();
                complex_near(response.field.currents_a[0], current, 3e-8);
                complex_near(
                    response.field.terminal_voltages_v[0],
                    expected.voltage,
                    3e-8,
                );
                complex_near(response.branches[0].current_a, -current, 3e-8);
                complex_near(
                    response.branches[branch].absorbed_power_va,
                    expected.terminal_supply,
                    3e-8,
                );
                assert!(response.branches[branch].absorbed_power_va.re < 0.0);
                assert!(response.sources.is_empty() && response.source_outputs.is_empty());
                complex_near(response.circuit_supplied_power_va, C64::ZERO, 3e-9);
                complex_near(
                    response.impressed_supplied_power_va,
                    expected.impressed_supply,
                    3e-8,
                );
                complex_near(response.supplied_power_va, expected.impressed_supply, 3e-8);
                near(response.mean_field_power_w, expected.heat, 3e-8);
                near(
                    response.mean_resistor_power_w,
                    0.2 * current.norm_sq(),
                    3e-8,
                );
                near(response.mean_field_mechanical_power_w, 0.0, 1e-15);
                let network_energy = if reactive {
                    0.25 * current.norm_sq() * (0.03 + 1.0 / (omega * omega * 0.7))
                } else {
                    0.0
                };
                near(response.mean_network_stored_energy_j, network_energy, 3e-8);
                complex_near(response.absorbed_power_va, expected.impressed_supply, 8e-8);
                complex_near(response.power_residual_va, C64::ZERO, 8e-8);
                assert!(response.maximum_kirchhoff_current_defect_a < 3e-8);
            }
        }
    }
}

#[test]
fn graph_and_impressed_sources_superpose_as_fields_with_original_interference_heat() {
    let (problem, sources) = cube(1, None);
    let source = current_source(&problem, sources, &electrodes(&problem));
    let j = C64::new(1.6, -0.4);
    let omega = 3.0;
    let velocity = [[0.0; 3]; 8];
    let graph = series_graph();
    let drive = C64::new(1.2, -0.7);
    let homogeneous = SpatialHarmonicTerminalFrequency::new(
        &source,
        omega,
        &velocity,
        terminal_controls(),
        &mut || false,
    )
    .unwrap();
    let graph_only = homogeneous
        .circuit_response(&graph, &[4], &[], &[drive], circuit_controls(), &mut || {
            false
        })
        .unwrap();
    let mut mixed = Vec::new();
    for sign in [1.0, -1.0] {
        let forcing = j.scale(0.5 * sign);
        let frequency = SpatialHarmonicTerminalFrequency::with_impressed_current(
            &source,
            &[[j.scale(sign); 3]; 6],
            omega,
            &velocity,
            terminal_controls(),
            &mut || false,
        )
        .unwrap();
        complex_near(
            frequency.impedance_ohm()[0],
            homogeneous.impedance_ohm()[0],
            1e-12,
        );
        let source_only = frequency
            .circuit_response(
                &graph,
                &[4],
                &[],
                &[C64::ZERO],
                circuit_controls(),
                &mut || false,
            )
            .unwrap();
        let response = frequency
            .circuit_response(&graph, &[4], &[], &[drive], circuit_controls(), &mut || {
                false
            })
            .unwrap();
        let open = cube_oracle(omega, forcing, C64::ZERO);
        let current = (drive - open.voltage)
            / (open.impedance + C64::new(0.4, 1.0 / (omega * 0.7) - omega * 0.03));
        let expected = cube_oracle(omega, forcing, current);
        complex_near(response.field.currents_a[0], current, 3e-8);
        complex_near(
            response.field.field.potential_wb()[diagonal(&problem)],
            expected.potential,
            3e-8,
        );
        complex_near(
            response.field.terminal_voltages_v[0],
            expected.voltage,
            3e-8,
        );
        complex_near(
            response.circuit_supplied_power_va,
            (drive * current.conj()).scale(0.5),
            3e-8,
        );
        complex_near(
            response.impressed_supplied_power_va,
            expected.impressed_supply,
            3e-8,
        );
        complex_near(
            response.supplied_power_va,
            response.circuit_supplied_power_va + response.impressed_supplied_power_va,
            3e-8,
        );
        near(response.mean_field_power_w, expected.heat, 3e-8);
        assert!(
            (response.mean_field_power_w
                - source_only.mean_field_power_w
                - graph_only.mean_field_power_w)
                .abs()
                > 1e-6
        );
        complex_near(
            response.field.currents_a[0],
            source_only.field.currents_a[0] + graph_only.field.currents_a[0],
            3e-8,
        );
        for ((&actual, &external), &induced) in response
            .field
            .field
            .potential_wb()
            .iter()
            .zip(graph_only.field.field.potential_wb())
            .zip(source_only.field.field.potential_wb())
        {
            complex_near(actual, external + induced, 3e-8);
        }
        complex_near(response.power_residual_va, C64::ZERO, 8e-8);
        mixed.push(response);
    }
    complex_near(
        mixed[0].field.currents_a[0] + mixed[1].field.currents_a[0],
        graph_only.field.currents_a[0].scale(2.0),
        3e-8,
    );
    complex_near(graph_only.impressed_supplied_power_va, C64::ZERO, 1e-15);
    assert_eq!(source.time_s(), 0.0);
}

#[test]
fn moving_three_electrode_affine_network_matches_independent_voltage_control() {
    let (problem, sources) = cube(2, None);
    let planes = electrodes(&problem);
    let first = planes[1][0];
    let second = *planes[1]
        .iter()
        .find(|&&face| {
            problem.mesh().faces[face]
                .iter()
                .all(|v| !problem.mesh().faces[first].contains(v))
        })
        .unwrap();
    let source = current_source(
        &problem,
        sources,
        &[planes[0].clone(), vec![first], vec![second]],
    );
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| 0.03 * x))
        .collect();
    let impressed = vec![[C64::new(0.4, 0.2); 3]; problem.mesh().tets.len()];
    let drive = [C64::new(0.9, -0.4), C64::new(-0.2, 0.3)];
    let voltage_owner = source
        .clone()
        .into_voltage_driven(preparation(), &mut || false)
        .unwrap();
    let direct = voltage_owner
        .solve_harmonic_voltage(
            &[C64::ZERO, drive[0], drive[1]],
            &impressed,
            &velocity,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let frequency = SpatialHarmonicTerminalFrequency::with_impressed_current(
        &source,
        &impressed,
        3.0,
        &velocity,
        terminal_controls(),
        &mut || false,
    )
    .unwrap();
    assert!(frequency.impedance_ohm()[1].abs() > 1e-6);
    assert!(frequency.impedance_ohm()[2].abs() > 1e-6);
    assert!(
        frequency
            .open_circuit_voltages_v()
            .iter()
            .any(|v| v.abs() > 1e-6)
    );
    // Both ports return to the same physical electrode. Original graph branch
    // order and external source indices deliberately differ from port order.
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::Inductor { henries: f64::NAN }),
            (1, 0, Branch::VoltageSource { port: 1 }),
            (1, 0, Branch::Inductor { henries: f64::NAN }),
            (2, 0, Branch::VoltageSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let response = frequency
        .circuit_response(
            &graph,
            &[2, 0],
            &[],
            &[drive[1], drive[0]],
            circuit_controls(),
            &mut || false,
        )
        .unwrap();
    let mut terminal_supply = C64::ZERO;
    for port in 0..2 {
        let expected_current = direct.terminals()[port + 1].incoming_current_a;
        complex_near(response.field.currents_a[port], expected_current, 3e-8);
        complex_near(response.field.terminal_voltages_v[port], drive[port], 3e-8);
        terminal_supply = terminal_supply + (drive[port] * expected_current.conj()).scale(0.5);
    }
    for (&actual, &expected) in response
        .field
        .field
        .potential_wb()
        .iter()
        .zip(direct.potential_wb())
    {
        complex_near(actual, expected, 3e-8);
    }
    for (actual, expected) in response.field.field.cells().iter().zip(direct.cells()) {
        for (a, e) in actual
            .current_density_a_m2
            .iter()
            .zip(&expected.current_density_a_m2)
        {
            for k in 0..3 {
                complex_near(a[k], e[k], 2e-7);
            }
        }
    }
    let balance = direct.balance();
    assert!(balance.mean_mechanical_power_w.abs() > 1e-9);
    assert!(balance.mean_impressed_power_w.abs() > 1e-9);
    near(
        response.mean_field_power_w,
        balance.mean_ohmic_power_w,
        8e-8,
    );
    near(
        response.mean_field_mechanical_power_w,
        balance.mean_mechanical_power_w,
        8e-8,
    );
    complex_near(response.circuit_supplied_power_va, terminal_supply, 8e-8);
    complex_near(
        response.impressed_supplied_power_va,
        C64::new(
            balance.mean_impressed_power_w,
            balance.reactive_impressed_power_var,
        ),
        8e-8,
    );
    near(
        response.supplied_power_va.re,
        balance.mean_ohmic_power_w + balance.mean_mechanical_power_w,
        8e-8,
    );
    near(response.mean_resistor_power_w, 0.0, 1e-15);
    complex_near(response.power_residual_va, C64::ZERO, 8e-8);
    assert!(response.maximum_kirchhoff_current_defect_a < 3e-8);
    assert_eq!(source.time_s(), 0.0);
}

#[test]
fn impressed_preparation_rejects_invalid_budgeted_and_cancelled_trials_without_mutating_source() {
    let (problem, sources) = cube(1, None);
    let source = current_source(&problem, sources, &electrodes(&problem));
    let before = source.potential_wb().to_vec();
    let velocity = [[0.0; 3]; 8];
    let impressed = [[C64::new(1.6, -0.4); 3]; 6];
    for invalid in [vec![], vec![[C64::new(f64::NAN, 0.0); 3]; 6]] {
        assert!(
            SpatialHarmonicTerminalFrequency::with_impressed_current(
                &source,
                &invalid,
                3.0,
                &velocity,
                terminal_controls(),
                &mut || false,
            )
            .is_err()
        );
    }
    assert!(
        SpatialHarmonicTerminalFrequency::with_impressed_current(
            &source,
            &impressed,
            3.0,
            &velocity,
            SpatialHarmonicTerminalControls {
                maximum_response_scalars: 0,
                ..terminal_controls()
            },
            &mut || false,
        )
        .is_err()
    );
    let mut polls = 0;
    let accepted = SpatialHarmonicTerminalFrequency::with_impressed_current(
        &source,
        &impressed,
        3.0,
        &velocity,
        terminal_controls(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 2);
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(
            SpatialHarmonicTerminalFrequency::with_impressed_current(
                &source,
                &impressed,
                3.0,
                &velocity,
                terminal_controls(),
                &mut || {
                    seen += 1;
                    seen >= stop
                },
            )
            .is_err()
        );
    }
    let replay = SpatialHarmonicTerminalFrequency::with_impressed_current(
        &source,
        &impressed,
        3.0,
        &velocity,
        terminal_controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        accepted.open_circuit_voltages_v(),
        replay.open_circuit_voltages_v()
    );
    assert_eq!(accepted.impedance_ohm(), replay.impedance_ohm());
    let zero = SpatialHarmonicTerminalFrequency::with_impressed_current(
        &source,
        &[[C64::ZERO; 3]; 6],
        3.0,
        &velocity,
        terminal_controls(),
        &mut || false,
    )
    .unwrap();
    assert!(zero.impressed_report().is_none());
    assert_eq!(zero.open_circuit_voltages_v(), &[C64::ZERO]);
    assert_eq!(source.potential_wb(), before);
    assert_eq!(source.time_s(), 0.0);
}
