use super::*;
use fs_material::gas::GasSpec;

fn duct(walls: bool, pressure: f64) -> DgDuct {
    duct_with_boundaries(walls, pressure, DgBoundary::Port(PortId(7)), DgBoundary::IdealOpen)
}

fn two_port_duct(walls: bool, pressure: f64) -> DgDuct {
    duct_with_boundaries(walls, pressure, DgBoundary::Port(PortId(7)), DgBoundary::Port(PortId(3)))
}

fn duct_with_boundaries(
    walls: bool, pressure: f64, inlet: DgBoundary, outlet: DgBoundary,
) -> DgDuct {
    let source = Duct {
        segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }],
    };
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let duct = DgDuct::new(
        &source, &gas, DgMesh { elements: 2, order: 2 }, inlet, outlet,
    ).unwrap();
    if walls {
        duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
    } else {
        duct
    }
}

#[test]
fn affine_port_predicts_the_actual_step_mean_and_integrated_work() {
    for walls in [false, true] {
        let duct = duct(walls, 101325.0);
        let mut state = duct.project_initial_state(|x| 0.3 * det::cos(30.0 * x), |_| 0.0);
        for dt in [1e-5, 2e-5] {
            let stepper = duct.propagator(dt, 1e-11).unwrap();
            for flow in [1e-8, -2e-8, 0.0, 3e-8] {
                let before = state.clone();
                let relation = stepper.port_response(&state, PortId(7)).unwrap();
                assert_eq!(
                    state, before,
                    "a boundary preview must not advance any memory"
                );
                assert!(relation.impedance_pa_s_m3 > 0.0);
                let expected = relation.history_pressure_pa + relation.impedance_pa_s_m3 * flow;
                let output = stepper
                    .step(
                        &mut state,
                        &[PortInput {
                            id: PortId(7),
                            effort: flow,
                        }],
                    )
                    .unwrap();
                let port = output.port_outputs[0];
                let scale =
                    relation.history_pressure_pa.abs() + (relation.impedance_pa_s_m3 * flow).abs();
                assert!((port.conjugate_rate - expected).abs() < 1e-11 * scale.max(1e-30));
                assert!(
                    (output.ledger.port_work_step_j - expected * flow * dt).abs()
                        < 1e-11 * (scale * flow.abs() * dt).max(1e-30)
                );
            }
        }
    }
}

#[test]
fn affine_port_refuses_unknown_ports_and_foreign_material_state() {
    let duct = duct(true, 101325.0);
    let stepper = duct.propagator(1e-5, 1e-11).unwrap();
    let state = duct.zero_state();
    assert!(stepper.port_response(&state, PortId(8)).is_err());
    let foreign = self::duct(true, 60000.0).zero_state();
    assert!(stepper.port_response(&foreign, PortId(7)).is_err());
    assert_eq!(state.step_count(), 0);
}

#[test]
fn multiport_predicts_simultaneous_pressure_and_work_in_both_kernels() {
    for walls in [false, true] {
        let duct = two_port_duct(walls, 101325.0);
        let mut state = duct.project_initial_state(
            |x| 0.3 * det::cos(30.0 * x),
            |x| 1e-9 * det::sin(20.0 * x),
        );
        for dt in [1e-5, 2e-5] {
            let stepper = duct.propagator(dt, 1e-11).unwrap();
            for flow in [
                [1e-8, 0.0], [0.0, -2e-8], [2e-8, -3e-8],
                [-1e-8, -2e-8], [0.0, 0.0],
            ] {
                let before = state.clone();
                let response = stepper.multiport_response(&state).unwrap();
                assert_eq!(state, before, "a multiport preview must not advance memory");
                assert_eq!(response.ports, vec![PortId(7), PortId(3)]);
                assert_eq!(response.history_pressure_pa.len(), 2);
                assert_eq!(response.impedance_pa_s_m3.len(), 4);
                let mut expected = [0.0; 2];
                let mut scales = [0.0; 2];
                for row in 0..2 {
                    let scalar = stepper.port_response(&state, response.ports[row]).unwrap();
                    assert_eq!(response.history_pressure_pa[row].to_bits(),
                        scalar.history_pressure_pa.to_bits());
                    assert_eq!(response.impedance_pa_s_m3[2 * row + row].to_bits(),
                        scalar.impedance_pa_s_m3.to_bits());
                    expected[row] = response.history_pressure_pa[row];
                    scales[row] = expected[row].abs();
                    for col in 0..2 {
                        let contribution = response.impedance_pa_s_m3[2 * row + col] * flow[col];
                        expected[row] += contribution;
                        scales[row] += contribution.abs();
                    }
                }
                // Reverse input order: PortId, not caller order, owns the column.
                let output = stepper.step(&mut state, &[
                    PortInput { id: PortId(3), effort: flow[1] },
                    PortInput { id: PortId(7), effort: flow[0] },
                ]).unwrap();
                for row in 0..2 {
                    assert_eq!(output.port_outputs[row].id, response.ports[row]);
                    assert!((output.port_outputs[row].conjugate_rate - expected[row]).abs()
                        <= 1e-11 * scales[row].max(1e-30));
                }
                let expected_work = dt * (flow[0] * expected[0] + flow[1] * expected[1]);
                let work_scale = dt * (flow[0].abs() * scales[0] + flow[1].abs() * scales[1]);
                assert!((output.ledger.port_work_step_j - expected_work).abs()
                    <= 1e-11 * work_scale.max(1e-30));
            }
        }
    }
}

#[test]
fn multiport_closes_a_passive_terminal_against_the_actual_same_step() {
    for walls in [false, true] {
        let duct = two_port_duct(walls, 101325.0);
        let stepper = duct.propagator(2e-5, 1e-11).unwrap();
        let mut state = duct.project_initial_state(|_| 0.2, |_| 0.0);
        let resistance = 1e8;
        for inlet_flow in [1e-8, 2e-8, 0.0, -1e-8] {
            let response = stepper.multiport_response(&state).unwrap();
            let z = &response.impedance_pa_s_m3;
            assert!(z[1] != 0.0 && z[2] != 0.0);
            // The terminal's outward flow is the negative of the DG inward flow.
            let outlet_flow =
                -(response.history_pressure_pa[1] + z[2] * inlet_flow) / (z[3] + resistance);
            let output = stepper.step(&mut state, &[
                PortInput { id: PortId(7), effort: inlet_flow },
                PortInput { id: PortId(3), effort: outlet_flow },
            ]).unwrap();
            let actual_outlet = output.port_outputs[1].conjugate_rate;
            let load_pressure = -resistance * outlet_flow;
            let scale = response.history_pressure_pa[1].abs()
                + (z[2] * inlet_flow).abs() + (z[3] * outlet_flow).abs() + load_pressure.abs();
            assert!((actual_outlet - load_pressure).abs() <= 1e-11 * scale.max(1e-30));
            let inlet_work = stepper.dt() * inlet_flow * output.port_outputs[0].conjugate_rate;
            let terminal_work = stepper.dt() * resistance * outlet_flow * outlet_flow;
            assert!((output.ledger.port_work_step_j - (inlet_work - terminal_work)).abs()
                <= 1e-11 * (inlet_work.abs() + terminal_work).max(1e-30));
        }
    }
}

#[test]
fn multiport_refuses_foreign_material_state_and_handles_no_ports() {
    for walls in [false, true] {
        let duct = two_port_duct(walls, 101325.0);
        let stepper = duct.propagator(1e-5, 1e-11).unwrap();
        let foreign = two_port_duct(walls, 60000.0).zero_state();
        assert!(stepper.multiport_response(&foreign).is_err());
        assert_eq!(foreign.step_count(), 0);
        let closed = duct_with_boundaries(
            walls, 101325.0, DgBoundary::Closed, DgBoundary::Closed,
        );
        let closed_stepper = closed.propagator(1e-5, 1e-11).unwrap();
        let state = closed.zero_state();
        let response = closed_stepper.multiport_response(&state).unwrap();
        assert!(response.ports.is_empty());
        assert!(response.history_pressure_pa.is_empty());
        assert!(response.impedance_pa_s_m3.is_empty());
        assert_eq!(state.step_count(), 0);
    }
}
