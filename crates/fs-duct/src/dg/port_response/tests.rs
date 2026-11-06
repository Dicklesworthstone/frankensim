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

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

/// The dense response exactly as it was before the cache: the product re-summed, the impedance
/// re-formed from the retained port columns, on every call.
fn reference_multiport(stepper: &DgPropagator, state: &DgState) -> (Vec<f64>, Vec<f64>) {
    let PropagatorKernel::DenseSeries { dt_phi1, port_cols_phi2, ports, .. } = &stepper.kernel else {
        panic!("a dense kernel");
    };
    let n = stepper.total_dofs;
    let count = ports.len();
    let mut integral = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            integral[i] += dt_phi1[i * n + j] * state.tilde[j];
        }
    }
    let mut history = vec![0.0; count];
    let mut impedance = vec![0.0; count * count];
    for (row, (_, b)) in ports.iter().enumerate() {
        for i in 0..n {
            history[row] += b[i] * integral[i];
            for col in 0..count {
                impedance[row * count + col] += b[i] * port_cols_phi2[col][i];
            }
        }
        history[row] /= stepper.dt;
        for col in 0..count {
            impedance[row * count + col] /= stepper.dt;
        }
    }
    (history, impedance)
}

/// The dense step as it was before the cache: one fused loop over the two matrices.
fn reference_step(stepper: &DgPropagator, state: &DgState, efforts: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let PropagatorKernel::DenseSeries { phi0, dt_phi1, port_cols_phi1, port_cols_phi2, .. } =
        &stepper.kernel
    else {
        panic!("a dense kernel");
    };
    let n = stepper.total_dofs;
    let mut next = vec![0.0; n];
    let mut int = vec![0.0; n];
    for i in 0..n {
        let mut acc_next = 0.0;
        let mut acc_int = 0.0;
        for j in 0..n {
            acc_next += phi0[i * n + j] * state.tilde[j];
            acc_int += dt_phi1[i * n + j] * state.tilde[j];
        }
        for (slot, &u) in efforts.iter().enumerate() {
            acc_next += port_cols_phi1[slot][i] * u;
            acc_int += port_cols_phi2[slot][i] * u;
        }
        next[i] = acc_next;
        int[i] = acc_int;
    }
    (next, int)
}

#[test]
fn the_dense_response_cache_delivers_the_uncached_bits_and_sums_each_state_once() {
    let duct = two_port_duct(true, 101325.0);
    let stepper = duct.propagator(1e-5, 1e-11).unwrap();
    assert!(matches!(stepper.kernel, PropagatorKernel::DenseSeries { .. }), "walls make it dense");
    let mut state = duct.project_initial_state(
        |x| 0.3 * det::cos(30.0 * x),
        |x| 1e-9 * det::sin(20.0 * x),
    );
    let ids = [PortId(7), PortId(3)];
    let flows = [[1e-8, 0.0], [0.0, -2e-8], [2e-8, -3e-8], [-1e-8, -2e-8], [0.0, 0.0], [3e-8, 1e-8]];
    assert_eq!(
        stepper.dense_response_report(),
        Some(DgDenseResponseReport { calls: 0, products: 0 })
    );
    let steps = 24usize;
    for k in 0..steps {
        let flow = flows[k % flows.len()];
        let inputs = [
            PortInput { id: PortId(7), effort: flow[0] },
            PortInput { id: PortId(3), effort: flow[1] },
        ];
        let (history_ref, impedance_ref) = reference_multiport(&stepper, &state);
        if k % 2 == 0 {
            // A valve iterating on one step: every response, repeatedly, then a trial step.
            for _ in 0..3 {
                let response = stepper.multiport_response(&state).unwrap();
                assert_eq!(bits(&response.history_pressure_pa), bits(&history_ref));
                assert_eq!(bits(&response.impedance_pa_s_m3), bits(&impedance_ref));
            }
            for (row, id) in ids.iter().enumerate() {
                let scalar = stepper.port_response(&state, *id).unwrap();
                assert_eq!(scalar.history_pressure_pa.to_bits(), history_ref[row].to_bits());
                assert_eq!(scalar.impedance_pa_s_m3.to_bits(), impedance_ref[2 * row + row].to_bits());
            }
            let mut trial = state.clone();
            stepper.step(&mut trial, &inputs).unwrap();
        }
        // Odd steps arrive cold: no response was asked for first.
        let (next_ref, int_ref) = reference_step(&stepper, &state, &flow);
        let output = stepper.step(&mut state, &inputs).unwrap();
        assert_eq!(bits(&state.tilde), bits(&next_ref), "step {k}");
        let PropagatorKernel::DenseSeries { ports, .. } = &stepper.kernel else { unreachable!() };
        for (slot, (_, b)) in ports.iter().enumerate() {
            let mut delta_x = 0.0;
            for i in 0..stepper.total_dofs {
                delta_x += b[i] * int_ref[i];
            }
            assert_eq!(output.port_outputs[slot].delta_x.to_bits(), delta_x.to_bits(), "step {k}");
        }
    }
    // One product per distinct state, however many times it was asked for: even steps ask 3 + 2 + 1
    // times before the real step, odd steps only at the step.
    let report = stepper.dense_response_report().unwrap();
    assert_eq!(report.products, steps as u64);
    assert_eq!(report.calls, 12 * 7 + 12);
    println!("DENSE-CACHE-MEASURE {steps} states: {report:?}");

    // A different dt is a different product, and a state last used at one dt is not served at another.
    let slow = duct.propagator(2e-5, 1e-11).unwrap();
    let first = stepper.multiport_response(&state).unwrap();
    let (history_fast, impedance_fast) = reference_multiport(&stepper, &state);
    assert_eq!(bits(&first.history_pressure_pa), bits(&history_fast));
    let on_slow = slow.multiport_response(&state).unwrap();
    let (history_slow, impedance_slow) = reference_multiport(&slow, &state);
    assert_eq!(bits(&on_slow.history_pressure_pa), bits(&history_slow));
    assert_eq!(bits(&on_slow.impedance_pa_s_m3), bits(&impedance_slow));
    assert_eq!(slow.dense_response_report().unwrap().products, 1);
    let again = stepper.multiport_response(&state).unwrap();
    assert_eq!(bits(&again.history_pressure_pa), bits(&history_fast));
    assert_eq!(bits(&again.impedance_pa_s_m3), bits(&impedance_fast));
    assert_eq!(stepper.dense_response_report().unwrap().products, steps as u64 + 2);
}

#[test]
fn the_sparse_response_forms_its_impedance_once_and_the_dense_cache_never_reaches_it() {
    let duct = two_port_duct(true, 101325.0);
    let dense = duct.propagator(1e-5, 1e-11).unwrap();
    // Two fresh sparse propagators: one solve of a propagator starts its convergence check from the
    // dimension of the previous solve, so only first solves are bit-comparable.
    let sparse_a = duct.build_krylov_propagator(1e-5, 1e-11).unwrap();
    let sparse_b = duct.build_krylov_propagator(1e-5, 1e-11).unwrap();
    let make = || duct.project_initial_state(|x| 0.3 * det::cos(30.0 * x), |x| 1e-9 * det::sin(20.0 * x));
    let warm = make();
    let cold = make();
    // The dense kernel fills the dense cache of `warm` and not that of `cold`.
    dense.multiport_response(&warm).unwrap();
    let on_warm = sparse_a.multiport_response(&warm).unwrap();
    let on_cold = sparse_b.multiport_response(&cold).unwrap();
    assert_eq!(bits(&on_warm.history_pressure_pa), bits(&on_cold.history_pressure_pa));
    assert_eq!(bits(&on_warm.impedance_pa_s_m3), bits(&on_cold.impedance_pa_s_m3));
    assert_eq!(sparse_a.dense_response_report(), None);
    assert_eq!(dense.dense_response_report().unwrap().products, 1);

    // The sparse response as it was before the build-time matrix: the impedance re-formed from the
    // retained port columns on every call, the history from the solve this state already holds.
    let PropagatorKernel::Krylov(kernel) = &sparse_a.kernel else { panic!("a sparse kernel") };
    let solved = kernel.solve_cached(sparse_a.dt, &warm).unwrap();
    let n = sparse_a.total_dofs;
    let count = kernel.ports.len();
    let mut history = vec![0.0; count];
    let mut impedance = vec![0.0; count * count];
    for (row, (_, b)) in kernel.ports.iter().enumerate() {
        for i in 0..n {
            history[row] += b[i] * solved.int[i];
            for col in 0..count {
                impedance[row * count + col] += b[i] * kernel.port_cols_phi2[col][i];
            }
        }
        history[row] /= sparse_a.dt;
        for col in 0..count {
            impedance[row * count + col] /= sparse_a.dt;
        }
    }
    assert_eq!(bits(&on_warm.history_pressure_pa), bits(&history));
    assert_eq!(bits(&on_warm.impedance_pa_s_m3), bits(&impedance));
    for (row, id) in [PortId(7), PortId(3)].iter().enumerate() {
        let scalar = sparse_a.port_response(&warm, *id).unwrap();
        assert_eq!(scalar.history_pressure_pa.to_bits(), history[row].to_bits());
        assert_eq!(scalar.impedance_pa_s_m3.to_bits(), impedance[2 * row + row].to_bits());
    }
}
