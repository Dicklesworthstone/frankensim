//! Prescribed volume ports exercise the original pressure solve, including the
//! pressure required by a flow source and its signed work at a reservoir.

use super::*;

fn injection(node: usize, flow: f64) -> FlowInjection {
    FlowInjection {
        node,
        flow: VolumetricFlowRate::new(flow),
    }
}

fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        a.is_finite()
            && b.is_finite()
            && (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

fn check_ports(solution: &CarreauGraphSolution, expected_injections: &[f64]) {
    assert_eq!(solution.pressures.len(), expected_injections.len());
    let mut is_reservoir = vec![false; solution.pressures.len()];
    let mut reservoir_power = 0.0;
    let mut reservoir_flow = 0.0;
    for (boundary, returned) in solution.boundaries.iter().zip(&solution.reservoir_outflows) {
        is_reservoir[boundary.node] = true;
        let expected =
            solution.node_outflows[boundary.node].value() - expected_injections[boundary.node];
        near(returned.value(), expected, 1e-12, 2e-18);
        reservoir_power += boundary.pressure.value() * expected;
        reservoir_flow += expected;
    }
    let mut injection_power = 0.0;
    let mut residual_power = 0.0;
    for (node, (&injection, pressure)) in expected_injections
        .iter()
        .zip(&solution.pressures)
        .enumerate()
    {
        injection_power += pressure.value() * injection;
        if !is_reservoir[node] {
            let residual = solution.node_outflows[node].value() - injection;
            assert!(residual.abs() < 2e-14);
            residual_power += pressure.value() * residual;
        }
    }
    near(solution.injection_power_w, injection_power, 1e-12, 2e-18);
    near(solution.boundary_power_w, reservoir_power, 1e-12, 2e-18);
    near(solution.node_residual_power_w, residual_power, 1e-10, 2e-18);
    assert!((reservoir_flow + expected_injections.iter().sum::<f64>()).abs() < 2e-14);
    check_power(solution);
}

#[test]
fn g1_newtonian_demands_resolve_unknown_pressures_and_remaining_reservoir_exchange() {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("01", 0, 1, pipe(0.5, true, 1e4)),
            branch("12", 1, 2, pipe(1.0, true, 1e4)),
            branch("02", 0, 2, pipe(1.5, true, 1e4)),
        ],
    )
    .unwrap();
    let sources = [injection(0, 3e-6), injection(1, -1e-6)];
    let reservoirs = [boundary(2, 2000.0)];
    let result =
        with_cx(|cx| graph.solve_with_injections(&reservoirs, &sources, controls(), cx)).unwrap();
    // Independent grounded 2x2 system in pressure relative to node 2:
    // g * [[4/3,-1],[-1,3/2]] * x = [3e-6,-1e-6], determinant g^2.
    let g = core::f64::consts::PI * 0.01_f64.powi(4) / (8.0 * 3.0 * 0.5);
    let x0 = (1.5 * 3e-6 - 1e-6) / g;
    let x1 = (3e-6 - 4.0 / 3.0 * 1e-6) / g;
    close(result.pressures[0].value(), 2000.0 + x0, 2e-12);
    close(result.pressures[1].value(), 2000.0 + x1, 2e-12);
    for (actual, expected) in
        result
            .branches
            .iter()
            .zip([g * (x0 - x1), 0.5 * g * x1, g / 3.0 * x0])
    {
        close(actual.profile.flow().volumetric_flow_m3_s, expected, 2e-12);
    }
    close(result.reservoir_outflows[0].value(), -2e-6, 2e-12);
    assert_eq!(result.injections, sources);
    assert_eq!(result.newton_iterations, 1);
    check_ports(&result, &[3e-6, -1e-6, 0.0]);

    // A separate supplied volume port at fixed pressure changes only the net
    // reservoir exchange; it must not change the physical pipe solve.
    let additional = [sources[0], sources[1], injection(2, 5e-6)];
    let extra = with_cx(|cx| graph.solve_with_injections(&reservoirs, &additional, controls(), cx))
        .unwrap();
    assert_eq!(extra.pressures, result.pressures);
    assert_eq!(extra.branches, result.branches);
    assert_eq!(extra.newton_iterations, result.newton_iterations);
    close(extra.reservoir_outflows[0].value(), -7e-6, 2e-12);
    near(
        extra.injection_power_w - result.injection_power_w,
        2000.0 * 5e-6,
        2e-12,
        0.0,
    );
    near(
        extra.boundary_power_w - result.boundary_power_w,
        -2000.0 * 5e-6,
        2e-12,
        0.0,
    );
    assert_eq!(extra.dissipated_power_w, result.dissipated_power_w);
    check_ports(&extra, &[3e-6, -1e-6, 5e-6]);
}

#[test]
fn g1_nonlinear_prescribed_flow_recovers_required_pressure_split_and_differential_response() {
    let pipes = [pipe(0.2, false, 1e4), pipe(1.0, false, 1e4)];
    let graph = CarreauPipeGraph::new(
        2,
        vec![
            branch("short", 0, 1, pipes[0]),
            branch("long", 0, 1, pipes[1]),
        ],
    )
    .unwrap();
    let profiles = pipes.map(|pipe| pipe.profile_at_pressure(4000.0, controls().pipe).unwrap());
    let flow = profiles
        .iter()
        .map(|p| p.flow().volumetric_flow_m3_s)
        .sum::<f64>();
    let reservoirs = [boundary(1, 700.0)];
    let solve = |flow| {
        with_cx(|cx| {
            graph.solve_with_injections(&reservoirs, &[injection(0, flow)], controls(), cx)
        })
        .unwrap()
    };
    let result = solve(flow);
    close(result.pressures[0].value(), 4700.0, 2e-8);
    for (actual, expected) in result.branches.iter().zip(&profiles) {
        close(
            actual.profile.flow().volumetric_flow_m3_s,
            expected.flow().volumetric_flow_m3_s,
            2e-8,
        );
        assert_eq!(actual.profile.pipe(), expected.pipe());
    }
    let short = result.branches[0].profile.flow().volumetric_flow_m3_s;
    let long = result.branches[1].profile.flow().volumetric_flow_m3_s;
    assert!(
        short / long > 5.5,
        "a frozen Newtonian viscosity predicts only the 5:1 geometry ratio"
    );
    check_ports(&result, &[flow, 0.0]);

    let step = flow * 1e-4;
    let plus = solve(flow + step);
    let minus = solve(flow - step);
    let observed = (plus.pressures[0].value() - minus.pressures[0].value()) / (2.0 * step);
    let conductance: f64 = result
        .branches
        .iter()
        .map(|branch| branch.differential_conductance_m3_s_pa)
        .sum();
    close(observed, 1.0 / conductance, 2e-5);
}

#[test]
fn g3_balanced_flow_transfer_drives_and_reverses_original_nonlinear_loop() {
    let graph = nonlinear_bridge();
    let reservoirs = [boundary(0, 1600.0)];
    let flow = 4e-5;
    let forward_sources = [injection(0, flow), injection(3, -flow)];
    let forward =
        with_cx(|cx| graph.solve_with_injections(&reservoirs, &forward_sources, controls(), cx))
            .unwrap();
    for (actual, expected) in forward.pressures[..3].iter().zip([1600.0, 1000.0, 400.0]) {
        close(actual.value(), expected, 2e-8);
    }
    assert!(forward.pressures[3].value().abs() < 2e-5);
    assert!(forward.reservoir_outflows[0].value().abs() < 2e-14);
    close(forward.injection_power_w, 0.064, 2e-8);
    close(forward.dissipated_power_w, 0.064, 2e-8);
    check_ports(&forward, &[flow, 0.0, 0.0, -flow]);

    let reverse_sources = [injection(0, -flow), injection(3, flow)];
    let reverse =
        with_cx(|cx| graph.solve_with_injections(&reservoirs, &reverse_sources, controls(), cx))
            .unwrap();
    for (node, pressure) in reverse.pressures.iter().enumerate() {
        near(
            pressure.value(),
            3200.0 - forward.pressures[node].value(),
            2e-8,
            2e-5,
        );
    }
    for (a, b) in forward.branches.iter().zip(&reverse.branches) {
        close(
            a.profile.flow().volumetric_flow_m3_s,
            -b.profile.flow().volumetric_flow_m3_s,
            2e-8,
        );
    }
    close(forward.dissipated_power_w, reverse.dissipated_power_w, 2e-8);
    close(forward.injection_power_w, reverse.injection_power_w, 2e-8);
    check_ports(&reverse, &[-flow, 0.0, 0.0, flow]);
}

#[test]
fn g3_prescribed_transfer_can_absorb_pressure_work_at_fixed_pressure_ports() {
    let graph =
        CarreauPipeGraph::new(2, vec![branch("pipe", 0, 1, pipe(0.5, false, 1e4))]).unwrap();
    let reservoirs = [boundary(0, 2000.0), boundary(1, 1000.0)];
    let flow = 2e-5;
    // Transfer fluid from the higher-pressure reservoir to the lower-pressure
    // reservoir. Its ideal pressure port absorbs work; it does not add heat.
    let inputs = [injection(0, -flow), injection(1, flow)];
    let result =
        with_cx(|cx| graph.solve_with_injections(&reservoirs, &inputs, controls(), cx)).unwrap();
    let pipe_flow = result.branches[0].profile.flow().volumetric_flow_m3_s;
    close(
        result.reservoir_outflows[0].value(),
        pipe_flow + flow,
        2e-13,
    );
    close(
        result.reservoir_outflows[1].value(),
        -pipe_flow - flow,
        2e-13,
    );
    close(result.injection_power_w, -1000.0 * flow, 2e-13);
    assert!(result.boundary_power_w > result.dissipated_power_w);
    assert_eq!(result.newton_iterations, 0);
    check_ports(&result, &[-flow, flow]);
}

#[test]
fn g3_compensated_port_combination_retains_small_net_demand_and_exact_rest() {
    let graph = CarreauPipeGraph::new(2, vec![branch("pipe", 0, 1, pipe(0.5, true, 1e4))]).unwrap();
    let flow = 1e-5;
    let inputs = [injection(0, 1e20), injection(0, flow), injection(0, -1e20)];
    let result =
        with_cx(|cx| graph.solve_with_injections(&[boundary(1, 0.0)], &inputs, controls(), cx))
            .unwrap();
    let resistance = 8.0 * 3.0 * 0.5 / (core::f64::consts::PI * 0.01_f64.powi(4));
    close(result.pressures[0].value(), resistance * flow, 2e-12);
    close(result.node_outflows[0].value(), flow, 2e-12);
    assert_eq!(result.injections, inputs);
    check_ports(&result, &[flow, 0.0]);

    let cancelled = [inputs[0], inputs[2], injection(0, 0.0)];
    let rest = with_cx(|cx| {
        graph.solve_with_injections(&[boundary(1, 2000.0)], &cancelled, controls(), cx)
    })
    .unwrap();
    assert_eq!(rest.pressures, [Pressure::new(2000.0); 2]);
    assert_eq!(rest.branches[0].profile.flow().volumetric_flow_m3_s, 0.0);
    assert_eq!(rest.injection_power_w, 0.0);
    assert_eq!(rest.dissipated_power_w, 0.0);
    assert_eq!(rest.newton_iterations, 0);
    assert_eq!(rest.pressure_attempts, 0);
    check_ports(&rest, &[0.0, 0.0]);
}

#[test]
fn g3_flow_forcing_does_not_disturb_independent_large_gauge_components() {
    let strict = pipe(0.5, false, 1e-6);
    let graph = CarreauPipeGraph::new(
        6,
        vec![
            branch("positive-gauge", 0, 1, strict),
            branch("negative-gauge", 2, 3, strict),
            branch("active", 4, 5, strict),
        ],
    )
    .unwrap();
    let flow = 1e-13;
    let expected = strict
        .profile_at_flow(flow, controls().pipe, 1e-11, 64)
        .unwrap();
    let reservoirs = [boundary(0, 1e120), boundary(2, -1e120), boundary(5, 0.0)];
    let result = with_cx(|cx| {
        graph.solve_with_injections(&reservoirs, &[injection(4, flow)], controls(), cx)
    })
    .unwrap();
    assert_eq!(result.pressures[1].value(), 1e120);
    assert_eq!(result.pressures[3].value(), -1e120);
    for branch in &result.branches[..2] {
        assert_eq!(branch.profile.flow().volumetric_flow_m3_s, 0.0);
    }
    close(
        result.pressures[4].value(),
        expected.profile().flow().pressure_drop_pa,
        2e-8,
    );
    close(
        result.branches[2].profile.flow().volumetric_flow_m3_s,
        flow,
        2e-8,
    );
    assert!(result.branches[2].profile.flow().wall_shear_rate_per_s < 1e-6);
    check_ports(&result, &[0.0, 0.0, 0.0, 0.0, flow, 0.0]);
}

#[test]
fn g4_invalid_or_unsupported_flow_ports_refuse_and_valid_retry_preserves_source() {
    let limited = pipe(0.5, false, 1.0);
    let graph = CarreauPipeGraph::new(2, vec![branch("limited", 0, 1, limited)]).unwrap();
    let reservoirs = [boundary(1, 0.0)];
    for bad in [
        injection(2, 1e-5),
        injection(0, f64::NAN),
        injection(0, f64::INFINITY),
    ] {
        assert!(matches!(
            with_cx(|cx| graph.solve_with_injections(&reservoirs, &[bad], controls(), cx)),
            Err(CarreauGraphError::Graph(GraphError::InvalidInput(_)))
        ));
    }
    let viscosity = law(false)
        .evaluate(300.0, 1.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    let ceiling_pressure = 2.0 * 0.5 * viscosity / 0.01;
    let ceiling_flow = limited
        .profile_at_pressure(ceiling_pressure, controls().pipe)
        .unwrap()
        .flow()
        .volumetric_flow_m3_s;
    let bounded = CarreauGraphControls {
        maximum_newton_iterations: 8,
        maximum_pressure_attempts: 8,
        maximum_backtracks: 12,
        ..controls()
    };
    assert!(
        with_cx(|cx| graph.solve_with_injections(
            &reservoirs,
            &[injection(0, 2.0 * ceiling_flow)],
            bounded,
            cx
        ))
        .is_err()
    );
    assert_eq!(graph.branches()[0].pipe, limited);
    for flow in [-0.25 * ceiling_flow, 0.25 * ceiling_flow] {
        let result = with_cx(|cx| {
            graph.solve_with_injections(&reservoirs, &[injection(0, flow)], controls(), cx)
        })
        .unwrap();
        let inverse = limited
            .profile_at_flow(flow, controls().pipe, 1e-11, 64)
            .unwrap();
        close(
            result.pressures[0].value(),
            inverse.profile().flow().pressure_drop_pa,
            2e-7,
        );
        assert!(
            result.branches[0]
                .profile
                .flow()
                .wall_shear_rate_per_s
                .abs()
                <= 1.0
        );
        check_ports(&result, &[flow, 0.0]);
    }

    let gate = CancelGate::new();
    gate.request();
    let cancelled = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 711,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        graph.solve_with_injections(&reservoirs, &[injection(0, 1e-6)], controls(), &cx)
    });
    assert!(matches!(
        cancelled,
        Err(CarreauGraphError::Graph(GraphError::Interrupted))
    ));
}
