//! Synthetic rheological fixtures, not material calibration. The Newtonian
//! bridge oracle is solved independently in closed form; nonlinear targets
//! use prescribed pipe flows to manufacture a different, looped pressure solve.

use super::super::tests::with_cx;
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn law(newtonian: bool) -> ThermalCarreauYasuda {
    let (zero, infinity) = if newtonian { (3.0, 3.0) } else { (12.0, 0.2) };
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(zero, infinity, 0.6, 0.45, 2.0, [0.0, 1e7]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::Vertical,
    )
}

fn pipe(length: f64, newtonian: bool, ceiling: f64) -> LaminarCarreauPipe {
    LaminarCarreauPipe::new(law(newtonian), 300.0, 0.01, length, ceiling).unwrap()
}

fn branch(name: &str, from: usize, to: usize, pipe: LaminarCarreauPipe) -> CarreauPipeBranch {
    CarreauPipeBranch {
        name: name.into(),
        from,
        to,
        pipe,
    }
}

fn boundary(node: usize, pressure: f64) -> FixedPressure {
    FixedPressure {
        node,
        pressure: Pressure::new(pressure),
    }
}

fn controls() -> CarreauGraphControls {
    CarreauGraphControls {
        absolute_flow_tolerance: VolumetricFlowRate::new(2e-15),
        relative_flow_tolerance: 1e-10,
        pipe: CarreauPipeControls {
            relative_quadrature_tolerance: 1e-10,
            ..CarreauPipeControls::default()
        },
        ..CarreauGraphControls::default()
    }
}

fn close(a: f64, b: f64, relative: f64) {
    let scale = a.abs().max(b.abs());
    assert!(
        a.is_finite() && b.is_finite() && (scale == 0.0 || ((a - b) / scale).abs() <= relative),
        "{a:e} != {b:e}"
    );
}

fn check_power(solution: &CarreauGraphSolution) {
    let requested_work: f64 = solution
        .branches
        .iter()
        .map(|branch| {
            let flow = branch.profile.flow();
            close(
                flow.requested_pressure_drop_pa,
                solution.pressures[branch.from].value() - solution.pressures[branch.to].value(),
                1e-13,
            );
            assert!(flow.dissipated_power_w >= 0.0);
            flow.requested_pressure_drop_pa * flow.volumetric_flow_m3_s
        })
        .sum();
    let all_node_work: f64 = solution
        .pressures
        .iter()
        .zip(&solution.node_outflows)
        .map(|(p, q)| p.value() * q.value())
        .sum();
    close(all_node_work, requested_work, 1e-12);
    close(
        solution.dissipated_power_w,
        requested_work + solution.constitutive_residual_power_w,
        1e-12,
    );
    assert!(
        solution.power_balance_residual_w.abs() <= 2e-13 * solution.dissipated_power_w.max(1e-12)
    );
}

#[test]
fn g1_newtonian_bridge_matches_independent_grounded_network_and_power() {
    let lengths = [1.0, 2.0, 3.0, 4.0, 1.5];
    let edges = [(0, 1), (1, 3), (2, 1), (0, 2), (2, 3)];
    let branches = edges
        .into_iter()
        .zip(lengths)
        .enumerate()
        .map(|(index, ((from, to), factor))| {
            branch(
                &format!("pipe-{index}"),
                from,
                to,
                pipe(0.5 * factor, true, 1e4),
            )
        })
        .collect();
    let graph = CarreauPipeGraph::new(4, branches).unwrap();
    let pressure = 1200.0;
    let g = core::f64::consts::PI * 0.01_f64.powi(4) / (8.0 * 3.0 * 0.5);
    let [g0, g1, g2, g3, g4] = lengths.map(|length| g / length);
    let (a, c) = (g0 + g1 + g2, g3 + g4 + g2);
    let determinant = a * c - g2 * g2;
    let p1 = pressure * (g0 * c + g2 * g3) / determinant;
    let p2 = pressure * (g3 * a + g2 * g0) / determinant;
    let expected_pressures = [pressure, p1, p2, 0.0];
    let expected_flows = [
        g0 * (pressure - p1),
        g1 * p1,
        g2 * (p2 - p1),
        g3 * (pressure - p2),
        g4 * p2,
    ];
    let result =
        with_cx(|cx| graph.solve(&[boundary(0, pressure), boundary(3, 0.0)], controls(), cx))
            .unwrap();
    for (actual, expected) in result.pressures.iter().zip(expected_pressures) {
        close(actual.value(), expected, 2e-12);
    }
    for (index, (actual, expected)) in result.branches.iter().zip(expected_flows).enumerate() {
        close(actual.profile.flow().volumetric_flow_m3_s, expected, 2e-12);
        close(
            actual.differential_conductance_m3_s_pa,
            g / lengths[index],
            2e-14,
        );
        assert_eq!(actual.profile.pipe(), &graph.branches()[index].pipe);
    }
    assert!(result.branches[2].profile.flow().volumetric_flow_m3_s < 0.0);
    close(
        result.boundary_power_w,
        pressure * (expected_flows[0] + expected_flows[3]),
        2e-12,
    );
    assert!(result.max_node_imbalance.value() < 1e-18);
    assert_eq!(
        result.newton_iterations, 1,
        "the true constant-viscosity Jacobian solves this bridge once"
    );
    assert_eq!(result.pressure_attempts, 1);
    check_power(&result);
}

fn nonlinear_bridge() -> CarreauPipeGraph {
    let base = pipe(1.0, false, 1e4);
    let specifications = [
        (0, 1, 600.0, 3e-5),
        (1, 3, 1000.0, 2e-5),
        (2, 1, -600.0, -1e-5),
        (0, 2, 1200.0, 1e-5),
        (2, 3, 400.0, 2e-5),
    ];
    let branches = specifications
        .into_iter()
        .enumerate()
        .map(|(index, (from, to, pressure, flow))| {
            let inverse = base
                .profile_at_flow(flow, controls().pipe, 1e-11, 96)
                .unwrap();
            let length = pressure / inverse.profile().flow().pressure_drop_pa;
            branch(
                &format!("nonlinear-{index}"),
                from,
                to,
                pipe(length, false, 1e4),
            )
        })
        .collect();
    CarreauPipeGraph::new(4, branches).unwrap()
}

#[test]
fn g1_g3_original_shear_thinning_material_closes_loop_continuity_and_signed_power() {
    let graph = nonlinear_bridge();
    let result =
        with_cx(|cx| graph.solve(&[boundary(0, 1600.0), boundary(3, 0.0)], controls(), cx))
            .unwrap();
    for (actual, expected) in result.pressures.iter().zip([1600.0, 1000.0, 400.0, 0.0]) {
        close(actual.value(), expected, 2e-8);
    }
    for (index, (actual, expected)) in result
        .branches
        .iter()
        .zip([3e-5, 2e-5, -1e-5, 1e-5, 2e-5])
        .enumerate()
    {
        let flow = actual.profile.flow();
        close(flow.volumetric_flow_m3_s, expected, 2e-8);
        assert_eq!(actual.profile.pipe(), &graph.branches()[index].pipe);
        // A frozen low-shear Poiseuille resistance and a secant Jacobian both
        // miss the actual differential response at this operating point.
        let secant = flow.volumetric_flow_m3_s / flow.pressure_drop_pa;
        assert!(actual.differential_conductance_m3_s_pa > 1.15 * secant);
        let original = actual.profile.pipe();
        let frozen = core::f64::consts::PI * original.radius_m().powi(4)
            / (8.0 * 12.0 * original.length_m());
        assert!(secant > 2.0 * frozen);
    }
    close(result.node_outflows[0].value(), 4e-5, 2e-8);
    close(result.node_outflows[3].value(), -4e-5, 2e-8);
    assert!(result.max_node_imbalance.value() < 1e-14);
    close(result.boundary_power_w, 0.064, 2e-8);
    check_power(&result);

    let reverse =
        with_cx(|cx| graph.solve(&[boundary(0, -1600.0), boundary(3, 0.0)], controls(), cx))
            .unwrap();
    for (forward, backward) in result.branches.iter().zip(&reverse.branches) {
        close(
            forward.profile.flow().volumetric_flow_m3_s,
            -backward.profile.flow().volumetric_flow_m3_s,
            2e-10,
        );
        close(
            forward.differential_conductance_m3_s_pa,
            backward.differential_conductance_m3_s_pa,
            2e-10,
        );
    }
    close(result.dissipated_power_w, reverse.dissipated_power_w, 2e-10);
    close(result.boundary_power_w, reverse.boundary_power_w, 2e-10);
    check_power(&reverse);
}

#[test]
fn g3_orientation_and_pipe_subdivision_preserve_actual_nonlinear_network() {
    let graph = nonlinear_bridge();
    let boundaries = [boundary(0, 1600.0), boundary(3, 0.0)];
    let original = with_cx(|cx| graph.solve(&boundaries, controls(), cx)).unwrap();
    let mut reversed = graph.branches().to_vec();
    (reversed[2].from, reversed[2].to) = (reversed[2].to, reversed[2].from);
    let reversed = CarreauPipeGraph::new(4, reversed).unwrap();
    let orientation = with_cx(|cx| reversed.solve(&boundaries, controls(), cx)).unwrap();
    for (index, (a, b)) in original
        .branches
        .iter()
        .zip(&orientation.branches)
        .enumerate()
    {
        close(
            a.profile.flow().volumetric_flow_m3_s,
            b.profile.flow().volumetric_flow_m3_s * if index == 2 { -1.0 } else { 1.0 },
            2e-10,
        );
    }

    let mut split = graph.branches()[1..].to_vec();
    let first = graph.branches()[0].pipe;
    let half = pipe(
        first.length_m() * 0.5,
        false,
        first.maximum_laminar_wall_rate_per_s(),
    );
    split.push(branch("first-half", 0, 4, half));
    split.push(branch("second-half", 4, 1, half));
    let split = CarreauPipeGraph::new(5, split).unwrap();
    let subdivided = with_cx(|cx| split.solve(&boundaries, controls(), cx)).unwrap();
    for node in 0..4 {
        close(
            original.pressures[node].value(),
            subdivided.pressures[node].value(),
            2e-8,
        );
    }
    close(subdivided.pressures[4].value(), 1300.0, 2e-8);
    for (a, b) in original.branches[1..].iter().zip(&subdivided.branches[..4]) {
        close(
            a.profile.flow().volumetric_flow_m3_s,
            b.profile.flow().volumetric_flow_m3_s,
            2e-8,
        );
    }
    for half in &subdivided.branches[4..] {
        close(
            half.profile.flow().volumetric_flow_m3_s,
            original.branches[0].profile.flow().volumetric_flow_m3_s,
            2e-8,
        );
    }
    close(
        original.dissipated_power_w,
        subdivided.dissipated_power_w,
        2e-8,
    );
    check_power(&subdivided);
}

#[test]
fn g3_pressure_continuation_recovers_supported_flow_after_unsupported_initial_guess() {
    let short = pipe(0.1, false, 10.0);
    let long = pipe(1.0, false, 10.0);
    let graph = CarreauPipeGraph::new(
        3,
        vec![branch("short", 0, 1, short), branch("long", 1, 2, long)],
    )
    .unwrap();
    // Both final pipes operate at exactly 5/s. The naive midpoint pressure
    // makes the short pipe exceed 10/s, even though the network is supported.
    let viscosity = law(false)
        .evaluate(300.0, 5.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    let total = 2.0 * 1.1 * viscosity * 5.0 / 0.01;
    assert!(
        short
            .profile_at_pressure(total * 0.5, controls().pipe)
            .is_err()
    );
    let result =
        with_cx(|cx| graph.solve(&[boundary(0, total), boundary(2, 0.0)], controls(), cx)).unwrap();
    close(result.pressures[1].value(), total / 1.1, 2e-9);
    for branch in &result.branches {
        close(branch.profile.flow().wall_shear_rate_per_s, 5.0, 2e-8);
    }
    assert!(result.pressure_attempts > 1);
    close(
        result.branches[0].profile.flow().volumetric_flow_m3_s,
        result.branches[1].profile.flow().volumetric_flow_m3_s,
        2e-9,
    );
    check_power(&result);
}

#[test]
fn g0_rest_and_anchored_disconnected_components_have_finite_true_conductance() {
    let fluid_pipe = pipe(0.5, false, 1e4);
    let graph = CarreauPipeGraph::new(
        4,
        vec![
            branch("left", 0, 1, fluid_pipe),
            branch("right", 2, 3, fluid_pipe),
        ],
    )
    .unwrap();
    let rest = with_cx(|cx| graph.solve(&[boundary(0, 500.0), boundary(2, 500.0)], controls(), cx))
        .unwrap();
    assert_eq!(rest.newton_iterations, 0);
    assert_eq!(rest.dissipated_power_w, 0.0);
    assert_eq!(rest.boundary_power_w, 0.0);
    for branch in &rest.branches {
        assert_eq!(branch.profile.flow().volumetric_flow_m3_s, 0.0);
        assert!(branch.differential_conductance_m3_s_pa > 0.0);
        assert!(branch.differential_conductance_m3_s_pa.is_finite());
    }
    let active = with_cx(|cx| {
        graph.solve(
            &[
                boundary(0, 100.0),
                boundary(1, 0.0),
                boundary(2, 200.0),
                boundary(3, 200.0),
            ],
            controls(),
            cx,
        )
    })
    .unwrap();
    assert!(active.branches[0].profile.flow().volumetric_flow_m3_s > 0.0);
    assert_eq!(active.branches[1].profile.flow().volumetric_flow_m3_s, 0.0);
    check_power(&active);
    assert!(matches!(
        with_cx(|cx| graph.solve(&[boundary(0, 0.0)], controls(), cx)),
        Err(CarreauGraphError::Graph(GraphError::UnreferencedNode(2)))
    ));
}

#[test]
fn g3_disconnected_pressure_gauges_do_not_consume_a_stagnant_pipes_rate_support() {
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
    let reservoirs = [
        boundary(0, 1e120),
        boundary(2, -1e120),
        boundary(4, 1e-5),
        boundary(5, 0.0),
    ];
    let result = with_cx(|cx| graph.solve(&reservoirs, controls(), cx)).unwrap();
    assert_eq!(result.pressures[1].value(), 1e120);
    assert_eq!(result.pressures[3].value(), -1e120);
    for branch in &result.branches[..2] {
        assert_eq!(branch.profile.flow().volumetric_flow_m3_s, 0.0);
        assert_eq!(branch.profile.flow().wall_shear_rate_per_s, 0.0);
    }
    assert!(result.branches[2].profile.flow().volumetric_flow_m3_s > 0.0);
    assert!(result.branches[2].profile.flow().wall_shear_rate_per_s < 1e-6);
    assert_eq!(result.newton_iterations, 0);
    assert_eq!(result.pressure_attempts, 1);
    check_power(&result);
}

#[test]
fn g4_source_ceiling_work_exhaustion_and_cancellation_publish_no_partial_network() {
    let limited = pipe(0.5, false, 1.0);
    let graph = CarreauPipeGraph::new(2, vec![branch("limited", 0, 1, limited)]).unwrap();
    let viscosity = law(false)
        .evaluate(300.0, 1.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    let ceiling = 2.0 * 0.5 * viscosity / 0.01;
    let bounded = CarreauGraphControls {
        maximum_pressure_attempts: 8,
        ..controls()
    };
    let failure =
        with_cx(|cx| graph.solve(&[boundary(0, 2.0 * ceiling), boundary(1, 0.0)], bounded, cx));
    match failure {
        Err(CarreauGraphError::PressureContinuationLimit {
            reached_fraction,
            attempts,
            ..
        }) => {
            assert!(reached_fraction < 1.0);
            assert_eq!(attempts, 8);
        }
        other => panic!("unsupported source must refuse: {other:?}"),
    }
    assert_eq!(graph.branches()[0].pipe, limited);
    let retry = with_cx(|cx| {
        graph.solve(
            &[boundary(0, 0.5 * ceiling), boundary(1, 0.0)],
            controls(),
            cx,
        )
    })
    .unwrap();
    assert!(retry.branches[0].profile.flow().wall_shear_rate_per_s < 1.0);
    check_power(&retry);

    let nonlinear = nonlinear_bridge();
    let short_budget = CarreauGraphControls {
        maximum_newton_iterations: 1,
        ..controls()
    };
    assert!(matches!(
        with_cx(|cx| nonlinear.solve(&[boundary(0, 1600.0), boundary(3, 0.0)], short_budget, cx)),
        Err(CarreauGraphError::NewtonIterationLimit { iterations: 1 })
    ));

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
        graph.solve(&[boundary(0, 100.0), boundary(1, 0.0)], controls(), &cx)
    });
    assert!(matches!(
        cancelled,
        Err(CarreauGraphError::Graph(GraphError::Interrupted))
    ));
}

#[test]
fn g0_invalid_topology_boundaries_controls_and_mixed_materials_refuse() {
    let physical = pipe(0.5, false, 1e4);
    assert!(CarreauPipeGraph::new(0, vec![]).is_err());
    assert!(CarreauPipeGraph::new(2, vec![]).is_err());
    assert!(CarreauPipeGraph::new(2, vec![branch("self", 0, 0, physical)]).is_err());
    assert!(CarreauPipeGraph::new(2, vec![branch("bad-endpoint", 0, 2, physical)]).is_err());
    assert!(
        CarreauPipeGraph::new(
            2,
            vec![
                branch("same", 0, 1, physical),
                branch("same", 1, 0, physical)
            ]
        )
        .is_err()
    );
    let warmer = LaminarCarreauPipe::new(law(false), 320.0, 0.01, 0.5, 1e4).unwrap();
    assert!(
        CarreauPipeGraph::new(
            3,
            vec![branch("cold", 0, 1, physical), branch("hot", 1, 2, warmer)]
        )
        .is_err()
    );
    let graph = CarreauPipeGraph::new(2, vec![branch("pipe", 0, 1, physical)]).unwrap();
    for boundaries in [
        vec![],
        vec![boundary(2, 0.0)],
        vec![boundary(0, f64::NAN)],
        vec![boundary(0, 1.0), boundary(0, 2.0)],
    ] {
        assert!(with_cx(|cx| graph.solve(&boundaries, controls(), cx)).is_err());
    }
    let malformed = CarreauGraphControls {
        absolute_flow_tolerance: VolumetricFlowRate::new(0.0),
        ..controls()
    };
    assert!(with_cx(|cx| graph.solve(&[boundary(0, 0.0)], malformed, cx)).is_err());
    let mut source_controls = controls();
    source_controls.pipe.maximum_quadrature_refinements = 0;
    assert!(matches!(
        with_cx(|cx| graph.solve(&[boundary(0, 0.0)], source_controls, cx)),
        Err(CarreauGraphError::Pipe { branch: 0, .. })
    ));
}
