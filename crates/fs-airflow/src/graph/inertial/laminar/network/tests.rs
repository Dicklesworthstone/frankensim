//! G1/G3 network limits: native pipe/series trajectories, Kirchhoff steady
//! flow, closed-loop source work, gauge invariance and original state transfer.
use super::*;
use crate::graph::inertial::laminar::{LaminarPipeSpec, series::LaminarSeries};
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::FluidTransportState;
use fs_phs::ZkDiffusiveBasis;
use fs_qty::{Length, VolumetricFlowRate};

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, mu).unwrap()
}
fn pipe(cx: &Cx<'_>, length: f64, mu: f64) -> LaminarPipe<FluidTransportState> {
    pipe_limit(cx, length, mu, 1000.0)
}
fn pipe_limit(cx: &Cx<'_>, length: f64, mu: f64, re: f64) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        fluid(mu),
        LaminarPipeSpec {
            length: Length::new(length),
            radius: Length::new(0.001),
            maximum_reynolds: re,
            band_rad_s: [0.01, 900.0],
            maximum_relative_wall_error: 1e-4,
        },
        &ZkDiffusiveBasis::new(16).unwrap(),
    )
    .unwrap()
}
fn topology(
    nodes: usize,
    branches: &[[usize; 2]],
    pressure_nodes: &[usize],
) -> LaminarNetworkTopology {
    LaminarNetworkTopology {
        node_count: nodes,
        branch_nodes: branches.to_vec(),
        pressure_nodes: pressure_nodes.to_vec(),
    }
}
fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scale {scale:e}"
    );
}
fn balance(report: LaminarNetworkStep, branches: &[LaminarNetworkBranchStep]) {
    let scale = report
        .native
        .supplied
        .abs()
        .max(report.native.dissipated)
        .max(report.initial_kinetic_energy_j)
        .max(report.final_kinetic_energy_j);
    near(
        report.native.delta_h + report.native.dissipated,
        report.native.supplied,
        scale,
        4e-10,
    );
    near(
        report.reservoir_work_j + report.branch_head_work_j,
        report.native.supplied,
        scale,
        4e-10,
    );
    near(
        branches.iter().map(|b| b.viscous_heat_j).sum(),
        report.native.dissipated,
        scale,
        4e-10,
    );
    near(
        branches
            .iter()
            .map(|b| b.final_mean_kinetic_energy_j + b.final_profile_memory_energy_j)
            .sum(),
        report.final_kinetic_energy_j,
        scale,
        4e-12,
    );
    assert!(branches.iter().all(|b| b.viscous_heat_j >= 0.0));
}

#[test]
fn single_branch_keeps_native_startup_reversal_and_heat() {
    with_cx(|cx| {
        let single = pipe(cx, 0.25, 0.001);
        let graph = LaminarNetwork::try_new(
            cx,
            topology(2, &[[0, 1]], &[0, 1]),
            vec![pipe(cx, 0.25, 0.001)],
        )
        .unwrap();
        let mut original = single
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        let mut state = graph.rest_state(0.0).unwrap();
        let mut one = single.prepare(cx, 0.025).unwrap();
        let mut prepared = graph.prepare(cx, 0.025).unwrap();
        for p in [80.0, 80.0, -50.0, -50.0, 0.0] {
            let pressure = [Pressure::new(100_000.0 + p), Pressure::new(100_000.0)];
            let reference = one.advance(cx, &mut original, pressure).unwrap();
            let report = prepared.advance(cx, &mut state, &pressure).unwrap();
            let branch = prepared.branch_steps()[0];
            near(
                branch.final_volume_flow_m3_s,
                reference.final_volume_flow_m3_s,
                1e-7,
                2e-9,
            );
            near(
                branch.average_volume_flow_m3_s,
                reference.average_volume_flow_m3_s,
                1e-7,
                2e-9,
            );
            near(
                branch.viscous_heat_j,
                reference.native.dissipated,
                reference.native.dissipated,
                5e-9,
            );
            near(
                state.kinetic_energy_j(),
                original.kinetic_energy_j(),
                original.kinetic_energy_j(),
                5e-9,
            );
            balance(report, prepared.branch_steps());
        }
    });
}

#[test]
fn internal_junction_matches_original_coupled_series_through_reversal() {
    with_cx(|cx| {
        let series =
            LaminarSeries::try_new(cx, vec![pipe(cx, 0.13, 0.001), pipe(cx, 0.27, 0.004)]).unwrap();
        let graph = LaminarNetwork::try_new(
            cx,
            topology(3, &[[0, 1], [1, 2]], &[0, 2]),
            vec![pipe(cx, 0.13, 0.001), pipe(cx, 0.27, 0.004)],
        )
        .unwrap();
        assert_eq!(graph.independent_mean_flows(), 1);
        let mut original = series
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        let mut state = graph.rest_state(0.0).unwrap();
        let mut serial = series.prepare(cx, 0.021).unwrap();
        let mut prepared = graph.prepare(cx, 0.021).unwrap();
        let mut heat = [0.0; 2];
        for p in [70.0, 70.0, -40.0, 0.0, 15.0] {
            let pressure = [Pressure::new(p), Pressure::new(0.0)];
            let reference = serial
                .advance(cx, &mut original, pressure, &mut heat)
                .unwrap();
            let report = prepared.advance(cx, &mut state, &pressure).unwrap();
            for (edge, branch) in prepared.branch_steps().iter().enumerate() {
                near(
                    branch.final_volume_flow_m3_s,
                    original.volume_flow_m3_s(),
                    1e-7,
                    3e-9,
                );
                near(
                    branch.average_volume_flow_m3_s,
                    reference.average_volume_flow_m3_s,
                    1e-7,
                    3e-9,
                );
                near(branch.viscous_heat_j, heat[edge], heat[edge], 6e-9);
            }
            assert!(report.maximum_internal_imbalance_m3_s < 1e-18);
            balance(report, prepared.branch_steps());
        }
    });
}

#[test]
fn mesh_steady_solution_matches_independent_poiseuille_kirchhoff_equations() {
    with_cx(|cx| {
        let lengths = [0.11, 0.25, 0.19, 0.13, 0.09];
        let viscosities = [0.001, 0.002, 0.0015, 0.0008, 0.003];
        let edges = [[0, 1], [1, 3], [0, 2], [2, 3], [1, 2]];
        let graph = LaminarNetwork::try_new(
            cx,
            topology(4, &edges, &[0, 3]),
            lengths
                .iter()
                .zip(viscosities)
                .map(|(&l, mu)| pipe(cx, l, mu))
                .collect(),
        )
        .unwrap();
        assert_eq!(graph.independent_mean_flows(), 3); // two through paths plus one cycle
        let resistance: Vec<_> = lengths
            .iter()
            .zip(viscosities)
            .map(|(&l, mu)| 8.0 * mu * l / (std::f64::consts::PI * 0.001_f64.powi(4)))
            .collect();
        let conductance: Vec<_> = resistance.iter().map(|r| 1.0 / r).collect();
        let a = conductance[0] + conductance[1] + conductance[4];
        let d = conductance[2] + conductance[3] + conductance[4];
        let c = -conductance[4];
        let determinant = a * d - c * c;
        let rhs1 = 60.0 * conductance[0];
        let rhs2 = 60.0 * conductance[2];
        let node = [
            60.0,
            (rhs1 * d - c * rhs2) / determinant,
            (a * rhs2 - c * rhs1) / determinant,
            0.0,
        ];
        let mut state = graph
            .steady_state(&[Pressure::new(60.0), Pressure::new(0.0)], 0.0)
            .unwrap();
        for (edge, &[from, to]) in edges.iter().enumerate() {
            let exact = (node[from] - node[to]) / resistance[edge];
            near(
                graph.branch_volume_flow_m3_s(&state, edge).unwrap(),
                exact,
                1e-7,
                5e-10,
            );
        }
        // A large interval at an already steady laminar state must not acquire
        // the unbounded dt*forcing Reynolds envelope of a zero-origin estimate.
        let mut prepared = graph.prepare(cx, 2.0).unwrap();
        let before = state.clone();
        let report = prepared
            .advance(cx, &mut state, &[Pressure::new(60.0), Pressure::new(0.0)])
            .unwrap();
        for edge in 0..edges.len() {
            near(
                graph.branch_volume_flow_m3_s(&state, edge).unwrap(),
                graph.branch_volume_flow_m3_s(&before, edge).unwrap(),
                1e-7,
                5e-9,
            );
        }
        balance(report, prepared.branch_steps());
    });
}

#[test]
fn closed_loop_pressure_head_drives_circulation_and_pays_the_original_heat() {
    with_cx(|cx| {
        let graph = LaminarNetwork::try_new(
            cx,
            topology(3, &[[0, 1], [1, 2], [2, 0]], &[0]),
            vec![
                pipe(cx, 0.1, 0.001),
                pipe(cx, 0.2, 0.001),
                pipe(cx, 0.15, 0.002),
            ],
        )
        .unwrap();
        assert_eq!(graph.independent_mean_flows(), 1);
        let pressure = [Pressure::new(1e9)];
        let heads = [Pressure::new(120.0), Pressure::new(0.0), Pressure::new(0.0)];
        let exact = 120.0
            / graph
                .branches()
                .iter()
                .map(|p| p.steady_resistance_pa_s_m3())
                .sum::<f64>();
        let mut state = graph
            .steady_state_with_heads(&pressure, &heads, 0.0)
            .unwrap();
        let mut prepared = graph.prepare(cx, 0.8).unwrap();
        let report = prepared
            .advance_with_heads(cx, &mut state, &pressure, &heads)
            .unwrap();
        assert_eq!(report.reservoir_work_j, 0.0);
        assert!(report.branch_head_work_j > 0.0);
        for branch in prepared.branch_steps() {
            near(branch.average_volume_flow_m3_s, exact, exact, 5e-9);
            near(branch.final_volume_flow_m3_s, exact, exact, 5e-9);
        }
        near(
            report.branch_head_work_j,
            120.0 * exact * 0.8,
            120.0 * exact * 0.8,
            5e-9,
        );
        assert!(prepared.boundary_volume_flows_m3_s()[0].abs() < 1e-18);
        balance(report, prepared.branch_steps());
    });
}

#[test]
fn separate_component_gauges_and_material_continuation_preserve_original_histories() {
    with_cx(|cx| {
        let graph = LaminarNetwork::try_new(
            cx,
            topology(4, &[[0, 1], [2, 3]], &[0, 1, 2, 3]),
            vec![pipe(cx, 0.2, 0.001), pipe(cx, 0.3, 0.003)],
        )
        .unwrap();
        let mut zero = graph.rest_state(0.0).unwrap();
        let mut prepared = graph.prepare(cx, 0.02).unwrap();
        let equal = [
            Pressure::new(1e12),
            Pressure::new(1e12),
            Pressure::new(-1e11),
            Pressure::new(-1e11),
        ];
        let rest = prepared.advance(cx, &mut zero, &equal).unwrap();
        assert_eq!(rest.native.supplied, 0.0);
        assert_eq!(zero.kinetic_energy_j(), 0.0);
        let mut a = graph.rest_state(0.0).unwrap();
        let mut b = a.clone();
        prepared
            .advance(
                cx,
                &mut a,
                &[
                    Pressure::new(80.0),
                    Pressure::new(0.0),
                    Pressure::new(25.0),
                    Pressure::new(0.0),
                ],
            )
            .unwrap();
        prepared
            .advance(
                cx,
                &mut b,
                &[
                    Pressure::new(1e12 + 80.0),
                    Pressure::new(1e12),
                    Pressure::new(-1e11 + 25.0),
                    Pressure::new(-1e11),
                ],
            )
            .unwrap();
        assert_eq!(a.coordinates(), b.coordinates());
        let originals: Vec<_> = (0..2).map(|e| graph.branch_state(&a, e).unwrap()).collect();
        let imported = graph.state_from_branches(&originals).unwrap();
        near(
            imported.kinetic_energy_j(),
            a.kinetic_energy_j(),
            a.kinetic_energy_j(),
            1e-14,
        );
        let updated = graph
            .with_sources(cx, vec![fluid(0.004), fluid(0.0007)])
            .unwrap();
        let before = a.clone();
        let mut next = updated.prepare(cx, 0.03).unwrap();
        next.advance(cx, &mut a, &[Pressure::new(0.0); 4]).unwrap();
        assert!(a.kinetic_energy_j() < before.kinetic_energy_j());
        assert_eq!(before.accepted_time_s(), 0.02);
    });
}

#[test]
fn incompatible_import_invalid_heads_regime_and_cancellation_refuse_atomically() {
    with_cx(|cx| {
        let graph = LaminarNetwork::try_new(
            cx,
            topology(3, &[[0, 1], [1, 2]], &[0, 2]),
            vec![pipe(cx, 0.1, 0.001), pipe_limit(cx, 0.2, 0.001, 0.1)],
        )
        .unwrap();
        let bad = graph
            .branches()
            .iter()
            .enumerate()
            .map(|(e, p)| {
                p.steady_state(
                    VolumetricFlowRate::new(if e == 0 { 1e-11 } else { 2e-11 }),
                    0.0,
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        assert!(graph.state_from_branches(&bad).is_err());
        let mut state = graph.rest_state(0.0).unwrap();
        let mut prepared = graph.prepare(cx, 0.02).unwrap();
        prepared
            .advance(cx, &mut state, &[Pressure::new(1e-4), Pressure::new(0.0)])
            .unwrap();
        let before = state.clone();
        let heat: Vec<_> = prepared
            .branch_steps()
            .iter()
            .map(|p| p.viscous_heat_j)
            .collect();
        let boundary = prepared.boundary_volume_flows_m3_s().to_vec();
        assert!(
            prepared
                .advance(cx, &mut state, &[Pressure::new(10.0), Pressure::new(0.0)])
                .is_err()
        );
        assert!(
            prepared
                .advance_with_heads(
                    cx,
                    &mut state,
                    &[Pressure::new(0.0); 2],
                    &[Pressure::new(0.0)]
                )
                .is_err()
        );
        assert!(
            prepared
                .advance(
                    cx,
                    &mut state,
                    &[Pressure::new(f64::NAN), Pressure::new(0.0)]
                )
                .is_err()
        );
        let gate = fs_exec::CancelGate::new();
        gate.request();
        let cancelled = Cx::new(&gate, cx.arena(), cx.stream_key(), cx.budget(), cx.mode());
        assert!(
            prepared
                .advance(&cancelled, &mut state, &[Pressure::new(0.0); 2])
                .is_err()
        );
        assert_eq!(state.coordinates(), before.coordinates());
        assert_eq!(state.accepted_time_s(), before.accepted_time_s());
        assert_eq!(prepared.boundary_volume_flows_m3_s(), boundary);
        assert_eq!(
            prepared
                .branch_steps()
                .iter()
                .map(|p| p.viscous_heat_j)
                .collect::<Vec<_>>(),
            heat
        );
    });
}
