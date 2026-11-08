//! Physical interval-mean pressure inverses retain the original trajectories.
use super::*;
use crate::graph::inertial::laminar::{LaminarPipe, LaminarPipeSpec};
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::FluidTransportState;
use fs_phs::ZkDiffusiveBasis;
use fs_qty::Length;

fn pipe(
    cx: &Cx<'_>,
    length: f64,
    mu: f64,
    maximum_reynolds: f64,
) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, mu).unwrap(),
        LaminarPipeSpec {
            length: Length::new(length),
            radius: Length::new(0.001),
            maximum_reynolds,
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

fn tolerance() -> LaminarNetworkFlowTolerance {
    LaminarNetworkFlowTolerance {
        absolute_m3_s: 1e-18,
        relative: 1e-8,
    }
}

fn near(a: f64, b: f64, scale: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * scale,
        "{a:e} != {b:e}; scale={scale:e}"
    );
}

fn same_coordinates(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().chain(b).fold(1e-30_f64, |s, x| s.max(x.abs()));
    for (&a, &b) in a.iter().zip(b) {
        near(a, b, scale, 1e-8);
    }
}

fn balance(report: LaminarNetworkStep, branches: &[LaminarNetworkBranchStep]) {
    let scale = report
        .native
        .supplied
        .abs()
        .max(report.native.dissipated)
        .max(report.initial_kinetic_energy_j)
        .max(report.final_kinetic_energy_j)
        .max(1e-30);
    near(
        report.native.delta_h + report.native.dissipated,
        report.native.supplied,
        scale,
        5e-9,
    );
    near(
        report.reservoir_work_j + report.branch_head_work_j,
        report.native.supplied,
        scale,
        5e-9,
    );
    near(
        branches.iter().map(|b| b.viscous_heat_j).sum(),
        report.native.dissipated,
        scale,
        5e-9,
    );
    near(
        branches
            .iter()
            .map(|b| b.final_mean_kinetic_energy_j + b.final_profile_memory_energy_j)
            .sum(),
        report.final_kinetic_energy_j,
        scale,
        5e-9,
    );
    assert!(branches.iter().all(|b| b.viscous_heat_j >= 0.0));
}

#[test]
fn scalar_inverse_matches_original_pipe_with_nonsteady_radial_memory() {
    with_cx(|cx| {
        let graph = LaminarNetwork::try_new(
            cx,
            topology(2, &[[0, 1]], &[0, 1]),
            vec![pipe(cx, 0.25, 0.001, 1000.0)],
        )
        .unwrap();
        let single = &graph.branches()[0];
        let dt = 0.025;
        let mut native = single.prepare(cx, dt).unwrap();
        let mut original = single
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        for pressure in [80.0, -50.0, 0.0] {
            native
                .advance(
                    cx,
                    &mut original,
                    [Pressure::new(pressure), Pressure::new(0.0)],
                )
                .unwrap();
        }
        let parabolic = single
            .steady_state(
                VolumetricFlowRate::new(original.volume_flow_m3_s()),
                original.accepted_time_s(),
            )
            .unwrap();
        assert!(original.profile_memory_energy_j() > 0.0);
        assert!(
            original.coordinates()[1..]
                .iter()
                .zip(&parabolic.coordinates()[1..])
                .any(|(a, b)| (a - b).abs() > 1e-8 * a.abs().max(b.abs()))
        );
        let mut state = graph
            .state_from_branches(std::slice::from_ref(&original))
            .unwrap();
        let mut controlled = graph
            .prepare_flow_control(cx, dt, &[LaminarNetworkFlowPort::Reservoir(0)])
            .unwrap();
        assert!(controlled.pressure_offsets_pa().is_empty());
        assert!(controlled.controlled_volume_flows_m3_s().is_empty());
        let mut differs_from_endpoint = false;
        for target in [2e-8, -1e-8, 0.0] {
            // Independent scalar inversion: two complete original single-pipe
            // trajectories identify Qbar(P)=Qbar(0)+P*dQbar/dP from this exact
            // initial history. No network inverse/output-map helper is used.
            let mut free_state = original.clone();
            let free = native
                .advance(cx, &mut free_state, [Pressure::new(0.0); 2])
                .unwrap()
                .average_volume_flow_m3_s;
            let mut probe_state = original.clone();
            let probe = native
                .advance(
                    cx,
                    &mut probe_state,
                    [Pressure::new(64.0), Pressure::new(0.0)],
                )
                .unwrap()
                .average_volume_flow_m3_s;
            let required = 64.0 * (target - free) / (probe - free);
            let reference = native
                .advance(
                    cx,
                    &mut original,
                    [Pressure::new(required), Pressure::new(0.0)],
                )
                .unwrap();
            // The reference reservoir is controlled beside a huge common
            // absolute pressure: the solved increment must stay in its gauge.
            let report = controlled
                .advance(
                    cx,
                    &mut state,
                    &[Pressure::new(1e30); 2],
                    &[Pressure::new(0.0)],
                    &[VolumetricFlowRate::new(target)],
                    tolerance(),
                )
                .unwrap();
            near(
                controlled.pressure_offsets_pa()[0],
                required,
                required.abs().max(64.0),
                2e-8,
            );
            near(
                controlled.controlled_volume_flows_m3_s()[0],
                target,
                1e-7,
                1e-8,
            );
            let branch = controlled.branch_steps()[0];
            near(
                branch.average_volume_flow_m3_s,
                reference.average_volume_flow_m3_s,
                1e-7,
                1e-8,
            );
            near(branch.displaced_volume_m3, target * dt, 1e-7 * dt, 1e-8);
            near(
                branch.final_volume_flow_m3_s,
                reference.final_volume_flow_m3_s,
                1e-7,
                1e-8,
            );
            near(
                branch.viscous_heat_j,
                reference.native.dissipated,
                reference.native.dissipated.max(1e-30),
                2e-8,
            );
            let recovered = graph.branch_state(&state, 0).unwrap();
            same_coordinates(recovered.coordinates(), original.coordinates());
            assert_eq!(state.accepted_time_s(), original.accepted_time_s());
            near(
                report.network.reservoir_work_j,
                required * target * dt,
                reference
                    .native
                    .dissipated
                    .max(reference.native.supplied.abs())
                    .max(1e-30),
                2e-8,
            );
            differs_from_endpoint |= (branch.final_volume_flow_m3_s - target).abs() > 1e-10;
            balance(report.network, controlled.branch_steps());
        }
        assert!(
            differs_from_endpoint,
            "the target is a true interval mean, not endpoint flow"
        );
    });
}

#[test]
fn two_port_fork_inverse_recovers_original_pressure_trajectory_and_each_branch_heat() {
    with_cx(|cx| {
        let graph = LaminarNetwork::try_new(
            cx,
            topology(4, &[[0, 1], [1, 2], [1, 3]], &[0, 2, 3]),
            vec![
                pipe(cx, 0.12, 0.001, 1000.0),
                pipe(cx, 0.18, 0.003, 1000.0),
                pipe(cx, 0.27, 0.0015, 1000.0),
            ],
        )
        .unwrap();
        assert_eq!(graph.independent_mean_flows(), 2);
        let mut original = graph.rest_state(0.0).unwrap();
        let mut state = original.clone();
        let mut native = graph.prepare(cx, 0.031).unwrap();
        let mut controlled = graph
            .prepare_flow_control(
                cx,
                0.031,
                &[
                    LaminarNetworkFlowPort::Reservoir(0),
                    LaminarNetworkFlowPort::Reservoir(1),
                ],
            )
            .unwrap();
        let heads = [Pressure::new(3.0), Pressure::new(-2.0), Pressure::new(1.0)];
        for pressures in [[45.0, -5.0, 0.0], [20.0, 15.0, 0.0], [-25.0, 10.0, 0.0]] {
            let reference = native
                .advance_with_heads(cx, &mut original, &pressures.map(Pressure::new), &heads)
                .unwrap();
            let desired = [
                VolumetricFlowRate::new(native.boundary_volume_flows_m3_s()[0]),
                VolumetricFlowRate::new(native.boundary_volume_flows_m3_s()[1]),
            ];
            let report = controlled
                .advance(
                    cx,
                    &mut state,
                    &[Pressure::new(1e30); 3],
                    &heads,
                    &desired,
                    tolerance(),
                )
                .unwrap();
            for port in 0..2 {
                near(
                    controlled.pressure_offsets_pa()[port],
                    pressures[port],
                    64.0,
                    3e-8,
                );
            }
            same_coordinates(state.coordinates(), original.coordinates());
            assert_eq!(state.accepted_time_s(), original.accepted_time_s());
            for (branch, expected) in controlled.branch_steps().iter().zip(native.branch_steps()) {
                near(
                    branch.average_volume_flow_m3_s,
                    expected.average_volume_flow_m3_s,
                    1e-7,
                    1e-8,
                );
                near(
                    branch.final_volume_flow_m3_s,
                    expected.final_volume_flow_m3_s,
                    1e-7,
                    1e-8,
                );
                near(
                    branch.viscous_heat_j,
                    expected.viscous_heat_j,
                    expected.viscous_heat_j.max(1e-30),
                    2e-8,
                );
                near(
                    branch.final_profile_memory_energy_j,
                    expected.final_profile_memory_energy_j,
                    expected.final_profile_memory_energy_j.max(1e-30),
                    2e-8,
                );
            }
            for (&actual, &expected) in controlled
                .boundary_volume_flows_m3_s()
                .iter()
                .zip(native.boundary_volume_flows_m3_s())
            {
                near(actual, expected, 1e-7, 1e-8);
            }
            near(
                report.network.native.supplied,
                reference.native.supplied,
                reference
                    .native
                    .dissipated
                    .max(reference.native.supplied.abs())
                    .max(1e-30),
                2e-8,
            );
            assert!(report.network.maximum_internal_imbalance_m3_s < 1e-17);
            balance(report.network, controlled.branch_steps());
        }
    });
}

#[test]
fn closed_loop_mean_flow_control_recovers_branch_head_work_and_circulation() {
    with_cx(|cx| {
        let graph = LaminarNetwork::try_new(
            cx,
            topology(3, &[[0, 1], [1, 2], [2, 0]], &[0]),
            vec![
                pipe(cx, 0.1, 0.001, 1000.0),
                pipe(cx, 0.2, 0.001, 1000.0),
                pipe(cx, 0.15, 0.002, 1000.0),
            ],
        )
        .unwrap();
        let dt = 0.04;
        let mut original = graph.rest_state(0.0).unwrap();
        let mut state = original.clone();
        let mut native = graph.prepare(cx, dt).unwrap();
        let mut controlled = graph
            .prepare_flow_control(cx, dt, &[LaminarNetworkFlowPort::BranchHead(1)])
            .unwrap();
        let reservoir = [Pressure::new(1e30)];
        let baseline = [Pressure::new(10.0), Pressure::new(-4.0), Pressure::new(3.0)];
        for required_head in [40.0, -35.0, 0.0] {
            let heads = [
                Pressure::new(10.0),
                Pressure::new(required_head),
                Pressure::new(3.0),
            ];
            let reference = native
                .advance_with_heads(cx, &mut original, &reservoir, &heads)
                .unwrap();
            let target = native.branch_steps()[1].average_volume_flow_m3_s;
            let report = controlled
                .advance(
                    cx,
                    &mut state,
                    &reservoir,
                    &baseline,
                    &[VolumetricFlowRate::new(target)],
                    tolerance(),
                )
                .unwrap();
            near(
                controlled.pressure_offsets_pa()[0],
                required_head + 4.0,
                64.0,
                3e-8,
            );
            same_coordinates(state.coordinates(), original.coordinates());
            for branch in controlled.branch_steps() {
                near(branch.average_volume_flow_m3_s, target, 1e-7, 1e-8);
                near(branch.displaced_volume_m3, target * dt, 1e-7 * dt, 1e-8);
            }
            assert_eq!(report.network.reservoir_work_j, 0.0);
            assert!(controlled.boundary_volume_flows_m3_s()[0].abs() < 1e-18);
            near(
                report.network.branch_head_work_j,
                (13.0 + required_head) * target * dt,
                reference
                    .native
                    .dissipated
                    .max(reference.native.supplied.abs())
                    .max(1e-30),
                2e-8,
            );
            balance(report.network, controlled.branch_steps());
        }
    });
}

#[test]
fn dependent_ports_and_failed_flow_proposals_preserve_accepted_history_and_reports() {
    with_cx(|cx| {
        let graph = LaminarNetwork::try_new(
            cx,
            topology(4, &[[0, 1], [1, 2], [1, 3]], &[0, 2, 3]),
            vec![
                pipe(cx, 0.12, 0.001, 0.1),
                pipe(cx, 0.18, 0.003, 1000.0),
                pipe(cx, 0.27, 0.0015, 1000.0),
            ],
        )
        .unwrap();
        let dt = 0.02;
        // These are different authored ports but the SAME physical flow:
        // reservoir 0 has only branch 0. Counting labels cannot establish rank.
        assert!(
            graph
                .prepare_flow_control(
                    cx,
                    dt,
                    &[
                        LaminarNetworkFlowPort::Reservoir(0),
                        LaminarNetworkFlowPort::BranchHead(0)
                    ],
                )
                .is_err()
        );
        assert!(
            graph
                .prepare_flow_control(
                    cx,
                    dt,
                    &[
                        LaminarNetworkFlowPort::Reservoir(0),
                        LaminarNetworkFlowPort::Reservoir(0)
                    ],
                )
                .is_err()
        );
        assert!(
            graph
                .prepare_flow_control(cx, dt, &[LaminarNetworkFlowPort::BranchHead(3)])
                .is_err()
        );
        let ports = [
            LaminarNetworkFlowPort::Reservoir(0),
            LaminarNetworkFlowPort::Reservoir(1),
        ];
        let mut controlled = graph.prepare_flow_control(cx, dt, &ports).unwrap();
        let mut state = graph.rest_state(0.0).unwrap();
        let reservoir = [Pressure::new(0.0); 3];
        let heads = [Pressure::new(0.0); 3];
        let targets = [
            VolumetricFlowRate::new(1e-12),
            VolumetricFlowRate::new(-4e-13),
        ];
        controlled
            .advance(cx, &mut state, &reservoir, &heads, &targets, tolerance())
            .unwrap();
        let before = state.clone();
        let offsets = controlled.pressure_offsets_pa().to_vec();
        let inputs = controlled.gauged_port_pressures_pa().to_vec();
        let flows = controlled.controlled_volume_flows_m3_s().to_vec();
        let boundary = controlled.boundary_volume_flows_m3_s().to_vec();
        let report = format!("{:?}", controlled.last_step());
        let branches = format!("{:?}", controlled.branch_steps());
        assert!(
            controlled
                .advance(
                    cx,
                    &mut state,
                    &reservoir,
                    &heads,
                    &[
                        VolumetricFlowRate::new(1e-6),
                        VolumetricFlowRate::new(-4e-7)
                    ],
                    tolerance(),
                )
                .is_err(),
            "the original whole-interval Reynolds ceiling still gates flow control"
        );
        assert!(
            controlled
                .advance(
                    cx,
                    &mut state,
                    &reservoir,
                    &heads,
                    &[VolumetricFlowRate::new(f64::NAN), targets[1]],
                    tolerance(),
                )
                .is_err()
        );
        assert!(
            controlled
                .advance(cx, &mut state, &reservoir, &[], &targets, tolerance())
                .is_err()
        );
        let gate = fs_exec::CancelGate::new();
        gate.request();
        let cancelled = Cx::new(&gate, cx.arena(), cx.stream_key(), cx.budget(), cx.mode());
        assert!(
            controlled
                .advance(
                    &cancelled,
                    &mut state,
                    &reservoir,
                    &heads,
                    &targets,
                    tolerance()
                )
                .is_err()
        );
        assert_eq!(state.coordinates(), before.coordinates());
        assert_eq!(state.accepted_time_s(), before.accepted_time_s());
        assert_eq!(controlled.pressure_offsets_pa(), offsets);
        assert_eq!(controlled.gauged_port_pressures_pa(), inputs);
        assert_eq!(controlled.controlled_volume_flows_m3_s(), flows);
        assert_eq!(controlled.boundary_volume_flows_m3_s(), boundary);
        assert_eq!(format!("{:?}", controlled.last_step()), report);
        assert_eq!(format!("{:?}", controlled.branch_steps()), branches);

        // Failed previews may alter scratch, but not the next accepted solve.
        let mut fresh = graph.prepare_flow_control(cx, dt, &ports).unwrap();
        let mut clean_state = before;
        fresh
            .advance(
                cx,
                &mut clean_state,
                &reservoir,
                &heads,
                &targets,
                tolerance(),
            )
            .unwrap();
        controlled
            .advance(cx, &mut state, &reservoir, &heads, &targets, tolerance())
            .unwrap();
        assert_eq!(state.coordinates(), clean_state.coordinates());
        assert_eq!(
            controlled.pressure_offsets_pa(),
            fresh.pressure_offsets_pa()
        );
    });
}
