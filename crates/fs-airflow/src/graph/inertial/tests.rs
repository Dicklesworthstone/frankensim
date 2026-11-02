//! Native branch-momentum stepping against independent small-network solutions.

use super::*;
use crate::graph::gravity::GraphGravitySection;
use crate::graph::porous::{LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::graph::tests::{boundary, edge, with_cx};
use crate::{SourceProvenance, ToleranceBasis};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};

fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101325.0).unwrap()
}
fn parallel(resistances: [f64; 3]) -> LossGraph {
    LossGraph::new_mixed(
        2,
        resistances
            .into_iter()
            .enumerate()
            .map(|(i, resistance)| PassiveBranch::Linear {
                from: 0,
                to: 1,
                loss: LinearLossElement::new(
                    format!("path-{i}"),
                    LinearLossResistance::new(resistance),
                    0.0,
                    SourceProvenance::new("authored linear network fixture", "inertial-graph-v1"),
                    ToleranceBasis::EngineeringAllowance,
                )
                .unwrap(),
            })
            .collect(),
    )
    .unwrap()
}
fn geometry(n: usize) -> Vec<InertialBranch> {
    vec![
        InertialBranch {
            length_m: 1.0,
            flow_area_m2: 1.0
        };
        n
    ]
}
fn controls() -> InertialGraphControls {
    InertialGraphControls {
        maximum_drag_step: 2.0,
        maximum_neglected_power_w: 1e6,
        maximum_neglected_energy_j: 1e6,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}
fn owner(
    cx: &Cx<'_>,
    graph: &LossGraph,
    flows: &[f64],
    boundaries: &[FixedPressure],
) -> InertialGraph {
    InertialGraph::try_new(
        cx,
        graph,
        Density::new(gas(300.0).density),
        &geometry(flows.len()),
        boundaries,
        flows,
        0.0,
        16,
    )
    .unwrap()
}

fn bent_loop(
    cx: &Cx<'_>,
    rising_density: f64,
    falling_density: f64,
    height_m: f64,
    maximum_sections: usize,
) -> LossGraph {
    let transport = |rho| FluidTransportState::try_new(rho, 1000.0, 0.1, 0.001).unwrap();
    let section = |start_m, end_m, rho| GraphGravitySection {
        start_m,
        end_m,
        transport: transport(rho),
    };
    let graph = LossGraph::new_mixed(
        2,
        [(0, 1), (1, 0)]
            .into_iter()
            .enumerate()
            .map(|(i, (from, to))| PassiveBranch::Linear {
                from,
                to,
                loss: LinearLossElement::new(
                    format!("bent-{i}"),
                    LinearLossResistance::new(1.0),
                    0.0,
                    SourceProvenance::new("analytic bent pipe", "integrated-gravity-v1"),
                    ToleranceBasis::Analytic,
                )
                .unwrap(),
            })
            .collect(),
    )
    .unwrap();
    graph
        .with_distributed_boussinesq_gravity(
            &transport(1.0),
            &[
                vec![
                    section([0.0; 3], [0.0, 0.0, height_m], rising_density),
                    section([0.0, 0.0, height_m], [1.0, 0.0, height_m], 1.0),
                    section([1.0, 0.0, height_m], [1.0, 0.0, 0.0], falling_density),
                ],
                vec![section([1.0, 0.0, 0.0], [0.0; 3], 1.0)],
            ],
            &[[0.0; 3], [1.0, 0.0, 0.0]],
            [0.0, 0.0, -8.0],
            0.3,
            maximum_sections,
            cx,
        )
        .unwrap()
}

#[test]
fn g1_bent_pipe_integral_accelerates_native_momentum_and_freezes_the_full_path() {
    with_cx(|cx| {
        let graph = bent_loop(cx, 0.75, 1.0, 2.0, 4);
        let boundaries = [boundary(0, 0.0)];
        // Actual paths have lengths 5 and 1. Areas 5 and 1 give L=1 each.
        let inertia = [
            InertialBranch {
                length_m: 5.0,
                flow_area_m2: 5.0,
            },
            InertialBranch {
                length_m: 1.0,
                flow_area_m2: 1.0,
            },
        ];
        let initial = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(1.0),
            &inertia,
            &boundaries,
            &[0.0; 2],
            0.0,
            16,
        )
        .unwrap();
        // Sum of native midpoint equations: 2L*dQ/dt + 2R*Q_average = 4 Pa.
        let average = 2.0 / 3.0;
        let (next, report) = initial
            .propose(
                cx,
                &graph,
                &boundaries,
                &[0.0; 2],
                &[average; 2],
                1.0,
                controls(),
            )
            .unwrap()
            .accept(cx, 1e-10, 1e-10)
            .unwrap();
        for flow in &report.average_volume_flows_m3_s {
            near(*flow, average, 1e-12);
        }
        for flow in next.volume_flows_m3_s() {
            near(flow, 2.0 * average, 1e-12);
        }
        assert_eq!(report.driving_pressure_pa, vec![0.0; 2]);
        assert_eq!(report.driving_work_j, 0.0);
        near(report.anomalous_gravity_work_j, 4.0 * average, 1e-12);
        near(
            next.kinetic_energy_j() + report.hydraulic_loss_j,
            report.anomalous_gravity_work_j,
            1e-12,
        );
        near(
            report.neglected_buoyancy_energy_bound_j,
            8.0 * average,
            1e-12,
        );

        let refreshed = bent_loop(cx, 0.875, 1.0, 2.0, 4);
        assert!(
            graph
                .gravity()
                .unwrap()
                .same_geometry_as(refreshed.gravity().unwrap())
        );
        let fresh = initial
            .propose(
                cx,
                &refreshed,
                &boundaries,
                &[0.0; 2],
                &[1.0 / 3.0; 2],
                1.0,
                controls(),
            )
            .unwrap();
        near(
            fresh.report().average_volume_flows_m3_s[0],
            1.0 / 3.0,
            1e-12,
        );
        for changed in [
            bent_loop(cx, 0.75, 1.0, 3.0, 4),
            bent_loop(cx, 0.75, 1.0, 2.0, 5),
        ] {
            assert!(matches!(
                initial.propose(
                    cx,
                    &changed,
                    &boundaries,
                    &[0.0; 2],
                    &[average; 2],
                    1.0,
                    controls()
                ),
                Err(InertialGraphError::Invalid(_))
            ));
        }
        assert_eq!(initial.accepted_time_s(), 0.0);
        assert_eq!(initial.volume_flows_m3_s(), vec![0.0; 2]);
    });
}

#[test]
fn g0_cancelling_bent_pipe_heads_cannot_hide_local_buoyancy_work() {
    with_cx(|cx| {
        let graph = bent_loop(cx, 0.75, 0.75, 2.0, 4);
        assert_eq!(graph.gravity().unwrap().departure_heads_pa(), &[0.0; 2]);
        assert_eq!(
            graph.gravity().unwrap().absolute_departure_heads_pa(),
            &[8.0, 0.0]
        );
        let boundaries = [boundary(0, 0.0)];
        let state = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(1.0),
            &geometry(2),
            &boundaries,
            &[1.0; 2],
            0.0,
            16,
        )
        .unwrap();
        let trial = state
            .propose(
                cx,
                &graph,
                &boundaries,
                &[1.0; 2],
                &[1.0; 2],
                1.0,
                controls(),
            )
            .unwrap();
        assert_eq!(trial.report().anomalous_gravity_work_j, 0.0);
        near(trial.report().neglected_buoyancy_energy_bound_j, 8.0, 1e-12);
        near(trial.report().neglected_energy_bound_j, 10.0, 1e-12);
        let tight = InertialGraphControls {
            maximum_neglected_power_w: 9.0,
            ..controls()
        };
        assert!(matches!(
            state.propose(cx, &graph, &boundaries, &[1.0; 2], &[1.0; 2], 1.0, tight),
            Err(InertialGraphError::NeglectedWork { .. })
        ));
        assert!(matches!(
            trial.recheck(cx, &graph, &boundaries, &[1.0; 2], tight),
            Err(InertialGraphError::NeglectedWork { .. })
        ));
        assert_eq!(state.accepted_time_s(), 0.0);
        assert_eq!(state.volume_flows_m3_s(), vec![1.0; 2]);
    });
}

#[test]
fn g0_constitutive_recheck_audits_the_same_trajectory_without_restepping() {
    with_cx(|cx| {
        let network = |resistance| {
            LossGraph::new_mixed(
                2,
                [(0, 1), (1, 0)]
                    .into_iter()
                    .enumerate()
                    .map(|(i, (from, to))| PassiveBranch::Linear {
                        from,
                        to,
                        loss: LinearLossElement::new(
                            format!("path-{i}"),
                            LinearLossResistance::new(resistance),
                            0.0,
                            SourceProvenance::new(
                                "authored linear network fixture",
                                "inertial-graph-v1",
                            ),
                            ToleranceBasis::EngineeringAllowance,
                        )
                        .unwrap(),
                    })
                    .collect(),
            )
            .unwrap()
        };
        let graph = network(1.0);
        let boundaries = [boundary(0, 0.0)];
        let initial = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(1.0),
            &geometry(2),
            &boundaries,
            &[0.0; 2],
            0.0,
            16,
        )
        .unwrap();
        let proposal = initial
            .propose(
                cx,
                &graph,
                &boundaries,
                &[3.0; 2],
                &[1.0; 2],
                1.0,
                controls(),
            )
            .unwrap();
        let old_x = proposal.report().momentum.x.clone();
        let old_y = proposal.report().momentum.y.clone();
        let refreshed = network(1.5);
        let checked = proposal
            .clone()
            .recheck(cx, &refreshed, &boundaries, &[3.0; 2], controls())
            .unwrap();
        assert_eq!(checked.report().momentum.x, old_x);
        assert_eq!(checked.report().momentum.y, old_y);
        near(checked.report().maximum_impulse_residual_pa_s, 0.5, 1e-12);
        near(checked.report().mechanical_balance_residual_j, 1.0, 1e-12);
        assert!(matches!(
            checked.accept(cx, 1e-8, 1e-8),
            Err(InertialGraphError::MomentumResidual { .. })
        ));

        // A compensating source makes the same trajectory physically valid.
        // The refreshed original-law ledger changes, not the native solve receipt.
        let (next, report) = proposal
            .recheck(cx, &refreshed, &boundaries, &[3.5; 2], controls())
            .unwrap()
            .accept(cx, 1e-12, 1e-12)
            .unwrap();
        assert_eq!(report.momentum.x, old_x);
        assert_eq!(report.momentum.y, old_y);
        near(report.hydraulic_loss_j, 3.0, 1e-12);
        near(report.driving_work_j, 7.0, 1e-12);
        near(report.momentum.dissipated, 2.0, 1e-12);
        assert_eq!(next.volume_flows_m3_s(), vec![2.0; 2]);
        assert_eq!(initial.volume_flows_m3_s(), vec![0.0; 2]);
        assert_eq!(initial.accepted_time_s(), 0.0);
    });
}

#[test]
fn g1_branched_start_stop_and_reversal_retain_two_independent_momenta() {
    with_cx(|cx| {
        // Three paths share two junctions: two independent circulations, not a
        // series loop. Equal L/R allows a separate closed-form projection.
        let graph = parallel([2.0; 3]);
        let boundaries = [boundary(0, 0.0)];
        let initial = owner(cx, &graph, &[0.0; 3], &boundaries);
        let mut state = initial.clone();
        let l = state.inertances_pa_s2_m3()[0];
        let dt = 0.125;
        let mut witnessed_reversal = false;
        for stage in 0..4 {
            let before = state.volume_flows_m3_s();
            let drive = match stage {
                0 => vec![3.0, 0.0, 0.0],
                // Exact desired Q_new=-Q_old/2, Q_average=Q_old/4.
                1 => before.iter().map(|q| -1.5 * l * q / dt + 0.5 * q).collect(),
                2 => vec![0.0, 2.0, -2.0],
                _ => vec![0.0; 3],
            };
            let mean = drive.iter().sum::<f64>() / 3.0;
            let expected: Vec<_> = before
                .iter()
                .zip(&drive)
                .map(|(q, p)| (q + 0.5 * dt * (p - mean) / l) / (1.0 + dt / l))
                .collect();
            let next_time = state.accepted_time_s() + dt;
            let old_energy = state.kinetic_energy_j();
            let trial = state
                .propose(
                    cx,
                    &graph,
                    &boundaries,
                    &drive,
                    &expected,
                    next_time,
                    controls(),
                )
                .unwrap();
            let (next, report) = trial.accept(cx, 1e-11, 1e-11).unwrap();
            for i in 0..3 {
                near(report.average_volume_flows_m3_s[i], expected[i], 2e-12);
                near(
                    report.final_flows_m3_s[i],
                    2.0 * expected[i] - before[i],
                    2e-12,
                );
            }
            near(report.average_pressures_pa[1], mean, 2e-12);
            near(report.average_volume_flows_m3_s.iter().sum(), 0.0, 2e-13);
            near(report.final_flows_m3_s.iter().sum(), 0.0, 2e-13);
            near(report.momentum.balance_residual(), 0.0, 1e-11);
            near(report.mechanical_balance_residual_j, 0.0, 1e-11);
            if stage == 1 {
                assert!(report.average_volume_flows_m3_s[0] > 0.0);
                assert!(report.final_flows_m3_s[0] < 0.0);
                witnessed_reversal = true;
            }
            if stage == 2 {
                assert!((report.final_flows_m3_s[1] - report.final_flows_m3_s[2]).abs() > 0.1);
            }
            if stage == 3 {
                assert!(next.kinetic_energy_j() < old_energy);
                near(
                    next.kinetic_energy_j() - old_energy,
                    -report.hydraulic_loss_j,
                    1e-11,
                );
            }
            near(state.accepted_time_s(), next_time - dt, 0.0);
            state = next;
        }
        assert!(witnessed_reversal);
        assert_eq!(initial.accepted_time_s(), 0.0);
        assert_eq!(initial.volume_flows_m3_s(), vec![0.0; 3]);
    });
}

#[test]
fn g1_original_quadratic_momentum_refuses_a_lagged_linearized_step() {
    with_cx(|cx| {
        let graph = LossGraph::new(
            2,
            (0..3)
                .map(|i| edge(&format!("path-{i}"), 0, 1, 10.0))
                .collect(),
        )
        .unwrap();
        let boundaries = [boundary(0, 0.0)];
        let state = owner(cx, &graph, &[0.0; 3], &boundaries);
        let dt = 0.125;
        let l = state.inertances_pa_s2_m3()[0];
        // Continuity gives [q,-q/2,-q/2]. Subtract the first return equation
        // from the driven equation: 3 L q/dt + 1.25 B q^2 = drive.
        let (mut low, mut high) = (0.0, 3.0 * dt / (3.0 * l));
        for _ in 0..80 {
            let q = 0.5 * (low + high);
            if 3.0 * l * q / dt + 12.5 * q * q > 3.0 {
                high = q;
            } else {
                low = q;
            }
        }
        let q = 0.5 * (low + high);
        let average = [q, -0.5 * q, -0.5 * q];
        let correct = state
            .propose(
                cx,
                &graph,
                &boundaries,
                &[3.0, 0.0, 0.0],
                &average,
                dt,
                controls(),
            )
            .unwrap();
        let (next, report) = correct.accept(cx, 1e-11, 1e-11).unwrap();
        for i in 0..3 {
            near(report.average_volume_flows_m3_s[i], average[i], 2e-12);
        }
        near(next.kinetic_energy_j(), 3.0 * l * q * q, 2e-12);
        let lagged = state
            .propose(
                cx,
                &graph,
                &boundaries,
                &[3.0, 0.0, 0.0],
                &[0.0; 3],
                dt,
                controls(),
            )
            .unwrap();
        assert!(matches!(
            lagged.accept(cx, 1e-11, 1e-11),
            Err(InertialGraphError::MomentumResidual { .. })
        ));
        assert_eq!(state.accepted_time_s(), 0.0);
    });
}

#[test]
fn g1_material_buoyancy_and_viscosity_drive_a_branched_network_from_rest() {
    with_cx(|cx| {
        let reference = gas(300.0);
        let positions = [[0.0; 3], [0.0, 0.0, 1.0]];
        let gravity = [0.0, 0.0, -9.81];
        let cold = parallel([2.0; 3])
            .with_boussinesq_gravity(&reference, &[reference; 3], &positions, gravity, 0.25, cx)
            .unwrap();
        let boundaries = [boundary(0, 0.0)];
        let state = owner(cx, &cold, &[0.0; 3], &boundaries);
        let dt = 0.125;
        let rest = state
            .propose(cx, &cold, &boundaries, &[0.0; 3], &[0.0; 3], dt, controls())
            .unwrap();
        let (_, rest_report) = rest.accept(cx, 0.0, 0.0).unwrap();
        assert_eq!(rest_report.average_volume_flows_m3_s, vec![0.0; 3]);
        near(
            rest_report.average_pressures_pa[1],
            -9.81 * reference.density,
            1e-13,
        );
        let sources = [gas(360.0), gas(300.0), gas(320.0)];
        let resistances =
            sources.map(|source| 2.0 * source.dynamic_viscosity / reference.dynamic_viscosity);
        let heated = parallel(resistances)
            .with_boussinesq_gravity(&reference, &sources, &positions, gravity, 0.25, cx)
            .unwrap();
        // Independent diagonal mobility plus the single continuity multiplier.
        let heads = sources.map(|source| (source.density - reference.density) * gravity[2]);
        let l = state.inertances_pa_s2_m3()[0];
        let mobility = resistances.map(|r| 1.0 / (2.0 * l / dt + r));
        let pressure = mobility.iter().zip(heads).map(|(a, h)| a * h).sum::<f64>()
            / mobility.iter().sum::<f64>();
        let average: Vec<_> = mobility
            .iter()
            .zip(heads)
            .map(|(a, h)| a * (h - pressure))
            .collect();
        let (next, report) = state
            .propose(
                cx,
                &heated,
                &boundaries,
                &[0.0; 3],
                &average,
                dt,
                controls(),
            )
            .unwrap()
            .accept(cx, 1e-11, 1e-11)
            .unwrap();
        for i in 0..3 {
            near(report.average_volume_flows_m3_s[i], average[i], 2e-12);
        }
        near(
            report.average_pressures_pa[1],
            -9.81 * reference.density + pressure,
            2e-12,
        );
        assert!(next.kinetic_energy_j() > 0.0);
        assert!(report.anomalous_gravity_work_j > 0.0);
        near(
            report.momentum.delta_h + report.hydraulic_loss_j,
            report.anomalous_gravity_work_j,
            1e-11,
        );
        assert!(report.neglected_buoyancy_energy_bound_j >= report.anomalous_gravity_work_j.abs());
    });
}

#[test]
fn g1_pressure_ports_supply_the_same_native_inertial_work() {
    with_cx(|cx| {
        let graph = LossGraph::new_mixed(
            3,
            (0..2)
                .map(|i| PassiveBranch::Linear {
                    from: i,
                    to: i + 1,
                    loss: LinearLossElement::new(
                        format!("path-{i}"),
                        LinearLossResistance::new(2.0),
                        0.0,
                        SourceProvenance::new("authored pressure path fixture", "inertial-open-v1"),
                        ToleranceBasis::EngineeringAllowance,
                    )
                    .unwrap(),
                })
                .collect(),
        )
        .unwrap();
        let boundaries = [boundary(0, 3.0), boundary(2, 0.0)];
        let state = owner(cx, &graph, &[0.0; 2], &boundaries);
        let dt = 0.125;
        let total_l = state.inertances_pa_s2_m3().iter().sum::<f64>();
        let q = 3.0 / (2.0 * total_l / dt + 4.0);
        let (next, report) = state
            .propose(cx, &graph, &boundaries, &[0.0; 2], &[q; 2], dt, controls())
            .unwrap()
            .accept(cx, 1e-11, 1e-11)
            .unwrap();
        near(report.average_pressures_pa[1], 1.5, 1e-12);
        for flow in report.average_volume_flows_m3_s {
            near(flow, q, 2e-12);
        }
        near(report.node_outflows_m3_s[0], q, 2e-12);
        near(report.node_outflows_m3_s[2], -q, 2e-12);
        near(report.reduced_boundary_work_j, dt * 3.0 * q, 2e-12);
        near(
            next.kinetic_energy_j() + report.hydraulic_loss_j,
            dt * 3.0 * q,
            1e-11,
        );
        assert_eq!(report.driving_work_j, 0.0);
        assert_eq!(report.anomalous_gravity_work_j, 0.0);
    });
}

#[test]
fn g0_model_and_continuity_refusals_preserve_retained_state_and_nonzero_speed_envelopes() {
    with_cx(|cx| {
        let graph = parallel([2.0; 3]);
        let boundaries = [boundary(0, 0.0)];
        assert!(matches!(
            InertialGraph::try_new(
                cx,
                &graph,
                Density::new(gas(300.0).density),
                &geometry(3),
                &boundaries,
                &[1.0, 0.0, 0.0],
                0.0,
                16
            ),
            Err(InertialGraphError::Continuity { .. })
        ));
        let state = owner(cx, &graph, &[1.0, -0.5, -0.5], &boundaries);
        let initial = state.volume_flows_m3_s();
        let dt = 0.125;
        let l = state.inertances_pa_s2_m3()[0];
        // An exactly reversed endpoint has zero native average but nonzero
        // unresolved speed/loss. Caps must not erase that interval.
        let drive: Vec<_> = initial.iter().map(|q| -2.0 * l * q / dt).collect();
        let mut strict = controls();
        strict.maximum_neglected_power_w = 0.0;
        assert!(matches!(
            state.propose(cx, &graph, &boundaries, &drive, &[0.0; 3], dt, strict),
            Err(InertialGraphError::NeglectedWork { .. })
        ));
        let long = 3.0 * l / 2.0;
        assert!(matches!(
            state.propose(
                cx,
                &graph,
                &boundaries,
                &[0.0; 3],
                &initial,
                long,
                controls()
            ),
            Err(InertialGraphError::DragStep { .. })
        ));
        assert!(matches!(
            state.propose(
                cx,
                &graph,
                &[boundary(1, 0.0)],
                &[0.0; 3],
                &initial,
                dt,
                controls()
            ),
            Err(InertialGraphError::Invalid(_))
        ));
        assert_eq!(state.volume_flows_m3_s(), initial);
        assert_eq!(state.accepted_time_s(), 0.0);
    });
}

#[test]
fn g4_cancellation_refuses_a_prepared_momentum_publication() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 712,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let graph = parallel([2.0; 3]);
        let boundaries = [boundary(0, 0.0)];
        let state = owner(&cx, &graph, &[0.0; 3], &boundaries);
        let trial = state
            .propose(
                &cx,
                &graph,
                &boundaries,
                &[3.0, 0.0, 0.0],
                &[0.0; 3],
                0.125,
                controls(),
            )
            .unwrap();
        gate.request();
        assert!(matches!(
            trial.accept(&cx, 1e-10, 1e-10),
            Err(InertialGraphError::Graph(GraphError::Interrupted))
        ));
        assert_eq!(state.accepted_time_s(), 0.0);
    });
}
