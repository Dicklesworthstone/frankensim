use super::*;
use crate::graph::porous::{LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::graph::tests::{boundary, with_cx};
use crate::{SourceProvenance, ToleranceBasis};
use fs_qty::Area;

fn series() -> LossGraph {
    LossGraph::new_mixed(
        3,
        [(0, 1), (1, 2)]
            .into_iter()
            .enumerate()
            .map(|(i, (from, to))| PassiveBranch::Linear {
                from,
                to,
                loss: LinearLossElement::new(
                    format!("total-pressure-{i}"),
                    LinearLossResistance::new(1.0),
                    0.0,
                    SourceProvenance::new(
                        "analytic reversible junction fixture",
                        "total-pressure-v1",
                    ),
                    ToleranceBasis::Analytic,
                )
                .unwrap(),
            })
            .collect(),
    )
    .unwrap()
}
fn controls() -> InertialGraphControls {
    InertialGraphControls {
        maximum_drag_step: 2.0,
        maximum_neglected_power_w: 1e6,
        maximum_neglected_energy_j: 1e6,
    }
}
fn inertia() -> [InertialBranch; 2] {
    [
        InertialBranch {
            length_m: 1.0,
            flow_area_m2: 1.0,
        },
        InertialBranch {
            length_m: 0.75,
            flow_area_m2: 0.5,
        },
    ]
}
fn pressure() -> [FixedTotalPressure; 2] {
    [
        FixedTotalPressure {
            node: 0,
            total_pressure: Pressure::new(4.0),
        },
        FixedTotalPressure {
            node: 2,
            total_pressure: Pressure::new(0.0),
        },
    ]
}
fn paths() -> [PipeLoadGeometry; 2] {
    [
        PipeLoadGeometry::try_new(vec![[0.0; 3], [1.0, 0.0, 0.0]], Area::new(1.0), 1, || false)
            .unwrap(),
        PipeLoadGeometry::try_new(
            vec![[1.0, 0.0, 0.0], [1.0, 0.75, 0.0]],
            Area::new(0.5),
            1,
            || false,
        )
        .unwrap(),
    ]
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-11, "{a} != {b}");
}
// Flows are stored as impulses p = L*q and read back as p/L (inertial.rs, volume_flows_m3_s), so a flow
// returns to two roundings, not bit for bit (0.4 reads back as 0.4000000000000001): 2 eps relative.
fn assert_flows_round_trip(flows: &[f64], expected: f64) {
    assert_eq!(flows.len(), 2);
    for q in flows {
        assert!(
            (q - expected).abs() <= 2.0 * f64::EPSILON * expected.abs(),
            "{q} != {expected}"
        );
    }
}

#[test]
fn typed_owner_retains_native_solve_and_uses_actual_static_port_loads() {
    with_cx(|cx| {
        let graph = series();
        let original = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(2.0),
            &inertia(),
            &[boundary(0, 4.0), boundary(2, 0.0)],
            &[0.4; 2],
            0.0,
            8,
        )
        .unwrap();
        let owner = TotalPressureGraph::try_new(
            cx,
            &graph,
            Density::new(2.0),
            &inertia(),
            &pressure(),
            &[0.4; 2],
            0.0,
            8,
        )
        .unwrap();
        let raw = original
            .propose(
                cx,
                &graph,
                &[boundary(0, 4.0), boundary(2, 0.0)],
                &[0.0; 2],
                &[0.4; 2],
                0.2,
                controls(),
            )
            .unwrap();
        let trial = owner
            .propose(
                cx,
                &graph,
                &pressure(),
                &[0.0; 2],
                &[0.4; 2],
                0.2,
                controls(),
            )
            .unwrap();
        assert_eq!(raw.report().momentum.x, trial.report().momentum().x);
        assert_eq!(raw.report().momentum.y, trial.report().momentum().y);
        let report = trial.report();
        near(report.average_total_pressures_pa()[1], 30.0 / 13.0);
        assert!(
            report.branch_streams[1].average_static_pressures_pa[0]
                < report.branch_streams[0].average_static_pressures_pa[1]
        );
        let native_static = raw
            .junction_load(cx, 1, &paths(), Pressure::new(0.0), [0.0; 3])
            .unwrap();
        let total_pressure = trial
            .junction_load(cx, 1, &paths(), Pressure::new(0.0), [0.0; 3])
            .unwrap();
        assert_ne!(native_static.load.force_n, total_pressure.load.force_n);
        near(
            total_pressure.kinetic_energy_transport_j,
            report.node_streams[1].kinetic_energy_out_j,
        );
        near(report.node_streams[1].zero_storage_balance_residual_j, 0.0);
        near(report.mechanical_balance_residual_j, 0.0);
        let (accepted, report, loads) = trial
            .accept_with_assembly_loads(
                cx,
                &paths(),
                &[],
                Pressure::new(0.0),
                [0.0; 3],
                1e-10,
                1e-10,
            )
            .unwrap();
        for axis in 0..3 {
            near(
                loads.total.impulse_ns[axis],
                loads.pipes.iter().map(|l| l.impulse_ns[axis]).sum::<f64>()
                    + loads
                        .junctions
                        .iter()
                        .map(|l| l.load.impulse_ns[axis])
                        .sum::<f64>(),
            );
            near(
                loads.total.angular_impulse_nms[axis],
                loads
                    .pipes
                    .iter()
                    .map(|l| l.angular_impulse_nms[axis])
                    .sum::<f64>()
                    + loads
                        .junctions
                        .iter()
                        .map(|l| l.load.angular_impulse_nms[axis])
                        .sum::<f64>(),
            );
        }
        near(accepted.accepted_time_s(), 0.2);
        near(
            report.boundary_static_pressure_work_j + report.boundary_kinetic_energy_inflow_j,
            report.momentum().supplied,
        );
        assert_flows_round_trip(&owner.volume_flows_m3_s(), 0.4);
    });
}

#[test]
fn total_pressure_full_gravity_balance_and_same_trajectory_recheck() {
    use fs_material::fluid_transport::FluidTransportState;
    with_cx(|cx| {
        let material = FluidTransportState::try_new(2.0, 1000.0, 0.1, 0.001).unwrap();
        let graph = series()
            .with_gravity(
                &material,
                &[[0.0; 3], [1.0, 0.0, 0.0], [1.0, 0.75, 0.0]],
                [0.0, -2.0, 0.0],
                cx,
            )
            .unwrap();
        let owner = TotalPressureGraph::try_new(
            cx,
            &graph,
            Density::new(2.0),
            &inertia(),
            &pressure(),
            &[0.4; 2],
            0.0,
            8,
        )
        .unwrap();
        let trial = owner
            .propose(
                cx,
                &graph,
                &pressure(),
                &[0.0; 2],
                &[0.4; 2],
                0.2,
                controls(),
            )
            .unwrap();
        let original_x = trial.report().momentum().x.clone();
        let checked = trial
            .recheck(cx, &graph, &pressure(), &[0.0; 2], controls())
            .unwrap();
        assert_eq!(checked.report().momentum().x, original_x);
        near(
            checked.report().gravity_work_j,
            -3.0 * 0.2 * checked.report().average_volume_flows_m3_s()[1],
        );
        near(checked.report().mechanical_balance_residual_j, 0.0);
        checked
            .clone()
            .accept_with_assembly_loads(
                cx,
                &paths(),
                &[],
                Pressure::new(0.0),
                [0.0; 3],
                1e-10,
                1e-10,
            )
            .unwrap();
        // Fresh drives audit the same trajectory; they do not quietly solve
        // a different momentum history in order to make acceptance pass.
        let changed = checked
            .recheck(cx, &graph, &pressure(), &[1.0; 2], controls())
            .unwrap();
        assert_eq!(changed.report().momentum().x, original_x);
        assert!(changed.accept(cx, 1e-10, 1e-10).is_err());
        assert_eq!(owner.accepted_time_s(), 0.0);
    });
}

#[test]
fn total_pressure_final_cancellation_never_publishes_history_or_loads() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
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
        let graph = series();
        let owner = TotalPressureGraph::try_new(
            &cx,
            &graph,
            Density::new(2.0),
            &inertia(),
            &pressure(),
            &[0.4; 2],
            0.0,
            8,
        )
        .unwrap();
        let trial = owner
            .propose(
                &cx,
                &graph,
                &pressure(),
                &[0.0; 2],
                &[0.4; 2],
                0.2,
                controls(),
            )
            .unwrap();
        gate.request();
        assert!(
            trial
                .accept_with_assembly_loads(
                    &cx,
                    &paths(),
                    &[],
                    Pressure::new(0.0),
                    [0.0; 3],
                    1e-10,
                    1e-10,
                )
                .is_err()
        );
        assert_eq!(owner.accepted_time_s(), 0.0);
        assert_flows_round_trip(&owner.volume_flows_m3_s(), 0.4);
    });
}
