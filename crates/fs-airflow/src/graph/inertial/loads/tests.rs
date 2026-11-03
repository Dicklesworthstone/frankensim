use super::*;
use crate::graph::LossGraph;
use crate::graph::gravity::GraphGravitySection;
use crate::graph::inertial::{InertialBranch, InertialGraph, InertialGraphControls};
use crate::graph::porous::{LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::graph::tests::{boundary, with_cx};
use crate::{SourceProvenance, ToleranceBasis};
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{Area, Density};

fn graph() -> LossGraph {
    LossGraph::new_mixed(
        2,
        vec![PassiveBranch::Linear {
            from: 0,
            to: 1,
            loss: LinearLossElement::new(
                "pipe",
                LinearLossResistance::new(4.0),
                0.0,
                SourceProvenance::new("analytic pipe fixture", "pipe-momentum-v1"),
                ToleranceBasis::Analytic,
            )
            .unwrap(),
        }],
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
fn near(a: [f64; 3], b: [f64; 3]) {
    for (a, b) in a.into_iter().zip(b) {
        assert!((a - b).abs() < 1e-9, "{a} != {b}");
    }
}

#[test]
fn g1_native_pressure_trajectory_loads_only_the_remaining_wall_impulse() {
    with_cx(|cx| {
        let graph = graph();
        let boundaries = [boundary(0, 61.0), boundary(1, 7.0)];
        let owner = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(2.0),
            &[InertialBranch {
                length_m: 4.0,
                flow_area_m2: 0.25,
            }],
            &boundaries,
            &[0.0],
            0.0,
            4,
        )
        .unwrap();
        let trial = owner
            .propose(cx, &graph, &boundaries, &[0.0], &[1.5], 2.0, controls())
            .unwrap();
        let geometry =
            PipeLoadGeometry::try_new(vec![[0.0; 3], [4.0, 0.0, 0.0]], Area::new(0.25), 1, || {
                false
            })
            .unwrap();
        let load = trial
            .pipe_load(cx, 0, &geometry, Pressure::new(2.0), [0.0; 3])
            .unwrap();
        near(load.force_n, [1.5, 0.0, 0.0]);
        near(load.impulse_ns, [3.0, 0.0, 0.0]);
        // Same L would also result from length=2, area=0.125, but its actual
        // fluid momentum and pressure traction are different. Do not accept it.
        let wrong =
            PipeLoadGeometry::try_new(vec![[0.0; 3], [2.0, 0.0, 0.0]], Area::new(0.125), 1, || {
                false
            })
            .unwrap();
        assert!(
            trial
                .pipe_load(cx, 0, &wrong, Pressure::new(2.0), [0.0; 3])
                .is_err()
        );
        assert!(
            trial
                .pipe_load(cx, 1, &geometry, Pressure::new(2.0), [0.0; 3])
                .is_err()
        );
        let (next, _) = trial.accept(cx, 1e-9, 1e-9).unwrap();
        assert!((next.volume_flows_m3_s()[0] - 3.0).abs() < 1e-10);
        assert_eq!(owner.volume_flows_m3_s(), vec![0.0]);
    });
}

#[test]
fn g1_rechecked_original_distributed_gravity_owns_pipe_weight_and_moment() {
    with_cx(|cx| {
        let path = vec![[0.0; 3], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0]];
        let material = |density| FluidTransportState::try_new(density, 1000.0, 0.1, 0.001).unwrap();
        let with_density = |density: [f64; 2]| {
            let sections = path
                .windows(2)
                .enumerate()
                .map(|(i, pair)| GraphGravitySection {
                    start_m: pair[0],
                    end_m: pair[1],
                    transport: material(density[i]),
                })
                .collect();
            graph()
                .with_distributed_boussinesq_gravity(
                    &material(1.0),
                    &[sections],
                    &[path[0], path[2]],
                    [0.0, 0.0, -10.0],
                    0.5,
                    2,
                    cx,
                )
                .unwrap()
        };
        let original = with_density([1.0; 2]);
        let refreshed = with_density([1.4, 0.6]);
        let boundaries = [boundary(0, 0.0), boundary(1, 0.0)];
        let owner = InertialGraph::try_new(
            cx,
            &original,
            Density::new(1.0),
            &[InertialBranch {
                length_m: 3.0,
                flow_area_m2: 0.5,
            }],
            &boundaries,
            &[0.0],
            0.0,
            4,
        )
        .unwrap();
        let trial = owner
            .propose(cx, &original, &boundaries, &[0.0], &[0.0], 1.0, controls())
            .unwrap()
            .recheck(cx, &refreshed, &boundaries, &[0.0], controls())
            .unwrap();
        let geometry = PipeLoadGeometry::try_new(path, Area::new(0.5), 2, || false).unwrap();
        let load = trial
            .pipe_load(cx, 0, &geometry, Pressure::new(0.0), [0.0; 3])
            .unwrap();
        near(load.force_n, [0.0, 0.0, -17.0]);
        near(load.torque_nm, [-1.5, 20.0, 0.0]);
        // Matching total length and area cannot substitute for original source sections.
        let other_path = PipeLoadGeometry::try_new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
            Area::new(0.5),
            3,
            || false,
        )
        .unwrap();
        assert!(
            trial
                .pipe_load(cx, 0, &other_path, Pressure::new(0.0), [0.0; 3])
                .is_err()
        );
        trial.accept(cx, 1e-9, 1e-9).unwrap();
    });
}
