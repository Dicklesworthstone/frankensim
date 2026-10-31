//! G0/G1/G3: hydrostatics, Bernoulli head, circulation, material and frame changes.

use super::*;
use crate::graph::sensitivity::{FlowWeight, GraphSensitivityConfig, PressureWeight};
use crate::graph::tests::{boundary, bridge, close, config, edge, with_cx};
use crate::graph::{GraphFanConfig, LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::{FanArrangement, FanBank, FanCurve, FanPoint, SourceProvenance, ToleranceBasis};
use fs_qty::{Density, Temperature, VolumetricFlowRate};
use fs_solver::CheckedCgConfig;

fn material(density: f64) -> FluidTransportState {
    FluidTransportState::try_new(density, 1000.0, 0.1, 0.001).unwrap()
}

#[test]
fn g1_uniform_hydrostatic_loop_is_exactly_quiescent_in_arbitrary_axes() {
    let positions = [
        [1.0, 2.0, 3.0],
        [2.0, 4.0, 6.0],
        [-1.0, 7.0, 2.0],
        [9.0, -2.0, -3.0],
    ];
    let acceleration = [2.0, -3.0, -4.0];
    for shift in [[0.0; 3], [1.0e8, -2.0e8, 3.0e8]] {
        let translated = positions.map(|p| core::array::from_fn(|i| p[i] + shift[i]));
        let result = with_cx(|cx| {
            bridge()
                .with_gravity(&material(3.0), &translated, acceleration, cx)
                .unwrap()
                .solve(&[boundary(0, 1.0e5)], config(), cx)
                .unwrap()
        });
        assert_eq!(result.sweeps, 0);
        for (node, p) in result.pressures.iter().enumerate() {
            let dot: f64 = (0..3)
                .map(|axis| acceleration[axis] * (positions[node][axis] - positions[0][axis]))
                .sum();
            assert_eq!(p.value(), 1.0e5 + 3.0 * dot);
        }
        for branch in &result.branches {
            assert_eq!(branch.flow.value(), 0.0);
            assert_eq!(branch.dissipative_pressure_drop.value(), 0.0);
            assert_eq!(
                branch.pressure_drop.value() + branch.gravity_pressure_head.value(),
                0.0
            );
        }
    }
}

#[test]
fn g1_gravity_and_density_drive_flow_and_dissipation_separately_from_static_pressure() {
    for (density, sign) in [(3.0, 1.0), (6.0, 1.0), (3.0, -1.0)] {
        let result = with_cx(|cx| {
            LossGraph::new(2, vec![edge("drop", 0, 1, 4.0)])
                .unwrap()
                .with_gravity(
                    &material(density),
                    &[[0.0, 0.0, 2.0], [0.0; 3]],
                    [0.0, 0.0, -10.0 * sign],
                    cx,
                )
                .unwrap()
                .solve(&[boundary(0, 0.0), boundary(1, 0.0)], config(), cx)
                .unwrap()
        });
        let branch = &result.branches[0];
        let head = sign * density * 20.0;
        close(branch.flow.value(), sign * (head.abs() / 4.0).sqrt());
        assert_eq!(branch.pressure_drop.value(), 0.0);
        assert_eq!(branch.gravity_pressure_head.value(), head);
        assert_eq!(branch.dissipative_pressure_drop.value(), head);
        assert!(branch.dissipative_pressure_drop.value() * branch.flow.value() > 0.0);
    }
}

#[test]
fn g1_material_density_difference_circulates_through_the_original_closed_graph() {
    let positions = [
        [0.0; 3],
        [0.0, 0.0, 10.0],
        [1.0, 0.0, 10.0],
        [1.0, 0.0, 0.0],
    ];
    let make = || {
        LossGraph::new(
            4,
            vec![
                edge("hot rise", 0, 1, 1.0),
                edge("top", 1, 2, 1.0),
                edge("cold fall", 2, 3, 1.0),
                edge("bottom", 3, 0, 1.0),
            ],
        )
        .unwrap()
    };
    with_cx(|cx| {
        let states = [material(0.9), material(1.0), material(1.0), material(1.0)];
        let graph = make()
            .with_boussinesq_gravity(
                &material(1.0),
                &states,
                &positions,
                [0.0, 0.0, -10.0],
                0.11,
                cx,
            )
            .unwrap();
        let solution = graph.solve(&[boundary(0, 0.0)], config(), cx).unwrap();
        // Loop head = (rho_cold-rho_hot)*g*height; four unit quadratic losses.
        let expected = 2.5_f64.sqrt();
        for branch in &solution.branches {
            close(branch.flow.value(), expected);
        }
        let gravity_power: f64 = solution
            .branches
            .iter()
            .map(|b| b.gravity_pressure_head.value() * b.flow.value())
            .sum();
        let dissipated: f64 = solution
            .branches
            .iter()
            .map(|b| b.dissipative_pressure_drop.value() * b.flow.value())
            .sum();
        close(gravity_power, 10.0 * expected);
        close(dissipated, gravity_power);
        assert!(solution.boundary_imbalance.value().abs() < 1e-9);
        let repeated = graph.solve(&[boundary(0, 0.0)], config(), cx).unwrap();
        assert_eq!(solution, repeated);
        let uniform = make()
            .with_gravity(&material(1.0), &positions, [0.0, 0.0, -10.0], cx)
            .unwrap()
            .solve(&[boundary(0, 0.0)], config(), cx)
            .unwrap();
        assert!(uniform.branches.iter().all(|b| b.flow.value() == 0.0));
    });
}

#[test]
fn g3_zero_gravity_preserves_original_flows_and_orientation_is_only_bookkeeping() {
    with_cx(|cx| {
        let plain = bridge();
        let boundaries = [boundary(0, 25.0), boundary(3, 0.0)];
        let before = plain.solve(&boundaries, config(), cx).unwrap();
        let zero = plain
            .with_gravity(&material(2.0), &[[0.0; 3]; 4], [0.0; 3], cx)
            .unwrap()
            .solve(&boundaries, config(), cx)
            .unwrap();
        assert_eq!(before.pressures, zero.pressures);
        assert_eq!(before.branches, zero.branches);
        let forward = LossGraph::new(2, vec![edge("drop", 0, 1, 4.0)]).unwrap();
        let backward = LossGraph::new(2, vec![edge("drop", 1, 0, 4.0)]).unwrap();
        let positions = [[0.0, 0.0, 1.0], [0.0; 3]];
        let a = forward
            .with_gravity(&material(3.0), &positions, [0.0, 0.0, -10.0], cx)
            .unwrap()
            .solve(&[boundary(0, 0.0), boundary(1, 0.0)], config(), cx)
            .unwrap();
        let b = backward
            .with_gravity(&material(3.0), &positions, [0.0, 0.0, -10.0], cx)
            .unwrap()
            .solve(&[boundary(0, 0.0), boundary(1, 0.0)], config(), cx)
            .unwrap();
        assert_eq!(a.branches[0].flow.value(), -b.branches[0].flow.value());
        assert_eq!(a.node_outflows, b.node_outflows);
    });
}

#[test]
fn g1_gravity_uses_real_fan_intersection_and_admits_a_zero_pressure_source() {
    with_cx(|cx| {
        for (intercept, slope, height) in [(100.0, 10.0, 1.0), (0.0, 0.0, -1.0)] {
            let fan = FanBank::new(
                FanCurve::new(
                    "analytic",
                    vec![
                        FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(intercept)),
                        FanPoint::new(
                            VolumetricFlowRate::new(10.0),
                            Pressure::new(intercept - slope * 10.0),
                        ),
                    ],
                    SourceProvenance::new("analytic", "gravity-fan"),
                    0.0,
                    ToleranceBasis::Analytic,
                    VolumetricFlowRate::new(0.0),
                    (0.5, 2.0),
                )
                .unwrap(),
                1,
                FanArrangement::Series,
                1.0,
            )
            .unwrap();
            let graph = LossGraph::new(2, vec![edge("duct", 0, 1, 4.0)])
                .unwrap()
                .with_gravity(
                    &material(3.0),
                    &[[0.0; 3], [0.0, 0.0, height]],
                    [0.0, 0.0, -10.0],
                    cx,
                )
                .unwrap();
            let point = graph
                .solve_with_fan(
                    0,
                    1,
                    &fan,
                    config(),
                    GraphFanConfig {
                        max_iterations: 128,
                        absolute_flow_tolerance: VolumetricFlowRate::new(1e-9),
                        absolute_pressure_tolerance: Pressure::new(1e-8),
                    },
                    cx,
                )
                .unwrap();
            let expected =
                (-slope + (slope * slope + 16.0 * (intercept - 30.0 * height)).sqrt()) / 8.0;
            close(point.flow.value(), expected);
            close(point.network_pressure.value(), intercept - slope * expected);
        }
    });
}

#[test]
fn g1_fixed_gravity_sensitivities_differentiate_dissipative_not_static_pressure() {
    with_cx(|cx| {
        for linear in [false, true] {
            let graph = if linear {
                LossGraph::new_mixed(
                    2,
                    vec![PassiveBranch::Linear {
                        from: 0,
                        to: 1,
                        loss: LinearLossElement::new(
                            "pipe",
                            LinearLossResistance::new(4.0),
                            0.0,
                            SourceProvenance::new("analytic", "gravity-linear"),
                            ToleranceBasis::Analytic,
                        )
                        .unwrap(),
                    }],
                )
                .unwrap()
            } else {
                LossGraph::new(2, vec![edge("pipe", 0, 1, 4.0)]).unwrap()
            };
            let graph = graph
                .with_gravity(
                    &material(1.0),
                    &[[0.0, 0.0, 1.0], [0.0; 3]],
                    [0.0, 0.0, -9.0],
                    cx,
                )
                .unwrap();
            let solved = graph
                .solve_linearized(
                    &[boundary(0, 0.0), boundary(1, 0.0)],
                    config(),
                    GraphSensitivityConfig {
                        minimum_absolute_flow: VolumetricFlowRate::new(1e-12),
                        linear: CheckedCgConfig {
                            tolerance: 1e-12,
                            max_iterations: 64,
                            max_corrections: 2,
                        },
                    },
                    cx,
                )
                .unwrap();
            let gradient = solved
                .adjoint(&[PressureWeight::new(0.0); 2], &[FlowWeight::new(1.0)], cx)
                .unwrap();
            let q = solved.solution().branches[0].flow.value();
            close(
                gradient.log_resistances[0],
                if linear { -q } else { -0.5 * q },
            );
            close(
                gradient.boundary_pressures[0].value(),
                if linear { 0.25 } else { 1.0 / 12.0 },
            );
        }
    });
}

#[test]
fn g1_downward_heat_transport_accepts_rising_static_pressure_and_shared_storage() {
    use crate::graph::thermal::transport::{
        BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork,
    };
    with_cx(|cx| {
        let state = material(2.0);
        let flow = LossGraph::new(2, vec![edge("down", 0, 1, 4.0)])
            .unwrap()
            .with_gravity(&state, &[[0.0, 0.0, 1.0], [0.0; 3]], [0.0, 0.0, -10.0], cx)
            .unwrap()
            .solve(&[boundary(0, 0.0), boundary(1, 4.0)], config(), cx)
            .unwrap();
        assert_eq!(flow.branches[0].flow.value(), 2.0);
        let air = TransportAir {
            density: Density::new(2.0),
            specific_heat_j_kg_k: 1000.0,
        };
        let budget = TransportConfig {
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-9),
            relative_flow_tolerance: 1e-9,
            absolute_heat_tolerance_w: 1e-7,
            relative_heat_tolerance: 1e-9,
        };
        let inlets = [TransportInlet {
            node: 0,
            temperature: Temperature::new(300.0),
        }];
        let net = TransportNetwork::new(
            cx,
            &flow,
            air,
            vec![BranchThermalModel::Adiabatic],
            &inlets,
            budget,
        )
        .unwrap();
        let result = net.march(cx, &[]).unwrap();
        assert_eq!(result.node_temperatures_k[1], Some(300.0));
        assert_eq!(result.wall_heat_rate_w, 0.0);
        assert!(
            TransportNetwork::new(
                cx,
                &flow,
                TransportAir {
                    density: Density::new(1.0),
                    ..air
                },
                vec![BranchThermalModel::Adiabatic],
                &inlets,
                budget
            )
            .is_err()
        );
    });
}

#[test]
fn g0_invalid_geometry_material_departures_and_nonfinite_heads_refuse() {
    with_cx(|cx| {
        let make = || LossGraph::new(2, vec![edge("pipe", 0, 1, 4.0)]).unwrap();
        assert!(
            make()
                .with_gravity(&material(1.0), &[[0.0; 3]], [0.0; 3], cx)
                .is_err()
        );
        assert!(
            make()
                .with_gravity(&material(1.0), &[[0.0; 3]; 2], [f64::NAN, 0.0, 0.0], cx)
                .is_err()
        );
        assert!(
            make()
                .with_gravity(
                    &material(2.0),
                    &[[0.0; 3], [1.0, 0.0, 0.0]],
                    [f64::MAX, 0.0, 0.0],
                    cx
                )
                .is_err()
        );
        for budget in [0.01, -1.0, 1.0, f64::NAN] {
            assert!(
                make()
                    .with_boussinesq_gravity(
                        &material(1.0),
                        &[material(0.9)],
                        &[[0.0; 3]; 2],
                        [0.0; 3],
                        budget,
                        cx
                    )
                    .is_err()
            );
        }
        let graph = make()
            .with_gravity(&material(1.0), &[[0.0; 3]; 2], [0.0; 3], cx)
            .unwrap();
        assert!(graph.solve(&[boundary(99, 0.0)], config(), cx).is_err());
    });
}
