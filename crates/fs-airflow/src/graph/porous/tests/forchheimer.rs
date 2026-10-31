//! G1/G3: calibrated pore inertia through the original signed graph and fan.
use super::*;
use fs_material::porous::permeability::{DarcyFlowConstraint, ForchheimerCalibration};

fn calibrated(permeability: f64, beta: f64, max_speed: f64) -> DarcyGasBarrier {
    let original = barrier(permeability, 0.0, max_speed);
    let calibration = ForchheimerCalibration::try_new(
        beta,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        original.support(),
        "synthetic superficial-velocity inertial calibration",
    )
    .unwrap();
    original
        .with_forchheimer_calibrations(&[calibration])
        .unwrap()
}

fn branch(barrier: DarcyGasBarrier, state: GasState, area: f64) -> PassiveBranch {
    PassiveBranch::from_darcy_forchheimer(
        "calibrated filter",
        0,
        1,
        barrier,
        state,
        area,
        0.01,
        0.0,
        source(),
        ToleranceBasis::Analytic,
    )
    .unwrap()
}

#[test]
fn original_material_and_gas_set_both_terms_and_signed_graph_flow() {
    for temperature in [300.0, 600.0] {
        let state = gas(temperature);
        for k in [1e-9, 2e-9] {
            for beta in [0.0, 1e4, 2e4] {
                let material = calibrated(k, beta, 1.0);
                let graph =
                    LossGraph::new_mixed(2, vec![branch(material.clone(), state, 0.02)]).unwrap();
                // Independent superficial-velocity pressure law, not the graph's
                // volume-flow coefficients or its quadratic inverse.
                let speed = 0.5;
                let dp = state.dynamic_viscosity * 0.01 / k * speed
                    + state.density * beta * 0.01 * speed * speed;
                for sign in [-1.0, 0.0, 1.0] {
                    let result = with_cx(|cx| {
                        graph.solve(&[boundary(0, sign * dp), boundary(1, 0.0)], config(), cx)
                    })
                    .unwrap();
                    close(result.branches[0].flow.value(), sign * speed * 0.02);
                    let retained = result.branches[0]
                        .linear_loss
                        .as_ref()
                        .unwrap()
                        .material()
                        .unwrap();
                    assert_eq!(retained.barrier, material);
                    assert_eq!(retained.gas, state);
                    assert!(
                        result.branches[0].flow.value() * result.branches[0].pressure_drop.value()
                            >= 0.0
                    );
                }
            }
        }
    }
}

#[test]
fn inertial_binding_cannot_be_dropped_or_changed_into_a_linear_edge() {
    let material = calibrated(1e-9, 1e4, 2.0);
    assert!(
        LinearLossElement::from_darcy(
            "not linear",
            material.clone(),
            gas(300.0),
            0.02,
            0.01,
            0.0,
            source(),
            ToleranceBasis::Analytic,
        )
        .is_err()
    );
    let PassiveBranch::LinearQuadratic {
        mut branch,
        linear_loss,
    } = branch(material, gas(300.0), 0.02)
    else {
        panic!("positive inertia needs both losses")
    };
    assert!(
        LossGraph::new_mixed(
            2,
            vec![PassiveBranch::Linear {
                from: 0,
                to: 1,
                loss: linear_loss.clone(),
            }]
        )
        .is_err()
    );
    branch.loss.resistance = LossResistance::new(2.0 * branch.loss.resistance.value());
    assert!(
        LossGraph::new_mixed(
            2,
            vec![PassiveBranch::LinearQuadratic {
                branch,
                linear_loss,
            }]
        )
        .is_err()
    );
}

#[test]
fn calibrated_inertia_changes_the_actual_fan_operating_point() {
    let state = gas(300.0);
    // With area=1, choose the independent calibration so dp=2Q+3Q².
    let k = state.dynamic_viscosity * 0.01 / 2.0;
    let beta = 3.0 / (state.density * 0.01);
    let graph =
        LossGraph::new_mixed(2, vec![branch(calibrated(k, beta, 5.0), state, 1.0)]).unwrap();
    let point =
        with_cx(|cx| graph.solve_with_fan(0, 1, &fan(), config(), fan_config(), cx)).unwrap();
    let expected = 200.0 / (12.0 + 1344.0_f64.sqrt());
    close(point.flow.value(), expected);
    close(
        point.network_pressure.value(),
        2.0 * expected + 3.0 * expected * expected,
    );
    assert_eq!(
        graph.linear_losses()[0]
            .as_ref()
            .unwrap()
            .material()
            .unwrap()
            .gas,
        state
    );
}

#[test]
fn nonlinear_pore_graph_keeps_actual_pressure_and_speed_support() {
    for (maximum_speed, from, to) in [(0.01, 100.0, 0.0), (5.0, 2000.0, 1900.0)] {
        let graph = LossGraph::new_mixed(
            2,
            vec![branch(
                calibrated(1e-9, 1e4, maximum_speed),
                gas(300.0),
                0.02,
            )],
        )
        .unwrap();
        assert!(matches!(
            with_cx(|cx| graph.solve(&[boundary(0, from), boundary(1, to)], config(), cx,)),
            Err(GraphError::InvalidBranch { branch: 0, .. })
        ));
    }
}
