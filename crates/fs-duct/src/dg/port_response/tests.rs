use super::*;
use fs_material::gas::GasSpec;

fn duct(walls: bool, pressure: f64) -> DgDuct {
    let source = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.05,
        }],
    };
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let duct = DgDuct::new(
        &source,
        &gas,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::IdealOpen,
    )
    .unwrap();
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
