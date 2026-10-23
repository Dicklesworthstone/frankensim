use super::*;
use fs_duct::dg::{DgMesh, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};

fn duct(inlet: DgBoundary, outlet: DgBoundary, lossy: bool) -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let source = Duct {
        segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }],
    };
    let source = DgDuct::new(
        &source, &gas, DgMesh { elements: 1, order: 2 }, inlet, outlet,
    ).unwrap();
    if lossy {
        source.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
    } else {
        source
    }
}

fn load(lossy: bool) -> DgPortLoad {
    DgPortLoad::new(
        duct(DgBoundary::Port(PortId(7)), DgBoundary::IdealOpen, lossy),
        1.0 / 96000.0, 1e-11, MAX_DG_PROPAGATOR_DOFS,
    ).unwrap()
}

#[test]
fn g0_transaction_retains_the_original_complete_propagation_and_work() {
    for lossy in [false, true] {
        let mut load = load(lossy);
        let propagator = load.propagator().clone();
        let mut reference = load.state().clone();
        let inlet = load.inlet();
        for sample in 0..32 {
            let flow = if sample < 16 { 2e-9 } else { 0.0 };
            let previous = load.state().clone();
            let relation = load.port_response().unwrap();
            let expected = propagator.step(
                &mut reference, &[PortInput { id: inlet, effort: flow }],
            ).unwrap();
            propagator.sync_nodal(&mut reference).unwrap();
            let trial = load.preview_step(flow).unwrap();
            assert_eq!(trial.previous_state(), &previous);
            assert_eq!(trial.state(), &reference);
            assert_eq!(trial.output(), &expected);
            let pressure = trial.output().port_outputs[0].conjugate_rate;
            let affine = relation.history_pressure_pa + relation.impedance_pa_s_m3 * flow;
            let scale = relation.history_pressure_pa.abs()
                + (relation.impedance_pa_s_m3 * flow).abs();
            assert!((pressure - affine).abs() <= 1e-10 * scale.max(1e-300));
            let accepted = trial.accept();
            assert_eq!(accepted, expected);
            assert_eq!(load.state(), &reference);
            assert_eq!(load.state().step_count(), sample + 1);
        }
    }
}

#[test]
fn g0_abandoned_or_invalid_trial_keeps_the_accepted_history() {
    let mut load = load(true);
    load.preview_step(1e-9).unwrap().accept();
    let previous = load.state().clone();
    let response = load.port_response().unwrap();
    {
        let trial = load.preview_step(2e-9).unwrap();
        assert_eq!(trial.state().step_count(), previous.step_count() + 1);
        assert_eq!(trial.previous_state(), &previous);
    }
    assert_eq!(load.state(), &previous);
    assert_eq!(load.port_response().unwrap(), response);
    for flow in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(load.preview_step(flow).is_err());
        assert_eq!(load.state(), &previous);
        assert_eq!(load.port_response().unwrap(), response);
    }
}

#[test]
fn g0_cold_admission_refuses_unowned_ports_and_unresolved_clock() {
    for (inlet, outlet) in [
        (DgBoundary::Closed, DgBoundary::IdealOpen),
        (DgBoundary::Port(PortId(7)), DgBoundary::Port(PortId(9))),
    ] {
        assert!(DgPortLoad::new(
            duct(inlet, outlet, false), 1.0 / 96000.0, 1e-11, MAX_DG_PROPAGATOR_DOFS,
        ).is_err());
    }
    let source = duct(DgBoundary::Port(PortId(7)), DgBoundary::IdealOpen, true);
    assert!(DgPortLoad::new(source.clone(), 1.0 / 96000.0, 1e-11, 1).is_err());
    assert!(DgPortLoad::new(source.clone(), 1.0 / 1000.0, 1e-11, MAX_DG_PROPAGATOR_DOFS).is_err());
    assert!(DgPortLoad::new(source, f64::NAN, 1e-11, MAX_DG_PROPAGATOR_DOFS).is_err());
}
