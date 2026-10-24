use super::*;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicApertureSpec};
use fs_dcontact::Obstacle;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};

fn load(pressure: f64) -> ApertureDgLoad {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let source = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.05,
        }],
    };
    let duct = DgDuct::new(
        &source,
        &gas,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_viscothermal(8, [50.0, 2000.0], 0.01)
    .unwrap();
    ApertureDgLoad::new(duct, 1.0 / 96000.0, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap()
}

fn aperture(load: &ApertureDgLoad, mode: BladePropagator, stiffness: f64) -> DynamicAperture {
    let rest = 0.0002;
    // An explicit no-contact fixture isolates the coupled pressure/flow and work
    // from independent contact quadrature. No instrument parameters are fitted.
    let lay = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        0.0,
        2.0,
        "DG reciprocal coupling test: explicit zero contact".into(),
    )
    .unwrap();
    DynamicAperture::new_with_propagator(
        DynamicApertureSpec {
            aperture: BernoulliAperture {
                rest_opening_m: rest,
                width_m: 0.001,
                closing_pressure_pa: stiffness * rest / 1e-5,
            },
            mass_kg: 1e-5,
            stiffness_n_m: stiffness,
            damping_ratio: 0.02,
            density_kg_m3: load.duct().inlet_gas().density,
            impedance_pa_s_m3: load.impedance_pa_s_m3(),
            time_step_s: load.propagator().dt(),
            max_steps: 256,
        },
        ApertureState {
            opening_m: rest,
            opening_velocity_m_s: 0.0,
        },
        lay,
        mode,
    )
    .unwrap()
}

fn model(pressure: f64, stiffness: f64) -> ApertureDg {
    let load = load(pressure);
    ApertureDg::new(aperture(&load, BladePropagator::ExactZoh, stiffness), load).unwrap()
}

#[test]
fn exact_dg_and_aperture_share_mean_pressure_flow_and_work() {
    let mut model = model(101325.0, 500.0);
    let mut inlet_energy = 0.0;
    for index in 0..128 {
        let drive = TubeDrive {
            upstream_pressure_pa: if index < 64 { 30.0 } else { 0.0 },
            body_flow_m3_s: if index < 16 { 2e-9 } else { 0.0 },
        };
        let force = if index < 32 { 1e-5 } else { 0.0 };
        let frame = model.step_with_force(drive, force).unwrap();
        let port = frame.duct.port_outputs[0];
        let scale = frame.stored_energy_j.abs()
            + frame.storage_change_j.abs()
            + frame.inferred_total_dissipation_j.abs()
            + frame.upstream_work_j.abs()
            + frame.body_work_j.abs()
            + frame.aperture.mechanical_work_j.abs();
        assert!(frame.balance_residual_j().abs() < 1e-8 * scale.max(1e-25));
        assert!(frame.port_work_residual_j.abs() < 1e-11 * scale.max(1e-25));
        assert!(
            (port.conjugate_rate - frame.aperture.bore_pressure_pa).abs()
                < 1e-9 * (1.0 + port.conjugate_rate.abs())
        );
        assert!(!frame.duct.ledger.balance_checked);
        assert!(frame.duct.ledger.cumulative_dissipated_j >= -1e-25);
        assert_eq!(frame.aperture.step, model.duct_state().step_count());
        inlet_energy += port.conjugate_rate * port.conjugate_rate;
    }
    assert!(inlet_energy > 0.0);
    assert!(model.duct_state().pressure().iter().any(|p| p.abs() > 1e-6));
}

#[test]
fn rejected_drive_preserves_both_accepted_histories() {
    let mut model = model(101325.0, 500.0);
    model
        .step(TubeDrive {
            upstream_pressure_pa: 20.0,
            body_flow_m3_s: 0.0,
        })
        .unwrap();
    let state = model.duct_state().clone();
    let opening = model.aperture().state();
    let steps = model.aperture().accepted_steps();
    for (drive, force) in [
        (
            TubeDrive {
                upstream_pressure_pa: f64::NAN,
                body_flow_m3_s: 0.0,
            },
            0.0,
        ),
        (TubeDrive::default(), f64::INFINITY),
    ] {
        assert!(model.step_with_force(drive, force).is_err());
        assert_eq!(model.duct_state(), &state);
        assert_eq!(model.aperture().state(), opening);
        assert_eq!(model.aperture().accepted_steps(), steps);
    }
    let replayed_load = model.load().clone();
    assert!(
        ApertureDg::new(
            aperture(&replayed_load, BladePropagator::ExactZoh, 500.0),
            replayed_load,
        )
        .is_err(),
        "an advanced acoustic history cannot bind to a fresh aperture clock"
    );
}

#[test]
fn gas_pressure_and_structural_stiffness_change_the_coupled_response() {
    let mut base = model(101325.0, 500.0);
    let mut gas = model(60000.0, 500.0);
    let mut structure = model(101325.0, 900.0);
    let mut gas_difference = 0.0;
    let mut structure_difference = 0.0;
    for index in 0..64 {
        let drive = TubeDrive {
            upstream_pressure_pa: if index < 32 { 30.0 } else { 0.0 },
            body_flow_m3_s: 0.0,
        };
        let reference = base.step(drive).unwrap().aperture.bore_pressure_pa;
        gas_difference += (reference - gas.step(drive).unwrap().aperture.bore_pressure_pa).abs();
        structure_difference +=
            (reference - structure.step(drive).unwrap().aperture.bore_pressure_pa).abs();
    }
    assert!(gas_difference > 1e-5);
    assert!(structure_difference > 1e-5);
}

#[test]
fn cold_binding_refuses_midpoint_mismatch_and_insufficient_state_budget() {
    let load = load(101325.0);
    assert!(
        ApertureDg::new(
            aperture(&load, BladePropagator::Midpoint, 5.0),
            load.clone()
        )
        .is_err()
    );
    assert!(ApertureDgLoad::new(load.duct().clone(), load.propagator().dt(), 1e-11, 1).is_err());
    let other = self::load(60000.0);
    assert!(ApertureDg::new(aperture(&other, BladePropagator::ExactZoh, 500.0), load).is_err());
}

#[test]
fn moving_aperture_retains_causal_terminal_history_and_signed_energy_exchange() {
    use fs_math::c64::C64;
    use fs_vfit::mimo::MimoModel;
    use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
    use fs_vfit::model::PoleTerm;
    for with_side in [false, true] {
    let dt = 1.0 / 96000.0;
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }] },
        &gas, DgMesh { elements: if with_side { 2 } else { 1 }, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Port(PortId(1)),
    ).unwrap().with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap();
    let duct = if with_side {
        duct.with_side_ports(&[fs_duct::dg::DgSidePort {
            id: PortId(2), downstream_element: 1,
        }]).unwrap()
    } else { duct };
    let mut model = MimoModel {
        ports: 1, poles: vec![PoleTerm::Real { pole: -2000.0, residue: 0.0 }],
        residues: vec![vec![C64::from_re(2e11)]], d: vec![1e8], e: vec![0.0],
        reciprocal: true,
    };
    if with_side {
        model.ports = 2;
        model.residues = vec![vec![C64::from_re(2e11), C64::from_re(4e10),
            C64::from_re(4e10), C64::from_re(1.5e11)]];
        model.d = vec![1e8, 2e7, 2e7, 1.2e8];
        model.e = vec![0.0; 4];
    }
    let terminal = PreparedMimoLoad::new(&model, dt, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / dt],
        sampled_passivity_points: 257, max_states: 4, max_memory_bytes: 1 << 20,
    }).unwrap();
    let load = if with_side {
        ApertureDgLoad::new_with_terminal_ports(duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS,
            vec![PortId(2), PortId(1)], terminal)
    } else {
        ApertureDgLoad::new_with_terminal(duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal)
    }.unwrap();
    let mut model = ApertureDg::new(
        aperture(&load, BladePropagator::ExactZoh, 500.0), load,
    ).unwrap();
    let mut transferred = 0.0;
    for sample in 0..128 {
        let frame = model.step(TubeDrive {
            upstream_pressure_pa: if sample < 64 { 30.0 } else { 0.0 },
            body_flow_m3_s: 0.0,
        }).unwrap();
        let terminal = frame.terminal.as_ref().unwrap();
        let scale = frame.storage_change_j.abs() + frame.inferred_total_dissipation_j.abs()
            + frame.upstream_work_j.abs() + frame.aperture.mechanical_work_j.abs()
            + frame.terminals.iter().map(|port| port.port_work_j.abs()).sum::<f64>();
        assert!(frame.balance_residual_j().abs() <= 1e-8 * scale.max(1e-25));
        assert!(frame.port_work_residual_j.abs() <= 1e-10 * scale.max(1e-25));
        assert!(terminal.port_work_residual_j.abs() <= 1e-10 * scale.max(1e-25));
        assert_eq!(frame.terminals.len(), if with_side { 2 } else { 1 });
        assert_eq!(model.load().last_terminal_frames(), frame.terminals);
        assert_eq!(model.load().last_terminal_frame(), frame.terminal.as_ref());
        transferred += sum_terminal_work(&frame.terminals);
    }
    assert!(transferred > 0.0);
    let state = model.duct_state().clone();
    let poles = model.load().terminal().unwrap().pole_states().to_vec();
    let terminals = model.load().last_terminal_frames().to_vec();
    let opening = model.aperture().state();
    assert!(model.step(TubeDrive {
        upstream_pressure_pa: f64::NAN, body_flow_m3_s: 0.0,
    }).is_err());
    assert_eq!(model.duct_state(), &state);
    assert_eq!(model.load().terminal().unwrap().pole_states(), poles);
    assert_eq!(model.load().last_terminal_frames(), terminals);
    assert_eq!(model.aperture().state(), opening);
    }
}
