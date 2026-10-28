//! Coupled port physics; these fixtures do not qualify an instrument.
use super::*;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::bernoulli_aperture::dg::LoadedApertureDg;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_dcontact::Obstacle;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::model::PoleTerm;

fn load(length: f64, pressure: f64, dt: f64) -> ApertureDgLoad {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length }] },
        &gas, DgMesh { elements: 1, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Port(PortId(1)),
    ).unwrap().with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap();
    let terminal = PreparedMimoLoad::new(&MimoModel {
        ports: 1, poles: vec![PoleTerm::Real { pole: -2000.0, residue: 0.0 }],
        residues: vec![vec![C64::from_re(2e11)]], d: vec![1e8], e: vec![0.0],
        reciprocal: true,
    }, dt, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / dt], sampled_passivity_points: 257,
        max_states: 4, max_memory_bytes: 1 << 20,
    }).unwrap();
    ApertureDgLoad::new_with_terminal(
        duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal,
    ).unwrap()
}

fn aperture(down: &ApertureDgLoad, impedance: f64, mode: BladePropagator) -> DynamicAperture {
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0],
        0.0, 2.0, "DG tract manufactured no-contact fixture".into()).unwrap();
    DynamicAperture::new_with_propagator(DynamicApertureSpec {
        aperture: BernoulliAperture {
            rest_opening_m: 0.0002, width_m: 0.001, closing_pressure_pa: 100.0,
        },
        mass_kg: 1e-5, stiffness_n_m: 5.0, damping_ratio: 0.02,
        density_kg_m3: down.duct().inlet_gas().density,
        impedance_pa_s_m3: impedance, time_step_s: down.propagator().dt(), max_steps: 256,
    }, ApertureState { opening_m: 0.0002, opening_velocity_m_s: 0.0 }, contact, mode).unwrap()
}

fn system(length: f64) -> LoadedApertureDg {
    let up = load(length, 101325.0, 1.0 / 96000.0);
    let down = load(0.05, 101325.0, 1.0 / 96000.0);
    let (_, z) = series_inlet_state(&up, &down).unwrap();
    LoadedApertureDg::new_loaded_with_upstream(
        aperture(&down, z, BladePropagator::ExactZoh), up, down,
    ).unwrap()
}

fn single(model: &LoadedApertureDg) -> &DynamicAperture {
    let LoadedAperture::Single(source) = model.aperture() else { unreachable!() };
    source
}

#[test]
fn both_exact_ducts_share_actual_face_pressure_flow_and_signed_work() {
    let mut model = system(0.035);
    let dt = model.load().propagator().dt();
    let mut nonzero_tract = false;
    let mut reverse_flow = false;
    for step in 0..128 {
        let drive = TubeDrive {
            upstream_pressure_pa: if step < 40 { 0.5 } else if step < 80 { -0.5 } else { 0.0 },
            body_flow_m3_s: if step < 16 { 2e-10 } else if step < 32 { -1e-10 } else { 0.0 },
        };
        let frame = model.step_with_forces(drive,
            [if step < 24 { 1e-6 } else { 0.0 }, 0.0]).unwrap();
        let tract = frame.upstream.as_ref().unwrap();
        let through = frame.aperture.bore_flow_m3_s - drive.body_flow_m3_s;
        let scale = frame.stored_energy_j.abs() + frame.storage_change_j.abs()
            + frame.inferred_total_dissipation_j.abs() + frame.upstream_work_j.abs()
            + frame.body_work_j.abs() + frame.aperture.mechanical_work_j.abs()
            + sum_terminal_work(&frame.terminals).abs() + sum_terminal_work(&tract.terminals).abs();
        assert!(frame.balance_residual_j().abs() <= 1e-8 * scale.max(1e-24));
        assert!(frame.port_work_residual_j.abs() <= 1e-10 * scale.max(1e-24));
        assert_eq!(tract.inlet_flow_m3_s, -through);
        assert_eq!(tract.upstream_face_pressure_pa,
            drive.upstream_pressure_pa + tract.inlet_pressure_pa);
        assert_eq!(frame.upstream_work_j, drive.upstream_pressure_pa * (through * dt));
        assert_eq!(frame.body_work_j, frame.aperture.bore_pressure_pa * (drive.body_flow_m3_s * dt));
        let pressure_work = tract.pressure_drop_pa * (through * dt);
        assert!((frame.aperture.pressure_work_j - pressure_work).abs()
            <= 1e-10 * (frame.aperture.pressure_work_j.abs() + pressure_work.abs()).max(1e-24));
        assert_eq!(frame.aperture.step, model.load().state().step_count());
        assert_eq!(frame.aperture.step, model.upstream_load().unwrap().state().step_count());
        assert_eq!(frame.aperture.bore_pressure_pa, frame.duct.port_outputs[0].conjugate_rate);
        nonzero_tract |= tract.duct.ledger.stored_energy_j > 0.0;
        reverse_flow |= through < 0.0;
    }
    assert!(nonzero_tract && reverse_flow);
    assert!(model.upstream_load().unwrap().terminal().unwrap()
        .pole_states().iter().any(|x| x.norm_sq() > 0.0));
}

#[test]
fn tract_geometry_changes_reed_motion_and_bore_without_changing_source_controls() {
    let mut first = system(0.035);
    let mut second = system(0.07);
    let mut pressure_difference = 0.0;
    let mut motion_difference = 0.0;
    for step in 0..128 {
        let drive = TubeDrive { upstream_pressure_pa: if step < 64 { 0.5 } else { 0.0 },
            body_flow_m3_s: 0.0 };
        let a = first.step(drive).unwrap();
        let b = second.step(drive).unwrap();
        pressure_difference += (a.aperture.bore_pressure_pa - b.aperture.bore_pressure_pa).abs();
        motion_difference += (single(&first).state().opening_m - single(&second).state().opening_m).abs();
    }
    assert!(pressure_difference > 1e-10 && motion_difference > 1e-16);
}

#[test]
fn refused_step_preserves_both_dg_terminal_histories_and_the_original_mechanics() {
    let mut model = system(0.035);
    let mut replay = system(0.035);
    let drive = TubeDrive { upstream_pressure_pa: 0.5, body_flow_m3_s: 2e-10 };
    for _ in 0..47 { assert_eq!(model.step(drive).unwrap(), replay.step(drive).unwrap()); }
    let before = [model.upstream_load().unwrap().state().clone(), model.load().state().clone()];
    let poles = [model.upstream_load().unwrap().terminal().unwrap().pole_states().to_vec(),
        model.load().terminal().unwrap().pole_states().to_vec()];
    let mechanics = single(&model).state();
    for (bad_drive, forces) in [
        (TubeDrive { upstream_pressure_pa: f64::NAN, ..drive }, [0.0; 2]),
        (TubeDrive { body_flow_m3_s: f64::INFINITY, ..drive }, [0.0; 2]),
        (drive, [0.0, 1.0]),
    ] {
        assert!(model.step_with_forces(bad_drive, forces).is_err());
        assert_eq!(model.upstream_load().unwrap().state(), &before[0]);
        assert_eq!(model.load().state(), &before[1]);
        assert_eq!(model.upstream_load().unwrap().terminal().unwrap().pole_states(), &poles[0]);
        assert_eq!(model.load().terminal().unwrap().pole_states(), &poles[1]);
        assert_eq!(single(&model).state(), mechanics);
        assert_eq!(single(&model).accepted_steps(), 47);
    }
    assert_eq!(model.step(drive).unwrap(), replay.step(drive).unwrap());
}

#[test]
fn cold_binding_refuses_wrong_summed_impedance_gas_clock_and_advanced_tract() {
    let up = load(0.035, 101325.0, 1.0 / 96000.0);
    let down = load(0.05, 101325.0, 1.0 / 96000.0);
    let (_, z) = series_inlet_state(&up, &down).unwrap();
    assert!(LoadedApertureDg::new_loaded_with_upstream(
        aperture(&down, down.impedance_pa_s_m3(), BladePropagator::ExactZoh),
        up.clone(), down.clone(),
    ).is_err());
    assert!(LoadedApertureDg::new_loaded_with_upstream(
        aperture(&down, z, BladePropagator::Midpoint), up.clone(), down.clone(),
    ).is_err());
    for other in [
        load(0.035, 81060.0, 1.0 / 96000.0),
        load(0.035, 101325.0, 1.0 / 192000.0),
    ] { assert!(series_inlet_state(&other, &down).is_err()); }
    let mut advanced = up;
    advanced.preview_step(1e-10).unwrap().accept();
    assert!(LoadedApertureDg::new_loaded_with_upstream(
        aperture(&down, z, BladePropagator::ExactZoh), advanced, down,
    ).is_err());
}
