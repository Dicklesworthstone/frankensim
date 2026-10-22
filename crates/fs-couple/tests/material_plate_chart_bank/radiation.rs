//! Synthetic positive-real poles exercise reciprocal composition; the Rayleigh
//! producer has separate geometry/material/gas tests, not fitted-timbre claims.
use super::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;

fn load() -> Model {
    Model { ports: 3, poles: vec![Pole { omega: 700.0, zeta: 0.1,
        coupling: vec![120.0,-80.0,60.0] }] }
}
fn base(cx: &Cx<'_>, gate: &CancelGate) -> ReciprocalStringBank {
    let mut settings=options();
    settings.network.max_modes=64;
    settings.network.max_connections=64;
    settings.network.max_setup_terms=1_000_000;
    let gas=gas(101_325.0);
    ReciprocalStringBank::from_material_plate_chart(cx,sources(),gas,body(),Some(cavity(gas)),
        law(),settings,gate).unwrap()
}
fn attach(bank: ReciprocalStringBank, gate: &CancelGate) -> ReciprocalStringBank {
    bank.with_body_radiation_load(load(),1,gate).unwrap()
}
fn extend(bank: ReciprocalStringBank, axial: bool, gate: &CancelGate) -> ReciprocalStringBank {
    if axial {
        bank.with_longitudinal_waves(vec![LongitudinalStringResolution {
            axial_modes: vec![1], maximum_slope: 0.2, maximum_strain: 0.01,
        };2],NonlinearModalConfig { max_iterations:64, maximum_generalized_force:1e5,
            force_absolute_tolerance:1e-11, force_relative_tolerance:1e-10 },gate).unwrap()
    } else {
        bank.with_cylinder_gas_loading(CylinderFlowResolution {
            minimum_wavenumber_radius:1e-3, maximum_wavenumber_radius:30.0, intervals:2,
        },gate).unwrap()
    }
}

#[test]
fn physical_body_radiation_changes_string_motion_and_pressure_in_the_same_work_balance() {
    with_cx(|cx,gate|{
        let mut unloaded=base(cx,gate);
        let mut loaded=attach(base(cx,gate),gate);
        let original=loaded.source().material_chart().unwrap().specimen.specimen_identity();
        assert_eq!(loaded.system().mode_count(),unloaded.system().mode_count());
        assert_eq!(loaded.source().body_radiation().unwrap().gas(),loaded.source().gas());
        let mut force=vec![0.0;loaded.system().mode_count()];
        force[7]=0.02;
        let released=FrictionDrive{speed_m_s:0.0,normal_force_n:0.0};
        let mut pressure_difference=0.0;
        for _ in 0..32 {
            let a=*unloaded.step(&force,released).unwrap();
            let b=*loaded.step(&force,released).unwrap();
            assert!(b.energy_residual_j.abs()<=b.energy_tolerance_j);
            pressure_difference+=(a.observer_pressure_pa-b.observer_pressure_pa).abs();
        }
        assert!(pressure_difference>1e-12);
        assert_ne!(loaded.system().components()[0].states(),unloaded.system().components()[0].states());
        assert!(loaded.system().radiation_energy_j().unwrap()>0.0);
        assert!(loaded.system().coupled_frame().unwrap().radiation_dissipation_j>0.0);
        assert_eq!(loaded.source().material_chart().unwrap().specimen.specimen_identity(),original);
    });
}

#[test]
fn body_basis_and_load_survive_both_orders_of_gas_and_axial_extension() {
    with_cx(|cx,gate|{
        for axial in [false,true] {
            let mut first=extend(attach(base(cx,gate),gate),axial,gate);
            let mut last=attach(extend(base(cx,gate),axial,gate),gate);
            assert_eq!(first.source().body_radiation().unwrap().basis_identity(),
                last.source().body_radiation().unwrap().basis_identity());
            assert_eq!(first.system().mode_count(),last.system().mode_count());
            assert_eq!(first.system().radiation_state().len(),2);
            let zero=vec![0.0;first.system().mode_count()];
            for _ in 0..8 {
                assert_eq!(first.step(&zero,drive()).unwrap(),last.step(&zero,drive()).unwrap());
                assert_eq!(first.system().radiation_state(),last.system().radiation_state());
            }
            assert!(first.system().radiation_state().iter().any(|v|*v!=0.0));
        }
    });
}

#[test]
fn a_vibrating_loaded_body_continues_through_multiple_bows_without_resetting_air() {
    with_cx(|cx,gate|{
        let mut reference=attach(base(cx,gate),gate);
        let mut bank=attach(base(cx,gate),gate);
        let zero=vec![0.0;bank.system().mode_count()];
        for _ in 0..4 {
            assert_eq!(bank.step(&zero,drive()).unwrap(),reference.step(&zero,drive()).unwrap());
        }
        let identity=bank.source().body_radiation().unwrap().basis_identity();
        let mut multi=bank.with_additional_bows(vec![],MultiFrictionConfig {
            max_ports:4,max_sweeps:128,max_setup_terms:16_384,
        },gate).unwrap();
        assert_eq!(multi.source().bank().body_radiation().unwrap().basis_identity(),identity);
        for _ in 0..8 {
            let a=*reference.step(&zero,drive()).unwrap();
            let b=multi.step(&zero,&[drive()]).unwrap();
            assert_eq!(a.observer_pressure_pa,b.observer_pressure_pa);
            assert_eq!(a.network_energy_j,b.network_energy_j);
            assert_eq!(a.network_dissipation_j,b.network_dissipation_j);
            assert_eq!(a.energy_residual_j,b.energy_residual_j);
        }
    });
}

#[test]
fn scheduled_pressure_is_the_same_loaded_physics_under_different_audio_blocks() {
    with_cx(|cx,gate|{
        let schedule=GestureSchedule::try_new(700,vec![GestureTrack { id:"bow".into(),
            target:GestureTarget::BowStroke{string:0},
            initial:GestureValue::Bow{velocity_m_per_s:0.01,normal_force_n:0.1,station:0.2},events:vec![],
        }]).unwrap();
        let config=FrictionGestureConfig {sample_rate_hz:48_000,samples:17,max_block:16,
            max_work:1000,max_events:16,station_fraction:0.2};
        let mut direct=attach(base(cx,gate),gate);
        let zero=vec![0.0;direct.system().mode_count()];
        let expected:Vec<_>=(0..17).map(|_|direct.step(&zero,drive()).unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1,7,16] {
            let mut renderer=attach(base(cx,gate),gate).into_scheduled(&schedule,"bow",config).unwrap();
            assert!(renderer.source().body_radiation().is_some());
            let mut output=vec![0.0;17];
            for chunk in output.chunks_mut(partition) {renderer.renderer_mut().block(chunk).unwrap();}
            assert_eq!(output.into_iter().map(f64::to_bits).collect::<Vec<_>>(),expected);
        }
    });
}
