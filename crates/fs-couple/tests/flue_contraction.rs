use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_couple::bernoulli_aperture::contraction::ContractionSpec;
use fs_couple::bernoulli_aperture::network::{NetworkNode,TubeNetworkSpec,TubeSection};
use fs_couple::flue_loop::{FlueControl,FlueGeometry,FlueVoice};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec,GasState};
const RATE:u32=48000;
fn spec(ratio:f64)->ContractionSpec {
    ContractionSpec {area_ratio:ratio,maximum_contracted_mach:0.2,
        provenance:"authored fixed-window low-Mach comparison".into()}
}
fn voice(pressure:f64)->FlueVoice {
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,pressure).unwrap();
    let geometry=FlueGeometry {flue_width_m:0.01,flue_height_m:0.001,cut_up_m:0.004,
        labium_offset_m:0.0,jet_angle_rad:0.0,
        bore_area_m2:core::f64::consts::PI*0.006*0.006,mouth_loss_rel:0.02};
    let source=FlueNetworkSpec {
        graph:TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,NetworkNode::Termination {reflection:-0.9}],
            sections:vec![TubeSection {nodes:[0,1],length_m:0.3,radius_m:0.006,max_length_error_m:0.005}],
            sound_speed_m_s:gas.sound_speed,max_wave_memory_bytes:1_000_000},
        gases:vec![gas],viscothermal:vec![],minimum_jet_speed_m_s:1.0,maximum_jet_mach:0.3,
        max_jet_history_samples:4096};
    FlueVoice::from_network(Some(&mint_tonal_interim_card().unwrap()),90.0,1e-6,
        geometry,source,RATE,3,0,&CancelGate::new_clock_free()).unwrap()
}

#[test]
fn separation_enters_the_solved_mouth_and_complete_network_work_account() {
    let mut v=voice(101325.0).with_mouth_contraction(spec(0.6)).unwrap();
    let mut linear=voice(101325.0);
    for p in [&mut v,&mut linear] {p.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();}
    let mut different=false;let mut dissipated=0.0;
    for n in 0..1024 {
        if n==256 {for p in [&mut v,&mut linear] {p.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();}}
        let mut a=[0.0];let mut b=[0.0];v.step_block(&mut a).unwrap();linear.step_block(&mut b).unwrap();
        different|=a!=b;
        let d=v.diagnostics().last().unwrap();let network=v.network().unwrap().last_frame().unwrap();
        let loss=v.mouth_contraction().unwrap().at_flow(v.mouth_flow_m3_s()).unwrap();
        assert_eq!(d.separation_loss_j,loss.dissipated_power_w*(1.0/f64::from(RATE)));
        assert!(d.separation_loss_j>=0.0);
        assert!((network.inlet_work_j-d.pipe_work_j).abs()<1e-12);
        let defect=d.source_work_j-d.stored_delta_j-d.numerical_dissipation_j
            -d.radiation_loss_j-d.separation_loss_j-network.storage_change_j
            -network.terminal_loss_j-network.interior_loss_j;
        assert!(defect.abs()<1e-12);
        assert!(d.ledger_defect_j.abs()<1e-12);
        assert_eq!(linear.diagnostics().last().unwrap().separation_loss_j,0.0);
        dissipated+=d.separation_loss_j;
    }
    assert!(different && dissipated>0.0);
    assert!(v.network().unwrap().stored_energy_j()>0.0,"release must retain traveling sound");
    assert!(v.diagnostics().last().unwrap().separation_loss_j>0.0,"acoustic release still loses energy");
}

#[test]
fn fixed_data_and_sample_controls_are_invariant_to_callback_partition() {
    let mut one=voice(101325.0).with_mouth_contraction(spec(0.6)).unwrap();
    let mut many=voice(101325.0).with_mouth_contraction(spec(0.6)).unwrap();
    for v in [&mut one,&mut many] {v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();}
    let mut a=[0.0;513];let mut b=[0.0;513];
    one.step_block(&mut a).unwrap();for chunk in b.chunks_mut(37) {many.step_block(chunk).unwrap();}
    assert_eq!(a,b);assert_eq!(one.mouth_flow_m3_s(),many.mouth_flow_m3_s());
    assert_eq!(one.network().unwrap().last_frame(),many.network().unwrap().last_frame());
    let x=one.diagnostics().iter().map(|d|d.separation_loss_j).sum::<f64>();
    let y=many.diagnostics().iter().map(|d|d.separation_loss_j).sum::<f64>();
    assert!((x-y).abs()<1e-12*x);
}

#[test]
fn actual_material_gas_and_window_define_the_contraction_not_the_bore_or_jet() {
    let a=voice(101325.0).with_mouth_contraction(spec(0.6)).unwrap();
    let b=voice(202650.0).with_mouth_contraction(spec(0.6)).unwrap();
    let la=a.mouth_contraction().unwrap();let lb=b.mouth_contraction().unwrap();
    assert_eq!(la.gas(),a.gas());assert_eq!(lb.gas(),b.gas());
    assert_eq!(la.width_m(),0.01);assert_eq!(la.opening_m(),0.004);
    let pa=la.at_flow(1e-5).unwrap().pressure_drop_pa;
    let pb=lb.at_flow(1e-5).unwrap().pressure_drop_pa;
    assert!((pb/pa-2.0).abs()<1e-12);
    assert_ne!(la.effective_area_m2(),core::f64::consts::PI*0.006*0.006);
}

#[test]
fn late_binding_and_solved_domain_excursions_do_not_become_clipped_success() {
    assert!(voice(101325.0).with_mouth_contraction(spec(0.0)).is_err());
    let v=voice(101325.0).with_mouth_contraction(spec(0.6)).unwrap();
    assert!(v.with_mouth_contraction(spec(0.5)).is_err());
    let mut v=voice(101325.0);v.step_block(&mut [0.0]).unwrap();
    assert!(v.with_mouth_contraction(spec(0.6)).is_err());
    let mut narrow=spec(0.6);narrow.maximum_contracted_mach=1e-12;
    let mut v=voice(101325.0).with_mouth_contraction(narrow).unwrap();
    v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let mut sentinel=[123.0];assert!(v.step_block(&mut sentinel).is_err());
    assert_eq!(sentinel,[123.0]);assert_eq!(v.samples_rendered(),0);
    assert_eq!(v.network().unwrap().stored_energy_j(),0.0);
    assert!(v.step_block(&mut sentinel).is_err());
}

#[test]
fn physical_receiver_observes_the_solved_loss_not_a_postprocessed_pressure_gain() {
    use fs_couple::flue_loop::network::performance::{FlueOutletConfig,FlueOutletPerformance};
    use fs_couple::pcm_wav::baffled::{BaffledPressure,CircularOutletReceiver,RayleighMedium};
    use fs_couple::pcm_wav::baffled::openings::CircularOpening;
    use fs_couple::pcm_wav::observation::PressureRenderer;
    use fs_scenario::gesture::GestureSchedule;
    let schedule=GestureSchedule::from_canonical_bytes(b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\nevent\t0e0\t0e0\tpressure-pa\t3e2\nevent\t5e-3\t0e0\tpressure-pa\t0e0\n").unwrap();
    let receiver=CircularOutletReceiver {position_m:[0.0,0.0,0.25],radial_rings:4,
        angular_points:16,maximum_frequency_hz:2000.0};
    let config=FlueOutletConfig {sample_rate_hz:RATE,samples:513,max_block:64,
        max_controls:128,max_compile_work:100000,outlet_node:1,receiver,maximum_abs_pressure_pa:1000.0};
    let mut actual=FlueOutletPerformance::new(voice(101325.0).with_mouth_contraction(spec(0.6)).unwrap(),
        schedule.clone(),config).unwrap();
    let mut direct=voice(101325.0).with_mouth_contraction(spec(0.6)).unwrap();
    let mut observed=BaffledPressure::circular_openings(&[CircularOpening {center_m:[0.0,0.0],radius_m:0.006}],
        RATE,receiver,RayleighMedium {density:direct.gas().density,sound_speed:direct.gas().sound_speed}).unwrap();
    let mut expected=vec![0.0;513];
    for (n,slot) in expected.iter_mut().enumerate() {
        if n==0 {direct.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();}
        if n==240 {direct.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();}
        direct.step_block(&mut [0.0]).unwrap();
        let q=direct.network().unwrap().node_frame(1).unwrap().net_flow_into_node_m3_s;
        *slot=observed.step(&[q]).unwrap();
    }
    let mut out=vec![0.0;513];for block in out.chunks_mut(37) {actual.block(block).unwrap();}
    assert_eq!(out,expected);assert!(out.iter().any(|p|*p!=0.0));
    assert_eq!(actual.voice().mouth_flow_m3_s(),direct.mouth_flow_m3_s());
    assert!(actual.voice().mouth_contraction().is_some());
    let mut other=FlueOutletPerformance::new(voice(101325.0).with_mouth_contraction(spec(1.0)).unwrap(),schedule,config).unwrap();
    let mut changed=vec![0.0;513];for block in changed.chunks_mut(64) {other.block(block).unwrap();}
    assert_ne!(out,changed,"same receiver; changed physical loss must reach the microphone");
}
