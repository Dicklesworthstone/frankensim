//! Frozen regional acoustics; neither thermal mixing nor measured instrument data.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, ApertureNetworkFrame, NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
const RATE: u32 = 96_000;
const DT: f64 = 1.0 / RATE as f64;
fn gas(t: f64, rh: f64) -> GasState { GasState::try_new_moist_air(t, 101_325.0, rh).unwrap() }
fn gases() -> Vec<GasState> { vec![gas(293.15, 0.0), gas(313.15, 0.5), gas(283.15, 0.0)] }
fn graph() -> TubeNetworkSpec {
    TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
            NetworkNode::Termination {reflection: -0.8}, NetworkNode::Termination {reflection: 0.2}],
        sections: [(0,1,0.12), (1,2,0.18), (1,3,0.09)].iter().map(|&(a,b,l)|
            TubeSection {nodes:[a,b], length_m:l, radius_m:0.007, max_length_error_m:0.002}).collect(),
        sound_speed_m_s:gas(293.15,0.0).sound_speed, max_wave_memory_bytes:1<<20,
    }
}
fn valve(density: f64, z: f64) -> DynamicAperture {
    DynamicAperture::new(DynamicApertureSpec {
        aperture: BernoulliAperture {rest_opening_m:0.0004,width_m:0.013,closing_pressure_pa:6000.0},
        mass_kg:1e-5,stiffness_n_m:500.0,damping_ratio:0.35,density_kg_m3:density,
        impedance_pa_s_m3:z,time_step_s:DT,max_steps:2048,
    }, ApertureState {opening_m:0.0004,opening_velocity_m_s:0.0},
        Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
            "synthetic regional-gas contact".into()).unwrap()).unwrap()
}
fn model(gases: Vec<GasState>) -> ApertureNetwork {
    let s=graph(); let z=s.inlet_impedance_with_gases(&gases).unwrap();
    ApertureNetwork::with_section_gases(valve(gases[0].density,z),s,gases).unwrap()
}
fn drive(n: usize) -> TubeDrive {
    TubeDrive {upstream_pressure_pa:if n<256 {100.0} else {0.0},body_flow_m3_s:0.0}
}

#[test]
fn regional_state_drives_the_original_scattering_and_complete_energy_balance() {
    use fs_vfit::waveguide::network::{NetworkSegment,WaveguideNetwork};
    let s=graph();let gas=gases();let mut actual=model(gas.clone());
    let segments:Vec<_>=s.sections.iter().zip(&gas).map(|(s,g)|NetworkSegment {
        nodes:s.nodes,one_way_samples:(s.length_m/(g.sound_speed*DT)).round() as usize,
        impedance_pa_s_m3:g.density*g.sound_speed/(core::f64::consts::PI*s.radius_m*s.radius_m),
    }).collect();
    let mut raw=WaveguideNetwork::new(&s.nodes,&segments,DT,1<<20).unwrap();
    let mut reed=valve(gas[0].density,segments[0].impedance_pa_s_m3);
    for n in 0..1024 {
        let before=actual.stored_energy_j();let d=drive(n);
        let r=reed.step(ApertureDrive {upstream_pressure_pa:d.upstream_pressure_pa,
            incoming_pressure_pa:raw.incoming_pressure_pa(),body_flow_m3_s:d.body_flow_m3_s}).unwrap();
        let wave=raw.step(r.outgoing_pressure_pa).unwrap();
        let f=actual.step(d).unwrap();assert_eq!(f.aperture,r);assert_eq!(f.network,wave);
        for node in 0..s.nodes.len() {assert_eq!(actual.node_frame(node),raw.node_frame(node));}
        let scale=before+f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs();
        assert!(f.balance_residual_j().abs()<=3e-10*scale.max(1e-30));
    }
    assert_eq!(actual.section_gases(),Some(gas.as_slice()));
    for (i,r) in actual.represented_sections().iter().enumerate() {
        assert_eq!(r.one_way_samples,segments[i].one_way_samples);
        assert_eq!(r.impedance_pa_s_m3,segments[i].impedance_pa_s_m3);
    }
}

#[test]
fn uniform_regional_path_is_bitwise_legacy_and_remote_gas_feedback_is_causal() {
    let g=gas(293.15,0.0);let s=graph();
    let mut legacy=ApertureNetwork::new(valve(g.density,s.inlet_impedance(g.density).unwrap()),s).unwrap();
    let mut uniform=model(vec![g;3]);let mut changed=model(gases());let mut differs=false;
    let first_return=2*uniform.represented_sections()[0].one_way_samples;
    for n in 0..1024 {
        let a=legacy.step(drive(n)).unwrap();let b=uniform.step(drive(n)).unwrap();
        assert_eq!(a,b);
        let c=changed.step(drive(n)).unwrap();
        if n<first_return {assert_eq!(b.aperture,c.aperture,"no remote feedthrough");}
        differs|=(b.aperture.state.opening_m-c.aperture.state.opening_m).abs()>1e-12;
    }
    assert!(differs,"remote temperature/density must change returning pressure and motion");
    assert_ne!(uniform.represented_sections()[1].one_way_samples,changed.represented_sections()[1].one_way_samples);
}

#[test]
fn cancellation_and_refused_input_preserve_regional_wave_history() {
    let mut a=model(gases());let mut b=model(gases());
    for n in 0..100 {assert_eq!(a.step(drive(n)).unwrap(),b.step(drive(n)).unwrap());}
    let state=b.aperture().state();let energy=b.stored_energy_j().to_bits();
    assert!(b.step(TubeDrive {upstream_pressure_pa:f64::NAN,body_flow_m3_s:0.0}).is_err());
    let stop=CancelGate::new();stop.request();
    let sentinel=ApertureNetworkFrame {stored_energy_j:-123.0,..ApertureNetworkFrame::default()};
    let mut out=[sentinel;37];
    assert_eq!(b.advance_block(&[drive(100);37],&mut out,&stop).unwrap().completed,0);
    assert_eq!(out,[sentinel;37]);assert_eq!(b.aperture().state(),state);assert_eq!(b.stored_energy_j().to_bits(),energy);
    for n in 100..1024 {assert_eq!(a.step(drive(n)).unwrap(),b.step(drive(n)).unwrap());}
}

#[test]
fn invalid_regional_maps_cannot_silently_select_uniform_air_or_steady_flow() {
    let s=graph();let g=gases();let z=s.inlet_impedance_with_gases(&g).unwrap();
    assert!(s.inlet_impedance_with_gases(&g[..2]).is_err());
    let mut bad=g.clone();bad[1].pressure=90_000.0;assert!(s.inlet_impedance_with_gases(&bad).is_err());
    let mut bad=g.clone();bad[1].sound_speed=f64::NAN;assert!(s.inlet_impedance_with_gases(&bad).is_err());
    let mut wrong=s.clone();wrong.sound_speed_m_s=343.0;assert!(wrong.validate_section_gases(&g).is_err());
    assert!(ApertureNetwork::with_section_gases(valve(g[0].density*2.0,z),s.clone(),g.clone()).is_err());
    let mut small=s;small.max_wave_memory_bytes=1;assert!(ApertureNetwork::with_section_gases(valve(g[0].density,z),small,g).is_err());
}

#[test]
fn exterior_receiver_uses_the_observed_outlet_medium_not_the_inlet() {
    use fs_couple::bernoulli_aperture::performance::{AperturePerformance,AperturePerformanceConfig,ApertureObservation,CoupledAperture};
    use fs_couple::pcm_wav::baffled::{BaffledPressure,CircularOutletReceiver,RayleighMedium};
    use fs_couple::pcm_wav::observation::PressureRenderer;
    use fs_scenario::gesture::{GestureSchedule,GestureTrack,GestureTarget,GestureValue};
    let loc=CircularOutletReceiver {position_m:[0.0,0.0,0.2],radial_rings:8,angular_points:32,maximum_frequency_hz:1000.0};
    let states=gases();let outlet=states[1];
    let mut receiver=BaffledPressure::circular_outlet(0.007,RATE,loc,
        RayleighMedium {density:outlet.density,sound_speed:outlet.sound_speed}).unwrap();
    let schedule=GestureSchedule::try_new(1,vec![GestureTrack {id:"mouth".into(),target:GestureTarget::BlowingPressure,
        initial:GestureValue::PressurePa(100.0),events:vec![]}]).unwrap();
    let mut r=AperturePerformance::new(CoupledAperture::Network(model(states.clone())),
        ApertureObservation::NetworkBaffled {node:2,receiver:loc},schedule,
        AperturePerformanceConfig {sample_rate_hz:RATE,samples:1024,max_block:37,max_compile_work:128,max_controls:2}).unwrap();
    let mut direct=model(states);let mut pressure=vec![0.0;1024];
    for chunk in pressure.chunks_mut(37) {r.block(chunk).unwrap();}
    let mut audible=false;
    for actual in pressure {
        direct.step(TubeDrive {upstream_pressure_pa:100.0,body_flow_m3_s:0.0}).unwrap();
        let expected=receiver.step(&[direct.node_frame(2).unwrap().net_flow_into_node_m3_s]).unwrap();
        assert_eq!(actual.to_bits(),expected.to_bits());audible|=expected.abs()>1e-12;
    }
    assert!(audible);assert_eq!(r.system().aperture().state(),direct.aperture().state());
}
