//! File selection must install the physical load, not merely name a microphone.
use super::*;
use fs_couple::bernoulli_aperture::network::{ApertureNetwork,NetworkNode,TubeNetworkSpec,TubeSection};
use fs_couple::bernoulli_aperture::performance::{ApertureObservation,CoupledAperture};
use fs_couple::pcm_wav::baffled::{BaffledPressure,CircularOutletReceiver,RayleighMedium};
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{RelaxationImpedanceSpec,RelaxationTerm};
const RADIATING:&str=include_str!("../../examples/plate-valve-radiating.performance");
const OBSERVER:&str="observation baffled-outlet 0 0 0.2 8 32 1500";

#[test]
fn file_load_and_receiver_match_independently_composed_reciprocal_physics() {
    let seed=load(INPUT);
    assert!(seed.info().radiation_load.is_none());
    let a=seed.renderer().system().aperture();let s=*a.spec();
    // Reuse the already-verified old specimen path, not the new load parser.
    // Supply the physical relaxation coefficients independently of the helper.
    let valve=DynamicAperture::from_plate_with_closure(
        a.plate_reduction().unwrap().clone(),a.plate_closure().unwrap().spec().clone(),
        s.density_kg_m3,s.impedance_pa_s_m3,s.time_step_s,s.max_steps,a.state()).unwrap()
        .with_plate_relaxation(a.relaxation().unwrap().spec().clone(),InitialApertureMemory::Relaxed).unwrap();
    let CoupledAperture::Tube(tube)=seed.renderer().system() else {panic!("old source path changed")};
    let t=*tube.spec();let alpha=8.0/(3.0*core::f64::consts::PI);
    let load_spec=RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
        resistance_pa_s_m3:0.0,inertance_pa_s2_m3:0.0,compliance_m3_pa:None,
    },&[RelaxationTerm {resistance_pa_s_m3:2.0*alpha*alpha*s.impedance_pa_s_m3,
        rate_per_s:2.0*alpha*t.sound_speed_m_s/t.radius_m}]).unwrap();
    let mut direct=ApertureNetwork::new(valve,TubeNetworkSpec {
        nodes:vec![NetworkNode::Inlet,NetworkNode::Relaxation {load:load_spec}],
        sections:vec![TubeSection {nodes:[0,1],length_m:t.length_m,radius_m:t.radius_m,max_length_error_m:t.max_length_error_m}],
        sound_speed_m_s:t.sound_speed_m_s,max_wave_memory_bytes:t.max_wave_memory_bytes,
    }).unwrap();
    let source=load(RADIATING);let info=source.info();
    assert_eq!(info.radiation_load.unwrap().load(),load_spec);
    assert_eq!(info.represented_tube_length_m,tube.represented_length_m());
    assert_eq!(info.memory_branches,seed.info().memory_branches);
    let mut actual=source.into_renderer();
    let mut terminal=load(&RADIATING.replace(OBSERVER,"observation terminal")).into_renderer();
    let receiver=CircularOutletReceiver {position_m:[0.0,0.0,0.2],radial_rings:8,
        angular_points:32,maximum_frequency_hz:1500.0};
    let mut mic=BaffledPressure::circular_outlet(t.radius_m,96000,receiver,
        RayleighMedium {density:s.density_kg_m3,sound_speed:t.sound_speed_m_s}).unwrap();
    let phrase=seed.renderer().schedule();let mut nonzero=false;
    for n in 0..info.samples {
        direct.step(TubeDrive {upstream_pressure_pa:phrase.sample("mouth",n*700/96000).unwrap(),body_flow_m3_s:0.0}).unwrap();
        let node=direct.node_frame(1).unwrap();
        let expected=mic.step(&[node.net_flow_into_node_m3_s]).unwrap();
        let mut value=[0.0];actual.block(&mut value).unwrap();
        assert_eq!(value[0].to_bits(),expected.to_bits());nonzero|=expected.abs()>1e-10;
        terminal.block(&mut value).unwrap();assert_eq!(value[0].to_bits(),node.pressure_pa.to_bits());
    }
    let CoupledAperture::Network(retained)=actual.system() else {panic!("radiation selection lost its load")};
    assert_eq!(retained.aperture().state(),direct.aperture().state());
    assert_eq!(retained.aperture().relaxation().unwrap().memory_sqrt_j(),direct.aperture().relaxation().unwrap().memory_sqrt_j());
    assert_eq!(retained.stored_energy_j().to_bits(),direct.stored_energy_j().to_bits());
    assert!(matches!(terminal.observation(),ApertureObservation::NetworkNode(1)));
    assert!(nonzero);
}

#[test]
fn receiver_is_observation_but_outlet_geometry_is_an_actual_mechanical_load() {
    let near=load(RADIATING);let original_load=near.info().radiation_load.unwrap();
    let far=load(&RADIATING.replace("baffled-outlet 0 0 0.2","baffled-outlet 0.1 0 0.4"));
    assert_eq!(far.info().radiation_load,Some(original_load));
    let wider=load(&RADIATING.replace("tube 0.25 0.007","tube 0.25 0.008"));
    assert_ne!(wider.info().radiation_load.unwrap().load(),original_load.load());
    let(mut near,mut far,mut wider)=(near.into_renderer(),far.into_renderer(),wider.into_renderer());
    assert!(far.baffled_receiver().unwrap().delay_samples.0>near.baffled_receiver().unwrap().delay_samples.1);
    let(mut sound_changed,mut motion_changed)=(false,false);
    for _ in 0..64 {
        let(mut a,mut b,mut c)=([0.0;37],[0.0;37],[0.0;37]);
        near.block(&mut a).unwrap();far.block(&mut b).unwrap();wider.block(&mut c).unwrap();
        sound_changed|=a!=b;
        assert_eq!(near.system().aperture().state(),far.system().aperture().state());
        assert_eq!(near.system().aperture().relaxation().unwrap().memory_sqrt_j(),far.system().aperture().relaxation().unwrap().memory_sqrt_j());
        motion_changed|=near.system().aperture().state()!=wider.system().aperture().state();
    }
    assert!(sound_changed && motion_changed);
}

#[test]
fn unsupported_or_contradictory_radiation_records_cannot_select_a_fallback() {
    for token in ["baffled-low-ka", "baffled-low-ka 0", "baffled-low-ka -1",
        "baffled-low-ka NaN", "baffled-low-ka 10000", "baffled-low-ka 1500 -0.8", "guessed 1500"] {
        let source=RADIATING.replace("baffled-low-ka 1500",token);
        assert!(PlateValvePerformance::from_bytes(source.as_bytes(),37,&CancelGate::new()).is_err(),"accepted {token}");
    }
    let broader=RADIATING.replace("8 32 1500","8 32 1600");
    assert!(PlateValvePerformance::from_bytes(broader.as_bytes(),37,&CancelGate::new()).is_err());
    let tiny_budget=RADIATING.replace("0.002 1048576","0.002 8");
    assert!(PlateValvePerformance::from_bytes(tiny_budget.as_bytes(),37,&CancelGate::new()).is_err());
    // The original explicitly numeric terminal remains the original tube path.
    let unchanged=load(INPUT);
    assert!(matches!(unchanged.renderer().system(),CoupledAperture::Tube(_)));
    assert!(unchanged.info().radiation_load.is_none());
}
