use super::*;
use fs_couple::bernoulli_aperture::network::{ApertureNetwork,NetworkNode,TubeNetworkSpec,TubeSection};
use fs_couple::bernoulli_aperture::viscothermal::{with_viscothermal_sections,ViscothermalSectionSpec};
use fs_couple::bernoulli_aperture::performance::{CoupledAperture,AperturePerformance};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::pcm_wav::baffled::{BaffledPressure,CircularOutletReceiver,RayleighMedium};
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::pcm_wav::stream::{Pcm16WavStream,ScheduledWavProgress,render_pressure_pcm16};
use std::io::Cursor;
const LOSSY:&str=include_str!("../../examples/plate-valve-viscothermal.performance");
fn network(p:&AperturePerformance)->&ApertureNetwork {
    let CoupledAperture::Network(n)=p.system() else {panic!("lossy sections were replaced by a bare tube")};n
}
fn fresh(seed:&PlateValvePerformance,z:f64)->DynamicAperture {
    let a=seed.renderer().system().aperture();let s=a.spec();
    DynamicAperture::from_plate_with_closure(a.plate_reduction().unwrap().clone(),
        a.plate_closure().unwrap().spec().clone(),s.density_kg_m3,z,s.time_step_s,s.max_steps,a.state()).unwrap()
        .with_plate_relaxation(a.relaxation().unwrap().spec().clone(),InitialApertureMemory::Relaxed).unwrap()
}

#[test]
fn supplied_gas_losses_reach_the_reciprocal_valve_and_actual_outlet_flow() {
    let seed=load(INPUT);let p=load(LOSSY);let gas=GasState::try_new_moist_air(293.15,101325.0,0.0).unwrap();
    let dt=1.0/96000.0;
    let radiation=BaffledRadiationLoad::new(0.007,gas.density,gas.sound_speed,dt,1000.0,&CancelGate::new()).unwrap();
    let original=TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,radiation.termination()],
        sections:vec![TubeSection {nodes:[0,1],length_m:0.25,radius_m:0.007,max_length_error_m:0.002}],
        sound_speed_m_s:gas.sound_speed,max_wave_memory_bytes:1<<20};
    let options=ViscothermalSectionSpec {section:0,minimum_frequency_hz:100.0,maximum_frequency_hz:1000.0,cells:12};
    let (graph,reports)=with_viscothermal_sections(original,&gas,dt,&[options],&CancelGate::new()).unwrap();
    assert_eq!(p.viscothermal_losses().len(),reports.len());
    assert_eq!(p.viscothermal_losses()[0].source,reports[0].source);
    assert_eq!(p.viscothermal_losses()[0].section_range,reports[0].section_range);
    assert_eq!(p.viscothermal_losses()[0].one_way_samples,reports[0].one_way_samples);
    assert_eq!(network(p.renderer()).spec(),&graph);
    assert!(load(INPUT).viscothermal_losses().is_empty());
    let mut direct=ApertureNetwork::new(fresh(&seed,graph.inlet_impedance(gas.density).unwrap()),graph).unwrap();
    let receiver=CircularOutletReceiver {position_m:[0.0,0.0,0.2],radial_rings:8,angular_points:32,maximum_frequency_hz:1000.0};
    let mut mic=BaffledPressure::circular_outlet(0.007,96000,receiver,
        RayleighMedium {density:gas.density,sound_speed:gas.sound_speed}).unwrap();
    let total=p.info().samples;let mut actual=p.into_renderer();let mut lossy_energy=0.0_f64;
    let mut reference=load(&LOSSY.replace(" viscothermal 100 1000 12 8","")).into_renderer();
    let mut different=false;
    for n in 0..total {
        let before=direct.stored_energy_j();
        let f=direct.step(TubeDrive {upstream_pressure_pa:seed.renderer().schedule().sample("mouth",n*700/96000).unwrap(),body_flow_m3_s:0.0}).unwrap();
        let scale=before+f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs();
        assert!(f.balance_residual_j().abs()<=3e-10*scale.max(f64::MIN_POSITIVE));
        lossy_energy+=f.network.interior_loss_j;
        let expected=mic.step(&[direct.node_frame(1).unwrap().net_flow_into_node_m3_s]).unwrap();
        let (mut out,mut baseline)=([0.0],[0.0]);actual.block(&mut out).unwrap();reference.block(&mut baseline).unwrap();
        assert_eq!(out[0].to_bits(),expected.to_bits());
        different|=actual.system().aperture().state()!=reference.system().aperture().state();
    }
    assert!(different && lossy_energy>0.0,"loss must change source mechanics, not only observation");
    assert_eq!(actual.system().aperture().relaxation().unwrap().memory_sqrt_j(),direct.aperture().relaxation().unwrap().memory_sqrt_j());
}

#[test]
fn selective_branch_loss_preserves_addresses_and_resumes_all_gas_solid_and_receiver_memory() {
    let graph="network 4 3 1048576\nduct_node inlet\nduct_node junction\nduct_node baffled-low-ka 1000\nduct_node cavity 0.00001 0.0015 0.02 200000\n\
        duct_section 0 1 0.125 0.007 0.002 viscothermal 100 1000 6 8\n\
        duct_section 1 2 0.125 0.007 0.002 viscothermal 100 1000 6 8\n\
        duct_section 1 3 0.043 0.002 0.002";
    let text=LOSSY.replace("tube 0.25 0.007 baffled-low-ka 1000 0.002 1048576 viscothermal 100 1000 12 8",graph)
        .replace("observation baffled-outlet","observation network-baffled 2");
    let p=load(&text);assert_eq!(p.viscothermal_losses().len(),2);
    assert_eq!(p.radiation_loads()[0].0,2);
    assert!(matches!(network(p.renderer()).spec().nodes[3],NetworkNode::Impedance{..}));
    let build=||DecimatedRenderer::new(load(&text).into_renderer(),96000,48000,37).unwrap();
    let stream=||Pcm16WavStream::new(Cursor::new(Vec::new()),48000,0.01,37).unwrap();
    let mut a=build();let mut expected=stream();let mut scratch=[0.0;37];
    render_pressure_pcm16(&mut a,&mut expected,&CancelGate::new(),&mut scratch,4801).unwrap();
    let mut b=build();let mut actual=stream();
    render_pressure_pcm16(&mut b,&mut actual,&CancelGate::new(),&mut scratch,1001).unwrap();
    let energy=network(b.source()).stored_energy_j().to_bits();
    let state=b.source().system().aperture().state();let pending=b.source().pending_controls().to_vec();
    let gate=CancelGate::new();gate.request();scratch.fill(-123.0);
    assert_eq!(render_pressure_pcm16(&mut b,&mut actual,&gate,&mut scratch,3800).unwrap(),ScheduledWavProgress::Cancelled{samples:0});
    assert_eq!(scratch,[-123.0;37]);assert_eq!(network(b.source()).stored_energy_j().to_bits(),energy);
    assert_eq!(b.source().system().aperture().state(),state);assert_eq!(b.source().pending_controls(),pending);
    render_pressure_pcm16(&mut b,&mut actual,&CancelGate::new(),&mut scratch,3800).unwrap();
    assert_eq!(actual.finish().unwrap().0.into_inner(),expected.finish().unwrap().0.into_inner());
}

#[test]
fn malformed_or_underresolved_loss_records_refuse_without_changing_other_source_paths() {
    for suffix in ["viscothermal", "viscothermal 100 1000 0 8", "viscothermal 100 1000 1 8",
        "viscothermal 100 1000 256 8", "viscothermal 100 1000 12 2", "viscothermal 0 1000 12 8",
        "viscothermal 100 1000 12 9", "viscothermal NaN 1000 12 8", "gain 0.5", "viscothermal 100 1000 12 8 ignored"] {
        let text=LOSSY.replace("viscothermal 100 1000 12 8",suffix);
        assert!(PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).is_err(),"accepted {suffix}");
    }
    let mismatched=LOSSY.replace("8 32 1000","8 32 1500");
    assert!(PlateValvePerformance::from_bytes(mismatched.as_bytes(),37,&CancelGate::new()).is_err());
    // A memoryless terminal can coexist with selected wall loss. It is not
    // silently converted to a geometry-derived radiation boundary.
    let memoryless=LOSSY.replace("baffled-low-ka 1000","-0.8");
    let p=load(&memoryless);assert!(p.info().radiation_load.is_none());assert_eq!(p.viscothermal_losses().len(),1);
}
