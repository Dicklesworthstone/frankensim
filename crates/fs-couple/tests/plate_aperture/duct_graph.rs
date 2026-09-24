//! Supplied branch geometry drives one reciprocal valve/duct system.
use super::*;
use fs_couple::bernoulli_aperture::network::{ApertureNetwork,NetworkNode,TubeNetworkSpec,TubeSection};
use fs_couple::bernoulli_aperture::performance::{AperturePerformance,ApertureObservation,CoupledAperture};
use fs_couple::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use fs_couple::bernoulli_aperture::wall::{WallPatch,WallPin};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::pcm_wav::baffled::{BaffledPressure,CircularOutletReceiver,RayleighMedium};
use fs_couple::pcm_wav::stream::{Pcm16WavStream,ScheduledWavProgress,render_pressure_pcm16};
use std::io::Cursor;

const GRAPH:&str = "network 5 4 1048576\n\
duct_node inlet\n\
duct_node junction\n\
duct_node wall 0.001 0.1 200000 40\n\
duct_node baffled-low-ka 1500\n\
duct_node cavity 0.00001 0.0015 0.02 200000\n\
duct_section 0 1 0.075 0.007 0.002\n\
duct_section 1 2 0.075 0.007 0.002\n\
duct_section 2 3 0.1 0.008 0.002\n\
duct_section 1 4 0.04 0.002 0.002";
const OBSERVATION:&str="observation network-baffled 3 0 0 0.2 8 32 1500";
fn graph_input() -> String {
    INPUT.replace("tube 0.25 0.007 -0.8 0.002 1048576",GRAPH)
        .replace("observation inlet",OBSERVATION)
}
fn wall() -> WallPatch {
    WallPatch { area_m2:0.001, wall:WallPin { surface_density:0.1, stiffness_per_area:200000.0, resistance:40.0 } }
}
fn cavity() -> HelmholtzLoadSpec {
    HelmholtzLoadSpec { volume_m3:1e-5, neck_radius_m:0.0015,
        effective_neck_length_m:0.02, resistance_pa_s_m3:200000.0 }
}
fn graph(speed:f64,density:f64) -> TubeNetworkSpec {
    let radiation=BaffledRadiationLoad::new(0.008,density,speed,1.0/96000.0,1500.0,&CancelGate::new()).unwrap();
    TubeNetworkSpec {
        nodes:vec![NetworkNode::Inlet,NetworkNode::Junction,wall().shunt().unwrap(),
            radiation.termination(),cavity().termination(density,speed).unwrap()],
        sections:vec![
            TubeSection {nodes:[0,1],length_m:0.075,radius_m:0.007,max_length_error_m:0.002},
            TubeSection {nodes:[1,2],length_m:0.075,radius_m:0.007,max_length_error_m:0.002},
            TubeSection {nodes:[2,3],length_m:0.1,radius_m:0.008,max_length_error_m:0.002},
            TubeSection {nodes:[1,4],length_m:0.04,radius_m:0.002,max_length_error_m:0.002},
        ],sound_speed_m_s:speed,max_wave_memory_bytes:1<<20,
    }
}
fn fresh_valve(seed:&PlateValvePerformance,z:f64) -> DynamicAperture {
    let a=seed.renderer().system().aperture();let s=*a.spec();
    DynamicAperture::from_plate_with_closure(a.plate_reduction().unwrap().clone(),
        a.plate_closure().unwrap().spec().clone(),s.density_kg_m3,z,s.time_step_s,s.max_steps,a.state()).unwrap()
        .with_plate_relaxation(a.relaxation().unwrap().spec().clone(),InitialApertureMemory::Relaxed).unwrap()
}
fn network(renderer:&AperturePerformance) -> &ApertureNetwork {
    let CoupledAperture::Network(n)=renderer.system() else {panic!("graph replaced with a uniform tube")};n
}
fn bits(values:&[f64]) -> Vec<u64> {values.iter().map(|v|v.to_bits()).collect()}

#[test]
fn branch_wall_and_cavity_file_matches_direct_physical_network_and_its_work() {
    let seed=load(INPUT);
    let CoupledAperture::Tube(t)=seed.renderer().system() else {unreachable!()};
    let density=t.aperture().spec().density_kg_m3;
    let topology=graph(t.spec().sound_speed_m_s,density);
    let mut direct=ApertureNetwork::new(fresh_valve(&seed,topology.inlet_impedance(density).unwrap()),topology.clone()).unwrap();
    let source=load(&graph_input());let info=source.info();
    assert_eq!((info.duct_nodes,info.duct_sections,info.radiation_terminals),(5,4,1));
    assert_eq!(source.radiation_loads()[0].0,3);
    assert_eq!(source.info().radiation_load.unwrap().radius_m(),0.008,"use selected outlet, not inlet radius");
    assert_eq!(network(source.renderer()).spec(),&topology);
    let mut actual=source.into_renderer();
    let receiver=CircularOutletReceiver {position_m:[0.0,0.0,0.2],radial_rings:8,angular_points:32,maximum_frequency_hz:1500.0};
    let mut mic=BaffledPressure::circular_outlet(0.008,96000,receiver,
        RayleighMedium {density,sound_speed:t.spec().sound_speed_m_s}).unwrap();
    let mut cavity_energy=0.0_f64;let mut wall_motion=0.0_f64;let mut acoustic_peak=0.0_f64;
    for sample in 0..info.samples {
        let before=direct.stored_energy_j();
        let pressure=seed.renderer().schedule().sample("mouth",sample*700/96000).unwrap();
        let f=direct.step(TubeDrive {upstream_pressure_pa:pressure,body_flow_m3_s:0.0}).unwrap();
        let scale=before+f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs();
        assert!(f.balance_residual_j().abs()<=3e-10*scale.max(f64::MIN_POSITIVE));
        cavity_energy=cavity_energy.max(direct.node_frame(4).unwrap().stored_energy_j);
        wall_motion=wall_motion.max(wall().observe(&direct,2).unwrap().displacement_m.abs());
        let expected=mic.step(&[direct.node_frame(3).unwrap().net_flow_into_node_m3_s]).unwrap();
        acoustic_peak=acoustic_peak.max(expected.abs());
        let mut output=[0.0];actual.block(&mut output).unwrap();assert_eq!(output[0].to_bits(),expected.to_bits());
    }
    assert!(cavity_energy>0.0 && wall_motion>0.0 && acoustic_peak>1e-12);
    let retained=network(&actual);
    assert_eq!(retained.aperture().state(),direct.aperture().state());
    assert_eq!(retained.aperture().relaxation().unwrap().memory_sqrt_j(),direct.aperture().relaxation().unwrap().memory_sqrt_j());
    for node in 0..5 {assert_eq!(retained.node_frame(node),direct.node_frame(node));assert_eq!(retained.load_state(node),direct.load_state(node));}
    assert_eq!(retained.stored_energy_j().to_bits(),direct.stored_energy_j().to_bits());
}

#[test]
fn side_chamber_wall_and_return_path_change_actual_valve_motion_not_just_output() {
    let text=graph_input();let mut baseline=load(&text).into_renderer();
    let mut changed=[
        load(&text.replace("cavity 0.00001","cavity 0.00002")).into_renderer(),
        load(&text.replace("wall 0.001 0.1 200000 40","wall 0.001 0.1 400000 40")).into_renderer(),
        load(&text.replace("1 4 0.04 0.002","1 4 0.08 0.002")).into_renderer(),
    ];
    let mut moved=[false;3];let mut sound=[false;3];
    for _ in 0..128 {
        let mut reference=[0.0;37];baseline.block(&mut reference).unwrap();
        for (i,source) in changed.iter_mut().enumerate() {
            let mut pressure=[0.0;37];source.block(&mut pressure).unwrap();
            sound[i]|=bits(&pressure)!=bits(&reference);
            moved[i]|=source.system().aperture().state()!=baseline.system().aperture().state();
        }
    }
    assert!(moved.into_iter().all(|v|v) && sound.into_iter().all(|v|v));
    // Relocation observes a different delayed field, but no physics is reloaded.
    let near=load(&text).into_renderer();
    let far=load(&text.replace("3 0 0 0.2","3 0.1 0 0.4")).into_renderer();
    assert_eq!(network(&near).spec(),network(&far).spec());
    assert!(far.baffled_receiver().unwrap().delay_samples.0>near.baffled_receiver().unwrap().delay_samples.1);
}

#[test]
fn graph_cancellation_and_callback_partitions_preserve_all_branch_and_load_history() {
    let text=graph_input();let total=load(&text).info().samples;
    let mut baseline=load(&text).into_renderer();let mut expected=vec![0.0;total as usize];
    for chunk in expected.chunks_mut(37) {baseline.block(chunk).unwrap();}
    for partition in [1,512] {
        let mut r=PlateValvePerformance::from_bytes(text.as_bytes(),partition,&CancelGate::new()).unwrap().into_renderer();
        let mut actual=vec![0.0;total as usize];for chunk in actual.chunks_mut(partition) {r.block(chunk).unwrap();}
        assert_eq!(bits(&actual),bits(&expected));
    }
    let mut r=load(&text).into_renderer();
    let mut pcm=Pcm16WavStream::new(Cursor::new(Vec::new()),96000,20.0,37).unwrap();let mut scratch=[0.0;37];
    render_pressure_pcm16(&mut r,&mut pcm,&CancelGate::new(),&mut scratch,1001).unwrap();
    let state=r.system().aperture().state();let energy=network(&r).stored_energy_j().to_bits();
    let nodes:Vec<_>=(0..5).map(|i|*network(&r).node_frame(i).unwrap()).collect();
    let pending=r.pending_controls().to_vec();let gate=CancelGate::new();gate.request();scratch.fill(-123.0);
    assert_eq!(render_pressure_pcm16(&mut r,&mut pcm,&gate,&mut scratch,total-1001).unwrap(),ScheduledWavProgress::Cancelled{samples:0});
    assert_eq!(scratch,[-123.0;37]);assert_eq!(r.system().aperture().state(),state);
    assert_eq!(network(&r).stored_energy_j().to_bits(),energy);assert_eq!(r.pending_controls(),pending);
    for (i,value) in nodes.iter().enumerate() {assert_eq!(network(&r).node_frame(i),Some(value));}
    render_pressure_pcm16(&mut r,&mut pcm,&CancelGate::new(),&mut scratch,total-1001).unwrap();
    let encoded=fs_couple::pcm_wav::encode_pcm16_wav(&expected,96000,20.0).unwrap().0;
    assert_eq!(pcm.finish().unwrap().0.into_inner(),encoded);
    assert!(r.validate_sample_count(1).is_err());
}

#[test]
fn graph_topology_physical_loads_and_observer_ambiguities_refuse_without_substitution() {
    let source=graph_input();
    for text in [
        source.replace("network 5 4","network 65 4"),source.replace("network 5 4","network 5 129"),
        source.replace("network 5 4 1048576","network 5 4 16"),
        source.replace("duct_node junction","duct_node inlet"),
        source.replace("duct_section 0 1","duct_section 0 0"),
        source.replace("duct_section 1 4","duct_section 1 99"),
        source.replace("duct_section 1 4","duct_section 2 3"),
        source.replace("0.04 0.002 0.002","0 0.002 0.002"),
        source.replace("0.04 0.002 0.002","0.04 0.002 0"),
        source.replace("cavity 0.00001","cavity -0.00001"),
        source.replace("wall 0.001 0.1","wall 0.001 -0.1"),
        source.replace("baffled-low-ka 1500","baffled-low-ka 10000"),
        source.replace(OBSERVATION,"observation terminal"),
        source.replace(OBSERVATION,"observation network-node 5"),
        source.replace("network-baffled 3","network-baffled 4"),
        source.replace("network-baffled 3","network-baffled 2"),
        source.replace("8 32 1500","8 32 1600"),
        source.replace("duct_node junction","duct_node reflection -0.8"),
        source.replace("duct_node wall 0.001 0.1 200000 40","duct_node series -1 1 none"),
    ] {assert!(PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).is_err(),"accepted {text}");}
    let old=load(INPUT);assert!(matches!(old.renderer().system(),CoupledAperture::Tube(_)));
    assert_eq!((old.info().duct_nodes,old.info().duct_sections),(2,1));
}

#[test]
fn more_than_one_radiating_terminal_is_retained_but_only_the_named_outlet_is_observed() {
    let two=graph_input().replace("network 5 4","network 6 5")
        .replace("duct_section 0 1","duct_node baffled-low-ka 1500\nduct_section 0 1")
        .replace(OBSERVATION,&format!("duct_section 1 5 0.06 0.003 0.002\n{OBSERVATION}"));
    let loaded=load(&two);
    assert_eq!(loaded.info().radiation_terminals,2);
    assert_eq!(loaded.radiation_loads().iter().map(|(n,_)|*n).collect::<Vec<_>>(),vec![3,5]);
    assert_eq!(loaded.info().radiation_load.unwrap().radius_m(),0.008);
    let selected=load(&two.replace("network-baffled 3","network-baffled 5"));
    assert_eq!(selected.info().radiation_load.unwrap().radius_m(),0.003);
    let(mut a,mut b)=(loaded.into_renderer(),selected.into_renderer());let mut distinct=false;
    for _ in 0..80 {
        let(mut p,mut q)=([0.0;37],[0.0;37]);a.block(&mut p).unwrap();b.block(&mut q).unwrap();distinct|=p!=q;
        assert_eq!(a.system().aperture().state(),b.system().aperture().state());
    }
    assert!(distinct);
    for node in [3,5] {assert!(network(&a).node_frame(node).unwrap().stored_energy_j>0.0);}
}
