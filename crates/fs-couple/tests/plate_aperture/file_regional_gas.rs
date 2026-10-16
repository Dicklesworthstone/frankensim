use super::*;
use fs_couple::bernoulli_aperture::network::{ApertureNetwork,NetworkNode,TubeNetworkSpec,TubeSection};
use fs_couple::bernoulli_aperture::performance::{AperturePerformance,CoupledAperture};
use fs_couple::bernoulli_aperture::viscothermal::{with_regional_viscothermal_sections,ViscothermalSectionSpec};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use fs_couple::pcm_wav::baffled::{BaffledPressure,CircularOutletReceiver,RayleighMedium};
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::pcm_wav::stream::{Pcm16WavStream,render_pressure_pcm16,ScheduledWavProgress};
use std::io::Cursor;
const REGIONAL:&str=include_str!("../../examples/plate-valve-regional-gas.performance");
fn gas(t:f64,rh:f64)->GasState {GasState::try_new_moist_air(t,101325.0,rh).unwrap()}
fn network(r:&AperturePerformance)->&ApertureNetwork {
    let CoupledAperture::Network(n)=r.system() else {panic!("regional source must retain a graph")};n
}
fn uniform()->String {REGIONAL.replace(" gas 313.15 0.5","").replace(" gas 283.15 0","")}

#[test]
fn regional_file_matches_direct_local_radiation_cavity_losses_and_physical_motion() {
    let seed=load(INPUT);let loaded=load(REGIONAL);let gases=vec![gas(293.15,0.0),gas(313.15,0.5),gas(283.15,0.0)];
    let dt=1.0/96000.0;
    let radiation=BaffledRadiationLoad::new(0.007,gases[1].density,gases[1].sound_speed,dt,1000.0,&CancelGate::new()).unwrap();
    let chamber=HelmholtzLoadSpec {volume_m3:1e-5,neck_radius_m:0.0015,effective_neck_length_m:0.02,resistance_pa_s_m3:200000.0};
    let original=TubeNetworkSpec {
        nodes:vec![NetworkNode::Inlet,NetworkNode::Junction,radiation.termination(),chamber.termination(gases[2].density,gases[2].sound_speed).unwrap()],
        sections:vec![TubeSection {nodes:[0,1],length_m:0.125,radius_m:0.007,max_length_error_m:0.002},
            TubeSection {nodes:[1,2],length_m:0.125,radius_m:0.007,max_length_error_m:0.002},
            TubeSection {nodes:[1,3],length_m:0.043,radius_m:0.002,max_length_error_m:0.002}],
        sound_speed_m_s:gases[0].sound_speed,max_wave_memory_bytes:1<<20,
    };
    let select=|section|ViscothermalSectionSpec {section,cells:6,minimum_frequency_hz:100.0,maximum_frequency_hz:1000.0};
    let (graph,mapped,reports)=with_regional_viscothermal_sections(original,gases.clone(),dt,&[select(0),select(1)],&CancelGate::new()).unwrap();
    assert_eq!(network(loaded.renderer()).spec(),&graph);
    assert_eq!(network(loaded.renderer()).section_gases(),Some(mapped.as_slice()));
    for (i,report) in reports.iter().enumerate() {
        assert_eq!(loaded.viscothermal_losses()[i].source,report.source);
        assert_eq!(report.loss.dynamic_viscosity_pa_s(),gases[i].dynamic_viscosity);
        assert_eq!(report.loss.prandtl(),gases[i].prandtl);
        assert_eq!(report.one_way_samples,(0.125/(gases[i].sound_speed*dt)).round() as usize);
        for index in report.section_range[0]..report.section_range[1] {assert_eq!(mapped[index],gases[i]);}
    }
    let a=seed.renderer().system().aperture();let s=a.spec();
    let valve=DynamicAperture::from_plate_with_closure(a.plate_reduction().unwrap().clone(),a.plate_closure().unwrap().spec().clone(),
        gases[0].density,graph.inlet_impedance_with_gases(&mapped).unwrap(),dt,s.max_steps,a.state()).unwrap()
        .with_plate_relaxation(a.relaxation().unwrap().spec().clone(),InitialApertureMemory::Relaxed).unwrap();
    let mut direct=ApertureNetwork::with_section_gases(valve,graph,mapped).unwrap();
    let mut mic=BaffledPressure::circular_outlet(0.007,96000,CircularOutletReceiver {
        position_m:[0.0,0.0,0.2],radial_rings:8,angular_points:32,maximum_frequency_hz:1000.0,
    },RayleighMedium {density:gases[1].density,sound_speed:gases[1].sound_speed}).unwrap();
    let mut actual=loaded.into_renderer();let mut baseline=load(&uniform()).into_renderer();let mut changed=false;
    for n in 0..1024 {
        let before=direct.stored_energy_j();
        let f=direct.step(TubeDrive {upstream_pressure_pa:seed.renderer().schedule().sample("mouth",n*700/96000).unwrap(),body_flow_m3_s:0.0}).unwrap();
        let scale=before+f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs();
        assert!(f.balance_residual_j().abs()<=3e-10*scale.max(1e-30));
        let expected=mic.step(&[direct.node_frame(2).unwrap().net_flow_into_node_m3_s]).unwrap();
        let mut out=[0.0];actual.block(&mut out).unwrap();baseline.block(&mut [0.0]).unwrap();
        assert_eq!(out[0].to_bits(),expected.to_bits());
        changed|=actual.system().aperture().state()!=baseline.system().aperture().state();
    }
    assert!(changed,"regional gas must affect the mechanical feedback, not only receiver metadata");
    assert_eq!(actual.system().aperture().state(),direct.aperture().state());
}

#[test]
fn explicit_uniform_gas_keeps_legacy_bits_and_inlet_override_sets_the_actual_flow_density() {
    let explicit=REGIONAL.replace("gas 313.15 0.5","gas 293.15 0").replace("gas 283.15 0","gas 293.15 0");
    let a=render(&explicit);let b=render(&uniform());
    assert_eq!(a.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),b.iter().map(|p|p.to_bits()).collect::<Vec<_>>());
    let warm_inlet=REGIONAL.replace("duct_section 0 1 0.125 0.007 0.002 viscothermal",
        "duct_section 0 1 0.125 0.007 0.002 gas 303.15 0.2 viscothermal");
    let p=load(&warm_inlet);let inlet=gas(303.15,0.2);
    assert_eq!(p.renderer().system().aperture().spec().density_kg_m3,inlet.density);
    assert_eq!(network(p.renderer()).spec().sound_speed_m_s,inlet.sound_speed);
}

#[test]
fn cancelled_regional_pcm_preserves_all_section_load_receiver_and_material_history() {
    let build=||DecimatedRenderer::new(load(REGIONAL).into_renderer(),96000,48000,37).unwrap();
    let stream=||Pcm16WavStream::new(Cursor::new(Vec::new()),48000,0.01,37).unwrap();
    let mut direct=build();let mut expected=stream();let mut scratch=[0.0;37];
    render_pressure_pcm16(&mut direct,&mut expected,&CancelGate::new(),&mut scratch,4801).unwrap();
    let mut resumed=build();let mut actual=stream();
    render_pressure_pcm16(&mut resumed,&mut actual,&CancelGate::new(),&mut scratch,1001).unwrap();
    let before=network(resumed.source()).stored_energy_j().to_bits();let state=resumed.source().system().aperture().state();
    let gate=CancelGate::new();gate.request();scratch.fill(-123.0);
    assert_eq!(render_pressure_pcm16(&mut resumed,&mut actual,&gate,&mut scratch,3800).unwrap(),ScheduledWavProgress::Cancelled {samples:0});
    assert_eq!(scratch,[-123.0;37]);assert_eq!(before,network(resumed.source()).stored_energy_j().to_bits());
    assert_eq!(state,resumed.source().system().aperture().state());
    render_pressure_pcm16(&mut resumed,&mut actual,&CancelGate::new(),&mut scratch,3800).unwrap();
    assert_eq!(actual.finish().unwrap().0.into_inner(),expected.finish().unwrap().0.into_inner());
}

#[test]
fn malformed_regional_inputs_never_fall_back_to_ambient_or_drop_selected_losses() {
    for replacement in ["gas","gas NaN 0.5","gas 313.15 2","gas 400 0.5","gas 313.15 0.5 gas 293.15 0",
        "gas 313.15 101325 0.5","temperature 313.15 0.5"] {
        let text=REGIONAL.replace("gas 313.15 0.5",replacement);
        assert!(PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).is_err(),"accepted {replacement}");
    }
    let unresolvable=REGIONAL.replace("gas 313.15 0.5 viscothermal 100 1000 6 8",
        "gas 313.15 0.5 viscothermal 100 1000 128 8");
    assert!(PlateValvePerformance::from_bytes(unresolvable.as_bytes(),37,&CancelGate::new()).is_err());
}
