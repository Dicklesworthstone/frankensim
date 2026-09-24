//! Numerical/physical-owner regressions; no measured-instrument claim.
use super::*;
use fs_couple::bernoulli_aperture::viscothermal::{WideTubeLoss, ViscothermalSectionSpec, with_viscothermal_sections};
use fs_couple::bernoulli_aperture::network::{TubeNetworkSpec, TubeSection, NetworkNode, ApertureNetwork, ApertureNetworkFrame};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_couple::bernoulli_aperture::dynamic::{DynamicAperture, ApertureState, ApertureTerminal};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_vfit::waveguide::network::{WaveguideNetwork, NetworkSegment};
const RATE:u32=96000;
fn gas()->GasState {GasState::try_new_moist_air(293.15,101325.0,0.0).unwrap()}
fn choice()->ViscothermalSectionSpec {ViscothermalSectionSpec {
    section:0,cells:16,minimum_frequency_hz:50.0,maximum_frequency_hz:1500.0,
}}
fn graph(g:&GasState)->TubeNetworkSpec {TubeNetworkSpec {
    nodes:vec![NetworkNode::Inlet,NetworkNode::Termination{reflection:0.0}],
    sections:vec![TubeSection{nodes:[0,1],length_m:0.25,radius_m:0.007,max_length_error_m:0.002}],
    sound_speed_m_s:g.sound_speed,max_wave_memory_bytes:4<<20,
}}
fn lowered()->TubeNetworkSpec {let g=gas();with_viscothermal_sections(graph(&g),&g,1.0/f64::from(RATE),&[choice()],&CancelGate::new()).unwrap().0}
fn waves(spec:&TubeNetworkSpec)->WaveguideNetwork {
    let dt=1.0/f64::from(RATE);let g=gas();
    let lines:Vec<_>=spec.sections.iter().map(|s|NetworkSegment{nodes:s.nodes,
        one_way_samples:(s.length_m/(g.sound_speed*dt)).round() as usize,
        impedance_pa_s_m3:g.density*g.sound_speed/(core::f64::consts::PI*s.radius_m*s.radius_m)}).collect();
    WaveguideNetwork::new(&spec.nodes,&lines,dt,spec.max_wave_memory_bytes).unwrap()
}

#[test]
fn complex_loss_follows_physical_viscous_and_thermal_coefficients_not_an_audio_gain() {
    let g=gas();let radius=0.007;let dt=1.0/f64::from(RATE);
    let loss=WideTubeLoss::new(radius,&g,dt,[50.0,1500.0],&CancelGate::new()).unwrap();
    let area=core::f64::consts::PI*radius*radius;
    for f in [53.0,97.0,233.0,601.0,1199.0,1497.0] {
        let w=core::f64::consts::TAU*f;
        let rv=radius*(w*g.density/g.dynamic_viscosity).sqrt();
        let r=w*(g.density/area)*core::f64::consts::SQRT_2/rv;
        let conductance=w*area/(g.density*g.sound_speed*g.sound_speed)*(g.gamma-1.0)
            *core::f64::consts::SQRT_2/(rv*g.prandtl.sqrt());
        for discrete in [false,true] {let (z,y)=loss.excess_at(f,discrete).unwrap();
            assert!((z-C64::new(r,-r)).abs()/r<0.071);
            assert!((y-C64::new(conductance,-conductance)).abs()/conductance<0.071);
            assert!(z.re>0.0 && y.re>0.0 && z.im<0.0 && y.im<0.0);
        }
    }
    // Synthetic transport scaling at fixed inertia/compliance, not a new gas card.
    let mut thicker=g;thicker.dynamic_viscosity*=4.0;thicker.thermal_conductivity*=4.0;
    let changed=WideTubeLoss::new(radius,&thicker,dt,[50.0,1500.0],&CancelGate::new()).unwrap();
    let (a,b)=loss.excess_at(601.0,true).unwrap();let (x,y)=changed.excess_at(601.0,true).unwrap();
    assert!((x-a.scale(2.0)).abs()<1e-12*x.abs());assert!((y-b.scale(2.0)).abs()<1e-12*y.abs());
    assert_eq!(loss.series(0.02).unwrap().base().inertance_pa_s2_m3,0.0);
    assert_eq!(loss.shunt(0.02).unwrap().compliance_m3_pa(),0.0);
    assert!(loss.excess_at(10.0,true).is_err());
}

#[test]
fn expansion_preserves_endpoints_inlet_impedance_and_every_original_transit_sample() {
    let g=gas();let original=graph(&g);let original_z=original.inlet_impedance(g.density).unwrap();
    let (expanded,report)=with_viscothermal_sections(original.clone(),&g,1.0/f64::from(RATE),&[choice()],&CancelGate::new()).unwrap();
    assert_eq!(&expanded.nodes[..2],original.nodes.as_slice());
    assert_eq!(expanded.inlet_impedance(g.density).unwrap().to_bits(),original_z.to_bits());
    assert_eq!(expanded.sections.first().unwrap().nodes[0],0);assert_eq!(expanded.sections.last().unwrap().nodes[1],1);
    let n:usize=expanded.sections.iter().map(|s|(s.length_m/g.sound_speed*f64::from(RATE)).round() as usize).sum();
    assert_eq!(n,report[0].one_way_samples);assert_eq!(expanded.sections.len(),64);
    assert_eq!(expanded.nodes.iter().filter(|x|matches!(x,NetworkNode::Series{..})).count(),32);
    assert_eq!(expanded.nodes.iter().filter(|x|matches!(x,NetworkNode::ShuntAdmittance{..})).count(),16);
    assert!(report[0].max_scattering_error<0.03);
    let mut w=waves(&expanded);
    for i in 0..n {w.step(if i==0 {1.0} else {0.0}).unwrap();assert_eq!(w.node_frame(1).unwrap().pressure_pa,0.0);}
    w.step(0.0).unwrap();assert!(w.node_frame(1).unwrap().pressure_pa>0.0);
    let (untouched,reports)=with_viscothermal_sections(original.clone(),&g,1.0/f64::from(RATE),&[],&CancelGate::new()).unwrap();
    assert_eq!(untouched,original);assert!(reports.is_empty());
}

#[test]
fn actual_discrete_waveguide_attenuates_and_closes_independent_storage_and_work() {
    let spec=lowered();let mut lossy=waves(&spec);let mut ideal=waves(&graph(&gas()));
    let mut dissipated=0.0;let mut return_energy=0.0;let mut stored=0.0;
    let mut pressure_difference=0.0;
    for n in 0..4096 {
        let drive=if n<512 {2.0*(core::f64::consts::TAU*700.0*n as f64/f64::from(RATE)).sin()} else {0.0};
        let a=lossy.step(drive).unwrap();let b=ideal.step(drive).unwrap();
        let scale=stored+a.stored_energy_j+a.inlet_work_j.abs()+a.interior_loss_j+a.terminal_loss_j;
        assert!(a.balance_residual_j().abs()<=3e-11*scale.max(1e-25));
        assert!(a.interior_loss_j>=0.0);dissipated+=a.interior_loss_j;stored=a.stored_energy_j;
        if n>700 {return_energy+=a.inlet_work_j.min(0.0).abs();}
        pressure_difference+=(lossy.node_frame(1).unwrap().pressure_pa-ideal.node_frame(1).unwrap().pressure_pa).abs();
        assert_eq!(b.interior_loss_j,0.0);
    }
    assert!(dissipated>0.0 && pressure_difference>0.1 && return_energy>0.0);
    let sum: f64=(2..spec.nodes.len()).map(|n|lossy.node_frame(n).unwrap().stored_energy_j).sum();
    assert!(sum>0.0,"thermal and viscous states survive release; they are not an output decay multiplier");
}

fn coupled()->ApertureNetwork {
    let g=gas();let spec=lowered();let r=reduction(4e9,900.0);let h=r.options().rest_opening_m;
    let lay=fs_dcontact::Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "authored viscothermal coupled test".into()).unwrap();
    let valve=DynamicAperture::from_plate(r,g.density,spec.inlet_impedance(g.density).unwrap(),
        1.0/f64::from(RATE),2048,ApertureState{opening_m:h,opening_velocity_m_s:0.0},lay).unwrap();
    ApertureNetwork::new(valve,spec).unwrap()
}
#[test]
fn viscous_thermal_and_valve_history_survive_cancellation_and_bad_input_retry() {
    let mut direct=coupled();let mut resumed=coupled();
    let inputs:Vec<_>=(0..1024).map(|n|TubeDrive{upstream_pressure_pa:if n<512 {5.0} else {0.0},body_flow_m3_s:0.0}).collect();
    let expected:Vec<_>=inputs.iter().map(|&d|direct.step(d).unwrap()).collect();
    let mut actual=vec![ApertureNetworkFrame::default();1024];
    resumed.advance_block(&inputs[..257],&mut actual[..257],&CancelGate::new()).unwrap();
    let before=resumed.stored_energy_j();let state=resumed.aperture().state();
    let cancelled=CancelGate::new();cancelled.request();
    let p=resumed.advance_block(&inputs[257..],&mut actual[257..],&cancelled).unwrap();
    assert_eq!(p.terminal,ApertureTerminal::Cancelled);assert_eq!(p.completed,0);
    assert!(resumed.step(TubeDrive{upstream_pressure_pa:f64::NAN,body_flow_m3_s:0.0}).is_err());
    assert_eq!(resumed.stored_energy_j().to_bits(),before.to_bits());assert_eq!(resumed.aperture().state(),state);
    for (source,out) in inputs[257..].chunks(37).zip(actual[257..].chunks_mut(37)) {
        resumed.advance_block(source,out,&CancelGate::new()).unwrap();
    }
    for (a,b) in actual.iter().zip(expected) {assert_eq!(a.aperture,b.aperture);assert_eq!(a.network,b.network);}
}

#[test]
fn unresolved_shear_clock_mesh_and_memory_refuse_without_a_lossless_fallback() {
    let g=gas();let dt=1.0/f64::from(RATE);
    assert!(WideTubeLoss::new(0.00001,&g,dt,[50.0,1500.0],&CancelGate::new()).is_err());
    assert!(WideTubeLoss::new(0.007,&g,dt,[1500.0,50.0],&CancelGate::new()).is_err());
    assert!(WideTubeLoss::new(0.007,&g,0.001,[50.0,1500.0],&CancelGate::new()).is_err());
    for selection in [ViscothermalSectionSpec{cells:1,..choice()},ViscothermalSectionSpec{cells:128,..choice()},
        ViscothermalSectionSpec{section:2,..choice()},ViscothermalSectionSpec{cells:0,..choice()}] {
        assert!(with_viscothermal_sections(graph(&g),&g,dt,&[selection],&CancelGate::new()).is_err());
    }
    assert!(with_viscothermal_sections(graph(&g),&g,dt,&[choice(),choice()],&CancelGate::new()).is_err());
    let mut small=graph(&g);small.max_wave_memory_bytes=64;
    assert!(with_viscothermal_sections(small,&g,dt,&[choice()],&CancelGate::new()).is_err());
    let cancelled=CancelGate::new();cancelled.request();
    assert!(with_viscothermal_sections(graph(&g),&g,dt,&[choice()],&cancelled).is_err());
}
