//! Existing physical owners are the playback oracle, not an invented waveform.
use super::*;
use crate::bernoulli_aperture::performance::force::{ApertureForceEvent, ApertureForceProgram};
use crate::bernoulli_aperture::dynamic::force::PlateForceFootprint;
use crate::bernoulli_aperture::dynamic::relaxation::{InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec};
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::paired::{BladeState, PairedAperture, PairedApertureLimits};
use crate::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::pcm_wav::observation::DecimatedRenderer;
use crate::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::{gas::{GasSpec,GasState},visco::GeneralizedMaxwell};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};
use fs_scenario::gesture::{GestureEvent, GestureTarget, GestureTrack, GestureValue};
use std::io::Cursor;

const RATE: u32=96000;
const SAMPLES: u64=768;
fn plate(e:f64,density:f64)->PlateApertureReduction {
    let mesh=PlateMesh::rectangle(0.025,0.01,6,2);
    let root=mesh.nodes.iter().enumerate().filter(|(_,p)|p.0==0.0).map(|(i,_)|i).collect();
    let edges=mesh.boundary_edges().into_iter().filter(|&(a,b)|
        (mesh.nodes[a].0-0.025).abs()<1e-12 && (mesh.nodes[b].0-0.025).abs()<1e-12)
        .map(|(a,b)|[a,b]).collect();
    let section=PlateSection::isotropic(e,0.3,0.0003,density).unwrap();
    let chart=PlateChart::with_boundary_and_regions(mesh,section,root,vec![]).unwrap();
    PlateApertureReduction::from_chart(chart,PlateApertureOptions {
        assembly:AssemblyOptions {pretension:0.0,support:EdgeSupport::Clamped},
        eigenvalue_window:(1.0,(core::f64::consts::TAU*450.0).powi(2)),mode_index:0,
        eigensolver:SliceOptions::default(),max_nodes:100,max_triangles:100,slit_edges:edges,
        rest_opening_m:0.0002,damping_ratio:0.0,max_slit_mode_variation:0.25,max_slope:0.1,
    },&CancelGate::new_clock_free()).unwrap()
}
fn model(second_tau:f64)->PairedApertureNetwork {
    let plates=[plate(4e9,900.0),plate(8e9,1300.0)];
    let law=|p:&PlateApertureReduction,e,tau|PlateRelaxationSpec {
        regions:vec![PlateRelaxationRegion {triangles:(0..p.chart().mesh.tris.len()).collect(),
            material:GeneralizedMaxwell::new(e,vec![(0.5*e,tau)]).unwrap(),poisson_ratio:0.3,
            band_hz:(0.0,1000.0),provenance:"synthetic independent material playback".into()}],
        max_branches:1,max_dt_over_tau:0.1,max_angular_step:0.1,
    };
    let laws=[law(&plates[0],4e9,0.001),law(&plates[1],8e9,second_tau)];
    let gases:Vec<_>=[293.15,323.15].into_iter().map(|t|
        GasState::try_new(&GasSpec::dry_air_ussa1976(),t,101325.0).unwrap()).collect();
    let graph=TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,NetworkNode::Junction,
        NetworkNode::Termination {reflection:-0.8}],
        sections:vec![TubeSection {nodes:[0,1],length_m:0.12,radius_m:0.006,max_length_error_m:0.005},
            TubeSection {nodes:[1,2],length_m:0.18,radius_m:0.007,max_length_error_m:0.005}],
        sound_speed_m_s:gases[0].sound_speed,max_wave_memory_bytes:1<<20};
    let contact=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],0.0,2.0,
        "explicit no-contact playback comparison".into()).unwrap();
    let source=PairedAperture::from_plates(plates,gases[0],graph.inlet_impedance_with_gases(&gases).unwrap(),
        1.0/f64::from(RATE),SAMPLES,PairedApertureLimits {
            maximum_abs_displacement_m:0.001,maximum_abs_velocity_m_s:10.0,maximum_energy_j:1.0,
            max_evaluations:256,velocity_tolerance_m_s:1e-15,relative_tolerance:1e-13,
        },[BladeState {displacement_m:-3e-6,velocity_m_s:0.0},
            BladeState {displacement_m:2e-6,velocity_m_s:0.0}],contact).unwrap()
        .with_blade_relaxation(0,laws[0].clone(),InitialApertureMemory::Unrelaxed).unwrap()
        .with_blade_relaxation(1,laws[1].clone(),InitialApertureMemory::Relaxed).unwrap();
    PairedApertureNetwork::new(source,graph,gases,&CancelGate::new_clock_free()).unwrap()
}
fn phrase()->GestureSchedule {
    GestureSchedule::try_new(700,vec![GestureTrack {
        id:"mouth".into(),target:GestureTarget::BlowingPressure,initial:GestureValue::PressurePa(0.0),
        events:vec![GestureEvent {time_s:0.0,transition_s:0.002,value:GestureValue::PressurePa(4.0)},
            GestureEvent {time_s:0.004,transition_s:0.0,value:GestureValue::PressurePa(0.0)}],
    }]).unwrap()
}
fn config(block:usize)->AperturePerformanceConfig {
    AperturePerformanceConfig {sample_rate_hz:RATE,samples:SAMPLES,max_block:block,
        max_compile_work:10000,max_controls:64}
}
fn program(p:&PlateApertureReduction,sign:f64)->ApertureForceProgram {
    let tip=p.chart().mesh.nodes.iter().position(|p|(p.0-0.025).abs()<1e-12).unwrap();
    ApertureForceProgram {footprints:vec![PlateForceFootprint::Node(tip),
        PlateForceFootprint::Patch((0..p.chart().mesh.tris.len()).collect())],
        events:vec![
            ApertureForceEvent {sample:96,port:0,force_n:-1e-4*sign},
            ApertureForceEvent {sample:0,port:0,force_n:2e-4*sign},
            ApertureForceEvent {sample:0,port:0,force_n:1e-4*sign},
            ApertureForceEvent {sample:65,port:1,force_n:-2e-5*sign},
            ApertureForceEvent {sample:192,port:0,force_n:0.0},
            ApertureForceEvent {sample:192,port:1,force_n:0.0},
        ]}
}
fn receiver()->CircularOutletReceiver {
    CircularOutletReceiver {position_m:[0.0,0.0,0.25],radial_rings:4,angular_points:16,maximum_frequency_hz:2000.0}
}
fn renderer(block:usize,observation:ApertureObservation,tau:f64)->PairedAperturePerformance {
    let model=model(tau);let plates=model.aperture().plates().unwrap();
    let a=program(&plates[0],1.0);let b=program(&plates[1],-0.7);
    PairedAperturePerformance::new(model,observation,phrase(),config(block)).unwrap()
        .with_plate_forces(0,a).unwrap().with_plate_forces(1,b).unwrap()
}

#[test]
fn pressure_forces_material_histories_and_outlet_delay_match_direct_physics_across_blocks() {
    for observation in [ApertureObservation::Inlet,ApertureObservation::NetworkNode(2),
        ApertureObservation::NetworkBaffled {node:2,receiver:receiver()}] {
        let mut direct=model(0.003);let plates=direct.aperture().plates().unwrap();
        let programs=[program(&plates[0],1.0),program(&plates[1],-0.7)];
        let ports:Vec<Vec<_>>=(0..2).map(|i|programs[i].footprints.iter().cloned()
            .map(|f|plates[i].force_port(f).unwrap()).collect()).collect();
        let mut outlet=direct.baffled_receiver(2,RATE,receiver()).unwrap();
        let schedule=phrase();let mut held=[[0.0;2];2];let mut expected=Vec::new();
        for sample in 0..SAMPLES {
            let mut forces=[0.0;2];
            for i in 0..2 {
                // Original input order, including same-sample last-write-wins.
                for e in &programs[i].events {if e.sample==sample {held[i][e.port]=e.force_n;}}
                for (p,&f) in ports[i].iter().zip(&held[i]) {forces[i]+=p.generalized_force_n(f).unwrap();}
            }
            let frame=direct.step_with_forces(TubeDrive {
                upstream_pressure_pa:schedule.sample("mouth",sample*700/u64::from(RATE)).unwrap(),
                body_flow_m3_s:0.0,
            },forces).unwrap();
            let node=direct.node_frame(2).unwrap();
            let exterior=outlet.step(&[node.net_flow_into_node_m3_s]).unwrap();
            expected.push(match observation {
                ApertureObservation::Inlet=>frame.mechanics.aperture.bore_pressure_pa,
                ApertureObservation::NetworkNode(_)=>node.pressure_pa,_=>exterior,
            });
        }
        assert!(expected[500..].iter().any(|p|p.abs()>1e-12),"release must not erase ringdown");
        for block in [1,37,128] {
            let mut playback=renderer(block,observation,0.003);let mut actual=vec![0.0;SAMPLES as usize];
            for out in actual.chunks_mut(block) {playback.block(out).unwrap();}
            assert_eq!(actual,expected);
            assert_eq!(playback.system().aperture().blade_states(),direct.aperture().blade_states());
            for i in 0..2 {
                assert_eq!(playback.system().aperture().blade_relaxation(i).unwrap().memory_sqrt_j(),
                    direct.aperture().blade_relaxation(i).unwrap().memory_sqrt_j());
                assert!(playback.pending_plate_force_controls(i).is_empty());
                assert_eq!(playback.applied_plate_force_controls(i).len(),6);
            }
            assert_eq!(playback.held_generalized_forces_n(),[0.0;2]);
            assert!(playback.pending_controls().is_empty());
        }
    }
}

#[test]
fn changing_only_one_blades_relaxation_changes_actual_radiated_pressure() {
    let observation=ApertureObservation::NetworkBaffled {node:2,receiver:receiver()};
    let mut a=renderer(128,observation,0.001);let mut b=renderer(128,observation,0.02);
    assert_eq!(a.system().aperture().spec(),b.system().aperture().spec());
    assert_eq!(a.system().source_gases(),b.system().source_gases());
    let mut difference=0.0;
    for _ in 0..6 {
        let (mut x,mut y)=([0.0;128],[0.0;128]);a.block(&mut x).unwrap();b.block(&mut y).unwrap();
        difference+=x.iter().zip(y).map(|(a,b)|(a-b).abs()).sum::<f64>();
    }
    assert!(difference>1e-12);
}

#[test]
fn existing_decimated_wav_cancellation_preserves_both_material_histories_and_force_clocks() {
    let observation=ApertureObservation::NetworkBaffled {node:2,receiver:receiver()};
    let make=||DecimatedRenderer::new(renderer(37,observation,0.003),RATE,48000,37).unwrap();
    let stream=||Pcm16WavStream::new(Cursor::new(Vec::new()),48000,20.0,37).unwrap();
    let (mut direct,mut resumed)=(make(),make());let (mut expected,mut actual)=(stream(),stream());
    let mut scratch=[0.0;37];let gate=CancelGate::new_clock_free();
    render_pressure_pcm16(&mut direct,&mut expected,&gate,&mut scratch,384).unwrap();
    render_pressure_pcm16(&mut resumed,&mut actual,&gate,&mut scratch,101).unwrap();
    let before=*resumed.source().system().aperture().blade_states();
    let memory:Vec<_>=(0..2).map(|i|resumed.source().system().aperture()
        .blade_relaxation(i).unwrap().memory_sqrt_j().to_vec()).collect();
    let controls=resumed.source().pending_controls().to_vec();
    let cancel=CancelGate::new_clock_free();cancel.request();scratch.fill(-123.0);
    assert_eq!(render_pressure_pcm16(&mut resumed,&mut actual,&cancel,&mut scratch,283).unwrap(),
        ScheduledWavProgress::Cancelled {samples:0});
    assert_eq!(scratch,[-123.0;37]);
    assert_eq!(resumed.source().system().aperture().blade_states(),&before);
    assert_eq!(resumed.source().pending_controls(),controls);
    for i in 0..2 {assert_eq!(resumed.source().system().aperture().blade_relaxation(i).unwrap().memory_sqrt_j(),memory[i]);}
    render_pressure_pcm16(&mut resumed,&mut actual,&gate,&mut scratch,283).unwrap();
    assert_eq!(actual.finish().unwrap().0.into_inner(),expected.finish().unwrap().0.into_inner());
}

#[test]
fn controls_share_budgets_and_bad_callbacks_do_not_advance_or_publish() {
    let mut playback=renderer(37,ApertureObservation::Inlet,0.003);
    let before=*playback.system().aperture().blade_states();let mut sentinel=[-123.0;38];
    assert!(playback.block(&mut sentinel).is_err());assert_eq!(sentinel,[-123.0;38]);
    assert_eq!(playback.system().aperture().blade_states(),&before);assert!(playback.applied_controls().is_empty());
    assert!(playback.validate_sample_rate(48000).is_err());assert!(playback.validate_sample_count(SAMPLES+1).is_err());
    let build=||PairedAperturePerformance::new(model(0.003),ApertureObservation::Inlet,phrase(),config(37)).unwrap();
    let get_program=|r:&PairedAperturePerformance,i|program(&r.system().aperture().plates().unwrap()[i],1.0);
    let r=build();let p=get_program(&r,0);assert!(r.with_plate_forces(2,p).is_err());
    for bad_work in [false,true] {
        let r=build();let p=get_program(&r,0);let mut one=r.with_plate_forces(0,p).unwrap();
        let next=get_program(&one,1);
        if bad_work {one.config.max_compile_work=one.plate_forces[0].as_ref().unwrap().work;}
        else {one.config.max_controls=one.controls.len()+one.plate_forces[0].as_ref().unwrap().events.len();}
        assert!(one.with_plate_forces(1,next).is_err());
    }
    let r=build();let p=get_program(&r,0);
    let raw=vec![BladeForceControl {sample:0,blade:0,force_n:1e-4}];
    assert!(r.with_blade_forces(raw.clone()).unwrap().with_plate_forces(0,p).is_err());
    let r=build();let p=get_program(&r,0);
    assert!(r.with_plate_forces(0,p).unwrap().with_blade_forces(raw).is_err());
    let r=build();let mut p=get_program(&r,0);p.events[0].sample=SAMPLES;
    assert!(r.with_plate_forces(0,p).is_err());
    assert!(PairedAperturePerformance::new(model(0.003),ApertureObservation::TubeTerminal,phrase(),config(37)).is_err());
    let r=build();let mut p=get_program(&r,0);
    p.events=vec![ApertureForceEvent {sample:0,port:0,force_n:1e100}];
    let mut huge=r.with_plate_forces(0,p).unwrap();let memory=huge.system().aperture().blade_relaxation(0).unwrap().memory_sqrt_j().to_vec();
    assert!(huge.block(&mut [0.0]).is_err());assert_eq!(huge.samples_rendered(),0);
    assert!(huge.applied_controls().is_empty());assert!(huge.applied_plate_force_controls(0).is_empty());
    assert_eq!(huge.system().aperture().blade_relaxation(0).unwrap().memory_sqrt_j(),memory);
    assert!(matches!(huge.validate_sample_count(0),Err(RenderError::Poisoned)));
}
