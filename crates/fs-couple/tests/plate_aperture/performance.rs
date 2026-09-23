use super::*;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use fs_couple::bernoulli_aperture::dynamic::relaxation::{InitialApertureMemory, PlateRelaxationSpec, PlateRelaxationRegion};
use fs_couple::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_couple::bernoulli_aperture::performance::{AperturePerformance, AperturePerformanceConfig, ApertureObservation, CoupledAperture};
use fs_couple::bernoulli_aperture::tube::{ApertureTube, TubeDrive, UniformTubeSpec};
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::pcm_wav::observation::{PressureRenderer, DecimatedRenderer};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_couple::render::RenderError;
use fs_material::visco::GeneralizedMaxwell;
use fs_scenario::gesture::{GestureSchedule, GestureTrack, GestureTarget, GestureEvent, GestureValue};
use std::io::Cursor;
const RATE: u32 = 96_000;
const SAMPLES: u64 = 9602;
fn phrase() -> GestureSchedule {
    GestureSchedule::try_new(700, vec![GestureTrack {
        id: "mouth".into(), target: GestureTarget::BlowingPressure,
        initial: GestureValue::PressurePa(0.0), events: vec![
            GestureEvent { time_s:0.0, transition_s:0.01, value:GestureValue::PressurePa(5.0) },
            GestureEvent { time_s:0.005, transition_s:0.01, value:GestureValue::PressurePa(6.0) },
            GestureEvent { time_s:0.02, transition_s:0.0, value:GestureValue::PressurePa(0.0) },
            GestureEvent { time_s:0.04, transition_s:0.005, value:GestureValue::PressurePa(5.0) },
            GestureEvent { time_s:0.06, transition_s:0.0, value:GestureValue::PressurePa(0.0) },
        ],
    }]).unwrap()
}
fn config(block: usize) -> AperturePerformanceConfig {
    AperturePerformanceConfig { sample_rate_hz:RATE, samples:SAMPLES, max_block:block,
        max_compile_work:10000, max_controls:128 }
}
fn uniform() -> UniformTubeSpec {
    UniformTubeSpec { length_m:0.25, radius_m:0.007, sound_speed_m_s:343.0,
        terminal_reflection:-0.8, max_length_error_m:0.004, max_wave_memory_bytes:1<<20 }
}
fn valve() -> DynamicAperture {
    let (chart, mut options) = fixture(4e9,900.0); options.damping_ratio=0.0;
    let p = PlateApertureReduction::from_chart(chart,options,&CancelGate::new()).unwrap();
    let triangles = p.chart().mesh.tris.len();
    let h = p.options().rest_opening_m;
    let closure = PlateClosureSpec {
        nodal_rest_gap_m:p.chart().mesh.nodes.iter().map(|n|h*(0.05+1.9*n.1/0.01)).collect(),
        lay_triangles:(0..triangles).collect(),stiffness_pa_per_m_alpha:1e12,
        alpha:2.0,internal_loss_s_per_m:0.5,provenance:"authored scheduler test".into(),max_penetration_m:0.0002,
    };
    DynamicAperture::from_plate_with_closure(p,closure,1.2,uniform().characteristic_impedance(1.2).unwrap(),
        1.0/f64::from(RATE),SAMPLES,ApertureState{opening_m:h,opening_velocity_m_s:0.0}).unwrap()
        .with_plate_relaxation(PlateRelaxationSpec {
            regions:vec![PlateRelaxationRegion {triangles:(0..triangles).collect(),
                material:GeneralizedMaxwell::new(4e9,vec![(2e9,0.001),(1e9,0.01)]).unwrap(),
                poisson_ratio:0.3,band_hz:(0.0,10000.0),provenance:"synthetic constitutive input".into()}],
            max_branches:2,max_dt_over_tau:1.0,max_angular_step:1.0,
        },InitialApertureMemory::Relaxed).unwrap()
}
fn model() -> ApertureTube { ApertureTube::new(valve(),uniform()).unwrap() }
fn renderer(block: usize, observation: ApertureObservation) -> AperturePerformance {
    AperturePerformance::new(CoupledAperture::Tube(model()),observation,phrase(),config(block)).unwrap()
}
fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|x|x.to_bits()).collect() }

#[test]
fn authored_phrases_drive_retained_contact_material_and_air_without_callback_dependence() {
    for observation in [ApertureObservation::Inlet,ApertureObservation::TubeTerminal] {
        let mut direct=model(); let source=phrase(); let mut expected=Vec::new();
        for i in 0..SAMPLES {
            let f=direct.step(TubeDrive {upstream_pressure_pa:source.sample("mouth",i*700/u64::from(RATE)).unwrap(),body_flow_m3_s:0.0}).unwrap();
            expected.push(if observation==ApertureObservation::Inlet {f.aperture.bore_pressure_pa} else {f.waveguide.terminal_pressure_pa});
        }
        assert!(expected[6000..].iter().any(|x|x.abs()>1e-12),"release must preserve ongoing physical ringing");
        for block in [1,37,512] {
            let mut scheduled=renderer(block,observation); let mut actual=vec![0.0;SAMPLES as usize];
            for chunk in actual.chunks_mut(block) {scheduled.block(chunk).unwrap();}
            assert_eq!(bits(&actual),bits(&expected));
            let retained=scheduled.system().aperture();
            assert_eq!(retained.state(),direct.aperture().state());
            assert_eq!(retained.relaxation().unwrap().memory_sqrt_j(),direct.aperture().relaxation().unwrap().memory_sqrt_j());
            assert!(retained.plate_closure().is_some());
            assert!(scheduled.pending_controls().is_empty());
        }
    }
}

#[test]
fn selected_network_pressure_is_the_existing_physical_node_not_an_extra_source() {
    let s=uniform();
    let topology=TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,NetworkNode::Termination{reflection:-0.8}],
        sections:vec![TubeSection{nodes:[0,1],length_m:s.length_m,radius_m:s.radius_m,max_length_error_m:s.max_length_error_m}],
        sound_speed_m_s:s.sound_speed_m_s,max_wave_memory_bytes:1<<20};
    let mut direct=ApertureNetwork::new(valve(),topology.clone()).unwrap();
    let network=ApertureNetwork::new(valve(),topology).unwrap();
    let mut scheduled=AperturePerformance::new(CoupledAperture::Network(network),ApertureObservation::NetworkNode(1),phrase(),config(37)).unwrap();
    let source=phrase();let mut actual=[0.0;37];
    for group in 0..8 {
        scheduled.block(&mut actual).unwrap();
        for (i,&pressure) in actual.iter().enumerate() {
            let n=group*37+i;
            direct.step(TubeDrive{upstream_pressure_pa:source.sample("mouth",n as u64*700/u64::from(RATE)).unwrap(),body_flow_m3_s:0.0}).unwrap();
            assert_eq!(pressure.to_bits(),direct.node_frame(1).unwrap().pressure_pa.to_bits());
        }
    }
}

#[test]
fn decimated_wav_cancellation_retains_complete_source_controls_and_material_memory() {
    fn render_source() -> DecimatedRenderer<AperturePerformance> {
        DecimatedRenderer::new(renderer(37,ApertureObservation::TubeTerminal),RATE,48000,37).unwrap()
    }
    fn stream() -> Pcm16WavStream<Cursor<Vec<u8>>> {
        Pcm16WavStream::new(Cursor::new(Vec::new()),48000,20.0,37).unwrap()
    }
    let mut direct=render_source();let mut expected=stream();let mut scratch=[0.0;37];
    render_pressure_pcm16(&mut direct,&mut expected,&CancelGate::new(),&mut scratch,4801).unwrap();
    let mut resumed=render_source();let mut actual=stream();
    render_pressure_pcm16(&mut resumed,&mut actual,&CancelGate::new(),&mut scratch,1001).unwrap();
    let before=resumed.source().system().aperture().state();
    let memory=resumed.source().system().aperture().relaxation().unwrap().memory_sqrt_j().to_vec();
    let pending=resumed.source().pending_controls().to_vec();
    let cancelled=CancelGate::new();cancelled.request();scratch.fill(-123.0);
    assert_eq!(render_pressure_pcm16(&mut resumed,&mut actual,&cancelled,&mut scratch,3800).unwrap(),ScheduledWavProgress::Cancelled{samples:0});
    assert_eq!(scratch,[-123.0;37]);assert_eq!(resumed.source().system().aperture().state(),before);
    assert_eq!(resumed.source().system().aperture().relaxation().unwrap().memory_sqrt_j(),memory);
    assert_eq!(resumed.source().pending_controls(),pending);
    render_pressure_pcm16(&mut resumed,&mut actual,&CancelGate::new(),&mut scratch,3800).unwrap();
    assert_eq!(actual.finish().unwrap().0.into_inner(),expected.finish().unwrap().0.into_inner());
}

#[test]
fn invalid_windows_and_observers_refuse_before_any_initial_control_or_state_changes() {
    let mut r=renderer(37,ApertureObservation::Inlet);let initial=r.system().aperture().state();
    assert!(r.validate_sample_count(SAMPLES+1).is_err());assert!(r.validate_sample_rate(48000).is_err());
    let mut sentinel=[-123.0;38];assert!(r.block(&mut sentinel).is_err());assert_eq!(sentinel,[-123.0;38]);
    assert!(r.applied_controls().is_empty());assert_eq!(r.system().aperture().state(),initial);
    let mut c=config(37);c.sample_rate_hz=48000;
    assert!(AperturePerformance::new(CoupledAperture::Tube(model()),ApertureObservation::Inlet,phrase(),c).is_err());
    assert!(AperturePerformance::new(CoupledAperture::Tube(model()),ApertureObservation::NetworkNode(0),phrase(),config(37)).is_err());
    let mut c=config(37);c.max_controls=1;
    assert!(AperturePerformance::new(CoupledAperture::Tube(model()),ApertureObservation::Inlet,phrase(),c).is_err());
    let mut c=config(37);c.samples=2;
    assert!(AperturePerformance::new(CoupledAperture::Tube(model()),ApertureObservation::Inlet,phrase(),c).is_err());
    let mut running=model();running.step(TubeDrive{upstream_pressure_pa:5.0,body_flow_m3_s:0.0}).unwrap();
    assert!(AperturePerformance::new(CoupledAperture::Tube(running),ApertureObservation::Inlet,phrase(),config(37)).is_err());
}

#[test]
fn failed_physical_callback_keeps_unapplied_control_and_refuses_continuation() {
    // Canonical sampler admits finite pressure; extreme input is a physical
    // refusal, not a reason for the scheduler to consume the command first.
    let track=GestureTrack{id:"mouth".into(),target:GestureTarget::BlowingPressure,
        initial:GestureValue::PressurePa(f64::MAX),events:vec![]};
    let source=GestureSchedule::try_new(700,vec![track]).unwrap();
    let mut r=AperturePerformance::new(CoupledAperture::Tube(model()),ApertureObservation::Inlet,source,config(37)).unwrap();
    let before=r.system().aperture().state();let mut out=[-123.0;1];
    assert!(r.block(&mut out).is_err());assert_eq!(r.system().aperture().state(),before);
    assert!(r.applied_controls().is_empty());assert_eq!(out,[-123.0]);
    assert!(matches!(r.block(&mut out),Err(RenderError::Poisoned)));
}
