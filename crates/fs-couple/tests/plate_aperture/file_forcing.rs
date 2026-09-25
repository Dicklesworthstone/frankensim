use super::*;
use fs_couple::bernoulli_aperture::dynamic::force::PlateForceFootprint;
use fs_couple::bernoulli_aperture::performance::{AperturePerformance, CoupledAperture};
use fs_couple::bernoulli_aperture::performance::force::{ApertureForceEvent as Event, ApertureForceProgram};
use fs_couple::bernoulli_aperture::network::ApertureNetwork;
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use std::io::Cursor;

const FORCED: &str = include_str!("../../examples/plate-valve-forced.performance");
const REGIONAL: &str = include_str!("../../examples/plate-valve-regional-gas.performance");
fn program() -> ApertureForceProgram {
    ApertureForceProgram { footprints: vec![PlateForceFootprint::Node(13), PlateForceFootprint::Patch(vec![10,11])],
        events: vec![
            Event {sample:0,port:0,force_n:-0.0002}, Event {sample:73,port:0,force_n:-0.001},
            Event {sample:512,port:1,force_n:-0.0008}, Event {sample:1024,port:0,force_n:0.0},
            Event {sample:1200,port:1,force_n:0.0}, Event {sample:4000,port:0,force_n:-0.0005},
            Event {sample:6000,port:0,force_n:0.0},
        ] }
}
fn runtime(base: &str, block: usize) -> AperturePerformance {
    PlateValvePerformance::from_bytes(base.as_bytes(),block,&CancelGate::new()).unwrap().into_renderer()
}
fn direct(source: &AperturePerformance) -> CoupledAperture {
    let a=source.system().aperture(); let spec=a.spec();
    let mut valve=DynamicAperture::from_plate_with_closure(a.plate_reduction().unwrap().clone(),
        a.plate_closure().unwrap().spec().clone(),spec.density_kg_m3,spec.impedance_pa_s_m3,
        spec.time_step_s,spec.max_steps,a.state()).unwrap();
    valve=valve.with_plate_relaxation(a.relaxation().unwrap().spec().clone(),InitialApertureMemory::Relaxed).unwrap();
    match source.system() {
        CoupledAperture::Tube(t)=>CoupledAperture::Tube(ApertureTube::new(valve,*t.spec()).unwrap()),
        CoupledAperture::Network(n)=>CoupledAperture::Network(if let Some(gas)=n.section_gases() {
            ApertureNetwork::with_section_gases(valve,n.spec().clone(),gas.to_vec()).unwrap()
        } else {ApertureNetwork::new(valve,n.spec().clone()).unwrap()}),
    }
}
fn trace(mut r: AperturePerformance, block: usize) -> (Vec<f64>,AperturePerformance) {
    let mut out=vec![0.0;r.config().samples as usize];
    for part in out.chunks_mut(block) {r.block(part).unwrap();}
    (out,r)
}

#[test]
fn scheduled_physical_forces_match_direct_pressure_contact_memory_and_regional_network_steps() {
    for base in [INPUT, REGIONAL] {
        // Ensure the reference samples the same internal port even if the example
        // chooses an exterior receiver by another spelling.
        let base=base.lines().map(|l|if l.starts_with("observation "){"observation inlet"}else{l}).collect::<Vec<_>>().join("\n")+"\n";
        let scheduled=runtime(&base,37).with_plate_forces(program()).unwrap();
        let mut system=direct(&scheduled);
        let ports=scheduled.force_ports().to_vec();let gesture=scheduled.schedule().clone();
        let events=program().events;let mut held=[0.0;2];let mut expected=Vec::new();let mut next=0;
        for n in 0..scheduled.config().samples {
            while next<events.len() && events[next].sample==n {let e=events[next];held[e.port]=e.force_n;next+=1;}
            let force=ports.iter().zip(held).fold(0.0,|sum,(p,f)|sum+p.generalized_force_n(f).unwrap());
            let drive=TubeDrive {upstream_pressure_pa:gesture.sample("mouth",n*700/96000).unwrap(),body_flow_m3_s:0.0};
            let (a,residual,total)=match &mut system {
                CoupledAperture::Tube(t)=>{let f=t.step_with_force(drive,force).unwrap();
                    (f.aperture,f.balance_residual_j(),f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs())},
                CoupledAperture::Network(t)=>{let f=t.step_with_force(drive,force).unwrap();
                    (f.aperture,f.balance_residual_j(),f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs())},
            };
            assert!(residual.abs()<1e-8*(total+a.mechanical_work_j.abs()).max(1e-25));
            expected.push(a.bore_pressure_pa);
        }
        let (actual,r)=trace(scheduled,37);assert_eq!(actual,expected);
        assert_eq!(r.system().aperture().state(),system.aperture().state());
        assert_eq!(r.system().aperture().relaxation().unwrap().memory_sqrt_j(),system.aperture().relaxation().unwrap().memory_sqrt_j());
        assert_eq!(r.applied_force_controls().len(),7);assert!(r.pending_force_controls().is_empty());
        assert!(actual.iter().zip(render(&base)).any(|(a,b)|(a-b).abs()>1e-6));
    }
}

#[test]
fn simultaneous_force_ports_and_half_open_callbacks_preserve_assignment_order_and_unforced_parity() {
    let mut data=program();data.events=vec![
        Event {sample:37,port:0,force_n:-0.0004}, Event {sample:0,port:1,force_n:-0.0001},
        Event {sample:37,port:0,force_n:-0.0003}, Event {sample:37,port:1,force_n:-0.0002},
        Event {sample:38,port:0,force_n:0.0},
    ];
    let mut r=runtime(INPUT,37).with_plate_forces(data.clone()).unwrap();
    let b=r.force_ports().iter().map(|p|p.coefficient()).collect::<Vec<_>>();
    let mut prefix=[0.0;37];r.block(&mut prefix).unwrap();
    assert_eq!(r.applied_force_controls().len(),1);assert_eq!(r.pending_force_controls()[0].sample,37);
    assert_eq!(r.held_generalized_force_n(),b[1]*-0.0001);
    r.block(&mut prefix[..1]).unwrap();assert_eq!(r.applied_force_controls().len(),4);
    assert_eq!(r.held_generalized_force_n(),b[0]*-0.0003+b[1]*-0.0002);
    r.block(&mut prefix[..1]).unwrap();assert_eq!(r.held_generalized_force_n(),b[1]*-0.0002);
    let expected=trace(runtime(INPUT,37).with_plate_forces(data.clone()).unwrap(),37).0;
    for block in [1,512] {assert_eq!(trace(runtime(INPUT,block).with_plate_forces(data.clone()).unwrap(),block).0,expected);}
    let mut zero=program();for e in &mut zero.events {e.force_n=0.0;}
    assert_eq!(trace(runtime(INPUT,37).with_plate_forces(zero).unwrap(),37).0,render(INPUT));
}

#[test]
fn force_file_and_decimated_cancellation_keep_contact_material_wave_and_control_history() {
    let source=load(FORCED);assert_eq!((source.info().force_ports,source.info().force_controls),(2,7));
    let mut source=DecimatedRenderer::new(source.into_renderer(),96000,48000,37).unwrap();
    let mut stream=Pcm16WavStream::new(Cursor::new(Vec::new()),48000,1.0,37).unwrap();
    let mut scratch=[0.0;37];
    render_pressure_pcm16(&mut source,&mut stream,&CancelGate::new(),&mut scratch,300).unwrap();
    let before=source.source().system().aperture().state();
    let memory=source.source().system().aperture().relaxation().unwrap().memory_sqrt_j().to_vec();
    let controls=source.source().applied_force_controls().to_vec();
    let gate=CancelGate::new();gate.request();scratch.fill(1234.0);
    assert_eq!(render_pressure_pcm16(&mut source,&mut stream,&gate,&mut scratch,4501).unwrap(),ScheduledWavProgress::Cancelled {samples:0});
    assert_eq!(scratch,[1234.0;37]);assert_eq!(source.source().system().aperture().state(),before);
    assert_eq!(source.source().applied_force_controls(),controls.as_slice());
    assert_eq!(source.source().system().aperture().relaxation().unwrap().memory_sqrt_j(),memory.as_slice());
    render_pressure_pcm16(&mut source,&mut stream,&CancelGate::new(),&mut scratch,4501).unwrap();
    let (out,_)=stream.finish().unwrap();
    let mut whole=DecimatedRenderer::new(runtime(FORCED,37),96000,48000,37).unwrap();
    let mut expected=vec![0.0;4801];for part in expected.chunks_mut(37){whole.block(part).unwrap();}
    let (wav,_)=fs_couple::pcm_wav::encode_pcm16_wav(&expected,48000,1.0).unwrap();
    assert_eq!(out.into_inner(),wav);assert_eq!(source.source().pending_force_controls().len(),0);
}

#[test]
fn invalid_force_input_refuses_and_late_physics_failure_never_consumes_pending_controls() {
    for text in [
        FORCED.replace("force_ports 2","force_ports 33"),FORCED.replace("force_port node 13","force_port node 999"),
        FORCED.replace("force_port patch 2 10 11","force_port patch 2 10 10"),
        FORCED.replace("force_event 73 0 -0.001","force_event 9602 0 -0.001"),
        FORCED.replace("force_event 73 0 -0.001","force_event 73 2 -0.001"),
        FORCED.replace("force_event 73 0 -0.001","force_event 73 0 NaN"),
        FORCED.replace("force_events 7","force_events 262145"),
        FORCED.replace("compile_limits 10000 128","compile_limits 10000 6"),
    ] {assert!(PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).is_err());}
    let mut data=program();data.events=vec![Event {sample:73,port:0,force_n:1e100}];
    let mut r=runtime(INPUT,37).with_plate_forces(data).unwrap();
    let mut out=[0.0;37];r.block(&mut out).unwrap();r.block(&mut out[..36]).unwrap();
    let state=r.system().aperture().state();out.fill(1234.0);
    assert!(r.block(&mut out[..1]).is_err());assert_eq!(r.system().aperture().state(),state);
    assert!(r.applied_force_controls().is_empty());assert_eq!(r.pending_force_controls().len(),1);
    assert!(matches!(r.block(&mut out[..1]),Err(fs_couple::render::RenderError::Poisoned)));
    assert_eq!(out,[1234.0;37]);
    let mut advanced=runtime(INPUT,37);advanced.block(&mut [0.0]).unwrap();
    assert!(advanced.with_plate_forces(program()).is_err());
    assert!(runtime(INPUT,37).with_plate_forces(program()).unwrap().with_plate_forces(program()).is_err());
}
