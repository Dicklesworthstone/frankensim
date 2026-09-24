//! Physical boundary feedback, not pressure-to-microphone gain tests.
use super::*;
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::pcm_wav::baffled::{BaffledPressure,RayleighMedium};
use fs_vfit::relaxation::RelaxationImpedance;

fn load() -> BaffledRadiationLoad {
    BaffledRadiationLoad::new(uniform().radius_m,1.2,uniform().sound_speed_m_s,
        1.0/f64::from(RATE),1500.0,&CancelGate::new()).unwrap()
}
fn network(load:BaffledRadiationLoad) -> ApertureNetwork {
    let s=uniform();
    ApertureNetwork::new(valve(),TubeNetworkSpec {
        nodes:vec![NetworkNode::Inlet,load.termination()],
        sections:vec![TubeSection{nodes:[0,1],length_m:s.length_m,radius_m:s.radius_m,
            max_length_error_m:s.max_length_error_m}],
        sound_speed_m_s:s.sound_speed_m_s,max_wave_memory_bytes:1<<20,
    }).unwrap()
}
fn observation() -> ApertureObservation {
    let mut receiver=receiver_location();receiver.maximum_frequency_hz=1500.0;
    ApertureObservation::NetworkBaffled{node:1,receiver}
}
fn scheduled(block:usize,load:BaffledRadiationLoad) -> AperturePerformance {
    AperturePerformance::new(CoupledAperture::Network(network(load)),observation(),phrase(),config(block)).unwrap()
}

#[test]
fn radiation_geometry_fixes_both_reactive_storage_and_frequency_dependent_loss() {
    let a=load();let branch=a.load();let term=branch.terms()[0];
    let alpha=8.0/(3.0*core::f64::consts::PI);
    let z=uniform().characteristic_impedance(1.2).unwrap();
    assert!((term.resistance_pa_s_m3/z-2.0*alpha*alpha).abs()<1e-14);
    let inertance=term.resistance_pa_s_m3/term.rate_per_s;
    assert!((inertance/(z*alpha*a.radius_m()/a.sound_speed_m_s())-1.0).abs()<1e-14);
    assert_eq!(a.impedance_at(0.0,true).unwrap(),(0.0,0.0));
    // Independent frequencies, neither construction endpoints nor its 32-row grid.
    for f in [71.0,271.0,823.0,1421.0] {
        let (r,x)=fs_phs::baffled_piston_impedance(1.2,a.sound_speed_m_s(),a.radius_m(),
            core::f64::consts::TAU*f,10).unwrap();
        for discrete in [true,false] {
            let (rr,xx)=a.impedance_at(f,discrete).unwrap();
            assert!((rr-r).hypot(xx-x)/r.hypot(x)<0.05);
            assert!((rr-r).abs()/r<0.05);
        }
    }
    let denser=BaffledRadiationLoad::new(a.radius_m(),2.4,a.sound_speed_m_s(),a.time_step_s(),1500.0,&CancelGate::new()).unwrap();
    let d=denser.load();assert_eq!(d.terms()[0].rate_per_s,term.rate_per_s);
    assert_eq!(d.terms()[0].resistance_pa_s_m3,2.0*term.resistance_pa_s_m3);
    assert!(a.impedance_at(1501.0,false).is_err());
    assert!(BaffledRadiationLoad::new(a.radius_m(),1.2,a.sound_speed_m_s(),a.time_step_s(),10000.0,&CancelGate::new()).is_err());
    assert!(BaffledRadiationLoad::new(a.radius_m(),1.2,a.sound_speed_m_s(),0.001,1500.0,&CancelGate::new()).is_err());
    let gate=CancelGate::new();gate.request();
    assert!(BaffledRadiationLoad::new(a.radius_m(),1.2,a.sound_speed_m_s(),a.time_step_s(),1500.0,&gate).is_err());
}

#[test]
fn actual_boundary_samples_match_the_checked_discrete_impedance_not_the_analogue_curve() {
    let l=load();let f=1000.0;let w=core::f64::consts::TAU*f;
    let mut owner=RelaxationImpedance::new(l.load(),uniform().characteristic_impedance(1.2).unwrap(),l.time_step_s()).unwrap();
    let(mut pc,mut ps,mut qc,mut qs)=(0.0,0.0,0.0,0.0);
    for n in 0..4096 {
        let phase=w*f64::from(n)*l.time_step_s();
        let (s,c)=phase.sin_cos();
        let frame=owner.step(c).unwrap();
        assert!(frame.port.dissipated_energy_j>=0.0);
        if n>=1024 {
            pc+=frame.port.pressure_pa*c;ps+=frame.port.pressure_pa*s;
            qc+=frame.port.flow_m3_s*c;qs+=frame.port.flow_m3_s*s;
        }
    }
    // 3072 retained samples = 32 complete cycles at the 96 kHz source clock.
    let den=qc*qc+qs*qs;
    let measured=((pc*qc+ps*qs)/den,(ps*qc-pc*qs)/den);
    let expected=l.impedance_at(f,true).unwrap();
    let err=(measured.0-expected.0).hypot(measured.1-expected.1)/expected.0.hypot(expected.1);
    assert!(err<1e-11,"actual discrete boundary response discrepancy {err:e}");
}

#[test]
fn radiation_reacts_after_the_round_trip_and_keeps_its_storage_out_of_the_loss_ledger() {
    let l=load();let mut loaded=network(l);let mut resistive=model();
    let delay=loaded.represented_sections()[0].one_way_samples;
    let source=phrase();let(mut changed,mut stored,mut dissipated)=(false,false,false);
    for n in 0..SAMPLES {
        let drive=TubeDrive {upstream_pressure_pa:source.sample("mouth",n*700/u64::from(RATE)).unwrap(),body_flow_m3_s:0.0};
        let before=loaded.stored_energy_j();
        let f=loaded.step(drive).unwrap();let old=resistive.step(drive).unwrap();
        let port=loaded.node_frame(1).unwrap();
        if n<2*delay as u64 {assert_eq!(f.aperture.state,old.aperture.state);}
        else {changed|=f.aperture.state!=old.aperture.state;}
        stored|=port.stored_energy_j>0.0;dissipated|=port.absorbed_energy_j>0.0;
        let work=port.pressure_pa*port.net_flow_into_node_m3_s*l.time_step_s();
        let scale=before+f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs()+f.body_work_j.abs();
        assert!(f.balance_residual_j().abs()<=3e-10*scale.max(f64::MIN_POSITIVE));
        assert!((port.storage_change_j+port.absorbed_energy_j-work).abs()<=3e-10*scale.max(f64::MIN_POSITIVE));
    }
    assert!(changed && stored && dissipated,"radiation must store, dissipate and feed back, not only change observations");
}

#[test]
fn loaded_outlet_receiver_uses_the_same_accepted_flow_without_changing_the_feedback() {
    let l=load();let mut direct=network(l);
    let ApertureObservation::NetworkBaffled{receiver,..}=observation() else {unreachable!()};
    let mut mic=BaffledPressure::circular_outlet(l.radius_m(),RATE,receiver,
        RayleighMedium {density:1.2,sound_speed:l.sound_speed_m_s()}).unwrap();
    let source=phrase();let mut expected=Vec::new();
    for n in 0..SAMPLES {
        direct.step(TubeDrive {upstream_pressure_pa:source.sample("mouth",n*700/u64::from(RATE)).unwrap(),body_flow_m3_s:0.0}).unwrap();
        expected.push(mic.step(&[direct.node_frame(1).unwrap().net_flow_into_node_m3_s]).unwrap());
    }
    assert!(expected.iter().any(|x|x.abs()>1e-10));
    for block in [1,37,512] {
        let mut source=scheduled(block,l);let mut output=vec![0.0;SAMPLES as usize];
        for chunk in output.chunks_mut(block) {source.block(chunk).unwrap();}
        assert_eq!(bits(&output),bits(&expected));
        let CoupledAperture::Network(m)=source.system() else {unreachable!()};
        assert_eq!(m.node_frame(1),direct.node_frame(1));
        assert_eq!(m.aperture().state(),direct.aperture().state());
        assert_eq!(m.stored_energy_j().to_bits(),direct.stored_energy_j().to_bits());
    }
    let bad=ApertureObservation::NetworkBaffled{node:0,receiver};
    assert!(AperturePerformance::new(CoupledAperture::Network(network(l)),bad,phrase(),config(37)).is_err());
}

#[test]
fn failed_physics_and_cancelled_audio_preserve_radiation_material_and_receiver_history() {
    let l=load();let mut a=network(l);let mut b=network(l);
    let drive=TubeDrive{upstream_pressure_pa:5.0,body_flow_m3_s:0.0};
    for _ in 0..256 {a.step(drive).unwrap();b.step(drive).unwrap();}
    let energy=a.stored_energy_j().to_bits();let node=*a.node_frame(1).unwrap();
    assert!(a.step(TubeDrive{upstream_pressure_pa:f64::NAN,..drive}).is_err());
    assert_eq!(a.stored_energy_j().to_bits(),energy);assert_eq!(*a.node_frame(1).unwrap(),node);
    for _ in 0..256 {assert_eq!(a.step(drive).unwrap(),b.step(drive).unwrap());}
    let mut expected=scheduled(37,l);let mut actual=scheduled(37,l);
    let mut uninterrupted=vec![0.0;SAMPLES as usize];
    for c in uninterrupted.chunks_mut(37) {expected.block(c).unwrap();}
    let mut pcm=Pcm16WavStream::new(Cursor::new(Vec::new()),RATE,0.01,37).unwrap();
    let mut scratch=[0.0;37];
    render_pressure_pcm16(&mut actual,&mut pcm,&CancelGate::new(),&mut scratch,271).unwrap();
    let before=actual.system().aperture().state();let pending=actual.pending_controls().to_vec();
    let gate=CancelGate::new();gate.request();scratch.fill(-123.0);
    assert_eq!(render_pressure_pcm16(&mut actual,&mut pcm,&gate,&mut scratch,SAMPLES as u64-271).unwrap(),ScheduledWavProgress::Cancelled{samples:0});
    assert_eq!(scratch,[-123.0;37]);assert_eq!(actual.system().aperture().state(),before);
    assert_eq!(actual.pending_controls(),pending);
    render_pressure_pcm16(&mut actual,&mut pcm,&CancelGate::new(),&mut scratch,SAMPLES as u64-271).unwrap();
    assert_eq!(pcm.finish().unwrap().0.into_inner(),fs_couple::pcm_wav::encode_pcm16_wav(&uninterrupted,RATE,0.01).unwrap().0);
}
