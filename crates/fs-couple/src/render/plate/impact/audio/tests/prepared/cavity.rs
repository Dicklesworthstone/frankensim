//! Integration fixture: a struck receiver drives the original analytic acoustic
//! observer. The cavity must affect sound through receiver motion, not a new gain.
use super::super::*;
use crate::render::plate::impact::{BodyPotential,cavity::CavityCoupling};
use crate::vibroacoustic::CavityModes;
use crate::modal_acoustic_time::ModalAcousticState;
use crate::pcm_wav::observation::DecimatedRenderer;
use crate::pcm_wav::stream::{Pcm16WavStream,render_pressure_pcm16,ScheduledWavProgress};
use fs_dcontact::Obstacle;
use std::io::Cursor;

#[test]
fn cavity_feedback_changes_observed_pressure_without_radiating_acoustic_coordinates_directly() {
    let radiation=artifact();
    let make=|standing:bool,connected:bool| {
        let (striker,weight)=ImpactBody::free_mass(0.04,-0.0001,0.2).unwrap();
        let receiver=ImpactBody {potential:BodyPotential::Linear(vec![core::f64::consts::TAU*300.0]),
            initial:vec![ModalAcousticState::default()],damping_per_s:vec![3.0]};
        let contact=Obstacle::new(vec![weight,-1.0],1,2,vec![0.0],vec![1.0],1e7,1.5,
            "synthetic cavity/audio integration; not a measured instrument".into()).unwrap();
        let count=if standing {2}else{1};
        let cavity=CavityModes {omegas:if standing {vec![0.0,core::f64::consts::TAU*700.0]}else{vec![0.0]},
            lambdas:vec![0.01;count],interface:vec![vec![1.0];count],loss_factor:0.0,rho0:1.2,c0:343.0};
        let mut coupling=vec![0.0;2*count];coupling[count..].fill(0.04);
        let config=ImpactConfig {dt_s:1.0/f64::from(RATE),max_steps:1028,maximum_energy_j:10.0,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1000.0};
        let (system,air)=CavityCoupling::new(&cavity,2,&coupling,&vec![0.0;count]).unwrap()
            .build(vec![striker,receiver],if connected {vec![contact]}else{vec![]},vec![],config,
                &CancelGate::new_clock_free()).unwrap();
        let n=air.total_modes();let mut weights=vec![0.0;n];weights[1]=1.0;
        let source=ImpactPressureRenderer::new(system.prepare().unwrap(),&radiation,
            vec![VelocityProjection {input_id:"surface_velocity_m_s".into(),weights}],
            vec![0.0;n],listener(1.0),RATE,97).unwrap();
        DecimatedRenderer::new(source,RATE,RATE/4,43).unwrap()
    };
    let mut baseline=make(true,true);let mut compact=make(false,true);
    let mut actual=vec![0.0;257];let mut old=vec![0.0;257];
    for (a,b) in actual.chunks_mut(17).zip(old.chunks_mut(17)) {baseline.block(a).unwrap();compact.block(b).unwrap();}
    assert!(actual.iter().any(|v|v.abs()>1e-6));
    assert!(actual.iter().zip(&old).any(|(a,b)|(a-b).abs()>1e-6),"standing-wave feedback must change receiver sound");
    assert_ne!(baseline.source().mechanics().state()[5],0.0,"cavity inertia stores history");
    assert!(baseline.source().last_mechanical_frame().unwrap().balance_residual_j.abs()<1e-8);
    let mut disconnected=make(true,false);
    for _ in 0..25 {let mut zero=[0.0;10];disconnected.block(&mut zero).unwrap();assert_eq!(zero,[0.0;10]);}
    let (expected,clips)=crate::pcm_wav::encode_pcm16_wav(&actual,RATE/4,20.0).unwrap();
    for block in [1,43] {
        let mut resumed=make(true,true);let gate=CancelGate::new_clock_free();
        let mut wav=Pcm16WavStream::new(Cursor::new(Vec::<u8>::new()),RATE/4,20.0,43).unwrap();
        let mut scratch=vec![0.0;block];
        render_pressure_pcm16(&mut resumed,&mut wav,&gate,&mut scratch,31).unwrap();
        let before=resumed.source().mechanics().state().to_vec();let stop=CancelGate::new_clock_free();stop.request();
        assert_eq!(render_pressure_pcm16(&mut resumed,&mut wav,&stop,&mut scratch,226).unwrap(),
            ScheduledWavProgress::Cancelled {samples:0});
        assert_eq!(resumed.source().mechanics().state(),before);
        render_pressure_pcm16(&mut resumed,&mut wav,&gate,&mut scratch,226).unwrap();
        let (bytes,report)=wav.finish().unwrap();assert_eq!(bytes.into_inner(),expected);
        assert_eq!(report.clipped_samples,clips as u64);
        assert_eq!(resumed.source().mechanics().state(),baseline.source().mechanics().state());
    }
}
