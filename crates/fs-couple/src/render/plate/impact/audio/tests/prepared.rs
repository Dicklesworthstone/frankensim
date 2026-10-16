use super::*;

#[test]
fn prepared_nonlinear_contact_reaches_pressure_decimation_pcm_and_exact_resume() {
    use super::super::super::BodyPotential;
    use crate::modal_acoustic_time::ModalAcousticState;
    use crate::pcm_wav::observation::DecimatedRenderer;
    use crate::pcm_wav::stream::{Pcm16WavStream, render_pressure_pcm16, ScheduledWavProgress};
    use fs_dcontact::Obstacle;
    use std::io::Cursor;
    let a = artifact();
    let make = |contact_enabled: bool| {
        let (striker, b) = ImpactBody::free_mass(0.04, -0.0001, 0.2).unwrap();
        let receiver = ImpactBody { potential: BodyPotential::Linear(vec![core::f64::consts::TAU*300.0]),
            initial: vec![ModalAcousticState::default()], damping_per_s: vec![3.0] };
        let contact = Obstacle::new(vec![b, -1.0], 1, 2, vec![0.0], vec![1.0],
            1e7, 1.5, "synthetic prepared impact, not a measured drum".into()).unwrap();
        let config = ImpactConfig { dt_s:1.0/f64::from(RATE), max_steps:1028,
            maximum_generalized_force:1000.0, maximum_energy_j:10.0,
            energy_absolute_tolerance_j:1e-10, energy_relative_tolerance:1e-7 };
        let mechanics = ImpactSystem::new(vec![striker, receiver],
            if contact_enabled {vec![contact]} else {vec![]}, vec![], vec![], config)
            .unwrap().prepare().unwrap();
        let source = ImpactPressureRenderer::new(mechanics, &a,
            vec![VelocityProjection { input_id: "surface_velocity_m_s".into(), weights: vec![0.0, 1.0] }],
            vec![0.0;2], listener(1.0), RATE, 97).unwrap();
        DecimatedRenderer::new(source, RATE, RATE/4, 43).unwrap()
    };
    let mut baseline = make(true);
    let mut expected = vec![0.0;257];
    for chunk in expected.chunks_mut(17) { baseline.block(chunk).unwrap(); }
    assert!(expected.iter().any(|p|p.abs()>1e-6));
    assert!(baseline.source().mechanics().state()[1] < 0.0, "physical rebound");
    assert_eq!(baseline.source().last_mechanical_frame().unwrap().sample, 1028);
    assert!(baseline.source().last_mechanical_frame().unwrap().balance_residual_j.abs() < 1e-8);
    let mut disconnected = make(false);
    for _ in 0..25 {
        let mut silent = [0.0;10]; disconnected.block(&mut silent).unwrap();
        assert!(silent.iter().all(|p|*p==0.0), "no contact must mean no receiver sound");
    }
    let (expected_wav, clips) = crate::pcm_wav::encode_pcm16_wav(&expected,RATE/4,20.0).unwrap();
    for block in [1,7,43] {
        let mut actual = make(true);
        let mut stream = Pcm16WavStream::new(Cursor::new(Vec::<u8>::new()),RATE/4,20.0,43).unwrap();
        let mut scratch = vec![0.0;block]; let gate=CancelGate::new_clock_free();
        render_pressure_pcm16(&mut actual,&mut stream,&gate,&mut scratch,31).unwrap();
        let saved=actual.source().mechanics().state().to_vec();
        let stopped=CancelGate::new_clock_free();stopped.request();
        assert_eq!(render_pressure_pcm16(&mut actual,&mut stream,&stopped,&mut scratch,226).unwrap(),
            ScheduledWavProgress::Cancelled{samples:0});
        assert_eq!(actual.source().mechanics().state(),saved);
        render_pressure_pcm16(&mut actual,&mut stream,&gate,&mut scratch,226).unwrap();
        let (bytes,summary)=stream.finish().unwrap();
        assert_eq!(bytes.into_inner(),expected_wav);
        assert_eq!(summary.samples,257);assert_eq!(summary.clipped_samples,clips as u64);
        assert_eq!(actual.source().mechanics().state(),baseline.source().mechanics().state());
    }
}

mod cavity;
