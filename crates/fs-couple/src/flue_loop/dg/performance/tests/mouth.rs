use super::*;
use crate::pcm_wav::baffled::openings::RectangularOutletReceiver;

#[test]
fn material_bound_mouth_reaction_and_signed_radiation_share_the_accepted_dg_step() {
    let input = baffled_input()
        .replace("flue 0.010 0.001 0.004 0 0 0.02", "flue 0.010 0.001 0.004 0 0 0")
        .replace("observation baffled 0 0 0.2 8 32 2000", concat!(
            "observation baffled 0 0 0.2 8 32 1000\n",
            "mouth_radiation -0.03 0 8 4\n",
            "mouth_reaction compact-baffled",
        ));
    // One computed outlet fit supplies every replica. Selecting a microphone
    // or host block size must not create a second radiation source.
    let source = FlueDgFilePerformance::from_bytes(
        input.as_bytes(), 37, &CancelGate::new_clock_free(),
    ).unwrap();
    let renderer = source.renderer();
    let receiver = renderer.baffled_receiver_config().unwrap();
    let mouth = renderer.mouth_radiation().unwrap();
    let load = renderer.voice().dg().unwrap().load();
    let exterior = *load.radiation_record().unwrap().gas_record().state();
    assert_eq!(mouth.center_m, [-0.03, 0.0]);
    assert_eq!(mouth.cells, [8, 4]);
    assert_ne!(exterior, *renderer.voice().gas());
    let reaction = renderer.voice().mouth_radiation_load().unwrap();
    assert_eq!(*reaction.exterior(), exterior);
    assert_eq!(reaction.size_m(), [source.geometry().flue_width_m, source.geometry().cut_up_m]);
    let mut expected_outlet = load.baffled_receiver(
        source.info().sample_rate_hz, receiver,
    ).unwrap();
    let mut expected_mouth = BaffledPressure::rectangular_outlet_multirate_differentiating(
        [source.geometry().flue_width_m, source.geometry().cut_up_m],
        source.info().sample_rate_hz, source.info().sample_rate_hz,
        RectangularOutletReceiver {
            center_m: mouth.center_m, position_m: receiver.position_m,
            cells: mouth.cells, maximum_frequency_hz: receiver.maximum_frequency_hz,
        }, (&exterior).into(),
    ).unwrap();
    assert_eq!(renderer.mouth_receiver().unwrap().delay_samples, expected_mouth.delay_samples);
    assert_eq!(renderer.mouth_receiver().unwrap().filter_delay_samples(),
        expected_mouth.filter_delay_samples());
    let support = (
        expected_outlet.delay_samples.0.min(expected_mouth.delay_samples.0),
        expected_outlet.delay_samples.1.max(expected_mouth.delay_samples.1),
    );
    assert_eq!(renderer.propagation_delay_samples(), Some(support));

    let make = |block| {
        let mut config = *renderer.config();
        config.max_block = block;
        FlueDgPerformance::with_baffled_outlet(
            replica_voice(&source), renderer.schedule().clone(), config, receiver,
        ).unwrap()
    };
    let mut whole = make(512).with_mouth_radiation(mouth).unwrap()
        .with_mouth_reaction(&CancelGate::new_clock_free()).unwrap();
    let mut partitioned = make(37).with_mouth_radiation(mouth).unwrap()
        .with_mouth_reaction(&CancelGate::new_clock_free()).unwrap();
    let mut observation_only = make(512).with_mouth_radiation(mouth).unwrap();
    assert!(observation_only.voice().mouth_radiation_load().is_none());

    // These are cold source/geometry refusals, reusing the same computed load.
    for refused in [
        FlueMouthRadiation { center_m: [0.0, 0.0], ..mouth },
        FlueMouthRadiation { cells: [0, 4], ..mouth },
        FlueMouthRadiation { center_m: [f64::NAN, 0.0], ..mouth },
    ] {
        assert!(make(37).with_mouth_radiation(refused).is_err());
    }
    assert!(make(37).with_mouth_reaction(&CancelGate::new_clock_free()).is_err());
    assert!(make(37).with_mouth_radiation(mouth).unwrap()
        .with_mouth_radiation(mouth).is_err());
    assert!(FlueDgPerformance::new(
        replica_voice(&source), renderer.schedule().clone(), *renderer.config(),
    ).unwrap().with_mouth_radiation(mouth).is_err());

    let mut single = source.into_renderer();
    let mut reference = Vec::with_capacity(512);
    let mut mouth_pressure_sq = 0.0;
    for _ in 0..512 {
        let mut pressure = [0.0];
        single.block(&mut pressure).unwrap();
        let duct = single.voice().dg().unwrap();
        let outlet_flow = duct.load().last_terminal_frame().unwrap().outward_flow_m3_s;
        let inward_flow = duct.last_frame().unwrap().inlet_flow_m3_s;
        let terminal = expected_outlet.step_trace(&[outlet_flow]).unwrap();
        let mouth_pressure = expected_mouth.step_trace(&[-inward_flow]).unwrap();
        assert_eq!(pressure[0], terminal + mouth_pressure);
        mouth_pressure_sq += mouth_pressure * mouth_pressure;
        reference.push(pressure[0]);
    }
    assert!(mouth_pressure_sq > 0.0, "the physical mouth contributes nonzero pressure");
    assert!(reference[..support.0].iter().all(|p| *p == 0.0));
    let mut together = vec![0.0; 512];
    let mut chunked = vec![0.0; 512];
    let mut unchanged_load = vec![0.0; 512];
    whole.block(&mut together).unwrap();
    for block in chunked.chunks_mut(37) { partitioned.block(block).unwrap(); }
    observation_only.block(&mut unchanged_load).unwrap();
    assert_eq!(reference, together);
    assert_eq!(reference, chunked);
    assert_eq!(single.voice().dg().unwrap().load().state(), whole.voice().dg().unwrap().load().state());
    assert_eq!(single.voice().dg().unwrap().load().state(), partitioned.voice().dg().unwrap().load().state());
    assert_ne!(whole.voice().dg().unwrap().load().state(), observation_only.voice().dg().unwrap().load().state(),
        "mouth reaction must feed the source/bore solve, not only the output metadata");
    for diagnostic in single.voice().diagnostics() {
        let scale = diagnostic.source_work_j.abs() + diagnostic.pipe_work_j.abs()
            + diagnostic.stored_delta_j.abs() + diagnostic.radiation_loss_j;
        assert!(diagnostic.ledger_defect_j.abs() <= 1e-10 * scale.max(1e-300));
        assert!(diagnostic.radiation_loss_j >= 0.0);
        assert_eq!(diagnostic.numerical_dissipation_j, 0.0);
    }
    assert!(single.with_mouth_radiation(mouth).is_err());
}
