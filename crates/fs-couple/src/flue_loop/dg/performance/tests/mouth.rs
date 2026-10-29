use super::*;
use crate::pcm_wav::baffled::openings::RectangularOutletReceiver;

#[test]
fn shared_observer_gas_covers_mouth_and_terminals_without_changing_the_source() {
    use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationLimits};
    let gate = CancelGate::new_clock_free();
    let input = baffled_input().replace("observation baffled 0 0 0.2 8 32 2000",
        concat!("observation baffled 0 0 0.2 8 32 2000\n",
            "observer_gas stokes-kirchhoff 4096 1000000 50000000\n",
            "mouth_radiation -0.03 0 4 2"));
    let source = FlueDgFilePerformance::from_bytes(input.as_bytes(), 37, &gate).unwrap();
    let renderer = source.renderer();
    let report = renderer.observer_gas().unwrap();
    let mouth_report = renderer.mouth_observer_gas().unwrap();
    let load = renderer.voice().dg().unwrap().load();
    assert_eq!(report.gas, *load.radiation_record().unwrap().gas_record().state());
    assert_eq!(report.gas, mouth_report.gas);
    assert_ne!(report.gas, *renderer.voice().gas());
    assert_eq!(report.additional_latency_samples, mouth_report.additional_latency_samples);
    assert_eq!(renderer.baffled_receiver().unwrap().filter_delay_samples(),
        renderer.mouth_receiver().unwrap().filter_delay_samples());
    let make = |block| {
        let mut config = *renderer.config();
        config.max_block = block;
        FlueDgPerformance::with_baffled_outlet(replica_voice(&source),
            renderer.schedule().clone(), config, renderer.baffled_receiver_config().unwrap())
            .unwrap().with_mouth_radiation(renderer.mouth_radiation().unwrap()).unwrap()
    };
    let mut original = make(512);
    let mut partitioned = make(37).with_observer_gas(
        &GasPropagation::classical(GasPropagationLimits::surface_reference()), &gate,
    ).unwrap();
    let mut authored = source.into_renderer();
    let (mut expected, mut whole, mut changed) = (vec![0.0;512], vec![0.0;512], vec![0.0;512]);
    original.block(&mut expected).unwrap();
    for block in whole.chunks_mut(37) { authored.block(block).unwrap(); }
    for block in changed.chunks_mut(37) { partitioned.block(block).unwrap(); }
    assert_eq!(whole, changed);
    assert!(whole.iter().zip(&expected).any(|(a,b)| a != b));
    assert_eq!(authored.voice().dg().unwrap().load().state(), original.voice().dg().unwrap().load().state());
    assert_eq!(authored.voice().dg().unwrap().load().state(), partitioned.voice().dg().unwrap().load().state());
    assert!(authored.with_observer_gas(
        &GasPropagation::classical(GasPropagationLimits::surface_reference()), &gate,
    ).is_err());
}

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


#[test]
fn joint_mouth_reaction_reaches_the_source_and_one_shared_physical_microphone() {
    use crate::flue_loop::dg::FlueDgSpec;
    let gate = CancelGate::new_clock_free();
    let text = baffled_input()
        .replace("flue 0.010 0.001 0.004 0 0 0.02", "flue 0.002 0.0005 0.001 0 0 0");
    // A small manufactured integration geometry. The original full-band fit,
    // passivity and source limits remain hard gates, pending native execution.
    let source = FlueDgFilePerformance::from_bytes(text.as_bytes(), 512, &gate).unwrap();
    let original = source.renderer();
    let receiver = original.baffled_receiver_config().unwrap();
    let mouth = FlueMouthRadiation { center_m: [-0.003, 0.0], cells: [8, 8], aperture_size_m: None };
    let base = || FlueDgPerformance::with_baffled_outlet(
        replica_voice(&source), original.schedule().clone(), *original.config(), receiver,
    ).unwrap().with_mouth_radiation(mouth).unwrap();
    let mut observation_only = base();
    let joint = base().with_mutual_mouth_reaction(&gate).unwrap();
    let load = joint.voice().dg().unwrap().load();
    let record = load.radiation_record().unwrap();
    // The selected fit is passive either unchanged or after an M2.5-W3b residue
    // repair inside a tenth of the allowance (order 4 alone, at origin
    // f9d059919, was -0.31 at 332 Hz against a data scale of 1.2e8).
    let report = record.report();
    let allowance = record.options().maximum_relative_error;
    let selected = report.rungs.iter().find(|r| r.order == report.selected_order).unwrap();
    assert!(selected.repair_relative_change.map_or(selected.continuous_passive,
        |change| !selected.continuous_passive && change <= 0.1 * allowance), "{selected:?}");
    assert!(report.continuous.passive && report.sampled.passive);
    assert!(report.maximum_fit_relative_error <= allowance);
    assert_eq!(record.driven_mouth().unwrap().center_m, mouth.center_m);
    assert_eq!(record.driven_mouth().unwrap().size_m, [0.002, 0.001]);
    assert_eq!(record.gas_record().state(), original.voice().dg().unwrap().load()
        .radiation_record().unwrap().gas_record().state());
    assert_ne!(record.realization_identity(), original.voice().dg().unwrap().load()
        .radiation_record().unwrap().realization_identity());
    assert_eq!(load.terminal().unwrap().ports(), load.terminal_ports().len() + 1);
    assert!(joint.voice().mouth_radiation_load().is_none());
    let exterior = *record.gas_record().state();
    let mut expected_outlet = load.baffled_receiver(source.info().sample_rate_hz, receiver).unwrap();
    let mut expected_mouth = BaffledPressure::rectangular_outlet_multirate_differentiating(
        [0.002, 0.001], source.info().sample_rate_hz, source.info().sample_rate_hz,
        RectangularOutletReceiver { center_m: mouth.center_m,
            position_m: receiver.position_m, cells: mouth.cells,
            maximum_frequency_hz: receiver.maximum_frequency_hz }, (&exterior).into(),
    ).unwrap();
    // Reconstructing a voice from the already joint-loaded cold source must
    // retain its physical operator and must not restore legacy mouth inertia.
    let make = |block| {
        let voice = FlueVoice::from_dg(Some(source.card()), 8.0, 1e-6, source.geometry(),
            FlueDgSpec { load: load.clone(), minimum_jet_speed_m_s: 1.0,
                maximum_jet_mach: 0.3, max_jet_history_samples: 8192 },
            source.info().sample_rate_hz, 3, 0, &gate,
        ).unwrap().with_viscous_windway(0.005).unwrap();
        let mut config = *joint.config();
        config.max_block = block;
        FlueDgPerformance::with_baffled_outlet(voice, joint.schedule().clone(), config, receiver)
            .unwrap().with_mouth_radiation(mouth).unwrap()
    };
    let mut single = make(1);
    let mut partitioned = make(37);
    assert!(make(1).with_mouth_reaction(&gate).is_err());
    assert!(make(1).with_mutual_mouth_reaction(&gate).is_err());
    let mut together = joint;
    let mut expected = Vec::with_capacity(512);
    let mut pressure_peak = 0.0_f64;
    let mut mouth_work = 0.0;
    for _ in 0..512 {
        let mut output = [0.0];
        single.block(&mut output).unwrap();
        let frame = single.voice().dg().unwrap().last_frame().unwrap();
        let outward: Vec<_> = frame.terminals.iter().map(|f| f.outward_flow_m3_s).collect();
        let at_outlet = expected_outlet.step_trace(&outward).unwrap();
        let at_mouth = expected_mouth.step_trace(&[-frame.inlet_flow_m3_s]).unwrap();
        assert_eq!(output[0], at_outlet + at_mouth);
        let exchange = frame.driven_mouth.unwrap();
        assert_eq!(exchange.outward_flow_m3_s, -frame.inlet_flow_m3_s);
        assert_eq!(exchange.source_pressure_pa,
            frame.inlet_pressure_pa - exchange.mean_pressure_pa);
        mouth_work += exchange.load_work_j;
        pressure_peak = pressure_peak.max(output[0].abs());
        expected.push(output[0]);
    }
    let mut whole = vec![0.0; 512];
    let mut chunks = vec![0.0; 512];
    together.block(&mut whole).unwrap();
    for chunk in chunks.chunks_mut(37) { partitioned.block(chunk).unwrap(); }
    observation_only.block(&mut vec![0.0; 512]).unwrap();
    assert!(pressure_peak > 0.0);
    assert_eq!(expected, whole);
    assert_eq!(expected, chunks);
    assert_eq!(single.voice().dg().unwrap().load().state(), together.voice().dg().unwrap().load().state());
    assert_eq!(single.voice().dg().unwrap().load().terminal().unwrap().pole_states(),
        partitioned.voice().dg().unwrap().load().terminal().unwrap().pole_states());
    assert_ne!(together.voice().dg().unwrap().load().state(),
        observation_only.voice().dg().unwrap().load().state());
    let reported_work: f64 = together.voice().diagnostics().iter()
        .map(|diagnostic| diagnostic.joint_mouth_work_j).sum();
    let work_scale: f64 = single.voice().diagnostics().iter()
        .map(|diagnostic| diagnostic.joint_mouth_work_j.abs()).sum();
    assert!((reported_work-mouth_work).abs() <= 2e-12*work_scale.max(1e-300));
    for diagnostic in together.voice().diagnostics() {
        assert_eq!(diagnostic.radiation_loss_j, 0.0);
        assert_eq!(diagnostic.numerical_dissipation_j, 0.0);
        let scale = diagnostic.source_work_j.abs() + diagnostic.pipe_work_j.abs()
            + diagnostic.joint_mouth_work_j.abs();
        assert!(diagnostic.ledger_defect_j.abs() <= 1e-10 * scale.max(1e-300));
    }
}
