use super::*;
use crate::bernoulli_aperture::contraction::ContractionSpec;
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_exec::CancelGate;

fn articulated() -> GestureSchedule {
    // Reversed source order deliberately exercises simultaneous assignments.
    GestureSchedule::try_new(700, vec![
        GestureTrack { id: "angle".into(), target: GestureTarget::JetAngle,
            initial: GestureValue::AngleRad(0.0), events: vec![
                GestureEvent { time_s: 1.0/700.0, transition_s: 4.0/700.0,
                    value: GestureValue::AngleRad(0.06) },
                GestureEvent { time_s: 3.0/700.0, transition_s: 2.0/700.0,
                    value: GestureValue::AngleRad(-0.02) },
            ] },
        GestureTrack { id: "blow".into(), target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(300.0), events: vec![
                GestureEvent { time_s: 2.0/700.0, transition_s: 2.0/700.0,
                    value: GestureValue::PressurePa(400.0) },
                GestureEvent { time_s: 6.0/700.0, transition_s: 0.0,
                    value: GestureValue::PressurePa(0.0) },
            ] },
    ]).unwrap()
}

#[test]
fn sampled_jet_angle_matches_physical_source_and_receiver_on_the_exact_clock() {
    let schedule = articulated();
    let g = gas();
    let spec = source(g, -0.9);
    let cfg = config(1);
    let mut original = voice(spec.clone());
    let mut receiver = BaffledPressure::circular_outlet(0.006, 48_000, cfg.receiver,
        RayleighMedium { density: g.density, sound_speed: g.sound_speed }).unwrap();
    let mut expected = vec![0.0; cfg.samples as usize];
    for (sample, pressure) in expected.iter_mut().enumerate() {
        let tick = (sample as u64 * 700)/48_000;
        original.apply(FlueControl::SetBlowingPressure(schedule.sample("blow", tick).unwrap())).unwrap();
        original.apply(FlueControl::SetJetAngle(schedule.sample("angle", tick).unwrap())).unwrap();
        original.step_block_observed(&mut [0.0], |line| {
            let FlueResonator::Network(n) = line else { panic!("network") };
            *pressure = receiver.step(&[n.node_frame(1).unwrap().net_flow_into_node_m3_s]).unwrap();
            Ok(())
        }).unwrap();
    }
    assert!(expected.iter().any(|p| p.abs() > 1.0e-12));
    for partition in [513, 37, 1] {
        let mut performance = FlueOutletPerformance::new(voice(spec.clone()), schedule.clone(), cfg).unwrap();
        let mut actual = vec![0.0; cfg.samples as usize];
        for block in actual.chunks_mut(partition) { performance.block(block).unwrap(); }
        assert_eq!(actual, expected, "source controls must not depend on callbacks");
        assert_eq!(performance.voice().jet_angle_rad(), -0.02);
        assert_eq!(performance.voice().jet_speed_m_s(), 0.0);
        assert!(performance.voice().network().unwrap().stored_energy_j() > 0.0);
        for frame in performance.voice().diagnostics() {
            assert!(frame.ledger_defect_j.abs() < 1.0e-12,
                "the changed jet split must enter the original source-work balance");
        }
    }
    let pressure_only = GestureSchedule::try_new(700, vec![schedule.tracks()[1].clone()]).unwrap();
    let mut fixed = FlueOutletPerformance::new(voice(spec), pressure_only, cfg).unwrap();
    let mut fixed_pressure = vec![0.0; cfg.samples as usize];
    fixed.block(&mut fixed_pressure).unwrap();
    assert_ne!(expected, fixed_pressure, "jet direction must reach actual microphone pressure");
}

#[test]
fn angle_control_is_static_geometry_equivalent_and_preserves_moving_gas() {
    let gate = CancelGate::new_clock_free();
    let spec = source(gas(), -0.9);
    let card = mint_tonal_interim_card().unwrap();
    let mut geometry = crate::flue_loop::network::tests::geometry();
    geometry.jet_angle_rad = 0.025;
    let mut geometric = FlueVoice::from_network(Some(&card), 90.0, 1.0e-6,
        geometry, spec.clone(), 48_000, 3, 0, &gate).unwrap();
    let mut controlled = voice(spec.clone());
    controlled.apply(FlueControl::SetJetAngle(0.025)).unwrap();
    for v in [&mut geometric, &mut controlled] {
        v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    }
    let mut a = [0.0; 256];
    let mut b = [0.0; 256];
    geometric.step_block(&mut a).unwrap();
    controlled.step_block(&mut b).unwrap();
    assert_eq!(a, b);

    let loss = ContractionSpec { area_ratio: 0.85, maximum_contracted_mach: 0.3,
        provenance: "manufactured-fixed-window".into() };
    let mut first = voice(spec.clone()).with_mouth_contraction(loss.clone()).unwrap()
        .with_unsteady_windway(0.015, 16, 32).unwrap();
    let mut second = voice(spec).with_mouth_contraction(loss).unwrap()
        .with_unsteady_windway(0.015, 16, 32).unwrap();
    for v in [&mut first, &mut second] {
        v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
        v.step_block(&mut [0.0; 256]).unwrap();
    }
    let velocity = first.unsteady_windway().unwrap().channel().velocities_m_s().to_vec();
    let stored = first.network().unwrap().stored_energy_j();
    let flow = first.mouth_flow_m3_s();
    let count = first.samples_rendered();
    assert!(velocity.iter().any(|v| *v > 0.0));
    assert!(stored > 0.0);
    first.apply(FlueControl::SetJetAngle(-0.03)).unwrap();
    assert_eq!(first.samples_rendered(), count);
    assert_eq!(first.network().unwrap().stored_energy_j(), stored);
    assert_eq!(first.mouth_flow_m3_s(), flow);
    assert_eq!(first.unsteady_windway().unwrap().channel().velocities_m_s(), velocity);
    assert_eq!(first.mouth_contraction(), second.mouth_contraction());
    // A pair of controls with no elapsed physical time cannot reset convection,
    // source filtering, stored waves, fluid momentum or the random stream.
    first.apply(FlueControl::SetJetAngle(0.0)).unwrap();
    first.step_block(&mut a).unwrap();
    second.step_block(&mut b).unwrap();
    assert_eq!(a, b);
    for angle in [f64::NAN, f64::INFINITY, core::f64::consts::FRAC_PI_2] {
        assert!(first.apply(FlueControl::SetJetAngle(angle)).is_err());
        assert_eq!(first.jet_angle_rad(), 0.0);
        assert_eq!(first.samples_rendered(), 512);
    }
}

#[test]
fn whole_angle_phrase_and_work_are_admitted_before_motion() {
    let schedule = articulated();
    let mut cfg = config(1);
    cfg.max_controls = 15; // eight ticks for each of the two declared tracks
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule.clone(), cfg).is_err());
    cfg.max_controls = 16;
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule.clone(), cfg).is_ok());
    cfg.max_compile_work = 1;
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule.clone(), cfg).is_err());
    cfg = config(1);
    let mut tracks = schedule.tracks().to_vec();
    tracks[0].events[1].value = GestureValue::AngleRad(core::f64::consts::FRAC_PI_2);
    let invalid_angle = GestureSchedule::try_new(700, tracks).unwrap();
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), invalid_angle, cfg).is_err());
    let mut tracks = schedule.tracks().to_vec();
    tracks[0].target = GestureTarget::JetSpeed;
    tracks[0].initial = GestureValue::VelocityMPerS(10.0);
    tracks[0].events.clear();
    let unsupported = GestureSchedule::try_new(700, tracks).unwrap();
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), unsupported, cfg).is_err());
    let mut tracks = schedule.tracks().to_vec();
    tracks[0].target = GestureTarget::BlowingPressure;
    tracks[0].initial = GestureValue::PressurePa(1.0);
    tracks[0].events.clear();
    let duplicate = GestureSchedule::try_new(700, tracks).unwrap();
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), duplicate, cfg).is_err());
}

#[test]
fn mouth_and_terminal_pressure_superpose_from_signed_physical_flows_without_changing_motion() {
    use crate::pcm_wav::baffled::{SurfaceSample, RayleighMedium};
    let g = gas();
    let spec = source(g, -0.9);
    let cfg = config(1);
    let input = FlueMouthRadiation { center_m: [-0.3, 0.0], cells: [4, 4], aperture_size_m: None };
    let mut original = voice(spec.clone());
    original.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let size = [original.geometry.flue_width_m, original.geometry.cut_up_m];
    let area = size[0] * size[1];
    let surface: Vec<_> = (0..4).flat_map(|y| (0..4).map(move |x| SurfaceSample {
        position_m: [input.center_m[0] + (x as f64 + 0.5) * size[0]/4.0 - size[0]/2.0,
            input.center_m[1] + (y as f64 + 0.5) * size[1]/4.0 - size[1]/2.0, 0.0],
        area_m2: area/16.0, mode_shape: vec![1.0/area],
    })).collect();
    let medium = RayleighMedium::from(&g);
    let mut mouth = BaffledPressure::new(&surface, 1, 48_000, cfg.receiver.position_m, medium).unwrap();
    let mut terminal = BaffledPressure::circular_outlet(0.006, 48_000, cfg.receiver, medium).unwrap();
    let mut expected = vec![0.0; cfg.samples as usize];
    let mut terminal_only = expected.clone();
    let mut mouth_only = expected.clone();
    for ((sum, duct), source_pressure) in expected.iter_mut().zip(&mut terminal_only).zip(&mut mouth_only) {
        original.step_block(&mut [0.0]).unwrap();
        let n = original.network().unwrap();
        let outward = n.node_frame(0).unwrap().net_flow_into_node_m3_s;
        assert!((outward + original.mouth_flow_m3_s()).abs() < 1e-16,
            "source-node sign must be outward, opposite to the solved q into the tube");
        *duct = terminal.step(&[n.node_frame(1).unwrap().net_flow_into_node_m3_s]).unwrap();
        *source_pressure = mouth.step(&[outward]).unwrap();
        *sum = *duct + *source_pressure;
    }
    assert!(mouth_only.iter().any(|p| p.abs() > 1e-10));
    assert!(terminal_only.iter().any(|p| p.abs() > 1e-10));
    assert!(terminal_only.iter().zip(&mouth_only).any(|(a,b)| a*b < 0.0),
        "signed aperture contributions must permit destructive interference");
    for partition in [513, 37, 1] {
        let mut performance = FlueOutletPerformance::new(voice(spec.clone()), schedule(300.0), cfg)
            .unwrap().with_mouth_radiation(input).unwrap();
        let mut actual = vec![0.0; cfg.samples as usize];
        for block in actual.chunks_mut(partition) { performance.block(block).unwrap(); }
        assert_eq!(actual, expected);
        assert_eq!(performance.voice().mouth_flow_m3_s(), original.mouth_flow_m3_s());
        assert_eq!(performance.voice().network().unwrap().stored_energy_j(),
            original.network().unwrap().stored_energy_j());
        assert_eq!(performance.mouth_receiver().unwrap().filter_delay_samples(),
            performance.receiver().filter_delay_samples());
        assert!(performance.with_mouth_radiation(input).is_err());
    }
}

#[test]
fn mouth_observation_shares_the_exterior_gas_and_rejects_duplicate_geometry() {
    let cfg = config(1);
    let input = FlueMouthRadiation { center_m: [-0.3, 0.0], cells: [4, 4], aperture_size_m: None };
    let build = |exterior| FlueOutletPerformance::with_radiation(voice(source(gas(), -0.9)),
        schedule(300.0), cfg, FlueRadiationInput { exterior, additional_outlets: vec![] })
        .unwrap().with_mouth_radiation(input).unwrap();
    let mut cold = build(gas());
    let mut hot = build(GasState::try_new(&GasSpec::dry_air_ussa1976(), 330.0, gas().pressure).unwrap());
    let mut a = [0.0; 513];
    let mut b = [0.0; 513];
    cold.block(&mut a).unwrap();
    hot.block(&mut b).unwrap();
    assert_ne!(a, b, "one exterior state must affect source and terminal radiation together");
    assert_eq!(cold.voice().mouth_flow_m3_s(), hot.voice().mouth_flow_m3_s(),
        "one-way observation must not alter source motion");
    for refused in [FlueMouthRadiation { center_m: [0.0, 0.0], ..input },
        FlueMouthRadiation { cells: [0, 4], ..input },
        FlueMouthRadiation { center_m: [f64::NAN, 0.0], ..input }] {
        let p = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(300.0), cfg).unwrap();
        assert!(p.with_mouth_radiation(refused).is_err());
    }
    assert!(build(gas()).with_mouth_radiation(input).is_err());
}

#[test]
fn independent_aperture_changes_observation_and_reaction_without_redefining_the_jet() {
    let gate = CancelGate::new_clock_free();
    let input = FlueMouthRadiation { center_m: [-0.3, 0.0], cells: [8, 8],
        aperture_size_m: Some([0.010, 0.012]) };
    let make = |mouth| FlueOutletPerformance::new(voice(source(gas(), -0.9)),
        schedule(300.0), config(1)).unwrap().with_mouth_radiation(mouth).unwrap();
    let mut explicit = make(input);
    let mut legacy = make(FlueMouthRadiation { aperture_size_m: None, ..input });
    let jet = explicit.voice.geometry;
    assert_ne!(input.size_m(&jet), [jet.flue_width_m, jet.cut_up_m]);
    let mut a = [0.0; 513];
    let mut b = [0.0; 513];
    explicit.block(&mut a).unwrap();
    legacy.block(&mut b).unwrap();
    assert_ne!(a, b, "the independent aperture must reach microphone pressure");
    assert_eq!(explicit.voice.mouth_flow_m3_s(), legacy.voice.mouth_flow_m3_s());
    assert_eq!(explicit.voice.network().unwrap().stored_energy_j(),
        legacy.voice.network().unwrap().stored_energy_j());
    assert_eq!(explicit.voice.geometry, jet);

    let mut cold = voice(source(gas(), -0.9));
    cold.geometry.mouth_loss_rel = 0.0;
    cold.resistance = 0.0;
    let reacted = FlueOutletPerformance::new(cold, schedule(300.0), config(1)).unwrap()
        .with_mouth_radiation(input).unwrap().with_mouth_reaction(&gate).unwrap();
    assert_eq!(reacted.voice.mouth_radiation_load().unwrap().size_m(), [0.010, 0.012]);
    assert_eq!(reacted.voice.geometry.flue_width_m, jet.flue_width_m);
    assert_eq!(reacted.voice.geometry.cut_up_m, jet.cut_up_m);
    for size in [[0.0, 0.012], [-0.01, 0.012], [f64::NAN, 0.012]] {
        assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)),
            schedule(300.0), config(1)).unwrap()
            .with_mouth_radiation(FlueMouthRadiation { aperture_size_m: Some(size), ..input }).is_err());
    }
}
