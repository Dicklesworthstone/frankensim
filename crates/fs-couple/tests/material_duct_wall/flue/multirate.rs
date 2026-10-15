//! End-to-end mechanics-frame gestures and material/gas output, no audio retiming.
use super::*;
use fs_couple::flue_loop::network::performance::{FlueRadiatingOutlet, FlueRadiationInput};
use fs_couple::pcm_wav::baffled::openings::CircularOpening;
use fs_couple::render::{ControlDelta, RenderError};
use fs_couple::render::schedule::{PressureGestureBinding, compile_pressure_gestures};
use fs_scenario::gesture::GestureEvent;

fn branched() -> FlueNetworkSpec {
    let mut spec = source();
    spec.graph.nodes = vec![NetworkNode::Inlet, NetworkNode::Junction, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.9 }, NetworkNode::Termination { reflection: -0.7 }];
    spec.graph.sections = vec![
        TubeSection { nodes: [0, 1], length_m: 0.1, radius_m: 0.006, max_length_error_m: 0.005 },
        TubeSection { nodes: [1, 2], length_m: 0.15, radius_m: 0.006, max_length_error_m: 0.005 },
        TubeSection { nodes: [2, 3], length_m: 0.2, radius_m: 0.006, max_length_error_m: 0.005 },
        TubeSection { nodes: [2, 4], length_m: 0.08, radius_m: 0.0035, max_length_error_m: 0.005 },
    ];
    spec.gases = vec![gas(); 4];
    spec.gases[3] = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, gas().pressure).unwrap();
    spec
}
fn high_rate_voice(rate: u32, young: f64) -> FlueVoice {
    let mut walls = input(young); walls.band_hz = [20.0, 22_000.0];
    let card = mint_tonal_interim_card().unwrap();
    FlueVoice::from_material_network(Some(&card), 90.0, 1e-6, geometry(), branched(), walls,
        rate, 3, 0, &CancelGate::new_clock_free()).unwrap()
}
fn radiation() -> FlueRadiationInput {
    FlueRadiationInput {
        exterior: GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas().pressure).unwrap(),
        additional_outlets: vec![FlueRadiatingOutlet { node: 4, center_m: [0.1, 0.0] }],
    }
}
fn cfg() -> FlueOutletConfig {
    let mut config = config(); config.outlet_node = 3;
    config.receiver.radial_rings = 64; config.receiver.angular_points = 64;
    config.receiver.maximum_frequency_hz = 18_000.0;
    config
}
fn gestures() -> GestureSchedule {
    GestureSchedule::try_new(700, vec![GestureTrack {
        id: "blow".into(), target: GestureTarget::BlowingPressure,
        initial: GestureValue::PressurePa(300.0), events: vec![
            GestureEvent { time_s: 1.0/700.0, transition_s: 4.0/700.0, value: GestureValue::PressurePa(450.0) },
            GestureEvent { time_s: 3.0/700.0, transition_s: 2.0/700.0, value: GestureValue::PressurePa(200.0) },
            GestureEvent { time_s: 6.0/700.0, transition_s: 0.0, value: GestureValue::PressurePa(0.0) },
        ],
    }]).unwrap()
}
fn observer(ratio: u32) -> BaffledPressure {
    let r = radiation();
    BaffledPressure::circular_openings_multirate_differentiating(&[
        CircularOpening { center_m: [0.0, 0.0], radius_m: 0.006 },
        CircularOpening { center_m: [0.1, 0.0], radius_m: 0.0035 },
    ], RATE * ratio, RATE, cfg().receiver,
        RayleighMedium { density: r.exterior.density, sound_speed: r.exterior.sound_speed }).unwrap()
}

#[test]
fn complete_substeps_and_subframe_controls_match_direct_material_physics_and_pressure() {
    let ratio = 4;
    let mechanics = RATE * ratio;
    let config = cfg();
    let schedule = gestures();
    let controls = compile_pressure_gestures(&schedule, &[PressureGestureBinding {
        track: "blow".into(), voice: 0,
    }], mechanics, config.samples * u64::from(ratio), config.max_compile_work).unwrap();
    assert!(controls.iter().any(|c| c.sample % u64::from(ratio) != 0));
    let mut direct = high_rate_voice(mechanics, 1e6);
    let mut receiver = observer(ratio);
    let mut endpoint_only = observer(ratio);
    let mut expected = Vec::new();
    let mut endpoint_signal = Vec::new();
    let mut next = 0;
    let mut trace = vec![0.0; 2 * ratio as usize];
    for frame in 0..config.samples {
        for sub in 0..ratio as usize {
            let sample = frame * u64::from(ratio) + sub as u64;
            if let Some(control) = controls.get(next) && control.sample == sample {
                let ControlDelta::SetBlowingPressure { pressure_pa, .. } = control.delta else { panic!("pressure") };
                direct.apply(FlueControl::SetBlowingPressure(pressure_pa)).unwrap();
                next += 1;
            }
            direct.step_block(&mut [0.0]).unwrap();
            for (i, node) in [3, 4].into_iter().enumerate() {
                trace[2 * sub + i] = direct.network().unwrap().node_frame(node).unwrap().net_flow_into_node_m3_s;
            }
        }
        expected.push(receiver.step_trace(&trace).unwrap());
        let last = [trace[trace.len()-2], trace[trace.len()-1]];
        endpoint_signal.push(endpoint_only.step_trace(&last.repeat(ratio as usize)).unwrap());
    }
    assert!(expected.iter().zip(&endpoint_signal).any(|(a,b)| (a-b).abs() > 1e-12));
    for partition in [1, 37, 1024] {
        let mut performance = FlueOutletPerformance::with_multirate_radiation(
            high_rate_voice(mechanics, 1e6), schedule.clone(), config, radiation(), mechanics).unwrap();
        let mut actual = vec![0.0; config.samples as usize];
        for block in actual.chunks_mut(partition) { performance.block(block).unwrap(); }
        assert_eq!(actual, expected);
        assert_eq!(performance.samples_rendered(), config.samples);
        assert_eq!(performance.voice().samples_rendered(), config.samples * u64::from(ratio));
        assert_eq!(performance.mechanics_rate_hz(), mechanics);
        assert_eq!(performance.mechanics_per_output(), ratio as usize);
        assert_eq!(performance.remaining_samples(), 0);
        let n = performance.voice().network().unwrap();
        assert_eq!(n.wall_state(1).unwrap(), direct.network().unwrap().wall_state(1).unwrap());
        assert_eq!(n.source().gases, branched().gases);
        assert_eq!(n.wall_input(), direct.network().unwrap().wall_input());
        assert_eq!(performance.radiation_input(), Some(&radiation()));
        assert_eq!(performance.voice().jet_speed_m_s(), 0.0);
        assert!(n.stored_energy_j() > 0.0); // release does not discard retained waves
    }
}

#[test]
fn material_only_change_reaches_high_rate_output_without_bypassing_the_body() {
    let mut outputs = Vec::new();
    for young in [1e6, 4e6] {
        let mut performance = FlueOutletPerformance::with_multirate_radiation(
            high_rate_voice(4 * RATE, young), gestures(), cfg(), radiation(), 4 * RATE).unwrap();
        let mut out = vec![0.0; cfg().samples as usize];
        performance.block(&mut out).unwrap();
        let n = performance.voice().network().unwrap();
        assert_eq!(n.material_walls()[0].source().material, input(young).assignments[0].material);
        assert!(n.wall_state(1).unwrap().stored_energy_j > 0.0);
        outputs.push(out);
    }
    assert_ne!(outputs[0], outputs[1]);
}

#[test]
fn clock_window_callback_and_pressure_failures_are_explicit_and_do_not_resume_partial_frames() {
    let build = |config, mechanics| FlueOutletPerformance::with_multirate_radiation(
        high_rate_voice(4 * RATE, 1e6), gestures(), config, radiation(), mechanics);
    for rate in [0, RATE, RATE + 1, 17 * RATE] { assert!(build(cfg(), rate).is_err()); }
    let mut config = cfg(); config.max_block = 1_048_576;
    assert!(build(config, 4 * RATE).is_err());
    config = cfg(); config.samples = 1_u64 << 53;
    assert!(build(config, 4 * RATE).is_err());
    config = cfg(); config.max_block = 64;
    let mut p = build(config, 4 * RATE).unwrap();
    let mut too_large = [7.0; 65];
    assert!(p.block(&mut too_large).is_err());
    assert_eq!(too_large, [7.0; 65]); assert_eq!(p.voice().samples_rendered(), 0);
    assert!(p.validate_sample_rate(4 * RATE).is_err());
    p.validate_sample_rate(RATE).unwrap();
    p.block(&mut [0.0; 64]).unwrap();
    assert_eq!(p.samples_rendered(), 64); assert_eq!(p.voice().samples_rendered(), 256);
    config = cfg(); config.maximum_abs_pressure_pa = f64::MIN_POSITIVE;
    let mut p = build(config, 4 * RATE).unwrap();
    assert!(p.block(&mut vec![0.0; config.samples as usize]).is_err());
    assert!(matches!(p.block(&mut [0.0]), Err(RenderError::Poisoned)));
    // A control clock above PCM rate is admissible when the physical source
    // clock resolves it. The control budget is still counted BEFORE coalescing.
    let schedule = GestureSchedule::try_new(64_000, vec![GestureTrack {
        id: "blow".into(), target: GestureTarget::BlowingPressure,
        initial: GestureValue::PressurePa(300.0), events: vec![],
    }]).unwrap();
    config = cfg(); config.max_controls = 2000;
    assert!(FlueOutletPerformance::with_multirate_radiation(high_rate_voice(4 * RATE, 1e6),
        schedule.clone(), config, radiation(), 4 * RATE).is_ok());
    config.max_controls = 64;
    assert!(FlueOutletPerformance::with_multirate_radiation(high_rate_voice(4 * RATE, 1e6),
        schedule, config, radiation(), 4 * RATE).is_err());
}
