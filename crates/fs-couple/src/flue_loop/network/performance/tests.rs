use super::*;
use crate::flue_loop::network::tests::{gas, source, voice};
use crate::bernoulli_aperture::network::TubeSection;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::gesture::{GestureEvent, GestureTarget, GestureTrack, GestureValue};

fn schedule(pressure: f64) -> GestureSchedule {
    GestureSchedule::try_new(700, vec![GestureTrack { id: "blow".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(pressure),
        events: vec![] }]).unwrap()
}
fn config(node: usize) -> FlueOutletConfig {
    FlueOutletConfig {
        sample_rate_hz: 48_000, samples: 513, max_block: 513, max_controls: 64,
        max_compile_work: 10_000, outlet_node: node,
        receiver: CircularOutletReceiver { position_m: [0.0, 0.0, 1.0], radial_rings: 4,
            angular_points: 16, maximum_frequency_hz: 1000.0 },
        maximum_abs_pressure_pa: 1.0e6,
    }
}

#[test]
fn exterior_output_matches_actual_terminal_flow_through_the_shared_receiver() {
    let g = gas();
    let spec = source(g, -0.9);
    let mut original = voice(spec.clone());
    original.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let cfg = config(1);
    let mut receiver = BaffledPressure::circular_outlet(0.006, cfg.sample_rate_hz, cfg.receiver,
        RayleighMedium { density: g.density, sound_speed: g.sound_speed }).unwrap();
    let tube_delay = original.network().unwrap().represented_sections()[0].one_way_samples;
    let first = tube_delay + receiver.delay_samples.0;
    let mut launched = vec![0.0; 513];
    let mut expected = Vec::new();
    original.step_block_observed(&mut launched, |line| {
        let FlueResonator::Network(n) = line else { panic!("network") };
        expected.push(receiver.step(&[n.node_frame(1).unwrap().net_flow_into_node_m3_s]).unwrap());
        Ok(())
    }).unwrap();
    let mut performance = FlueOutletPerformance::new(voice(spec), schedule(300.0), cfg).unwrap();
    let mut actual = vec![0.0; 513];
    performance.block(&mut actual).unwrap();
    assert_eq!(actual, expected);
    assert!(actual[..first].iter().all(|p| *p == 0.0));
    assert!(actual[first..].iter().any(|p| p.abs() > 1.0e-12));
    assert_ne!(actual, launched); // internal waves must not masquerade as microphone pressure
    assert_eq!(performance.voice().diagnostics().len(), 1);
}

#[test]
fn nondivisor_pressure_gestures_and_release_preserve_physics_across_callbacks() {
    let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
        id: "blow".into(), target: GestureTarget::BlowingPressure,
        initial: GestureValue::PressurePa(300.0), events: vec![
            GestureEvent { time_s: 1.0/700.0, transition_s: 4.0/700.0, value: GestureValue::PressurePa(450.0) },
            GestureEvent { time_s: 3.0/700.0, transition_s: 2.0/700.0, value: GestureValue::PressurePa(200.0) },
            GestureEvent { time_s: 6.0/700.0, transition_s: 0.0, value: GestureValue::PressurePa(0.0) },
        ],
    }]).unwrap();
    let mut expected = Vec::new();
    for partition in [513, 1, 37, 256] {
        let mut performance = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule.clone(), config(1)).unwrap();
        let mut out = vec![0.0; 513];
        for block in out.chunks_mut(partition) { performance.block(block).unwrap(); }
        let bits: Vec<_> = out.iter().map(|p| p.to_bits()).collect();
        if expected.is_empty() { expected = bits; } else { assert_eq!(bits, expected); }
        assert_eq!(performance.samples_rendered(), 513);
        assert_eq!(performance.remaining_samples(), 0);
        assert_eq!(performance.voice().jet_speed_m_s(), 0.0);
        assert!(performance.voice().network().unwrap().stored_energy_j() > 0.0);
        assert!(out[480..].iter().any(|p| p.abs() > 1.0e-12));
    }
}

#[test]
fn microphone_uses_the_outlets_local_gas_not_the_jets_inlet_gas() {
    let inlet = gas();
    let outlet = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, inlet.pressure).unwrap();
    let mut spec = source(inlet, -0.9);
    spec.graph.nodes = vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.9 }];
    spec.graph.sections = vec![
        TubeSection { nodes: [0,1], length_m: 0.15, radius_m: 0.006, max_length_error_m: 0.005 },
        TubeSection { nodes: [1,2], length_m: 0.15, radius_m: 0.006, max_length_error_m: 0.005 },
    ];
    spec.gases = vec![inlet, outlet];
    let cfg = config(2);
    let mut direct = voice(spec.clone());
    direct.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let mut correct = BaffledPressure::circular_outlet(0.006, 48_000, cfg.receiver,
        RayleighMedium { density: outlet.density, sound_speed: outlet.sound_speed }).unwrap();
    let mut wrong = BaffledPressure::circular_outlet(0.006, 48_000, cfg.receiver,
        RayleighMedium { density: inlet.density, sound_speed: inlet.sound_speed }).unwrap();
    let mut expected = Vec::new();
    let mut incorrect = Vec::new();
    direct.step_block_observed(&mut [0.0; 513], |line| {
        let FlueResonator::Network(n) = line else { panic!("network") };
        let q = n.node_frame(2).unwrap().net_flow_into_node_m3_s;
        expected.push(correct.step(&[q]).unwrap());
        incorrect.push(wrong.step(&[q]).unwrap());
        Ok(())
    }).unwrap();
    let mut performance = FlueOutletPerformance::new(voice(spec), schedule(300.0), cfg).unwrap();
    let mut actual = vec![0.0; 513];
    performance.block(&mut actual).unwrap();
    assert_eq!(actual, expected);
    assert!(actual.iter().zip(&incorrect).any(|(a,b)| (a-b).abs() > 1.0e-9));
}

#[test]
fn zero_pressure_silence_and_failures_do_not_become_fake_audio() {
    let mut quiet = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(0.0), config(1)).unwrap();
    let mut out = [7.0; 513];
    quiet.block(&mut out).unwrap();
    assert!(out.iter().all(|p| *p == 0.0));
    let mut cfg = config(1);
    cfg.max_block = 64;
    let mut model = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(300.0), cfg).unwrap();
    let mut oversized = [7.0; 65];
    assert!(model.block(&mut oversized).is_err());
    assert_eq!(oversized, [7.0; 65]);
    assert_eq!(model.samples_rendered(), 0);
    assert!(model.validate_sample_rate(44_100).is_err());
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(1.0e9), config(1)).is_err());
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(300.0), config(0)).is_err());
    cfg = config(1);
    cfg.maximum_abs_pressure_pa = f64::MIN_POSITIVE;
    let mut limited = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(300.0), cfg).unwrap();
    assert!(limited.block(&mut out).is_err());
    assert!(limited.samples_rendered() > 0);
    assert!(matches!(limited.block(&mut [0.0]), Err(RenderError::Poisoned)));
}

mod embouchure;
