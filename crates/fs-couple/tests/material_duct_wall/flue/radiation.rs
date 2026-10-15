//! Independent opening superposition and the complete material/jet/output path.
use super::*;
use fs_couple::pcm_wav::baffled::openings::CircularOpening;
use fs_couple::flue_loop::network::performance::{FlueRadiatingOutlet, FlueRadiationInput};

fn medium(g: GasState) -> RayleighMedium {
    RayleighMedium { density: g.density, sound_speed: g.sound_speed }
}
fn shifted(mut receiver: CircularOutletReceiver, center: [f64; 2]) -> CircularOutletReceiver {
    receiver.position_m[0] -= center[0];
    receiver.position_m[1] -= center[1];
    receiver
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12 * (1.0 + a.abs() + b.abs()), "{a} != {b}");
}

#[test]
fn fused_openings_equal_original_single_disk_arithmetic_and_signed_superposition() {
    let cfg = config().receiver;
    let openings = [CircularOpening { center_m: [-0.1, 0.0], radius_m: 0.006 },
        CircularOpening { center_m: [0.1, 0.0], radius_m: 0.0035 }];
    let mut a = BaffledPressure::circular_outlet(openings[0].radius_m, RATE,
        shifted(cfg, openings[0].center_m), medium(gas())).unwrap();
    let mut one = BaffledPressure::circular_openings(&openings[..1], RATE, cfg, medium(gas())).unwrap();
    let mut b = BaffledPressure::circular_outlet(openings[1].radius_m, RATE,
        shifted(cfg, openings[1].center_m), medium(gas())).unwrap();
    let mut combined = BaffledPressure::circular_openings(&openings, RATE, cfg, medium(gas())).unwrap();
    let mut heard = false;
    for n in 0..1024 {
        let q = [1e-5 * (n as f64 * 0.043).sin(), -2e-5 * (n as f64 * 0.057).cos()];
        let pa = a.step(&q[..1]).unwrap();
        assert_eq!(pa, one.step(&q[..1]).unwrap());
        let pb = b.step(&q[1..]).unwrap();
        let p = combined.step(&q).unwrap();
        near(p, pa + pb);
        let mut parts = [0.0; 2];
        combined.mode_pressures(&mut parts).unwrap();
        near(parts[0], pa); near(parts[1], pb);
        heard |= p.abs() > 1e-6;
    }
    assert!(heard);
}

#[test]
fn opposite_flows_cancel_on_axis_but_not_at_an_off_axis_microphone() {
    let openings = [CircularOpening { center_m: [-0.1, 0.0], radius_m: 0.006 },
        CircularOpening { center_m: [0.1, 0.0], radius_m: 0.006 }];
    let cfg = config().receiver;
    let mut center = BaffledPressure::circular_openings(&openings, RATE, cfg, medium(gas())).unwrap();
    let mut side_cfg = cfg; side_cfg.position_m[0] = 0.4;
    let mut side = BaffledPressure::circular_openings(&openings, RATE, side_cfg, medium(gas())).unwrap();
    let mut peak = 0.0_f64;
    for n in 0..1024 {
        let q = 1e-5 * (n as f64 * 0.09).sin();
        assert!(center.step(&[q, -q]).unwrap().abs() < 1e-12);
        peak = peak.max(side.step(&[q, -q]).unwrap().abs());
    }
    assert!(peak > 1e-5); // Geometry and signed phase, not a positive amplitude mix.
}

#[test]
fn opening_admission_and_invalid_samples_preserve_the_shared_history() {
    let cfg = config().receiver;
    let opening = CircularOpening { center_m: [0.0, 0.0], radius_m: 0.006 };
    assert!(BaffledPressure::circular_openings(&[], RATE, cfg, medium(gas())).is_err());
    assert!(BaffledPressure::circular_openings(&[opening, opening], RATE, cfg, medium(gas())).is_err());
    let far = CircularOpening { center_m: [0.2, 0.0], ..opening };
    let mut excessive = cfg; excessive.radial_rings = 128; excessive.angular_points = 512;
    assert!(BaffledPressure::circular_openings(&[opening, far], RATE, excessive, medium(gas())).is_err());
    let mut bad = far; bad.center_m[1] = f64::NAN;
    assert!(BaffledPressure::circular_openings(&[opening, bad], RATE, cfg, medium(gas())).is_err());
    let mut a = BaffledPressure::circular_openings(&[opening, far], RATE, cfg, medium(gas())).unwrap();
    let mut b = BaffledPressure::circular_openings(&[opening, far], RATE, cfg, medium(gas())).unwrap();
    for n in 0..600 {
        assert!(a.step(&[0.0, f64::NAN]).is_err());
        assert!(a.step(&[0.0]).is_err());
        let q = [1e-5 * (n as f64 * 0.07).sin(), 2e-5 * (n as f64 * 0.05).sin()];
        assert_eq!(a.step(&q).unwrap(), b.step(&q).unwrap());
    }
}

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
    spec.gases[3] = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas().pressure).unwrap();
    spec
}
fn radiation() -> FlueRadiationInput {
    FlueRadiationInput {
        exterior: GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, gas().pressure).unwrap(),
        additional_outlets: vec![FlueRadiatingOutlet { node: 4, center_m: [0.2, 0.0] }],
    }
}
fn multi_config() -> FlueOutletConfig { FlueOutletConfig { outlet_node: 3, ..config() } }

#[test]
fn complete_material_jet_and_regional_bore_feed_both_physical_outlets_on_one_clock() {
    let r = radiation();
    let cfg = multi_config();
    let mut reference = voice(branched(), input(1e6));
    reference.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let mut a = BaffledPressure::circular_outlet(0.006, RATE, cfg.receiver, medium(r.exterior)).unwrap();
    let mut b = BaffledPressure::circular_outlet(0.0035, RATE,
        shifted(cfg.receiver, r.additional_outlets[0].center_m), medium(r.exterior)).unwrap();
    let mut expected = Vec::new();
    for _ in 0..1024 {
        reference.step_block(&mut [0.0]).unwrap();
        let n = reference.network().unwrap();
        let pa = a.step(&[n.node_frame(3).unwrap().net_flow_into_node_m3_s]).unwrap();
        let pb = b.step(&[n.node_frame(4).unwrap().net_flow_into_node_m3_s]).unwrap();
        expected.push(pa + pb);
    }
    let mut exact = None;
    for partition in [1024, 1, 37] {
        let mut performance = FlueOutletPerformance::with_radiation(voice(branched(), input(1e6)),
            schedule(), cfg, r.clone()).unwrap();
        assert_eq!(performance.outlet_nodes(), &[3, 4]);
        assert_eq!(performance.radiation_input(), Some(&r));
        let mut output = vec![0.0; 1024];
        for block in output.chunks_mut(partition) { performance.block(block).unwrap(); }
        for (&actual, &expected) in output.iter().zip(&expected) { near(actual, expected); }
        if let Some(previous) = &exact { assert_eq!(&output, previous); } else { exact = Some(output); }
        assert_eq!(performance.samples_rendered(), 1024);
        let n = performance.voice().network().unwrap();
        assert_eq!(n.wall_state(1).unwrap(), reference.network().unwrap().wall_state(1).unwrap());
        assert_eq!(n.source().gases, branched().gases);
    }
    assert!(expected.iter().any(|p| p.abs() > 1e-8));
}

#[test]
fn exterior_only_changes_observation_and_single_port_legacy_remains_identical() {
    let cfg = multi_config();
    let mut outputs = Vec::new();
    let mut states = Vec::new();
    for t in [293.15, 350.0] {
        let r = FlueRadiationInput { exterior: GasState::try_new(&GasSpec::dry_air_ussa1976(),
            t, gas().pressure).unwrap(), additional_outlets: vec![] };
        let mut p = FlueOutletPerformance::with_radiation(voice(branched(), input(1e6)), schedule(), cfg, r).unwrap();
        let mut output = vec![0.0; 1024]; p.block(&mut output).unwrap();
        states.push(p.voice().network().unwrap().wall_state(1).unwrap());
        outputs.push(output);
    }
    assert_eq!(states[0], states[1]); // An observer must not reinitialize the source.
    assert_ne!(outputs[0], outputs[1]); // Actual exterior travel/density changed.
    let mut legacy = FlueOutletPerformance::new(voice(branched(), input(1e6)), schedule(), cfg).unwrap();
    let mut output = vec![0.0; 1024]; legacy.block(&mut output).unwrap();
    assert_eq!(outputs[0], output);
}

#[test]
fn duplicate_ports_incoherent_gases_and_hidden_mean_flow_are_refused() {
    let check = |r| assert!(FlueOutletPerformance::with_radiation(voice(branched(), input(1e6)),
        schedule(), multi_config(), r).is_err());
    let mut r = radiation(); r.additional_outlets[0].node = 3; check(r);
    let mut r = radiation(); r.additional_outlets[0].node = 1; check(r);
    let mut r = radiation(); r.additional_outlets[0].center_m = [0.001, 0.0]; check(r);
    let mut r = radiation(); r.exterior.density *= 2.0; check(r);
    let mut r = radiation();
    r.exterior = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, 80_000.0).unwrap(); check(r);
}
