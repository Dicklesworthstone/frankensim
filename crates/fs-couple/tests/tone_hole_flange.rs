//! Per-aperture unflanged open-end reaction (MB 92): the card word `outlet ... unflanged`, the
//! per-opening binding it reaches, and the guarantee that selecting nothing builds the same
//! network as before. Reaction unflanged, observation baffled half-space (CONTRACT.md).
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::radiation::environment::{
    bind_baffled_environment, bind_baffled_environment_with_flanges,
};
use fs_couple::bernoulli_aperture::radiation::{BaffledRadiationLoad, RadiationFlange};
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctRadiationSpec, DuctReconstruction,
};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const RATE: u32 = 48_000;
const HOLE: f64 = 0.003;
fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101_325.0).unwrap()
}
/// A bore with two open side-hole chimneys (nodes 2 and 3), both ending in a memoryless open end.
fn graph(g: GasState) -> TubeNetworkSpec {
    let section = |nodes, length_m, radius_m| TubeSection { nodes, length_m, radius_m, max_length_error_m: 0.005 };
    TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -1.0 }, NetworkNode::Termination { reflection: -1.0 }],
        sections: vec![section([0, 1], 0.12, 0.009), section([1, 2], 0.012, HOLE), section([1, 3], 0.012, HOLE)],
        sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20 }
}
fn radiation(exterior: GasState) -> DuctRadiationSpec {
    DuctRadiationSpec {
        openings: vec![DuctOpening { network: 0, node: 2, center_m: [0.0, 0.0] },
            DuctOpening { network: 0, node: 3, center_m: [0.05, 0.0] }],
        exterior: Some(exterior),
        receiver: CircularOutletReceiver { position_m: [0.0, 0.0, 0.5], radial_rings: 16,
            angular_points: 64, maximum_frequency_hz: 400.0 },
        reconstruction: DuctReconstruction::Differentiating {
            mechanics_rate_hz: RATE, output_rate_hz: RATE / 2 },
    }
}
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn load(g: GasState, flange: RadiationFlange) -> BaffledRadiationLoad {
    BaffledRadiationLoad::with_flange(HOLE, g.density, g.sound_speed, 1.0 / f64::from(RATE), 400.0,
        flange, &gate()).unwrap()
}

#[test]
fn an_empty_or_all_baffled_list_is_the_original_binding_bit_for_bit() {
    let g = gas(300.15);
    let gases = [g; 3];
    let run = |flanges: Option<&[RadiationFlange]>| {
        let mut graphs = [graph(g)];
        let environment = match flanges {
            None => bind_baffled_environment(&mut graphs, &[&gases], radiation(g), &gate()),
            Some(list) => bind_baffled_environment_with_flanges(&mut graphs, &[&gases], radiation(g), list, &gate()),
        }.unwrap();
        (format!("{graphs:?}"), environment.bindings().to_vec(), environment.radiation().clone())
    };
    let original = run(None);
    assert_eq!(run(Some(&[])), original);
    assert_eq!(run(Some(&[RadiationFlange::Baffled, RadiationFlange::Baffled])), original);
    // The load is exactly the one the original constructor builds.
    for binding in &original.1 {
        assert_eq!(*binding.load(), load(g, RadiationFlange::Baffled));
        assert_eq!(binding.load().flange(), RadiationFlange::Baffled);
    }
}

#[test]
fn an_unflanged_opening_gets_the_unflanged_load_and_only_that_opening() {
    let g = gas(300.15);
    let gases = [g; 3];
    let mut graphs = [graph(g)];
    let before = graph(g);
    let environment = bind_baffled_environment_with_flanges(&mut graphs, &[&gases], radiation(g),
        &[RadiationFlange::Baffled, RadiationFlange::Unflanged], &gate()).unwrap();
    assert_eq!(graphs[0].nodes[2], load(g, RadiationFlange::Baffled).termination());
    assert_eq!(graphs[0].nodes[3], load(g, RadiationFlange::Unflanged).termination());
    assert_eq!(graphs[0].nodes[..2], before.nodes[..2]);
    assert_eq!(format!("{:?}", graphs[0].sections), format!("{:?}", before.sections));
    let bindings = environment.bindings();
    assert_eq!(bindings[0].load().flange(), RadiationFlange::Baffled);
    assert_eq!(bindings[1].load().flange(), RadiationFlange::Unflanged);
    for binding in bindings {
        binding.validate_receiver(&graphs[0], environment.radiation(), g).unwrap();
    }
}

#[test]
fn planted_negative_unflanged_must_differ_from_baffled_by_the_source_law() {
    let g = gas(300.15);
    let gases = [g; 3];
    let bind = |list: &[RadiationFlange]| {
        let mut graphs = [graph(g)];
        let environment = bind_baffled_environment_with_flanges(&mut graphs, &[&gases], radiation(g),
            list, &gate()).unwrap();
        (graphs, environment.bindings().to_vec())
    };
    let (baffled, baffled_bindings) = bind(&[]);
    let (unflanged, unflanged_bindings) = bind(&[RadiationFlange::Unflanged, RadiationFlange::Unflanged]);
    assert_ne!(baffled[0].nodes, unflanged[0].nodes);
    for (b, u) in baffled_bindings.iter().zip(&unflanged_bindings) {
        // At low ka: resistance (ka)^2/4 against (ka)^2/2, end correction 0.6133 a against 8/(3 pi) a.
        let (rb, xb) = b.load().impedance_at(100.0, false).unwrap();
        let (ru, xu) = u.load().impedance_at(100.0, false).unwrap();
        assert!((ru / rb - 0.5).abs() < 0.01, "resistance ratio {}", ru / rb);
        let expected = 0.6133 / (8.0 / (3.0 * core::f64::consts::PI));
        assert!((xu / xb - expected).abs() < 0.01, "reactance ratio {}", xu / xb);
    }
}

#[test]
fn a_flange_list_of_the_wrong_length_refuses_and_changes_nothing() {
    let g = gas(300.15);
    let gases = [g; 3];
    let mut graphs = [graph(g)];
    let original = format!("{graphs:?}");
    for list in [&[RadiationFlange::Unflanged][..], &[RadiationFlange::Unflanged; 3][..]] {
        assert!(bind_baffled_environment_with_flanges(&mut graphs, &[&gases], radiation(g), list, &gate()).is_err());
        assert_eq!(format!("{graphs:?}"), original);
    }
}

const RECORDER_A1: &str =
    include_str!("../../../data/instruments/descriptions/recorder/recorder-a1.performance");
const HOLE_ROWS: [&str; 2] = ["outlet 19 -0.1199000 0\n", "outlet 21 -0.0863000 0\n"];

/// The alto recorder's a1 card cut to a 1024-sample window with its onset only (as flue_window_mass.rs).
fn recorder_a1() -> String {
    assert!(RECORDER_A1.contains("\noutput 120000 5 "));
    let events = "events\t2\nevent\t0e0\t5e-2\tpressure-pa\t2.1e2\nevent\t2.2e0\t1.05e-2\tpressure-pa\t0e0";
    assert!(RECORDER_A1.contains(events), "the a1 gesture changed; update this test");
    for row in HOLE_ROWS { assert!(RECORDER_A1.contains(row), "the a1 outlet rows changed; update this test"); }
    RECORDER_A1
        .replace("\noutput 120000 5 ", "\noutput 1024 1 ")
        .replace(events, "events\t1\nevent\t0e0\t5e-2\tpressure-pa\t2.1e2")
}
fn unflanged_holes(text: &str) -> String {
    HOLE_ROWS.iter().fold(text.to_string(), |t, row|
        t.replace(row, &format!("{} unflanged\n", row.trim_end())))
}
fn parse(text: &str) -> Result<FlueFilePerformance, impl core::fmt::Debug> {
    FlueFilePerformance::from_bytes(text.as_bytes(), 256, &gate())
}

#[test]
fn the_card_word_binds_the_two_tone_holes_and_keeps_every_chimney_term() {
    let base = recorder_a1();
    let plain = parse(&base).unwrap();
    let option = parse(&unflanged_holes(&base)).unwrap();
    let flanges = |file: &FlueFilePerformance| file.radiation_bindings().iter()
        .map(|b| (b.opening().node, b.load().flange())).collect::<Vec<_>>();
    assert_eq!(flanges(&plain), [(24, RadiationFlange::Baffled), (19, RadiationFlange::Baffled),
        (21, RadiationFlange::Baffled)]);
    assert_eq!(flanges(&option), [(24, RadiationFlange::Baffled), (19, RadiationFlange::Unflanged),
        (21, RadiationFlange::Unflanged)]);
    // Absent word: every installed load equals the original constructor's at the binding's own radius, gas and band.
    for binding in plain.radiation_bindings() {
        let l = binding.load();
        assert_eq!(*l, BaffledRadiationLoad::new(l.radius_m(), l.density_kg_m3(), l.sound_speed_m_s(),
            l.time_step_s(), l.maximum_frequency_hz(), &gate()).unwrap());
    }
    // Only nodes 19 and 21 change; chimney lengths with t_i, t_m, t_a (the `inner ... matching series` rows) stay.
    let nodes = |file: &FlueFilePerformance| file.renderer().voice().network().unwrap().source().graph.nodes.clone();
    let (a, b) = (nodes(&plain), nodes(&option));
    assert_eq!(a.len(), b.len());
    for (index, (x, y)) in a.iter().zip(&b).enumerate() {
        assert_eq!(x == y, index != 19 && index != 21, "node {index}");
    }
    let sections = |file: &FlueFilePerformance|
        format!("{:?}", file.renderer().voice().network().unwrap().source().graph.sections);
    assert_eq!(sections(&plain), sections(&option));
    for node in [19, 21] {
        let l = option.radiation_bindings().iter().find(|b| b.opening().node == node).unwrap().load();
        assert_eq!(b[node], l.termination());
    }
}

#[test]
fn the_card_word_is_refused_outside_its_scope() {
    let base = recorder_a1();
    let with_word = unflanged_holes(&base);
    assert!(parse(&with_word).is_ok());
    // Needs the compact per-outlet reaction.
    assert!(parse(&with_word.replace("radiation compact-baffled", "radiation prescribed")).is_err());
    // The primary row takes no word; unknown or doubled words refuse.
    assert!(parse(&base.replace("outlets 24 2\n", "outlets 24 2 unflanged\n")).is_err());
    assert!(parse(&base.replace(HOLE_ROWS[0], "outlet 19 -0.1199000 0 flanged\n")).is_err());
    assert!(parse(&base.replace(HOLE_ROWS[0], "outlet 19 -0.1199000 0 unflanged unflanged\n")).is_err());
}

#[test]
fn planted_negative_the_option_changes_the_rendered_pressure() {
    let base = recorder_a1();
    let render = |text: &str| {
        let file = parse(text).unwrap();
        let mut renderer = file.into_renderer();
        let mut pressure = vec![0.0; 1024];
        for block in pressure.chunks_mut(256) { renderer.block(block).unwrap(); }
        assert!(pressure.iter().all(|p| p.is_finite()));
        assert!(pressure.iter().any(|p| p.abs() > 1e-10));
        pressure
    };
    let plain = render(&base);
    assert_eq!(plain, render(&base), "the control is deterministic");
    assert_ne!(plain, render(&unflanged_holes(&base)));
}
