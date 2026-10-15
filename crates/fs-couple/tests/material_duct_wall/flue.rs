//! End-to-end material wall -> original jet -> actual outlet pressure.
use super::*;
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_couple::bernoulli_aperture::wall::material::MaterialWallInput;
use fs_couple::bernoulli_aperture::viscothermal::ViscothermalSectionSpec;
use fs_couple::flue_loop::{FlueControl, FlueError, FlueGeometry, FlueVoice};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::flue_loop::network::performance::{FlueOutletConfig, FlueOutletPerformance};
use fs_couple::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

fn geometry() -> FlueGeometry {
    FlueGeometry { flue_width_m: 0.010, flue_height_m: 0.001, cut_up_m: 0.004,
        labium_offset_m: 0.0, jet_angle_rad: 0.0,
        bore_area_m2: core::f64::consts::PI * 0.006 * 0.006, mouth_loss_rel: 0.02 }
}
fn source() -> FlueNetworkSpec {
    FlueNetworkSpec { graph: graph(), gases: vec![gas(); 2], viscothermal: vec![],
        minimum_jet_speed_m_s: 1.0, maximum_jet_mach: 0.3, max_jet_history_samples: 4096 }
}
fn input(young: f64) -> MaterialWallInput {
    MaterialWallInput { assignments: vec![assignment(1000.0, young, 2.0)], band_hz: BAND, max_walls: 1 }
}
fn voice(source: FlueNetworkSpec, walls: MaterialWallInput) -> FlueVoice {
    let card = mint_tonal_interim_card().unwrap();
    FlueVoice::from_material_network(Some(&card), 90.0, 1e-6, geometry(), source, walls,
        RATE, 3, 0, &CancelGate::new_clock_free()).unwrap()
}
fn schedule() -> GestureSchedule {
    GestureSchedule::try_new(700, vec![GestureTrack { id: "blow".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(300.0),
        events: vec![] }]).unwrap()
}
fn config() -> FlueOutletConfig {
    FlueOutletConfig { sample_rate_hz: RATE, samples: 1024, max_block: 1024,
        max_controls: 64, max_compile_work: 10000, outlet_node: 2,
        receiver: CircularOutletReceiver { position_m: [0.0, 0.0, 1.0], radial_rings: 4,
            angular_points: 16, maximum_frequency_hz: 1000.0 },
        maximum_abs_pressure_pa: 1e6 }
}

#[test]
fn no_wall_selection_preserves_the_existing_jet_path_bit_for_bit() {
    let card = mint_tonal_interim_card().unwrap();
    let mut old = FlueVoice::from_network(Some(&card), 90.0, 1e-6, geometry(), source(),
        RATE, 3, 0, &CancelGate::new_clock_free()).unwrap();
    let mut empty = input(1e6); empty.assignments.clear(); empty.max_walls = 0;
    let mut new = voice(source(), empty.clone());
    assert!(old.network().unwrap().wall_input().is_none());
    assert_eq!(new.network().unwrap().wall_input(), Some(&empty));
    for v in [&mut old, &mut new] { v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap(); }
    let mut a = [0.0; 1024]; let mut b = a;
    old.step_block(&mut a).unwrap(); new.step_block(&mut b).unwrap();
    assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    assert_eq!(old.network().unwrap().last_frame(), new.network().unwrap().last_frame());
}

#[test]
fn independent_material_changes_reach_both_source_feedback_and_outlet_flow() {
    let mut a = voice(source(), input(1e6));
    let mut b = voice(source(), input(4e6));
    let delay = a.network().unwrap().represented_sections()[0].one_way_samples;
    for v in [&mut a, &mut b] { v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap(); }
    let mut waves = [Vec::new(), Vec::new()];
    let mut flows = [Vec::new(), Vec::new()];
    for (index, v) in [&mut a, &mut b].into_iter().enumerate() {
        for _ in 0..1024 {
            let mut one = [0.0]; v.step_block(&mut one).unwrap(); waves[index].push(one[0]);
            let n = v.network().unwrap();
            let state = n.wall_state(1).unwrap();
            close(state.stored_energy_j, n.node_frame(1).unwrap().stored_energy_j);
            let frame = n.last_frame().unwrap();
            let scale = (frame.stored_energy_j + frame.storage_change_j.abs()
                + frame.inlet_work_j.abs() + frame.interior_loss_j + frame.terminal_loss_j).max(1e-30);
            assert!(frame.balance_residual_j().abs() < 1e-10 * scale);
            flows[index].push(n.node_frame(2).unwrap().net_flow_into_node_m3_s);
        }
        assert_eq!(n_original(v).sections, graph().sections);
        assert!(v.network().unwrap().wall_state(0).is_err());
    }
    // No instantaneous material-dependent gain: the changed wall is downstream.
    assert_eq!(waves[0][..2 * delay], waves[1][..2 * delay]);
    assert_ne!(waves[0][2 * delay..], waves[1][2 * delay..]);
    assert_ne!(flows[0], flows[1]);
}
fn n_original(v: &FlueVoice) -> &TubeNetworkSpec { &v.network().unwrap().source().graph }

#[test]
fn solid_sources_and_regional_gases_survive_viscothermal_subdivision() {
    let mut spec = source();
    spec.gases[1] = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, gas().pressure).unwrap();
    // The existing admitted 300..500 Hz wide-tube lane; no altered numerical gate.
    spec.viscothermal.push(ViscothermalSectionSpec { section: 0, cells: 4,
        minimum_frequency_hz: 300.0, maximum_frequency_hz: 500.0 });
    let original = spec.clone();
    let walls = input(1e6);
    let mut v = voice(spec, walls.clone());
    let n = v.network().unwrap();
    assert_eq!(n.source().graph, original.graph);
    assert_eq!(n.source().gases, original.gases);
    assert_eq!(n.source().viscothermal, original.viscothermal);
    assert_eq!(n.wall_input(), Some(&walls));
    assert_eq!(n.material_walls()[0].source(), &walls.assignments[0]);
    assert_eq!(n.material_walls()[0].length_m(), 0.15);
    assert_eq!(n.graph().nodes[1], n.material_walls()[0].patch().shunt().unwrap());
    assert!(n.graph().nodes.len() > original.graph.nodes.len());
    assert_eq!(n.viscothermal_sections()[0].source.section, 0);
    for (section, g) in n.graph().sections.iter().zip(n.section_gases()) {
        if section.nodes.contains(&2) { assert_eq!(*g, original.gases[1]); }
    }
    v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    v.step_block(&mut [0.0; 256]).unwrap();
    assert!(v.network().unwrap().wall_state(1).unwrap().stored_energy_j > 0.0);
}

#[test]
fn actual_microphone_uses_material_feedback_and_keeps_it_across_callbacks() {
    let mut spec = source();
    spec.gases[1] = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, gas().pressure).unwrap();
    let mut signals = Vec::new();
    for young in [1e6, 4e6] {
        let walls = input(young);
        let mut direct = voice(spec.clone(), walls.clone());
        direct.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
        let outlet_gas = spec.gases[1];
        let mut receiver = BaffledPressure::circular_outlet(0.006, RATE, config().receiver,
            RayleighMedium { density: outlet_gas.density, sound_speed: outlet_gas.sound_speed }).unwrap();
        let mut expected = Vec::new();
        for _ in 0..1024 {
            direct.step_block(&mut [0.0]).unwrap();
            let q = direct.network().unwrap().node_frame(2).unwrap().net_flow_into_node_m3_s;
            expected.push(receiver.step(&[q]).unwrap());
        }
        for partition in [1, 37, 1024] {
            let mut performance = FlueOutletPerformance::new(voice(spec.clone(), walls.clone()),
                schedule(), config()).unwrap();
            let mut output = vec![0.0; 1024];
            for block in output.chunks_mut(partition) { performance.block(block).unwrap(); }
            assert_eq!(output, expected);
            let n = performance.voice().network().unwrap();
            assert_eq!(n.wall_input(), Some(&walls));
            assert_eq!(n.wall_state(1).unwrap(), direct.network().unwrap().wall_state(1).unwrap());
            assert_eq!(performance.samples_rendered(), 1024);
        }
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        signals.push(expected);
    }
    assert_ne!(signals[0], signals[1]);
}

#[test]
fn bad_wall_inputs_keep_the_original_cold_refusal_boundaries() {
    let gate = CancelGate::new_clock_free();
    let card = mint_tonal_interim_card().unwrap();
    let mut wrong = input(1e6); wrong.assignments[0].node = 0;
    assert!(matches!(FlueVoice::from_material_network(None, 90.0, 1e-6, geometry(),
        source(), wrong.clone(), RATE, 3, 0, &gate), Err(FlueError::NoCard)));
    assert!(FlueVoice::from_material_network(Some(&card), 90.0, 1e-6, geometry(),
        source(), wrong, RATE, 3, 0, &gate).is_err());
    let mut nyquist = input(1e6); nyquist.band_hz[1] = f64::from(RATE) / 2.0;
    assert!(FlueVoice::from_material_network(Some(&card), 90.0, 1e-6, geometry(),
        source(), nyquist, RATE, 3, 0, &gate).is_err());
    gate.request();
    assert!(FlueVoice::from_material_network(Some(&card), 90.0, 1e-6, geometry(),
        source(), input(1e6), RATE, 3, 0, &gate).is_err());
}
