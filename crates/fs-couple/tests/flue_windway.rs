use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_couple::bernoulli_aperture::contraction::ContractionSpec;
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::flue_loop::{FlueControl, FlueGeometry, FlueVoice};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::flue_loop::network::performance::{FlueOutletConfig, FlueOutletPerformance};
use fs_couple::pcm_wav::{baffled::CircularOutletReceiver, observation::PressureRenderer};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

const RATE: u32 = 48_000;
const HEIGHT: f64 = 0.001;
const LENGTH: f64 = 0.015;

fn gas(viscosity_scale: f64) -> GasState {
    // Synthetic transport substitution isolates the causal viscosity path:
    // density and sound speed remain identical, not a measured named gas.
    let mut spec = GasSpec::dry_air_ussa1976();
    spec.sutherland_beta *= viscosity_scale;
    GasState::try_new(&spec, 293.15, 101_325.0).unwrap()
}
fn voice(gas: GasState, maximum_mach: f64) -> FlueVoice {
    let geometry = FlueGeometry {
        flue_width_m: 0.01, flue_height_m: HEIGHT, cut_up_m: 0.004,
        labium_offset_m: 0.0, jet_angle_rad: 0.0,
        bore_area_m2: core::f64::consts::PI * 0.006 * 0.006, mouth_loss_rel: 0.02,
    };
    let source = FlueNetworkSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.9 }],
            sections: vec![TubeSection { nodes: [0, 1], length_m: 0.3, radius_m: 0.006,
                max_length_error_m: 0.005 }],
            sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
        },
        gases: vec![gas], viscothermal: vec![], minimum_jet_speed_m_s: 1.0,
        maximum_jet_mach: maximum_mach, max_jet_history_samples: 4096,
    };
    FlueVoice::from_network(Some(&mint_tonal_interim_card().unwrap()), 90.0, 1e-6,
        geometry, source, RATE, 3, 0, &CancelGate::new_clock_free()).unwrap()
}
fn schedule(pressure: f64) -> GestureSchedule {
    GestureSchedule::try_new(700, vec![GestureTrack { id: "blow".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(pressure),
        events: vec![] }]).unwrap()
}
fn config() -> FlueOutletConfig {
    FlueOutletConfig { sample_rate_hz: RATE, samples: 513, max_block: 513,
        max_controls: 64, max_compile_work: 10_000, outlet_node: 1,
        receiver: CircularOutletReceiver { position_m: [0.0, 0.0, 0.25],
            radial_rings: 4, angular_points: 16, maximum_frequency_hz: 2000.0 },
        maximum_abs_pressure_pa: 1e6 }
}
fn render(voice: FlueVoice, partition: usize) -> (FlueOutletPerformance, Vec<f64>) {
    let mut performance = FlueOutletPerformance::new(voice, schedule(300.0), config()).unwrap();
    let mut pressure = vec![0.0; 513];
    for chunk in pressure.chunks_mut(partition) { performance.block(chunk).unwrap(); }
    (performance, pressure)
}

#[test]
fn viscosity_only_data_swap_changes_jet_timing_and_physical_microphone_pressure() {
    let a = gas(1.0);
    let b = gas(2.0);
    assert_eq!(a.density, b.density);
    assert_eq!(a.sound_speed, b.sound_speed);
    let (ideal_a, original) = render(voice(a, 0.3), 513);
    let (_, other_ideal) = render(voice(b, 0.3), 37);
    assert_eq!(original, other_ideal, "the omitted windway keeps the original ideal jet");
    assert_eq!(ideal_a.voice().jet_speed_m_s(), (2.0 * 300.0 / a.density).sqrt());
    let (low, low_pressure) = render(voice(a, 0.3).with_viscous_windway(LENGTH).unwrap(), 513);
    let (high, high_pressure) = render(voice(b, 0.3).with_viscous_windway(LENGTH).unwrap(), 37);
    for (performance, gas) in [(&low, a), (&high, b)] {
        let v = performance.voice();
        assert_eq!(v.viscous_windway().unwrap().gas(), Some(v.gas()));
        assert_eq!(v.viscous_windway().unwrap().dynamic_viscosity_pa_s(), gas.dynamic_viscosity);
        let speed = v.jet_speed_m_s();
        let viscous_drop = 12.0 * gas.dynamic_viscosity * LENGTH * speed / HEIGHT.powi(2);
        let kinetic_head = 0.5 * gas.density * speed.powi(2);
        assert!((viscous_drop + kinetic_head - 300.0).abs() < 1e-11);
        assert_eq!(v.reynolds(), speed * HEIGHT / (gas.dynamic_viscosity / gas.density));
        assert!((v.transit_s() * speed * v.island().convection_ratio - 0.004).abs() < 1e-15);
        assert!(v.diagnostics().iter().all(|d| d.ledger_defect_j.abs() < 1e-12));
    }
    assert!(low.voice().jet_speed_m_s() < ideal_a.voice().jet_speed_m_s());
    assert!(high.voice().jet_speed_m_s() < low.voice().jet_speed_m_s());
    assert!(high.voice().transit_s() > low.voice().transit_s());
    assert!(low_pressure.iter().any(|p| p.abs() > 1e-12));
    assert_ne!(original, low_pressure);
    assert_ne!(low_pressure, high_pressure, "transport data must reach actual outlet pressure");
}

#[test]
fn pressure_controls_and_schedule_admission_use_the_actual_viscous_speed() {
    let g = gas(1.0);
    // The ideal 200 Pa speed exceeds Mach 0.05, but the finite windway
    // stays below that SAME physically declared speed limit.
    let ideal = voice(g, 0.05);
    assert!(ideal.validate_blowing_pressure(200.0).is_err());
    let mut viscous = voice(g, 0.05).with_viscous_windway(LENGTH).unwrap();
    viscous.validate_blowing_pressure(200.0).unwrap();
    viscous.apply(FlueControl::SetBlowingPressure(200.0)).unwrap();
    assert!(viscous.jet_speed_m_s() < 0.05 * g.sound_speed);
    assert!(FlueOutletPerformance::new(voice(g, 0.05).with_viscous_windway(LENGTH).unwrap(),
        schedule(200.0), config()).is_ok());
    // Conversely, an ideal 1 Pa jet fits the old history while viscosity
    // makes its real transit too long. Do not silently clamp that delay.
    assert!(voice(g, 0.3).validate_blowing_pressure(1.0).is_ok());
    let old_speed = viscous.jet_speed_m_s();
    assert!(viscous.apply(FlueControl::SetBlowingPressure(1.0)).is_err());
    assert_eq!(viscous.jet_speed_m_s(), old_speed);
    assert_eq!(viscous.samples_rendered(), 0);
    assert_eq!(viscous.network().unwrap().stored_energy_j(), 0.0);
    assert!(FlueOutletPerformance::new(voice(g, 0.3).with_viscous_windway(LENGTH).unwrap(),
        schedule(1.0), config()).is_err());
    viscous.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
    assert_eq!(viscous.jet_speed_m_s(), 0.0);
    assert!(viscous.transit_s().is_infinite());
}

#[test]
fn windway_and_mouth_losses_share_the_same_partition_independent_receiver() {
    let build = || voice(gas(1.0), 0.3).with_viscous_windway(LENGTH).unwrap()
        .with_mouth_contraction(ContractionSpec { area_ratio: 0.6,
            maximum_contracted_mach: 0.2, provenance: "authored finite-channel fixture".into() }).unwrap();
    let (one, a) = render(build(), 513);
    let (many, b) = render(build(), 37);
    assert_eq!(a, b);
    assert_eq!(one.voice().mouth_flow_m3_s(), many.voice().mouth_flow_m3_s());
    assert_eq!(one.voice().network().unwrap().last_frame(), many.voice().network().unwrap().last_frame());
    assert!(one.voice().diagnostics().iter().any(|d| d.separation_loss_j > 0.0));
}

#[test]
fn cold_attachment_validates_length_and_already_held_pressure() {
    for length in [0.0, -0.001, f64::NAN, f64::INFINITY] {
        assert!(voice(gas(1.0), 0.3).with_viscous_windway(length).is_err());
    }
    let mut held = voice(gas(1.0), 0.3);
    held.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let ideal_speed = held.jet_speed_m_s();
    let held = held.with_viscous_windway(LENGTH).unwrap();
    assert!(held.jet_speed_m_s() < ideal_speed);
    assert!(held.with_viscous_windway(LENGTH).is_err());
    let mut late = voice(gas(1.0), 0.3);
    late.step_block(&mut [0.0]).unwrap();
    assert!(late.with_viscous_windway(LENGTH).is_err());
    let mut too_slow = voice(gas(1.0), 0.3);
    too_slow.apply(FlueControl::SetBlowingPressure(1.0)).unwrap();
    assert!(too_slow.with_viscous_windway(LENGTH).is_err());
}
