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
    voice_geometry(gas, maximum_mach, 0.004, 90.0)
}
fn voice_geometry(gas: GasState, maximum_mach: f64, cut_up_m: f64, gain: f64) -> FlueVoice {
    let geometry = FlueGeometry {
        flue_width_m: 0.01, flue_height_m: HEIGHT, cut_up_m,
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
    FlueVoice::from_network(Some(&mint_tonal_interim_card().unwrap()), gain, 1e-6,
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
    // makes its real transit too long. Do not silently clamp that delay: the
    // jet keeps its real transit, below the declared 1 m/s minimum speed, where
    // the voice reads its onset's rest state or takes the declared state (ruling 8132).
    assert!(voice(g, 0.3).validate_blowing_pressure(1.0).is_ok());
    viscous.apply(FlueControl::SetBlowingPressure(1.0)).unwrap();
    let speed = viscous.jet_speed_m_s();
    assert!(speed < viscous.minimum_jet_speed_m_s());
    assert!((viscous.transit_s() * speed * viscous.island().convection_ratio - 0.004).abs() < 1e-15);
    assert_eq!(viscous.samples_rendered(), 0);
    assert_eq!(viscous.network().unwrap().stored_energy_j(), 0.0);
    assert!(FlueOutletPerformance::new(voice(g, 0.3).with_viscous_windway(LENGTH).unwrap(),
        schedule(1.0), config()).is_ok());
    viscous.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
    assert_eq!(viscous.jet_speed_m_s(), 0.0);
    assert!(viscous.transit_s().is_infinite());
}

#[test]
fn centreline_pressure_admission_uses_the_card_mach_limit_without_changing_bulk_flow() {
    let g = gas(1.0);
    let length = 0.06;
    let mach = 0.025;
    let limit = mach * g.sound_speed;
    let build = || voice(g, mach).with_viscous_windway(length).unwrap();

    // At 100 Pa the bulk flow is admissible, while its physical Poiseuille centre
    // exceeds the declared Mach limit (and still lies below the inviscid bound).
    let mut bulk = build();
    bulk.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
    assert!(bulk.jet_speed_m_s() < limit);
    assert!(1.5 * bulk.jet_speed_m_s() > limit);
    assert!(1.5 * bulk.jet_speed_m_s() < (2.0 * 100.0 / g.density).sqrt());
    assert!(build().with_card_centreline_speed(100.0).is_err());
    assert!(bulk.with_card_centreline_speed(40.0).is_err(),
        "a valid receipt must not admit an invalid already-held pressure");

    let centred = || build().with_card_centreline_speed(40.0).unwrap();
    assert!(FlueOutletPerformance::new(centred(), schedule(40.0), config()).is_ok());
    assert!(FlueOutletPerformance::new(centred(), schedule(100.0), config()).is_err());
    let mut v = centred();
    v.apply(FlueControl::SetBlowingPressure(40.0)).unwrap();
    let mut bulk = build();
    bulk.apply(FlueControl::SetBlowingPressure(40.0)).unwrap();
    assert_eq!(v.jet_speed_m_s(), bulk.jet_speed_m_s());
    assert_eq!(v.card_speed_m_s(), 1.5 * v.jet_speed_m_s());

    // Invert the physical slit law at the declared centre-speed boundary. Nearby
    // pressures straddle that boundary without clipping either speed or the flux.
    let bulk_at_limit = limit / 1.5;
    let pressure_at_limit = 12.0 * g.dynamic_viscosity * length * bulk_at_limit
        / HEIGHT.powi(2) + 0.5 * g.density * bulk_at_limit.powi(2);
    let below = pressure_at_limit * (1.0 - 1e-6);
    let above = pressure_at_limit * (1.0 + 1e-6);
    v.apply(FlueControl::SetBlowingPressure(below)).unwrap();
    assert!(v.card_speed_m_s() < limit);
    let accepted = (v.jet_speed_m_s(), v.card_speed_m_s());
    assert!(v.validate_blowing_pressure(above).is_err());
    assert!(v.apply(FlueControl::SetBlowingPressure(above)).is_err());
    assert_eq!((v.jet_speed_m_s(), v.card_speed_m_s()), accepted);
    bulk.apply(FlueControl::SetBlowingPressure(above)).unwrap();
    assert!(bulk.jet_speed_m_s() < limit);
    assert!(1.5 * bulk.jet_speed_m_s() > limit);
}

#[test]
fn centreline_controls_preserve_the_inviscid_bound_and_refuse_atomically() {
    let g = gas(1.0);
    let build = || voice(g, 0.3).with_viscous_windway(0.06).unwrap();
    let mut bulk = build();
    bulk.apply(FlueControl::SetBlowingPressure(1000.0)).unwrap();
    assert!(bulk.jet_speed_m_s() < 0.3 * g.sound_speed);
    assert!(1.5 * bulk.jet_speed_m_s() < 0.3 * g.sound_speed);
    assert!(1.5 * bulk.jet_speed_m_s() > (2.0 * 1000.0 / g.density).sqrt());
    assert!(build().with_card_centreline_speed(1000.0).is_err());
    assert!(bulk.with_card_centreline_speed(100.0).is_err());

    let centred = || build().with_card_centreline_speed(100.0).unwrap();
    assert!(FlueOutletPerformance::new(centred(), schedule(1000.0), config()).is_err());
    let mut v = centred();
    let mut unchanged = centred();
    for source in [&mut v, &mut unchanged] {
        source.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        source.step_block(&mut [0.0; 32]).unwrap();
    }
    let before = (v.jet_speed_m_s(), v.card_speed_m_s(), v.samples_rendered(),
        v.mouth_flow_m3_s(), v.network().unwrap().stored_energy_j());
    assert!(v.validate_blowing_pressure(1000.0).is_err());
    assert!(v.apply(FlueControl::SetBlowingPressure(1000.0)).is_err());
    assert_eq!((v.jet_speed_m_s(), v.card_speed_m_s(), v.samples_rendered(),
        v.mouth_flow_m3_s(), v.network().unwrap().stored_energy_j()), before);
    let (mut actual, mut expected) = ([0.0; 64], [0.0; 64]);
    v.step_block(&mut actual).unwrap();
    unchanged.step_block(&mut expected).unwrap();
    assert_eq!(actual, expected, "refusal must also preserve the accepted pressure and source history");
    v.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
    assert_eq!(v.jet_speed_m_s(), 0.0);
    assert_eq!(v.card_speed_m_s(), 0.0);
}

#[test]
fn developing_card_preflight_uses_the_proposed_pressure_for_its_inviscid_core() {
    let g = gas(1.0);
    let build = || voice(g, 0.05).with_developing_windway(LENGTH, 100.0).unwrap();
    let mut bulk = build();
    bulk.apply(FlueControl::SetBlowingPressure(200.0)).unwrap();
    let core = (2.0 * 200.0 / g.density).sqrt();
    assert!(bulk.jet_speed_m_s() < 0.05 * g.sound_speed);
    assert!(1.5 * bulk.jet_speed_m_s() > core);
    assert!(core > 0.05 * g.sound_speed);

    // The accepted pressure is initially zero. Reading it during preflight
    // would report a zero core speed and incorrectly admit the 200 Pa input.
    let mut v = build().with_card_centreline_speed(100.0).unwrap();
    assert_eq!(v.card_speed_m_s(), 0.0);
    assert!(v.validate_blowing_pressure(200.0).is_err());
    assert!(v.apply(FlueControl::SetBlowingPressure(200.0)).is_err());
    assert_eq!(v.card_speed_m_s(), 0.0);
    assert!(FlueOutletPerformance::new(
        build().with_card_centreline_speed(100.0).unwrap(), schedule(200.0), config()).is_err());
    v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
    assert_eq!(v.card_speed_m_s(), (2.0 * 100.0 / g.density).sqrt());
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
    // A held pressure whose viscous jet is below the declared minimum speed attaches: that
    // speed is a validity edge, not a refusal (ruling 8132).
    let mut slow = voice(gas(1.0), 0.3);
    slow.apply(FlueControl::SetBlowingPressure(1.0)).unwrap();
    let slow = slow.with_viscous_windway(LENGTH).unwrap();
    assert!(slow.jet_speed_m_s() > 0.0 && slow.jet_speed_m_s() < slow.minimum_jet_speed_m_s());
}

#[test]
fn resolved_windway_changes_real_receiver_and_preserves_callback_independence() {
    let build = |g| voice(g, 0.3).with_unsteady_windway(LENGTH, 32, 96).unwrap()
        .with_mouth_contraction(ContractionSpec { area_ratio: 0.6,
            maximum_contracted_mach: 0.2, provenance: "authored unsteady-channel fixture".into() }).unwrap();
    let (one, a) = render(build(gas(1.0)), 513);
    let (many, b) = render(build(gas(1.0)), 37);
    assert_eq!(a, b);
    let first = one.voice().unsteady_windway().unwrap();
    let second = many.voice().unsteady_windway().unwrap();
    assert_eq!(first.channel().velocities_m_s(), second.channel().velocities_m_s());
    assert_eq!(first.channel().last_step(), second.channel().last_step());
    assert_eq!(first.convected_velocity_m_s(), second.convected_velocity_m_s());
    assert_eq!(first.boundary_velocity_m_s(), second.boundary_velocity_m_s());
    assert_eq!(first.channel().gas(), one.voice().gas());
    assert_eq!(first.channel().spec().transverse_cells, 32);
    assert_eq!(first.channel().spec().length_m, LENGTH);
    assert_eq!(first.convection_cells(), 96);
    assert_eq!(one.voice().network().unwrap().last_frame(),
        many.voice().network().unwrap().last_frame());
    let (_, high_viscosity) = render(build(gas(2.0)), 37);
    let (_, quasisteady) = render(voice(gas(1.0), 0.3).with_viscous_windway(LENGTH).unwrap(), 513);
    assert!(a.iter().any(|p| p.abs() > 1e-12));
    assert_ne!(a, quasisteady, "source momentum must reach physical outlet pressure");
    assert_ne!(a, high_viscosity, "material transport must reach physical outlet pressure");
    let d = one.voice().diagnostics().last().unwrap();
    assert!(d.windway_work_j > 0.0 && d.windway_stored_delta_j > 0.0);
    assert!(d.windway_viscous_loss_j > 0.0 && d.windway_exit_work_j > 0.0);
    assert!(d.windway_numerical_dissipation_j > 0.0);
    assert!(d.windway_ledger_defect_j.abs() < 1e-12 * d.windway_work_j);
    assert!(d.ledger_defect_j.abs() < 1e-12, "the mouth still owns its separate ledger");
    assert!(d.reynolds_range.0 < d.reynolds_range.1, "diagnostics must retain the whole attack");
}

#[test]
fn pressure_release_preserves_and_dissipates_the_actual_source_gas() {
    let mut v = voice(gas(1.0), 0.3).with_unsteady_windway(LENGTH, 32, 96).unwrap();
    v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    assert_eq!(v.jet_speed_m_s(), 0.0, "applying force cannot teleport stationary gas");
    v.step_block(&mut [0.0]).unwrap();
    assert!(v.jet_speed_m_s() > 0.0 && v.jet_speed_m_s() < 1.0,
        "startup below the former history-speed floor must be physically admitted");
    v.step_block(&mut [0.0; 512]).unwrap();
    let before = v.unsteady_windway().unwrap().channel().last_step();
    let speed = v.jet_speed_m_s();
    v.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
    assert_eq!(v.jet_speed_m_s(), speed);
    v.step_block(&mut [0.0]).unwrap();
    let after = v.unsteady_windway().unwrap().channel().last_step();
    assert!(after.bulk_speed_m_s > 0.0 && after.bulk_speed_m_s < speed);
    assert!(after.kinetic_energy_j < before.kinetic_energy_j);
    assert_eq!(after.pressure_work_j, 0.0);
    assert!(after.viscous_dissipation_j > 0.0 && after.exit_work_j > 0.0);
    // A slow release may last longer than the old allocated time-history ring.
    // Spatial transport carries the existing field without a clamped delay.
    v.step_block(&mut [0.0; 8192]).unwrap();
    assert!(v.jet_speed_m_s() < 1e-4);
    assert!(v.diagnostics().last().unwrap().windway_stored_delta_j < 0.0);
}

#[test]
fn unsteady_admission_retains_mach_cap_and_refuses_ambiguous_or_late_attachments() {
    let mut v = voice(gas(1.0), 0.05).with_unsteady_windway(LENGTH, 32, 64).unwrap();
    v.validate_blowing_pressure(200.0).unwrap();
    v.validate_blowing_pressure(1.0).unwrap();
    assert!(v.validate_blowing_pressure(300.0).is_err());
    v.apply(FlueControl::SetBlowingPressure(200.0)).unwrap();
    let before = v.unsteady_windway().unwrap().channel().last_step();
    assert!(v.apply(FlueControl::SetBlowingPressure(300.0)).is_err());
    assert_eq!(v.unsteady_windway().unwrap().channel().last_step(), before);
    for cells in [0, 1, 4097] {
        assert!(voice(gas(1.0), 0.3).with_unsteady_windway(LENGTH, cells, 64).is_err());
        assert!(voice(gas(1.0), 0.3).with_unsteady_windway(LENGTH, 32, cells).is_err());
    }
    assert!(voice(gas(1.0), 0.3).with_viscous_windway(LENGTH).unwrap()
        .with_unsteady_windway(LENGTH, 32, 64).is_err());
    assert!(voice(gas(1.0), 0.3).with_unsteady_windway(LENGTH, 32, 64).unwrap()
        .with_viscous_windway(LENGTH).is_err());
    assert!(voice(gas(1.0), 0.3).with_unsteady_windway(LENGTH, 32, 64).unwrap()
        .with_unsteady_windway(LENGTH, 32, 64).is_err());
    v.step_block(&mut [0.0]).unwrap();
    assert!(v.with_unsteady_windway(LENGTH, 32, 64).is_err());
}

#[test]
fn transported_feedback_uses_the_solved_current_mouth_at_physical_entry_time() {
    // The small manufactured cut-up crosses within one source sample. The
    // ordinary cut-up instead exercises boundary entry into interior cells.
    for cut_up in [0.004, 0.00002] {
        let cells = 64;
        let mut v = voice_geometry(gas(1.0), 0.3, cut_up, 0.01)
            .with_unsteady_windway(LENGTH, 32, cells).unwrap();
        v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
        let mut current_outlet_steps = 0;
        let mut live_boundary_steps = 0;
        for _ in 0..256 {
            let before = v.unsteady_windway().unwrap();
            let old = before.convected_velocity_m_s().to_vec();
            let previous_boundary = before.boundary_velocity_m_s();
            v.step_block(&mut [0.0]).unwrap();
            let after = v.unsteady_windway().unwrap();
            let current_boundary = v.mouth_flow_m3_s() / (0.01 * cut_up);
            assert_eq!(after.boundary_velocity_m_s(), current_boundary);
            let travel = after.channel().last_step().bulk_speed_m_s
                * v.island().convection_ratio / f64::from(RATE);
            let dx = cut_up / cells as f64;
            if travel > cut_up { current_outlet_steps += 1; }
            for (index, &actual) in after.convected_velocity_m_s().iter().enumerate() {
                let x = (index + 1) as f64 * dx;
                let departure = x - travel;
                let expected = if departure <= 0.0 {
                    // The exact characteristic entered at t_new-x/c. Its
                    // linear boundary trace interpolates accepted endpoints.
                    live_boundary_steps += usize::from(current_boundary != previous_boundary);
                    current_boundary + (previous_boundary-current_boundary)*(x/travel)
                } else {
                    let position = departure/dx;
                    let lower = position.floor() as usize;
                    let left = if lower == 0 { previous_boundary } else { old[lower-1] };
                    let fraction = position-lower as f64;
                    if fraction == 0.0 { left }
                    else { left+(old[lower]-left)*fraction }
                };
                assert!((actual-expected).abs() <= 2e-12*(1.0+expected.abs()),
                    "cut_up={cut_up}, cell={index}, actual={actual}, expected={expected}");
            }
        }
        assert!(live_boundary_steps > 0);
        assert_eq!(current_outlet_steps > 0, cut_up < 0.0001);
    }
}

#[test]
fn transported_feedback_is_read_at_the_frozen_band_lag_delay() {
    // 1718d79fb reads the frozen history early by the band filters' low-frequency lag,
    // 2/omega_c + 1/omega_h = 16 theta/U, so the loop delay is the convective transit. The
    // transported field is read at the same delay. Through the attack, with the speed new every
    // sample, tap delay plus band lag is W / c_p to the frozen read's 64 epsilon (flue_loop's
    // g1_retarded_affine_signal_has_no_extra_sample_at_any_ring_position), and the tap stays
    // 16 theta c_p/U upstream of the labium.
    let cut_up = 0.004;
    let mut v = voice(gas(1.0), 0.3).with_unsteady_windway(LENGTH, 32, 64).unwrap();
    assert_eq!(v.unsteady_windway().unwrap().feedback_tap_m(), cut_up);
    let theta = v.island().theta_over_b * 0.5 * HEIGHT;
    let ratio = v.island().convection_ratio;
    v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let mut speeds = Vec::new();
    for _ in 0..512 {
        v.step_block(&mut [0.0]).unwrap();
        let windway = v.unsteady_windway().unwrap();
        let u = windway.channel().last_step().bulk_speed_m_s;
        let c = ratio * u;
        let omega_c = u / (4.0 * theta);
        let band_lag = 2.0 / omega_c + 1.0 / (0.5 * omega_c);
        let transit = cut_up / c;
        let tap = windway.feedback_tap_m();
        assert!((tap / c + band_lag - transit).abs() < 64.0 * f64::EPSILON * transit,
            "{u} m/s: tap {tap} m, band lag {band_lag} s, transit {transit} s");
        assert!((cut_up - tap - 16.0 * theta * ratio).abs() < 64.0 * f64::EPSILON * cut_up,
            "{u} m/s: tap {tap} m");
        speeds.push(u);
    }
    assert!(speeds[0] > 0.0 && speeds[0] < 0.5 * speeds[511], "{} {}", speeds[0], speeds[511]);
}

#[test]
fn rejected_mouth_solve_keeps_the_accepted_source_momentum_and_feedback() {
    let mut v = voice(gas(1.0), 0.3).with_unsteady_windway(LENGTH, 32, 64).unwrap()
        .with_mouth_contraction(ContractionSpec { area_ratio: 0.6,
            maximum_contracted_mach: 1e-15, provenance: "manufactured mouth refusal".into() }).unwrap();
    v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let before = v.unsteady_windway().unwrap().clone();
    let speed = v.jet_speed_m_s();
    let wave_energy = v.network().unwrap().stored_energy_j();
    assert!(v.step_block(&mut [0.0]).is_err());
    let after = v.unsteady_windway().unwrap();
    assert_eq!(after.channel().velocities_m_s(), before.channel().velocities_m_s());
    assert_eq!(after.channel().last_step(), before.channel().last_step());
    assert_eq!(after.convected_velocity_m_s(), before.convected_velocity_m_s());
    assert_eq!(after.boundary_velocity_m_s(), before.boundary_velocity_m_s());
    assert_eq!(v.jet_speed_m_s(), speed);
    assert_eq!(v.samples_rendered(), 0);
    assert_eq!(v.network().unwrap().stored_energy_j(), wave_energy);
    // Preserving physical source state does not relax the voice's existing
    // poison policy or promise resumability of unrelated source participants.
    assert!(v.apply(FlueControl::SetBlowingPressure(0.0)).is_err());
    assert!(v.step_block(&mut [0.0]).is_err());
}

#[test]
fn resolved_centre_speed_drives_the_same_sample_jet_and_convected_field() {
    let g = gas(1.0);
    for cells in [31, 32] {
        let mut v = voice_geometry(g, 0.3, 0.004, 0.01)
            .with_unsteady_windway(0.06, cells, 64).unwrap()
            .with_card_centreline_speed(300.0).unwrap();
        let steady = v.unsteady_windway().unwrap().channel()
            .steady_centre_speed_m_s(300.0).unwrap();
        assert_eq!(v.card_centreline_speed().unwrap().held_card_speed_m_s, steady);
        v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
        assert_eq!(v.card_speed_m_s(), 0.0);
        let mut first_ratio = 0.0;
        let mut last_ratio = 0.0;
        for step in 0..320 {
            if step == 256 {
                let speed = v.card_speed_m_s();
                v.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
                assert_eq!(v.card_speed_m_s(), speed, "release cannot reset the momentum field");
            }
            let before = v.unsteady_windway().unwrap();
            let old = before.convected_velocity_m_s().to_vec();
            let previous_boundary = before.boundary_velocity_m_s();
            v.step_block(&mut [0.0]).unwrap();
            let after = v.unsteady_windway().unwrap();
            let frame = after.channel().last_step();
            let centre = after.channel().centre_speed_m_s();
            assert_eq!(v.card_speed_m_s(), centre);
            assert_eq!(v.jet_speed_m_s(), frame.bulk_speed_m_s);
            let expected_flow = 0.01 * HEIGHT * frame.bulk_speed_m_s;
            assert!((frame.flow_m3_s - expected_flow).abs()
                <= 64.0 * f64::EPSILON * expected_flow.abs());
            assert_eq!(v.reynolds(), centre * HEIGHT / (g.dynamic_viscosity / g.density));
            assert!((v.transit_s() * centre * v.island().convection_ratio - 0.004).abs() < 1e-15);
            let ratio = centre / frame.bulk_speed_m_s;
            if step == 0 { first_ratio = ratio; }
            if step == 255 { last_ratio = ratio; }
            let travel = centre * v.island().convection_ratio / f64::from(RATE);
            let current_boundary = after.boundary_velocity_m_s();
            let dx = 0.004 / 64.0;
            for (index, &actual) in after.convected_velocity_m_s().iter().enumerate() {
                let x = (index + 1) as f64 * dx;
                let from = (x - travel) / dx;
                let expected = if from <= 0.0 {
                    current_boundary + (previous_boundary - current_boundary) * (x / travel)
                } else {
                    let lower = from.floor() as usize;
                    let fraction = from - lower as f64;
                    let left = if lower == 0 { previous_boundary } else { old[lower - 1] };
                    if fraction == 0.0 { left } else { left + fraction * (old[lower] - left) }
                };
                assert!((actual - expected).abs() < 2e-12 * (1.0 + expected.abs()),
                    "centre-speed characteristic at cell {index}");
            }
            if step >= 256 {
                assert!(centre > 0.0);
                assert_eq!(frame.pressure_work_j, 0.0);
                assert!(frame.stored_delta_j < 0.0);
            }
        }
        assert!(first_ratio < 1.15);
        assert!(last_ratio > first_ratio + 0.1, "the profile must develop, not use a fixed ratio");
        assert!(v.diagnostics().iter().all(|d| d.windway_ledger_defect_j.abs() < 1e-12));
    }
}

#[test]
fn centreline_material_profile_reaches_microphone_and_preserves_block_partition() {
    let build = |g, centred| {
        let v = voice_geometry(g, 0.3, 0.004, 0.01)
            .with_unsteady_windway(0.06, 32, 96).unwrap();
        if centred { v.with_card_centreline_speed(300.0).unwrap() } else { v }
    };
    let (one, a) = render(build(gas(1.0), true), 513);
    let (many, b) = render(build(gas(1.0), true), 37);
    assert_eq!(a, b);
    assert_eq!(one.voice().card_speed_m_s(), many.voice().card_speed_m_s());
    let (bulk, bulk_audio) = render(build(gas(1.0), false), 37);
    assert_eq!(one.voice().unsteady_windway().unwrap().channel().velocities_m_s(),
        bulk.voice().unsteady_windway().unwrap().channel().velocities_m_s(),
        "centre observation cannot change the mass or energy owner");
    assert_ne!(a, bulk_audio, "source feedback must use the selected physical speed");
    let (viscous, changed_audio) = render(build(gas(2.0), true), 37);
    assert!(viscous.voice().card_speed_m_s() < one.voice().card_speed_m_s());
    assert_ne!(a, changed_audio, "material transport must reach real receiver pressure");
}

#[test]
fn unsteady_centreline_admission_uses_steady_mesh_but_release_keeps_live_state() {
    let build = || voice(gas(1.0), 0.05).with_unsteady_windway(0.06, 32, 64).unwrap();
    let mut bulk = build();
    bulk.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    assert!(bulk.with_card_centreline_speed(100.0).is_err(),
        "a cold receipt cannot admit an already-held centre speed above Mach");
    let mut v = build().with_card_centreline_speed(100.0).unwrap();
    v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
    v.step_block(&mut [0.0; 64]).unwrap();
    let state = v.unsteady_windway().unwrap().channel().velocities_m_s().to_vec();
    assert!(v.apply(FlueControl::SetBlowingPressure(300.0)).is_err());
    assert_eq!(v.unsteady_windway().unwrap().channel().velocities_m_s(), state);
    let centre = v.card_speed_m_s();
    v.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
    assert_eq!(v.card_speed_m_s(), centre);
    v.step_block(&mut [0.0]).unwrap();
    assert!(v.card_speed_m_s() > 0.0);
}
