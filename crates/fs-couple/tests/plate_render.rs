//! Geometry reduction -> the existing body stepper -> scheduled block pressure.
use fs_couple::render::{CompactPlateVoice, ControlDelta, PlateVoiceConfig, RenderContext, RenderError, RenderVoice};
use fs_couple::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_couple::thin_plate::{CompactBody, certified_radiators};
use fs_scenario::{RadiatingPlate, ThinPlate};
use fs_exec::CancelGate;
use fs_material::gas::GasState;

fn air() -> GasState { GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap() }
fn config() -> PlateVoiceConfig { config_in(air()) }
fn config_in(gas: GasState) -> PlateVoiceConfig {
    PlateVoiceConfig { sample_rate_hz: 48_000, max_modes: 8, nyquist_guard_fraction: 0.9,
        gas, listener_m: 1.0, maximum_abs_force_n: 100.0,
        maximum_abs_pressure_pa: 1e6 }
}
fn bodies() -> Vec<CompactBody> {
    let mut a = CompactBody::from_radiator(RadiatingPlate {
        area_m2: 0.03, mass_kg: 0.2, frequency_hz: 220.0, damping_ratio: 0.02,
    }).unwrap();
    let mut b = CompactBody::from_radiator(RadiatingPlate {
        area_m2: 0.01, mass_kg: 0.1, frequency_hz: 371.0, damping_ratio: 0.03,
    }).unwrap();
    a.drive_participation = 0.7;
    b.drive_participation = -0.4;
    b.area_m2 = -0.01;
    vec![a,b]
}
fn plate() -> ThinPlate {
    let e = 7e9;
    let nu = 0.3;
    ThinPlate { length_m: 0.4, width_m: 0.3, thickness_m: 0.003,
        density_kg_m3: 700.0, e1_pa: e, e2_pa: e, nu12: nu, g12_pa: e/(2.0*(1.0+nu)),
        material_angle_rad: 0.0, damping_ratio: 0.02, thermoelastic: None,
        kelvin_voigt_bending: None, n_modes: 3, geometric_nonlinearity: false,
        pretension_n_m: 0.0, clamped: false }
}
fn event(sample: u64, force_n: f64) -> ScheduledControl {
    ScheduledControl { sample, delta: ControlDelta::SetPlateForce { voice: 0, force_n } }
}
fn renderer(bodies: Vec<CompactBody>) -> ScheduledRenderer {
    let voice = CompactPlateVoice::from_radiators(bodies, 1.0, config()).unwrap();
    ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::CompactPlate(voice)], 571),
        vec![event(13, 0.0), event(91, -0.4), event(170, 0.0)], 3).unwrap()
}
fn scheduled(i: usize) -> f64 { if i < 13 { 1.0 } else if i < 91 { 0.0 } else if i < 170 { -0.4 } else { 0.0 } }
fn direct(bodies: Vec<CompactBody>, n: usize) -> Vec<f64> { retarded(bodies, n, &air(), scheduled) }
fn direct_in(bodies: Vec<CompactBody>, n: usize, gas: &GasState) -> Vec<f64> { retarded(bodies, n, gas, scheduled) }
/// Samples of the 1 m time of flight in the reference air.
fn onset() -> usize { (48_000.0 / air().sound_speed).floor() as usize }
/// Independent retarded reference through the public stepper: the listener at
/// 1 m hears output i = t_{i+1} - r/c, which lies (1 - f) dt into step i - m of
/// a delay of m + f samples. Step i - m is re-driven for that partial duration
/// from a clone of its start state; output before the first wavefront is zero.
fn retarded(mut bodies: Vec<CompactBody>, n: usize, gas: &GasState, force: fn(usize) -> f64) -> Vec<f64> {
    let delay = 1.0 * 48_000.0 / gas.sound_speed;
    let whole = delay.floor() as usize;
    let elapsed = (1.0 - (delay - delay.floor())) * 48_000.0_f64.recip();
    let mut starts = Vec::with_capacity(n);
    (0..n).map(|i| {
        let f = force(i);
        starts.push((bodies.clone(), f));
        for body in &mut bodies { body.drive(f*body.drive_participation, 1.0/48_000.0).unwrap(); }
        if i < whole { return 0.0; }
        let (start, f) = &starts[i - whole];
        let mut pressure = 0.0;
        for body in start {
            let mut probe = body.clone();
            let acc = probe.drive(f*probe.drive_participation, elapsed).unwrap();
            pressure += probe.radiate(acc, gas.density, 1.0);
        }
        pressure
    }).collect()
}
#[test]
fn signed_footprint_and_radiation_projections_match_the_existing_mechanics_bitwise() {
    let expected = direct(bodies(), 577);
    assert!(expected[onset() + 13..onset() + 91].iter().any(|p| p.abs() > 1e-8));
    for block in [1, 7, 37, 64, 571] {
        let mut r = renderer(bodies());
        let mut actual = vec![0.0; 577];
        for chunk in actual.chunks_mut(block) { r.block(chunk).unwrap(); }
        assert_eq!(actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        assert_eq!(r.samples_rendered(), 577);
    }
}
#[test]
fn the_declared_gas_sets_the_pressure_scale_and_the_time_of_flight() {
    // Only the gas changes. Without a reaction load the mechanics ignore it, so
    // p is exactly rho-proportional at equal c, and arrives r/c after the start.
    let thin = GasState::try_new_moist_air(293.15, 0.7 * 101_325.0, 0.0).unwrap();
    let co2 = fs_material::gas_species::input::read_gas_state(
        &mut "mixture 293.15 101325 mole 1 CO2 1".split_ascii_whitespace(), None).unwrap();
    let render = |gas: GasState| {
        let voice = CompactPlateVoice::from_radiators(bodies(), 1.0, config_in(gas)).unwrap();
        let mut r = ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::CompactPlate(voice)], 571),
            vec![event(13, 0.0), event(91, -0.4), event(170, 0.0)], 3).unwrap();
        let mut out = vec![0.0; 2000];
        for chunk in out.chunks_mut(571) { r.block(chunk).unwrap(); }
        out
    };
    let (a, b, c) = (render(air()), render(thin), render(co2));
    assert_eq!(thin.sound_speed.to_bits(), air().sound_speed.to_bits(), "ideal-gas c is pressure-independent");
    let peak = a.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    assert!(peak > 1e-6);
    let scale = thin.density / air().density;
    assert!((scale - 0.7).abs() < 1e-12);
    for (x, y) in a.iter().zip(&b) { assert!((y - scale * x).abs() <= 1e-12 * peak, "{x} {y}"); }
    for (gas, out) in [(air(), &a), (thin, &b), (co2, &c)] {
        let onset = (48_000.0 / gas.sound_speed).floor() as usize;
        assert!(out[..onset].iter().all(|p| *p == 0.0), "no sound before r/c");
        assert!(out[onset] != 0.0, "the wavefront arrives at r/c");
        assert_eq!(out.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
            direct_in(bodies(), 2000, &gas).iter().map(|p| p.to_bits()).collect::<Vec<_>>());
    }
    // CO2 is slower than air: 1 m takes about 179 rather than 140 samples.
    assert!((48_000.0 / co2.sound_speed).floor() > (48_000.0 / air().sound_speed).floor() + 30.0);
}
#[test]
fn a_listener_beyond_the_retarded_history_budget_refuses() {
    let mut far = config();
    far.listener_m = 1e6;
    assert!(CompactPlateVoice::from_radiators(bodies(), 0.0, far).is_err());
    let mut silent = config();
    silent.gas.sound_speed = f64::NAN;
    assert!(CompactPlateVoice::from_radiators(bodies(), 0.0, silent).is_err());
}
#[test]
fn actual_plate_geometry_reaches_the_renderer_without_hand_authored_frequencies() {
    let geometry = plate();
    let reduced = certified_radiators(geometry).unwrap();
    let mut voice = CompactPlateVoice::from_plate(geometry, 1.0, config()).unwrap();
    assert_eq!(voice.bodies().len(), reduced.len());
    for (a,b) in voice.bodies().iter().zip(&reduced) {
        assert_eq!(a.omega.to_bits(), b.omega.to_bits());
        assert_eq!(a.mass_kg.to_bits(), b.mass_kg.to_bits());
        assert_eq!(a.drive_participation.to_bits(), b.drive_participation.to_bits());
    }
    // Constant initial force; 13 samples of sound after the 1 m time of flight.
    let expected = retarded(reduced, onset() + 13, &air(), |_| 1.0);
    let mut actual = vec![0.0; onset() + 13];
    voice.step_block(&mut actual).unwrap();
    assert!(actual[onset()..].iter().any(|x| x.abs() > 1e-10));
    assert_eq!(actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    let mut thicker = geometry;
    thicker.thickness_m *= 1.5;
    let other = CompactPlateVoice::from_plate(thicker, 1.0, config()).unwrap();
    assert!((other.bodies()[0].omega / voice.bodies()[0].omega - 1.5).abs() < 0.02);
}
#[test]
fn invalid_control_batch_preserves_force_state_clock_and_log() {
    for invalid in [event(0, f64::NAN).delta, event(0, 101.0).delta,
        ControlDelta::SetPlateForce { voice: 9, force_n: 1.0 },
        ControlDelta::SetBlowingPressure { voice: 0, pressure_pa: 1.0 },
        ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 1.0 }]
    {
        let mut context = renderer(bodies()).into_context();
        assert!(context.apply_controls(&[event(0, -3.0).delta, invalid]).is_err());
        assert!(context.control_log().is_empty());
        let mut actual = vec![0.0; onset() + 13];
        context.block(&mut actual).unwrap();
        let expected = retarded(bodies(), onset() + 13, &air(), |_| 1.0);
        assert!(actual[onset()..].iter().any(|x| *x != 0.0));
        assert_eq!(actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    }
}
#[test]
fn cancellation_resumes_plate_ringdown_and_keeps_boundary_events_pending() {
    let expected = direct(bodies(), 577);
    let mut r = renderer(bodies());
    let mut actual = vec![0.0; 577];
    r.block(&mut actual[..13]).unwrap();
    let cancelled = CancelGate::new();
    cancelled.request();
    let mut untouched = [123.0; 64];
    let outcome = r.render_under_gate(&cancelled, &mut untouched, 64, 1).unwrap();
    assert!(matches!(outcome, fs_couple::render::GatedRenderOutcome::Cancelled { blocks: 0 }));
    assert_eq!(untouched, [123.0; 64]);
    assert_eq!(r.pending_controls()[0].sample, 13);
    for chunk in actual[13..].chunks_mut(37) { r.block(chunk).unwrap(); }
    assert_eq!(actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
}
#[test]
fn plate_clock_nyquist_and_nonlinear_admission_refuse_without_silently_retuning() {
    let mut r = renderer(bodies());
    assert!(r.validate_sample_rate(48_000).is_ok());
    assert!(r.validate_sample_rate(44_100).is_err());
    assert_eq!(r.samples_rendered(), 0);
    let mut aliased = bodies();
    aliased[0].omega = core::f64::consts::TAU * 24_000.0;
    assert!(CompactPlateVoice::from_radiators(aliased, 0.0, config()).is_err());
    let mut nonlinear = plate();
    nonlinear.geometric_nonlinearity = true;
    assert!(CompactPlateVoice::from_plate(nonlinear, 0.0, config()).is_err());
    let mut bad = config();
    bad.max_modes = 0;
    assert!(CompactPlateVoice::from_radiators(bodies(), 0.0, bad).is_err());
    // A refused sink clock does not poison a valid renderer.
    r.block(&mut [0.0; 1]).unwrap();
}

#[test]
fn retained_reaction_must_match_the_observer_gas_and_clock() {
    let loaded = CompactPlateVoice::from_radiators_with_piston_loads(bodies(), 1.0, config()).unwrap();
    let bodies = loaded.bodies().to_vec();
    let mut different_gas = config();
    different_gas.gas = GasState::try_new_moist_air(293.15, 202_650.0, 0.0).unwrap();
    assert!(CompactPlateVoice::from_radiators(bodies.clone(), 0.0, different_gas).is_err());
    let mut different_clock = config();
    different_clock.sample_rate_hz = 96_000;
    assert!(CompactPlateVoice::from_radiators(bodies.clone(), 0.0, different_clock).is_err());
    let mut incoherent = config();
    incoherent.gas.sound_speed *= 2.0;
    assert!(CompactPlateVoice::from_radiators(bodies.clone(), 0.0, incoherent).is_err());
    assert!(CompactPlateVoice::from_radiators(bodies, 0.0, config()).is_ok());
}

#[test]
fn retained_piston_load_keeps_its_gas_clock_and_memory() {
    let gas = air();
    let other = GasState::try_new_moist_air(293.15, 202_650.0, 0.0).unwrap();
    let make = || {
        let mut bank = fs_couple::thin_plate::PlateBank::default();
        bank.linear.push(CompactBody::from_radiator(RadiatingPlate {
            area_m2: 0.01, mass_kg: 0.1, frequency_hz: 220.0, damping_ratio: 0.001,
        }).unwrap());
        bank.attach_radiation_loads(&gas, 48_000).unwrap();
        bank
    };
    let mut bank = make();
    let mut reference = make();
    let dt = 48_000.0_f64.recip();
    for _ in 0..19 {
        bank.advance_controlled(1.0, dt, gas.density, 1.0, || false).unwrap();
        reference.advance_controlled(1.0, dt, gas.density, 1.0, || false).unwrap();
    }
    assert!(bank.attach_radiation_loads(&other, 48_000).is_err());
    assert!(bank.attach_radiation_loads(&gas, 96_000).is_err());
    bank.attach_radiation_loads(&gas, 48_000).unwrap();
    for _ in 0..20 {
        let full = bank.advance_controlled(0.0, dt, gas.density, 1.0, || false).unwrap();
        let first = reference.advance_controlled(0.0, dt*0.5, gas.density, 1.0, || false).unwrap();
        let second = reference.advance_controlled(0.0, dt*0.5, gas.density, 1.0, || false).unwrap();
        let scale = full.volume_velocity_m3_s.abs().max(1e-12);
        assert!((full.volume_velocity_m3_s-second.volume_velocity_m3_s).abs() < 1e-10*scale);
        let heat = first.work.dissipated+second.work.dissipated;
        assert!((full.work.dissipated-heat).abs() < 1e-10*heat.max(1e-18));
    }
    assert!((bank.stored_energy_j()-reference.stored_energy_j()).abs()
        < 1e-10*bank.stored_energy_j());
}

#[test]
fn selected_gas_reaction_changes_body_kinetic_energy_not_only_output() {
    let dense = GasState::try_new_moist_air(293.15, 202_650.0, 0.0).unwrap();
    let energy = |gas, loaded| {
        let mut voice = if loaded {
            CompactPlateVoice::from_radiators_with_piston_loads(bodies(), 1.0, config_in(gas))
        } else {
            CompactPlateVoice::from_radiators(bodies(), 1.0, config_in(gas))
        }.unwrap();
        voice.step_block(&mut [0.0; 2048]).unwrap();
        // The retained body's actual velocity, independent of the observer's
        // density scale and time-of-flight: Q = modal_area * q_dot.
        voice.bodies().iter().map(|body|
            0.5*body.mass_kg*(body.volume_velocity()/body.area_m2).powi(2)
        ).sum::<f64>()
    };
    assert_eq!(energy(air(), false), energy(dense, false));
    let (ordinary, doubled) = (energy(air(), true), energy(dense, true));
    assert!(ordinary > 0.0 && doubled > 0.0);
    assert!((ordinary-doubled).abs()/ordinary.max(doubled) > 1e-5,
        "the actual body's kinetic energy must respond to the selected gas: {ordinary}, {doubled}");
}
#[test]
fn pressure_refusal_poisoning_prevents_reuse_of_partial_plate_state() {
    let mut c = config();
    c.maximum_abs_pressure_pa = 1e-30;
    let plate = CompactPlateVoice::from_radiators(bodies(), 1.0, c).unwrap();
    // The first nonzero pressure arrives after the 1 m time of flight.
    let mut context = RenderContext::new(vec![RenderVoice::CompactPlate(plate)], 256);
    assert!(onset() < 256);
    assert!(context.block(&mut [0.0; 256]).is_err());
    let mut untouched = [17.0; 32];
    assert!(matches!(context.block(&mut untouched), Err(RenderError::Poisoned)));
    assert_eq!(untouched, [17.0; 32]);
    assert!(matches!(context.apply_controls(&[event(0, 0.0).delta]), Err(RenderError::Poisoned)));
}
