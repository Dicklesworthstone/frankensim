//! Geometry reduction -> the existing body stepper -> scheduled block pressure.
use fs_couple::render::{CompactPlateVoice, ControlDelta, PlateVoiceConfig, RenderContext, RenderError, RenderVoice};
use fs_couple::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_couple::thin_plate::{CompactBody, certified_radiators};
use fs_scenario::{RadiatingPlate, ThinPlate};
use fs_exec::CancelGate;

fn config() -> PlateVoiceConfig {
    PlateVoiceConfig { sample_rate_hz: 48_000, max_modes: 8, nyquist_guard_fraction: 0.9,
        density_kg_m3: 1.2, listener_m: 1.0, maximum_abs_force_n: 100.0,
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
fn direct(mut bodies: Vec<CompactBody>, n: usize) -> Vec<f64> {
    (0..n).map(|i| {
        let f = if i < 13 { 1.0 } else if i < 91 { 0.0 } else if i < 170 { -0.4 } else { 0.0 };
        let mut pressure = 0.0;
        for body in &mut bodies {
            pressure += body.drive_and_radiate(f*body.drive_participation, 1.0/48_000.0, 1.2, 1.0).unwrap();
        }
        pressure
    }).collect()
}
#[test]
fn signed_footprint_and_radiation_projections_match_the_existing_mechanics_bitwise() {
    let expected = direct(bodies(), 577);
    assert!(expected[13..91].iter().any(|p| p.abs() > 1e-8));
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
    let expected = direct(reduced, 13);
    let mut actual = [0.0; 13];
    voice.step_block(&mut actual).unwrap();
    assert!(actual.iter().any(|x| x.abs() > 1e-10));
    assert_eq!(actual.map(f64::to_bits).as_slice(), expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
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
        let mut actual = [0.0; 13];
        context.block(&mut actual).unwrap();
        let expected = direct(bodies(), 13);
        assert_eq!(actual.map(f64::to_bits).as_slice(), expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
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
fn pressure_refusal_poisoning_prevents_reuse_of_partial_plate_state() {
    let mut c = config();
    c.maximum_abs_pressure_pa = 1e-30;
    let plate = CompactPlateVoice::from_radiators(bodies(), 1.0, c).unwrap();
    let mut context = RenderContext::new(vec![RenderVoice::CompactPlate(plate)], 32);
    assert!(context.block(&mut [0.0; 32]).is_err());
    let mut untouched = [17.0; 32];
    assert!(matches!(context.block(&mut untouched), Err(RenderError::Poisoned)));
    assert_eq!(untouched, [17.0; 32]);
    assert!(matches!(context.apply_controls(&[event(0, 0.0).delta]), Err(RenderError::Poisoned)));
}
