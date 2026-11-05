use super::*;
use crate::air_path::GasAbsorptionModel;
use crate::pcm_wav::baffled::SurfaceSample;
use crate::render::schedule::{ScheduledControl, ScheduledRenderer};
use crate::render::{ControlDelta, RenderContext, RenderVoice};
use crate::thin_plate::CompactBody;
use fs_material::gas::GasState;
use fs_scenario::RadiatingPlate;

fn air() -> GasState { GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap() }

fn config() -> PlateVoiceConfig {
    PlateVoiceConfig { sample_rate_hz: 48_000, max_modes: 8, nyquist_guard_fraction: 0.9,
        gas: air(), listener_m: 1.0, maximum_abs_force_n: 100.0,
        maximum_abs_pressure_pa: 1e6 }
}

// A genuine zero-area mode has separated moving lobes. A second independent
// mode exercises coordinate ordering, signed force participation and mass.
fn surface(two_modes: bool) -> PlateSurfaceRadiators {
    let mut dark = CompactBody::from_radiator(RadiatingPlate {
        area_m2: 0.02, mass_kg: 0.2, frequency_hz: 1000.0, damping_ratio: 0.01,
    }).unwrap();
    dark.area_m2 = 0.0;
    dark.drive_participation = 0.7;
    let mut bodies = vec![dark];
    let mut left = vec![1.0];
    let mut right = vec![-1.0];
    if two_modes {
        let mut bright = CompactBody::from_radiator(RadiatingPlate {
            area_m2: 0.02, mass_kg: 0.37, frequency_hz: 371.0, damping_ratio: 0.03,
        }).unwrap();
        bright.drive_participation = -0.4;
        bodies.push(bright);
        left.push(0.8);
        right.push(1.2);
    }
    PlateSurfaceRadiators {
        bodies,
        surface: vec![
            SurfaceSample { position_m: [-0.12, 0.0, 0.0], area_m2: 0.01, mode_shape: left },
            SurfaceSample { position_m: [0.12, 0.0, 0.0], area_m2: 0.01, mode_shape: right },
        ],
        maximum_cell_diameter_m: 0.02,
    }
}

fn force(i: usize) -> f64 {
    if i < 13 { 1.0 } else if i < 91 { 0.0 } else if i < 170 { -0.4 } else { 0.0 }
}

fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|v| v.to_bits()).collect() }

fn render(config: PlateVoiceConfig, block: usize) -> Vec<f64> {
    let voice = CompactPlateVoice::from_surface_radiators(
        surface(true), 1.0, config, [0.6, 0.15, 0.5], 1500.0,
    ).unwrap();
    let events = [(13, 0.0), (91, -0.4), (170, 0.0)].into_iter()
        .map(|(sample, force_n)| ScheduledControl {
            sample, delta: ControlDelta::SetPlateForce { voice: 0, force_n },
        }).collect();
    let mut renderer = ScheduledRenderer::new(
        RenderContext::new(vec![RenderVoice::CompactPlate(voice)], 571), events, 3,
    ).unwrap();
    let mut out = vec![0.0; 2048];
    for chunk in out.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    assert_eq!(renderer.samples_rendered(), out.len() as u64);
    out
}

#[test]
fn zero_monopole_area_keeps_its_off_axis_radiation_and_exact_mechanics() {
    let reduction = surface(false);
    let mut compact = CompactPlateVoice::from_radiators(
        reduction.bodies().to_vec(), 1.0, config(),
    ).unwrap();
    let mut axial = CompactPlateVoice::from_surface_radiators(
        reduction.clone(), 1.0, config(), [0.0, 0.0, 0.5], 1500.0,
    ).unwrap();
    let mut side = CompactPlateVoice::from_surface_radiators(
        reduction, 1.0, config(), [0.6, 0.0, 0.5], 1500.0,
    ).unwrap();
    assert!(axial.history.is_empty() && side.history.is_empty());
    let mut peak = 0.0_f64;
    for i in 0..2048 {
        let mut values = [0.0; 3];
        for (voice, slot) in [&mut compact, &mut axial, &mut side].into_iter().zip(&mut values) {
            voice.validate_force(force(i)).unwrap();
            voice.set_force_admitted(force(i));
            voice.step_block(core::slice::from_mut(slot)).unwrap();
        }
        assert_eq!(values[0], 0.0, "compact signed area cancels exactly");
        assert_eq!(values[1], 0.0, "equal axial paths cancel exactly");
        peak = peak.max(values[2].abs());
        let velocity = compact.bodies()[0].modal_velocity_m_s().to_bits();
        assert_eq!(axial.bodies()[0].modal_velocity_m_s().to_bits(), velocity);
        assert_eq!(side.bodies()[0].modal_velocity_m_s().to_bits(), velocity);
    }
    assert!(peak > 1e-5, "separated opposite lobes must radiate off-axis: {peak}");
}

#[test]
fn original_modal_velocities_match_the_shared_rayleigh_owner_bitwise() {
    let reduction = surface(true);
    let position = [0.6, 0.15, 0.5];
    let cfg = config();
    let mut receiver = BaffledPressure::new_differentiating(
        reduction.surface(), 2, cfg.sample_rate_hz, position, RayleighMedium::from(&cfg.gas),
    ).unwrap();
    assert!(receiver.differentiate_pressure().unwrap());
    let mut bodies = reduction.bodies().to_vec();
    let mut voice = CompactPlateVoice::from_surface_radiators(
        reduction, 1.0, cfg, position, 1500.0,
    ).unwrap();
    let receipt = *voice.spatial_observer_receipt().unwrap();
    assert_eq!(receipt.filter_delay_samples, receiver.filter_delay_samples());
    assert_eq!(receipt.filter_delay_samples, 80.0);
    assert_eq!(receipt.multiply_adds_per_sample, receiver.multiply_adds_per_sample());
    assert!(receipt.multiply_adds_per_sample > 0);
    for i in 0..2048 {
        let mut velocities = [0.0; 2];
        for (body, velocity) in bodies.iter_mut().zip(&mut velocities) {
            body.drive(force(i) * body.drive_participation, 1.0 / 48_000.0).unwrap();
            *velocity = body.modal_velocity_m_s();
        }
        let expected = receiver.step_trace(&velocities).unwrap();
        voice.validate_force(force(i)).unwrap();
        voice.set_force_admitted(force(i));
        let mut actual = [0.0];
        voice.step_block(&mut actual).unwrap();
        assert_eq!(actual[0].to_bits(), expected.to_bits(), "sample {i}");
        for (actual, expected) in voice.bodies().iter().zip(&bodies) {
            assert_eq!(actual.modal_velocity_m_s().to_bits(), expected.modal_velocity_m_s().to_bits());
        }
    }
}

#[test]
fn scheduled_surface_pressure_is_block_invariant_without_a_second_compact_delay() {
    let expected = render(config(), 1);
    assert!(expected.iter().any(|v| v.abs() > 1e-5));
    for block in [7, 37, 64, 571] {
        assert_eq!(bits(&render(config(), block)), bits(&expected), "block {block}");
    }
    // A reference distance that would exceed compact history cannot allocate
    // that history or introduce another travel delay into a spatial image.
    let mut distant_reference = config();
    distant_reference.listener_m = 1e6;
    assert_eq!(bits(&render(distant_reference, 37)), bits(&expected));
}

#[test]
fn spatial_receipt_retains_each_original_path_and_the_actual_gas() {
    let position = [0.6, 0.15, 0.5];
    let mut warm = config();
    warm.gas = GasState::try_new_moist_air(333.15, 101_325.0, 0.0).unwrap();
    for cfg in [config(), warm] {
        let reduction = surface(true);
        let mut expected = reduction.surface().iter().map(|s| {
            fs_math::det::sqrt((0..3).map(|i| (position[i] - s.position_m[i]).powi(2)).sum())
                / cfg.gas.sound_speed
        }).collect::<Vec<_>>();
        expected.sort_by(f64::total_cmp);
        let voice = CompactPlateVoice::from_surface_radiators(
            reduction, 1.0, cfg, position, 1500.0,
        ).unwrap();
        let receipt = voice.spatial_observer_receipt().unwrap();
        assert_eq!(receipt.position_m, position);
        assert_eq!(receipt.surface_samples, 2);
        assert_eq!(receipt.maximum_cell_diameter_m, 0.02);
        assert_eq!(receipt.maximum_frequency_hz, 1500.0);
        assert_eq!(receipt.propagation_delay_s, [expected[0], expected[1]]);
        assert!(receipt.propagation_delay_s[1] > receipt.propagation_delay_s[0]);
        assert!(voice.observer_receipt().is_none());
    }
    // Same sound speed but different density: the original Rayleigh pressure
    // changes linearly with rho, while its mechanics and path delays do not.
    let mut thin = config();
    thin.gas = GasState::try_new_moist_air(293.15, 0.7 * 101_325.0, 0.0).unwrap();
    assert_eq!(thin.gas.sound_speed.to_bits(), config().gas.sound_speed.to_bits());
    let scale = thin.gas.density / config().gas.density;
    let original = render(config(), 37);
    let attenuated = render(thin, 37);
    let peak = original.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    for (a, b) in original.iter().zip(&attenuated) {
        assert!((b - scale * a).abs() <= 2e-12 * peak, "{a} {b}");
    }
    assert_ne!(bits(&render(warm, 37)), bits(&original), "changed c must change path phase");
}

#[test]
fn spatial_admission_requires_modal_band_wavelength_and_near_field_resolution() {
    let position = [0.6, 0.15, 0.5];
    for band in [0.0, f64::NAN, 999.0, 0.45 * 48_000.0 + 1.0] {
        assert!(CompactPlateVoice::from_surface_radiators(
            surface(true), 0.0, config(), position, band,
        ).is_err(), "band {band}");
    }
    let mut phase_unresolved = surface(true);
    phase_unresolved.maximum_cell_diameter_m = 0.04;
    assert!(CompactPlateVoice::from_surface_radiators(
        phase_unresolved, 0.0, config(), position, 1500.0,
    ).is_err());
    // This cell resolves the acoustic phase but not the receiver clearance.
    assert!(CompactPlateVoice::from_surface_radiators(
        surface(true), 0.0, config(), [0.6, 0.15, 0.05], 1500.0,
    ).is_err());
    for diameter in [0.0, f64::NAN, f64::INFINITY] {
        let mut invalid = surface(true);
        invalid.maximum_cell_diameter_m = diameter;
        assert!(CompactPlateVoice::from_surface_radiators(
            invalid, 0.0, config(), position, 1500.0,
        ).is_err());
    }
    for position in [[0.6, 0.15, -0.5], [f64::NAN, 0.0, 0.5]] {
        assert!(CompactPlateVoice::from_surface_radiators(
            surface(true), 0.0, config(), position, 1500.0,
        ).is_err());
    }
    assert!(CompactPlateVoice::from_surface_radiators(
        surface(true), 0.0, config(), position, 1000.0,
    ).is_ok(), "every retained natural frequency fits the band");
}

#[test]
fn surface_observer_refuses_single_range_atmosphere_and_late_load_attachment() {
    let make = || CompactPlateVoice::from_surface_radiators(
        surface(true), 1.0, config(), [0.6, 0.15, 0.5], 1500.0,
    ).unwrap();
    for law in [GasAbsorptionModel::Lossless, GasAbsorptionModel::StokesKirchhoff] {
        assert!(make().with_observer_law(law, None).is_err());
    }
    let mut voice = make();
    voice.step_block(&mut [0.0; 1]).unwrap();
    let velocities: Vec<_> = voice.bodies().iter().map(|b| b.modal_velocity_m_s()).collect();
    assert!(voice.attach_piston_loads().is_err());
    assert_eq!(voice.samples_rendered(), 1);
    assert_eq!(velocities, voice.bodies().iter().map(|b| b.modal_velocity_m_s()).collect::<Vec<_>>());
    assert!(voice.step_block(&mut [0.0; 1]).is_ok(), "cold-only refusal does not poison rendering");
}
