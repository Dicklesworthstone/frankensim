use super::*;
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use crate::render::plate::impact::{BodyPotential, linear::LinearImpactConfig};
use crate::render::schedule::force::coupled::{ModalCouplingConfig,
    contact::{ModalContactConfig, multiple::MultiContactConfig}};

const RHO: f64 = 1.2;
const AREA: f64 = 1.0e-4;

fn monopole(maximum_hz: f64) -> ImpactRadiation {
    // Independent point-monopole Green function: F/Q = -i*omega*rho/(4*pi).
    // Q is an actual volume velocity; no frequency-dependent amplitude is
    // supplied by the mechanical source or inferred from its mode frequency.
    let diagnostic = RadiationSampleDiagnostics { captured_fraction: 1.0,
        panels_per_wavelength: 100.0, condition_lower_bound: 1.0 };
    let root = (4.0 * core::f64::consts::PI).sqrt();
    let transfer = |omega: f64| C64::new(0.0, -omega * RHO / (4.0 * core::f64::consts::PI));
    let directions = vec![[1.0,0.0,0.0], [-1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,-1.0,0.0],
        [0.0,0.0,1.0], [0.0,0.0,-1.0], [1.0,1.0,1.0], [-1.0,1.0,-1.0]];
    let samples = SampledRadiationData {
        source_id: "analytic compact monopole, volume-flow input".into(),
        harmonic_time_convention: HarmonicTimeConvention::ExpNegativeIOmegaT,
        l_max: 0, input_ids: vec!["volume_velocity_m3_s".into()],
        training: (1..=24).map(|i| {
            let omega = core::f64::consts::TAU * (maximum_hz * (f64::from(i) / 24.0));
            ComplexShTrainingSample { omega_rad_s: omega,
                coefficients_by_input: vec![vec![transfer(omega).scale(root)]], diagnostics: diagnostic }
        }).collect(),
        held_out: (0..24).map(|i| {
            let omega = core::f64::consts::TAU * maximum_hz * (f64::from(i) + 0.5) / 24.0;
            DirectFarFieldHeldOutSample { omega_rad_s: omega,
                directions: directions.clone(),
                far_field_by_input: vec![vec![transfer(omega); directions.len()]], diagnostics: diagnostic }
        }).collect(),
    };
    let mut controls = controls();
    controls.maximum_normalized_error = 1e-8;
    controls.rms_normalized_error = 1e-8;
    ImpactRadiation::from_velocity_samples(&samples, controls, 0.0, 343.0).unwrap()
}

fn projection() -> Vec<VelocityProjection> {
    // This test's mass-normalized coordinate has mass 1 kg. A physical area
    // converts its normal velocity [m/s] to the acoustic input Q [m^3/s].
    vec![VelocityProjection { input_id: "volume_velocity_m3_s".into(), weights: vec![AREA] }]
}

fn oscillator(frequency_hz: f64, velocity_amplitude: f64, steps: u64) -> LinearImpactSystem {
    let omega = core::f64::consts::TAU * frequency_hz;
    let body = ImpactBody { potential: BodyPotential::Linear(vec![omega]),
        initial: vec![ModalAcousticState { displacement_m_sqrt_kg: velocity_amplitude / omega,
            velocity_m_sqrt_kg_per_s: 0.0 }], damping_per_s: vec![0.0] };
    let mut component = ModalAcousticTimeBudget::audible_reference();
    component.nyquist_guard_fraction = 0.95;
    let config = LinearImpactConfig { sample_rate_hz: RATE, max_steps: steps,
        maximum_generalized_force: 1000.0, component,
        coupling: ModalCouplingConfig { max_modes: 4, max_connections: 0, max_setup_terms: 10_000,
            nyquist_guard_fraction: 0.95, maximum_total_energy_j: 10.0,
            maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e4,
            solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-10,
            energy_relative_tolerance: 1e-8 },
        contact: ModalContactConfig { max_iterations: 100, maximum_force_n: 1e4,
            maximum_penetration_m: 0.01, force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9 },
        multiple: MultiContactConfig { max_contacts: 4, max_sweeps: 100, max_setup_terms: 10_000 },
    };
    LinearImpactSystem::new(vec![body], vec![], vec![], config, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn physical_high_harmonics_keep_derivative_gain_phase_and_original_mechanics() {
    let radiation = monopole(0.45 * f64::from(RATE));
    let range = 343.0 * 32.0 / f64::from(RATE);
    let amplitude = 1e-3;
    let gate = CancelGate::new_clock_free();
    for fraction in [0.05, 0.2, 0.4, 0.45] {
        let frequency = fraction * f64::from(RATE);
        let omega = core::f64::consts::TAU * frequency;
        let mut direct = oscillator(frequency, amplitude, 2048);
        let mut observed = ImpactPressureRenderer::new(oscillator(frequency, amplitude, 2048),
            &radiation, projection(), vec![0.0], listener(range), RATE, 97).unwrap();
        let timing = observed.timing();
        assert!((timing.velocity_difference_delay_s - 48.0 / f64::from(RATE)).abs() < 1e-15);
        assert_eq!(timing.source_time_shift_s, 0.0);
        let delay = (timing.remaining_travel_time_s + timing.velocity_difference_delay_s
            + timing.propagation_interpolation_delay_s) * f64::from(RATE);
        let scale = RHO * AREA / (4.0 * core::f64::consts::PI * range);
        let (mut error, mut signal) = (0.0, 0.0);
        for n in 0..2048 {
            let direct_frame = direct.step(&[0.0], &gate).unwrap();
            let mut pressure = [0.0];
            observed.block(&mut pressure).unwrap();
            assert_eq!(observed.mechanics().state(), direct.state());
            assert_eq!(observed.last_mechanical_frame().unwrap(), &direct_frame);
            if n >= 256 {
                let time = ((n + 1) as f64 - delay) / f64::from(RATE);
                let expected = -scale * omega * amplitude * (omega * time).cos();
                error += (pressure[0] - expected).powi(2);
                signal += expected * expected;
            }
        }
        let relative = (error / signal).sqrt();
        assert!(relative < 2.5e-4, "fraction {fraction}: complex waveform error {relative}");
        // The former first difference loses 24-30% at these upper harmonics;
        // a low-frequency-only or magnitude-only test would conceal that loss.
        if fraction >= 0.4 {
            let old_gain = (core::f64::consts::PI * fraction).sin()
                / (core::f64::consts::PI * fraction);
            assert!(1.0 - old_gain > 0.24);
        }
    }
}

#[test]
fn fractional_listener_distance_keeps_physical_phase_and_separate_numerical_latency() {
    let radiation = monopole(0.45 * f64::from(RATE));
    let amplitude = 1e-3;
    for travel in [0.125, 32.125, 32.5, 32.875] {
        let range = 343.0 * travel / f64::from(RATE);
        for fraction in [0.05, 0.2, 0.4, 0.45] {
            let omega = core::f64::consts::TAU * fraction * f64::from(RATE);
            let mut receiver = listener(range);
            receiver.maximum_delay_error_s = 0.0;
            let mut observed = ImpactPressureRenderer::new(
                oscillator(fraction * f64::from(RATE), amplitude, 1024),
                &radiation, projection(), vec![0.0], receiver, RATE, 97).unwrap();
            let timing = observed.timing();
            assert_eq!(timing.propagation_rounding_s, 0.0);
            assert!((timing.propagation_samples as f64 + timing.propagation_fractional_samples
                - travel).abs() < 1e-13);
            assert!((timing.propagation_interpolation_delay_s * f64::from(RATE) - 32.0).abs() < 1e-13);
            let scale = RHO * AREA / (4.0 * core::f64::consts::PI * range);
            let (mut error, mut signal) = (0.0, 0.0);
            for n in 0..1024 {
                let mut pressure = [0.0];
                observed.block(&mut pressure).unwrap();
                if n >= 256 {
                    // Physical time uses the original geometry; the separate
                    // derivative/interpolation delays are retained, never removed.
                    let time = (n + 1) as f64 / f64::from(RATE)
                        - range / 343.0 - 80.0 / f64::from(RATE);
                    let expected = -scale * omega * amplitude * (omega * time).cos();
                    error += (pressure[0] - expected).powi(2);
                    signal += expected * expected;
                }
            }
            let relative = (error / signal).sqrt();
            assert!(relative < 5e-4, "travel {travel}, fraction {fraction}: {relative}");
        }
    }
}

#[test]
fn coasting_velocity_has_no_pressure_and_reconstruction_never_changes_motion() {
    let radiation = monopole(0.45 * f64::from(RATE));
    let mut direct = system(512);
    let mut observed = ImpactPressureRenderer::new(system(512), &radiation,
        projection(), vec![0.0], listener(1.0), RATE, 97).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut peak = 0.0_f64;
    for n in 0..512 {
        let force = if n < 16 { 100.0 } else { 0.0 };
        direct.step(&[force], &gate).unwrap();
        observed.set_forces(&[force]).unwrap();
        let mut pressure = [0.0];
        observed.block(&mut pressure).unwrap();
        assert_eq!(observed.mechanics().state(), direct.state());
        peak = peak.max(pressure[0].abs());
        if n > 400 { assert!(pressure[0].abs() < 1e-12); }
    }
    assert!(peak > 1e-7);
    assert!(observed.mechanics().state()[1] > 0.01);
}

#[test]
fn fitted_bands_above_the_reconstruction_limit_refuse() {
    let radiation = monopole(0.46 * f64::from(RATE));
    let result = ImpactPressureRenderer::new(system(8), &radiation,
        projection(), vec![0.0], listener(1.0), RATE, 97);
    assert!(matches!(result, Err(RenderError::Control { what })
        if what.contains("0.45-rate")));
}
