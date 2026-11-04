//! Actual loaded-response peaks, including coherent reflected-wave feedthrough.
use fs_math::c64::C64;
use fs_phs::harmonic::resonance::{DrivenResonanceError, ResonanceControls};
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel, ScatteringObservation};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};
use fs_phs::{PortHamiltonian, QuadraticStorage};

fn prepare(system: &PortHamiltonian, ports: &[(usize, f64)]) -> ScatteringModel {
    ScatteringModel::new(
        HarmonicModel::from_ode(
            system,
            HarmonicBudget {
                max_dimension: 32,
                max_dense_bytes: 2_000_000,
                max_dense_work: 4_000_000,
                relative_residual_tolerance: 2e-12,
            },
            || false,
        )
        .unwrap(),
        ports,
        ScatteringBudget {
            maximum_ports: 8,
            maximum_payload_bytes: 500_000,
            maximum_work: 1_000_000,
        },
        || false,
    )
    .unwrap()
}

fn controls(bandwidth: bool) -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-10,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: bandwidth,
    }
}

fn two_ports() -> ScatteringModel {
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.2],
        vec![0.0, 0.0, 1.0, 1.0],
        Box::new(QuadraticStorage::new(vec![100.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap();
    // Deliberately reverse wave order: outgoing weights use selected wave order,
    // while input/output weights always use original physical source order.
    prepare(&system, &[(1, 1.0), (0, 1.0)])
}

fn close(value: C64, expected: C64) {
    assert!(
        (value - expected).abs() < 3e-11 * expected.abs().max(1e-30),
        "{value:?} != {expected:?}"
    );
}

#[test]
fn matched_boundary_resolves_a_lossless_bare_resonance_and_its_physical_width() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 0.0).unwrap();
    let model = prepare(&system, &[(0, 2.0)]);
    let input = [C64::new(0.7, -0.3)];
    let result = model
        .driven_resonance(
            &input,
            ScatteringObservation::Output(&[1.0]),
            (7.0, 13.0),
            controls(true),
            || false,
        )
        .unwrap();
    assert!((result.peak.angular_frequency_rad_s - 10.0).abs() < 1e-10);
    let band = result.half_power.unwrap();
    let gamma = 0.5_f64;
    let root = (gamma * gamma + 400.0).sqrt();
    assert!((band.lower.angular_frequency_rad_s - (root - gamma) / 2.0).abs() < 1e-10);
    assert!((band.upper.angular_frequency_rad_s - (root + gamma) / 2.0).abs() < 1e-10);
    assert!((band.bandwidth_rad_s - gamma).abs() < 1e-10);
    assert!((band.response_quality_factor - 20.0).abs() < 1e-8);
    close(result.observation, input[0].scale(2.0));
    assert!(result.maximum_relative_residual <= 2e-12);
    let peak = model
        .factor(result.peak.angular_frequency_rad_s, || false)
        .unwrap()
        .response(&input, || false)
        .unwrap();
    assert_eq!(peak.physical.dissipated_power_w, 0.0);
    assert!(peak.wave_power_residual_w.abs() < 1e-12);
    // A finite loaded width does not imply material heat in a lossless body.
    assert!((peak.waves[0].incident_power_w - peak.waves[0].outgoing_power_w).abs() < 1e-12);
}

#[test]
fn coherent_incident_feedthrough_shifts_the_outgoing_peak_from_the_motion_peak() {
    let model = two_ports();
    let input = [C64::ONE, C64::new(0.0, 1.0)];
    let gamma = 2.2_f64;
    // With x=100/omega-omega, b0=(2-gamma+i(2-x))/(gamma+i*x).
    // Its squared-amplitude derivative vanishes at x²+0.4*x-gamma²=0.
    let x = (-0.4 - (0.16 + 4.0 * gamma * gamma).sqrt()) / 2.0;
    let expected = ((x * x + 400.0).sqrt() - x) / 2.0;
    let result = model
        .driven_resonance(
            &input,
            ScatteringObservation::Outgoing(&[0.0, 1.0]),
            (10.0, 16.0),
            controls(false),
            || false,
        )
        .unwrap();
    assert!((result.peak.angular_frequency_rad_s - expected).abs() < 1e-10);
    assert!(result.peak.angular_frequency_rad_s > 11.0);
    let omega = result.peak.angular_frequency_rad_s;
    let velocity = C64::new(2.0, 2.0) / C64::new(gamma, 100.0 / omega - omega);
    close(result.observation, velocity - C64::ONE);
    let motion = model
        .driven_resonance(
            &input,
            ScatteringObservation::Output(&[1.0, 0.0]),
            (7.0, 13.0),
            controls(false),
            || false,
        )
        .unwrap();
    assert!((motion.peak.angular_frequency_rad_s - 10.0).abs() < 1e-10);
    // Common phase and excitation strength must not alter the physical peak.
    let phase = C64::new(-0.3, 0.7);
    let scaled = model
        .driven_resonance(
            &[input[0] * phase, input[1] * phase],
            ScatteringObservation::Outgoing(&[0.0, 1.0]),
            (10.0, 16.0),
            controls(false),
            || false,
        )
        .unwrap();
    assert!((scaled.peak.angular_frequency_rad_s - expected).abs() < 1e-10);
    close(scaled.observation, result.observation * phase);
}

#[test]
fn actual_boundary_force_peak_includes_the_reciprocal_load_reaction() {
    let model = two_ports();
    let input = [C64::ONE, C64::new(0.0, 1.0)];
    // |u0|²=4*((gamma-1)²+(x-1)²)/(gamma²+x²), gamma=2.2.
    let gamma = 2.2_f64;
    let a = gamma * gamma - (gamma - 1.0).powi(2) - 1.0;
    let x = (-a - (a * a + 4.0 * gamma * gamma).sqrt()) / 2.0;
    let expected = ((x * x + 400.0).sqrt() - x) / 2.0;
    let result = model
        .driven_resonance(
            &input,
            ScatteringObservation::Input(&[1.0, 0.0]),
            (10.0, 16.0),
            controls(false),
            || false,
        )
        .unwrap();
    assert!((result.peak.angular_frequency_rad_s - expected).abs() < 1e-10);
    let omega = result.peak.angular_frequency_rad_s;
    let velocity = C64::new(2.0, 2.0) / C64::new(gamma, 100.0 / omega - omega);
    close(result.observation, C64::new(2.0, 0.0) - velocity);
    assert!(result.maximum_relative_residual < 2e-12);
}

#[test]
fn loaded_displacement_and_storage_effort_share_the_correct_shifted_peak() {
    let model = two_ports();
    let input = [C64::ONE, C64::new(0.0, 1.0)];
    let expected = (100.0_f64 - 2.2_f64.powi(2) / 2.0).sqrt();
    let displacement = model
        .driven_resonance(
            &input,
            ScatteringObservation::State(&[1.0, 0.0]),
            (7.0, 13.0),
            controls(false),
            || false,
        )
        .unwrap();
    let effort = model
        .driven_resonance(
            &input,
            ScatteringObservation::Effort(&[1.0, 0.0]),
            (7.0, 13.0),
            controls(false),
            || false,
        )
        .unwrap();
    assert!((displacement.peak.angular_frequency_rad_s - expected).abs() < 1e-10);
    assert!((effort.peak.angular_frequency_rad_s - expected).abs() < 1e-10);
    close(effort.observation, displacement.observation.scale(100.0));
}

#[test]
fn reflection_dips_invalid_observations_and_interrupted_searches_refuse() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 0.2).unwrap();
    let model = prepare(&system, &[(0, 1.0)]);
    assert!(matches!(
        model.driven_resonance(
            &[C64::ONE],
            ScatteringObservation::Outgoing(&[1.0]),
            (7.0, 13.0),
            controls(false),
            || false,
        ),
        Err(DrivenResonanceError::NotBracketed(_))
    ));
    assert!(matches!(
        model.driven_resonance(
            &[C64::ONE],
            ScatteringObservation::Input(&[1.0, 0.0]),
            (7.0, 13.0),
            controls(false),
            || false,
        ),
        Err(DrivenResonanceError::Invalid(_))
    ));
    let mut polls = 0;
    assert_eq!(
        model.driven_resonance(
            &[C64::ONE],
            ScatteringObservation::Output(&[1.0]),
            (7.0, 13.0),
            controls(true),
            || {
                polls += 1;
                polls > 120
            },
        ),
        Err(DrivenResonanceError::Harmonic(HarmonicError::Cancelled))
    );
    let result = model
        .driven_resonance(
            &[C64::ONE],
            ScatteringObservation::Output(&[1.0]),
            (7.0, 13.0),
            controls(true),
            || false,
        )
        .unwrap();
    let replay = model
        .driven_resonance(
            &[C64::ONE],
            ScatteringObservation::Output(&[1.0]),
            (7.0, 13.0),
            controls(true),
            || false,
        )
        .unwrap();
    assert_eq!(result, replay);
}
