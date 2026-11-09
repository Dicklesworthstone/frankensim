//! Physical notch tuning and absolute recovery-edge derivatives, including exact zeros.
use fs_math::c64::C64;
use fs_phs::{PortHamiltonian, QuadraticStorage};
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection as Direction,
};
use fs_phs::harmonic::resonance::{
    AntiresonanceControls, DrivenResonanceError, HarmonicObservation,
};
use fs_phs::harmonic::scattering::{
    ScatteringBudget, ScatteringModel, ScatteringObservation,
};

const OBSERVATION: [f64; 4] = [1.0, 0.0, 0.0, 0.0];

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 8,
        max_dense_bytes: 2_000_000,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 2e-12,
    }
}

fn controls(level: Option<f64>) -> AntiresonanceControls {
    AntiresonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-11,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        recovery_power: level,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}"
    );
}

// Original [q_host, p_host, q_absorber, p_absorber] coordinates.
// M_host=2, M_absorber=.5, K_host=200, K_absorber=50, C_host=.4.
// Only the original host momentum is forced; the absorber remains a full state.
fn absorber() -> PortHamiltonian {
    let q = vec![
        250.0, 0.0, -50.0, 0.0,
        0.0, 0.5, 0.0, 0.0,
        -50.0, 0.0, 50.0, 0.0,
        0.0, 0.0, 0.0, 2.0,
    ];
    let j = vec![
        0.0, 1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, -1.0, 0.0,
    ];
    let mut r = vec![0.0; 16];
    r[5] = 0.4;
    PortHamiltonian::new(
        4, 1, j, r, vec![0.0, 1.0, 0.0, 0.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    ).unwrap()
}

// Independent inverse of the original 2x2 mechanical dynamic stiffness.
// The native pHS owner uses exp(-i omega t), so damping is -i omega C.
fn pressure(omega: f64, force: C64) -> (C64, C64) {
    let host = C64::new(250.0 - 2.0 * omega * omega, -0.4 * omega);
    let host_slope = C64::new(-4.0 * omega, -0.4);
    let absorber = 50.0 - 0.5 * omega * omega;
    let absorber_slope = -omega;
    let determinant = host.scale(absorber) - C64::from_re(2500.0);
    let determinant_slope = host_slope.scale(absorber) + host.scale(absorber_slope);
    let value = force.scale(absorber) / determinant;
    let slope = (force.scale(absorber_slope) - value * determinant_slope) / determinant;
    (value, slope)
}

fn recovery_root(mut low: f64, mut high: f64, force: C64, level: f64) -> f64 {
    let below = pressure(low, force).0.norm_sq() < level;
    assert_ne!(below, pressure(high, force).0.norm_sq() < level);
    for _ in 0..80 {
        let middle = f64::midpoint(low, high);
        if (pressure(middle, force).0.norm_sq() < level) == below {
            low = middle;
        } else {
            high = middle;
        }
    }
    f64::midpoint(low, high)
}

fn loaded(c: f64) -> ScatteringModel {
    let source = fs_phs::mass_spring_damper(1.0, 100.0, c).unwrap();
    ScatteringModel::new(
        HarmonicModel::from_ode(&source, budget(), || false).unwrap(),
        &[(0, 2.0)],
        ScatteringBudget {
            maximum_ports: 1,
            maximum_payload_bytes: 500_000,
            maximum_work: 1_000_000,
        },
        || false,
    ).unwrap()
}

#[test]
fn simple_zero_tuning_retains_original_storage_and_damping_at_extreme_source_scales() {
    let model = HarmonicModel::from_ode(&absorber(), budget(), || false).unwrap();
    let mut stiffness = [0.0; 16];
    stiffness[0] = 1.0;
    stiffness[2] = -1.0;
    stiffness[8] = -1.0;
    stiffness[10] = 1.0;
    let mut mass = [0.0; 16];
    mass[15] = -4.0; // d(1/m_absorber)/dm = -1/m_absorber^2.
    let mut absorber_damping = [0.0; 16];
    absorber_damping[15] = 1.0;
    let unit_force = C64::new(0.6, -0.8);
    for scale in [1e-120, 1.0, 1e120] {
        let force = unit_force.scale(scale);
        for (direction, frequency_slope, depth_slope) in [
            (Direction { storage: Some(&stiffness), ..Default::default() }, 0.1, C64::ZERO),
            (Direction { storage: Some(&mass), ..Default::default() }, -10.0, C64::ZERO),
            (Direction { dissipation: Some(&absorber_damping), ..Default::default() },
                0.0, C64::new(0.0, 0.004) * unit_force),
        ] {
            let result = model.driven_antiresonance_with_parameter_derivative(
                &[force], HarmonicObservation::State(&OBSERVATION), direction, None,
                (9.0, 11.0), controls(Some(1e-6 * scale * scale)), || false,
            ).unwrap();
            // omega0=sqrt(k_absorber/m_absorber); changing absorber loss
            // opens a complex notch depth to first order without detuning it.
            near(result.antiresonance.minimum.angular_frequency_rad_s, 10.0, 3e-11);
            near(result.minimum_frequency_derivative_rad_s, frequency_slope, 3e-9);
            complex_near(result.antiresonance.observation.scale(1.0 / scale), C64::ZERO, 2e-10);
            complex_near(result.observation_derivative.scale(1.0 / scale), depth_slope, 3e-9);
            // Recover physical curvature after cancelling source units FIRST,
            // rather than squaring a 1e120 or 1e-120 conditioning scale.
            let a = result.curvature_amplitude_scale / scale;
            near(result.scaled_power_curvature_s2 * a * a, 2.0 * 0.004_f64.powi(2), 3e-10);
            assert!(result.scaled_power_curvature_s2 > 0.0);
            assert!(result.scaled_stationarity_residual_s.abs() < 1e-8);
            assert_eq!(result.evaluations, result.antiresonance.evaluations + 3);
            assert!(result.evaluations <= controls(Some(1.0)).maximum_evaluations);
            assert!(result.maximum_relative_residual <= 2e-12);
        }
    }
}

#[test]
fn source_gain_moves_absolute_recovery_edges_while_phase_and_exact_zero_remain_resolved() {
    let model = HarmonicModel::from_ode(&absorber(), budget(), || false).unwrap();
    let unit_force = C64::new(0.6, -0.8);
    let level = 1e-6;
    let lower = recovery_root(9.0, 10.0, unit_force, level);
    let upper = recovery_root(10.0, 11.0, unit_force, level);
    for alpha in [C64::new(0.3, -0.4), C64::new(0.0, 0.7)] {
        let edge_slope = |omega| {
            let (value, slope) = pressure(omega, unit_force);
            -alpha.re * value.norm_sq() / (value.conj() * slope).re
        };
        let dl = edge_slope(lower);
        let du = edge_slope(upper);
        if alpha.re != 0.0 {
            assert!(dl > 0.0 && du < 0.0); // A fixed absolute threshold narrows.
        }
        for scale in [1e-120, 1.0, 1e120] {
            let force = unit_force.scale(scale);
            let result = model.driven_antiresonance_with_parameter_derivative(
                &[force], HarmonicObservation::State(&OBSERVATION), Direction::default(),
                Some(&[alpha * force]), (9.0, 11.0),
                controls(Some(level * scale * scale)), || false,
            ).unwrap();
            near(result.minimum_frequency_derivative_rad_s, 0.0, 3e-9);
            complex_near(result.observation_derivative.scale(1.0 / scale), C64::ZERO, 3e-9);
            let band = result.antiresonance.recovery.unwrap();
            near(band.lower.angular_frequency_rad_s, lower, 3e-10);
            near(band.upper.angular_frequency_rad_s, upper, 3e-10);
            let derivative = result.recovery_derivative.unwrap();
            near(derivative.lower_frequency_derivative_rad_s, dl, 3e-9);
            near(derivative.upper_frequency_derivative_rad_s, du, 3e-9);
            near(derivative.bandwidth_derivative_rad_s, du - dl, 5e-9);
        }
    }
}

#[test]
fn nonunit_wave_reference_retains_exact_zero_depth_and_absolute_band_tuning() {
    // Z is velocity/force admittance: matching requires c=1/Z, not c=Z.
    let (m, k, c, z) = (1.0_f64, 100.0_f64, 0.5_f64, 2.0_f64);
    let model = loaded(c);
    let input = C64::new(0.6, -0.8);
    let rho = 0.25_f64;
    let s = (rho / (1.0 - rho)).sqrt();
    let x = 2.0 * s / z;
    let root = (x * x + 4.0 * m * k).sqrt();
    let lower = (root - x) / (2.0 * m);
    let upper = (root + x) / (2.0 * m);
    for (dm, dk, dc, dz, alpha) in [
        (1.0, 0.0, 0.0, 0.0, C64::ZERO),
        (0.0, 1.0, 0.0, 0.0, C64::ZERO),
        (0.0, 0.0, 1.0, 0.0, C64::ZERO),
        (0.0, 0.0, 0.0, 1.0, C64::ZERO),
        (0.0, 0.0, 0.0, 0.0, C64::new(0.3, -0.4)),
        (0.2, 0.7, -0.1, 0.3, C64::new(-0.2, 0.5)),
    ] {
        let dq = [dk, 0.0, 0.0, -dm / (m * m)];
        let dr = [0.0, 0.0, 0.0, dc];
        let result = model.driven_antiresonance_with_parameter_derivative(
            &[input], ScatteringObservation::Outgoing(&[1.0]),
            Direction { storage: Some(&dq), dissipation: Some(&dr), ..Default::default() },
            Some(&[alpha * input]), Some(&[dz]), (8.0, 12.0),
            controls(Some(rho)), || false,
        ).unwrap();
        // |b/a|²=((1-Zc)²+Z²x²)/((1+Zc)²+Z²x²), x=k/omega-m*omega.
        // At c=1/Z, X=2sqrt(rho/(1-rho))/Z. The reference derivative
        // here holds physical c fixed; it is not a rematched boundary path.
        let dx = s * (dc - dz / (z * z)) - x * alpha.re / (1.0 - rho);
        let droot = (x * dx + 2.0 * (k * dm + m * dk)) / root;
        let dl = (droot - dx) / (2.0 * m) - lower * dm / m;
        let du = (droot + dx) / (2.0 * m) - upper * dm / m;
        near(result.antiresonance.minimum.angular_frequency_rad_s, 10.0, 3e-11);
        near(result.minimum_frequency_derivative_rad_s, 5.0 * (dk / k - dm / m), 3e-9);
        complex_near(result.antiresonance.observation, C64::ZERO, 2e-10);
        complex_near(result.observation_derivative, input.scale(-0.5 * (z * dc + c * dz)), 3e-9);
        let band = result.antiresonance.recovery.unwrap();
        near(band.lower.angular_frequency_rad_s, lower, 3e-10);
        near(band.upper.angular_frequency_rad_s, upper, 3e-10);
        let derivative = result.recovery_derivative.unwrap();
        near(derivative.lower_frequency_derivative_rad_s, dl, 3e-9);
        near(derivative.upper_frequency_derivative_rad_s, du, 3e-9);
        near(derivative.bandwidth_derivative_rad_s, dx / m - x * dm / (m * m), 5e-9);
        assert_eq!(result.evaluations, result.antiresonance.evaluations + 3);
        assert!(result.maximum_relative_residual <= 2e-12);
    }
}

#[test]
fn zero_directions_silent_or_flat_experiments_and_native_caps_keep_their_refusals() {
    let system = absorber();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let run = || model.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
        (9.0, 11.0), controls(None), || false,
    ).unwrap();
    let zero = run();
    assert_eq!(zero.minimum_frequency_derivative_rad_s, 0.0);
    assert_eq!(zero.observation_derivative, C64::ZERO);
    assert_eq!(zero.recovery_derivative, None);
    assert_eq!(zero.evaluations, zero.antiresonance.evaluations + 1);
    assert_eq!(zero, run());
    assert!(matches!(model.driven_antiresonance_with_parameter_derivative(
        &[C64::ZERO], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
        (9.0, 11.0), controls(None), || false,
    ), Err(DrivenResonanceError::NotBracketed(_))));
    assert!(model.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], HarmonicObservation::State(&[0.0; 4]), Direction::default(), None,
        (9.0, 11.0), controls(None), || false,
    ).is_err());
    assert!(model.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), Some(&[]),
        (9.0, 11.0), controls(None), || false,
    ).is_err());
    for (bytes, work) in [(4000, 4_000_000), (2_000_000, 1500)] {
        // The complete 4-state constructor fits (2592 bytes, 400 work),
        // but the simultaneous second/mixed jet query does not (4896, 2400).
        let limited = HarmonicModel::from_ode(&system,
            HarmonicBudget { max_dense_bytes: bytes, max_dense_work: work, ..budget() },
            || false,
        ).unwrap();
        assert!(matches!(limited.driven_antiresonance_with_parameter_derivative(
            &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
            (9.0, 11.0), controls(None), || false,
        ), Err(DrivenResonanceError::Invalid(_))));
    }
    let flat = loaded(0.0); // A lossless one-port reflector has constant |b/a|=1.
    assert!(flat.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], ScatteringObservation::Outgoing(&[1.0]), Direction::default(),
        None, None, (8.0, 12.0), controls(None), || false,
    ).is_err());
    assert!(loaded(0.5).driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], ScatteringObservation::Outgoing(&[1.0]), Direction::default(),
        None, Some(&[]), (8.0, 12.0), controls(None), || false,
    ).is_err());
}

#[test]
fn total_query_budget_and_late_cancellation_never_publish_partial_sensitivities() {
    let model = HarmonicModel::from_ode(&absorber(), budget(), || false).unwrap();
    let mut capped = controls(Some(1e-6));
    capped.maximum_evaluations = 5;
    assert!(matches!(model.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
        (9.0, 11.0), capped, || false,
    ), Err(DrivenResonanceError::Invalid(_))));
    capped.maximum_evaluations = 8;
    assert!(matches!(model.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
        (9.0, 11.0), capped, || false,
    ), Err(DrivenResonanceError::EvaluationLimit { .. })));
    let mut polls = 0;
    let completed = model.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
        (9.0, 11.0), controls(Some(1e-6)), || { polls += 1; false },
    ).unwrap();
    for stop in [1, 120, polls] {
        let mut checks = 0;
        assert!(matches!(model.driven_antiresonance_with_parameter_derivative(
            &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
            (9.0, 11.0), controls(Some(1e-6)), || { checks += 1; checks == stop },
        ), Err(DrivenResonanceError::Harmonic(HarmonicError::Cancelled))));
    }
    let replay = model.driven_antiresonance_with_parameter_derivative(
        &[C64::ONE], HarmonicObservation::State(&OBSERVATION), Direction::default(), None,
        (9.0, 11.0), controls(Some(1e-6)), || false,
    ).unwrap();
    assert_eq!(completed, replay);
}

#[test]
fn finite_depth_wave_notches_retain_mixed_curvature_and_material_tuning() {
    // A finite reflection minimum needs both |H_omega|^2 and
    // Re(conj(H) H_omegaomega); a simple-zero-only formula is insufficient.
    let (m, k, z, rho) = (1.0_f64, 100.0_f64, 2.0_f64, 0.4_f64);
    let input = C64::new(0.6, -0.8);
    for c in [0.2_f64, 1.3_f64] {
        let model = loaded(c);
        let r = (1.0 - z * c) / (1.0 + z * c);
        let sum = c + 1.0 / z;
        let difference = c - 1.0 / z;
        let x = ((rho * sum * sum - difference * difference) / (1.0 - rho)).sqrt();
        let root = (x * x + 4.0 * m * k).sqrt();
        let lower = (root - x) / (2.0 * m);
        let upper = (root + x) / (2.0 * m);
        assert!(r * r < rho);
        for (dm, dk, dc, dz, alpha) in [
            (1.0, 0.0, 0.0, 0.0, C64::ZERO),
            (0.0, 1.0, 0.0, 0.0, C64::ZERO),
            (0.0, 0.0, 1.0, 0.0, C64::ZERO),
            (0.0, 0.0, 0.0, 1.0, C64::ZERO),
            (0.0, 0.0, 0.0, 0.0, C64::new(0.3, -0.4)),
            (0.2, 0.7, -0.1, 0.3, C64::new(-0.2, 0.5)),
        ] {
            let dq = [dk, 0.0, 0.0, -dm / (m * m)];
            let dr = [0.0, 0.0, 0.0, dc];
            let result = model.driven_antiresonance_with_parameter_derivative(
                &[input], ScatteringObservation::Outgoing(&[1.0]),
                Direction { storage: Some(&dq), dissipation: Some(&dr), ..Default::default() },
                Some(&[alpha * input]), Some(&[dz]), (8.0, 12.0),
                controls(Some(rho)), || false,
            ).unwrap();
            // Solve the independent original scattering rational equation:
            // rho = ((c-1/Z)^2 + X^2) / ((c+1/Z)^2 + X^2).
            // The fixed physical recovery level means rho changes when |a| does.
            let d_rho = -2.0 * rho * alpha.re;
            let d_sum = dc - dz / (z * z);
            let d_difference = dc + dz / (z * z);
            let dx = (2.0 * rho * sum * d_sum
                - 2.0 * difference * d_difference + d_rho * (sum * sum + x * x))
                / (2.0 * x * (1.0 - rho));
            let droot = (x * dx + 2.0 * (k * dm + m * dk)) / root;
            let dl = (droot - dx) / (2.0 * m) - lower * dm / m;
            let du = (droot + dx) / (2.0 * m) - upper * dm / m;
            near(result.antiresonance.minimum.angular_frequency_rad_s, 10.0, 3e-11);
            near(result.minimum_frequency_derivative_rad_s, 5.0 * (dk / k - dm / m), 3e-9);
            complex_near(result.antiresonance.observation, input.scale(r), 2e-10);
            let d_reflection = -2.0 * (z * dc + c * dz) / (1.0 + z * c).powi(2);
            complex_near(result.observation_derivative,
                input * (alpha.scale(r) + C64::from_re(d_reflection)), 3e-9);
            let curvature = 32.0 * m * m * z.powi(3) * c / (1.0 + z * c).powi(4);
            near(result.scaled_power_curvature_s2 * result.curvature_amplitude_scale.powi(2),
                curvature, 3e-8);
            let band = result.antiresonance.recovery.unwrap();
            near(band.lower.angular_frequency_rad_s, lower, 3e-10);
            near(band.upper.angular_frequency_rad_s, upper, 3e-10);
            let derivative = result.recovery_derivative.unwrap();
            near(derivative.lower_frequency_derivative_rad_s, dl, 3e-9);
            near(derivative.upper_frequency_derivative_rad_s, du, 3e-9);
            near(derivative.bandwidth_derivative_rad_s, du - dl, 5e-9);
            assert!(result.maximum_relative_residual <= 2e-12);
        }
    }
}
