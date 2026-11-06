use super::{controls, derivative_controls, near, orbit, scales};
use fs_couple::phs_orbit::{
    AutonomousPortFeedback,
    sideband::{
        PhsPortSidebandError, PhsPortSidebandSource,
        solve_autonomous_feedback_port_sideband_response,
        solve_autonomous_feedback_port_sideband_response_analytic,
    },
};
use fs_math::{c64::C64, det};
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_orbit::{
    HbOrbit,
    sideband::{SidebandError, SidebandScales, refinement::SidebandRefinementControls},
};
use fs_phs::{PortHamiltonian, RelaxationBranch, Storage};

// An authored active source maintains the exact mechanical limit cycle
// r'=gain*(1-r^2)*r, theta'=-omega while preserving the original H/R/G and
// Maxwell state. Radial decay is -2*gain, phase is neutral, memory decays -1/tau.
struct CycleDrive {
    omega: f64,
    gain: f64,
    resistance: f64,
    cubic: f64,
    memory: f64,
}
impl AutonomousPortFeedback for CycleDrive {
    fn port_count(&self) -> usize {
        2
    }
    fn input(&self, x: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        let radial = self.gain * (1.0 - x[0] * x[0] - x[1] * x[1]);
        out[0] = self.resistance * flows[0] + radial * x[0];
        out[1] = self.omega * (flows[0] - x[0]) + self.resistance * flows[1] + radial * x[1];
        true
    }
}
impl CycleDrive {
    fn total_derivative(&self, x: &[f64], _: &[f64], out: &mut [f64]) -> bool {
        let n = x.len();
        let radial = self.gain * (1.0 - x[0] * x[0] - x[1] * x[1]);
        let (rq, rp) = (-2.0 * self.gain * x[0], -2.0 * self.gain * x[1]);
        let hqq = 1.0 + 3.0 * self.cubic * x[0] * x[0] + self.memory;
        out.fill(0.0);
        out[0] = self.resistance * hqq + radial + x[0] * rq;
        out[1] = x[0] * rp;
        out[n] = self.omega * (hqq - 1.0) + x[1] * rq;
        out[n + 1] = self.resistance + radial + x[1] * rp;
        if n == 3 {
            // These are the original effort/memory cross derivatives, not a
            // feedback derivative at fixed physical flow.
            out[2] = -self.resistance * det::sqrt(self.memory);
            out[n + 2] = -self.omega * det::sqrt(self.memory);
        }
        true
    }
}

fn system(source: &CycleDrive, tau: f64, storage: Option<Box<dyn Storage>>) -> PortHamiltonian {
    let storage = storage.unwrap_or_else(|| {
        Box::new(SosModalStorage {
            omegas: vec![1.0],
            channels: vec![StressChannel {
                coefficient: source.cubic,
                coupling: vec![1.0],
            }],
        })
    });
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, source.omega, -source.omega, 0.0],
        vec![source.resistance, 0.0, 0.0, source.resistance],
        vec![1.0, 0.0, 0.0, 1.0],
        storage,
    )
    .unwrap();
    if source.memory == 0.0 {
        system
    } else {
        system
            .with_relaxation_branches(vec![RelaxationBranch {
                projection: vec![1.0, 0.0],
                stiffness: source.memory,
                relaxation_time_s: tau,
            }])
            .unwrap()
    }
}
fn candidate(source: &CycleDrive, tau: f64, phase: f64) -> HbOrbit {
    let q = C64::new(det::cos(phase), det::sin(phase)).scale(0.5);
    let mut coefficients = vec![vec![C64::ZERO, q], vec![C64::ZERO, C64::new(0.0, 1.0) * q]];
    if source.memory != 0.0 {
        coefficients.push(vec![
            C64::ZERO,
            q.scale(det::sqrt(source.memory)) / C64::new(1.0, tau * source.omega),
        ]);
    }
    orbit(source.omega, coefficients)
}

#[test]
fn fixed_base_autonomous_stability_keeps_memory_phase_and_extra_neutral_family() {
    use fs_spectral::floquet::FloquetClassification;

    // These all have the same exact finite-amplitude original periodic base.
    // Its radial dynamics may decay, grow, or retain an amplitude family. The
    // original Maxwell memory stays in both the phase flow and complete map.
    for (gain, expected_class) in [
        (0.2, FloquetClassification::Stable),
        (-0.1, FloquetClassification::Unstable),
        (0.0, FloquetClassification::Unresolved),
    ] {
        let drive = CycleDrive {
            omega: 1.2,
            gain,
            resistance: 0.2,
            cubic: 0.6,
            memory: 0.3,
        };
        let tau = 0.7;
        let system = system(&drive, tau, None);
        let base = candidate(&drive, tau, 0.0);
        let original = base.coeffs.clone();
        let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| drive.total_derivative(x, y, out);
        let physical = controls(3, 2, 2, 4);
        let mut source = PhsPortSidebandSource::autonomous_feedback(
            &system,
            &drive,
            &derivative,
            &base,
            &physical,
        )
        .unwrap();
        let mut prepared = source.prepare(&scales(3), || false).unwrap();
        let settings = super::base_stability_controls(&[0.5, 1.5, 2.5]);
        let report = prepared
            .analyze_base_stability(&settings, || false)
            .unwrap();
        assert!(report.autonomous);
        assert_eq!(report.classification, expected_class);
        assert_eq!(report.state_scales, settings.variational.state_scales);
        let phase_flow = report.scaled_phase_flow.as_ref().unwrap();
        assert_eq!(phase_flow.len(), 3);
        near(phase_flow[0], C64::ZERO, 2e-14);
        near(phase_flow[1], C64::from_re(-drive.omega / 1.5), 2e-14);
        let memory_flow = (det::sqrt(drive.memory) - 2.0 * base.coeffs[2][1].re) / tau;
        near(phase_flow[2], C64::from_re(memory_flow / 2.5), 2e-14);
        assert!(phase_flow[2].abs() > 0.05);
        let radial = det::exp(-2.0 * gain * report.period_s);
        let memory = det::exp(-report.period_s / tau);
        for (index, level) in report.levels.iter().enumerate() {
            assert_eq!(level.monodromy.len(), 9);
            assert_eq!(level.spectrum.spectrum.values.len(), 3);
            assert_eq!(level.spectrum.spectrum.right_vectors.len(), 3);
            assert!(level.maximum_scaled_base_excursion > 3.0);
            let accuracy = if index == 2 { 3e-8 } else { 3e-6 };
            for expected in [1.0, radial, memory] {
                assert!(
                    level
                        .spectrum
                        .spectrum
                        .values
                        .iter()
                        .any(|&value| (value - C64::from_re(expected)).abs() < accuracy)
                );
            }
            assert!(
                level
                    .spectrum
                    .spectrum
                    .relative_residuals
                    .iter()
                    .all(|&value| value <= settings.floquet.spectrum.maximum_relative_residual)
            );
            if gain != 0.0 {
                let phase_index = level.spectrum.phase_multiplier_index.unwrap();
                assert!(
                    (level.spectrum.spectrum.values[phase_index] - C64::ONE).abs()
                        < settings.floquet.phase_multiplier_tolerance
                );
                assert!(
                    level.spectrum.phase_residual.unwrap()
                        <= settings.floquet.maximum_phase_residual
                );
                assert!(
                    (level.spectrum.effective_spectral_radius.unwrap() - radial.max(memory)).abs()
                        < accuracy
                );
            } else {
                assert_eq!(
                    level.spectrum.classification,
                    FloquetClassification::Unresolved
                );
                assert!(level.spectrum.phase_multiplier_index.is_none());
                assert!(level.spectrum.effective_spectral_radius.is_none());
                assert!(level.spectrum.dominant_growth_rate_per_s.is_none());
                assert_eq!(
                    level
                        .spectrum
                        .spectrum
                        .values
                        .iter()
                        .filter(|&&value| (value - C64::ONE).abs()
                            < settings.floquet.phase_multiplier_tolerance)
                        .count(),
                    2
                );
            }
        }
        if gain != 0.0 {
            assert!((report.dominant_growth_rate_per_s.unwrap() + 2.0 * gain).abs() < 3e-8);
        }
        assert_eq!(base.coeffs, original);
    }
}

// Independent rotating-coordinate solution. Q(theta)^T turns the original
// perturbation into a scalar radial filter (s+2*gain)^-1 and a neutral phase
// integrator s^-1. Rotate back to the SAME Cartesian state and phase clock.
fn exact_mechanical(
    source: &CycleDrive,
    phase: f64,
    carrier: f64,
    force: [C64; 2],
) -> Vec<[C64; 2]> {
    let mut answer = vec![[C64::ZERO; 2]; 5];
    for sigma in [-1, 1] {
        let angle = sigma as f64 * phase;
        let clock = C64::new(det::cos(angle), det::sin(angle));
        let i = C64::new(0.0, 0.5 * sigma as f64);
        let rate = C64::new(0.0, carrier + sigma as f64 * source.omega);
        let radial =
            (force[0].scale(0.5) + i * force[1]) * clock / (rate + C64::from_re(2.0 * source.gain));
        let tangent = (-i * force[0] + force[1].scale(0.5)) * clock / rate;
        for tau in [-1, 1] {
            let angle = tau as f64 * phase;
            let clock = C64::new(det::cos(angle), det::sin(angle));
            let i = C64::new(0.0, 0.5 * tau as f64);
            let band = (sigma + tau + 2) as usize;
            answer[band][0] = answer[band][0] + clock * (radial.scale(0.5) - i * tangent);
            answer[band][1] = answer[band][1] + clock * (i * radial + tangent.scale(0.5));
        }
    }
    answer
}

// Differentiate the independent rotating-coordinate filters analytically.
// The original base phase/rotation is fixed; only the signed carrier changes.
fn exact_mechanical_derivative(
    source: &CycleDrive,
    phase: f64,
    carrier: f64,
    force: [C64; 2],
) -> Vec<[C64; 2]> {
    let mut answer = vec![[C64::ZERO; 2]; 5];
    for sigma in [-1, 1] {
        let angle = sigma as f64 * phase;
        let clock = C64::new(det::cos(angle), det::sin(angle));
        let i = C64::new(0.0, 0.5 * sigma as f64);
        let rate = C64::new(0.0, carrier + sigma as f64 * source.omega);
        let radial_rate = rate + C64::from_re(2.0 * source.gain);
        let radial = C64::new(0.0, -1.0) * (force[0].scale(0.5) + i * force[1]) * clock
            / (radial_rate * radial_rate);
        let tangent =
            C64::new(0.0, -1.0) * (-i * force[0] + force[1].scale(0.5)) * clock / (rate * rate);
        for rotate_sign in [-1, 1] {
            let angle = rotate_sign as f64 * phase;
            let clock = C64::new(det::cos(angle), det::sin(angle));
            let i = C64::new(0.0, 0.5 * rotate_sign as f64);
            let band = (sigma + rotate_sign + 2) as usize;
            answer[band][0] = answer[band][0] + clock * (radial.scale(0.5) - i * tangent);
            answer[band][1] = answer[band][1] + clock * (i * radial + tangent.scale(0.5));
        }
    }
    answer
}

#[test]
fn autonomous_probe_matches_radial_phase_solution_with_nonlinear_output_and_maxwell_memory() {
    let source = CycleDrive {
        omega: 1.2,
        gain: 0.45,
        resistance: 0.2,
        cubic: 0.6,
        memory: 0.3,
    };
    let (tau, phase, carrier) = (0.7, 0.31, 0.47);
    let system = system(&source, tau, None);
    let base = candidate(&source, tau, phase);
    let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| source.total_derivative(x, y, out);
    let force = [C64::new(0.12, -0.04), C64::new(-0.03, 0.05)];
    let mut forcing = vec![C64::ZERO; 10];
    forcing[4..6].copy_from_slice(&force);
    let result = solve_autonomous_feedback_port_sideband_response(
        &system,
        &source,
        &derivative,
        &base,
        carrier,
        &forcing,
        &controls(3, 2, 2, 4),
        &scales(3),
        || false,
    )
    .unwrap();
    let mechanics = exact_mechanical(&source, phase, carrier, force);
    let mut expected = vec![[C64::ZERO; 3]; 5];
    for band in 0..5 {
        let w = carrier + (band as f64 - 2.0) * source.omega;
        expected[band] = [
            mechanics[band][0],
            mechanics[band][1],
            mechanics[band][0].scale(det::sqrt(source.memory)) / C64::new(1.0, tau * w),
        ];
        for j in 0..3 {
            near(
                result.state_response.coefficients[band * 3 + j],
                expected[band][j],
                4e-11,
            );
        }
    }
    let state = |k: isize, j: usize| {
        if k.unsigned_abs() > 2 {
            C64::ZERO
        } else {
            expected[(k + 2) as usize][j]
        }
    };
    let modulation = (base.coeffs[0][1] * base.coeffs[0][1]).scale(3.0 * source.cubic);
    for k in -4..=4 {
        let effort = state(k, 0).scale(1.0 + source.memory + 1.5 * source.cubic)
            + modulation * state(k - 2, 0)
            + modulation.conj() * state(k + 2, 0)
            - state(k, 2).scale(det::sqrt(source.memory));
        near(
            result.port_coefficients[(k + 4) as usize * 2],
            effort,
            4e-11,
        );
        near(
            result.port_coefficients[(k + 4) as usize * 2 + 1],
            state(k, 1),
            4e-11,
        );
    }
    assert_eq!(result.state_dimension, 3);
    assert!(result.state_response.coefficients[0].abs() > 1e-3);
    assert!(result.state_response.coefficients[12].abs() > 1e-3);
    assert!(result.base_residual < 1e-11);
    assert!(result.state_response.omitted_residual < 1e-11);
}

struct UndeclaredQuadratic;
impl Storage for UndeclaredQuadratic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (x[0] * x[0] + x[1] * x[1])
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(x);
    }
}

#[test]
fn autonomous_probe_requires_complete_storage_and_total_feedback_derivatives() {
    let source = CycleDrive {
        omega: 1.0,
        gain: 0.4,
        resistance: 0.2,
        cubic: 0.0,
        memory: 0.0,
    };
    let base = candidate(&source, 1.0, 0.0);
    let system = system(&source, 1.0, Some(Box::new(UndeclaredQuadratic)));
    let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| source.total_derivative(x, y, out);
    let mut forcing = vec![C64::ZERO; 10];
    forcing[4] = C64::from_re(0.2);
    let c = controls(2, 2, 2, 2);
    assert!(matches!(
        solve_autonomous_feedback_port_sideband_response(
            &system,
            &source,
            &derivative,
            &base,
            0.3,
            &forcing,
            &c,
            &scales(2),
            || false,
        ),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::TangentUnavailable
        ))
    ));
    let hessian = |_: &[f64], v: &[f64], out: &mut [f64]| {
        out.copy_from_slice(v);
        true
    };
    let result = solve_autonomous_feedback_port_sideband_response_analytic(
        &system,
        &source,
        &derivative,
        &base,
        0.3,
        &forcing,
        &c,
        &scales(2),
        &hessian,
        || false,
    )
    .unwrap();
    let expected = exact_mechanical(&source, 0.0, 0.3, [C64::from_re(0.2), C64::ZERO]);
    for (actual, expected) in result
        .state_response
        .coefficients
        .iter()
        .zip(expected.iter().flatten())
    {
        near(*actual, *expected, 3e-11);
    }
    for incomplete in [false, true] {
        let refused = |_: &[f64], _: &[f64], out: &mut [f64]| {
            if incomplete {
                out[0] = 0.0;
            }
            incomplete
        };
        assert!(matches!(
            solve_autonomous_feedback_port_sideband_response_analytic(
                &system,
                &source,
                &refused,
                &base,
                0.3,
                &forcing,
                &c,
                &scales(2),
                &hessian,
                || false,
            ),
            Err(PhsPortSidebandError::Sideband(
                SidebandError::TangentRefused { .. }
            ))
        ));
    }
}

#[test]
fn autonomous_neutral_carrier_is_not_phase_pinned_or_given_artificial_damping() {
    let source = CycleDrive {
        omega: 1.0,
        gain: 0.0,
        resistance: 0.2,
        cubic: 0.0,
        memory: 0.0,
    };
    let base = candidate(&source, 1.0, 0.0);
    let system = system(&source, 1.0, None);
    let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| source.total_derivative(x, y, out);
    let mut forcing = vec![C64::ZERO; 6];
    forcing[2] = C64::ONE;
    let scales = SidebandScales {
        states: vec![1.0; 2],
        equations: vec![1.0; 2],
    };
    assert!(matches!(
        solve_autonomous_feedback_port_sideband_response(
            &system,
            &source,
            &derivative,
            &base,
            0.0,
            &forcing,
            &controls(2, 2, 1, 1),
            &scales,
            || false,
        ),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::Singular { .. }
        ))
    ));
}

#[test]
fn prepared_autonomous_sweep_reuses_one_complete_tangent_and_keeps_physical_memory_outputs() {
    let drive = CycleDrive {
        omega: 1.2,
        gain: 0.45,
        resistance: 0.2,
        cubic: 0.6,
        memory: 0.3,
    };
    let (tau, phase) = (0.7, 0.31);
    let system = system(&drive, tau, None);
    let base = candidate(&drive, tau, phase);
    let original = base.coeffs.clone();
    let calls = core::cell::Cell::new(0_usize);
    let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        drive.total_derivative(x, y, out)
    };
    let c = controls(3, 2, 2, 4);
    let mut source =
        PhsPortSidebandSource::autonomous_feedback(&system, &drive, &derivative, &base, &c)
            .unwrap();
    let mut prepared = source.prepare(&scales(3), || false).unwrap();
    let preparations = 3 * c.sideband.samples_per_period + 1;
    assert_eq!(calls.get(), preparations);
    assert!(prepared.base_residual() < 1e-11);
    assert!(prepared.tangent_sampling_defect() < 1e-11);
    let forces = [
        [C64::new(0.12, -0.04), C64::new(-0.03, 0.05)],
        [C64::new(-0.07, 0.09), C64::new(0.11, -0.02)],
    ];
    for carrier in [0.47, 1.93, -0.37] {
        let mut factor = prepared.factor(carrier, || false).unwrap();
        assert_eq!(factor.carrier_angular_frequency(), carrier);
        for force in forces {
            let mut forcing = vec![C64::ZERO; 10];
            forcing[4..6].copy_from_slice(&force);
            let result = factor.solve(&forcing, || false).unwrap();
            // Factorization and physical observation must not resample the
            // feedback tangent, despite changing carrier and complex port drive.
            assert_eq!(calls.get(), preparations);
            let mechanics = exact_mechanical(&drive, phase, carrier, force);
            let expected: Vec<_> = mechanics
                .iter()
                .enumerate()
                .map(|(band, response)| {
                    let w = carrier + (band as f64 - 2.0) * drive.omega;
                    [
                        response[0],
                        response[1],
                        response[0].scale(det::sqrt(drive.memory)) / C64::new(1.0, tau * w),
                    ]
                })
                .collect();
            for (actual, expected) in result
                .state_response
                .coefficients
                .iter()
                .zip(expected.iter().flatten())
            {
                near(*actual, *expected, 5e-11);
            }
            let state = |k: isize, j: usize| {
                if k.unsigned_abs() > 2 {
                    C64::ZERO
                } else {
                    expected[(k + 2) as usize][j]
                }
            };
            let modulation = (base.coeffs[0][1] * base.coeffs[0][1]).scale(3.0 * drive.cubic);
            for k in -4..=4 {
                let effort = state(k, 0).scale(1.0 + drive.memory + 1.5 * drive.cubic)
                    + modulation * state(k - 2, 0)
                    + modulation.conj() * state(k + 2, 0)
                    - state(k, 2).scale(det::sqrt(drive.memory));
                near(
                    result.port_coefficients[(k + 4) as usize * 2],
                    effort,
                    5e-11,
                );
                near(
                    result.port_coefficients[(k + 4) as usize * 2 + 1],
                    state(k, 1),
                    5e-11,
                );
            }
            assert_eq!(result.observation_trace.len(), 3);
            assert_eq!(result.hessian_actions, 2 * (64 + 128 + 256));
        }
        assert!(matches!(
            factor.solve(&[C64::ZERO], || false),
            Err(PhsPortSidebandError::Invalid(_))
        ));
        assert!(matches!(
            factor.solve(&[C64::ZERO; 10], || true),
            Err(PhsPortSidebandError::Cancelled)
        ));
        // A refused/cancelled solve leaves the same cached factor usable.
        let zero = factor.solve(&[C64::ZERO; 10], || false).unwrap();
        assert!(
            zero.state_response
                .coefficients
                .iter()
                .all(|z| *z == C64::ZERO)
        );
        assert!(zero.port_coefficients.iter().all(|z| *z == C64::ZERO));
    }
    assert_eq!(calls.get(), preparations);
    assert_eq!(base.coeffs, original);
}

#[test]
fn adaptive_autonomous_probe_resolves_generated_bands_before_original_memory_and_port_observation()
{
    let drive = CycleDrive {
        omega: 1.2,
        gain: 0.45,
        resistance: 0.2,
        cubic: 0.6,
        memory: 0.3,
    };
    let (tau, phase, carrier) = (0.7, 0.31, 0.47);
    let system = system(&drive, tau, None);
    let base = candidate(&drive, tau, phase);
    let original = base.coeffs.clone();
    let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| drive.total_derivative(x, y, out);
    let mut c = controls(3, 2, 1, 4);
    c.sideband.maximum_unknowns = 3 * 15;
    let refinement = SidebandRefinementControls {
        maximum_harmonics: 7,
        equation_tolerances: vec![2e-10; 3],
        state_tolerances: vec![2e-10; 3],
    };
    let force = [C64::new(0.12, -0.04), C64::new(-0.03, 0.05)];
    let mut forcing = vec![C64::ZERO; 6];
    forcing[2..4].copy_from_slice(&force);
    let mut source =
        PhsPortSidebandSource::autonomous_feedback(&system, &drive, &derivative, &base, &c)
            .unwrap();
    let adaptive = source
        .solve_adaptive(carrier, &forcing, &scales(3), &refinement, || false)
        .unwrap();
    assert_eq!(
        adaptive
            .refinement_trace
            .iter()
            .map(|check| check.harmonics)
            .collect::<Vec<_>>(),
        [1, 3, 7]
    );
    assert!(adaptive.refinement_trace[0].equation_residual_ratio > 1.0);
    assert!(adaptive.refinement_trace[1].state_change_ratio.unwrap() > 1.0);
    let result = adaptive.response;
    assert_eq!(result.state_harmonics, 7);
    assert_eq!(result.output_harmonics, 4);
    assert_eq!(result.state_dimension, 3);
    assert_eq!(result.port_count, 2);
    let mechanics = exact_mechanical(&drive, phase, carrier, force);
    let state = |k: isize, j: usize| {
        if k.unsigned_abs() > 2 {
            return C64::ZERO;
        }
        let mechanical = mechanics[(k + 2) as usize];
        if j < 2 {
            mechanical[j]
        } else {
            mechanical[0].scale(det::sqrt(drive.memory))
                / C64::new(1.0, tau * (carrier + k as f64 * drive.omega))
        }
    };
    for k in -7..=7 {
        for j in 0..3 {
            near(
                result.state_response.coefficients[(k + 7) as usize * 3 + j],
                state(k, j),
                5e-11,
            );
        }
    }
    let modulation = (base.coeffs[0][1] * base.coeffs[0][1]).scale(3.0 * drive.cubic);
    for k in -4..=4 {
        let effort = state(k, 0).scale(1.0 + drive.memory + 1.5 * drive.cubic)
            + modulation * state(k - 2, 0)
            + modulation.conj() * state(k + 2, 0)
            - state(k, 2).scale(det::sqrt(drive.memory));
        near(
            result.port_coefficients[(k + 4) as usize * 2],
            effort,
            5e-11,
        );
        near(
            result.port_coefficients[(k + 4) as usize * 2 + 1],
            state(k, 1),
            5e-11,
        );
    }
    assert!(result.state_response.coefficients[(7 - 2) * 3].abs() > 1e-3);
    assert!(result.state_response.coefficients[(7 + 2) * 3].abs() > 1e-3);
    assert!(result.port_coefficients[0].abs() > 1e-4);
    assert!(result.port_coefficients[16].abs() > 1e-4);
    assert_eq!(base.coeffs, original);
}

#[test]
fn autonomous_carrier_derivative_matches_rotating_filters_original_memory_and_nonlinear_outputs() {
    let drive = CycleDrive {
        omega: 1.2,
        gain: 0.45,
        resistance: 0.2,
        cubic: 0.6,
        memory: 0.3,
    };
    let (tau, phase, carrier) = (0.7, 0.31, 0.47);
    let system = system(&drive, tau, None);
    let base = candidate(&drive, tau, phase);
    let original = base.coeffs.clone();
    let calls = core::cell::Cell::new(0_usize);
    let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        drive.total_derivative(x, y, out)
    };
    let mut c = controls(3, 2, 2, 4);
    c.max_hessian_actions = 1792;
    let force = [C64::new(0.12, -0.04), C64::new(-0.03, 0.05)];
    let mut forcing = vec![C64::ZERO; 10];
    forcing[4..6].copy_from_slice(&force);
    let mut source =
        PhsPortSidebandSource::autonomous_feedback(&system, &drive, &derivative, &base, &c)
            .unwrap();
    let mut prepared = source.prepare(&scales(3), || false).unwrap();
    assert_eq!(calls.get(), 385);
    let mut factor = prepared.factor(carrier, || false).unwrap();
    let pair = factor
        .solve_with_carrier_derivative(&forcing, None, derivative_controls(), &[2e-10; 2], || false)
        .unwrap();
    assert_eq!(calls.get(), 385);
    assert_eq!(pair.response.hessian_actions, 896);
    assert_eq!(pair.total_hessian_actions, 1792);
    let mechanics = exact_mechanical(&drive, phase, carrier, force);
    let mechanical_derivative = exact_mechanical_derivative(&drive, phase, carrier, force);
    let state = |k: isize, j: usize, tangent: bool| {
        if k.unsigned_abs() > 2 {
            return C64::ZERO;
        }
        let band = (k + 2) as usize;
        if j < 2 {
            if tangent {
                mechanical_derivative[band][j]
            } else {
                mechanics[band][j]
            }
        } else {
            let relaxation = C64::new(1.0, tau * (carrier + k as f64 * drive.omega));
            let memory = mechanics[band][0].scale(det::sqrt(drive.memory)) / relaxation;
            if tangent {
                mechanical_derivative[band][0].scale(det::sqrt(drive.memory)) / relaxation
                    - C64::new(0.0, tau) * memory / relaxation
            } else {
                memory
            }
        }
    };
    for k in -2..=2 {
        for j in 0..3 {
            near(
                pair.response.state_response.coefficients[(k + 2) as usize * 3 + j],
                state(k, j, false),
                5e-11,
            );
            near(
                pair.state_angular_derivative.coefficients[(k + 2) as usize * 3 + j],
                state(k, j, true),
                6e-11,
            );
        }
    }
    let modulation = (base.coeffs[0][1] * base.coeffs[0][1]).scale(3.0 * drive.cubic);
    for k in -4..=4 {
        let effort = state(k, 0, true).scale(1.0 + drive.memory + 1.5 * drive.cubic)
            + modulation * state(k - 2, 0, true)
            + modulation.conj() * state(k + 2, 0, true)
            - state(k, 2, true).scale(det::sqrt(drive.memory));
        near(
            pair.port_angular_derivative[(k + 4) as usize * 2],
            effort,
            6e-11,
        );
        near(
            pair.port_angular_derivative[(k + 4) as usize * 2 + 1],
            state(k, 1, true),
            6e-11,
        );
    }
    assert!(pair.port_angular_derivative[0].abs() > 1e-6);
    assert!(pair.port_angular_derivative[16].abs() > 1e-6);
    assert_eq!(base.coeffs, original);
}

#[test]
fn joint_adaptive_autonomous_pair_resolves_complete_rotating_and_memory_fields_before_ports() {
    use fs_orbit::sideband::sensitivity::refinement::{
        SidebandCarrierRefinementControls, SidebandCarrierRefinementTolerances,
    };
    let drive = CycleDrive {
        omega: 1.2,
        gain: 0.45,
        resistance: 0.2,
        cubic: 0.6,
        memory: 0.3,
    };
    let (tau, phase, carrier) = (0.7, 0.31, 0.47);
    let system = system(&drive, tau, None);
    let base = candidate(&drive, tau, phase);
    let derivative = |x: &[f64], y: &[f64], out: &mut [f64]| drive.total_derivative(x, y, out);
    let mut c = controls(3, 2, 1, 4);
    c.sideband.maximum_unknowns = 45;
    c.max_hessian_actions = 1792;
    let refinement = SidebandCarrierRefinementControls {
        maximum_harmonics: 7,
        tolerances: SidebandCarrierRefinementTolerances {
            equation_tolerances: vec![2e-9; 3],
            state_tolerances: vec![2e-9; 3],
            derivative_equation_tolerances: vec![2e-9; 3],
            derivative_state_tolerances: vec![2e-9; 3],
        },
    };
    let force = [C64::new(0.08, -0.03), C64::new(-0.04, 0.06)];
    let forcing = [
        C64::ZERO,
        C64::ZERO,
        force[0],
        force[1],
        C64::ZERO,
        C64::ZERO,
    ];
    let mut source =
        PhsPortSidebandSource::autonomous_feedback(&system, &drive, &derivative, &base, &c)
            .unwrap();
    let resolved = source
        .solve_adaptive_with_carrier_derivative(
            carrier,
            &forcing,
            None,
            &scales(3),
            derivative_controls(),
            &[2e-10; 2],
            &refinement,
            || false,
        )
        .unwrap();
    assert_eq!(
        resolved
            .refinement_trace
            .iter()
            .map(|c| c.response.harmonics)
            .collect::<Vec<_>>(),
        [1, 3, 7]
    );
    assert!(!resolved.refinement_trace[0].equations_resolved());
    assert!(
        resolved.refinement_trace[1]
            .angular_derivative
            .state_change_ratio
            .unwrap()
            > 1.0
    );
    assert!(resolved.refinement_trace[2].is_resolved());
    let pair = resolved.response;
    assert_eq!(pair.response.state_harmonics, 7);
    assert_eq!(pair.response.output_harmonics, 4);
    assert_eq!(pair.total_hessian_actions, 1792);
    let mechanics = exact_mechanical(&drive, phase, carrier, force);
    let slope = exact_mechanical_derivative(&drive, phase, carrier, force);
    let exact = |k: isize, j: usize, derivative: bool| {
        if k.unsigned_abs() > 2 {
            return C64::ZERO;
        }
        let band = (k + 2) as usize;
        if j < 2 {
            return if derivative {
                slope[band][j]
            } else {
                mechanics[band][j]
            };
        }
        let r = C64::new(1.0, tau * (carrier + k as f64 * drive.omega));
        let memory = mechanics[band][0].scale(det::sqrt(drive.memory)) / r;
        if derivative {
            slope[band][0].scale(det::sqrt(drive.memory)) / r - C64::new(0.0, tau) * memory / r
        } else {
            memory
        }
    };
    for k in -7..=7 {
        for j in 0..3 {
            near(
                pair.response.state_response.coefficients[(k + 7) as usize * 3 + j],
                exact(k, j, false),
                8e-11,
            );
            near(
                pair.state_angular_derivative.coefficients[(k + 7) as usize * 3 + j],
                exact(k, j, true),
                1e-10,
            );
        }
    }
    let modulation = (base.coeffs[0][1] * base.coeffs[0][1]).scale(3.0 * drive.cubic);
    for (field, differentiated) in [
        (&pair.response.port_coefficients, false),
        (&pair.port_angular_derivative, true),
    ] {
        for k in -4..=4 {
            let effort = exact(k, 0, differentiated).scale(1.0 + drive.memory + 1.5 * drive.cubic)
                + modulation * exact(k - 2, 0, differentiated)
                + modulation.conj() * exact(k + 2, 0, differentiated)
                - exact(k, 2, differentiated).scale(det::sqrt(drive.memory));
            near(field[(k + 4) as usize * 2], effort, 1e-10);
            near(
                field[(k + 4) as usize * 2 + 1],
                exact(k, 1, differentiated),
                1e-10,
            );
        }
    }
}
