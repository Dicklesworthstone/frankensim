//! Full original forced response, physical phase and local-branch peak queries.

use fs_couple::phs_orbit::{
    PeriodicPortInput,
    forced_resonance::{
        PhsForcedResonanceControls, PhsForcedResonanceError, survey_forced_port_resonances,
    },
    sideband::PhsPortSidebandControls,
    solve_forced_harmonic_balance_scaled_seeded,
    spectrum::{PortSpectrumControls, observe_port_spectrum},
};
use fs_math::{c64::C64, det};
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbScales,
    sideband::{SidebandControls, SidebandScales},
};
use fs_phs::{
    PortHamiltonian, RelaxationBranch,
    harmonic::resonance::{ResonanceSearchError, ResonanceSurveyControls, ResonanceSurveyError},
};

fn controls(n: usize, m: usize, h: usize, output: usize) -> PhsForcedResonanceControls {
    PhsForcedResonanceControls {
        harmonic_balance: HbBudget {
            harmonics: h,
            max_newton: 40,
            tolerance: 2e-12,
            deflate_equilibrium: false,
        },
        harmonic_balance_scales: HbScales {
            states: vec![1.0; n],
            equations: vec![1.0; n],
            angular_frequency: 2.0,
        },
        tangent: PhsPortSidebandControls {
            sideband: SidebandControls {
                harmonics: h,
                samples_per_period: 128,
                maximum_unknowns: n * (2 * h + 1),
                maximum_payload_bytes: 32 << 20,
                maximum_work: 1_000_000_000,
                maximum_base_residual: 2e-10,
                maximum_tangent_sampling_defect: 2e-10,
                maximum_relative_residual: 1e-12,
                maximum_omitted_residual: 2e-8,
            },
            output_harmonics: output,
            initial_samples: 64,
            max_samples: 256,
            max_hessian_actions: 1024,
            max_workspace_bytes: 8 << 20,
            max_observation_work: 100_000_000,
            envelope_tolerances: vec![2e-9; m],
        },
        tangent_scales: SidebandScales {
            states: vec![1.0; n],
            equations: vec![1.0; n],
        },
        spectrum: PortSpectrumControls {
            harmonics: output,
            initial_samples: 64,
            max_samples: 256,
            max_source_evaluations: 512,
            max_workspace_bytes: 8 << 20,
            rms_tolerances: vec![2e-9; m],
        },
        maximum_log_frequency_distance: 0.3,
        predictor_rms_tolerances: vec![1.0; n],
        survey: ResonanceSurveyControls {
            absolute_frequency_tolerance_rad_s: 2e-8,
            relative_frequency_tolerance: 1e-9,
            maximum_evaluations: 128,
            maximum_peaks: 4,
        },
    }
}

fn orbit(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs,
        residual_trace: vec![0.0],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual {actual:?}, expected {expected:?}"
    );
}

const FORCE: [C64; 4] = [
    C64::new(0.02, 0.0),
    C64::new(0.04, -0.02),
    C64::new(0.0, 0.0),
    C64::new(-0.013, 0.007),
];

// Independent original mass/spring/damper + Maxwell dynamic stiffness, in the
// positive-time Fourier chart. The memory coordinate is retained explicitly.
fn exact_memory(w: f64, force: C64) -> ([C64; 3], C64) {
    let (mass, stiffness, damping, arm, tau) = (2.0, 4.0, 0.3, 3.0, 0.4);
    let s = C64::new(0.0, w);
    let relaxation = C64::ONE + s.scale(tau);
    let dynamic =
        C64::new(stiffness - mass * w * w, damping * w) + (s.scale(tau) / relaxation).scale(arm);
    let log_dynamic = C64::new(-2.0 * mass * w * w, damping * w)
        + (s.scale(tau) / (relaxation * relaxation)).scale(arm);
    let q = force / dynamic;
    let log_q = -(q * log_dynamic / dynamic);
    (
        [q, s * q.scale(mass), q.scale(det::sqrt(arm)) / relaxation],
        s * (q + log_q),
    )
}

fn memory_fixture(omega: f64) -> (PortHamiltonian, PeriodicPortInput, HbOrbit) {
    let system = fs_phs::mass_spring_damper(2.0, 4.0, 0.3)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 3.0,
            relaxation_time_s: 0.4,
        }])
        .unwrap();
    let mut states = vec![vec![C64::ZERO; 4]; 3];
    for (k, &force) in FORCE.iter().enumerate() {
        let (state, _) = exact_memory(k as f64 * omega, force);
        for j in 0..3 {
            states[j][k] = state[j];
        }
    }
    (
        system,
        PeriodicPortInput::new(omega, vec![FORCE.to_vec()]).unwrap(),
        orbit(omega, states),
    )
}

#[test]
fn original_memory_phase_and_higher_harmonic_use_drive_frequency_derivatives() {
    for harmonic in [1, 3] {
        let factor = harmonic as f64;
        let (system, input, anchor) = memory_fixture(1.6 / factor);
        let grid: Vec<_> = (0..=7).map(|i| (1.3 + 0.1 * i as f64) / factor).collect();
        let controls = controls(3, 1, 3, 3);
        let result = survey_forced_port_resonances(
            &system,
            &input,
            &anchor,
            0,
            harmonic,
            &grid,
            &controls,
            || false,
        )
        .unwrap();
        assert_eq!(result.resonances.peaks.len(), 1);
        assert!(result.resonances.evaluations <= controls.survey.maximum_evaluations);
        let peak = &result.resonances.peaks[0];
        let omega = peak.peak.angular_frequency_rad_s;
        let physical = factor * omega;
        let (states, log_velocity) = exact_memory(physical, FORCE[harmonic]);
        let value = states[1].scale(0.5).conj().scale(2.0);
        // Divide by the drive fundamental: k=3 must retain its factor of three.
        let derivative = log_velocity.conj().scale(2.0 / omega);
        near(peak.observation, value, 3e-10);
        near(peak.observation_derivative, derivative, 3e-9);
        assert!((derivative / value).re.abs() < 3e-6);
        assert!(physical > 1.45 && physical < 1.9);
        assert!(peak.observation.im.abs() > 1e-4, "authored phase was lost");
        assert!(
            result.maximum_predictor_correction_rms[2] > 1e-5,
            "material-memory prediction was omitted from branch admission"
        );
        for side in [-1.0, 1.0] {
            let nearby = physical * (1.0 + 0.01 * side);
            let response = exact_memory(nearby, FORCE[harmonic]).0[1].abs();
            assert!(response < states[1].abs());
        }
    }
}

fn nonlinear_fixture(amplitude: f64) -> (PortHamiltonian, PeriodicPortInput, HbOrbit) {
    let (stiffness, cubic, damping, phase) = (4.0, 2.0, 0.3, 0.37);
    let omega = det::sqrt(stiffness + 0.75 * cubic * amplitude * amplitude);
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(SosModalStorage {
            omegas: vec![det::sqrt(stiffness)],
            channels: vec![StressChannel {
                coefficient: cubic,
                coupling: vec![1.0],
            }],
        }),
    )
    .unwrap();
    let q = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let mut states = vec![vec![C64::ZERO; 10]; 2];
    states[0][1] = q;
    states[1][1] = C64::new(0.0, omega) * q;
    let mut force = vec![vec![C64::ZERO; 10]; 2];
    force[1][1] = q * C64::new(0.0, damping * omega);
    force[1][3] = (q * q * q).scale(cubic);
    (
        system,
        PeriodicPortInput::new(omega, force).unwrap(),
        orbit(omega, states),
    )
}

#[test]
fn finite_amplitude_peaks_follow_the_original_nonlinear_output_and_waveform() {
    let mut peaks = Vec::new();
    for amplitude in [0.15, 0.5] {
        let (system, input, anchor) = nonlinear_fixture(amplitude);
        let controls = controls(2, 2, 9, 11);
        let grid: Vec<_> = [0.97, 0.985, 1.0, 1.015, 1.03]
            .map(|s| s * anchor.omega)
            .to_vec();
        // Port 0 is the actual nonlinear restoring effort k*q + cubic*q^3,
        // not displacement or a fixed scalar multiple of that coordinate.
        let result =
            survey_forced_port_resonances(&system, &input, &anchor, 0, 1, &grid, &controls, || {
                false
            })
            .unwrap();
        assert_eq!(result.resonances.peaks.len(), 1);
        let peak = &result.resonances.peaks[0];
        let omega = peak.peak.angular_frequency_rad_s;
        peaks.push(omega);
        let solve = |w: f64, seed: &HbOrbit| {
            solve_forced_harmonic_balance_scaled_seeded(
                &system,
                &input.at_frequency(w).unwrap(),
                seed.packed_seed(&HbAnchor::Forced { omega: w }),
                &controls.harmonic_balance,
                &controls.harmonic_balance_scales,
            )
            .unwrap()
            .orbit
        };
        let corrected = solve(omega, &anchor);
        assert!(
            corrected.coeffs[0][3].abs() > 1e-9,
            "nonlinear waveform was frozen"
        );
        let observe = |candidate: &HbOrbit| {
            observe_port_spectrum(&system, candidate, &controls.spectrum, || false)
                .unwrap()
                .peak_phasors[0][0]
        };
        near(peak.observation, observe(&corrected), 2e-9);
        let linear_proxy = corrected.coeffs[0][1].conj().scale(8.0);
        assert!(
            (peak.observation - linear_proxy).abs() > 1e-3,
            "physical nonlinear output was replaced by a state proxy"
        );
        // Independent native solves verify the analytic output derivative;
        // finite differences belong only to this regression, never the survey.
        let step = 1e-5;
        let minus = observe(&solve(omega - step, &corrected));
        let plus = observe(&solve(omega + step, &corrected));
        near(
            peak.observation_derivative,
            (plus - minus).scale(0.5 / step),
            2e-5,
        );
        assert!(minus.abs() < peak.observation.abs() + 1e-9);
        assert!(plus.abs() < peak.observation.abs() + 1e-9);
    }
    assert!(
        peaks[1] > peaks[0] + 0.04,
        "finite-amplitude hardening was replaced by a small-signal peak"
    );
}

#[test]
fn branch_distance_source_clock_resolution_and_cancellation_are_not_bypassed() {
    let (system, input, anchor) = memory_fixture(1.6);
    let grid = [1.3, 1.5, 1.7, 1.9];
    let mut c = controls(3, 1, 3, 3);
    c.predictor_rms_tolerances = vec![1e-14; 3];
    assert!(matches!(
        survey_forced_port_resonances(&system, &input, &anchor, 0, 1, &grid, &c, || false,),
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
            PhsForcedResonanceError::PredictorDisagreement { .. }
        )))
    ));
    let c = controls(3, 1, 3, 3);
    assert!(
        survey_forced_port_resonances(
            &system,
            &input.at_frequency(1.61).unwrap(),
            &anchor,
            0,
            1,
            &grid,
            &c,
            || false,
        )
        .is_err()
    );
    let mut corrupted = anchor.clone();
    corrupted.coeffs[2][1].re += 0.1;
    assert!(matches!(
        survey_forced_port_resonances(&system, &input, &corrupted, 0, 1, &grid, &c, || false,),
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
            PhsForcedResonanceError::Tangent(_)
        )))
    ));
    assert!(
        survey_forced_port_resonances(&system, &input, &anchor, 0, 1, &grid, &c, || true,).is_err()
    );
    let (nonlinear, forcing, base) = nonlinear_fixture(0.5);
    let mut unresolved = controls(2, 2, 1, 11);
    // Input harmonic 3 and the exact nonlinear tangent must not be discarded.
    unresolved.predictor_rms_tolerances = vec![10.0; 2];
    assert!(
        survey_forced_port_resonances(
            &nonlinear,
            &forcing,
            &base,
            0,
            1,
            &[0.99 * base.omega, 1.01 * base.omega],
            &unresolved,
            || false,
        )
        .is_err()
    );
}
