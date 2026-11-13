use super::{base_stability_controls, controls, near, orbit, scales};
use core::cell::Cell;
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PeriodicPortInput,
    sideband::{PhsBaseStabilityError, PhsPortSidebandSource},
};
use fs_math::{c64::C64, det};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};
use fs_spectral::floquet::FloquetClassification;

struct VelocityFeedback(f64);
impl AutonomousPortFeedback for VelocityFeedback {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, _: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        out[0] = -self.0 * flows[0];
        true
    }
}

fn rk4_multiplier(pole: C64, period: f64, steps: usize) -> C64 {
    let z = pole.scale(period / steps as f64);
    let z2 = z * z;
    let step =
        C64::ONE + z + z2.scale(0.5) + (z2 * z).scale(1.0 / 6.0) + (z2 * z2).scale(1.0 / 24.0);
    (0..steps).fold(C64::ONE, |value, _| value * step)
}

#[test]
fn same_velocity_magnitude_can_have_opposite_full_base_stability_without_changing_the_base() {
    let (mass, stiffness, passive_damping, omega) = (2.0, 8.0, 0.4, 2.7);
    let system = fs_phs::mass_spring_damper(mass, stiffness, passive_damping).unwrap();
    let input_force = C64::new(0.03, -0.01);
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO, input_force]]).unwrap();
    let probe_force = C64::new(0.7, -0.2);
    let carrier = 1.9;
    let mut first_magnitude: Option<f64> = None;
    for (gain, expected_class) in [
        (0.4, FloquetClassification::Stable),
        (-1.2, FloquetClassification::Unstable),
    ] {
        let damping = passive_damping + gain;
        let feedback = VelocityFeedback(gain);
        let q = input_force / C64::new(stiffness - mass * omega * omega, damping * omega);
        let base = orbit(
            omega,
            vec![
                vec![C64::ZERO, q],
                vec![C64::ZERO, C64::new(0.0, mass * omega) * q],
            ],
        );
        let original = base.coeffs.clone();
        let calls = Cell::new(0_usize);
        let refused = Cell::new(false);
        let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
            calls.set(calls.get() + 1);
            out.fill(0.0);
            out[1] = -gain / mass;
            !refused.get()
        };
        let physical = controls(2, 1, 1, 1);
        let mut source = PhsPortSidebandSource::forced_feedback(
            &system,
            &input,
            &feedback,
            &derivative,
            &base,
            &physical,
        )
        .unwrap();
        let mut prepared = source.prepare(&scales(2), || false).unwrap();
        let settings = base_stability_controls(&[0.5, 1.5]);
        let before = calls.get();
        let report = prepared
            .analyze_base_stability(&settings, || false)
            .unwrap();
        assert_eq!(report.classification, expected_class);
        assert!(!report.autonomous);
        assert!(report.scaled_phase_flow.is_none());
        assert_eq!(report.state_scales, settings.variational.state_scales);
        assert_eq!(
            report.period_s.to_bits(),
            (core::f64::consts::TAU / omega).to_bits()
        );
        assert_eq!(calls.get() - before, report.original_jacobian_evaluations);
        assert_eq!(
            report
                .levels
                .iter()
                .map(|l| l.steps_per_period)
                .collect::<Vec<_>>(),
            [128, 256, 512]
        );
        let real = -damping / (2.0 * mass);
        for level in &report.levels {
            assert_eq!(level.monodromy.len(), 4);
            assert_eq!(level.spectrum.spectrum.values.len(), 2);
            assert_eq!(level.spectrum.spectrum.right_vectors.len(), 2);
            assert!(level.spectrum.phase_multiplier_index.is_none());
            assert!(level.spectrum.phase_residual.is_none());
            assert_eq!(level.linear_rhs_evaluations, 4 * level.steps_per_period);
            assert_eq!(
                level.linear_jacobian_evaluations,
                4 * level.steps_per_period
            );
            // The midpoint and identical RHS/Jacobian times reuse the SAME
            // complete original tangent. No original callback is spent twice.
            assert!(level.original_jacobian_evaluations >= 2 * level.steps_per_period + 1);
            assert!(level.original_jacobian_evaluations <= 3 * level.steps_per_period);
            for sign in [-1.0, 1.0] {
                let pole = C64::new(real, sign * det::sqrt(stiffness / mass - real * real));
                let expected = rk4_multiplier(pole, report.period_s, level.steps_per_period);
                assert!(
                    level
                        .spectrum
                        .spectrum
                        .values
                        .iter()
                        .any(|&z| (z - expected).abs() < 3e-11)
                );
            }
            assert!(
                level
                    .spectrum
                    .spectrum
                    .relative_residuals
                    .iter()
                    .all(|&r| r <= settings.floquet.spectrum.maximum_relative_residual)
            );
        }
        assert!((report.dominant_growth_rate_per_s.unwrap() - real).abs() < 1e-8);
        if gain > 0.0 {
            // A complete map comparison is an independent gate, even when all
            // guarded radii and their deliberately looser changes look stable.
            let mut strict = settings.clone();
            strict.monodromy_change_tolerance = 1e-16;
            strict.spectral_radius_change_tolerance = 1e-3;
            let unresolved = prepared.analyze_base_stability(&strict, || false).unwrap();
            assert_eq!(unresolved.classification, FloquetClassification::Unresolved);
            assert!(unresolved.levels[1].monodromy_change_ratio.unwrap() > 1.0);
            assert!(unresolved.levels[1].spectral_radius_change_ratio.unwrap() <= 1.0);
            for payload in [false, true] {
                let mut limited = settings.clone();
                if payload {
                    limited.maximum_workspace_bytes = report.admitted_workspace_bytes - 1;
                } else {
                    limited.maximum_work = report.admitted_work - 1;
                }
                let before = calls.get();
                assert!(matches!(
                    prepared.analyze_base_stability(&limited, || false),
                    Err(PhsBaseStabilityError::Budget(_))
                ));
                assert_eq!(calls.get(), before);
            }
            let before = calls.get();
            assert!(matches!(
                prepared.analyze_base_stability(&settings, || calls.get() > before + 2),
                Err(PhsBaseStabilityError::Cancelled)
            ));
            refused.set(true);
            assert!(matches!(
                prepared.analyze_base_stability(&settings, || false),
                Err(PhsBaseStabilityError::SourceRefused { .. })
            ));
            refused.set(false);
            assert_eq!(
                prepared
                    .analyze_base_stability(&settings, || false)
                    .unwrap()
                    .classification,
                FloquetClassification::Stable
            );
        }
        // Stability analysis does not change the original cached probe model.
        let response = prepared
            .factor(carrier, || false)
            .unwrap()
            .solve(&[C64::ZERO, probe_force, C64::ZERO], || false)
            .unwrap();
        let expected = C64::new(0.0, carrier) * probe_force
            / C64::new(stiffness - mass * carrier * carrier, damping * carrier);
        near(response.port_coefficients[1], expected, 3e-12);
        let magnitude = response.port_coefficients[1].abs();
        if let Some(first) = first_magnitude {
            assert!((magnitude - first).abs() < 3e-12);
        } else {
            first_magnitude = Some(magnitude);
        }
        assert_eq!(base.coeffs, original);
    }
}

#[test]
fn driven_maxwell_spectrum_preserves_three_original_modes_and_forced_unit_mode_is_not_removed() {
    // The original closed-loop polynomial is (s+1)(s+2)(s+3), with physical
    // passive damping1.5 and separate feedback damping0.5, not a replacement R.
    let (stiffness, arm, tau, omega) = (1.5, 1.5_f64, 0.25, core::f64::consts::TAU);
    let system = fs_phs::mass_spring_damper(1.0, stiffness, 1.5)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: arm,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let feedback = VelocityFeedback(0.5);
    let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        out[1] = -0.5;
        true
    };
    let force = C64::new(1e-4, 3e-5);
    let s = C64::new(0.0, omega);
    let relaxation = C64::ONE + s.scale(tau);
    let q = force
        / (C64::new(stiffness - omega * omega, 2.0 * omega)
            + (s.scale(tau) / relaxation).scale(arm));
    let base = orbit(
        omega,
        vec![
            vec![C64::ZERO, q],
            vec![C64::ZERO, s * q],
            vec![C64::ZERO, q.scale(det::sqrt(arm)) / relaxation],
        ],
    );
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO, force]]).unwrap();
    let physical = controls(3, 1, 1, 1);
    let mut source = PhsPortSidebandSource::forced_feedback(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        &physical,
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(3), || false).unwrap();
    let mut reference: Option<Vec<Vec<f64>>> = None;
    for unit in [1e-30, 1.0, 1e30] {
        let settings = base_stability_controls(&[0.1 * unit, unit, 0.2 * unit]);
        let report = prepared
            .analyze_base_stability(&settings, || false)
            .unwrap();
        assert_eq!(report.classification, FloquetClassification::Stable);
        assert_eq!(report.period_s, 1.0);
        for (index, level) in report.levels.iter().enumerate() {
            assert_eq!(level.monodromy.len(), 9);
            assert_eq!(level.spectrum.spectrum.values.len(), 3);
            for pole in [-1.0, -2.0, -3.0] {
                let expected = rk4_multiplier(C64::from_re(pole), 1.0, level.steps_per_period);
                assert!(
                    level
                        .spectrum
                        .spectrum
                        .values
                        .iter()
                        .any(|&z| (z - expected).abs() < 3e-11)
                );
            }
            if let Some(reference) = &reference {
                for (&a, &b) in level.monodromy.iter().zip(&reference[index]) {
                    assert!((a - b).abs() < 3e-11);
                }
            }
        }
        reference = Some(
            report
                .levels
                .iter()
                .map(|level| level.monodromy.clone())
                .collect(),
        );
    }
    // Even zero imposed input preserves the externally authored clock. A unit
    // mode is not an autonomous phase direction and cannot be excluded.
    let integrator = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let base = orbit(omega, vec![vec![C64::ZERO]]);
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO]]).unwrap();
    let mut source =
        PhsPortSidebandSource::forced(&integrator, &input, &base, &controls(1, 1, 0, 0)).unwrap();
    let report = source
        .analyze_base_stability(&scales(1), &base_stability_controls(&[1.0]), || false)
        .unwrap();
    assert_eq!(report.classification, FloquetClassification::Unresolved);
    assert!(report.scaled_phase_flow.is_none());
    for level in &report.levels {
        assert_eq!(level.spectrum.spectrum.values, [C64::ONE]);
        assert_eq!(level.spectrum.spectral_radius, 1.0);
        assert_eq!(level.spectrum.effective_spectral_radius, Some(1.0));
        assert!(level.spectrum.phase_multiplier_index.is_none());
    }
}
