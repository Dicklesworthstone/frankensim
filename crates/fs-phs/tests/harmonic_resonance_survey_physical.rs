//! Native band surveys retain observed resonance shifts, mutual damping and wave interference.
use fs_math::c64::C64;
use fs_phs::harmonic::resonance::{
    DrivenResonanceError, HarmonicObservation, ResonanceSearchError, ResonanceSurveyControls,
    ResonanceSurveyError,
};
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel, ScatteringObservation};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};
use fs_phs::{PortHamiltonian, QuadraticStorage};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 16,
        max_dense_bytes: 1_000_000,
        max_dense_work: 1_000_000,
        relative_residual_tolerance: 2e-12,
    }
}

fn controls() -> ResonanceSurveyControls {
    ResonanceSurveyControls {
        absolute_frequency_tolerance_rad_s: 1e-10,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        maximum_peaks: 4,
    }
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-30),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn two_physical_peaks_survive_mutual_damping_and_an_intervening_antiresonance() {
    let c = 0.4;
    // q1,p1,q2,p2, with one real force/velocity port and rank-one C=c*b*b^T.
    // Retaining the mutual loss is essential: Y=-iw(D1+D2)/(D1*D2-iwc(D1+D2)).
    let system = PortHamiltonian::new(
        4,
        1,
        vec![
            0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0,
        ],
        vec![
            0.0, 0.0, 0.0, 0.0, 0.0, c, 0.0, c, 0.0, 0.0, 0.0, 0.0, 0.0, c, 0.0, c,
        ],
        vec![0.0, 1.0, 0.0, 1.0],
        Box::new(
            QuadraticStorage::new(
                vec![
                    25.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 175.0, 0.0, 0.0, 0.0, 0.0,
                    1.0,
                ],
                4,
            )
            .unwrap(),
        ),
    )
    .unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let input = [C64::new(0.8, -0.3)];
    let grid = [2.0, 4.0, 5.0, 7.0, 10.0, 12.0, 14.0, 16.0];
    let result = model
        .driven_resonances(
            &input,
            HarmonicObservation::Output(&[1.0]),
            &grid,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result.peaks.len(), 2);
    assert_eq!(result.grid_points, grid.len());
    assert_eq!(result.interval_rad_s, (2.0, 16.0));
    for (peak, expected) in result.peaks.iter().zip([5.0, 175.0_f64.sqrt()]) {
        let w = peak.peak.angular_frequency_rad_s;
        assert!((w - expected).abs() <= 1e-10, "{w} != {expected}");
        let (d1, d2) = (25.0 - w * w, 175.0 - w * w);
        let exact =
            input[0] * C64::new(0.0, -w * (d1 + d2)) / C64::new(d1 * d2, -c * w * (d1 + d2));
        close(peak.observation, exact);
        assert!(peak.half_power.is_none());
    }
    assert!(result.maximum_relative_residual <= budget().relative_residual_tolerance);
    assert!(result.evaluations <= controls().maximum_evaluations);
    let silent = model
        .driven_resonances(
            &[C64::ZERO],
            HarmonicObservation::Output(&[1.0]),
            &grid,
            controls(),
            || false,
        )
        .unwrap();
    assert!(silent.peaks.is_empty());
}

#[test]
fn state_effort_and_port_observations_have_their_own_physical_peak() {
    let (m, k, c) = (2.0_f64, 200.0_f64, 8.0_f64);
    let system = fs_phs::mass_spring_damper(m, k, c).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let input = [C64::new(-0.3, 0.8)];
    for (observation, expected, kind) in [
        (
            HarmonicObservation::State(&[1.0, 0.0]),
            (k / m - c * c / (2.0 * m * m)).sqrt(),
            0,
        ),
        (
            HarmonicObservation::Effort(&[1.0, 0.0]),
            (k / m - c * c / (2.0 * m * m)).sqrt(),
            1,
        ),
        (HarmonicObservation::Output(&[-2.0]), (k / m).sqrt(), 2),
    ] {
        let result = model
            .driven_resonances(
                &input,
                observation,
                &[4.0, 8.0, 10.0, 12.0, 15.0],
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(result.peaks.len(), 1);
        let peak = &result.peaks[0];
        let w = peak.peak.angular_frequency_rad_s;
        assert!((w - expected).abs() <= 1e-10);
        let q = input[0] / C64::new(k - m * w * w, -c * w);
        let exact = match kind {
            0 => q,
            1 => q.scale(k),
            _ => q * C64::new(0.0, 2.0 * w),
        };
        close(peak.observation, exact);
    }
}

#[test]
fn wave_surveys_retain_boundary_reaction_and_coherent_incident_feedthrough() {
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.2],
        vec![0.0, 0.0, 1.0, 1.0],
        Box::new(QuadraticStorage::new(vec![100.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap();
    let model = ScatteringModel::new(
        HarmonicModel::from_ode(&system, budget(), || false).unwrap(),
        &[(1, 1.0), (0, 1.0)], // Wave order intentionally differs from source order.
        ScatteringBudget {
            maximum_ports: 2,
            maximum_payload_bytes: 500_000,
            maximum_work: 1_000_000,
        },
        || false,
    )
    .unwrap();
    let gamma = 2.2_f64;
    let outgoing_x = (-0.4 - (0.16 + 4.0 * gamma * gamma).sqrt()) / 2.0;
    let a = gamma * gamma - (gamma - 1.0).powi(2) - 1.0;
    let input_x = (-a - (a * a + 4.0 * gamma * gamma).sqrt()) / 2.0;
    let from_x = |x: f64| ((x * x + 400.0).sqrt() - x) / 2.0;
    let displacement_peak = (100.0 - gamma * gamma / 2.0).sqrt();
    for (observation, expected, kind) in [
        (
            ScatteringObservation::State(&[1.0, 0.0]),
            displacement_peak,
            0,
        ),
        (
            ScatteringObservation::Effort(&[1.0, 0.0]),
            displacement_peak,
            1,
        ),
        (
            ScatteringObservation::Input(&[1.0, 0.0]),
            from_x(input_x),
            2,
        ),
        (ScatteringObservation::Output(&[1.0, 0.0]), 10.0, 3),
        (
            ScatteringObservation::Outgoing(&[0.0, 1.0]),
            from_x(outgoing_x),
            4,
        ),
    ] {
        let result = model
            .driven_resonances(
                &[C64::ONE, C64::new(0.0, 1.0)],
                observation,
                &[7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 16.0],
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(result.peaks.len(), 1);
        let peak = &result.peaks[0];
        let w = peak.peak.angular_frequency_rad_s;
        assert!(
            (w - expected).abs() <= 1e-10,
            "observation {kind}: {w} != {expected}"
        );
        let v = C64::new(2.0, 2.0) / C64::new(gamma, 100.0 / w - w);
        let exact = match kind {
            0 => v * C64::new(0.0, 1.0 / w),
            1 => v * C64::new(0.0, 100.0 / w),
            2 => C64::new(2.0, 0.0) - v,
            3 => v,
            _ => v - C64::ONE,
        };
        close(peak.observation, exact);
        assert!(result.maximum_relative_residual <= budget().relative_residual_tolerance);
    }
}

#[test]
fn invalid_physical_observations_cancellation_and_true_source_poles_refuse() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 0.0).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let survey = |observation, stop| {
        model.driven_resonances(
            &[C64::ONE],
            observation,
            &[9.0, 10.0, 11.0],
            controls(),
            || stop,
        )
    };
    assert!(matches!(
        survey(HarmonicObservation::State(&[1.0]), false),
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Invalid(
            _
        )))
    ));
    assert!(matches!(
        survey(HarmonicObservation::Output(&[1.0]), true),
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
            DrivenResonanceError::Harmonic(HarmonicError::Cancelled)
        )))
    ));
    assert!(matches!(
        survey(HarmonicObservation::Output(&[1.0]), false),
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
            DrivenResonanceError::Harmonic(HarmonicError::Singular { .. })
        )))
    ));
}
