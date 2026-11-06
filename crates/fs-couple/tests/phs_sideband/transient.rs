use super::{base_stability_controls, controls, near, orbit, scales};
use core::cell::Cell;
use fs_couple::phs_orbit::{
    PeriodicPortInput,
    sideband::{PhsBaseTransientControls, PhsBaseTransientError, PhsPortSidebandSource},
};
use fs_math::{c64::C64, det};
use fs_spectral::floquet::{FloquetClassification, transient::StroboscopicGrowthControls};

fn transient_controls() -> PhsBaseTransientControls {
    PhsBaseTransientControls {
        growth: StroboscopicGrowthControls {
            maximum_dimension: 16,
            maximum_horizons: 16,
            maximum_periods: 64,
            maximum_payload_bytes: 1 << 20,
            maximum_work: 100_000_000,
            svd_max_sweeps: 100,
            svd_tolerance: 2e-14,
            maximum_relative_residual: 1e-10,
        },
        log_gain_change_tolerance: 2e-8,
        maximum_workspace_bytes: 8 << 20,
        maximum_work: 1_000_000_000,
    }
}

// Independent exact continuous physical oscillator transition, then its
// explicitly authored S^-1*M*S chart. No numerical map power or SVD oracle.
fn exact_map(time: f64) -> [f64; 4] {
    let frequency = det::sqrt(9.0 - 0.2 * 0.2);
    let envelope = det::exp(-0.2 * time);
    let (cosine, sine) = (
        det::cos(frequency * time),
        det::sin(frequency * time) / frequency,
    );
    [
        envelope * (cosine + 0.2 * sine),
        envelope * sine / 0.02,
        -9.0 * envelope * sine * 0.02,
        envelope * (cosine - 0.2 * sine),
    ]
}

fn exact_gain(matrix: [f64; 4]) -> f64 {
    let squared = matrix.iter().map(|value| value * value).sum::<f64>();
    let determinant = (matrix[0] * matrix[3] - matrix[1] * matrix[2]).abs();
    det::sqrt(
        0.5 * (squared
            + det::sqrt((squared - 2.0 * determinant).max(0.0) * (squared + 2.0 * determinant))),
    )
}

#[test]
fn stable_physical_oscillator_has_resolved_nonnormal_growth_in_its_disclosed_chart() {
    let period = 0.3;
    let omega = core::f64::consts::TAU / period;
    let system = fs_phs::mass_spring_damper(1.0, 9.0, 0.4).unwrap();
    let force = C64::new(0.01, -0.003);
    let q = force / C64::new(9.0 - omega * omega, 0.4 * omega);
    let base = orbit(
        omega,
        vec![
            vec![C64::ZERO, q],
            vec![C64::ZERO, C64::new(0.0, omega) * q],
        ],
    );
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO, force]]).unwrap();
    let calls = Cell::new(0_usize);
    let hessian = |_: &[f64], direction: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out[0] = 9.0 * direction[0];
        out[1] = direction[1];
        true
    };
    let physical = controls(2, 1, 1, 1);
    let mut source =
        PhsPortSidebandSource::forced_analytic(&system, &input, &base, &physical, &hessian)
            .unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let horizons = [0, 1, 2, 4, 8, 16];
    let settings = transient_controls();
    let mut reference: Option<Vec<Vec<f64>>> = None;
    for unit in [1e-30, 1.0, 1e30] {
        let stability_controls = base_stability_controls(&[0.02 * unit, unit]);
        let stability = prepared
            .analyze_base_stability(&stability_controls, || false)
            .unwrap();
        assert_eq!(stability.classification, FloquetClassification::Stable);
        let original_maps = stability
            .levels
            .iter()
            .map(|level| level.monodromy.clone())
            .collect::<Vec<_>>();
        let before = calls.get();
        let result = stability
            .analyze_stroboscopic_growth(&horizons, &settings, || false)
            .unwrap();
        assert_eq!(calls.get(), before);
        assert_eq!(
            result.asymptotic_classification,
            FloquetClassification::Stable
        );
        assert!(result.full_maps_resolved && result.temporally_resolved);
        assert_eq!(
            result.state_scales,
            stability_controls.variational.state_scales
        );
        assert_eq!(result.period_s, stability.period_s);
        assert_eq!(result.omega_rad_s, omega);
        assert!(!result.autonomous);
        assert_eq!(result.levels.len(), 3);
        for (index, level) in result.levels.iter().enumerate() {
            assert_eq!(level.steps_per_period, 128 << index);
            assert_eq!(level.growth.dimension, 2);
            assert_eq!(
                level
                    .growth
                    .samples
                    .iter()
                    .map(|sample| sample.periods)
                    .collect::<Vec<_>>(),
                horizons
            );
            assert!(level.growth.matrix_products > 0);
            for (sample_index, sample) in level.growth.samples.iter().enumerate() {
                let matrix = exact_map(sample.periods as f64 * result.period_s);
                let gain = sample.gain.unwrap();
                let expected = exact_gain(matrix);
                assert!((gain - expected).abs() < 2e-7 * expected);
                assert!((sample.log_gain.unwrap() - det::ln(expected)).abs() < 2e-7);
                assert_eq!(sample.initial_direction.len(), 2);
                let final_direction = sample.final_direction.as_ref().unwrap();
                assert_eq!(final_direction.len(), 2);
                assert!(
                    (sample
                        .initial_direction
                        .iter()
                        .map(|value| value.norm_sq())
                        .sum::<f64>()
                        - 1.0)
                        .abs()
                        < 1e-11
                );
                assert!(
                    (final_direction
                        .iter()
                        .map(|value| value.norm_sq())
                        .sum::<f64>()
                        - 1.0)
                        .abs()
                        < 1e-11
                );
                // The complete reported directions attain the gain in THIS
                // physical chart, not in an inferred energy or modal norm.
                for row in 0..2 {
                    let evolved = sample.initial_direction[0].scale(matrix[2 * row])
                        + sample.initial_direction[1].scale(matrix[2 * row + 1]);
                    near(evolved, final_direction[row].scale(gain), 3e-7);
                }
                assert!(
                    sample.forward_relative_residual <= settings.growth.maximum_relative_residual
                );
                assert!(
                    sample.adjoint_relative_residual <= settings.growth.maximum_relative_residual
                );
                if let Some(reference) = &reference {
                    assert!(
                        (sample.log_gain.unwrap() - reference[index][sample_index]).abs() < 2e-10
                    );
                }
            }
            assert_eq!(level.growth.samples[0].gain, Some(1.0));
            assert_eq!(level.growth.samples[0].log_gain, Some(0.0));
            assert!(level.growth.samples[1].gain.unwrap() > 10.0);
            if index != 0 {
                assert!(
                    level
                        .log_gain_change_ratios
                        .as_ref()
                        .unwrap()
                        .iter()
                        .all(|&value| value <= 1.0)
                );
            }
        }
        reference = Some(
            result
                .levels
                .iter()
                .map(|level| {
                    level
                        .growth
                        .samples
                        .iter()
                        .map(|sample| sample.log_gain.unwrap())
                        .collect()
                })
                .collect(),
        );
        for (old, level) in original_maps.iter().zip(&stability.levels) {
            assert_eq!(old, &level.monodromy);
            assert_eq!(level.spectrum.spectrum.values.len(), 2);
            assert_eq!(level.spectrum.spectrum.right_vectors.len(), 2);
        }
        if unit == 1.0 {
            // Gain convergence has its own meaningful accuracy allowance, and
            // cannot be inferred just from stable, converged Floquet radii.
            let mut strict = settings.clone();
            strict.log_gain_change_tolerance = 1e-15;
            let unresolved = stability
                .analyze_stroboscopic_growth(&horizons, &strict, || false)
                .unwrap();
            assert!(unresolved.full_maps_resolved);
            assert!(!unresolved.temporally_resolved);
            assert_eq!(
                unresolved.asymptotic_classification,
                FloquetClassification::Stable
            );
            assert!(unresolved.levels.iter().skip(1).any(|level| {
                level
                    .log_gain_change_ratios
                    .as_ref()
                    .unwrap()
                    .iter()
                    .any(|&value| value > 1.0)
            }));
            for payload in [false, true] {
                let mut limited = settings.clone();
                if payload {
                    limited.maximum_workspace_bytes = result.admitted_workspace_bytes - 1;
                } else {
                    limited.maximum_work = result.admitted_work - 1;
                }
                assert!(matches!(
                    stability.analyze_stroboscopic_growth(&horizons, &limited, || false),
                    Err(PhsBaseTransientError::Budget(_))
                ));
            }
            let mut polls = 0;
            assert!(matches!(
                stability.analyze_stroboscopic_growth(&horizons, &settings, || {
                    polls += 1;
                    polls > 3
                }),
                Err(PhsBaseTransientError::Cancelled)
            ));
            assert!(
                stability
                    .analyze_stroboscopic_growth(&horizons, &settings, || false)
                    .unwrap()
                    .temporally_resolved
            );
            assert_eq!(calls.get(), before);

            // Even an adequate gain ladder cannot override its unresolved
            // original full-map comparison. Recompute through the same source.
            let mut strict_map = stability_controls.clone();
            strict_map.monodromy_change_tolerance = 1e-16;
            let parent = prepared
                .analyze_base_stability(&strict_map, || false)
                .unwrap();
            let unresolved = parent
                .analyze_stroboscopic_growth(&horizons, &settings, || false)
                .unwrap();
            assert!(!unresolved.full_maps_resolved && !unresolved.temporally_resolved);
            assert!(unresolved.levels.iter().skip(1).all(|level| {
                level
                    .log_gain_change_ratios
                    .as_ref()
                    .unwrap()
                    .iter()
                    .all(|&value| value <= 1.0)
            }));
        }
    }
}
