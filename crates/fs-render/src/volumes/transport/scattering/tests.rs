//! Analytic radiometry, spectral-mixture expectations and phase normalization.

use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::AbsorbingGridBudget;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn budget(steps: u32) -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 4096,
        max_steps_per_segment: steps,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-100),
        "{actual:e} != {expected:e} (relative tolerance {relative:e})"
    );
}

// Independent direct SI equation, not the production log-domain Planck helper.
fn planck_si(temperature: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let lambda = wavelength_nm * 1e-9;
    2.0 * h * c * c / lambda.powi(5) / ((h * c / (lambda * k * temperature)).exp() - 1.0) * 1e-9
}

#[test]
fn g1_spectral_free_flight_matches_escape_and_scattering_integrals() {
    with_cx(|cx, _| {
        let density = [2.0, 0.0, 1.0, 3.0];
        let wavelengths = [380.0, 580.0, 780.0];
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 4], &density, [0.0; 3], [1.0, 1.0, 0.25]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (580.0, 0.3), (780.0, 1.0)]).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::table(&[(380.0, 0.0), (580.0, 0.7), (780.0, 2.0)]).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        assert!(core::ptr::eq(
            medium.absorption.grid.data.as_ptr(),
            density.as_ptr()
        ));
        let origin = Point3::new(0.5, 0.5, -0.5);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let transmitted = medium
            .transmittance(cx, origin, direction, f64::INFINITY, &wavelengths)
            .unwrap();
        let expected_escape = [1.0, (-1.5_f64).exp(), (-4.5_f64).exp()];
        let expected_scatter = [
            0.0,
            0.7 * (1.0 - expected_escape[1]),
            (2.0 / 3.0) * (1.0 - expected_escape[2]),
        ];
        for lane in 0..3 {
            close(transmitted[lane], expected_escape[lane], 2e-13);
        }
        let mut escaped = [0.0; 3];
        let mut scattered = [0.0; 3];
        const SAMPLES: usize = 8192;
        for hero in 0..3 {
            for sample in 0..SAMPLES {
                let result = medium
                    .sample_segment(
                        cx,
                        origin,
                        direction,
                        f64::INFINITY,
                        &wavelengths,
                        (hero as f64 + 0.5) / 3.0,
                        (sample as f64 + 0.5) / SAMPLES as f64,
                    )
                    .unwrap();
                let accumulator = if let Some(event) = result.event {
                    assert!((0.0..1.0).contains(&event.position.z));
                    assert!(
                        !(0.25..0.5).contains(&event.position.z),
                        "event in vacuum gap"
                    );
                    close(event.distance_m, event.position.z + 0.5, 2e-13);
                    &mut scattered
                } else {
                    close(
                        result.escape_probability,
                        expected_escape.iter().sum::<f64>() / 3.0,
                        2e-13,
                    );
                    &mut escaped
                };
                for lane in 0..3 {
                    accumulator[lane] += result.weight[lane] / (3 * SAMPLES) as f64;
                }
            }
        }
        for lane in 0..3 {
            assert!((escaped[lane] - expected_escape[lane]).abs() < 2e-4);
            assert!((scattered[lane] - expected_scatter[lane]).abs() < 2e-4);
        }
    });
}

#[test]
fn g1_lte_emission_uses_absorption_and_the_complete_segment() {
    with_cx(|cx, _| {
        let density = [2.0, 1.0];
        let temperature = [1600.0, 2200.0];
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0, 1.0, 0.5]),
            MassAbsorptionSpectrum::constant(0.25).unwrap(),
            Some(&temperature),
            budget(16),
        )
        .unwrap();
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(0.75).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        assert!(medium.has_thermal_emission());
        let origin = Point3::new(0.5, 0.5, -1.0);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let collision = medium
            .sample_segment(
                cx,
                origin,
                direction,
                f64::INFINITY,
                &[550.0],
                0.5,
                (-0.4_f64).exp(),
            )
            .unwrap();
        let escape = medium
            .sample_segment(
                cx,
                origin,
                direction,
                f64::INFINITY,
                &[550.0],
                0.5,
                (-5.0_f64).exp(),
            )
            .unwrap();
        close(collision.event.unwrap().distance_m, 1.2, 2e-13);
        close(collision.weight[0], 0.75, 2e-13);
        assert!(escape.event.is_none());
        close(escape.weight[0], 1.0, 2e-13);
        close(escape.escape_probability, (-1.5_f64).exp(), 2e-13);
        let source = 0.25
            * (planck_si(1600.0, 550.0) * (1.0 - (-1.0_f64).exp())
                + (-1.0_f64).exp() * planck_si(2200.0, 550.0) * (1.0 - (-0.5_f64).exp()));
        close(collision.emission_w_m2_sr_nm[0], source, 2e-12);
        assert_eq!(collision.emission_w_m2_sr_nm, escape.emission_w_m2_sr_nm);
        let deterministic = medium
            .transport(cx, origin, direction, f64::INFINITY, &[550.0])
            .unwrap();
        close(deterministic.transmittance[0], (-1.5_f64).exp(), 2e-13);
        assert_eq!(
            deterministic.emission_w_m2_sr_nm, collision.emission_w_m2_sr_nm,
            "source-light transport must share the complete free-flight source integral"
        );
    });
}

#[test]
fn g0_zero_scattering_and_vacuum_preserve_deterministic_absorption() {
    with_cx(|cx, _| {
        for rho in [0.0, 2.0] {
            let density = [rho];
            let absorption = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::constant(0.3).unwrap(),
                Some(&[1800.0]),
                budget(4),
            )
            .unwrap();
            let medium = ScatteringGrid::try_new(
                cx,
                absorption,
                MassScatteringSpectrum::constant(0.0).unwrap(),
                PhaseFunction::Isotropic,
            )
            .unwrap();
            assert!(!medium.has_scattering());
            assert_eq!(medium.is_vacuum(), rho == 0.0);
            let origin = Point3::new(0.5, 0.5, -1.0);
            let direction = Vec3::new(0.0, 0.0, 1.0);
            let wavelengths = [400.0, 500.0, 600.0, 700.0];
            let reference = medium
                .absorption_grid()
                .transport(cx, origin, direction, f64::INFINITY, &wavelengths)
                .unwrap();
            let sample = medium
                .sample_segment(cx, origin, direction, f64::INFINITY, &wavelengths, 0.3, 0.7)
                .unwrap();
            assert_eq!(sample.weight, reference.transmittance);
            assert_eq!(sample.emission_w_m2_sr_nm, reference.emission_w_m2_sr_nm);
            assert_eq!(sample.event, None);
            assert_eq!(sample.escape_probability, 1.0);
        }
    });
}

#[test]
fn g0_absorption_only_spectral_lanes_retain_sampled_survival_density() {
    with_cx(|cx, _| {
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &[1.0], [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            None,
            budget(4),
        )
        .unwrap();
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::table(&[(380.0, 0.0), (580.0, 0.0), (780.0, 1.0)]).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let trace = |u| {
            medium
                .sample_segment(
                    cx,
                    Point3::new(0.5, 0.5, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    1.0,
                    &[400.0, 500.0],
                    0.5,
                    u,
                )
                .unwrap()
        };
        let escaped = trace(0.1);
        assert!(escaped.event.is_none());
        assert_eq!(escaped.weight, [1.0; 2]);
        close(escaped.escape_probability, (-1.0_f64).exp(), 2e-13);
        let absorbed = trace(0.9);
        assert!(absorbed.event.is_some());
        assert_eq!(absorbed.weight, [0.0; 2]);
    });
}

#[test]
fn g1_log_mixture_weights_survive_underflow_and_pure_scatter_has_no_lte_source() {
    // A subnormal intermediate product must not lose percent-level accuracy
    // when the final multiplication/division restores a normal result.
    close(positive_product(1e-200, 5e-124, 1e300), 5e-24, 2e-12);
    close(divide_product(1e-16, 1e-200, 5e-124), 2e307, 2e-12);
    with_cx(|cx, _| {
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &[1e300], [0.0; 3], [1000.0; 3]),
            MassAbsorptionSpectrum::constant(0.0).unwrap(),
            Some(&[2000.0]),
            budget(4),
        )
        .unwrap();
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::table(&[(380.0, 1e-300), (580.0, 2e-300), (780.0, 3e-300)])
                .unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        assert!(!medium.has_thermal_emission());
        let sample = medium
            .sample_segment(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                1000.0,
                &[380.0, 580.0, 780.0],
                0.0,
                (-600.0_f64).exp(),
            )
            .unwrap();
        close(sample.event.unwrap().distance_m, 600.0, 2e-13);
        close(sample.weight[0], 3.0, 2e-13);
        assert!(sample.weight.iter().all(|x| x.is_finite() && *x >= 0.0));
        assert_eq!(sample.emission_w_m2_sr_nm, [0.0; 3]);
    });
}

#[test]
fn g1_phase_normalization_sampling_moments_and_direction_convention() {
    let cases = [
        (PhaseFunction::Isotropic, 0.0, 1.0 / 3.0),
        (
            PhaseFunction::HenyeyGreenstein { g: -0.6 },
            -0.6,
            (1.0 + 2.0 * 0.36) / 3.0,
        ),
        (
            PhaseFunction::HenyeyGreenstein { g: 1e-12 },
            1e-12,
            1.0 / 3.0,
        ),
        (
            PhaseFunction::HenyeyGreenstein { g: 0.6 },
            0.6,
            (1.0 + 2.0 * 0.36) / 3.0,
        ),
        (PhaseFunction::Rayleigh, 0.0, 0.4),
    ];
    let direction = Vec3::new(0.0, 0.0, 1.0);
    const SAMPLES: usize = 8192;
    for (phase, expected_mean, expected_second) in cases {
        let (mut integral, mut mean, mut second) = (0.0, 0.0, 0.0);
        for i in 0..SAMPLES {
            let u = (i as f64 + 0.5) / SAMPLES as f64;
            integral += phase.eval_cos(2.0 * u - 1.0).unwrap() * 4.0 * core::f64::consts::PI
                / SAMPLES as f64;
            let sample = phase
                .sample(direction, u, (u * 1.618_033_988_749_895).fract())
                .unwrap();
            let c = direction.dot(sample.direction);
            mean += c / SAMPLES as f64;
            second += c * c / SAMPLES as f64;
            close(sample.direction.dot(sample.direction), 1.0, 1e-14);
            assert_eq!(sample.pdf, phase.eval_cos(c).unwrap());
            if phase == PhaseFunction::Rayleigh {
                close(0.5 + 0.375 * (c + c * c * c / 3.0), u, 2e-11);
            }
        }
        assert!((integral - 1.0).abs() < 5e-7);
        assert!((mean - expected_mean).abs() < 2e-7);
        assert!((second - expected_second).abs() < 3e-7);
    }
    for axis in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0)] {
        let phase = PhaseFunction::HenyeyGreenstein { g: 0.9 };
        let sample = phase.sample(axis, 0.5, 0.25).unwrap();
        assert!(
            axis.dot(sample.direction) > 0.9,
            "positive g must be forward"
        );
        close(
            sample.pdf,
            phase.eval_cos(axis.dot(sample.direction)).unwrap(),
            2e-13,
        );
    }
}

#[test]
fn g0_g4_admission_identity_sampling_budgets_and_cancellation() {
    assert_eq!(
        MassScatteringSpectrum::constant(-1.0),
        Err(VolumeTransportError::ScatteringSpectrum)
    );
    assert_eq!(
        PhaseFunction::henyey_greenstein(1.0),
        Err(VolumeTransportError::Phase)
    );
    assert_eq!(
        PhaseFunction::Isotropic.eval_cos(2.0),
        Err(VolumeTransportError::Phase)
    );
    with_cx(|cx, gate| {
        let make = |phase, steps, temperature| {
            let absorption = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], &[1.0, 1.0], [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::constant(0.5).unwrap(),
                temperature,
                budget(steps),
            )
            .unwrap();
            ScatteringGrid::try_new(
                cx,
                absorption,
                MassScatteringSpectrum::constant(0.5).unwrap(),
                phase,
            )
            .unwrap()
        };
        let medium = make(PhaseFunction::Isotropic, 1, None);
        assert_eq!(
            medium.content_identity(),
            make(PhaseFunction::Isotropic, 5, None).content_identity()
        );
        assert_ne!(
            medium.content_identity(),
            make(PhaseFunction::Rayleigh, 1, None).content_identity()
        );
        let origin = Point3::new(0.5, 0.5, 0.0);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        assert_eq!(
            medium.sample_segment(cx, origin, direction, 2.0, &[550.0], 0.5, 0.0),
            Err(VolumeTransportError::Budget)
        );
        assert!(
            medium
                .sample_segment(cx, origin, direction, 2.0, &[550.0], 0.5, 0.9)
                .is_ok(),
            "a nonemissive event needs only the cells preceding it"
        );
        let emissive = make(PhaseFunction::Isotropic, 1, Some(&[1800.0; 2]));
        assert_eq!(
            emissive.sample_segment(cx, origin, direction, 2.0, &[550.0], 0.5, 0.9),
            Err(VolumeTransportError::Budget),
            "complete emission still visits the second cell"
        );
        assert_eq!(
            medium.sample_segment(cx, origin, direction, 2.0, &[550.0], 1.0, 0.5),
            Err(VolumeTransportError::Sample)
        );
        assert_eq!(
            medium.sample_segment(cx, origin, direction, 2.0, &[], 0.5, 0.5),
            Err(VolumeTransportError::Sample)
        );
        gate.request();
        assert_eq!(
            medium.transmittance(cx, origin, direction, 2.0, &[550.0]),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
