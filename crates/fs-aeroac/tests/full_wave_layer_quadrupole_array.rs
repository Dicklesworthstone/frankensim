//! G1/G3/G4: complete original stress arrays and acquired factors around a layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::recorded::{
    QuadrupoleRecordingBudget, QuadrupoleSampleClock, QuadrupoleWelchStream,
};
use fs_aeroac::curle3d::quadrupole::{QuadrupolePhasor, resolved_reynolds_moment_n_m};
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::PreparedInteriorFullWaveQuadrupole;
use fs_aeroac::curle3d::refraction::full_wave::layer::reflection::PreparedLayeredReflectedFullWaveQuadrupole;
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    LayerQuadrupoleField, PlanarFluidLayer, PreparedLayeredFullWaveQuadrupole,
    PreparedLayeredFullWaveQuadrupoleObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::FullWaveQuadrupoleSourceLimits;
use fs_aeroac::curle3d::{CoherenceBudget, CompactFieldStatistics};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

const SOURCES: [[f64; 3]; 2] = [[-0.01, 0.0, -0.025], [0.02, -0.015, -0.04]];
const KINDS: [LayerQuadrupoleField; 4] = [
    LayerQuadrupoleField::Transmitted,
    LayerQuadrupoleField::Reflected,
    LayerQuadrupoleField::IncidentTotal,
    LayerQuadrupoleField::Interior,
];
fn receiver(kind: LayerQuadrupoleField) -> [f64; 3] {
    [
        0.035,
        0.01,
        match kind {
            LayerQuadrupoleField::Transmitted => 0.075,
            LayerQuadrupoleField::Interior => 0.007,
            _ => -0.06,
        },
    ]
}
fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn layer(same: bool) -> PlanarFluidLayer {
    let a = fluid(1.2, 300.0);
    PlanarFluidLayer::new(
        a,
        if same { a } else { fluid(2.4, 600.0) },
        if same { a } else { fluid(1.8, 400.0) },
        [0.0; 3],
        [0.0, 0.0, 1.0],
        0.018,
    )
    .unwrap()
}
fn limits() -> FullWaveQuadrupoleSourceLimits {
    FullWaveQuadrupoleSourceLimits {
        max_sources: 4,
        maximum_total_evaluations: 1_048_576,
        quadrupole: FullWaveQuadrupoleLimits {
            absolute_pressure_tolerance_pa_per_n_m: 2e-6,
            absolute_velocity_tolerance_m_s_per_n_m: 2e-9,
            relative_tolerance: 1e-9,
            maximum_evaluations: 262_144,
            maximum_depth: 40,
        },
    }
}
fn coherence() -> CoherenceBudget {
    CoherenceBudget {
        max_modes: 8,
        max_source_mode_pairs: 32,
    }
}
fn moment(cell: usize) -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: core::array::from_fn(|i| {
            C64::new(1e-4 * (i + 1) as f64, 2e-5 * (cell + i + 1) as f64)
        }),
    }
}
fn near(a: C64, b: C64) {
    assert!((a - b).abs() <= 3e-12 * b.abs() + 1e-13, "{a:?} != {b:?}");
}
fn statistics_equal(a: &CompactFieldStatistics, b: &CompactFieldStatistics) {
    near(
        C64::from_re(a.mean_square_pressure_pa2),
        C64::from_re(b.mean_square_pressure_pa2),
    );
    for i in 0..3 {
        near(
            C64::new(a.active_intensity_w_m2[i], a.reactive_intensity_w_m2[i]),
            C64::new(b.active_intensity_w_m2[i], b.reactive_intensity_w_m2[i]),
        );
        for j in 0..3 {
            near(
                a.velocity_covariance_m2_s2[i][j],
                b.velocity_covariance_m2_s2[i][j],
            );
        }
    }
}

#[test]
fn all_four_domains_sum_original_tensor_fields_before_flux_and_statistics() {
    let slab = layer(false);
    let moments = [moment(0), moment(1)];
    for kind in KINDS {
        let position = receiver(kind);
        let array = PreparedLayeredFullWaveQuadrupoleObserver::new(
            slab,
            kind,
            700.0,
            position,
            &SOURCES,
            limits(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(array.field_kind(), kind);
        assert_eq!(array.source_count(), 2);
        assert_eq!(array.source_position_m(1), Some(SOURCES[1]));
        assert_eq!(array.source_position_m(2), None);
        assert_eq!(array.observer_position_m(), position);
        let mut expected = [C64::ZERO; 4];
        let mut evaluations = 0;
        for (index, (&source, &moment)) in SOURCES.iter().zip(&moments).enumerate() {
            let (line, diagnostics) = match kind {
                LayerQuadrupoleField::Transmitted => {
                    let native = PreparedLayeredFullWaveQuadrupole::new(
                        slab,
                        700.0,
                        source,
                        position,
                        limits().quadrupole,
                        &mut || false,
                    )
                    .unwrap();
                    (
                        native.evaluate_harmonic_line(moment).unwrap(),
                        *native.diagnostics(),
                    )
                }
                LayerQuadrupoleField::Interior => {
                    let native = PreparedInteriorFullWaveQuadrupole::new(
                        slab,
                        700.0,
                        source,
                        position,
                        limits().quadrupole,
                        &mut || false,
                    )
                    .unwrap();
                    (
                        native.evaluate_harmonic_line(moment).unwrap(),
                        *native.diagnostics(),
                    )
                }
                _ => {
                    let native = PreparedLayeredReflectedFullWaveQuadrupole::new(
                        slab,
                        700.0,
                        source,
                        position,
                        limits().quadrupole,
                        &mut || false,
                    )
                    .unwrap();
                    let value = if kind == LayerQuadrupoleField::Reflected {
                        native.evaluate_reflected_harmonic_line(moment)
                    } else {
                        native.evaluate_total_harmonic_line(moment)
                    };
                    (value.unwrap(), *native.diagnostics())
                }
            };
            evaluations += diagnostics.evaluations;
            assert_eq!(
                array.source_diagnostics(index).unwrap().evaluations,
                diagnostics.evaluations
            );
            expected[0] = expected[0] + line.pressure_pa;
            for i in 0..3 {
                expected[i + 1] = expected[i + 1] + line.particle_velocity_m_s[i];
            }
        }
        assert_eq!(array.evaluations(), evaluations);
        let actual = array.evaluate(&moments).unwrap();
        let line = array.evaluate_harmonic_line(&moments).unwrap();
        assert_eq!(line.frequency_hz, 700.0);
        assert_eq!(line.pressure_pa, actual.pressure_pa);
        near(actual.pressure_pa, expected[0]);
        let stats = array
            .evaluate_coherence_modes(&[&moments], coherence(), || false)
            .unwrap();
        near(
            C64::from_re(stats.mean_square_pressure_pa2),
            C64::from_re(0.5 * expected[0].abs().powi(2)),
        );
        for i in 0..3 {
            near(actual.particle_velocity_m_s[i], expected[i + 1]);
            let flux = expected[0].scale(0.5) * expected[i + 1].conj();
            near(
                C64::new(
                    actual.active_intensity_w_m2[i],
                    actual.reactive_intensity_w_m2[i],
                ),
                flux,
            );
            near(
                C64::new(
                    stats.active_intensity_w_m2[i],
                    stats.reactive_intensity_w_m2[i],
                ),
                flux,
            );
        }
        let pair = PreparedLayeredFullWaveQuadrupoleObserver::new(
            slab,
            kind,
            700.0,
            position,
            &[SOURCES[0]; 2],
            limits(),
            &mut || false,
        )
        .unwrap();
        let opposite = QuadrupolePhasor {
            stress_moment_n_m: moments[0].stress_moment_n_m.map(|m| m.scale(-1.0)),
        };
        assert_eq!(
            pair.evaluate(&[moments[0], opposite]).unwrap().pressure_pa,
            C64::ZERO
        );
        assert_eq!(
            pair.evaluate(&[moments[0], opposite])
                .unwrap()
                .particle_velocity_m_s,
            [C64::ZERO; 3]
        );
        let independent = pair
            .evaluate_independent(&[moments[0], opposite], || false)
            .unwrap();
        assert!(independent.mean_square_pressure_pa2 > 0.0);
        let a = [moments[0], QuadrupolePhasor::default()];
        let b = [QuadrupolePhasor::default(), opposite];
        let modes = pair
            .evaluate_coherence_modes(&[&a, &b], coherence(), || false)
            .unwrap();
        statistics_equal(&modes, &independent);
    }
}

#[test]
fn recorded_reynolds_stresses_keep_original_clock_geometry_and_native_bin_power() {
    let clock = QuadrupoleSampleClock::new(0.125, 1.0 / 8192.0, 0.0).unwrap();
    let config = WelchConfig {
        segment_samples: 64,
        hop_samples: 32,
        window: WelchWindow::PeriodicHann,
        detrend: WelchDetrend::SegmentMean,
    };
    let budget = QuadrupoleRecordingBudget {
        maximum_sources: 2,
        welch: WelchBudget {
            maximum_channels: 12,
            maximum_samples_per_channel: 96,
            maximum_transforms: 128,
            maximum_workspace_bytes: 1_000_000,
        },
    };
    let mut recording =
        QuadrupoleWelchStream::new(&SOURCES, clock, 96, config, budget, || false).unwrap();
    for sample in 0..96 {
        let t = core::f64::consts::TAU * sample as f64 / 64.0;
        let moments: [[f64; 6]; 2] = core::array::from_fn(|cell| {
            let phase = t + cell as f64 * 0.2;
            resolved_reynolds_moment_n_m(
                1.2,
                0.002,
                [
                    1.0 + (3.0 * phase).cos(),
                    (2.0 * phase).sin(),
                    (7.0 * phase).cos(),
                ],
            )
            .unwrap()
        });
        recording
            .push_sample(&SOURCES, clock.time_at(sample).unwrap(), &moments, || false)
            .unwrap();
    }
    let recorded = recording.finish().unwrap();
    let slab = layer(true);
    let mut exact_limits = limits();
    exact_limits.maximum_total_evaluations = 1;
    exact_limits.quadrupole.maximum_evaluations = 1;
    for kind in KINDS {
        let position = receiver(kind);
        let observed = recorded
            .prepare_layered_observer(slab, kind, position, 6, exact_limits, || false)
            .unwrap();
        assert_eq!(observed.evaluations(), 0);
        assert_eq!(observed.frequency_hz(), 768.0);
        assert_eq!(recorded.clock(), clock);
        let actual = recorded
            .evaluate_layered_bin(&observed, 6, coherence(), || false)
            .unwrap();
        if kind == LayerQuadrupoleField::Reflected {
            assert_eq!(actual.mean_square_pressure_pa2, 0.0);
            assert_eq!(actual.velocity_covariance_m2_s2, [[C64::ZERO; 3]; 3]);
        } else {
            let native = recorded
                .prepare_observer(slab.incident_state(), position, 6, || false)
                .unwrap();
            let expected = recorded
                .evaluate_bin(&native, 6, coherence(), || false)
                .unwrap();
            assert!(expected.mean_square_pressure_pa2 > 0.0);
            statistics_equal(&actual, &expected);
        }
        assert!(
            recorded
                .evaluate_layered_bin(&observed, 7, coherence(), || false)
                .is_err()
        );
        let swapped = PreparedLayeredFullWaveQuadrupoleObserver::new(
            slab,
            kind,
            768.0,
            position,
            &[SOURCES[1], SOURCES[0]],
            exact_limits,
            &mut || false,
        )
        .unwrap();
        assert!(
            recorded
                .evaluate_layered_bin(&swapped, 6, coherence(), || false)
                .is_err()
        );
        assert!(matches!(
            recorded.evaluate_layered_bin(&observed, 6, coherence(), || true),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    assert!(
        recorded
            .prepare_layered_observer(
                slab,
                LayerQuadrupoleField::Interior,
                receiver(LayerQuadrupoleField::Interior),
                0,
                exact_limits,
                || false
            )
            .is_err()
    );
}

#[test]
fn observer_domains_self_fields_and_shared_work_refuse_before_partial_publication() {
    let slab = layer(false);
    for kind in KINDS {
        let position = receiver(kind);
        let one = PreparedLayeredFullWaveQuadrupoleObserver::new(
            slab,
            kind,
            700.0,
            position,
            &SOURCES[..1],
            limits(),
            &mut || false,
        )
        .unwrap();
        assert!(one.evaluations() > 0);
        let mut capped = limits();
        capped.maximum_total_evaluations = 2 * one.evaluations() - 1;
        assert!(
            PreparedLayeredFullWaveQuadrupoleObserver::new(
                slab,
                kind,
                700.0,
                position,
                &[SOURCES[0]; 2],
                capped,
                &mut || false
            )
            .is_err()
        );
        capped = limits();
        capped.max_sources = 0;
        let empty = PreparedLayeredFullWaveQuadrupoleObserver::new(
            slab,
            kind,
            700.0,
            position,
            &[],
            capped,
            &mut || false,
        )
        .unwrap();
        assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
        assert!(
            PreparedLayeredFullWaveQuadrupoleObserver::new(
                slab,
                kind,
                700.0,
                [0.0; 3],
                &[],
                capped,
                &mut || false
            )
            .is_err()
        );
        assert!(
            PreparedLayeredFullWaveQuadrupoleObserver::new(
                slab,
                kind,
                0.0,
                position,
                &[],
                capped,
                &mut || false
            )
            .is_err()
        );
        assert!(
            PreparedLayeredFullWaveQuadrupoleObserver::new(
                slab,
                kind,
                700.0,
                position,
                &[[0.0, 0.0, 0.007]],
                limits(),
                &mut || false
            )
            .is_err()
        );
        assert!(
            PreparedLayeredFullWaveQuadrupoleObserver::new(
                slab,
                kind,
                700.0,
                position,
                &SOURCES[..1],
                capped,
                &mut || false
            )
            .is_err()
        );
        assert!(matches!(
            one.evaluate_cancellable(&[moment(0)], || true),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let returned = PreparedLayeredFullWaveQuadrupoleObserver::new(
        slab,
        LayerQuadrupoleField::Reflected,
        700.0,
        SOURCES[0],
        &SOURCES[..1],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        returned
            .evaluate(&[moment(0)])
            .unwrap()
            .pressure_pa
            .abs()
            .is_finite()
    );
    assert!(
        PreparedLayeredFullWaveQuadrupoleObserver::new(
            slab,
            LayerQuadrupoleField::IncidentTotal,
            700.0,
            SOURCES[0],
            &SOURCES[..1],
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayeredFullWaveQuadrupoleObserver::new(
            slab,
            LayerQuadrupoleField::Interior,
            700.0,
            receiver(LayerQuadrupoleField::Interior),
            &SOURCES,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}
