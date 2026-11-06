//! G1/G3/G4: original inside-layer stress arrays, exterior domains and recorded clocks.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::recorded::{
    QuadrupoleRecordingBudget, QuadrupoleSampleClock, QuadrupoleWelchStream,
};
use fs_aeroac::curle3d::quadrupole::{QuadrupolePhasor, resolved_reynolds_moment_n_m};
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::PreparedInteriorSourceFullWaveQuadrupole;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::quadrupole::array::PreparedInteriorSourceFullWaveQuadrupoleObserver;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::FullWaveQuadrupoleSourceLimits;
use fs_aeroac::curle3d::{CoherenceBudget, CompactFieldStatistics};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

const H: f64 = 0.04;
const SOURCES: [[f64; 3]; 2] = [[-0.01, 0.0, 0.012], [0.02, -0.015, 0.027]];
const OBSERVER: [f64; 3] = [0.035, 0.01, -0.03];

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn slab(homogeneous: bool) -> PlanarFluidLayer {
    let a = fluid(1.2, 300.0);
    PlanarFluidLayer::new(
        a,
        if homogeneous { a } else { fluid(2.4, 600.0) },
        if homogeneous { a } else { fluid(1.8, 400.0) },
        [0.0; 3],
        [0.0, 0.0, 1.0],
        H,
    )
    .unwrap()
}
fn reverse(layer: PlanarFluidLayer) -> PlanarFluidLayer {
    PlanarFluidLayer::new(
        layer.receiving_state(),
        layer.layer_state(),
        layer.incident_state(),
        [0.0, 0.0, H],
        [0.0, 0.0, -1.0],
        H,
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
fn both_exteriors_keep_native_stress_fields_phases_and_statistics() {
    let moments = [moment(0), moment(1)];
    let original = slab(false);
    for (layer, observer) in [
        (original, OBSERVER),
        (reverse(original), [0.035, 0.01, H + 0.03]),
    ] {
        let array = PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            layer,
            700.0,
            observer,
            &SOURCES,
            limits(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(array.source_count(), 2);
        assert_eq!(array.source_position_m(0), Some(SOURCES[0]));
        assert_eq!(array.source_position_m(1), Some(SOURCES[1]));
        assert_eq!(array.source_position_m(2), None);
        assert_eq!(array.observer_position_m(), observer);
        let mut expected = [C64::ZERO; 4];
        let mut work = 0;
        for (index, (&source, &stress)) in SOURCES.iter().zip(&moments).enumerate() {
            let native = PreparedInteriorSourceFullWaveQuadrupole::new(
                layer,
                700.0,
                source,
                observer,
                limits().quadrupole,
                &mut || false,
            )
            .unwrap();
            let line = native.evaluate_harmonic_line(stress).unwrap();
            work += native.diagnostics().evaluations;
            assert_eq!(
                array.source_diagnostics(index).unwrap().evaluations,
                native.diagnostics().evaluations
            );
            expected[0] = expected[0] + line.pressure_pa;
            for i in 0..3 {
                expected[i + 1] = expected[i + 1] + line.particle_velocity_m_s[i];
            }
        }
        assert_eq!(array.evaluations(), work);
        let result = array.evaluate(&moments).unwrap();
        let line = array.evaluate_harmonic_line(&moments).unwrap();
        assert_eq!(line.frequency_hz, 700.0);
        assert_eq!(line.pressure_pa, result.pressure_pa);
        near(result.pressure_pa, expected[0]);
        let statistics = array
            .evaluate_coherence_modes(&[&moments], coherence(), || false)
            .unwrap();
        near(
            C64::from_re(statistics.mean_square_pressure_pa2),
            C64::from_re(0.5 * expected[0].abs().powi(2)),
        );
        for i in 0..3 {
            near(result.particle_velocity_m_s[i], expected[i + 1]);
            let flux = expected[0].scale(0.5) * expected[i + 1].conj();
            near(
                C64::new(
                    result.active_intensity_w_m2[i],
                    result.reactive_intensity_w_m2[i],
                ),
                flux,
            );
            near(
                C64::new(
                    statistics.active_intensity_w_m2[i],
                    statistics.reactive_intensity_w_m2[i],
                ),
                flux,
            );
        }
        let first = [moments[0], QuadrupolePhasor::default()];
        let second = [QuadrupolePhasor::default(), moments[1]];
        statistics_equal(
            &array.evaluate_independent(&moments, || false).unwrap(),
            &array
                .evaluate_coherence_modes(&[&first, &second], coherence(), || false)
                .unwrap(),
        );
    }
}

#[test]
fn cell_cancellation_is_coherent_but_independent_cells_keep_power() {
    let array = PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
        slab(false),
        700.0,
        OBSERVER,
        &[SOURCES[0]; 2],
        limits(),
        &mut || false,
    )
    .unwrap();
    let a = moment(0);
    let b = QuadrupolePhasor {
        stress_moment_n_m: a.stress_moment_n_m.map(|m| m.scale(-1.0)),
    };
    let field = array.evaluate(&[a, b]).unwrap();
    assert_eq!(field.pressure_pa, C64::ZERO);
    assert_eq!(field.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(
        array
            .evaluate_coherence_modes(&[&[a, b]], coherence(), || false)
            .unwrap()
            .mean_square_pressure_pa2,
        0.0
    );
    assert!(
        array
            .evaluate_independent(&[a, b], || false)
            .unwrap()
            .mean_square_pressure_pa2
            > 0.0
    );
}

#[test]
fn original_inside_source_domain_and_shared_work_apply_even_to_empty_arrays() {
    let layer = slab(false);
    let one = PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
        layer,
        700.0,
        OBSERVER,
        &SOURCES[..1],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(one.evaluations() > 0);
    let mut cap = limits();
    cap.maximum_total_evaluations = 2 * one.evaluations() - 1;
    assert!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            layer,
            700.0,
            OBSERVER,
            &[SOURCES[0]; 2],
            cap,
            &mut || false,
        )
        .is_err()
    );
    cap = limits();
    cap.max_sources = 0;
    let empty = PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
        layer,
        700.0,
        OBSERVER,
        &[],
        cap,
        &mut || false,
    )
    .unwrap();
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
    for observer in [
        [0.0; 3],
        [0.0, 0.0, 0.02],
        [0.0, 0.0, 0.08],
        [f64::NAN, 0.0, -0.03],
    ] {
        assert!(
            PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
                layer,
                700.0,
                observer,
                &[],
                cap,
                &mut || false,
            )
            .is_err()
        );
    }
    for source in [
        [0.0; 3],
        [0.0, 0.0, H],
        [0.0, 0.0, -0.01],
        [0.0, 0.0, H + 0.01],
        [f64::NAN, 0.0, 0.02],
    ] {
        assert!(
            PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
                layer,
                700.0,
                OBSERVER,
                &[source],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    let zero = PlanarFluidLayer::new(
        layer.incident_state(),
        layer.layer_state(),
        layer.receiving_state(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        0.0,
    )
    .unwrap();
    assert!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            zero,
            700.0,
            OBSERVER,
            &[],
            cap,
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            layer,
            0.0,
            OBSERVER,
            &[],
            cap,
            &mut || false
        )
        .is_err()
    );
    cap.maximum_total_evaluations = 0;
    assert!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            layer,
            700.0,
            OBSERVER,
            &[],
            cap,
            &mut || false
        )
        .is_err()
    );
    cap = limits();
    cap.quadrupole.relative_tolerance = f64::NAN;
    assert!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            layer,
            700.0,
            OBSERVER,
            &[],
            cap,
            &mut || false
        )
        .is_err()
    );
    assert!(one.evaluate(&[]).is_err());
    assert!(
        one.evaluate_coherence_modes(
            &[&[moment(0)]],
            CoherenceBudget {
                max_modes: 0,
                max_source_mode_pairs: 0
            },
            || false
        )
        .is_err()
    );
    assert!(matches!(
        one.evaluate_cancellable(&[moment(0)], || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            layer,
            700.0,
            OBSERVER,
            &[],
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));

    let a = fluid(1.2, 300.0);
    let b = fluid(2.4, 300.0);
    let exact = PlanarFluidLayer::new(a, b, b, [0.0; 3], [0.0, 0.0, 1.0], H).unwrap();
    let density_layer = PlanarFluidLayer::new(a, b, a, [0.0; 3], [0.0, 0.0, 1.0], H).unwrap();
    cap = limits();
    cap.maximum_total_evaluations = 1;
    cap.quadrupole.maximum_evaluations = 1;
    assert_eq!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            exact,
            700.0,
            OBSERVER,
            &SOURCES,
            cap,
            &mut || false
        )
        .unwrap()
        .evaluations(),
        0
    );
    assert!(
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            density_layer,
            700.0,
            OBSERVER,
            &SOURCES[..1],
            cap,
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn recorded_stresses_preserve_original_clock_geometry_and_joint_bin_power() {
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
    let mut stream =
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
        stream
            .push_sample(&SOURCES, clock.time_at(sample).unwrap(), &moments, || false)
            .unwrap();
    }
    let recorded = stream.finish().unwrap();
    let layer = slab(true);
    let mut exact_limits = limits();
    exact_limits.maximum_total_evaluations = 1;
    exact_limits.quadrupole.maximum_evaluations = 1;
    let observer = recorded
        .prepare_interior_source_layered_observer(layer, OBSERVER, 6, exact_limits, || false)
        .unwrap();
    assert_eq!(observer.evaluations(), 0);
    assert_eq!(observer.frequency_hz(), 768.0);
    assert_eq!(recorded.clock(), clock);
    let actual = recorded
        .evaluate_interior_source_layered_bin(&observer, 6, coherence(), || false)
        .unwrap();
    let native = recorded
        .prepare_observer(layer.layer_state(), OBSERVER, 6, || false)
        .unwrap();
    let expected = recorded
        .evaluate_bin(&native, 6, coherence(), || false)
        .unwrap();
    assert!(expected.mean_square_pressure_pa2 > 0.0);
    statistics_equal(&actual, &expected);
    assert!(
        recorded
            .evaluate_interior_source_layered_bin(&observer, 7, coherence(), || false)
            .is_err()
    );
    assert!(
        recorded
            .prepare_interior_source_layered_observer(layer, OBSERVER, 0, exact_limits, || false)
            .is_err()
    );
    let swapped = PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
        layer,
        768.0,
        OBSERVER,
        &[SOURCES[1], SOURCES[0]],
        exact_limits,
        &mut || false,
    )
    .unwrap();
    assert!(
        recorded
            .evaluate_interior_source_layered_bin(&swapped, 6, coherence(), || false)
            .is_err()
    );
    assert!(matches!(
        recorded.evaluate_interior_source_layered_bin(&observer, 6, coherence(), || true),
        Err(AeroacError::Cancelled { .. })
    ));
}
