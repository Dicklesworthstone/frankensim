//! Original Q/F/stress recordings through coherent finite-layer reflections.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::refraction::full_wave::array::{
    FullWaveCompactFieldKind as Kind, FullWaveCompactSourceLimits,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::array::{
    PreparedLayeredFullWaveCompactSourceObserver as CompactObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    LayerQuadrupoleField as StressKind, PlanarFluidLayer,
    PreparedLayeredFullWaveQuadrupoleObserver as StressObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::FullWaveQuadrupoleSourceLimits;
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
    PreparedCompactSourceObserver,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

const N: usize = 64;
const BIN: usize = 5;
const HZ: f64 = 640.0;
const QF_POSITIONS: [[f64; 3]; 2] = [[-0.02, 0.01, -0.16], [0.05, -0.03, -0.13]];
const STRESS_POSITIONS: [[f64; 3]; 1] = [[0.01, -0.015, -0.15]];
const ZERO: CompactSourcePhasor = CompactSourcePhasor {
    volume_flow_m3_s: C64::ZERO, force_on_fluid_n: [C64::ZERO; 3],
};
const WORK: CoherenceBudget = CoherenceBudget {
    max_modes: 3, max_source_mode_pairs: 9,
};

fn layer(rho: [f64; 3], speed: [f64; 3], point: [f64; 3], normal: [f64; 3], h: f64)
    -> PlanarFluidLayer
{
    let fluid = |i: usize| FluidAcousticState::try_new(rho[i], speed[i]).unwrap();
    PlanarFluidLayer::new(fluid(0), fluid(1), fluid(2), point, normal, h).unwrap()
}
fn uniform(h: f64) -> PlanarFluidLayer {
    layer([1.2; 3], [343.0; 3], [0.0; 3], [0.0, 0.0, 1.0], h)
}
fn contrast() -> PlanarFluidLayer {
    layer([1.2, 2.1, 1.6], [343.0, 420.0, 365.0], [0.0; 3], [0.0, 0.0, 1.0], 0.05)
}
fn compact_limits() -> FullWaveCompactSourceLimits {
    FullWaveCompactSourceLimits {
        max_sources: 2, maximum_total_evaluations: 1_048_576,
        monopole: FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
            absolute_velocity_tolerance_per_m2: 1e-10,
            relative_tolerance: 1e-9, maximum_evaluations: 131_072, maximum_depth: 40,
        },
        force: FullWaveForceLimits {
            absolute_pressure_tolerance_pa_per_n: 1e-9,
            absolute_velocity_tolerance_m_s_per_n: 1e-12,
            relative_tolerance: 1e-9, maximum_evaluations: 131_072, maximum_depth: 40,
        },
    }
}
fn stress_limits() -> FullWaveQuadrupoleSourceLimits {
    FullWaveQuadrupoleSourceLimits {
        max_sources: 1, maximum_total_evaluations: 262_144,
        quadrupole: FullWaveQuadrupoleLimits {
            absolute_pressure_tolerance_pa_per_n_m: 1e-8,
            absolute_velocity_tolerance_m_s_per_n_m: 1e-11,
            relative_tolerance: 1e-10, maximum_evaluations: 262_144, maximum_depth: 40,
        },
    }
}
fn receiver(kind: Kind) -> [f64; 3] {
    [0.03, 0.02, if kind == Kind::Transmitted { 0.13 } else { -0.07 }]
}
fn observers(boundary: PlanarFluidLayer, kind: Kind) -> (CompactObserver, StressObserver) {
    let compact = match kind {
        Kind::Transmitted => CompactObserver::new(
            boundary, HZ, receiver(kind), &QF_POSITIONS, compact_limits(), &mut || false,
        ),
        Kind::Reflected => CompactObserver::new_reflected(
            boundary, HZ, receiver(kind), &QF_POSITIONS, compact_limits(), &mut || false,
        ),
        Kind::IncidentTotal => CompactObserver::new_total(
            boundary, HZ, receiver(kind), &QF_POSITIONS, compact_limits(), &mut || false,
        ),
    }.unwrap();
    let stress_kind = match kind {
        Kind::Transmitted => StressKind::Transmitted,
        Kind::Reflected => StressKind::Reflected,
        Kind::IncidentTotal => StressKind::IncidentTotal,
    };
    let stress = StressObserver::new(
        boundary, stress_kind, HZ, receiver(kind), &STRESS_POSITIONS,
        stress_limits(), &mut || false,
    ).unwrap();
    (compact, stress)
}
fn estimate<const D: usize>(segments: &[[C64; D]]) -> WelchCrossSpectrum {
    let record: Vec<_> = segments.iter().flat_map(|amplitudes| (0..N).map(move |j| {
        let phase = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
        amplitudes.map(|a| a.re * phase.cos() + a.im * phase.sin())
    })).collect();
    WelchCrossSpectrum::estimate_interleaved(
        &record, 1.0 / 8192.0,
        WelchConfig {
            segment_samples: N, hop_samples: N,
            window: WelchWindow::Rectangular, detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 14, maximum_samples_per_channel: 3 * N,
            maximum_transforms: 42, maximum_workspace_bytes: 1_048_576,
        },
        || false,
    ).unwrap()
}
fn sources() -> [([CompactSourcePhasor; 2], QuadrupolePhasor); 3] {
    core::array::from_fn(|segment| {
        let phases = [C64::ONE, C64::new(0.0, 1.0), C64::new(-0.4, 0.3)];
        let compact = core::array::from_fn(|cell| {
            let phase = phases[(segment + cell) % 3];
            CompactSourcePhasor {
                volume_flow_m3_s: C64::new(2e-5, -1e-5) * phase,
                force_on_fluid_n: [
                    C64::new(0.012, -0.006) * phase, C64::new(-0.003, 0.008) * phase,
                    C64::new(0.002, 0.005) * phase,
                ],
            }
        });
        let stress = QuadrupolePhasor {
            stress_moment_n_m: core::array::from_fn(|entry| {
                C64::new(2e-4 * (entry + 1) as f64, 3e-4 * (3.0 - entry as f64))
                    * phases[segment].conj()
            }),
        };
        (compact, stress)
    })
}
fn channels(compact: [CompactSourcePhasor; 2], stress: QuadrupolePhasor) -> [C64; 14] {
    core::array::from_fn(|i| {
        if i < 8 {
            let source = compact[i / 4];
            if i % 4 == 0 { source.volume_flow_m3_s } else { source.force_on_fluid_n[i % 4 - 1] }
        } else {
            stress.stress_moment_n_m[i - 8]
        }
    })
}
fn values(field: CoherentCompactField) -> [C64; 4] {
    [field.pressure_pa, field.particle_velocity_m_s[0],
        field.particle_velocity_m_s[1], field.particle_velocity_m_s[2]]
}
fn add(a: CoherentCompactField, b: CoherentCompactField) -> [C64; 4] {
    let a = values(a);
    let b = values(b);
    core::array::from_fn(|i| a[i] + b[i])
}
fn close(actual: C64, expected: C64) {
    assert!((actual - expected).abs() < 1e-23 + 8e-10 * expected.abs(),
        "actual={actual:?}, expected={expected:?}");
}
fn moments(actual: &CompactFieldStatistics, fields: &[[C64; 4]]) {
    let moment = |i: usize, j: usize| fields.iter()
        .fold(C64::ZERO, |sum, f| sum + f[i] * f[j].conj())
        .scale(0.5 / fields.len() as f64);
    close(C64::from_re(actual.mean_square_pressure_pa2), moment(0, 0));
    for i in 0..3 {
        close(C64::new(actual.active_intensity_w_m2[i], actual.reactive_intensity_w_m2[i]),
            moment(0, i + 1));
        for j in 0..3 { close(actual.velocity_covariance_m2_s2[i][j], moment(i + 1, j + 1)); }
    }
}

#[test]
fn original_segment_phase_and_full_tensor_moments_survive_all_layer_views() {
    let input = sources();
    let joint = estimate(&input.map(|(qf, m)| channels(qf, m)));
    let compact_only = estimate(&input.map(|(qf, m)| {
        let all = channels(qf, m);
        core::array::from_fn::<_, 8, _>(|i| all[i])
    }));
    for kind in [Kind::Transmitted, Kind::Reflected, Kind::IncidentTotal] {
        let (compact, stress) = observers(contrast(), kind);
        assert!(compact.evaluations() > 0);
        assert!(stress.evaluations() > 0);
        let qf_fields: Vec<_> = input.iter().map(|(qf, _)| {
            values(compact.evaluate(qf).unwrap())
        }).collect();
        let mixed_fields: Vec<_> = input.iter().map(|(qf, m)| {
            add(compact.evaluate(qf).unwrap(), stress.evaluate(&[*m]).unwrap())
        }).collect();
        moments(&compact.evaluate_welch_bin(&compact_only, BIN, WORK, || false).unwrap(),
            &qf_fields);
        moments(&compact.evaluate_welch_bin_with_quadrupoles(
            &stress, &joint, BIN, WORK, || false,
        ).unwrap(), &mixed_fields);
    }
}

#[test]
fn absent_and_homogeneous_layers_keep_the_original_free_field_record() {
    let input = sources();
    let joint = estimate(&input.map(|(qf, m)| channels(qf, m)));
    for h in [0.0, 0.08] {
        for kind in [Kind::Transmitted, Kind::Reflected, Kind::IncidentTotal] {
            let boundary = uniform(h);
            let (compact, stress) = observers(boundary, kind);
            let free_compact = PreparedCompactSourceObserver::new(
                boundary.incident_state(), HZ, receiver(kind), &QF_POSITIONS, 2,
            ).unwrap();
            let free_stress = PreparedQuadrupoleObserver::new(
                boundary.incident_state(), HZ, receiver(kind), &STRESS_POSITIONS, 1,
            ).unwrap();
            let expected: Vec<_> = input.iter().map(|(qf, m)| {
                if kind == Kind::Reflected { [C64::ZERO; 4] } else {
                    add(free_compact.evaluate(qf).unwrap(), free_stress.evaluate(&[*m]).unwrap())
                }
            }).collect();
            assert_eq!(compact.evaluations(), 0);
            assert_eq!(stress.evaluations(), 0);
            moments(&compact.evaluate_welch_bin_with_quadrupoles(
                &stress, &joint, BIN, WORK, || false,
            ).unwrap(), &expected);
        }
    }
}

#[test]
fn measured_cross_mechanism_phase_preserves_pressure_null_velocity() {
    let (compact, stress) = observers(contrast(), Kind::IncidentTotal);
    let mut m = QuadrupolePhasor::default();
    m.stress_moment_n_m[2] = C64::new(0.001, -0.0002);
    let returned = stress.evaluate(&[m]).unwrap();
    let unit = compact.evaluate(&[CompactSourcePhasor::monopole(C64::ONE), ZERO]).unwrap();
    let qf = [CompactSourcePhasor::monopole(-returned.pressure_pa / unit.pressure_pa), ZERO];
    let q = compact.evaluate(&qf).unwrap();
    let independent_power = 0.5 * (q.pressure_pa.abs().powi(2) + returned.pressure_pa.abs().powi(2));
    assert!(independent_power > 0.0);
    let coherent = estimate(&[channels(qf, m)]);
    let actual = compact.evaluate_welch_bin_with_quadrupoles(
        &stress, &coherent, BIN, WORK, || false,
    ).unwrap();
    assert!(actual.mean_square_pressure_pa2 < 1e-20 * independent_power);
    assert!((0..3).map(|i| actual.velocity_covariance_m2_s2[i][i].re).sum::<f64>() > 0.0);
    let opposite = QuadrupolePhasor {
        stress_moment_n_m: m.stress_moment_n_m.map(|value| -value),
    };
    let uncorrelated = estimate(&[channels(qf, m), channels(qf, opposite)]);
    let actual = compact.evaluate_welch_bin_with_quadrupoles(
        &stress, &uncorrelated, BIN, WORK, || false,
    ).unwrap();
    close(C64::from_re(actual.mean_square_pressure_pa2), C64::from_re(independent_power));
}

#[test]
fn all_three_fluids_frame_thickness_clock_and_total_work_are_bound() {
    let input = sources();
    let joint = estimate(&input.map(|(qf, m)| channels(qf, m)));
    let (compact, stress) = observers(uniform(0.0), Kind::IncidentTotal);
    let original_receiver = receiver(Kind::IncidentTotal);
    let changed_layers = [
        layer([1.3, 1.2, 1.2], [343.0; 3], [0.0; 3], [0.0, 0.0, 1.0], 0.0),
        layer([1.2, 1.3, 1.2], [343.0; 3], [0.0; 3], [0.0, 0.0, 1.0], 0.0),
        layer([1.2, 1.2, 1.3], [343.0; 3], [0.0; 3], [0.0, 0.0, 1.0], 0.0),
        layer([1.2; 3], [343.0; 3], [0.0, 0.0, 0.02], [0.0, 0.0, 1.0], 0.0),
        layer([1.2; 3], [343.0; 3], [0.0; 3], [0.1, 0.0, 1.0], 0.0),
        uniform(0.05),
    ];
    for boundary in changed_layers {
        let other = StressObserver::new(
            boundary, StressKind::IncidentTotal, HZ, original_receiver, &STRESS_POSITIONS,
            stress_limits(), &mut || false,
        ).unwrap();
        assert!(compact.evaluate_welch_bin_with_quadrupoles(
            &other, &joint, BIN, WORK, || false,
        ).is_err());
    }
    for (point, hz, kind) in [
        ([0.04, 0.02, -0.07], HZ, StressKind::IncidentTotal),
        (original_receiver, HZ + 128.0, StressKind::IncidentTotal),
        (original_receiver, HZ, StressKind::Reflected),
    ] {
        let other = StressObserver::new(
            uniform(0.0), kind, hz, point, &STRESS_POSITIONS, stress_limits(), &mut || false,
        ).unwrap();
        assert!(compact.evaluate_welch_bin_with_quadrupoles(
            &other, &joint, BIN, WORK, || false,
        ).is_err());
    }
    assert!(compact.evaluate_welch_bin_with_quadrupoles(
        &stress, &joint, BIN + 1, WORK, || false,
    ).is_err());
    for work in [
        CoherenceBudget { max_modes: 3, max_source_mode_pairs: 6 },
        CoherenceBudget { max_modes: 2, max_source_mode_pairs: 9 },
    ] {
        assert!(compact.evaluate_welch_bin_with_quadrupoles(
            &stress, &joint, BIN, work, || false,
        ).is_err());
    }
    let compact_only = estimate(&[core::array::from_fn::<_, 8, _>(|i| {
        channels(input[0].0, input[0].1)[i]
    })]);
    assert!(compact.evaluate_welch_bin_with_quadrupoles(
        &stress, &compact_only, BIN, WORK, || false,
    ).is_err());
    assert!(compact.evaluate_welch_bin(&joint, BIN, WORK, || false).is_err());
    assert!(matches!(
        compact.evaluate_welch_bin(&compact_only, BIN, WORK, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        compact.evaluate_welch_bin_with_quadrupoles(&stress, &joint, BIN, WORK, || true),
        Err(AeroacError::Cancelled { .. })
    ));
}
