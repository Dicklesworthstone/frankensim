//! G1/G3/G4: original joint compact/stress recordings through a refracting plane.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::QuadrupolePhasor;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::array::{
    FullWaveCompactFieldKind, FullWaveCompactSourceLimits,
    PreparedFullWaveCompactSourceObserver as CompactObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::{
    FullWaveQuadrupoleSourceLimits, PreparedFullWaveQuadrupoleObserver as StressObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const N: usize = 64;
const BIN: usize = 7;
const DT: f64 = 1.0 / 8192.0;
const FREQUENCY: f64 = BIN as f64 / (N as f64 * DT);
const COMPACT_POSITIONS: [[f64; 3]; 2] = [[-0.02, 0.01, -0.08], [0.035, -0.02, -0.06]];
const STRESS_POSITIONS: [[f64; 3]; 1] = [[0.01, -0.015, -0.07]];
const RECEIVER: [f64; 3] = [0.03, 0.02, 0.13];
const WORK: CoherenceBudget = CoherenceBudget {
    max_modes: 3,
    max_source_mode_pairs: 9,
};

fn boundary(
    rho1: f64,
    rho2: f64,
    c2: f64,
    point: [f64; 3],
    normal: [f64; 3],
) -> PlanarFluidInterface {
    PlanarFluidInterface::new(
        FluidAcousticState::try_new(rho1, 340.0).unwrap(),
        FluidAcousticState::try_new(rho2, c2).unwrap(),
        point,
        normal,
    )
    .unwrap()
}

fn interface(c2: f64) -> PlanarFluidInterface {
    boundary(1.2, 2.4, c2, [0.0; 3], [0.0, 0.0, 1.0])
}

fn compact_limits() -> FullWaveCompactSourceLimits {
    FullWaveCompactSourceLimits {
        max_sources: 2,
        maximum_total_evaluations: 1_048_576,
        monopole: FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-6,
            absolute_velocity_tolerance_per_m2: 1e-9,
            relative_tolerance: 1e-8,
            maximum_evaluations: 131_072,
            maximum_depth: 40,
        },
        force: FullWaveForceLimits {
            absolute_pressure_tolerance_pa_per_n: 1e-7,
            absolute_velocity_tolerance_m_s_per_n: 1e-10,
            relative_tolerance: 1e-8,
            maximum_evaluations: 131_072,
            maximum_depth: 40,
        },
    }
}

fn stress_limits() -> FullWaveQuadrupoleSourceLimits {
    FullWaveQuadrupoleSourceLimits {
        max_sources: 1,
        maximum_total_evaluations: 262_144,
        quadrupole: FullWaveQuadrupoleLimits {
            absolute_pressure_tolerance_pa_per_n_m: 1e-6,
            absolute_velocity_tolerance_m_s_per_n_m: 1e-9,
            relative_tolerance: 1e-8,
            maximum_evaluations: 262_144,
            maximum_depth: 40,
        },
    }
}

fn observers(kind: FullWaveCompactFieldKind, c2: f64) -> (CompactObserver, StressObserver) {
    let plane = interface(c2);
    let mut receiver = RECEIVER;
    if kind != FullWaveCompactFieldKind::Transmitted {
        receiver[2] = -receiver[2];
    }
    let compact = match kind {
        FullWaveCompactFieldKind::Transmitted => CompactObserver::new(
            plane,
            FREQUENCY,
            receiver,
            &COMPACT_POSITIONS,
            compact_limits(),
            &mut || false,
        ),
        FullWaveCompactFieldKind::Reflected => CompactObserver::new_reflected(
            plane,
            FREQUENCY,
            receiver,
            &COMPACT_POSITIONS,
            compact_limits(),
            &mut || false,
        ),
        FullWaveCompactFieldKind::IncidentTotal => CompactObserver::new_total(
            plane,
            FREQUENCY,
            receiver,
            &COMPACT_POSITIONS,
            compact_limits(),
            &mut || false,
        ),
    }
    .unwrap();
    let stress = match kind {
        FullWaveCompactFieldKind::Transmitted => StressObserver::new(
            plane,
            FREQUENCY,
            receiver,
            &STRESS_POSITIONS,
            stress_limits(),
            &mut || false,
        ),
        FullWaveCompactFieldKind::Reflected => StressObserver::new_reflected(
            plane,
            FREQUENCY,
            receiver,
            &STRESS_POSITIONS,
            stress_limits(),
            &mut || false,
        ),
        FullWaveCompactFieldKind::IncidentTotal => StressObserver::new_total(
            plane,
            FREQUENCY,
            receiver,
            &STRESS_POSITIONS,
            stress_limits(),
            &mut || false,
        ),
    }
    .unwrap();
    (compact, stress)
}

fn estimate<const D: usize>(segments: &[[C64; D]]) -> WelchCrossSpectrum {
    let record: Vec<_> = segments
        .iter()
        .flat_map(|amplitudes| {
            (0..N).map(move |j| {
                let theta = TAU * BIN as f64 * j as f64 / N as f64;
                amplitudes.map(|a| a.re * theta.cos() + a.im * theta.sin())
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &record,
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 14,
            maximum_samples_per_channel: 3 * N,
            maximum_transforms: 42,
            maximum_workspace_bytes: 1024 * 1024,
        },
        || false,
    )
    .unwrap()
}

fn sources() -> [([CompactSourcePhasor; 2], QuadrupolePhasor); 3] {
    core::array::from_fn(|segment| {
        let phases = [C64::ONE, C64::new(0.0, 1.0), C64::new(-0.4, 0.3)];
        let compact = core::array::from_fn(|source| {
            let phase = phases[(segment + source) % 3];
            CompactSourcePhasor {
                volume_flow_m3_s: C64::new(2e-5, -1e-5) * phase,
                force_on_fluid_n: [
                    C64::new(0.012, -0.006) * phase,
                    C64::new(-0.003, 0.008) * phase,
                    C64::new(0.002, 0.005) * phase,
                ],
            }
        });
        let stress = QuadrupolePhasor {
            stress_moment_n_m: core::array::from_fn(|i| {
                C64::new((i + 1) as f64 * 2e-4, (3.0 - i as f64) * 3e-4) * phases[segment].conj()
            }),
        };
        (compact, stress)
    })
}

fn channels(input: ([CompactSourcePhasor; 2], QuadrupolePhasor)) -> [C64; 14] {
    core::array::from_fn(|i| {
        if i < 8 {
            let source = input.0[i / 4];
            if i % 4 == 0 {
                source.volume_flow_m3_s
            } else {
                source.force_on_fluid_n[i % 4 - 1]
            }
        } else {
            input.1.stress_moment_n_m[i - 8]
        }
    })
}

fn field(value: CoherentCompactField) -> [C64; 4] {
    [
        value.pressure_pa,
        value.particle_velocity_m_s[0],
        value.particle_velocity_m_s[1],
        value.particle_velocity_m_s[2],
    ]
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 6e-11 * expected.abs() + 1e-24,
        "actual={actual:?}, expected={expected:?}"
    );
}

fn check_moments(actual: &CompactFieldStatistics, fields: &[[C64; 4]]) {
    let moment = |i: usize, j: usize| {
        fields
            .iter()
            .fold(C64::ZERO, |sum, value| sum + value[i] * value[j].conj())
            .scale(0.5 / fields.len() as f64)
    };
    close(C64::from_re(actual.mean_square_pressure_pa2), moment(0, 0));
    for i in 0..3 {
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            moment(0, i + 1),
        );
        for j in 0..3 {
            close(actual.velocity_covariance_m2_s2[i][j], moment(i + 1, j + 1));
        }
    }
}

#[test]
fn g1_original_joint_factors_match_known_phasors_in_all_three_refracting_fields() {
    let inputs = sources();
    let spectrum = estimate(&inputs.map(channels));
    let compact_record =
        estimate(&inputs.map(|input| core::array::from_fn::<_, 8, _>(|i| channels(input)[i])));
    let stress_record = estimate(&inputs.map(|input| input.1.stress_moment_n_m));
    for c2 in [340.0, 510.0] {
        for kind in [
            FullWaveCompactFieldKind::Transmitted,
            FullWaveCompactFieldKind::Reflected,
            FullWaveCompactFieldKind::IncidentTotal,
        ] {
            let (compact, stress) = observers(kind, c2);
            let compact_fields = inputs.map(|(c, _)| field(compact.evaluate(&c).unwrap()));
            let stress_fields = inputs.map(|(_, q)| field(stress.evaluate(&[q]).unwrap()));
            let expected: [[C64; 4]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|j| compact_fields[i][j] + stress_fields[i][j])
            });
            let result = compact
                .evaluate_welch_bin_with_quadrupoles(&stress, &spectrum, BIN, WORK, || false)
                .unwrap();
            check_moments(&result, &expected);
            assert!(
                result
                    .scope
                    .contains("without permeable-surface double counting")
            );
            let pure_compact = compact
                .evaluate_welch_bin(&compact_record, BIN, WORK, || false)
                .unwrap();
            check_moments(&pure_compact, &compact_fields);
            let pure_stress = stress
                .evaluate_welch_bin(&stress_record, BIN, WORK, || false)
                .unwrap();
            check_moments(&pure_stress, &stress_fields);
            let repeated = compact
                .evaluate_welch_bin_with_quadrupoles(&stress, &spectrum, BIN, WORK, || false)
                .unwrap();
            assert_eq!(
                result.mean_square_pressure_pa2,
                repeated.mean_square_pressure_pa2
            );
            assert_eq!(
                result.velocity_covariance_m2_s2,
                repeated.velocity_covariance_m2_s2
            );
            assert_eq!(compact.evaluations() == 0, c2 == 340.0);
            assert_eq!(stress.evaluations() == 0, c2 == 340.0);
        }
    }
}

#[test]
fn g3_transmitted_cross_mechanism_null_keeps_velocity_and_measured_correlation() {
    let plane = interface(340.0);
    let source = [[0.0, 0.0, -0.08]];
    let receiver = [0.0, 0.0, 0.13];
    let compact = CompactObserver::new(
        plane,
        FREQUENCY,
        receiver,
        &source,
        compact_limits(),
        &mut || false,
    )
    .unwrap();
    let stress = StressObserver::new(
        plane,
        FREQUENCY,
        receiver,
        &source,
        stress_limits(),
        &mut || false,
    )
    .unwrap();
    let flow = C64::from_re(2e-5);
    let omega = TAU * FREQUENCY;
    let k = omega / plane.incident_state().sound_speed_m_s();
    let r = receiver[2] - source[0][2];
    // Exact equal-speed transmission multiplies Q and zz pressure by the same
    // factor. Cancel their analytic spherical pressures; their gradients differ.
    let zz = C64::new(0.0, omega * plane.incident_state().density_kg_m3()) * flow
        / C64::new(2.0 / (r * r) - k * k, -2.0 * k / r);
    let mut cancelling = [C64::ZERO; 10];
    cancelling[0] = flow;
    cancelling[6] = zz;
    let mut reinforcing = cancelling;
    reinforcing[6] = -zz;
    let single = compact
        .evaluate(&[CompactSourcePhasor::monopole(flow)])
        .unwrap();
    let single_power = 0.5 * single.pressure_pa.norm_sq();
    let quiet = compact
        .evaluate_welch_bin_with_quadrupoles(
            &stress,
            &estimate(&[cancelling, cancelling]),
            BIN,
            WORK,
            || false,
        )
        .unwrap();
    assert!(quiet.mean_square_pressure_pa2 < single_power * 1e-24);
    assert!(
        quiet.velocity_covariance_m2_s2[2][2].re > 1e-5 * single.particle_velocity_m_s[2].norm_sq()
    );
    let decorrelated = compact
        .evaluate_welch_bin_with_quadrupoles(
            &stress,
            &estimate(&[cancelling, reinforcing]),
            BIN,
            WORK,
            || false,
        )
        .unwrap();
    close(
        C64::from_re(decorrelated.mean_square_pressure_pa2),
        C64::from_re(2.0 * single_power),
    );
    let loud = compact
        .evaluate_welch_bin_with_quadrupoles(
            &stress,
            &estimate(&[reinforcing, reinforcing]),
            BIN,
            WORK,
            || false,
        )
        .unwrap();
    close(
        C64::from_re(loud.mean_square_pressure_pa2),
        C64::from_re(4.0 * single_power),
    );
}

#[test]
fn g4_joint_domains_clock_shape_combined_work_and_cancellation_are_admitted_first() {
    let spectrum = estimate(&sources().map(channels));
    let (compact, stress) = observers(FullWaveCompactFieldKind::Transmitted, 340.0);
    let observe = |q: &StressObserver, bin, budget| {
        compact.evaluate_welch_bin_with_quadrupoles(q, &spectrum, bin, budget, || false)
    };
    for (plane, frequency, receiver) in [
        (
            boundary(1.3, 2.4, 340.0, [0.0; 3], [0.0, 0.0, 1.0]),
            FREQUENCY,
            RECEIVER,
        ),
        (
            boundary(1.2, 2.5, 340.0, [0.0; 3], [0.0, 0.0, 1.0]),
            FREQUENCY,
            RECEIVER,
        ),
        (
            boundary(1.2, 2.4, 340.0, [0.0, 0.0, 0.001], [0.0, 0.0, 1.0]),
            FREQUENCY,
            RECEIVER,
        ),
        (
            boundary(1.2, 2.4, 340.0, [0.0; 3], [0.05, 0.0, 1.0]),
            FREQUENCY,
            RECEIVER,
        ),
        (
            interface(340.0),
            f64::from_bits(FREQUENCY.to_bits() + 1),
            RECEIVER,
        ),
        (interface(340.0), FREQUENCY, [0.0301, 0.02, 0.13]),
    ] {
        let wrong = StressObserver::new(
            plane,
            frequency,
            receiver,
            &STRESS_POSITIONS,
            stress_limits(),
            &mut || false,
        )
        .unwrap();
        assert!(matches!(
            observe(&wrong, BIN, WORK),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    let (total, _) = observers(FullWaveCompactFieldKind::IncidentTotal, 340.0);
    let (_, returning) = observers(FullWaveCompactFieldKind::Reflected, 340.0);
    assert!(matches!(
        total.evaluate_welch_bin_with_quadrupoles(&returning, &spectrum, BIN, WORK, || false,),
        Err(AeroacError::InvalidParameter { .. })
    ));
    for bin in [0, N / 2, BIN + 1, usize::MAX] {
        assert!(matches!(
            observe(&stress, bin, WORK),
            Err(AeroacError::Welch(_))
        ));
    }
    let wrong_shape = estimate(&[[C64::ONE; 13]; 3]);
    assert!(
        compact
            .evaluate_welch_bin_with_quadrupoles(&stress, &wrong_shape, BIN, WORK, || false)
            .is_err()
    );
    for budget in [
        CoherenceBudget {
            max_modes: 2,
            ..WORK
        },
        CoherenceBudget {
            max_source_mode_pairs: 8,
            ..WORK
        },
    ] {
        assert!(matches!(
            observe(&stress, BIN, budget),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    for stop in [1, 4, 7] {
        let mut polls = 0;
        assert!(matches!(
            compact.evaluate_welch_bin_with_quadrupoles(&stress, &spectrum, BIN, WORK, || {
                polls += 1;
                polls == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let compact_record =
        estimate(&sources().map(|input| core::array::from_fn::<_, 8, _>(|i| channels(input)[i])));
    assert!(
        compact
            .evaluate_welch_bin(
                &compact_record,
                BIN,
                CoherenceBudget {
                    max_source_mode_pairs: 5,
                    ..WORK
                },
                || false
            )
            .is_err()
    );
    assert!(
        compact
            .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .is_err()
    );
    assert!(matches!(
        compact.evaluate_welch_bin(&compact_record, BIN, WORK, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        observe(&stress, BIN, WORK)
            .unwrap()
            .mean_square_pressure_pa2
            > 0.0
    );
}
