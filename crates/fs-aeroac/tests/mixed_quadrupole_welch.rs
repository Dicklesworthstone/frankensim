//! Joint surface/volume observations retain recorded cross-mechanism sound.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::{
    CoherenceBudget, CompactFieldStatistics, CompactSourcePhasor, IdealAcousticPlane,
    IdealPlaneBoundary, PreparedCompactSourceObserver,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium, SpeciesRelaxation};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;
use std::f64::consts::TAU;

const N: usize = 64;
const BIN: usize = 7;
const DT: f64 = 1.0 / 8192.0;
const FREQUENCY: f64 = BIN as f64 / (N as f64 * DT);
const COMPACT_POSITIONS: [[f64; 3]; 2] = [[0.05, 0.03, 0.04], [0.01, 0.07, -0.03]];
const STRESS_POSITION: [[f64; 3]; 1] = [[0.02, -0.01, 0.03]];
const RECEIVER: [f64; 3] = [0.15, 0.12, 0.09];
const WORK: CoherenceBudget = CoherenceBudget {
    max_modes: 3,
    max_source_mode_pairs: 9,
};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.19, 344.0).unwrap()
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
    let compact = [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(2e-5, -1e-5),
            force_on_fluid_n: [
                C64::new(0.012, -0.006),
                C64::new(-0.003, 0.008),
                C64::new(0.002, 0.005),
            ],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-1e-5, 3e-5),
            force_on_fluid_n: [
                C64::new(0.004, 0.007),
                C64::new(0.009, -0.002),
                C64::new(-0.006, 0.003),
            ],
        },
    ];
    let stress = [
        C64::new(0.002, -0.001),
        C64::new(-0.001, 0.0007),
        C64::new(0.0008, 0.0003),
        C64::new(0.0015, -0.0004),
        C64::new(-0.0006, 0.0009),
        C64::new(0.0005, 0.0002),
    ];
    [
        (C64::ONE, C64::ONE),
        (C64::new(0.0, 1.0), C64::new(0.0, -1.0)),
        (C64::new(-0.4, 0.3), C64::new(0.6, -0.2)),
    ]
    .map(|(compact_phase, stress_phase)| {
        (
            [
                compact[0],
                CompactSourcePhasor {
                    volume_flow_m3_s: compact[1].volume_flow_m3_s * compact_phase,
                    force_on_fluid_n: compact[1].force_on_fluid_n.map(|v| v * compact_phase),
                },
            ],
            QuadrupolePhasor {
                stress_moment_n_m: stress.map(|q| q * stress_phase),
            },
        )
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

fn observers() -> (PreparedCompactSourceObserver, PreparedQuadrupoleObserver) {
    (
        PreparedCompactSourceObserver::new(medium(), FREQUENCY, RECEIVER, &COMPACT_POSITIONS, 2)
            .unwrap(),
        PreparedQuadrupoleObserver::new(medium(), FREQUENCY, RECEIVER, &STRESS_POSITION, 1)
            .unwrap(),
    )
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 4e-11 * expected.abs().max(1e-18),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn check_moments(actual: &CompactFieldStatistics, fields: &[[C64; 4]]) {
    let moment = |i: usize, j: usize| {
        fields
            .iter()
            .fold(C64::ZERO, |sum, field| sum + field[i] * field[j].conj())
            .scale(0.5 / fields.len() as f64)
    };
    close(C64::from_re(actual.mean_square_pressure_pa2), moment(0, 0));
    for axis in 0..3 {
        close(
            C64::new(
                actual.active_intensity_w_m2[axis],
                actual.reactive_intensity_w_m2[axis],
            ),
            moment(0, axis + 1),
        );
        for other in 0..3 {
            close(
                actual.velocity_covariance_m2_s2[axis][other],
                moment(axis + 1, other + 1),
            );
        }
    }
}

#[test]
fn g1_recorded_mixed_mechanisms_match_coherent_known_phasors_and_all_field_moments() {
    let inputs = sources();
    let spectrum = estimate(&inputs.map(channels));
    let (compact, quadrupoles) = observers();
    // Reference from independently supplied physical source amplitudes, with
    // mechanism fields added BEFORE the p/u products. No Welch factors enter.
    let fields: Vec<_> = inputs
        .iter()
        .map(|(c, q)| {
            let a = compact.evaluate(c).unwrap();
            let b = quadrupoles.evaluate(&[*q]).unwrap();
            core::array::from_fn(|i| {
                if i == 0 {
                    a.pressure_pa + b.pressure_pa
                } else {
                    a.particle_velocity_m_s[i - 1] + b.particle_velocity_m_s[i - 1]
                }
            })
        })
        .collect();
    let actual = compact
        .evaluate_welch_bin_with_quadrupoles(&quadrupoles, &spectrum, BIN, WORK, || false)
        .unwrap();
    check_moments(&actual, &fields);
    assert!(actual.scope.contains("disjoint physical mechanisms"));
    assert!(
        actual
            .scope
            .contains("no permeable-surface double counting")
    );
    let repeated = compact
        .evaluate_welch_bin_with_quadrupoles(&quadrupoles, &spectrum, BIN, WORK, || false)
        .unwrap();
    assert_eq!(
        repeated.mean_square_pressure_pa2,
        actual.mean_square_pressure_pa2
    );
    assert_eq!(
        repeated.velocity_covariance_m2_s2,
        actual.velocity_covariance_m2_s2
    );
}

#[test]
fn g1_matching_reflecting_domains_preserve_recorded_cross_mechanism_coherence() {
    let inputs = sources();
    let spectrum = estimate(&inputs.map(channels));
    let plane = IdealAcousticPlane::new(
        [0.0, 0.0, -0.1],
        [0.0, 0.0, 1.0],
        IdealPlaneBoundary::PressureRelease,
    )
    .unwrap();
    let wedge = [
        IdealAcousticPlane::new([-0.1, 0.0, 0.0], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid)
            .unwrap(),
        IdealAcousticPlane::new(
            [0.0, -0.1, 0.0],
            [0.0, 1.0, 0.0],
            IdealPlaneBoundary::PressureRelease,
        )
        .unwrap(),
    ];
    for use_wedge in [false, true] {
        let (compact, quadrupoles) = observers();
        let (compact, quadrupoles) = if use_wedge {
            (
                compact.with_wedge(wedge, 2, 8).unwrap(),
                quadrupoles.with_wedge(wedge, 2, 4).unwrap(),
            )
        } else {
            (
                compact.with_orthogonal_planes(&[plane], 4).unwrap(),
                quadrupoles.with_orthogonal_planes(&[plane], 2).unwrap(),
            )
        };
        let expected: Vec<_> = inputs
            .iter()
            .map(|(c, q)| {
                let a = compact.evaluate(c).unwrap();
                let b = quadrupoles.evaluate(&[*q]).unwrap();
                core::array::from_fn(|i| {
                    if i == 0 {
                        a.pressure_pa + b.pressure_pa
                    } else {
                        a.particle_velocity_m_s[i - 1] + b.particle_velocity_m_s[i - 1]
                    }
                })
            })
            .collect();
        let result = compact
            .evaluate_welch_bin_with_quadrupoles(&quadrupoles, &spectrum, BIN, WORK, || false)
            .unwrap();
        check_moments(&result, &expected);
        assert!(result.scope.contains("ideal"));
        let (_, free_space) = observers();
        assert!(
            compact
                .evaluate_welch_bin_with_quadrupoles(&free_space, &spectrum, BIN, WORK, || false,)
                .is_err()
        );
    }
}

#[test]
fn g3_cross_mechanism_cancellation_retains_velocity_and_changes_with_measured_correlation() {
    let acoustic = medium();
    let r = 0.21;
    let omega = TAU * FREQUENCY;
    let k = omega / acoustic.sound_speed_m_s();
    let flow = C64::from_re(2e-5);
    // On the x axis, p_Q=-i omega rho G Q and p_xx=(2/r²-k²-2ik/r) G Qxx.
    // Choose physical Qxx so the two pressures cancel at this receiver; their
    // radial gradients differ, so their velocities must not be discarded.
    let stress = C64::new(0.0, omega * acoustic.density_kg_m3()) * flow
        / C64::new(2.0 / (r * r) - k * k, -2.0 * k / r);
    let compact =
        PreparedCompactSourceObserver::new(acoustic, FREQUENCY, [r, 0.0, 0.0], &[[0.0; 3]], 1)
            .unwrap();
    let quadrupoles =
        PreparedQuadrupoleObserver::new(acoustic, FREQUENCY, [r, 0.0, 0.0], &[[0.0; 3]], 1)
            .unwrap();
    let a = [
        flow,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        stress,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
    ];
    let mut b = a;
    b[4] = -stress;
    let single = compact
        .evaluate(&[CompactSourcePhasor::monopole(flow)])
        .unwrap();
    let single_power = 0.5 * single.pressure_pa.norm_sq();
    let coherent = estimate(&[a, a]);
    let quiet = compact
        .evaluate_welch_bin_with_quadrupoles(&quadrupoles, &coherent, BIN, WORK, || false)
        .unwrap();
    assert!(quiet.mean_square_pressure_pa2 < single_power * 1e-24);
    assert!(quiet.velocity_covariance_m2_s2[0][0].re > 0.0);

    // Reversing the stress in half the measured segments cancels the cross
    // covariance and restores the sum of separate powers, 2 P. Keeping the
    // constructive phase in both segments gives 4 P, not 2 P.
    let decorrelated = estimate(&[a, b]);
    let restored = compact
        .evaluate_welch_bin_with_quadrupoles(&quadrupoles, &decorrelated, BIN, WORK, || false)
        .unwrap();
    close(
        C64::from_re(restored.mean_square_pressure_pa2),
        C64::from_re(2.0 * single_power),
    );
    let constructive = estimate(&[b, b]);
    let loud = compact
        .evaluate_welch_bin_with_quadrupoles(&quadrupoles, &constructive, BIN, WORK, || false)
        .unwrap();
    close(
        C64::from_re(loud.mean_square_pressure_pa2),
        C64::from_re(4.0 * single_power),
    );
}

#[test]
fn g4_mixed_source_medium_geometry_path_shape_work_and_cancellation_refusals_are_explicit() {
    let spectrum = estimate(&sources().map(channels));
    let (compact, quadrupoles) = observers();
    let observe =
        |c: &PreparedCompactSourceObserver, q: &PreparedQuadrupoleObserver, bin, budget| {
            c.evaluate_welch_bin_with_quadrupoles(q, &spectrum, bin, budget, || false)
        };
    for wrong in [
        PreparedQuadrupoleObserver::new(
            FluidAcousticState::try_new(1.2, 344.0).unwrap(),
            FREQUENCY,
            RECEIVER,
            &STRESS_POSITION,
            1,
        )
        .unwrap(),
        PreparedQuadrupoleObserver::new(
            medium(),
            FREQUENCY,
            [0.16, 0.12, 0.09],
            &STRESS_POSITION,
            1,
        )
        .unwrap(),
        PreparedQuadrupoleObserver::new(
            medium(),
            f64::from_bits(FREQUENCY.to_bits() + 1),
            RECEIVER,
            &STRESS_POSITION,
            1,
        )
        .unwrap(),
    ] {
        assert!(matches!(
            observe(&compact, &wrong, BIN, WORK),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    for bin in [0, N / 2, BIN + 1, usize::MAX] {
        assert!(matches!(
            observe(&compact, &quadrupoles, bin, WORK),
            Err(AeroacError::Welch(_))
        ));
    }
    let wrong_shape = estimate(&[[C64::ONE; 13]; 3]);
    assert!(
        compact
            .evaluate_welch_bin_with_quadrupoles(&quadrupoles, &wrong_shape, BIN, WORK, || false)
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
            observe(&compact, &quadrupoles, BIN, budget),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    let plane =
        IdealAcousticPlane::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], IdealPlaneBoundary::Rigid)
            .unwrap();
    let reflected = compact.clone().with_orthogonal_planes(&[plane], 4).unwrap();
    assert!(observe(&reflected, &quadrupoles, BIN, WORK).is_err());
    let walls = [
        IdealAcousticPlane::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid)
            .unwrap(),
        IdealAcousticPlane::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid)
            .unwrap(),
    ];
    let wedge = compact.clone().with_wedge(walls, 2, 8).unwrap();
    assert!(observe(&wedge, &quadrupoles, BIN, WORK).is_err());

    // Even a zero-reservoir molecular owner remains an explicitly different
    // propagation path; match its acoustic state so that this guard is tested.
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let species = SpeciesId::new("synthetic").unwrap();
    let relaxation = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: Vec::new(),
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap();
    let molecular = PreparedCompactSourceObserver::new_molecular(
        relaxation,
        FREQUENCY,
        RECEIVER,
        &COMPACT_POSITIONS,
        2,
    )
    .unwrap();
    let matched = PreparedQuadrupoleObserver::new(
        molecular.acoustic(),
        FREQUENCY,
        RECEIVER,
        &STRESS_POSITION,
        1,
    )
    .unwrap();
    assert!(observe(&molecular, &matched, BIN, WORK).is_err());

    for stop in [1, 6] {
        let mut polls = 0;
        assert!(matches!(
            compact.evaluate_welch_bin_with_quadrupoles(&quadrupoles, &spectrum, BIN, WORK, || {
                polls += 1;
                polls == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let repeated = observe(&compact, &quadrupoles, BIN, WORK).unwrap();
    assert!(repeated.mean_square_pressure_pa2 > 0.0);
}
