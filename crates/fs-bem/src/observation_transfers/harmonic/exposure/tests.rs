use super::*;
use crate::helmholtz::bc::{
    Admittance,
    incident::{self, IncidentField, Monopole},
};
use crate::helmholtz::{self, Medium};
use crate::near_field::Options;
use crate::observation_transfers::{HarmonicField, PointSourceSpec, ScatteringObservation};
use crate::panel3d::SpherePanels;

fn limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 64,
        maximum_factors: 8,
        maximum_transfer_products: 1_000_000,
        maximum_receiver_frequency_pairs: 1_000_000,
        maximum_workspace_bytes: 4_000_000,
    }
}

fn medium() -> Medium {
    Medium {
        density: 1.2,
        sound_speed: 343.0,
    }
}

fn phase(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}

fn analytic(
    omega: f64,
    fields: &[[C64; 4]],
    sources: usize,
    receivers: usize,
) -> PreparedHarmonicObservation {
    let medium = medium();
    let ids: Vec<_> = (0..sources)
        .map(|i| SourceId::StructuralMode(format!("coordinate-{i}")))
        .collect();
    let receivers: Vec<_> = (0..receivers)
        .map(|i| ReceiverSpec {
            id: format!("receiver-{i}"),
            position_m: [i as f64, 0.0, 0.0],
        })
        .collect();
    let mut samples = Vec::new();
    for field in fields {
        samples.push(field[0]);
        samples.extend(
            field[1..]
                .iter()
                .map(|&u| C64::new(0.0, omega * medium.density) * u),
        );
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        medium,
        &ids,
        &receivers,
        ContentHash([9; 32]),
        samples,
    )
    .unwrap()
}

fn pressure_only(omega: f64) -> PreparedHarmonicObservation {
    analytic(omega, &[[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]], 1, 1)
}

fn line<'a>(
    observation: &'a PreparedHarmonicObservation,
    coefficients: &'a [C64],
) -> HarmonicExposureLine<'a> {
    HarmonicExposureLine {
        observation,
        coefficients,
    }
}

fn values(point: HarmonicExposureMoments) -> [f64; 14] {
    [
        point.mean_pressure_pa,
        point.mean_particle_velocity_m_s[0],
        point.mean_particle_velocity_m_s[1],
        point.mean_particle_velocity_m_s[2],
        point.mean_square_pressure_pa2,
        point.mean_intensity_w_m2[0],
        point.mean_intensity_w_m2[1],
        point.mean_intensity_w_m2[2],
        point.velocity_second_m2_s2[0][0],
        point.velocity_second_m2_s2[0][1],
        point.velocity_second_m2_s2[0][2],
        point.velocity_second_m2_s2[1][1],
        point.velocity_second_m2_s2[1][2],
        point.velocity_second_m2_s2[2][2],
    ]
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

#[derive(Clone, Copy)]
struct Tone {
    omega: f64,
    field: [C64; 4],
}

fn tones(lines: &[HarmonicExposureLine<'_>], receiver: usize) -> Vec<Tone> {
    lines
        .iter()
        .map(|line| {
            let mut fields = vec![HarmonicField::default(); line.observation.receivers.len()];
            line.observation
                .evaluate_coherent(line.coefficients, &mut fields, || false)
                .unwrap();
            let point = fields[receiver];
            Tone {
                omega: line.observation.omega,
                field: [
                    point.pressure_pa,
                    point.particle_velocity_m_s[0],
                    point.particle_velocity_m_s[1],
                    point.particle_velocity_m_s[2],
                ],
            }
        })
        .collect()
}

// Independent real-signal Simpson integration, with an exact finite ensemble
// for the declared joint proper second moments. Each factor independently
// takes phases 1, i, -1, -i; no production temporal kernel is called here.
fn real_time_oracle(mean: &[Tone], factors: &[Vec<Tone>], window: HarmonicWindow) -> [f64; 14] {
    let n = 8192;
    let states = 4usize.pow(factors.len() as u32);
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let mut integral = [0.0; 14];
    for sample in 0..=n {
        let time = window.start_s() + window.duration_s() * sample as f64 / n as f64;
        let weight = if sample == 0 || sample == n {
            1.0
        } else if sample % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for state in 0..states {
            let mut real = [0.0; 4];
            let mut add = |lines: &[Tone], phase: C64| {
                for tone in lines {
                    let angle = tone.omega * time;
                    let (sin, cos) = angle.sin_cos();
                    for (out, amplitude) in real.iter_mut().zip(tone.field) {
                        let amplitude = amplitude * phase;
                        *out += amplitude.re * cos + amplitude.im * sin;
                    }
                }
            };
            add(mean, C64::ONE);
            let mut code = state;
            for factor in factors {
                add(factor, phases[code % 4]);
                code /= 4;
            }
            let mut raw = [0.0; 14];
            raw[..4].copy_from_slice(&real);
            let mut index = 4;
            for row in 0..4 {
                for col in row..4 {
                    raw[index] = real[row] * real[col];
                    index += 1;
                }
            }
            for (sum, value) in integral.iter_mut().zip(raw) {
                *sum += weight * value;
            }
        }
    }
    integral.map(|sum| sum / (3.0 * n as f64 * states as f64))
}

/// G1: the two-frequency observation uses the original incident PLUS
/// scattered BEM field. Independently solve simultaneous physical sources at
/// each frequency, then integrate the resulting real p/u signals in time.
#[test]
fn physical_diffraction_exposure_matches_independent_solves_and_time_integration() {
    let medium = medium();
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let sources = [
        PointSourceSpec {
            id: "left".into(),
            position_m: [0.02, 0.01, 0.09],
        },
        PointSourceSpec {
            id: "right".into(),
            position_m: [-0.04, 0.03, 0.07],
        },
    ];
    let receivers = [
        ReceiverSpec {
            id: "shadow".into(),
            position_m: [-0.03, 0.02, -0.08],
        },
        ReceiverSpec {
            id: "side".into(),
            position_m: [0.07, -0.03, 0.04],
        },
    ];
    let frequencies = [2100.0, 2600.0];
    let owner = ScatteringObservation::new(
        &surface,
        medium,
        &sources,
        &receivers,
        frequencies[1],
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap();
    let prepared = [
        owner.prepare_harmonic(frequencies[0], || false).unwrap(),
        owner.prepare_harmonic(frequencies[1], || false).unwrap(),
    ];
    let coefficients = [
        [C64::new(0.002, 0.001), C64::new(-0.003, 0.002)],
        [C64::new(-0.001, 0.002), C64::new(0.002, -0.0007)],
    ];
    let mean = [
        line(&prepared[0], &coefficients[0]),
        line(&prepared[1], &coefficients[1]),
    ];
    let exposure =
        PreparedHarmonicExposure::new(&prepared[0], "source/emission-clock", &mean, &[], limits())
            .unwrap();
    let window = HarmonicWindow::new(0.00017, 0.00083).unwrap();
    let observation = exposure.observe_window(window).unwrap();
    assert_eq!(observation.producer_hash(), owner.operator_hash());
    assert_eq!(observation.receivers(), &receivers);
    assert_eq!(observation.sources(), owner.sources());
    assert_eq!(observation.source_clock(), "source/emission-clock");
    assert!(matches!(
        observation.producer_uncertainty(),
        ProducerUncertainty::Unresolved { .. }
    ));
    let mut independent: [Vec<Tone>; 2] = core::array::from_fn(|_| Vec::new());
    for f in 0..2 {
        let monopoles: [Monopole; 2] = core::array::from_fn(|s| Monopole {
            position_m: sources[s].position_m,
            volume_velocity_m3_s: coefficients[f][s],
        });
        let solved = incident::solve_scattering(
            &surface,
            frequencies[f] / medium.sound_speed,
            medium,
            Admittance::Rigid,
            IncidentField::Monopoles(&monopoles),
            1e-3,
            helmholtz::Formulation::BurtonMiller,
        )
        .unwrap();
        let points = receivers.clone().map(|r| r.position_m);
        let field = solved
            .total_pressure_velocity_at(&points, 1e-3, Options::default())
            .unwrap();
        for r in 0..2 {
            let velocity = field.particle_velocity_m_s[r];
            independent[r].push(Tone {
                omega: frequencies[f],
                field: [
                    field.scalar.pressure[r],
                    velocity[0],
                    velocity[1],
                    velocity[2],
                ],
            });
        }
    }
    for r in 0..2 {
        let expected = real_time_oracle(&independent[r], &[], window);
        for (actual, expected) in values(observation.points()[r]).into_iter().zip(expected) {
            close(actual, expected, expected.abs().max(1e-4));
        }
    }
}

/// G0/G1/G3: proper factors preserve joint source/frequency coherence, full
/// Cartesian second moments and zero ensemble fluctuation means. Orthogonal
/// factor mixing leaves the physical covariance and finite window unchanged.
#[test]
fn mean_and_joint_proper_factors_match_finite_phase_ensemble_and_unitary_mixing() {
    let fields = [
        [
            C64::new(1.3, 0.4),
            C64::new(0.03, -0.01),
            C64::new(0.02, 0.04),
            C64::new(-0.01, 0.03),
        ],
        [
            C64::new(-0.2, 0.7),
            C64::new(0.01, 0.04),
            C64::new(-0.03, 0.02),
            C64::new(0.01, -0.02),
        ],
    ];
    let owners = [
        analytic(1200.0, &fields, 2, 1),
        analytic(1470.0, &fields, 2, 1),
    ];
    let mean_coeff = [
        [C64::new(0.3, 0.2), C64::new(-0.1, 0.4)],
        [C64::new(-0.2, 0.1), C64::new(0.2, -0.3)],
    ];
    let a = [
        [C64::new(0.7, -0.2), C64::new(0.1, 0.3)],
        [C64::new(-0.4, 0.6), C64::new(0.2, -0.1)],
    ];
    let b = [
        [C64::new(-0.2, 0.5), C64::new(0.4, -0.3)],
        [C64::new(0.1, 0.3), C64::new(-0.2, 0.7)],
    ];
    let mean = [
        line(&owners[0], &mean_coeff[0]),
        line(&owners[1], &mean_coeff[1]),
    ];
    let first = [line(&owners[0], &a[0]), line(&owners[1], &a[1])];
    let second = [line(&owners[0], &b[0]), line(&owners[1], &b[1])];
    let exposure = PreparedHarmonicExposure::new(
        &owners[0],
        "joint-clock",
        &mean,
        &[&first, &second],
        limits(),
    )
    .unwrap();
    let window = HarmonicWindow::new(0.00013, 0.0021).unwrap();
    let observed = exposure.observe_window(window).unwrap();
    let expected = real_time_oracle(
        &tones(&mean, 0),
        &[tones(&first, 0), tones(&second, 0)],
        window,
    );
    for (actual, expected) in values(observed.points()[0]).into_iter().zip(expected) {
        close(actual, expected, expected.abs().max(1e-3));
    }
    let root = core::f64::consts::FRAC_1_SQRT_2;
    let plus: [[C64; 2]; 2] =
        core::array::from_fn(|f| core::array::from_fn(|s| (a[f][s] + b[f][s]).scale(root)));
    let minus: [[C64; 2]; 2] =
        core::array::from_fn(|f| core::array::from_fn(|s| (a[f][s] - b[f][s]).scale(root)));
    let p = [line(&owners[0], &plus[0]), line(&owners[1], &plus[1])];
    let m = [line(&owners[0], &minus[0]), line(&owners[1], &minus[1])];
    let mixed =
        PreparedHarmonicExposure::new(&owners[0], "joint-clock", &mean, &[&p, &m], limits())
            .unwrap()
            .observe_window(window)
            .unwrap()
            .points()[0];
    for (actual, expected) in values(mixed).into_iter().zip(values(observed.points()[0])) {
        close(actual, expected, expected.abs().max(1e-3));
    }
    assert_eq!(
        observed.points()[0].velocity_second_m2_s2[0][1],
        observed.points()[0].velocity_second_m2_s2[1][0]
    );
    assert!(observed.points()[0].velocity_second_m2_s2[0][1].abs() > 1e-5);
}

/// G1: a reactive standing wave transports signed energy during a short
/// window even though its cycle-average active intensity is zero.
#[test]
fn subperiod_standing_wave_retains_signed_energy_flux() {
    let impedance = medium().density * medium().sound_speed;
    let root = core::f64::consts::SQRT_2;
    let owner = analytic(
        1.0,
        &[[
            C64::from_re(root),
            C64::new(0.0, root / impedance),
            C64::ZERO,
            C64::ZERO,
        ]],
        1,
        1,
    );
    let coefficients = [C64::ONE];
    let lines = [line(&owner, &coefficients)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "standing-clock", &lines, &[], limits()).unwrap();
    let window =
        HarmonicWindow::new(core::f64::consts::FRAC_PI_2, core::f64::consts::FRAC_PI_4).unwrap();
    let point = exposure.observe_window(window).unwrap().points()[0];
    close(
        point.mean_square_pressure_pa2,
        1.0 - 2.0 / core::f64::consts::PI,
        1.0,
    );
    close(
        point.velocity_second_m2_s2[0][0],
        (1.0 + 2.0 / core::f64::consts::PI) / impedance.powi(2),
        1e-5,
    );
    close(
        point.mean_intensity_w_m2[0],
        -2.0 / (core::f64::consts::PI * impedance),
        0.01,
    );
    assert!(point.mean_intensity_w_m2[0] < 0.0);
}

/// G0: beats below one f64 ulp of the individual line powers survive. Making
/// the lines independent would instead give unit pressure mean square.
#[test]
fn joint_cross_frequency_proper_factors_retain_a_sub_ulp_beat() {
    let delta = 2.0_f64.powi(-40);
    let first = pressure_only(1.0);
    let second = pressure_only(1.0 + delta);
    let positive = [C64::ONE];
    let negative = [-C64::ONE];
    let a = [line(&first, &positive)];
    let b = [line(&second, &negative)];
    let joint = [a[0], b[0]];
    let exposure =
        PreparedHarmonicExposure::new(&first, "beat-clock", &[], &[&joint], limits()).unwrap();
    let window = HarmonicWindow::new(0.0, 1.0).unwrap();
    let value = exposure.observe_window(window).unwrap().points()[0].mean_square_pressure_pa2;
    let expected = delta * delta / 6.0;
    assert!(
        value > 0.0 && (value / expected - 1.0).abs() < 2e-5,
        "tiny beat {value:e}, expected {expected:e}"
    );
    let independent =
        PreparedHarmonicExposure::new(&first, "beat-clock", &[], &[&a, &b], limits()).unwrap();
    assert_eq!(
        independent.observe_window(window).unwrap().points()[0].mean_square_pressure_pa2,
        1.0
    );
}

/// G0/G3: all original source products and repeated equal-frequency lines are
/// field-summed before squaring. Weak residuals survive strong cancellation.
#[test]
fn source_products_and_equal_frequency_groups_retain_weak_fields() {
    let field = [
        C64::ONE,
        C64::new(0.01, 0.02),
        C64::new(-0.03, 0.01),
        C64::ZERO,
    ];
    let owner = analytic(7.0, &[field; 3], 3, 1);
    let large = [C64::from_re(1e16), C64::ZERO, C64::ZERO];
    let weak = [C64::ZERO, C64::ONE, C64::ZERO];
    let opposite = [C64::from_re(-1e16), C64::ZERO, C64::ZERO];
    let combined = [C64::from_re(1e16), C64::ONE, C64::from_re(-1e16)];
    let lines = [
        line(&owner, &large),
        line(&owner, &weak),
        line(&owner, &opposite),
    ];
    let single = [line(&owner, &combined)];
    let only_weak = [line(&owner, &weak)];
    let grouped = PreparedHarmonicExposure::new(&owner, "residual", &lines, &[], limits()).unwrap();
    assert_eq!(grouped.receiver_frequency_pairs(), 1);
    let window = HarmonicWindow::new(0.17, 0.31).unwrap();
    let actual = grouped.observe_window(window).unwrap().points()[0];
    let direct = PreparedHarmonicExposure::new(&owner, "residual", &single, &[], limits())
        .unwrap()
        .observe_window(window)
        .unwrap()
        .points()[0];
    let expected = PreparedHarmonicExposure::new(&owner, "residual", &only_weak, &[], limits())
        .unwrap()
        .observe_window(window)
        .unwrap()
        .points()[0];
    for field in [actual, direct] {
        for (a, b) in values(field).into_iter().zip(values(expected)) {
            close(a, b, b.abs().max(1e-4));
        }
    }
}

/// G3: shifting the observation clock is the same physical operation as
/// applying its negative-time phase shift to every source/frequency amplitude.
#[test]
fn clock_shift_covariance_uses_every_original_frequency() {
    let owners = [pressure_only(17.0), pressure_only(23.0)];
    let coefficients = [C64::new(0.7, -0.2), C64::new(-0.4, 0.5)];
    let a = [coefficients[0]];
    let b = [coefficients[1]];
    let lines = [line(&owners[0], &a), line(&owners[1], &b)];
    let original =
        PreparedHarmonicExposure::new(&owners[0], "original", &lines, &[&lines], limits()).unwrap();
    let shift = 0.13;
    let ar = [a[0] * phase(-owners[0].omega * shift)];
    let br = [b[0] * phase(-owners[1].omega * shift)];
    let rotated = [line(&owners[0], &ar), line(&owners[1], &br)];
    let shifted =
        PreparedHarmonicExposure::new(&owners[0], "shifted", &rotated, &[&rotated], limits())
            .unwrap();
    let expected = original
        .observe_window(HarmonicWindow::new(0.07 + shift, 0.11).unwrap())
        .unwrap();
    let actual = shifted
        .observe_window(HarmonicWindow::new(0.07, 0.11).unwrap())
        .unwrap();
    for (a, b) in values(actual.points()[0])
        .into_iter()
        .zip(values(expected.points()[0]))
    {
        close(a, b, b.abs().max(1e-3));
    }
}

/// G0/G4: exact original scope matters; matching labels cannot make a different
/// geometry, source-coordinate kind, medium or producer into the same owner.
#[test]
fn owner_identity_shapes_and_resource_limits_refuse_before_publication() {
    let owner = pressure_only(1.0);
    let coefficients = [C64::ONE];
    let good = [line(&owner, &coefficients)];
    for mutation in 0..5 {
        let mut other = pressure_only(2.0);
        match mutation {
            0 => other.producer_hash = ContentHash([8; 32]),
            1 => other.receivers[0].position_m[0] += 1.0,
            2 => other.sources[0] = SourceId::PointVolumeFlow("coordinate-0".into()),
            3 => other.medium.density *= 2.0,
            _ => other.medium.sound_speed *= 2.0,
        }
        let bad = [line(&other, &coefficients)];
        assert!(PreparedHarmonicExposure::new(&owner, "clock", &bad, &[], limits()).is_err());
    }
    let empty: [C64; 0] = [];
    let nan = [C64::new(f64::NAN, 0.0)];
    for coefficients in [&empty[..], &nan[..]] {
        assert!(
            PreparedHarmonicExposure::new(
                &owner,
                "clock",
                &[line(&owner, coefficients)],
                &[],
                limits()
            )
            .is_err()
        );
    }
    assert!(PreparedHarmonicExposure::new(&owner, " ", &good, &[], limits()).is_err());
    let accepted =
        PreparedHarmonicExposure::new(&owner, "clock", &good, &[&good], limits()).unwrap();
    for budget in [
        HarmonicExposureLimits {
            maximum_lines: 1,
            ..limits()
        },
        HarmonicExposureLimits {
            maximum_factors: 0,
            ..limits()
        },
        HarmonicExposureLimits {
            maximum_transfer_products: 7,
            ..limits()
        },
        HarmonicExposureLimits {
            maximum_receiver_frequency_pairs: 1,
            ..limits()
        },
        HarmonicExposureLimits {
            maximum_workspace_bytes: accepted.workspace_bytes() - 1,
            ..limits()
        },
    ] {
        assert!(matches!(
            PreparedHarmonicExposure::new(&owner, "clock", &good, &[&good], budget),
            Err(TransferError::Cap { .. })
        ));
    }
    let too_many = vec![&good[..]; MAX_EXPOSURE_FACTORS + 1];
    let roomy = HarmonicExposureLimits {
        maximum_factors: usize::MAX,
        ..limits()
    };
    assert!(PreparedHarmonicExposure::new(&owner, "clock", &[], &too_many, roomy).is_err());
}

/// G4/G5: cancellation cannot publish a partial field and does not poison a
/// reusable owner. Empty authored banks retain their actual physical scope.
#[test]
fn cancellation_and_explicit_dark_fields_preserve_reusable_owner() {
    let owners = [pressure_only(3.0), pressure_only(5.0)];
    let a = [C64::new(0.7, 0.4)];
    let b = [C64::new(-0.2, 0.8)];
    let lines = [line(&owners[0], &a), line(&owners[1], &b)];
    let mut preparation_polls = 0;
    let prepared = PreparedHarmonicExposure::new_cancellable(
        &owners[0],
        "clock",
        &lines,
        &[&lines],
        limits(),
        || {
            preparation_polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, preparation_polls / 2, preparation_polls] {
        let mut calls = 0;
        assert!(matches!(
            PreparedHarmonicExposure::new_cancellable(
                &owners[0],
                "clock",
                &lines,
                &[&lines],
                limits(),
                || {
                    calls += 1;
                    calls == stop
                },
            ),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(calls, stop);
    }
    let window = HarmonicWindow::new(0.1, 0.2).unwrap();
    let mut query_polls = 0;
    let expected = prepared
        .observe_window_cancellable(window, || {
            query_polls += 1;
            false
        })
        .unwrap()
        .points()[0];
    for stop in [1, query_polls / 2, query_polls] {
        let mut calls = 0;
        assert!(matches!(
            prepared.observe_window_cancellable(window, || {
                calls += 1;
                calls == stop
            }),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(calls, stop);
        assert_eq!(
            prepared.observe_window(window).unwrap().points()[0],
            expected
        );
    }
    let dark = PreparedHarmonicExposure::new(&owners[0], "dark-clock", &[], &[], limits()).unwrap();
    let observation = dark.observe_window(window).unwrap();
    assert_eq!(observation.points()[0], HarmonicExposureMoments::default());
    assert_eq!(observation.producer_hash(), owners[0].producer_hash());
    assert_eq!(observation.receivers(), owners[0].receivers());
    assert_eq!(dark.receiver_frequency_pairs(), 0);
}
