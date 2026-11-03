use super::*;
use crate::helmholtz::Medium;
use crate::observation_transfers::HarmonicWorkBudget;

fn medium() -> Medium {
    Medium {
        density: 2.0,
        sound_speed: 4.0,
    }
}

fn exposure_limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 100_000,
        maximum_receiver_frequency_pairs: 100_000,
        maximum_workspace_bytes: 4_000_000,
    }
}

fn limits() -> HarmonicMomentumLimits {
    HarmonicMomentumLimits {
        maximum_quadrature_points: 32,
        maximum_quadrature_frequency_pairs: 100_000,
        maximum_generated_terms: 100_000,
        maximum_workspace_bytes: 4_000_000,
    }
}

fn observer(
    omega: f64,
    positions: &[[f64; 3]],
    columns: &[Vec<[C64; 4]>],
) -> PreparedHarmonicObservation {
    let medium = medium();
    let sources: Vec<_> = (0..columns.len())
        .map(|s| SourceId::StructuralMode(format!("source-{s}")))
        .collect();
    let receivers: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(r, &position_m)| ReceiverSpec {
            id: format!("receiver-{r}"),
            position_m,
        })
        .collect();
    let mut samples = Vec::new();
    for column in columns {
        assert_eq!(column.len(), positions.len());
        for field in column {
            samples.push(field[0]);
            samples.extend(
                field[1..]
                    .iter()
                    .map(|&u| C64::new(0.0, omega * medium.density) * u),
            );
        }
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        medium,
        &sources,
        &receivers,
        ContentHash([17; 32]),
        samples,
    )
    .unwrap()
}

fn line<'a>(
    owner: &'a PreparedHarmonicObservation,
    coefficients: &'a [C64],
) -> HarmonicExposureLine<'a> {
    HarmonicExposureLine {
        observation: owner,
        coefficients,
    }
}

fn values(v: HarmonicWaveMomentumFlux) -> [f64; 6] {
    [
        v.linear_momentum_rate_n[0],
        v.linear_momentum_rate_n[1],
        v.linear_momentum_rate_n[2],
        v.angular_momentum_rate_nm[0],
        v.angular_momentum_rate_nm[1],
        v.angular_momentum_rate_nm[2],
    ]
}

fn transfers(v: HarmonicWaveMomentumTransfer) -> [Dd; 6] {
    [
        v.linear_momentum_ns[0],
        v.linear_momentum_ns[1],
        v.linear_momentum_ns[2],
        v.angular_momentum_nms[0],
        v.angular_momentum_nms[1],
        v.angular_momentum_nms[2],
    ]
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance={tolerance:e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Independent real-time sum of the original transfer columns; no exposure
// contraction, product spectrum, temporal kernel or production stress is used.
fn real_bank(
    lines: &[HarmonicExposureLine<'_>],
    receiver: usize,
    time: f64,
    factor: C64,
) -> [f64; 4] {
    let mut field = [0.0; 4];
    for line in lines {
        let (sin, cos) = (line.observation.omega * time).sin_cos();
        for (s, coefficient) in line.coefficients.iter().enumerate() {
            let transfer =
                line.observation.transfers[s * line.observation.receivers.len() + receiver];
            for c in 0..4 {
                let amplitude = transfer[c] * *coefficient * factor;
                field[c] += amplitude.re * cos + amplitude.im * sin;
            }
        }
    }
    field
}

fn oracle(
    mean: &[HarmonicExposureLine<'_>],
    factors: &[(
        &[HarmonicExposureLine<'_>],
        HarmonicExposureFactorStatistics,
    )],
    reference: &PreparedHarmonicObservation,
    quadrature: &[HarmonicSurfaceQuadrature],
    origin: [f64; 3],
    time: f64,
) -> [f64; 6] {
    let states: usize = factors
        .iter()
        .map(|(_, law)| {
            if *law == HarmonicExposureFactorStatistics::Real {
                2
            } else {
                4
            }
        })
        .product();
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    let mut answer = [0.0; 6];
    for state in 0..states {
        for q in quadrature {
            let mut f = real_bank(mean, q.receiver, time, C64::ONE);
            let mut code = state;
            for (bank, law) in factors {
                let radix = if *law == HarmonicExposureFactorStatistics::Real {
                    2
                } else {
                    4
                };
                let add = real_bank(bank, q.receiver, time, phases[code % radix]);
                code /= radix;
                for i in 0..4 {
                    f[i] += add[i];
                }
            }
            let rho = reference.medium.density;
            let speed = reference.medium.sound_speed;
            let iso = (f[0] * f[0] / (rho * speed * speed)
                - rho * (f[1] * f[1] + f[2] * f[2] + f[3] * f[3]))
                / 2.0;
            let mut force = [0.0; 3];
            for i in 0..3 {
                for j in 0..3 {
                    force[i] += (rho * f[i + 1] * f[j + 1] + if i == j { iso } else { 0.0 })
                        * q.oriented_area_m2[j];
                }
            }
            let arm =
                core::array::from_fn(|i| reference.receivers[q.receiver].position_m[i] - origin[i]);
            let torque = cross(arm, force);
            for i in 0..3 {
                answer[i] += force[i] / states as f64;
                answer[i + 3] += torque[i] / states as f64;
            }
        }
    }
    answer
}

fn oracle_mean(
    mean: &[HarmonicExposureLine<'_>],
    factors: &[(
        &[HarmonicExposureLine<'_>],
        HarmonicExposureFactorStatistics,
    )],
    reference: &PreparedHarmonicObservation,
    quadrature: &[HarmonicSurfaceQuadrature],
    origin: [f64; 3],
    window: HarmonicWindow,
) -> [f64; 6] {
    let n = 4096;
    let mut answer = [0.0; 6];
    for i in 0..=n {
        let time = window.start_s() + window.duration_s() * i as f64 / n as f64;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, x) in answer
            .iter_mut()
            .zip(oracle(mean, factors, reference, quadrature, origin, time))
        {
            *sum += weight * x / (3.0 * n as f64);
        }
    }
    answer
}

fn fixture(omega: f64) -> PreparedHarmonicObservation {
    observer(
        omega,
        &[[1.0, 0.3, -0.2], [-0.4, 0.7, 0.6]],
        &[
            vec![
                [
                    C64::new(1.3, 0.4),
                    C64::new(0.1, -0.2),
                    C64::new(-0.3, 0.04),
                    C64::new(0.07, 0.02),
                ],
                [
                    C64::new(-0.2, 0.7),
                    C64::new(0.04, 0.2),
                    C64::new(0.1, 0.04),
                    C64::new(-0.03, 0.1),
                ],
            ],
            vec![
                [
                    C64::new(-0.5, 0.6),
                    C64::new(0.1, 0.02),
                    C64::new(0.3, -0.04),
                    C64::new(-0.1, 0.2),
                ],
                [
                    C64::new(0.3, 0.4),
                    C64::new(0.2, -0.1),
                    C64::new(-0.1, 0.2),
                    C64::new(0.3, 0.03),
                ],
            ],
        ],
    )
}

fn quadrature() -> [HarmonicSurfaceQuadrature; 2] {
    [
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [0.4, -0.2, 0.1],
        },
        HarmonicSurfaceQuadrature {
            receiver: 1,
            oriented_area_m2: [-0.3, 0.5, 0.2],
        },
    ]
}

/// G1: all six rates and finite-interval transfers agree with independent
/// Simpson time integration and a finite real/proper ensemble. Cross-frequency
/// and cross-source coherence stay inside each authored factor.
#[test]
fn mixed_fields_match_independent_real_time_stress_and_torque() {
    let owners = [fixture(7.0), fixture(11.0)];
    let a = [C64::new(0.7, -0.2), C64::new(-0.3, 0.1)];
    let b = [C64::new(-0.1, 0.4), C64::new(0.2, 0.3)];
    let c = [C64::new(0.2, 0.6), C64::new(-0.5, 0.1)];
    let mean = [line(&owners[0], &a), line(&owners[1], &b)];
    let real = [line(&owners[0], &c), line(&owners[1], &a)];
    let proper = [line(&owners[0], &b), line(&owners[1], &c)];
    let laws = [
        HarmonicExposureFactorStatistics::Real,
        HarmonicExposureFactorStatistics::ProperComplex,
    ];
    let exposure = PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "same-emission-clock",
        &mean,
        &[&real, &proper],
        &laws,
        exposure_limits(),
    )
    .unwrap();
    let q = quadrature();
    let origin = [0.2, -0.1, 0.3];
    let prepared =
        PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &q, origin, limits()).unwrap();
    assert_eq!(prepared.source_clock(), "same-emission-clock");
    assert_eq!(prepared.factor_statistics(), &laws);
    assert_eq!(prepared.producer_hash(), owners[0].producer_hash());
    let factors = [(&real[..], laws[0]), (&proper[..], laws[1])];
    for time in [0.0, 0.17, 0.83] {
        let expected = oracle(&mean, &factors, &owners[0], &q, origin, time);
        for (actual, expected) in values(prepared.evaluate(time).unwrap())
            .into_iter()
            .zip(expected)
        {
            close(actual, expected, 3e-13);
        }
    }
    let window = HarmonicWindow::new(0.17, 0.43).unwrap();
    let expected = oracle_mean(&mean, &factors, &owners[0], &q, origin, window);
    let integrated = prepared
        .transported_momentum_dd(window.start_s(), window.duration_s())
        .unwrap();
    for ((actual, integral), expected) in values(prepared.mean_flux(window).unwrap())
        .into_iter()
        .zip(transfers(integrated))
        .zip(expected)
    {
        close(actual, expected, 3e-12);
        close(integral.to_f64(), expected * window.duration_s(), 2e-12);
    }
}

/// G1/G3: a complete-period deterministic observation recovers the original
/// cycle-mean radiation-stress owner, including full moment transport.
#[test]
fn complete_cycles_match_original_momentum_and_observation_geometry() {
    let owner = fixture(7.0);
    let coefficients = [C64::new(0.7, -0.2), C64::new(-0.1, 0.4)];
    let lines = [line(&owner, &coefficients)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &lines, &[], exposure_limits()).unwrap();
    let q = quadrature();
    let origin = [0.0; 3];
    let wave = PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &q, origin, limits()).unwrap();
    let cycle = owner
        .momentum_surface(&q, origin, || false)
        .unwrap()
        .evaluate_coherent(
            &coefficients,
            HarmonicWorkBudget {
                maximum_products: 100_000,
            },
            || false,
        )
        .unwrap();
    let actual = wave
        .mean_flux(HarmonicWindow::new(0.17, core::f64::consts::TAU / 7.0).unwrap())
        .unwrap();
    for i in 0..3 {
        close(
            actual.linear_momentum_rate_n[i],
            cycle.linear_momentum_rate_n[i],
            3e-14,
        );
        close(
            actual.angular_momentum_rate_nm[i],
            cycle.angular_momentum_rate_nm[i],
            3e-14,
        );
    }
    let shift = [0.3, -0.5, 0.1];
    let shifted = PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &q, shift, limits()).unwrap();
    let f = wave.evaluate(0.19).unwrap();
    let g = shifted.evaluate(0.19).unwrap();
    let orbital = cross(shift, f.linear_momentum_rate_n);
    for i in 0..3 {
        close(
            g.linear_momentum_rate_n[i],
            f.linear_momentum_rate_n[i],
            2e-14,
        );
        close(
            g.angular_momentum_rate_nm[i],
            f.angular_momentum_rate_nm[i] - orbital[i],
            2e-14,
        );
    }
    let reversed = q.map(|q| HarmonicSurfaceQuadrature {
        receiver: q.receiver,
        oriented_area_m2: q.oriented_area_m2.map(|x| -x),
    });
    let reverse =
        PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &reversed, origin, limits()).unwrap();
    for (a, b) in values(reverse.evaluate(0.19).unwrap())
        .into_iter()
        .zip(values(f))
    {
        close(a, -b, 2e-14);
    }
}

/// G1: signed subperiod transverse wave stress survives a standing field with
/// zero cycle-mean energy flux. A travelling field has its own instantaneous
/// quadratic flux, while a proper single line is constant in expectation.
#[test]
fn travelling_and_standing_fields_keep_short_window_stress() {
    let z = medium().density * medium().sound_speed;
    let p = C64::from_re(z);
    let travelling = observer(
        1.0,
        &[[2.0, 0.0, 0.0]],
        &[vec![[p, C64::ONE, C64::ZERO, C64::ZERO]]],
    );
    let standing = observer(
        1.0,
        &[[2.0, 0.0, 0.0]],
        &[vec![[p, C64::new(0.0, 1.0), C64::ZERO, C64::ZERO]]],
    );
    let coefficient = [C64::ONE];
    let t = [line(&travelling, &coefficient)];
    let s = [line(&standing, &coefficient)];
    let et =
        PreparedHarmonicExposure::new(&travelling, "clock", &t, &[], exposure_limits()).unwrap();
    let es = PreparedHarmonicExposure::new(&standing, "clock", &s, &[], exposure_limits()).unwrap();
    let x = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [1.0, 0.0, 0.0],
    }];
    let y = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [0.0, 1.0, 0.0],
    }];
    let ts = PreparedHarmonicSurfaceWaveMomentum::new(&et, &x, [0.0; 3], limits()).unwrap();
    let ss = PreparedHarmonicSurfaceWaveMomentum::new(&es, &y, [0.0; 3], limits()).unwrap();
    for time in [0.0, 0.2, 1.4] {
        close(
            ts.evaluate(time).unwrap().linear_momentum_rate_n[0],
            2.0 * time.cos().powi(2),
            2e-14,
        );
        let observed = ss.evaluate(time).unwrap();
        close(
            observed.linear_momentum_rate_n[1],
            (2.0 * time).cos(),
            2e-14,
        );
        close(
            observed.angular_momentum_rate_nm[2],
            2.0 * (2.0 * time).cos(),
            4e-14,
        );
    }
    let w =
        HarmonicWindow::new(core::f64::consts::FRAC_PI_2, core::f64::consts::FRAC_PI_4).unwrap();
    let value = ss.mean_flux(w).unwrap();
    close(
        value.linear_momentum_rate_n[1],
        -2.0 / core::f64::consts::PI,
        2e-14,
    );
    assert!(value.linear_momentum_rate_n[1] < 0.0);
    let proper =
        PreparedHarmonicExposure::new(&travelling, "clock", &[], &[&t], exposure_limits()).unwrap();
    let ps = PreparedHarmonicSurfaceWaveMomentum::new(&proper, &x, [0.0; 3], limits()).unwrap();
    for time in [0.0, 0.17, 8.4] {
        close(
            ps.evaluate(time).unwrap().linear_momentum_rate_n[0],
            1.0,
            1e-14,
        );
    }
}

/// G3: finite transport subdivides on the original DD clock without rounding
/// away elapsed time. A coherent phase change implements a source-clock shift.
#[test]
fn finite_transport_subdivides_and_retains_clock_phase() {
    let owners = [fixture(7.0), fixture(11.0)];
    let a = [C64::new(0.7, -0.2), C64::new(0.1, 0.3)];
    let b = [C64::new(-0.3, 0.4), C64::new(0.2, -0.5)];
    let lines = [line(&owners[0], &a), line(&owners[1], &b)];
    let exposure =
        PreparedHarmonicExposure::new(&owners[0], "clock", &lines, &[&lines], exposure_limits())
            .unwrap();
    let q = quadrature();
    let surface =
        PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &q, [0.1, 0.2, 0.3], limits()).unwrap();
    for start in [
        Dd::from_f64(0.25),
        Dd::from_pair(4_503_599_627_370_496.0, 0.125),
    ] {
        let whole = transfers(
            surface
                .transported_momentum_from_clock_dd(start, 0.75)
                .unwrap(),
        );
        let first = transfers(
            surface
                .transported_momentum_from_clock_dd(start, 0.25)
                .unwrap(),
        );
        let last = transfers(
            surface
                .transported_momentum_from_clock_dd(start + Dd::from_f64(0.25), 0.5)
                .unwrap(),
        );
        for i in 0..6 {
            close((whole[i] - first[i] - last[i]).to_f64(), 0.0, 2e-14);
        }
    }
    let shift = 0.13;
    let phase = |w: f64| {
        let (s, c) = (-w * shift).sin_cos();
        C64::new(c, s)
    };
    let ar = a.map(|c| c * phase(owners[0].omega));
    let br = b.map(|c| c * phase(owners[1].omega));
    let shifted_lines = [line(&owners[0], &ar), line(&owners[1], &br)];
    let shifted = PreparedHarmonicExposure::new(
        &owners[0],
        "shifted",
        &shifted_lines,
        &[&shifted_lines],
        exposure_limits(),
    )
    .unwrap();
    let shifted_surface =
        PreparedHarmonicSurfaceWaveMomentum::new(&shifted, &q, [0.1, 0.2, 0.3], limits()).unwrap();
    for (a, b) in values(surface.evaluate(0.19 + shift).unwrap())
        .into_iter()
        .zip(values(shifted_surface.evaluate(0.19).unwrap()))
    {
        close(a, b, 3e-14);
    }
}

/// G0: field coalescing and signed-area accumulation precede tone rounding,
/// preserving a weak field/patch residual and a jointly coherent sub-ulp beat.
#[test]
fn weak_fields_signed_patches_and_slow_beats_survive() {
    let field = [C64::from_re(8.0), C64::ONE, C64::ZERO, C64::ZERO];
    let owner = observer(
        1.0,
        &[[0.0, 2.0, 0.0]],
        &[vec![field], vec![field], vec![field]],
    );
    let coefficients = [C64::from_re(1e16), C64::ONE, C64::from_re(-1e16)];
    let lines = [line(&owner, &coefficients)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &lines, &[], exposure_limits()).unwrap();
    let q = [1e16, 1.0, -1e16].map(|a| HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [a, 0.0, 0.0],
    });
    let surface =
        PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &q, [0.0; 3], limits()).unwrap();
    close(
        surface.evaluate(0.0).unwrap().linear_momentum_rate_n[0],
        2.0,
        1e-14,
    );
    close(
        surface.evaluate(0.0).unwrap().angular_momentum_rate_nm[2],
        -4.0,
        2e-14,
    );
    let delta = 2.0_f64.powi(-40);
    let owners = [
        observer(1.0, &[[0.0, 2.0, 0.0]], &[vec![field]]),
        observer(1.0 + delta, &[[0.0, 2.0, 0.0]], &[vec![field]]),
    ];
    let plus = [C64::ONE];
    let minus = [-C64::ONE];
    let joint = [line(&owners[0], &plus), line(&owners[1], &minus)];
    let beat = PreparedHarmonicExposure::new(&owners[0], "beat", &[], &[&joint], exposure_limits())
        .unwrap();
    let q = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [1.0, 0.0, 0.0],
    }];
    let wave = PreparedHarmonicSurfaceWaveMomentum::new(&beat, &q, [0.0; 3], limits()).unwrap();
    let value = wave
        .transported_momentum_dd(0.0, 1.0)
        .unwrap()
        .linear_momentum_ns[0]
        .to_f64();
    let expected = delta * delta / 3.0;
    assert!(
        value > 0.0 && (value / expected - 1.0).abs() < 3e-5,
        "{value:e}, expected {expected:e}"
    );
}

/// G0/G4/G5: admission includes every original medium, even silent lines;
/// failures and late cancellation preserve immutable preparation for replay.
#[test]
fn scope_budgets_dark_fields_and_cancellation_are_atomic() {
    let owner = fixture(7.0);
    let coefficients = [C64::ONE, C64::ZERO];
    let lines = [line(&owner, &coefficients)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &lines, &[], exposure_limits()).unwrap();
    let q = quadrature();
    let prepare = |budget, cancelled: fn() -> bool| {
        PreparedHarmonicSurfaceWaveMomentum::new_cancellable(
            &exposure, &q, [0.0; 3], budget, cancelled,
        )
    };
    let good = prepare(limits(), || false).unwrap();
    for budget in [
        HarmonicMomentumLimits {
            maximum_quadrature_points: q.len() - 1,
            ..limits()
        },
        HarmonicMomentumLimits {
            maximum_quadrature_frequency_pairs: good.quadrature_frequency_pairs() - 1,
            ..limits()
        },
        HarmonicMomentumLimits {
            maximum_generated_terms: good.generated_terms() - 1,
            ..limits()
        },
        HarmonicMomentumLimits {
            maximum_workspace_bytes: good.workspace_bytes() - 1,
            ..limits()
        },
    ] {
        assert!(matches!(
            prepare(budget, || false),
            Err(TransferError::Cap { .. })
        ));
    }
    let invalid = [HarmonicSurfaceQuadrature {
        receiver: usize::MAX,
        oriented_area_m2: [0.0; 3],
    }];
    assert!(
        PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &invalid, [0.0; 3], limits()).is_err()
    );
    assert!(
        PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &q, [f64::NAN, 0.0, 0.0], limits())
            .is_err()
    );
    let bad_area = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [f64::NAN, 0.0, 0.0],
    }];
    assert!(
        PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &bad_area, [0.0; 3], limits()).is_err()
    );
    let dark = PreparedHarmonicExposure::new(&owner, "dark", &[], &[], exposure_limits()).unwrap();
    let dark_surface =
        PreparedHarmonicSurfaceWaveMomentum::new(&dark, &q, [0.0; 3], limits()).unwrap();
    assert_eq!(
        dark_surface.evaluate(0.1).unwrap(),
        HarmonicWaveMomentumFlux::default()
    );
    assert_eq!(
        good.transported_momentum_dd(0.1, 0.0).unwrap(),
        HarmonicWaveMomentumTransfer::default()
    );
    assert!(good.evaluate(f64::NAN).is_err());
    assert!(good.transported_momentum_dd(0.0, -0.1).is_err());
    let no_stress = fixture(11.0).without_radiation_momentum();
    let silent = [C64::ZERO; 2];
    let invalid_lines = [line(&no_stress, &silent)];
    let mixed =
        PreparedHarmonicExposure::new(&owner, "clock", &invalid_lines, &[], exposure_limits())
            .unwrap();
    assert!(PreparedHarmonicSurfaceWaveMomentum::new(&mixed, &[], [0.0; 3], limits()).is_err());
    let invalid_dark =
        PreparedHarmonicExposure::new(&no_stress, "dark", &[], &[], exposure_limits()).unwrap();
    assert!(
        PreparedHarmonicSurfaceWaveMomentum::new(&invalid_dark, &[], [0.0; 3], limits()).is_err()
    );
    let distant = observer(
        1.0,
        &[[f64::MAX, 0.0, 0.0]],
        &[vec![[C64::from_re(8.0), C64::ONE, C64::ZERO, C64::ZERO]]],
    );
    let one = [C64::ONE];
    let distant_lines = [line(&distant, &one)];
    let distant_exposure =
        PreparedHarmonicExposure::new(&distant, "clock", &distant_lines, &[], exposure_limits())
            .unwrap();
    let zero = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [0.0; 3],
    }];
    let zero_surface = PreparedHarmonicSurfaceWaveMomentum::new(
        &distant_exposure,
        &zero,
        [-f64::MAX, 0.0, 0.0],
        limits(),
    )
    .unwrap();
    assert_eq!(zero_surface.quadrature_frequency_pairs(), 1);
    assert_eq!(
        zero_surface.evaluate(0.2).unwrap(),
        HarmonicWaveMomentumFlux::default()
    );
    assert_eq!(
        zero_surface.transported_momentum_dd(0.2, 0.1).unwrap(),
        HarmonicWaveMomentumTransfer::default()
    );
    let mut preparation_polls = 0;
    PreparedHarmonicSurfaceWaveMomentum::new_cancellable(&exposure, &q, [0.0; 3], limits(), || {
        preparation_polls += 1;
        false
    })
    .unwrap();
    for stop in [1, preparation_polls / 2, preparation_polls] {
        let mut calls = 0;
        assert!(matches!(
            PreparedHarmonicSurfaceWaveMomentum::new_cancellable(
                &exposure,
                &q,
                [0.0; 3],
                limits(),
                || {
                    calls += 1;
                    calls == stop
                }
            ),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(calls, stop);
    }
    let expected = good.transported_momentum_dd(0.13, 0.37).unwrap();
    let mut polls = 0;
    good.transported_momentum_dd_cancellable(0.13, 0.37, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        assert!(matches!(
            good.transported_momentum_dd_cancellable(0.13, 0.37, || {
                calls += 1;
                calls == stop
            }),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(calls, stop);
        assert_eq!(good.transported_momentum_dd(0.13, 0.37).unwrap(), expected);
    }
}
