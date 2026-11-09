//! G0/G1/G3: original-field energy, reservoir work and retained source clocks.
use super::*;

fn append(
    banks: &mut Vec<EnergyBank>,
    lines: &mut Vec<EnergyLine>,
    proper: bool,
    account: usize,
    coefficient: f64,
    values: impl IntoIterator<Item = (f64, [C64; 3])>,
) {
    let start = lines.len();
    lines.extend(values.into_iter().map(|(omega, fields)| EnergyLine {
        omega,
        field: fields.map(|value| filtered(value, C64::ONE)),
    }));
    banks.push(EnergyBank {
        start,
        end: lines.len(),
        proper,
        account,
        scale: ScaledDd::from_f64(0.5) * ScaledDd::from_f64(coefficient),
    });
}

fn near(actual: Dd, expected: f64, tolerance: f64) {
    assert!(
        (actual.to_f64() - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:e}"
    );
}

fn observe(
    banks: &[EnergyBank],
    lines: &[EnergyLine],
    time: Dd,
    duration: Option<f64>,
) -> AngularSpectrumExposureLayerEnergyDensity {
    energy(banks, lines, time, duration, &mut || false).unwrap()
}

#[test]
fn dark_real_and_proper_stores_keep_the_coherent_field_before_the_square() {
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let mut banks = Vec::new();
        let mut lines = Vec::new();
        append(
            &mut banks,
            &mut lines,
            proper,
            0,
            1.0,
            [1.0, -2.0, 1.0]
                .into_iter()
                .enumerate()
                .map(|(i, a)| (carrier + i as f64, [C64::from_re(a), C64::ZERO, C64::ZERO])),
        );
        assert_eq!(observe(&banks, &lines, Dd::ZERO, None).total_j_m3, Dd::ZERO);
        for time in [1e-6_f64, 1e-9, 1e-12] {
            let amplitude = 4.0 * (0.5 * time).sin().powi(2);
            let expected = if proper {
                0.25 * amplitude * amplitude
            } else {
                0.5 * (amplitude * (2.0 * time).cos()).powi(2)
            };
            near(
                observe(&banks, &lines, Dd::from_f64(time), None).total_j_m3,
                expected,
                3e-7,
            );
            near(
                observe(&banks, &lines, Dd::ZERO, Some(time)).total_j_m3,
                expected,
                3e-7,
            );
        }
    }
}

#[test]
fn increments_retain_finite_storage_changes_below_dd_endpoint_resolution() {
    for proper in [false, true] {
        let mut banks = Vec::new();
        let mut lines = Vec::new();
        let amplitude = if proper { 1e200 } else { 1e100 };
        append(
            &mut banks,
            &mut lines,
            proper,
            0,
            1.0,
            [
                (1.0, [C64::from_re(amplitude), C64::ZERO, C64::ZERO]),
                (2.0, [C64::from_re(amplitude), C64::ZERO, C64::ZERO]),
            ],
        );
        let start = Dd::from_pair(0.3, 1e-18);
        let dt = if proper { 1e-120 } else { 1e-40 };
        assert_eq!(start, start + Dd::from_f64(dt));
        let expected = if proper {
            -0.5 * 0.3_f64.sin() * 1e280
        } else {
            -(0.3_f64.cos() + 0.6_f64.cos()) * (0.3_f64.sin() + 2.0 * 0.6_f64.sin()) * 1e160
        };
        let change = observe(&banks, &lines, start, Some(dt));
        near(change.total_j_m3, expected, 3e-14);
        assert!(change.total_j_m3.to_f64() < 0.0);
        if proper {
            // Neither bright endpoint energy fits f64/DD, but its physical
            // interval change does: query the original increments directly.
            assert!(energy(&banks, &lines, start, None, &mut || false).is_err());
        }
    }
}

#[test]
fn precise_epochs_and_partitioned_changes_keep_real_carrier_and_proper_beats() {
    for proper in [false, true] {
        let mut banks = Vec::new();
        let mut lines = Vec::new();
        append(
            &mut banks,
            &mut lines,
            proper,
            0,
            0.7,
            [
                (1.0, [C64::new(0.7, -0.2), C64::new(0.1, 0.3), C64::ZERO]),
                (2.0, [C64::new(-0.4, 0.6), C64::new(0.5, -0.1), C64::ZERO]),
                (3.0, [C64::new(0.2, 0.8), C64::new(-0.3, 0.2), C64::ZERO]),
            ],
        );
        let start = Dd::from_pair(2.0_f64.powi(49), 0.001);
        let first = observe(&banks, &lines, start, Some(0.02));
        let second = observe(&banks, &lines, start + Dd::from_f64(0.02), Some(0.02));
        let whole = observe(&banks, &lines, start, Some(0.04));
        near(
            first.total_j_m3 + second.total_j_m3,
            whole.total_j_m3.to_f64(),
            4e-14,
        );
        let initial = observe(&banks, &lines, start, None);
        let final_state = observe(&banks, &lines, start + Dd::from_f64(0.04), None);
        near(
            initial.total_j_m3 + whole.total_j_m3,
            final_state.total_j_m3.to_f64(),
            4e-14,
        );
        let rounded = observe(&banks, &lines, Dd::from_f64(start.to_f64()), None);
        assert!((rounded.total_j_m3 - initial.total_j_m3).to_f64().abs() > 1e-5);
    }
    let mut result = Vec::new();
    for carrier in [0.0, 2.0_f64.powi(49)] {
        let mut banks = Vec::new();
        let mut lines = Vec::new();
        append(
            &mut banks,
            &mut lines,
            true,
            0,
            1.0,
            [
                (carrier + 1.0, [C64::new(0.8, 0.2), C64::ZERO, C64::ZERO]),
                (carrier + 2.0, [C64::new(-0.3, 0.5), C64::ZERO, C64::ZERO]),
            ],
        );
        result.push(observe(&banks, &lines, Dd::from_f64(0.17), Some(0.013)));
    }
    assert_eq!(result[0], result[1]);
}

fn parameters() -> (f64, f64, f64, f64, f64) {
    (0.7, 0.4, 0.2, 3.0, 0.5)
}

// Independently specified original polarization phasors P and J=i omega P.
fn polarization(omega: f64, e: C64) -> [C64; 3] {
    let (delta_d, tau, delta_l, natural, gamma) = parameters();
    let pd = e.scale(delta_d) / C64::new(1.0, omega * tau);
    let pl = e.scale(delta_l * natural * natural)
        / C64::new(natural * natural - omega * omega, gamma * omega);
    [pd, pl, C64::new(-pl.im, pl.re).scale(omega)]
}

fn physical(time: f64, proper_phase: C64, real_sign: f64) -> ([f64; 3], f64, f64) {
    let (delta_d, tau, delta_l, natural, gamma) = parameters();
    let mut e = [0.0; 3];
    let mut h = [0.0; 3];
    let mut pd = [0.0; 3];
    let mut jd = [0.0; 3];
    let mut pl = [0.0; 3];
    let mut jl = [0.0; 3];
    for source in 0..3 {
        let factor = [C64::ONE, C64::from_re(real_sign), proper_phase][source];
        for (i, omega) in [2.0_f64, 2.3].into_iter().enumerate() {
            let phase = C64::new((omega * time).cos(), (omega * time).sin());
            for axis in 0..3 {
                let a = original(source, i, axis) * factor;
                let magnetic = original(source, i, axis + 3) * factor;
                let [d, l, j] = polarization(omega, a);
                e[axis] += (a * phase).re;
                h[axis] += (magnetic * phase).re;
                pd[axis] += (d * phase).re;
                jd[axis] += (C64::new(-d.im, d.re).scale(omega) * phase).re;
                pl[axis] += (l * phase).re;
                jl[axis] += (j * phase).re;
            }
        }
    }
    let norm = |v: [f64; 3]| v.into_iter().map(|x| x * x).sum::<f64>();
    let stores = [
        0.5 * (0.9 * norm(e) + 0.6 * norm(h)),
        0.5 * (norm(pd) / delta_d + norm(pl) / delta_l),
        0.5 * norm(jl) / (delta_l * natural * natural),
    ];
    let heat = tau * norm(jd) / delta_d + gamma * norm(jl) / (delta_l * natural * natural);
    let work = (0..3).map(|i| e[i] * (jd[i] + jl[i])).sum();
    (stores, heat, work)
}

fn original(source: usize, frequency: usize, axis: usize) -> C64 {
    C64::new(
        0.2 + 0.3 * axis as f64 - 0.13 * frequency as f64 - 0.17 * source as f64,
        0.12 * (1 + frequency + source) as f64 - 0.11 * axis as f64,
    )
}

fn mixed() -> (Vec<EnergyBank>, Vec<EnergyLine>) {
    let mut banks = Vec::new();
    let mut lines = Vec::new();
    let (delta_d, _, delta_l, natural, _) = parameters();
    for source in 0..3 {
        for coordinate in 0..5 {
            let (account, coefficient) = match coordinate {
                0 => (0, 0.9),
                1 => (0, 0.6),
                2 | 3 => (1, 1.0),
                _ => (2, 1.0),
            };
            append(
                &mut banks,
                &mut lines,
                source == 2,
                account,
                coefficient,
                [2.0, 2.3].into_iter().enumerate().map(|(i, omega)| {
                    let field = core::array::from_fn(|axis| match coordinate {
                        0 => original(source, i, axis),
                        1 => original(source, i, axis + 3),
                        // The retained owner stores the original normalized
                        // coordinates P/sqrt(Delta), J/(omega0 sqrt(Delta)).
                        2 => polarization(omega, original(source, i, axis))[0]
                            .scale(1.0 / delta_d.sqrt()),
                        3 => polarization(omega, original(source, i, axis))[1]
                            .scale(1.0 / delta_l.sqrt()),
                        _ => polarization(omega, original(source, i, axis))[2]
                            .scale(1.0 / (natural * delta_l.sqrt())),
                    });
                    (omega, field)
                }),
            );
        }
    }
    (banks, lines)
}

fn ensemble(time: f64) -> ([f64; 3], f64, f64) {
    let mut stores = [0.0; 3];
    let mut heat = 0.0;
    let mut work = 0.0;
    for sign in [-1.0, 1.0] {
        for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
            let (u, q, w) = physical(time, z, sign);
            for i in 0..3 {
                stores[i] += u[i] / 8.0;
            }
            heat += q / 8.0;
            work += w / 8.0;
        }
    }
    (stores, heat, work)
}

#[test]
fn independent_real_and_proper_ensemble_retains_every_physical_store() {
    let (banks, lines) = mixed();
    for time in [-0.3, 0.0, 0.17, 0.49, 1.7] {
        let actual = observe(&banks, &lines, Dd::from_f64(time), None);
        let (expected, _, _) = ensemble(time);
        near(actual.electromagnetic_j_m3, expected[0], 3e-14);
        near(actual.polarization_restoring_j_m3, expected[1], 3e-14);
        near(actual.polarization_kinetic_j_m3, expected[2], 3e-14);
        near(actual.total_j_m3, expected.into_iter().sum(), 3e-14);
    }
}

#[test]
fn reservoir_storage_change_and_irreversible_heat_close_original_polarization_work() {
    let (banks, lines) = mixed();
    for (start, duration) in [(0.17, 0.013), (0.41, 0.17), (0.7, 0.3)] {
        let delta = observe(&banks, &lines, Dd::from_f64(start), Some(duration));
        let n = 2048;
        let mut q = 0.0;
        let mut w = 0.0;
        for i in 0..=n {
            let time = start + duration * i as f64 / n as f64;
            let (_, heat, work) = ensemble(time);
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            q += weight * heat;
            w += weight * work;
        }
        q *= duration / (3 * n) as f64;
        w *= duration / (3 * n) as f64;
        near(
            delta.polarization_restoring_j_m3 + delta.polarization_kinetic_j_m3 + Dd::from_f64(q),
            w,
            2e-12,
        );
        let (before, _, _) = ensemble(start);
        let (after, _, _) = ensemble(start + duration);
        near(
            delta.total_j_m3,
            after.into_iter().sum::<f64>() - before.into_iter().sum::<f64>(),
            2e-12,
        );
    }
}

#[test]
fn filtering_and_material_scaling_precede_large_or_tiny_energy_products() {
    let mut real_banks = Vec::new();
    let mut real_lines = Vec::new();
    append(
        &mut real_banks,
        &mut real_lines,
        false,
        0,
        1.0,
        [(1.0, [C64::from_re(1e200), C64::ZERO, C64::ZERO])],
    );
    // cos(dt)-1 itself is below DD range. Its squared-field storage change
    // is finite once the original amplitude meets the retained increment.
    near(
        observe(&real_banks, &real_lines, Dd::ZERO, Some(1e-200)).total_j_m3,
        -0.5,
        3e-14,
    );
    assert!(
        energy(
            &real_banks,
            &real_lines,
            Dd::ZERO,
            Some(f64::from_bits(1)),
            &mut || false
        )
        .is_err()
    );
    for (amplitude, coefficient, expected) in [(1e200, 1e-300, 2.5e99), (1e-200, 1e300, 2.5e-101)] {
        let mut banks = Vec::new();
        let mut lines = Vec::new();
        append(
            &mut banks,
            &mut lines,
            true,
            0,
            coefficient,
            [(1.0, [C64::from_re(amplitude), C64::ZERO, C64::ZERO])],
        );
        near(
            observe(&banks, &lines, Dd::from_f64(0.4), None).total_j_m3,
            expected,
            3e-14,
        );
    }
    // The original filter product is larger than f64, yet the authored
    // physical normalization restores a representable reversible energy.
    let lines = [EnergyLine {
        omega: 1.0,
        field: [
            filtered(C64::from_re(1e200), C64::from_re(1e100)),
            [ScaledDd::ZERO; 2],
            [ScaledDd::ZERO; 2],
        ],
    }];
    let banks = [EnergyBank {
        start: 0,
        end: 1,
        proper: true,
        account: 1,
        scale: ScaledDd::from_f64(0.5) * ScaledDd::from_f64(1e-300),
    }];
    near(
        observe(&banks, &lines, Dd::ZERO, None).total_j_m3,
        2.5e299,
        3e-14,
    );
}

#[test]
fn invalid_clock_and_late_cancellation_preserve_reusable_original_fields() {
    let (banks, lines) = mixed();
    let start = Dd::from_f64(0.17);
    let mut polls = 0;
    let before = energy(&banks, &lines, start, Some(0.013), &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        energy(&banks, &lines, start, Some(0.013), &mut || {
            seen += 1;
            seen == polls
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(observe(&banks, &lines, start, Some(0.013)), before);
    for time in [
        Dd::from_f64(f64::NAN),
        Dd::from_pair(f64::MAX, f64::MAX),
        Dd {
            hi: 0.0,
            lo: f64::INFINITY,
        },
    ] {
        assert!(energy(&banks, &lines, time, None, &mut || false).is_err());
        assert!(energy(&[], &[], time, None, &mut || false).is_err());
    }
    assert_eq!(
        observe(&[], &[], Dd::ZERO, None),
        AngularSpectrumExposureLayerEnergyDensity::default()
    );
}
