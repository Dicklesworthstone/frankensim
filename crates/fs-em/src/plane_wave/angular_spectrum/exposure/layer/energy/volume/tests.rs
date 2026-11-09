//! Original sampled-volume products, physical units and stable increments.
use super::*;

fn add(
    banks: &mut Vec<EnergyBank>,
    lines: &mut Vec<EnergyLine>,
    proper: bool,
    account: usize,
    volume: f64,
    coefficient: f64,
    values: impl IntoIterator<Item = (f64, C64)>,
) {
    let start = lines.len();
    lines.extend(values.into_iter().map(|(omega, a)| EnergyLine {
        omega,
        field: [
            filtered(a, C64::ONE),
            [ScaledDd::ZERO; 2],
            [ScaledDd::ZERO; 2],
        ],
    }));
    banks.push(EnergyBank {
        start,
        end: lines.len(),
        proper,
        account,
        scale: ScaledDd::from_f64(volume)
            * ScaledDd::from_f64(coefficient)
            * ScaledDd::from_f64(0.5),
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
) -> AngularSpectrumExposureLayerVolumeEnergy {
    volume_accounts(banks, lines, time, duration, &mut || false).unwrap()
}

#[test]
fn nonuniform_volume_keeps_dark_fields_and_separate_material_accounts() {
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let mut banks = Vec::new();
        let mut lines = Vec::new();
        let mut weights = [0.0; 3];
        for (sample, volume) in [0.2, 0.03, 1.1].into_iter().enumerate() {
            for (account, coefficient) in [(0, 0.7), (0, 0.4), (1, 1.0), (2, 1.0)] {
                let amplitude = (sample + 1) as f64 * (account + 1) as f64;
                weights[account] += volume * coefficient * amplitude * amplitude;
                add(
                    &mut banks,
                    &mut lines,
                    proper,
                    account,
                    volume,
                    coefficient,
                    [1.0, -2.0, 1.0]
                        .into_iter()
                        .enumerate()
                        .map(|(i, a)| (carrier + i as f64, C64::from_re(a * amplitude))),
                );
            }
        }
        assert_eq!(observe(&banks, &lines, Dd::ZERO, None).total_j, Dd::ZERO);
        for time in [1e-6_f64, 1e-9, 1e-12] {
            let field = 4.0 * (time * 0.5).sin().powi(2);
            let product = if proper {
                0.25 * field * field
            } else {
                0.5 * (field * (2.0 * time).cos()).powi(2)
            };
            for actual in [
                observe(&banks, &lines, Dd::from_f64(time), None),
                observe(&banks, &lines, Dd::ZERO, Some(time)),
            ] {
                near(actual.electromagnetic_j, weights[0] * product, 3e-7);
                near(actual.polarization_restoring_j, weights[1] * product, 3e-7);
                near(actual.polarization_kinetic_j, weights[2] * product, 3e-7);
                near(
                    actual.total_j,
                    weights.into_iter().sum::<f64>() * product,
                    3e-7,
                );
            }
        }
    }
}

fn amplitude(sample: usize, source: usize, frequency: usize) -> C64 {
    C64::new(
        0.7 + 0.3 * sample as f64 - 0.2 * source as f64,
        0.1 + 0.2 * frequency as f64 - 0.15 * source as f64,
    )
}

#[test]
fn independent_mixed_ensemble_uses_physical_unnormalized_volume_weights() {
    let mut banks = Vec::new();
    let mut lines = Vec::new();
    let weights = [0.7, 0.03, 1.2];
    for (sample, volume) in weights.into_iter().enumerate() {
        for source in 0..3 {
            add(
                &mut banks,
                &mut lines,
                source == 2,
                0,
                volume,
                0.8,
                [1.0, 2.3]
                    .into_iter()
                    .enumerate()
                    .map(|(frequency, omega)| (omega, amplitude(sample, source, frequency))),
            );
        }
    }
    let direct = |time: f64| {
        let mut energy = 0.0;
        for (sample, volume) in weights.into_iter().enumerate() {
            for real in [-1.0, 1.0] {
                for proper in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
                    let mut field = 0.0;
                    for (source, factor) in [C64::ONE, C64::from_re(real), proper]
                        .into_iter()
                        .enumerate()
                    {
                        for (frequency, omega) in [1.0, 2.3].into_iter().enumerate() {
                            let angle = omega * time;
                            field += (amplitude(sample, source, frequency)
                                * factor
                                * C64::new(angle.cos(), angle.sin()))
                            .re;
                        }
                    }
                    energy += volume * 0.4 * field * field / 8.0;
                }
            }
        }
        energy
    };
    for (start, duration) in [(-0.3, 0.1), (0.17, 0.03), (1.2, 0.2)] {
        near(
            observe(&banks, &lines, Dd::from_f64(start), None).total_j,
            direct(start),
            3e-14,
        );
        near(
            observe(&banks, &lines, Dd::from_f64(start), Some(duration)).total_j,
            direct(start + duration) - direct(start),
            2e-12,
        );
    }
}

#[test]
fn physical_volume_precedes_products_and_unrepresentable_point_densities() {
    for (field, volume, expected) in [(1e200, 1e-200, 2.5e199), (1e-200, 1e200, 2.5e-201)] {
        let mut banks = Vec::new();
        let mut lines = Vec::new();
        add(
            &mut banks,
            &mut lines,
            true,
            0,
            volume,
            1.0,
            [(1.0, C64::from_re(field))],
        );
        near(
            observe(&banks, &lines, Dd::ZERO, None).total_j,
            expected,
            3e-14,
        );
    }
    let mut banks = Vec::new();
    let mut lines = Vec::new();
    add(
        &mut banks,
        &mut lines,
        false,
        0,
        1e100,
        1.0,
        [(1.0, C64::from_re(1e200))],
    );
    assert!(volume_accounts(&banks, &lines, Dd::ZERO, None, &mut || false).is_err());
    // Both endpoint energies overflow; the weighted analytic increment is finite.
    near(
        observe(&banks, &lines, Dd::ZERO, Some(1e-200)).total_j,
        -5e99,
        3e-14,
    );
}

#[test]
fn spatial_storage_exchange_cancels_and_precise_clock_changes_partition() {
    let mut banks = Vec::new();
    let mut lines = Vec::new();
    for value in [C64::ONE, C64::new(0.0, 1.0)] {
        add(&mut banks, &mut lines, false, 0, 0.7, 1.0, [(1.0, value)]);
    }
    let start = Dd::from_pair(2.0_f64.powi(49), 0.001);
    let first = observe(&banks[..1], &lines, start, Some(0.013)).total_j;
    let second = observe(&banks[1..], &lines, start, Some(0.013)).total_j;
    assert!(first.to_f64().abs() > 1e-5);
    assert!(second.to_f64().abs() > 1e-5);
    assert!(
        observe(&banks, &lines, start, Some(0.013))
            .total_j
            .to_f64()
            .abs()
            < 1e-30
    );
    near(observe(&banks, &lines, start, None).total_j, 0.35, 3e-14);
    add(
        &mut banks,
        &mut lines,
        true,
        1,
        0.03,
        1.0,
        [(1.0, C64::new(0.8, 0.2)), (1.2, C64::new(-0.3, 0.4))],
    );
    let part = observe(&banks, &lines, start, Some(0.02)).total_j
        + observe(&banks, &lines, start + Dd::from_f64(0.02), Some(0.02)).total_j;
    near(
        part,
        observe(&banks, &lines, start, Some(0.04)).total_j.to_f64(),
        4e-14,
    );
}

#[test]
fn cancellation_and_invalid_clocks_do_not_expose_partial_volume_accounts() {
    let mut banks = Vec::new();
    let mut lines = Vec::new();
    add(
        &mut banks,
        &mut lines,
        true,
        2,
        0.2,
        1.0,
        [(1.0, C64::ONE), (2.0, C64::ONE)],
    );
    let expected = observe(&banks, &lines, Dd::ZERO, Some(0.1));
    let mut polls = 0;
    volume_accounts(&banks, &lines, Dd::ZERO, Some(0.1), &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        volume_accounts(&banks, &lines, Dd::ZERO, Some(0.1), &mut || {
            seen += 1;
            seen == polls
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(observe(&banks, &lines, Dd::ZERO, Some(0.1)), expected);
    assert!(volume_accounts(&[], &[], Dd::from_f64(f64::NAN), None, &mut || false).is_err());
    assert_eq!(
        observe(&[], &[], Dd::ZERO, None),
        AngularSpectrumExposureLayerVolumeEnergy::default()
    );
}
