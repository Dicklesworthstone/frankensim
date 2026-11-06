//! G0/G1 format and unit checks, including two pinned evaluated-data records.

use super::*;

fn selection() -> EndfDecaySelection {
    EndfDecaySelection {
        material_number: 777,
        nuclide: NuclideId::new(27, 60, 1).unwrap(),
    }
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("fixture value does not fit ENDF");
}

fn record(fields: [String; 6], mt: i32) -> String {
    assert!(fields.iter().all(|field| field.len() == 11));
    format!("{}{:4}{:2}{mt:3}{:5}", fields.concat(), 777, 8, 1)
}

fn values(values: [f64; 6]) -> String {
    record(values.map(real), 457)
}

fn control(c1: f64, c2: f64, integers: [i32; 4]) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        457,
    )
}

fn send() -> String {
    record(std::array::from_fn(|_| format!("{:11}", 0)), 0)
}

fn fixture(energies_ev: &[(f64, f64)], modes: &[[f64; 6]]) -> Vec<String> {
    let mut lines = vec![
        control(27060.0, 59.4, [2, 1, 0, 0]),
        control(8.0, 0.25, [0, 0, (energies_ev.len() * 2) as i32, 0]),
    ];
    for chunk in energies_ev.chunks(3) {
        let mut row = [0.0; 6];
        for (index, &(energy, uncertainty)) in chunk.iter().enumerate() {
            row[2 * index] = energy;
            row[2 * index + 1] = uncertainty;
        }
        lines.push(values(row));
    }
    lines.push(control(
        0.5,
        -1.0,
        [0, 0, (modes.len() * 6) as i32, modes.len() as i32],
    ));
    for &mode in modes {
        lines.push(values(mode));
    }
    lines.push(send());
    lines
}

fn simple() -> Vec<String> {
    fixture(
        &[(50.0, 0.5), (100.0, 1.0), (0.0, 0.0)],
        &[[1.5, 2.0, 300.0, 3.0, 1.0, 0.01]],
    )
}

fn parse(lines: &[String]) -> Result<DecayEvaluation, EndfError> {
    DecayEvaluation::parse(
        "manufactured MF8".into(),
        selection(),
        &lines.join("\n"),
        EndfDecayReadLimits::default(),
        || Ok(()),
    )
}

fn field(line: &mut String, index: usize, value: &str) {
    assert!(value.len() <= 11);
    line.replace_range(index * 11..(index + 1) * 11, &format!("{value:>11}"));
}

#[test]
fn pinned_endf_b_viii_tritium_and_stable_helium_preserve_evaluated_values() {
    // Original file bytes and source attribution are retained with the fixtures.
    let h3 = DecayEvaluation::parse(
        "ENDF/B-VIII.0/H-3".into(),
        EndfDecaySelection {
            material_number: 4,
            nuclide: NuclideId::new(1, 3, 0).unwrap(),
        },
        include_str!("../../../tests/data/endf_b_viii_0/dec-001_H_003.endf"),
        EndfDecayReadLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(h3.law().half_life().unwrap().value(), 3.887_896e8);
    assert_eq!(h3.half_life_uncertainty_s(), 6.311_520e5);
    assert_eq!(h3.average_energies().len(), 3);
    assert_eq!(h3.average_energies()[0].value_j, Some(5690.0 * EV_J));
    assert_eq!(h3.average_energies()[0].uncertainty_j, 40.0 * EV_J);
    assert_eq!(h3.modes()[0].q_j, 18590.0 * EV_J);
    assert_eq!(h3.modes()[0].probability, 1.0);
    assert_eq!(h3.metadata().spectrum_count, 1);
    let beta = &h3.spectra()[0];
    assert_eq!(beta.radiation, spectra::DecayRadiation::BetaMinus);
    assert_eq!(beta.discrete_normalization.value, 1.0);
    assert_eq!(beta.discrete.len(), 1);
    assert_eq!(beta.discrete[0].energy_j.value, 18590.0 * EV_J);
    assert_eq!(beta.discrete[0].energy_j.uncertainty, 2.0 * EV_J);
    assert_eq!(beta.discrete[0].intensity.value, 1.0);
    assert_eq!(beta.discrete[0].transition_type, 1.0);
    assert_eq!(beta.average_energy_j.value_j, Some(5690.0 * EV_J));
    assert!(beta.continuum.is_none());
    assert_eq!(h3.metadata().spin, Some(0.5));
    assert_eq!(h3.metadata().parity, Some(1));
    // Neither the Q remainder nor the descriptive MF1 text becomes a guessed
    // neutrino component when the actual MF8 list contains only three values.
    assert!(h3.average_energies().get(16).is_none());

    let he3 = DecayEvaluation::parse(
        "ENDF/B-VIII.0/He-3".into(),
        EndfDecaySelection {
            material_number: 9,
            nuclide: NuclideId::new(2, 3, 0).unwrap(),
        },
        include_str!("../../../tests/data/endf_b_viii_0/dec-002_He_003.endf"),
        EndfDecayReadLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(he3.law(), DecayLawV1::stable());
    assert!(he3.modes().is_empty());
    assert_eq!(he3.branch_fraction_sum(), 0.0);
    assert_eq!(he3.metadata().spectrum_count, 0);
    assert!(he3.spectra().is_empty());
    assert!(
        he3.average_energies()
            .iter()
            .all(|energy| { energy.value_j == Some(0.0) && energy.uncertainty_j == 0.0 })
    );
}

#[test]
fn seventeen_components_unknown_values_and_mode_probabilities_remain_explicit() {
    let mut energies: Vec<_> = (0..17).map(|i| (i as f64 * 10.0, 0.25)).collect();
    energies[6] = (-1.0, 0.0);
    energies[16] = (18000.0, 20.0);
    let lines = fixture(
        &energies,
        &[
            [1.5, 2.0, 30000.0, 30.0, 0.7, 0.01],
            [6.0, 0.0, 200000.0, 0.0, 0.2, 0.02],
            [10.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    );
    let data = parse(&lines).unwrap();
    assert_eq!(data.data_id(), "manufactured MF8");
    assert_eq!(data.metadata().nuclide, selection().nuclide);
    assert_eq!(data.metadata().excited_state, 2);
    assert_eq!(data.law().half_life().unwrap().value(), 8.0);
    assert_eq!(data.half_life_uncertainty_s(), 0.25);
    assert_eq!(data.average_energies()[6].value_j, None);
    assert_eq!(data.average_energies()[16].value_j, Some(18000.0 * EV_J));
    assert_eq!(data.modes()[0].rtyp, 1.5);
    assert_eq!(data.modes()[0].daughter_isomer, 2);
    assert_eq!(data.modes()[1].rtyp, 6.0);
    assert_eq!(data.modes()[2].rtyp, 10.0);
    assert_eq!(data.modes()[2].probability, 0.0);
    assert!((data.branch_fraction_sum() - 0.9).abs() < 1.0e-15);
    assert_ne!(data.branch_fraction_sum(), 1.0);
}

#[test]
fn published_cobalt_gamma_lines_keep_multiplicity_and_conversion_coefficients() {
    let data = DecayEvaluation::parse(
        "ENDF/B-VIII.0/Co-60".into(),
        EndfDecaySelection {
            material_number: 558,
            nuclide: NuclideId::new(27, 60, 0).unwrap(),
        },
        include_str!("../../../tests/data/endf_b_viii_0/dec-027_Co_060.endf"),
        EndfDecayReadLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let gamma = &data.spectra()[0];
    assert_eq!(gamma.radiation, spectra::DecayRadiation::Gamma);
    assert_eq!(gamma.discrete.len(), 6);
    assert_eq!(gamma.discrete[2].energy_j.value, 1.173228e6 * EV_J);
    assert_eq!(gamma.discrete[2].intensity.value, 0.9985);
    assert_eq!(gamma.discrete[3].energy_j.value, 1.332492e6 * EV_J);
    assert_eq!(gamma.discrete[3].intensity.value, 0.999826);
    assert_eq!(
        gamma.discrete[2].conversion_total.unwrap().value,
        1.667988e-4
    );
    assert_eq!(gamma.discrete[2].conversion_k.unwrap().value, 1.5e-4);
    assert_eq!(gamma.discrete[2].conversion_l.unwrap().value, 1.465e-5);
    let photons: f64 = gamma.discrete.iter().map(|line| line.intensity.value).sum();
    assert!(photons > 1.99 && photons < 2.0);
    let emitted_j: f64 = gamma
        .discrete
        .iter()
        .map(|line| line.energy_j.value * line.intensity.value)
        .sum();
    // Agreement at the precision of the published rounded mean; conversion
    // coefficients and the separate electron spectrum are not extra photons.
    assert!((emitted_j / EV_J - gamma.average_energy_j.value_j.unwrap() / EV_J).abs() < 1.0);
    assert_eq!(
        data.spectra()[2].radiation,
        spectra::DecayRadiation::Electron
    );
    assert_eq!(data.spectra()[2].discrete.len(), 20);
    assert_eq!(data.spectra()[3].radiation, spectra::DecayRadiation::XRay);
}

#[test]
fn spectral_values_reject_negative_data_and_preserve_unknown_means() {
    let original = spectral_fixture();
    for (row, column, value) in [(6, 0, "-1"), (8, 2, "-1"), (13, 3, "-0.1"), (7, 4, "10")] {
        let mut broken = original.clone();
        field(&mut broken[row], column, value);
        assert!(
            parse(&broken).is_err(),
            "invalid spectral field {row}:{column}"
        );
    }
    let mut unknown = original;
    field(&mut unknown[6], 2, "-1");
    assert_eq!(
        parse(&unknown).unwrap().spectra()[0]
            .average_energy_j
            .value_j,
        None
    );
}

fn spectral_fixture() -> Vec<String> {
    let mut lines = simple();
    field(&mut lines[0], 5, "3");
    lines.pop();
    for (styp, lcon, lcov, ner) in [(0.0, 0, 0, 1), (1.0, 1, 1, 0), (2.0, 2, 3, 1)] {
        lines.push(control(0.0, styp, [lcon, lcov, 6, ner]));
        lines.push(values([1.0, 0.0, 10.0, 0.0, 1.0, 0.0]));
        if lcon != 1 {
            lines.push(control(9.0, 0.0, [0, 0, 6, 0]));
            lines.push(values([1.0, 0.0, 1.0, 0.0, 0.0, 0.0]));
        }
        if lcon != 0 {
            lines.push(control(1.0, 0.0, [0, 0, 1, 2]));
            lines.push(record([2, 2, 0, 0, 0, 0].map(|x| format!("{x:11}")), 457));
            lines.push(values([0.0, 0.0, 20.0, 0.1, 0.0, 0.0]));
            lines.push(control(0.0, 0.0, [0, 2, 4, 2]));
            lines.push(values([0.0, -0.1, 20.0, 0.0, 0.0, 0.0]));
        }
        if lcov == 3 {
            lines.push(control(0.0, 0.0, [1, 5, 3, 2]));
            lines.push(values([9.0, -1.0, 0.01, 0.0, 0.0, 0.0]));
        }
    }
    lines.push(send());
    lines
}

#[test]
fn retained_discrete_and_continuous_spectra_preserve_units_and_covariance_presence() {
    let lines = spectral_fixture();
    let data = parse(&lines).unwrap();
    let simple = parse(&simple()).unwrap();
    assert_eq!(data.metadata().spectrum_count, 3);
    assert_eq!(data.spectra().len(), 3);
    let gamma = &data.spectra()[0];
    assert_eq!(gamma.radiation, spectra::DecayRadiation::Gamma);
    assert_eq!(gamma.discrete[0].energy_j.value, 9.0 * EV_J);
    assert_eq!(gamma.discrete[0].intensity.value, 1.0);
    assert_eq!(gamma.discrete[0].decay_mode, 1.0);
    assert_eq!(gamma.covariance_flag, 0);
    assert!(gamma.continuum.is_none());
    let beta = &data.spectra()[1];
    assert_eq!(beta.covariance_flag, 1);
    let continuous = beta.continuum.as_ref().unwrap();
    assert_eq!(continuous.decay_mode, 1.0);
    assert_eq!(
        continuous.interpolation[0],
        spectra::DecayInterpolation {
            end_point: 2,
            law: 2
        }
    );
    assert_eq!(continuous.points[1].energy_j, 20.0 * EV_J);
    assert_eq!(continuous.points[1].density_per_j, 0.1 * EV_J.recip());
    // The triangular probability integral is invariant under eV -> J conversion.
    assert!(
        (continuous.points[1].energy_j * continuous.points[1].density_per_j * 0.5 - 1.0).abs()
            < 3e-16
    );
    assert_eq!(data.spectra()[2].covariance_flag, 3);
    assert_eq!(data.average_energies(), simple.average_energies());
    assert_eq!(data.modes(), simple.modes());
    assert!(parse(&lines[..lines.len() - 3]).is_err());
    let mut interrupted = lines.clone();
    interrupted[7].replace_range(70..72, " 3");
    assert!(matches!(parse(&interrupted), Err(EndfError::Format { .. })));
    let mut bad_range = lines;
    let interpolation = bad_range
        .iter()
        .position(|line| line.starts_with("          2          2"))
        .unwrap();
    field(&mut bad_range[interpolation], 0, "3");
    assert!(matches!(parse(&bad_range), Err(EndfError::Format { .. })));
}

#[test]
fn selected_identity_counts_padding_and_section_terminators_are_checked() {
    let lines = simple();
    let text = lines.join("\n");
    let duplicated = format!("{text}\n{text}");
    assert_eq!(
        DecayEvaluation::parse(
            "duplicate".into(),
            selection(),
            &duplicated,
            EndfDecayReadLimits::default(),
            || Ok(())
        ),
        Err(EndfError::AmbiguousSelection)
    );
    for id in [
        EndfDecaySelection {
            material_number: 778,
            ..selection()
        },
        EndfDecaySelection {
            nuclide: NuclideId::new(27, 60, 0).unwrap(),
            ..selection()
        },
    ] {
        assert!(
            DecayEvaluation::parse(
                "wrong selection".into(),
                id,
                &text,
                EndfDecayReadLimits::default(),
                || Ok(())
            )
            .is_err()
        );
    }
    for (row, column, value) in [(0, 4, "2"), (1, 4, "8"), (3, 4, "12"), (5, 2, "1")] {
        let mut invalid = lines.clone();
        field(&mut invalid[row], column, value);
        assert!(matches!(parse(&invalid), Err(EndfError::Format { .. })));
    }
    let mut padded = fixture(&[(0.0, 0.0); 17], &[[1.0, 0.0, 1.0, 0.0, 1.0, 0.0]]);
    field(&mut padded[7], 4, "1");
    assert!(matches!(parse(&padded), Err(EndfError::Format { .. })));
    assert!(parse(&lines[..lines.len() - 1]).is_err());
}

#[test]
fn limits_and_original_checkpoint_errors_apply_to_entire_scanned_input() {
    let text = spectral_fixture().join("\n");
    let admitted = EndfDecayReadLimits {
        max_bytes: text.len(),
        max_records: text.lines().count(),
        max_decay_modes: 1,
        max_spectra: 3,
        max_discrete_lines: 2,
        max_continuum_points: 4,
        max_interpolation_regions: 2,
    };
    assert!(
        DecayEvaluation::parse("exact limits".into(), selection(), &text, admitted, || {
            Ok(())
        })
        .is_ok()
    );
    for limited in [
        EndfDecayReadLimits {
            max_bytes: text.len() - 1,
            ..admitted
        },
        EndfDecayReadLimits {
            max_records: text.lines().count() - 1,
            ..admitted
        },
        EndfDecayReadLimits {
            max_decay_modes: 0,
            ..admitted
        },
        EndfDecayReadLimits {
            max_spectra: 2,
            ..admitted
        },
        EndfDecayReadLimits {
            max_discrete_lines: 1,
            ..admitted
        },
        EndfDecayReadLimits {
            max_continuum_points: 3,
            ..admitted
        },
        EndfDecayReadLimits {
            max_interpolation_regions: 1,
            ..admitted
        },
    ] {
        assert!(matches!(
            DecayEvaluation::parse("limited".into(), selection(), &text, limited, || Ok(())),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    for stop in [1, 12, text.lines().count()] {
        let mut calls = 0;
        assert_eq!(
            DecayEvaluation::parse("cancel".into(), selection(), &text, admitted, || {
                calls += 1;
                if calls == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(EndfError::Physics(FissionError::Cancelled))
        );
        assert_eq!(calls, stop);
    }
    let mut oversized = simple();
    field(&mut oversized[3], 4, "2147483646");
    field(&mut oversized[3], 5, "357913941");
    assert!(matches!(
        parse(&oversized),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
}

#[test]
fn unknown_sentinels_do_not_hide_invalid_or_unrepresentable_quantities() {
    let lines = simple();
    for (row, column, value) in [
        (1, 0, "0"),
        (1, 0, "1.00000-310"),
        (1, 0, "1.00000-999"),
        (1, 1, "-1"),
        (2, 0, "-2"),
        (2, 0, "1.00000-310"),
        (2, 1, "-1"),
        (4, 0, "0"),
        (4, 1, "0.5"),
        (4, 3, "-1"),
        (4, 4, "1.1"),
        (4, 5, "-1"),
    ] {
        let mut invalid = lines.clone();
        field(&mut invalid[row], column, value);
        assert!(
            parse(&invalid).is_err(),
            "accepted row {row} field {column} = {value}"
        );
    }
}

#[test]
fn stable_records_require_explicit_zero_decay_payloads() {
    let mut lines = simple();
    field(&mut lines[0], 4, "1");
    lines[1] = control(0.0, 0.0, [0, 0, 6, 0]);
    lines[2] = values([0.0; 6]);
    lines[3] = control(-77.777, 0.0, [0, 0, 6, 0]);
    lines[4] = values([0.0; 6]);
    let data = parse(&lines).unwrap();
    assert_eq!(data.law(), DecayLawV1::stable());
    assert_eq!(data.metadata().spin, None);
    assert_eq!(data.metadata().parity, None);
    assert_eq!(data.average_energies().len(), 3);
    for (row, column) in [(0, 5), (1, 0), (2, 0), (3, 5), (4, 0)] {
        let mut invalid = lines.clone();
        field(&mut invalid[row], column, "1");
        assert!(parse(&invalid).is_err());
    }
}
