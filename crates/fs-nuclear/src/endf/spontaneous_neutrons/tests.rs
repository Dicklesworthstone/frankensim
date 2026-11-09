//! Manufactured records exercise input semantics, not empirical isotope data.

use super::*;

fn selection() -> EndfSpontaneousNeutronSelection {
    EndfSpontaneousNeutronSelection {
        material_number: 777,
        nuclide: NuclideId::new(80, 200, 1).unwrap(),
    }
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let text = format!("{mantissa}{exponent:+}");
        if text.len() <= 11 {
            return format!("{text:>11}");
        }
    }
    panic!("manufactured number too wide");
}

fn row(
    selection: EndfSpontaneousNeutronSelection,
    fields: [String; 6],
    mf: i32,
    mt: i32,
) -> String {
    format!(
        "{}{:4}{mf:2}{mt:3}{:5}",
        fields.concat(),
        selection.material_number,
        1
    )
}

fn control(
    selection: EndfSpontaneousNeutronSelection,
    c1: f64,
    c2: f64,
    integers: [i32; 4],
    mf: i32,
    mt: i32,
) -> String {
    row(
        selection,
        [
            real(c1),
            real(c2),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        mf,
        mt,
    )
}

fn ints(selection: EndfSpontaneousNeutronSelection, values: [i32; 6], mf: i32, mt: i32) -> String {
    row(selection, values.map(|value| format!("{value:11}")), mf, mt)
}

fn append_values(
    lines: &mut Vec<String>,
    selection: EndfSpontaneousNeutronSelection,
    values: &[f64],
    mf: i32,
    mt: i32,
) {
    for chunk in values.chunks(6) {
        lines.push(row(
            selection,
            std::array::from_fn(|i| real(chunk.get(i).copied().unwrap_or(0.0))),
            mf,
            mt,
        ));
    }
}

fn za(selection: EndfSpontaneousNeutronSelection) -> f64 {
    1000.0 * f64::from(selection.nuclide.atomic_number())
        + f64::from(selection.nuclide.mass_number())
}

fn metadata(selection: EndfSpontaneousNeutronSelection, sublibrary: i32) -> Vec<String> {
    vec![
        control(
            selection,
            za(selection),
            f64::from(selection.nuclide.mass_number()),
            [-1, 0, 1, 3],
            1,
            451,
        ),
        control(
            selection,
            0.0,
            0.0,
            [
                i32::from(selection.nuclide.isomer()),
                i32::from(selection.nuclide.isomer()),
                0,
                6,
            ],
            1,
            451,
        ),
        control(selection, 0.0, 0.0, [2, 0, sublibrary, 8], 1, 451),
        control(selection, 0.0, 0.0, [0; 4], 1, 451),
        control(selection, 0.0, 0.0, [0; 4], 1, 0),
    ]
}

fn yield_section(
    selection: EndfSpontaneousNeutronSelection,
    kind: FissionNeutronKind,
    value: f64,
    families: &[f64],
) -> Vec<String> {
    let mt = kind.mt();
    let mut lines = vec![control(
        selection,
        za(selection),
        f64::from(selection.nuclide.mass_number()),
        [0, 1, 0, 0],
        1,
        mt,
    )];
    if kind == FissionNeutronKind::Delayed {
        lines.push(control(
            selection,
            0.0,
            0.0,
            [0, 0, families.len() as i32, 0],
            1,
            mt,
        ));
        append_values(&mut lines, selection, families, 1, mt);
    }
    lines.push(control(selection, 0.0, 0.0, [0, 0, 1, 0], 1, mt));
    append_values(&mut lines, selection, &[value], 1, mt);
    lines.push(control(selection, 0.0, 0.0, [0; 4], 1, 0));
    lines
}

fn spectrum_section(
    selection: EndfSpontaneousNeutronSelection,
    mt: i32,
    components: &[(f64, Vec<(f64, f64)>)],
    incident_rows: usize,
    outgoing_law: i32,
) -> Vec<String> {
    let mut lines = vec![control(
        selection,
        za(selection),
        f64::from(selection.nuclide.mass_number()),
        [0, 0, components.len() as i32, 0],
        5,
        mt,
    )];
    for (probability, points) in components {
        lines.push(control(selection, 0.0, 0.0, [0, 1, 1, 2], 5, mt));
        lines.push(ints(selection, [2, 2, 0, 0, 0, 0], 5, mt));
        append_values(
            &mut lines,
            selection,
            &[0.0, *probability, 1.0, *probability],
            5,
            mt,
        );
        lines.push(control(
            selection,
            0.0,
            0.0,
            [0, 0, 1, incident_rows as i32],
            5,
            mt,
        ));
        lines.push(ints(
            selection,
            [incident_rows as i32, 2, 0, 0, 0, 0],
            5,
            mt,
        ));
        for incident in 0..incident_rows {
            lines.push(control(
                selection,
                0.0,
                incident as f64,
                [0, 0, 1, points.len() as i32],
                5,
                mt,
            ));
            lines.push(ints(
                selection,
                [points.len() as i32, outgoing_law, 0, 0, 0, 0],
                5,
                mt,
            ));
            let values: Vec<_> = points.iter().flat_map(|&(e, p)| [e, p]).collect();
            append_values(&mut lines, selection, &values, 5, mt);
        }
    }
    lines.push(control(selection, 0.0, 0.0, [0; 4], 5, 0));
    lines
}

/// Shared manufactured prompt source for native binding regressions. Energies
/// and densities are eV and eV^-1 in the source; no production constructor is
/// exposed. Each outgoing component is supplied in its original point order.
pub(in crate::endf) fn prompt_fixture(
    selection: EndfSpontaneousNeutronSelection,
    yield_per_fission: f64,
    components: &[(f64, Vec<(f64, f64)>)],
) -> String {
    let mut lines = metadata(selection, 4);
    lines.extend(yield_section(
        selection,
        FissionNeutronKind::Prompt,
        yield_per_fission,
        &[],
    ));
    lines.extend(spectrum_section(selection, 18, components, 1, 2));
    lines.join("\n")
}

fn triangle() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (2.0, 0.5), (4.0, 0.0)]
}

fn parse_yield(
    text: &str,
    kind: FissionNeutronKind,
) -> Result<SpontaneousFissionNeutronYield, EndfError> {
    SpontaneousFissionNeutronYield::parse(
        "manufactured-SF-yield".into(),
        selection(),
        kind,
        text,
        SpontaneousNeutronReadLimits::default(),
        || Ok(()),
    )
}

fn parse_spectrum(
    text: &str,
    kind: FissionNeutronKind,
) -> Result<SpontaneousFissionNeutronSpectrum, EndfError> {
    SpontaneousFissionNeutronSpectrum::parse(
        "manufactured-SF-spectrum".into(),
        selection(),
        kind,
        text,
        SpontaneousNeutronReadLimits::default(),
        || Ok(()),
    )
}

fn changed(text: &str, line: usize, field: usize, value: &str) -> String {
    let mut lines: Vec<_> = text.lines().map(str::to_owned).collect();
    lines[line].replace_range(field * 11..(field + 1) * 11, value);
    lines.join("\n")
}

#[test]
fn g0_constant_yield_kinds_lifetimes_and_source_identity_are_distinct() {
    let mut lines = metadata(selection(), 4);
    for (kind, value, families) in [
        (FissionNeutronKind::Prompt, 2.25, vec![]),
        (FissionNeutronKind::Total, 3.0, vec![]),
        (FissionNeutronKind::Delayed, 0.75, vec![0.7, 0.1]),
    ] {
        lines.extend(yield_section(selection(), kind, value, &families));
    }
    let text = lines.join("\n");
    let prompt = parse_yield(&text, FissionNeutronKind::Prompt).unwrap();
    assert_eq!(prompt.yield_per_fission(), 2.25);
    assert_eq!(prompt.kind(), FissionNeutronKind::Prompt);
    assert!(prompt.decay_constants_s_inv().is_empty());
    assert_eq!(prompt.metadata().nuclide, selection().nuclide);
    assert_eq!(prompt.metadata().library_version, 8);
    assert_eq!(prompt.metadata().library_release, 2);
    assert_eq!(prompt.metadata().material_revision, 3);
    assert_eq!(prompt.metadata().resonance_flag, -1);
    assert_eq!(
        parse_yield(&text, FissionNeutronKind::Total)
            .unwrap()
            .yield_per_fission(),
        3.0
    );
    let delayed = parse_yield(&text, FissionNeutronKind::Delayed).unwrap();
    assert_eq!(delayed.yield_per_fission(), 0.75);
    assert_eq!(delayed.decay_constants_s_inv(), &[0.7, 0.1]);
    let zero = prompt_fixture(selection(), 0.0, &[(1.0, triangle())]);
    assert_eq!(
        parse_yield(&zero, FissionNeutronKind::Prompt)
            .unwrap()
            .yield_per_fission(),
        0.0
    );
    assert_eq!(
        parse_yield(&zero, FissionNeutronKind::Total),
        Err(EndfError::MissingSelection)
    );
}

#[test]
fn g0_lf1_keeps_component_order_source_grids_and_unrepaired_values() {
    let text = prompt_fixture(
        selection(),
        2.25,
        &[
            (0.25, triangle()),
            (0.5, vec![(1.0, 0.0), (2.0, 1.0), (3.0, 0.0)]),
        ],
    );
    let source = parse_spectrum(&text, FissionNeutronKind::Prompt).unwrap();
    assert_eq!(
        source.metadata(),
        parse_yield(&text, FissionNeutronKind::Prompt)
            .unwrap()
            .metadata()
    );
    assert_eq!(source.kind(), FissionNeutronKind::Prompt);
    assert_eq!(source.data_id(), "manufactured-SF-spectrum");
    assert_eq!(source.components().len(), 2);
    assert_eq!(
        source
            .components()
            .iter()
            .map(|c| c.probability)
            .sum::<f64>(),
        0.75
    );
    let first = &source.components()[0];
    assert_eq!(
        first.probability_energies_j.as_ref(),
        &[0.0, super::super::EV_J]
    );
    assert_eq!(first.incident_energies_j.as_ref(), &[0.0]);
    assert_eq!(first.points[1].energy_j, 2.0 * super::super::EV_J);
    assert_eq!(
        first.points[1].density_per_j,
        0.5 * super::super::EV_J.recip()
    );
    assert_eq!(first.interpolation[0].law, 2);
    assert_eq!(source, source.clone());

    let mut delayed = metadata(selection(), 4);
    delayed.extend(spectrum_section(
        selection(),
        455,
        &[(0.25, triangle()), (0.75, triangle())],
        2,
        1,
    ));
    let delayed = parse_spectrum(&delayed.join("\n"), FissionNeutronKind::Delayed).unwrap();
    assert_eq!(delayed.components()[0].probability, 0.25);
    assert_eq!(delayed.components()[1].probability, 0.75);
    assert_eq!(delayed.components()[0].incident_energies_j.len(), 2);
    assert_eq!(delayed.components()[0].interpolation[0].law, 1);
}

#[test]
fn g0_identity_source_kind_and_energy_dependence_cannot_be_silently_reinterpreted() {
    let text = prompt_fixture(selection(), 2.25, &[(1.0, triangle())]);
    for nsub in [5, 10, 11] {
        let other = changed(&text, 2, 4, &format!("{nsub:11}"));
        assert!(matches!(
            parse_yield(&other, FissionNeutronKind::Prompt),
            Err(EndfError::Unsupported(_))
        ));
        assert!(matches!(
            parse_spectrum(&other, FissionNeutronKind::Prompt),
            Err(EndfError::Unsupported(_))
        ));
    }
    assert!(
        SpontaneousFissionNeutronYield::parse(
            "wrong-isomer".into(),
            EndfSpontaneousNeutronSelection {
                nuclide: NuclideId::new(80, 200, 0).unwrap(),
                ..selection()
            },
            FissionNeutronKind::Prompt,
            &text,
            SpontaneousNeutronReadLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    let duplicate = format!("{text}\n{text}");
    assert_eq!(
        parse_yield(&duplicate, FissionNeutronKind::Prompt),
        Err(EndfError::AmbiguousSelection)
    );
    assert_eq!(
        parse_spectrum(&duplicate, FissionNeutronKind::Prompt),
        Err(EndfError::AmbiguousSelection)
    );
    assert!(parse_spectrum(&text, FissionNeutronKind::Total).is_err());
    assert_eq!(
        parse_spectrum(&text, FissionNeutronKind::Delayed),
        Err(EndfError::MissingSelection)
    );
    // MF1 selected HEAD, polynomial LIST and MF5 LF declaration respectively.
    assert!(
        parse_yield(
            &changed(&text, 5, 3, &format!("{:11}", 2)),
            FissionNeutronKind::Prompt
        )
        .is_err()
    );
    assert!(
        parse_yield(
            &changed(&text, 6, 4, &format!("{:11}", 2)),
            FissionNeutronKind::Prompt
        )
        .is_err()
    );
    assert!(
        parse_spectrum(
            &changed(&text, 10, 3, &format!("{:11}", 9)),
            FissionNeutronKind::Prompt
        )
        .is_err()
    );
    assert!(
        parse_spectrum(
            &changed(&text, 12, 3, &real(0.5)),
            FissionNeutronKind::Prompt
        )
        .is_err()
    );
    let mut repeated = metadata(selection(), 4);
    repeated.extend(spectrum_section(
        selection(),
        18,
        &[(1.0, triangle())],
        2,
        2,
    ));
    let repeated = repeated.join("\n");
    assert!(parse_spectrum(&repeated, FissionNeutronKind::Prompt).is_ok());
    // Second outgoing row's center density differs; no incident-dependent fit.
    let changed_row = changed(&repeated, 16, 3, &real(0.75));
    assert!(matches!(
        parse_spectrum(&changed_row, FissionNeutronKind::Prompt),
        Err(EndfError::Unsupported(_))
    ));
}

#[test]
fn g0_lf1_domains_and_declared_source_limits_are_enforced() {
    let text = prompt_fixture(selection(), 2.25, &[(1.0, triangle())]);
    let limits = SpontaneousNeutronReadLimits::default();
    for bounded in [
        SpontaneousNeutronReadLimits {
            records: EndfReadLimits {
                max_bytes: text.len() - 1,
                ..limits.records
            },
            ..limits
        },
        SpontaneousNeutronReadLimits {
            records: EndfReadLimits {
                max_records: 12,
                ..limits.records
            },
            ..limits
        },
        SpontaneousNeutronReadLimits {
            max_components: 0,
            ..limits
        },
        SpontaneousNeutronReadLimits {
            max_incident_points: 2,
            ..limits
        },
        SpontaneousNeutronReadLimits {
            max_spectrum_points: 2,
            ..limits
        },
        SpontaneousNeutronReadLimits {
            max_interpolation_regions: 2,
            ..limits
        },
    ] {
        assert!(matches!(
            SpontaneousFissionNeutronSpectrum::parse(
                "bounded".into(),
                selection(),
                FissionNeutronKind::Prompt,
                &text,
                bounded,
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    assert!(
        parse_spectrum(
            &prompt_fixture(selection(), 1.0, &[(1.0, vec![(0.0, 0.5), (4.0, 0.0)])]),
            FissionNeutronKind::Prompt
        )
        .is_err()
    );
    let mut log = metadata(selection(), 4);
    log.extend(spectrum_section(
        selection(),
        18,
        &[(1.0, vec![(1.0, 0.0), (2.0, 0.5), (4.0, 0.0)])],
        1,
        3,
    ));
    assert!(parse_spectrum(&log.join("\n"), FissionNeutronKind::Prompt).is_ok());
    assert!(
        parse_spectrum(
            &changed(&text, 16, 1, &format!("{:11}", 4)),
            FissionNeutronKind::Prompt
        )
        .is_err()
    );
    let truncated = text.lines().take(17).collect::<Vec<_>>().join("\n");
    assert!(parse_spectrum(&truncated, FissionNeutronKind::Prompt).is_err());
    let mut delayed = metadata(selection(), 4);
    delayed.extend(yield_section(
        selection(),
        FissionNeutronKind::Delayed,
        0.5,
        &[0.7, 0.1],
    ));
    assert!(matches!(
        SpontaneousFissionNeutronYield::parse(
            "bounded-delayed".into(),
            selection(),
            FissionNeutronKind::Delayed,
            &delayed.join("\n"),
            SpontaneousNeutronReadLimits {
                max_delayed_families: 1,
                ..limits
            },
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
}

#[test]
fn g4_final_cancellation_does_not_publish_a_partial_source() {
    let text = prompt_fixture(selection(), 2.25, &[(1.0, triangle())]);
    let mut polls = 0;
    SpontaneousFissionNeutronSpectrum::parse(
        "poll-count".into(),
        selection(),
        FissionNeutronKind::Prompt,
        &text,
        SpontaneousNeutronReadLimits::default(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        SpontaneousFissionNeutronSpectrum::parse(
            "cancelled".into(),
            selection(),
            FissionNeutronKind::Prompt,
            &text,
            SpontaneousNeutronReadLimits::default(),
            || {
                seen += 1;
                if seen == polls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(EndfError::Physics(FissionError::Cancelled))
    );
    assert!(parse_spectrum(&text, FissionNeutronKind::Prompt).is_ok());
}
