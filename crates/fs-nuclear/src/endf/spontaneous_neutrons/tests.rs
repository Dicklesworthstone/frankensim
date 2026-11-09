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

/// LF5 input with one explicit zero-energy incident row. Shape coordinates and
/// densities are dimensionless; theta and U are source eV values. Singleton
/// INT5 metadata deliberately exercises the absence of an interpolation interval.
pub(in crate::endf) fn prompt_general_evaporation_fixture(
    selection: EndfSpontaneousNeutronSelection,
    yield_per_fission: f64,
    probability: f64,
    temperature_ev: f64,
    restriction_ev: f64,
    shape: &[(f64, f64)],
    interpolation_law: i32,
) -> String {
    let mut lines = metadata(selection, 4);
    lines.extend(yield_section(
        selection,
        FissionNeutronKind::Prompt,
        yield_per_fission,
        &[],
    ));
    lines.extend(evaporation_prefix(
        selection,
        probability,
        0.0,
        temperature_ev,
        restriction_ev,
        5,
    ));
    lines.push(control(
        selection,
        0.0,
        0.0,
        [0, 0, 1, shape.len() as i32],
        5,
        18,
    ));
    lines.push(ints(
        selection,
        [shape.len() as i32, interpolation_law, 0, 0, 0, 0],
        5,
        18,
    ));
    let values: Vec<_> = shape.iter().flat_map(|&(x, g)| [x, g]).collect();
    append_values(&mut lines, selection, &values, 5, 18);
    lines.push(control(selection, 0.0, 0.0, [0; 4], 5, 0));
    lines.join("\n")
}

/// LF7/9 input with one explicitly supplied reference E and finite E-U cutoff.
pub(in crate::endf) fn prompt_analytic_fixture(
    selection: EndfSpontaneousNeutronSelection,
    yield_per_fission: f64,
    probability: f64,
    incident_ev: f64,
    temperature_ev: f64,
    restriction_ev: f64,
    law: i32,
) -> String {
    let mut lines = metadata(selection, 4);
    lines.extend(yield_section(
        selection,
        FissionNeutronKind::Prompt,
        yield_per_fission,
        &[],
    ));
    lines.extend(evaporation_prefix(
        selection,
        probability,
        incident_ev,
        temperature_ev,
        restriction_ev,
        law,
    ));
    lines.push(control(selection, 0.0, 0.0, [0; 4], 5, 0));
    lines.join("\n")
}

fn evaporation_prefix(
    selection: EndfSpontaneousNeutronSelection,
    probability: f64,
    incident_ev: f64,
    temperature_ev: f64,
    restriction_ev: f64,
    law: i32,
) -> Vec<String> {
    let mut lines = Vec::new();
    lines.push(control(
        selection,
        za(selection),
        f64::from(selection.nuclide.mass_number()),
        [0, 0, 1, 0],
        5,
        18,
    ));
    lines.push(control(
        selection,
        restriction_ev,
        0.0,
        [0, law, 1, 1],
        5,
        18,
    ));
    lines.push(ints(selection, [1, 5, 0, 0, 0, 0], 5, 18));
    append_values(&mut lines, selection, &[incident_ev, probability], 5, 18);
    lines.push(control(selection, 0.0, 0.0, [0, 0, 1, 1], 5, 18));
    lines.push(ints(selection, [1, 5, 0, 0, 0, 0], 5, 18));
    append_values(&mut lines, selection, &[incident_ev, temperature_ev], 5, 18);
    lines
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
    assert_eq!(first.restriction_energy_j, None);
    let SpontaneousFissionSpectrumLaw::Tabulated {
        points,
        interpolation,
    } = &first.law
    else {
        panic!("LF1 must retain its tabulated law");
    };
    assert_eq!(points[1].energy_j, 2.0 * super::super::EV_J);
    assert_eq!(points[1].density_per_j, 0.5 * super::super::EV_J.recip());
    assert_eq!(interpolation[0].law, 2);
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
    assert!(matches!(
        &delayed.components()[0].law,
        SpontaneousFissionSpectrumLaw::Tabulated { interpolation, .. }
            if interpolation[0].law == 1
    ));
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
fn g0_lf5_retains_dimensionless_shape_temperature_restriction_and_singletons() {
    let shape = [(0.5, 0.5), (1.0, 3.0), (1.5, 0.25)];
    let text = prompt_general_evaporation_fixture(selection(), 2.25, 0.75, 2.0, -3.0, &shape, 2);
    let source = parse_spectrum(&text, FissionNeutronKind::Prompt).unwrap();
    let component = &source.components()[0];
    assert_eq!(component.probability, 0.75);
    assert_eq!(component.probability_energies_j.as_ref(), &[0.0]);
    assert_eq!(component.incident_energies_j.as_ref(), &[0.0]);
    assert_eq!(
        component.restriction_energy_j,
        Some(-3.0 * super::super::EV_J)
    );
    let SpontaneousFissionSpectrumLaw::GeneralEvaporation {
        temperature_j,
        shape: original_shape,
        interpolation,
    } = &component.law
    else {
        panic!("LF5 must retain its general evaporation law");
    };
    assert_eq!(*temperature_j, 2.0 * super::super::EV_J);
    assert_eq!(interpolation[0].law, 2);
    assert_eq!(original_shape.len(), shape.len());
    for (point, &(x, g)) in original_shape.iter().zip(&shape) {
        assert_eq!(point.coordinate, x);
        assert_eq!(point.density, g);
    }
    assert_eq!(source, source.clone());
    // Every ENDF interpolation law is valid metadata for a singleton: there
    // is no log(0) interval to evaluate at its supplied zero incident energy.
    for law in 1..=5 {
        let changed_law = changed(&text, 11, 1, &format!("{law:11}"));
        let changed_law = changed(&changed_law, 14, 1, &format!("{law:11}"));
        let changed_source = parse_spectrum(&changed_law, FissionNeutronKind::Prompt).unwrap();
        assert_eq!(
            changed_source.components()[0].probability_interpolation[0].law,
            law as u8
        );
        assert_eq!(
            changed_source.components()[0].incident_interpolation[0].law,
            law as u8
        );
    }
    // The incident value is retained; the parser does not impose E=0 as a
    // spontaneous-fission convention or clip the original outgoing curve.
    let shifted = changed(&text, 15, 0, &real(2.0));
    let shifted = changed(&shifted, 10, 0, &real(0.5));
    let shifted = parse_spectrum(&shifted, FissionNeutronKind::Prompt).unwrap();
    assert_eq!(
        shifted.components()[0].incident_energies_j.as_ref(),
        &[2.0 * super::super::EV_J]
    );
    assert_eq!(
        shifted.components()[0].restriction_energy_j,
        Some(0.5 * super::super::EV_J)
    );
    assert_eq!(shifted.components()[0].law, component.law);
}

#[test]
fn g0_lf5_interpolation_domains_apply_to_actual_shape_intervals() {
    let fixture = |shape: &[(f64, f64)], law| {
        prompt_general_evaporation_fixture(selection(), 2.25, 1.0, 2.0, -10.0, shape, law)
    };
    for law in 1..=5 {
        let text = fixture(&[(0.5, 0.25), (1.0, 0.5), (2.0, 0.25)], law);
        assert!(parse_spectrum(&text, FissionNeutronKind::Prompt).is_ok());
    }
    // Positive endpoint densities are valid for LF5. A zero-width jump does
    // not apply log-density interpolation to either side of the jump.
    for shape in [
        vec![(0.0, 0.5), (1.0, 0.5)],
        vec![(0.0, 0.0), (0.0, 0.5), (1.0, 0.5)],
    ] {
        assert!(parse_spectrum(&fixture(&shape, 4), FissionNeutronKind::Prompt).is_ok());
    }
    for (shape, law) in [
        (vec![(0.0, 0.5), (1.0, 0.5)], 3),
        (vec![(0.0, 0.5), (1.0, 0.5)], 5),
        (vec![(1.0, 0.0), (2.0, 0.5)], 4),
        (vec![(1.0, 0.5), (2.0, 0.0)], 5),
        (vec![(1.0, 0.5), (1.0, 0.0)], 1),
        (vec![(0.0, 0.0), (1.0, 0.5), (1.0, 0.4), (1.0, 0.0)], 2),
        (vec![(0.0, 0.5), (1.0, 0.5)], 6),
    ] {
        assert!(parse_spectrum(&fixture(&shape, law), FissionNeutronKind::Prompt).is_err());
    }
}

#[test]
fn g0_lf5_refuses_ambiguous_temperature_rows_invalid_cutoffs_and_excess_work() {
    let text = prompt_general_evaporation_fixture(
        selection(),
        2.25,
        1.0,
        2.0,
        -3.0,
        &[(0.5, 0.0), (1.0, 2.0), (1.5, 0.0)],
        2,
    );
    for (line, field, value) in [
        (15, 1, real(0.0)),
        (15, 1, real(-2.0)),
        (15, 0, real(-1.0)),
        (10, 0, real(0.0)),
        (10, 0, real(1.0)),
        (13, 5, format!("{:11}", 0)),
        (13, 5, format!("{:11}", 2)),
        (15, 2, real(1.0)),
    ] {
        assert!(
            parse_spectrum(
                &changed(&text, line, field, &value),
                FissionNeutronKind::Prompt
            )
            .is_err()
        );
    }
    let limits = SpontaneousNeutronReadLimits::default();
    for bounded in [
        SpontaneousNeutronReadLimits {
            max_incident_points: 1,
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
                "bounded-LF5".into(),
                selection(),
                FissionNeutronKind::Prompt,
                &text,
                bounded,
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    let truncated = text.lines().take(18).collect::<Vec<_>>().join("\n");
    assert!(parse_spectrum(&truncated, FissionNeutronKind::Prompt).is_err());
}

#[test]
fn g0_lf7_lf9_retain_explicit_reference_temperature_and_finite_cutoff() {
    for law in [7, 9] {
        let text = prompt_analytic_fixture(selection(), 2.25, 0.75, 4.0, 2.0, 1.0, law);
        let source = parse_spectrum(&text, FissionNeutronKind::Prompt).unwrap();
        let component = &source.components()[0];
        assert_eq!(component.probability, 0.75);
        assert_eq!(
            component.probability_energies_j.as_ref(),
            &[4.0 * super::super::EV_J]
        );
        assert_eq!(
            component.incident_energies_j.as_ref(),
            &[4.0 * super::super::EV_J]
        );
        assert_eq!(component.restriction_energy_j, Some(super::super::EV_J));
        assert_eq!(component.probability_interpolation[0].law, 5);
        assert_eq!(component.incident_interpolation[0].law, 5);
        let expected = if law == 7 {
            SpontaneousFissionSpectrumLaw::Maxwellian {
                temperature_j: 2.0 * super::super::EV_J,
            }
        } else {
            SpontaneousFissionSpectrumLaw::Evaporation {
                temperature_j: 2.0 * super::super::EV_J,
            }
        };
        assert_eq!(component.law, expected);
        // No sampled outgoing points or infinite-support substitution are
        // introduced. The zero-energy reference remains a valid singleton.
        let zero = prompt_analytic_fixture(selection(), 2.25, 1.0, 0.0, 2.0, -3.0, law);
        let zero = SpontaneousFissionNeutronSpectrum::parse(
            "analytic-without-tabulation".into(),
            selection(),
            FissionNeutronKind::Prompt,
            &zero,
            SpontaneousNeutronReadLimits {
                max_spectrum_points: 0,
                ..SpontaneousNeutronReadLimits::default()
            },
            || Ok(()),
        )
        .unwrap();
        assert_eq!(zero.components()[0].incident_energies_j.as_ref(), &[0.0]);
        assert_eq!(
            zero.components()[0].restriction_energy_j,
            Some(-3.0 * super::super::EV_J)
        );
    }
}

#[test]
fn g0_lf7_lf9_refuse_undefined_parameters_incident_dependence_and_excess_work() {
    for law in [7, 9] {
        let text = prompt_analytic_fixture(selection(), 2.25, 1.0, 0.0, 2.0, -3.0, law);
        for (line, field, value) in [
            (15, 1, real(0.0)),
            (15, 1, real(-2.0)),
            (15, 0, real(-1.0)),
            (10, 0, real(0.0)),
            (10, 0, real(1.0)),
            (13, 5, format!("{:11}", 0)),
            (15, 2, real(1.0)),
        ] {
            assert!(
                parse_spectrum(
                    &changed(&text, line, field, &value),
                    FissionNeutronKind::Prompt
                )
                .is_err()
            );
        }
        // A second incident energy changes E-U even when theta is constant.
        // This is a complete two-point TAB1, not a truncated malformed source.
        let repeated = changed(&text, 13, 5, &format!("{:11}", 2));
        let repeated = changed(&repeated, 14, 0, &format!("{:11}", 2));
        let repeated = changed(&repeated, 14, 1, &format!("{:11}", 2));
        let repeated = changed(&repeated, 15, 2, &real(1.0));
        let repeated = changed(&repeated, 15, 3, &real(2.0));
        assert!(matches!(
            parse_spectrum(&repeated, FissionNeutronKind::Prompt),
            Err(EndfError::Unsupported(_))
        ));
        for bounded in [
            SpontaneousNeutronReadLimits {
                max_incident_points: 1,
                ..SpontaneousNeutronReadLimits::default()
            },
            SpontaneousNeutronReadLimits {
                max_interpolation_regions: 1,
                ..SpontaneousNeutronReadLimits::default()
            },
        ] {
            assert!(matches!(
                SpontaneousFissionNeutronSpectrum::parse(
                    "bounded-analytic".into(),
                    selection(),
                    FissionNeutronKind::Prompt,
                    &text,
                    bounded,
                    || Ok(())
                ),
                Err(EndfError::Physics(FissionError::Budget(_)))
            ));
        }
        let unsupported = changed(&text, 10, 3, &format!("{:11}", 11));
        assert!(matches!(
            parse_spectrum(&unsupported, FissionNeutronKind::Prompt),
            Err(EndfError::Unsupported(_))
        ));
    }
}

#[test]
fn g4_final_cancellation_does_not_publish_a_partial_source() {
    for text in [
        prompt_fixture(selection(), 2.25, &[(1.0, triangle())]),
        prompt_general_evaporation_fixture(
            selection(),
            2.25,
            1.0,
            2.0,
            -3.0,
            &[(0.5, 0.0), (1.0, 2.0), (1.5, 0.0)],
            2,
        ),
        prompt_analytic_fixture(selection(), 2.25, 1.0, 0.0, 2.0, -3.0, 7),
        prompt_analytic_fixture(selection(), 2.25, 1.0, 0.0, 2.0, -3.0, 9),
    ] {
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
}
