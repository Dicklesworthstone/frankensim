//! Manufactured format and equation fixtures, not evaluated isotope data.

use super::*;

pub(super) fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("manufactured value exceeds ENDF field width");
}

pub(super) fn record(fields: [String; 6], mt: i32) -> String {
    assert!(fields.iter().all(|field| field.len() == 11));
    format!("{}{:4}{:2}{mt:3}{:5}", fields.concat(), 101, 1, 1)
}

pub(super) fn control(c1: f64, c2: f64, integers: [i32; 4], mt: i32) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        mt,
    )
}

pub(super) fn selection() -> EndfFissionEnergySelection {
    EndfFissionEnergySelection {
        material_number: 101,
        nuclide: NuclideId::new(10, 20, 1).unwrap(),
        temperature_k: 0.0,
    }
}

pub(super) fn coefficients(values: [f64; 9]) -> [FissionEnergyCoefficient; 9] {
    values.map(|value| FissionEnergyCoefficient {
        value,
        uncertainty: 0.0,
    })
}

#[derive(Clone)]
pub(super) struct TableFixture {
    pub component: FissionEnergyComponent,
    pub derivation: u8,
    /// Original incident energies and component energies, both numerically eV.
    pub points: Vec<(f64, f64)>,
    pub regions: Vec<(usize, i32)>,
}

pub(super) fn fixture(orders: &[[FissionEnergyCoefficient; 9]], tables: &[TableFixture]) -> String {
    assert!(!orders.is_empty());
    let mut lines = vec![
        control(10020.0, 19.8, [1, 1, 0, 3], 451),
        control(1.0, 0.0, [1, 1, 0, 6], 451),
        control(1.0, 2.0e7, [2, 0, 10, 8], 451),
        control(0.0, 0.0, [0; 4], 451),
        control(0.0, 0.0, [0; 4], 0),
        control(
            10020.0,
            19.8,
            [0, i32::from(!tables.is_empty()), 0, tables.len() as i32],
            458,
        ),
        control(
            0.0,
            0.0,
            [
                0,
                orders.len() as i32 - 1,
                18 * orders.len() as i32,
                9 * orders.len() as i32,
            ],
            458,
        ),
    ];
    for order in orders {
        for row in order.chunks(3) {
            lines.push(record(
                std::array::from_fn(|index| {
                    let coefficient = row[index / 2];
                    real(if index % 2 == 0 {
                        coefficient.value
                    } else {
                        coefficient.uncertainty
                    })
                }),
                458,
            ));
        }
    }
    for table in tables {
        lines.push(control(
            0.0,
            0.0,
            [
                i32::from(table.derivation),
                table.component.index() as i32 + 1,
                table.regions.len() as i32,
                table.points.len() as i32,
            ],
            458,
        ));
        for chunk in table.regions.chunks(3) {
            let mut fields = std::array::from_fn(|_| format!("{:11}", 0));
            for (index, &(end, interpolation)) in chunk.iter().enumerate() {
                fields[2 * index] = format!("{end:11}");
                fields[2 * index + 1] = format!("{interpolation:11}");
            }
            lines.push(record(fields, 458));
        }
        for chunk in table.points.chunks(3) {
            let mut fields = std::array::from_fn(|_| real(0.0));
            for (index, &(incident, energy)) in chunk.iter().enumerate() {
                fields[2 * index] = real(incident);
                fields[2 * index + 1] = real(energy);
            }
            lines.push(record(fields, 458));
        }
    }
    lines.push(control(0.0, 0.0, [0; 4], 0));
    lines.join("\n") + "\n"
}

pub(super) fn parse(source: &str) -> Result<FissionEnergyRelease, EndfError> {
    FissionEnergyRelease::parse(
        "manufactured-energy-release".into(),
        selection(),
        source,
        EndfFissionEnergyReadLimits::default(),
        || Ok(()),
    )
}

pub(super) fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(EV_J),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn replace_line(source: &str, index: usize, replacement: String) -> String {
    let mut lines: Vec<_> = source.lines().map(str::to_owned).collect();
    lines[index] = replacement;
    lines.join("\n") + "\n"
}

fn table() -> TableFixture {
    TableFixture {
        component: FissionEnergyComponent::Fragments,
        derivation: 2,
        points: vec![(1.0, 10.0), (2.0, 20.0), (2.0, 30.0), (4.0, 50.0)],
        regions: vec![(4, 2)],
    }
}

#[test]
fn g0_reference_components_and_uncertainties_keep_their_original_meanings() {
    let reference = std::array::from_fn(|index| FissionEnergyCoefficient {
        value: (index + 1) as f64 * 1.0e6,
        uncertainty: (index + 1) as f64 * 100.0,
    });
    let law = parse(&fixture(&[reference], &[])).unwrap();
    assert_eq!(law.data_id(), "manufactured-energy-release");
    assert_eq!(law.polynomial_order(), 0);
    assert_eq!(law.metadata().nuclide, selection().nuclide);
    assert_eq!(law.metadata().temperature_k, 0.0);
    assert_eq!(law.metadata().resonance_flag, 1);
    for (index, component) in FissionEnergyComponent::ALL.into_iter().enumerate() {
        assert_eq!(component.index(), index);
        assert_eq!(law.coefficients(component), &[reference[index]]);
        near(
            law.reference_value_j(component),
            reference[index].value * EV_J,
        );
        assert!(law.tabulation(component).is_none());
        for energy in [0.0, EV_J, 1.0e6 * EV_J] {
            assert!(matches!(
                law.value_at_energy_j(component, energy),
                Err(EndfError::Unsupported(_))
            ));
        }
    }
    // Deliberately unequal component sums and source totals remain unchanged;
    // neither an energy destination nor a normalization is inferred here.
    near(
        law.reference_value_j(FissionEnergyComponent::Recoverable),
        8.0e6 * EV_J,
    );
    near(
        law.reference_value_j(FissionEnergyComponent::Total),
        9.0e6 * EV_J,
    );
}

#[test]
fn g1_signed_polynomials_are_evaluated_in_source_units_without_unit_repair() {
    let reference = coefficients(std::array::from_fn(|i| 100.0 * (i + 1) as f64));
    let linear = std::array::from_fn(|i| FissionEnergyCoefficient {
        value: -(i as f64 + 1.0),
        uncertainty: 0.125 * (i + 1) as f64,
    });
    let quadratic = coefficients(std::array::from_fn(|i| 0.25 * (i + 1) as f64));
    let law = parse(&fixture(&[reference, linear, quadratic], &[])).unwrap();
    assert_eq!(law.polynomial_order(), 2);
    for component in FissionEnergyComponent::ALL {
        let i = component.index();
        assert_eq!(
            law.coefficients(component),
            &[reference[i], linear[i], quadratic[i]]
        );
        for incident_ev in [0.0, 2.0, 3.0, 10.0] {
            let expected_ev =
                (i + 1) as f64 * (100.0 - incident_ev + 0.25 * incident_ev * incident_ev);
            near(
                law.value_at_energy_j(component, incident_ev * EV_J)
                    .unwrap(),
                expected_ev * EV_J,
            );
        }
    }
    // A large native quadratic is retained exactly; no library-specific
    // heuristic divides it by 1e6 as a speculative historical unit correction.
    assert_eq!(
        law.coefficients(FissionEnergyComponent::Fragments)[2].value,
        0.25
    );
}

#[test]
fn g1_component_tables_override_references_and_retain_jumps_and_histograms() {
    let mut photons = table();
    photons.component = FissionEnergyComponent::PromptPhotons;
    photons.derivation = 1;
    photons.points = vec![(1.0, 3.0), (4.0, 6.0)];
    photons.regions = vec![(2, 1)];
    let law = parse(&fixture(&[coefficients([100.0; 9])], &[table(), photons])).unwrap();
    let fragments = law.tabulation(FissionEnergyComponent::Fragments).unwrap();
    assert_eq!(fragments.derivation(), 2);
    assert_eq!(
        fragments.energies_j(),
        &[EV_J, 2.0 * EV_J, 2.0 * EV_J, 4.0 * EV_J]
    );
    for (energy_ev, value_ev) in [
        (1.0, 10.0),
        (1.5, 15.0),
        (2.0, 30.0),
        (3.0, 40.0),
        (4.0, 50.0),
    ] {
        near(
            law.value_at_energy_j(FissionEnergyComponent::Fragments, energy_ev * EV_J)
                .unwrap(),
            value_ev * EV_J,
        );
    }
    let photons = law
        .tabulation(FissionEnergyComponent::PromptPhotons)
        .unwrap();
    assert_eq!(photons.derivation(), 1);
    near(photons.value_at_energy_j(3.0 * EV_J).unwrap(), 3.0 * EV_J);
    near(photons.value_at_energy_j(4.0 * EV_J).unwrap(), 6.0 * EV_J);
    near(
        law.reference_value_j(FissionEnergyComponent::Fragments),
        100.0 * EV_J,
    );
    assert!(
        law.value_at_energy_j(FissionEnergyComponent::PromptNeutrons, 2.0 * EV_J)
            .is_err()
    );
    assert!(fragments.value_at_energy_j(0.5 * EV_J).is_err());
    assert!(fragments.value_at_energy_j(5.0 * EV_J).is_err());
}

#[test]
fn g0_selection_is_exact_and_ambiguous_or_truncated_records_refuse() {
    let source = fixture(&[coefficients([100.0; 9])], &[]);
    for changed in [
        EndfFissionEnergySelection {
            temperature_k: 300.0,
            ..selection()
        },
        EndfFissionEnergySelection {
            material_number: 102,
            ..selection()
        },
        EndfFissionEnergySelection {
            nuclide: NuclideId::new(10, 20, 0).unwrap(),
            ..selection()
        },
    ] {
        assert!(
            FissionEnergyRelease::parse(
                "fixture".into(),
                changed,
                &source,
                EndfFissionEnergyReadLimits::default(),
                || Ok(())
            )
            .is_err()
        );
    }
    assert!(matches!(
        parse(&(source.clone() + &source)),
        Err(EndfError::AmbiguousSelection)
    ));
    let truncated = source.lines().take(9).collect::<Vec<_>>().join("\n");
    assert!(parse(&truncated).is_err());
    let last = source.lines().count() - 1;
    assert!(parse(&replace_line(&source, last, control(1.0, 0.0, [0; 4], 0))).is_err());
    assert!(
        parse(&replace_line(
            &source,
            5,
            control(10021.0, 19.8, [0; 4], 458)
        ))
        .is_err()
    );
    assert!(
        parse(&replace_line(
            &source,
            5,
            control(10020.0, 20.0, [0; 4], 458)
        ))
        .is_err()
    );
}

#[test]
fn g0_malformed_component_shapes_flags_and_uncertainties_refuse() {
    let reference = coefficients([100.0; 9]);
    let source = fixture(&[reference], &[]);
    for head in [
        control(10020.0, 19.8, [0, 2, 0, 0], 458),
        control(10020.0, 19.8, [1, 0, 0, 0], 458),
        control(10020.0, 19.8, [0, 0, 0, 1], 458),
        control(10020.0, 19.8, [0, 1, 0, 0], 458),
    ] {
        assert!(parse(&replace_line(&source, 5, head)).is_err());
    }
    for list in [
        control(0.0, 0.0, [0, 0, 17, 9], 458),
        control(0.0, 0.0, [0, 0, 18, 8], 458),
        control(1.0, 0.0, [0, 0, 18, 9], 458),
    ] {
        assert!(parse(&replace_line(&source, 6, list)).is_err());
    }
    let mut invalid = reference;
    invalid[0].value = -1.0;
    assert!(parse(&fixture(&[invalid], &[])).is_err());
    invalid = reference;
    invalid[4].uncertainty = -1.0;
    assert!(parse(&fixture(&[reference, invalid], &[])).is_err());
    assert!(parse(&fixture(&[reference, coefficients([0.0; 9])], &[table()])).is_err());
}

#[test]
fn g0_duplicate_ifc_invalid_derivation_and_unsupported_table_laws_refuse() {
    let reference = coefficients([100.0; 9]);
    assert!(parse(&fixture(&[reference], &[table(), table()])).is_err());
    for derivation in [0, 3] {
        let mut invalid = table();
        invalid.derivation = derivation;
        assert!(parse(&fixture(&[reference], &[invalid])).is_err());
    }
    let source = fixture(&[reference], &[table()]);
    for ifc in [0, 10] {
        assert!(
            parse(&replace_line(
                &source,
                10,
                control(0.0, 0.0, [2, ifc, 1, 4], 458)
            ))
            .is_err()
        );
    }
    let mut invalid = table();
    invalid.regions = vec![(4, 5)];
    assert!(matches!(
        parse(&fixture(&[reference], &[invalid])),
        Err(EndfError::Unsupported(_))
    ));
    let mut invalid = table();
    invalid.points[0].1 = -1.0;
    assert!(parse(&fixture(&[reference], &[invalid])).is_err());
}

#[test]
fn g4_input_storage_and_work_limits_apply_before_retaining_the_source() {
    let source = fixture(&[coefficients([100.0; 9])], &[table()]);
    for limits in [
        EndfFissionEnergyReadLimits {
            records: EndfReadLimits {
                max_bytes: source.len() - 1,
                ..EndfReadLimits::default()
            },
            ..EndfFissionEnergyReadLimits::default()
        },
        EndfFissionEnergyReadLimits {
            records: EndfReadLimits {
                max_records: 9,
                ..EndfReadLimits::default()
            },
            ..EndfFissionEnergyReadLimits::default()
        },
        EndfFissionEnergyReadLimits {
            records: EndfReadLimits {
                max_points: 3,
                ..EndfReadLimits::default()
            },
            ..EndfFissionEnergyReadLimits::default()
        },
        EndfFissionEnergyReadLimits {
            records: EndfReadLimits {
                max_interpolation_regions: 0,
                ..EndfReadLimits::default()
            },
            ..EndfFissionEnergyReadLimits::default()
        },
        EndfFissionEnergyReadLimits {
            max_table_values: 7,
            ..EndfFissionEnergyReadLimits::default()
        },
    ] {
        assert!(
            FissionEnergyRelease::parse("fixture".into(), selection(), &source, limits, || Ok(()))
                .is_err()
        );
    }
    let linear = fixture(&[coefficients([100.0; 9]), coefficients([-1.0; 9])], &[]);
    assert!(
        FissionEnergyRelease::parse(
            "fixture".into(),
            selection(),
            &linear,
            EndfFissionEnergyReadLimits {
                max_polynomial_order: 0,
                ..EndfFissionEnergyReadLimits::default()
            },
            || Ok(())
        )
        .is_err()
    );
    let mut second = table();
    second.component = FissionEnergyComponent::PromptPhotons;
    let shared = fixture(&[coefficients([100.0; 9])], &[table(), second]);
    assert!(
        FissionEnergyRelease::parse(
            "fixture".into(),
            selection(),
            &shared,
            EndfFissionEnergyReadLimits {
                max_table_values: 15,
                ..EndfFissionEnergyReadLimits::default()
            },
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn g4_callbacks_include_early_mid_and_final_publication_boundaries() {
    let source = fixture(&[coefficients([100.0; 9])], &[table()]);
    let limits = EndfFissionEnergyReadLimits::default();
    let mut polls = 0;
    FissionEnergyRelease::parse("fixture".into(), selection(), &source, limits, || {
        polls += 1;
        Ok(())
    })
    .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut calls = 0;
        let result =
            FissionEnergyRelease::parse("fixture".into(), selection(), &source, limits, || {
                calls += 1;
                if calls == cutoff {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            });
        assert!(matches!(
            result,
            Err(EndfError::Physics(FissionError::Cancelled))
        ));
    }
}

#[test]
fn g0_negative_nonfinite_and_unrepresentable_physical_outputs_refuse() {
    let law = parse(&fixture(
        &[coefficients([1.0; 9]), coefficients([-2.0; 9])],
        &[],
    ))
    .unwrap();
    for energy in [-EV_J, f64::NAN, f64::INFINITY, EV_J] {
        assert!(
            law.value_at_energy_j(FissionEnergyComponent::Fragments, energy)
                .is_err()
        );
    }
    let overflow = parse(&fixture(
        &[coefficients([1.0; 9]), coefficients([1e308; 9])],
        &[],
    ))
    .unwrap();
    assert!(
        overflow
            .value_at_energy_j(FissionEnergyComponent::Fragments, 2.0 * EV_J)
            .is_err()
    );
    assert!(parse(&fixture(&[coefficients([1e-310; 9])], &[])).is_err());
    let mut too_small = table();
    too_small.points[0].1 = 1e-310;
    assert!(parse(&fixture(&[coefficients([1.0; 9])], &[too_small])).is_err());
}
