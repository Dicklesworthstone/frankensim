//! Manufactured ENDF records and analytic laws; no evaluated isotope data.

use super::*;

mod absorption;
mod neutron_production;

fn selection(temperature_k: f64, reaction: NeutronReaction) -> EndfSelection {
    EndfSelection {
        material_number: 9237,
        nuclide: NuclideId::new(92, 235, 0).unwrap(),
        reaction,
        temperature_k,
    }
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let endf = format!("{mantissa}{exponent:+}");
        if endf.len() <= 11 {
            return format!("{endf:>11}");
        }
    }
    panic!("fixture number exceeds ENDF field");
}

fn record(fields: [String; 6], (mat, mf, mt): (i32, i32, i32)) -> String {
    assert!(fields.iter().all(|s| s.len() == 11));
    format!("{}{mat:4}{mf:2}{mt:3}{:5}", fields.concat(), 1)
}

fn control(c1: f64, c2: f64, integers: [i32; 4], section: (i32, i32, i32)) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        section,
    )
}

fn fixture(id: EndfSelection, lrp: i32, points: &[(f64, f64)], regions: &[(usize, i32)]) -> String {
    let mat = i32::from(id.material_number);
    let mt = reaction_mt(id.reaction);
    let za = f64::from(id.nuclide.atomic_number()) * 1000.0 + f64::from(id.nuclide.mass_number());
    let isomer = i32::from(id.nuclide.isomer());
    let mut lines = vec![
        control(za, 233.0, [lrp, 1, 0, 7], (mat, 1, 451)),
        control(0.0, 0.0, [isomer, isomer, 0, 6], (mat, 1, 451)),
        control(1.0, 2.0e7, [3, 0, 10, 8], (mat, 1, 451)),
        control(id.temperature_k, 0.0, [1, 0, 5, 2], (mat, 1, 451)),
    ];
    for description in [
        "MANUFACTURED FORMAT FIXTURE",
        "NO EVALUATED COEFFICIENTS",
        "TEST ANALYTIC POINTWISE LAW",
        "NO ACCURACY CLAIM",
        "NO PROCESSING CLAIM",
    ] {
        lines.push(format!("{description:<66}{mat:4}{:2}{:3}{:5}", 1, 451, 1));
    }
    lines.push(control(0.0, 0.0, [1, 451, 11, 0], (mat, 1, 451)));
    let section_count = 2 + regions.len().div_ceil(3) + points.len().div_ceil(3);
    lines.push(control(
        0.0,
        0.0,
        [3, mt, section_count as i32, 0],
        (mat, 1, 451),
    ));
    lines.push(control(0.0, 0.0, [0; 4], (mat, 1, 0)));
    lines.push(control(0.0, 0.0, [0; 4], (mat, 0, 0)));
    lines.push(control(za, 233.0, [0; 4], (mat, 3, mt)));
    lines.push(control(
        2.0e8,
        -1.0e6,
        [0, 0, regions.len() as i32, points.len() as i32],
        (mat, 3, mt),
    ));
    for chunk in regions.chunks(3) {
        let mut fields = std::array::from_fn(|_| format!("{:11}", 0));
        for (index, &(end, law)) in chunk.iter().enumerate() {
            fields[index * 2] = format!("{end:11}");
            fields[index * 2 + 1] = format!("{law:11}");
        }
        lines.push(record(fields, (mat, 3, mt)));
    }
    for chunk in points.chunks(3) {
        let mut fields = std::array::from_fn(|_| real(0.0));
        for (index, &(energy, sigma)) in chunk.iter().enumerate() {
            fields[index * 2] = real(energy);
            fields[index * 2 + 1] = real(sigma);
        }
        lines.push(record(fields, (mat, 3, mt)));
    }
    lines.push(control(0.0, 0.0, [0; 4], (mat, 3, 0)));
    lines.push(control(0.0, 0.0, [0; 4], (mat, 0, 0)));
    lines.push(control(0.0, 0.0, [0; 4], (0, 0, 0)));
    lines.join("\n") + "\n"
}

fn parse(source: &str, id: EndfSelection) -> Result<PointwiseCrossSection, EndfError> {
    PointwiseCrossSection::parse(
        "manufactured-processed-source".into(),
        id,
        source,
        EndfReadLimits::default(),
        || Ok(()),
    )
}

fn groups(edges_ev: &[f64]) -> Arc<NeutronGroups> {
    Arc::new(
        NeutronGroups::new(
            "manufactured-energy-grid".into(),
            edges_ev.iter().map(|&e| e * EV_J).collect(),
        )
        .unwrap(),
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs().max(1.0e-300) * 3.0e-14,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn fixed_width_numbers_accept_endf_implicit_exponents_and_refuse_numeric_loss() {
    for (input, expected) in [
        (" 1.234567+5", 123456.7),
        ("-1.234567-5", -0.00001234567),
        (" 9.223500+4", 92235.0),
        (" 1.25 E + 2", 125.0),
        ("           ", 0.0),
        ("-0.00000+99", -0.0),
    ] {
        assert_eq!(number(input), Some(expected));
    }
    for input in ["1.0000D+02", "NaN", "1.0-999", "1.0+999", "1-+2", "1.2.3"] {
        assert_eq!(number(input), None, "{input}");
    }
    assert_eq!(integer(" 1 2 3 "), Some(123));
    assert_eq!(integer("2147483648"), None);
    assert_eq!(
        convert_si(1e-310, BARN_M2),
        Err(EndfError::Physics(FissionError::Arithmetic))
    );
}

#[test]
fn energy_interpolation_retains_near_endpoint_tails_and_refuses_underflow() {
    let id = selection(300.0, NeutronReaction::Fission);
    let curve = parse(&fixture(id, 2, &[(0.25, 1e20), (1.0, 0.0)], &[(2, 2)]), id).unwrap();
    let upper = curve.energies_j()[1];
    let probe = f64::from_bits(upper.to_bits() - 1);
    let expected =
        curve.cross_sections_m2()[0] * ((upper - probe) / (upper - curve.energies_j()[0]));
    close(curve.cross_section_at_energy(probe).unwrap(), expected);

    let tiny = parse(&fixture(id, 2, &[(0.0, 0.0), (1.0, 3e-296)], &[(2, 2)]), id).unwrap();
    assert!(tiny.cross_sections_m2()[1] > 0.0);
    assert_eq!(
        tiny.cross_section_at_energy(0.25 * EV_J),
        Err(EndfError::Physics(FissionError::Arithmetic))
    );
    assert_eq!(
        tiny.collapse_flat_in_energy(&groups(&[0.0, 1.0]), EndfCollapseLimits::default(), || Ok(
            ()
        )),
        Err(EndfError::Physics(FissionError::Arithmetic))
    );
}

#[test]
fn metadata_material_temperature_reaction_and_isomer_stay_bound() {
    let mut selected = selection(600.0, NeutronReaction::Capture);
    selected.nuclide = NuclideId::new(92, 235, 1).unwrap();
    let mut other = selected;
    other.material_number = 9238;
    other.nuclide = NuclideId::new(92, 238, 0).unwrap();
    let source = fixture(other, 2, &[(0.0, 99.0), (4.0, 99.0)], &[(2, 2)])
        + &fixture(
            EndfSelection {
                temperature_k: 300.0,
                ..selected
            },
            2,
            &[(0.0, 7.0), (4.0, 7.0)],
            &[(2, 2)],
        )
        + &fixture(selected, 2, &[(0.0, 2.0), (4.0, 2.0)], &[(2, 2)]);
    let curve = parse(&source.replace('\n', "\r\n"), selected).unwrap();
    assert_eq!(curve.metadata().nuclide, selected.nuclide);
    assert_eq!(curve.metadata().temperature_k, 600.0);
    assert_eq!(curve.metadata().material_number, 9237);
    assert_eq!(curve.metadata().library_version, 8);
    assert_eq!(curve.metadata().library_release, 3);
    assert_eq!(curve.metadata().material_revision, 7);
    assert_eq!(curve.reaction(), NeutronReaction::Capture);
    close(
        curve.cross_section_at_energy(2.0 * EV_J).unwrap(),
        2.0 * BARN_M2,
    );
    close(curve.mass_difference_q_j(), 2.0e8 * EV_J);
    close(curve.reaction_q_j(), -1.0e6 * EV_J);
    assert_eq!(
        parse(
            &source,
            EndfSelection {
                temperature_k: 450.0,
                ..selected
            }
        ),
        Err(EndfError::MissingSelection)
    );
    assert!(matches!(
        parse(
            &source,
            EndfSelection {
                nuclide: NuclideId::new(92, 235, 0).unwrap(),
                ..selected
            }
        ),
        Err(EndfError::Format { .. })
    ));
    assert_eq!(
        parse(&(source.clone() + &source), selected),
        Err(EndfError::AmbiguousSelection)
    );
}

#[test]
fn collapse_integrates_thresholds_mixed_laws_and_discontinuities() {
    let id = selection(300.0, NeutronReaction::Fission);
    // Histogram zero on [0,1), then linear to a left limit 4 at E=2;
    // jump to 10 at E=2 and linearly fall to zero at E=4.
    let source = fixture(
        id,
        2,
        &[(0.0, 0.0), (1.0, 0.0), (2.0, 4.0), (2.0, 10.0), (4.0, 0.0)],
        &[(2, 1), (5, 2)],
    );
    let curve = parse(&source, id).unwrap();
    assert_eq!(curve.cross_section_at_energy(0.5 * EV_J).unwrap(), 0.0);
    close(
        curve.cross_section_at_energy(1.5 * EV_J).unwrap(),
        2.0 * BARN_M2,
    );
    close(
        curve.cross_section_at_energy(2.0 * EV_J).unwrap(),
        10.0 * BARN_M2,
    );
    let values = curve
        .collapse_flat_in_energy(
            &groups(&[0.0, 1.5, 2.0, 3.0, 4.0]),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    for (&value, barns) in values.iter().zip([1.0 / 3.0, 3.0, 7.5, 2.5]) {
        close(value, barns * BARN_M2);
    }
    let whole = curve
        .collapse_flat_in_energy(&groups(&[0.0, 4.0]), EndfCollapseLimits::default(), || {
            Ok(())
        })
        .unwrap();
    close(whole[0], 3.0 * BARN_M2);
    assert!(curve.cross_section_at_energy(5.0 * EV_J).is_err());
    assert!(
        curve
            .collapse_flat_in_energy(&groups(&[0.0, 5.0]), EndfCollapseLimits::default(), || Ok(
                ()
            ))
            .is_err()
    );
}

#[test]
fn histogram_uses_left_interval_value_and_final_knot_value() {
    let id = selection(300.0, NeutronReaction::Capture);
    let curve = parse(
        &fixture(id, 0, &[(0.0, 2.0), (1.0, 8.0), (4.0, 100.0)], &[(3, 1)]),
        id,
    )
    .unwrap();
    close(
        curve.cross_section_at_energy(1.0 * EV_J).unwrap(),
        8.0 * BARN_M2,
    );
    close(
        curve.cross_section_at_energy(4.0 * EV_J).unwrap(),
        100.0 * BARN_M2,
    );
    let row = curve
        .collapse_flat_in_energy(&groups(&[0.0, 4.0]), EndfCollapseLimits::default(), || {
            Ok(())
        })
        .unwrap();
    close(row[0], 6.5 * BARN_M2);
}

#[test]
fn backgrounds_unsupported_laws_and_malformed_sections_are_refused() {
    let id = selection(300.0, NeutronReaction::Fission);
    let points = [(0.0, 0.0), (1.0, 2.0), (4.0, 3.0)];
    assert!(matches!(
        parse(&fixture(id, 1, &points, &[(3, 2)]), id),
        Err(EndfError::Unsupported(_))
    ));
    assert!(matches!(
        parse(&fixture(id, 2, &points, &[(3, 5)]), id),
        Err(EndfError::Unsupported(_))
    ));
    assert!(matches!(
        parse(&fixture(id, 2, &points, &[(2, 2)]), id),
        Err(EndfError::Format { .. })
    ));
    assert!(parse(&fixture(id, 2, &[(1.0, 0.0), (0.0, 2.0)], &[(2, 2)]), id).is_err());
    assert!(parse(&fixture(id, 2, &[(0.0, -1.0), (1.0, 2.0)], &[(2, 2)]), id).is_err());
    assert!(
        parse(
            &fixture(
                id,
                2,
                &[(0.0, 0.0), (1.0, 1.0), (1.0, 2.0), (1.0, 3.0)],
                &[(4, 2)]
            ),
            id
        )
        .is_err()
    );
    let source = fixture(id, 2, &points, &[(3, 2)]);
    let mut lines: Vec<String> = source.lines().map(str::to_owned).collect();
    // Pointwise payload must not switch to another material inside the TAB1.
    lines[16].replace_range(66..70, "9238");
    assert!(matches!(
        parse(&lines.join("\n"), id),
        Err(EndfError::Format { .. })
    ));
    lines = source.lines().map(str::to_owned).collect();
    // MF3 HEAD must carry the same isotope as MF1, even if MAT agrees.
    lines[13].replace_range(0..11, &real(92238.0));
    assert!(matches!(
        parse(&lines.join("\n"), id),
        Err(EndfError::Format { .. })
    ));
    assert!(matches!(
        parse("short record", id),
        Err(EndfError::Format { .. })
    ));
    assert!(matches!(
        parse(&source.lines().take(17).collect::<Vec<_>>().join("\n"), id),
        Err(EndfError::Format { .. })
    ));
}

#[test]
fn point_record_region_and_cancellation_budgets_apply_before_publication() {
    let id = selection(300.0, NeutronReaction::Fission);
    let source = fixture(id, 2, &[(0.0, 1.0), (4.0, 1.0)], &[(2, 2)]);
    for limits in [
        EndfReadLimits {
            max_bytes: source.len() - 1,
            ..EndfReadLimits::default()
        },
        EndfReadLimits {
            max_records: 10,
            ..EndfReadLimits::default()
        },
        EndfReadLimits {
            max_points: 1,
            ..EndfReadLimits::default()
        },
        EndfReadLimits {
            max_interpolation_regions: 0,
            ..EndfReadLimits::default()
        },
    ] {
        assert!(matches!(
            PointwiseCrossSection::parse("source".into(), id, &source, limits, || Ok(())),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    let mut calls = 0;
    let result = PointwiseCrossSection::parse(
        "source".into(),
        id,
        &source,
        EndfReadLimits::default(),
        || {
            calls += 1;
            if calls == 7 {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(result, Err(EndfError::Physics(FissionError::Cancelled)));
    let curve = parse(&source, id).unwrap();
    assert!(matches!(
        curve.collapse_flat_in_energy(
            &groups(&[0.0, 4.0]),
            EndfCollapseLimits {
                max_segment_visits: 0,
                ..EndfCollapseLimits::default()
            },
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    assert_eq!(
        curve.collapse_flat_in_energy(&groups(&[0.0, 4.0]), EndfCollapseLimits::default(), || Err(
            FissionError::Cancelled
        )),
        Err(EndfError::Physics(FissionError::Cancelled))
    );
}

#[test]
fn temperature_assembly_preserves_domain_and_shares_integration_budget() {
    let curves: Vec<_> = [300.0, 600.0]
        .into_iter()
        .map(|temperature| {
            let id = selection(temperature, NeutronReaction::Fission);
            parse(
                &fixture(
                    id,
                    2,
                    &[(0.0, temperature / 300.0), (4.0, temperature / 300.0)],
                    &[(2, 2)],
                ),
                id,
            )
            .unwrap()
        })
        .collect();
    let table = collapse_temperature_table(
        "manufactured-flat-energy-table".into(),
        &curves,
        groups(&[0.0, 4.0]),
        EndfCollapseLimits::default(),
        CrossSectionLimits::default(),
        || Ok(()),
    )
    .unwrap();
    close(
        table.cross_section_at_group(450.0, 0).unwrap(),
        1.5 * BARN_M2,
    );
    assert!(table.cross_section_at_group(700.0, 0).is_err());
    assert!(matches!(
        collapse_temperature_table(
            "source".into(),
            &curves,
            groups(&[0.0, 4.0]),
            EndfCollapseLimits {
                max_segment_visits: 1,
                ..EndfCollapseLimits::default()
            },
            CrossSectionLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    assert!(
        collapse_temperature_table(
            "source".into(),
            &[curves[1].clone(), curves[0].clone()],
            groups(&[0.0, 4.0]),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn imported_fission_and_capture_tables_drive_competing_burnup_events() {
    use crate::EnergyChannels;
    use crate::burnup::{
        BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
        BurnupProduct, BurnupRate, BurnupReaction,
    };
    use crate::irradiation::GroupFlux;
    use fs_qty::Frequency;
    use fs_thermochem::decay::DecayLawV1;

    let grid = groups(&[0.0, 4.0]);
    let mut reactions = Vec::new();
    for (kind, scale, product, yield_per_event) in [
        (NeutronReaction::Fission, 1.0, 1, 2.0),
        (NeutronReaction::Capture, 0.5, 2, 1.0),
    ] {
        let curves: Vec<_> = [(300.0, scale), (600.0, 3.0 * scale)]
            .into_iter()
            .map(|(t, sigma)| {
                let id = selection(t, kind);
                parse(
                    &fixture(id, 2, &[(0.0, sigma), (4.0, sigma)], &[(2, 2)]),
                    id,
                )
                .unwrap()
            })
            .collect();
        let table = Arc::new(
            collapse_temperature_table(
                "manufactured-row-assembly".into(),
                &curves,
                Arc::clone(&grid),
                EndfCollapseLimits::default(),
                CrossSectionLimits::default(),
                || Ok(()),
            )
            .unwrap(),
        );
        reactions.push(BurnupReaction {
            data_id: "manufactured-energy-account".into(),
            parent: 0,
            rate: match kind {
                NeutronReaction::Fission => BurnupRate::TabulatedFission {
                    cross_sections: table,
                },
                NeutronReaction::Capture
                | NeutronReaction::ProtonProduction
                | NeutronReaction::DeuteronProduction
                | NeutronReaction::TritonProduction
                | NeutronReaction::Helium3Production
                | NeutronReaction::AlphaProduction => BurnupRate::TabulatedCapture {
                    cross_sections: table,
                },
                NeutronReaction::TwoNeutronProduction
                | NeutronReaction::ThreeNeutronProduction
                | NeutronReaction::FourNeutronProduction => {
                    BurnupRate::TabulatedNeutronProduction {
                        cross_sections: table,
                    }
                }
            },
            products: vec![BurnupProduct {
                nuclide: product,
                yield_per_event,
            }],
            energy_per_event: EnergyChannels {
                deposited_j: 1e-12,
                radiation_j: 0.0,
                escaped_j: 0.0,
            },
        });
    }
    let limits = BurnupLimits::default();
    let stable = DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap();
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-import-to-reaction".into(),
                groups: Arc::clone(&grid),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![
                    BurnupNuclide {
                        id: selection(300.0, NeutronReaction::Fission).nuclide,
                        law: stable,
                    },
                    BurnupNuclide {
                        id: NuclideId::new(40, 95, 0).unwrap(),
                        law: stable,
                    },
                    BurnupNuclide {
                        id: NuclideId::new(92, 236, 0).unwrap(),
                        law: stable,
                    },
                ],
                reactions,
            },
            limits,
        )
        .unwrap(),
    );
    let mut inventory = BurnupInventory::new(model, vec![100.0, 0.0, 0.0]).unwrap();
    let flux = GroupFlux::new(grid, vec![1e27]).unwrap();
    let report = inventory
        .advance(
            1.0,
            &flux,
            450.0,
            BurnupBudget {
                limits,
                nuclei_tolerance: 1e-9,
                energy_tolerance_j: 1e-18,
            },
            || Ok(()),
        )
        .unwrap();
    let removals = 100.0 * -(-0.3_f64).exp_m1();
    close(report.fissions, removals * 2.0 / 3.0);
    close(report.captures, removals / 3.0);
    close(inventory.populations()[0], 100.0 * (-0.3_f64).exp());
    close(inventory.populations()[1], 2.0 * report.fissions);
    close(inventory.populations()[2], report.captures);
    // Imported QM/QI have not replaced the explicit prompt-energy owner.
    close(report.released.deposited_j, removals * 1e-12);
}
