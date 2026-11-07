//! Manufactured coefficient rows test interpolation, units and admission;
//! these values are not evaluated cross sections for the named nuclide.

use super::*;

fn limits() -> CrossSectionLimits {
    CrossSectionLimits {
        max_temperatures: 8,
        max_neutron_groups: 8,
        max_values: 64,
    }
}

fn data(temperatures: Vec<f64>, values: Vec<Vec<f64>>) -> CrossSectionData {
    let groups = values.first().map_or(3, Vec::len);
    CrossSectionData {
        data_id: "manufactured-temperature-data-v1".into(),
        nuclide: NuclideId::new(92, 235, 0).unwrap(),
        reaction: NeutronReaction::Fission,
        groups: Arc::new(
            NeutronGroups::new(
                "manufactured-grid-v1".into(),
                (0..=groups).map(|i| i as f64).collect(),
            )
            .unwrap(),
        ),
        temperatures_k: temperatures,
        cross_sections_m2: values,
    }
}

fn fixture() -> CrossSectionData {
    data(
        vec![300.0, 600.0, 900.0],
        vec![
            vec![2.0, 8.0, 0.0],
            vec![4.0, 3.0, 1.0],
            vec![10.0, 1.0, 5.0],
        ],
    )
}

#[test]
fn g3_temperature_components_recover_original_rates_with_shared_storage() {
    let table = CrossSectionTable::new(fixture(), limits()).unwrap();
    let mut parts = Vec::new();
    for row in 0..3 {
        for group in 0..3 {
            let part = table
                .temperature_component(format!("knot {row}, group {group}"), row, group)
                .unwrap();
            assert!(Arc::ptr_eq(&part.values_m2, &table.values_m2));
            assert!(Arc::ptr_eq(&part.temperatures_k, &table.temperatures_k));
            assert!(Arc::ptr_eq(&part.groups, &table.groups));
            for knot in 0..3 {
                for output_group in 0..3 {
                    let expected = if knot == row && output_group == group {
                        fixture().cross_sections_m2[row][group]
                    } else {
                        0.0
                    };
                    assert_eq!(
                        part.cross_section_at_group(table.temperatures_k()[knot], output_group)
                            .unwrap(),
                        expected
                    );
                }
            }
            parts.push(part);
        }
    }
    let incident = flux(&table, vec![5.0, 7.0, 11.0]);
    for temperature in [300.0, 379.0, 450.0, 600.0, 742.0, 750.0, 900.0] {
        let original = table
            .rate_constant(&incident, temperature, 3, || Ok(()))
            .unwrap();
        let mut recovered = Accumulator::new(0.0);
        for part in &parts {
            recovered
                .add(
                    part.rate_constant(&incident, temperature, 3, || Ok(()))
                        .unwrap(),
                )
                .unwrap();
        }
        assert!((recovered.value() - original).abs() <= 8.0 * f64::EPSILON * original);
    }
}

#[test]
fn g1_temperature_components_preserve_cross_section_weighted_observables() {
    let table = CrossSectionTable::new(
        data(
            vec![300.0, 600.0, 900.0],
            vec![vec![2.0], vec![6.0], vec![0.0]],
        ),
        limits(),
    )
    .unwrap();
    let low = table.temperature_component("low".into(), 0, 0).unwrap();
    let high = table.temperature_component("high".into(), 1, 0).unwrap();
    // At 450 K the contributions are 1 and 3: for yields 0.5 and 1.5,
    // the production coefficient is 5, not sigma(450)*mean(yield)=4.
    let a = low.cross_section_at_group(450.0, 0).unwrap();
    let b = high.cross_section_at_group(450.0, 0).unwrap();
    assert_eq!((a, b), (1.0, 3.0));
    assert_eq!(a * 0.5 + b * 1.5, 5.0);
    assert_eq!(a * 4.0 + b * 8.0, 28.0);
    // A zero final knot remains an anchor; the preceding channel decays to
    // zero there instead of extrapolating its nonzero source coefficient.
    assert_eq!(high.cross_section_at_group(750.0, 0).unwrap(), 3.0);
    assert_eq!(high.cross_section_at_group(900.0, 0).unwrap(), 0.0);
    assert_eq!(low.cross_section_at_group(750.0, 0).unwrap(), 0.0);
}

#[test]
fn g0_temperature_components_keep_domain_identity_and_positive_range_admission() {
    let table = CrossSectionTable::new(fixture(), limits()).unwrap();
    for (name, row, group) in [("", 0, 0), ("bad row", 3, 0), ("bad group", 0, 3)] {
        assert!(
            table
                .temperature_component(name.into(), row, group)
                .is_err()
        );
    }
    let part = table
        .temperature_component("selected".into(), 1, 0)
        .unwrap();
    assert!(part.temperature_component("nested".into(), 1, 0).is_err());
    for temperature in [299.0, 901.0, f64::NAN] {
        assert!(part.cross_section_at_group(temperature, 2).is_err());
        assert!(
            part.rate_constant(&flux(&table, vec![0.0; 3]), temperature, 3, || Ok(()))
                .is_err()
        );
    }
    let foreign = GroupFlux::new(
        Arc::new(NeutronGroups::new("foreign".into(), vec![0.0, 1.0, 2.0, 3.0]).unwrap()),
        vec![0.0; 3],
    )
    .unwrap();
    assert!(part.rate_constant(&foreign, 600.0, 3, || Ok(())).is_err());
    let tiny = CrossSectionTable::new(
        data(vec![1.0, 3.0], vec![vec![0.0], vec![f64::from_bits(1)]]),
        limits(),
    )
    .unwrap()
    .temperature_component("tiny".into(), 1, 0)
    .unwrap();
    assert_eq!(
        tiny.cross_section_at_group(1.0_f64.next_up(), 0),
        Err(FissionError::Arithmetic)
    );
    let wide = CrossSectionTable::new(
        data(vec![1e-300, 1e308], vec![vec![0.0], vec![1e308]]),
        limits(),
    )
    .unwrap()
    .temperature_component("wide".into(), 1, 0)
    .unwrap();
    assert_eq!(
        wide.cross_section_at_group(2e-300, 0),
        Err(FissionError::Arithmetic)
    );
}

#[test]
fn g4_temperature_component_hazards_visit_one_group_and_preserve_cancellation() {
    let table = CrossSectionTable::new(
        data(vec![300.0, 600.0], vec![vec![2.0; 256], vec![4.0; 256]]),
        CrossSectionLimits {
            max_temperatures: 2,
            max_neutron_groups: 256,
            max_values: 512,
        },
    )
    .unwrap();
    let part = table
        .temperature_component("one group".into(), 0, 173)
        .unwrap();
    let incident = flux(&table, vec![3.0; 256]);
    let mut polls = 0;
    let rate = part
        .rate_constant(&incident, 450.0, 256, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(rate, 3.0);
    assert_eq!(polls, 3, "start, selected group and final publication only");
    for stop in 1..=polls {
        let mut calls = 0;
        assert_eq!(
            part.rate_constant(&incident, 450.0, 256, || {
                calls += 1;
                if calls == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(FissionError::Cancelled)
        );
    }
    assert_eq!(
        part.rate_constant(&incident, 450.0, 255, || Ok(())),
        Err(FissionError::Budget("neutron groups"))
    );
    assert_eq!(
        table
            .rate_constant(&incident, 450.0, 256, || Ok(()))
            .unwrap(),
        2304.0
    );
}

fn flux(table: &CrossSectionTable, values: Vec<f64>) -> GroupFlux {
    GroupFlux::new(Arc::clone(table.groups()), values).unwrap()
}

#[test]
fn g1_exact_knots_piecewise_interpolation_and_absolute_reaction_rates() {
    let input = fixture();
    let table = CrossSectionTable::new(input.clone(), limits()).unwrap();
    assert_eq!(table.data_id(), input.data_id);
    assert_eq!(table.nuclide(), input.nuclide);
    assert_eq!(table.reaction(), NeutronReaction::Fission);
    assert_eq!(table.temperatures_k(), &[300.0, 600.0, 900.0]);
    assert_eq!(table.temperature_min_k(), 300.0);
    assert_eq!(table.temperature_max_k(), 900.0);
    for (row, &temperature) in input.temperatures_k.iter().enumerate() {
        for group in 0..3 {
            assert_eq!(
                table
                    .cross_section_at_group(temperature, group)
                    .unwrap()
                    .to_bits(),
                input.cross_sections_m2[row][group].to_bits()
            );
        }
    }
    let spectrum = flux(&table, vec![5.0, 7.0, 11.0]);
    for (temperature, expected) in [
        (300.0, 66.0),
        (450.0, 59.0),
        (600.0, 52.0),
        (750.0, 82.0),
        (900.0, 112.0),
    ] {
        assert_eq!(
            table
                .rate_constant(&spectrum, temperature, 3, || Ok(()))
                .unwrap(),
            expected
        );
    }
    assert_eq!(table.cross_section_at_group(450.0, 0).unwrap(), 3.0);
    assert_eq!(table.cross_section_at_group(450.0, 1).unwrap(), 5.5);
    assert_eq!(table.cross_section_at_group(750.0, 2).unwrap(), 3.0);
}

#[test]
fn g3_interpolation_is_positive_and_bounded_in_each_temperature_interval() {
    let table = CrossSectionTable::new(fixture(), limits()).unwrap();
    for interval in table.temperatures_k().windows(2) {
        for k in 0..=32 {
            let temperature = interval[0] + (interval[1] - interval[0]) * (k as f64 / 32.0);
            for group in 0..3 {
                let low = table.cross_section_at_group(interval[0], group).unwrap();
                let high = table.cross_section_at_group(interval[1], group).unwrap();
                let value = table.cross_section_at_group(temperature, group).unwrap();
                assert!((low.min(high)..=low.max(high)).contains(&value));
                assert!(value >= 0.0);
            }
        }
    }
}

#[test]
fn g1_near_endpoint_interpolation_retains_a_small_coefficient_at_large_scales() {
    for scale in [2.0_f64.powi(-800), 1.0, 2.0_f64.powi(800)] {
        let small = scale * (f64::EPSILON / 4.0);
        let table = CrossSectionTable::new(
            data(
                vec![1.0, 2.0],
                vec![vec![scale, small, scale], vec![small, scale, scale]],
            ),
            limits(),
        )
        .unwrap();
        let expected = scale * (1.25 * f64::EPSILON);
        // Using (b-a)*w+a only from the lower endpoint loses the small b
        // here, producing roughly epsilon*scale instead of 1.25 times that.
        assert_eq!(
            table
                .cross_section_at_group(2.0_f64.next_down(), 0)
                .unwrap(),
            expected
        );
        assert_eq!(
            table.cross_section_at_group(1.0_f64.next_up(), 1).unwrap(),
            expected
        );
        assert_eq!(table.cross_section_at_group(1.5, 2).unwrap(), scale);
    }
    let huge = CrossSectionTable::new(
        data(vec![300.0, 900.0], vec![vec![f64::MAX], vec![f64::MAX]]),
        limits(),
    )
    .unwrap();
    assert_eq!(huge.cross_section_at_group(600.0, 0).unwrap(), f64::MAX);
}

#[test]
fn g0_single_temperature_rows_have_an_exact_non_extrapolating_domain() {
    let table = CrossSectionTable::new(data(vec![450.0], vec![vec![0.0, 2.5]]), limits()).unwrap();
    assert_eq!(table.cross_section_at_group(450.0, 0).unwrap(), 0.0);
    assert_eq!(
        table
            .rate_constant(&flux(&table, vec![7.0, 4.0]), 450.0, 2, || Ok(()))
            .unwrap(),
        10.0
    );
    for temperature in [
        450.0_f64.next_down(),
        450.0_f64.next_up(),
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
    ] {
        assert!(matches!(
            table.cross_section_at_group(temperature, 0),
            Err(FissionError::InvalidInput(_))
        ));
    }
    assert!(matches!(
        table.cross_section_at_group(450.0, 2),
        Err(FissionError::InvalidInput(_))
    ));
}

#[test]
fn g0_zero_flux_keeps_grid_identity_temperature_and_work_admission() {
    let table = CrossSectionTable::new(fixture(), limits()).unwrap();
    let zero = flux(&table, vec![0.0; 3]);
    assert_eq!(
        table.rate_constant(&zero, 450.0, 3, || Ok(())).unwrap(),
        0.0
    );
    let equivalent_grid = Arc::new(
        NeutronGroups::new(
            table.groups().data_id().into(),
            table.groups().edges_j().to_vec(),
        )
        .unwrap(),
    );
    let equivalent = GroupFlux::new(equivalent_grid, vec![0.0; 3]).unwrap();
    assert_eq!(
        table
            .rate_constant(&equivalent, 450.0, 3, || Ok(()))
            .unwrap(),
        0.0
    );
    for grid in [
        NeutronGroups::new("different-grid-identity".into(), vec![0.0, 1.0, 2.0, 3.0]).unwrap(),
        NeutronGroups::new(table.groups().data_id().into(), vec![0.0, 1.0, 2.0, 4.0]).unwrap(),
    ] {
        let zero = GroupFlux::new(Arc::new(grid), vec![0.0; 3]).unwrap();
        assert!(matches!(
            table.rate_constant(&zero, 450.0, 3, || Ok(())),
            Err(FissionError::InvalidInput(_))
        ));
    }
    assert!(matches!(
        table.rate_constant(&zero, 300.0_f64.next_down(), 3, || Ok(())),
        Err(FissionError::InvalidInput(_))
    ));
    assert!(matches!(
        table.rate_constant(&zero, 450.0, 2, || Ok(())),
        Err(FissionError::Budget(_))
    ));
}

#[test]
fn g0_invalid_rows_and_explicit_payload_budgets_refuse() {
    for budget in [
        CrossSectionLimits {
            max_temperatures: 2,
            ..limits()
        },
        CrossSectionLimits {
            max_neutron_groups: 2,
            ..limits()
        },
        CrossSectionLimits {
            max_values: 8,
            ..limits()
        },
    ] {
        assert!(matches!(
            CrossSectionTable::new(fixture(), budget),
            Err(FissionError::Budget(_))
        ));
    }
    let mut invalid = Vec::new();
    let mut input = fixture();
    input.data_id.clear();
    invalid.push(input);
    let mut input = fixture();
    input.temperatures_k.clear();
    invalid.push(input);
    let mut input = fixture();
    input.temperatures_k[1] = 300.0;
    invalid.push(input);
    let mut input = fixture();
    input.temperatures_k[0] = 0.0;
    invalid.push(input);
    let mut input = fixture();
    input.temperatures_k[1] = f64::NAN;
    invalid.push(input);
    let mut input = fixture();
    input.cross_sections_m2[1].pop();
    invalid.push(input);
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut input = fixture();
        input.cross_sections_m2[0][0] = value;
        invalid.push(input);
    }
    for input in invalid {
        assert!(matches!(
            CrossSectionTable::new(input, limits()),
            Err(FissionError::InvalidInput(_))
        ));
    }
}

#[test]
fn g0_positive_numerical_range_failures_do_not_become_zero_hazards() {
    let table = CrossSectionTable::new(
        data(vec![1e-300, 1e308], vec![vec![0.0], vec![1e308]]),
        limits(),
    )
    .unwrap();
    // The temperature weight underflows even though its product with the
    // distant cross section could be resolved. Refusal is conservative.
    assert_eq!(
        table.cross_section_at_group(2e-300, 0),
        Err(FissionError::Arithmetic)
    );
    assert_eq!(
        table
            .rate_constant(&flux(&table, vec![0.0]), 2e-300, 1, || Ok(()))
            .unwrap(),
        0.0
    );
    let constant = CrossSectionTable::new(
        data(vec![1e-300, 1e308], vec![vec![1e-28], vec![1e-28]]),
        limits(),
    )
    .unwrap();
    assert_eq!(constant.cross_section_at_group(2e-300, 0).unwrap(), 1e-28);
    let tiny = CrossSectionTable::new(
        data(vec![1.0, 3.0], vec![vec![0.0], vec![f64::from_bits(1)]]),
        limits(),
    )
    .unwrap();
    assert_eq!(
        tiny.cross_section_at_group(1.0_f64.next_up(), 0),
        Err(FissionError::Arithmetic)
    );
    for (sigma, flux_value) in [(1e-300, 1e-300), (f64::MAX, 2.0)] {
        let table = CrossSectionTable::new(data(vec![300.0], vec![vec![sigma]]), limits()).unwrap();
        assert_eq!(
            table.rate_constant(&flux(&table, vec![flux_value]), 300.0, 1, || Ok(())),
            Err(FissionError::Arithmetic)
        );
    }
    let table =
        CrossSectionTable::new(data(vec![300.0], vec![vec![f64::MAX; 2]]), limits()).unwrap();
    assert_eq!(
        table.rate_constant(&flux(&table, vec![1.0; 2]), 300.0, 2, || Ok(())),
        Err(FissionError::Arithmetic)
    );
}

#[test]
fn g3_compensated_group_collapse_and_reciprocal_scaling_preserve_exposure() {
    let table = CrossSectionTable::new(
        data(vec![300.0, 900.0], vec![vec![1e16, 1.0, 1.0]; 2]),
        limits(),
    )
    .unwrap();
    let spectrum = flux(&table, vec![1.0; 3]);
    assert_eq!(
        table.rate_constant(&spectrum, 500.0, 3, || Ok(())).unwrap(),
        1e16 + 2.0
    );
    let input = fixture();
    let original = CrossSectionTable::new(input.clone(), limits()).unwrap();
    let factor = 2.0_f64.powi(100);
    let mut scaled_input = input;
    for row in &mut scaled_input.cross_sections_m2 {
        for value in row {
            *value *= factor;
        }
    }
    let scaled = CrossSectionTable::new(scaled_input, limits()).unwrap();
    let expected = original
        .rate_constant(&flux(&original, vec![5.0, 7.0, 11.0]), 437.0, 3, || Ok(()))
        .unwrap();
    let actual = scaled
        .rate_constant(
            &flux(&scaled, vec![5.0 / factor, 7.0 / factor, 11.0 / factor]),
            437.0,
            3,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(actual.to_bits(), expected.to_bits());
}

#[test]
fn g4_all_observed_checkpoints_refuse_and_immutable_tables_replay() {
    let mut construction_polls = 0;
    let table = CrossSectionTable::new_with_checkpoint(fixture(), limits(), || {
        construction_polls += 1;
        Ok(())
    })
    .unwrap();
    for stop in 1..=construction_polls {
        let mut polls = 0;
        assert_eq!(
            CrossSectionTable::new_with_checkpoint(fixture(), limits(), || {
                polls += 1;
                if polls == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(FissionError::Cancelled)
        );
    }
    let spectrum = flux(&table, vec![5.0, 7.0, 11.0]);
    let mut evaluation_polls = 0;
    let expected = table
        .rate_constant(&spectrum, 437.0, 3, || {
            evaluation_polls += 1;
            Ok(())
        })
        .unwrap();
    for stop in 1..=evaluation_polls {
        let mut polls = 0;
        assert_eq!(
            table.rate_constant(&spectrum, 437.0, 3, || {
                polls += 1;
                if polls == stop {
                    Err(FissionError::Budget("caller checkpoint"))
                } else {
                    Ok(())
                }
            }),
            Err(FissionError::Budget("caller checkpoint"))
        );
    }
    let cloned = table.clone();
    assert!(Arc::ptr_eq(&cloned.values_m2, &table.values_m2));
    assert_eq!(
        cloned
            .rate_constant(&spectrum, 437.0, 3, || Ok(()))
            .unwrap()
            .to_bits(),
        expected.to_bits()
    );
}
