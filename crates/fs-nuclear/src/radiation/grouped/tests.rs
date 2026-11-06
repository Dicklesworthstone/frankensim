use std::collections::BTreeMap;

use super::*;

fn limits() -> GroupedRadiationResponseLimits {
    GroupedRadiationResponseLimits {
        max_sources: 4,
        max_groups: 8,
        max_receivers: 8,
        max_source_groups: 16,
        max_transfers: 32,
    }
}

fn group(deposits: &[(u64, f64)], escape_fraction: f64) -> RadiationGroupResponse {
    RadiationGroupResponse {
        deposits: deposits
            .iter()
            .map(|&(receiver_id, fraction)| RadiationDeposit {
                receiver_id,
                fraction,
            })
            .collect(),
        escape_fraction,
    }
}

fn row(source_id: u64, groups: Vec<RadiationGroupResponse>) -> GroupedRadiationSourceResponse {
    GroupedRadiationSourceResponse { source_id, groups }
}

fn response(rows: Vec<GroupedRadiationSourceResponse>) -> GroupedRadiationResponse {
    GroupedRadiationResponse::new(
        "manufactured two-group transport response".into(),
        DecayRadiation::Gamma,
        vec![0.0, 1e-13, 4e-13],
        rows,
        limits(),
    )
    .unwrap()
}

#[test]
fn g1_changing_energy_mix_changes_receiver_heating_at_fixed_emission() {
    let response = response(vec![row(
        1,
        vec![group(&[(10, 0.75)], 0.25), group(&[(20, 0.25)], 0.75)],
    )]);
    let early = response.apply(&[(1, &[6.0, 2.0])], 0.0, || Ok(())).unwrap();
    let late = response.apply(&[(1, &[2.0, 6.0])], 0.0, || Ok(())).unwrap();
    assert_eq!(
        early.receiver_energy_j,
        BTreeMap::from([(10, 4.5), (20, 0.5)])
    );
    assert_eq!(
        late.receiver_energy_j,
        BTreeMap::from([(10, 1.5), (20, 1.5)])
    );
    assert_eq!((early.emitted_j, late.emitted_j), (8.0, 8.0));
    assert_eq!((early.deposited_j, early.escaped_j), (5.0, 3.0));
    assert_eq!((late.deposited_j, late.escaped_j), (3.0, 5.0));

    // A fixed linear transport response applies to an integrated history:
    // summing the spectra before applying it equals summing accepted intervals.
    let total = response.apply(&[(1, &[8.0, 8.0])], 0.0, || Ok(())).unwrap();
    for id in [10, 20] {
        assert_eq!(
            total.receiver_energy_j[&id],
            early.receiver_energy_j[&id] + late.receiver_energy_j[&id]
        );
    }
    assert_eq!(total.escaped_j, early.escaped_j + late.escaped_j);
    assert_eq!((total.residual_j, total.roundoff_j), (0.0, 0.0));
    assert_eq!(response.radiation(), DecayRadiation::Gamma);
    assert_eq!(response.group_edges_j(), &[0.0, 1e-13, 4e-13]);
    assert_eq!(
        response.data_id(),
        "manufactured two-group transport response"
    );
}

#[test]
fn g1_sparse_sources_self_return_and_receiver_overlap_preserve_physical_identities() {
    let response = response(vec![
        row(
            u64::MAX,
            vec![
                group(&[(u64::MAX, 0.5), (99, 0.25)], 0.25),
                group(&[(99, 0.75)], 0.25),
            ],
        ),
        row(
            2,
            vec![
                group(&[(99, 0.125)], 0.875),
                group(&[(2, 0.5), (u64::MAX, 0.25)], 0.25),
            ],
        ),
    ]);
    let result = response
        .apply(
            &[(u64::MAX, &[8.0, 4.0]), (2, &[16.0, 8.0])],
            0.0,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(
        result.receiver_energy_j,
        BTreeMap::from([(2, 4.0), (99, 7.0), (u64::MAX, 6.0)])
    );
    assert_eq!(
        (result.emitted_j, result.deposited_j, result.escaped_j),
        (36.0, 17.0, 19.0)
    );
    assert_eq!(response.source_ids().collect::<Vec<_>>(), [2, u64::MAX]);
    assert_eq!(
        response.receiver_ids().collect::<Vec<_>>(),
        [2, 99, u64::MAX]
    );
    assert_eq!((result.residual_j, result.roundoff_j), (0.0, 0.0));
}

#[test]
fn g5_source_and_receiver_order_do_not_change_grouped_results() {
    let rows = vec![
        row(
            20,
            vec![group(&[(9, 0.3), (7, 0.2)], 0.5), group(&[(7, 0.1)], 0.9)],
        ),
        row(
            10,
            vec![group(&[(7, 0.1)], 0.9), group(&[(9, 0.3), (7, 0.2)], 0.5)],
        ),
    ];
    let a = response(rows.clone());
    let mut reordered = rows;
    reordered.reverse();
    for row in &mut reordered {
        for group in &mut row.groups {
            group.deposits.reverse();
        }
    }
    let b = response(reordered);
    assert_eq!(
        a.apply(&[(20, &[3.1, 2.7]), (10, &[7.3, 0.2])], 1e-12, || Ok(()))
            .unwrap(),
        b.apply(&[(10, &[7.3, 0.2]), (20, &[3.1, 2.7])], 1e-12, || Ok(()))
            .unwrap()
    );
}

#[test]
fn g0_empty_all_escape_and_zero_emitters_are_supported_without_omitting_groups() {
    let empty = response(vec![]).apply(&[], 0.0, || Ok(())).unwrap();
    assert!(empty.receiver_energy_j.is_empty());
    assert_eq!(
        (empty.emitted_j, empty.deposited_j, empty.escaped_j),
        (0.0, 0.0, 0.0)
    );
    let response = response(vec![
        row(1, vec![group(&[], 1.0), group(&[], 1.0)]),
        row(2, vec![group(&[(9, 1.0)], 0.0), group(&[(9, 1.0)], 0.0)]),
    ]);
    let result = response
        .apply(&[(1, &[3.0, 5.0]), (2, &[0.0, 0.0])], 0.0, || Ok(()))
        .unwrap();
    assert_eq!(result.receiver_energy_j[&9], 0.0);
    assert_eq!(
        (result.emitted_j, result.deposited_j, result.escaped_j),
        (8.0, 0.0, 8.0)
    );
    assert!(response.apply(&[(1, &[3.0, 5.0])], 0.0, || Ok(())).is_err());
}

#[test]
fn g0_grid_source_group_coverage_and_individual_row_closure_are_admitted() {
    for edges in [
        vec![],
        vec![0.0],
        vec![0.0, 0.0],
        vec![2.0, 1.0],
        vec![-1.0, 1.0],
        vec![0.0, f64::INFINITY],
        vec![0.0, f64::NAN],
    ] {
        assert!(
            GroupedRadiationResponse::new(
                "synthetic".into(),
                DecayRadiation::Gamma,
                edges,
                vec![],
                limits()
            )
            .is_err()
        );
    }
    let valid = vec![group(&[(9, 0.5)], 0.5), group(&[(9, 0.5)], 0.5)];
    for rows in [
        vec![row(1, vec![])],
        vec![row(1, vec![group(&[], 1.0)])],
        vec![row(1, vec![group(&[], 1.0); 3])],
        vec![row(1, valid.clone()), row(1, valid.clone())],
        vec![row(
            1,
            vec![group(&[(9, 0.4)], 0.5), group(&[(9, 0.6)], 0.5)],
        )],
        vec![row(1, vec![group(&[(9, -0.5)], 1.0), group(&[], 1.0)])],
        vec![row(1, vec![group(&[(9, f64::NAN)], 1.0), group(&[], 1.0)])],
        vec![row(
            1,
            vec![group(&[(9, 0.5), (9, 0.5)], 0.0), group(&[], 1.0)],
        )],
    ] {
        assert!(
            GroupedRadiationResponse::new(
                "synthetic".into(),
                DecayRadiation::Gamma,
                vec![0.0, 1.0, 2.0],
                rows,
                limits()
            )
            .is_err()
        );
    }
    assert!(
        GroupedRadiationResponse::new(
            " ".into(),
            DecayRadiation::Gamma,
            vec![0.0, 1.0, 2.0],
            vec![row(1, valid)],
            limits()
        )
        .is_err()
    );
}

#[test]
fn g0_every_expanded_storage_and_sparse_work_limit_is_enforced() {
    let rows = vec![row(
        1,
        vec![group(&[(9, 1.0)], 0.0), group(&[(7, 1.0)], 0.0)],
    )];
    for constrained in [
        GroupedRadiationResponseLimits {
            max_sources: 0,
            ..limits()
        },
        GroupedRadiationResponseLimits {
            max_groups: 1,
            ..limits()
        },
        GroupedRadiationResponseLimits {
            max_receivers: 1,
            ..limits()
        },
        GroupedRadiationResponseLimits {
            max_source_groups: 1,
            ..limits()
        },
        GroupedRadiationResponseLimits {
            max_transfers: 1,
            ..limits()
        },
    ] {
        assert!(matches!(
            GroupedRadiationResponse::new(
                "synthetic".into(),
                DecayRadiation::Gamma,
                vec![0.0, 1.0, 2.0],
                rows.clone(),
                constrained
            ),
            Err(RadiationError::Budget(_))
        ));
    }
}

#[test]
fn g0_application_requires_finite_exact_source_and_group_coverage() {
    let response = response(vec![
        row(1, vec![group(&[], 1.0), group(&[], 1.0)]),
        row(2, vec![group(&[], 1.0), group(&[], 1.0)]),
    ]);
    let good = [1.0, 2.0];
    let short = [1.0];
    let long = [1.0, 2.0, 3.0];
    for inputs in [
        vec![(1, good.as_slice())],
        vec![(1, good.as_slice()), (1, good.as_slice())],
        vec![(1, good.as_slice()), (3, good.as_slice())],
        vec![(1, short.as_slice()), (2, good.as_slice())],
        vec![(1, long.as_slice()), (2, good.as_slice())],
    ] {
        assert!(matches!(
            response.apply(&inputs, 1e-9, || Ok(())),
            Err(RadiationError::InvalidInput(_))
        ));
    }
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            response.apply(&[(1, &[invalid, 1.0]), (2, &good)], 1e-9, || Ok(())),
            Err(RadiationError::InvalidInput(_))
        ));
        assert!(matches!(
            response.apply(&[(1, &good), (2, &good)], invalid, || Ok(())),
            Err(RadiationError::InvalidInput(_))
        ));
    }
}

#[test]
fn g0_original_roundoff_and_positive_underflow_refusals_survive_group_mapping() {
    let mixed = response(vec![row(
        1,
        vec![group(&[(10, 1.0)], 0.0), group(&[(20, 1.0)], 0.0)],
    )]);
    assert!(matches!(
        mixed.apply(&[(1, &[1e300, 1.0])], 0.5, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    let accepted = mixed.apply(&[(1, &[1e300, 1.0])], 10.0, || Ok(())).unwrap();
    assert_eq!(accepted.receiver_energy_j[&20], 1.0);
    assert!(accepted.roundoff_j >= 2.0);
    let fractional = response(vec![row(
        1,
        vec![group(&[(10, 0.5)], 0.5), group(&[], 1.0)],
    )]);
    assert_eq!(
        fractional.apply(&[(1, &[f64::from_bits(1), 0.0])], 1.0, || Ok(())),
        Err(RadiationError::Arithmetic)
    );
}

#[test]
fn g4_cancellation_preserves_exact_callback_error_and_reusable_response() {
    let response = response(vec![row(
        1,
        vec![group(&[(10, 0.5)], 0.5), group(&[(20, 1.0)], 0.0)],
    )]);
    let mut visits = 0;
    let expected = response
        .apply(&[(1, &[4.0, 8.0])], 0.0, || {
            visits += 1;
            Ok(())
        })
        .unwrap();
    assert!(visits > 8);
    for stop in [1, 4, visits / 2, visits] {
        let mut seen = 0;
        let error = RadiationError::InvalidInput("caller checkpoint stopped the transaction");
        assert_eq!(
            response.apply(&[(1, &[4.0, 8.0])], 0.0, || {
                seen += 1;
                if seen == stop {
                    Err(error.clone())
                } else {
                    Ok(())
                }
            }),
            Err(error)
        );
        assert_eq!(
            response.apply(&[(1, &[4.0, 8.0])], 0.0, || Ok(())).unwrap(),
            expected
        );
    }
}
