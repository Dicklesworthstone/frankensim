use super::*;

fn response(rows: Vec<RadiationSourceResponse>) -> RadiationResponse {
    RadiationResponse::new(
        "manufactured energy response v1".into(),
        rows,
        RadiationResponseLimits {
            max_sources: 8,
            max_receivers: 8,
            max_transfers: 16,
        },
    )
    .unwrap()
}

fn row(source_id: u64, deposits: &[(u64, f64)], escape_fraction: f64) -> RadiationSourceResponse {
    RadiationSourceResponse {
        source_id,
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

#[test]
fn g1_unequal_fractions_and_overlapping_receivers_conserve_energy() {
    let response = response(vec![
        row(2, &[(20, 0.25), (10, 0.5)], 0.25),
        row(1, &[(20, 0.375), (10, 0.125)], 0.5),
    ]);
    let result = response
        .apply(&[(2, 16.0), (1, 8.0)], 0.0, || Ok(()))
        .unwrap();
    assert_eq!(
        result.receiver_energy_j,
        BTreeMap::from([(10, 9.0), (20, 7.0)])
    );
    assert_eq!(
        (result.emitted_j, result.deposited_j, result.escaped_j),
        (24.0, 16.0, 8.0)
    );
    assert_eq!((result.residual_j, result.roundoff_j), (0.0, 0.0));
    assert_eq!(response.source_ids().collect::<Vec<_>>(), [1, 2]);
    assert_eq!(response.receiver_ids().collect::<Vec<_>>(), [10, 20]);
    assert_eq!(response.data_id(), "manufactured energy response v1");
}

#[test]
fn g5_input_and_authored_row_order_do_not_change_results() {
    let a = response(vec![
        row(2, &[(20, 0.3), (10, 0.2)], 0.5),
        row(1, &[(10, 0.1), (20, 0.3)], 0.6),
    ]);
    let b = response(vec![
        row(1, &[(20, 0.3), (10, 0.1)], 0.6),
        row(2, &[(10, 0.2), (20, 0.3)], 0.5),
    ]);
    assert_eq!(
        a.apply(&[(1, 3.1), (2, 7.3)], 1e-12, || Ok(())).unwrap(),
        b.apply(&[(2, 7.3), (1, 3.1)], 1e-12, || Ok(())).unwrap()
    );
    assert_eq!(a.rows(), b.rows());
}

#[test]
fn g0_all_escape_full_deposition_and_no_sources_are_valid() {
    let escape = response(vec![row(1, &[], 1.0)])
        .apply(&[(1, 3.0)], 0.0, || Ok(()))
        .unwrap();
    assert!(escape.receiver_energy_j.is_empty());
    assert_eq!((escape.deposited_j, escape.escaped_j), (0.0, 3.0));
    let full = response(vec![row(1, &[(10, 1.0)], 0.0)])
        .apply(&[(1, 3.0)], 0.0, || Ok(()))
        .unwrap();
    assert_eq!((full.deposited_j, full.escaped_j), (3.0, 0.0));
    let empty = response(vec![]).apply(&[], 0.0, || Ok(())).unwrap();
    assert_eq!(
        (empty.emitted_j, empty.deposited_j, empty.escaped_j),
        (0.0, 0.0, 0.0)
    );
}

#[test]
fn g0_invalid_coefficients_closure_and_duplicate_identities_are_refused() {
    let limits = RadiationResponseLimits {
        max_sources: 8,
        max_receivers: 8,
        max_transfers: 16,
    };
    for rows in [
        vec![row(1, &[(10, -0.1)], 1.0)],
        vec![row(1, &[(10, f64::NAN)], 0.0)],
        vec![row(1, &[(10, f64::INFINITY)], 0.0)],
        vec![row(1, &[], -0.1)],
        vec![row(1, &[], f64::NAN)],
        vec![row(1, &[], 1.1)],
        vec![row(1, &[(10, 0.4)], 0.5)],
        vec![row(1, &[(10, 0.5), (10, 0.5)], 0.0)],
        vec![row(1, &[], 1.0), row(1, &[], 1.0)],
    ] {
        assert!(RadiationResponse::new("synthetic".into(), rows, limits).is_err());
    }
    assert!(RadiationResponse::new(" ".into(), vec![], limits).is_err());
}

#[test]
fn g0_each_sparse_work_limit_is_enforced() {
    let rows = vec![row(1, &[(10, 0.5), (20, 0.25)], 0.25)];
    for limits in [
        RadiationResponseLimits {
            max_sources: 0,
            max_receivers: 2,
            max_transfers: 2,
        },
        RadiationResponseLimits {
            max_sources: 1,
            max_receivers: 1,
            max_transfers: 2,
        },
        RadiationResponseLimits {
            max_sources: 1,
            max_receivers: 2,
            max_transfers: 1,
        },
    ] {
        assert!(matches!(
            RadiationResponse::new("synthetic".into(), rows.clone(), limits),
            Err(RadiationError::Budget(_))
        ));
    }
}

#[test]
fn g0_exact_source_coverage_and_finite_nonnegative_energy_are_required() {
    let response = response(vec![row(1, &[], 1.0), row(2, &[], 1.0)]);
    for values in [
        vec![(1, 1.0)],
        vec![(1, 1.0), (2, 0.0), (3, 0.0)],
        vec![(1, 1.0), (1, 2.0)],
        vec![(1, 1.0), (3, 0.0)],
        vec![(1, -1.0), (2, 0.0)],
        vec![(1, f64::NAN), (2, 0.0)],
        vec![(1, f64::INFINITY), (2, 0.0)],
    ] {
        assert!(matches!(
            response.apply(&values, 1e-6, || Ok(())),
            Err(RadiationError::InvalidInput(_))
        ));
    }
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            response
                .apply(&[(1, 0.0), (2, 0.0)], tolerance, || Ok(()))
                .is_err()
        );
    }
}

#[test]
fn g0_small_admitted_row_defect_is_not_normalized_or_hidden() {
    let fraction = 0.5 + f64::EPSILON;
    let response = response(vec![row(1, &[(10, fraction)], 0.5)]);
    assert_eq!(response.rows()[0].deposits[0].fraction, fraction);
    assert!(matches!(
        response.apply(&[(1, 1e15)], 0.01, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    let accepted = response.apply(&[(1, 1e15)], 10.0, || Ok(())).unwrap();
    assert_eq!(accepted.receiver_energy_j[&10], 1e15 * fraction);
    assert!(accepted.roundoff_j > 0.2);
}

#[test]
fn g0_tiny_fraction_times_huge_emission_retains_real_receiver_heat() {
    let response = response(vec![row(1, &[(10, 1e-300)], 1.0)]);
    assert!(matches!(
        response.apply(&[(1, 1e300)], 0.5, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    let accepted = response.apply(&[(1, 1e300)], 10.0, || Ok(())).unwrap();
    assert_eq!(accepted.receiver_energy_j[&10], 1e300 * 1e-300);
    assert!(accepted.roundoff_j >= 1.0);
}

#[test]
fn g0_scalar_exports_charge_loss_even_when_global_residual_cancels() {
    let response = response(vec![row(1, &[(10, 1.0)], 0.0), row(2, &[(20, 1.0)], 0.0)]);
    assert!(matches!(
        response.apply(&[(1, 1e300), (2, 1.0)], 0.5, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    let result = response
        .apply(&[(1, 1e300), (2, 1.0)], 10.0, || Ok(()))
        .unwrap();
    assert_eq!(result.receiver_energy_j[&20], 1.0);
    assert_eq!(result.residual_j, 0.0);
    assert!(result.roundoff_j >= 2.0);
}

#[test]
fn g0_underflowing_positive_heat_is_refused_but_exact_subnormals_survive() {
    let response = response(vec![row(1, &[(10, 0.5)], 0.5)]);
    let tiny = f64::from_bits(1);
    assert_eq!(
        response.apply(&[(1, tiny)], 1.0, || Ok(())),
        Err(RadiationError::Arithmetic)
    );
    let exact = response.apply(&[(1, tiny * 2.0)], 0.0, || Ok(())).unwrap();
    assert_eq!(exact.receiver_energy_j[&10], tiny);
    assert_eq!(exact.escaped_j, tiny);
    let inexact = response
        .apply(&[(1, tiny * 3.0)], tiny * 8.0, || Ok(()))
        .unwrap();
    assert!(inexact.roundoff_j >= tiny * 2.0);
}

#[test]
fn g4_checkpoint_refusal_is_preserved_and_response_can_be_reused() {
    let response = response(vec![row(1, &[(10, 0.5), (20, 0.25)], 0.25)]);
    let mut calls = 0;
    let refused = response.apply(&[(1, 8.0)], 0.0, || {
        calls += 1;
        if calls == 7 {
            Err(RadiationError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert_eq!(refused, Err(RadiationError::Cancelled));
    let marker = RadiationError::InvalidInput("caller checkpoint refusal");
    assert_eq!(
        response.apply(&[(1, 8.0)], 0.0, || Err(marker.clone())),
        Err(marker)
    );
    assert_eq!(
        response
            .clone()
            .apply(&[(1, 8.0)], 0.0, || Ok(()))
            .unwrap()
            .deposited_j,
        6.0
    );
}
