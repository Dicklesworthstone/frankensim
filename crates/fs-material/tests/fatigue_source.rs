//! G0/G1/G3 original-source admission and analytic S-N convention checks.
//! Synthetic fixtures establish the implementation, not measured service life.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::source::{
    SnLifeMeasure, SnSourceRequest, SnStressMeasure, resolve_log_log_sn,
};
use fs_material::fatigue::{FatigueError, MeanStressCorrection, RainflowLimits, RainflowMiner};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn request(stress_measure: SnStressMeasure, life_measure: SnLifeMeasure) -> SnSourceRequest {
    SnSourceRequest {
        stress_measure,
        life_measure,
        temperature_axis: "temperature".into(),
        stress_ratio_axis: "stress-ratio".into(),
        surface_and_geometry: "smooth coupon".into(),
        environment: "synthetic dry air".into(),
        load_basis: "fully reversed axial normal stress".into(),
        maximum_points: 3,
    }
}
fn point(request: &SnSourceRequest) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(Temperature::DIMS),
            293.15,
        )
        .unwrap()
        .with_quantity("stress-ratio", QuantitySpec::dimensional(Dims::NONE), -1.0)
        .unwrap()
        .with_quantity(
            request.life_measure.axis_name(),
            QuantitySpec::dimensional(Dims::NONE),
            if request.life_measure == SnLifeMeasure::Cycles {
                10.0
            } else {
                20.0
            },
        )
        .unwrap()
}
fn claim(request: &SnSourceRequest, interpolation: InterpolationPolicy) -> PropertyClaim {
    let count_scale = if request.life_measure == SnLifeMeasure::Cycles {
        1.0
    } else {
        2.0
    };
    let stress_scale = if request.stress_measure == SnStressMeasure::Amplitude {
        1.0
    } else {
        2.0
    };
    PropertyClaim {
        key: request.stress_measure.property_key(),
        value: PropertyValue::Curve {
            abscissa: request.life_measure.axis_name().into(),
            abscissa_dims: Dims::NONE,
            knots: [(10.0, 80e6), (100.0, 40e6), (1000.0, 20e6)]
                .map(|(n, s)| (n * count_scale, s * stress_scale))
                .into(),
            dims: Pressure::DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                QuantitySpec::dimensional(Temperature::DIMS),
                293.15,
                293.15,
            )
            .with_quantity(
                "stress-ratio",
                QuantitySpec::dimensional(Dims::NONE),
                -1.0,
                -1.0,
            )
            .with_quantity(
                request.life_measure.axis_name(),
                QuantitySpec::dimensional(Dims::NONE),
                10.0 * count_scale,
                1000.0 * count_scale,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic S-N convention fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for row in rows {
        claims.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "synthetic machining".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn g1_range_and_reversal_sources_give_identical_amplitudes_lives_and_damage() {
    let amplitude = request(SnStressMeasure::Amplitude, SnLifeMeasure::Cycles);
    let range = request(SnStressMeasure::Range, SnLifeMeasure::Reversals);
    let a = card(vec![claim(&amplitude, InterpolationPolicy::TabulatedOnly)]);
    let b = card(vec![claim(&range, InterpolationPolicy::TabulatedOnly)]);
    let curve_a = resolve_log_log_sn(&a, &point(&amplitude), &amplitude, SINGLE).unwrap();
    let curve_b = resolve_log_log_sn(&b, &point(&range), &range, SINGLE).unwrap();
    assert_eq!(curve_a.points(), curve_b.points());
    assert_eq!(curve_a.cycles_to_failure(40e6).unwrap(), Some(100.0));
    let interpolated = curve_a
        .cycles_to_failure((20e6_f64 * 40e6).sqrt())
        .unwrap()
        .unwrap();
    assert!((interpolated / (1000.0_f64 * 100.0).sqrt() - 1.0).abs() < 1e-12);
    let limits = RainflowLimits {
        maximum_samples: 20,
        maximum_residue: 10,
    };
    let mut ca =
        RainflowMiner::new(Arc::new(curve_a), MeanStressCorrection::Uncorrected, limits).unwrap();
    let mut cb =
        RainflowMiner::new(Arc::new(curve_b), MeanStressCorrection::Uncorrected, limits).unwrap();
    for s in [40e6, -40e6, 40e6, -40e6, 40e6] {
        ca.push(s).unwrap();
        cb.push(s).unwrap();
    }
    let result = ca.finish().unwrap();
    assert_eq!(result, cb.finish().unwrap());
    assert!((result.miner_damage - 2.0 / 100.0).abs() < 1e-15);
    assert_eq!(
        ca.curve().source().unwrap().card_identity(),
        a.content_hash()
    );
    assert_eq!(
        cb.curve().source().unwrap().card_identity(),
        b.content_hash()
    );
}

#[test]
fn g0_original_tabulated_policy_uncertainty_and_all_receipts_survive() {
    let req = request(SnStressMeasure::Amplitude, SnLifeMeasure::Cycles);
    let original = claim(&req, InterpolationPolicy::TabulatedOnly);
    let material = card(vec![original.clone()]);
    let at = point(&req);
    let curve = resolve_log_log_sn(&material, &at, &req, SINGLE).unwrap();
    let retained = curve.source().unwrap();
    assert_eq!(retained.claim(), &original);
    assert_eq!(retained.query(), &at);
    assert_eq!(retained.selection(), SINGLE);
    assert_eq!(retained.material_state(), material.id());
    assert_eq!(retained.observations().len(), 3);
    for (life, answer) in retained.observations() {
        let expected = material
            .claims()
            .query_typed(
                &req.stress_measure.property_key(),
                &at.clone()
                    .with(req.life_measure.axis_name(), *life)
                    .unwrap(),
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(*answer, expected);
        material.claims().verify_receipt(&answer.receipt).unwrap();
    }
    assert!(
        material
            .claims()
            .query_typed(
                &req.stress_measure.property_key(),
                &at.with(req.life_measure.axis_name(), 50.0).unwrap(),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    assert!(curve.endurance_limit_pa().is_none());
    assert!(matches!(
        curve.cycles_to_failure(10e6),
        Err(FatigueError::AmplitudeOutsideCurve { .. })
    ));
    assert!(matches!(
        curve.cycles_to_failure(90e6),
        Err(FatigueError::AmplitudeOutsideCurve { .. })
    ));
}

#[test]
fn g0_missing_or_incompatible_mean_stress_temperature_and_count_conventions_refuse() {
    let req = request(SnStressMeasure::Amplitude, SnLifeMeasure::Cycles);
    let material = card(vec![claim(&req, InterpolationPolicy::TabulatedOnly)]);
    for at in [
        point(&req).with("stress-ratio", 0.0).unwrap(),
        point(&req).with("temperature", 300.0).unwrap(),
        point(&req)
            .with_quantity("temperature", QuantitySpec::dimensional(Dims::NONE), 293.15)
            .unwrap(),
    ] {
        assert!(resolve_log_log_sn(&material, &at, &req, SINGLE).is_err());
    }
    let mut wrong = req.clone();
    wrong.life_measure = SnLifeMeasure::Reversals;
    assert!(resolve_log_log_sn(&material, &point(&wrong), &wrong, SINGLE).is_err());
    wrong = req.clone();
    wrong.stress_measure = SnStressMeasure::Range;
    assert!(resolve_log_log_sn(&material, &point(&wrong), &wrong, SINGLE).is_err());
    let mut non_reversed = claim(&req, InterpolationPolicy::TabulatedOnly);
    non_reversed.validity = non_reversed.validity.with_quantity(
        "stress-ratio",
        QuantitySpec::dimensional(Dims::NONE),
        -1.0,
        0.0,
    );
    assert!(resolve_log_log_sn(&card(vec![non_reversed]), &point(&req), &req, SINGLE).is_err());
}

#[test]
fn g0_no_source_knot_is_discarded_to_repair_a_curve_or_escape_its_validity() {
    let req = request(SnStressMeasure::Amplitude, SnLifeMeasure::Cycles);
    let mut bounded = claim(&req, InterpolationPolicy::TabulatedOnly);
    bounded.validity = bounded.validity.with_quantity(
        req.life_measure.axis_name(),
        QuantitySpec::dimensional(Dims::NONE),
        10.0,
        100.0,
    );
    assert!(resolve_log_log_sn(&card(vec![bounded]), &point(&req), &req, SINGLE).is_err());
    let mut nonmonotone = claim(&req, InterpolationPolicy::TabulatedOnly);
    if let PropertyValue::Curve { knots, .. } = &mut nonmonotone.value {
        knots[1].1 = 10e6;
    }
    assert!(resolve_log_log_sn(&card(vec![nonmonotone]), &point(&req), &req, SINGLE).is_err());
    let material = card(vec![claim(&req, InterpolationPolicy::LinearInside)]);
    let mut budget = req.clone();
    budget.maximum_points = 2;
    assert!(resolve_log_log_sn(&material, &point(&req), &budget, SINGLE).is_err());
    assert!(
        resolve_log_log_sn(
            &material,
            &point(&req)
                .with(req.life_measure.axis_name(), 20.0)
                .unwrap(),
            &req,
            SINGLE
        )
        .is_err()
    );
}

#[test]
fn g0_conflicting_sources_require_the_same_explicit_pin_at_every_knot() {
    let req = request(SnStressMeasure::Amplitude, SnLifeMeasure::Cycles);
    let first = claim(&req, InterpolationPolicy::TabulatedOnly);
    let mut conflict = first.clone();
    conflict.provenance.source = "different synthetic fatigue test".into();
    conflict.validity = conflict.validity.with_quantity(
        req.life_measure.axis_name(),
        QuantitySpec::dimensional(Dims::NONE),
        100.0,
        1000.0,
    );
    let material = card(vec![first, conflict]);
    let anchor = material
        .claims()
        .query_typed(
            &req.stress_measure.property_key(),
            &point(&req),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!(resolve_log_log_sn(&material, &point(&req), &req, SINGLE).is_err());
    let pin = ClaimSelection::Pinned(anchor.receipt.selected);
    let resolved = resolve_log_log_sn(&material, &point(&req), &req, pin).unwrap();
    assert!(
        resolved
            .source()
            .unwrap()
            .observations()
            .iter()
            .all(|(_, a)| a.receipt.selected == anchor.receipt.selected)
    );
}

#[test]
fn g0_static_strength_is_not_a_fatigue_curve() {
    let req = request(SnStressMeasure::Amplitude, SnLifeMeasure::Cycles);
    let mut row = claim(&req, InterpolationPolicy::ConstantWithinValidity);
    row.value = PropertyValue::Scalar {
        value: 80e6,
        dims: Pressure::DIMS,
    };
    assert!(resolve_log_log_sn(&card(vec![row]), &point(&req), &req, SINGLE).is_err());
}
