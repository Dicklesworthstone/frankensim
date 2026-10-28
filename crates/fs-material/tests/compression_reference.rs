//! G0/G1/G3 source-bound compression references; synthetic numerical fixtures.
//! These prove the bounded reference calculation, not specimen calibration.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, EvaluationDecision, InterpolationPolicy, MaterialCard,
    MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::compression_reference::{
    CompressionReferenceInterpolation, CompressionReferenceRequest, ResolvedCompressionReference,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn context() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        "source-grade",
        QuantitySpec::dimensional(Dims::NONE),
        1.0,
        1.0,
    )
}
fn density_claim() -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("density", Density::DIMS),
        value: PropertyValue::Scalar {
            value: 1000.0,
            dims: Density::DIMS,
        },
        validity: context(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic compression reference density".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn stress_claim(interpolation: InterpolationPolicy) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("compressive-stress", Pressure::DIMS),
        value: PropertyValue::Curve {
            abscissa: "strain".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.01, 10e6), (0.02, 30e6), (0.05, 60e6)],
            dims: Pressure::DIMS,
        },
        validity: context()
            .with_quantity("strain", QuantitySpec::dimensional(Dims::NONE), 0.01, 0.05)
            .with_quantity("test-speed", QuantitySpec::dimensional(SPEED), 0.001, 0.001),
        uncertainty: UncertaintyModel::Unstated,
        interpolation,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic compression reference stress".into(),
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
            chemistry: "synthetic compression reference".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("source-grade", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
        .with_quantity("strain", QuantitySpec::dimensional(Dims::NONE), 0.01)
        .unwrap()
        .with_quantity("test-speed", QuantitySpec::dimensional(SPEED), 0.001)
        .unwrap()
}
fn request(interpolation: CompressionReferenceInterpolation) -> CompressionReferenceRequest {
    CompressionReferenceRequest {
        density: PropertyKey::new("density", Density::DIMS),
        stress: PropertyKey::new("compressive-stress", Pressure::DIMS),
        strain_axis: "strain".into(),
        crosshead_speed_axis: Some("test-speed".into()),
        interpolation,
    }
}
fn resolve(card: &MaterialCard) -> ResolvedCompressionReference {
    ResolvedCompressionReference::resolve(
        card,
        &point(),
        &request(CompressionReferenceInterpolation::SourceInterpolationOnly),
        SINGLE,
        SINGLE,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_nonuniform_slopes_give_independent_force_tangent_and_segment_work() {
    let source = card(vec![
        density_claim(),
        stress_claim(InterpolationPolicy::LinearInside),
    ]);
    let law = resolve(&source);
    let mut coupon = law.coupon(0.1, 0.002).unwrap();
    close(coupon.mass_kg(), 0.2);
    close(coupon.initial_shortening_m(), 0.001);
    assert!(coupon.last_response().is_none());
    let first = coupon.advance_shortening(0.0015, 0.5).unwrap();
    close(first.nominal_stress_pa, 20e6);
    close(first.force_n, 40_000.0);
    close(first.tangent_n_m, 40e6);
    close(first.model_loading_work_j, 15.0);
    let crossed = coupon.advance_shortening(0.004, 2.5).unwrap();
    close(crossed.nominal_stress_pa, 50e6);
    close(crossed.force_n, 100_000.0);
    close(crossed.tangent_n_m, 20e6);
    // Two distinct trapezoids: 0.005*(20+30)/2 + 0.02*(30+50)/2 MPa.
    close(crossed.incremental_model_loading_work_j, 185.0);
    close(crossed.model_loading_work_j, 200.0);
    let end = coupon.advance_shortening(0.005, 1.0).unwrap();
    close(end.force_n, 120_000.0);
    close(end.tangent_n_m, 20e6);
    close(end.incremental_model_loading_work_j, 110.0);
    close(end.model_loading_work_j, 310.0);
    let knot = law
        .coupon(0.1, 0.002)
        .unwrap()
        .advance_shortening(0.002, 1.0)
        .unwrap();
    close(knot.tangent_n_m, 20e6); // Right-hand tangent at the interior knot.
    close(knot.model_loading_work_j, 40.0);
}

#[test]
fn g0_discrete_observations_require_opt_in_and_keep_original_receipts() {
    let source = card(vec![
        density_claim(),
        stress_claim(InterpolationPolicy::TabulatedOnly),
    ]);
    let identity = source.content_hash();
    let mut req = request(CompressionReferenceInterpolation::SourceInterpolationOnly);
    assert!(
        ResolvedCompressionReference::resolve(&source, &point(), &req, SINGLE, SINGLE).is_err()
    );
    req.interpolation = CompressionReferenceInterpolation::PiecewiseLinearBetweenObservations;
    let law =
        ResolvedCompressionReference::resolve(&source, &point(), &req, SINGLE, SINGLE).unwrap();
    for (strain, retained) in law.observations() {
        let expected = source
            .claims()
            .query_typed(
                &req.stress,
                &point().with("strain", *strain).unwrap(),
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(*retained, expected);
        assert!(matches!(
            retained.receipt.decision,
            EvaluationDecision::ExactTabulated { .. }
        ));
        assert_eq!(
            retained.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        source.claims().verify_receipt(&retained.receipt).unwrap();
    }
    assert_eq!(source.content_hash(), identity);
    assert!(
        source
            .claims()
            .query_typed(
                &req.stress,
                &point().with("strain", 0.015).unwrap(),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    close(
        law.coupon(0.1, 0.002)
            .unwrap()
            .advance_shortening(0.0015, 0.5)
            .unwrap()
            .force_n,
        40_000.0,
    );
    assert_eq!(
        source
            .claims()
            .claim(law.observations()[0].1.receipt.selected)
            .unwrap()
            .interpolation,
        InterpolationPolicy::TabulatedOnly
    );
}

#[test]
fn g0_origin_extrapolation_and_trimmed_source_domains_refuse() {
    let source = card(vec![
        density_claim(),
        stress_claim(InterpolationPolicy::LinearInside),
    ]);
    let law = resolve(&source);
    assert_eq!(law.strain_domain(), [0.01, 0.05]);
    let mut coupon = law.coupon(0.1, 0.002).unwrap();
    let identity = coupon.identity();
    for (shortening, duration) in [(0.0, 1.0), (0.0005, 0.5), (0.0051, 4.1)] {
        assert!(coupon.advance_shortening(shortening, duration).is_err());
        assert_eq!(coupon.identity(), identity);
        assert!(coupon.last_response().is_none());
    }
    let mut stress = stress_claim(InterpolationPolicy::LinearInside);
    stress.validity =
        stress
            .validity
            .with_quantity("strain", QuantitySpec::dimensional(Dims::NONE), 0.02, 0.05);
    let narrowed = card(vec![density_claim(), stress]);
    // Selecting an admitted point cannot silently remove the original 1% knot.
    assert!(
        ResolvedCompressionReference::resolve(
            &narrowed,
            &point().with("strain", 0.02).unwrap(),
            &request(CompressionReferenceInterpolation::SourceInterpolationOnly),
            SINGLE,
            SINGLE
        )
        .is_err()
    );
}

#[test]
fn g0_rejected_rate_hold_unloading_and_domain_steps_are_atomic() {
    let source = card(vec![
        density_claim(),
        stress_claim(InterpolationPolicy::LinearInside),
    ]);
    let mut coupon = resolve(&source).coupon(0.1, 0.002).unwrap();
    let accepted = coupon.advance_shortening(0.002, 1.0).unwrap();
    for (shortening, duration) in [
        (0.002, 0.5),
        (0.0015, 0.5),
        (-0.001, 0.5),
        (0.003, 2.0),
        (0.003, 0.0),
        (0.003, f64::NAN),
        (0.003, f64::INFINITY),
        (0.006, 4.0),
    ] {
        assert!(coupon.advance_shortening(shortening, duration).is_err());
        assert_eq!(coupon.identity(), accepted.identity);
        assert_eq!(coupon.last_response(), Some(&accepted));
    }
    let next = coupon.advance_shortening(0.003, 1.0).unwrap();
    close(next.elapsed_s, 2.0);
    close(next.force_n, 80_000.0);
    close(next.model_loading_work_j, 110.0);
    let mut no_speed_owner = request(CompressionReferenceInterpolation::SourceInterpolationOnly);
    no_speed_owner.crosshead_speed_axis = None;
    assert!(
        ResolvedCompressionReference::resolve(&source, &point(), &no_speed_owner, SINGLE, SINGLE)
            .is_err()
    );
}

#[test]
fn g3_step_subdivision_preserves_response_and_work_but_retains_path_identity() {
    let source = card(vec![
        density_claim(),
        stress_claim(InterpolationPolicy::LinearInside),
    ]);
    let law = resolve(&source);
    let direct = law
        .coupon(0.1, 0.002)
        .unwrap()
        .advance_shortening(0.005, 4.0)
        .unwrap();
    let mut split = law.coupon(0.1, 0.002).unwrap();
    let mut increments = 0.0;
    for (shortening, duration) in [(0.0015, 0.5), (0.002, 0.5), (0.003, 1.0), (0.005, 2.0)] {
        increments += split
            .advance_shortening(shortening, duration)
            .unwrap()
            .incremental_model_loading_work_j;
    }
    let end = split.last_response().unwrap();
    close(end.elapsed_s, direct.elapsed_s);
    close(end.force_n, direct.force_n);
    close(end.tangent_n_m, direct.tangent_n_m);
    close(end.model_loading_work_j, direct.model_loading_work_j);
    close(increments, 310.0);
    assert_ne!(end.identity, direct.identity);
    let replay = law
        .coupon(0.1, 0.002)
        .unwrap()
        .advance_shortening(0.005, 4.0)
        .unwrap();
    assert_eq!(direct, replay);
}

#[test]
fn g0_scalar_strength_wrong_units_and_sampled_density_cannot_be_substituted() {
    let req = request(CompressionReferenceInterpolation::SourceInterpolationOnly);
    let mut scalar_stress = stress_claim(InterpolationPolicy::ConstantWithinValidity);
    scalar_stress.value = PropertyValue::Scalar {
        value: 30e6,
        dims: Pressure::DIMS,
    };
    let mut sampled_density = density_claim();
    sampled_density.interpolation = InterpolationPolicy::TabulatedOnly;
    for rows in [
        vec![density_claim(), scalar_stress],
        vec![
            sampled_density,
            stress_claim(InterpolationPolicy::LinearInside),
        ],
    ] {
        assert!(
            ResolvedCompressionReference::resolve(&card(rows), &point(), &req, SINGLE, SINGLE)
                .is_err()
        );
    }
    let source = card(vec![
        density_claim(),
        stress_claim(InterpolationPolicy::LinearInside),
    ]);
    let mut wrong_units = req.clone();
    wrong_units.stress = PropertyKey::new("compressive-stress", Density::DIMS);
    assert!(
        ResolvedCompressionReference::resolve(&source, &point(), &wrong_units, SINGLE, SINGLE)
            .is_err()
    );
    let mut strength = stress_claim(InterpolationPolicy::LinearInside);
    strength.key = PropertyKey::new("compressive-strength", Pressure::DIMS);
    let mut wrong_role = req.clone();
    wrong_role.stress = strength.key.clone();
    assert!(
        ResolvedCompressionReference::resolve(
            &card(vec![density_claim(), strength]),
            &point(),
            &wrong_role,
            SINGLE,
            SINGLE
        )
        .is_err()
    );
    assert!(
        ResolvedCompressionReference::resolve(
            &source,
            &point().with("source-grade", 2.0).unwrap(),
            &req,
            SINGLE,
            SINGLE
        )
        .is_err()
    );
}

#[test]
fn g0_explicit_claim_pins_preserve_one_source_over_all_knots() {
    let first = stress_claim(InterpolationPolicy::LinearInside);
    let mut second = first.clone();
    second.provenance.source = "independent synthetic stress observation".into();
    if let PropertyValue::Curve { knots, .. } = &mut second.value {
        for (_, value) in knots {
            *value *= 2.0;
        }
    }
    let source = card(vec![density_claim(), first, second]);
    let req = request(CompressionReferenceInterpolation::SourceInterpolationOnly);
    assert!(
        ResolvedCompressionReference::resolve(&source, &point(), &req, SINGLE, SINGLE).is_err()
    );
    let mut work = Vec::new();
    for (id, _) in source.claims().claims_for("compressive-stress") {
        let law = ResolvedCompressionReference::resolve(
            &source,
            &point(),
            &req,
            SINGLE,
            ClaimSelection::Pinned(id),
        )
        .unwrap();
        assert!(
            law.observations()
                .iter()
                .all(|(_, answer)| answer.receipt.selected == id)
        );
        work.push(
            law.coupon(0.1, 0.002)
                .unwrap()
                .advance_shortening(0.005, 4.0)
                .unwrap()
                .model_loading_work_j,
        );
    }
    work.sort_by(f64::total_cmp);
    close(work[0], 310.0);
    close(work[1], 620.0);
    let density_pin = source.claims().claims_for("density")[0].0;
    assert!(
        ResolvedCompressionReference::resolve(
            &source,
            &point(),
            &req,
            SINGLE,
            ClaimSelection::Pinned(density_pin)
        )
        .is_err()
    );
}
