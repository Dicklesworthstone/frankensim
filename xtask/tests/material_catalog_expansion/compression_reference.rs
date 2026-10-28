use super::*;

use fs_matdb::{ClaimSelection, MaterialStateId, NormalizedMaterialCardPack};
use fs_material::compression_reference::{
    CompressionReferenceInterpolation, CompressionReferenceRequest, ResolvedCompressionReference,
};

fn source_reference(
    pack: NormalizedPack,
    point: &QueryPoint,
    speed: Option<&str>,
) -> ResolvedCompressionReference {
    let card = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid".into(),
            process: "source-bound compression reference".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let request = CompressionReferenceRequest {
        density: key("density", RHO),
        stress: key("compressive-stress", PA),
        strain_axis: "compressive-strain".into(),
        crosshead_speed_axis: speed.map(str::to_owned),
        interpolation: CompressionReferenceInterpolation::PiecewiseLinearBetweenObservations,
    };
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let reference =
        ResolvedCompressionReference::resolve(card.card(), point, &request, selection, selection)
            .unwrap();
    for (_, answer) in reference.observations() {
        card.card()
            .claims()
            .verify_receipt(&answer.receipt)
            .unwrap();
        assert!(matches!(
            answer.receipt.decision,
            fs_matdb::EvaluationDecision::ExactTabulated { .. }
        ));
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    reference
}

#[test]
fn g1_g3_cushioning_packs_drive_distinct_bounded_compression_responses() {
    for (grade, count, density, stress25, stress50) in [
        ("plastazote-ld33", 8, 33.0, 66000.0, 134000.0),
        ("plastazote-ld45", 8, 45.0, 91000.0, 164000.0),
        ("evazote-ev30", 7, 30.0, 39000.0, 100000.0),
        ("evazote-ev50", 7, 50.0, 50000.0, 116000.0),
    ] {
        let selector = format!("source-grade-{grade}-2026");
        let point = flags(&[
            (&selector, 1.0),
            ("density-includes-skins", 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-orientation-known", 0.0),
            ("source-test-rate-known", 0.0),
            ("source-specimen-geometry-known", 0.0),
            ("compressive-strain", 0.25),
        ]);
        let pack = compile(&format!("zotefoams-{grade}-2026"), count);
        // Authoring a reference does not grant the source an interpolation policy.
        refuses(
            &pack,
            "compressive-stress",
            PA,
            &point.clone().with("compressive-strain", 0.375).unwrap(),
        );
        let reference = source_reference(pack, &point, None);
        let mut coupon = reference.coupon(0.01, 0.0001).unwrap();
        assert_eq!(coupon.initial_shortening_m(), 0.0025);
        assert!((coupon.mass_kg() - density * 1e-6).abs() < 1e-15);
        let middle = coupon.advance_shortening(0.00375, 1.0).unwrap();
        assert!((middle.force_n - 0.0001 * (stress25 + stress50) / 2.0).abs() < 1e-12);
        let end = coupon.advance_shortening(0.005, 1.0).unwrap();
        assert!((end.force_n - 0.0001 * stress50).abs() < 1e-12);
        assert!((end.tangent_n_m - (stress50 - stress25) / 0.25 * 0.0001 / 0.01).abs() < 1e-9);
        assert!(
            (end.model_loading_work_j - 1e-6 * 0.25 * (stress25 + stress50) / 2.0).abs() < 1e-15
        );
        assert!(coupon.advance_shortening(0.0051, 1.0).is_err());
        assert_eq!(coupon.last_response(), Some(&end));
    }
}

#[test]
fn g1_g3_reinforced_peek_uses_all_compression_segments_and_source_speed() {
    let point = flags(&[
        ("source-grade-tecapeek-gf30-natural", 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-orientation-known", 0.0),
        ("compressive-strain", 0.01),
    ])
    .with_quantity(
        "test-speed",
        QuantitySpec::dimensional(SPEED),
        5.0 / 60000.0,
    )
    .unwrap();
    let reference = source_reference(
        compile("ensinger-tecapeek-gf30-natural", 6),
        &point,
        Some("test-speed"),
    );
    let mut coupon = reference.coupon(0.01, 0.0001).unwrap();
    assert!((coupon.mass_kg() - 0.00153).abs() < 1e-15);
    assert!(coupon.advance_shortening(0.0002, 1.0).is_err());
    assert!(coupon.last_response().is_none());
    let at_two_percent = coupon.advance_shortening(0.0002, 1.2).unwrap();
    assert_eq!(at_two_percent.force_n, 5200.0);
    assert!((at_two_percent.model_loading_work_j - 0.405).abs() < 1e-14);
    let at_five_percent = coupon.advance_shortening(0.0005, 3.6).unwrap();
    assert_eq!(at_five_percent.force_n, 12000.0);
    assert!((at_five_percent.model_loading_work_j - 2.985).abs() < 1e-13);
    assert!((at_five_percent.incremental_model_loading_work_j - 2.58).abs() < 1e-13);
}
