//! Synthetic affine-profile integration references; no material calibration claim.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::depth::{EXPOSED_MATERIAL_DEPTH_AXIS as Z, ResolvedDepthHardness};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec};

fn profile(knots: Vec<(f64, f64)>) -> ResolvedDepthHardness {
    let maximum = knots.last().unwrap().0;
    let length = Dims([1, 0, 0, 0, 0, 0]);
    let provenance = Provenance {
        source: "synthetic hardness-depth integral fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic graded specimen".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"depth integral"),
            caveats: "numerical references only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic depth profile revision 1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: Z.into(),
                abscissa_dims: length,
                knots,
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                Z,
                QuantitySpec::dimensional(length),
                0.0,
                maximum,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "synthetic depth profile".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedDepthHardness::new(
        &card,
        &QueryPoint::new(),
        key,
        MaterialPropertySelection::SingleClaimOnly,
        maximum,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn constant_hardening_and_softening_profiles_match_independent_affine_integrals() {
    for (h0, h1) in [(1.0e9, 1.0e9), (1.0e8, 1.0e9), (1.0e9, 1.0e8)] {
        let source = profile(vec![(0.0, h0), (1.0, h1)]);
        for (start, end) in [(0.0, 0.25), (0.125, 0.5), (0.25, 1.0)] {
            let expected = h0 * (end - start) + 0.5 * (h1 - h0) * (end * end - start * start);
            close(
                source.hardness_integral_between(start, end).unwrap(),
                expected,
                3.0e-15,
            );
            let next = source.advance_by_integral(start, expected).unwrap();
            close(next.depth_m, end, 3.0e-15);
            close(next.increment_m, end - start, 5.0e-15);
            close(
                next.path_mean_hardness_pa,
                expected / (end - start),
                5.0e-15,
            );
        }
    }
}

#[test]
fn crossed_knots_and_path_partition_retain_the_same_original_profile() {
    let source = profile(vec![
        (0.0, 100.0),
        (1.0, 100.0),
        (2.0, 300.0),
        (3.0, 50.0),
        (4.0, 50.0),
    ]);
    let refined = profile(vec![
        (0.0, 100.0),
        (0.5, 100.0),
        (1.0, 100.0),
        (1.5, 200.0),
        (2.0, 300.0),
        (2.5, 175.0),
        (3.0, 50.0),
        (3.5, 50.0),
        (4.0, 50.0),
    ]);
    assert_eq!(source.hardness_integral_between(0.0, 4.0).unwrap(), 525.0);
    let advance = source.advance_by_integral(0.5, 450.0).unwrap();
    assert_eq!(advance.depth_m, 3.5);
    assert_eq!(advance.path_mean_hardness_pa, 150.0);
    for endpoint in [1.0, 2.0, 3.0, 4.0] {
        let integral = source.hardness_integral_between(0.0, endpoint).unwrap();
        assert_eq!(
            source.advance_by_integral(0.0, integral).unwrap().depth_m,
            endpoint
        );
    }
    let first = source.hardness_integral_between(0.37, 1.2).unwrap();
    let second = source.hardness_integral_between(1.2, 3.4).unwrap();
    let split = source
        .advance_by_integral(
            source.advance_by_integral(0.37, first).unwrap().depth_m,
            second,
        )
        .unwrap();
    let whole = source.advance_by_integral(0.37, first + second).unwrap();
    let refined_whole = refined.advance_by_integral(0.37, first + second).unwrap();
    close(split.depth_m, 3.4, 3.0e-15);
    close(whole.depth_m, split.depth_m, 3.0e-15);
    close(refined_whole.depth_m, whole.depth_m, 3.0e-15);
}

#[test]
fn small_late_increments_do_not_subtract_large_reference_integrals() {
    let source = profile(vec![(0.0, 1.0e12), (1.0, 2.0e12)]);
    let start = 0.875;
    let increment = 2.0_f64.powi(-40);
    let end = start + increment;
    let expected = (1.0e12 * (1.0 + start) + 0.5e12 * increment) * increment;
    close(
        source.hardness_integral_between(start, end).unwrap(),
        expected,
        2.0e-15,
    );
    let next = source.advance_by_integral(start, expected).unwrap();
    assert_eq!(next.depth_m.to_bits(), end.to_bits());
    assert_eq!(next.increment_m.to_bits(), increment.to_bits());
}

#[test]
fn scaled_inverse_avoids_squaring_dimensional_hardness() {
    for scale in [1.0e-200, 1.0, 1.0e200] {
        let source = profile(vec![(0.0, scale), (0.1, 2.0 * scale)]);
        let expected = scale * 0.0725;
        close(
            source.hardness_integral_between(0.02, 0.07).unwrap(),
            expected,
            4.0e-15,
        );
        close(
            source.advance_by_integral(0.02, expected).unwrap().depth_m,
            0.07,
            4.0e-15,
        );
    }
}

#[test]
fn source_limits_invalid_measures_and_lost_depth_resolution_refuse() {
    let source = profile(vec![(0.0, 100.0), (1.0, 200.0)]);
    let before = source.clone();
    for measure in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(source.advance_by_integral(0.0, measure).is_err());
    }
    for depth in [-1.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(source.advance_by_integral(depth, 0.0).is_err());
    }
    assert!(source.advance_by_integral(0.0, 150.000001).is_err());
    assert!(source.advance_by_integral(1.0, 1.0).is_err());
    assert!(source.advance_by_integral(0.5, f64::MIN_POSITIVE).is_err());
    assert!(source.hardness_integral_between(0.5, 0.25).is_err());
    assert!(source.hardness_integral_between(0.0, 1.01).is_err());
    assert_eq!(source.hardness_integral_between(1.0, 1.0).unwrap(), 0.0);
    let unchanged = source.advance_by_integral(1.0, 0.0).unwrap();
    assert_eq!(unchanged.depth_m, 1.0);
    assert_eq!(unchanged.increment_m, 0.0);
    assert_eq!(unchanged.path_mean_hardness_pa, 200.0);
    assert_eq!(source, before);
}
