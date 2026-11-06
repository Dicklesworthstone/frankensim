//! G0/G3: directional tensile limits retain material axes and source conditions.
//! All coefficients are manufactured; these tests do not validate a material.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::strength::OrthotropicTensileLimit;
use fs_material::tensor;
use fs_qty::{Dims, Pressure, Temperature};

const NAMES: [&str; 3] = ["tensile_along_1", "tensile_along_2", "tensile_along_3"];
const SELECTION: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;

fn card(limits: [f64; 3], dims: Dims, extra_axis: Option<&str>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value) in NAMES.into_iter().zip(limits) {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 500.0);
        if let Some(axis) = extra_axis {
            validity = validity.with(axis, 1.0, 2.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, value), (500.0, value * 0.5)],
                    dims,
                },
                validity,
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured directional tensile test".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "declared orthotropic specimen axes".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
}

fn limit(limits: [f64; 3]) -> OrthotropicTensileLimit {
    OrthotropicTensileLimit::resolve(
        &card(limits, Pressure::DIMS, None),
        &point(300.0),
        NAMES,
        SELECTION,
    )
    .unwrap()
}

#[test]
fn g0_each_material_axis_uses_its_own_tensile_strength() {
    let limit = limit([100.0, 10.0, 2.0]);
    for (axis, strength) in limit.limits_pa().into_iter().enumerate() {
        let mut stress = [0.0; 6];
        stress[axis] = strength;
        let at_limit = limit.assess_material_stress(&stress).unwrap();
        assert_eq!(at_limit.governing_axis(), axis);
        assert_eq!(at_limit.utilization(), 1.0);
        assert!(at_limit.reaches_limit());
        stress[axis] = strength * 0.99;
        assert!(
            !limit
                .assess_material_stress(&stress)
                .unwrap()
                .reaches_limit()
        );
        stress[axis] = strength * 1.01;
        assert!(
            limit
                .assess_material_stress(&stress)
                .unwrap()
                .reaches_limit()
        );
    }
    let assessment = limit
        .assess_material_stress(&[80.0, 9.0, 3.0, 4.0, 5.0, 6.0])
        .unwrap();
    assert_eq!(assessment.material_normal_stresses_pa(), [80.0, 9.0, 3.0]);
    assert_eq!(assessment.tensile_utilizations(), [0.8, 0.9, 1.5]);
    assert_eq!(assessment.governing_axis(), 2);
    assert_eq!(
        assessment.material_stress_pa(),
        [80.0, 9.0, 3.0, 4.0, 5.0, 6.0]
    );
}

#[test]
fn g0_tension_only_observation_retains_unassessed_compression_and_shear() {
    let limit = limit([100.0, 10.0, 2.0]);
    for stress in [
        [-1e6, -2e6, -3e6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1e6, -2e6, 3e6],
    ] {
        let result = limit.assess_material_stress(&stress).unwrap();
        assert_eq!(result.material_stress_pa(), stress);
        assert_eq!(result.tensile_utilizations(), [0.0; 3]);
        assert_eq!(result.governing_axis(), 0);
        assert!(!result.reaches_limit());
    }
    let mixed = limit
        .assess_material_stress(&[-1e6, 20.0, -3e6, 0.0, 0.0, 0.0])
        .unwrap();
    assert_eq!(mixed.utilization(), 2.0);
    assert_eq!(mixed.governing_axis(), 1);
    let tied = limit
        .assess_material_stress(&[50.0, 5.0, 1.0, 0.0, 0.0, 0.0])
        .unwrap();
    assert_eq!(tied.governing_axis(), 0);
}

#[test]
fn g3_material_axes_do_not_follow_the_principal_stresses() {
    let limit = limit([100.0, 10.0, 2.0]);
    let along = [20.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let across = [0.0, 20.0, 0.0, 0.0, 0.0, 0.0];
    assert_eq!(
        limit.assess_material_stress(&along).unwrap().utilization(),
        0.2
    );
    assert_eq!(
        limit.assess_material_stress(&across).unwrap().utilization(),
        2.0
    );
    let rotated = tensor::rotate(
        &along,
        &tensor::rotation([0.0, 0.0, 1.0], core::f64::consts::FRAC_PI_4),
    );
    let result = limit.assess_material_stress(&rotated).unwrap();
    assert!((result.utilization() - 1.0).abs() < 1e-14);
    assert_eq!(result.governing_axis(), 1);
}

#[test]
fn g3_common_stress_unit_scaling_preserves_utilization() {
    for scale in [1e-280, 1.0, 1e280] {
        let limit = limit([100.0 * scale, 10.0 * scale, 2.0 * scale]);
        let result = limit
            .assess_material_stress(&[50.0 * scale, 15.0 * scale, scale, 0.0, 0.0, 0.0])
            .unwrap();
        assert!((result.utilization() - 1.5).abs() < 1e-14);
        assert_eq!(result.governing_axis(), 1);
    }
}

#[test]
fn g0_refresh_preserves_axis_selection_and_source_receipts() {
    let card = card([100.0, 10.0, 2.0], Pressure::DIMS, None);
    let cold = OrthotropicTensileLimit::resolve(&card, &point(300.0), NAMES, SELECTION).unwrap();
    let warm = cold.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(cold.limits_pa(), [100.0, 10.0, 2.0]);
    assert_eq!(warm.limits_pa(), [75.0, 7.5, 1.5]);
    assert_eq!(warm.limit_properties(), NAMES);
    assert_eq!(
        cold.resolved().card_identity(),
        warm.resolved().card_identity()
    );
    for (a, b) in cold.sources().into_iter().zip(warm.sources()) {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    let stress = [80.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(
        !cold
            .assess_material_stress(&stress)
            .unwrap()
            .reaches_limit()
    );
    assert!(
        warm.assess_material_stress(&stress)
            .unwrap()
            .reaches_limit()
    );
    assert!(warm.at_query_point(&card, &point(501.0)).is_err());
    assert_eq!(warm.limits_pa(), [75.0, 7.5, 1.5]);
}

#[test]
fn g0_missing_axes_wrong_units_and_invalid_strengths_refuse() {
    let good = card([100.0, 10.0, 2.0], Pressure::DIMS, None);
    for names in [
        [NAMES[0], NAMES[1], "missing"],
        [NAMES[0], NAMES[0], NAMES[2]],
    ] {
        assert!(OrthotropicTensileLimit::resolve(&good, &point(300.0), names, SELECTION).is_err());
    }
    for axis in 0..3 {
        for invalid in [0.0, -1.0] {
            let mut values = [100.0, 10.0, 2.0];
            values[axis] = invalid;
            assert!(
                OrthotropicTensileLimit::resolve(
                    &card(values, Pressure::DIMS, None),
                    &point(300.0),
                    NAMES,
                    SELECTION,
                )
                .is_err()
            );
        }
    }
    let wrong = card([100.0, 10.0, 2.0], Dims::NONE, None);
    let requirements = NAMES.map(|name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    let resolved =
        resolve_material_state_point(&wrong, &point(300.0), &requirements, SELECTION).unwrap();
    assert!(OrthotropicTensileLimit::try_from_resolved(resolved, NAMES).is_err());
}

#[test]
fn g0_unowned_rate_frequency_and_foreign_card_refuse() {
    for axis in ["strain_rate", "omega"] {
        let card = card([100.0, 10.0, 2.0], Pressure::DIMS, Some(axis));
        let query = point(300.0).with(axis, 1.5).unwrap();
        assert!(OrthotropicTensileLimit::resolve(&card, &query, NAMES, SELECTION).is_err());
    }
    let limit = limit([100.0, 10.0, 2.0]);
    assert!(
        limit
            .at_query_point(
                &card([101.0, 10.0, 2.0], Pressure::DIMS, None),
                &point(300.0),
            )
            .is_err()
    );
}

#[test]
fn g0_invalid_stresses_and_unrepresentable_utilization_refuse() {
    let limit = limit([100.0, 10.0, 2.0]);
    for axis in 0..6 {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut stress = [0.0; 6];
            stress[axis] = invalid;
            assert!(limit.assess_material_stress(&stress).is_err());
        }
    }
    assert!(
        limit
            .assess_material_stress(&[f64::from_bits(1), 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
    let tiny = OrthotropicTensileLimit::resolve(
        &card([1e-300; 3], Pressure::DIMS, None),
        &point(300.0),
        NAMES,
        SELECTION,
    )
    .unwrap();
    assert!(
        tiny.assess_material_stress(&[f64::MAX, 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
}
