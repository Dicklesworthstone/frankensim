//! G0/G3: local tensile limits from actual stresses and original source states.
//! Synthetic coefficients verify the model; no measured material is implied.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::strength::{IsotropicTensileLimit, principal_stresses};
use fs_material::tensor;
use fs_qty::{Dims, Pressure};

const LIMIT: &str = "specified_tensile_strength";

fn card(value: PropertyValue, validity: ValidityDomain) -> MaterialCard {
    let dims = match &value {
        PropertyValue::Scalar { dims, .. } | PropertyValue::Curve { dims, .. } => *dims,
        _ => unreachable!(),
    };
    let interpolation = if matches!(value, PropertyValue::Curve { .. }) {
        InterpolationPolicy::LinearInside
    } else {
        InterpolationPolicy::ConstantWithinValidity
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(LIMIT, dims),
        value,
        validity,
        interpolation,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic tensile assessment".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "test specimen".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap()
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
}

fn constant(strength: f64) -> MaterialCard {
    card(
        PropertyValue::Scalar { value: strength, dims: Pressure::DIMS },
        ValidityDomain::unconstrained().with("T", 280.0, 500.0),
    )
}

fn limit(strength: f64) -> IsotropicTensileLimit {
    IsotropicTensileLimit::resolve(
        &constant(strength),
        &point(300.0),
        LIMIT,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 5e-14 * scale,
        "actual {actual:e}, expected {expected:e}");
}

#[test]
fn g0_tension_cutoff_includes_hydrostatic_tension() {
    let limit = limit(4.0e8);
    for stress in [
        [4.0e8, 0.0, 0.0, 0.0, 0.0, 0.0],
        [4.0e8, 4.0e8, 4.0e8, 0.0, 0.0, 0.0],
    ] {
        let result = limit.assess(&stress).unwrap();
        assert_eq!(result.tensile_stress_pa(), 4.0e8);
        assert_eq!(result.utilization(), 1.0);
        assert!(result.reaches_limit());
    }
    assert_eq!(tensor::von_mises(&[4.0e8, 4.0e8, 4.0e8, 0.0, 0.0, 0.0]), 0.0);
    assert!(!limit.assess(&[3.9e8, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap().reaches_limit());
    assert!(limit.assess(&[4.1e8, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap().reaches_limit());
}

#[test]
fn g0_compression_does_not_consume_tension_and_shear_uses_tensor_convention() {
    let limit = limit(4.0e8);
    for stress in [
        [-1e9, -2e9, -3e9, 0.0, 0.0, 0.0],
        [-1e9, 0.0, 0.0, 0.0, 0.0, 0.0],
    ] {
        let result = limit.assess(&stress).unwrap();
        assert_eq!(result.tensile_stress_pa(), 0.0);
        assert_eq!(result.utilization(), 0.0);
        assert!(!result.reaches_limit());
    }
    for component in 3..6 {
        for sign in [-1.0, 1.0] {
            let mut stress = [0.0; 6];
            stress[component] = sign * 2.0e8;
            let result = limit.assess(&stress).unwrap();
            assert_eq!(result.principal_stresses_pa(), [-2.0e8, 0.0, 2.0e8]);
            assert_eq!(result.utilization(), 0.5);
        }
    }
}

#[test]
fn g3_rotation_and_repeated_principal_stresses_preserve_tensile_demand() {
    let limit = limit(4.0e8);
    let axis = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    for diagonal in [[-3e8, 2e8, 3e8], [-3e8, 2e8, 2e8], [2e8, 2e8, 2e8]] {
        let stress = [diagonal[0], diagonal[1], diagonal[2], 0.0, 0.0, 0.0];
        for angle in [-2.1, -0.25, 0.7, 2.3] {
            let rotated = tensor::rotate(&stress, &tensor::rotation(axis, angle));
            let result = limit.assess(&rotated).unwrap();
            for (actual, expected) in result.principal_stresses_pa().into_iter().zip(diagonal) {
                near(actual, expected, 4e8);
            }
            near(result.utilization(), diagonal[2] / 4e8, 1.0);
        }
    }
}

#[test]
fn g0_scaled_eigenvalues_do_not_square_dimensional_stress() {
    for scale in [1e-300, 1.0, 1e300] {
        let result = limit(scale).assess(&[0.0, 0.0, 0.0, scale, 0.0, 0.0]).unwrap();
        assert_eq!(result.principal_stresses_pa(), [-scale, 0.0, scale]);
        assert_eq!(result.utilization(), 1.0);
    }
    assert_eq!(
        principal_stresses(&[f64::MAX, f64::MAX, f64::MAX, 0.0, 0.0, 0.0]).unwrap(),
        [f64::MAX; 3],
    );
}

#[test]
fn g0_temperature_refresh_preserves_selected_claim_and_changes_capacity() {
    let card = card(
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, 4e8), (500.0, 2e8)],
            dims: Pressure::DIMS,
        },
        ValidityDomain::unconstrained().with("T", 300.0, 500.0),
    );
    let cold = IsotropicTensileLimit::resolve(
        &card, &point(300.0), LIMIT, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let warm = cold.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(cold.limit_pa(), 4e8);
    assert_eq!(warm.limit_pa(), 3e8);
    assert_eq!(cold.limit_property(), LIMIT);
    assert_eq!(cold.resolved().card_identity(), warm.resolved().card_identity());
    assert_eq!(cold.source().answer().receipt.selected, warm.source().answer().receipt.selected);
    assert!(!cold.assess(&[3.5e8, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap().reaches_limit());
    assert!(warm.assess(&[3.5e8, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap().reaches_limit());
    assert!(warm.at_query_point(&card, &point(501.0)).is_err());
    assert!(warm.at_query_point(&constant(3e8), &point(400.0)).is_err());
    assert_eq!(warm.limit_pa(), 3e8);
}

#[test]
fn g0_unowned_rate_and_frequency_conditions_refuse() {
    for axis in ["strain_rate", "omega"] {
        let card = card(
            PropertyValue::Scalar { value: 4e8, dims: Pressure::DIMS },
            ValidityDomain::unconstrained().with("T", 280.0, 500.0).with(axis, 1.0, 2.0),
        );
        let query = point(300.0).with(axis, 1.5).unwrap();
        assert!(IsotropicTensileLimit::resolve(
            &card, &query, LIMIT, MaterialPropertySelection::SingleClaimOnly,
        ).is_err());
    }
}

#[test]
fn g0_source_admission_never_fills_missing_or_invalid_strength() {
    for strength in [0.0, -1.0] {
        assert!(IsotropicTensileLimit::resolve(
            &constant(strength), &point(300.0), LIMIT, MaterialPropertySelection::SingleClaimOnly,
        ).is_err());
    }
    assert!(IsotropicTensileLimit::resolve(
        &constant(4e8), &point(300.0), "missing", MaterialPropertySelection::SingleClaimOnly,
    ).is_err());
    let wrong = card(
        PropertyValue::Scalar { value: 1.0, dims: Dims::NONE },
        ValidityDomain::unconstrained().with("T", 280.0, 500.0),
    );
    let requirement = ScalarPropertyRequirement::try_new(
        LIMIT, Dims::NONE, ScalarAdmissibility::StrictlyPositive,
    ).unwrap();
    let resolved = resolve_material_state_point(
        &wrong, &point(300.0), &[requirement], MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    assert!(IsotropicTensileLimit::try_from_resolved(resolved, LIMIT).is_err());
}

#[test]
fn g0_nonfinite_stress_and_unrepresentable_results_refuse() {
    let limit = limit(1e-300);
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(limit.assess(&[invalid, 0.0, 0.0, 0.0, 0.0, 0.0]).is_err());
    }
    assert!(limit.assess(&[f64::MAX, 0.0, 0.0, 0.0, 0.0, 0.0]).is_err());
    assert!(principal_stresses(&[f64::MAX; 6]).is_err());
    assert_eq!(limit.assess(&[0.0; 6]).unwrap().utilization(), 0.0);
}

#[test]
fn g0_small_tensile_stress_is_not_hidden_by_large_compression() {
    for (compression, shear) in [(1.0e16, 1.0), (1.0e300, 1.0e-300)] {
        for (normal, transverse_shear) in [(0, 4), (1, 5), (2, 3)] {
            let mut stress = [0.0; 6];
            stress[normal] = -compression;
            stress[transverse_shear] = shear;
            let result = limit(shear / 2.0).assess(&stress).unwrap();
            assert_eq!(result.principal_stresses_pa(), [-compression, -shear, shear]);
            assert_eq!(result.utilization(), 2.0);
            assert!(result.reaches_limit());
        }
    }
    let result = limit(5.0e-201)
        .assess(&[-1.0e300, 0.0, 0.0, 1.0e50, 0.0, 0.0])
        .unwrap();
    near(result.utilization(), 2.0, 2.0);
    assert!(result.reaches_limit());
}
