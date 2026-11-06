//! Original r-value provenance follows calibration into the same solid law.
//! All observations here are synthetic numerical fixtures.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::plastic::hill::Hill48Coefficients;
use fs_material::plastic::hill::source::planar::{
    HILL_PLANAR_R_PROPERTIES as R, HILL_TRANSVERSE_SHEAR_RATIO_PROPERTIES as S,
    planar_r_requirements, resolve_hill48_from_planar_r_values,
};
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY as H,
    HILL_REFERENCE_YIELD_STRESS_PROPERTY as Y, ResolvedHill48,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure, Temperature};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const PLANAR: [f64; 3] = [1.8, 1.1, 2.2];
const SHEAR: [f64; 2] = [0.8, 1.3];

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: format!("synthetic original measurement: {name}"),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn properties() -> Vec<PropertyClaim> {
    let mut values = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 210e9),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.3),
        scalar(Y, Pressure::DIMS, 250e6),
        scalar(H, Pressure::DIMS, 1.2e9),
    ];
    values.extend(
        R.into_iter()
            .zip(PLANAR)
            .chain(S.into_iter().zip(SHEAR))
            .map(|(name, value)| scalar(name, Dims([0; 6]), value)),
    );
    values
}

fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic sheet metal".into(),
            phase: "solid".into(),
            process: "original planar strain-ratio fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn calibrated_native_coefficients_retain_original_r_value_and_shear_receipts() {
    let source = card(properties());
    let law = resolve_hill48_from_planar_r_values(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.04,
        IDENTITY,
    )
    .unwrap();
    assert_eq!(law.material().card_identity(), source.content_hash());
    assert_eq!(law.material().properties().len(), 9);
    assert_eq!(
        law.local_law().coefficients(),
        &Hill48Coefficients::from_planar_r_values(PLANAR, SHEAR).unwrap()
    );
    assert_eq!(law.reference_yield_source().value_si(), 250e6);
    assert_eq!(law.hardening_source().value_si(), 1.2e9);
    for property in law.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    for (name, value) in R.into_iter().zip(PLANAR).chain(S.into_iter().zip(SHEAR)) {
        assert_eq!(law.material().property(name).unwrap().value_si(), value);
    }
    for name in HILL_COEFFICIENT_PROPERTIES {
        assert!(law.material().property(name).is_none());
    }
    let strain = [-0.006, 0.002, 0.001, 0.001, 0.0, 0.0];
    let response = law.response(&strain, &law.initial_state()).unwrap();
    assert_eq!(
        response,
        law.local_law()
            .response(&strain, &law.initial_state())
            .unwrap()
    );
    assert!(response.state.alpha > 0.0);
}

#[test]
fn fixed_query_orientation_and_extra_physics_receipts_survive_calibration() {
    let mut data = properties();
    data[4].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1.6), (400.0, 2.0)],
        dims: Dims([0; 6]),
    };
    data[4].validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
    data[4].interpolation = InterpolationPolicy::LinearInside;
    data.push(scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0));
    let source = card(data);
    let mut requirements = planar_r_requirements().unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let query = QueryPoint::new().with("T", 350.0).unwrap();
    let material = resolve_material_state_point(
        &source,
        &query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let rotation = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.61);
    let law = ResolvedHill48::try_from_planar_r_material(material.clone(), 0.04, rotation).unwrap();
    assert_eq!(law.material(), &material);
    assert_eq!(law.law().material_to_spatial(), &rotation);
    assert_eq!(law.material().properties().len(), 10);
    assert_eq!(law.material().property(R[0]).unwrap().value_si(), 1.8);
    let strain = [-0.006, 0.002, 0.001, 0.001, 0.0, 0.0];
    assert_eq!(
        law.response(&strain, &law.initial_state()),
        law.law().response(&strain, &law.initial_state())
    );
    for property in law.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
}

#[test]
fn missing_shear_wrong_units_and_evolving_strain_ratios_refuse() {
    for missing in 4..9 {
        let mut data = properties();
        data.remove(missing);
        assert!(
            resolve_hill48_from_planar_r_values(
                &card(data),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                0.04,
                IDENTITY,
            )
            .is_err()
        );
    }
    let mut wrong_units = properties();
    wrong_units[4] = scalar(R[0], Pressure::DIMS, 1.8);
    assert!(
        resolve_hill48_from_planar_r_values(
            &card(wrong_units),
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            0.04,
            IDENTITY,
        )
        .is_err()
    );
    let mut evolving = properties();
    evolving[5].value = PropertyValue::Curve {
        abscissa: "strain".into(),
        abscissa_dims: Dims([0; 6]),
        knots: vec![(0.0, 1.0), (0.02, 1.2)],
        dims: Dims([0; 6]),
    };
    evolving[5].validity = ValidityDomain::unconstrained().with("strain", 0.0, 0.02);
    evolving[5].interpolation = InterpolationPolicy::LinearInside;
    let material = resolve_material_state_point(
        &card(evolving),
        &QueryPoint::new().with("strain", 0.01).unwrap(),
        &planar_r_requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(ResolvedHill48::try_from_planar_r_material(material, 0.04, IDENTITY).is_err());
}
