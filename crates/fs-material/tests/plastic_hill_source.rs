//! G0/G3: explicit directional coefficients drive the original material law.
//! These are authored numerical fixtures, not measured rolling/formability data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES as C, HILL_HARDENING_MODULUS_PROPERTY as H,
    HILL_REFERENCE_YIELD_STRESS_PROPERTY as Y, ResolvedHill48, hill48_requirements, resolve_hill48,
};
use fs_material::plastic::hill::{Hill48, Hill48Coefficients};
use fs_material::plastic::source::ISOTROPIC_HARDENING_MODULUS_PROPERTY;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_qty::{Density, Dims, Pressure, Temperature};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const VALUES: [f64; 6] = [0.35, 0.8, 0.45, 1.5, 1.1, 2.0];

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic Hill48 coefficients".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn properties() -> Vec<PropertyClaim> {
    let mut properties = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.3),
        scalar(Y, Pressure::DIMS, 1.0),
        scalar(H, Pressure::DIMS, 60.0),
    ];
    properties.extend(
        C.into_iter()
            .zip(VALUES)
            .map(|(name, value)| scalar(name, Dims([0; 6]), value)),
    );
    properties
}
fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored rolled metal".into(),
            phase: "solid".into(),
            process: "numerical Hill48 fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(card: &MaterialCard) -> ResolvedHill48 {
    resolve_hill48(
        card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.04,
        IDENTITY,
    )
    .unwrap()
}

#[test]
fn original_source_bundle_drives_directional_compression_shear_and_unload() {
    let source = card(properties());
    let resolved = resolve(&source);
    let [f, g, h, l, m, n] = VALUES;
    let native = Hill48::new(
        IsotropicElastic::new(1000.0, 0.3, 0.04).unwrap(),
        1.0,
        60.0,
        Hill48Coefficients { f, g, h, l, m, n },
    )
    .unwrap();
    assert_eq!(resolved.local_law(), &native);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    assert_eq!(resolved.material().properties().len(), 10);
    assert_eq!(resolved.reference_yield_source().value_si(), 1.0);
    assert_eq!(resolved.hardening_source().value_si(), 60.0);
    for property in resolved.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let mut state = resolved.initial_state();
    for strain in [
        [-0.01, 0.003, 0.002, 0.001, 0.0, 0.0],
        [-0.008, 0.002, 0.001, 0.003, 0.001, 0.0],
        [0.002, -0.003, 0.001, 0.0, 0.001, 0.0],
    ] {
        let response = resolved.response(&strain, &state).unwrap();
        assert_eq!(response, native.response(&strain, &state).unwrap());
        assert_eq!(
            resolved.tangent(&strain, &state),
            native.tangent(&strain, &state)
        );
        state = response.state;
    }
    assert!(state.alpha > 0.0);
}

#[test]
fn fixed_temperature_and_orientation_keep_all_other_physics_receipts() {
    let mut data = properties();
    data[2].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1.0), (400.0, 0.8)],
        dims: Pressure::DIMS,
    };
    data[2].validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
    data[2].interpolation = InterpolationPolicy::LinearInside;
    data.push(scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0));
    let source = card(data);
    let mut requirements = hill48_requirements().unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let point = QueryPoint::new().with("T", 350.0).unwrap();
    let bundle = resolve_material_state_point(
        &source,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let rotation = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.47);
    let resolved = ResolvedHill48::try_from_material(bundle.clone(), 0.04, rotation).unwrap();
    assert_eq!(resolved.material(), &bundle);
    assert_eq!(resolved.local_law().yield_stress_pa(), 0.9);
    assert_eq!(resolved.material().properties().len(), 11);
    assert_eq!(resolved.law().material_to_spatial(), &rotation);
    let independent = OrientedHill48::new(resolved.local_law().clone(), rotation).unwrap();
    let strain = [0.007, -0.002, -0.001, 0.001, 0.0, 0.0];
    assert_eq!(
        resolved.response(&strain, &resolved.initial_state()),
        independent.response(&strain, &independent.initial_state())
    );
    assert!(
        resolve_hill48(
            &source,
            &QueryPoint::new().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            0.04,
            rotation
        )
        .is_err()
    );
}

#[test]
fn uniaxial_j2_or_incomplete_data_cannot_replace_hill_calibration() {
    for index in [2, 3, 4, 9] {
        let mut data = properties();
        data.remove(index);
        assert!(
            resolve_hill48(
                &card(data),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                0.04,
                IDENTITY
            )
            .is_err()
        );
    }
    for (index, name, dims, value) in [
        (2, YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0),
        (
            3,
            ISOTROPIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            60.0,
        ),
        (4, C[0], Pressure::DIMS, 0.35),
        (9, C[5], Dims([0; 6]), -1.0),
        (4, C[0], Dims([0; 6]), -2.0),
    ] {
        let mut data = properties();
        data[index] = scalar(name, dims, value);
        assert!(
            resolve_hill48(
                &card(data),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                0.04,
                IDENTITY
            )
            .is_err()
        );
    }
    // Negative F is valid when the complete deviatoric quadratic stays positive.
    let mut data = properties();
    data[4] = scalar(C[0], Dims([0; 6]), -0.1);
    assert!(resolve(&card(data)).local_law().coefficients().f < 0.0);
}

#[test]
fn evolving_anisotropy_and_ambiguous_sources_refuse_without_hidden_defaults() {
    let mut data = properties();
    data[4].value = PropertyValue::Curve {
        abscissa: "strain".into(),
        abscissa_dims: Dims([0; 6]),
        knots: vec![(0.0, 0.35), (0.02, 0.6)],
        dims: Dims([0; 6]),
    };
    data[4].validity = ValidityDomain::unconstrained().with("strain", 0.0, 0.02);
    data[4].interpolation = InterpolationPolicy::LinearInside;
    let bundle = resolve_material_state_point(
        &card(data),
        &QueryPoint::new().with("strain", 0.01).unwrap(),
        &hill48_requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(ResolvedHill48::try_from_material(bundle, 0.04, IDENTITY).is_err());
    let mut ambiguous = properties();
    ambiguous.push(scalar(C[0], Dims([0; 6]), 0.6));
    assert!(
        resolve_hill48(
            &card(ambiguous),
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            0.04,
            IDENTITY
        )
        .is_err()
    );
}
