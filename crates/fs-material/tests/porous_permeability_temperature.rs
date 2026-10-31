//! Source-backed permeability fixtures; these values are not material measurements.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::porous::permeability::{
    INTRINSIC_PERMEABILITY_DIMS, temperature::PermeabilityTemperatureLaw,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};

const K_NAMES: [&str; 3] = ["intrinsic-k1", "intrinsic-k2", "intrinsic-k3"];
const AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), 300.0)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100000.0,
        )
        .unwrap()
}

fn keys() -> [PropertyKey; 3] {
    K_NAMES.map(|name| PropertyKey::new(name, INTRINSIC_PERMEABILITY_DIMS))
}

fn card(split_source: Option<&str>, first_values: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (K_NAMES[0], INTRINSIC_PERMEABILITY_DIMS, first_values),
        (K_NAMES[1], INTRINSIC_PERMEABILITY_DIMS, [8e-6, 16e-6]),
        (K_NAMES[2], INTRINSIC_PERMEABILITY_DIMS, [18e-6, 36e-6]),
        (DENSITY_PROPERTY, Density::DIMS, [2000.0, 2000.0]),
    ] {
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, values[0]), (400.0, values[1])],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", absolute_temperature(), 300.0, 400.0)
                .with_quantity(
                    "pressure",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Pressure,
                        ValueForm::Static,
                    )),
                    100000.0,
                    100000.0,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic directional permeability, not measurements".into(),
                license: "test-only".into(),
                artifact: None,
            },
        };
        if split_source == Some(name) {
            claim.validity = claim.validity.with("T", 300.0, 350.0);
            claims.insert_claim(claim.clone()).unwrap();
            claim.validity = claim.validity.with("T", 351.0, 400.0);
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic porous skeleton".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(card: &MaterialCard, point: &QueryPoint) -> ResolvedMaterialStatePoint {
    let mut requirements = K_NAMES
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                INTRINSIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap()
        })
        .to_vec();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    resolve_material_state_point(
        card,
        point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn temperature_refresh_preserves_full_source_bundle_and_rotated_tensor() {
    let card = card(None, [2e-6, 4e-6]);
    let material = resolve(&card, &point());
    let c = std::f64::consts::FRAC_1_SQRT_2;
    let axes = [[c, c, 0.0], [-c, c, 0.0], [0.0, 0.0, 1.0]];
    let law =
        PermeabilityTemperatureLaw::try_new(&card, material.clone(), keys(), axes, "T").unwrap();
    let hot = law.at_temperature(Temperature::new(350.0)).unwrap();

    for (actual, expected) in hot
        .permeability()
        .principal_m2()
        .iter()
        .zip([3e-6, 12e-6, 27e-6])
    {
        assert!((actual - expected).abs() < 1e-18);
    }
    // K_world = R diag(k1,k2,k3) R^T: a 45-degree normal sees (k1+k2)/2.
    for normal in [AXES[0], AXES[1]] {
        assert!((hot.permeability().normal_permeability(normal).unwrap() - 7.5e-6).abs() < 1e-18);
    }
    assert_eq!(hot.permeability().material_axes_in_world(), &axes);
    assert_eq!(hot.material().query(), point().with("T", 350.0).unwrap());
    assert_eq!(
        hot.material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        2000.0
    );
    for name in K_NAMES.into_iter().chain([DENSITY_PROPERTY]) {
        let receipt = &hot.material().property(name).unwrap().answer().receipt;
        assert_eq!(
            receipt.selected,
            material.property(name).unwrap().answer().receipt.selected
        );
        assert_eq!(receipt.axis_quantities, *point().axis_quantities());
        assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
    }
    assert_eq!(law.reference().material(), &material);
    assert_eq!(law.source_keys(), &keys());
    assert_eq!(law.temperature_axis(), "T");
    assert_eq!(hot, law.at_temperature(Temperature::new(350.0)).unwrap());
}

#[test]
fn every_original_claim_stays_pinned_including_other_material_properties() {
    for split_source in [K_NAMES[0], DENSITY_PROPERTY] {
        let card = card(Some(split_source), [2e-6, 4e-6]);
        let law =
            PermeabilityTemperatureLaw::try_new(&card, resolve(&card, &point()), keys(), AXES, "T")
                .unwrap();
        let reference = law.reference().clone();
        // Fresh selection can find another claim, but continuation must not switch to it.
        let fresh = resolve(&card, &point().with("T", 375.0).unwrap());
        assert_ne!(
            fresh
                .property(split_source)
                .unwrap()
                .answer()
                .receipt
                .selected,
            reference
                .material()
                .property(split_source)
                .unwrap()
                .answer()
                .receipt
                .selected,
        );
        assert!(law.at_temperature(Temperature::new(375.0)).is_err());
        assert!(law.at_temperature(Temperature::new(325.0)).is_ok());
        assert_eq!(law.reference(), &reference);
    }
}

#[test]
fn equal_principal_keys_and_exact_seals_survive_temperature_updates() {
    let card = card(None, [0.0, 0.0]);
    let material = resolve(&card, &point());
    let k = keys()[0].clone();
    let law =
        PermeabilityTemperatureLaw::try_new(&card, material, [k.clone(), k.clone(), k], AXES, "T")
            .unwrap();
    let hot = law.at_temperature(Temperature::new(400.0)).unwrap();
    assert_eq!(hot.permeability().principal_m2(), &[0.0; 3]);
    assert_eq!(
        hot.permeability().normal_permeability(AXES[0]).unwrap(),
        0.0
    );
}

#[test]
fn wrong_card_quantity_frame_and_unsupported_temperature_refuse() {
    let card = card(None, [2e-6, 4e-6]);
    let material = resolve(&card, &point());
    let other_card = self::card(None, [1e-6, 2e-6]);
    assert!(
        PermeabilityTemperatureLaw::try_new(&other_card, material.clone(), keys(), AXES, "T")
            .is_err()
    );
    let mut wrong_keys = keys();
    wrong_keys[0] = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    assert!(
        PermeabilityTemperatureLaw::try_new(&card, material.clone(), wrong_keys, AXES, "T")
            .is_err()
    );
    assert!(
        PermeabilityTemperatureLaw::try_new(&card, material.clone(), keys(), [AXES[0]; 3], "T")
            .is_err()
    );
    assert!(
        PermeabilityTemperatureLaw::try_new(
            &card,
            material.clone(),
            keys(),
            AXES,
            "missing-temperature"
        )
        .is_err()
    );

    let law =
        PermeabilityTemperatureLaw::try_new(&card, material.clone(), keys(), AXES, "T").unwrap();
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY, 299.0, 401.0] {
        assert!(law.at_temperature(Temperature::new(temperature)).is_err());
    }
    assert_eq!(law.reference().material(), &material);
    assert!(law.at_temperature(Temperature::new(350.0)).is_ok());
}
