//! G0/G3: ordered interface data reach the native fracture law without defaults.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    ResolvedElectricalInterface,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, interface_stiffness_per_area_dims, resolve_interface_state_point,
};
use fs_qty::{Dims, Length, Pressure, Temperature};

fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured cohesive interface; not measured adhesive data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn properties() -> Vec<PropertyClaim> {
    vec![
        claim(
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        claim(COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 1e6),
        claim(
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            2000.0,
        ),
        claim(
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ]
}

fn card(properties: Vec<PropertyClaim>, history: &str) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    let surface = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.into(),
            phase: "solid".into(),
            process: "manufactured surface".into(),
            revision: 1,
        },
        texture_frame: format!("{chemistry} reference surface"),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit cohesive bond".into(),
            third_body: None,
            environment: "fixed external conditions".into(),
            history: history.into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(card: &InterfaceSystemCard) -> ResolvedBilinearCohesive {
    resolve_bilinear_cohesive(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap()
}

fn bundle(card: &InterfaceSystemCard, properties: &[PropertyClaim]) -> ResolvedInterfaceStatePoint {
    let requirements: Vec<_> = properties
        .iter()
        .map(|p| {
            ScalarPropertyRequirement::try_with_key(&p.key, ScalarAdmissibility::Finite).unwrap()
        })
        .collect();
    resolve_interface_state_point(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn ordered_source_bundle_is_shared_by_failure_and_electrical_response() {
    let mut properties = properties();
    properties.push(claim(
        ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
        ELECTRICAL_CONTACT_RESISTANCE_DIMS,
        2e-8,
    ));
    let card = card(properties.clone(), "original bond");
    let resolved = bundle(&card, &properties);
    let electrical =
        ResolvedElectricalInterface::try_from_resolved(&card, resolved.clone()).unwrap();
    let mechanical =
        ResolvedBilinearCohesive::try_from_resolved(&card, resolved.clone(), 0.02).unwrap();
    assert_eq!(mechanical.resolved(), electrical.resolved());
    assert_eq!(mechanical.resolved(), &resolved);
    assert_eq!(
        mechanical.surfaces(),
        &[card.surface_a().clone(), card.surface_b().clone()]
    );
    mechanical
        .require_surface_materials(&card.surface_a().material, &card.surface_b().material)
        .unwrap();
    assert!(
        mechanical
            .require_surface_materials(&card.surface_b().material, &card.surface_a().material)
            .is_err()
    );
    let failed = mechanical
        .law()
        .bonded_point()
        .trial([0.005, 0.0, 0.0], 1e-8)
        .unwrap();
    assert!(failed.point.failed());
    assert!((failed.energy.dissipated - 2000.0).abs() < 1e-9);
    assert_eq!(mechanical.resolved(), &resolved);
}

#[test]
fn fixed_temperature_data_change_damage_onset_and_failure_work() {
    let mut properties = properties();
    properties[1].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 1.2e6), (340.0, 0.6e6)],
        dims: Pressure::DIMS,
    };
    properties[1].interpolation = InterpolationPolicy::LinearInside;
    let card = card(properties, "unchanged preparation");
    let cold = resolve(&card);
    let warm = resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new().with("T", 320.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap();
    assert_eq!(
        cold.resolved().card_identity(),
        warm.resolved().card_identity()
    );
    assert_ne!(cold.resolved().identity(), warm.resolved().identity());
    assert!(warm.law().onset_separation_m() < cold.law().onset_separation_m());
    assert!(warm.law().failure_separation_m() > cold.law().failure_separation_m());
    let cold_point = cold
        .law()
        .bonded_point()
        .trial([0.0009, 0.0, 0.0], 1e-8)
        .unwrap();
    let warm_point = warm
        .law()
        .bonded_point()
        .trial([0.0009, 0.0, 0.0], 1e-8)
        .unwrap();
    assert_eq!(cold_point.point.damage(), 0.0);
    assert!(warm_point.point.damage() > 0.0);
    for source in [&cold, &warm] {
        let failed = source
            .law()
            .bonded_point()
            .trial([0.01, 0.0, 0.0], 1e-8)
            .unwrap();
        assert!((failed.energy.dissipated - 2000.0).abs() < 1e-9);
    }
}

#[test]
fn wrong_cards_missing_units_and_invalid_envelopes_refuse() {
    let original = card(properties(), "original");
    let other = card(properties(), "different surface preparation history");
    let resolved = bundle(&original, &properties());
    assert!(ResolvedBilinearCohesive::try_from_resolved(&other, resolved.clone(), 0.02).is_err());
    for limit in [0.0, f64::NAN, 0.001] {
        assert!(
            ResolvedBilinearCohesive::try_from_resolved(&original, resolved.clone(), limit)
                .is_err()
        );
    }
    for kind in 0..4 {
        let mut p = properties();
        match kind {
            0 => {
                p.pop();
            }
            1 => {
                p[1] = claim(COHESIVE_STRENGTH_PROPERTY, Dims::NONE, 1e6);
            }
            2 => {
                p[0] = claim(
                    COHESIVE_STIFFNESS_PROPERTY,
                    interface_stiffness_per_area_dims(),
                    0.0,
                );
            }
            _ => {
                p[2] = claim(
                    COHESIVE_FRACTURE_ENERGY_PROPERTY,
                    FRACTURE_ENERGY_PER_AREA_DIMS,
                    100.0,
                );
            }
        }
        let changed = card(p, "invalid input fixture");
        assert!(
            resolve_bilinear_cohesive(
                &changed,
                &QueryPoint::new().with("T", 300.0).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
                0.02
            )
            .is_err()
        );
    }
}

#[test]
fn evolving_gap_and_rate_applicability_cannot_be_frozen_into_damage_law() {
    for axis in ["gap", "damage", "sliding_rate", "omega"] {
        let mut p = properties();
        p[0].validity = p[0].validity.clone().with(axis, 0.0, 1.0);
        let card = card(p, "bounded evolving condition");
        let point = QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with(axis, 0.1)
            .unwrap();
        assert!(
            resolve_bilinear_cohesive(
                &card,
                &point,
                MaterialPropertySelection::SingleClaimOnly,
                0.02
            )
            .is_err(),
            "{axis}"
        );
    }
    let mut p = properties();
    p[0].value = PropertyValue::Curve {
        abscissa: "gap".into(),
        abscissa_dims: Length::DIMS,
        knots: vec![(0.0, 1e9), (0.01, 2e9)],
        dims: interface_stiffness_per_area_dims(),
    };
    p[0].interpolation = InterpolationPolicy::LinearInside;
    p[0].validity = p[0].validity.clone().with("gap", 0.0, 0.01);
    let card = card(p, "separation-dependent coefficients");
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("gap", 0.0)
        .unwrap();
    assert!(
        resolve_bilinear_cohesive(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.02
        )
        .is_err()
    );
}
