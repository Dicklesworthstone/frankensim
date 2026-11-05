//! Source meaning and actual attenuation in the native ductile J2 update.
//! All coefficients are synthetic verification data, not an alloy calibration.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::ductile::evolution::DuctileJ2Conditions;
use fs_material::ductile::evolution::source::{
    DuctileJ2SourceError, DuctileJ2SourceSettings, PLASTIC_WORK_PER_AREA_DIMS,
    POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY, ResolvedDuctileJ2, ductile_j2_requirements,
    resolve_ductile_j2,
};
use fs_material::ductile::source::{
    DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY, PLASTIC_STRAIN_RATE_AXIS,
    ResolvedDuctileLocus, STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{HardeningPropertyBasis, resolve_j2_plasticity};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const TEMPERATURE: &str = "specimen-temperature";
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const WORK_INDEX: usize = 4;

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn fixed_point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(TEMPERATURE, temperature_quantity(), 300.0)
        .unwrap()
}

fn fixed_domain() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(TEMPERATURE, temperature_quantity(), 290.0, 310.0)
}

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: fixed_domain(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic ductile J2 source integration fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn claims(work_j_m2: f64) -> Vec<PropertyClaim> {
    let mut locus = scalar(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE, 0.002);
    locus.value = PropertyValue::Curve {
        abscissa: STRESS_TRIAXIALITY_AXIS.into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(-1.0, 0.002), (0.0, 0.002), (1.0, 0.002)],
        dims: Dims::NONE,
    };
    locus.interpolation = InterpolationPolicy::TabulatedOnly;
    locus.validity = fixed_domain()
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
            1.0,
        )
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE),
            1e-8,
            1.0,
        );
    vec![
        scalar("young_modulus", Pressure::DIMS, 1200.0),
        scalar("poisson_ratio", Dims::NONE, 0.2),
        scalar("yield_stress", Pressure::DIMS, 3.0),
        scalar("isotropic_hardening_modulus", Pressure::DIMS, 100.0),
        scalar(
            POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY,
            PLASTIC_WORK_PER_AREA_DIMS,
            work_j_m2,
        ),
        locus,
    ]
}

fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for row in rows {
        set.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic ductile coupon".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 1,
        },
        set,
        vec![],
    )
    .unwrap()
}

fn initiation(source: &MaterialCard) -> ResolvedDuctileLocus {
    let point = fixed_point()
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE),
            0.01,
        )
        .unwrap();
    resolve_ductile_locus(
        source,
        &point,
        &DuctileSourceRequest {
            temperature_axis: TEMPERATURE.into(),
            maximum_points: 3,
        },
        SINGLE,
    )
    .unwrap()
}

fn settings(width: f64) -> DuctileJ2SourceSettings {
    DuctileJ2SourceSettings {
        characteristic_length_m: width,
        maximum_damage_increment: 0.95,
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-10,
            yield_pa: 1e-8,
        },
    }
}

fn conditions() -> DuctileJ2Conditions {
    DuctileJ2Conditions {
        duration_s: 1.0,
        temperature_k: 300.0,
    }
}

fn resolve(source: &MaterialCard, width: f64) -> Result<ResolvedDuctileJ2, DuctileJ2SourceError> {
    resolve_ductile_j2(
        source,
        &fixed_point(),
        MaterialPropertySelection::SingleClaimOnly,
        BASIS,
        0.05,
        initiation(source),
        settings(width),
    )
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}",
    );
}

fn shear_retention(resolved: &ResolvedDuctileJ2) -> f64 {
    let strain = [0.0, 0.0, 0.0, 0.01, 0.0, 0.0];
    let law = resolved.small_strain_law(conditions()).unwrap();
    let stress = law.stress(&strain, &law.initial_state());
    let native = resolved.plastic();
    let undamaged = native.stress(&strain, &native.initial_state());
    assert!(stress.iter().all(|value| value.is_finite()));
    stress[3] / undamaged[3]
}

#[test]
fn original_bundle_drives_native_attenuation_and_retains_every_receipt() {
    let mut rows = claims(0.001);
    rows.push(scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0));
    let source = card(rows);
    let mut requirements = ductile_j2_requirements(BASIS).unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &source,
        &fixed_point(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = bundle.clone();
    let resolved = ResolvedDuctileJ2::try_from_material(
        bundle,
        BASIS,
        0.05,
        initiation(&source),
        settings(0.01),
    )
    .unwrap();
    assert_eq!(resolved.material(), &original);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    assert_eq!(resolved.law().native(), resolved.plastic().law());
    assert_eq!(
        resolved.law().initiation_law().source().card_identity(),
        source.content_hash(),
    );
    for property in resolved.material().properties() {
        source
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert_eq!(
        resolved.law().settings().post_initiation_plastic_work_j_m2,
        0.001,
    );
    // Independent pure-shear radial return: mu=500, H=100, sigma_y=3.
    // Only the post-initiation segment alpha>0.002 enters the work exponent.
    let alpha = (3.0_f64.sqrt() * 10.0 - 3.0) / 1600.0;
    let work_density = 3.0 * (alpha - 0.002) + 50.0 * (alpha * alpha - 0.002_f64.powi(2));
    near(
        shear_retention(&resolved),
        (-0.01 * work_density / 0.001).exp(),
    );
}

#[test]
fn source_work_and_original_band_width_change_response_without_changing_j2() {
    let low = resolve(&card(claims(0.001)), 0.01).unwrap();
    let high = resolve(&card(claims(0.002)), 0.01).unwrap();
    let wider = resolve(&card(claims(0.001)), 0.02).unwrap();
    assert_eq!(low.law().native(), high.law().native());
    assert_eq!(low.law().native(), wider.law().native());
    let retention = shear_retention(&low);
    assert!(retention > 0.0 && retention < 1.0);
    assert!(shear_retention(&high) > retention);
    near(shear_retention(&high).powi(2), retention);
    near(shear_retention(&wider), retention.powi(2));
}

#[test]
fn absent_source_wrong_dimensions_and_other_fracture_properties_never_substitute() {
    let mut absent = claims(0.001);
    absent.remove(WORK_INDEX);
    assert!(resolve(&card(absent), 0.01).is_err());
    for name in [
        "cohesive_fracture_energy",
        "brittle_mode_i_fracture_energy",
        "ductile_fracture_energy",
    ] {
        let mut rows = claims(0.001);
        rows[WORK_INDEX] = scalar(name, PLASTIC_WORK_PER_AREA_DIMS, 0.001);
        assert!(resolve(&card(rows), 0.01).is_err());
    }
    for dims in [Pressure::DIMS, Dims([2, 1, -2, 0, 0, 0]), Dims::NONE] {
        let mut rows = claims(0.001);
        rows[WORK_INDEX] = scalar(
            POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY,
            dims,
            0.001,
        );
        assert!(resolve(&card(rows), 0.01).is_err());
    }
    assert!(resolve(&card(claims(0.0)), 0.01).is_err());
    // A J2 resolution that omitted Gp cannot borrow an unretained card claim.
    let source = card(claims(0.001));
    let only_j2 = resolve_j2_plasticity(
        &source,
        &fixed_point(),
        MaterialPropertySelection::SingleClaimOnly,
        BASIS,
        0.05,
    )
    .unwrap();
    assert!(matches!(
        ResolvedDuctileJ2::try_from_j2(only_j2, initiation(&source), settings(0.01)),
        Err(DuctileJ2SourceError::Invalid(_)),
    ));
}

#[test]
fn card_temperature_and_extra_fixed_context_must_match_before_preparation() {
    let source = card(claims(0.001));
    let foreign = card(claims(0.002));
    let pressure = fixed_point()
        .with_quantity(
            "ambient-pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            101325.0,
        )
        .unwrap();
    for (mechanical_card, point) in [
        (&foreign, fixed_point()),
        (&source, fixed_point().with(TEMPERATURE, 301.0).unwrap()),
        (&source, pressure),
    ] {
        assert!(matches!(
            resolve_ductile_j2(
                mechanical_card,
                &point,
                MaterialPropertySelection::SingleClaimOnly,
                BASIS,
                0.05,
                initiation(&source),
                settings(0.01),
            ),
            Err(DuctileJ2SourceError::Invalid(_)),
        ));
    }
    let resolved = resolve(&source, 0.01).unwrap();
    assert!(
        resolved
            .small_strain_law(DuctileJ2Conditions {
                temperature_k: 301.0,
                ..conditions()
            })
            .is_err()
    );
}

#[test]
fn evolving_calibration_restrictions_and_invalid_geometry_refuse() {
    let mut rows = claims(0.001);
    rows[WORK_INDEX].validity = fixed_domain().with_quantity(
        PLASTIC_STRAIN_RATE_AXIS,
        QuantitySpec::dimensional(RATE),
        1e-8,
        1.0,
    );
    let rate_bounded = card(rows);
    let point = fixed_point()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE),
            0.01,
        )
        .unwrap();
    assert!(matches!(
        resolve_ductile_j2(
            &rate_bounded,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            BASIS,
            0.05,
            initiation(&rate_bounded),
            settings(0.01),
        ),
        Err(DuctileJ2SourceError::Invalid(_)),
    ));
    let source = card(claims(0.001));
    for width in [0.0, -0.01, f64::NAN] {
        assert!(matches!(
            resolve(&source, width),
            Err(DuctileJ2SourceError::Native(_))
        ));
    }
}
