//! Free-sided sourced thermocreep checks for fiber/member consumption.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::Uniaxial;
use fs_material::creep::CreepTolerances;
use fs_material::creep::plastic::source::NortonJ2Source;
use fs_material::creep::source::*;
use fs_material::creep::thermal::expansion::ExpandingCreepSource;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY, POISSON_RATIO_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const E: f64 = 1.0e9;
const NU: f64 = 0.25;
const ALPHA: f64 = 2.0e-5;

fn temp() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", temp(), t).unwrap()
}
fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity("T", temp(), 300.0, 500.0),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured uniaxial thermo-creep fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims::NONE, NU),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, 2.0e6),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 1.0e-5),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims::NONE, 1.0),
        (CREEP_ACTIVATION_ENERGY_PROPERTY, CREEP_ACTIVATION_ENERGY_DIMS, 0.0),
        (CREEP_REFERENCE_TEMPERATURE_PROPERTY, Temperature::DIMS, 300.0),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0e6),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 2.0e8),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ALPHA,
        ),
    ] {
        claims.insert_claim(claim(name, dims, value)).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal fiber".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn expansion(card: &MaterialCard) -> fs_material::state_point::ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(300.0),
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn norton(card: &MaterialCard) -> ExpandingCreepSource<NortonCreepSource<'_>> {
    let source = NortonCreepSource::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        || false,
    )
    .unwrap();
    ExpandingCreepSource::new(card, source, &expansion(card), || false).unwrap()
}
fn joint(card: &MaterialCard) -> ExpandingCreepSource<NortonJ2Source<'_>> {
    let source = NortonJ2Source::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        HardeningPropertyBasis::PlasticStrainModulus,
        || false,
    )
    .unwrap();
    ExpandingCreepSource::new(card, source, &expansion(card), || false).unwrap()
}
fn tol() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 1e-4,
        energy_j_m3: 1e-5,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn free_expansion_is_stress_free_and_does_not_age_norton_history() {
    let card = card();
    let source = norton(&card);
    let cold = source.at_uniaxial_temperature(300.0, 0.5).unwrap();
    let hot = source.at_uniaxial_temperature(400.0, 0.5).unwrap();
    let old = cold.initial_state();
    let free = ALPHA * 100.0;
    let response = hot.response(free, &old).unwrap();
    close(response.stress_pa, 0.0, 1e-8);
    close(response.state.mechanical.creep_strain(), 0.0, 1e-16);
    assert_eq!(response.state.temperature_k(), 400.0);
    for value in &response.strain[..3] {
        close(*value, free, 1e-16);
    }
    let checked = hot.advance_with_energy(0.0, &old, free, tol()).unwrap();
    assert_eq!(checked.state, response.state);
    close(checked.energy.physical_dissipation_j_m3, 0.0, 1e-12);
    close(checked.energy.reversible_heat_j_m3, 0.0, 1e-8);
    close(checked.energy.residual_j_m3, 0.0, tol().energy_j_m3);
}

#[test]
fn restrained_joint_thermal_fiber_yields_creeps_and_closes_reciprocal_work() {
    let card = card();
    let source = joint(&card);
    let cold = source.at_uniaxial_temperature(300.0, 0.5).unwrap();
    let hot = source.at_uniaxial_temperature(500.0, 0.5).unwrap();
    let old = cold.initial_state();
    let response = hot.response(0.0, &old).unwrap();
    assert!(response.stress_pa < -1.0e6);
    assert!(response.state.mechanical.creep_strain() < 0.0);
    assert!(response.state.mechanical.plastic_strain() < 0.0);
    assert!(response.state.mechanical.equivalent_creep_strain() > 0.0);
    assert!(response.state.mechanical.equivalent_plastic_strain() > 0.0);
    let checked = hot.advance_with_energy(0.0, &old, 0.0, tol()).unwrap();
    assert_eq!(checked.state, response.state);
    close(checked.stress_pa[0], response.stress_pa, 1e-5);
    assert!(checked.stress_pa[1..].iter().all(|s| s.abs() < 1e-5));
    assert!(checked.energy.creep_dissipation_j_m3 > 0.0);
    assert!(checked.energy.plastic_dissipation_j_m3 > 0.0);
    assert!(checked.energy.thermoelastic_storage_change_j_m3 != 0.0);
    assert!(checked.energy.reversible_heat_j_m3 != 0.0);
    close(checked.energy.residual_j_m3, 0.0, tol().energy_j_m3);
    assert_eq!(old.temperature_k(), 300.0);
}
