//! G0/G1/G3: source data reach the original brittle/cohesive owner without defaults.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::brittle::source::*;
use fs_material::brittle::{BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings};
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::elastic::IsotropicElastic;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure, Temperature};

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured brittle law; not measured glass or ceramic data".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn properties() -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1e9),
        scalar(POISSON_RATIO_PROPERTY, Dims::NONE, 0.2),
        scalar(BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, 1e6),
        scalar(BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 0.8e6),
        scalar(BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY, FRACTURE_ENERGY_PER_AREA_DIMS, 100.0),
        scalar(BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY, FRACTURE_ENERGY_PER_AREA_DIMS, 200.0),
        scalar(BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5),
    ]
}

fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle coupon".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 1,
        },
        claims,
        vec![],
    ).unwrap()
}

fn settings(width: f64) -> BrittleCrackBandSettings {
    BrittleCrackBandSettings {
        band_width_m: width,
        normal: [1.0, 0.0, 0.0],
        maximum_strain: 0.05,
        energy_tolerance_j_m3: 1e-6,
    }
}

fn resolve(source: &MaterialCard, width: f64) -> ResolvedBrittleCrackBand {
    resolve_brittle_crack_band(
        source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        settings(width),
    ).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!((actual-expected).abs() <= 2e-10*expected.abs().max(1.0),
        "{actual:e} != {expected:e}");
}

#[test]
fn original_source_coefficients_drive_native_damage_tangent_and_work() {
    let source = card(properties());
    let resolved = resolve(&source, 0.01);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    assert_eq!(resolved.material().properties().len(), 7);
    let native = BrittleCrackBand::new(
        IsotropicElastic::new(1e9, 0.2, 0.05).unwrap(),
        BrittleCrackBandParameters {
            tensile_strength_pa: 1e6,
            shear_strength_pa: 0.8e6,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 200.0,
            bk_exponent: 1.5,
        },
        settings(0.01),
    ).unwrap();
    let mut sourced = resolved.law().initial_state();
    let mut direct = native.initial_state();
    for strain in [
        [0.0003, 0.0, 0.0, 0.0001, 0.0, 0.0],
        [0.004, 0.0, 0.0, 0.001, 0.0, 0.0],
        [0.001, 0.0, 0.0, 0.0001, 0.0, 0.0],
        [-0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
    ] {
        let a = resolved.law().trial(&strain, &sourced).unwrap();
        let b = native.trial(&strain, &direct).unwrap();
        assert_eq!(a.stress_pa, b.stress_pa);
        assert_eq!(a.tangent_pa, b.tangent_pa);
        assert_eq!(a.branch_tangent_pa, b.branch_tangent_pa);
        near(a.state.damage(), b.state.damage());
        near(a.energy.work, b.energy.work);
        near(a.energy.dissipated, b.energy.dissipated);
        let solver_law = resolved.small_strain_law();
        assert_eq!(solver_law.stress(&strain, &sourced), a.stress_pa);
        assert_eq!(solver_law.tangent(&strain, &sourced), a.branch_tangent_pa);
        sourced = a.state;
        direct = b.state;
    }
    assert!(sourced.damage() > 0.0);
}

#[test]
fn original_multiphysics_bundle_retains_density_and_each_source_receipt() {
    let mut claims = properties();
    claims.push(scalar(DENSITY_PROPERTY, Density::DIMS, 2300.0));
    let source = card(claims);
    let mut requirements = brittle_crack_band_requirements().unwrap();
    requirements.push(ScalarPropertyRequirement::try_new(
        DENSITY_PROPERTY,
        Density::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    ).unwrap());
    let bundle = resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let original = bundle.clone();
    let resolved = ResolvedBrittleCrackBand::try_from_material(bundle, settings(0.01)).unwrap();
    assert_eq!(resolved.material(), &original);
    near(resolved.material().property(DENSITY_PROPERTY).unwrap().value_si(), 2300.0);
}

#[test]
fn every_missing_coefficient_and_wrong_fracture_units_refuse() {
    for omitted in 0..7 {
        let mut claims = properties();
        claims.remove(omitted);
        assert!(resolve_brittle_crack_band(
            &card(claims),
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        ).is_err());
    }
    // Bulk energy density, total work and dimensionless data do not become Gc.
    for wrong in [Pressure::DIMS, Dims([2, 1, -2, 0, 0, 0]), Dims::NONE] {
        let mut claims = properties();
        claims[4] = scalar(BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY, wrong, 100.0);
        assert!(resolve_brittle_crack_band(
            &card(claims),
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        ).is_err());
    }
    let mut claims = properties();
    claims[2] = scalar("modulus_of_rupture", Pressure::DIMS, 1e6);
    assert!(resolve_brittle_crack_band(
        &card(claims),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
    ).is_err());
}

#[test]
fn fixed_temperature_changes_fracture_response_without_resetting_original_history() {
    let mut claims = properties();
    claims[4].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 100.0), (400.0, 200.0)],
        dims: FRACTURE_ENERGY_PER_AREA_DIMS,
    };
    claims[4].validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
    claims[4].interpolation = InterpolationPolicy::LinearInside;
    let source = card(claims);
    let at = |temperature| resolve_brittle_crack_band(
        &source,
        &QueryPoint::new().with("T", temperature, Temperature::DIMS),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
    );
    let cold = at(300.0).unwrap();
    let hot = at(400.0).unwrap();
    let strain = [0.025, 0.0, 0.0, 0.0, 0.0, 0.0];
    let failed = cold.law().trial(&strain, &cold.law().initial_state()).unwrap();
    let surviving = hot.law().trial(&strain, &hot.law().initial_state()).unwrap();
    assert!(failed.state.failed());
    assert!(!surviving.state.failed());
    assert!(surviving.stress_pa[0] > 0.0);
    let saved = failed.state;
    assert!(hot.law().trial(&strain, &saved).is_err());
    assert!(saved.failed());
    assert!(at(401.0).is_err());
    assert_ne!(cold.material().identity(), hot.material().identity());
}

#[test]
fn mechanical_and_rate_applicability_cannot_be_frozen_as_brittle_constants() {
    for (axis, dims, value) in [
        ("strain", Dims::NONE, 0.001),
        ("strain_rate", Dims([0, 0, -1, 0, 0, 0]), 1.0),
        ("frequency", Dims([0, 0, -1, 0, 0, 0]), 1.0),
    ] {
        let mut claims = properties();
        claims[2].validity = ValidityDomain::unconstrained().with(axis, 0.0, 2.0);
        assert!(resolve_brittle_crack_band(
            &card(claims),
            &QueryPoint::new().with(axis, value, dims),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        ).is_err());
    }
}

#[test]
fn geometric_band_width_preserves_original_source_fracture_energy() {
    let source = card(properties());
    for width in [0.01, 0.02, 0.04] {
        let resolved = resolve(&source, width);
        let strain = [0.025, 0.0, 0.0, 0.0, 0.0, 0.0];
        let step = resolved.law().trial(&strain, &resolved.law().initial_state()).unwrap();
        assert!(step.state.failed());
        near(step.state.dissipated_j_m3()*width, 100.0);
        near(step.energy.work, step.energy.stored_change+step.energy.dissipated);
        assert_eq!(resolved.material().card_identity(), source.content_hash());
    }
}
