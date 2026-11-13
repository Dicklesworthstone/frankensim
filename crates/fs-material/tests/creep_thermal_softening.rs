//! Manufactured source curves test thermodynamic consistency, not calibrated data.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::source::*;
use fs_material::creep::thermal::expansion::ExpandingCreepSource;
use fs_material::creep::{CreepTolerances, GAS_CONSTANT_J_MOL_K};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_material::{SmallStrainLaw, Uniaxial, Voigt};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const E0: f64 = 1.2e9;
const E_SLOPE: f64 = -3.0e6;
const NU0: f64 = 0.2;
const NU_SLOPE: f64 = 5.0e-4;
const H: f64 = 5.0e7;
const RATE: f64 = 1.0e-4;
const REFERENCE_STRESS: f64 = 2.0e6;
const ACTIVATION: f64 = 25_000.0;

fn youngs(t: f64) -> f64 {
    E0 + E_SLOPE * (t - 300.0)
}
fn poisson(t: f64) -> f64 {
    NU0 + NU_SLOPE * (t - 300.0)
}
fn bulk(t: f64) -> f64 {
    youngs(t) / (3.0 * (1.0 - 2.0 * poisson(t)))
}
fn bulk_derivative(t: f64) -> f64 {
    let d = 1.0 - 2.0 * poisson(t);
    (E_SLOPE + 2.0 * youngs(t) * NU_SLOPE / d) / (3.0 * d)
}
fn shear(t: f64) -> f64 {
    youngs(t) / (2.0 * (1.0 + poisson(t)))
}
fn shear_derivative(t: f64) -> f64 {
    let d = 1.0 + poisson(t);
    (E_SLOPE - youngs(t) * NU_SLOPE / d) / (2.0 * d)
}
fn rate(t: f64) -> f64 {
    RATE * (ACTIVATION / GAS_CONSTANT_J_MOL_K * (1.0 / 300.0 - 1.0 / t)).exp()
}
fn dot(a: &Voigt, b: &Voigt) -> f64 {
    (0..6).map(|i| a[i] * b[i] * if i < 3 { 1.0 } else { 2.0 }).sum()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; tolerance={tolerance}"
    );
}
fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 2.0e-4,
        energy_j_m3: 2.0e-6,
    }
}
fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", temperature_quantity(), t).unwrap()
}
fn query<'a>(lower: &'a QueryPoint, upper: &'a QueryPoint) -> NortonSourceQuery<'a> {
    NortonSourceQuery {
        lower,
        upper,
        temperature_axis: "T",
        selection: MaterialPropertySelection::SingleClaimOnly,
        strain_limit: 0.05,
    }
}

// One optional interior perturbation checks complete-path admission; equal
// endpoints must not conceal changing kinetics/hardening or a nonsmooth potential.
fn card(alpha: f64, perturbation: Option<(&str, f64)>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, lower, upper) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, youngs(300.0), youngs(500.0)),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), poisson(300.0), poisson(500.0)),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, REFERENCE_STRESS, REFERENCE_STRESS),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, RATE, RATE),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 1.0, 1.0),
        (CREEP_ACTIVATION_ENERGY_PROPERTY, CREEP_ACTIVATION_ENERGY_DIMS, ACTIVATION, ACTIVATION),
        (CREEP_REFERENCE_TEMPERATURE_PROPERTY, Temperature::DIMS, 300.0, 300.0),
        (CREEP_HARDENING_REFERENCE_STRAIN_PROPERTY, Dims([0; 6]), 0.001, 0.001),
        (CREEP_HARDENING_EXPONENT_PROPERTY, Dims([0; 6]), 1.5, 1.5),
        (CREEP_RESIDUAL_RATE_FRACTION_PROPERTY, Dims([0; 6]), 0.2, 0.2),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 6.0e6, 1.0e6),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H, H),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, Dims([0, 0, 0, -1, 0, 0]), alpha, alpha),
    ] {
        let offset = perturbation.filter(|(property, _)| *property == name)
            .map_or(0.0, |(_, offset)| offset);
        let variable = lower != upper || offset != 0.0;
        let value = if variable {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, lower),
                    (400.0, 0.5 * (lower + upper) + offset),
                    (500.0, upper),
                ],
                dims,
            }
        } else {
            PropertyValue::Scalar { value: lower, dims }
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value,
            validity: ValidityDomain::unconstrained().with_quantity(
                "T", temperature_quantity(), 300.0, 500.0,
            ),
            interpolation: if variable {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured thermo-creep softening curves".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured elastic/creep reference".into(),
            phase: "solid".into(),
            process: "analytic test".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap()
}
fn expansion_reference(card: &MaterialCard) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(300.0),
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ScalarAdmissibility::Finite,
        ).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}
fn source(card: &MaterialCard) -> ExpandingCreepSource<NortonCreepSource<'_>> {
    ExpandingCreepSource::from_norton(
        card,
        query(&point(300.0), &point(500.0)),
        &expansion_reference(card),
        || false,
    ).unwrap()
}

#[test]
fn restrained_hydrostatic_heating_and_cooling_close_the_actual_internal_energy_cycle() {
    let card = card(0.0, None); // Explicitly sourced zero expansion isolates stiffness work.
    let source = source(&card);
    let strain = [0.002, 0.002, 0.002, 0.0, 0.0, 0.0];
    let trace: f64 = 0.006;
    let mut old_temperature = 300.0;
    let mut history = source.at_temperature(old_temperature, 1.0).unwrap().initial_state();
    let mut heat = 0.0;
    for temperature in [450.0, 320.0, 500.0, 300.0] {
        let law = source.at_temperature(temperature, 2.0).unwrap();
        let step = law.advance_with_energy(&strain, &history, &strain, tolerances()).unwrap();
        let initial_energy = 0.5 * bulk(old_temperature) * trace.powi(2);
        let final_energy = 0.5 * bulk(temperature) * trace.powi(2);
        let old_entropy = -0.5 * bulk_derivative(old_temperature) * trace.powi(2);
        let entropy = -0.5 * bulk_derivative(temperature) * trace.powi(2);
        let entropy_storage = temperature * entropy - old_temperature * old_entropy;
        for stress in &step.stress_pa[..3] {
            close(*stress, bulk(temperature) * trace, 2.0e-7);
        }
        close(step.energy.initial_stored_j_m3, initial_energy, 1.0e-7);
        close(step.energy.final_stored_j_m3, final_energy, 1.0e-7);
        close(step.energy.thermoelastic_storage_change_j_m3, entropy_storage, 1.0e-7);
        close(step.energy.reversible_heat_j_m3, initial_energy - final_energy - entropy_storage, 1.0e-7);
        close(step.energy.physical_dissipation_j_m3, 0.0, 1.0e-10);
        close(step.energy.numerical_dissipation_j_m3, 0.0, 1.0e-10);
        close(step.energy.mechanical_work_j_m3, 0.0, 1.0e-10);
        close(step.energy.residual_j_m3, 0.0, 1.0e-7);
        assert!(step.state.mechanical.equivalent_creep_strain < 1.0e-15);
        heat += step.energy.caloric_heat_j_m3();
        history = step.state;
        old_temperature = temperature;
    }
    close(heat, 0.0, 1.0e-7);
    assert!(source.minimum_thermal_capacity_j_m3_k().unwrap().is_finite());
}

#[test]
fn heated_nonvirgin_tensor_creep_matches_the_linear_norton_solution_and_tangent() {
    let card = card(0.0, None);
    let source = source(&card);
    let strain: Voigt = [0.004, -0.002, -0.002, 0.0006, 0.0003, -0.0002];
    let cold = source.at_temperature(300.0, 2.0).unwrap();
    let loaded = cold.advance_with_energy(
        &[0.0; 6], &cold.initial_state(), &strain, tolerances(),
    ).unwrap();
    assert!(loaded.state.mechanical.equivalent_creep_strain > 0.0);
    let accepted = loaded.state.clone();
    let temperature = 450.0;
    let dt = 0.5;
    let hot = source.at_temperature(temperature, dt).unwrap();
    let step = hot.advance_with_energy(&strain, &loaded.state, &strain, tolerances()).unwrap();
    let old_elastic: Voigt = core::array::from_fn(|i| {
        strain[i] - loaded.state.mechanical.creep_strain[i]
    });
    let ratio = (1.0 + 3.0 * shear(temperature) * dt * rate(temperature) / REFERENCE_STRESS).recip();
    let final_elastic: Voigt = core::array::from_fn(|i| ratio * old_elastic[i]);
    let creep_increment: Voigt = core::array::from_fn(|i| (1.0 - ratio) * old_elastic[i]);
    let expected_stress: Voigt = core::array::from_fn(|i| 2.0 * shear(temperature) * final_elastic[i]);
    for i in 0..6 {
        close(step.stress_pa[i], expected_stress[i], 1.0e-5);
        close(
            step.state.mechanical.creep_strain[i],
            loaded.state.mechanical.creep_strain[i] + creep_increment[i],
            1.0e-13,
        );
        for j in 0..6 {
            let expected = if i < 3 && j < 3 {
                bulk(temperature) + 2.0 * shear(temperature) * ratio
                    * (if i == j { 1.0 } else { 0.0 } - 1.0 / 3.0)
            } else if i == j {
                2.0 * shear(temperature) * ratio
            } else {
                0.0
            };
            close(step.tangent_pa[i][j], expected, 2.0e-3);
        }
    }
    let initial_energy = shear(300.0) * dot(&old_elastic, &old_elastic);
    let final_energy = shear(temperature) * dot(&final_elastic, &final_elastic);
    let entropy_storage = -temperature * shear_derivative(temperature) * dot(&final_elastic, &final_elastic)
        + 300.0 * shear_derivative(300.0) * dot(&old_elastic, &old_elastic);
    let stiffness_work = (shear(temperature) - shear(300.0)) * dot(&old_elastic, &old_elastic);
    close(step.energy.initial_stored_j_m3, initial_energy, 1.0e-6);
    close(step.energy.final_stored_j_m3, final_energy, 1.0e-6);
    close(step.energy.physical_dissipation_j_m3, dot(&expected_stress, &creep_increment), 1.0e-6);
    close(step.energy.thermoelastic_storage_change_j_m3, entropy_storage, 1.0e-6);
    close(step.energy.reversible_heat_j_m3, -stiffness_work - entropy_storage, 1.0e-6);
    close(step.energy.residual_j_m3, 0.0, 1.0e-6);
    assert!(step.state.mechanical.equivalent_creep_strain > accepted.mechanical.equivalent_creep_strain);
    assert_eq!(loaded.state, accepted);
}

#[test]
fn free_sided_continuation_reconstructs_old_poisson_ratio_and_old_stored_energy() {
    let card = card(0.0, None);
    let source = source(&card);
    let cold = source.at_uniaxial_temperature(300.0, 1.0).unwrap();
    let axial = 0.005;
    let loaded = cold.advance_with_energy(0.0, &cold.initial_state(), axial, tolerances()).unwrap();
    let hot = source.at_uniaxial_temperature(450.0, 0.5).unwrap();
    let reconstructed = hot.tensor_strain(axial, &loaded.state).unwrap();
    for (actual, original) in reconstructed.iter().zip(loaded.strain) {
        close(*actual, original, 1.0e-17);
    }
    let old_elastic = axial - loaded.state.mechanical.creep_strain();
    close(reconstructed[1], -poisson(300.0) * old_elastic - 0.5 * loaded.state.mechanical.creep_strain(), 1.0e-17);
    let step = hot.advance_with_energy(axial, &loaded.state, axial, tolerances()).unwrap();
    let expected = youngs(450.0) * old_elastic
        / (1.0 + youngs(450.0) * 0.5 * rate(450.0) / REFERENCE_STRESS);
    close(step.stress_pa[0], expected, 1.0e-5);
    for stress in &step.stress_pa[1..] {
        close(*stress, 0.0, 1.0e-5);
    }
    close(step.energy.initial_stored_j_m3, 0.5 * youngs(300.0) * old_elastic.powi(2), 1.0e-6);
    let final_elastic = axial - step.state.mechanical.creep_strain();
    close(step.energy.final_stored_j_m3, 0.5 * youngs(450.0) * final_elastic.powi(2), 1.0e-6);
    close(step.energy.residual_j_m3, 0.0, 1.0e-6);
    let property = hot.material().unwrap().property(YOUNG_MODULUS_PROPERTY).unwrap();
    close(property.value_si(), youngs(450.0), 1.0e-6);
    assert_eq!(
        property.answer().receipt.selected,
        source.mechanical_source().reference().property(YOUNG_MODULUS_PROPERTY).unwrap().answer().receipt.selected,
    );
}

#[test]
fn primary_creep_and_yield_softening_retain_both_histories_and_constant_hardening_storage() {
    let alpha = 2.0e-6;
    let card = card(alpha, None);
    let source = ExpandingCreepSource::from_norton_j2_with_strain_hardening(
        &card,
        query(&point(300.0), &point(500.0)),
        HardeningPropertyBasis::PlasticStrainModulus,
        &expansion_reference(&card),
        || false,
    ).unwrap();
    let cold = source.at_uniaxial_temperature(300.0, 1.0).unwrap();
    let loaded = cold.advance_with_energy(0.0, &cold.initial_state(), 0.008, tolerances()).unwrap();
    assert!(loaded.state.mechanical.equivalent_creep_strain() > 0.0);
    assert!(loaded.state.mechanical.equivalent_plastic_strain() > 0.0);
    let accepted = loaded.state.clone();
    let hot = source.at_uniaxial_temperature(450.0, 0.4).unwrap();
    let step = hot.advance_with_energy(0.008, &loaded.state, 0.010, tolerances()).unwrap();
    let old = &loaded.state.mechanical;
    let next = &step.state.mechanical;
    assert!(next.equivalent_creep_strain() > old.equivalent_creep_strain());
    assert!(next.equivalent_plastic_strain() > old.equivalent_plastic_strain());
    close(hot.hardening(), H, 0.0);
    close(hot.yield_stress(), 2.25e6, 1.0e-8);
    close(step.stress_pa[0], hot.yield_stress() + H * next.equivalent_plastic_strain(), 1.0e-4);
    let old_elastic = 0.008 - old.creep_strain() - old.plastic_strain();
    let final_elastic = 0.010 - alpha * 150.0 - next.creep_strain() - next.plastic_strain();
    close(
        step.energy.initial_stored_j_m3,
        0.5 * youngs(300.0) * old_elastic.powi(2) + 0.5 * H * old.equivalent_plastic_strain().powi(2),
        1.0e-6,
    );
    close(
        step.energy.final_stored_j_m3,
        0.5 * youngs(450.0) * final_elastic.powi(2) + 0.5 * H * next.equivalent_plastic_strain().powi(2),
        1.0e-6,
    );
    let h = 1.0e-7;
    let derivative = (hot.stress(0.010 + h, &loaded.state) - hot.stress(0.010 - h, &loaded.state)) / (2.0 * h);
    let tangent = hot.tangent(0.010, &loaded.state);
    close(tangent, derivative, 1.0e-5 * tangent.abs().max(1.0));
    assert!(step.energy.physical_dissipation_j_m3 > 0.0);
    assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
    close(step.energy.residual_j_m3, 0.0, 1.0e-6);
    assert_eq!(loaded.state, accepted);
}

#[test]
fn varying_elasticity_requires_its_thermal_owner_and_full_smooth_source_support() {
    let lower = point(300.0);
    let upper = point(500.0);
    let card = card(0.0, None);
    assert!(NortonCreepSource::new(&card, query(&lower, &upper), || false).is_err());
    let bound = source(&card);
    assert!(bound.mechanical_source().at_temperature(400.0, 1.0).is_err());
    assert!(bound.at_temperature(299.0, 1.0).is_err());
    assert!(bound.at_temperature(400.0, 0.0).is_err());
    assert!(ExpandingCreepSource::from_norton_j2(
        &card, query(&lower, &upper), HardeningPropertyBasis::UniaxialTangentModulus,
        &expansion_reference(&card), || false,
    ).is_err());
    for (property, offset) in [
        (YOUNG_MODULUS_PROPERTY, 1.0e8),
        (POISSON_RATIO_PROPERTY, 0.01),
        (CREEP_REFERENCE_RATE_PROPERTY, 1.0e-5),
        (CREEP_HARDENING_EXPONENT_PROPERTY, 0.2),
    ] {
        let altered = self::card(0.0, Some((property, offset)));
        assert!(ExpandingCreepSource::from_norton_with_strain_hardening(
            &altered, query(&lower, &upper), &expansion_reference(&altered), || false,
        ).is_err(), "unexpected admission of {property}");
    }
    // E and nu must stay affine, but a continuous piecewise-affine expansion coefficient is
    // admitted (tabulated CTE through the shared thermoelastic potential, since 273f716cd3).
    let kinked_alpha = self::card(
        0.0,
        Some((LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, 1.0e-6)),
    );
    assert!(
        ExpandingCreepSource::from_norton_with_strain_hardening(
            &kinked_alpha,
            query(&lower, &upper),
            &expansion_reference(&kinked_alpha),
            || false,
        )
        .is_ok()
    );
    let altered = self::card(0.0, Some((ISOTROPIC_HARDENING_MODULUS_PROPERTY, 1.0e6)));
    assert!(ExpandingCreepSource::from_norton_j2(
        &altered, query(&lower, &upper), HardeningPropertyBasis::PlasticStrainModulus,
        &expansion_reference(&altered), || false,
    ).is_err());
}
