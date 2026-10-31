//! Manufactured thermo-creep references, not measured material qualification.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::CreepTolerances;
use fs_material::creep::plastic::source::NortonJ2Source;
use fs_material::creep::source::*;
use fs_material::creep::thermal::expansion::ExpandingCreepSource;
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_material::{SmallStrainLaw, Voigt};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const H: f64 = 4.0e8;
const ALPHA: f64 = 1.0e-5;
const SLOPE: f64 = 1.0e-8;
const LIMIT: f64 = 0.02;
fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), t)
        .unwrap()
}
fn source_card(slope: f64, nonaffine: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), NU),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, 2e6),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 1e-5),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 1.0),
        (
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            0.0,
        ),
        (
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            300.0,
        ),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 4e6),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ALPHA,
        ),
    ] {
        let knots = if name == YIELD_STRESS_PROPERTY {
            Some(vec![(300.0, 4e6), (500.0, 1e6)])
        } else if name == LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY {
            Some(vec![
                (300.0, ALPHA),
                (
                    400.0,
                    ALPHA + 100.0 * slope + if nonaffine { 1e-6 } else { 0.0 },
                ),
                (500.0, ALPHA + 200.0 * slope),
            ])
        } else {
            None
        };
        let interpolation = if knots.is_some() {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::ConstantWithinValidity
        };
        let value = match knots {
            Some(knots) => PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims,
            },
            None => PropertyValue::Scalar { value, dims },
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    absolute_temperature(),
                    300.0,
                    500.0,
                ),
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured affine thermo-creep".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal solid".into(),
            phase: "solid".into(),
            process: "manufactured reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn reference(card: &MaterialCard) -> ResolvedMaterialStatePoint {
    let alpha = ScalarPropertyRequirement::try_new(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        Dims([0, 0, 0, -1, 0, 0]),
        ScalarAdmissibility::Finite,
    )
    .unwrap();
    resolve_material_state_point(
        card,
        &point(300.0),
        &[alpha],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn mechanical(card: &MaterialCard) -> NortonCreepSource<'_> {
    NortonCreepSource::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: LIMIT,
        },
        || false,
    )
    .unwrap()
}
fn source(card: &MaterialCard) -> ExpandingCreepSource<NortonCreepSource<'_>> {
    ExpandingCreepSource::new(card, mechanical(card), &reference(card), || false).unwrap()
}
fn joint_source(card: &MaterialCard) -> ExpandingCreepSource<NortonJ2Source<'_>> {
    let joint = NortonJ2Source::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: LIMIT,
        },
        HardeningPropertyBasis::PlasticStrainModulus,
        || false,
    )
    .unwrap();
    ExpandingCreepSource::new(card, joint, &reference(card), || false).unwrap()
}
fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 1e-4,
        energy_j_m3: 1e-6,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance={tolerance}"
    );
}
fn free(t: f64) -> f64 {
    ALPHA * (t - 300.0) + 0.5 * SLOPE * (t - 300.0).powi(2)
}
fn alpha(t: f64) -> f64 {
    ALPHA + SLOPE * (t - 300.0)
}
fn bulk() -> f64 {
    E / (3.0 * (1.0 - 2.0 * NU))
}

#[test]
fn free_and_fully_restrained_affine_expansion_have_independent_entropy_and_work() {
    let card = source_card(SLOPE, false);
    let source = source(&card);
    let cold = source.at_temperature(300.0, 1.0).unwrap();
    let old = cold.initial_state();
    let t = 450.0;
    let hot = source.at_temperature(t, 1.0).unwrap();
    let f = free(t);
    let expanded = [f, f, f, 0.0, 0.0, 0.0];
    let next = hot
        .advance_with_energy(&[0.0; 6], &old, &expanded, tolerances())
        .unwrap();
    assert_eq!(next.strain, expanded);
    for s in next.stress_pa {
        close(s, 0.0, 2e-8);
    }
    close(next.energy.physical_dissipation_j_m3, 0.0, 1e-14);
    close(next.energy.reversible_heat_j_m3, 0.0, 1e-7);
    let fixed = hot
        .advance_with_energy(&[0.0; 6], &old, &[0.0; 6], tolerances())
        .unwrap();
    let stress = -3.0 * bulk() * f;
    for s in &fixed.stress_pa[..3] {
        close(*s, stress, 1e-7);
    }
    assert_eq!(fixed.state.mechanical, old.mechanical);
    let entropy = 3.0 * alpha(t) * stress;
    let stored = 4.5 * bulk() * f * f;
    close(fixed.energy.final_stored_j_m3, stored, 1e-7);
    close(fixed.energy.numerical_dissipation_j_m3, stored, 1e-7);
    close(
        fixed.energy.thermoelastic_storage_change_j_m3,
        t * entropy,
        1e-6,
    );
    close(
        fixed.energy.reversible_heat_j_m3,
        3.0 * stress * f - t * entropy,
        1e-6,
    );
    close(fixed.energy.mechanical_work_j_m3, 0.0, 0.0);
    assert_eq!(fixed.energy.physical_dissipation_j_m3, 0.0);
    assert_eq!(fixed.state.temperature_k(), t);
    assert_eq!(old.temperature_k(), 300.0);
}

#[test]
fn mixed_shear_uses_tensor_tangent_and_mandel_thermodynamics_consistently() {
    let card = source_card(SLOPE, false);
    let source = source(&card);
    let t = 450.0;
    let dt = 0.4;
    let at = source.at_temperature(t, dt).unwrap();
    let old = source.at_temperature(300.0, dt).unwrap().initial_state();
    let mechanical = [0.002, -0.001, 0.0004, 0.0007, -0.0005, 0.0009];
    let total: Voigt = core::array::from_fn(|i| mechanical[i] + if i < 3 { free(t) } else { 0.0 });
    let next = at
        .advance_with_energy(&[0.0; 6], &old, &total, tolerances())
        .unwrap();
    let g = E / (2.0 * (1.0 + NU));
    let ratio = 1.0 / (1.0 + 3.0 * g * dt * 1e-5 / 2e6);
    let trace = mechanical[..3].iter().sum::<f64>();
    for i in 0..6 {
        let expected = if i < 3 {
            bulk() * trace + 2.0 * g * ratio * (mechanical[i] - trace / 3.0)
        } else {
            2.0 * g * ratio * mechanical[i]
        };
        close(next.stress_pa[i], expected, 1e-7);
        for j in 0..6 {
            let volume = if i < 3 && j < 3 { 1.0 } else { 0.0 };
            let identity = if i == j { 1.0 } else { 0.0 };
            close(
                next.tangent_pa[i][j],
                bulk() * volume + 2.0 * g * ratio * (identity - volume / 3.0),
                2e-6,
            );
        }
    }
    let mut recoverable: Voigt =
        core::array::from_fn(|i| total[i] - next.state.mechanical.creep_strain[i]);
    for v in &mut recoverable[3..] {
        *v *= 2.0_f64.sqrt();
    }
    let potential = at.thermoelastic().state(recoverable, 1.0).unwrap();
    close(potential.free_energy_j, next.energy.final_stored_j_m3, 1e-7);
    close(
        potential.entropy_j_k,
        alpha(t) * next.stress_pa[..3].iter().sum::<f64>(),
        1e-9,
    );
    close(
        next.energy.thermoelastic_storage_change_j_m3,
        t * potential.entropy_j_k,
        1e-6,
    );
}

#[test]
fn axial_thermal_restraint_drives_joint_creep_yield_and_preserves_cooling_history() {
    let card = source_card(SLOPE, false);
    let source = joint_source(&card);
    let cold = source.at_temperature(300.0, 1.0).unwrap();
    let hot = source.at_temperature(500.0, 1.0).unwrap();
    let old = cold.initial_state();
    let f = free(500.0);
    let a = 1e-5 / 2e6;
    let yield_hot = 1e6;
    let q = (H * E * f + E * yield_hot) / (E + H + E * H * a);
    let dc = a * q;
    let dp = (q - yield_hot) / H;
    let lateral = f + NU * q / E + 0.5 * (dc + dp);
    let total = [0.0, lateral, lateral, 0.0, 0.0, 0.0];
    let next = hot
        .advance_with_energy(&[0.0; 6], &old, &total, tolerances())
        .unwrap();
    close(next.stress_pa[0], -q, 1e-7);
    for s in &next.stress_pa[1..] {
        close(*s, 0.0, 1e-7);
    }
    close(next.state.mechanical.creep.creep_strain[0], -dc, 2e-17);
    close(next.state.mechanical.plastic.plastic_strain[0], -dp, 2e-17);
    close(next.state.mechanical.plastic.alpha, dp, 2e-17);
    close(next.energy.creep_dissipation_j_m3, q * dc, 1e-7);
    close(next.energy.plastic_dissipation_j_m3, yield_hot * dp, 1e-7);
    let stored = q * q / (2.0 * E) + 0.5 * H * dp * dp;
    close(next.energy.final_stored_j_m3, stored, 1e-7);
    close(next.energy.numerical_dissipation_j_m3, stored, 1e-7);
    let ts = -500.0 * alpha(500.0) * q;
    close(next.energy.thermoelastic_storage_change_j_m3, ts, 1e-6);
    close(next.energy.reversible_heat_j_m3, -q * f - ts, 1e-6);
    close(next.energy.mechanical_work_j_m3, 0.0, 1e-8);
    let before = next.state.clone();
    let cooled = cold
        .advance_with_energy(&total, &before, &[0.0; 6], tolerances())
        .unwrap();
    assert_eq!(cooled.state.temperature_k(), 300.0);
    assert_eq!(before.temperature_k(), 500.0);
    close(
        cooled.state.mechanical.plastic.alpha,
        before.mechanical.plastic.alpha,
        2e-17,
    );
    assert!(
        cooled.state.mechanical.creep.equivalent_creep_strain
            > before.mechanical.creep.equivalent_creep_strain
    );
    assert_eq!(next.state, before);
}

#[test]
fn capacity_floor_covers_separate_histories_and_original_source_admission() {
    let card = source_card(SLOPE, false);
    let source = source(&card);
    let bound = source.minimum_thermal_capacity_j_m3_k().unwrap();
    let expected = -500.0
        * 9.0
        * bulk()
        * (alpha(500.0).powi(2) + SLOPE * (2.0 * LIMIT + 200.0 * alpha(500.0)));
    assert!(bound <= expected);
    close(bound, expected, 1e-8);
    assert!(
        joint_source(&card)
            .minimum_thermal_capacity_j_m3_k()
            .unwrap()
            < bound
    );
    for t in [300.0, 400.0, 500.0] {
        for sign in [-1.0, 1.0] {
            let strain = [
                sign * 0.015,
                sign * 0.015,
                sign * 0.015,
                0.002,
                -0.003,
                0.004,
            ];
            let state = source
                .potential()
                .at_temperature(t)
                .unwrap()
                .state(strain, 1.0)
                .unwrap();
            let trace = 3.0 * sign * 0.015 - 3.0 * free(t);
            let capacity = t * (3.0 * bulk() * SLOPE * trace - 9.0 * bulk() * alpha(t).powi(2));
            close(state.heat_capacity_j_k, capacity, 1e-8);
            assert!(capacity >= bound);
        }
    }
    let foreign = source_card(SLOPE * 2.0, false);
    assert!(
        ExpandingCreepSource::new(&card, mechanical(&card), &reference(&foreign), || false)
            .is_err()
    );
    let missing = mechanical(&card);
    assert!(
        ExpandingCreepSource::new(&card, missing.clone(), missing.reference(), || false).is_err()
    );
    assert!(matches!(
        ExpandingCreepSource::new(&card, mechanical(&card), &reference(&card), || true),
        Err(CreepSourceError::Cancelled)
    ));
    let kinked = source_card(SLOPE, true);
    assert!(
        ExpandingCreepSource::new(&kinked, mechanical(&kinked), &reference(&kinked), || false)
            .is_err()
    );
}

#[test]
fn frozen_trials_temperature_reversal_and_refusals_do_not_double_advance() {
    let card = source_card(SLOPE, false);
    let source = source(&card);
    assert!(source.admit_stress_free_temperature(300.0).is_ok());
    assert!(source.admit_stress_free_temperature(400.0).is_err());
    assert!(source.at_temperature(500.01, 1.0).is_err());
    let cold = source.at_temperature(300.0, 0.2).unwrap();
    let hot = source.at_temperature(450.0, 0.2).unwrap();
    let old = cold.initial_state();
    let total = [0.0014, -0.0004, 0.0001, 0.0005, 0.0, 0.0];
    let next = hot
        .advance_with_energy(&[0.0; 6], &old, &total, tolerances())
        .unwrap();
    assert_eq!(
        next,
        hot.advance_with_energy(&[0.0; 6], &old, &total, tolerances())
            .unwrap()
    );
    assert_eq!(
        hot.update_state(&total, &old),
        hot.update_state(&total, &old)
    );
    assert_eq!(next.strain, total);
    assert_eq!(old, cold.initial_state());
    assert!(
        hot.advance_with_energy(&[0.0; 6], &old, &[f64::NAN; 6], tolerances())
            .is_err()
    );
    assert!(hot.stress(&[f64::NAN; 6], &old).iter().all(|s| s.is_nan()));
    let fixed = hot
        .advance_with_energy(&[0.0; 6], &old, &[0.0; 6], tolerances())
        .unwrap();
    let returned = cold
        .advance_with_energy(&[0.0; 6], &fixed.state, &[0.0; 6], tolerances())
        .unwrap();
    assert_eq!(returned.state.mechanical, old.mechanical);
    assert_eq!(returned.state.temperature_k(), 300.0);
    assert!(returned.energy.reversible_heat_j_m3 < 0.0);
    assert_eq!(returned.energy.physical_dissipation_j_m3, 0.0);
    assert!(returned.energy.residual_j_m3.abs() < tolerances().energy_j_m3);
}
