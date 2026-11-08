//! Manufactured thermal fibers: permanent creep, yield and recoverable resistance.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::recovery::plastic::source::RecoveringNortonJ2Source;
use fs_material::creep::recovery::plastic::uniaxial::UniaxialRecoveringNortonJ2State;
use fs_material::creep::recovery::source::{
    CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY, CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
    CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY, RecoveringNortonSource,
};
use fs_material::creep::recovery::uniaxial::UniaxialRecoveringNortonState;
use fs_material::creep::source::*;
use fs_material::creep::thermal::expansion::{ExpandingCreepSource, ThermalCreepState};
use fs_material::creep::{CreepTolerances, GAS_CONSTANT_J_MOL_K};
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_material::{SmallStrainLaw, Uniaxial};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const E: f64 = 1e9;
const NU: f64 = 0.25;
const EXPANSION: f64 = 2e-5;
const STRESS: f64 = 2e6;
const RATE: f64 = 1e-4;
const PRIMARY: f64 = 2e-4;
const H: f64 = 2e8;

fn soft_youngs(t: f64) -> f64 {
    E - 2e6 * (t - 300.0)
}
fn soft_poisson(t: f64) -> f64 {
    NU + 5e-4 * (t - 300.0)
}
fn soft_alpha(t: f64) -> f64 {
    EXPANSION + 5e-8 * (t - 300.0)
}
fn soft_free(t: f64) -> f64 {
    let delta = t - 300.0;
    EXPANSION * delta + 0.5 * 5e-8 * delta * delta
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn yield_stress(t: f64) -> f64 {
    3e6 - 1e4 * (t - 300.0)
}
fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 500.0),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured recovering thermal fiber".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card() -> MaterialCard {
    source_card(false, 0.8, None)
}
fn source_card(
    thermal_elasticity: bool,
    recovery_rate: f64,
    perturbation: Option<(&str, f64)>,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims::NONE, NU),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, STRESS),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, RATE),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims::NONE, 1.0),
        (
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            20_000.0,
        ),
        (
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            400.0,
        ),
        (
            CREEP_HARDENING_REFERENCE_STRAIN_PROPERTY,
            Dims::NONE,
            PRIMARY,
        ),
        (CREEP_HARDENING_EXPONENT_PROPERTY, Dims::NONE, 1.0),
        (CREEP_RESIDUAL_RATE_FRACTION_PROPERTY, Dims::NONE, 0.0),
        (
            CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
            CREEP_RATE_DIMS,
            recovery_rate,
        ),
        (
            CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            15_000.0,
        ),
        (
            CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            400.0,
        ),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            EXPANSION,
        ),
    ] {
        let mut property = claim(name, dims, value);
        let (lower, upper) = if thermal_elasticity {
            match name {
                YOUNG_MODULUS_PROPERTY => (soft_youngs(300.0), soft_youngs(500.0)),
                POISSON_RATIO_PROPERTY => (soft_poisson(300.0), soft_poisson(500.0)),
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY => {
                    (soft_alpha(300.0), soft_alpha(500.0))
                }
                _ => (value, value),
            }
        } else {
            (value, value)
        };
        let offset = perturbation
            .filter(|(property, _)| *property == name)
            .map_or(0.0, |(_, offset)| offset);
        if lower != upper || offset != 0.0 {
            property.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, lower),
                    (400.0, 0.5 * (lower + upper) + offset),
                    (500.0, upper),
                ],
                dims,
            };
            property.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(property).unwrap();
    }
    let mut yield_claim = claim(YIELD_STRESS_PROPERTY, Pressure::DIMS, yield_stress(300.0));
    yield_claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, yield_stress(300.0)), (500.0, yield_stress(500.0))],
        dims: Pressure::DIMS,
    };
    yield_claim.interpolation = InterpolationPolicy::LinearInside;
    claims.insert_claim(yield_claim).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic recovering thermal fiber".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn expansion(card: &MaterialCard) -> ResolvedMaterialStatePoint {
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
fn pure_source(card: &MaterialCard) -> ExpandingCreepSource<RecoveringNortonSource<'_>> {
    let mechanical = RecoveringNortonSource::new(
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
    ExpandingCreepSource::new(card, mechanical, &expansion(card), || false).unwrap()
}
fn joint_source(card: &MaterialCard) -> ExpandingCreepSource<RecoveringNortonJ2Source<'_>> {
    let mechanical = RecoveringNortonJ2Source::new(
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
    ExpandingCreepSource::new(card, mechanical, &expansion(card), || false).unwrap()
}
fn soft_pure_source(card: &MaterialCard) -> ExpandingCreepSource<RecoveringNortonSource<'_>> {
    ExpandingCreepSource::from_recovering_norton(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        &expansion(card),
        || false,
    )
    .unwrap()
}
fn soft_joint_source(card: &MaterialCard) -> ExpandingCreepSource<RecoveringNortonJ2Source<'_>> {
    ExpandingCreepSource::from_recovering_norton_j2(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        HardeningPropertyBasis::PlasticStrainModulus,
        &expansion(card),
        || false,
    )
    .unwrap()
}
fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 1e-4,
        energy_j_m3: 1e-5,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; allowance={tolerance:e}"
    );
}

#[derive(Clone, Copy, Default)]
struct History {
    creep: f64,
    creep_path: f64,
    plastic: f64,
    plastic_path: f64,
    resistance: f64,
}
fn pure_history(state: &UniaxialRecoveringNortonState) -> History {
    History {
        creep: state.creep_strain(),
        creep_path: state.equivalent_creep_strain(),
        resistance: state.resistance_strain(),
        ..History::default()
    }
}
fn joint_history(state: &UniaxialRecoveringNortonJ2State) -> History {
    History {
        creep: state.creep_strain(),
        creep_path: state.equivalent_creep_strain(),
        plastic: state.plastic_strain(),
        plastic_path: state.equivalent_plastic_strain(),
        resistance: state.resistance_strain(),
    }
}

// Independent endpoint: eliminate primary flow by its quadratic at a fixed q,
// then solve q + E*(dc+dp) = |trial|. No production constitutive owner is used.
fn reference(joint: bool, t: f64, dt: f64, strain: f64, old: History) -> (History, f64, f64) {
    reference_at_elasticity(joint, t, dt, strain, old, E, EXPANSION * (t - 300.0))
}
fn reference_at_elasticity(
    joint: bool,
    t: f64,
    dt: f64,
    strain: f64,
    old: History,
    youngs: f64,
    free: f64,
) -> (History, f64, f64) {
    let rate = RATE * (20_000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / 400.0 - 1.0 / t)).exp();
    let recovery = 0.8 * (15_000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / 400.0 - 1.0 / t)).exp();
    let g = 1.0 + recovery * dt;
    let a = PRIMARY * g;
    let coefficient = dt * rate * a / STRESS;
    let trial = youngs * (strain - free - old.creep - old.plastic);
    let increments = |q: f64| {
        let b = a + old.resistance;
        let c = coefficient * q;
        let dc = 2.0 * c / (b + (b * b + 4.0 * c).sqrt());
        let dp = if joint {
            ((q - yield_stress(t)) / H - old.plastic_path).max(0.0)
        } else {
            0.0
        };
        (dc, dp)
    };
    let (mut lo, mut hi) = (0.0, trial.abs());
    for _ in 0..100 {
        let q = 0.5 * (lo + hi);
        let (dc, dp) = increments(q);
        if q + youngs * (dc + dp) > trial.abs() {
            hi = q;
        } else {
            lo = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let (dc, dp) = increments(q);
    let tangent = youngs
        / (1.0
            + youngs
                * (coefficient / (a + old.resistance + 2.0 * dc)
                    + if dp > 0.0 { 1.0 / H } else { 0.0 }));
    (
        History {
            creep: old.creep + trial.signum() * dc,
            creep_path: old.creep_path + dc,
            plastic: old.plastic + trial.signum() * dp,
            plastic_path: old.plastic_path + dp,
            resistance: (old.resistance + dc) / g,
        },
        trial.signum() * q,
        tangent,
    )
}

macro_rules! trajectory {
    ($name:ident, $source:ident, $history:ident, $joint:expr, $claims:expr, $soft:expr) => {
        #[test]
        fn $name() {
            let card = if $soft {
                source_card(true, 0.8, None)
            } else {
                card()
            };
            let source = $source(&card);
            let youngs = |t| if $soft { soft_youngs(t) } else { E };
            let poisson = |t| if $soft { soft_poisson(t) } else { NU };
            let free = |t| {
                if $soft {
                    soft_free(t)
                } else {
                    EXPANSION * (t - 300.0)
                }
            };
            let alpha = |t| if $soft { soft_alpha(t) } else { EXPANSION };
            let youngs_slope = if $soft { -2e6 } else { 0.0 };
            let cold = source.at_uniaxial_temperature(300.0, 0.25).unwrap();
            let mut old = cold.initial_state();
            let mut expected = History::default();
            let (mut old_t, mut old_strain, mut old_stress) = (300.0, 0.0, 0.0);
            let mut old_tensor = [0.0; 6];
            for (step, (t, dt)) in [
                (400.0, 0.25),
                (500.0, 0.25),
                (500.0, 0.5),
                (if $soft { 450.0 } else { 500.0 }, 2.0),
                (500.0, 0.25),
            ]
            .into_iter()
            .enumerate()
            {
                let law = source.at_uniaxial_temperature(t, dt).unwrap();
                let strain = match step {
                    2 | 3 => free(t) + expected.creep + expected.plastic,
                    4 => -0.001,
                    _ => 0.0,
                };
                let (next, stress, tangent) = if $soft {
                    reference_at_elasticity($joint, t, dt, strain, expected, youngs(t), free(t))
                } else {
                    reference($joint, t, dt, strain, expected)
                };
                let frozen = old.clone();
                // A hot law must not reconstruct old history with its new nu.
                for (actual, wanted) in law
                    .tensor_strain(old_strain, &old)
                    .unwrap()
                    .iter()
                    .zip(old_tensor)
                {
                    close(*actual, wanted, 1e-16);
                }
                let response = law.response(strain, &old).unwrap();
                let checked = law
                    .advance_with_energy(old_strain, &old, strain, tolerances())
                    .unwrap();
                assert_eq!(old, frozen);
                assert_eq!(checked.state, response.state);
                assert_eq!(checked.state.temperature_k(), t);
                close(response.stress_pa, stress, 2e-6);
                close(response.tangent_pa, tangent, 3e-5);
                let actual = $history(&checked.state.mechanical);
                for (a, b) in [
                    (actual.creep, next.creep),
                    (actual.creep_path, next.creep_path),
                    (actual.plastic, next.plastic),
                    (actual.plastic_path, next.plastic_path),
                    (actual.resistance, next.resistance),
                ] {
                    close(a, b, 5e-14);
                }
                close(response.strain[0], strain, 1e-16);
                close(
                    response.strain[1],
                    free(t) - poisson(t) * stress / youngs(t) - 0.5 * (next.creep + next.plastic),
                    5e-14,
                );
                assert_eq!(response.strain[1], response.strain[2]);
                assert!(checked.stress_pa[1..].iter().all(|s| s.abs() < 2e-6));

                let dc = next.creep_path - expected.creep_path;
                let dp = next.plastic_path - expected.plastic_path;
                let energy = checked.energy;
                let old_elastic = old_stress / youngs(old_t);
                let elastic = stress / youngs(t);
                let old_entropy =
                    alpha(old_t) * old_stress - 0.5 * youngs_slope * old_elastic * old_elastic;
                let entropy = alpha(t) * stress - 0.5 * youngs_slope * elastic * elastic;
                let entropy_storage = t * entropy - old_t * old_entropy;
                let nu = poisson(t);
                let lambda = youngs(t) * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
                let mu = youngs(t) / (2.0 * (1.0 + nu));
                // At frozen old tensor strain, current elasticity need not be
                // free-sided. Independently evaluate its Lamé quadratic form.
                let old_trace = (1.0 - 2.0 * poisson(old_t)) * old_elastic;
                let stiffness_work = if $soft {
                    0.5 * lambda * old_trace * old_trace
                        + mu * (1.0 + 2.0 * poisson(old_t).powi(2)) * old_elastic.powi(2)
                        - 0.5 * youngs(old_t) * old_elastic.powi(2)
                } else {
                    0.0
                };
                close(energy.creep_dissipation_j_m3, stress.abs() * dc, 1e-7);
                close(energy.plastic_dissipation_j_m3, yield_stress(t) * dp, 1e-7);
                close(
                    energy.initial_stored_j_m3,
                    0.5 * old_stress * old_stress / youngs(old_t)
                        + 0.5 * H * expected.plastic_path.powi(2),
                    1e-7,
                );
                close(
                    energy.final_stored_j_m3,
                    0.5 * stress * stress / youngs(t) + 0.5 * H * next.plastic_path.powi(2),
                    1e-7,
                );
                close(
                    energy.reversible_heat_j_m3,
                    stress * (free(t) - free(old_t)) - stiffness_work - entropy_storage,
                    1e-7,
                );
                close(
                    energy.thermoelastic_storage_change_j_m3,
                    entropy_storage,
                    1e-7,
                );
                close(
                    energy.mechanical_work_j_m3,
                    stress * (strain - old_strain),
                    1e-7,
                );
                let delta_axial = elastic - old_elastic;
                let delta_lateral = -nu * elastic + poisson(old_t) * old_elastic;
                close(
                    energy.numerical_dissipation_j_m3,
                    0.5 * lambda * (delta_axial + 2.0 * delta_lateral).powi(2)
                        + mu * (delta_axial.powi(2) + 2.0 * delta_lateral.powi(2))
                        + 0.5 * H * dp * dp,
                    1e-7,
                );
                close(energy.residual_j_m3, 0.0, tolerances().energy_j_m3);
                let native = law
                    .tensor_law()
                    .advance_with_energy(
                        &law.tensor_strain(old_strain, &old).unwrap(),
                        &ThermalCreepState::from_parts(old.mechanical.tensor_history(), old_t),
                        &response.strain,
                        tolerances(),
                    )
                    .unwrap();
                assert_eq!(energy, native.energy);
                let c = checked.tangent_pa;
                let determinant = c[1][1] * c[2][2] - c[1][2] * c[2][1];
                let dy = (c[2][2] * c[1][0] - c[1][2] * c[2][0]) / determinant;
                let dz = (c[1][1] * c[2][0] - c[2][1] * c[1][0]) / determinant;
                close(
                    response.tangent_pa,
                    c[0][0] - c[0][1] * dy - c[0][2] * dz,
                    3e-5,
                );
                let h = 1e-8;
                close(
                    response.tangent_pa,
                    (law.stress(strain + h, &old) - law.stress(strain - h, &old)) / (2.0 * h),
                    1e-5 * E,
                );

                assert_eq!(law.material().unwrap().properties().len(), $claims);
                assert_eq!(
                    law.material().unwrap(),
                    law.tensor_law().mechanical().material()
                );
                assert_eq!(law.kinetics().duration_s(), dt);
                close(law.kinetics().elastic().youngs, youngs(t), 1e-6);
                close(law.kinetics().elastic().poisson, poisson(t), 1e-15);
                close(
                    law.mechanical_law().backbone().recovery_rate_s(),
                    0.8 * (15_000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / 400.0 - 1.0 / t)).exp(),
                    1e-14,
                );
                assert_eq!(
                    law.mechanical_law()
                        .backbone()
                        .recovery_parameters()
                        .reference_rate_s,
                    0.8
                );
                if step == 0 || step == 1 {
                    assert_eq!(energy.mechanical_work_j_m3, 0.0);
                    if !$soft && $joint && step == 1 {
                        // Thermal yield softening releases compression faster
                        // than expansion rebuilds it, so reciprocal heat is negative.
                        assert!(energy.reversible_heat_j_m3 < 0.0);
                    } else if !$soft {
                        assert!(energy.reversible_heat_j_m3 > 0.0);
                    }
                }
                if step == 2 || step == 3 {
                    assert!(next.resistance < expected.resistance);
                    close(energy.physical_dissipation_j_m3, 0.0, 1e-9);
                    close(actual.creep, expected.creep, 5e-14);
                    close(actual.plastic, expected.plastic, 5e-14);
                    if step == 3 {
                        close(energy.caloric_heat_j_m3(), 0.0, 1e-7);
                    }
                } else {
                    assert!(dc > 0.0);
                    if $joint && step > 0 {
                        assert!(dp > 0.0);
                    }
                }
                old = checked.state;
                expected = next;
                old_t = t;
                old_strain = strain;
                old_stress = stress;
                old_tensor = response.strain;
            }
            if $soft {
                assert!(
                    source
                        .mechanical_source()
                        .at_temperature(400.0, 1.0)
                        .is_err()
                );
                assert!(
                    source
                        .minimum_thermal_capacity_j_m3_k()
                        .unwrap()
                        .is_finite()
                );
            }
            assert!(source.at_uniaxial_temperature(500.01, 0.25).is_err());
            assert!(source.at_uniaxial_temperature(299.99, 0.25).is_err());
            assert!(source.at_uniaxial_temperature(500.0, 0.0).is_err());
            let law = source.at_uniaxial_temperature(500.0, 0.25).unwrap();
            let poisoned = law.update_state(f64::NAN, &old);
            assert!(poisoned.temperature_k().is_nan());
            assert!(law.response(old_strain, &poisoned).is_err());
            assert_eq!(poisoned.mechanical, old.mechanical);
        }
    };
}

trajectory!(
    heated_recovering_fiber_matches_native_thermal_energy_and_axial_reference,
    pure_source,
    pure_history,
    false,
    13,
    false
);
trajectory!(
    heated_recovering_joint_fiber_yields_then_preserves_both_memories_during_dwell,
    joint_source,
    joint_history,
    true,
    15,
    false
);

trajectory!(
    g1_softening_recovering_fiber_retains_old_poisson_history_and_recovery_kinetics,
    soft_pure_source,
    pure_history,
    false,
    13,
    true
);
trajectory!(
    g1_softening_recovering_joint_fiber_closes_native_storage_heat_and_dwell,
    soft_joint_source,
    joint_history,
    true,
    15,
    true
);

#[test]
fn g0_zero_recovery_preserves_exact_primary_thermocreep_under_nonproportional_loading() {
    let card = source_card(true, 0.0, None);
    let recovering = soft_pure_source(&card);
    let primary = ExpandingCreepSource::from_norton_with_strain_hardening(
        &card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        &expansion(&card),
        || false,
    )
    .unwrap();
    let mut primary_state = primary.at_temperature(300.0, 1.0).unwrap().initial_state();
    let mut recovering_state = recovering
        .at_temperature(300.0, 1.0)
        .unwrap()
        .initial_state();
    let mut old_strain = [0.0; 6];
    for (t, dt, strain) in [
        (300.0, 0.25, [0.004, -0.001, -0.002, 0.0008, 0.0, 0.0002]),
        (450.0, 0.5, [0.002, -0.0005, -0.001, -0.0004, 0.0, 0.0007]),
        (350.0, 0.75, [-0.001, 0.003, -0.001, 0.0, 0.0006, -0.0003]),
    ] {
        let expected = primary
            .at_temperature(t, dt)
            .unwrap()
            .advance_with_energy(&old_strain, &primary_state, &strain, tolerances())
            .unwrap();
        let actual = recovering
            .at_temperature(t, dt)
            .unwrap()
            .advance_with_energy(&old_strain, &recovering_state, &strain, tolerances())
            .unwrap();
        assert_eq!(actual.stress_pa, expected.stress_pa);
        assert_eq!(actual.tangent_pa, expected.tangent_pa);
        assert_eq!(actual.energy, expected.energy);
        assert_eq!(actual.state.mechanical.native, expected.state.mechanical);
        assert_eq!(
            actual.state.mechanical.resistance_strain,
            expected.state.mechanical.equivalent_creep_strain,
        );
        old_strain = strain;
        primary_state = expected.state;
        recovering_state = actual.state;
    }
}

#[test]
fn g0_recovering_thermal_sources_keep_recovery_claims_hardening_and_storage_admission() {
    let lower = point(300.0);
    let upper = point(500.0);
    let query = || NortonSourceQuery {
        lower: &lower,
        upper: &upper,
        temperature_axis: "T",
        selection: MaterialPropertySelection::SingleClaimOnly,
        strain_limit: 0.03,
    };
    let card = source_card(true, 0.8, None);
    assert!(RecoveringNortonSource::new(&card, query(), || false).is_err());
    assert!(
        RecoveringNortonJ2Source::new(
            &card,
            query(),
            HardeningPropertyBasis::PlasticStrainModulus,
            || false,
        )
        .is_err()
    );
    let joint = soft_joint_source(&card);
    assert!(
        joint
            .mechanical_source()
            .recovery_source()
            .at_temperature(400.0, 1.0)
            .is_err()
    );
    assert!(joint.creep_source().at_temperature(400.0, 1.0).is_err());
    assert!(matches!(
        ExpandingCreepSource::from_recovering_norton_j2(
            &card,
            query(),
            HardeningPropertyBasis::UniaxialTangentModulus,
            &expansion(&card),
            || false,
        ),
        Err(CreepSourceError::Invalid(_))
    ));
    for (property, offset) in [
        (CREEP_RECOVERY_REFERENCE_RATE_PROPERTY, 0.1),
        (CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY, 1000.0),
        (CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY, 10.0),
        (YOUNG_MODULUS_PROPERTY, 1e8),
    ] {
        let changed = source_card(true, 0.8, Some((property, offset)));
        assert!(
            ExpandingCreepSource::from_recovering_norton(
                &changed,
                query(),
                &expansion(&changed),
                || false,
            )
            .is_err(),
            "interior {property} change was admitted"
        );
    }
    let changed = source_card(true, 0.8, Some((ISOTROPIC_HARDENING_MODULUS_PROPERTY, 1e7)));
    assert!(matches!(
        ExpandingCreepSource::from_recovering_norton_j2(
            &changed,
            query(),
            HardeningPropertyBasis::PlasticStrainModulus,
            &expansion(&changed),
            || false,
        ),
        Err(CreepSourceError::Invalid(_))
    ));
    assert!(matches!(
        ExpandingCreepSource::from_recovering_norton(&card, query(), &expansion(&card), || true,),
        Err(CreepSourceError::Cancelled)
    ));
    assert!(matches!(
        ExpandingCreepSource::from_recovering_norton_j2(
            &card,
            query(),
            HardeningPropertyBasis::PlasticStrainModulus,
            &expansion(&card),
            || true,
        ),
        Err(CreepSourceError::Cancelled)
    ));
}
