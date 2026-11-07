//! Manufactured thermal fibers: permanent creep, yield and recoverable resistance.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::Uniaxial;
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
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const E: f64 = 1e9;
const NU: f64 = 0.25;
const EXPANSION: f64 = 2e-5;
const STRESS: f64 = 2e6;
const RATE: f64 = 1e-4;
const PRIMARY: f64 = 2e-4;
const H: f64 = 2e8;

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
        (CREEP_RECOVERY_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 0.8),
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
        claims.insert_claim(claim(name, dims, value)).unwrap();
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
    let rate = RATE * (20_000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / 400.0 - 1.0 / t)).exp();
    let recovery = 0.8 * (15_000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / 400.0 - 1.0 / t)).exp();
    let g = 1.0 + recovery * dt;
    let a = PRIMARY * g;
    let coefficient = dt * rate * a / STRESS;
    let trial = E * (strain - EXPANSION * (t - 300.0) - old.creep - old.plastic);
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
        if q + E * (dc + dp) > trial.abs() {
            hi = q;
        } else {
            lo = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let (dc, dp) = increments(q);
    let tangent = E
        / (1.0
            + E * (coefficient / (a + old.resistance + 2.0 * dc)
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
    ($name:ident, $source:ident, $history:ident, $joint:expr, $claims:expr) => {
        #[test]
        fn $name() {
            let card = card();
            let source = $source(&card);
            let cold = source.at_uniaxial_temperature(300.0, 0.25).unwrap();
            let mut old = cold.initial_state();
            let mut expected = History::default();
            let (mut old_t, mut old_strain, mut old_stress) = (300.0, 0.0, 0.0);
            for (step, (t, dt)) in [
                (400.0, 0.25),
                (500.0, 0.25),
                (500.0, 0.5),
                (500.0, 2.0),
                (500.0, 0.25),
            ]
            .into_iter()
            .enumerate()
            {
                let law = source.at_uniaxial_temperature(t, dt).unwrap();
                let strain = match step {
                    2 | 3 => EXPANSION * (t - 300.0) + expected.creep + expected.plastic,
                    4 => -0.001,
                    _ => 0.0,
                };
                let (next, stress, tangent) = reference($joint, t, dt, strain, expected);
                let frozen = old.clone();
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
                    EXPANSION * (t - 300.0) - NU * stress / E - 0.5 * (next.creep + next.plastic),
                    5e-14,
                );
                assert_eq!(response.strain[1], response.strain[2]);
                assert!(checked.stress_pa[1..].iter().all(|s| s.abs() < 2e-6));

                let dc = next.creep_path - expected.creep_path;
                let dp = next.plastic_path - expected.plastic_path;
                let energy = checked.energy;
                close(energy.creep_dissipation_j_m3, stress.abs() * dc, 1e-7);
                close(energy.plastic_dissipation_j_m3, yield_stress(t) * dp, 1e-7);
                close(
                    energy.final_stored_j_m3,
                    0.5 * stress * stress / E + 0.5 * H * next.plastic_path.powi(2),
                    1e-7,
                );
                close(
                    energy.reversible_heat_j_m3,
                    -EXPANSION * old_t * (stress - old_stress),
                    1e-7,
                );
                close(
                    energy.thermoelastic_storage_change_j_m3,
                    EXPANSION * (t * stress - old_t * old_stress),
                    1e-7,
                );
                close(
                    energy.mechanical_work_j_m3,
                    stress * (strain - old_strain),
                    1e-7,
                );
                close(
                    energy.numerical_dissipation_j_m3,
                    0.5 * (stress - old_stress).powi(2) / E + 0.5 * H * dp * dp,
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
                assert_eq!(
                    law.mechanical_law()
                        .backbone()
                        .recovery_parameters()
                        .reference_rate_s,
                    0.8
                );
                if step == 0 || step == 1 {
                    assert_eq!(energy.mechanical_work_j_m3, 0.0);
                    if $joint && step == 1 {
                        // Thermal yield softening releases compression faster
                        // than expansion rebuilds it, so reciprocal heat is negative.
                        assert!(energy.reversible_heat_j_m3 < 0.0);
                    } else {
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
    13
);
trajectory!(
    heated_recovering_joint_fiber_yields_then_preserves_both_memories_during_dwell,
    joint_source,
    joint_history,
    true,
    15
);
