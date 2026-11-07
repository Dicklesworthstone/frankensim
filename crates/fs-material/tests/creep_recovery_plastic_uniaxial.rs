//! Exact free-sided recovering creep/J2 against independent scalar solutions,
//! retained load/dwell/reload histories and original source/energy ownership.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::plastic::uniaxial::UniaxialNortonJ2;
use fs_material::creep::recovery::plastic::RecoveringNortonJ2;
use fs_material::creep::recovery::plastic::source::RecoveringNortonJ2Source;
use fs_material::creep::recovery::plastic::uniaxial::UniaxialRecoveringNortonJ2;
use fs_material::creep::recovery::source::{
    CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY, CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
    CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
};
use fs_material::creep::recovery::uniaxial::UniaxialRecoveringNorton;
use fs_material::creep::recovery::{CreepRecoveryParameters, RecoveringNortonCreep};
use fs_material::creep::source::*;
use fs_material::creep::{
    CreepTolerances, GAS_CONSTANT_J_MOL_K, NortonCreep, NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_material::{IsotropicElastic, Uniaxial};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const STRESS: f64 = 3e6;
const RATE: f64 = 2e-3;
const A: f64 = 8e-4;
const YIELD: f64 = 2e6;
const H: f64 = 2e8;

fn primary(n: f64, dt: f64) -> NortonCreep {
    NortonCreep::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: n,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 500.0,
            temperature_range_k: [300.0, 900.0],
        },
        500.0,
        dt,
    )
    .unwrap()
    .with_strain_hardening(StrainHardeningCreep {
        reference_strain: A,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    })
    .unwrap()
}

fn recovery(rate: f64) -> CreepRecoveryParameters {
    CreepRecoveryParameters {
        reference_rate_s: rate,
        activation_energy_j_mol: 0.0,
        reference_temperature_k: 500.0,
    }
}

fn law(n: f64, dt: f64, hardening: f64, rate: f64) -> UniaxialRecoveringNortonJ2 {
    let creep = primary(n, dt);
    let plastic = J2Plasticity::new(creep.elastic().clone(), YIELD, hardening).unwrap();
    UniaxialRecoveringNortonJ2::new(
        RecoveringNortonJ2::new(NortonJ2::new(creep, plastic).unwrap(), recovery(rate)).unwrap(),
    )
}

fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 3e-6,
        energy_j_m3: 2e-7,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tolerance {tolerance})"
    );
}

fn quadratic(b: f64, c: f64) -> f64 {
    2.0 * c / (b + (b * b + 4.0 * c).sqrt())
}

#[test]
fn axial_joint_recovery_matches_independent_quadratic_in_tension_and_reversal() {
    let dt = 0.25;
    let rate = 1.3;
    let divisor = 1.0 + rate * dt;
    let effective_a = A * divisor;
    let coefficient = dt * RATE * effective_a / STRESS;
    for hardening in [0.0, H] {
        let axial = law(1.0, dt, hardening, rate);
        let old_strain = 0.0025;
        let old = axial
            .response(old_strain, &axial.initial_state())
            .unwrap()
            .state;
        let frozen = old.clone();
        for strain in [-0.006, 0.008] {
            let trial = E * (strain - old.creep_strain() - old.plastic_strain());
            let old_yield = YIELD + hardening * old.equivalent_plastic_strain();
            let weight = hardening / (E + hardening);
            let effective_trial = old_yield + weight * (trial.abs() - old_yield);
            let modulus = E * weight;
            let dc = quadratic(
                effective_a + old.resistance_strain() + coefficient * modulus,
                coefficient * effective_trial,
            );
            let q = effective_trial - modulus * dc;
            let dp = (trial.abs() - old_yield - E * dc) / (E + hardening);
            assert!(dc > 0.0 && dp > 0.0);
            let response = axial.response(strain, &old).unwrap();
            let checked = axial
                .advance_with_energy(old_strain, &old, strain, tolerances())
                .unwrap();
            close(response.stress_pa, trial.signum() * q, 5e-8);
            let tangent = E * weight
                / (1.0
                    + modulus * coefficient / (effective_a + old.resistance_strain() + 2.0 * dc));
            close(response.tangent_pa, tangent, 3e-5);
            close(
                response.state.creep_strain(),
                old.creep_strain() + trial.signum() * dc,
                3e-17,
            );
            close(
                response.state.plastic_strain(),
                old.plastic_strain() + trial.signum() * dp,
                3e-17,
            );
            close(
                response.state.equivalent_creep_strain(),
                old.equivalent_creep_strain() + dc,
                3e-17,
            );
            close(
                response.state.equivalent_plastic_strain(),
                old.equivalent_plastic_strain() + dp,
                3e-17,
            );
            close(
                response.state.resistance_strain(),
                (old.resistance_strain() + dc) / divisor,
                3e-17,
            );
            close(
                response.strain[1],
                -NU * response.stress_pa / E
                    - 0.5 * (response.state.creep_strain() + response.state.plastic_strain()),
                3e-17,
            );
            assert_eq!(response.strain[1], response.strain[2]);
            assert_eq!(checked.state, response.state);
            close(checked.energy.creep_dissipation_j_m3, q * dc, 2e-8);
            close(checked.energy.plastic_dissipation_j_m3, YIELD * dp, 2e-8);
            close(
                checked.energy.physical_dissipation_j_m3,
                q * dc + YIELD * dp,
                2e-8,
            );
            close(checked.energy.residual_j_m3, 0.0, tolerances().energy_j_m3);
            let native = axial
                .backbone()
                .advance_with_energy(
                    &axial.tensor_strain(old_strain, &old).unwrap(),
                    &old.tensor_history(),
                    &response.strain,
                    tolerances(),
                )
                .unwrap();
            assert_eq!(checked.energy, native.energy);
            assert_eq!(checked.stress_pa, native.stress_pa);
            for stress in &checked.stress_pa[1..] {
                close(*stress, 0.0, 5e-8);
            }
        }
        assert_eq!(old, frozen);
        assert_eq!(axial.kinetics(), axial.backbone().norton_j2().creep());
        assert_eq!(
            axial
                .kinetics()
                .strain_hardening()
                .unwrap()
                .reference_strain,
            A
        );
        assert!(axial.material().is_none());
    }
}

#[test]
fn subyield_joint_response_reduces_to_pure_recovering_axial_flow() {
    let axial = law(1.0, 0.25, H, 1.3);
    let pure = UniaxialRecoveringNorton::new(
        RecoveringNortonCreep::new(axial.kinetics().clone(), recovery(1.3)).unwrap(),
    );
    for strain in [0.0, 0.001, -0.001] {
        let next = axial.response(strain, &axial.initial_state()).unwrap();
        let expected = pure.response(strain, &pure.initial_state()).unwrap();
        assert_eq!(next.stress_pa, expected.stress_pa);
        assert_eq!(next.tangent_pa, expected.tangent_pa);
        assert_eq!(next.strain, expected.strain);
        assert_eq!(next.state.creep_strain(), expected.state.creep_strain());
        assert_eq!(
            next.state.equivalent_creep_strain(),
            expected.state.equivalent_creep_strain()
        );
        assert_eq!(
            next.state.resistance_strain(),
            expected.state.resistance_strain()
        );
        assert_eq!(next.state.plastic_strain(), 0.0);
        assert_eq!(next.state.equivalent_plastic_strain(), 0.0);
    }
}

#[test]
fn nonlinear_joint_axial_tangent_matches_finite_difference_and_tensor_condensation() {
    for hardening in [0.0, H] {
        let axial = law(3.2, 0.125, hardening, 1.7);
        let old = axial.response(0.004, &axial.initial_state()).unwrap().state;
        let frozen = old.clone();
        let free = old.creep_strain() + old.plastic_strain();
        for strain in [0.007, -0.007, free + 0.0002] {
            let next = axial.response(strain, &old).unwrap();
            let step = 1e-8;
            let fd = (axial.stress(strain + step, &old) - axial.stress(strain - step, &old))
                / (2.0 * step);
            close(next.tangent_pa, fd, 2e-6 * next.tangent_pa.abs().max(1e6));
            let full = axial
                .backbone()
                .response(&next.strain, &old.tensor_history())
                .unwrap();
            close(full.stress_pa[0], next.stress_pa, 5e-8);
            let c = full.tangent_pa;
            let determinant = c[1][1] * c[2][2] - c[1][2] * c[2][1];
            let dy = (c[2][2] * c[1][0] - c[1][2] * c[2][0]) / determinant;
            let dz = (c[1][1] * c[2][0] - c[2][1] * c[1][0]) / determinant;
            close(next.tangent_pa, c[0][0] - c[0][1] * dy - c[0][2] * dz, 3e-5);
            assert_eq!(next, axial.response(strain, &old).unwrap());
            axial
                .advance_with_energy(0.004, &old, strain, tolerances())
                .unwrap();
        }
        assert_eq!(old, frozen);
    }
}

#[test]
fn zero_recovery_is_exact_native_joint_axial_behavior_through_holds_and_reversals() {
    for hardening in [0.0, H] {
        let axial = law(3.2, 0.25, hardening, 0.0);
        let original = UniaxialNortonJ2::new(axial.backbone().norton_j2().clone());
        let mut state = axial.initial_state();
        let mut old_strain = 0.0;
        for strain in [0.004, 0.004, -0.006, 0.007] {
            let expected = original.response(strain, state.native()).unwrap();
            let next = axial.response(strain, &state).unwrap();
            assert_eq!(next.strain, expected.strain);
            assert_eq!(next.stress_pa, expected.stress_pa);
            assert_eq!(next.tangent_pa, expected.tangent_pa);
            assert_eq!(next.state.native(), &expected.state);
            assert_eq!(
                next.state.resistance_strain(),
                expected.state.equivalent_creep_strain()
            );
            let expected = original
                .advance_with_energy(old_strain, state.native(), strain, tolerances())
                .unwrap();
            let checked = axial
                .advance_with_energy(old_strain, &state, strain, tolerances())
                .unwrap();
            assert_eq!(checked.energy, expected.energy);
            assert_eq!(checked.stress_pa, expected.stress_pa);
            assert_eq!(checked.tangent_pa, expected.tangent_pa);
            state = next.state;
            old_strain = strain;
        }
    }
}

#[test]
fn constant_stress_hold_dwell_and_reload_keep_both_permanent_paths() {
    let dt = 0.2;
    let rate = 1.0;
    let axial = law(1.0, dt, H, rate);
    let sigma = 3e6;
    let plastic = (sigma - YIELD) / H;
    let divisor = 1.0 + rate * dt;
    let effective_a = A * divisor;
    let mut old = axial.initial_state();
    let mut old_strain = 0.0;
    for step in 0..5 {
        let dc = quadratic(
            effective_a + old.resistance_strain(),
            dt * RATE * (sigma / STRESS) * effective_a,
        );
        let strain = sigma / E + plastic + old.creep_strain() + dc;
        let next = axial
            .advance_with_energy(old_strain, &old, strain, tolerances())
            .unwrap();
        close(next.stress_pa[0], sigma, 1e-7);
        close(
            next.state.equivalent_creep_strain(),
            old.equivalent_creep_strain() + dc,
            3e-17,
        );
        close(next.state.plastic_strain(), plastic, 3e-17);
        close(next.state.equivalent_plastic_strain(), plastic, 3e-17);
        close(
            next.state.resistance_strain(),
            (old.resistance_strain() + dc) / divisor,
            3e-17,
        );
        close(next.energy.creep_dissipation_j_m3, sigma * dc, 2e-8);
        close(
            next.energy.plastic_dissipation_j_m3,
            if step == 0 { YIELD * plastic } else { 0.0 },
            2e-8,
        );
        old = next.state;
        old_strain = strain;
    }
    let free = old.creep_strain() + old.plastic_strain();
    let unloaded = axial
        .advance_with_energy(old_strain, &old, free, tolerances())
        .unwrap();
    close(unloaded.energy.physical_dissipation_j_m3, 0.0, 1e-12);
    let mut rested = unloaded.state.clone();
    for step in 1..=8 {
        let dwell = axial
            .advance_with_energy(free, &rested, free, tolerances())
            .unwrap();
        close(dwell.state.creep_strain(), old.creep_strain(), 3e-17);
        close(dwell.state.plastic_strain(), old.plastic_strain(), 3e-17);
        close(
            dwell.state.equivalent_creep_strain(),
            old.equivalent_creep_strain(),
            3e-17,
        );
        close(
            dwell.state.equivalent_plastic_strain(),
            old.equivalent_plastic_strain(),
            3e-17,
        );
        close(
            dwell.state.resistance_strain(),
            unloaded.state.resistance_strain() / divisor.powi(step),
            3e-17,
        );
        close(dwell.energy.physical_dissipation_j_m3, 0.0, 1e-12);
        rested = dwell.state;
    }
    let reloaded_strain = free + 0.006;
    let immediate = axial.response(reloaded_strain, &unloaded.state).unwrap();
    let delayed = axial
        .advance_with_energy(free, &rested, reloaded_strain, tolerances())
        .unwrap();
    assert!(delayed.state.equivalent_creep_strain() > immediate.state.equivalent_creep_strain());
    assert!(
        delayed.state.equivalent_plastic_strain() < immediate.state.equivalent_plastic_strain()
    );
    assert!(delayed.state.equivalent_plastic_strain() > rested.equivalent_plastic_strain());
    assert!(delayed.stress_pa[0] < immediate.stress_pa);
    let reversed = axial
        .advance_with_energy(reloaded_strain, &delayed.state, -0.010, tolerances())
        .unwrap();
    assert!(reversed.state.creep_strain() < delayed.state.creep_strain());
    assert!(reversed.state.plastic_strain() < delayed.state.plastic_strain());
    assert!(reversed.state.equivalent_creep_strain() > delayed.state.equivalent_creep_strain());
    assert!(reversed.state.equivalent_plastic_strain() > delayed.state.equivalent_plastic_strain());
    assert!(reversed.energy.creep_dissipation_j_m3 > 0.0);
    assert!(reversed.energy.plastic_dissipation_j_m3 > 0.0);
}

#[test]
fn exact_creep_increment_survives_signed_and_accumulated_path_roundoff() {
    let loading = law(1.0, 0.2, H, 1.0);
    let loaded = loading
        .response(0.006, &loading.initial_state())
        .unwrap()
        .state;
    let free = loaded.creep_strain() + loaded.plastic_strain();
    let rested = law(1.0, 1e3, H, 1e15)
        .response(free, &loaded)
        .unwrap()
        .state;
    let short = law(1.0, 1e-17, H, 1e17);
    let strain = free + 1e6 / E;
    let next = short
        .advance_with_energy(strain, &rested, strain, tolerances())
        .unwrap();
    assert_eq!(next.state.creep_strain(), rested.creep_strain());
    assert_eq!(
        next.state.equivalent_creep_strain(),
        rested.equivalent_creep_strain()
    );
    assert_eq!(next.state.plastic_strain(), rested.plastic_strain());
    assert_eq!(
        next.state.equivalent_plastic_strain(),
        rested.equivalent_plastic_strain()
    );
    assert!(next.state.resistance_strain() > rested.resistance_strain());
    assert!(next.energy.creep_dissipation_j_m3 > 0.0);
    assert_eq!(next.energy.plastic_dissipation_j_m3, 0.0);
    let dc = 2.0 * next.state.resistance_strain() - rested.resistance_strain();
    let sigma = next.stress_pa[0].abs();
    let rate = short
        .kinetics()
        .equivalent_creep_rate_at_strain_s(sigma, next.state.resistance_strain())
        .unwrap();
    close(dc, 1e-17 * rate, 2e-32);
    close(next.energy.creep_dissipation_j_m3, sigma * dc, 2e-26);
    let frozen = rested.clone();
    for invalid in [f64::NAN, f64::INFINITY, 1.0] {
        assert!(short.response(invalid, &rested).is_err());
        assert!(short.stress(invalid, &rested).is_nan());
    }
    assert_eq!(rested, frozen);
}

fn sourced_card() -> MaterialCard {
    let coefficients = [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), NU),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, STRESS),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, RATE),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
        (
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            30_000.0,
        ),
        (
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            500.0,
        ),
        (CREEP_HARDENING_REFERENCE_STRAIN_PROPERTY, Dims([0; 6]), A),
        (CREEP_HARDENING_EXPONENT_PROPERTY, Dims([0; 6]), 1.0),
        (CREEP_RESIDUAL_RATE_FRACTION_PROPERTY, Dims([0; 6]), 0.0),
        (CREEP_RECOVERY_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 0.6),
        (
            CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            40_000.0,
        ),
        (
            CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            500.0,
        ),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
    ];
    let mut claims = ClaimSet::new();
    for (name, dims, value) in coefficients {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    700.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured recovering axial creep/plasticity fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic recovering creep/plastic solid".into(),
            phase: "solid".into(),
            process: "manufactured reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn source_axial_refresh_keeps_fifteen_receipts_both_histories_and_independent_recovery() {
    let card = sourced_card();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
            .unwrap()
    };
    let source = RecoveringNortonJ2Source::new(
        &card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(700.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        },
        HardeningPropertyBasis::PlasticStrainModulus,
        || false,
    )
    .unwrap();
    let resolved = source.at_temperature(500.0, 0.1).unwrap();
    let cold = UniaxialRecoveringNortonJ2::from_source(resolved.clone());
    assert_eq!(cold.backbone(), resolved.law());
    assert_eq!(cold.material().unwrap(), resolved.material());
    assert_eq!(
        cold.material().unwrap().card_identity(),
        card.content_hash()
    );
    assert_eq!(cold.material().unwrap().properties().len(), 15);
    let cloned = cold.clone();
    assert!(std::ptr::eq(
        cold.material().unwrap(),
        cloned.material().unwrap()
    ));
    let strain = 0.006;
    let first = cold
        .advance_with_energy(0.0, &cold.initial_state(), strain, tolerances())
        .unwrap();
    let frozen = first.state.clone();
    let hot = source.at_uniaxial_temperature(700.0, 0.2).unwrap();
    let expected_recovery =
        0.6 * f64::exp(40_000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / 500.0 - 1.0 / 700.0));
    close(
        hot.backbone().recovery_rate_s(),
        expected_recovery,
        3e-13 * expected_recovery,
    );
    let next = hot
        .advance_with_energy(strain, &first.state, strain, tolerances())
        .unwrap();
    assert!(next.state.equivalent_creep_strain() > first.state.equivalent_creep_strain());
    assert_eq!(
        next.state.equivalent_plastic_strain(),
        first.state.equivalent_plastic_strain()
    );
    assert!(next.stress_pa[0] < first.stress_pa[0]);
    assert!(next.energy.creep_dissipation_j_m3 > 0.0);
    close(next.energy.plastic_dissipation_j_m3, 0.0, 1e-12);
    close(next.energy.mechanical_work_j_m3, 0.0, 2e-7);
    let dc = next.state.equivalent_creep_strain() - first.state.equivalent_creep_strain();
    close(
        next.state.resistance_strain(),
        (first.state.resistance_strain() + dc) / (1.0 + 0.2 * expected_recovery),
        3e-17,
    );
    for property in hot.material().unwrap().properties() {
        let original = cold
            .material()
            .unwrap()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            original.answer().receipt.selected
        );
    }
    assert_eq!(first.state, frozen);
    assert!(source.at_uniaxial_temperature(701.0, 0.1).is_err());
    assert!(source.at_uniaxial_temperature(500.0, 0.0).is_err());
}
