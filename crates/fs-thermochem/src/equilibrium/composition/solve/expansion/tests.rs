use super::*;
use crate::{ElementId, ElementalReferenceIdV1, Nasa9RegionV1,
    ReferenceEquationOfStateV1, StandardStatePhaseV1, required_nasa9_parameter_dims};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use std::collections::BTreeMap;

const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
fn id(s: &str) -> SpeciesId { SpeciesId::new(s).unwrap() }
fn model(name: &str, mass: f64, cp_r: f64, h_r: f64, s_r: f64) -> Nasa9StandardStateModelV1 {
    let coefficients = [0.0, 0.0, cp_r, 0.0, 0.0, 0.0, 0.0, h_r, s_r];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard {
        law: LawId(crate::NASA9_LAW_ID_V1.into()), law_version: crate::NASA9_LAW_VERSION_V1,
        parameters, state_schema_version: crate::NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0), sources: Vec::new(),
        provenance: Provenance { source: "manufactured reacting expansion".into(), license: "test fixture".into(), artifact: None },
    };
    Nasa9StandardStateModelV1::new(id(name), MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("expansion-reference").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}
fn owner(models: Vec<Nasa9StandardStateModelV1>, counts: Vec<i128>) -> IdealGasEquilibriumV1 {
    let species: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], species.clone(), vec![counts]).unwrap(),
        ChargeVector::new(species.clone(), vec![0; species.len()]).unwrap(), models).unwrap()
}
fn inert() -> IdealGasEquilibriumV1 { owner(vec![model("H2", 0.002, 3.5, 0.0, 0.0)], vec![2]) }
fn dissociation(energy_shift: f64, entropy_shift: f64) -> IdealGasEquilibriumV1 {
    owner(vec![model("H", 0.001, 2.5, 12000.0 + energy_shift, 1.0 + entropy_shift),
        model("H2", 0.002, 3.5, 2.0 * energy_shift, 2.0 * entropy_shift)], vec![1, 2])
}
fn controls(scale: f64) -> EquilibriumEntropySolveV1 {
    EquilibriumEntropySolveV1 { equilibrium: EquilibriumSolveV1::default(),
        minimum_temperature_k: 200.0, maximum_temperature_k: 6000.0,
        entropy_tolerance_j_k: 1e-8 * scale, maximum_temperature_iterations: 100 }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:.17e} != {expected:.17e}");
}
fn tp(o: &IdealGasEquilibriumV1, t: f64, p: f64, scale: f64) -> EquilibriumThermodynamicStateV1 {
    o.equilibrate_tp_with_response(Temperature::new(t), Pressure::new(p),
        &[(id("H2"), Amount::new(scale))], EquilibriumSolveV1::default(), || false).unwrap()
}

#[test]
fn g1_inert_gas_response_and_sp_reduce_to_independent_perfect_gas_laws() {
    let o = inert();
    let initial = tp(&o, 1600.0, 800000.0, 2.0);
    let r = initial.response();
    close(r.heat_capacity_p_j_k, 7.0 * R, 1e-10);
    close(r.heat_capacity_v_j_k, 5.0 * R, 1e-10);
    close(r.thermal_expansion_per_k, 1.0 / 1600.0, 1e-10);
    close(r.isothermal_compressibility_per_pa, 1.0 / 800000.0, 1e-10);
    close(r.sound_speed_m_s, (1.4 * R / 0.002 * 1600.0).sqrt(), 1e-10);
    close(r.isentropic_pressure_exponent, 1.4, 1e-10);
    for p in [100000.0, 400000.0, 800000.0, 1600000.0] {
        let result = o.equilibrate_sp(Pressure::new(p), &[(id("H2"), Amount::new(2.0))],
            initial.entropy_j_k(), controls(2.0), || false).unwrap();
        close(result.equilibrium().temperature().value(), 1600.0 * (p / 800000.0_f64).powf(2.0 / 7.0), 2e-9);
        assert!(result.entropy_residual_j_k().unwrap().abs() <= controls(2.0).entropy_tolerance_j_k);
        assert!(result.equilibrium().species()[0].source.is_some());
    }
}

#[test]
fn g1_dissociation_response_matches_closed_form_not_frozen_gamma() {
    let o = dissociation(0.0, 0.0);
    for t in [1600.0_f64, 2000.0, 2600.0] {
        for p in [10000.0_f64, 100000.0, 1000000.0] {
            let state = tp(&o, t, p, 1.0);
            // H2 <=> 2H: Kp=exp(1/2 - 24000/T + 3/2 ln(T)).
            let k = (0.5 - 24000.0 / t + 1.5 * t.ln()).exp();
            let d = (k / (k + 4.0 * p / 100000.0)).sqrt();
            let dt = 0.5 * d * (1.0 - d * d) * (1.5 + 24000.0 / t);
            let dp = -0.5 * d * (1.0 - d * d);
            let cp = R * (3.5 + 1.5 * d + (1.5 + 24000.0 / t) * dt);
            let a = 1.0 + dt / (1.0 + d);
            let b = 1.0 - dp / (1.0 + d);
            let density = 0.002 * p / ((1.0 + d) * R * t);
            let cv = cp - (1.0 + d) * R * a * a / b;
            let exponent = 1.0 / (b - (1.0 + d) * R * a * a / cp);
            let r = state.response();
            close(state.equilibrium().total_amount().value(), 1.0 + d, 2e-8);
            close(r.heat_capacity_p_j_k, cp, 2e-8);
            close(r.heat_capacity_v_j_k, cv, 2e-8);
            close(r.thermal_expansion_per_k, a / t, 2e-8);
            close(r.isothermal_compressibility_per_pa, b / p, 2e-8);
            close(r.density_kg_m3, density, 2e-8);
            close(r.sound_speed_m_s, (exponent * p / density).sqrt(), 2e-8);
            close(r.isentropic_pressure_exponent, exponent, 2e-8);
            assert!((cp - cv - (1.0 + d) * R).abs() > 0.01 * R);
            assert!((r.isentropic_pressure_exponent - cp / cv).abs() > 1e-4);
        }
    }
}

#[test]
fn g1_sp_sound_speed_and_maxwell_relation_match_independent_state_differences() {
    let o = dissociation(0.0, 0.0);
    let (t, p) = (2000.0, 100000.0);
    let center = tp(&o, t, p, 1.0);
    let dp = 10.0;
    let sp = |pressure| o.equilibrate_sp(Pressure::new(pressure), &[(id("H2"), Amount::new(1.0))],
        center.entropy_j_k(), controls(1.0), || false).unwrap();
    let low = sp(p - dp);
    let high = sp(p + dp);
    let c2_fd = 2.0 * dp / (high.response().density_kg_m3 - low.response().density_kg_m3);
    close(center.response().sound_speed_m_s.powi(2), c2_fd, 3e-6);
    let dt = 0.1;
    let hot = tp(&o, t + dt, p, 1.0);
    let cold = tp(&o, t - dt, p, 1.0);
    close((hot.entropy_j_k() - cold.entropy_j_k()) / (2.0 * dt),
        center.response().heat_capacity_p_j_k / t, 3e-6);
    let ds_dp = (tp(&o, t, p + dp, 1.0).entropy_j_k() - tp(&o, t, p - dp, 1.0).entropy_j_k()) / (2.0 * dp);
    let volume = 0.002 / center.response().density_kg_m3;
    close(ds_dp, -volume * center.response().thermal_expansion_per_k, 3e-6);
}

#[test]
fn g3_amount_scaling_and_element_reference_shifts_preserve_reacting_response() {
    let ordinary = dissociation(0.0, 0.0);
    let base = tp(&ordinary, 2000.0, 100000.0, 1.0);
    for scale in [1e-8, 1.0, 1e8] {
        let state = tp(&ordinary, 2000.0, 100000.0, scale);
        close(state.response().sound_speed_m_s, base.response().sound_speed_m_s, 2e-8);
        close(state.response().heat_capacity_p_j_k, scale * base.response().heat_capacity_p_j_k, 2e-8);
        close(state.entropy_j_k(), scale * base.entropy_j_k(), 2e-8);
    }
    let shifted = dissociation(100000.0, 100.0);
    let state = tp(&shifted, 2000.0, 100000.0, 1.0);
    close(state.response().sound_speed_m_s, base.response().sound_speed_m_s, 2e-8);
    close(state.entropy_j_k() - base.entropy_j_k(), 200.0 * R, 2e-8);
    let expand = |o: &IdealGasEquilibriumV1, s| o.equilibrate_sp(Pressure::new(30000.0),
        &[(id("H2"), Amount::new(1.0))], s, controls(1.0), || false).unwrap();
    close(expand(&ordinary, base.entropy_j_k()).equilibrium().temperature().value(),
        expand(&shifted, state.entropy_j_k()).equilibrium().temperature().value(), 2e-8);
}

#[test]
fn g3_g4_bad_entropy_domains_and_cancelled_trials_publish_nothing() {
    let o = inert();
    let snapshot = o.clone();
    let feed = [(id("H2"), Amount::new(1.0))];
    let target = tp(&o, 1600.0, 100000.0, 1.0).entropy_j_k();
    assert!(matches!(o.equilibrate_sp(Pressure::new(100000.0), &feed, 1e9, controls(1.0), || false),
        Err(EquilibriumExpansionErrorV1::EntropyNotBracketed { .. })));
    assert!(o.equilibrate_sp(Pressure::new(100000.0), &feed, f64::NAN, controls(1.0), || false).is_err());
    let mut bad = controls(1.0); bad.entropy_tolerance_j_k = 0.0;
    assert!(o.equilibrate_sp(Pressure::new(100000.0), &feed, target, bad, || false).is_err());
    for limit in [0, 3, 20] {
        let mut calls = 0;
        let result = o.equilibrate_sp(Pressure::new(30000.0), &feed, target, controls(1.0), || {
            calls += 1; calls > limit
        });
        assert!(matches!(result, Err(EquilibriumExpansionErrorV1::Equilibrium(E::Cancelled))));
    }
    assert_eq!(o, snapshot);
}

#[test]
fn g3_multiregion_sp_requires_the_original_explicit_caloric_policy() {
    let raw = model("H2", 0.002, 3.5, 0.0, 0.0);
    let mut lower = raw.regions()[0].card().clone();
    lower.validity = ValidityDomain::unconstrained().with("T", 200.0, 1000.0);
    let mut upper = lower.clone();
    upper.validity = ValidityDomain::unconstrained().with("T", 1000.0, 6000.0);
    let split = Nasa9StandardStateModelV1::new(id("H2"), MolarMass::new(0.002), raw.convention().clone(),
        vec![Nasa9RegionV1::from_card(lower).unwrap(), Nasa9RegionV1::from_card(upper).unwrap()]).unwrap();
    let corrected = split.with_caloric_continuity(Temperature::new(500.0),
        crate::MolarEnergyQuantityV1::new(0.0), crate::MolarThermalQuantityV1::new(0.0)).unwrap();
    let raw_owner = owner(vec![split], vec![2]);
    let corrected_owner = owner(vec![corrected], vec![2]);
    let target = tp(&raw_owner, 1600.0, 100000.0, 1.0).entropy_j_k();
    let solve = |o: &IdealGasEquilibriumV1| o.equilibrate_sp(Pressure::new(30000.0),
        &[(id("H2"), Amount::new(1.0))], target, controls(1.0), || false);
    assert!(matches!(solve(&raw_owner), Err(EquilibriumExpansionErrorV1::Invalid(_))));
    close(solve(&corrected_owner).unwrap().equilibrium().temperature().value(),
        1600.0 * 0.3_f64.powf(2.0 / 7.0), 2e-9);
}
