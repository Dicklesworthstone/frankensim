use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, MolarEnergyQuantityV1,
    NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReactionId, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

#[path = "tests/hydrogen.rs"]
mod hydrogen;
#[path = "tests/methane.rs"]
mod methane;
#[path = "tests/prepared.rs"]
mod prepared;
#[path = "tests/reduced.rs"]
mod reduced;
#[path = "tests/sdirk.rs"]
mod sdirk;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn source(
    name: &str,
    mass: f64,
    cp_r: f64,
    formation: f64,
    entropy_offset: f64,
) -> Nasa9StandardStateModelV1 {
    let coefficients = [
        0.0,
        0.0,
        cp_r,
        0.0,
        0.0,
        0.0,
        0.0,
        formation / R,
        entropy_offset,
    ];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0),
        sources: vec![],
        provenance: Provenance {
            source: "manufactured finite-rate reactor".into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("kinetic-reactor-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn mechanism(
    models: Vec<Nasa9StandardStateModelV1>,
    atoms: Vec<i128>,
    reactions: Vec<ElementaryReactionV1>,
) -> IdealGasKineticsV1 {
    let species: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            species.clone(),
            vec![atoms],
        )
        .unwrap(),
        ChargeVector::new(species.clone(), vec![0; species.len()]).unwrap(),
        models,
    )
    .unwrap();
    IdealGasKineticsV1::new(source, "manufactured-reactor", 1, reactions, 1e-12).unwrap()
}
fn reaction(
    name: &str,
    lhs: &[(&str, u32)],
    rhs: &[(&str, u32)],
    rate: f64,
) -> ElementaryReactionV1 {
    ElementaryReactionV1::new(
        ReactionId::new(name).unwrap(),
        lhs.iter().map(|(s, n)| (id(s), *n)).collect(),
        rhs.iter().map(|(s, n)| (id(s), *n)).collect(),
        ArrheniusRateV1::new(
            Temperature::new(500.0),
            ReactionProgressRateQuantityV1::new(rate),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap()
}
fn isomers(q: f64, rate: f64) -> IdealGasKineticsV1 {
    mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, -q, 0.0),
        ],
        vec![1, 1],
        vec![reaction("A-B", &[("A", 1)], &[("B", 1)], rate)],
    )
}
fn controls(relative: f64) -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        maximum_temperature_k: 2000.0,
        energy_tolerance_j: 1e-8,
        relative_tolerance: relative,
        absolute_amount_tolerance_mol: 1e-13,
        entropy_tolerance_j_k: 1e-10,
        initial_step_s: 0.01,
        maximum_step_s: 0.5,
        maximum_steps: 30000,
        ..IsochoricKineticsSolve::default()
    }
}
fn run(
    m: &IdealGasKineticsV1,
    t: f64,
    duration: f64,
    control: IsochoricKineticsSolve,
) -> IsochoricKineticsResult {
    advance_adiabatic(
        m,
        &[(id("A"), Amount::new(1.0))],
        Volume::new(0.01),
        Temperature::new(t),
        Time::new(duration),
        control,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_finite_elapsed_chemistry_matches_analytic_relaxation_and_refines() {
    let m = isomers(0.0, 3.0);
    let expected = 0.5 * (1.0 - (-0.6_f64).exp());
    let coarse = run(&m, 500.0, 0.1, controls(1e-3));
    let fine = run(&m, 500.0, 0.1, controls(1e-5));
    let coarse_error = (coarse.amounts()[1].value() - expected).abs();
    let fine_error = (fine.amounts()[1].value() - expected).abs();
    assert!(coarse_error < 0.004, "{coarse_error:e}");
    assert!(fine_error < 0.0004, "{fine_error:e}");
    assert!(
        fine_error < coarse_error * 0.2,
        "coarse={coarse_error:e},fine={fine_error:e}"
    );
    assert!(fine.amounts()[1].value() < 0.3); // Does not jump to equilibrium0.5.
    assert_eq!(fine.elapsed_time(), Time::new(0.1));
    assert!((fine.temperature().value() - 500.0).abs() < 1e-9);
    assert!(fine.balance().entropy_change_j_k > 0.0);
    assert!(fine.balance().energy_residual_j.abs() <= controls(1e-5).energy_tolerance_j);
    assert!(fine.balance().maximum_estimated_local_error <= 1.0);
    assert_eq!(fine.sources().len(), 2);
}

#[test]
fn g1_stiff_relaxation_accepts_macroscopic_time_without_explicit_stability_steps() {
    let m = isomers(0.0, 1e8);
    let result = run(&m, 500.0, 0.1, controls(1e-5));
    assert!((result.amounts()[0].value() - 0.5).abs() < 2e-6);
    assert!((result.amounts()[1].value() - 0.5).abs() < 2e-6);
    assert!(
        result.balance().accepted_steps < 100,
        "{:?}",
        result.balance()
    );
    assert_eq!(result.elapsed_time(), Time::new(0.1));
    assert!(result.amounts().iter().all(|n| n.value() >= 0.0));
}

#[test]
fn g1_coupled_three_species_reactor_matches_independent_linear_eigenmodes() {
    let rate = 5.0;
    let m = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, 0.0, 0.0),
            source("C", 0.01, 4.0, 0.0, 0.0),
        ],
        vec![1, 1, 1],
        vec![
            reaction("A-B", &[("A", 1)], &[("B", 1)], rate),
            reaction("B-C", &[("B", 1)], &[("C", 1)], rate),
        ],
    );
    let duration = 0.2;
    let result = run(&m, 500.0, duration, controls(1e-5));
    // The symmetric three-state generator has eigenvalues0,-k,-3k.
    let slow = (-rate * duration).exp();
    let fast = (-3.0 * rate * duration).exp();
    let expected = [
        1.0 / 3.0 + slow / 2.0 + fast / 6.0,
        (1.0 - fast) / 3.0,
        1.0 / 3.0 - slow / 2.0 + fast / 6.0,
    ];
    for (amount, expected) in result.amounts().iter().zip(expected) {
        assert!(
            (amount.value() - expected).abs() < 3e-4,
            "actual={}, expected={expected}",
            amount.value()
        );
    }
    assert!(result.balance().mass_residual_kg.abs() < 1e-15);
    assert!(result.balance().maximum_relative_element_residual < 1e-14);
    assert!(result.balance().entropy_change_j_k > 0.0);
}

#[test]
fn g1_exothermic_uv_reactor_conserves_energy_and_reaches_independent_equilibrium() {
    let q = 12000.0;
    let m = isomers(q, 20.0);
    let result = run(&m, 300.0, 5.0, controls(2e-5));
    let b = result.amounts()[1].value();
    let expected_temperature = 300.0 + q * b / (3.0 * R);
    assert!((result.temperature().value() - expected_temperature).abs() < 2e-8);
    let (mut lo, mut hi) = (300.0, 300.0 + q / (3.0 * R));
    for _ in 0..100 {
        let t = (lo + hi) * 0.5;
        let x = 1.0 / (1.0 + (-q / (R * t)).exp());
        if 3.0 * R * (t - 300.0) - q * x < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    let equilibrium = (lo + hi) * 0.5;
    assert!((result.temperature().value() - equilibrium).abs() < 2e-4);
    assert!(
        (result.pressure().value() / (R * 300.0 / 0.01) - result.temperature().value() / 300.0)
            .abs()
            < 1e-12
    );
    assert!(result.balance().mass_residual_kg.abs() < 1e-15);
    assert!(result.balance().maximum_relative_element_residual < 1e-14);
    assert!(
        result.balance().minimum_step_entropy_change_j_k >= -controls(2e-5).entropy_tolerance_j_k
    );
}

#[test]
fn g1_dissociation_advances_molecular_count_and_coupled_pressure_temperature() {
    let tr = 1000.0_f64;
    let m = mechanism(
        vec![
            source("A", 0.001, 2.5, 0.0, 2.5 * (1.0 - tr.ln())),
            source("A2", 0.002, 3.5, 0.0, 3.5 * (1.0 - tr.ln())),
        ],
        vec![1, 2],
        vec![reaction("A2-2A", &[("A2", 1)], &[("A", 2)], 50.0)],
    );
    let v = Volume::new(0.01);
    let feed = [(id("A2"), Amount::new(1.0))];
    let result = advance_adiabatic(
        &m,
        &feed,
        v,
        Temperature::new(tr),
        Time::new(2.0),
        controls(2e-5),
        || false,
    )
    .unwrap();
    let a = result.amounts()[0].value();
    let a2 = result.amounts()[1].value();
    assert!((a + 2.0 * a2 - 2.0).abs() < 2e-14);
    assert!(a > 0.0 && a2 > 0.0);
    assert!((result.temperature().value() - 2.5 * tr / (1.5 * a + 2.5 * a2)).abs() < 2e-8);
    assert!(
        (result.pressure().value() - (a + a2) * R * result.temperature().value() / v.value()).abs()
            < 1e-7
    );
    let equilibrium = m
        .thermodynamics()
        .equilibrate_uv(
            v,
            &feed,
            Temperature::new(tr),
            fs_thermochem::equilibrium::composition::volume::EquilibriumInternalEnergySolveV1 {
                volume: Default::default(),
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 2000.0,
                energy_tolerance_j: 1e-8,
                maximum_temperature_iterations: 80,
            },
            || false,
        )
        .unwrap();
    assert!(
        (result.temperature().value() - equilibrium.equilibrium().temperature().value()).abs()
            < 1e-4
    );
}

#[test]
fn g3_extensive_rescaling_and_explicit_energy_replay_the_same_physics() {
    let m = isomers(5000.0, 5.0);
    let base = run(&m, 500.0, 0.2, controls(2e-5));
    for scale in [0.125, 8.0] {
        let mut ctl = controls(2e-5);
        ctl.absolute_amount_tolerance_mol *= scale;
        ctl.energy_tolerance_j *= scale;
        ctl.entropy_tolerance_j_k *= scale;
        let result = advance_adiabatic_energy(
            &m,
            &[(id("A"), Amount::new(scale))],
            Energy::new(base.target_internal_energy().value() * scale),
            Volume::new(0.01 * scale),
            Temperature::new(500.0),
            Time::new(0.2),
            ctl,
            || false,
        )
        .unwrap();
        assert_eq!(
            result.temperature().value().to_bits(),
            base.temperature().value().to_bits()
        );
        assert_eq!(
            result.pressure().value().to_bits(),
            base.pressure().value().to_bits()
        );
        assert_eq!(
            result.amounts()[1].value() / scale,
            base.amounts()[1].value()
        );
        assert_eq!(
            result.balance().accepted_steps,
            base.balance().accepted_steps
        );
    }
}

#[test]
fn g0_one_newton_correction_zero_duration_and_failures_have_explicit_semantics() {
    let m = isomers(0.0, 3.0);
    let mut ctl = controls(1e-3);
    ctl.maximum_newton_iterations = 1;
    ctl.initial_step_s = 0.001;
    let result = run(&m, 500.0, 0.001, ctl);
    assert!(result.amounts()[1].value() > 0.0);
    let zero = run(&m, 500.0, 0.0, ctl);
    assert_eq!(zero.amounts(), &[Amount::new(1.0), Amount::new(0.0)]);
    assert_eq!(zero.balance().accepted_steps, 0);
    let feed = [(id("A"), Amount::new(1.0))];
    let m_before = m.clone();
    assert!(matches!(
        advance_adiabatic_energy(
            &m,
            &feed,
            Energy::new(1e9),
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(0.1),
            ctl,
            || false
        ),
        Err(IsochoricKineticsError::EnergyNotBracketed { .. })
    ));
    ctl.maximum_steps = 1;
    assert!(matches!(
        advance_adiabatic(
            &m,
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(0.1),
            ctl,
            || false
        ),
        Err(IsochoricKineticsError::WorkLimit { .. })
    ));
    assert_eq!(m, m_before);
    assert!(matches!(
        advance_adiabatic(
            &m,
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(0.1),
            ctl,
            || true
        ),
        Err(IsochoricKineticsError::Cancelled)
    ));
}

#[test]
fn g4_late_cancellation_and_extensive_overflow_publish_no_partial_state() {
    let m = isomers(0.0, 3.0);
    let feed = [(id("A"), Amount::new(1.0))];
    let mut count = 0;
    let first = advance_adiabatic(
        &m,
        &feed,
        Volume::new(0.01),
        Temperature::new(500.0),
        Time::new(0.01),
        controls(1e-3),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    let mut polls = 0;
    let stopped = advance_adiabatic(
        &m,
        &feed,
        Volume::new(0.01),
        Temperature::new(500.0),
        Time::new(0.01),
        controls(1e-3),
        || {
            polls += 1;
            polls == count
        },
    );
    assert_eq!(stopped, Err(IsochoricKineticsError::Cancelled));
    assert_eq!(first, run(&m, 500.0, 0.01, controls(1e-3)));
    let huge = mechanism(
        vec![
            source("A", 1000.0, 4.0, -3.0 * R * 500.0, 0.0),
            source("B", 1000.0, 4.0, -3.0 * R * 500.0, 0.0),
        ],
        vec![1, 1],
        vec![reaction("heavy-isomers", &[("A", 1)], &[("B", 1)], 1.0)],
    );
    assert!(matches!(
        advance_adiabatic_energy(
            &huge,
            &[(id("A"), Amount::new(1e308))],
            Energy::new(0.0),
            Volume::new(1e308),
            Temperature::new(500.0),
            Time::new(0.0),
            controls(1e-3),
            || false
        ),
        Err(IsochoricKineticsError::Invalid(_))
    ));
}
