use super::*;
use crate::{
    ElementId, ElementalReferenceIdV1, InitialStatePolicy, MolarEnergyQuantityV1, NASA9_LAW_ID_V1,
    NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, ReferenceEquationOfStateV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use std::collections::BTreeMap;

const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
const REFERENCE_T: f64 = 1000.0;

fn id(value: &str) -> SpeciesId {
    SpeciesId::new(value).unwrap()
}

// Manufactured standard states. Formation constants and entropy offsets are
// explicit; g_at_reference fixes an independently chosen dimensionless g/RT.
fn model(
    name: &str,
    mass: f64,
    cp_r: f64,
    formation_j_mol: f64,
    g_at_reference: f64,
) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[7] = formation_j_mol / R;
    coefficients[8] =
        cp_r + formation_j_mol / (R * REFERENCE_T) - cp_r * REFERENCE_T.ln() - g_at_reference;
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.to_string(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".to_string(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.to_string()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "manufactured equilibrium conformance".into(),
            license: "test fixture".into(),
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
            ElementalReferenceIdV1::new("equilibrium-test-reference").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn owner(
    models: Vec<Nasa9StandardStateModelV1>,
    elements: &[&str],
    rows: Vec<Vec<i128>>,
) -> IdealGasEquilibriumV1 {
    let species: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let matrix = ElementalMatrix::new(
        elements
            .iter()
            .map(|e| ElementId::new(*e).unwrap())
            .collect(),
        species.clone(),
        rows,
    )
    .unwrap();
    let charge = ChargeVector::new(species.clone(), vec![0; species.len()]).unwrap();
    IdealGasEquilibriumV1::new(matrix, charge, models).unwrap()
}
fn dissociation() -> IdealGasEquilibriumV1 {
    owner(
        vec![
            model("A", 0.001, 2.5, 0.0, 0.0),
            model("A2", 0.002, 3.5, 0.0, 0.0),
        ],
        &["H"],
        vec![vec![1, 2]],
    )
}
fn isomers(q: f64) -> IdealGasEquilibriumV1 {
    owner(
        vec![
            model("A", 0.01, 4.0, 0.0, 0.0),
            model("B", 0.01, 4.0, -q, -q / (R * REFERENCE_T)),
        ],
        &["H"],
        vec![vec![1, 1]],
    )
}
fn species<'a>(
    solution: &'a IdealGasEquilibriumSolutionV1,
    name: &str,
) -> &'a EquilibriumSpeciesAmountV1 {
    solution
        .species()
        .iter()
        .find(|s| s.species.as_str() == name)
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn hp_controls() -> EquilibriumEnthalpySolveV1 {
    EquilibriumEnthalpySolveV1 {
        equilibrium: EquilibriumSolveV1::default(),
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 2000.0,
        enthalpy_tolerance_j: 1e-6,
        maximum_temperature_iterations: 80,
    }
}

#[test]
fn g1_dissociation_matches_closed_form_across_six_decades_of_pressure() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    for pressure_ratio in [0.001_f64, 0.01, 1.0, 10.0, 1000.0, 1000000.0] {
        let solution = chemistry
            .equilibrate_tp(
                Temperature::new(REFERENCE_T),
                Pressure::new(100000.0 * pressure_ratio),
                &feed,
                EquilibriumSolveV1::default(),
                || false,
            )
            .unwrap();
        // A2 <=> 2A, Kp=1: alpha=sqrt(1/(1+4p/p0)).
        let alpha = (1.0 / (1.0 + 4.0 * pressure_ratio)).sqrt();
        close(species(&solution, "A").amount.value(), 2.0 * alpha, 2e-9);
        close(species(&solution, "A2").amount.value(), 1.0 - alpha, 2e-9);
        close(solution.mass().value(), 0.002, 1e-12);
        assert!(solution.residuals().maximum_relative_element_residual <= 1e-10);
        assert!(
            solution
                .residuals()
                .maximum_dimensionless_potential_residual
                <= 1e-9
        );
        for product in solution.species() {
            assert!(product.amount.value() > 0.0);
            assert_eq!(product.source.as_ref().unwrap().species(), &product.species);
        }
    }
}

#[test]
fn g0_multielement_equilibrium_retains_minor_element_inventory_and_mass_action() {
    let chemistry = owner(
        vec![
            model("H", 0.001, 2.5, 0.0, 0.0),
            model("O", 0.016, 2.5, 0.0, 0.0),
            model("HO", 0.017, 3.5, 0.0, -5.0),
        ],
        &["H", "O"],
        vec![vec![1, 0, 1], vec![0, 1, 1]],
    );
    for hydrogen in [2.0_f64, 1.0, 1e-8] {
        let feed = [
            (id("H"), Amount::new(hydrogen)),
            (id("O"), Amount::new(1.0)),
        ];
        let state = chemistry
            .equilibrate_tp(
                Temperature::new(REFERENCE_T),
                Pressure::new(100000.0),
                &feed,
                EquilibriumSolveV1::default(),
                || false,
            )
            .unwrap();
        let h = species(&state, "H");
        let o = species(&state, "O");
        let ho = species(&state, "HO");
        assert!(((h.amount.value() + ho.amount.value()) / hydrogen - 1.0).abs() < 1.1e-10);
        close(o.amount.value() + ho.amount.value(), 1.0, 1.1e-10);
        let log_k = ho.log_mole_fraction.unwrap()
            - h.log_mole_fraction.unwrap()
            - o.log_mole_fraction.unwrap();
        close(log_k, 5.0, 1e-9);
        assert!(state.relative_mass_residual().abs() <= 1e-10);
    }
}

#[test]
fn g0_trace_underflow_is_not_structural_species_absence() {
    let chemistry = owner(
        vec![
            model("A", 0.01, 4.0, 0.0, 0.0),
            model("B", 0.01, 4.0, 0.0, -1000.0),
        ],
        &["H"],
        vec![vec![1, 1]],
    );
    let state = chemistry
        .equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &[(id("A"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap();
    let trace = species(&state, "A");
    assert_eq!(trace.mole_fraction, 0.0);
    close(trace.log_mole_fraction.unwrap(), -1000.0, 1e-9);
    assert!(trace.source.is_some());
    assert!(
        state.composition_iterations() < 32,
        "large downward log steps need no arbitrary trace floor"
    );
}

#[test]
fn g0_absent_elements_and_redundant_elemental_rows_have_explicit_semantics() {
    let chemistry = owner(
        vec![
            model("H", 0.001, 2.5, 0.0, 0.0),
            model("H2", 0.002, 3.5, 0.0, 0.0),
            model("O2", 0.032, 3.5, 0.0, 0.0),
        ],
        &["H", "O"],
        vec![vec![1, 2, 0], vec![0, 0, 2]],
    );
    let state = chemistry
        .equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &[(id("H2"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap();
    let absent = species(&state, "O2");
    assert_eq!(absent.amount.value(), 0.0);
    assert!(absent.log_mole_fraction.is_none());
    assert!(absent.source.is_none());
    let water = owner(
        vec![model("H2O", 0.018, 4.0, 0.0, 0.0)],
        &["H", "O"],
        vec![vec![2], vec![1]],
    );
    let pure = water
        .equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &[(id("H2O"), Amount::new(2.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap();
    close(pure.total_amount().value(), 2.0, 2e-10);
    close(
        pure.equilibrium_heat_capacity_p_j_k().unwrap(),
        8.0 * R,
        2e-9,
    );
}

#[test]
fn g0_strong_multielement_products_use_full_kkt_when_reduced_gram_loses_rank() {
    let chemistry = owner(
        vec![
            model("H", 0.001, 2.5, 0.0, 0.0),
            model("O", 0.016, 2.5, 0.0, 0.0),
            model("HO", 0.017, 3.5, 0.0, -100.0),
        ],
        &["H", "O"],
        vec![vec![1, 0, 1], vec![0, 1, 1]],
    );
    let state = chemistry
        .equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &[(id("H"), Amount::new(1.0)), (id("O"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap();
    close(species(&state, "HO").amount.value(), 1.0, 2e-10);
    let log_q = species(&state, "HO").log_mole_fraction.unwrap()
        - species(&state, "H").log_mole_fraction.unwrap()
        - species(&state, "O").log_mole_fraction.unwrap();
    close(log_q, 100.0, 2e-9);
    assert!(state.residuals().maximum_relative_element_residual <= 1e-10);
    // Residual bounds do not certify relative accuracy of sub-tolerance traces.
}

#[test]
fn g1_adiabatic_hp_matches_independent_isomer_energy_root() {
    let q = 24000.0;
    let chemistry = isomers(q);
    let inlet = 300.0;
    let state = chemistry
        .equilibrate_hp(
            Pressure::new(400000.0),
            &[(id("A"), Amount::new(1.0))],
            Temperature::new(inlet),
            hp_controls(),
            || false,
        )
        .unwrap();
    // Independent scalar oracle: x_B=1/(1+exp(-Q/(RT))) and
    // 4R(T-T_in)-Q*x_B=0, with no additional heat release in H.
    let (mut low, mut high) = (300.0_f64, 1500.0_f64);
    for _ in 0..80 {
        let t = (low + high) * 0.5;
        let f = 1.0 / (1.0 + (-q / (R * t)).exp());
        if 4.0 * R * (t - inlet) - q * f < 0.0 {
            low = t;
        } else {
            high = t;
        }
    }
    close(state.temperature().value(), (low + high) * 0.5, 2e-9);
    close(state.enthalpy().value(), 4.0 * R * inlet, 2e-10);
    assert!(state.enthalpy_residual().unwrap().value().abs() <= hp_controls().enthalpy_tolerance_j);
    assert!(state.temperature_iterations() < 40);
}

#[test]
fn g1_hp_core_combines_distinct_inlet_temperatures_in_one_enthalpy_target() {
    let chemistry = isomers(24000.0);
    let a = [(id("A"), Amount::new(0.8))];
    let b = [(id("B"), Amount::new(0.2))];
    let ha = chemistry
        .feed_enthalpy(Temperature::new(300.0), &a, || false)
        .unwrap();
    let hb = chemistry
        .feed_enthalpy(Temperature::new(700.0), &b, || false)
        .unwrap();
    let feed = [a[0].clone(), b[0].clone()];
    let target = Energy::new(ha.value() + hb.value());
    let mixed = chemistry
        .equilibrate_hp_enthalpy(
            Pressure::new(100000.0),
            &feed,
            target,
            hp_controls(),
            || false,
        )
        .unwrap();
    // Equal species heat capacities make the common-temperature reference
    // exactly 0.8*300+0.2*700 K, retaining B's formation energy in both routes.
    let common = chemistry
        .equilibrate_hp(
            Pressure::new(100000.0),
            &feed,
            Temperature::new(380.0),
            hp_controls(),
            || false,
        )
        .unwrap();
    close(
        mixed.temperature().value(),
        common.temperature().value(),
        2e-10,
    );
    close(mixed.enthalpy().value(), target.value(), 2e-10);
}

#[test]
fn g1_equilibrium_heat_capacity_contains_composition_response() {
    let q = 24000.0;
    let chemistry = isomers(q);
    let feed = [(id("A"), Amount::new(1.0))];
    let t = 1000.0;
    let solve = |temperature| {
        chemistry
            .equilibrate_tp(
                Temperature::new(temperature),
                Pressure::new(100000.0),
                &feed,
                EquilibriumSolveV1::default(),
                || false,
            )
            .unwrap()
    };
    let state = solve(t);
    let f = 1.0 / (1.0 + (-q / (R * t)).exp());
    let expected = 4.0 * R + f * (1.0 - f) * q * q / (R * t * t);
    close(
        state.equilibrium_heat_capacity_p_j_k().unwrap(),
        expected,
        2e-9,
    );
    assert!(state.equilibrium_heat_capacity_p_j_k().unwrap() > 4.0 * R);
    let derivative = (solve(t + 0.1).enthalpy().value() - solve(t - 0.1).enthalpy().value()) / 0.2;
    close(
        state.equilibrium_heat_capacity_p_j_k().unwrap(),
        derivative,
        2e-7,
    );
}

#[test]
fn g3_feed_scale_order_and_common_element_reference_preserve_equilibrium() {
    let base = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let solve = |owner: &IdealGasEquilibriumV1, feed: &[(SpeciesId, Amount)]| {
        owner
            .equilibrate_tp(
                Temperature::new(REFERENCE_T),
                Pressure::new(200000.0),
                feed,
                EquilibriumSolveV1::default(),
                || false,
            )
            .unwrap()
    };
    let ordinary = solve(&base, &feed);
    let shifted = owner(
        vec![
            model("A2", 0.002, 3.5, 2.0 * R * 100000.0, 200.0),
            model("A", 0.001, 2.5, R * 100000.0, 100.0),
        ],
        &["H"],
        vec![vec![2, 1]],
    );
    let altered = solve(
        &shifted,
        &[(id("A"), Amount::new(0.0)), (id("A2"), Amount::new(1e7))],
    );
    for name in ["A", "A2"] {
        close(
            species(&ordinary, name).mole_fraction,
            species(&altered, name).mole_fraction,
            2e-9,
        );
        close(
            species(&altered, name).amount.value() / 1e7,
            species(&ordinary, name).amount.value(),
            2e-9,
        );
    }
}

#[test]
fn g0_source_domain_mass_inconsistency_and_work_limits_refuse_without_mutation() {
    let chemistry = dissociation();
    let original = chemistry.clone();
    let feed = [(id("A2"), Amount::new(1.0))];
    assert!(matches!(
        chemistry.equilibrate_tp(
            Temperature::new(100.0),
            Pressure::new(100000.0),
            &feed,
            EquilibriumSolveV1::default(),
            || false
        ),
        Err(CompositionEquilibriumErrorV1::SpeciesEvaluation { .. })
    ));
    let mut budget = EquilibriumSolveV1::default();
    budget.maximum_iterations = 1;
    assert!(matches!(
        chemistry.equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &feed,
            budget,
            || false
        ),
        Err(CompositionEquilibriumErrorV1::NoConvergence { .. })
    ));
    assert_eq!(chemistry, original);
    let wrong = owner(
        vec![
            model("A", 0.01, 4.0, 0.0, 0.0),
            model("B", 0.02, 4.0, 0.0, 0.0),
        ],
        &["H"],
        vec![vec![1, 1]],
    );
    assert!(matches!(
        wrong.equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &[(id("A"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false
        ),
        Err(CompositionEquilibriumErrorV1::MassBalance { .. })
    ));
    assert!(matches!(
        isomers(24000.0).equilibrate_hp_enthalpy(
            Pressure::new(100000.0),
            &[(id("A"), Amount::new(1.0))],
            Energy::new(1e9),
            hp_controls(),
            || false
        ),
        Err(CompositionEquilibriumErrorV1::EnthalpyNotBracketed { .. })
    ));
}

#[test]
fn g4_late_cancellation_never_publishes_a_partial_composition() {
    let chemistry = dissociation();
    let original = chemistry.clone();
    let feed = [(id("A2"), Amount::new(1.0))];
    let mut calls = 0;
    chemistry
        .equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &feed,
            EquilibriumSolveV1::default(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        chemistry.equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &feed,
            EquilibriumSolveV1::default(),
            || {
                count += 1;
                count == calls
            }
        ),
        Err(CompositionEquilibriumErrorV1::Cancelled)
    ));
    assert_eq!(chemistry, original);
}

#[test]
fn g0_neutrality_identity_and_invalid_feed_are_admitted_explicitly() {
    let chemistry = dissociation();
    let charges =
        ChargeVector::new(chemistry.elemental_matrix().species().to_vec(), vec![1, 0]).unwrap();
    assert!(matches!(
        IdealGasEquilibriumV1::new(
            chemistry.elemental_matrix().clone(),
            charges,
            chemistry.models().to_vec()
        ),
        Err(CompositionEquilibriumErrorV1::Invalid(_))
    ));
    let bad = [(id("A2"), Amount::new(1.0)), (id("A2"), Amount::new(0.0))];
    assert!(matches!(
        chemistry.equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(100000.0),
            &bad,
            EquilibriumSolveV1::default(),
            || false
        ),
        Err(CompositionEquilibriumErrorV1::DuplicateFeedSpecies(_))
    ));
    assert!(matches!(
        chemistry.equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(0.0),
            &[(id("A2"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false
        ),
        Err(CompositionEquilibriumErrorV1::Invalid(_))
    ));
    let _typed_total: Energy = Energy::new(0.0);
    let _typed_molar: MolarEnergyQuantityV1 = MolarEnergyQuantityV1::new(0.0);
}
