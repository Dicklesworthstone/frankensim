use super::*;
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};

fn chart(phase: SolidLiquidPhase, values: &[(f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    let liquid = match phase {
        SolidLiquidPhase::Solid => 0.0,
        SolidLiquidPhase::Liquid => 1.0,
        SolidLiquidPhase::SolidLiquid => panic!("fixture requires a pure caloric branch"),
    };
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([71; 32]),
        phase,
        values
            .iter()
            .map(
                |&(temperature_k, specific_enthalpy_j_kg)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg,
                    temperature_k,
                    liquid_mass_fraction: liquid,
                    bulk_density_kg_m3: if liquid == 0.0 { 1000.0 } else { 800.0 },
                },
            )
            .collect(),
    )
    .unwrap()
}

fn branch(phase: SolidLiquidPhase, cp: f64, h0: f64) -> EquilibriumEnthalpyPhaseCurve {
    chart(phase, &[(280.0, h0), (600.0, h0 + 320.0 * cp)])
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_unequal_original_capacities_define_phase_entropy_and_gibbs_gaps() {
    let parent = branch(SolidLiquidPhase::Liquid, 2000.0, 100_000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1000.0, 0.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 600.0]).unwrap();
    let entropy = CrystallizationEntropy::new(&pair, 500.0).unwrap();
    assert_eq!(entropy.caloric_pair().identity(), pair.identity());
    assert_eq!(
        entropy.caloric_pair().parent().source_receipts(),
        pair.parent().source_receipts()
    );
    assert_eq!(
        entropy.caloric_pair().crystal().source_receipts(),
        pair.crystal().source_receipts()
    );
    assert_eq!(entropy.coexistence_temperature_k(), 500.0);
    assert_eq!(entropy.phase_gibbs_gap_j_kg(500.0).unwrap(), 0.0);
    for temperature in [280.0_f64, 300.0, 400.0, 500.0, 600.0] {
        let gap_s = 640.0 + 1000.0 * (temperature / 500.0).ln();
        let gap_h = 100_000.0 + 1000.0 * (temperature - 280.0);
        close(
            entropy.phase_entropy_gap_j_kg_k(temperature).unwrap(),
            gap_s,
            5e-13,
        );
        close(
            entropy.phase_gibbs_gap_j_kg(temperature).unwrap(),
            gap_h - temperature * gap_s,
            2e-10,
        );
    }
    assert!(entropy.phase_gibbs_gap_j_kg(400.0).unwrap() > 0.0);
    assert!(entropy.phase_gibbs_gap_j_kg(600.0).unwrap() < 0.0);
    // Independent state function. The phase densities differ, so replacing X
    // by crystal VOLUME fraction would fail this caloric oracle.
    let primitive = |t: f64, x: f64| (2000.0 - 1000.0 * x) * (t / 500.0).ln() - 640.0 * x;
    close(
        entropy
            .specific_entropy_change_j_kg_k(300.0, 0.2, 420.0, 0.65)
            .unwrap(),
        primitive(420.0, 0.65) - primitive(300.0, 0.2),
        1e-12,
    );
    // A warm amorphous state still uses its original single-phase entropy.
    close(
        entropy
            .specific_entropy_change_j_kg_k(500.0, 0.0, 550.0, 0.0)
            .unwrap(),
        2000.0 * 1.1_f64.ln(),
        5e-13,
    );
}

#[test]
fn g1_entropy_retains_unequal_original_temperature_breakpoints() {
    let parent = chart(
        SolidLiquidPhase::Liquid,
        &[
            (280.0, 120_000.0),
            (350.0, 246_000.0),
            (500.0, 576_000.0),
            (600.0, 816_000.0),
        ],
    );
    let crystal = chart(
        SolidLiquidPhase::Solid,
        &[
            (280.0, 0.0),
            (400.0, 108_000.0),
            (500.0, 218_000.0),
            (600.0, 348_000.0),
        ],
    );
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 600.0]).unwrap();
    let entropy = CrystallizationEntropy::new(&pair, 500.0).unwrap();
    let parent_initial = -1800.0 * (350.0_f64 / 330.0).ln() - 2200.0 * (500.0_f64 / 350.0).ln();
    let crystal_initial =
        -900.0 * (400.0_f64 / 330.0).ln() - 1100.0 * (500.0_f64 / 400.0).ln() - 716.0;
    let parent_final = 2400.0 * 1.1_f64.ln();
    let crystal_final = 1300.0 * 1.1_f64.ln() - 716.0;
    let expected =
        0.9 * parent_final + 0.1 * crystal_final - (0.75 * parent_initial + 0.25 * crystal_initial);
    close(
        entropy
            .specific_entropy_change_j_kg_k(330.0, 0.25, 550.0, 0.1)
            .unwrap(),
        expected,
        2e-12,
    );
}

#[test]
fn g1_accepted_plateau_remelting_has_the_original_latent_entropy() {
    let parent = branch(SolidLiquidPhase::Liquid, 1000.0, 100_000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1000.0, 0.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 600.0]).unwrap();
    let entropy = CrystallizationEntropy::new(&pair, 500.0).unwrap();
    let rate = |t: f64| Ok(0.002 * (500.0 - t).max(0.0));
    let law = pair
        .law(&rate, ContentHash([72; 32]), 2.0, 0.8, 500.0)
        .unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(500.0, 0.8).unwrap();
    let policy = law.equivalent_age_remelting().unwrap();
    let partial = law
        .advance_specific_enthalpy_with_remelting(initial, 310_000.0, 1.0, solve, policy, || false)
        .unwrap();
    assert_eq!(partial.state.temperature_k(), 500.0);
    close(partial.state.crystalline_mass_fraction(), 0.1, 1e-13);
    close(
        entropy
            .step_entropy_change_j_kg_k(&law, &partial, solve)
            .unwrap(),
        partial.specific_heat_j_kg / 500.0,
        3e-11,
    );
    let melted = law
        .advance_specific_enthalpy_with_remelting(
            partial.state,
            370_000.0,
            1.0,
            solve,
            policy,
            || false,
        )
        .unwrap();
    assert_eq!(melted.state.crystalline_mass_fraction(), 0.0);
    close(melted.state.temperature_k(), 550.0, 1e-9);
    close(
        entropy
            .step_entropy_change_j_kg_k(&law, &melted, solve)
            .unwrap(),
        1000.0 * 1.1_f64.ln() + 200.0 * partial.state.crystalline_mass_fraction(),
        3e-10,
    );
}

#[test]
fn g1_adiabatic_recalescence_and_isothermal_growth_keep_distinct_entropy_accounts() {
    let parent = branch(SolidLiquidPhase::Liquid, 1000.0, 100_000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1000.0, 0.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 600.0]).unwrap();
    let entropy = CrystallizationEntropy::new(&pair, 500.0).unwrap();
    let rate = |t: f64| Ok(0.002 * (500.0 - t).max(0.0));
    let law = pair
        .law(&rate, ContentHash([72; 32]), 2.0, 0.8, 500.0)
        .unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(400.0, 0.0).unwrap();
    let adiabatic = law
        .advance_specific_enthalpy(
            initial,
            initial.specific_enthalpy_j_kg(),
            0.5,
            solve,
            || false,
        )
        .unwrap();
    assert_eq!(adiabatic.specific_heat_j_kg, 0.0);
    assert!(adiabatic.state.temperature_k() > initial.temperature_k());
    let entropy_change = entropy
        .step_entropy_change_j_kg_k(&law, &adiabatic, solve)
        .unwrap();
    let x = adiabatic.state.crystalline_mass_fraction();
    close(
        entropy_change,
        1000.0 * (adiabatic.state.temperature_k() / 400.0).ln() - 200.0 * x,
        5e-13,
    );
    assert!(entropy_change > 0.0);

    let isothermal = law.advance_isothermal(initial, 0.5).unwrap();
    let material_entropy = entropy
        .step_entropy_change_j_kg_k(&law, &isothermal, solve)
        .unwrap();
    assert!(material_entropy < 0.0);
    let reservoir_entropy = -isothermal.specific_heat_j_kg / 400.0;
    let produced = material_entropy + reservoir_entropy;
    close(
        produced,
        50.0 * isothermal.state.crystalline_mass_fraction(),
        1e-12,
    );
    assert!(produced > 0.0);
}

#[test]
fn g3_entropy_is_path_additive_reversible_and_independent_of_common_enthalpy_datum() {
    let evaluate = |datum: f64| {
        let parent = branch(SolidLiquidPhase::Liquid, 1024.0, datum + 65_536.0);
        let crystal = branch(SolidLiquidPhase::Solid, 1024.0, datum);
        let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 600.0]).unwrap();
        let entropy = CrystallizationEntropy::new(&pair, 500.0).unwrap();
        let change = |a: (f64, f64), b: (f64, f64)| {
            entropy
                .specific_entropy_change_j_kg_k(a.0, a.1, b.0, b.1)
                .unwrap()
        };
        let (a, b, c) = ((300.0, 0.1), (410.0, 0.4), (550.0, 0.0));
        let ab = change(a, b);
        let bc = change(b, c);
        let ac = change(a, c);
        assert_eq!(change(b, a).to_bits(), (-ab).to_bits());
        assert_eq!(change(a, a), 0.0);
        close(ab + bc, ac, 5e-13);
        (ab, bc, ac, entropy.phase_gibbs_gap_j_kg(410.0).unwrap())
    };
    assert_eq!(evaluate(0.0), evaluate(1.0e12));
}

#[test]
fn g0_entropy_refuses_uncovered_coexistence_invalid_states_and_wrong_bindings() {
    let parent = branch(SolidLiquidPhase::Liquid, 1000.0, 100_000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1000.0, 0.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 600.0]).unwrap();
    for tm in [f64::NAN, f64::INFINITY, 0.0, 279.0, 601.0] {
        assert!(CrystallizationEntropy::new(&pair, tm).is_err());
    }
    let entropy = CrystallizationEntropy::new(&pair, 500.0).unwrap();
    for (t, x) in [
        (279.0, 0.0),
        (601.0, 0.0),
        (f64::NAN, 0.0),
        (400.0, -0.1),
        (400.0, 1.01),
        (400.0, f64::NAN),
    ] {
        assert!(
            entropy
                .specific_entropy_change_j_kg_k(400.0, 0.0, t, x)
                .is_err()
        );
        assert!(
            entropy
                .specific_entropy_change_j_kg_k(t, x, 400.0, 0.0)
                .is_err()
        );
    }
    let rate = |_: f64| Ok(0.0);
    let law = pair
        .law(&rate, ContentHash([72; 32]), 2.0, 0.8, 500.0)
        .unwrap();
    let state = law.state_at_temperature(400.0, 0.5).unwrap();
    let solve = CrystallizationSolve::default();
    let wrong_tm = CrystallizationEntropy::new(&pair, 510.0).unwrap();
    assert!(matches!(
        wrong_tm.state_entropy_change_j_kg_k(&law, state, state, solve),
        Err(CrystallizationEntropyError::LawMismatch)
    ));
    let wrong_crystal = branch(SolidLiquidPhase::Solid, 1000.0, 1.0);
    let wrong_pair =
        CrystallizationCaloricPair::new(&parent, &wrong_crystal, [280.0, 600.0]).unwrap();
    let wrong_sources = CrystallizationEntropy::new(&wrong_pair, 500.0).unwrap();
    assert!(matches!(
        wrong_sources.state_entropy_change_j_kg_k(&law, state, state, solve),
        Err(CrystallizationEntropyError::LawMismatch)
    ));
    let narrow_pair = CrystallizationCaloricPair::new(&parent, &crystal, [300.0, 550.0]).unwrap();
    let narrow = CrystallizationEntropy::new(&narrow_pair, 500.0).unwrap();
    assert!(matches!(
        narrow.state_entropy_change_j_kg_k(&law, state, state, solve),
        Err(CrystallizationEntropyError::LawMismatch)
    ));
    let other = pair
        .law(&rate, ContentHash([73; 32]), 2.0, 0.8, 500.0)
        .unwrap();
    assert!(matches!(
        entropy.state_entropy_change_j_kg_k(&other, state, state, solve),
        Err(CrystallizationEntropyError::Kinetics(
            CrystallizationError::StateMismatch
        ))
    ));
    assert_ne!(entropy.identity(), wrong_tm.identity());
    assert_ne!(entropy.identity(), wrong_sources.identity());
}
