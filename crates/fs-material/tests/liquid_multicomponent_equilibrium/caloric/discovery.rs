//! Original material receipts after finite discovered phase-count transitions.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::ResolvedMulticomponentLiquidEquilibriumCaloricState;
use fs_thermochem::solution::multicomponent::caloric::equilibrium::discovery::MulticomponentDiscoveredCaloricSolveV1;
use fs_thermochem::solution::multicomponent::discovery::MulticomponentPhaseDiscoveryControlV1;

fn discovered_controls() -> MulticomponentDiscoveredCaloricSolveV1 {
    let mut thermal = controls().thermal;
    thermal.temperature_bracket = [Temperature::new(380.0), Temperature::new(595.0)];
    MulticomponentDiscoveredCaloricSolveV1 {
        thermal,
        response: Default::default(),
        discovery: MulticomponentPhaseDiscoveryControlV1 {
            equilibrium: MulticomponentEquilibriumControlV1 {
                chemical_potential_tolerance: 2e-13,
                ..Default::default()
            },
            ..Default::default()
        },
        maximum_candidate_solves: 16_384,
    }
}

fn separated(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
    temperature: f64,
    amounts: [f64; 3],
) -> ResolvedMulticomponentLiquidEquilibriumCaloricState {
    caloric
        .equilibrium_caloric_state(
            caloric
                .liquid_law()
                .resolve_equilibrium(
                    candidate(caloric.model().thermodynamics(), temperature, amounts, true),
                    || false,
                )
                .unwrap(),
            Default::default(),
            || false,
        )
        .unwrap()
}

fn homogeneous_oracle(temperature: f64) -> [f64; 2] {
    let d = temperature - 450.0;
    [
        840.0 * (135.0 + 1.2 * d + 0.003 * d * d) + 10.0 * GAS_CONSTANT * 1600.0 * 0.31,
        840.0 * (-0.6 * 1.5_f64.ln() + 0.6 - 1.5 * (temperature / 450.0).ln() + 0.006 * d)
            - 10.0 * GAS_CONSTANT * (0.2 * 0.2_f64.ln() + 0.3 * 0.3_f64.ln() + 0.5 * 0.5_f64.ln()),
    ]
}

fn assert_original_sources(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
    final_state: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
    card: &MaterialCard,
) {
    assert_eq!(final_state.sources(), caloric.sources());
    assert_eq!(
        final_state.caloric_state().candidate(),
        final_state.equilibrium().equilibrium_candidate()
    );
    let mut volume = 0.0;
    for phase in final_state.equilibrium().phases() {
        let liquid = phase.liquid();
        assert_eq!(liquid.material().card_identity(), card.content_hash());
        assert_eq!(
            liquid.material().state_coordinate("T"),
            Some(final_state.temperature().value())
        );
        assert_eq!(liquid.material().state_coordinate("p"), Some(100_000.0));
        near(
            liquid.density_kg_m3(),
            source_value(phase.mole_fractions(), DENSITIES),
        );
        near(
            liquid.thermal_conductivity_w_m_k(),
            source_value(phase.mole_fractions(), CONDUCTIVITIES),
        );
        volume += phase.mass_kg() / liquid.density_kg_m3();
    }
    near(final_state.equilibrium().total_volume_m3(), volume);
    near(
        final_state.specific_heat_cp_j_kg_k(),
        final_state.caloric_state().isobaric_heat_capacity_j_per_k()
            / final_state.equilibrium().total_mass_kg(),
    );
}

#[test]
fn g1_three_to_one_heating_and_one_to_three_cooling_refresh_actual_sources() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    // The equimolar feed has a homogeneous equilibrium at 590 K, above its
    // chi=4 ln(2) transition and inside the original 600 K material/Cp support.
    // The unequal AMOUNTS fixture remains spinodal there and cannot serve as
    // an independently admitted one-phase target.
    let amounts = [1.0; 3];
    let initial = separated(&caloric, 400.0, amounts);
    let retained = initial.clone();
    let (r, q) = phase_composition(400.0);
    let original_h = 220.0 * 80.0 + 3.0 * GAS_CONSTANT * 1600.0 * (2.0 * r * q + q * q);
    let original_s = 220.0 * (0.4 - 0.6 * (4.0_f64 / 3.0).ln())
        - 3.0 * GAS_CONSTANT * (r * r.ln() + 2.0 * q * q.ln());
    let target_h = 220.0 * (135.0 + 1.2 * 140.0 + 0.003 * 140.0 * 140.0) + GAS_CONSTANT * 1600.0;
    let target_s = 220.0
        * (-0.6 * 1.5_f64.ln() + 0.6 - 1.5 * (590.0_f64 / 450.0).ln() + 0.006 * 140.0)
        + 3.0 * GAS_CONSTANT * 3.0_f64.ln();
    let heat = target_h - original_h;
    let entropy = target_s - original_s;
    let heated = caloric
        .add_heat_discovered(&initial, heat, discovered_controls(), || false)
        .unwrap();
    assert_eq!(
        heated.final_state().equilibrium().phases().len(),
        1,
        "T={}, Q={}, dH={}, H_initial={}, H_final={}, residual={}",
        heated.final_state().temperature().value(),
        heat,
        heated.balance().balance().enthalpy_change_j(),
        initial.caloric_state().enthalpy_j(),
        heated.final_state().caloric_state().enthalpy_j(),
        heated.balance().balance().residual_j()
    );
    assert_eq!(
        heated.balance().discovery_report().selected_phase_count(),
        1
    );
    assert!((heated.final_state().temperature().value() - 590.0).abs() < 2e-7);
    assert!((heated.balance().balance().enthalpy_change_j() - heat).abs() < 1e-6);
    assert!((heated.entropy_change_j_per_k() - entropy).abs() < 1e-7);
    assert_original_sources(&caloric, heated.final_state(), &card);
    near(
        heated.final_state().equilibrium().total_mass_kg(),
        initial.equilibrium().total_mass_kg(),
    );
    assert_eq!(
        heated.final_state().caloric_state().component_moles(),
        amounts
    );

    let cooled = caloric
        .add_heat_discovered(heated.final_state(), -heat, discovered_controls(), || false)
        .unwrap();
    assert_eq!(cooled.final_state().equilibrium().phases().len(), 3);
    assert_eq!(
        cooled.balance().discovery_report().selected_phase_count(),
        3
    );
    assert!((cooled.final_state().temperature().value() - 400.0).abs() < 2e-7);
    assert!((cooled.entropy_change_j_per_k() + entropy).abs() < 1e-7);
    assert_original_sources(&caloric, cooled.final_state(), &card);
    near(
        cooled.final_state().equilibrium().total_mass_kg(),
        initial.equilibrium().total_mass_kg(),
    );
    assert_eq!(initial, retained);
}

#[test]
fn g1_unequal_temperature_discovered_mixing_retains_initial_phase_storage() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let left = separated(&caloric, 400.0, AMOUNTS);
    let right = caloric
        .homogeneous(Temperature::new(590.0), &[1.0; 3], || false)
        .unwrap();
    let result = caloric
        .mix_discovered(&left, &right, 0.0, discovered_controls(), || false)
        .unwrap();
    assert_original_sources(&caloric, result.final_state(), &card);
    assert_eq!(
        result.final_state().caloric_state().component_moles(),
        &[3.0, 4.0, 6.0]
    );
    let energy = left.caloric_state().enthalpy_j() + right.caloric_state().enthalpy_j();
    assert!((result.final_state().caloric_state().enthalpy_j() - energy).abs() < 1e-6);
    assert!(result.entropy_change_j_per_k() > 0.0);
    near(
        result.final_state().equilibrium().total_mass_kg(),
        3.0 * MASSES[0] + 4.0 * MASSES[1] + 6.0 * MASSES[2],
    );
}

#[test]
fn g3_newly_discovered_phase_outside_original_material_chart_refuses_atomically() {
    let card = card(0.9);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = caloric
        .homogeneous(Temperature::new(590.0), &AMOUNTS, || false)
        .unwrap();
    let retained = initial.clone();
    let heat = enthalpy_oracle(400.0) - homogeneous_oracle(590.0)[0];
    let numeric = caloric
        .model()
        .add_heat_discovered(initial.caloric_state(), heat, discovered_controls(), || {
            false
        })
        .unwrap();
    assert_eq!(numeric.final_state().candidate().phases().len(), 3);
    assert!(matches!(
        caloric.add_heat_discovered(&initial, heat, discovered_controls(), || false),
        Err(MulticomponentLiquidCaloricError::Material(
            MulticomponentLiquidError::Source(_)
        ))
    ));
    assert_eq!(initial, retained);
}

#[test]
fn g3_foreign_source_binding_and_last_poll_cancellation_publish_no_endpoint() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = caloric
        .homogeneous(Temperature::new(590.0), &AMOUNTS, || false)
        .unwrap();
    let retained = initial.clone();
    let foreign = MulticomponentLiquidCaloricLaw::try_new(
        &liquid,
        (0..3)
            .map(|i| (species(i), cp_source(i, 100_000.0, 310.0)))
            .collect(),
        "p",
    )
    .unwrap();
    assert!(matches!(
        foreign.add_heat_discovered(&initial, 0.0, discovered_controls(), || false),
        Err(MulticomponentLiquidCaloricError::Invalid(_))
    ));
    let mut controls = discovered_controls();
    controls.thermal.temperature_bracket = [Temperature::new(585.0), Temperature::new(595.0)];
    let mut complete_polls = 0;
    caloric
        .add_heat_discovered(&initial, 0.0, controls, || {
            complete_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert_eq!(
        caloric.add_heat_discovered(&initial, 0.0, controls, || {
            polls += 1;
            polls == complete_polls
        }),
        Err(MulticomponentLiquidCaloricError::Cancelled)
    );
    assert_eq!(polls, complete_polls);
    assert_eq!(initial, retained);
}
