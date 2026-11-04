//! G0/G3: finite multicomponent inventories, independent analytic kinetics,
//! shared calorics, permutation/subdivision and atomic source refusal.
use super::*;
use crate::cavitation::dissolved_gas::mixture::advance_mixture;

fn species_law(
    catalog: &ThermochemicalGasModel,
    name: &str,
    h: f64,
    d: f64,
    heat_j_mol: f64,
    pressures: [f64; 2],
) -> HenryDissolution {
    let species = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|s| s.species() == &id(name))
        .unwrap();
    let mut caloric = calorics(catalog);
    caloric.enthalpy_offset_j_kg = catalog
        .species_specific_enthalpy_j_kg(&id(name), T)
        .unwrap()
        - caloric
            .partial_enthalpy
            .state_at_temperature(T)
            .unwrap()
            .specific_enthalpy_j_kg()
        - heat_j_mol / species.molar_mass().value();
    HenryDissolution::new(
        catalog,
        id(name),
        HenrySource {
            material: card(
                &format!("{name} in the same declared solvent"),
                vec![
                    claim(
                        HENRY_SOLUBILITY_PROPERTY,
                        HENRY_SOLUBILITY_DIMS,
                        h,
                        pressures,
                    ),
                    claim(
                        DISSOLVED_DIFFUSIVITY_PROPERTY,
                        DIFFUSIVITY_DIMS,
                        d,
                        pressures,
                    ),
                ],
            ),
            point: context(T),
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        caloric,
    )
    .unwrap()
}

fn laws(catalog: &ThermochemicalGasModel) -> [HenryDissolution; 2] {
    [
        species_law(catalog, "A", H, D, 6000.0, [1000.0, 1e7]),
        species_law(catalog, "B", 0.5 * H, 3.0 * D, 2500.0, [1000.0, 1e7]),
    ]
}
fn inventories() -> [DissolvedGasInventory; 2] {
    [
        DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap(),
        DissolvedGasInventory::new(VL, 0.006, 0.0).unwrap(),
    ]
}

#[test]
fn unequal_opposing_films_refuse_an_interior_mixture_concentration_peak_atomically() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let laws = laws(&catalog);
    let mut gas = gas(&catalog);
    let before_gas = gas.clone();
    let mut liquid = inventories();
    let before_liquid = liquid;
    let films = [
        DissolutionFilm {
            area_m2: 0.01,
            thickness_m: 1e-3,
            gas_side_transfer: None,
        },
        DissolutionFilm {
            area_m2: 0.01,
            thickness_m: 1e-5,
            gas_side_transfer: None,
        },
    ];
    let mut limits = controls();
    limits.maximum_concentration_mol_m3 = 3.1;
    let bulk = [0.0, 3.0];
    let interface = [H * 60_000.0, 0.5 * H * 40_000.0];
    assert!(bulk.iter().sum::<f64>() < limits.maximum_concentration_mol_m3);
    assert!(interface.iter().sum::<f64>() < limits.maximum_concentration_mol_m3);
    // At B's film edge, B has reached bulk concentration while A is still
    // near its high interfacial value. This independently evaluates the
    // declared linear film profiles, not the implementation's envelope.
    let x = films[1].thickness_m;
    let interior_a = interface[0] + (bulk[0] - interface[0]) * x / films[0].thickness_m;
    near(interior_a + bulk[1], 4.188, 1e-12);
    assert!(interior_a + bulk[1] > limits.maximum_concentration_mol_m3);
    assert!(matches!(
        advance_mixture(&laws, &mut gas, &mut liquid, &films, 0.01, limits, || false),
        Err(DissolvedGasError::Concentration { .. })
    ));
    assert_eq!(gas, before_gas);
    assert_eq!(liquid, before_liquid);
}

#[test]
fn opposite_species_exchange_matches_finite_pair_oracles_and_one_heat_account() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let laws = laws(&catalog);
    let mut gas = gas(&catalog);
    let before_gas = gas.clone();
    let mut liquid = inventories();
    let before_liquid = liquid;
    let duration = 200.0;
    let report = advance_mixture(
        &laws,
        &mut gas,
        &mut liquid,
        &[film(); 2],
        duration,
        controls(),
        || false,
    )
    .unwrap();
    assert!(report.exchanges[0].dissolved_moles_change > 0.0);
    assert!(report.exchanges[1].dissolved_moles_change < 0.0);
    let mut expected_gas_heat = 0.0;
    let mut expected_liquid_heat = 0.0;
    for (index, law) in laws.iter().enumerate() {
        let mw = law.molar_mass_kg_mol;
        let total = before_gas.species_mass_kg(law.species()).unwrap() / mw
            + before_liquid[index].dissolved_moles();
        let fraction = law.henry_mol_m3_pa() * R * T * VL / VG;
        let equilibrium = total * fraction / (1.0 + fraction);
        let rate = law.diffusivity_m2_s() * film().area_m2 / film().thickness_m
            * (1.0 / VL + law.henry_mol_m3_pa() * R * T / VG);
        let expected = equilibrium
            + (before_liquid[index].dissolved_moles() - equilibrium) * (-rate * duration).exp();
        near(liquid[index].dissolved_moles(), expected, 2e-15);
        near(
            gas.species_mass_kg(law.species()).unwrap() / mw + liquid[index].dissolved_moles(),
            total,
            1e-14,
        );
        let dn = expected - before_liquid[index].dissolved_moles();
        expected_gas_heat += dn * R * T;
        expected_liquid_heat -= dn * [6000.0, 2500.0][index];
        assert_eq!(liquid[index].time_s(), duration);
        assert_eq!(
            report.exchanges[index].sources[1].material[0]
                .query_point()
                .iter()
                .find(|(axis, _)| axis == "pressure")
                .unwrap()
                .1,
            gas.state().state().pressure
        );
    }
    near(report.thermostat_heat_into_gas_j, expected_gas_heat, 1e-11);
    near(
        report.thermostat_heat_into_liquid_j,
        expected_liquid_heat,
        1e-11,
    );
    near(
        report.thermostat_heat_into_system_j,
        expected_gas_heat + expected_liquid_heat,
        1e-11,
    );
    near(
        report.energy_residual_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
    assert_eq!(gas.state().state().temperature, T);
    assert_eq!(gas.volume_m3(), VG);
    assert_ne!(
        gas.state().state().sound_speed,
        before_gas.state().state().sound_speed
    );
}

#[test]
fn input_permutation_is_exact_and_subdivision_reuses_the_same_physical_clock() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let laws = laws(&catalog);
    let mut whole_gas = gas(&catalog);
    let mut whole_liquid = inventories();
    let whole = advance_mixture(
        &laws,
        &mut whole_gas,
        &mut whole_liquid,
        &[film(); 2],
        200.0,
        controls(),
        || false,
    )
    .unwrap();
    let mut reversed_gas = gas(&catalog);
    let mut reversed_liquid = [inventories()[1], inventories()[0]];
    let reversed = advance_mixture(
        &[laws[1].clone(), laws[0].clone()],
        &mut reversed_gas,
        &mut reversed_liquid,
        &[film(); 2],
        200.0,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(whole, reversed);
    assert_eq!(whole_gas, reversed_gas);
    assert_eq!(whole_liquid, [reversed_liquid[1], reversed_liquid[0]]);

    let mut split_gas = gas(&catalog);
    let mut split_liquid = inventories();
    let mut heat = 0.0;
    for end in [20.0, 70.0, 130.0, 200.0] {
        heat += advance_mixture(
            &laws,
            &mut split_gas,
            &mut split_liquid,
            &[film(); 2],
            end,
            controls(),
            || false,
        )
        .unwrap()
        .thermostat_heat_into_system_j;
    }
    for index in 0..2 {
        near(
            whole_liquid[index].dissolved_moles(),
            split_liquid[index].dissolved_moles(),
            3e-15,
        );
        near(
            whole_gas.species_mass_kg(laws[index].species()).unwrap(),
            split_gas.species_mass_kg(laws[index].species()).unwrap(),
            1e-15,
        );
    }
    near(whole.thermostat_heat_into_system_j, heat, 1e-10);
}

#[test]
fn one_species_is_the_original_scalar_owner_without_another_gas_update() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let law = law(&catalog);
    let mut scalar_gas = gas(&catalog);
    let mut scalar_liquid = inventories()[0];
    let scalar = law
        .advance(
            &mut scalar_gas,
            &mut scalar_liquid,
            film(),
            200.0,
            controls(),
            || false,
        )
        .unwrap();
    let mut mixture_gas = gas(&catalog);
    let mut mixture_liquid = [inventories()[0]];
    let mixture = advance_mixture(
        &[law],
        &mut mixture_gas,
        &mut mixture_liquid,
        &[film()],
        200.0,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(mixture_gas, scalar_gas);
    assert_eq!(mixture_liquid[0], scalar_liquid);
    assert_eq!(mixture.gas, scalar.gas);
    assert_eq!(
        mixture.thermostat_heat_into_system_j,
        scalar.thermostat_heat_into_system_j
    );
}

#[test]
fn combined_dilution_limits_and_layout_refusals_are_atomic() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let laws = laws(&catalog);
    let original_gas = gas(&catalog);
    let original_liquid = inventories();
    let mut gas = original_gas.clone();
    let mut liquid = original_liquid;
    let mut dilute = controls();
    // Both individual initial interfaces fit 1.3 mol/m3, their sum does not.
    // Use unloaded liquid so the test isolates combined interface concentration.
    let mut unloaded = [original_liquid[0]; 2];
    let before = unloaded;
    dilute.maximum_concentration_mol_m3 = 1.3;
    assert!(matches!(
        advance_mixture(
            &laws,
            &mut gas,
            &mut unloaded,
            &[film(); 2],
            1.0,
            dilute,
            || false
        ),
        Err(DissolvedGasError::Concentration { .. })
    ));
    assert_eq!(gas, original_gas);
    assert_eq!(unloaded, before);
    assert!(
        advance_mixture(
            &[laws[0].clone(), laws[0].clone()],
            &mut gas,
            &mut liquid,
            &[film(); 2],
            200.0,
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(gas, original_gas);
    assert_eq!(liquid, original_liquid);
    for changed in [
        DissolvedGasInventory::new(2.0 * VL, 0.006, 0.0).unwrap(),
        DissolvedGasInventory::new(VL, 0.006, 1.0).unwrap(),
    ] {
        liquid[1] = changed;
        let before = liquid;
        assert!(
            advance_mixture(
                &laws,
                &mut gas,
                &mut liquid,
                &[film(); 2],
                200.0,
                controls(),
                || false
            )
            .is_err()
        );
        assert_eq!(gas, original_gas);
        assert_eq!(liquid, before);
    }
}

#[test]
fn final_cancellation_rolls_back_all_species_and_retry_is_identical() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let laws = laws(&catalog);
    let mut accepted_gas = gas(&catalog);
    let mut accepted_liquid = inventories();
    let mut polls = 0;
    let expected = advance_mixture(
        &laws,
        &mut accepted_gas,
        &mut accepted_liquid,
        &[film(); 2],
        200.0,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut gas = gas(&catalog);
    let before_gas = gas.clone();
    let mut liquid = inventories();
    let before_liquid = liquid;
    let mut seen = 0;
    assert_eq!(
        advance_mixture(
            &laws,
            &mut gas,
            &mut liquid,
            &[film(); 2],
            200.0,
            controls(),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(DissolvedGasError::Cancelled)
    );
    assert_eq!(gas, before_gas);
    assert_eq!(liquid, before_liquid);
    assert_eq!(
        advance_mixture(
            &laws,
            &mut gas,
            &mut liquid,
            &[film(); 2],
            200.0,
            controls(),
            || false
        )
        .unwrap(),
        expected
    );
}

#[test]
fn nonmonotone_total_pressure_cannot_hide_an_unsupported_interior_state() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let make_laws = |pressures| {
        [
            species_law(&catalog, "A", H, D, 6000.0, pressures),
            species_law(&catalog, "B", 0.5 * H, 100.0 * D, 2500.0, pressures),
        ]
    };
    let mut broad_gas = gas(&catalog);
    let mut broad_liquid = inventories();
    let broad = advance_mixture(
        &make_laws([1000.0, 1e7]),
        &mut broad_gas,
        &mut broad_liquid,
        &[film(); 2],
        1e6,
        controls(),
        || false,
    )
    .unwrap();
    // Initial and final total pressure fit the narrow sources, but B exsolves
    // rapidly before A has time to absorb its equilibrium share.
    assert!(broad_gas.state().state().pressure < 110_000.0);
    assert!(broad.pressure_envelope_pa[1] > 110_000.0);
    let mut interior_gas = gas(&catalog);
    let mut interior_liquid = inventories();
    advance_mixture(
        &make_laws([1000.0, 1e7]),
        &mut interior_gas,
        &mut interior_liquid,
        &[film(); 2],
        100.0,
        controls(),
        || false,
    )
    .unwrap();
    assert!(interior_gas.state().state().pressure > 110_000.0);

    let mut narrow_gas = gas(&catalog);
    let before_gas = narrow_gas.clone();
    let mut narrow_liquid = inventories();
    let before_liquid = narrow_liquid;
    assert!(matches!(
        advance_mixture(
            &make_laws([90_000.0, 110_000.0]),
            &mut narrow_gas,
            &mut narrow_liquid,
            &[film(); 2],
            1e6,
            controls(),
            || false
        ),
        Err(DissolvedGasError::Material(_))
    ));
    assert_eq!(narrow_gas, before_gas);
    assert_eq!(narrow_liquid, before_liquid);

    // The same excursion must respect the retained gas transport sources,
    // even when all Henry laws themselves admit the wider pressure envelope.
    let narrow_catalog = crate::thermochemical_gas::tests::model_in_domains(
        vec![0.6, 0.4],
        0.0,
        [200.0, 2000.0],
        [200.0, 2000.0],
        [90_000.0, 110_000.0],
    );
    let wide_henry = [
        species_law(&narrow_catalog, "A", H, D, 6000.0, [1000.0, 1e7]),
        species_law(
            &narrow_catalog,
            "B",
            0.5 * H,
            100.0 * D,
            2500.0,
            [1000.0, 1e7],
        ),
    ];
    let mut gas = gas(&narrow_catalog);
    let before_gas = gas.clone();
    let mut liquid = inventories();
    let before_liquid = liquid;
    assert!(matches!(
        advance_mixture(
            &wide_henry,
            &mut gas,
            &mut liquid,
            &[film(); 2],
            1e6,
            controls(),
            || false
        ),
        Err(DissolvedGasError::Gas(_))
    ));
    assert_eq!(gas, before_gas);
    assert_eq!(liquid, before_liquid);
}

#[test]
fn two_film_mixture_retains_opposing_species_and_one_original_gas_update() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let laws = laws(&catalog);
    let films = [
        DissolutionFilm {
            gas_side_transfer: Some(MolarPermeance::try_new(H * D / film().thickness_m).unwrap()),
            ..film()
        },
        DissolutionFilm {
            gas_side_transfer: Some(
                MolarPermeance::try_new(9.0 * (0.5 * H) * (3.0 * D) / film().thickness_m).unwrap(),
            ),
            ..film()
        },
    ];
    let initial_gas = gas(&catalog);
    let mut gas = initial_gas.clone();
    let initial = inventories();
    let mut liquid = initial;
    let step = advance_mixture(
        &laws,
        &mut gas,
        &mut liquid,
        &films,
        200.0,
        controls(),
        || false,
    )
    .unwrap();
    assert!(step.exchanges[0].dissolved_moles_change > 0.0);
    assert!(step.exchanges[1].dissolved_moles_change < 0.0);
    for i in 0..2 {
        let law = &laws[i];
        let mw = law.molar_mass_kg_mol;
        let total =
            initial_gas.species_mass_kg(law.species()).unwrap() / mw + initial[i].dissolved_moles();
        let h = law.henry_mol_m3_pa();
        let kl = law.diffusivity_m2_s() / films[i].thickness_m;
        let kg = films[i].gas_side_transfer.unwrap().value_si();
        let overall = 1.0 / (1.0 / kl + h / kg);
        let lambda = films[i].area_m2 * overall * (1.0 / VL + h * R * T / VG);
        let equilibrium = total * (h * R * T * VL) / (VG + h * R * T * VL);
        let expected =
            equilibrium + (initial[i].dissolved_moles() - equilibrium) * (-lambda * 200.0).exp();
        near(liquid[i].dissolved_moles(), expected, 3e-15);
        near(
            gas.species_mass_kg(law.species()).unwrap() / mw + liquid[i].dissolved_moles(),
            total,
            1e-14,
        );
        near(step.exchanges[i].relaxation_rate_per_s, lambda, 2e-18);
    }
    near(step.energy_residual_j, 0.0, 1e-8);
    for endpoint in 0..2 {
        let interface_total: f64 = step
            .exchanges
            .iter()
            .map(|e| e.sources[endpoint].interface_partial_pressure_pa)
            .sum();
        let bulk_total = if endpoint == 0 {
            initial_gas.state().state().pressure
        } else {
            gas.state().state().pressure
        };
        assert!(interface_total <= bulk_total);
    }
    let mut reversed_gas = initial_gas;
    let mut reversed_liquid = [initial[1], initial[0]];
    let reversed = advance_mixture(
        &[laws[1].clone(), laws[0].clone()],
        &mut reversed_gas,
        &mut reversed_liquid,
        &[films[1], films[0]],
        200.0,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(step, reversed);
    assert_eq!(gas, reversed_gas);
    assert_eq!(liquid, [reversed_liquid[1], reversed_liquid[0]]);
}

#[test]
fn two_film_mixture_refuses_incompatible_shared_interface_atomically() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let laws = laws(&catalog);
    let films = [
        DissolutionFilm {
            gas_side_transfer: Some(MolarPermeance::try_new(H * D / film().thickness_m).unwrap()),
            ..film()
        },
        DissolutionFilm {
            // p_B,interface is below the total gas pressure by itself, but
            // together with A it would require a negative carrier fraction.
            gas_side_transfer: Some(
                MolarPermeance::try_new(4.0 * (0.5 * H) * (3.0 * D) / film().thickness_m).unwrap(),
            ),
            ..film()
        },
    ];
    let mut gas = gas(&catalog);
    let original_gas = gas.clone();
    let mut liquid = inventories();
    let original_liquid = liquid;
    let a = 0.5 * 60_000.0;
    let b = 0.8 * 40_000.0 + 0.2 * (3.0 / (0.5 * H));
    assert!(b < 100_000.0 && a + b > 100_000.0);
    assert!(matches!(
        advance_mixture(
            &laws,
            &mut gas,
            &mut liquid,
            &films,
            1.0,
            controls(),
            || false
        ),
        Err(DissolvedGasError::Invalid(
            "mixture gas-film interfaces exceed the common total pressure"
        ))
    ));
    assert_eq!(gas, original_gas);
    assert_eq!(liquid, original_liquid);
}
