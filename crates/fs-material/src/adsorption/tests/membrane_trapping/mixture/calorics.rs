//! G0/G3 caloric closure, source limits and reference-invariant thermostat heat.
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::mixture::{
    adaptive::MembraneMixtureAdvanceControls,
    competitive::calorics::*,
};

const BYTES: usize = 4 * 1024 * 1024;

fn caloric_card(q: f64, pressure_upper: f64) -> MaterialCard {
    card(
        "source-qualified condensed mobile calorics",
        vec![
            source_claim(
                MEMBRANE_DIFFUSIVITY_PROPERTY,
                MEMBRANE_DIFFUSIVITY_DIMS,
                1e-10,
                1e7,
            ),
            source_claim(
                MOLAR_HENRY_SOLUBILITY_PROPERTY,
                MOLAR_HENRY_SOLUBILITY_DIMS,
                2e-5,
                1e7,
            ),
            source_claim(
                MOBILE_RELATIVE_INTERNAL_ENERGY_PROPERTY,
                MOLAR_INTERNAL_ENERGY_DIMS,
                q,
                pressure_upper,
            ),
        ],
    )
}
fn controls() -> MembraneMixtureAdvanceControls {
    MembraneMixtureAdvanceControls {
        initial_step_s: 0.5,
        minimum_step_s: 1e-6,
        maximum_step_s: 0.5,
        maximum_evaluations: 20000,
        maximum_total_cells: 100,
        absolute_pressure_tolerance_pa: 200.0,
        absolute_amount_tolerance_mol: 2e-9,
        relative_tolerance: 0.0,
        species_balance_tolerance_mol: 1e-13,
        mass_balance_tolerance_kg: 1e-14,
        gas_energy_tolerance_j: 1e-9,
    }
}
fn reference_catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.0, reference),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, reference),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap()
}

#[test]
fn fixed_and_adaptive_mobile_heat_survive_huge_formation_reference_changes() {
    let material = caloric_card(-4000.0, 1e7);
    let context = point(T0, PRESSURE);
    let mut fixed = None;
    let mut adaptive = None;
    for reference in [0.0, 1e24] {
        let catalog = reference_catalog(reference);
        let raw = mixture(
            &material,
            &context,
            gases(&catalog, PRESSURE),
            [0.0; 2],
            false,
        )
        .with_competitive_sites(
            &[None],
            &[(id("A"), vec![0.0; 8]), (id("B"), vec![0.0; 8])],
            solve(),
            || false,
        )
        .unwrap();
        let original = raw
            .step_candidate(1.0, &mut raw.workspace(BYTES).unwrap(), || false)
            .unwrap();
        let state = raw.with_caloric_closure(1e-10, || false).unwrap();
        let result = state
            .step_candidate(1.0, &mut state.workspace(BYTES).unwrap(), || false)
            .unwrap();
        assert_eq!(
            result.state().state().finite_state(),
            original.state().finite_state()
        );
        assert_eq!(result.state().state().gases(), original.state().gases());
        let account = result.calorics();
        let amount: f64 = result
            .transport()
            .diffusion
            .species
            .iter()
            .map(|(_, row)| row.mobile_change_mol)
            .sum();
        near(
            account.thermostat_heat_into_membrane_j,
            -4000.0 * amount,
            1e-13,
        );
        assert!(account.thermostat_heat_into_membrane_j < 0.0);
        assert_eq!(account.trapping_contribution_j, 0.0);
        near(
            account.total_relative_storage_change_j,
            account.total_thermostat_heat_into_system_j,
            1e-10,
        );
        assert!(account.energy_residual_j.abs() <= 1e-10);
        if let Some(expected) = &fixed {
            assert_eq!(account, expected);
        } else {
            fixed = Some(account.clone());
        }
        for name in ["A", "B"] {
            let source = result.state().state().source_support(&id(name)).unwrap();
            for endpoint in &source.layers[0] {
                assert_eq!(
                    endpoint
                        .property(MOBILE_RELATIVE_INTERNAL_ENERGY_PROPERTY)
                        .unwrap()
                        .value_si(),
                    -4000.0
                );
                assert_eq!(endpoint.card_identity(), material.content_hash());
            }
        }
        let window = state
            .advance_candidate(
                1.0,
                controls(),
                &mut state.adaptive_workspace(BYTES).unwrap(),
                || false,
            )
            .unwrap();
        assert_eq!(window.calorics().interval_s, [0.0, 1.0]);
        near(
            window.calorics().membrane_relative_storage_change_j,
            -4000.0
                * window
                    .transport()
                    .species
                    .iter()
                    .map(|r| r.actual_sorbed_amount_change_mol)
                    .sum::<f64>(),
            1e-13,
        );
        assert_eq!(
            window.calorics().thermostat_heat_into_gas_j,
            window.transport().thermostat_heat_into_gas_j
        );
        if let Some(expected) = &adaptive {
            assert_eq!(window.calorics(), expected);
        } else {
            adaptive = Some(window.calorics().clone());
        }
        if reference > 0.0 {
            assert!(window.transport().enthalpy_into_gas_j[0].abs() > 1e10);
        }
    }
}

#[test]
fn layered_displacement_closes_independent_internal_energy_without_double_counting_binding_heat() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let cards = [caloric_card(-2000.0, 1e7), caloric_card(-5000.0, 1e7)];
    let context = point(T0, PRESSURE);
    let sources: Vec<_> = ["A", "B"]
        .into_iter()
        .map(|name| MembraneMixtureChannelSource {
            layers: cards
                .iter()
                .map(|material| MembraneDiffusionLayerSource {
                    thickness_m: 5e-5,
                    cells: 4,
                    material: MembranePermeabilityMaterial {
                        species: id(name),
                        material,
                        context: &context,
                        selection: MaterialPropertySelection::SingleClaimOnly,
                    },
                })
                .collect(),
            interfaces: vec![None],
            initial_partial_pressures_pa: vec![50_000.0; 8],
        })
        .collect();
    let site = sites(&catalog, [0.0, 1e7]);
    let raw = IsothermalGasMembraneMixture::new(
        T0,
        0.01,
        &sources,
        gases(&catalog, PRESSURE),
        gas_solve(),
        || false,
    )
    .unwrap()
    .with_competitive_sites(
        &[Some(site.clone()), Some(site)],
        &[(id("A"), vec![2e-7; 8]), (id("B"), vec![0.0; 8])],
        solve(),
        || false,
    )
    .unwrap();
    let state = raw.with_caloric_closure(1e-10, || false).unwrap();
    let result = state
        .step_candidate(1.0, &mut state.workspace(BYTES).unwrap(), || false)
        .unwrap();
    let mut independent_energy = 0.0;
    for (s, row) in result.calorics().species.iter().enumerate() {
        let molar_mass = [0.028, 0.004][s];
        let h = catalog
            .species_specific_enthalpy_j_kg(&row.species, T0)
            .unwrap();
        for side in 0..2 {
            independent_energy += (result.state().state().gases()[side]
                .species_mass_kg(&row.species)
                .unwrap()
                - state.state().gases()[side]
                    .species_mass_kg(&row.species)
                    .unwrap())
                * (h - R * T0 / molar_mass);
        }
        for (layer, &qm) in row.layers.iter().zip(&[-2000.0, -5000.0]) {
            let qb = layer.trapped_relative_internal_energy_j_mol.unwrap();
            near(qb, -HEAT, 1e-10);
            independent_energy += layer.mobile_amount_change_mol * (h * molar_mass + qm)
                + layer.trapped_amount_change_mol * (h * molar_mass + qb);
            near(
                layer.dissolution_contribution_j + layer.trapping_contribution_j,
                layer.mobile_relative_storage_change_j + layer.trapped_relative_storage_change_j,
                1e-15,
            );
            near(
                layer.trapping_contribution_j,
                (qb - qm) * layer.trapped_amount_change_mol,
                1e-15,
            );
            if s == 0 {
                assert!(layer.trapped_amount_change_mol < 0.0);
                assert!(layer.trapping_contribution_j > 0.0);
            } else {
                assert!(layer.trapped_amount_change_mol > 0.0);
                assert!(layer.trapping_contribution_j < 0.0);
            }
        }
    }
    near(
        independent_energy,
        result.calorics().total_thermostat_heat_into_system_j,
        1e-10,
    );
    assert!(result.calorics().energy_residual_j.abs() <= 1e-10);
}

#[test]
fn missing_or_late_invalid_caloric_sources_and_late_cancellation_publish_no_state() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let context = point(T0, PRESSURE);
    let build = |material| {
        mixture(
            material,
            &context,
            gases(&catalog, PRESSURE),
            [50_000.0; 2],
            false,
        )
        .with_competitive_sites(
            &[Some(sites(&catalog, [0.0, 1e7]))],
            &[(id("A"), vec![2.5e-7; 8]), (id("B"), vec![0.0; 8])],
            solve(),
            || false,
        )
        .unwrap()
    };
    let absent = membrane_card();
    assert!(matches!(
        build(&absent).with_caloric_closure(1e-10, || false),
        Err(MembraneCaloricError::Source(_))
    ));
    let limited = caloric_card(-4000.0, PRESSURE);
    let state = build(&limited)
        .with_caloric_closure(1e-10, || false)
        .unwrap();
    let saved = state.state().finite_state().clone();
    let saved_gases = state.state().gases().clone();
    assert!(matches!(
        state.step_candidate(10.0, &mut state.workspace(BYTES).unwrap(), || false),
        Err(MembraneCaloricError::Native(
            IsothermalCompetitiveGasError::Source(_)
        ))
    ));
    assert!(matches!(
        state.advance_candidate(
            10.0,
            controls(),
            &mut state.adaptive_workspace(BYTES).unwrap(),
            || false
        ),
        Err(MembraneCaloricError::Adaptive(
            crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::mixture::competitive::adaptive::IsothermalCompetitiveGasAdvanceError::Native(
                IsothermalCompetitiveGasError::Source(_)
            )
        ))
    ));
    assert_eq!(state.state().finite_state(), &saved);
    assert_eq!(state.state().gases(), &saved_gases);
    let valid = caloric_card(-4000.0, 1e7);
    let state = build(&valid).with_caloric_closure(1e-10, || false).unwrap();
    let mut work = state.workspace(BYTES).unwrap();
    let mut polls = 0;
    let success = state
        .step_candidate(1.0, &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        state.step_candidate(1.0, &mut work, || {
            seen += 1;
            seen >= polls
        }),
        Err(MembraneCaloricError::Cancelled)
    ));
    assert_eq!(state.state().elapsed_s(), 0.0);
    assert_eq!(state.state().finite_state(), &saved);
    assert_eq!(
        success.calorics(),
        state
            .step_candidate(1.0, &mut work, || false)
            .unwrap()
            .calorics()
    );
}
