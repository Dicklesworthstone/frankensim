//! G0/G1/G3 actual finite-mixture uptake and shared-site displacement.
#[path = "mixture/adaptive.rs"]
mod adaptive;
use super::*;
use crate::adsorption::competitive::CompetitiveAdsorbent;
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::GasSpeciesBoundarySolve;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::{
    mixture::{IsothermalGasMembraneMixture, MembraneMixtureChannelSource, competitive::*},
    reservoir::IsothermalGasMembraneSolve,
    slab::{MembraneDiffusionSolve, competitive::CompetitiveMembraneSolve},
};

fn id(s: &str) -> SpeciesId {
    SpeciesId::new(s).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}
fn solve() -> CompetitiveMembraneSolve {
    CompetitiveMembraneSolve {
        balance: MembraneDiffusionSolve {
            absolute_amount_tolerance_mol: 1e-18,
            relative_amount_tolerance: 1e-12,
        },
        ..CompetitiveMembraneSolve::default()
    }
}
fn gas_solve() -> IsothermalGasMembraneSolve {
    IsothermalGasMembraneSolve {
        membrane: solve().balance,
        gas: [GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 250.0,
                maximum_temperature_k: 500.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-16,
            energy_tolerance_j: 1e-10,
        }; 2],
        mass_tolerance_kg: 1e-16,
    }
}
fn site_law(
    catalog: &ThermochemicalGasModel,
    name: &str,
    affinity: f64,
    bounds: [f64; 2],
) -> LangmuirAdsorbent {
    let mut source = vacuum_source(bounds[1]);
    let mut claims = Vec::new();
    for (key, dims, value) in [
        (LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_CAPACITY_DIMS, 0.5),
        (
            LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
            LANGMUIR_AFFINITY_DIMS,
            affinity,
        ),
        (LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS, 0.2),
    ] {
        let mut claim = source_claim(key, dims, value, bounds[1]);
        claim.validity = claim.validity.with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            bounds[0],
            bounds[1],
        );
        claims.push(claim);
    }
    source.material = card(&format!("explicit mixed-gas lattice {name}"), claims);
    let mut caloric = calorics(catalog, 0.0, HEAT);
    if name != "A" {
        caloric.bound_adsorbate = chart("bound B", 2.5 * R / 0.004, 0.0);
        caloric.bound_enthalpy_offset_j_kg = catalog
            .species_specific_enthalpy_j_kg(&id(name), T0)
            .unwrap()
            - caloric
                .bound_adsorbate
                .state_at_temperature(T0)
                .unwrap()
                .specific_enthalpy_j_kg()
            - HEAT / 0.004;
    }
    LangmuirAdsorbent::try_new(catalog, id(name), source, caloric).unwrap()
}
fn sites(catalog: &ThermochemicalGasModel, bounds: [f64; 2]) -> CompetitiveMembraneTrapSource {
    CompetitiveMembraneTrapSource {
        adsorbent: Arc::new(
            CompetitiveAdsorbent::try_new(vec![
                site_law(catalog, "A", 2e-5, bounds),
                site_law(catalog, "B", 2e-3, bounds),
            ])
            .unwrap(),
        ),
        dry_skeleton_density_kg_m3: 4.0,
    }
}
fn mixture<'a>(
    card: &'a MaterialCard,
    context: &'a QueryPoint,
    gases: [ClosedGasParcel; 2],
    profile: [f64; 2],
    reverse: bool,
) -> IsothermalGasMembraneMixture<'a> {
    let mut sources = Vec::new();
    for (s, name) in ["A", "B"].iter().enumerate() {
        sources.push(MembraneMixtureChannelSource {
            layers: vec![MembraneDiffusionLayerSource {
                thickness_m: 1e-4,
                cells: 8,
                material: MembranePermeabilityMaterial {
                    species: id(name),
                    material: card,
                    context,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                },
            }],
            interfaces: vec![],
            initial_partial_pressures_pa: vec![profile[s]; 8],
        });
    }
    if reverse {
        sources.reverse();
    }
    IsothermalGasMembraneMixture::new(T0, 0.01, &sources, gases, gas_solve(), || false).unwrap()
}
fn gases(catalog: &ThermochemicalGasModel, pressure: f64) -> [ClosedGasParcel; 2] {
    [1e-6, 2e-6].map(|v| ClosedGasParcel::new(catalog.clone(), point(T0, pressure), v).unwrap())
}
fn audit(state: &IsothermalCompetitiveGasMembrane<'_>, step: &IsothermalCompetitiveGasStep) {
    assert!(step.mass_residual_kg.abs() <= gas_solve().mass_tolerance_kg);
    assert!(step.total_mass_residual_kg.abs() <= gas_solve().mass_tolerance_kg);
    for (s, row) in step.species.iter().enumerate() {
        assert!(row.mass_residual_kg.abs() <= gas_solve().mass_tolerance_kg);
        assert!(row.total_mass_residual_kg.abs() <= gas_solve().mass_tolerance_kg);
        for side in 0..2 {
            assert_eq!(
                row.actual_partial_pressures_pa[side],
                state.gases()[side]
                    .species_partial_pressure_pa(&row.species)
                    .unwrap()
            );
            assert_eq!(
                state.finite_state().reservoir_partial_pressures_pa()[s][side],
                row.actual_partial_pressures_pa[side]
            );
        }
    }
    for gas in state.gases() {
        assert_eq!(gas.state().state().temperature, T0);
        assert_eq!(&gas.model().evaluate(gas.point()).unwrap(), gas.state());
        near(
            gas.state().state().density * gas.volume_m3(),
            gas.mass_kg(),
            1e-18,
        );
    }
    assert_eq!(step.gas_energy_residual_j, 0.0);
}

#[test]
fn incoming_gas_displaces_loaded_species_and_actual_gases_conserve_each_inventory() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let original_gases = gases(&catalog, PRESSURE);
    let source = sites(&catalog, [60_000.0, 200_000.0]);
    let loadings = [(id("A"), vec![2e-7; 8]), (id("B"), vec![0.0; 8])];
    let state = mixture(
        &card,
        &context,
        original_gases.clone(),
        [50_000.0; 2],
        false,
    )
    .with_competitive_sites(&[Some(source.clone())], &loadings, solve(), || false)
    .unwrap();
    // Both partial pressures are below the source's 60 kPa floor. Source
    // admission correctly uses summed/actual totals, while kinetics use q_i.
    for name in ["A", "B"] {
        assert_eq!(
            state.trap_source_support(&id(name)).unwrap()[0]
                .as_ref()
                .unwrap()[0]
                .source
                .state_coordinate("pressure"),
            Some(PRESSURE)
        );
    }
    let saved = state.finite_state().clone();
    let mut work = state.workspace(1_000_000).unwrap();
    let next = state.step_candidate(1.0, &mut work, || false).unwrap();
    audit(next.state(), next.step());
    assert!(next.step().species[0].trapped_mass_change_kg < 0.0);
    assert!(next.step().species[1].trapped_mass_change_kg > 0.0);
    assert!(
        next.step().species[0]
            .actual_gas_mass_changes_kg
            .iter()
            .sum::<f64>()
            > 0.0
    );
    assert!(
        next.step().species[1]
            .actual_gas_mass_changes_kg
            .iter()
            .sum::<f64>()
            < 0.0
    );
    for cell in 0..8 {
        let sum: f64 = (0..2)
            .map(|s| next.state().membrane().trapped_moles(s).unwrap()[cell])
            .sum();
        assert!(sum <= next.state().membrane().site_capacities_mol()[cell]);
    }
    assert_eq!(state.finite_state(), &saved);
    assert_eq!(state.gases(), &original_gases);
    let reversed = mixture(&card, &context, original_gases, [50_000.0; 2], true)
        .with_competitive_sites(
            &[Some(source)],
            &[loadings[1].clone(), loadings[0].clone()],
            solve(),
            || false,
        )
        .unwrap();
    let repeated = reversed
        .step_candidate(1.0, &mut reversed.workspace(1_000_000).unwrap(), || false)
        .unwrap();
    assert_eq!(next.state().finite_state(), repeated.state().finite_state());
    assert_eq!(next.step(), repeated.step());
}

#[test]
fn absent_gas_species_reenters_from_traps_and_no_site_limit_matches_original_mixture() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let pure_b = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap())
        .unwrap();
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let mut endpoints = gases(&catalog, PRESSURE);
    endpoints[1] = ClosedGasParcel::new(pure_b, point(T0, PRESSURE), 2e-6).unwrap();
    let loaded = mixture(&card, &context, endpoints, [0.0, PRESSURE], false)
        .with_competitive_sites(
            &[Some(sites(&catalog, [0.0, 1e7]))],
            &[(id("A"), vec![2e-7; 8]), (id("B"), vec![0.0; 8])],
            solve(),
            || false,
        )
        .unwrap();
    assert_eq!(loaded.gases()[1].species_mass_kg(&id("A")).unwrap(), 0.0);
    let next = loaded
        .step_candidate(1.0, &mut loaded.workspace(1_000_000).unwrap(), || false)
        .unwrap();
    audit(next.state(), next.step());
    assert!(next.state().gases()[0].species_mass_kg(&id("A")).unwrap() > 0.0);
    assert!(next.state().gases()[1].species_mass_kg(&id("A")).unwrap() > 0.0);

    let original = mixture(
        &card,
        &context,
        gases(&catalog, PRESSURE),
        [20_000.0, 80_000.0],
        false,
    );
    let free = original
        .clone()
        .with_competitive_sites(
            &[None],
            &[(id("A"), vec![0.0; 8]), (id("B"), vec![0.0; 8])],
            solve(),
            || false,
        )
        .unwrap();
    let plain = original
        .step_candidate(0.5, &mut original.workspace().unwrap(), || false)
        .unwrap();
    let empty = free
        .step_candidate(0.5, &mut free.workspace(1_000_000).unwrap(), || false)
        .unwrap();
    audit(empty.state(), empty.step());
    for (s, name) in ["A", "B"].iter().enumerate() {
        let a = empty
            .state()
            .membrane()
            .mobile(s)
            .unwrap()
            .partial_pressures_pa();
        let b = plain
            .state()
            .finite_state(&id(name))
            .unwrap()
            .membrane()
            .partial_pressures_pa();
        for (&a, &b) in a.iter().zip(b) {
            near(a, b, 1e-7);
        }
        for side in 0..2 {
            near(
                empty.state().gases()[side]
                    .species_mass_kg(&id(name))
                    .unwrap(),
                plain.state().gases()[side]
                    .species_mass_kg(&id(name))
                    .unwrap(),
                1e-18,
            );
        }
    }
}

#[test]
fn source_domain_failure_foreign_calorics_and_late_cancellation_publish_nothing() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    // Loaded sites release enough adsorbate to exceed the original
    // total-pressure support. The raw numerical endpoint remains valid.
    let state = mixture(
        &card,
        &context,
        gases(&catalog, PRESSURE),
        [50_000.0; 2],
        false,
    )
    .with_competitive_sites(
        &[Some(sites(&catalog, [0.0, PRESSURE]))],
        &[(id("A"), vec![2.5e-7; 8]), (id("B"), vec![0.0; 8])],
        solve(),
        || false,
    )
    .unwrap();
    let saved = state.finite_state().clone();
    let saved_gases = state.gases().clone();
    let receipts = state.trap_source_support(&id("A")).unwrap().to_vec();
    let mut work = state.workspace(1_000_000).unwrap();
    let numerical = state
        .finite_state()
        .step_candidate(10.0, solve(), &mut work, || false)
        .unwrap();
    assert!((0..2).any(|side| {
        numerical
            .state()
            .reservoir_partial_pressures_pa()
            .iter()
            .map(|p| p[side])
            .sum::<f64>()
            > PRESSURE
    }));
    assert!(matches!(
        state.step_candidate(10.0, &mut work, || false),
        Err(IsothermalCompetitiveGasError::Source(
            ResolvedTrappingMembraneError::Adsorption(AdsorptionError::Material(_))
        ))
    ));
    assert_eq!(state.finite_state(), &saved);
    assert_eq!(state.gases(), &saved_gases);
    assert_eq!(state.trap_source_support(&id("A")).unwrap(), receipts);

    let valid = mixture(
        &card,
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
    .unwrap();
    let mut polls = 0;
    valid
        .step_candidate(1.0, &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut interrupted = 0;
    assert!(matches!(
        valid.step_candidate(1.0, &mut work, || {
            interrupted += 1;
            interrupted >= polls
        }),
        Err(IsothermalCompetitiveGasError::Cancelled)
    ));
    assert_eq!(valid.elapsed_s(), 0.0);
    assert_eq!(valid.gases(), &saved_gases);
    assert!(valid.workspace(1).is_err());

    let foreign = model(vec![0.5, 0.5], 0.01);
    assert!(matches!(
        mixture(
            &card,
            &context,
            gases(&catalog, PRESSURE),
            [50_000.0; 2],
            false
        )
        .with_competitive_sites(
            &[Some(sites(&foreign, [0.0, 1e7]))],
            &[(id("A"), vec![0.0; 8]), (id("B"), vec![0.0; 8])],
            solve(),
            || false
        ),
        Err(IsothermalCompetitiveGasError::Source(
            ResolvedTrappingMembraneError::Adsorption(AdsorptionError::Invalid(_))
        ))
    ));
}
