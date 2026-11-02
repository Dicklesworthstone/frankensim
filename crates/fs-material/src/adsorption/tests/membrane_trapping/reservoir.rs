//! Original gas sources, finite inventories and source-relative thermal closure.
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::GasSpeciesBoundarySolve;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::{
    reservoir::{IsothermalGasMembraneError,IsothermalGasMembraneSolve},
    slab::{MembraneDiffusionSolve,finite::FiniteMembraneDiffusionWorkspace,
        adaptive::{MembraneDiffusionAdvanceControls,MembraneDiffusionAdvanceError},
        trapping::adaptive::TrappingMembraneAdvanceControls},
    trapping::reservoir::{IsothermalTrappingGasSolve,IsothermalTrappingGasError},
};

fn gas_model(reference: f64, maximum_pressure: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species_with_energy_reference(
                "A", MASS, 3.5, 0.0, reference,
            )],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![crate::thermochemical_gas::tests::transport_in_domain(
            "A",
            1.7e-5,
            0.025,
            [200.0, 2000.0],
            [1000.0, maximum_pressure],
        )],
        None,
    )
    .unwrap()
}
fn gases(model: &ThermochemicalGasModel, pressures: [f64; 2]) -> [ClosedGasParcel; 2] {
    pressures.map(|p| ClosedGasParcel::new(model.clone(), point(T0, p), 2e-11 * R * T0).unwrap())
}
fn controls() -> IsothermalTrappingGasSolve {
    IsothermalTrappingGasSolve {
        membrane: TrappingMembraneSolve {
            balance: MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: 1e-17,
                relative_amount_tolerance: 1e-12,
            },
            ..TrappingMembraneSolve::default()
        },
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
fn temporal() -> TrappingMembraneAdvanceControls {
    TrappingMembraneAdvanceControls {
        temporal: MembraneDiffusionAdvanceControls {
            step: controls().membrane.balance,
            initial_step_s: 2.0,
            minimum_step_s: 1e-8,
            maximum_step_s: 5.0,
            maximum_evaluations: 20_000,
            absolute_pressure_tolerance_pa: 30.0,
            absolute_amount_tolerance_mol: 1e-10,
            relative_tolerance: 0.0,
            amount_balance_tolerance_mol: 1e-13,
        },
        maximum_iterations: 64,
        maximum_backtracks: 32,
    }
}
fn film<'a>(
    catalog: &ThermochemicalGasModel,
    card: &'a MaterialCard,
    context: &'a QueryPoint,
    traps: bool,
    preloaded: bool,
) -> ResolvedTrappingMembraneDiffusion<'a> {
    let source = traps.then(|| LangmuirMembraneTrapSource {
        adsorbent: law(catalog, vacuum_source(PRESSURE)),
        dry_skeleton_density_kg_m3: 4.0,
    });
    let amount = if preloaded {
        2.0 * ((1e-4 / 8.0) * 0.01)
    } else {
        0.0
    };
    ResolvedTrappingMembraneDiffusion::new(
        mobile(card, context),
        &[source],
        &[amount; 8],
        [PRESSURE, 0.0],
        || false,
    )
    .unwrap()
}

#[test]
fn no_trap_parcels_recover_original_finite_gas_membrane() {
    let model = gas_model(0.0, 1e7);
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let controls = controls();
    let parcels = gases(&model, [PRESSURE, 20_000.0]);
    let original = mobile(&card, &context)
        .with_finite_gas_reservoirs(
            parcels.clone(),
            IsothermalGasMembraneSolve {
                membrane: controls.membrane.balance,
                gas: controls.gas,
                mass_tolerance_kg: controls.mass_tolerance_kg,
            },
            || false,
        )
        .unwrap();
    let trapped = film(&model, &card, &context, false, false)
        .with_finite_gas_reservoirs(parcels, controls, || false)
        .unwrap();
    let legacy = original
        .step_candidate(
            10.0,
            &mut FiniteMembraneDiffusionWorkspace::new(
                original.membrane().partial_pressures_pa().len(),
            )
            .unwrap(),
            || false,
        )
        .unwrap();
    let next = trapped
        .step_candidate(10.0, &mut trapped.workspace().unwrap(), || false)
        .unwrap();
    for (a, b) in next
        .state()
        .membrane()
        .mobile()
        .partial_pressures_pa()
        .iter()
        .zip(legacy.state().membrane().partial_pressures_pa())
    {
        assert!((a - b).abs() < 1e-7);
    }
    for i in 0..2 {
        assert!(
            (next.state().gases()[i].mass_kg() - legacy.state().gases()[i].mass_kg()).abs() < 1e-19
        );
        assert!(
            (next.state().gases()[i].state().state().pressure
                - legacy.state().gases()[i].state().state().pressure)
                .abs()
                < 1e-7
        );
    }
    assert_eq!(next.report().accounting.trapped_mass_change_kg, 0.0);
    assert!(next.report().accounting.total_mass_residual_kg.abs() < 1e-16);
}

#[test]
fn actual_gas_uptake_and_adaptive_windows_preserve_mass_temperature_and_caloric_datum_invariance() {
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let mut reference_pressures = None;
    for reference in [0.0, 1e12] {
        let model = gas_model(reference, 1e7);
        let state = film(&model, &card, &context, true, false)
            .with_finite_gas_reservoirs(gases(&model, [PRESSURE, 20_000.0]), controls(), || false)
            .unwrap();
        let mut work = state.workspace().unwrap();
        let step = state.step_candidate(10.0, &mut work, || false).unwrap();
        assert!(step.report().diffusion.trapped_change_mol > 0.0);
        assert!(step.report().accounting.mass_residual_kg.abs() < 1e-16);
        for i in 0..2 {
            assert_eq!(step.state().gases()[i].state().state().temperature, T0);
            assert_eq!(step.state().gases()[i].model(), &model);
            assert_eq!(
                step.state().gases()[i].volume_m3(),
                state.gases()[i].volume_m3()
            );
            assert_eq!(
                step.state().gases()[i].state().state().pressure,
                step.state().finite_state().reservoir_pressures_pa()[i]
            );
            assert!(step.report().gas[i].energy_residual_j.abs() < 1e-10);
            assert_eq!(
                step.report().accounting.thermostat_heat_into_gas_j[i],
                -step.report().gas[i].mass_into_gas_kg * (R * T0 / MASS)
            );
        }
        let window = state
            .advance_candidate(
                20.0,
                temporal(),
                &mut state.adaptive_workspace().unwrap(),
                || false,
            )
            .unwrap();
        assert!(window.report().diffusion.transport.accepted_substeps > 2);
        assert_eq!(window.state().membrane().mobile().elapsed_s(), 20.0);
        assert!(window.report().accounting.total_mass_residual_kg.abs() < 1e-16);
        assert!(window.report().accounting.mass_residual_kg.abs() < 1e-16);
        let pressures = window.state().finite_state().reservoir_pressures_pa();
        if let Some(previous) = reference_pressures {
            assert_eq!(pressures, previous);
        }
        reference_pressures = Some(pressures);
    }
}

#[test]
fn gas_source_failure_caloric_substitution_and_late_cancellation_are_atomic() {
    let model = gas_model(0.0, 3_500.0);
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let state = film(&model, &card, &context, true, true)
        .with_finite_gas_reservoirs(gases(&model, [2_000.0; 2]), controls(), || false)
        .unwrap();
    let saved_gases = state.gases().clone();
    let saved_film = state.finite_state().state().clone();
    let saved_support = state.finite_state().trap_source_support().to_vec();
    let mut work = state.workspace().unwrap();
    let native = state
        .finite_state()
        .step_candidate(100.0, controls().membrane, &mut work, || false)
        .unwrap();
    assert!(native.state().reservoir_pressures_pa()[0] > 3_500.0);
    assert!(matches!(
        state.step_candidate(100.0, &mut work, || false),
        Err(IsothermalTrappingGasError::Gas(
            IsothermalGasMembraneError::Gas { .. }
        ))
    ));
    let mut limits = temporal();
    limits.temporal.initial_step_s = 100.0;
    limits.temporal.maximum_step_s = 100.0;
    assert!(matches!(
        state.advance_candidate(
            100.0,
            limits,
            &mut state.adaptive_workspace().unwrap(),
            || false
        ),
        Err(MembraneDiffusionAdvanceError::Native(
            IsothermalTrappingGasError::Gas(IsothermalGasMembraneError::Gas { .. })
        ))
    ));
    assert_eq!(state.gases(), &saved_gases);
    assert_eq!(state.finite_state().state(), &saved_film);
    assert_eq!(state.finite_state().trap_source_support(), saved_support);

    let normal = gas_model(0.0, 1e7);
    let foreign = gas_model(1e12, 1e7);
    assert!(matches!(
        film(&normal, &card, &context, true, false).with_finite_gas_reservoirs(
            gases(&foreign, [PRESSURE, 20_000.0]),
            controls(),
            || false
        ),
        Err(IsothermalTrappingGasError::Source(
            ResolvedTrappingMembraneError::Adsorption(AdsorptionError::Invalid(_))
        ))
    ));
    let state = film(&normal, &card, &context, true, false)
        .with_finite_gas_reservoirs(gases(&normal, [PRESSURE, 20_000.0]), controls(), || false)
        .unwrap();
    let saved_gases = state.gases().clone();
    let saved_film = state.finite_state().state().clone();
    let mut polls = 0;
    let mut work = state.adaptive_workspace().unwrap();
    state
        .advance_candidate(2.0, temporal(), &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        state.advance_candidate(2.0, temporal(), &mut work, || {
            calls += 1;
            calls == polls
        }),
        Err(MembraneDiffusionAdvanceError::Cancelled)
    ));
    assert_eq!(state.gases(), &saved_gases);
    assert_eq!(state.finite_state().state(), &saved_film);
}
