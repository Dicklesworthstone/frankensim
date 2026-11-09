//! Material handoff: live transported stocks, thermal aqueous state and refusal.

use super::*;
use fs_qty::{Amount, AmountConcentration, VolumetricFlowRate};
use fs_thermochem::{
    MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::{
        caloric::{MolarHeatCapacityCurveV1, PureLiquidHeatCapacityV1},
        solubility::{
            IdealSolubilityV1,
            caloric::{
                SolubilityCaloricComponentV1, SolubilityCaloricModelV1, SolubilityCaloricSolveV1,
            },
            kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
        },
    },
};

const MOLAR_VOLUME: f64 = 2e-5;
const DISSOLUTION_ENTHALPY: f64 = 20_000.0;
const SOLVENT_CP: f64 = 50.0;
const SOLUTE_CP: f64 = 30.0;

fn knots(capacity: f64) -> Vec<(Temperature, MolarThermalQuantityV1)> {
    [280.0, 340.0]
        .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(capacity)))
        .to_vec()
}

fn model(
    conductance: f64,
    policy: PrecipitationPolicyV1,
    molar_volume: f64,
) -> SolubilityThermalKineticsV1 {
    let law = IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("A").unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(600.0),
        DISSOLUTION_ENTHALPY,
        [Temperature::new(280.0), Temperature::new(340.0)],
    )
    .unwrap();
    let kinetic =
        SolubilityKineticsV1::new(law.clone(), VolumetricFlowRate::new(conductance), policy)
            .unwrap();
    let solvent = PureLiquidHeatCapacityV1::new(
        law.solvent().clone(),
        Temperature::new(300.0),
        knots(SOLVENT_CP),
    )
    .unwrap();
    let solute = MolarHeatCapacityCurveV1::new(
        law.solute().clone(),
        Temperature::new(300.0),
        knots(SOLUTE_CP),
    )
    .unwrap();
    SolubilityThermalKineticsV1::new(
        SolubilityCaloricModelV1::new(
            molar_volume,
            solvent,
            vec![SolubilityCaloricComponentV1::new(law, solute).unwrap()],
        )
        .unwrap(),
        vec![kinetic],
    )
    .unwrap()
}

fn material(
    model: &SolubilityThermalKineticsV1,
    dissolved: f64,
    solid: f64,
) -> KineticSolubilitySolution {
    KineticSolubilitySolution::new(
        Temperature::new(300.0),
        model.caloric().solvent_molar_volume_m3_mol(),
        Volume::new(1e-3),
        vec![(
            model.kinetics()[0].clone(),
            SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))
                .unwrap(),
        )],
    )
    .unwrap()
}

fn controls() -> SolubilityThermalKineticControlsV1 {
    SolubilityThermalKineticControlsV1 {
        root: SolubilityCaloricSolveV1 {
            temperature_bracket: [Temperature::new(280.0), Temperature::new(340.0)],
            temperature_tolerance_k: 1e-10,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        },
        initial_step: Time::new(0.25),
        minimum_step: Time::new(1e-6),
        maximum_step: Time::new(0.5),
        relative_tolerance: 2e-5,
        amount_absolute_tolerance_mol: 1e-10,
        temperature_absolute_tolerance_k: 1e-7,
        maximum_attempts: 8192,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}

#[test]
fn transported_inventories_cool_adiabatically_and_drive_the_actual_aqueous_endpoint() {
    let model = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let mut transported = material(&model, 0.15, 1.5);
    transported
        .replace_dissolved_inventories(7e-4, &[0.03])
        .unwrap();
    let before = transported.clone();
    let step = transported
        .advance_with_heat(&model, Time::new(3.0), 0.0, controls(), || false)
        .unwrap();
    let (endpoint, balance) = step.into_parts();
    let aqueous = endpoint.solution();
    let phase = endpoint.last_steps()[0].state();
    let initial = before.last_steps()[0].state();
    let temperature = aqueous.temperature_k();
    let volume = before.solution().solvent_volume_m3();
    let dissolved = aqueous.solute_moles()[0];

    assert_eq!(transported, before);
    assert_eq!(aqueous.solvent_volume_m3(), volume);
    assert_eq!(balance.final_state().solvent_volume().value(), volume);
    assert_eq!(balance.final_state().temperature().value(), temperature);
    assert_eq!(phase, &balance.final_state().phases()[0]);
    assert_eq!(phase.total_solute(), initial.total_solute());
    assert_eq!(dissolved, phase.dissolved().value());
    assert!(temperature < before.solution().temperature_k());
    assert!(dissolved > initial.dissolved().value());
    assert!(dissolved < endpoint.last_steps()[0].equilibrium().dissolved.value());
    assert!(phase.solid().value() > 0.0);

    // Independent constant-Cp balance on the CURRENT transported volume and
    // total. Neither the creation-time inventory nor the equilibrium target
    // supplies this bridge's sensible or dissolution heat.
    let capacity = volume / MOLAR_VOLUME * SOLVENT_CP + initial.total_solute().value() * SOLUTE_CP;
    let energy = capacity * (temperature - before.solution().temperature_k())
        + DISSOLUTION_ENTHALPY * (dissolved - initial.dissolved().value());
    near(energy, 0.0, 1e-7);
    near(balance.energy_residual_j(), energy, 2e-9);
    near(
        aqueous.osmotic_pressure_pa().unwrap(),
        R * temperature * dissolved / volume,
        1e-7,
    );
    near(
        aqueous.water_activity().unwrap(),
        (-MOLAR_VOLUME * dissolved / volume).exp(),
        2e-14,
    );
}

#[test]
fn equal_equilibrium_laws_cannot_substitute_rates_policy_or_solvent_and_late_cancel_is_atomic() {
    let original_model = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let owner = material(&original_model, 0.05, 1.5);
    let before = owner.clone();
    for replacement in [
        model(2e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME),
        model(1e-4, PrecipitationPolicyV1::SeededSurface, MOLAR_VOLUME),
        model(
            1e-4,
            PrecipitationPolicyV1::ExistingSolidOnly,
            1.1 * MOLAR_VOLUME,
        ),
    ] {
        assert_eq!(replacement.kinetics()[0].law(), owner.models()[0].law());
        assert!(matches!(
            owner.advance_with_heat(&replacement, Time::new(1.0), 0.0, controls(), || false),
            Err(KineticSolubilityHeatError::SourceMismatch)
        ));
        assert_eq!(owner, before);
    }

    let mut calls = 0;
    let accepted = owner
        .advance_with_heat(&original_model, Time::new(1.0), 0.0, controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 1);
    let mut observed = 0;
    assert!(matches!(
        owner.advance_with_heat(&original_model, Time::new(1.0), 0.0, controls(), || {
            observed += 1;
            observed == calls
        }),
        Err(KineticSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::Cancelled
        ))
    ));
    assert_eq!(observed, calls);
    assert_eq!(owner, before);
    assert_eq!(
        owner
            .advance_with_heat(&original_model, Time::new(1.0), 0.0, controls(), || false)
            .unwrap(),
        accepted,
    );
}

#[test]
fn closed_surface_heats_sensibly_and_preserves_explicit_trace_crystals() {
    let model = model(0.0, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let owner = material(&model, 1.0, 1e-100);
    let before = owner.clone();
    let capacity = owner.solution().solvent_volume_m3() / MOLAR_VOLUME * SOLVENT_CP
        + owner.last_steps()[0].state().total_solute().value() * SOLUTE_CP;
    let heat = 10.0 * capacity;
    let step = owner
        .advance_with_heat(&model, Time::new(2.0), heat, controls(), || false)
        .unwrap();
    let endpoint = step.solution();

    assert_eq!(owner, before);
    near(endpoint.solution().temperature_k(), 310.0, 1e-9);
    assert_eq!(
        endpoint.last_steps()[0].state(),
        owner.last_steps()[0].state()
    );
    assert_eq!(endpoint.last_steps()[0].state().solid().value(), 1e-100);
    assert_eq!(
        endpoint.solution().solute_moles(),
        owner.solution().solute_moles()
    );
    assert!(endpoint.last_steps()[0].equilibrium().dissolved.value() < 1.0);
    near(step.balance().enthalpy_change_j(), heat, 1e-7);
    near(step.balance().energy_residual_j(), 0.0, 1e-7);
    near(
        endpoint.solution().water_activity().unwrap(),
        owner.solution().water_activity().unwrap(),
        2e-14,
    );
}

#[test]
fn shrinking_particles_share_the_original_heat_balance_and_geometry_source_identity() {
    let fixed = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let particles = SolubilityThermalKineticsV1::new(
        fixed.caloric().clone(),
        vec![
            fixed.kinetics()[0]
                .clone()
                .with_particle_surface(Amount::new(0.2))
                .unwrap(),
        ],
    )
    .unwrap();
    let owner = material(&particles, 0.05, 0.2);
    let before = owner.clone();
    let result = owner
        .advance_with_heat(&particles, Time::new(2.0), 0.0, controls(), || false)
        .unwrap();
    let endpoint = result.solution();
    let phase = endpoint.last_steps()[0].state();
    let temperature = endpoint.solution().temperature_k();
    assert_eq!(owner, before);
    assert_eq!(endpoint.models(), owner.models());
    assert_eq!(
        phase.total_solute(),
        before.last_steps()[0].state().total_solute()
    );
    assert!(phase.solid().value() > 0.0 && phase.solid().value() < 0.2);
    assert!(
        endpoint.models()[0]
            .conductance_at(phase.solid())
            .unwrap()
            .value()
            < 1e-4
    );
    assert!(temperature < 300.0);

    let capacity = 1e-3 / MOLAR_VOLUME * SOLVENT_CP + 0.25 * SOLUTE_CP;
    let energy = capacity * (temperature - 300.0)
        + DISSOLUTION_ENTHALPY * (phase.dissolved().value() - 0.05);
    near(energy, 0.0, 1e-7);
    near(result.balance().energy_residual_j(), energy, 2e-9);
    near(
        endpoint.solution().osmotic_pressure_pa().unwrap(),
        R * temperature * phase.dissolved().value() / 1e-3,
        1e-7,
    );

    let fixed_owner = material(&fixed, 0.05, 0.2);
    let fixed_result = fixed_owner
        .advance_with_heat(&fixed, Time::new(2.0), 0.0, controls(), || false)
        .unwrap();
    assert!(phase.dissolved().value() < fixed_result.solution().solution().solute_moles()[0]);
    assert!(temperature > fixed_result.solution().solution().temperature_k());

    let other_geometry = SolubilityThermalKineticsV1::new(
        fixed.caloric().clone(),
        vec![
            fixed.kinetics()[0]
                .clone()
                .with_particle_surface(Amount::new(0.4))
                .unwrap(),
        ],
    )
    .unwrap();
    assert!(matches!(
        owner.advance_with_heat(&other_geometry, Time::new(2.0), 0.0, controls(), || false),
        Err(KineticSolubilityHeatError::SourceMismatch)
    ));
    assert_eq!(owner, before);
}
