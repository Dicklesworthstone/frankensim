//! Material handoff: live transported stocks, thermal aqueous state and refusal.

use super::*;
use crate::osmosis::{
    nonideal::VirialSolutionModel, solubility::kinetics::VirialSolubilityThermalModel,
};
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

#[test]
fn virial_phase_binding_cannot_silently_drop_excess_calorics() {
    let thermal = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let original = material(&thermal, 0.15, 1.5)
        .with_nonideal_thermodynamics(Arc::new(
            crate::osmosis::nonideal::VirialSolutionModel::new(
                300.0,
                MOLAR_VOLUME,
                1.0,
                vec![2000.0],
                vec![0.001],
            )
            .unwrap(),
        ))
        .unwrap();
    assert!(matches!(
        original.advance_with_heat(&thermal, Time::new(1.0), 0.0, controls(), || false),
        Err(KineticSolubilityHeatError::Material(_))
    ));
}

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

fn thermal_family(
    thermal: &SolubilityThermalKineticsV1,
    inverse_temperature: f64,
) -> Arc<VirialSolubilityThermalModel> {
    Arc::new(
        VirialSolubilityThermalModel::new(
            thermal.caloric().clone(),
            Arc::new(
                VirialSolutionModel::new(300.0, MOLAR_VOLUME, 1.0, vec![2000.0], vec![1e-5])
                    .unwrap(),
            ),
            vec![inverse_temperature],
            [Temperature::new(280.0), Temperature::new(340.0)],
        )
        .unwrap(),
    )
}

#[test]
fn explicit_thermal_family_survives_temperature_transport_and_source_checks() {
    let thermal = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let family = thermal_family(&thermal, 0.003);
    let fixed = material(&thermal, 0.15, 1.5)
        .with_nonideal_thermodynamics(family.at_temperature(Temperature::new(300.0)).unwrap())
        .unwrap();
    let bound = fixed
        .clone()
        .with_thermal_nonideal_thermodynamics(family.clone())
        .unwrap();
    assert_eq!(bound.last_steps(), fixed.last_steps());
    let warmer = bound.with_temperature(Temperature::new(320.0)).unwrap();
    assert!(Arc::ptr_eq(
        warmer.thermal_nonideal_source().unwrap(),
        &family
    ));
    assert!(std::ptr::eq(warmer.models(), bound.models()));
    assert_eq!(
        warmer.last_steps()[0].state(),
        bound.last_steps()[0].state()
    );
    assert_eq!(
        warmer.solution().nonideal_model(),
        Some(
            family
                .at_temperature(Temperature::new(320.0))
                .unwrap()
                .as_ref()
        )
    );
    assert_ne!(
        warmer.solution().nonideal_model(),
        bound.solution().nonideal_model()
    );
    assert!(bound.with_temperature(Temperature::new(350.0)).is_err());
    assert!(fixed.with_temperature(Temperature::new(320.0)).is_err());
    let restricted = Arc::new(
        VirialSolubilityThermalModel::new(
            thermal.caloric().clone(),
            family.reference_model().clone(),
            family.inverse_temperature_m3_k_mol().to_vec(),
            [Temperature::new(290.0), Temperature::new(330.0)],
        )
        .unwrap(),
    );
    let restricted = material(&thermal, 0.15, 1.5)
        .with_thermal_nonideal_thermodynamics(restricted)
        .unwrap();
    // Inside the original pair-law domain but outside the thermal fit.
    assert!(
        restricted
            .with_temperature(Temperature::new(335.0))
            .is_err()
    );
    assert!(
        restricted
            .advance_with_heat(&thermal, Time::new(0.1), 0.0, controls(), || false)
            .is_err()
    );

    let mut diluted = warmer.with_solvent_volume_m3(2e-3).unwrap();
    diluted
        .replace_dissolved_inventories(1.5e-3, &[0.1])
        .unwrap();
    assert_eq!(diluted.last_steps()[0].state().solid(), Amount::new(1.5));
    assert!(Arc::ptr_eq(
        diluted.thermal_nonideal_source().unwrap(),
        &family
    ));
    let mut copy = warmer.clone();
    copy.copy_state_from(&diluted).unwrap();
    assert_eq!(copy, diluted);
    assert!(Arc::ptr_eq(
        copy.thermal_nonideal_source().unwrap(),
        &family
    ));

    // Equal current B(T_ref) does not make a different thermal slope/source
    // compatible. Copying also cannot discard a declared thermal family.
    let different = thermal_family(&thermal, 0.006);
    assert!(
        bound
            .clone()
            .with_thermal_nonideal_thermodynamics(different.clone())
            .is_err()
    );
    let other = material(&thermal, 0.15, 1.5)
        .with_thermal_nonideal_thermodynamics(different)
        .unwrap();
    let mut destination = bound.clone();
    assert!(destination.copy_state_from(&other).is_err());
    assert!(destination.copy_state_from(&fixed).is_err());
    assert_eq!(destination, bound);
    let equal = Arc::new((*family).clone());
    let rebound = bound
        .clone()
        .with_thermal_nonideal_thermodynamics(equal)
        .unwrap();
    assert!(Arc::ptr_eq(
        rebound.thermal_nonideal_source().unwrap(),
        &family
    ));
}

#[test]
fn transported_particles_cool_with_original_excess_and_solid_reference_heat() {
    let base = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let thermal = SolubilityThermalKineticsV1::new(
        base.caloric().clone(),
        vec![
            base.kinetics()[0]
                .clone()
                .with_particle_surface(Amount::new(1.5))
                .unwrap(),
        ],
    )
    .unwrap();
    let family = thermal_family(&thermal, 0.003);
    let mut original = material(&thermal, 0.15, 1.5)
        .with_thermal_nonideal_thermodynamics(family.clone())
        .unwrap();
    original
        .replace_dissolved_inventories(7e-4, &[0.03])
        .unwrap();
    let before = original.clone();
    let step = original
        .advance_with_heat(&thermal, Time::new(3.0), 0.0, controls(), || false)
        .unwrap();
    let accepted = step.solution();
    let after = accepted.last_steps()[0].state();
    let temperature = accepted.solution().temperature_k();
    assert!(temperature < 300.0 && after.dissolved().value() > 0.03);
    assert!(after.solid().value() < 1.5);
    assert_eq!(
        after.total_solute(),
        original.last_steps()[0].state().total_solute()
    );
    assert_eq!(
        accepted.models()[0].particle_reference_solid(),
        Some(Amount::new(1.5))
    );
    assert!(Arc::ptr_eq(
        accepted.thermal_nonideal_source().unwrap(),
        &family
    ));
    assert_eq!(
        accepted.solution().nonideal_model(),
        Some(
            family
                .at_temperature(Temperature::new(temperature))
                .unwrap()
                .as_ref()
        )
    );
    assert_eq!(
        accepted.solution().solute_moles()[0],
        after.dissolved().value()
    );
    assert_eq!(original, before);

    // Independent H_extra for B(T)=B_ref+C(1/T-1/T_ref), retaining BOTH
    // the quadratic mixing term and the pure-solid reference correction.
    let extra = |t: f64, dissolved: f64| {
        let c = 0.003;
        let b = 1e-5 + c * (1.0 / t - 1.0 / 300.0);
        let saturation = 600.0 * (DISSOLUTION_ENTHALPY / R * (1.0 / 300.0 - 1.0 / t)).exp();
        0.5 * R * c * dissolved * dissolved / 7e-4
            + dissolved * saturation * (b * DISSOLUTION_ENTHALPY - R * c)
    };
    let base_heat = (7e-4 / MOLAR_VOLUME * SOLVENT_CP + after.total_solute().value() * SOLUTE_CP)
        * (temperature - 300.0)
        + (after.dissolved().value() - 0.03) * DISSOLUTION_ENTHALPY;
    let extra_change = extra(temperature, after.dissolved().value()) - extra(300.0, 0.03);
    assert!(extra_change.abs() > 1e-3);
    near(base_heat + extra_change, 0.0, 5e-7);
    near(
        step.balance().enthalpy_change_j(),
        base_heat + extra_change,
        5e-7,
    );
    assert!(step.balance().energy_residual_j().abs() <= controls().root.energy_tolerance_j);

    // A later heat step continues with this same family and original particle
    // geometry rather than silently falling back to ideal calorics.
    let next = accepted
        .advance_with_heat(&thermal, Time::new(0.1), 5.0, controls(), || false)
        .unwrap();
    assert!(Arc::ptr_eq(
        next.solution().thermal_nonideal_source().unwrap(),
        &family
    ));
    near(
        next.balance().enthalpy_change_j(),
        5.0,
        controls().root.energy_tolerance_j,
    );
}

#[test]
fn nonideal_heat_requires_exact_core_source_identity_and_rolls_back_on_refusal() {
    let thermal = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let family = thermal_family(&thermal, 0.003);
    let original = material(&thermal, 0.15, 1.5)
        .with_thermal_nonideal_thermodynamics(family.clone())
        .unwrap();
    let checkpoint = original.clone();
    let equal_but_distinct: Arc<dyn SolubilityThermalInteractionV1> = Arc::new((*family).clone());
    let foreign = thermal.with_interactions(equal_but_distinct).unwrap();
    assert!(matches!(
        original.advance_with_heat(&foreign, Time::new(1.0), 0.0, controls(), || false),
        Err(KineticSolubilityHeatError::SourceMismatch)
    ));
    let declared: Arc<dyn SolubilityThermalInteractionV1> = family.clone();
    let bound_core = thermal.with_interactions(declared).unwrap();
    assert!(matches!(
        material(&thermal, 0.15, 1.5).advance_with_heat(
            &bound_core,
            Time::new(1.0),
            0.0,
            controls(),
            || false
        ),
        Err(KineticSolubilityHeatError::SourceMismatch)
    ));
    let mut polls = 0;
    assert!(matches!(
        original.advance_with_heat(&bound_core, Time::new(3.0), 0.0, controls(), || {
            polls += 1;
            polls == 16
        }),
        Err(KineticSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::Cancelled
        ))
    ));
    assert_eq!(original, checkpoint);
    let mut limited = controls();
    limited.maximum_attempts = 1;
    assert!(
        original
            .advance_with_heat(&bound_core, Time::new(3.0), 0.0, limited, || false)
            .is_err()
    );
    assert_eq!(original, checkpoint);
    let accepted = original
        .advance_with_heat(&bound_core, Time::new(0.1), 0.0, controls(), || false)
        .unwrap();
    assert!(Arc::ptr_eq(
        accepted.solution().thermal_nonideal_source().unwrap(),
        &family
    ));
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

#[test]
fn reservoir_cooling_tracks_temperature_feedback_and_preserves_closed_phase_stocks() {
    let thermal = model(0.0, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let original = material(&thermal, 0.5, 1e-100);
    let before = original.clone();
    let capacity = 1e-3 / MOLAR_VOLUME * SOLVENT_CP + 0.5 * SOLUTE_CP;
    let duration = 40.0;
    let reservoir = 285.0;
    let conductance = 100.0;
    let exact = reservoir + (300.0 - reservoir) * (-conductance * duration / capacity).exp();
    let run = |maximum_step| {
        let mut control = controls();
        control.initial_step = Time::new(maximum_step);
        control.maximum_step = Time::new(maximum_step);
        control.relative_tolerance = 1e-3;
        original
            .advance_with_heat_transfer(
                &thermal,
                Time::new(duration),
                Temperature::new(reservoir),
                conductance,
                control,
                || false,
            )
            .unwrap()
    };
    let coarse = run(4.0);
    let fine = run(0.25);
    let temperature = fine.solution().solution().temperature_k();
    let coarse_error = (coarse.solution().solution().temperature_k() - exact).abs();
    let fine_error = (temperature - exact).abs();
    assert!(fine_error < 0.03 && 4.0 * fine_error < coarse_error);
    assert!(temperature > reservoir && temperature < 300.0);
    assert_eq!(original, before);
    assert_eq!(
        fine.solution().last_steps()[0].state(),
        original.last_steps()[0].state()
    );
    assert_eq!(
        fine.solution().last_steps()[0].state().solid().value(),
        1e-100
    );
    assert_eq!(
        fine.solution().solution().solute_moles(),
        original.solution().solute_moles()
    );
    let heat = capacity * (temperature - 300.0);
    near(fine.balance().heat_input_j(), heat, 1e-7);
    near(fine.balance().enthalpy_change_j(), heat, 1e-7);
    assert!(fine.balance().energy_residual_j().abs() <= 1e-7);
    // A frozen initial cooling power would remove twice as much heat here.
    assert!(fine.balance().heat_input_j().abs() < 0.75 * conductance * duration * 15.0);
}

#[test]
fn jacket_cooling_grows_original_particles_and_keeps_excess_heat_in_the_balance() {
    let fixed = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let thermal = SolubilityThermalKineticsV1::new(
        fixed.caloric().clone(),
        vec![
            fixed.kinetics()[0]
                .clone()
                .with_particle_surface(Amount::new(0.2))
                .unwrap(),
        ],
    )
    .unwrap();
    let family = thermal_family(&thermal, 0.003);
    let original = material(&thermal, 0.6, 0.2)
        .with_thermal_nonideal_thermodynamics(family.clone())
        .unwrap();
    let before = original.clone();
    let cooled = original
        .advance_with_heat_transfer(
            &thermal,
            Time::new(40.0),
            Temperature::new(285.0),
            80.0,
            controls(),
            || false,
        )
        .unwrap();
    let endpoint = cooled.solution();
    let phase = endpoint.last_steps()[0].state();
    let temperature = endpoint.solution().temperature_k();
    assert!(temperature > 285.0 && temperature < 300.0);
    assert!(phase.dissolved().value() < 0.6 && phase.solid().value() > 0.2);
    assert_eq!(
        phase.total_solute(),
        original.last_steps()[0].state().total_solute()
    );
    assert_eq!(
        endpoint.models()[0].particle_reference_solid(),
        Some(Amount::new(0.2))
    );
    assert!(
        endpoint.models()[0]
            .conductance_at(phase.solid())
            .unwrap()
            .value()
            > 1e-4
    );
    assert!(Arc::ptr_eq(
        endpoint.thermal_nonideal_source().unwrap(),
        &family
    ));
    assert_eq!(
        endpoint.solution().solute_moles()[0],
        phase.dissolved().value()
    );
    assert_eq!(original, before);

    // Independent sensible + dissolution + virial/reference enthalpy; the
    // jacket heat comes from the native accepted boundary powers, not this H.
    let extra = |t: f64, dissolved: f64| {
        let c = 0.003;
        let b = 1e-5 + c * (1.0 / t - 1.0 / 300.0);
        let saturation = 600.0 * (DISSOLUTION_ENTHALPY / R * (1.0 / 300.0 - 1.0 / t)).exp();
        0.5 * R * c * dissolved * dissolved / 1e-3
            + dissolved * saturation * (b * DISSOLUTION_ENTHALPY - R * c)
    };
    let sensible = (1e-3 / MOLAR_VOLUME * SOLVENT_CP + phase.total_solute().value() * SOLUTE_CP)
        * (temperature - 300.0);
    let phase_heat = DISSOLUTION_ENTHALPY * (phase.dissolved().value() - 0.6);
    let extra_change = extra(temperature, phase.dissolved().value()) - extra(300.0, 0.6);
    assert!(phase_heat < -10.0 && extra_change.abs() > 1e-3);
    let heat = sensible + phase_heat + extra_change;
    assert!(cooled.balance().heat_input_j() < 0.0);
    near(cooled.balance().heat_input_j(), heat, 5e-7);
    near(cooled.balance().enthalpy_change_j(), heat, 5e-7);
    assert!(cooled.balance().energy_residual_j().abs() <= controls().root.energy_tolerance_j);

    // A subsequent physical dilution retains the actual crystals and source.
    // Its warmer jacket then dissolves those stocks on the same particle law.
    let diluted = endpoint.with_solvent_volume_m3(1.2e-3).unwrap();
    let warmed = diluted
        .advance_with_heat_transfer(
            &thermal,
            Time::new(1.0),
            Temperature::new(320.0),
            80.0,
            controls(),
            || false,
        )
        .unwrap();
    assert!(warmed.balance().heat_input_j() > 0.0);
    assert!(warmed.solution().solution().solute_moles()[0] > phase.dissolved().value());
    assert_eq!(warmed.solution().solution().solvent_volume_m3(), 1.2e-3);
    assert_eq!(warmed.solution().models(), original.models());
    assert!(Arc::ptr_eq(
        warmed.solution().thermal_nonideal_source().unwrap(),
        &family
    ));
}

#[test]
fn heat_transfer_retains_exact_sources_zero_contact_limit_and_atomic_cancellation() {
    let thermal = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let family = thermal_family(&thermal, 0.003);
    let original = material(&thermal, 0.15, 1.5)
        .with_thermal_nonideal_thermodynamics(family.clone())
        .unwrap();
    let before = original.clone();
    let foreign_source: Arc<dyn SolubilityThermalInteractionV1> = Arc::new((*family).clone());
    let foreign = thermal.with_interactions(foreign_source).unwrap();
    assert!(matches!(
        original.advance_with_heat_transfer(
            &foreign,
            Time::new(0.5),
            Temperature::new(290.0),
            10.0,
            controls(),
            || false,
        ),
        Err(KineticSolubilityHeatError::SourceMismatch)
    ));
    let fixed_only = material(&thermal, 0.15, 1.5)
        .with_nonideal_thermodynamics(family.at_temperature(Temperature::new(300.0)).unwrap())
        .unwrap();
    assert!(matches!(
        fixed_only.advance_with_heat_transfer(
            &thermal,
            Time::new(0.5),
            Temperature::new(290.0),
            10.0,
            controls(),
            || false,
        ),
        Err(KineticSolubilityHeatError::Material(_))
    ));
    let adiabatic = original
        .advance_with_heat(&thermal, Time::new(0.5), 0.0, controls(), || false)
        .unwrap();
    let isolated = original
        .advance_with_heat_transfer(
            &thermal,
            Time::new(0.5),
            Temperature::new(290.0),
            0.0,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(isolated, adiabatic);

    let mut calls = 0;
    let accepted = original
        .advance_with_heat_transfer(
            &thermal,
            Time::new(0.5),
            Temperature::new(290.0),
            10.0,
            controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(calls > 1);
    let mut observed = 0;
    assert!(matches!(
        original.advance_with_heat_transfer(
            &thermal,
            Time::new(0.5),
            Temperature::new(290.0),
            10.0,
            controls(),
            || {
                observed += 1;
                observed == calls
            },
        ),
        Err(KineticSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::Cancelled
        ))
    ));
    assert_eq!(observed, calls);
    assert_eq!(original, before);
    assert_eq!(
        original
            .advance_with_heat_transfer(
                &thermal,
                Time::new(0.5),
                Temperature::new(290.0),
                10.0,
                controls(),
                || false,
            )
            .unwrap(),
        accepted
    );
}

#[test]
fn finite_reservoir_tracks_two_body_cooling_and_retains_heat_across_windows() {
    let thermal = model(0.0, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let original = material(&thermal, 0.5, 1e-100);
    let reservoir = SensibleHeatReservoirV1::new(Temperature::new(280.0), 1000.0).unwrap();
    let capacity = 1e-3 / MOLAR_VOLUME * SOLVENT_CP + 0.5 * SOLUTE_CP;
    let duration = 60.0;
    let conductance = 40.0;
    let total_capacity = capacity + reservoir.heat_capacity_j_k();
    let equilibrium = (capacity * 300.0 + reservoir.heat_capacity_j_k() * 280.0) / total_capacity;
    let difference = 20.0
        * (-conductance * (1.0 / capacity + 1.0 / reservoir.heat_capacity_j_k()) * duration).exp();
    let expected_solution =
        equilibrium + reservoir.heat_capacity_j_k() / total_capacity * difference;
    let expected_reservoir = equilibrium - capacity / total_capacity * difference;
    let mut control = controls();
    control.initial_step = Time::new(0.25);
    control.maximum_step = Time::new(0.25);
    control.relative_tolerance = 1e-3;
    let whole = original
        .advance_with_finite_heat_reservoir(
            &thermal,
            Time::new(duration),
            &reservoir,
            conductance,
            control,
            || false,
        )
        .unwrap();
    let temperature = whole.solution().solution().temperature_k();
    let final_reservoir = whole.balance().reservoir();
    near(temperature, expected_solution, 0.02);
    near(
        final_reservoir.temperature().value(),
        expected_reservoir,
        0.02,
    );
    assert!(temperature > final_reservoir.temperature().value());
    assert!(final_reservoir.temperature().value() > 280.0);
    assert_eq!(
        whole.solution().last_steps()[0].state(),
        original.last_steps()[0].state()
    );
    assert_eq!(
        whole.solution().last_steps()[0].state().solid().value(),
        1e-100
    );
    let solution_heat = capacity * (temperature - 300.0);
    let reservoir_heat =
        reservoir.heat_capacity_j_k() * (final_reservoir.temperature().value() - 280.0);
    near(solution_heat + reservoir_heat, 0.0, 2e-7);
    near(
        whole.balance().balance().heat_input_j(),
        solution_heat,
        1e-7,
    );
    near(
        whole.balance().reservoir_enthalpy_change_j(),
        reservoir_heat,
        1e-7,
    );
    assert!(
        whole.balance().closed_system_energy_residual_j().abs() <= control.root.energy_tolerance_j
    );

    // Pass the ORIGINAL accepted reservoir value onward. Reinitializing from
    // its displayed temperature would discard its compensated energy datum.
    let mut continued = original.clone();
    let mut finite = reservoir.clone();
    let mut delivered = 0.0;
    for _ in 0..6 {
        let next = continued
            .advance_with_finite_heat_reservoir(
                &thermal,
                Time::new(10.0),
                &finite,
                conductance,
                control,
                || false,
            )
            .unwrap();
        let (material, balance) = next.into_parts();
        delivered += balance.balance().heat_input_j();
        finite = balance.reservoir().clone();
        continued = material;
    }
    near(continued.solution().temperature_k(), temperature, 1e-7);
    near(
        finite.temperature().value(),
        final_reservoir.temperature().value(),
        1e-7,
    );
    assert_eq!(
        finite.reference_temperature(),
        reservoir.reference_temperature()
    );
    near(finite.energy_change_j(), -delivered, 2e-7);
    near(
        reservoir.enthalpy_change_j(&finite).unwrap(),
        finite.energy_change_j(),
        1e-10,
    );
    assert_eq!(reservoir.temperature(), Temperature::new(280.0));
    assert_eq!(reservoir.energy_change_j(), 0.0);
    assert_eq!(original.solution().temperature_k(), 300.0);
}

#[test]
fn finite_reservoir_grows_particles_with_the_complete_original_nonideal_heat() {
    let fixed = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let thermal = SolubilityThermalKineticsV1::new(
        fixed.caloric().clone(),
        vec![
            fixed.kinetics()[0]
                .clone()
                .with_particle_surface(Amount::new(0.2))
                .unwrap(),
        ],
    )
    .unwrap();
    let family = thermal_family(&thermal, 0.003);
    let original = material(&thermal, 0.6, 0.2)
        .with_thermal_nonideal_thermodynamics(family.clone())
        .unwrap();
    let reservoir = SensibleHeatReservoirV1::new(Temperature::new(285.0), 2000.0).unwrap();
    let first = original
        .advance_with_finite_heat_reservoir(
            &thermal,
            Time::new(40.0),
            &reservoir,
            80.0,
            controls(),
            || false,
        )
        .unwrap();
    let endpoint = first.solution();
    let phase = endpoint.last_steps()[0].state();
    let temperature = endpoint.solution().temperature_k();
    assert!(temperature < 300.0 && temperature > first.balance().reservoir().temperature().value());
    assert!(first.balance().reservoir().temperature().value() > 285.0);
    assert!(phase.dissolved().value() < 0.6 && phase.solid().value() > 0.2);
    assert_eq!(
        phase.total_solute(),
        original.last_steps()[0].state().total_solute()
    );
    assert_eq!(endpoint.models(), original.models());
    assert!(
        endpoint.models()[0]
            .conductance_at(phase.solid())
            .unwrap()
            .value()
            > 1e-4
    );
    assert!(Arc::ptr_eq(
        endpoint.thermal_nonideal_source().unwrap(),
        &family
    ));
    assert_eq!(
        endpoint.solution().solute_moles()[0],
        phase.dissolved().value()
    );
    assert_eq!(phase, &first.balance().balance().final_state().phases()[0]);

    let extra = |t: f64, dissolved: f64| {
        let c = 0.003;
        let b = 1e-5 + c * (1.0 / t - 1.0 / 300.0);
        let saturation = 600.0 * (DISSOLUTION_ENTHALPY / R * (1.0 / 300.0 - 1.0 / t)).exp();
        0.5 * R * c * dissolved * dissolved / 1e-3
            + dissolved * saturation * (b * DISSOLUTION_ENTHALPY - R * c)
    };
    let solution_heat = (1e-3 / MOLAR_VOLUME * SOLVENT_CP
        + phase.total_solute().value() * SOLUTE_CP)
        * (temperature - 300.0)
        + DISSOLUTION_ENTHALPY * (phase.dissolved().value() - 0.6)
        + extra(temperature, phase.dissolved().value())
        - extra(300.0, 0.6);
    let reservoir_heat = 2000.0 * (first.balance().reservoir().temperature().value() - 285.0);
    near(solution_heat + reservoir_heat, 0.0, 5e-7);
    near(
        first.balance().balance().enthalpy_change_j(),
        solution_heat,
        5e-7,
    );
    near(
        first.balance().reservoir_enthalpy_change_j(),
        reservoir_heat,
        2e-7,
    );
    assert!(
        first.balance().closed_system_energy_residual_j().abs()
            <= controls().root.energy_tolerance_j
    );

    let second = endpoint
        .advance_with_finite_heat_reservoir(
            &thermal,
            Time::new(10.0),
            first.balance().reservoir(),
            80.0,
            controls(),
            || false,
        )
        .unwrap();
    assert!(Arc::ptr_eq(
        second.solution().thermal_nonideal_source().unwrap(),
        &family
    ));
    assert_eq!(second.solution().models(), original.models());
    assert_eq!(
        second.balance().reservoir().reference_temperature(),
        Temperature::new(285.0)
    );
    let all_solution_heat =
        first.balance().balance().heat_input_j() + second.balance().balance().heat_input_j();
    near(
        reservoir
            .enthalpy_change_j(second.balance().reservoir())
            .unwrap(),
        -all_solution_heat,
        2e-7,
    );
    assert_eq!(original.last_steps()[0].state().solid(), Amount::new(0.2));
    assert_eq!(reservoir.energy_change_j(), 0.0);
}

#[test]
fn finite_reservoir_keeps_source_zero_contact_and_late_cancellation_atomic() {
    let thermal = model(1e-4, PrecipitationPolicyV1::ExistingSolidOnly, MOLAR_VOLUME);
    let family = thermal_family(&thermal, 0.003);
    let original = material(&thermal, 0.15, 1.5)
        .with_thermal_nonideal_thermodynamics(family.clone())
        .unwrap();
    let reservoir = SensibleHeatReservoirV1::new(Temperature::new(293.0), 5000.0).unwrap();
    let before = original.clone();
    let reservoir_before = reservoir.clone();
    let foreign_source: Arc<dyn SolubilityThermalInteractionV1> = Arc::new((*family).clone());
    let foreign = thermal.with_interactions(foreign_source).unwrap();
    assert!(matches!(
        original.advance_with_finite_heat_reservoir(
            &foreign,
            Time::new(0.5),
            &reservoir,
            30.0,
            controls(),
            || false,
        ),
        Err(KineticSolubilityHeatError::SourceMismatch)
    ));
    let adiabatic = original
        .advance_with_heat(&thermal, Time::new(0.5), 0.0, controls(), || false)
        .unwrap();
    let isolated = original
        .advance_with_finite_heat_reservoir(
            &thermal,
            Time::new(0.5),
            &reservoir,
            0.0,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(isolated.solution(), adiabatic.solution());
    assert_eq!(isolated.balance().balance(), adiabatic.balance());
    assert_eq!(isolated.balance().reservoir(), &reservoir);
    assert_eq!(isolated.balance().reservoir_enthalpy_change_j(), 0.0);

    let mut calls = 0;
    let accepted = original
        .advance_with_finite_heat_reservoir(
            &thermal,
            Time::new(0.5),
            &reservoir,
            30.0,
            controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(calls > 1);
    let mut observed = 0;
    assert!(matches!(
        original.advance_with_finite_heat_reservoir(
            &thermal,
            Time::new(0.5),
            &reservoir,
            30.0,
            controls(),
            || {
                observed += 1;
                observed == calls
            },
        ),
        Err(KineticSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::Cancelled
        ))
    ));
    assert_eq!(observed, calls);
    assert_eq!(original, before);
    assert_eq!(reservoir, reservoir_before);
    assert_eq!(
        original
            .advance_with_finite_heat_reservoir(
                &thermal,
                Time::new(0.5),
                &reservoir,
                30.0,
                controls(),
                || false,
            )
            .unwrap(),
        accepted
    );
}
