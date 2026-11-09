//! One resolved liquid/phase trajectory supplies the actual acoustic source.
//! These check numerical coupling and conservation within the admitted model.

#![cfg(feature = "cavitation")]

use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits};
use fs_couple::cavitation::bubble::vapor::conducting::{
    MaterialConductingVaporBubble, MaterialConductingVaporBubbleError,
};
use fs_flux::cavitation::vapor::conducting::{
    ConductingVaporBubble, ConductingVaporBubbleError, ConductingVaporBubbleNumerics,
};
use fs_flux::cavitation::vapor::liquid::{
    LiquidThermalLimits, SphericalLiquidField, SphericalLiquidGrid,
};
use fs_flux::cavitation::vapor::{
    SaturatedVaporBubble, VaporBubbleLimits, VaporBubbleNumerics, VaporBubbleState,
};
use fs_flux::cavitation::{BubbleLimits, BubbleNumerics};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

const RADIUS: f64 = 100e-6;
const TEMPERATURE: f64 = 373.15;
const END_TIME: f64 = 2e-6;
const CELLS: usize = 16;

fn water() -> SaturatedWaterCavitation {
    let pair = If97SaturationPair::at_temperature(Temperature::new(TEMPERATURE)).unwrap();
    SaturatedWaterCavitation::try_new(
        pair,
        pair.surface_tension().unwrap().surface_tension().value(),
    )
    .unwrap()
}
fn initial() -> VaporBubbleState {
    VaporBubbleState::try_new(0.0, RADIUS, 0.0, TEMPERATURE).unwrap()
}
fn limits() -> VaporBubbleLimits {
    VaporBubbleLimits {
        bubble: BubbleLimits {
            minimum_radius_m: 0.8 * RADIUS,
            maximum_radius_m: 1.3 * RADIUS,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
        minimum_temperature_k: TEMPERATURE - 2.0,
        maximum_temperature_k: TEMPERATURE + 2.0,
        maximum_pressure_ratio: 0.001,
        maximum_vapor_mach_number: 0.01,
        maximum_frozen_liquid_relative_change: 0.02,
        maximum_neglected_enthalpy_fraction: 0.001,
    }
}
fn field(water: SaturatedWaterCavitation) -> SphericalLiquidField {
    SphericalLiquidField::new(
        water,
        SphericalLiquidGrid::equal_volume(RADIUS, RADIUS + 30e-6, CELLS).unwrap(),
        LiquidThermalLimits {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            maximum_source_relative_change: 0.02,
            maximum_pressure_path_capacity_fraction: 0.01,
            maximum_pressure_heating_temperature_k: 0.01,
            maximum_outflow_peclet: 0.5,
        },
    )
    .unwrap()
}
fn numerics() -> ConductingVaporBubbleNumerics {
    ConductingVaporBubbleNumerics {
        bubble: VaporBubbleNumerics {
            radial: BubbleNumerics {
                relative_tolerance: 1e-9,
                radius_absolute_tolerance_m: 1e-14,
                velocity_absolute_tolerance_m_s: 1e-9,
                initial_step_s: 1e-8,
                minimum_step_s: 1e-13,
                maximum_step_s: 5e-8,
                maximum_attempts: 20_000,
            },
            temperature_absolute_tolerance_k: 1e-9,
            energy_absolute_tolerance_j: 1e-16,
            mass_absolute_tolerance_kg: 1e-22,
            energy_balance_tolerance_j: 1e-13,
            mass_balance_tolerance_kg: 1e-19,
        },
        liquid_temperature_absolute_tolerance_k: 1e-9,
        liquid_energy_balance_tolerance_j: 1e-13,
        total_energy_balance_tolerance_j: 1e-13,
    }
}
fn equilibrium(water: &SaturatedWaterCavitation) -> f64 {
    water.state().saturation_pressure_pa() - 2.0 * water.state().surface_tension_n_m() / RADIUS
}
fn bubble(
    water: &SaturatedWaterCavitation,
    liquid_temperature: f64,
) -> MaterialConductingVaporBubble<'_> {
    MaterialConductingVaporBubble::new(
        water,
        initial(),
        limits(),
        field(*water),
        &[liquid_temperature; CELLS],
    )
    .unwrap()
}
fn observer_limits() -> MonopoleLimits {
    MonopoleLimits {
        maximum_ka: 0.1,
        minimum_distance_radii: 10.0,
        maximum_pressure_ratio: 0.01,
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-13 * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn native_temperature_field_phase_history_and_causal_sound_share_one_trajectory() {
    let water = water();
    let p = equilibrium(&water) - 1000.0;
    let mut direct = ConductingVaporBubble::new(
        SaturatedVaporBubble::new(water, initial(), limits()).unwrap(),
        field(water),
        &[TEMPERATURE; CELLS],
    )
    .unwrap();
    let mut coupled = bubble(&water, TEMPERATURE);
    let observer = coupled
        .monopole_observer(0.1, 100_000.0, observer_limits())
        .unwrap();
    let expected = direct
        .advance_controlled(END_TIME, &numerics(), &|_| p, &|_| TEMPERATURE, &mut || {
            false
        })
        .unwrap();
    let actual = coupled
        .advance_observed(
            END_TIME,
            &numerics(),
            &|_| p,
            &|_| TEMPERATURE,
            &observer,
            &mut || false,
        )
        .unwrap();
    assert_eq!(actual.bubble, expected);
    assert_eq!(coupled.state(), direct.state());
    assert_eq!(
        coupled.liquid_temperatures_k(),
        direct.liquid_temperatures_k()
    );
    assert_eq!(coupled.source(), &water);
    assert_eq!(
        actual.acoustic.volume_acceleration_m3_s2,
        expected.bubble.observation.liquid_volume_flow_rate_m3_s2
    );
    near(
        actual.acoustic.arrival_time_s,
        END_TIME + 0.1 / water.state().liquid_sound_speed_m_s(),
    );
    near(
        actual.acoustic.pressure_pa,
        water.state().liquid_density_kg_m3()
            * expected.bubble.observation.liquid_volume_flow_rate_m3_s2
            / (4.0 * core::f64::consts::PI * 0.1),
    );
    assert!(expected.total_energy_residual_j.abs() <= numerics().total_energy_balance_tolerance_j);
}

#[test]
fn finite_liquid_heat_changes_phase_inventory_and_sound_without_a_heat_coefficient() {
    let water = water();
    let p = equilibrium(&water);
    let mut cold = bubble(&water, TEMPERATURE);
    let mut warm = bubble(&water, TEMPERATURE + 0.2);
    let cold_step = cold
        .advance_compact_source(END_TIME, &numerics(), &|_| p, &|_| TEMPERATURE, &mut || {
            false
        })
        .unwrap();
    let warm_step = warm
        .advance_compact_source(END_TIME, &numerics(), &|_| p, &|_| TEMPERATURE, &mut || {
            false
        })
        .unwrap();
    let phase = warm_step.bubble.bubble;
    assert!(phase.supplied_heat_j > 0.0);
    assert!(phase.observation.vapor_mass_kg > cold_step.bubble.bubble.observation.vapor_mass_kg);
    assert!(warm.state().temperature_k() > cold.state().temperature_k());
    assert!(warm.liquid_temperatures_k()[0] < TEMPERATURE + 0.2);
    assert_ne!(
        warm_step.source.strength.volume_flow_rate_m3_s2,
        cold_step.source.strength.volume_flow_rate_m3_s2
    );
    assert_eq!(
        warm_step.source,
        warm.compact_source_sample(p, TEMPERATURE).unwrap()
    );
    assert!(
        warm_step.bubble.liquid_energy_residual_j.abs()
            <= numerics().liquid_energy_balance_tolerance_j
    );
    assert!(
        warm_step.bubble.total_energy_residual_j.abs()
            <= numerics().total_energy_balance_tolerance_j
    );
}

#[test]
fn wrong_medium_and_late_observer_failure_preserve_every_liquid_cell() {
    let water = water();
    let mut bubble = bubble(&water, TEMPERATURE);
    let wrong = CompactMonopoleObserver::new(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        0.1,
        100_000.0,
        observer_limits(),
    )
    .unwrap();
    let result = bubble.advance_observed(
        END_TIME,
        &numerics(),
        &|_| equilibrium(&water) - 1000.0,
        &|_| TEMPERATURE,
        &wrong,
        &mut || false,
    );
    assert_eq!(
        result,
        Err(MaterialConductingVaporBubbleError::ObserverMediumMismatch)
    );
    assert_eq!(bubble.state(), initial());
    assert_eq!(bubble.liquid_temperatures_k(), &[TEMPERATURE; CELLS]);
    let close = bubble
        .monopole_observer(10.0 * RADIUS, 100_000.0, observer_limits())
        .unwrap();
    let result = bubble.advance_observed(
        END_TIME,
        &numerics(),
        &|_| equilibrium(&water) - 1000.0,
        &|_| TEMPERATURE,
        &close,
        &mut || false,
    );
    assert!(matches!(
        result,
        Err(MaterialConductingVaporBubbleError::Acoustic(_))
    ));
    assert_eq!(bubble.state(), initial());
    assert_eq!(bubble.liquid_temperatures_k(), &[TEMPERATURE; CELLS]);
}

#[test]
fn final_cancellation_discards_phase_and_the_entire_resolved_liquid_history() {
    let water = water();
    let p = equilibrium(&water) - 1000.0;
    let mut reference = bubble(&water, TEMPERATURE);
    let mut polls = 0;
    reference
        .advance_compact_source(END_TIME, &numerics(), &|_| p, &|_| TEMPERATURE, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut candidate = bubble(&water, TEMPERATURE);
    let mut calls = 0;
    let result = candidate.advance_compact_source(
        END_TIME,
        &numerics(),
        &|_| p,
        &|_| TEMPERATURE,
        &mut || {
            calls += 1;
            calls == polls
        },
    );
    assert_eq!(
        result,
        Err(MaterialConductingVaporBubbleError::Dynamics(
            ConductingVaporBubbleError::Cancelled
        ))
    );
    assert_eq!(candidate.state(), initial());
    assert_eq!(candidate.liquid_temperatures_k(), &[TEMPERATURE; CELLS]);
}
