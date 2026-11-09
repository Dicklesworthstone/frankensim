//! Phase mass, liquid-displacement acoustics and atomic observation on the
//! original saturated-water trajectory. These are numerical model checks,
//! not an experimental validation of rapid phase equilibrium or heat transfer.

#![cfg(feature = "cavitation")]

use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits};
use fs_couple::cavitation::bubble::vapor::{MaterialVaporBubble, MaterialVaporBubbleError};
use fs_flux::cavitation::vapor::{
    SaturatedVaporBubble, VaporBubbleError, VaporBubbleLimits, VaporBubbleNumerics,
    VaporBubbleState, VaporBubbleThermalState,
};
use fs_flux::cavitation::{BubbleLimits, BubbleNumerics};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

const RADIUS: f64 = 100e-6;
const TEMPERATURE: f64 = 373.15;
const END_TIME: f64 = 2e-6;

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
        maximum_pressure_ratio: 0.01,
        maximum_vapor_mach_number: 0.01,
        maximum_frozen_liquid_relative_change: 0.02,
        maximum_neglected_enthalpy_fraction: 0.001,
    }
}
fn numerics() -> VaporBubbleNumerics {
    VaporBubbleNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-14,
            velocity_absolute_tolerance_m_s: 1e-9,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-13,
            maximum_step_s: 5e-8,
            maximum_attempts: 10_000,
        },
        temperature_absolute_tolerance_k: 1e-9,
        energy_absolute_tolerance_j: 1e-16,
        mass_absolute_tolerance_kg: 1e-22,
        energy_balance_tolerance_j: 1e-12,
        mass_balance_tolerance_kg: 1e-19,
    }
}
fn pressure(water: &SaturatedWaterCavitation) -> f64 {
    water.state().saturation_pressure_pa()
        - 2.0 * water.state().surface_tension_n_m() / RADIUS
        - 1000.0
}
fn no_heat(_: f64, _: &VaporBubbleThermalState) -> Result<f64, VaporBubbleError> {
    Ok(0.0)
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
fn native_phase_history_and_causal_liquid_sound_remain_one_trajectory() {
    let water = water();
    let p = pressure(&water);
    let mut direct = SaturatedVaporBubble::new(water, initial(), limits()).unwrap();
    let mut coupled = MaterialVaporBubble::new(&water, initial(), limits()).unwrap();
    let observer = coupled
        .monopole_observer(0.1, 100_000.0, observer_limits())
        .unwrap();
    let expected = direct
        .advance_controlled(END_TIME, &numerics(), &|_| p, &no_heat, &mut || false)
        .unwrap();
    let actual = coupled
        .advance_observed(
            END_TIME,
            &numerics(),
            &|_| p,
            &no_heat,
            &observer,
            &mut || false,
        )
        .unwrap();
    assert_eq!(actual.bubble, expected);
    assert_eq!(coupled.state(), direct.state());
    assert_eq!(coupled.source(), &water);
    assert_eq!(
        actual.acoustic.volume_acceleration_m3_s2,
        expected.observation.liquid_volume_flow_rate_m3_s2
    );
    near(
        actual.acoustic.arrival_time_s,
        END_TIME + 0.1 / water.state().liquid_sound_speed_m_s(),
    );
    near(
        actual.acoustic.pressure_pa,
        water.state().liquid_density_kg_m3() * expected.observation.liquid_volume_flow_rate_m3_s2
            / (4.0 * core::f64::consts::PI * 0.1),
    );
    assert!(expected.mass_residual_kg.abs() <= numerics().mass_balance_tolerance_kg);
    assert!(expected.energy_residual_j.abs() <= numerics().energy_balance_tolerance_j);
}

#[test]
fn phase_exchange_uses_liquid_displacement_instead_of_geometric_volume_rate() {
    let water = water();
    let bubble = MaterialVaporBubble::new(&water, initial(), limits()).unwrap();
    let observation = bubble.observe(pressure(&water), 0.001).unwrap();
    let source = bubble
        .compact_source_sample(pressure(&water), 0.001)
        .unwrap();
    assert!(observation.vapor_mass_rate_kg_s > 0.0);
    assert!(observation.volume_rate_m3_s > 0.0);
    assert_eq!(source.strength.volume_flow_m3_s, 0.0);
    near(
        observation.volume_rate_m3_s - source.strength.volume_flow_m3_s,
        observation.vapor_mass_rate_kg_s / water.state().liquid_density_kg_m3(),
    );
    assert_eq!(
        source.strength.volume_flow_rate_m3_s2,
        observation.liquid_volume_flow_rate_m3_s2
    );
    assert_eq!(source.strength.force_on_fluid_n, [0.0; 3]);
    assert_eq!(bubble.state(), initial());
}

#[test]
fn heat_changes_phase_mass_and_the_computed_acoustic_source() {
    let water = water();
    let p = pressure(&water);
    let mut adiabatic = MaterialVaporBubble::new(&water, initial(), limits()).unwrap();
    let mut heated = MaterialVaporBubble::new(&water, initial(), limits()).unwrap();
    let cold = adiabatic
        .advance_compact_source(END_TIME, &numerics(), &|_| p, &no_heat, &mut || false)
        .unwrap();
    let hot = heated
        .advance_compact_source(
            END_TIME,
            &numerics(),
            &|_| p,
            &|_, _| Ok(0.001),
            &mut || false,
        )
        .unwrap();
    assert!(
        hot.bubble.observation.state.temperature_k()
            > cold.bubble.observation.state.temperature_k()
    );
    assert!(hot.bubble.observation.vapor_mass_kg > cold.bubble.observation.vapor_mass_kg);
    assert_ne!(
        hot.source.strength.volume_flow_rate_m3_s2,
        cold.source.strength.volume_flow_rate_m3_s2
    );
    near(hot.bubble.supplied_heat_j, 0.001 * END_TIME);
    assert_eq!(hot.source, heated.compact_source_sample(p, 0.001).unwrap());
}

#[test]
fn wrong_medium_and_late_observer_refusal_preserve_phase_and_motion() {
    let water = water();
    let mut bubble = MaterialVaporBubble::new(&water, initial(), limits()).unwrap();
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
        &|_| pressure(&water),
        &no_heat,
        &wrong,
        &mut || false,
    );
    assert_eq!(
        result,
        Err(MaterialVaporBubbleError::ObserverMediumMismatch)
    );
    assert_eq!(bubble.state(), initial());
    let close = bubble
        .monopole_observer(10.0 * RADIUS, 100_000.0, observer_limits())
        .unwrap();
    let result = bubble.advance_observed(
        END_TIME,
        &numerics(),
        &|_| pressure(&water),
        &no_heat,
        &close,
        &mut || false,
    );
    assert!(matches!(result, Err(MaterialVaporBubbleError::Acoustic(_))));
    assert_eq!(bubble.state(), initial());
}

#[test]
fn cancellation_after_source_admission_discards_the_complete_physical_window() {
    let water = water();
    let mut reference = MaterialVaporBubble::new(&water, initial(), limits()).unwrap();
    let mut polls = 0;
    reference
        .advance_compact_source(
            END_TIME,
            &numerics(),
            &|_| pressure(&water),
            &no_heat,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut bubble = MaterialVaporBubble::new(&water, initial(), limits()).unwrap();
    let mut calls = 0;
    let result = bubble.advance_compact_source(
        END_TIME,
        &numerics(),
        &|_| pressure(&water),
        &no_heat,
        &mut || {
            calls += 1;
            calls == polls
        },
    );
    assert_eq!(
        result,
        Err(MaterialVaporBubbleError::Dynamics(
            VaporBubbleError::Cancelled
        ))
    );
    assert_eq!(bubble.state(), initial());
}
