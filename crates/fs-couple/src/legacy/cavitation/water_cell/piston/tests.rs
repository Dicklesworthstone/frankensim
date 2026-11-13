use super::*;
use crate::cavitation::water_cell::WaterCellControls;
use fs_material::phase::water::WaterPhase;
use fs_material::phase::water::acoustics::HomogeneousEquilibriumWater;
use fs_material::phase::water::isochoric::IsochoricWaterControls;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97SaturationPair, If97State};

fn cell_controls(mass: f64, lower: f64, upper: f64) -> WaterCellControls {
    WaterCellControls {
        equilibrium: IsochoricWaterControls {
            minimum_temperature_k: lower,
            maximum_temperature_k: upper,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 100e6,
            volume_tolerance_m3_kg: 1e-13,
            internal_energy_tolerance_j_kg: 1e-7,
            maximum_pressure_iterations: 96,
            maximum_temperature_iterations: 96,
        },
        maximum_relative_mass_exchange: 0.1,
        maximum_inflow_parcels: 4,
        mass_tolerance_kg: mass * 1e-12,
        energy_tolerance_j: mass * 1e-6,
        volume_tolerance_m3: mass * 1e-13,
    }
}

fn mixed(temperature: f64, quality: f64) -> EquilibriumWaterCell {
    let mass = 0.05;
    let pair = If97SaturationPair::at_temperature(Temperature::new(temperature)).unwrap();
    let volume = mass
        * ((1.0 - quality) * pair.liquid().specific_volume().value()
            + quality * pair.vapor().specific_volume().value());
    EquilibriumWaterCell::at_temperature(
        volume,
        mass,
        Temperature::new(temperature),
        0.0,
        cell_controls(mass, temperature - 20.0, temperature + 20.0),
        || false,
    )
    .unwrap()
}

fn controls(piston: &WaterPiston, relative_bracket: f64) -> WaterPistonControls {
    let cell = piston.water();
    let volume = cell.state().volume_m3();
    WaterPistonControls {
        water: WaterVolumeControls {
            maximum_step_s: 1.0,
            maximum_relative_volume_change: 0.1,
            maximum_temperature_increment_k: 2.0,
            energy_tolerance_j: cell.controls().energy_tolerance_j,
            maximum_iterations: 96,
        },
        volume_bracket_m3: [
            volume * (1.0 - relative_bracket),
            volume * (1.0 + relative_bracket),
        ],
        momentum_tolerance_n_s: 1e-11,
        energy_tolerance_j: 10.0 * cell.controls().energy_tolerance_j,
        maximum_iterations: 96,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn stationary_and_fixed_support_retain_the_original_water_heat_owner() {
    let cell = mixed(400.0, 0.02);
    let pressure = cell.state().source_state().pressure_pa();
    let input = WaterPistonInput {
        heat_j: 0.0,
        exterior_pressure_pa: pressure,
        external_force_n: 0.0,
    };
    let mut free = WaterPiston::new(cell, 1.0, 0.01, 0.0).unwrap();
    let before = free;
    let c = controls(&free, 0.01);
    let stationary = free.step(0.01, input, c, || false).unwrap();
    assert_eq!(free.water().state().volume_m3(), cell.state().volume_m3());
    assert_eq!(
        free.water().state().internal_energy_j(),
        cell.state().internal_energy_j()
    );
    assert_eq!(
        free.water().state().phase_masses_kg(),
        cell.state().phase_masses_kg()
    );
    assert_eq!(free.velocity_m_s(), 0.0);
    assert_eq!(stationary.water.water.boundary_work_j, 0.0);
    assert_eq!(stationary.total_energy_residual_j, 0.0);
    assert_eq!(stationary.mechanical_trials, 1);
    assert_eq!(free.water().state().time_s(), 0.01);

    let mut fixed = before.with_fixed_support().unwrap();
    let mut expected = cell;
    let heat = 10.0;
    let expected_step = expected.advance(&[], 0.0, heat, 0.1, || false).unwrap();
    let input = WaterPistonInput {
        heat_j: heat,
        external_force_n: 7.0,
        ..input
    };
    let mut c = controls(&fixed, 0.0);
    c.volume_bracket_m3 = [cell.state().volume_m3(); 2];
    let step = fixed.step(0.1, input, c, || false).unwrap();
    assert_eq!(fixed.water(), &expected);
    assert_eq!(step.water.water, expected_step);
    assert_eq!(fixed.velocity_m_s(), 0.0);
    assert_eq!(step.displacement_m, 0.0);
    assert_eq!(step.support_work_j, 0.0);
    near(
        step.support_impulse_n_s,
        -0.1 * (fixed.area_m2() * (step.water.mean_boundary_pressure_pa - pressure) + 7.0),
        1e-14,
    );
    near(step.total_energy_residual_j, 0.0, c.energy_tolerance_j);
}

#[test]
fn isothermal_evaporation_and_condensation_supply_exact_reciprocal_work() {
    for sign in [-1.0, 1.0] {
        let cell = mixed(400.0, 0.01);
        let initial = cell.state();
        let pair = initial.source_state().saturation_pair().unwrap();
        let delta_mass = sign * initial.mass_kg() * 0.0005;
        let delta_volume = delta_mass
            * (pair.vapor().specific_volume().value() - pair.liquid().specific_volume().value());
        let heat = delta_mass * pair.latent_enthalpy().value();
        let duration = 0.1;
        let area = 0.01;
        let speed = delta_volume / (area * duration);
        let mut piston = WaterPiston::new(cell, 2.0, area, speed).unwrap();
        let c = controls(&piston, 0.09);
        let input = WaterPistonInput {
            heat_j: heat,
            exterior_pressure_pa: pair.pressure().value(),
            external_force_n: 0.0,
        };
        let step = piston.step(duration, input, c, || false).unwrap();
        near(
            piston.water().state().source_state().temperature_k(),
            400.0,
            2e-7,
        );
        near(
            piston.water().state().volume_m3() - initial.volume_m3(),
            delta_volume,
            1e-12,
        );
        near(piston.velocity_m_s(), speed, 2e-9);
        near(step.vapor_mass_change_kg, delta_mass, 2e-12);
        near(
            step.water.water.boundary_work_j,
            pair.pressure().value() * delta_volume,
            2e-7,
        );
        near(
            step.exterior_pressure_work_j,
            step.water.water.boundary_work_j,
            2e-7,
        );
        near(step.kinetic_energy_change_j, 0.0, 1e-9);
        near(
            step.water.water.internal_energy_change_j,
            delta_mass
                * (pair.vapor().specific_internal_energy().value()
                    - pair.liquid().specific_internal_energy().value()),
            2e-7,
        );
        near(step.total_energy_residual_j, 0.0, c.energy_tolerance_j);
        near(
            step.source_total_energy_residual_j,
            0.0,
            c.energy_tolerance_j,
        );
        assert_eq!(piston.water().state().mass_kg(), initial.mass_kg());
    }
}

#[test]
fn pressure_load_creates_and_collapses_vapor_without_prescribing_the_motion() {
    let mass = 0.001;
    let source = If97State::at_pressure_temperature(
        Pressure::new(120_000.0),
        Temperature::new(373.15),
        If97Region::Region1,
    )
    .unwrap();
    let cell = EquilibriumWaterCell::new(
        mass * source.specific_volume().value(),
        mass,
        mass * source.specific_internal_energy().value(),
        0.0,
        cell_controls(mass, 371.15, 375.15),
        || false,
    )
    .unwrap();
    assert_eq!(
        cell.state().source_state().phase(),
        WaterPhase::CompressedLiquid
    );
    let mut piston = WaterPiston::new(cell, 0.1, 1e-4, 0.0).unwrap();
    let initial = piston;
    let duration = 2e-6;
    let mut external_work = 0.0;
    let mut mechanical_work = 0.0;
    let mut largest_vapor_mass = 0.0_f64;
    let mut condensed = false;
    let mut expanded = None;
    for index in 1..=300 {
        let input = WaterPistonInput {
            heat_j: 0.0,
            exterior_pressure_pa: if index <= 100 { 80_000.0 } else { 140_000.0 },
            external_force_n: 0.0,
        };
        let c = controls(&piston, 0.001);
        let step = piston
            .step(f64::from(index) * duration, input, c, || false)
            .unwrap();
        external_work += step.exterior_pressure_work_j;
        mechanical_work += step.water.water.boundary_work_j;
        largest_vapor_mass = largest_vapor_mass.max(piston.water().state().phase_masses_kg()[1]);
        condensed |= step.vapor_mass_change_kg < 0.0;
        near(step.mechanical_energy_residual_j, 0.0, c.energy_tolerance_j);
        near(
            step.source_total_energy_residual_j,
            0.0,
            c.energy_tolerance_j,
        );
        assert_eq!(piston.water().state().mass_kg(), mass);
        if index == 100 {
            expanded = Some(piston);
        }
    }
    let expanded = expanded.unwrap();
    assert!(expanded.water().state().volume_m3() > initial.water().state().volume_m3());
    assert!(expanded.velocity_m_s() > 0.0);
    assert_eq!(
        expanded.water().state().source_state().phase(),
        WaterPhase::LiquidVapor
    );
    assert!(largest_vapor_mass > 1e-12);
    assert!(condensed);
    assert_eq!(
        piston.water().state().source_state().phase(),
        WaterPhase::CompressedLiquid
    );
    assert_eq!(piston.water().state().phase_masses_kg()[1], 0.0);
    near(
        piston.kinetic_energy_j() - initial.kinetic_energy_j() - mechanical_work + external_work,
        0.0,
        2e-9,
    );
    near(
        piston.water().state().internal_energy_j() - initial.water().state().internal_energy_j()
            + piston.kinetic_energy_j()
            - initial.kinetic_energy_j()
            + external_work,
        0.0,
        2e-9,
    );
}

#[test]
fn oscillation_timing_converges_to_the_original_equilibrium_acoustic_tangent() {
    let cell = mixed(400.0, 0.02);
    let initial = cell.state();
    let pair = initial.source_state().saturation_pair().unwrap();
    let equilibrium = HomogeneousEquilibriumWater::try_new(pair, 0.02).unwrap();
    let bulk = equilibrium.state().adiabatic_bulk_modulus_pa();
    let area = 0.01;
    let wall_mass = 1.0;
    let omega = (bulk * area * area / (wall_mass * initial.volume_m3())).sqrt();
    let amplitude = 1e-5 * initial.volume_m3();
    let speed = amplitude * omega / area;
    let duration = 1.0 / omega;
    let input = WaterPistonInput {
        heat_j: 0.0,
        exterior_pressure_pa: initial.source_state().pressure_pa(),
        external_force_n: 0.0,
    };
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let mut piston = WaterPiston::new(cell, wall_mass, area, speed).unwrap();
        for index in 1..=count {
            let mut c = controls(&piston, 1e-4);
            c.momentum_tolerance_n_s = 1e-12;
            piston
                .step(
                    duration * f64::from(index) / f64::from(count),
                    input,
                    c,
                    || false,
                )
                .unwrap();
        }
        let displacement = (piston.water().state().volume_m3() - initial.volume_m3()) / amplitude;
        let velocity = piston.velocity_m_s() / speed;
        let discrete_angle = 2.0 * f64::from(count) * (0.5 / f64::from(count)).atan();
        near(displacement, discrete_angle.sin(), 1e-4);
        near(velocity, discrete_angle.cos(), 1e-4);
        errors.push((displacement - 1.0_f64.sin()).abs() + (velocity - 1.0_f64.cos()).abs());
    }
    assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.4 * errors[1], "{errors:?}");
}

#[test]
fn final_cancellation_unbracketed_motion_and_dropped_proposals_are_atomic() {
    let cell = mixed(400.0, 0.02);
    let pressure = cell.state().source_state().pressure_pa();
    let mut piston = WaterPiston::new(cell, 1.0, 0.01, 0.001).unwrap();
    let initial = piston;
    let input = WaterPistonInput {
        heat_j: 0.0,
        exterior_pressure_pa: pressure,
        external_force_n: 0.0,
    };
    let c = controls(&piston, 0.001);
    let mut polls = 0;
    {
        let proposal = piston
            .propose_step(0.001, input, c, || {
                polls += 1;
                false
            })
            .unwrap();
        assert!(proposal.endpoint().position_m() > initial.position_m());
        assert_eq!(proposal.endpoint().water().state().time_s(), 0.001);
    }
    assert_eq!(piston, initial);
    let mut visited = 0;
    assert!(matches!(
        piston.propose_step(0.001, input, c, || {
            visited += 1;
            visited == polls
        }),
        Err(WaterPistonError::Cancelled)
    ));
    assert_eq!(piston, initial);
    assert_eq!(visited, polls);
    assert!(piston.with_fixed_support().is_err());

    let mut still = WaterPiston::new(cell, 1.0, 0.01, 0.0).unwrap();
    let before = still;
    let narrow = controls(&still, 1e-9);
    let input = WaterPistonInput {
        external_force_n: 100.0,
        ..input
    };
    assert!(matches!(
        still.step(0.01, input, narrow, || false),
        Err(WaterPistonError::OutsideBracket { .. })
    ));
    assert_eq!(still, before);
}
