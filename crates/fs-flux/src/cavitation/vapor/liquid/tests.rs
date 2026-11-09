use super::*;
use crate::cavitation::BubbleLimits;
use crate::cavitation::vapor::{SaturatedVaporBubble, VaporBubbleLimits};
use fs_material::cavitation::vapor::SaturatedWaterVapor;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

const RADIUS: f64 = 1.0e-4;
const TEMPERATURE: f64 = 373.15;
const CELLS: usize = 8;

fn source() -> SaturatedWaterCavitation {
    let pair = If97SaturationPair::at_temperature(Temperature::new(TEMPERATURE)).unwrap();
    SaturatedWaterCavitation::try_new(
        pair,
        pair.surface_tension().unwrap().surface_tension().value(),
    )
    .unwrap()
}

fn limits() -> LiquidThermalLimits {
    LiquidThermalLimits {
        minimum_temperature_k: 371.15,
        maximum_temperature_k: 375.15,
        maximum_source_relative_change: 0.02,
        maximum_pressure_path_capacity_fraction: 0.01,
        maximum_pressure_heating_temperature_k: 0.01,
        maximum_outflow_peclet: 0.5,
    }
}

fn field() -> SphericalLiquidField {
    SphericalLiquidField::new(
        source(),
        SphericalLiquidGrid::equal_volume(RADIUS, 1.5 * RADIUS, CELLS).unwrap(),
        limits(),
    )
    .unwrap()
}

fn state(temperature: f64, velocity: f64) -> VaporBubbleState {
    VaporBubbleState::try_new(0.0, RADIUS, velocity, temperature).unwrap()
}

fn thermal(state: VaporBubbleState) -> VaporBubbleThermalState {
    let vapor = SaturatedWaterVapor::at_temperature_k(state.temperature_k()).unwrap();
    let volume = (4.0 * PI / 3.0) * state.radius_m().powi(3);
    VaporBubbleThermalState {
        state,
        vapor,
        volume_m3: volume,
        vapor_mass_kg: vapor.vapor_density_kg_m3() * volume,
    }
}

fn bubble(state: VaporBubbleState) -> SaturatedVaporBubble {
    SaturatedVaporBubble::new(
        source(),
        state,
        VaporBubbleLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.5 * RADIUS,
                maximum_radius_m: 2.0 * RADIUS,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            minimum_temperature_k: 371.15,
            maximum_temperature_k: 375.15,
            maximum_pressure_ratio: 0.001,
            maximum_vapor_mach_number: 0.05,
            maximum_frozen_liquid_relative_change: 0.02,
            maximum_neglected_enthalpy_fraction: 0.001,
        },
    )
    .unwrap()
}

fn pressure() -> f64 {
    let source = source();
    let point = source.state();
    point.saturation_pressure_pa() - 2.0 * point.surface_tension_n_m() / RADIUS
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn volume_coordinates_preserve_each_liquid_inventory_during_radial_motion() {
    let grid = SphericalLiquidGrid::equal_volume(RADIUS, 1.7 * RADIUS, CELLS).unwrap();
    for radius in [0.6 * RADIUS, RADIUS, 1.6 * RADIUS] {
        let geometry = grid.geometry(radius).unwrap();
        assert_eq!(geometry.faces[0], radius);
        for i in 0..CELLS {
            let inner = geometry.faces[i];
            let outer = geometry.faces[i + 1];
            let volume = (4.0 * PI / 3.0) * (outer.powi(3) - inner.powi(3));
            near(volume, grid.volumes_m3()[i], 1.0e-13, 0.0);
            assert!(geometry.centres[i] > inner);
            assert!(geometry.centres[i] < outer);
        }
    }
    assert!(SphericalLiquidGrid::from_face_radii(&[RADIUS, RADIUS]).is_err());
    assert!(SphericalLiquidGrid::equal_volume(RADIUS, 2.0 * RADIUS, 0).is_err());
    assert!(SphericalLiquidGrid::equal_volume(RADIUS, 2.0 * RADIUS, 65).is_err());
    assert!(SphericalLiquidGrid::equal_volume(RADIUS, 2.0 * RADIUS, 64).is_ok());
}

#[test]
fn resolved_conduction_recovers_the_spherical_laplace_profile() {
    let field = field();
    let geometry = field.grid.geometry(RADIUS).unwrap();
    let outer = geometry.faces[CELLS];
    let increment = 1.0e-3;
    let interface = thermal(state(TEMPERATURE + increment, 0.0));
    let mut temperatures = [TEMPERATURE; MAX_LIQUID_THERMAL_CELLS];
    for (i, temperature) in temperatures.iter_mut().enumerate().take(CELLS) {
        *temperature = TEMPERATURE
            + increment * (1.0 / geometry.centres[i] - 1.0 / outer) / (1.0 / RADIUS - 1.0 / outer);
    }
    let stage = field
        .prepare(&interface, &temperatures, TEMPERATURE)
        .unwrap();
    let expected =
        4.0 * PI * field.reference.conductivity_w_m_k() * increment / (1.0 / RADIUS - 1.0 / outer);
    // The actual source conductivity varies with T; its small-signal solution
    // converges to the exact spherical Laplace profile.
    for &power in &stage.conductive_power[..=CELLS] {
        near(power, expected, 1.0e-5, 0.0);
    }
    near(stage.interface_heat_power_w(), -expected, 1.0e-5, 0.0);
    let observation = bubble(interface.state)
        .observe(pressure(), stage.interface_heat_power_w())
        .unwrap();
    let rates = field.finish(stage, &observation, pressure()).unwrap();
    assert!(rates.interface_heat_power_w < 0.0);
    assert!(rates.outer_conductive_power_w > 0.0);
    assert!(rates.energy_rate_residual_w.abs() < 1.0e-18);
}

#[test]
fn both_phase_directions_have_reciprocal_original_enthalpy_and_heat() {
    let field = field();
    for offset in [-0.05, 0.05] {
        let interface = thermal(state(TEMPERATURE, 0.0));
        let temperatures = [TEMPERATURE + offset; MAX_LIQUID_THERMAL_CELLS];
        let stage = field
            .prepare(&interface, &temperatures, TEMPERATURE + offset)
            .unwrap();
        let observation = bubble(interface.state)
            .observe(pressure(), stage.interface_heat_power_w())
            .unwrap();
        let rates = field.finish(stage, &observation, pressure()).unwrap();
        assert_eq!(observation.vapor_mass_rate_kg_s.signum(), offset.signum());

        // Reconstruct the unshifted total first law from original material
        // endpoints rather than reusing the operator's reported balance.
        let point = SaturatedWaterLiquidThermal::at_temperature_k(TEMPERATURE + offset).unwrap();
        let density = field.source.state().liquid_density_kg_m3();
        let liquid_power: f64 = field
            .grid
            .volumes_m3()
            .iter()
            .zip(rates.temperature_rate_k_s.iter())
            .map(|(volume, rate)| {
                density * volume * point.enthalpy_temperature_derivative_j_kg_k() * rate
            })
            .sum();
        let vapor = observation.vapor;
        let vapor_power = observation.vapor_mass_kg
            * vapor.vapor_internal_energy_temperature_derivative_j_kg_k()
            * observation.temperature_rate_k_s
            + vapor.vapor_internal_energy_j_kg() * observation.vapor_mass_rate_kg_s;
        let outer_enthalpy_power = -observation.vapor_mass_rate_kg_s * point.enthalpy_j_kg();
        let residual = vapor_power
            + liquid_power
            + vapor.pressure_pa() * observation.volume_rate_m3_s
            + rates.outer_conductive_power_w
            + outer_enthalpy_power
            - rates.viscous_heat_power_w;
        assert!(
            residual.abs() < 1.0e-14,
            "physical first-law residual={residual:e}"
        );
        assert!(rates.energy_rate_residual_w.abs() < 1.0e-16);
        assert!(rates.minimum_transport_time_s > 0.0);
    }
}

#[test]
fn original_radial_viscous_loss_is_returned_to_its_resolved_liquid_cells() {
    let field = field();
    for velocity in [-0.03, 0.03] {
        let interface = thermal(state(TEMPERATURE, velocity));
        let temperatures = [TEMPERATURE; MAX_LIQUID_THERMAL_CELLS];
        let stage = field
            .prepare(&interface, &temperatures, TEMPERATURE)
            .unwrap();
        let observation = bubble(interface.state)
            .observe(pressure(), stage.interface_heat_power_w())
            .unwrap();
        let rates = field.finish(stage, &observation, pressure()).unwrap();
        let expected = 16.0
            * PI
            * field.source.state().dynamic_viscosity_pa_s()
            * RADIUS
            * velocity
            * velocity;
        near(
            rates.viscous_heat_power_w + rates.exterior_viscous_heat_power_w,
            expected,
            2.0e-14,
            0.0,
        );
        assert!(rates.viscous_heat_power_w > 0.0);
        assert!(rates.exterior_viscous_heat_power_w > 0.0);
        assert!(
            rates.temperature_rate_k_s[..CELLS]
                .iter()
                .all(|rate| *rate > 0.0)
        );
        assert_eq!(rates.interface_heat_power_w, 0.0);
        assert_eq!(rates.outer_advected_power_w, 0.0);
        assert!(rates.energy_rate_residual_w.abs() < 1.0e-20);
    }
}

#[test]
fn thermal_approximation_and_unresolved_outflow_refuse_without_changing_the_field() {
    let mut thermal_limits = limits();
    thermal_limits.maximum_outflow_peclet = 1.0e-12;
    let field = SphericalLiquidField::new(source(), field().grid.clone(), thermal_limits).unwrap();
    let before = field.clone();
    for offset in [-0.05, 0.05] {
        let interface = thermal(state(TEMPERATURE, 0.0));
        let temperatures = [TEMPERATURE + offset; MAX_LIQUID_THERMAL_CELLS];
        let stage = field
            .prepare(&interface, &temperatures, TEMPERATURE + offset)
            .unwrap();
        let observation = bubble(interface.state)
            .observe(pressure(), stage.interface_heat_power_w())
            .unwrap();
        assert!(matches!(
            field.finish(stage, &observation, pressure()),
            Err(VaporBubbleError::OutsideModel {
                quantity: "outflow Dirichlet liquid half-cell Peclet number",
                ..
            })
        ));
        assert_eq!(field, before);
    }

    thermal_limits.maximum_outflow_peclet = 0.5;
    thermal_limits.maximum_pressure_heating_temperature_k = 1.0e-10;
    let pressure_limited =
        SphericalLiquidField::new(source(), field.grid.clone(), thermal_limits).unwrap();
    let interface = thermal(state(TEMPERATURE, 0.0));
    let temperatures = [TEMPERATURE; MAX_LIQUID_THERMAL_CELLS];
    let stage = pressure_limited
        .prepare(&interface, &temperatures, TEMPERATURE)
        .unwrap();
    let observation = bubble(interface.state)
        .observe(pressure(), stage.interface_heat_power_w())
        .unwrap();
    assert!(matches!(
        pressure_limited.finish(stage, &observation, pressure()),
        Err(VaporBubbleError::OutsideModel {
            quantity: "omitted liquid pressure-heating temperature scale",
            ..
        })
    ));
    let unsupported = [400.0; MAX_LIQUID_THERMAL_CELLS];
    assert!(matches!(
        field.prepare(&interface, &unsupported, TEMPERATURE),
        Err(VaporBubbleError::TemperatureOutsideDomain { .. })
    ));
    thermal_limits.minimum_temperature_k = 273.15;
    assert!(thermal_limits.validate().is_err());
}

#[test]
fn energy_changes_use_original_liquid_calorics_and_fixed_inventory() {
    let field = field();
    let observe = |temperature| {
        let interface = thermal(state(TEMPERATURE, 0.0));
        let temperatures = [temperature; MAX_LIQUID_THERMAL_CELLS];
        let stage = field
            .prepare(&interface, &temperatures, temperature)
            .unwrap();
        let observation = bubble(interface.state)
            .observe(pressure(), stage.interface_heat_power_w())
            .unwrap();
        field.finish(stage, &observation, pressure()).unwrap()
    };
    let initial = observe(TEMPERATURE);
    let final_state = observe(TEMPERATURE + 0.1);
    let final_material = SaturatedWaterLiquidThermal::at_temperature_k(TEMPERATURE + 0.1).unwrap();
    let expected = field.source.state().liquid_density_kg_m3()
        * field.grid.volumes_m3().iter().sum::<f64>()
        * (final_material.enthalpy_j_kg() - field.reference.enthalpy_j_kg());
    near(
        final_state.energy_change_from(&initial).unwrap(),
        expected,
        2.0e-15,
        0.0,
    );
    near(
        initial.energy_change_from(&final_state).unwrap(),
        -expected,
        2.0e-15,
        0.0,
    );
    assert_eq!(initial.energy_change_from(&initial).unwrap(), 0.0);
}
