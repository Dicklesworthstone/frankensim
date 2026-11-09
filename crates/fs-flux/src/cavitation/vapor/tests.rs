use super::*;
use crate::cavitation::{BubbleGas, BubbleState, RayleighPlessetBubble};
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

const RADIUS: f64 = 1e-4;

fn source(temperature: f64) -> SaturatedWaterCavitation {
    SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(Temperature::new(temperature)).unwrap(),
        0.072,
    )
    .unwrap()
}

fn limits() -> VaporBubbleLimits {
    VaporBubbleLimits {
        bubble: BubbleLimits {
            minimum_radius_m: 0.5 * RADIUS,
            maximum_radius_m: 2.0 * RADIUS,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
        minimum_temperature_k: 273.15,
        maximum_temperature_k: 500.0,
        maximum_pressure_ratio: 0.001,
        maximum_vapor_mach_number: 0.05,
        maximum_frozen_liquid_relative_change: 0.05,
        maximum_neglected_enthalpy_fraction: 0.001,
    }
}

fn numerics() -> VaporBubbleNumerics {
    VaporBubbleNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-15,
            velocity_absolute_tolerance_m_s: 1e-10,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-13,
            maximum_step_s: 5e-7,
            maximum_attempts: 20_000,
        },
        temperature_absolute_tolerance_k: 1e-10,
        energy_absolute_tolerance_j: 1e-16,
        mass_absolute_tolerance_kg: 1e-22,
        energy_balance_tolerance_j: 2e-13,
        mass_balance_tolerance_kg: 1e-20,
    }
}

fn model(temperature: f64, liquid_velocity: f64) -> SaturatedVaporBubble {
    SaturatedVaporBubble::new(
        source(temperature),
        VaporBubbleState::try_new(0.0, RADIUS, liquid_velocity, temperature).unwrap(),
        limits(),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn no_heat(_: f64, _: &VaporBubbleThermalState) -> Result<f64, VaporBubbleError> {
    Ok(0.0)
}

#[test]
fn stationary_saturated_nucleus_preserves_native_radial_equilibrium() {
    let mut bubble = model(293.15, 0.0);
    let liquid = bubble.source().state();
    let pressure = liquid.saturation_pressure_pa() - 2.0 * liquid.surface_tension_n_m() / RADIUS;
    let original = RayleighPlessetBubble::try_new(
        liquid,
        BubbleGas::try_new(RADIUS, 0.0, 1.0).unwrap(),
        limits().bubble,
    )
    .unwrap();
    let native = original
        .observe(&BubbleState::try_new(0.0, RADIUS, 0.0).unwrap(), pressure)
        .unwrap();
    let observation = bubble.observe(pressure, 0.0).unwrap();
    assert_eq!(
        observation.liquid_acceleration_m_s2,
        native.wall_acceleration_m_s2
    );
    assert_eq!(observation.wall_velocity_m_s, 0.0);
    assert_eq!(observation.temperature_rate_k_s, 0.0);
    assert_eq!(observation.vapor_mass_rate_kg_s, 0.0);
    let step = bubble
        .advance_controlled(1e-5, &numerics(), &|_| pressure, &no_heat, &mut || false)
        .unwrap();
    assert_eq!(bubble.state().radius_m(), RADIUS);
    assert_eq!(bubble.state().temperature_k(), 293.15);
    assert_eq!(step.mass_residual_kg, 0.0);
    assert_eq!(step.energy_residual_j, 0.0);
    assert_eq!(step.supplied_heat_j, 0.0);
    assert_eq!(step.transferred_vapor_mass_kg, 0.0);
}

#[test]
fn phase_mass_energy_and_generalized_rp_close_at_each_actual_interface_state() {
    for temperature in [293.15, 450.0] {
        for velocity in [-0.03, 0.0, 0.03] {
            let bubble = model(temperature, velocity);
            let liquid = bubble.source().state();
            for heat in [-1e-4, 0.0, 1e-4] {
                let pressure = liquid.saturation_pressure_pa() - 1000.0;
                let o = bubble.observe(pressure, heat).unwrap();
                let rho = liquid.liquid_density_kg_m3();
                let rv = o.vapor.vapor_density_kg_m3();
                let j = o.vapor_mass_flux_kg_m2_s;
                let area = 4.0 * PI * RADIUS * RADIUS;
                near(o.wall_velocity_m_s, velocity + j / rho, 2e-15, 1e-20);
                let mass_rate = rv * o.volume_rate_m3_s
                    + o.volume_m3
                        * o.vapor.vapor_density_temperature_derivative_kg_m3_k()
                        * o.temperature_rate_k_s;
                near(o.vapor_mass_rate_kg_s, mass_rate, 2e-13, 1e-26);
                let energy_rate = o.vapor.vapor_internal_energy_j_kg() * mass_rate
                    + o.vapor_mass_kg
                        * o.vapor
                            .vapor_internal_energy_temperature_derivative_j_kg_k()
                        * o.temperature_rate_k_s;
                let first_law = heat + o.vapor.liquid_enthalpy_j_kg() * mass_rate
                    - o.vapor.pressure_pa() * o.volume_rate_m3_s;
                near(energy_rate, first_law, 3e-12, 1e-17);
                // Published generalized RP, written in J=-j convention,
                // independently of the native shared acceleration helper.
                let outward_flux = -j;
                let left = RADIUS * o.liquid_acceleration_m_s2 + 1.5 * velocity * velocity
                    - outward_flux / rho * (2.0 * velocity + outward_flux * (1.0 / rv - 1.0 / rho));
                let right = (o.vapor.pressure_pa()
                    - pressure
                    - 2.0 * liquid.surface_tension_n_m() / RADIUS
                    - 4.0 * liquid.dynamic_viscosity_pa_s() * velocity / RADIUS)
                    / rho;
                near(left, right, 2e-13, 1e-15);
                near(o.liquid_volume_flow_m3_s, area * velocity, 2e-15, 0.0);
                near(
                    o.volume_rate_m3_s - o.liquid_volume_flow_m3_s,
                    o.vapor_mass_rate_kg_s / rho,
                    2e-9,
                    1e-24,
                );
            }
        }
    }
}

#[test]
fn manufactured_isothermal_evaporation_and_condensation_have_exact_latent_heat_and_liquid_flux() {
    for temperature in [293.15, 450.0] {
        for wall_velocity in [-0.1, 0.1] {
            let source = source(temperature);
            let liquid = source.state();
            let vapor = SaturatedWaterVapor::try_from_saturation(*source.saturation()).unwrap();
            let rho = liquid.liquid_density_kg_m3();
            let rv = vapor.vapor_density_kg_m3();
            let ratio = rv / rho;
            let velocity = (1.0 - ratio) * wall_velocity;
            let mut bubble = model(temperature, velocity);
            let end = 2e-5;
            let pressure = |time: f64| {
                let radius = RADIUS + wall_velocity * time;
                let flux = rv * wall_velocity;
                vapor.pressure_pa()
                    - 2.0 * liquid.surface_tension_n_m() / radius
                    - 4.0 * liquid.dynamic_viscosity_pa_s() * velocity / radius
                    + flux * flux * (1.0 / rv - 1.0 / rho)
                    - rho * (2.0 * wall_velocity * velocity - 0.5 * velocity * velocity)
            };
            let heat = |_: f64, s: &VaporBubbleThermalState| {
                Ok(4.0
                    * PI
                    * s.state.radius_m().powi(2)
                    * s.vapor.vapor_density_kg_m3()
                    * s.vapor.latent_enthalpy_j_kg()
                    * s.state.liquid_velocity_m_s()
                    / (1.0 - s.vapor.vapor_density_kg_m3() / rho))
            };
            let step = bubble
                .advance_controlled(end, &numerics(), &pressure, &heat, &mut || false)
                .unwrap();
            let expected_radius = RADIUS + wall_velocity * end;
            near(bubble.state().radius_m(), expected_radius, 2e-12, 1e-16);
            near(bubble.state().liquid_velocity_m_s(), velocity, 2e-11, 1e-12);
            near(bubble.state().temperature_k(), temperature, 0.0, 2e-10);
            let delta_volume = (4.0 * PI / 3.0)
                * (expected_radius - RADIUS)
                * (expected_radius.powi(2) + expected_radius * RADIUS + RADIUS.powi(2));
            let delta_mass = rv * delta_volume;
            near(step.vapor_mass_change_kg, delta_mass, 2e-10, 1e-24);
            near(step.transferred_vapor_mass_kg, delta_mass, 2e-10, 1e-24);
            near(
                step.supplied_heat_j,
                vapor.latent_enthalpy_j_kg() * delta_mass,
                2e-10,
                1e-18,
            );
            near(
                step.pressure_work_j,
                vapor.pressure_pa() * delta_volume,
                2e-10,
                1e-18,
            );
            near(
                step.vapor_energy_change_j,
                vapor.vapor_internal_energy_j_kg() * delta_mass,
                2e-10,
                1e-18,
            );
            near(
                step.transferred_enthalpy_j,
                vapor.liquid_enthalpy_j_kg() * delta_mass,
                2e-10,
                1e-18,
            );
            let o = step.observation;
            near(
                o.liquid_volume_flow_m3_s,
                (1.0 - ratio) * o.volume_rate_m3_s,
                2e-13,
                1e-25,
            );
            near(
                o.liquid_volume_flow_rate_m3_s2,
                (1.0 - ratio) * 8.0 * PI * expected_radius * wall_velocity.powi(2),
                2e-8,
                1e-14,
            );
            assert!((o.liquid_volume_flow_m3_s - o.volume_rate_m3_s).abs() > 0.0);
            assert!(step.energy_residual_j.abs() <= numerics().energy_balance_tolerance_j);
            assert!(step.mass_residual_kg.abs() <= numerics().mass_balance_tolerance_kg);
            assert!(step.neglected_interface_energy_scale_j > 0.0);
        }
    }
}

#[test]
fn signed_heat_drives_real_nonlinear_phase_mass_and_saturation_pressure_feedback() {
    for heat in [-1e-4, 1e-4] {
        let mut bubble = model(293.15, 0.0);
        let original = bubble;
        let liquid = bubble.source().state();
        let pressure =
            liquid.saturation_pressure_pa() - 2.0 * liquid.surface_tension_n_m() / RADIUS;
        let end = 5e-5;
        let step = bubble
            .advance_controlled(
                end,
                &numerics(),
                &|_| pressure,
                &|_, _| Ok(heat),
                &mut || false,
            )
            .unwrap();
        assert_eq!(
            step.vapor_mass_change_kg.is_sign_positive(),
            heat.is_sign_positive()
        );
        assert_eq!(
            (bubble.state().temperature_k() - 293.15).is_sign_positive(),
            heat.is_sign_positive()
        );
        assert_eq!(
            (step.observation.vapor.pressure_pa() - liquid.saturation_pressure_pa())
                .is_sign_positive(),
            heat.is_sign_positive()
        );
        assert_eq!(
            (bubble.state().radius_m() - RADIUS).is_sign_positive(),
            heat.is_sign_positive()
        );
        assert_ne!(bubble.state(), original.state());
        assert_eq!(bubble.source(), original.source());
        near(step.supplied_heat_j, heat * end, 2e-14, 1e-24);
        near(
            step.reference_invariant_energy_change_j + step.pressure_work_j,
            step.supplied_heat_j + step.reference_invariant_transferred_enthalpy_j,
            0.0,
            numerics().energy_balance_tolerance_j,
        );
        near(
            step.original_datum_energy_residual_j,
            step.energy_residual_j + step.reference_liquid_enthalpy_j_kg * step.mass_residual_kg,
            2e-15,
            1e-30,
        );
    }
}

#[test]
fn physical_refusals_late_cancellation_and_dropped_proposals_preserve_the_complete_state() {
    let mut bubble = model(293.15, 0.0);
    let original = bubble;
    let liquid = bubble.source().state();
    let pressure = liquid.saturation_pressure_pa() - 2.0 * liquid.surface_tension_n_m() / RADIUS;
    {
        let proposal = bubble
            .prepare_advance_controlled(
                1e-5,
                &numerics(),
                &|_| pressure,
                &|_, _| Ok(1e-4),
                &mut || false,
            )
            .unwrap();
        assert!(proposal.step().observation.state.temperature_k() > 293.15);
    }
    assert_eq!(bubble, original);
    let mut polls = 0;
    assert!(matches!(
        bubble.advance_controlled(
            1e-5,
            &numerics(),
            &|_| pressure,
            &|_, _| Ok(1e-4),
            &mut || {
                polls += 1;
                polls >= 20
            }
        ),
        Err(VaporBubbleError::Cancelled)
    ));
    assert_eq!(bubble, original);
    assert!(
        bubble
            .advance_controlled(
                1e-5,
                &numerics(),
                &|_| pressure,
                &|time, _| if time > 2e-6 {
                    Err(VaporBubbleError::InvalidInput("heat port refused"))
                } else {
                    Ok(1e-4)
                },
                &mut || false
            )
            .is_err()
    );
    assert_eq!(bubble, original);
    let restricted = VaporBubbleLimits {
        maximum_temperature_k: 293.15001,
        ..limits()
    };
    let mut bounded =
        SaturatedVaporBubble::new(*bubble.source(), bubble.state(), restricted).unwrap();
    let before = bounded;
    assert!(
        bounded
            .advance_controlled(
                1e-5,
                &numerics(),
                &|_| pressure,
                &|_, _| Ok(1e-4),
                &mut || false
            )
            .is_err()
    );
    assert_eq!(bounded, before);
    assert!(matches!(
        bubble.observe(pressure, 1000.0),
        Err(VaporBubbleError::OutsideModel { .. })
    ));
    assert!(matches!(
        bubble.observe(1e12, 0.0),
        Err(VaporBubbleError::OutsideModel { .. })
    ));
    assert!(bubble.observe(pressure, f64::NAN).is_err());
    let strict = VaporBubbleLimits {
        maximum_neglected_enthalpy_fraction: 1e-12,
        ..limits()
    };
    let bounded = SaturatedVaporBubble::new(*bubble.source(), bubble.state(), strict).unwrap();
    assert!(matches!(
        bounded.observe(pressure, 0.0),
        Err(VaporBubbleError::OutsideModel { .. })
    ));
    let strict = VaporBubbleLimits {
        maximum_frozen_liquid_relative_change: 1e-5,
        ..limits()
    };
    assert!(matches!(
        SaturatedVaporBubble::new(
            *bubble.source(),
            VaporBubbleState::try_new(0.0, RADIUS, 0.0, 294.15).unwrap(),
            strict
        ),
        Err(VaporBubbleError::OutsideModel {
            quantity: "frozen liquid coefficient variation",
            ..
        })
    ));
    assert_eq!(bubble, original);
}
