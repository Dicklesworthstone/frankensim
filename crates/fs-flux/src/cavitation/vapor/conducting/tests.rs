use super::*;
use crate::cavitation::vapor::VaporBubbleLimits;
use crate::cavitation::vapor::liquid::{LiquidThermalLimits, SphericalLiquidGrid};
use crate::cavitation::{BubbleLimits, BubbleNumerics};
use core::f64::consts::PI;
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_material::cavitation::vapor::liquid::SaturatedWaterLiquidThermal;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

const RADIUS: f64 = 1e-4;
const TEMPERATURE: f64 = 373.15;
const CELLS: usize = 8;

fn source(surface_tension: f64) -> SaturatedWaterCavitation {
    SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(Temperature::new(TEMPERATURE)).unwrap(),
        surface_tension,
    )
    .unwrap()
}

fn phase_limits() -> VaporBubbleLimits {
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
        maximum_vapor_mach_number: 0.05,
        maximum_frozen_liquid_relative_change: 0.05,
        maximum_neglected_enthalpy_fraction: 0.001,
    }
}

fn field(source: SaturatedWaterCavitation) -> SphericalLiquidField {
    field_with_cells(source, CELLS)
}

fn field_with_cells(source: SaturatedWaterCavitation, cells: usize) -> SphericalLiquidField {
    SphericalLiquidField::new(
        source,
        SphericalLiquidGrid::equal_volume(RADIUS, 1.1 * RADIUS, cells).unwrap(),
        LiquidThermalLimits {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            maximum_source_relative_change: 0.05,
            maximum_pressure_path_capacity_fraction: 0.01,
            maximum_pressure_heating_temperature_k: 0.01,
            maximum_outflow_peclet: 0.9,
        },
    )
    .unwrap()
}

fn model(offset: f64, velocity: f64) -> ConductingVaporBubble {
    model_with_cells(offset, velocity, CELLS)
}

fn model_with_cells(offset: f64, velocity: f64, cells: usize) -> ConductingVaporBubble {
    let source = source(0.06);
    let phase = SaturatedVaporBubble::new(
        source,
        VaporBubbleState::try_new(0.0, RADIUS, velocity, TEMPERATURE).unwrap(),
        phase_limits(),
    )
    .unwrap();
    ConductingVaporBubble::new(
        phase,
        field_with_cells(source, cells),
        &vec![TEMPERATURE + offset; cells],
    )
    .unwrap()
}

fn equilibrium_pressure() -> f64 {
    let liquid = source(0.06).state();
    liquid.saturation_pressure_pa() - 2.0 * liquid.surface_tension_n_m() / RADIUS
}

fn numerics() -> ConductingVaporBubbleNumerics {
    ConductingVaporBubbleNumerics {
        bubble: VaporBubbleNumerics {
            radial: BubbleNumerics {
                relative_tolerance: 1e-10,
                radius_absolute_tolerance_m: 1e-14,
                velocity_absolute_tolerance_m_s: 1e-9,
                initial_step_s: 1e-9,
                minimum_step_s: 1e-14,
                maximum_step_s: 1e-7,
                maximum_attempts: 20_000,
            },
            temperature_absolute_tolerance_k: 1e-10,
            energy_absolute_tolerance_j: 1e-16,
            mass_absolute_tolerance_kg: 1e-22,
            energy_balance_tolerance_j: 1e-13,
            mass_balance_tolerance_kg: 1e-19,
        },
        liquid_temperature_absolute_tolerance_k: 1e-9,
        liquid_energy_balance_tolerance_j: 1e-13,
        total_energy_balance_tolerance_j: 2e-13,
    }
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn uniform_liquid_preserves_saturated_radial_equilibrium() {
    let mut bubble = model(0.0, 0.0);
    let step = bubble
        .advance_controlled(
            5e-6,
            &numerics(),
            &|_| equilibrium_pressure(),
            &|_| TEMPERATURE,
            &mut || false,
        )
        .unwrap();
    near(bubble.state().radius_m(), RADIUS, 0.0, 1e-17);
    near(bubble.state().temperature_k(), TEMPERATURE, 0.0, 1e-11);
    near(bubble.state().liquid_velocity_m_s(), 0.0, 0.0, 1e-12);
    for temperature in bubble.liquid_temperatures_k() {
        near(*temperature, TEMPERATURE, 0.0, 1e-11);
    }
    near(step.bubble.vapor_mass_change_kg, 0.0, 0.0, 1e-24);
    near(step.bubble.supplied_heat_j, 0.0, 0.0, 1e-18);
    near(step.liquid_energy_residual_j, 0.0, 0.0, 1e-17);
    near(step.total_energy_residual_j, 0.0, 0.0, 1e-17);
    near(step.viscous_heat_j, 0.0, 0.0, 1e-24);
    assert_eq!(bubble.state().time_s(), 5e-6);
}

#[test]
fn resolved_temperature_reciprocally_drives_phase_change_and_retains_its_energy() {
    for offset in [-0.05, 0.05] {
        let mut bubble = model(offset, 0.0);
        let bath = TEMPERATURE + offset;
        let pressure = equilibrium_pressure();
        let initial = bubble.observe(pressure, bath).unwrap();
        assert_eq!(
            initial.bubble.heat_power_w,
            initial.liquid.interface_heat_power_w
        );
        assert_eq!(initial.bubble.heat_power_w.is_sign_positive(), offset > 0.0);
        assert_eq!(
            initial.bubble.vapor_mass_rate_kg_s.is_sign_positive(),
            offset > 0.0
        );
        let anchor = bubble.liquid().reference_enthalpy_j_kg();
        let rho = bubble.model().source().state().liquid_density_kg_m3();
        let initial_enthalpy = SaturatedWaterLiquidThermal::at_temperature_k(bath)
            .unwrap()
            .enthalpy_j_kg();
        let mut external_work_and_heat = 0.0;
        let mut transferred_mass = 0.0;
        for end in [1e-6, 2e-6] {
            let step = bubble
                .advance_controlled(end, &numerics(), &|_| pressure, &|_| bath, &mut || false)
                .unwrap();
            external_work_and_heat += step.bubble.pressure_work_j
                + step.outer_conductive_heat_j
                + step.outer_advected_enthalpy_j
                - step.viscous_heat_j;
            transferred_mass += step.bubble.transferred_vapor_mass_kg;
            assert!(
                step.liquid_energy_residual_j.abs() <= numerics().liquid_energy_balance_tolerance_j
            );
            assert!(
                step.total_energy_residual_j.abs() <= numerics().total_energy_balance_tolerance_j
            );
        }
        let endpoint = bubble.observe(pressure, bath).unwrap();
        let liquid_energy_change: f64 = bubble
            .liquid_temperatures_k()
            .iter()
            .zip(bubble.liquid().grid().volumes_m3())
            .map(|(temperature, volume)| {
                let enthalpy = SaturatedWaterLiquidThermal::at_temperature_k(*temperature)
                    .unwrap()
                    .enthalpy_j_kg();
                rho * volume * (enthalpy - initial_enthalpy)
            })
            .sum();
        let phase_energy_change = endpoint.bubble.vapor_mass_kg
            * (endpoint.bubble.vapor.vapor_internal_energy_j_kg() - anchor)
            - initial.bubble.vapor_mass_kg
                * (initial.bubble.vapor.vapor_internal_energy_j_kg() - anchor);
        near(
            phase_energy_change + liquid_energy_change + external_work_and_heat,
            0.0,
            0.0,
            2.0 * numerics().total_energy_balance_tolerance_j,
        );
        near(
            endpoint.bubble.vapor_mass_kg - initial.bubble.vapor_mass_kg,
            transferred_mass,
            0.0,
            2.0 * numerics().bubble.mass_balance_tolerance_kg,
        );
        assert_eq!(
            (endpoint.bubble.vapor_mass_kg - initial.bubble.vapor_mass_kg).is_sign_positive(),
            offset > 0.0
        );
        assert_eq!(
            (bubble.liquid_temperatures_k()[0] - bath).is_sign_negative(),
            offset > 0.0
        );
        assert!(bubble.liquid_temperatures_k()[0] != bubble.liquid_temperatures_k()[CELLS - 1]);
        assert_eq!(bubble.state().time_s(), 2e-6);
    }
}

#[test]
fn native_viscous_dissipation_returns_only_its_resolved_share_to_the_liquid() {
    let velocity = 0.02;
    let mut bubble = model(0.0, velocity);
    let pressure = equilibrium_pressure();
    let observation = bubble.observe(pressure, TEMPERATURE).unwrap();
    let viscosity = bubble.model().source().state().dynamic_viscosity_pa_s();
    let all_liquid = 16.0 * PI * viscosity * RADIUS * velocity * velocity;
    near(
        observation.liquid.viscous_heat_power_w + observation.liquid.exterior_viscous_heat_power_w,
        all_liquid,
        3e-15,
        1e-25,
    );
    assert!(observation.liquid.viscous_heat_power_w > 0.0);
    assert!(observation.liquid.exterior_viscous_heat_power_w > 0.0);
    let step = bubble
        .advance_controlled(
            2e-7,
            &numerics(),
            &|_| pressure,
            &|_| TEMPERATURE,
            &mut || false,
        )
        .unwrap();
    assert!(step.viscous_heat_j > 0.0);
    assert!(step.total_energy_residual_j.abs() <= numerics().total_energy_balance_tolerance_j);
    assert!(step.liquid_energy_residual_j.abs() <= numerics().liquid_energy_balance_tolerance_j);
}

#[test]
fn dropped_proposals_late_cancellation_and_work_refusal_preserve_both_phases() {
    let mut bubble = model(0.05, 0.0);
    let original_state = bubble.state();
    let original_profile = bubble.liquid_temperatures_k().to_vec();
    let pressure = equilibrium_pressure();
    let bath = TEMPERATURE + 0.05;
    let mut complete_poll_count = 0;
    {
        let proposal = bubble
            .prepare_advance_controlled(1e-6, &numerics(), &|_| pressure, &|_| bath, &mut || {
                complete_poll_count += 1;
                false
            })
            .unwrap();
        assert_ne!(proposal.step().bubble.observation.state, original_state);
        assert_ne!(
            proposal.liquid_temperatures_k(),
            original_profile.as_slice()
        );
    }
    assert_eq!(bubble.state(), original_state);
    assert_eq!(bubble.liquid_temperatures_k(), original_profile.as_slice());
    let mut polls = 0;
    assert!(matches!(
        bubble.advance_controlled(1e-6, &numerics(), &|_| pressure, &|_| bath, &mut || {
            polls += 1;
            polls >= complete_poll_count
        }),
        Err(ConductingVaporBubbleError::Cancelled)
    ));
    assert_eq!(polls, complete_poll_count);
    assert_eq!(bubble.state(), original_state);
    assert_eq!(bubble.liquid_temperatures_k(), original_profile.as_slice());
    let mut bounded = numerics();
    bounded.bubble.radial.maximum_attempts = 1;
    assert!(
        bubble
            .advance_controlled(1e-6, &bounded, &|_| pressure, &|_| bath, &mut || false)
            .is_err()
    );
    assert_eq!(bubble.state(), original_state);
    assert_eq!(bubble.liquid_temperatures_k(), original_profile.as_slice());
    assert!(
        bubble
            .advance_controlled(1e-6, &numerics(), &|_| pressure, &|_| f64::NAN, &mut || {
                false
            })
            .is_err()
    );
    assert_eq!(bubble.state(), original_state);
    assert_eq!(bubble.liquid_temperatures_k(), original_profile.as_slice());
}

#[test]
fn independent_source_data_and_incomplete_profiles_are_refused() {
    let reference = model(0.0, 0.0);
    assert!(matches!(
        ConductingVaporBubble::new(
            *reference.model(),
            field(source(0.061)),
            &[TEMPERATURE; CELLS],
        ),
        Err(ConductingVaporBubbleError::SourceMismatch)
    ));
    assert!(matches!(
        ConductingVaporBubble::new(
            *reference.model(),
            reference.liquid().clone(),
            &[TEMPERATURE; CELLS - 1],
        ),
        Err(ConductingVaporBubbleError::InvalidInput(_))
    ));
}

#[test]
fn coupled_temperature_and_phase_exchange_converge_under_space_and_time_refinement() {
    let pressure = equilibrium_pressure();
    let bath = TEMPERATURE + 0.05;
    let solve = |cells: usize, controls: &ConductingVaporBubbleNumerics| {
        let mut bubble = model_with_cells(0.05, 0.0, cells);
        let step = bubble
            .advance_controlled(2e-6, controls, &|_| pressure, &|_| bath, &mut || false)
            .unwrap();
        (
            bubble.state().temperature_k(),
            step.bubble.supplied_heat_j,
            step.bubble.vapor_mass_change_kg,
        )
    };

    // This exercises the moving field and its reciprocal nonlinear phase
    // response, including the initially unresolved liquid/interface jump.
    // The test checks convergence of physical observables, not just a small
    // conserved residual at one spatial resolution.
    let spatial = [16, 32, 64].map(|cells| solve(cells, &numerics()));
    for (coarse, fine, floor) in [
        (
            (spatial[0].0 - spatial[1].0).abs(),
            (spatial[1].0 - spatial[2].0).abs(),
            1e-8,
        ),
        (
            (spatial[0].1 - spatial[1].1).abs(),
            (spatial[1].1 - spatial[2].1).abs(),
            1e-16,
        ),
        (
            (spatial[0].2 - spatial[1].2).abs(),
            (spatial[1].2 - spatial[2].2).abs(),
            1e-22,
        ),
    ] {
        assert!(
            fine < 0.8 * coarse + floor,
            "spatial refinement: {fine:e} >= {coarse:e}"
        );
    }

    // Fix a single spatial operator and reduce the authored maximum step.
    // Loose local scales keep the step ceiling in control; the independent
    // first-law budgets still check every complete physical window.
    let temporal = [8e-8, 4e-8, 2e-8].map(|step_s| {
        let mut controls = numerics();
        controls.bubble.radial.relative_tolerance = 1e-4;
        controls.bubble.radial.initial_step_s = step_s;
        controls.bubble.radial.maximum_step_s = step_s;
        controls.bubble.temperature_absolute_tolerance_k = 1e-6;
        controls.liquid_temperature_absolute_tolerance_k = 1e-6;
        controls.bubble.energy_absolute_tolerance_j = 1e-12;
        controls.bubble.mass_absolute_tolerance_kg = 1e-17;
        controls.bubble.energy_balance_tolerance_j = 1e-10;
        controls.bubble.mass_balance_tolerance_kg = 1e-17;
        controls.liquid_energy_balance_tolerance_j = 1e-10;
        controls.total_energy_balance_tolerance_j = 2e-10;
        solve(32, &controls)
    });
    for (coarse, fine, floor) in [
        (
            (temporal[0].0 - temporal[1].0).abs(),
            (temporal[1].0 - temporal[2].0).abs(),
            1e-10,
        ),
        (
            (temporal[0].1 - temporal[1].1).abs(),
            (temporal[1].1 - temporal[2].1).abs(),
            1e-18,
        ),
        (
            (temporal[0].2 - temporal[1].2).abs(),
            (temporal[1].2 - temporal[2].2).abs(),
            1e-24,
        ),
    ] {
        assert!(
            fine < 0.8 * coarse + floor,
            "temporal refinement: {fine:e} >= {coarse:e}"
        );
    }
}
