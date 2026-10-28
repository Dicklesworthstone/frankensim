use super::*;
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

fn law(temperature: f64) -> IsobaricWaterVaporization {
    IsobaricWaterVaporization::from_saturation(
        SaturatedWaterState::at_temperature(Temperature::new(temperature)).unwrap(),
    )
    .unwrap()
}

fn enthalpy(law: IsobaricWaterVaporization, quality: f64) -> f64 {
    let source = law.saturation();
    source.liquid_specific_enthalpy().value()
        + quality
            * (source.vapor_specific_enthalpy().value() - source.liquid_specific_enthalpy().value())
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (budget {tolerance})"
    );
}

#[test]
fn coexistence_uses_real_vapor_inventory_and_additive_specific_volume() {
    let law = law(373.15);
    let source = law.saturation();
    for quality in [0.0, 0.1, 0.5, 0.9, 1.0] {
        let state = law
            .state_at_specific_enthalpy(enthalpy(law, quality))
            .unwrap();
        near(state.vapor_mass_fraction(), quality, 2e-15);
        near(state.liquid_mass_fraction(), 1.0 - quality, 2e-15);
        assert_eq!(state.temperature_k(), 373.15);
        assert_eq!(state.pressure_pa(), source.pressure().value());
        let volume = (1.0 - quality) / source.liquid_density().value()
            + quality / source.vapor_density().value();
        near(state.specific_volume_m3_kg(), volume, 2e-15);
        near(state.bulk_density_kg_m3(), 1.0 / volume, 2e-11);
        near(
            state.specific_internal_energy_j_kg(),
            (1.0 - quality) * source.liquid_specific_internal_energy().value()
                + quality * source.vapor_specific_internal_energy().value(),
            1e-9,
        );
        near(
            state.specific_entropy_j_kg_k(),
            (1.0 - quality) * source.liquid_specific_entropy().value()
                + quality * source.vapor_specific_entropy().value(),
            4e-12,
        );
        near(
            state.specific_enthalpy_j_kg() - state.specific_internal_energy_j_kg(),
            state.pressure_pa() * state.specific_volume_m3_kg(),
            1e-9,
        );
    }
    let liquid = law
        .state_at_specific_enthalpy(source.liquid_specific_enthalpy().value())
        .unwrap();
    let vapor = law
        .state_at_specific_enthalpy(source.vapor_specific_enthalpy().value())
        .unwrap();
    assert_eq!(liquid.vapor_mass_fraction(), 0.0);
    assert_eq!(vapor.vapor_mass_fraction(), 1.0);
    assert_eq!(
        liquid.specific_entropy_j_kg_k(),
        source.liquid_specific_entropy().value()
    );
    assert_eq!(
        vapor.specific_entropy_j_kg_k(),
        source.vapor_specific_entropy().value()
    );
    assert_eq!(
        vapor.specific_internal_energy_j_kg(),
        source.vapor_specific_internal_energy().value()
    );
    assert_eq!(ThermalLaw::maximum_temperature_slope(&law), 0.0);
    near(
        ThermalLaw::specific_entropy_change_j_kg_k(
            &law,
            liquid.specific_enthalpy_j_kg(),
            vapor.specific_enthalpy_j_kg(),
        )
        .unwrap()
        .unwrap(),
        source.latent_entropy().value(),
        4e-12,
    );
}

#[test]
fn heating_and_condensation_pay_real_expansion_work() {
    let law = law(373.15);
    let source = law.saturation();
    let mass = 0.5;
    let initial = law
        .state_at_specific_enthalpy(enthalpy(law, 0.125))
        .unwrap();
    let heat = mass * (enthalpy(law, 0.75) - initial.specific_enthalpy_j_kg());
    let step = law
        .advance_heat(initial, mass, heat, 1e-8, || false)
        .unwrap();
    let vaporized = mass * (0.75 - 0.125);
    near(step.state.vapor_mass_fraction(), 0.75, 2e-15);
    near(step.initial_phase_masses_kg[0], mass * 0.875, 1e-15);
    near(step.phase_masses_kg[1], mass * 0.75, 1e-15);
    near(step.phase_masses_kg.iter().sum(), mass, 1e-15);
    let volume_change =
        vaporized * (1.0 / source.vapor_density().value() - 1.0 / source.liquid_density().value());
    near(step.volume_change_m3, volume_change, 2e-15);
    near(
        step.boundary_work_j,
        source.pressure().value() * volume_change,
        1e-8,
    );
    near(
        step.internal_energy_change_j,
        vaporized
            * (source.vapor_specific_internal_energy().value()
                - source.liquid_specific_internal_energy().value()),
        1e-8,
    );
    near(
        step.entropy_change_j_k,
        heat / initial.temperature_k(),
        1e-10,
    );
    near(step.enthalpy_change_j, heat, 1e-8);
    assert!(step.boundary_work_j > 0.05 * heat);
    assert!(step.internal_energy_change_j < 0.95 * heat);
    near(
        step.internal_energy_change_j + step.boundary_work_j,
        heat,
        1e-8,
    );
    let cooled = law
        .advance_heat(step.state, mass, -heat, 1e-8, || false)
        .unwrap();
    near(
        cooled.state.specific_enthalpy_j_kg(),
        initial.specific_enthalpy_j_kg(),
        1e-9,
    );
    near(
        cooled.state.vapor_mass_fraction(),
        initial.vapor_mass_fraction(),
        2e-15,
    );
    near(cooled.boundary_work_j, -step.boundary_work_j, 1e-8);
    near(cooled.entropy_change_j_k, -step.entropy_change_j_k, 1e-10);
    assert!(cooled.internal_energy_change_j < 0.0);
}

#[test]
fn complete_boiling_condensation_and_phase_reentry_keep_exact_source_endpoints() {
    for temperature in [273.16, 300.0, 373.15, 450.0, 600.0, 647.095] {
        let law = law(temperature);
        let source = law.saturation();
        let liquid_h = source.liquid_specific_enthalpy().value();
        let vapor_h = source.vapor_specific_enthalpy().value();
        let liquid = law.state_at_specific_enthalpy(liquid_h).unwrap();
        for mass in [0.5, 0.7, 1.0] {
            let heat = mass * (vapor_h - liquid_h);
            let boiled = law
                .advance_heat(liquid, mass, heat, 1e-7, || false)
                .unwrap();
            assert_eq!(boiled.state.specific_enthalpy_j_kg(), vapor_h);
            assert_eq!(boiled.phase_masses_kg, [0.0, mass]);
            let condensed = law
                .advance_heat(boiled.state, mass, -heat, 1e-7, || false)
                .unwrap();
            assert_eq!(condensed.state, liquid);
            assert_eq!(condensed.phase_masses_kg, [mass, 0.0]);
            assert!(
                law.advance_heat(condensed.state, mass, 1.0, 1e-7, || false)
                    .unwrap()
                    .state
                    .vapor_mass_fraction()
                    > 0.0
            );
            assert_eq!(
                law.advance_heat(liquid, mass, heat + 1.0, heat, || false),
                Err(WaterVaporizationError::OutsideCoexistence)
            );
            assert_eq!(
                law.advance_heat(boiled.state, mass, -heat - 1.0, heat, || false),
                Err(WaterVaporizationError::OutsideCoexistence)
            );
        }
    }
}

#[test]
fn pressure_binding_domain_accuracy_and_cancellation_are_explicit() {
    let controls = WaterSaturationSolve {
        pressure_tolerance_pa: 1e-6,
        maximum_iterations: 128,
    };
    let law = IsobaricWaterVaporization::at_pressure(Pressure::new(101_325.0), controls, || false)
        .unwrap();
    near(law.saturation().pressure().value(), 101_325.0, 1e-6);
    assert!(law.saturation().temperature().value() > 373.0);
    assert!(law.saturation().temperature().value() < 374.0);
    let source = law.saturation();
    assert_eq!(
        law.state_at_specific_enthalpy(source.liquid_specific_enthalpy().value() - 1.0),
        Err(WaterVaporizationError::OutsideCoexistence)
    );
    assert_eq!(
        law.state_at_specific_enthalpy(source.vapor_specific_enthalpy().value() + 1.0),
        Err(WaterVaporizationError::OutsideCoexistence)
    );
    assert!(law.state_at_specific_enthalpy(f64::NAN).is_err());
    assert_eq!(
        IsobaricWaterVaporization::from_saturation(
            SaturatedWaterState::at_temperature(Temperature::new(647.096)).unwrap()
        ),
        Err(WaterVaporizationError::CollapsedCoexistence)
    );
    let current = law.state_at_specific_enthalpy(enthalpy(law, 0.5)).unwrap();
    assert!(matches!(
        law.advance_heat(current, 0.7, 0.1, 0.0, || false),
        Err(WaterVaporizationError::EnergyBalance { .. })
    ));
    let represented = law.advance_heat(current, 0.7, 0.1, 1e-8, || false).unwrap();
    assert!(represented.heat_residual_j != 0.0);
    assert!(matches!(
        law.advance_heat(current, 0.7, 1e-20, 1.0, || false),
        Err(WaterVaporizationError::Unrepresentable(_))
    ));
    assert_eq!(
        super::tests::law(400.0).advance_heat(current, 0.7, 1.0, 1e-8, || false),
        Err(WaterVaporizationError::StateMismatch)
    );
    let mut polls = 0;
    assert_eq!(
        law.advance_heat(current, 0.7, 100.0, 1e-8, || {
            polls += 1;
            polls == 2
        }),
        Err(WaterVaporizationError::Cancelled)
    );
    assert_eq!(
        law.advance_heat(current, 0.7, 0.0, 0.0, || false)
            .unwrap()
            .state,
        current
    );
}

fn boundary(ambient: f64, htc: f64, source: f64) -> Boundary {
    Boundary {
        ambient_k: ambient,
        htc_w_m2_k: htc,
        reservoir: None,
        source_w: source,
    }
}

fn birth(
    id: u64,
    law: IsobaricWaterVaporization,
    quality: f64,
    mass: f64,
    boundary: Boundary,
) -> TimedEvent {
    let h = enthalpy(law, quality);
    let state = law.state_at_specific_enthalpy(h).unwrap();
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: id,
            material_id: MaterialId(id),
            mass_kg: mass,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64, 0.0, 1.0],
                volume_m3: mass / state.bulk_density_kg_m3(),
                surface_area_m2: 1.0,
            },
            boundary,
        }),
    }
}

fn budget() -> DepositionBudget {
    DepositionBudget {
        max_cells: 2,
        max_links: 1,
        max_events: 4,
        max_substeps: 1000,
    }
}

fn step_config() -> StepConfig {
    StepConfig {
        max_dt_s: 2.0,
        energy_tolerance_j: 1e-7,
    }
}

#[test]
fn real_thermal_history_boils_and_condenses_with_ambient_loss() {
    let law = law(373.15);
    let mass = 0.5;
    let initial_h = enthalpy(law, 0.25);
    let loss_w = 2.0 * (373.15 - 300.0);
    let initial = DepositionSimulation::new(
        vec![
            birth(1, law, 0.25, mass, boundary(300.0, 2.0, 1000.0)),
            TimedEvent {
                time_s: 20.0,
                event: Event::Boundary(CellId(1), boundary(300.0, 2.0, -1000.0 + 2.0 * loss_w)),
            },
        ],
        &[(MaterialId(1), &law)],
        budget(),
    )
    .unwrap();
    let (active, _) = initial.advance_to(0.0, step_config(), || Ok(())).unwrap();
    let (heated, heat_audit, heat_entropy) = active
        .advance_to_with_entropy(20.0, step_config(), || Ok(()))
        .unwrap();
    let heat = 20.0 * (1000.0 - loss_w);
    near(heat_audit.external_heat_j, heat, 1e-10);
    near(heat_audit.residual_j, 0.0, 1e-8);
    let heated_h = heated.cells()[&CellId(1)].specific_h_j_kg;
    near(mass * (heated_h - initial_h), heat, 1e-8);
    assert_eq!(heated.temperature_k(CellId(1)).unwrap(), 373.15);
    assert!(
        law.state_at_specific_enthalpy(heated_h)
            .unwrap()
            .vapor_mass_fraction()
            > 0.25
    );
    near(
        heat_entropy.thermal_storage_change_j_k,
        heat / 373.15,
        1e-10,
    );
    assert!(heat_entropy.reservoir_entropy_input_j_k < 0.0);
    assert!(heat_entropy.conduction_generation_j_k > 0.0);
    near(heat_entropy.temporal_defect_j_k, 0.0, 1e-10);
    let (cooled, cool_audit, _) = heated
        .advance_to_with_entropy(40.0, step_config(), || Ok(()))
        .unwrap();
    near(cool_audit.external_heat_j, -heat, 1e-10);
    near(cooled.cells()[&CellId(1)].specific_h_j_kg, initial_h, 1e-8);
    assert_eq!(cooled.temperature_k(CellId(1)).unwrap(), 373.15);
    let (single, _, _) = active
        .advance_to_with_entropy(40.0, step_config(), || Ok(()))
        .unwrap();
    assert_eq!(single.cells(), cooled.cells());
    assert_eq!(active.cells()[&CellId(1)].specific_h_j_kg, initial_h);
    assert_eq!(active.time_s(), 0.0);
}

#[test]
fn finite_contact_condenses_hot_steam_and_boils_cooler_water_conservatively() {
    let cold = law(373.15);
    let hot = law(450.0);
    let mass = 0.5;
    let cold_h = enthalpy(cold, 0.2);
    let hot_h = enthalpy(hot, 0.8);
    let initial = DepositionSimulation::new(
        vec![
            birth(1, cold, 0.2, mass, boundary(300.0, 0.0, 0.0)),
            birth(2, hot, 0.8, mass, boundary(300.0, 0.0, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.5,
                    conductance_w_k: 3.0,
                }),
            },
        ],
        &[(MaterialId(1), &cold), (MaterialId(2), &hot)],
        budget(),
    )
    .unwrap();
    let (active, _) = initial.advance_to(0.0, step_config(), || Ok(())).unwrap();
    let (next, audit, entropy) = active
        .advance_to_with_entropy(20.0, step_config(), || Ok(()))
        .unwrap();
    let heat = 20.0 * 3.0 * (450.0 - 373.15);
    let cold_final = cold
        .state_at_specific_enthalpy(next.cells()[&CellId(1)].specific_h_j_kg)
        .unwrap();
    let hot_final = hot
        .state_at_specific_enthalpy(next.cells()[&CellId(2)].specific_h_j_kg)
        .unwrap();
    near(
        mass * (cold_final.specific_enthalpy_j_kg() - cold_h),
        heat,
        1e-8,
    );
    near(
        mass * (hot_final.specific_enthalpy_j_kg() - hot_h),
        -heat,
        1e-8,
    );
    assert!(cold_final.vapor_mass_fraction() > 0.2);
    assert!(hot_final.vapor_mass_fraction() < 0.8);
    assert_eq!(audit.external_heat_j, 0.0);
    near(audit.residual_j, 0.0, 1e-8);
    let produced = heat * (1.0 / 373.15 - 1.0 / 450.0);
    near(entropy.thermal_storage_change_j_k, produced, 1e-10);
    near(entropy.conduction_generation_j_k, produced, 1e-10);
    near(entropy.temporal_defect_j_k, 0.0, 1e-10);
    assert!(!entropy.has_unclassified_sources);
    assert!(produced > 0.0);
    for (law, h0, h1, signed_heat) in [
        (cold, cold_h, cold_final.specific_enthalpy_j_kg(), heat),
        (hot, hot_h, hot_final.specific_enthalpy_j_kg(), -heat),
    ] {
        let calorimetry = law
            .advance_heat(
                law.state_at_specific_enthalpy(h0).unwrap(),
                mass,
                signed_heat,
                1e-8,
                || false,
            )
            .unwrap();
        near(calorimetry.state.specific_enthalpy_j_kg(), h1, 1e-8);
        near(
            calorimetry.internal_energy_change_j + calorimetry.boundary_work_j,
            signed_heat,
            1e-8,
        );
        assert_eq!(
            calorimetry.boundary_work_j.is_sign_positive(),
            signed_heat.is_sign_positive()
        );
    }
}

#[test]
fn thermal_domain_exit_and_cancellation_publish_no_history() {
    let law = law(373.15);
    let initial = DepositionSimulation::new(
        vec![birth(1, law, 0.99, 0.5, boundary(300.0, 0.0, 10_000.0))],
        &[(MaterialId(1), &law)],
        budget(),
    )
    .unwrap();
    let (active, _) = initial.advance_to(0.0, step_config(), || Ok(())).unwrap();
    let before = active.cells().clone();
    assert!(matches!(
        active.advance_to(4.0, step_config(), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert!(matches!(
        active.advance_to(0.01, step_config(), || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(active.cells(), &before);
    assert_eq!(active.time_s(), 0.0);
}
