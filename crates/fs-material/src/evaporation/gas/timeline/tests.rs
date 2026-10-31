//! G0 original-owner coupling and G1 full heat/finite-humidity refinement.
use super::*;
use crate::thermochemical_gas::energy::{
    GasEnergySolve, parcel::transfer::species::boundary::GasSpeciesBoundarySolve,
};
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, MaterialId, Reservoir, ThermalLaw,
    TimedEvent,
};
use fs_thermal::evaporation::EvaporationProperties;
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, mixture::FrozenIdealGasMixtureModelV1,
};

const M: f64 = 0.018015;
const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / M;
const LATENT: f64 = 200000.0;
fn h(t: f64) -> f64 {
    R * (4.0 * t + 0.0005 * t * t) - LATENT
}
struct Liquid;
impl ThermalLaw for Liquid {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let c = (h + LATENT) / R;
        let t = 2.0 * c / (4.0 + (16.0 + 0.002 * c).sqrt());
        if (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / (4.2 * R)
    }
}
impl EvaporationLaw for Liquid {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 900.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let p = 1500.0 * fs_math::det::exp(LATENT / R * (1.0 / 300.0 - 1.0 / t));
        Ok(EvaporationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p * LATENT / (R * t * t),
            latent_heat_j_kg: LATENT,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: M,
        })
    }
}
static LIQUID: Liquid = Liquid;
fn id() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", M, 4.0, 0.001, reference),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.03)],
        Some(&id()),
    )
    .unwrap()
}
fn gas(catalog: &ThermochemicalGasModel, volume: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        catalog
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap(),
            )
            .unwrap(),
        point(300.0, 100000.0),
        volume,
    )
    .unwrap()
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn birth(time_s: f64, source_w: f64) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(DepositedCell {
            id: CellId(1),
            segment_id: 19,
            body_id: 2,
            material_id: MaterialId(1),
            mass_kg: 0.01,
            specific_h_j_kg: h(300.0),
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [0.01, 0.0, 0.0],
                volume_m3: 1e-5,
                surface_area_m2: 1e-3,
            },
            boundary: boundary(source_w),
        }),
    }
}
fn thermal(events: Vec<TimedEvent>, steps: usize) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &LIQUID)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 4,
            max_substeps: steps,
        },
    )
    .unwrap()
}
fn heat(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-6,
    }
}
fn controls() -> GasEvaporationSolve {
    GasEvaporationSolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 900.0,
                energy_tolerance_j_kg: 1e-8,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        maximum_inventory_fraction: 0.25,
        maximum_temperature_increment_k: 10.0,
        latent_compatibility_tolerance_j_kg: 1e-5,
    }
}
fn film_controls() -> GasFilmHeatSolve {
    GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [h(250.0), h(350.0)],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}
fn surface(binding: &GasLiquidBinding) -> GasEvaporationTimelineSurface<'_> {
    GasEvaporationTimelineSurface {
        cell_id: CellId(1),
        law: &LIQUID,
        binding,
        pure_liquid_enthalpy_range_j_kg: [h(250.0), h(350.0)],
        mass_transfer_coefficient_m_s: 0.02,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e} ({tolerance:e})");
}

#[test]
fn heat_evaporation_and_finite_humidity_share_original_energy_and_mass() {
    let catalog = catalog(0.0);
    let initial_gas = gas(&catalog, 1e-4);
    let thermal = thermal(vec![birth(0.0, 2.0)], 1000)
        .advance_to(0.0, heat(0.05), || Ok(()))
        .unwrap()
        .0;
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let result = advance_deposition_with_gas(
        &thermal,
        surface(&binding),
        &initial_gas,
        &catalog,
        2.0,
        heat(0.05),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = result.audit();
    near(audit.thermal.external_heat_j, 4.0, 1e-12);
    assert!(audit.transferred_mass_kg > 0.0);
    assert!(result.thermal().temperature_k(CellId(1)).unwrap() > 300.0);
    assert!(result.gas().species_partial_pressure_pa(&id()).unwrap() > 0.0);
    near(
        result.thermal().mass_kg() - thermal.mass_kg() + result.gas().mass_kg()
            - initial_gas.mass_kg(),
        0.0,
        1e-14,
    );
    let initial_u = initial_gas.mass_kg() * initial_gas.state().specific_internal_energy_j_kg();
    let final_u = result.gas().mass_kg() * result.gas().state().specific_internal_energy_j_kg();
    near(
        result.thermal().energy_j() - thermal.energy_j() + final_u - initial_u,
        4.0,
        1e-6,
    );
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j,
        4.0,
        1e-6,
    );
    assert!(audit.energy_residual_j.abs() < 1e-6);
    assert_eq!(
        result.thermal().cells()[&CellId(1)].geometry,
        thermal.cells()[&CellId(1)].geometry
    );
    assert!(std::ptr::eq(result.catalog(), &catalog));
    // A larger physical gas inventory supplies less back-pressure.
    let large_gas = gas(&catalog, 0.1);
    let large = advance_deposition_with_gas(
        &thermal,
        surface(&binding),
        &large_gas,
        &catalog,
        2.0,
        heat(0.05),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(large.audit().transferred_mass_kg > audit.transferred_mass_kg);
}

#[test]
fn births_and_endpoint_covering_preserve_phase_timing_and_reference_offsets() {
    let mut endpoints = Vec::new();
    for reference in [0.0, 1e8] {
        let catalog = catalog(reference);
        let initial_gas = gas(&catalog, 1e-3);
        let binding = GasLiquidBinding {
            species: id(),
            liquid_enthalpy_offset_j_kg: reference,
        };
        let mut covered = boundary(0.0);
        covered.reservoir = Some(Reservoir {
            area_m2: 1e-3,
            conductance_w_k: 0.0,
            temperature_k: 300.0,
        });
        let original = thermal(
            vec![
                birth(0.1, 0.0),
                TimedEvent {
                    time_s: 0.2,
                    event: Event::Boundary(CellId(1), covered),
                },
            ],
            100,
        );
        let born = advance_deposition_with_gas(
            &original,
            surface(&binding),
            &initial_gas,
            &catalog,
            0.1,
            heat(0.1),
            controls(),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(born.audit().phase_steps, 0);
        assert_eq!(born.gas(), &initial_gas);
        let result = advance_deposition_with_gas(
            &original,
            surface(&binding),
            &initial_gas,
            &catalog,
            0.4,
            heat(0.1),
            controls(),
            || Ok(()),
        )
        .unwrap();
        assert!(result.audit().transferred_mass_kg > 0.0);
        assert_eq!(
            result
                .audit()
                .last_step
                .as_ref()
                .unwrap()
                .transferred_mass_kg,
            0.0
        );
        assert_eq!(result.thermal().exposed_area_m2(CellId(1)), Ok(0.0));
        near(
            result.audit().incoming_energy_j,
            0.01 * (h(300.0) + reference),
            1e-8,
        );
        assert!(result.audit().energy_residual_j.abs() < 1e-6);
        endpoints.push((
            result.thermal().temperature_k(CellId(1)).unwrap(),
            result.gas().state().state().temperature,
            result.audit().transferred_mass_kg,
        ));
    }
    near(endpoints[0].0, endpoints[1].0, 1e-10);
    near(endpoints[0].1, endpoints[1].1, 1e-8);
    near(endpoints[0].2, endpoints[1].2, 1e-14);
}

#[test]
fn full_timeline_cancellation_and_work_refusal_publish_neither_owner() {
    let catalog = catalog(0.0);
    let initial_gas = gas(&catalog, 1e-3);
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let original = thermal(vec![birth(0.0, 2.0)], 100);
    let before = initial_gas.clone();
    let mut count = 0;
    advance_deposition_with_gas(
        &original,
        surface(&binding),
        &initial_gas,
        &catalog,
        0.2,
        heat(0.05),
        controls(),
        || {
            count += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut seen = 0;
    let result = advance_deposition_with_gas(
        &original,
        surface(&binding),
        &initial_gas,
        &catalog,
        0.2,
        heat(0.05),
        controls(),
        || {
            seen += 1;
            if seen == count {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        result,
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(initial_gas, before);
    assert!(original.cells().is_empty());
    assert_eq!(original.time_s(), 0.0);
    let limited = thermal(vec![birth(0.0, 2.0)], 1);
    assert!(matches!(
        advance_deposition_with_gas(
            &limited,
            surface(&binding),
            &initial_gas,
            &catalog,
            0.2,
            heat(0.05),
            controls(),
            || Ok(())
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert!(limited.cells().is_empty());
    assert_eq!(initial_gas, before);
    let future = thermal(vec![birth(1.0, 2.0)], 100);
    let mut invalid = controls();
    invalid.maximum_inventory_fraction = 0.0;
    assert!(
        advance_deposition_with_gas(
            &future,
            surface(&binding),
            &initial_gas,
            &catalog,
            0.5,
            heat(0.05),
            invalid,
            || Ok(())
        )
        .is_err()
    );
    let changed_catalog = self::catalog(1e6);
    // A dry parcel has no original vapor card to compare. Retain actual vapor
    // inventory so this tests substitution of an owned source.
    let humid_gas = ClosedGasParcel::new(catalog.clone(), point(300.0, 100000.0), 1e-3).unwrap();
    let before_humid = humid_gas.clone();
    assert!(
        advance_deposition_with_gas(
            &future,
            surface(&binding),
            &humid_gas,
            &changed_catalog,
            0.0,
            heat(0.05),
            controls(),
            || Ok(())
        )
        .is_err()
    );
    assert!(future.cells().is_empty());
    assert_eq!(initial_gas, before);
    assert_eq!(humid_gas, before_humid);

    let mut invalid_heat = film_controls();
    invalid_heat.maximum_iterations = 0;
    let mut different_branch = film_controls();
    different_branch.pure_liquid_enthalpy_range_j_kg[0] += 1.0;
    for (coefficient, heat_controls) in [
        (-1.0, film_controls()),
        (0.0, invalid_heat),
        (1.0, different_branch),
    ] {
        assert!(
            advance_deposition_with_gas_and_heat(
                &future,
                surface(&binding),
                &initial_gas,
                &catalog,
                0.0,
                heat(0.05),
                controls(),
                coefficient,
                heat_controls,
                || Ok(()),
            )
            .is_err()
        );
    }
}

#[test]
fn complete_heat_and_finite_gas_feedback_converges_under_time_refinement() {
    let catalog = catalog(0.0);
    let initial_gas = gas(&catalog, 1e-4);
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let thermal = thermal(vec![birth(0.0, 2.0)], 1000);
    let run = |dt| {
        advance_deposition_with_gas(
            &thermal,
            surface(&binding),
            &initial_gas,
            &catalog,
            2.0,
            heat(dt),
            controls(),
            || Ok(()),
        )
        .unwrap()
        .gas()
        .species_partial_pressure_pa(&id())
        .unwrap()
    };
    let reference = run(0.00625);
    let a = (run(0.2) - reference).abs();
    let b = (run(0.1) - reference).abs();
    let c = (run(0.05) - reference).abs();
    assert!(a > 1e-5);
    assert!(b < 0.7 * a, "{a}, {b}");
    assert!(c < 0.7 * b, "{b}, {c}");
}

#[test]
fn bounded_local_retry_recovers_finite_backpressure_without_replaying_heat_or_events() {
    let catalog = catalog(0.0);
    let initial_gas = gas(&catalog, 1e-5);
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let mut covered = boundary(2.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 1e-3,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let original = thermal(
        vec![
            birth(0.0, 2.0),
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        1000,
    );
    let result = advance_deposition_with_gas(
        &original,
        surface(&binding),
        &initial_gas,
        &catalog,
        2.0,
        heat(2.0),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = result.audit();
    assert!(audit.rejected_intervals > 0);
    assert!(audit.attempted_thermal_substeps >= audit.thermal.substeps + audit.rejected_intervals);
    assert!(audit.attempted_thermal_substeps <= original.work_budget().max_substeps);
    assert!(audit.transferred_mass_kg > 0.0);
    assert_eq!(result.thermal().time_s(), 2.0);
    assert_eq!(audit.thermal.events, 2);
    assert_eq!(result.thermal().exposed_area_m2(CellId(1)), Ok(0.0));
    assert_eq!(audit.last_step.as_ref().unwrap().transferred_mass_kg, 0.0);
    near(audit.thermal.external_heat_j, 4.0, 1e-12);
    near(audit.incoming_energy_j, 0.01 * h(300.0), 1e-10);
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j - audit.incoming_energy_j - 4.0,
        0.0,
        1e-6,
    );
    near(audit.mass_residual_kg, 0.0, 1e-14);
    assert!(audit.energy_residual_j.abs() <= 1e-6);
    assert!(original.cells().is_empty());
}

#[test]
fn rejected_intervals_spend_work_and_late_retry_cancellation_preserves_both_inputs() {
    let catalog = catalog(0.0);
    let initial_gas = gas(&catalog, 1e-5);
    let before_gas = initial_gas.clone();
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let limited = thermal(vec![birth(0.0, 2.0)], 1);
    assert!(matches!(
        advance_deposition_with_gas(
            &limited,
            surface(&binding),
            &initial_gas,
            &catalog,
            2.0,
            heat(2.0),
            controls(),
            || Ok(())
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert!(limited.cells().is_empty());
    assert_eq!(limited.time_s(), 0.0);
    assert_eq!(initial_gas, before_gas);

    let original = thermal(vec![birth(0.0, 2.0)], 1000);
    let mut polls = 0;
    let accepted = advance_deposition_with_gas(
        &original,
        surface(&binding),
        &initial_gas,
        &catalog,
        2.0,
        heat(2.0),
        controls(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    assert!(accepted.audit().rejected_intervals > 0);
    let mut seen = 0;
    assert!(matches!(
        advance_deposition_with_gas(
            &original,
            surface(&binding),
            &initial_gas,
            &catalog,
            2.0,
            heat(2.0),
            controls(),
            || {
                seen += 1;
                if seen == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert!(original.cells().is_empty());
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(initial_gas, before_gas);
}

#[test]
fn sensible_heat_retries_are_reciprocal_and_stop_at_the_covering_event() {
    let catalog = catalog(0.0);
    let initial_gas = ClosedGasParcel::new(
        gas(&catalog, 1e-3).model().clone(),
        point(340.0, 100000.0),
        1e-3,
    )
    .unwrap();
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 1e-3,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let original = thermal(
        vec![
            birth(0.0, 0.0),
            TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        1000,
    )
    .advance_to(0.0, heat(1.0), || Ok(()))
    .unwrap()
    .0;
    let mut selected = surface(&binding);
    selected.mass_transfer_coefficient_m_s = 0.0;
    let result = advance_deposition_with_gas_and_heat(
        &original,
        selected,
        &initial_gas,
        &catalog,
        1.0,
        heat(1.0),
        controls(),
        1000.0,
        film_controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = result.audit();
    assert!(audit.rejected_intervals > 0);
    assert!(audit.attempted_thermal_substeps >= audit.thermal.substeps + audit.rejected_intervals);
    assert!(audit.sensible_heat_into_liquid_j > 0.0);
    assert_eq!(audit.transferred_mass_kg, 0.0);
    assert_eq!(audit.energy_into_gas_j, 0.0);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.thermal.events, 1);
    assert_eq!(result.thermal().exposed_area_m2(CellId(1)), Ok(0.0));
    assert_eq!(result.thermal().mass_kg(), original.mass_kg());
    assert_eq!(result.gas().mass_kg(), initial_gas.mass_kg());
    assert!(result.thermal().temperature_k(CellId(1)).unwrap() > 300.0);
    assert!(result.gas().state().state().temperature < 340.0);
    near(
        result.thermal().energy_j() - original.energy_j(),
        audit.sensible_heat_into_liquid_j,
        1e-8,
    );
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j,
        0.0,
        1e-6,
    );
    near(
        audit.gas_energy_change_j,
        -audit.sensible_heat_into_liquid_j,
        1e-6,
    );
    assert!(audit.energy_residual_j.abs() <= 1e-6);
    let last = audit.last_step.as_ref().unwrap();
    assert_eq!(last.surface_temperature_k[0], last.surface_temperature_k[1]);
    assert_eq!(last.gas_temperature_k[0], last.gas_temperature_k[1]);
    assert_eq!(original.exposed_area_m2(CellId(1)), Ok(1e-3));
}

#[test]
fn complete_sensible_heat_phase_feedback_refines_and_zero_heat_preserves_old_path() {
    let catalog = catalog(0.0);
    let initial_gas = ClosedGasParcel::new(
        gas(&catalog, 1e-3).model().clone(),
        point(340.0, 100000.0),
        1e-3,
    )
    .unwrap();
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let original = thermal(vec![birth(0.0, 2.0)], 1000)
        .advance_to(0.0, heat(0.05), || Ok(()))
        .unwrap()
        .0;
    let run = |dt, coefficient| {
        let result = advance_deposition_with_gas_and_heat(
            &original,
            surface(&binding),
            &initial_gas,
            &catalog,
            2.0,
            heat(dt),
            controls(),
            coefficient,
            film_controls(),
            || Ok(()),
        )
        .unwrap();
        let audit = result.audit();
        near(audit.thermal.external_heat_j, 4.0, 1e-12);
        near(
            audit.thermal_energy_change_j + audit.gas_energy_change_j,
            4.0,
            1e-6,
        );
        near(audit.mass_residual_kg, 0.0, 1e-14);
        assert!(audit.energy_residual_j.abs() <= 1e-6);
        assert!(audit.transferred_mass_kg > 0.0);
        if coefficient > 0.0 {
            assert!(audit.sensible_heat_into_liquid_j > 0.0);
        }
        result
    };
    let reference = run(0.0125, 100.0).gas().state().state().temperature;
    let a = (run(0.2, 100.0).gas().state().state().temperature - reference).abs();
    let b = (run(0.1, 100.0).gas().state().state().temperature - reference).abs();
    let c = (run(0.05, 100.0).gas().state().state().temperature - reference).abs();
    assert!(a > 1e-5);
    assert!(b < 0.75 * a, "{a}, {b}");
    assert!(c < 0.75 * b, "{b}, {c}");

    let old = advance_deposition_with_gas(
        &original,
        surface(&binding),
        &initial_gas,
        &catalog,
        2.0,
        heat(0.05),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let disabled = run(0.05, 0.0);
    assert_eq!(disabled.gas(), old.gas());
    assert_eq!(disabled.thermal().cells(), old.thermal().cells());
    assert_eq!(disabled.audit(), old.audit());
    assert_eq!(disabled.audit().sensible_heat_into_liquid_j, 0.0);
}
