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
    assert!(
        advance_deposition_with_gas(
            &future,
            surface(&binding),
            &initial_gas,
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
