//! The actual heat timeline owns storage, heat and immutable failed advances.

use super::*;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn birth(id: u64, mass_kg: f64, specific_h_j_kg: f64, source_w: f64) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: 1,
            material_id: MaterialId(1),
            mass_kg,
            specific_h_j_kg,
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64, 1.0, 0.0],
                volume_m3: mass_kg / 1000.0,
                surface_area_m2: 1.0,
            },
            boundary: boundary(source_w),
        }),
    }
}
fn budget() -> DepositionBudget {
    DepositionBudget {
        max_cells: 2,
        max_links: 1,
        max_events: 3,
        max_substeps: 256,
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-6,
    }
}

#[test]
fn thermal_storage_retains_bounds_and_heat_timeline_crosses_miscibility_both_directions() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let storage = thermal
        .thermal_storage([2.0, 3.0], controls(), || false)
        .unwrap();
    assert_eq!(storage.sources(), thermal.sources());
    assert_eq!(storage.temperature_domain(), controls().temperature_bracket);
    assert_eq!(storage.caloric_identity(), None);
    assert_eq!(
        storage.specific_entropy_change_j_kg_k(0.0, 1.0).unwrap(),
        None
    );
    assert!(storage.advected_composition().is_none());
    let [minimum, maximum] = storage.specific_enthalpy_domain_j_kg();
    assert_eq!(storage.temperature_k(minimum).unwrap(), 290.0);
    assert_eq!(storage.temperature_k(maximum).unwrap(), 600.0);
    for invalid in [minimum.next_down(), maximum.next_up(), f64::NAN] {
        assert_eq!(
            storage.temperature_k(invalid),
            Err(DepositionError::MaterialDomain)
        );
    }
    assert!(
        storage
            .specific_enthalpy_j_kg(Temperature::new(289.0), || false)
            .is_err()
    );
    assert!(matches!(
        thermal.thermal_storage([2.0, 3.0], controls(), || true),
        Err(BinaryLiquidCaloricError::Cancelled)
    ));
    let bound = storage.maximum_temperature_slope();
    assert!(bound >= 0.336 / 500.0 && bound < (0.336 / 500.0) * (1.0 + 1e-13));
    for (a, b) in [
        (300.0, 301.0),
        (430.0, 431.0),
        (449.0, 451.0),
        (500.0, 501.0),
    ] {
        let ha = storage
            .specific_enthalpy_j_kg(Temperature::new(a), || false)
            .unwrap();
        let hb = storage
            .specific_enthalpy_j_kg(Temperature::new(b), || false)
            .unwrap();
        assert!((b - a) / (hb - ha) <= bound * (1.0 + 1e-10));
    }

    // Independent chi=3 binodal root; Cp is 100 J/(mol K) for both liquids.
    // Double the reference inventory's mass to exercise an intensive law.
    let a = 0.070_720_181_679_944_8;
    let heat_j = 10.0 * (100.0 * 200.0 + R * 900.0 * (0.4 * 0.6 - a * (1.0 - a)));
    let initial_h = storage
        .specific_enthalpy_j_kg(Temperature::new(300.0), || false)
        .unwrap();
    let initial = DepositionSimulation::new(
        vec![
            birth(1, 0.672, initial_h, heat_j),
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), boundary(-heat_j)),
            },
        ],
        &[(MaterialId(1), &storage)],
        budget(),
    )
    .unwrap();
    let (initial, _) = initial.advance_to(0.0, config(0.25), || Ok(())).unwrap();
    let (warm, heating) = initial.advance_to(1.0, config(0.25), || Ok(())).unwrap();
    assert!((warm.temperature_k(CellId(1)).unwrap() - 500.0).abs() < 1e-8);
    assert!((heating.external_heat_j - heat_j).abs() < 1e-7);
    assert!(heating.residual_j.abs() < 1e-7);
    let state = storage
        .caloric_state_at_specific_enthalpy(warm.cells()[&CellId(1)].specific_h_j_kg, || false)
        .unwrap();
    assert_eq!(state.equilibrium().unwrap().phases().len(), 1);
    let (cold, cooling) = warm.advance_to(2.0, config(0.25), || Ok(())).unwrap();
    assert!((cold.temperature_k(CellId(1)).unwrap() - 300.0).abs() < 1e-8);
    assert!(cooling.residual_j.abs() < 1e-7);
    let state = storage
        .caloric_state_at_specific_enthalpy(cold.cells()[&CellId(1)].specific_h_j_kg, || false)
        .unwrap();
    assert_eq!(state.equilibrium().unwrap().phases().len(), 2);
    assert_eq!(cold.mass_kg(), initial.mass_kg());
    assert_eq!(
        cold.cells()[&CellId(1)].geometry,
        initial.cells()[&CellId(1)].geometry
    );
    let before = cold.cells().clone();
    assert!(matches!(
        cold.advance_to(4.0, config(0.25), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert_eq!(cold.cells(), &before);
    assert_eq!(cold.time_s(), 2.0);
}

#[test]
fn original_thermal_links_exchange_miscibility_energy_without_external_heat() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let storage = thermal
        .thermal_storage([2.0, 3.0], controls(), || false)
        .unwrap();
    let h = [350.0, 590.0].map(|t| {
        storage
            .specific_enthalpy_j_kg(Temperature::new(t), || false)
            .unwrap()
    });
    let expected = storage.temperature_k(0.5 * h[0] + 0.5 * h[1]).unwrap();
    assert!(expected > 450.0);
    let simulation = DepositionSimulation::new(
        vec![
            birth(1, 0.336, h[0], 0.0),
            birth(2, 0.336, h[1], 0.0),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.1,
                    conductance_w_k: 10.0,
                }),
            },
        ],
        &[(MaterialId(1), &storage)],
        budget(),
    )
    .unwrap();
    let (initial, _) = simulation.advance_to(0.0, config(5.0), || Ok(())).unwrap();
    let (final_state, audit) = initial.advance_to(500.0, config(5.0), || Ok(())).unwrap();
    assert_eq!(audit.external_heat_j, 0.0);
    assert!(audit.residual_j.abs() < 1e-7);
    assert!((final_state.energy_j() - initial.energy_j()).abs() < 1e-7);
    for id in [CellId(1), CellId(2)] {
        assert!((final_state.temperature_k(id).unwrap() - expected).abs() < 1e-6);
        let state = storage
            .caloric_state_at_specific_enthalpy(final_state.cells()[&id].specific_h_j_kg, || false)
            .unwrap();
        assert_eq!(state.equilibrium().unwrap().phases().len(), 1);
        assert_eq!(state.component_moles(), [2.0, 3.0]);
    }
    assert!((initial.temperature_k(CellId(1)).unwrap() - 350.0).abs() < 1e-8);
}
