//! Batched coverage preserves the original scalar area and actual heat equations.
use fs_thermal::enthalpy::*;

struct Law;
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > 0.0 { Ok(h) } else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
}
static LAW: Law = Law;
fn config() -> StepConfig { StepConfig { max_dt_s: 0.01, energy_tolerance_j: 1e-8 } }
fn initial(count: u64) -> DepositionSimulation<'static> {
    let events = (1..=count).map(|id| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: CellId(id), segment_id: id, body_id: 1, material_id: MaterialId(1),
        mass_kg: 1.0, specific_h_j_kg: 300.0,
        geometry: CellGeometry { start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 1.0, 0.0], volume_m3: 0.001, surface_area_m2: 1.0 },
        boundary: Boundary { ambient_k: 200.0, htc_w_m2_k: 2.0, source_w: 0.0,
            reservoir: Some(Reservoir { area_m2: 0.125, conductance_w_k: 0.0, temperature_k: 300.0 }) },
    }) }).collect();
    DepositionSimulation::new(events, &[(MaterialId(1), &LAW)], DepositionBudget {
        max_cells: count as usize, max_links: 3 * count as usize,
        max_events: 4 * count as usize, max_substeps: 100,
    }).unwrap().advance_to(0.0, config(), || Ok(())).unwrap().0
}
fn patch(id: u64, a: u64, b: u64, area: [f64; 2]) -> ContactPatch {
    ContactPatch { link: ThermalLink { id, a: CellId(a), b: CellId(b),
        area_m2: 0.02, conductance_w_k: if id % 2 == 0 { 0.0 } else { 0.125 } },
        covered_area_m2: area }
}
fn install(state: &DepositionSimulation<'static>, patches: &[ContactPatch]) -> DepositionSimulation<'static> {
    let updates: Vec<_> = patches.iter().map(|&p| ContactPatchUpdate {
        id: p.link.id, expected: state.contact_patch(p.link.id), replacement: Some(p),
    }).collect();
    state.replace_contact_patches(&updates, || Ok(())).unwrap()
}
fn compare_scalar(state: &DepositionSimulation<'_>, ids: &[CellId]) {
    let areas = state.exposed_areas_m2(ids, || Ok(())).unwrap();
    for &id in ids {
        assert_eq!(areas[&id].to_bits(), state.exposed_area_m2(id).unwrap().to_bits());
    }
    let reverse: Vec<_> = ids.iter().rev().copied().collect();
    assert_eq!(areas, state.exposed_areas_m2(&reverse, || Ok(())).unwrap());
}

#[test]
fn selected_areas_match_scalar_bits_with_reservoirs_asymmetry_and_insulating_contacts() {
    for count in [3, 17, 257] {
        let patches: Vec<_> = (1..=count).map(|id|
            patch(id, id, id % count + 1, [0.031 + id as f64 * 1e-5, 0.047])).collect();
        let state = install(&initial(count), &patches);
        let ids: Vec<_> = (1..=count).rev().map(CellId).collect();
        compare_scalar(&state, &ids);
        compare_scalar(&state, &ids.iter().step_by(2).copied().collect::<Vec<_>>());
        let swapped: Vec<_> = patches.iter().map(|p| ContactPatch {
            link: ThermalLink { a: p.link.b, b: p.link.a, ..p.link },
            covered_area_m2: [p.covered_area_m2[1], p.covered_area_m2[0]],
        }).collect();
        let other = install(&initial(count), &swapped);
        assert_eq!(state.exposed_areas_m2(&ids, || Ok(())).unwrap(),
            other.exposed_areas_m2(&ids, || Ok(())).unwrap());
        let mut polls = 0;
        assert!(state.exposed_areas_m2(&[], || { polls += 1; Ok(()) }).unwrap().is_empty());
        assert_eq!(polls, 2, "empty selection must not walk contacts");
    }
}

#[test]
fn final_contact_repartition_owns_actual_ambient_heat_without_changing_material_inventory() {
    let original = install(&initial(3), &[patch(10, 1, 2, [0.6, 0.2]), patch(20, 1, 3, [0.2, 0.4])]);
    let replacements = [patch(20, 1, 3, [0.6, 0.3]), patch(10, 1, 2, [0.2, 0.5])];
    let next = install(&original, &replacements);
    assert_eq!(original.cells(), next.cells());
    assert_eq!(original.mass_kg(), next.mass_kg());
    assert_eq!(original.energy_j(), next.energy_j());
    assert_eq!(next.time_s(), 0.0);
    compare_scalar(&next, &[CellId(1), CellId(2), CellId(3)]);
    let areas = next.exposed_areas_m2(&[CellId(1), CellId(2), CellId(3)], || Ok(())).unwrap();
    let (cooled, audit) = next.advance_to(0.01, config(), || Ok(())).unwrap();
    let mut external = 0.0;
    for (&id, &area) in &areas {
        let heat = 0.01 * area * 2.0 * (200.0 - 300.0);
        external += heat;
        assert!((cooled.cells()[&id].specific_h_j_kg - (300.0 + heat)).abs() < 1e-12);
    }
    assert!((audit.external_heat_j - external).abs() < 1e-12);
    assert!(audit.residual_j.abs() < 1e-10);
    let too_large = [ContactPatchUpdate { id: 20, expected: next.contact_patch(20),
        replacement: Some(patch(20, 1, 3, [0.9, 0.3])) }];
    assert!(matches!(next.replace_contact_patches(&too_large, || Ok(())),
        Err(DepositionError::InvalidInput("contacts exceed surface area"))));
    assert_eq!(next.contact_patch(20), Some(replacements[0]));
}

#[test]
fn selection_and_late_cancellation_publish_no_partial_area_map_or_contact_change() {
    let original = initial(3);
    for bad in [vec![CellId(1), CellId(1)], vec![CellId(99)]] {
        assert!(original.exposed_areas_m2(&bad, || Ok(())).is_err());
    }
    assert!(matches!(original.exposed_areas_m2(&[CellId(1); 4], || Ok(())),
        Err(DepositionError::Budget("exposed surface selection"))));
    let updates = [ContactPatchUpdate { id: 10, expected: None,
        replacement: Some(patch(10, 1, 2, [0.3, 0.4])) }];
    let mut polls = 0;
    let expected = original.replace_contact_patches(&updates, || { polls += 1; Ok(()) }).unwrap();
    let mut calls = 0;
    assert!(matches!(original.replace_contact_patches(&updates, || {
        calls += 1;
        if calls == polls { Err(DepositionError::Cancelled) } else { Ok(()) }
    }), Err(DepositionError::Cancelled)));
    assert!(original.links().is_empty());
    let retried = original.replace_contact_patches(&updates, || Ok(())).unwrap();
    assert_eq!(retried.contact_patch(10), expected.contact_patch(10));
    let ids = [CellId(1), CellId(2), CellId(3)];
    polls = 0;
    expected.exposed_areas_m2(&ids, || { polls += 1; Ok(()) }).unwrap();
    calls = 0;
    assert!(matches!(expected.exposed_areas_m2(&ids, || {
        calls += 1;
        if calls == polls { Err(DepositionError::Cancelled) } else { Ok(()) }
    }), Err(DepositionError::Cancelled)));
    compare_scalar(&expected, &ids);
}

#[test]
fn a_late_overcovered_cell_refuses_the_complete_repartition() {
    let original = initial(3);
    let updates = [
        ContactPatchUpdate { id: 10, expected: None, replacement: Some(patch(10, 1, 2, [0.3, 0.4])) },
        ContactPatchUpdate { id: 20, expected: None, replacement: Some(patch(20, 2, 3, [0.3, 0.9])) },
    ];
    assert!(matches!(original.replace_contact_patches(&updates, || Ok(())),
        Err(DepositionError::InvalidInput("contacts exceed surface area"))));
    assert!(original.links().is_empty());
    for id in 1..=3 { assert_eq!(original.exposed_area_m2(CellId(id)), Ok(0.875)); }
    original.replace_contact_patches(&updates[..1], || Ok(())).unwrap();
}
