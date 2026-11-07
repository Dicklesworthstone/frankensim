use super::*;

fn original(thermal: &DepositionSimulation<'_>) -> ConductionObserver<'static> {
    ConductionObserver::new(
        thermal,
        &[(MaterialId(1), &UNIT_K)],
        &[ConductiveLink {
            link_id: 1,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.125,
        }],
        &[face(1.0, 0.5)],
        || Ok(()),
    )
    .unwrap()
    .with_ambient_faces(
        thermal,
        &[ConductiveAmbient {
            cell_id: CellId(1),
            half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.25,
        }],
        || Ok(()),
    )
    .unwrap()
}

fn metrics() -> CurrentConductionGeometry {
    CurrentConductionGeometry {
        links: [(
            1,
            CurrentConductiveLink {
                area_m2: 2.0,
                half_distances_m: [0.25, 0.75],
            },
        )]
        .into(),
        reservoirs: [(
            CellId(0),
            CurrentConductiveBoundary {
                area_m2: 1.5,
                half_distance_m: 0.25,
            },
        )]
        .into(),
        ambient_faces: [(
            CellId(1),
            CurrentConductiveBoundary {
                area_m2: 0.5,
                half_distance_m: 1.0,
            },
        )]
        .into(),
    }
}

#[test]
fn current_faces_drive_actual_heat_and_cfl_without_replacing_reference_inventory() {
    let initial = simulation(
        cells(&[300.0, 400.0], 1.0, Some(600.0)),
        &SENSIBLE,
        vec![],
        100,
    );
    let original = original(&initial);
    let current = original
        .with_current_geometry(&initial, metrics(), || Ok(()))
        .unwrap();
    let rates = initial.power_balance_observed(&current, || Ok(())).unwrap();
    // Independent normal resistance calculation. Only original pairwise heat
    // contributes internally; the two finite-bath exchanges remain external.
    let pair = 100.0 * 2.0 / (0.25 + 0.125 + 0.75);
    let hot_bath = 300.0 * 1.5 / (0.25 + 0.5);
    let cool_bath = -100.0 * 0.5 / (1.0 + 0.25);
    near(rates.net_power_w[&CellId(0)], pair + hot_bath, 1e-12);
    near(rates.net_power_w[&CellId(1)], -pair + cool_bath, 1e-12);
    near(rates.external_power_w, hot_bath + cool_bath, 1e-12);
    let dt = 0.025;
    let (next, retained, audit, entropy) = initial
        .advance_observed_with_entropy(&current, dt, config(dt), || Ok(()))
        .unwrap();
    assert_eq!(audit.substeps, 1);
    near(
        next.temperature_k(CellId(0)).unwrap(),
        300.0 + dt * (pair + hot_bath),
        1e-12,
    );
    near(
        next.temperature_k(CellId(1)).unwrap(),
        400.0 + dt * (-pair + cool_bath),
        1e-12,
    );
    near(audit.external_heat_j, dt * (hot_bath + cool_bath), 1e-12);
    near(audit.stored_energy_change_j, audit.external_heat_j, 1e-10);
    assert!(entropy.conduction_generation_j_k > 0.0);
    assert!(!entropy.has_unclassified_sources);
    for (&id, old) in initial.cells() {
        let new = &next.cells()[&id];
        assert_eq!(new.geometry, old.geometry);
        assert_eq!(new.mass_kg.to_bits(), old.mass_kg.to_bits());
        assert_eq!(new.boundary, old.boundary);
        assert_eq!(new.material_id, old.material_id);
    }
    assert_eq!(next.links(), initial.links());
    assert!(original.current_geometry().is_none());
    assert_eq!(retained.current_geometry(), Some(&metrics()));

    // Native CFL must see current G as well as its flux. A much larger current
    // area forces multiple steps while leaving the declared maximum unchanged.
    let mut fast = metrics();
    fast.links.get_mut(&1).unwrap().area_m2 *= 1000.0;
    let fast = original
        .with_current_geometry(&initial, fast, || Ok(()))
        .unwrap();
    let (_, _, audit) = initial
        .advance_observed(&fast, dt, config(dt), || Ok(()))
        .unwrap();
    assert!(audit.substeps > 1);
}

#[test]
fn incomplete_foreign_or_conflicting_metrics_refuse_and_retry_preserves_original_owner() {
    let initial = simulation(
        cells(&[300.0, 400.0], 1.0, Some(600.0)),
        &SENSIBLE,
        vec![],
        100,
    );
    let original = original(&initial);
    let mut missing = metrics();
    missing.ambient_faces.clear();
    assert!(
        original
            .with_current_geometry(&initial, missing, || Ok(()))
            .is_err()
    );
    let mut foreign = metrics();
    let face = foreign.links.remove(&1).unwrap();
    foreign.links.insert(99, face);
    assert!(
        original
            .with_current_geometry(&initial, foreign, || Ok(()))
            .is_err()
    );
    let mut invalid = metrics();
    invalid
        .reservoirs
        .get_mut(&CellId(0))
        .unwrap()
        .half_distance_m = f64::NAN;
    assert!(
        original
            .with_current_geometry(&initial, invalid, || Ok(()))
            .is_err()
    );
    assert!(matches!(
        original.with_current_geometry(&initial, metrics(), || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    let layered = original
        .clone()
        .with_planar_phase_layers(
            &initial,
            &[PlanarPhaseLayer {
                cell_id: CellId(0),
                footprint_area_m2: 1.0,
                density_kg_m3: 1.0,
            }],
            || Ok(()),
        )
        .unwrap();
    assert!(
        layered
            .with_current_geometry(&initial, metrics(), || Ok(()))
            .is_err()
    );
    let current = original
        .with_current_geometry(&initial, metrics(), || Ok(()))
        .unwrap();
    assert!(
        current
            .clone()
            .with_planar_phase_layers(
                &initial,
                &[PlanarPhaseLayer {
                    cell_id: CellId(0),
                    footprint_area_m2: 1.0,
                    density_kg_m3: 1.0,
                }],
                || Ok(())
            )
            .is_err()
    );
    let mut changed = metrics();
    changed.links.get_mut(&1).unwrap().half_distances_m[0] *= 2.0;
    near(
        metrics().relative_difference(&changed, || Ok(())).unwrap(),
        0.5,
        0.0,
    );
    assert_eq!(
        metrics()
            .relative_difference(&metrics(), || Ok(()))
            .unwrap(),
        0.0
    );
    assert!(original.current_geometry().is_none());
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 300.0);
    assert!(
        initial
            .advance_observed(&current, 0.025, config(0.025), || Ok(()))
            .is_ok()
    );
}
