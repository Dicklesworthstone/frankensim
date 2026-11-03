//! G0/G1/G3/G4: deposited mass drives the original resistance and heat storage.
use super::*;
use crate::enthalpy::sublimation::external::SolidCellExchange;
use crate::sublimation::SolidSurfaceState;

fn layer(id: u64) -> PlanarPhaseLayer {
    PlanarPhaseLayer {
        cell_id: CellId(id),
        footprint_area_m2: 1.0,
        density_kg_m3: 1.0,
    }
}

// Prescribed conservative mass ingress is already solved by the phase owner.
// This exercises its existing atomic thermal boundary, not private field edits.
fn transfer(
    thermal: &DepositionSimulation<'static>,
    id: u64,
    mass: f64,
    h: f64,
) -> DepositionSimulation<'static> {
    let mut next = thermal.clone();
    let expected = thermal.cells()[&CellId(id)].clone();
    let area = thermal.exposed_area_m2(CellId(id)).unwrap();
    next.apply_solid_exchanges(
        &[SolidCellExchange {
            expected,
            interval_exposed_area_m2: area,
            replacement: SolidSurfaceState {
                mass_kg: mass,
                specific_enthalpy_j_kg: h,
                area_m2: area,
            },
            pure_solid_enthalpy_range_j_kg: [1.0, 1000.0],
        }],
        || Ok(()),
    )
    .unwrap();
    next
}

#[test]
fn deposited_and_sublimated_mass_change_resistance_and_storage_together() {
    let thermal = simulation(cells(&[300.0, 400.0], 1.0, None), &SENSIBLE, vec![], 100);
    let fixed = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &UNIT_K)],
        &link_geometry(2, 1.0),
        &[],
        || Ok(()),
    )
    .unwrap();
    let observer = fixed
        .clone()
        .with_planar_phase_layers(&thermal, &[layer(0)], || Ok(()))
        .unwrap();
    let grown = transfer(&thermal, 0, 2.0, 300.0);
    assert!(grown.power_balance_observed(&fixed, || Ok(())).is_err());
    let rates = grown.power_balance_observed(&observer, || Ok(())).unwrap();
    // Thickness doubles: old half-cell resistance .5 becomes 1, while the
    // fixed neighbour retains .5. The heat capacities are now 2 and 1 J/K.
    let power = 100.0 / 1.5;
    near(rates.net_power_w[&CellId(0)], power, 1e-12);
    near(rates.net_power_w[&CellId(1)], -power, 1e-12);
    let (heated, prepared, audit) = grown
        .advance_observed(&observer, 0.1, config(0.1), || Ok(()))
        .unwrap();
    near(
        heated.temperature_k(CellId(0)).unwrap(),
        300.0 + 0.1 * power / 2.0,
        1e-12,
    );
    near(
        heated.temperature_k(CellId(1)).unwrap(),
        400.0 - 0.1 * power,
        1e-12,
    );
    near(heated.energy_j(), grown.energy_j(), 1e-12);
    near(audit.residual_j, 0.0, 1e-12);
    assert_eq!(prepared.last_planar_half_distance_m(CellId(0)), Some(1.0));
    assert_eq!(observer.last_planar_half_distance_m(CellId(0)), Some(0.5));
    assert_eq!(
        grown.cells()[&CellId(0)].geometry,
        thermal.cells()[&CellId(0)].geometry
    );
    assert_eq!(thermal.cells()[&CellId(0)].mass_kg, 1.0);

    let shrunk = transfer(&thermal, 0, 0.5, 300.0);
    let rates = shrunk.power_balance_observed(&observer, || Ok(())).unwrap();
    near(rates.net_power_w[&CellId(0)], 100.0 / 0.75, 1e-12);
    let both_layers = fixed
        .with_planar_phase_layers(&thermal, &[layer(0), layer(1)], || Ok(()))
        .unwrap();
    let both_grown = transfer(&grown, 1, 2.0, 400.0);
    let rates = both_grown
        .power_balance_observed(&both_layers, || Ok(()))
        .unwrap();
    near(rates.net_power_w[&CellId(0)], 50.0, 1e-12);
    // Retrying the exact original state retains its original heat coefficients.
    let rates = thermal
        .power_balance_observed(&observer, || Ok(()))
        .unwrap();
    near(rates.net_power_w[&CellId(0)], 100.0, 1e-12);
}

#[test]
fn layer_distance_reaches_reservoir_and_ambient_with_one_material_query() {
    let count = Cell::new(0);
    let conductivity = LinearK(&count);
    let mut cell = cells(&[300.0], 1.0, Some(500.0));
    cell[0].boundary.ambient_k = 450.0;
    let thermal = simulation(cell, &SENSIBLE, vec![], 100);
    let observer = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &conductivity)],
        &[],
        &[face(1.0, 0.25)],
        || Ok(()),
    )
    .unwrap()
    .with_planar_phase_layers(&thermal, &[layer(0)], || Ok(()))
    .unwrap()
    .with_ambient_faces(
        &thermal,
        &[ConductiveAmbient {
            cell_id: CellId(0),
            half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.5,
        }],
        || Ok(()),
    )
    .unwrap();
    let grown = transfer(&thermal, 0, 2.0, 350.0);
    count.set(0);
    let rates = grown.power_balance_observed(&observer, || Ok(())).unwrap();
    assert_eq!(
        count.get(),
        1,
        "growth and multiple faces must not duplicate conductivity queries"
    );
    let expected = (500.0 - 350.0) / (1.0 / 3.5 + 0.25) + (450.0 - 350.0) / (1.0 / 3.5 + 0.5);
    near(rates.net_power_w[&CellId(0)], expected, 1e-12);
    near(rates.external_power_w, expected, 1e-12);
}

#[test]
fn planar_reference_and_per_side_contact_measures_are_admitted_explicitly() {
    let thermal = simulation(cells(&[300.0, 400.0], 1.0, None), &SENSIBLE, vec![], 100);
    let build = || {
        ConductionObserver::new(
            &thermal,
            &[(MaterialId(1), &UNIT_K)],
            &link_geometry(2, 1.0),
            &[],
            || Ok(()),
        )
        .unwrap()
    };
    for invalid in [
        PlanarPhaseLayer {
            density_kg_m3: 2.0,
            ..layer(0)
        },
        PlanarPhaseLayer {
            footprint_area_m2: 0.5,
            ..layer(0)
        },
        PlanarPhaseLayer {
            footprint_area_m2: f64::NAN,
            ..layer(0)
        },
        layer(99),
    ] {
        assert!(
            build()
                .with_planar_phase_layers(&thermal, &[invalid], || Ok(()))
                .is_err()
        );
    }
    assert!(
        build()
            .with_planar_phase_layers(&thermal, &[layer(0), layer(0)], || Ok(()))
            .is_err()
    );
    let wrong_distance = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &UNIT_K)],
        &link_geometry(2, 0.5),
        &[],
        || Ok(()),
    )
    .unwrap();
    assert!(
        wrong_distance
            .with_planar_phase_layers(&thermal, &[layer(0)], || Ok(()))
            .is_err()
    );

    let patch = thermal.contact_patch(1).unwrap();
    let asymmetric = thermal
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: 1,
                expected: Some(patch),
                replacement: Some(ContactPatch {
                    covered_area_m2: [1.5, 0.5],
                    ..patch
                }),
            }],
            || Ok(()),
        )
        .unwrap();
    let observer = ConductionObserver::new(
        &asymmetric,
        &[(MaterialId(1), &UNIT_K)],
        &link_geometry(2, 1.0),
        &[],
        || Ok(()),
    )
    .unwrap();
    assert!(
        observer
            .with_planar_phase_layers(&asymmetric, &[layer(0)], || Ok(()))
            .is_err()
    );
    let observer = build()
        .with_planar_phase_layers(&thermal, &[layer(0)], || Ok(()))
        .unwrap();
    assert!(
        asymmetric
            .power_balance_observed(&observer, || Ok(()))
            .is_err()
    );
    let partial = thermal
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: 1,
                expected: Some(patch),
                replacement: Some(ContactPatch {
                    covered_area_m2: [0.5, 1.0],
                    ..patch
                }),
            }],
            || Ok(()),
        )
        .unwrap();
    let observer = ConductionObserver::new(
        &partial,
        &[(MaterialId(1), &UNIT_K)],
        &link_geometry(2, 1.0),
        &[],
        || Ok(()),
    )
    .unwrap()
    .with_planar_phase_layers(&partial, &[layer(0)], || Ok(()))
    .unwrap();
    let grown = transfer(&partial, 0, 2.0, 300.0);
    assert!(grown.power_balance_observed(&observer, || Ok(())).is_ok());
    // A contact nominally larger than the selected side's footprint is also
    // outside this explicitly planar normal-conduction model.
    assert!(layer(0).half_distance(0.0).is_err());
}

#[test]
fn changed_unselected_mass_cancellation_and_fixed_implicit_paths_refuse() {
    let thermal = simulation(cells(&[300.0, 400.0], 1.0, None), &SENSIBLE, vec![], 100);
    let observer = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &UNIT_K)],
        &link_geometry(2, 1.0),
        &[],
        || Ok(()),
    )
    .unwrap()
    .with_planar_phase_layers(&thermal, &[layer(0)], || Ok(()))
    .unwrap();
    let unsupported = transfer(&thermal, 1, 2.0, 400.0);
    assert!(
        unsupported
            .power_balance_observed(&observer, || Ok(()))
            .is_err()
    );
    let grown = transfer(&thermal, 0, 2.0, 300.0);
    let original = grown.cells().clone();
    let mut calls = 0;
    assert!(
        grown
            .advance_observed(&observer, 0.1, config(0.05), || {
                calls += 1;
                if calls > 10 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(grown.cells(), &original);
    assert_eq!(grown.time_s(), 0.0);
    assert_eq!(observer.last_planar_half_distance_m(CellId(0)), Some(0.5));
    assert!(matches!(
        implicit::ImplicitConductionChain::new(&thermal, &observer, &[], || Ok(())),
        Err(implicit::ImplicitConductionError::Invalid(
            "mass-evolving planar phase layers require the explicit conduction observer"
        ))
    ));
}
