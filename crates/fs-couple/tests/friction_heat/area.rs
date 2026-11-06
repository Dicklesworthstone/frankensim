//! Changing contact coverage must enter the original conduction AND ambient
//! heat rates. These two-cell references solve in temperature coordinates;
//! the production residual still receives cell-energy increments in joules.
use super::*;
use fs_thermal::enthalpy::backward_euler::BackwardEulerStep;

struct Reference {
    energy_j: [f64; 2],
    temperature_k: [f64; 2],
    conductances_w_k: [f64; 2],
    body_heat_j: [f64; 2],
    external_heat_j: f64,
}

fn reference(
    old: &DepositionSimulation<'_>,
    route: &FrictionHeatRoute,
    area: f64,
    k: [f64; 2],
    dt: f64,
    interface_k: f64,
) -> Reference {
    let mut result = Reference {
        energy_j: [0.0; 2],
        temperature_k: [0.0; 2],
        conductances_w_k: [0.0; 2],
        body_heat_j: [0.0; 2],
        external_heat_j: 0.0,
    };
    for (side, id) in route.cells().into_iter().enumerate() {
        let cell = &old.cells()[&id];
        // SENSIBLE has exactly one joule per kilogram-kelvin. Each body's
        // full surface is unchanged, while its actual covered part changes.
        let capacity = cell.mass_kg;
        let theta_old = old.temperature_k(id).unwrap() - 300.0;
        let g = k[side] * area / route.distances_m()[side];
        let ambient = cell.boundary.htc_w_m2_k * (cell.geometry.surface_area_m2 - area);
        let theta = (capacity * theta_old + dt * g * (interface_k - 300.0))
            / (capacity + dt * (g + ambient));
        result.energy_j[side] = capacity * (theta - theta_old);
        result.temperature_k[side] = 300.0 + theta;
        result.conductances_w_k[side] = g;
        result.body_heat_j[side] = dt * g * (interface_k - 300.0 - theta);
        result.external_heat_j -= dt * ambient * theta;
    }
    result
}

#[test]
fn changing_area_uses_actual_conduction_exposure_and_entropy_then_rebinds_the_route() {
    let mut thermal = network(1, [1.0, 3.0], [2.0, 9.0], [2.0; 2], [1.5, 2.5], [8.0, 3.0]);
    let mut route = FrictionHeatRoute::from_link(&thermal, CONTACT, [2.0, 9.0]).unwrap();
    let dt = 0.2;
    let interface_k = 330.0;
    // Shrink, then grow from the newly accepted coverage. The endpoint source
    // conductivity also changes, so the successor cannot retain a stale link.
    for (area, k) in [(0.25, [3.0, 12.0]), (0.75, [5.0, 6.0])] {
        let original_cells = thermal.cells().clone();
        let original_patch = thermal.contact_patch(CONTACT).unwrap();
        let original_time = thermal.time_s();
        let expected = reference(&thermal, &route, area, k, dt, interface_k);
        let dissipated = expected.body_heat_j.iter().sum::<f64>();
        assert!(dissipated > 0.0);
        let step = FrictionHeatStep::new(&thermal, original_time + dt, &route, || Ok(())).unwrap();
        let trial = step
            .evaluate_with_area_and_entropy(
                &expected.energy_j,
                interface_k,
                |endpoint, t| {
                    assert_eq!(endpoint.contact_patch(CONTACT), Some(original_patch));
                    close(t, interface_k, 0.0);
                    Ok::<_, DepositionError>((area, Box::new((area, t))))
                },
                |endpoint, t, staged| {
                    assert_eq!(*staged, (area, t));
                    let patch = endpoint.contact_patch(CONTACT).unwrap();
                    assert_eq!(patch.covered_area_m2, [area; 2]);
                    assert_eq!(patch.link.area_m2, area);
                    // The massless interface owns this exchange during the
                    // residual; a second direct conductive path must be off.
                    assert_eq!(patch.link.conductance_w_k, 0.0);
                    for (side, id) in route.cells().into_iter().enumerate() {
                        let cell = &endpoint.cells()[&id];
                        close(
                            endpoint.exposed_area_m2(id).unwrap(),
                            cell.geometry.surface_area_m2 - area,
                            0.0,
                        );
                        close(
                            endpoint.temperature_k(id).unwrap(),
                            expected.temperature_k[side],
                            1e-12,
                        );
                        assert_eq!(cell.mass_kg, original_cells[&id].mass_kg);
                        assert_eq!(cell.geometry, original_cells[&id].geometry);
                    }
                    Ok::<_, DepositionError>(source(dissipated, k, staged))
                },
                || Ok(()),
            )
            .unwrap();
        assert!(trial.residual_j().iter().all(|r| r.abs() < 1e-11));
        assert_eq!(trial.route().contact_area_m2(), area);
        assert_eq!(**trial.payload(), (area, interface_k));
        let (next, staged, report, successor) = trial.accept(tolerances(), || Ok(())).unwrap();
        assert_eq!(*staged, (area, interface_k));
        for side in 0..2 {
            close(
                report.conductances_w_k[side],
                expected.conductances_w_k[side],
                1e-13,
            );
            close(report.body_heat_j[side], expected.body_heat_j[side], 1e-12);
        }
        close(report.external_heat_j, expected.external_heat_j, 1e-12);
        close(
            report.thermal.stored_energy_change_j,
            expected.energy_j.iter().sum(),
            1e-12,
        );
        close(report.total_energy_residual_j, 0.0, 1e-11);
        // The old footprint produces a measurably different ambient account;
        // endpoint cooling cannot be evaluated before installing the patch.
        let stale_external = route
            .cells()
            .into_iter()
            .enumerate()
            .map(|(side, id)| {
                -dt * original_cells[&id].boundary.htc_w_m2_k
                    * thermal.exposed_area_m2(id).unwrap()
                    * (expected.temperature_k[side] - 300.0)
            })
            .sum::<f64>();
        assert!((stale_external - expected.external_heat_j).abs() > 0.1);
        let entropy = report.entropy.unwrap();
        let storage = route
            .cells()
            .into_iter()
            .enumerate()
            .map(|(side, id)| {
                original_cells[&id].mass_kg
                    * (expected.temperature_k[side] / thermal.temperature_k(id).unwrap()).ln()
            })
            .sum::<f64>();
        close(entropy.thermal_storage_change_j_k, storage, 1e-12);
        close(
            entropy.reservoir_entropy_input_j_k,
            expected.external_heat_j / 300.0,
            1e-12,
        );
        close(
            entropy.friction_generation_j_k,
            dissipated / interface_k,
            1e-12,
        );
        let ambient_generation = route
            .cells()
            .into_iter()
            .enumerate()
            .map(|(side, id)| {
                let body_k = expected.temperature_k[side];
                let heat = -dt
                    * original_cells[&id].boundary.htc_w_m2_k
                    * (original_cells[&id].geometry.surface_area_m2 - area)
                    * (body_k - 300.0);
                heat * (1.0 / body_k - 1.0 / 300.0)
            })
            .sum::<f64>();
        assert!(ambient_generation > 0.0);
        close(
            entropy.network_conduction_generation_j_k,
            ambient_generation,
            1e-12,
        );
        assert!(!entropy.has_unclassified_sources);
        let patch = next.contact_patch(CONTACT).unwrap();
        assert_eq!(patch.covered_area_m2, [area; 2]);
        close(
            patch.link.conductance_w_k,
            series_conductance_w_k(area, route.distances_m(), k).unwrap(),
            1e-13,
        );
        assert_eq!(successor.contact_area_m2(), area);
        assert_eq!(successor.cells(), route.cells());
        assert_eq!(successor.distances_m(), route.distances_m());
        assert!(FrictionHeatRoute::from_link(&next, CONTACT, k).is_ok());
        assert!(FrictionHeatStep::new(&next, next.time_s() + dt, &successor, || Ok(())).is_ok());
        assert!(FrictionHeatStep::new(&next, next.time_s() + dt, &route, || Ok(())).is_err());
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(thermal.contact_patch(CONTACT), Some(original_patch));
        assert_eq!(thermal.time_s(), original_time);
        thermal = next;
        route = successor;
    }
    close(thermal.time_s(), 0.4, 0.0);
}

#[test]
fn endpoint_area_change_preserves_exact_target_event_ordering() {
    let seed = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let mut events: Vec<_> = seed
        .cells()
        .values()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(seed.links()[&CONTACT]),
    });
    events.push(TimedEvent {
        time_s: 0.2,
        event: Event::Separation(CONTACT),
    });
    let config = StepConfig {
        max_dt_s: 0.1,
        energy_tolerance_j: 1e-12,
    };
    let old = DepositionSimulation::new(events, &[(MaterialId(1), &SENSIBLE)], seed.budget())
        .unwrap()
        .advance_to(0.0, config, || Ok(()))
        .unwrap()
        .0;
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [1.0; 2]).unwrap();
    let (left, (), report, successor) = FrictionHeatStep::new(&old, 0.2, &route, || Ok(()))
        .unwrap()
        .evaluate_with_area(
            &[0.0; 2],
            300.0,
            |_, _| Ok::<_, DepositionError>((0.5, ())),
            |endpoint, _, ()| {
                assert_eq!(
                    endpoint.contact_patch(CONTACT).unwrap().covered_area_m2,
                    [0.5; 2]
                );
                Ok::<_, DepositionError>(source(0.0, [1.0; 2], ()))
            },
            || Ok(()),
        )
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    assert_eq!(left.time_s(), 0.2);
    assert_eq!(
        left.contact_patch(CONTACT).unwrap().covered_area_m2,
        [0.5; 2]
    );
    assert_eq!(report.thermal.events, 0);
    assert!(matches!(
        FrictionHeatStep::new(&left, 0.3, &successor, || Ok(())),
        Err(DepositionError::InvalidInput(_))
    ));
    // Processing the due separation changes topology without advancing time
    // or replaying the accepted interval's heat/storage.
    let (right, audit) = left.advance_to(0.2, config, || Ok(())).unwrap();
    assert_eq!(audit.events, 1);
    assert_eq!(audit.substeps, 0);
    assert_eq!(right.cells(), left.cells());
    assert_eq!(right.time_s(), 0.2);
    assert_eq!(right.contact_patch(CONTACT), None);
    assert!(FrictionHeatStep::new(&right, 0.3, &successor, || Ok(())).is_err());
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(
        old.contact_patch(CONTACT).unwrap().covered_area_m2,
        [1.0; 2]
    );
}

#[test]
fn refused_or_cancelled_endpoint_patch_never_publishes_partial_coverage() {
    let old = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let cells = old.cells().clone();
    let patch = old.contact_patch(CONTACT).unwrap();
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [1.0; 2]).unwrap();
    let step = FrictionHeatStep::new(&old, 0.2, &route, || Ok(())).unwrap();
    let prepare = |_: &DepositionSimulation<'_>, _| Ok::<_, DepositionError>((0.5, ()));
    assert!(matches!(
        step.evaluate_with_area(
            &[0.0; 2],
            300.0,
            prepare,
            |endpoint, _, ()| {
                assert_eq!(
                    endpoint.contact_patch(CONTACT).unwrap().covered_area_m2,
                    [0.5; 2]
                );
                Err::<FrictionHeatInput<()>, _>(DepositionError::MaterialDomain)
            },
            || Ok(())
        ),
        Err(DepositionError::MaterialDomain)
    ));
    let resolve =
        |_: &DepositionSimulation<'_>, _, ()| Ok::<_, DepositionError>(source(0.0, [1.0; 2], ()));
    let trial = step
        .evaluate_with_area(&[0.0; 2], 300.0, prepare, resolve, || Ok(()))
        .unwrap();
    assert_eq!(trial.route().contact_area_m2(), 0.5);
    assert!(matches!(
        trial.accept(tolerances(), || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    let unbalanced = step
        .evaluate_with_area(&[0.0; 2], 301.0, prepare, resolve, || Ok(()))
        .unwrap();
    assert!(matches!(
        unbalanced.accept(tolerances(), || Ok(())),
        Err(DepositionError::EnergyBalance(_))
    ));

    // An ignored error must poison the narrow thermal ingress: retaining the
    // first successful patch after a stale second update would publish a
    // partially completed coupled constitutive operation.
    let native = BackwardEulerStep::new(&old, 0.2, || Ok(())).unwrap();
    let half = ContactPatch::symmetric(ThermalLink {
        area_m2: 0.5,
        ..patch.link
    });
    let change = ContactPatchUpdate {
        id: CONTACT,
        expected: Some(patch),
        replacement: Some(half),
    };
    assert!(matches!(
        native.evaluate_with_contact_patches(
            &[0.0; 2],
            |endpoint| {
                endpoint.replace_contact_patches(&[change])?;
                assert_eq!(endpoint.thermal().contact_patch(CONTACT), Some(half));
                assert!(endpoint.replace_contact_patches(&[change]).is_err());
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || Ok(())
        ),
        Err(DepositionError::InvalidInput(_))
    ));
    assert_eq!(old.cells(), &cells);
    assert_eq!(old.contact_patch(CONTACT), Some(patch));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(route.contact_area_m2(), 1.0);
    // A fresh, fully admitted retry remains possible from the same original.
    let (accepted, (), _, rebound) = step
        .evaluate_with_area(&[0.0; 2], 300.0, prepare, resolve, || Ok(()))
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    assert_eq!(
        accepted.contact_patch(CONTACT).unwrap().covered_area_m2,
        [0.5; 2]
    );
    assert_eq!(rebound.contact_area_m2(), 0.5);
    assert_eq!(accepted.cells(), &cells);
}
