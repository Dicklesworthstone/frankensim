use super::*;
use fs_thermal::enthalpy::{
    DepositionBudget, DepositionSimulation, StepConfig, ThermalLaw, ThermalStorageBasis,
    conduction::{ConductionObserver, ConductivityLaw},
};

fn boundary() -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}
fn cell(id: u64, mass: f64) -> TransientFinCell {
    TransientFinCell {
        id: CellId(id),
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: 300.0,
        boundary: boundary(),
    }
}
fn input<'a>(
    edges: &'a [f64],
    cells: &'a [TransientFinCell],
    axial: &'a [u64],
    sides: &'a [u64],
    tip: Option<u64>,
) -> TransientFinInput<'a> {
    TransientFinInput {
        geometry: StraightFinGeometry::try_new(1.0, 1.0, 2.0).unwrap(),
        axial_edges_m: edges,
        inventory: cells,
        placement: FinPlacement {
            root_m: [1.0, 2.0, 3.0],
            direction: [0.0, 0.0, 1.0],
            body_id: 8,
            segment_id: 9,
        },
        base_cell: CellId(1),
        fluid_cell: CellId(2),
        root_link: 100,
        axial_links: axial,
        side_links: sides,
        tip_link: tip,
        base_half_distance_m: 0.5,
        root_contact_resistance_m2_k_w: 0.0,
    }
}
fn budget() -> TransientFinBudget {
    TransientFinBudget {
        maximum_cells: 8,
        maximum_links: 20,
    }
}

// Authored analytic test material, with one joule per kilogram per kelvin.
// The original native timeline below owns every actual heat/storage update.
struct UnitCapacity;
impl ThermalLaw for UnitCapacity {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::DeclaredSensibleStorage
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > 0.0 {
            Ok(h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
struct UnitConductivity;
impl ConductivityLaw for UnitConductivity {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        UnitCapacity.temperature_k(h)?;
        Ok(1.0)
    }
}

fn endpoint(id: CellId, temperature: f64) -> DepositedCell {
    DepositedCell {
        id,
        segment_id: id.0,
        body_id: id.0,
        material_id: MaterialId(1),
        mass_kg: 1.0,
        specific_h_j_kg: temperature,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, 0.0],
            end_m: [0.0, 1.0, 0.0],
            volume_m3: 1.0,
            surface_area_m2: 8.0,
        },
        boundary: boundary(),
    }
}
fn activate<'a>(
    topology: &TransientFinTopology,
    law: &'a UnitCapacity,
    base_k: f64,
    fluid_k: f64,
    film_h: f64,
) -> DepositionSimulation<'a> {
    let mut events = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(endpoint(topology.base_cell(), base_k)),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(endpoint(topology.fluid_cell(), fluid_k)),
        },
    ];
    let mut fin_events = topology.events_at(0.0, || Ok(())).unwrap();
    for event in &mut fin_events {
        if let Event::Contact(link) = &mut event.event {
            if topology.side_links().contains(&link.id) || topology.tip_link() == Some(link.id) {
                // Authored constant film for this independent analytic test.
                // Its actual native link supplies both heat signs and CFL.
                link.conductance_w_k = film_h * link.area_m2;
            }
        }
    }
    events.extend(fin_events);
    let simulation = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 10,
            max_links: 20,
            max_events: 30,
            max_substeps: 10_000,
        },
    )
    .unwrap();
    simulation
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0
}

#[test]
fn nonuniform_native_cells_retain_inventory_pose_and_exact_surface_ownership() {
    let inventory = [cell(42, 0.25), cell(17, 0.75)];
    let mut authored = input(
        &[0.0, 0.25, 1.0],
        &inventory,
        &[101],
        &[102, 103],
        Some(104),
    );
    authored.root_contact_resistance_m2_k_w = 0.02;
    let topology = TransientFinTopology::new(authored, budget(), || Ok(())).unwrap();
    assert_eq!(
        topology.cells().iter().map(|c| c.id).collect::<Vec<_>>(),
        [CellId(42), CellId(17)]
    );
    assert_eq!(topology.cells()[0].geometry.start_m, [1.0, 2.0, 3.0]);
    assert_eq!(topology.cells()[0].geometry.end_m, [1.0, 2.0, 3.25]);
    assert_eq!(topology.cells()[1].geometry.end_m, [1.0, 2.0, 4.0]);
    assert_eq!(topology.cells()[0].geometry.volume_m3, 0.25);
    assert_eq!(topology.cells()[1].geometry.volume_m3, 0.75);
    assert_eq!(topology.cells().iter().map(|c| c.mass_kg).sum::<f64>(), 1.0);
    assert!(
        topology
            .cells()
            .iter()
            .all(|c| c.body_id == 8 && c.segment_id == 9)
    );
    assert_eq!(
        topology.conductive_links()[0].half_distances_m,
        [0.5, 0.125]
    );
    assert_eq!(
        topology.conductive_links()[0].contact_resistance_m2_k_w,
        0.02
    );
    assert_eq!(
        topology.conductive_links()[1].half_distances_m,
        [0.125, 0.375]
    );
    assert_eq!(topology.links().len(), 5);
    assert!(
        topology
            .links()
            .iter()
            .all(|link| link.conductance_w_k == 0.0)
    );
    let law = UnitCapacity;
    let thermal = activate(&topology, &law, 300.0, 300.0, 0.0);
    assert_eq!(thermal.exposed_area_m2(CellId(42)).unwrap(), 0.0);
    assert_eq!(thermal.exposed_area_m2(CellId(17)).unwrap(), 0.0);
    assert_eq!(thermal.exposed_area_m2(topology.base_cell()).unwrap(), 7.0);
    assert_eq!(thermal.exposed_area_m2(topology.fluid_cell()).unwrap(), 5.0);
    assert_eq!(thermal.energy_j(), 900.0);
    assert_eq!(topology.events_at(2.0, || Ok(())).unwrap().len(), 7);
}

#[test]
fn native_three_store_startup_retains_fin_lag_and_refines_to_independent_modes() {
    let inventory = [cell(3, 1.0)];
    let topology = TransientFinTopology::new(
        input(&[0.0, 1.0], &inventory, &[], &[102], None),
        budget(),
        || Ok(()),
    )
    .unwrap();
    let law = UnitCapacity;
    let conductivity = UnitConductivity;
    // Root G=A/(0.5/k+0.5/k)=1, side G=h*P*L=1, all three capacities=1.
    let thermal = activate(&topology, &law, 400.0, 300.0, 0.5);
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &conductivity as &dyn ConductivityLaw)],
        topology.conductive_links(),
        &[],
        || Ok(()),
    )
    .unwrap();
    let step = StepConfig {
        max_dt_s: 0.01,
        energy_tolerance_j: 1e-9,
    };
    let (first, _, audit) = thermal
        .advance_observed(&conduction, 0.01, step, || Ok(()))
        .unwrap();
    assert_eq!(first.temperature_k(CellId(1)).unwrap(), 399.0);
    assert_eq!(first.temperature_k(CellId(3)).unwrap(), 301.0);
    assert_eq!(first.temperature_k(CellId(2)).unwrap(), 300.0);
    assert_eq!(first.energy_j(), thermal.energy_j());
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(3)).unwrap(), 300.0);

    // Independent eigenmodes of a three-capacity path: rates 0, G/C, 3G/C.
    let t: f64 = 0.7;
    let a = (-t).exp();
    let b = (-3.0 * t).exp();
    let expected = [
        300.0 + 100.0 * (1.0 / 3.0 + a / 2.0 + b / 6.0),
        300.0 + 100.0 * (1.0 / 3.0 - b / 3.0),
        300.0 + 100.0 * (1.0 / 3.0 - a / 2.0 + b / 6.0),
    ];
    let mut errors = Vec::new();
    for dt in [0.02, 0.01, 0.005] {
        let (next, _, report) = thermal
            .advance_observed(
                &conduction,
                t,
                StepConfig {
                    max_dt_s: dt,
                    energy_tolerance_j: 1e-8,
                },
                || Ok(()),
            )
            .unwrap();
        let error = [CellId(1), CellId(3), CellId(2)]
            .iter()
            .zip(expected)
            .map(|(&id, exact)| (next.temperature_k(id).unwrap() - exact).abs())
            .fold(0.0, f64::max);
        errors.push(error);
        assert!((next.energy_j() - thermal.energy_j()).abs() < 1e-8);
        assert_eq!(report.external_heat_j, 0.0);
        assert!(next.temperature_k(CellId(1)).unwrap() > next.temperature_k(CellId(3)).unwrap());
        assert!(next.temperature_k(CellId(3)).unwrap() > next.temperature_k(CellId(2)).unwrap());
    }
    assert!(errors[0] / errors[1] > 1.98, "{errors:?}");
    assert!(errors[1] / errors[2] > 1.98, "{errors:?}");
    assert!(errors[2] < 0.08, "{errors:?}");

    // Fin mass is independent of base mass and geometric volume.
    let double_mass = [cell(3, 2.0)];
    let heavier = TransientFinTopology::new(
        input(&[0.0, 1.0], &double_mass, &[], &[102], None),
        budget(),
        || Ok(()),
    )
    .unwrap();
    let heavy = activate(&heavier, &law, 400.0, 300.0, 0.5);
    let heavy_conduction = ConductionObserver::new(
        &heavy,
        &[(MaterialId(1), &conductivity as &dyn ConductivityLaw)],
        heavier.conductive_links(),
        &[],
        || Ok(()),
    )
    .unwrap();
    let (next, _, _) = heavy
        .advance_observed(&heavy_conduction, 0.01, step, || Ok(()))
        .unwrap();
    assert_eq!(next.temperature_k(CellId(3)).unwrap(), 300.5);
    assert_eq!(next.temperature_k(CellId(2)).unwrap(), 300.0);
    assert_eq!(next.energy_j(), heavy.energy_j());
}

#[test]
fn insulated_tip_is_uncovered_and_authored_heat_is_not_relabelled_as_convection() {
    let mut inventory = [cell(3, 1.0)];
    inventory[0].boundary.source_w = 5.0;
    let topology = TransientFinTopology::new(
        input(&[0.0, 1.0], &inventory, &[], &[102], None),
        budget(),
        || Ok(()),
    )
    .unwrap();
    let law = UnitCapacity;
    let thermal = activate(&topology, &law, 300.0, 300.0, 0.0);
    assert_eq!(thermal.exposed_area_m2(CellId(3)).unwrap(), 1.0);
    assert_eq!(topology.links().len(), 2);
    let (next, audit) = thermal
        .advance_to(
            0.2,
            StepConfig {
                max_dt_s: 0.2,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.temperature_k(CellId(3)).unwrap(), 301.0);
    assert_eq!(audit.external_heat_j, 1.0);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(next.temperature_k(CellId(2)).unwrap(), 300.0);
}

#[test]
fn malformed_geometry_aliases_occupied_boundaries_and_budgets_refuse() {
    let inventory = [cell(3, 0.5), cell(4, 0.5)];
    let edges = [0.0, 0.5, 1.0];
    for kind in 0..12 {
        let mut cells = inventory;
        match kind {
            7 => cells[1].id = cells[0].id,
            8 => cells[0].id = CellId(1),
            9 => cells[0].mass_kg = f64::NAN,
            10 => cells[0].boundary.htc_w_m2_k = 1.0,
            11 => cells[0].specific_h_j_kg = f64::INFINITY,
            _ => {}
        }
        let mut authored = input(&edges, &cells, &[101], &[102, 103], Some(104));
        match kind {
            0 => authored.base_cell = authored.fluid_cell,
            1 => authored.root_link = 102,
            2 => authored.tip_link = Some(101),
            3 => authored.base_half_distance_m = 0.0,
            4 => authored.root_contact_resistance_m2_k_w = -1.0,
            5 => authored.placement.direction = [0.0, 0.0, 2.0],
            6 => authored.placement.root_m = [0.0, 0.0, 1e30],
            7..=11 => {}
            _ => unreachable!(),
        }
        assert!(
            TransientFinTopology::new(authored, budget(), || Ok(())).is_err(),
            "case {kind}"
        );
    }
    for edges in [&[0.0, 0.0, 1.0][..], &[0.0, 1.2, 1.0], &[0.0, 0.5, 0.9]] {
        assert!(
            TransientFinTopology::new(
                input(edges, &inventory, &[101], &[102, 103], Some(104)),
                budget(),
                || Ok(()),
            )
            .is_err()
        );
    }
    for limited in [
        TransientFinBudget {
            maximum_cells: 1,
            maximum_links: 10,
        },
        TransientFinBudget {
            maximum_cells: 2,
            maximum_links: 4,
        },
    ] {
        assert!(matches!(
            TransientFinTopology::new(
                input(&edges, &inventory, &[101], &[102, 103], Some(104)),
                limited,
                || Ok(()),
            ),
            Err(DepositionError::Budget(_))
        ));
    }
}

#[test]
fn cancellation_leaves_inputs_and_native_candidate_unchanged_and_retry_succeeds() {
    let inventory = [cell(3, 1.0)];
    let mut calls = 0;
    assert_eq!(
        TransientFinTopology::new(
            input(&[0.0, 1.0], &inventory, &[], &[102], None),
            budget(),
            || {
                calls += 1;
                if calls == 5 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .err()
        .unwrap(),
        DepositionError::Cancelled
    );
    assert_eq!(inventory, [cell(3, 1.0)]);
    let topology = TransientFinTopology::new(
        input(&[0.0, 1.0], &inventory, &[], &[102], None),
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(
        topology
            .events_at(0.0, || Err(DepositionError::Cancelled))
            .unwrap_err(),
        DepositionError::Cancelled
    );
    let law = UnitCapacity;
    let thermal = activate(&topology, &law, 400.0, 300.0, 0.5);
    let saved = thermal.cells().clone();
    let mut calls = 0;
    assert_eq!(
        thermal
            .advance_to(
                1.0,
                StepConfig {
                    max_dt_s: 0.01,
                    energy_tolerance_j: 1e-8
                },
                || {
                    calls += 1;
                    if calls == 20 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                },
            )
            .err()
            .unwrap(),
        DepositionError::Cancelled
    );
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert!(
        thermal
            .advance_to(
                1.0,
                StepConfig {
                    max_dt_s: 0.01,
                    energy_tolerance_j: 1e-8
                },
                || Ok(()),
            )
            .is_ok()
    );
}
