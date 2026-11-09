//! Actual coupled wall states, independent triangle loss integrals, and native
//! thermal events. Authored coefficients exercise equations, not calibration.

use super::*;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::MechanicalBoundary;
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};

const WALL_IDS: [CellId; 2] = [CellId(40), CellId(12)];
const AREAS: [f64; 2] = [0.5, 1.0];
const SURFACE_MASSES: [f64; 2] = [0.6, 0.8];

fn boundaries() -> [MechanicalBoundary; 2] {
    core::array::from_fn(|wall| MechanicalBoundary {
        faces: vec![[0, 1, [3, 4][wall]]],
        wall: WallPin {
            surface_density: SURFACE_MASSES[wall],
            stiffness_per_area: [4.0, 6.0][wall],
            resistance: [0.4, 0.7][wall],
        },
        relaxation: vec![
            AcousticWallRelaxation {
                stiffness_per_area: 2.5,
                relaxation_time_s: 0.3,
            },
            AcousticWallRelaxation {
                stiffness_per_area: 0.7,
                relaxation_time_s: 0.08,
            },
        ],
    })
}

fn mechanical_prepared<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> PreparedEquivalentFluid<'f, 'c> {
    PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
        geometry,
        FREQUENCY_HZ,
        &[],
        &boundaries(),
        &[],
        EquivalentFluidControls::default(),
        &mut |element, _| Ok::<_, HeterogeneousError>(coefficients(element)),
        &mut || false,
    )
    .unwrap()
}

fn wall_routes(law: &dyn ThermalLaw) -> [InteriorWallHeatCell<'_>; 2] {
    core::array::from_fn(|wall| InteriorWallHeatCell {
        cell: WALL_IDS[wall],
        caloric: law,
        acoustic_area_m2: AREAS[wall],
    })
}

fn all_events() -> Vec<TimedEvent> {
    let mut events = schedule([0.2, -0.05]);
    for wall in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: WALL_IDS[wall],
                segment_id: 200 + wall as u64,
                body_id: 9,
                material_id: MaterialId(6),
                mass_kg: AREAS[wall] * SURFACE_MASSES[wall],
                specific_h_j_kg: 10.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, 0.0],
                    end_m: [0.0, 0.01, 0.0],
                    volume_m3: AREAS[wall] * 0.01,
                    // Includes both thermal faces and edges; deliberately NOT
                    // the one-sided acoustic contact area.
                    surface_area_m2: 2.0 * AREAS[wall] + 0.1,
                },
                boundary: boundary([0.1, -0.02][wall]),
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.375,
        event: Event::Boundary(WALL_IDS[0], boundary(0.3)),
    });
    events
}

// Independent exact integral of a linear phasor on a flat triangle:
// integral |v|² = A (sum |v_i|² + |sum v_i|²) / 12. The original
// Maxwell law gives e_j = sqrt(k_j)*v / (1/tau_j+i*omega), so neither
// the retained e heat matrix nor a pressure-work residual supplies this oracle.
fn analytic_wall_powers(states: &[MechanicalWallPhasors]) -> [[f64; 2]; 2] {
    let omega = core::f64::consts::TAU * FREQUENCY_HZ;
    let boundaries = boundaries();
    core::array::from_fn(|wall| {
        let velocity = &states[wall].normal_velocity_m_s;
        assert_eq!(velocity.len(), 3);
        let sum = velocity.iter().copied().fold(C64::ZERO, |a, b| a + b);
        let integral = AREAS[wall]
            * (velocity.iter().map(|v| v.norm_sq()).sum::<f64>() + sum.norm_sq())
            / 12.0;
        let relaxation = boundaries[wall]
            .relaxation
            .iter()
            .map(|arm| {
                let rate = 1.0 / arm.relaxation_time_s;
                arm.stiffness_per_area * rate / (rate * rate + omega * omega)
            })
            .sum::<f64>();
        [
            0.5 * boundaries[wall].wall.resistance * integral,
            0.5 * relaxation * integral,
        ]
    })
}

#[test]
fn original_volume_and_wall_heat_enter_distinct_cells_with_separate_native_receipts() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = mechanical_prepared(&geometry);
    let field = prepared
        .solve(&[([0, 2, 3], [C64::new(0.2, -0.05); 3])], &mut || false)
        .unwrap();
    let expected = analytic_wall_powers(&field.walls);
    let law = caloric();
    let thermal = initial(all_events(), &law);
    let heating = InteriorAcousticHeating::try_new_with_mechanical_walls(
        &prepared,
        &field.pressure_pa,
        &field.walls,
        &thermal,
        &routes(&law),
        &wall_routes(&law),
        4,
        || Ok(()),
    )
    .unwrap();
    assert!(core::ptr::eq(
        heating.wall_phasors().unwrap(),
        field.walls.as_slice()
    ));
    assert_eq!(heating.cells().len(), 2);
    assert_eq!(heating.walls().len(), 2);
    for (wall, heat) in heating.walls().enumerate() {
        assert_eq!(heat.cell, WALL_IDS[wall]);
        assert_eq!(heat.wall, wall);
        close(heat.area_m2, AREAS[wall], 2e-14);
        close(heat.mass_kg, AREAS[wall] * SURFACE_MASSES[wall], 2e-14);
        assert!(expected[wall].iter().all(|&power| power > 0.0));
        close(heat.resistance_dissipation_w, expected[wall][0], 2e-13);
        close(heat.relaxation_dissipation_w, expected[wall][1], 2e-13);
    }
    let step = heating
        .advance_to(&thermal, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let bulk = field
        .element_power
        .iter()
        .map(|power| power.total_dissipation_w())
        .sum::<f64>();
    let resistance = expected.iter().map(|power| power[0]).sum::<f64>();
    let relaxation = expected.iter().map(|power| power[1]).sum::<f64>();
    close(
        step.acoustic.wall_resistance_dissipation_j,
        resistance,
        2e-13,
    );
    close(
        step.acoustic.wall_relaxation_dissipation_j,
        relaxation,
        2e-13,
    );
    close(
        step.acoustic.deposited_j,
        bulk + resistance + relaxation,
        2e-13,
    );
    for wall in 0..2 {
        let acoustic_j = expected[wall].iter().sum::<f64>();
        let authored_j = [0.1 * 0.375 + 0.3 * 0.625, -0.02][wall];
        let original = &thermal.cells()[&WALL_IDS[wall]];
        let final_cell = &step.thermal.cells()[&WALL_IDS[wall]];
        assert_eq!(step.wall_energy_j[wall].0, WALL_IDS[wall]);
        close(step.wall_energy_j[wall].1, acoustic_j, 2e-13);
        close(
            (final_cell.specific_h_j_kg - original.specific_h_j_kg) * original.mass_kg,
            acoustic_j + authored_j,
            2e-12,
        );
        assert_eq!(final_cell.geometry, original.geometry);
        assert_eq!(final_cell.mass_kg, original.mass_kg);
    }
    close(
        step.heat_report.external_heat_j,
        step.acoustic.deposited_j + 0.2 - 0.05 + 0.1 * 0.375 + 0.3 * 0.625 - 0.02,
        2e-12,
    );
    close(step.heat_report.residual_j, 0.0, 2e-12);
    assert_eq!(step.heat_report.events, 1);
    assert_eq!(thermal.time_s(), 0.0);

    // Existing bulk-only API remains explicit: it never silently sends wall
    // work into the gas or deposits wall heat without a physical receiver.
    let bulk_only = InteriorAcousticHeating::try_new(
        &prepared,
        &field.pressure_pa,
        &thermal,
        &routes(&law),
        2,
        || Ok(()),
    )
    .unwrap();
    assert!(bulk_only.wall_phasors().is_none());
    assert_eq!(bulk_only.walls().len(), 0);
    let bulk_step = bulk_only
        .advance_to(&thermal, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert!(bulk_step.wall_energy_j.is_empty());
    assert_eq!(bulk_step.acoustic.wall_resistance_dissipation_j, 0.0);
    close(bulk_step.acoustic.deposited_j, bulk, 2e-13);
}

#[test]
fn wall_heating_composes_with_original_contacts_and_preserves_split_observers() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = mechanical_prepared(&geometry);
    let field = prepared
        .solve(&[([0, 2, 3], [C64::from_re(0.2); 3])], &mut || false)
        .unwrap();
    let law = caloric();
    let mut events = all_events();
    events.insert(
        events.len() - 1,
        TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 5,
                a: IDS[0],
                b: WALL_IDS[0],
                area_m2: 0.25,
                conductance_w_k: 384.0,
            }),
        },
    );
    let thermal = initial(events.clone(), &law);
    let heating = InteriorAcousticHeating::try_new_with_mechanical_walls(
        &prepared,
        &field.pressure_pa,
        &field.walls,
        &thermal,
        &routes(&law),
        &wall_routes(&law),
        4,
        || Ok(()),
    )
    .unwrap();
    let observers = CoupledObservers::new((), heating.observer());
    let (middle, middle_observers, _) = thermal
        .advance_observed(&observers, 0.5, config(10.0), || Ok(()))
        .unwrap();
    let (end, accepted, _) = middle
        .advance_observed(&middle_observers, 1.0, config(10.0), || Ok(()))
        .unwrap();
    let (whole, whole_observers, _) = thermal
        .advance_observed(&observers, 1.0, config(10.0), || Ok(()))
        .unwrap();

    // Original native reference uses independent analytic wall watts and the
    // separately evaluated volume forms as authored thermal sources. Boundary
    // replacement events carry their same acoustic source exactly once.
    let wall_power = analytic_wall_powers(&field.walls);
    let extra: BTreeMap<_, _> = IDS
        .into_iter()
        .zip(&field.element_power)
        .map(|(id, power)| (id, power.total_dissipation_w()))
        .chain(
            WALL_IDS
                .into_iter()
                .zip(wall_power)
                .map(|(id, power)| (id, power[0] + power[1])),
        )
        .collect();
    for event in &mut events {
        match &mut event.event {
            Event::Birth(cell) => cell.boundary.source_w += extra[&cell.id],
            Event::Boundary(id, boundary) => boundary.source_w += extra[id],
            _ => {}
        }
    }
    let reference = initial(events, &law)
        .advance_to(1.0, config(10.0), || Ok(()))
        .unwrap()
        .0;
    for id in IDS.into_iter().chain(WALL_IDS) {
        close(
            end.cells()[&id].specific_h_j_kg,
            reference.cells()[&id].specific_h_j_kg,
            2e-11,
        );
        close(
            end.cells()[&id].specific_h_j_kg,
            whole.cells()[&id].specific_h_j_kg,
            2e-11,
        );
    }
    assert_eq!(end.links(), thermal.links());
    close(
        accepted.second.energy().unwrap().deposited_j,
        whole_observers.second.energy().unwrap().deposited_j,
        2e-13,
    );
    for (wall, &(_, energy)) in accepted.second.wall_energy_j().unwrap().iter().enumerate() {
        close(energy, wall_power[wall].iter().sum::<f64>(), 2e-13);
    }
    assert_eq!(observers.second.energy().unwrap().deposited_j, 0.0);

    let zero = prepared.solve(&[], &mut || false).unwrap();
    let zero_heat = InteriorAcousticHeating::try_new_with_mechanical_walls(
        &prepared,
        &zero.pressure_pa,
        &zero.walls,
        &thermal,
        &routes(&law),
        &wall_routes(&law),
        4,
        || Ok(()),
    )
    .unwrap();
    let zero_step = zero_heat
        .advance_to(&thermal, 1.0, config(10.0), || Ok(()))
        .unwrap();
    let native = thermal.advance_to(1.0, config(10.0), || Ok(())).unwrap();
    assert_eq!(zero_step.thermal.cells(), native.0.cells());
    assert_eq!(zero_step.heat_report, native.1);
    assert_eq!(zero_step.acoustic.deposited_j, 0.0);
}

#[test]
fn original_wall_mass_area_caloric_and_complete_state_are_required() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = mechanical_prepared(&geometry);
    let field = prepared
        .solve(&[([0, 2, 3], [C64::ONE; 3])], &mut || false)
        .unwrap();
    let law = caloric();
    let twin = law;
    let thermal = initial(all_events(), &law);
    let mapping = wall_routes(&law);
    let make = |wall_phasors,
                thermal: &DepositionSimulation<'_>,
                wall_destinations: &[InteriorWallHeatCell<'_>],
                cap| {
        InteriorAcousticHeating::try_new_with_mechanical_walls(
            &prepared,
            &field.pressure_pa,
            wall_phasors,
            thermal,
            &routes(&law),
            wall_destinations,
            cap,
            || Ok(()),
        )
        .map(|_| ())
    };
    assert!(make(&field.walls, &thermal, &mapping, 3).is_err());
    assert!(make(&field.walls, &thermal, &mapping[..1], 4).is_err());
    assert!(make(&field.walls[..1], &thermal, &mapping, 4).is_err());
    for bad in [
        [mapping[1], mapping[0]],
        [mapping[0], mapping[0]],
        [
            InteriorWallHeatCell {
                cell: IDS[0],
                ..mapping[0]
            },
            mapping[1],
        ],
        [
            InteriorWallHeatCell {
                acoustic_area_m2: 0.51,
                ..mapping[0]
            },
            mapping[1],
        ],
        [
            InteriorWallHeatCell {
                acoustic_area_m2: f64::NAN,
                ..mapping[0]
            },
            mapping[1],
        ],
        wall_routes(&twin),
    ] {
        assert!(matches!(
            make(&field.walls, &thermal, &bad, 4),
            Err(InteriorHeatingError::Invalid(_))
        ));
    }
    let mut incomplete = field.walls.clone();
    incomplete[0].relaxation_extension.clear();
    assert!(matches!(
        make(&incomplete, &thermal, &mapping, 4),
        Err(InteriorHeatingError::Acoustic(_))
    ));

    let heating = InteriorAcousticHeating::try_new_with_mechanical_walls(
        &prepared,
        &field.pressure_pa,
        &field.walls,
        &thermal,
        &routes(&law),
        &mapping,
        4,
        || Ok(()),
    )
    .unwrap();
    for mutation in 0..6 {
        let mut changed = all_events();
        let Event::Birth(cell) = &mut changed[2].event else {
            panic!("wall birth")
        };
        match mutation {
            0 => cell.mass_kg *= 1.01,
            1 => cell.geometry.surface_area_m2 *= 1.01,
            2 => cell.geometry.volume_m3 *= 1.01,
            3 => cell.body_id += 1,
            4 => cell.segment_id += 1,
            _ => cell.material_id = MaterialId(5),
        }
        let changed = initial(changed, &law);
        assert!(matches!(
            heating.advance_to(&changed, 0.1, config(0.1), || Ok(())),
            Err(InteriorHeatingError::Invalid(_))
        ));
        if mutation == 0 {
            assert!(matches!(
                make(&field.walls, &changed, &mapping, 4),
                Err(InteriorHeatingError::Invalid(_))
            ));
        }
    }
}

#[test]
fn late_wall_heating_refusals_preserve_native_cells_events_and_all_four_accounts() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = mechanical_prepared(&geometry);
    let field = prepared
        .solve(&[([0, 2, 3], [C64::ONE; 3])], &mut || false)
        .unwrap();
    let law = caloric();
    let mut events = all_events();
    events.push(TimedEvent {
        time_s: 0.5,
        event: Event::Boundary(CellId(999), boundary(0.0)),
    });
    let thermal = initial(events, &law);
    let original = thermal.cells().clone();
    let heating = InteriorAcousticHeating::try_new_with_mechanical_walls(
        &prepared,
        &field.pressure_pa,
        &field.walls,
        &thermal,
        &routes(&law),
        &wall_routes(&law),
        4,
        || Ok(()),
    )
    .unwrap();
    let observer = heating.observer();
    let mut polls = 0;
    assert!(matches!(
        thermal.advance_observed(&observer, 0.5, config(0.1), || {
            polls += 1;
            Ok(())
        }),
        Err(InteriorHeatingError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    assert!(polls > 20);
    let mut calls = 0;
    assert!(matches!(
        thermal.advance_observed(&observer, 0.5, config(0.1), || {
            calls += 1;
            if calls >= polls * 3 / 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(InteriorHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
    let energy = observer.energy().unwrap();
    assert_eq!(energy.deposited_j, 0.0);
    assert_eq!(energy.wall_resistance_dissipation_j, 0.0);
    assert_eq!(energy.wall_relaxation_dissipation_j, 0.0);
    assert!(
        observer
            .wall_energy_j()
            .unwrap()
            .iter()
            .all(|&(_, value)| value == 0.0)
    );
}
