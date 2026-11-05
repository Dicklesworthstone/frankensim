//! Native moving bulk, source-refreshed stationary joints and one thermal store.
use super::*;
use fs_couple::circuit_enthalpy::conduction::{
    heat::step_spatial_moving_conduction_candidate,
    spatial::motion::SpatialMovingConductionEnthalpyCandidate,
};
use fs_em::spatial::conduction::motion::SpatialConductionMotion;

fn moving_specimen(
    cards: &[MaterialCard; 2],
    interface: &InterfaceSystemCard,
    context: &QueryPoint,
    subdivisions: usize,
    speed: f64,
) -> (
    SpatialConduction,
    Vec<SpatialConductionRegion>,
    Vec<[f64; 3]>,
) {
    assert!(subdivisions > 0 && subdivisions % 2 == 0);
    let sources = cards.each_ref().map(|card| {
        resolve_conductor_state_with_basis(
            card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap()
    });
    let source = resolve_electrical_interface_state(
        interface,
        context,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let xs = (0..=subdivisions)
        .map(|i| i as f64 / subdivisions as f64)
        .collect::<Vec<_>>();
    let (mut positions, mut tets) = slab(&xs);
    let offset = positions.len() as u32;
    let (right, right_tets) = slab(&xs.iter().map(|x| 1.0 + x).collect::<Vec<_>>());
    positions.extend(right);
    tets.extend(right_tets.into_iter().map(|t| t.map(|v| v + offset)));
    let left = offset - 4;
    let contacts = [[0, 1, 3], [0, 2, 3]]
        .map(|face| SpatialConductionContact {
            first_vertices: face.map(|v| left + v),
            second_vertices: face.map(|v| offset + v),
            source: source.clone(),
        })
        .to_vec();
    let cells = 12 * subdivisions;
    let problem = SpatialConduction::new_with_contacts(
        positions,
        tets,
        (0..cells)
            .map(|i| sources[i / (6 * subdivisions)].clone())
            .collect(),
        contacts,
        limits(),
        &mut || false,
    )
    .unwrap();
    let regions = (0..cells)
        .map(|tetrahedron| SpatialConductionRegion {
            tetrahedron,
            cell: if tetrahedron < 6 * subdivisions {
                CellId(30)
            } else {
                CellId(10)
            },
        })
        .collect();
    // Motion occupies exactly the first half metre. Both original materials
    // adjacent to the interface at x=1 are stationary; no surface work is inferred.
    let velocities = (0..cells)
        .map(|i| {
            if i < 3 * subdivisions {
                [0.0, speed, 0.0]
            } else {
                [0.0; 3]
            }
        })
        .collect();
    (problem, regions, velocities)
}

fn assert_moving_balance<S>(
    next: &SpatialMovingConductionEnthalpyCandidate<'_, ResolvedConductorState, S>,
    original: &DepositionSimulation<'_>,
    r0: f64,
    conductance: f64,
    speed: f64,
) {
    let [a, b] = expected(r0, conductance);
    let surface_resistance = r0 * (1.0 + 0.003 * a + 0.007 * b);
    let bulk_resistance = 1.0 + 0.002 * a + 3.0;
    let (dt, current) = (0.5, 2.0);
    let heat = dt * current * current * (bulk_resistance + surface_resistance);
    let mechanical = dt * current * speed * 0.5;
    let report = next.report();
    let step = &report.electrical;
    near(next.thermal().temperature_k(CellId(30)).unwrap(), 300.0 + a);
    near(next.thermal().temperature_k(CellId(10)).unwrap(), 300.0 + b);
    near(
        step.field.voltage_v().last().copied().unwrap(),
        current * (bulk_resistance + surface_resistance) + 0.5 * speed,
    );
    near(
        step.field.cell_joule_power_w().iter().sum(),
        current * current * bulk_resistance,
    );
    near(
        step.field.contact_power_w(),
        current * current * surface_resistance,
    );
    near(step.heat_j, heat);
    near(step.mechanical_work_j, mechanical);
    near(step.terminal_work_j, heat + mechanical);
    near(step.net_supplied_work_j, heat);
    near(report.deposited_joule_j, heat);
    near(next.thermal().energy_j() - original.energy_j(), heat);
    near(report.total_energy_residual_j, 0.0);
    assert_eq!(step.contacts.len(), 2);
    for (contact, source) in step.contacts.iter().zip(next.problem().contacts()) {
        assert_eq!(contact.cells, [CellId(30), CellId(10)]);
        near(contact.heat_fraction_to_first, 0.25);
        near(contact.response.current_a, -1.0);
        near(contact.response.power_w, 2.0 * surface_resistance);
        near(source.source.resistance_ohm_m2(), surface_resistance);
        near(
            source.source.resolved().state_coordinate("Ta").unwrap(),
            300.0 + a,
        );
        near(
            source.source.resolved().state_coordinate("Tb").unwrap(),
            300.0 + b,
        );
    }
    if speed < -20.0 {
        assert!(step.terminal_work_j < 0.0);
        assert!(step.mechanical_work_j < step.terminal_work_j);
    }
}

#[test]
fn stationary_joints_and_moving_bulk_share_source_feedback_and_one_heat_deposit() {
    for (subdivisions, speed) in [(2, 3.0), (4, -30.0), (2, 0.0)] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let interface = interface(&cards, 0.5, 1000.0);
        let ip = interface_point(300.0, 1e5);
        let surfaces = surface_bindings(&cards, &interface, &ip);
        let (problem, regions, velocities) =
            moving_specimen(&cards, &interface, &ip, subdivisions, speed);
        let fields = vec![[0.0, 0.0, 1.0]; velocities.len()];
        let motion = SpatialConductionMotion {
            velocity_m_s: &velocities,
            flux_density_t: &fields,
        };
        let thermal = initial(&curves, 999.0, None);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let drive = boundary(&problem, 2.0, 2.0, None);
        let next = layout
            .step_moving_contact_candidate(
                &thermal,
                0.5,
                &problem,
                &drive,
                motion,
                &surfaces,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        assert_moving_balance(&next, &thermal, 0.5, 1.0, speed);
        // Motion changes terminal/mechanical work, while fixed current keeps
        // both original bulk I²R and contact heat at the same thermal endpoint.
        if speed == 0.0 {
            let stationary = layout
                .step_contact_candidate(
                    &thermal,
                    0.5,
                    &problem,
                    &drive,
                    &surfaces,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            near(stationary.thermal().energy_j(), next.thermal().energy_j());
            near(
                stationary.report().electrical.heat_j,
                next.report().electrical.heat_j,
            );
        }
        assert!(
            layout
                .step_moving_candidate(
                    &thermal,
                    0.5,
                    &problem,
                    &drive,
                    motion,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .is_err()
        );
        assert!(matches!(
            layout.step_moving_contact_candidate(
                &thermal,
                0.5,
                &problem,
                &drive,
                motion,
                &surfaces,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || true,
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        let moving_interfaces = vec![[0.0, 1.0, 0.0]; velocities.len()];
        assert!(
            layout
                .step_moving_contact_candidate(
                    &thermal,
                    0.5,
                    &problem,
                    &drive,
                    SpatialConductionMotion {
                        velocity_m_s: &moving_interfaces,
                        flux_density_t: &fields
                    },
                    &surfaces,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
        near(problem.contacts()[0].source.resistance_ohm_m2(), 0.5);
    }
}

#[test]
fn complete_heat_assembly_keeps_thermal_resistance_independent_of_contact_joule_loss() {
    for (pressure, speed) in [(1e5, 3.0), (2e5, -30.0)] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let interface = interface(&cards, 0.5, 1000.0);
        let ip = interface_point(300.0, pressure);
        let surfaces = surface_bindings(&cards, &interface, &ip);
        let (problem, regions, velocities) = moving_specimen(&cards, &interface, &ip, 2, speed);
        let fields = vec![[0.0, 0.0, 1.0]; velocities.len()];
        let thermal = initial(&curves, 999.0, None);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &LINKS,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let thermal_contacts = [MaterialContactLink {
            link: 7,
            point: surfaces[0].point,
            sides: [0, 1].map(|i| ContactCell {
                cell: cells[i].cell,
                material_id: cells[i].material_id,
                curve: curves[i].curve(),
            }),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }];
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &thermal_contacts,
            radiation: None,
        };
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 73,
                    kernel_id: 52,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat.prepare(&cx, &thermal).unwrap();
            let drive = boundary(&problem, 2.0, 2.0, None);
            let next = step_spatial_moving_conduction_candidate(
                &prepared,
                &cx,
                &thermal,
                0.5,
                &problem,
                &drive,
                SpatialConductionMotion {
                    velocity_m_s: &velocities,
                    flux_density_t: &fields,
                },
                &surfaces,
                &regions,
                1e-12,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
            )
            .unwrap();
            let thermal_resistance = 2.0 - 1.8 * (pressure - 1e5) / 1e5;
            assert_moving_balance(
                &next,
                &thermal,
                0.5,
                1.0 / (1.0 + thermal_resistance),
                speed,
            );
            near(
                next.report().source.thermal_contacts[0]
                    .interface
                    .value_m2_k_per_w(),
                thermal_resistance,
            );
            near(next.report().non_electrical_heat_j, 0.0);
        });
        assert_eq!(thermal.time_s(), 0.0);
    }
}
