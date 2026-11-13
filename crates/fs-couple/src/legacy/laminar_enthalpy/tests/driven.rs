//! Original endpoint mechanical feedback and independently measured work.

use super::super::boussinesq::BoussinesqConfig;
use super::super::drive::*;
use super::*;
use fs_flux::OpenBoundaryPolicy;

struct Actuator {
    thermal_gain: f64,
    transverse_force: f64,
    pressure_gradient: [f64; 2],
}
struct Snapshot<'a> {
    original: &'a Actuator,
    temperature: f64,
    boundaries: Vec<LaminarDrivenBoundary>,
}
impl LaminarMechanicalDrive for Snapshot<'_> {
    fn boundaries(&self) -> &[LaminarDrivenBoundary] {
        &self.boundaries
    }
    fn force_density(&self, _: [f64; 2]) -> [f64; 2] {
        [
            self.original.thermal_gain * (self.temperature - 300.0),
            self.original.transverse_force,
        ]
    }
    fn velocity(&self, _: usize, _: [f64; 2]) -> [f64; 2] {
        [0.0; 2]
    }
    fn traction(&self, _: usize, point: [f64; 2], normal: [f64; 2]) -> [f64; 2] {
        let pressure = 17.0
            + self.original.pressure_gradient[0] * point[0]
            + self.original.pressure_gradient[1] * point[1];
        normal.map(|n| -pressure * n)
    }
}
fn snapshot<'a>(
    original: &'a Actuator,
    mesh: &TriMesh,
    thermal: &DepositionSimulation<'_>,
    cells: &[CellId],
) -> Snapshot<'a> {
    let temperature = cells
        .iter()
        .map(|&cell| thermal.temperature_k(cell).unwrap())
        .sum::<f64>()
        / cells.len() as f64;
    Snapshot {
        original,
        temperature,
        boundaries: (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| LaminarDrivenBoundary::Open {
                edge,
                backflow: OpenBoundaryPolicy::NaturalTraction,
            })
            .collect(),
    }
}
fn heating(
    mesh: &TriMesh,
    law: &Caloric,
    old: &DepositionSimulation<'_>,
    cells: &[CellId],
) -> LaminarHeatExchange<()> {
    let mut response = exchange(mesh, law.h(310.0));
    for &cell in cells {
        response
            .heat
            .bulk
            .sources_w
            .insert(cell, old.cells()[&cell].mass_kg * CP * 10.0 / 0.2);
    }
    response
}
fn verify<G>(
    system: &FluxSystem<'_>,
    cells: &[CellId],
    original: &Actuator,
    solved: &LaminarDrivenCandidate<'_, Snapshot<'_>, G, ()>,
    buoyancy_y: f64,
) {
    let force = [2.0, 0.3];
    let velocity = [
        0.2 * (force[0] - 0.4) / RHO,
        0.2 * (force[1] + buoyancy_y) / RHO,
    ];
    let (error, divergence) =
        system.velocity_error(&solved.flow().step.solution.discrete().x, &|_| velocity);
    assert!(error < 2e-8 && divergence < 2e-8, "{error}, {divergence}");
    for (cell, &id) in cells.iter().enumerate() {
        near(solved.thermal().temperature_k(id).unwrap(), 310.0, 2e-8);
        near(
            solved.flow().samples()[cell]
                .material()
                .thermal_state()
                .temperature_k(),
            310.0,
            2e-8,
        );
        near(
            solved.flow().step.solution.pressure_pa()[cell],
            17.0 + 0.4 * system.mesh.centroids[cell][0],
            2e-8,
        );
    }
    let receipt = &solved.report().exchange;
    assert_eq!(
        receipt.query,
        LaminarDriveQuery {
            time_s: 0.2,
            segment: 7
        }
    );
    assert!(core::ptr::eq(receipt.drive.original, original));
    near(receipt.drive.temperature, 310.0, 2e-8);
    let mechanical = solved.report().mechanical.unwrap();
    near(
        mechanical.body_work_j,
        0.2 * DEPTH * (force[0] * velocity[0] + force[1] * velocity[1]),
        2e-8,
    );
    near(
        mechanical.imposed_traction_work_j,
        -0.2 * DEPTH * 0.4 * velocity[0],
        2e-8,
    );
    near(mechanical.outward_kinetic_energy_j, 0.0, 2e-10);
    near(mechanical.boundary_upwind_loss_j, 0.0, 2e-10);
    assert_eq!(mechanical.velocity_boundary_work_j, Some(0.0));
    assert_eq!(mechanical.backflow_resistance_work_j, Some(0.0));
    assert!(mechanical.body_work_j > 0.0 && mechanical.imposed_traction_work_j < 0.0);
    near(solved.report().thermal.stored_energy_change_j, 24.0, 2e-7);
    near(solved.report().other_external_heat_j, 24.0, 2e-7);
    near(solved.report().external_stream_enthalpy_j, 0.0, 2e-7);
    near(solved.report().bulk_viscous_heat_j, 0.0, 1e-12);
    near(solved.report().thermal_balance_residual_j, 0.0, 1e-8);
    match solved.report().buoyancy.as_ref() {
        Some(report) => near(
            report.endpoint_work_j,
            0.2 * DEPTH * buoyancy_y * velocity[1],
            2e-8,
        ),
        None => assert_eq!(buoyancy_y, 0.0),
    }
}

#[test]
fn actual_caloric_endpoint_drives_original_force_pressure_and_buoyancy_without_heat_duplication() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 1e5 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let before = old.cells().clone();
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let original = Actuator {
            thermal_gain: 0.2,
            transverse_force: 0.3,
            pressure_gradient: [0.4, 0.0],
        };
        let previous = vec![0.0; system.n_u];
        let history = LaminarEnthalpyHistory {
            velocity_m_s: &previous,
            time_s: 0.0,
        };
        let query = LaminarDriveQuery {
            time_s: 0.2,
            segment: 7,
        };
        let drive = |_: &Cx<'_>, _: LaminarDriveQuery, state: &DepositionSimulation<'_>| {
            Ok::<_, ()>(snapshot(&original, &mesh, state, &cells))
        };
        let solved = owner
            .advance_driven(
                cx,
                &old,
                query,
                history,
                controls(),
                &[0.0; 2],
                &drive,
                |_, t, _, _| Ok(source(t, 0.4)),
                |_, _, _, _| Ok::<_, ()>(heating(&mesh, &law, &old, &cells)),
            )
            .unwrap();
        verify(&system, &cells, &original, &solved, 0.0);
        // Same original caloric inventory, with bounded SOURCE density used
        // only for the independently requested buoyancy contribution.
        let buoyant = owner
            .advance_boussinesq_driven(
                cx,
                &old,
                query,
                history,
                BoussinesqConfig {
                    gravity_m_s2: [0.0, -10.0],
                    max_relative_density_departure: 0.05,
                    max_relative_specific_heat_departure: 0.0,
                },
                controls(),
                &[0.0; 2],
                &drive,
                |_, t, _, _| {
                    ThermalFluidTransportState::try_new(
                        FluidTransportState::try_new(
                            RHO * (1.0 - 0.001 * (t - 300.0)),
                            CP,
                            0.5,
                            0.4,
                        )?,
                        t,
                        None,
                    )
                },
                |_, _, _, _| Ok::<_, ()>(heating(&mesh, &law, &old, &cells)),
            )
            .unwrap();
        verify(&system, &cells, &original, &buoyant, 0.2);
        for sample in buoyant.flow().samples() {
            near(
                sample
                    .material()
                    .source()
                    .resolved_transport()
                    .density_kg_m3(),
                0.99 * RHO,
                2e-9,
            );
            assert_eq!(sample.transport().density_kg_m3(), RHO);
        }
        assert_eq!(old.cells(), &before);
        assert_eq!(old.time_s(), 0.0);
        assert!(previous.iter().all(|&v| v == 0.0));
    });
}

#[test]
fn moving_slip_source_work_uses_original_heat_partition_once_in_physical_depth_and_time() {
    struct Wall {
        boundaries: Vec<LaminarDrivenBoundary>,
        routes: Vec<StokesSlipHeatRoute>,
    }
    impl LaminarMechanicalDrive for Wall {
        fn boundaries(&self) -> &[LaminarDrivenBoundary] {
            &self.boundaries
        }
        fn force_density(&self, _: [f64; 2]) -> [f64; 2] {
            [0.0; 2]
        }
        fn velocity(&self, edge: usize, point: [f64; 2]) -> [f64; 2] {
            let slip = if self.routes.iter().any(|route| route.edge == edge) {
                0.08
            } else {
                0.0
            };
            [0.1 + 0.2 * point[1] + slip, 0.0]
        }
        fn traction(&self, _: usize, _: [f64; 2], normal: [f64; 2]) -> [f64; 2] {
            [
                0.08 * normal[1] - 3.0 * normal[0],
                0.08 * normal[0] - 3.0 * normal[1],
            ]
        }
        fn slip_heat_routes(&self) -> &[StokesSlipHeatRoute] {
            &self.routes
        }
    }
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], true);
        let before = old.cells().clone();
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let velocity = |p: [f64; 2]| [0.1 + 0.2 * p[1], 0.0];
        let previous = moments(&mesh, &velocity);
        let dt = 0.04;
        let solved = owner
            .advance_driven(
                cx,
                &old,
                LaminarDriveQuery {
                    time_s: dt,
                    segment: 0,
                },
                LaminarEnthalpyHistory {
                    velocity_m_s: &previous,
                    time_s: 0.0,
                },
                controls(),
                &[0.0; 3],
                |_, _, _| {
                    let mut routes = Vec::new();
                    let boundaries = (0..mesh.edges.len())
                        .filter(|&e| mesh.is_boundary(e))
                        .map(|edge| {
                            let normal = mesh.edges[edge].normal;
                            if normal[1] > 0.5 {
                                routes.push(StokesSlipHeatRoute {
                                    edge,
                                    wall_cell: WALL,
                                    fluid_heat_fraction: 0.25,
                                });
                                LaminarDrivenBoundary::NavierSlip {
                                    edge,
                                    tangential_resistance_pa_s_per_m: 1.0,
                                }
                            } else if normal[1] < -0.5 {
                                LaminarDrivenBoundary::Velocity { edge }
                            } else {
                                LaminarDrivenBoundary::Open {
                                    edge,
                                    backflow: OpenBoundaryPolicy::NaturalTraction,
                                }
                            }
                        })
                        .collect();
                    Ok::<_, ()>(Wall { boundaries, routes })
                },
                |_, t, _, _| Ok(source(t, 0.4)),
                |_, _, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .unwrap();
        let (error, divergence) =
            system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
        assert!(error < 2e-8 && divergence < 2e-8);
        let bulk = 0.4 * 0.2_f64.powi(2) * DEPTH * dt;
        let slip = 0.4_f64.powi(2) * 0.2_f64.powi(2) * DEPTH * dt;
        let report = solved.report();
        let mechanical = report.mechanical.unwrap();
        near(
            mechanical.velocity_boundary_work_j.unwrap(),
            0.4 * 0.2 * (0.38 - 0.1) * DEPTH * dt,
            2e-10,
        );
        near(
            mechanical.velocity_boundary_work_j.unwrap(),
            bulk + slip,
            2e-10,
        );
        near(mechanical.imposed_traction_work_j, 0.0, 2e-10);
        near(mechanical.backflow_resistance_work_j.unwrap(), 0.0, 1e-12);
        near(report.bulk_viscous_heat_j, bulk, 2e-10);
        near(report.slip_heat_j, slip, 2e-10);
        near(report.other_external_heat_j, 0.0, 1e-12);
        let wall_gain = old.cells()[&WALL].mass_kg
            * (solved.thermal().cells()[&WALL].specific_h_j_kg
                - old.cells()[&WALL].specific_h_j_kg);
        near(wall_gain, 0.75 * slip, 1e-8);
        let stored: f64 = old
            .cells()
            .iter()
            .map(|(&id, cell)| {
                cell.mass_kg
                    * (solved.thermal().cells()[&id].specific_h_j_kg - cell.specific_h_j_kg)
            })
            .sum();
        near(
            stored,
            bulk + slip + report.external_stream_enthalpy_j,
            1e-8,
        );
        near(report.thermal_balance_residual_j, 0.0, 1e-8);
        assert_eq!(old.cells(), &before);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn mechanical_source_and_late_exchange_cancellation_leave_both_native_histories_intact() {
    with_cx(|cx, gate| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let before = old.cells().clone();
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics::default(),
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let original = Actuator {
            thermal_gain: 0.2,
            transverse_force: 0.3,
            pressure_gradient: [0.4, 0.0],
        };
        let previous = vec![0.0; system.n_u];
        let history = LaminarEnthalpyHistory {
            velocity_m_s: &previous,
            time_s: 0.0,
        };
        let query = LaminarDriveQuery {
            time_s: 0.2,
            segment: 0,
        };
        let failed = owner.advance_driven(
            cx,
            &old,
            query,
            history,
            controls(),
            &[0.0; 2],
            |_, _, _| Err::<Snapshot<'_>, _>("mechanical source refused"),
            |_, t, _, _| Ok(source(t, 0.4)),
            |_, _, _, _| Ok::<_, &str>(heating(&mesh, &law, &old, &cells)),
        );
        assert!(matches!(
            failed,
            Err(LaminarEnthalpyError::Drive("mechanical source refused"))
        ));
        let reached_flow = Cell::new(false);
        let cancelled = owner.advance_driven(
            cx,
            &old,
            query,
            history,
            controls(),
            &[0.0; 2],
            |_, _, state| Ok::<_, ()>(snapshot(&original, &mesh, state, &cells)),
            |_, t, _, _| Ok(source(t, 0.4)),
            |_, _, flow, _| {
                assert_eq!(flow.samples().len(), cells.len());
                reached_flow.set(true);
                gate.request();
                Ok::<_, ()>(heating(&mesh, &law, &old, &cells))
            },
        );
        assert!(reached_flow.get());
        assert!(matches!(cancelled, Err(LaminarEnthalpyError::Cancelled)));
        assert_eq!(old.cells(), &before);
        assert_eq!(old.time_s(), 0.0);
        assert!(previous.iter().all(|&v| v == 0.0));
    });
}
