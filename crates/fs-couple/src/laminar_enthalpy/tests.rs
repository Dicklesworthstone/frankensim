//! Independent control-volume, source-rheology and physical heat regressions.

mod adaptive;
mod buoyancy;
mod driven;
mod scheduled;
mod windkessel;

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_flux::TriMesh;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_solid::Mesh2;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, StepConfig, ThermalLaw,
    TimedEvent,
};
use std::collections::BTreeMap;

const RHO: f64 = 2.0;
const CP: f64 = 4.0;
const DEPTH: f64 = 0.3;
const Z: f64 = -0.2;
const MATERIAL: MaterialId = MaterialId(9);
const WALL: CellId = CellId(999);

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 511,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

struct Caloric {
    datum: f64,
}
impl Caloric {
    fn h(&self, t: f64) -> f64 {
        self.datum + CP * (t - 300.0)
    }
}
impl ThermalLaw for Caloric {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut id = [91; 32];
        id[..8].copy_from_slice(&self.datum.to_bits().to_le_bytes());
        Some(id)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + (h - self.datum) / CP;
        if temperature.is_finite() && (280.0..=360.0).contains(&temperature) {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        CP.recip()
    }
}

fn reference(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(RHO, CP, 0.5, mu).unwrap()
}
fn source(t: f64, mu: f64) -> ThermalFluidTransportState {
    ThermalFluidTransportState::try_new(reference(mu), t, None).unwrap()
}
fn zero(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual} versus {expected}"
    );
}
fn moments(mesh: &TriMesh, velocity: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    mesh.edges
        .iter()
        .flat_map(|edge| {
            let a = velocity(mesh.verts[edge.verts.0]);
            let b = velocity(mesh.verts[edge.verts.1]);
            [
                0.5 * ((a[0] + b[0]) * edge.normal[0] + (a[1] + b[1]) * edge.normal[1]),
                ((b[0] - a[0]) * edge.normal[0] + (b[1] - a[1]) * edge.normal[1]) / 12.0,
            ]
        })
        .collect()
}
fn controls() -> LaminarEnthalpyControls {
    LaminarEnthalpyControls {
        momentum: NavierStokesControls {
            maximum_iterations: 120,
            relative_tolerance: 1e-11,
            relaxation: 0.8,
        },
        maximum_active_cells: 16,
        maximum_endpoint_evaluations: 600,
        maximum_thermal_iterations: 20,
        difference_step_j: 1e-4,
        cell_tolerance_j: 1e-8,
        total_tolerance_j: 1e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 0.0,
            forcing_minimum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn inventory<'m>(
    mesh: &TriMesh,
    law: &'m Caloric,
    temperatures: &[f64],
    wall: bool,
) -> (DepositionSimulation<'m>, Vec<CellId>) {
    let cells: Vec<_> = (0..mesh.tris.len())
        .map(|i| CellId(100 + 3 * i as u64))
        .collect();
    let mut events: Vec<_> = cells
        .iter()
        .enumerate()
        .map(|(i, &id)| {
            let [x, y] = mesh.centroids[i];
            let volume = mesh.areas[i] * DEPTH;
            let perimeter: f64 = mesh.tri_edges[i]
                .iter()
                .map(|&(e, _)| mesh.edges[e].len)
                .sum();
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: i as u64,
                    body_id: 1,
                    material_id: MATERIAL,
                    mass_kg: RHO * volume,
                    specific_h_j_kg: law.h(temperatures[i]),
                    geometry: CellGeometry {
                        start_m: [x, y, Z],
                        end_m: [x, y, Z + DEPTH],
                        volume_m3: volume,
                        surface_area_m2: 2.0 * mesh.areas[i] + perimeter * DEPTH,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect();
    if wall {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: WALL,
                segment_id: 99,
                body_id: 2,
                material_id: MATERIAL,
                mass_kg: 1.0,
                specific_h_j_kg: law.h(300.0),
                geometry: CellGeometry {
                    start_m: [0.0, 2.0, 0.0],
                    end_m: [1.0, 2.0, 0.0],
                    volume_m3: 0.5,
                    surface_area_m2: 3.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    let old = DepositionSimulation::new(
        events,
        &[(MATERIAL, law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 16,
            max_links: 100,
            max_events: 100,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    (old, cells)
}
fn map(cells: &[CellId]) -> StokesThermalCells<'_> {
    StokesThermalCells {
        triangle_cells: cells,
        plane_z_m: Z,
    }
}
fn exchange(mesh: &TriMesh, inlet_h: f64) -> LaminarHeatExchange<()> {
    LaminarHeatExchange {
        inlet_enthalpies: (0..mesh.edges.len())
            .map(|edge| {
                mesh.is_boundary(edge).then_some(LaminarInletEnthalpy {
                    material_id: MATERIAL,
                    specific_h_j_kg: inlet_h,
                })
            })
            .collect(),
        heat: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        }
        .into(),
    }
}

#[test]
fn original_bidirectional_face_transport_keeps_zero_net_counterflow() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (thermal, cells) = inventory(&mesh, &law, &[300.0, 310.0], false);
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics::default(),
            &thermal,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let velocity = moments(&mesh, &|p| [p[1] - 0.5, 0.0]);
        let incoming = exchange(&mesh, law.h(320.0)).inlet_enthalpies;
        let advected = owner
            .binding
            .advection::<()>(cx, &system, &velocity, &incoming)
            .unwrap();
        let interior = mesh
            .edges
            .iter()
            .position(|g| g.tris.1 != usize::MAX)
            .unwrap();
        let pairs: Vec<_> = advected
            .flows()
            .iter()
            .filter(|f| f.id / 2 == interior as u64)
            .collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!((pairs[0].a, pairs[0].b), (pairs[1].b, pairs[1].a));
        near(pairs[0].kg_s, RHO * DEPTH / 8.0, 1e-14);
        near(pairs[1].kg_s, RHO * DEPTH / 8.0, 1e-14);
        let original = system.edge_volume_transport(&velocity).unwrap();
        near(original[interior].net_m2_s(), 0.0, 1e-15);
        assert!(advected.endpoint_transport(&thermal, || Ok(())).is_ok());
        let mut missing = incoming;
        let incoming_edge = original
            .iter()
            .enumerate()
            .find(|(edge, rate)| mesh.is_boundary(*edge) && rate.against_normal_m2_s() > 0.0)
            .unwrap()
            .0;
        missing[incoming_edge] = None;
        assert!(matches!(
            owner
                .binding
                .advection::<()>(cx, &system, &velocity, &missing),
            Err(LaminarEnthalpyError::Invalid(_)),
        ));
    });
}

#[test]
fn reversible_resolved_throughflow_matches_two_original_cell_enthalpy_balances() {
    with_cx(|cx, _| {
        for datum in [0.0, 1e6] {
            for speed in [-0.7_f64, 0.7] {
                let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
                let system = FluxSystem::new(&mesh);
                let law = Caloric { datum };
                let (old, cells) = inventory(&mesh, &law, &[298.0, 308.0], false);
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
                let velocity = |_| [speed, 0.0];
                let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
                let boundaries: Vec<_> = mesh
                    .edges
                    .iter()
                    .enumerate()
                    .filter(|(e, _)| mesh.is_boundary(*e))
                    .map(|(edge, geometry)| {
                        if geometry.normal[0].abs() > 0.5 {
                            MaterialNavierStokesBoundary::OpenReservoir {
                                edge,
                                traction: &traction,
                                velocity: &velocity,
                            }
                        } else {
                            MaterialNavierStokesBoundary::Velocity {
                                edge,
                                value: &velocity,
                            }
                        }
                    })
                    .collect();
                let previous = moments(&mesh, &velocity);
                let before = old.cells().clone();
                let dt = 0.04;
                let solved = owner
                    .advance(
                        cx,
                        &old,
                        dt,
                        LaminarEnthalpyInput {
                            history: LaminarEnthalpyHistory {
                                velocity_m_s: &previous,
                                time_s: 0.0,
                            },
                            force_density: &zero,
                            boundaries: &boundaries,
                            slip_heat_routes: &[],
                        },
                        controls(),
                        &[0.0; 2],
                        |_, t, _| Ok(source(t, 0.4)),
                        |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(330.0))),
                    )
                    .unwrap();
                let first = mesh
                    .edges
                    .iter()
                    .find(|e| e.tris.1 == usize::MAX && speed * e.normal[0] < -0.1)
                    .unwrap()
                    .tris
                    .0;
                let second = 1 - first;
                let alpha = 2.0 * speed.abs() * dt;
                let old_t: Vec<_> = cells
                    .iter()
                    .map(|&id| old.temperature_k(id).unwrap())
                    .collect();
                let expected_first = (old_t[first] + alpha * 330.0) / (1.0 + alpha);
                let expected_second = (old_t[second] + alpha * expected_first) / (1.0 + alpha);
                near(
                    solved.thermal().temperature_k(cells[first]).unwrap(),
                    expected_first,
                    2e-8,
                );
                near(
                    solved.thermal().temperature_k(cells[second]).unwrap(),
                    expected_second,
                    2e-8,
                );
                let (error, divergence) =
                    system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
                assert!(error < 1e-9 && divergence < 1e-9);
                near(solved.report().bulk_viscous_heat_j, 0.0, 1e-15);
                near(solved.report().thermal_balance_residual_j, 0.0, 1e-8);
                let expected_heat = RHO * DEPTH * speed.abs() * dt * CP * (330.0 - expected_second);
                near(
                    solved.report().external_stream_enthalpy_j,
                    expected_heat,
                    2e-8,
                );
                for (i, sample) in solved.flow().samples().iter().enumerate() {
                    near(
                        sample.material().source().temperature_k(),
                        solved.thermal().temperature_k(cells[i]).unwrap(),
                        1e-12,
                    );
                }
                assert_eq!(old.cells(), &before);
                assert_eq!(old.time_s(), 0.0);
            }
        }
    });
}

#[test]
fn original_thermal_carreau_closes_convection_viscosity_strain_heat_and_storage() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let caloric = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &caloric, &[300.0; 2], false);
        let material = reference(2.0);
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
        let curve = ThermalCarreauYasuda::new(
            CarreauYasuda::new(2.0, 0.2, 1.2, 0.45, 2.0, [0.0, 1e4]).unwrap(),
            TemperatureShift::new(
                TemperatureShiftModel::Arrhenius {
                    activation_energy_j_mol: 12_000.0,
                },
                300.0,
                [280.0, 360.0],
            )
            .unwrap(),
            ViscosityShiftMode::VerticalAndRate,
        );
        let (gamma, transverse, dt) = (3.0, 0.1, 0.7);
        // Independent declared Carreau/Arrhenius equation, not solver output.
        let viscosity = |t: f64| {
            let a = (12_000.0 / fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * (1.0 / t - 1.0 / 300.0))
                .exp();
            a * (0.2 + 1.8 * (1.0 + (1.2 * a * gamma).powi(2)).powf((0.45 - 1.0) / 2.0))
        };
        let (mut lower, mut upper) = (300.0, 310.0);
        for _ in 0..70 {
            let middle = 0.5 * (lower + upper);
            if RHO * CP * (middle - 300.0) / dt > viscosity(middle) * gamma * gamma {
                upper = middle;
            } else {
                lower = middle;
            }
        }
        let expected_t = 0.5 * (lower + upper);
        let expected_mu = viscosity(expected_t);
        let velocity = |p: [f64; 2]| [0.2 + gamma * p[1], transverse];
        // Nonzero original u.grad(u) = (gamma*transverse,0).
        let force = |_| [RHO * gamma * transverse, 0.0];
        let previous = moments(&mesh, &velocity);
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| MaterialNavierStokesBoundary::Velocity {
                edge,
                value: &velocity,
            })
            .collect();
        let solved = owner
            .advance(
                cx,
                &old,
                dt,
                LaminarEnthalpyInput {
                    history: LaminarEnthalpyHistory {
                        velocity_m_s: &previous,
                        time_s: 0.0,
                    },
                    force_density: &force,
                    boundaries: &boundaries,
                    slip_heat_routes: &[],
                },
                controls(),
                &[0.0; 2],
                |point, t, _| {
                    curve
                        .evaluate(t, point.shear_rate_per_s)?
                        .bind_thermal(source(t, 2.0))
                },
                |_, _, _| Ok::<_, ()>(exchange(&mesh, caloric.h(expected_t))),
            )
            .unwrap();
        assert!(expected_t > 300.1 && expected_mu < viscosity(300.0));
        let (error, divergence) =
            system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
        assert!(error < 2e-8 && divergence < 2e-8, "{error}, {divergence}");
        for (i, sample) in solved.flow().samples().iter().enumerate() {
            near(
                solved.thermal().temperature_k(cells[i]).unwrap(),
                expected_t,
                2e-8,
            );
            near(
                sample.transport().dynamic_viscosity_pa_s(),
                expected_mu,
                2e-8,
            );
            let original = sample.material().source();
            near(
                original.rheology().temperature().temperature_k(),
                expected_t,
                2e-8,
            );
            near(original.rheology().shear_rate_per_s(), gamma, 2e-8);
            near(original.source().temperature_k(), expected_t, 2e-8);
        }
        let expected_heat = expected_mu * gamma * gamma * DEPTH * dt;
        near(solved.report().bulk_viscous_heat_j, expected_heat, 2e-8);
        near(
            solved.report().thermal.stored_energy_change_j,
            expected_heat,
            2e-8,
        );
        near(solved.report().external_stream_enthalpy_j, 0.0, 2e-8);
        near(solved.report().thermal_balance_residual_j, 0.0, 1e-8);
        let top: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].normal[1] > 0.5)
            .collect();
        let reaction = solved
            .flow()
            .step
            .solution
            .boundary_loads(&top, [0.0; 2], &[Some(0.0)])
            .unwrap();
        near(reaction.stress.force_n_per_m[0], -expected_mu * gamma, 3e-8);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn native_slip_heat_is_partitioned_once_without_depositing_numerical_losses() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], true);
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
        let (gamma, mu, beta, dt, fraction) = (0.2, 0.4, 1.0, 0.04, 0.25);
        let velocity = |p: [f64; 2]| [0.1 + gamma * p[1], 0.0];
        let wall = |_| [0.1 + gamma + mu * gamma / beta, 0.0];
        let mut routes = Vec::new();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| {
                if mesh.edges[edge].normal[1] > 0.5 {
                    routes.push(StokesSlipHeatRoute {
                        edge,
                        wall_cell: WALL,
                        fluid_heat_fraction: fraction,
                    });
                    MaterialNavierStokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_pa_s_per_m: beta,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                }
            })
            .collect();
        let previous = moments(&mesh, &velocity);
        let solved = owner
            .advance(
                cx,
                &old,
                dt,
                LaminarEnthalpyInput {
                    history: LaminarEnthalpyHistory {
                        velocity_m_s: &previous,
                        time_s: 0.0,
                    },
                    force_density: &zero,
                    boundaries: &boundaries,
                    slip_heat_routes: &routes,
                },
                controls(),
                &[0.0; 3],
                |_, t, _| Ok(source(t, mu)),
                |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .unwrap();
        let slip = mu * mu * gamma * gamma / beta * DEPTH * dt;
        near(
            solved.report().bulk_viscous_heat_j,
            mu * gamma * gamma * DEPTH * dt,
            1e-10,
        );
        near(solved.report().slip_heat_j, slip, 1e-10);
        let wall_gain =
            solved.thermal().cells()[&WALL].specific_h_j_kg - old.cells()[&WALL].specific_h_j_kg;
        near(wall_gain, (1.0 - fraction) * slip, 1e-9);
        near(solved.report().thermal_balance_residual_j, 0.0, 1e-8);
    });
}

#[test]
fn source_exchange_work_and_cancel_refusals_preserve_both_physical_histories() {
    with_cx(|cx, gate| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
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
        let velocity = |_| [0.2, 0.0];
        let previous = moments(&mesh, &velocity);
        let original_velocity = previous.clone();
        let original_cells = old.cells().clone();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| MaterialNavierStokesBoundary::Velocity {
                edge,
                value: &velocity,
            })
            .collect();
        let input = LaminarEnthalpyInput {
            history: LaminarEnthalpyHistory {
                velocity_m_s: &previous,
                time_s: 0.0,
            },
            force_density: &zero,
            boundaries: &boundaries,
            slip_heat_routes: &[],
        };
        for mode in 0..5 {
            let mut policy = controls();
            if mode == 3 {
                policy.maximum_endpoint_evaluations = 1;
            }
            let attempted = owner.advance(
                cx,
                &old,
                0.04,
                input,
                policy,
                &[0.0; 2],
                |_, temperature, _| {
                    if mode == 0 {
                        return Err(material_error("original material domain refusal"));
                    }
                    let transport = FluidTransportState::try_new(
                        if mode == 1 { RHO * 2.0 } else { RHO },
                        CP,
                        0.5,
                        0.4,
                    )?;
                    ThermalFluidTransportState::try_new(transport, temperature, None)
                },
                |_, _, _| {
                    if mode == 2 {
                        return Err("original heat source refused");
                    }
                    let mut heat = exchange(&mesh, law.h(320.0));
                    if mode == 4 {
                        heat.inlet_enthalpies.fill(None);
                    }
                    Ok(heat)
                },
            );
            match mode {
                0 | 1 => assert!(matches!(attempted, Err(LaminarEnthalpyError::Momentum(_)))),
                2 => assert!(matches!(attempted, Err(LaminarEnthalpyError::Exchange(_)))),
                3 => assert!(matches!(
                    attempted,
                    Err(LaminarEnthalpyError::EndpointBudget)
                )),
                _ => assert!(matches!(attempted, Err(LaminarEnthalpyError::Invalid(_)))),
            }
            assert_eq!(previous, original_velocity);
            assert_eq!(old.cells(), &original_cells);
            assert_eq!(old.time_s(), 0.0);
        }
        let attempted = owner.advance(
            cx,
            &old,
            0.04,
            input,
            controls(),
            &[0.0; 2],
            |_, t, _| Ok(source(t, 0.4)),
            |_, _, _| {
                gate.request();
                Ok::<_, ()>(exchange(&mesh, law.h(320.0)))
            },
        );
        assert!(matches!(attempted, Err(LaminarEnthalpyError::Cancelled)));
        assert_eq!(previous, original_velocity);
        assert_eq!(old.cells(), &original_cells);
        assert_eq!(old.time_s(), 0.0);
    });
}
