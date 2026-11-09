//! Physical Stokes strain heat reaches the same finite thermal state atomically.

#[path = "stokes_motion_heat/planar.rs"]
mod planar;

#[path = "stokes_motion_heat/feedback.rs"]
mod feedback;

#[path = "stokes_motion_heat/planar_feedback.rs"]
mod planar_feedback;

#[path = "stokes_motion_heat/thermal_dilation.rs"]
mod thermal_dilation;

use std::collections::BTreeMap;

use fs_couple::stokes_motion::heating::{
    GuidedStokesHeatError, GuidedStokesObservedHeatError, StokesSlipHeatRoute, StokesThermalCells,
};
use fs_couple::stokes_motion::{
    GuidedStokesConfig, GuidedStokesControls, GuidedStokesError, GuidedStokesMotion,
};
use fs_flux::{FluxSystem, MaterialFluxNumerics, MaterialStokesBoundary, TriMesh};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{MassProperties, Pose, RigidBodyState, Vec3};
use fs_solid::Mesh2;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw,
    ThermalLink, TimedEvent,
};

const RHO: f64 = 2.0;
const CP: f64 = 4.0;
const MU: f64 = 0.4;
const DEPTH: f64 = 0.3;
const Z: f64 = -0.2;
const SPEED: f64 = 0.2;

struct Caloric {
    maximum_h: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (0.0..=self.maximum_h).contains(&h) {
            Ok(300.0 + h / CP)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        CP.recip()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} vs {expected}"
    );
}
fn zero(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn shear(x: [f64; 2]) -> [f64; 2] {
    [SPEED * x[1], 0.0]
}
fn traction(_: [f64; 2], n: [f64; 2]) -> [f64; 2] {
    [MU * SPEED * n[1], MU * SPEED * n[0]]
}
fn config() -> GuidedStokesConfig {
    GuidedStokesConfig {
        direction: [1.0, 0.0],
        depth_m: DEPTH,
        max_displacement_m: 0.05,
        stiffness_n_m: 0.0,
        damping_n_s_m: 0.0,
        velocity_tolerance_m_s: 1e-8,
    }
}
fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.01,
        energy_tolerance_j: 1e-10,
    }
}
fn fluid() -> FluidTransportState {
    // Explicit coherent numerical fixture, not a named-fluid material claim.
    FluidTransportState::try_new(RHO, CP, 0.5, MU).unwrap()
}
fn top_edges(mesh: &TriMesh) -> Vec<usize> {
    (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].normal[1] > 0.5)
        .collect()
}
fn exterior<'a>(mesh: &TriMesh, body: &[usize]) -> Vec<MaterialStokesBoundary<'a>> {
    (0..mesh.edges.len())
        .filter(|e| mesh.is_boundary(*e) && !body.contains(e))
        .map(|edge| {
            if mesh.edges[edge].normal[1] < -0.5 {
                MaterialStokesBoundary::Velocity { edge, value: &zero }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect()
}
fn motion<'s, 'm, 'a, F: FluidTransport + ?Sized>(
    system: &'s FluxSystem<'m>,
    fluid: &'a F,
    edges: &[usize],
    cfg: GuidedStokesConfig,
    field: impl Fn([f64; 2]) -> [f64; 2],
    speed: f64,
) -> GuidedStokesMotion<'s, 'm, 'a, F> {
    let mass = MassProperties::new(0.7, Vec3::ZERO, Vec3::new(0.1, 0.1, 0.1)).unwrap();
    let body = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(mass.mass() * speed, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    let mut history = vec![0.0; system.n_u];
    for (e, g) in system.mesh.edges.iter().enumerate() {
        let a = field(system.mesh.verts[g.verts.0]);
        let b = field(system.mesh.verts[g.verts.1]);
        for i in 0..2 {
            history[2 * e] += 0.5 * (a[i] + b[i]) * g.normal[i];
            history[2 * e + 1] += (b[i] - a[i]) * g.normal[i] / 12.0;
        }
    }
    GuidedStokesMotion::new(
        system,
        fluid,
        MaterialFluxNumerics::default(),
        mass,
        body,
        history,
        edges,
        cfg,
    )
    .unwrap()
}
fn initial_thermal<'a>(
    mesh: &TriMesh,
    law: &'a dyn ThermalLaw,
    linked: bool,
    edit: impl FnMut(usize, &mut DepositedCell),
) -> (DepositionSimulation<'a>, Vec<CellId>) {
    initial_thermal_with_events(mesh, law, linked, edit, Vec::new())
}
fn initial_thermal_with_events<'a>(
    mesh: &TriMesh,
    law: &'a dyn ThermalLaw,
    linked: bool,
    mut edit: impl FnMut(usize, &mut DepositedCell),
    extra_events: Vec<TimedEvent>,
) -> (DepositionSimulation<'a>, Vec<CellId>) {
    let ids: Vec<_> = (0..mesh.tris.len())
        .map(|i| CellId(100 + i as u64))
        .collect();
    let mut events = Vec::new();
    for (i, &id) in ids.iter().enumerate() {
        let [x, y] = mesh.centroids[i];
        let volume = mesh.areas[i] * DEPTH;
        let perimeter: f64 = mesh.tri_edges[i]
            .iter()
            .map(|&(e, _)| mesh.edges[e].len)
            .sum();
        let mut cell = DepositedCell {
            id,
            segment_id: i as u64,
            body_id: 7,
            material_id: MaterialId(9),
            mass_kg: RHO * volume,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [x, y, Z],
                end_m: [x, y, Z + DEPTH],
                volume_m3: volume,
                surface_area_m2: 2.0 * mesh.areas[i] + DEPTH * perimeter,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        };
        edit(i, &mut cell);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        });
    }
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 71,
                a: ids[0],
                b: ids[1],
                area_m2: 0.01,
                conductance_w_k: 0.2,
            }),
        });
        events.push(TimedEvent {
            time_s: 0.02,
            event: Event::Boundary(
                ids[0],
                Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.8,
                },
            ),
        });
    }
    events.extend(extra_events);
    events.sort_by(|a, b| a.time_s.total_cmp(&b.time_s));
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(9), law)],
        DepositionBudget {
            max_cells: 100,
            max_links: 10,
            max_events: 110,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0;
    (thermal, ids)
}
fn binding(ids: &[CellId]) -> StokesThermalCells<'_> {
    StokesThermalCells {
        triangle_cells: ids,
        plane_z_m: Z,
    }
}

#[test]
fn couette_wall_work_heats_original_finite_cells_at_both_mesh_resolutions() {
    let law = Caloric { maximum_h: 100.0 };
    for n in [1, 2] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let fluid = fluid();
        let edges = top_edges(&mesh);
        let exterior = exterior(&mesh, &edges);
        let mut motion = motion(&system, &fluid, &edges, config(), shear, SPEED);
        let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, _| {});
        let step = motion
            .step_with_heating(
                &mut thermal,
                binding(&ids),
                0.05,
                &zero,
                &exterior,
                &[None],
                DEPTH * MU * SPEED,
                GuidedStokesControls::default(),
                thermal_config(),
                || false,
            )
            .unwrap();
        let expected_heat = 0.05 * DEPTH * MU * SPEED * SPEED;
        near(step.motion.speed_m_s, SPEED, 2e-8);
        near(step.deposited_viscous_heat_j, expected_heat, 2e-10);
        assert_eq!(step.deposited_slip_heat_j, 0.0);
        near(
            step.deposited_viscous_heat_j,
            -step.motion.interface_work_j,
            2e-10,
        );
        near(step.thermal.stored_energy_change_j, expected_heat, 2e-10);
        near(step.thermal.external_heat_j, expected_heat, 2e-10);
        near(step.thermal.residual_j, 0.0, 1e-12);
        assert_eq!(thermal.time_s(), motion.state().time_s());
        assert_eq!(thermal.time_s(), 0.05);
        for id in ids {
            near(
                thermal.temperature_k(id).unwrap(),
                300.0 + 0.05 * MU * SPEED * SPEED / (RHO * CP),
                2e-9,
            );
            assert_eq!(thermal.cells()[&id].boundary.source_w, 0.0);
        }
    }
}

#[test]
fn authored_heat_contacts_and_events_survive_and_final_cancellation_rolls_back_both_owners() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let fluid = fluid();
    let law = Caloric { maximum_h: 100.0 };
    let edges = top_edges(&mesh);
    let exterior = exterior(&mesh, &edges);
    let mut control_motion = motion(&system, &fluid, &edges, config(), shear, SPEED);
    let (original, ids) = initial_thermal(&mesh, &law, true, |i, c| {
        if i == 0 {
            c.specific_h_j_kg = 1.0;
            c.boundary.source_w = 0.5;
        }
    });
    let mut control_thermal = original.clone();
    let mut polls = 0;
    let step = control_motion
        .step_with_heating(
            &mut control_thermal,
            binding(&ids),
            0.05,
            &zero,
            &exterior,
            &[None],
            DEPTH * MU * SPEED,
            GuidedStokesControls::default(),
            thermal_config(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let viscous_w = MU * SPEED * SPEED * mesh.areas[0] * DEPTH;
    let mass = original.cells()[&ids[0]].mass_kg;
    let mut h = [1.0, 0.0];
    for k in 0..5 {
        let conduction_w = 0.2 * (h[0] - h[1]) / CP;
        let authored_w = if k < 2 { 0.5 } else { 0.8 };
        h[0] += 0.01 * (viscous_w + authored_w - conduction_w) / mass;
        h[1] += 0.01 * (viscous_w + conduction_w) / mass;
    }
    for (i, &id) in ids.iter().enumerate() {
        near(control_thermal.cells()[&id].specific_h_j_kg, h[i], 2e-10);
    }
    near(
        step.thermal.stored_energy_change_j,
        0.02 * 0.5 + 0.03 * 0.8 + 0.05 * DEPTH * MU * SPEED * SPEED,
        2e-10,
    );
    assert_eq!(step.thermal.events, 1);
    assert_eq!(control_thermal.links(), original.links());
    assert_eq!(control_thermal.cells()[&ids[0]].boundary.source_w, 0.8);

    let mut retry_motion = motion(&system, &fluid, &edges, config(), shear, SPEED);
    let old_motion = retry_motion.state().clone();
    let mut retry_thermal = original.clone();
    let mut attempts = 0;
    let failure = retry_motion.step_with_heating(
        &mut retry_thermal,
        binding(&ids),
        0.05,
        &zero,
        &exterior,
        &[None],
        DEPTH * MU * SPEED,
        GuidedStokesControls::default(),
        thermal_config(),
        || {
            attempts += 1;
            attempts == polls
        },
    );
    assert!(matches!(
        failure,
        Err(GuidedStokesHeatError::Motion(GuidedStokesError::Cancelled))
    ));
    assert_eq!(retry_motion.state(), &old_motion);
    assert_eq!(retry_thermal.cells(), original.cells());
    assert_eq!(retry_thermal.links(), original.links());
    assert_eq!(retry_thermal.time_s(), 0.0);
    retry_motion
        .step_with_heating(
            &mut retry_thermal,
            binding(&ids),
            0.05,
            &zero,
            &exterior,
            &[None],
            DEPTH * MU * SPEED,
            GuidedStokesControls::default(),
            thermal_config(),
            || false,
        )
        .unwrap();
    assert_eq!(retry_motion.state(), control_motion.state());
    assert_eq!(retry_thermal.cells(), control_thermal.cells());
    assert_eq!(retry_thermal.time_s(), control_thermal.time_s());
}

#[test]
fn late_caloric_refusal_preserves_motion_and_thermal_history_then_short_retry_succeeds() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let fluid = fluid();
    let law = Caloric { maximum_h: 0.0002 };
    let edges = top_edges(&mesh);
    let exterior = exterior(&mesh, &edges);
    let mut motion = motion(&system, &fluid, &edges, config(), shear, SPEED);
    let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, _| {});
    let old_motion = motion.state().clone();
    let old_cells = thermal.cells().clone();
    let refusal = motion.step_with_heating(
        &mut thermal,
        binding(&ids),
        0.05,
        &zero,
        &exterior,
        &[None],
        DEPTH * MU * SPEED,
        GuidedStokesControls::default(),
        thermal_config(),
        || false,
    );
    assert!(matches!(
        refusal,
        Err(GuidedStokesHeatError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(motion.state(), &old_motion);
    assert_eq!(thermal.cells(), &old_cells);
    assert_eq!(thermal.time_s(), 0.0);
    motion
        .step_with_heating(
            &mut thermal,
            binding(&ids),
            0.01,
            &zero,
            &exterior,
            &[None],
            DEPTH * MU * SPEED,
            GuidedStokesControls::default(),
            thermal_config(),
            || false,
        )
        .unwrap();
    for id in ids {
        near(
            thermal.cells()[&id].specific_h_j_kg,
            0.01 * MU * SPEED * SPEED / RHO,
            2e-10,
        );
    }
}

#[test]
fn rigid_translation_does_not_deposit_body_damping_or_backward_euler_loss() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let fluid = fluid();
    let law = Caloric { maximum_h: 100.0 };
    let edges: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    let cfg = GuidedStokesConfig {
        stiffness_n_m: 7.0,
        damping_n_s_m: 0.4,
        ..config()
    };
    let mut motion = motion(&system, &fluid, &edges, cfg, |_| [0.03, 0.0], 0.03);
    let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, _| {});
    let step = motion
        .step_with_heating(
            &mut thermal,
            binding(&ids),
            0.002,
            &zero,
            &[],
            &[Some(91.0)],
            1.3,
            GuidedStokesControls::default(),
            thermal_config(),
            || false,
        )
        .unwrap();
    assert!(step.motion.damping_loss_j > 1e-8);
    assert!(step.motion.mechanical_time_discretization_loss_j > 1e-8);
    near(step.deposited_viscous_heat_j, 0.0, 1e-18);
    near(step.thermal.stored_energy_change_j, 0.0, 1e-18);
    for id in ids {
        near(thermal.temperature_k(id).unwrap(), 300.0, 1e-12);
    }
}

#[test]
fn mismatched_or_ambiguous_receivers_refuse_before_flow_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let fluid = fluid();
    let law = Caloric { maximum_h: 100.0 };
    let edges = top_edges(&mesh);
    let exterior = exterior(&mesh, &edges);
    let mut motion = motion(&system, &fluid, &edges, config(), shear, SPEED);
    let original = motion.state().clone();
    for case in 0..8 {
        let (mut thermal, mut ids) = initial_thermal(&mesh, &law, false, |i, c| {
            if i == 0 {
                match case {
                    3 => c.mass_kg *= 1.1,
                    4 => c.geometry.volume_m3 *= 1.1,
                    5 => c.geometry.surface_area_m2 *= 1.1,
                    6 => c.geometry.end_m[0] += 0.01,
                    _ => {}
                }
            }
        });
        match case {
            0 => {
                ids.pop();
            }
            1 => ids[1] = ids[0],
            2 => ids.swap(0, 1),
            7 => {
                thermal = thermal
                    .advance_to(0.001, thermal_config(), || Ok(()))
                    .unwrap()
                    .0;
            }
            _ => {}
        }
        let before = thermal.cells().clone();
        let queries = std::cell::Cell::new(0);
        let force = |_| {
            queries.set(queries.get() + 1);
            [0.0; 2]
        };
        let result = motion.step_with_heating(
            &mut thermal,
            binding(&ids),
            0.05,
            &force,
            &exterior,
            &[None],
            DEPTH * MU * SPEED,
            GuidedStokesControls::default(),
            thermal_config(),
            || false,
        );
        assert!(
            matches!(result, Err(GuidedStokesHeatError::Binding(_))),
            "case {case}"
        );
        assert_eq!(queries.get(), 0);
        assert_eq!(motion.state(), &original);
        assert_eq!(thermal.cells(), &before);
    }
}

const WALL: CellId = CellId(900);

fn wall_cell(depth: f64) -> DepositedCell {
    DepositedCell {
        id: WALL,
        segment_id: 900,
        body_id: 8,
        material_id: MaterialId(9),
        mass_kg: 3.0 * 0.1 * depth,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.5, -0.05, Z],
            end_m: [0.5, -0.05, Z + depth],
            volume_m3: 0.1 * depth,
            surface_area_m2: 0.2 + 2.2 * depth,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn slip_exterior<'a>(
    mesh: &TriMesh,
    body: &[usize],
    beta: f64,
    wall: &'a dyn Fn([f64; 2]) -> [f64; 2],
    traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
) -> Vec<MaterialStokesBoundary<'a>> {
    (0..mesh.edges.len())
        .filter(|e| mesh.is_boundary(*e) && !body.contains(e))
        .map(|edge| {
            if mesh.edges[edge].normal[1] < -0.5 {
                MaterialStokesBoundary::NavierSlip {
                    edge,
                    value: wall,
                    tangential_resistance_pa_s_per_m: beta,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: traction,
                }
            }
        })
        .collect()
}

fn slip_routes(exterior: &[MaterialStokesBoundary<'_>], fraction: f64) -> Vec<StokesSlipHeatRoute> {
    exterior
        .iter()
        .filter_map(|boundary| match *boundary {
            MaterialStokesBoundary::NavierSlip { edge, .. } => Some(StokesSlipHeatRoute {
                edge,
                wall_cell: WALL,
                fluid_heat_fraction: fraction,
            }),
            _ => None,
        })
        .collect()
}

struct SingleRead(std::cell::Cell<usize>);
impl FluidTransport for SingleRead {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        assert_eq!(self.0.replace(self.0.get() + 1), 0, "source requeried");
        Ok(fluid())
    }
}

#[test]
fn navier_slip_wall_work_reaches_fluid_and_wall_with_explicit_partition() {
    let law = Caloric { maximum_h: 100.0 };
    // Partition endpoints, shared wall receiver, depth scaling, uniform frame
    // translation and perfect slip each retain the same independent solution.
    for (n, beta, fraction, boost, depth) in [
        (1, 0.4, 0.0, 0.0, 0.3),
        (1, 0.4, 0.25, 0.0, 0.3),
        (1, 0.4, 1.0, 0.0, 0.3),
        (2, 0.4, 0.25, 0.3, 0.6),
        (2, 0.0, 0.25, 0.3, 0.6),
    ] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let source = SingleRead(std::cell::Cell::new(0));
        let edges = top_edges(&mesh);
        let gamma = SPEED * beta / (MU + beta);
        let slip = SPEED - gamma;
        let stress = MU * gamma;
        let wall_velocity = |_| [boost, 0.0];
        let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
        let exterior = slip_exterior(&mesh, &edges, beta, &wall_velocity, &traction);
        let routes = slip_routes(&exterior, fraction);
        let cfg = GuidedStokesConfig {
            depth_m: depth,
            ..config()
        };
        let mut motion = motion(
            &system,
            &source,
            &edges,
            cfg,
            |x| [boost + slip + gamma * x[1], 0.0],
            boost + SPEED,
        );
        let mut wall = wall_cell(depth);
        wall.boundary.source_w = 0.002;
        let (mut thermal, ids) = initial_thermal_with_events(
            &mesh,
            &law,
            false,
            |i, c| {
                c.mass_kg *= depth / DEPTH;
                c.geometry.volume_m3 *= depth / DEPTH;
                c.geometry.end_m[2] = Z + depth;
                let perimeter: f64 = mesh.tri_edges[i]
                    .iter()
                    .map(|&(e, _)| mesh.edges[e].len)
                    .sum();
                c.geometry.surface_area_m2 = 2.0 * mesh.areas[i] + perimeter * depth;
                if i == 0 {
                    c.boundary.source_w = 0.003;
                }
            },
            vec![
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(wall.clone()),
                },
                TimedEvent {
                    time_s: 0.02,
                    event: Event::Boundary(
                        WALL,
                        Boundary {
                            source_w: 0.004,
                            ..wall.boundary
                        },
                    ),
                },
            ],
        );
        let dt = 0.05;
        let step = motion
            .step_with_slip_heating(
                &mut thermal,
                binding(&ids),
                &routes,
                dt,
                &zero,
                &exterior,
                &[None],
                depth * stress,
                GuidedStokesControls::default(),
                thermal_config(),
                || false,
            )
            .unwrap();
        let bulk_heat = dt * depth * MU * gamma * gamma;
        let slip_heat = dt * depth * beta * slip * slip;
        near(step.motion.speed_m_s, boost + SPEED, 2e-8);
        near(step.deposited_viscous_heat_j, bulk_heat, 2e-10);
        near(step.deposited_slip_heat_j, slip_heat, 2e-10);
        let bottom: Vec<_> = routes.iter().map(|route| route.edge).collect();
        let load = step
            .motion
            .flow
            .boundary_loads(&bottom, [0.0; 2], &[None])
            .unwrap();
        near(
            -step.motion.interface_work_j - dt * depth * boost * load.stress.force_n_per_m[0],
            bulk_heat + slip_heat,
            2e-10,
        );
        let authored_wall = 0.02 * 0.002 + 0.03 * 0.004;
        let total_heat = bulk_heat + slip_heat + dt * 0.003 + authored_wall;
        near(step.thermal.stored_energy_change_j, total_heat, 2e-10);
        near(step.thermal.external_heat_j, total_heat, 2e-10);
        near(step.thermal.residual_j, 0.0, 1e-12);
        for (i, id) in ids.iter().enumerate() {
            let edge_length: f64 = routes
                .iter()
                .filter(|route| mesh.edges[route.edge].tris.0 == i)
                .map(|route| mesh.edges[route.edge].len)
                .sum();
            let cell_heat = bulk_heat * mesh.areas[i]
                + slip_heat * fraction * edge_length
                + if i == 0 { dt * 0.003 } else { 0.0 };
            near(
                thermal.cells()[id].specific_h_j_kg,
                cell_heat / thermal.cells()[id].mass_kg,
                2e-9,
            );
        }
        near(
            thermal.cells()[&WALL].specific_h_j_kg,
            (slip_heat * (1.0 - fraction) + authored_wall) / wall.mass_kg,
            2e-9,
        );
        assert_eq!(thermal.cells()[&WALL].boundary.source_w, 0.004);
        assert_eq!(thermal.cells()[&ids[0]].boundary.source_w, 0.003);
        assert_eq!(step.thermal.events, 1);
        assert_eq!(thermal.time_s(), motion.state().time_s());
        assert_eq!(source.0.get(), 1);
    }
}

#[test]
fn slip_routes_refuse_before_callbacks_and_late_cancellation_rolls_back_every_receiver() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let fluid = fluid();
    let law = Caloric { maximum_h: 100.0 };
    let edges = top_edges(&mesh);
    let gamma = SPEED / 2.0;
    let field = |x: [f64; 2]| [gamma * (1.0 + x[1]), 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [MU * gamma * n[1], MU * gamma * n[0]];
    let exterior = slip_exterior(&mesh, &edges, MU, &zero, &traction);
    let routes = slip_routes(&exterior, 0.25);
    let (original, ids) = initial_thermal_with_events(
        &mesh,
        &law,
        false,
        |_, _| {},
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(wall_cell(DEPTH)),
        }],
    );
    let mut owner = motion(&system, &fluid, &edges, config(), field, SPEED);
    let original_motion = owner.state().clone();
    for case in 0..9 {
        let mut invalid = routes.clone();
        match case {
            0 => invalid.clear(),
            1 => invalid.push(invalid[0]),
            2 => invalid[0].edge = edges[0],
            3 => invalid[0].wall_cell = ids[0],
            4 => invalid[0].wall_cell = CellId(999),
            5 => invalid[0].fluid_heat_fraction = f64::NAN,
            6 => invalid[0].fluid_heat_fraction = -0.1,
            7 => invalid[0].fluid_heat_fraction = 1.1,
            _ => invalid[0].edge = usize::MAX,
        }
        let mut thermal = original.clone();
        let queries = std::cell::Cell::new(0);
        let force = |_| {
            queries.set(queries.get() + 1);
            [0.0; 2]
        };
        let result = owner.step_with_slip_heating(
            &mut thermal,
            binding(&ids),
            &invalid,
            0.05,
            &force,
            &exterior,
            &[None],
            DEPTH * MU * gamma,
            GuidedStokesControls::default(),
            thermal_config(),
            || false,
        );
        assert!(
            matches!(result, Err(GuidedStokesHeatError::Binding(_))),
            "case {case}"
        );
        assert_eq!(queries.get(), 0);
        assert_eq!(owner.state(), &original_motion);
        assert_eq!(thermal.cells(), original.cells());
        assert_eq!(thermal.time_s(), 0.0);
    }
    let mut control_thermal = original.clone();
    let mut polls = 0;
    owner
        .step_with_slip_heating(
            &mut control_thermal,
            binding(&ids),
            &routes,
            0.05,
            &zero,
            &exterior,
            &[None],
            DEPTH * MU * gamma,
            GuidedStokesControls::default(),
            thermal_config(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(control_thermal.cells()[&WALL].specific_h_j_kg > 0.0);
    let mut retry = motion(&system, &fluid, &edges, config(), field, SPEED);
    let mut retry_thermal = original.clone();
    let mut attempts = 0;
    let result = retry.step_with_slip_heating(
        &mut retry_thermal,
        binding(&ids),
        &routes,
        0.05,
        &zero,
        &exterior,
        &[None],
        DEPTH * MU * gamma,
        GuidedStokesControls::default(),
        thermal_config(),
        || {
            attempts += 1;
            attempts == polls
        },
    );
    assert!(matches!(
        result,
        Err(GuidedStokesHeatError::Motion(GuidedStokesError::Cancelled))
    ));
    assert_eq!(retry.state(), &original_motion);
    assert_eq!(retry_thermal.cells(), original.cells());
    assert_eq!(retry_thermal.time_s(), 0.0);
    retry
        .step_with_slip_heating(
            &mut retry_thermal,
            binding(&ids),
            &routes,
            0.05,
            &zero,
            &exterior,
            &[None],
            DEPTH * MU * gamma,
            GuidedStokesControls::default(),
            thermal_config(),
            || false,
        )
        .unwrap();
    assert_eq!(retry.state(), owner.state());
    assert_eq!(retry_thermal.cells(), control_thermal.cells());
}


#[derive(Clone, Debug, PartialEq)]
enum CoupledHeatRefusal {
    Thermal(DepositionError),
    Endpoint,
}
impl From<DepositionError> for CoupledHeatRefusal {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}

#[derive(Clone, Debug, PartialEq)]
struct ClockedHeat {
    cell: CellId,
    starts_and_durations: Vec<(f64, f64)>,
    prepared: Option<(f64, f64, f64)>,
    heat_j: f64,
    begins: usize,
    events: usize,
    refuse_at_s: Option<f64>,
}
impl ClockedHeat {
    fn new(cell: CellId) -> Self {
        Self {
            cell, starts_and_durations: Vec::new(), prepared: None,
            heat_j: 0.0, begins: 0, events: 0, refuse_at_s: None,
        }
    }
    fn source_w(time_s: f64) -> f64 {
        if time_s < 0.015 { 0.01 } else { 0.03 }
    }
}
impl DepositionObserver for ClockedHeat {
    type Error = CoupledHeatRefusal;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.begins += 1;
        Ok(())
    }
    fn next_transition_time_s(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        checkpoint()?;
        Ok((thermal.time_s() < 0.015).then_some(0.015))
    }
    fn sources(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        Ok(BTreeMap::from([(self.cell, Self::source_w(thermal.time_s()))]))
    }
    fn conductances(
        &mut self, _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        checkpoint()?;
        Ok(BTreeMap::from([(71, 0.4)]))
    }
    fn limit_timestep_s(
        &mut self, _: &DepositionSimulation<'_>, proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        checkpoint()?;
        Ok(proposed_s.min(0.007))
    }
    fn requires_step_interval(&self) -> bool { true }
    fn step_sources(
        &mut self, thermal: &DepositionSimulation<'_>, start_s: f64, dt_s: f64, end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        assert_eq!(thermal.time_s(), start_s);
        assert!(self.prepared.is_none());
        self.prepared = Some((start_s, dt_s, end_s));
        // Exact average of the additional ramp q(t) = 0.4 t watts.
        Ok(BTreeMap::from([(self.cell, 0.2 * (start_s + end_s))]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.events += 1;
        Ok(())
    }
    fn step(
        &mut self, start_s: f64, dt_s: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        let (start, duration, end) = self.prepared.take().expect("interval hook was forwarded");
        assert_eq!((start, duration, end), (start_s, dt_s, thermal.time_s()));
        self.heat_j += dt_s * (Self::source_w(start_s) + 0.2 * (start_s + end));
        self.starts_and_durations.push((start_s, dt_s));
        if self.refuse_at_s.is_some_and(|time| end >= time) {
            return Err(CoupledHeatRefusal::Endpoint);
        }
        Ok(())
    }
}

#[test]
fn observed_stokes_heat_shares_physical_clock_and_commits_all_histories_atomically() {
    let law = Caloric { maximum_h: 100.0 };
    for with_slip in [false, true] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let fluid = fluid();
        let edges = top_edges(&mesh);
        let gamma = if with_slip { SPEED / 2.0 } else { SPEED };
        let slip = SPEED - gamma;
        let field = |x: [f64; 2]| [slip + gamma * x[1], 0.0];
        let stress = MU * gamma;
        let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
        let boundaries = if with_slip {
            slip_exterior(&mesh, &edges, MU, &zero, &traction)
        } else {
            exterior(&mesh, &edges)
        };
        let routes = slip_routes(&boundaries, 0.25);
        let (original, ids) = initial_thermal_with_events(
            &mesh, &law, true,
            |i, cell| { if i == 0 { cell.specific_h_j_kg = 1.0; } },
            vec![TimedEvent { time_s: 0.0, event: Event::Birth(wall_cell(DEPTH)) }],
        );
        let advance = |
            owner: &mut GuidedStokesMotion<'_, '_, '_, FluidTransportState>,
            thermal: &mut DepositionSimulation<'_>,
            observer: &mut ClockedHeat,
            cancelled: &mut dyn FnMut() -> bool,
        | {
            if with_slip {
                owner.step_with_slip_heating_observed(
                    thermal, observer, binding(&ids), &routes, 0.05, &zero, &boundaries,
                    &[None], DEPTH * stress, GuidedStokesControls::default(),
                    thermal_config(), cancelled,
                )
            } else {
                owner.step_with_heating_observed(
                    thermal, observer, binding(&ids), 0.05, &zero, &boundaries,
                    &[None], DEPTH * stress, GuidedStokesControls::default(),
                    thermal_config(), cancelled,
                )
            }
        };

        let mut owner = motion(&system, &fluid, &edges, config(), field, SPEED);
        let mut thermal = original.clone();
        let mut observer = ClockedHeat::new(ids[0]);
        let mut polls = 0;
        let step = advance(&mut owner, &mut thermal, &mut observer, &mut || {
            polls += 1;
            false
        }).unwrap();
        let bulk_heat = 0.05 * DEPTH * MU * gamma * gamma;
        let slip_heat = 0.05 * DEPTH * MU * slip * slip;
        let added_heat = 0.015 * 0.01 + 0.035 * 0.03 + 0.2 * 0.05 * 0.05;
        let total_heat = bulk_heat + slip_heat + added_heat + 0.03 * 0.8;
        near(step.deposited_viscous_heat_j, bulk_heat, 2e-10);
        near(step.deposited_slip_heat_j, slip_heat, 2e-10);
        near(observer.heat_j, added_heat, 1e-14);
        near(step.thermal.external_heat_j, total_heat, 2e-10);
        near(step.thermal.stored_energy_change_j, total_heat, 2e-10);
        near(step.thermal.residual_j, 0.0, 1e-12);
        assert_eq!(observer.begins, 1);
        assert_eq!(observer.events, 1);
        assert_eq!(step.thermal.events, 1);
        assert_eq!(step.thermal.substeps, observer.starts_and_durations.len());
        assert_eq!(thermal.time_s(), 0.05);
        assert_eq!(thermal.time_s(), owner.state().time_s());
        assert!(observer.starts_and_durations.iter().any(|&(start, _)| start == 0.015));
        assert!(observer.starts_and_durations.iter().any(|&(start, _)| start == 0.02));
        assert!(observer.starts_and_durations.iter().all(|&(_, dt)| dt <= 0.007 + 1e-15));
        assert_eq!(thermal.links(), original.links());

        // Independent explicit energy recurrence checks the overridden link,
        // pulse, ramp and fluid-side slip heat on the actual accepted intervals.
        let mut h = [1.0, 0.0];
        for &(start, dt) in &observer.starts_and_durations {
            let conduction_w = 0.4 * (h[0] - h[1]) / CP;
            for i in 0..2 {
                let edge_length: f64 = routes.iter()
                    .filter(|route| mesh.edges[route.edge].tris.0 == i)
                    .map(|route| mesh.edges[route.edge].len).sum();
                let mut watts = DEPTH * MU * gamma * gamma * mesh.areas[i]
                    + DEPTH * MU * slip * slip * 0.25 * edge_length;
                watts += if i == 0 { -conduction_w } else { conduction_w };
                if i == 0 {
                    watts += ClockedHeat::source_w(start) + 0.2 * (2.0 * start + dt);
                    if start >= 0.02 { watts += 0.8; }
                }
                h[i] += dt * watts / original.cells()[&ids[i]].mass_kg;
            }
        }
        for (i, id) in ids.iter().enumerate() {
            near(thermal.cells()[id].specific_h_j_kg, h[i], 2e-9);
        }
        near(
            thermal.cells()[&WALL].specific_h_j_kg,
            0.75 * slip_heat / thermal.cells()[&WALL].mass_kg,
            2e-9,
        );

        // Both a late typed refusal and the final cancellation must discard
        // already completed thermal substeps and all tentative observer history.
        for refuse in [true, false] {
            let mut retry = motion(&system, &fluid, &edges, config(), field, SPEED);
            let old_motion = retry.state().clone();
            let mut retry_thermal = original.clone();
            let mut retry_observer = ClockedHeat::new(ids[0]);
            retry_observer.refuse_at_s = refuse.then_some(0.03);
            let old_observer = retry_observer.clone();
            let mut attempts = 0;
            let result = advance(
                &mut retry, &mut retry_thermal, &mut retry_observer, &mut || {
                    attempts += 1;
                    !refuse && attempts == polls
                },
            );
            if refuse {
                assert!(matches!(result, Err(GuidedStokesObservedHeatError::Observer(
                    CoupledHeatRefusal::Endpoint
                ))));
            } else {
                assert!(matches!(result, Err(GuidedStokesObservedHeatError::Heat(
                    GuidedStokesHeatError::Motion(GuidedStokesError::Cancelled)
                ))));
            }
            assert_eq!(retry.state(), &old_motion);
            assert_eq!(retry_thermal.time_s(), original.time_s());
            assert_eq!(retry_thermal.cells(), original.cells());
            assert_eq!(retry_thermal.links(), original.links());
            assert_eq!(retry_observer, old_observer);

            retry_observer.refuse_at_s = None;
            advance(&mut retry, &mut retry_thermal, &mut retry_observer, &mut || false).unwrap();
            assert_eq!(retry.state(), owner.state());
            assert_eq!(retry_thermal.cells(), thermal.cells());
            assert_eq!(retry_thermal.time_s(), thermal.time_s());
            assert_eq!(retry_observer, observer);
        }
    }
}
