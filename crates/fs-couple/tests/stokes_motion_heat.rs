//! Physical Stokes strain heat reaches the same finite thermal state atomically.

use fs_couple::stokes_motion::heating::{GuidedStokesHeatError, StokesThermalCells};
use fs_couple::stokes_motion::{
    GuidedStokesConfig, GuidedStokesControls, GuidedStokesError, GuidedStokesMotion,
};
use fs_flux::{FluxSystem, MaterialFluxNumerics, MaterialStokesBoundary, TriMesh};
use fs_material::fluid_transport::FluidTransportState;
use fs_mbd::{MassProperties, Pose, RigidBodyState, Vec3};
use fs_solid::Mesh2;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
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
fn motion<'s, 'm, 'a>(
    system: &'s FluxSystem<'m>,
    fluid: &'a FluidTransportState,
    edges: &[usize],
    cfg: GuidedStokesConfig,
    field: impl Fn([f64; 2]) -> [f64; 2],
    speed: f64,
) -> GuidedStokesMotion<'s, 'm, 'a, FluidTransportState> {
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
    mut edit: impl FnMut(usize, &mut DepositedCell),
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
