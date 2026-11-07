//! G1 analytic heat, G3 owner parity and G4 atomic planar heat composition.

use super::{
    CP, DEPTH, MU, RHO, SPEED, WALL, Caloric, ClockedHeat, CoupledHeatRefusal, SingleRead,
    binding, config, exterior, fluid,
    initial_thermal, initial_thermal_with_events, motion, near, shear, slip_exterior,
    slip_routes, thermal_config, top_edges, wall_cell, zero,
};
use fs_couple::stokes_motion::GuidedStokesControls;
use fs_couple::stokes_motion::heating::{
    PlanarStokesHeatError, PlanarStokesHeatStep, PlanarStokesObservedHeatError,
};
use fs_couple::stokes_motion::planar::{
    PlanarStokesConfig, PlanarStokesControls, PlanarStokesError, PlanarStokesMotion,
};
use fs_flux::{FluxSystem, MaterialFluxNumerics, MaterialStokesBoundary, TriMesh};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solid::Mesh2;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, Event, TimedEvent};
use std::cell::Cell;

fn planar_config() -> PlanarStokesConfig {
    PlanarStokesConfig {
        depth_m: DEPTH,
        max_translation_m: 0.05,
        max_rotation_rad: 0.05,
        translation_stiffness_n_m: [0.0; 2],
        rotation_stiffness_n_m_rad: 0.0,
        translation_damping_n_s_m: [0.0; 2],
        rotation_damping_n_m_s_rad: 0.0,
        velocity_tolerance_m_s: 1e-8,
        angular_velocity_tolerance_rad_s: 1e-8,
    }
}

pub(super) fn planar_motion<'s, 'm, 'a, F: FluidTransport + ?Sized>(
    system: &'s FluxSystem<'m>,
    fluid: &'a F,
    edges: &[usize],
    field: impl Fn([f64; 2]) -> [f64; 2],
    rates: [f64; 3],
) -> PlanarStokesMotion<'s, 'm, 'a, F> {
    let mass = MassProperties::new(0.7, Vec3::ZERO, Vec3::new(0.1, 0.1, 0.1)).unwrap();
    // COM lies on the centre of the top edge, so uniform shear exerts no
    // moment there. Rotation still prescribes the edge's linear normal moment.
    let body = RigidBodyState::new(
        Pose::new(Vec3::new(0.5, 1.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(mass.mass() * rates[0], mass.mass() * rates[1], 0.0),
        Vec3::new(0.0, 0.0, 0.1 * rates[2]),
    )
    .unwrap();
    let mut history = vec![0.0; system.n_u];
    for (edge, geometry) in system.mesh.edges.iter().enumerate() {
        let a = field(system.mesh.verts[geometry.verts.0]);
        let b = field(system.mesh.verts[geometry.verts.1]);
        for axis in 0..2 {
            history[2 * edge] += 0.5 * (a[axis] + b[axis]) * geometry.normal[axis];
            history[2 * edge + 1] += (b[axis] - a[axis]) * geometry.normal[axis] / 12.0;
        }
    }
    PlanarStokesMotion::new(
        system, fluid, MaterialFluxNumerics::default(), mass, body, history, edges,
        planar_config(),
    )
    .unwrap()
}

fn shear_heat(
    motion: &mut PlanarStokesMotion<'_, '_, '_, FluidTransportState>,
    thermal: &mut DepositionSimulation<'_>,
    ids: &[CellId],
    exterior: &[MaterialStokesBoundary<'_>],
    dt: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<PlanarStokesHeatStep, PlanarStokesHeatError> {
    motion.step_with_heating(
        thermal, binding(ids), dt, &zero, exterior, &[None],
        [DEPTH * MU * SPEED, 0.0], 0.0, PlanarStokesControls::default(),
        thermal_config(), cancelled,
    )
}

#[test]
fn g3_planar_translation_preserves_guided_heat_and_native_motion() {
    let law = Caloric { maximum_h: 100.0 };
    for n in [1, 2] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let fluid = fluid();
        let edges = top_edges(&mesh);
        let exterior = exterior(&mesh, &edges);
        let mut planar = planar_motion(&system, &fluid, &edges, shear, [SPEED, 0.0, 0.0]);
        let mut plain = planar_motion(&system, &fluid, &edges, shear, [SPEED, 0.0, 0.0]);
        let mut guided = motion(&system, &fluid, &edges, config(), shear, SPEED);
        let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, _| {});
        let mut guided_thermal = thermal.clone();
        let step = shear_heat(&mut planar, &mut thermal, &ids, &exterior, 0.05, || false)
            .unwrap();
        let guided_step = guided.step_with_heating(
            &mut guided_thermal, binding(&ids), 0.05, &zero, &exterior, &[None],
            DEPTH * MU * SPEED, GuidedStokesControls::default(), thermal_config(), || false,
        ).unwrap();
        plain.step(
            0.05, &zero, &exterior, &[None], [DEPTH * MU * SPEED, 0.0], 0.0,
            PlanarStokesControls::default(), || false,
        ).unwrap();
        assert_eq!(planar.state(), plain.state(), "thermal staging must not step motion twice");
        near(step.motion.velocity_m_s[0], guided_step.motion.speed_m_s, 2e-8);
        near(step.motion.angular_velocity_rad_s, 0.0, 2e-8);
        for (actual, expected) in step.motion.flow.velocity_dofs().iter()
            .zip(guided_step.motion.flow.velocity_dofs())
        {
            near(*actual, *expected, 2e-8);
        }
        near(step.deposited_viscous_heat_j, 0.05 * DEPTH * MU * SPEED * SPEED, 2e-10);
        near(step.deposited_viscous_heat_j, guided_step.deposited_viscous_heat_j, 2e-10);
        assert_eq!(step.deposited_slip_heat_j, 0.0);
        near(step.thermal.stored_energy_change_j, step.deposited_viscous_heat_j, 2e-10);
        for id in ids {
            near(thermal.temperature_k(id).unwrap(), guided_thermal.temperature_k(id).unwrap(), 2e-9);
        }
        assert_eq!(thermal.time_s(), planar.state().time_s());
    }
}

#[test]
fn g1_rotating_boundary_heats_by_strain_while_rigid_rotation_produces_no_heat() {
    let law = Caloric { maximum_h: 100.0 };
    for omega in [-SPEED, 0.12] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
        let system = FluxSystem::new(&mesh);
        let fluid = fluid();
        let edges = top_edges(&mesh);
        // Steady divergence-free affine field. At the top it is exactly the
        // body's translation plus rotation about (0.5, 1). Its symmetric
        // strain has 2 mu D:D = mu * (SPEED + omega)^2 everywhere.
        let field = |x: [f64; 2]| [SPEED * x[1], omega * (x[0] - 0.5)];
        let stress = MU * (SPEED + omega);
        let traction = |_: [f64; 2], n: [f64; 2]| [stress * n[1], stress * n[0]];
        // Side tractions admit all three rigid impedance probes. Sealing every
        // exterior normal velocity would forbid the unit vertical response.
        let exterior: Vec<_> = (0..mesh.edges.len())
            .filter(|edge| mesh.is_boundary(*edge) && !edges.contains(edge))
            .map(|edge| {
                if mesh.edges[edge].normal[1] < -0.5 {
                    MaterialStokesBoundary::Velocity { edge, value: &field }
                } else {
                    MaterialStokesBoundary::TotalTraction { edge, value: &traction }
                }
            })
            .collect();
        let mut motion = planar_motion(
            &system, &fluid, &edges, field, [SPEED, 0.0, omega],
        );
        let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, _| {});
        let step = motion.step_with_heating(
            &mut thermal, binding(&ids), 0.05, &zero, &exterior, &[None],
            [DEPTH * MU * (SPEED + omega), 0.0], 0.0,
            PlanarStokesControls::default(), thermal_config(), || false,
        ).unwrap();
        let density_heat = 0.05 * MU * (SPEED + omega).powi(2);
        near(step.motion.velocity_m_s[0], SPEED, 2e-8);
        near(step.motion.angular_velocity_rad_s, omega, 2e-8);
        near(motion.state().rotation_rad(), 0.05 * omega, 2e-9);
        near(step.deposited_viscous_heat_j, DEPTH * density_heat, 2e-10);
        near(step.thermal.stored_energy_change_j, DEPTH * density_heat, 2e-10);
        near(step.thermal.residual_j, 0.0, 1e-12);
        for id in ids {
            near(thermal.temperature_k(id).unwrap(), 300.0 + density_heat / (RHO * CP), 2e-9);
        }
        assert_eq!(thermal.time_s(), motion.state().time_s());
    }
}

#[test]
fn g4_late_caloric_refusal_preserves_planar_state_and_short_retry_deposits_once() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let fluid = fluid();
    let law = Caloric { maximum_h: 0.0002 };
    let edges = top_edges(&mesh);
    let exterior = exterior(&mesh, &edges);
    let mut motion = planar_motion(&system, &fluid, &edges, shear, [SPEED, 0.0, 0.0]);
    let before = motion.state().clone();
    let (mut thermal, ids) = initial_thermal(&mesh, &law, false, |_, _| {});
    let original_cells = thermal.cells().clone();
    let refusal = shear_heat(&mut motion, &mut thermal, &ids, &exterior, 0.05, || false);
    assert!(matches!(refusal, Err(PlanarStokesHeatError::Thermal(
        DepositionError::MaterialDomain
    ))));
    assert_eq!(motion.state(), &before);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    shear_heat(&mut motion, &mut thermal, &ids, &exterior, 0.01, || false).unwrap();
    for id in ids {
        near(
            thermal.cells()[&id].specific_h_j_kg,
            0.01 * MU * SPEED * SPEED / RHO, 2e-10,
        );
    }
    assert_eq!(motion.state().time_s(), 0.01);
    assert_eq!(thermal.time_s(), 0.01);
}

#[test]
fn g4_final_cancellation_preserves_sources_events_and_both_planar_clocks() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let fluid = fluid();
    let law = Caloric { maximum_h: 100.0 };
    let edges = top_edges(&mesh);
    let exterior = exterior(&mesh, &edges);
    let make = || planar_motion(&system, &fluid, &edges, shear, [SPEED, 0.0, 0.0]);
    let mut control = make();
    let mut retry = make();
    let original_motion = retry.state().clone();
    let (original, ids) = initial_thermal(&mesh, &law, true, |i, cell| {
        if i == 0 {
            cell.specific_h_j_kg = 1.0;
            cell.boundary.source_w = 0.5;
        }
    });
    let mut control_thermal = original.clone();
    let mut retry_thermal = original.clone();
    let mut polls = 0;
    let control_step = shear_heat(
        &mut control, &mut control_thermal, &ids, &exterior, 0.05,
        || { polls += 1; false },
    ).unwrap();
    let mut attempts = 0;
    let refusal = shear_heat(
        &mut retry, &mut retry_thermal, &ids, &exterior, 0.05,
        || { attempts += 1; attempts == polls },
    );
    assert!(matches!(refusal, Err(PlanarStokesHeatError::Motion(
        PlanarStokesError::Cancelled
    ))));
    assert_eq!(retry.state(), &original_motion);
    assert_eq!(retry_thermal.cells(), original.cells());
    assert_eq!(retry_thermal.links(), original.links());
    assert_eq!(retry_thermal.time_s(), 0.0);
    let step = shear_heat(
        &mut retry, &mut retry_thermal, &ids, &exterior, 0.05, || false,
    ).unwrap();
    assert_eq!(retry.state(), control.state());
    assert_eq!(retry_thermal.cells(), control_thermal.cells());
    assert_eq!(retry_thermal.links(), control_thermal.links());
    assert_eq!(step.thermal, control_step.thermal);
    assert_eq!(step.thermal.events, 1);
    near(
        step.thermal.external_heat_j,
        0.02 * 0.5 + 0.03 * 0.8 + 0.05 * DEPTH * MU * SPEED * SPEED, 2e-10,
    );
    assert_eq!(retry_thermal.time_s(), retry.state().time_s());
}

#[test]
fn g1_planar_navier_slip_heat_is_partitioned_once_after_complete_route_admission() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let source = SingleRead(Cell::new(0));
    let law = Caloric { maximum_h: 100.0 };
    let edges = top_edges(&mesh);
    let gamma = SPEED / 2.0;
    let field = |x: [f64; 2]| [gamma * (1.0 + x[1]), 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [MU * gamma * n[1], MU * gamma * n[0]];
    let exterior = slip_exterior(&mesh, &edges, MU, &zero, &traction);
    let routes = slip_routes(&exterior, 0.25);
    let mut motion = planar_motion(&system, &source, &edges, field, [SPEED, 0.0, 0.0]);
    let original_motion = motion.state().clone();
    let (mut thermal, ids) = initial_thermal_with_events(
        &mesh, &law, false, |_, _| {},
        vec![TimedEvent { time_s: 0.0, event: Event::Birth(wall_cell(DEPTH)) }],
    );
    let original_cells = thermal.cells().clone();
    let queries = Cell::new(0);
    let force = |_| { queries.set(queries.get() + 1); [0.0; 2] };
    let refusal = motion.step_with_slip_heating(
        &mut thermal, binding(&ids), &[], 0.05, &force, &exterior, &[None],
        [DEPTH * MU * gamma, 0.0], 0.0, PlanarStokesControls::default(),
        thermal_config(), || false,
    );
    assert!(matches!(refusal, Err(PlanarStokesHeatError::Binding(_))));
    assert_eq!(queries.get(), 0, "missing heat ownership must refuse before flow work");
    assert_eq!(motion.state(), &original_motion);
    assert_eq!(thermal.cells(), &original_cells);
    let step = motion.step_with_slip_heating(
        &mut thermal, binding(&ids), &routes, 0.05, &zero, &exterior, &[None],
        [DEPTH * MU * gamma, 0.0], 0.0, PlanarStokesControls::default(),
        thermal_config(), || false,
    ).unwrap();
    let bulk = 0.05 * DEPTH * MU * gamma * gamma;
    let slip = 0.05 * DEPTH * MU * (SPEED - gamma).powi(2);
    near(step.deposited_viscous_heat_j, bulk, 2e-10);
    near(step.deposited_slip_heat_j, slip, 2e-10);
    near(step.thermal.stored_energy_change_j, bulk + slip, 2e-10);
    near(step.thermal.external_heat_j, bulk + slip, 2e-10);
    near(
        thermal.cells()[&WALL].specific_h_j_kg * thermal.cells()[&WALL].mass_kg,
        0.75 * slip, 2e-10,
    );
    let fluid_heat: f64 = ids.iter().map(|id| {
        let cell = &thermal.cells()[id];
        cell.mass_kg * cell.specific_h_j_kg
    }).sum();
    near(fluid_heat, bulk + 0.25 * slip, 2e-10);
    assert_eq!(source.0.get(), 1, "thermalization must not resolve material twice");
    assert_eq!(thermal.time_s(), motion.state().time_s());
}


#[test]
fn observed_planar_heat_preserves_rotation_slip_and_atomic_observer_history() {
    let law = Caloric { maximum_h: 100.0 };
    for with_slip in [false, true] {
        let n = if with_slip { 1 } else { 2 };
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let fluid = fluid();
        let edges = top_edges(&mesh);
        let gamma = if with_slip { SPEED / 2.0 } else { SPEED };
        let omega = if with_slip { 0.0 } else { 0.12 };
        let slip = SPEED - gamma;
        let field = |x: [f64; 2]| [slip + gamma * x[1], omega * (x[0] - 0.5)];
        let stress = MU * (gamma + omega);
        let traction = |_: [f64; 2], normal: [f64; 2]| [stress * normal[1], stress * normal[0]];
        let boundaries: Vec<_> = if with_slip {
            slip_exterior(&mesh, &edges, MU, &zero, &traction)
        } else {
            (0..mesh.edges.len())
                .filter(|edge| mesh.is_boundary(*edge) && !edges.contains(edge))
                .map(|edge| {
                    if mesh.edges[edge].normal[1] < -0.5 {
                        MaterialStokesBoundary::Velocity { edge, value: &field }
                    } else {
                        MaterialStokesBoundary::TotalTraction { edge, value: &traction }
                    }
                })
                .collect()
        };
        let routes = slip_routes(&boundaries, 0.25);
        let pressure = [None];
        let (original, ids) = initial_thermal_with_events(
            &mesh, &law, true,
            |i, cell| { if i == 0 { cell.specific_h_j_kg = 1.0; } },
            vec![TimedEvent { time_s: 0.0, event: Event::Birth(wall_cell(DEPTH)) }],
        );
        let advance = |
            owner: &mut PlanarStokesMotion<'_, '_, '_, FluidTransportState>,
            thermal: &mut DepositionSimulation<'_>,
            observer: &mut ClockedHeat,
            cancelled: &mut dyn FnMut() -> bool,
        | {
            if with_slip {
                owner.step_with_slip_heating_observed(
                    thermal, observer, binding(&ids), &routes, 0.05, &zero, &boundaries,
                    &pressure, [DEPTH * stress, 0.0], 0.0, PlanarStokesControls::default(),
                    thermal_config(), cancelled,
                )
            } else {
                owner.step_with_heating_observed(
                    thermal, observer, binding(&ids), 0.05, &zero, &boundaries,
                    &pressure, [DEPTH * stress, 0.0], 0.0, PlanarStokesControls::default(),
                    thermal_config(), cancelled,
                )
            }
        };
        let mut owner = planar_motion(&system, &fluid, &edges, field, [SPEED, 0.0, omega]);
        let mut thermal = original.clone();
        let mut observer = ClockedHeat::new(ids[0]);
        let mut polls = 0;
        let step = advance(&mut owner, &mut thermal, &mut observer, &mut || {
            polls += 1;
            false
        }).unwrap();
        let bulk_heat = 0.05 * DEPTH * MU * (gamma + omega).powi(2);
        let slip_heat = 0.05 * DEPTH * MU * slip * slip;
        let observed_heat = 0.015 * 0.01 + 0.035 * 0.03 + 0.2 * 0.05 * 0.05;
        let total_heat = bulk_heat + slip_heat + observed_heat + 0.03 * 0.8;
        near(step.motion.velocity_m_s[0], SPEED, 2e-8);
        near(step.motion.angular_velocity_rad_s, omega, 2e-8);
        near(owner.state().rotation_rad(), 0.05 * omega, 2e-9);
        near(step.deposited_viscous_heat_j, bulk_heat, 2e-10);
        near(step.deposited_slip_heat_j, slip_heat, 2e-10);
        near(observer.heat_j, observed_heat, 1e-14);
        near(step.thermal.external_heat_j, total_heat, 2e-10);
        near(step.thermal.stored_energy_change_j, total_heat, 2e-10);
        near(step.thermal.residual_j, 0.0, 1e-12);
        near(
            thermal.cells()[&WALL].specific_h_j_kg,
            0.75 * slip_heat / thermal.cells()[&WALL].mass_kg,
            2e-9,
        );
        assert_eq!(observer.begins, 1);
        assert_eq!(observer.events, 1);
        assert_eq!(step.thermal.events, 1);
        assert_eq!(step.thermal.substeps, observer.starts_and_durations.len());
        assert!(observer.starts_and_durations.iter().any(|&(start, _)| start == 0.015));
        assert!(observer.starts_and_durations.iter().any(|&(start, _)| start == 0.02));
        assert_eq!(thermal.links(), original.links());
        assert_eq!(thermal.time_s(), 0.05);
        assert_eq!(owner.state().time_s(), thermal.time_s());

        for refuse in [true, false] {
            let mut retry = planar_motion(&system, &fluid, &edges, field, [SPEED, 0.0, omega]);
            let old_motion = retry.state().clone();
            let mut retry_thermal = original.clone();
            let mut retry_observer = ClockedHeat::new(ids[0]);
            retry_observer.refuse_at_s = refuse.then_some(0.03);
            let old_observer = retry_observer.clone();
            let mut attempts = 0;
            let refusal = advance(
                &mut retry, &mut retry_thermal, &mut retry_observer, &mut || {
                    attempts += 1;
                    !refuse && attempts == polls
                },
            );
            if refuse {
                assert!(matches!(refusal, Err(PlanarStokesObservedHeatError::Observer(
                    CoupledHeatRefusal::Endpoint
                ))));
            } else {
                assert!(matches!(refusal, Err(PlanarStokesObservedHeatError::Heat(
                    PlanarStokesHeatError::Motion(PlanarStokesError::Cancelled)
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
            assert_eq!(retry_thermal.time_s(), thermal.time_s());
            assert_eq!(retry_thermal.cells(), thermal.cells());
            assert_eq!(retry_observer, observer);
        }
    }
}
