//! G1/G3/G4: native pressure-driven body motion, added mass and atomic refusal.

use fs_couple::stokes_motion::{
    GuidedStokesConfig, GuidedStokesControls, GuidedStokesError, GuidedStokesMotion,
};
use fs_flux::{FluxSystem, MaterialFluxNumerics, MaterialStokesBoundary, TriMesh};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{MassProperties, Pose, RigidBodyState, Vec3};
use fs_solid::Mesh2;
use std::cell::Cell;

#[path = "stokes_motion/planar.rs"]
mod planar;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}"
    );
}
fn material(density: f64, viscosity: f64) -> FluidTransportState {
    // Explicit numerical material fixture; no sourced-liquid authority inferred.
    FluidTransportState::try_new(density, 1000.0, 0.5, viscosity).unwrap()
}
fn config() -> GuidedStokesConfig {
    GuidedStokesConfig {
        direction: [1.0, 0.0],
        depth_m: 0.3,
        max_displacement_m: 0.05,
        stiffness_n_m: 7.0,
        damping_n_s_m: 0.4,
        velocity_tolerance_m_s: 1e-8,
    }
}
fn mass() -> MassProperties {
    MassProperties::new(0.7, Vec3::ZERO, Vec3::new(0.1, 0.1, 0.1)).unwrap()
}
fn body(speed: f64) -> RigidBodyState {
    RigidBodyState::new(
        Pose::identity(),
        Vec3::new(mass().mass() * speed, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap()
}
fn zero(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn boundary_edges(mesh: &TriMesh) -> Vec<usize> {
    (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect()
}
fn uniform_history(mesh: &TriMesh, speed: f64) -> Vec<f64> {
    let mut values = vec![0.0; 2 * mesh.edges.len()];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        values[2 * edge] = speed * geometry.normal[0];
    }
    values
}

#[test]
fn closed_translating_box_resolves_fluid_added_mass_and_complete_be_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.2, 0.8, 3, 2));
    let system = FluxSystem::new(&mesh);
    let area: f64 = mesh.areas.iter().sum();
    let edges = boundary_edges(&mesh);
    let dt = 0.002;
    let applied = 1.3;
    let initial_speed = 0.03;
    for (density, depth) in [(2.0, 0.3), (8.0, 0.6)] {
        let fluid = material(density, 0.8);
        let cfg = GuidedStokesConfig {
            depth_m: depth,
            ..config()
        };
        let mut motion = GuidedStokesMotion::new(
            &system,
            &fluid,
            MaterialFluxNumerics::default(),
            mass(),
            body(initial_speed),
            uniform_history(&mesh, initial_speed),
            &edges,
            cfg,
        )
        .unwrap();
        let fluid_mass = density * area * depth;
        let total_mass = mass().mass() + fluid_mass;
        let mut v0 = initial_speed;
        let mut q0 = 0.0;
        for _ in 0..4 {
            let expected = (total_mass * v0 + dt * (applied - cfg.stiffness_n_m * q0))
                / (total_mass + cfg.damping_n_s_m * dt + cfg.stiffness_n_m * dt * dt);
            let momentum_before = motion.state().body().linear_momentum_world();
            let step = motion
                .step(
                    dt,
                    &zero,
                    &[],
                    &[Some(91.0)],
                    applied,
                    GuidedStokesControls::default(),
                    || false,
                )
                .unwrap();
            // The conservation gate must describe the ACTUALLY stored body
            // momentum. With nonbinary mass and nonzero initial motion, m*dv
            // can round differently from this independently observed change.
            let inertia_force = step
                .body_after
                .linear_momentum_world()
                .sub(momentum_before)
                .dot(Vec3::new(cfg.direction[0], cfg.direction[1], 0.0))
                / dt;
            let force_residual = inertia_force
                - (step.fluid_force_n + applied
                    - cfg.stiffness_n_m * motion.state().displacement_m()
                    - cfg.damping_n_s_m * step.speed_m_s);
            assert_eq!(step.force_residual_n, force_residual);
            near(step.speed_m_s, expected, 3e-9);
            near(step.fluid_force_n, -fluid_mass * (expected - v0) / dt, 3e-7);
            near(step.fluid_force_slope_n_s_m, -fluid_mass / dt, 3e-6);
            near(motion.state().displacement_m(), q0 + dt * expected, 3e-11);
            let (velocity_error, divergence) =
                system.velocity_error(&step.flow.discrete().x, &|_| [expected, 0.0]);
            assert!(velocity_error < 3e-9 && divergence < 3e-9);
            // Independent full-box identity: fluid moves uniformly, so physical
            // viscous loss vanishes; its BE loss is known from its actual mass.
            let q1 = motion.state().displacement_m();
            let dv = step.speed_m_s - v0;
            let dx = q1 - q0;
            let storage = 0.5 * total_mass * (step.speed_m_s + v0) * dv
                + 0.5 * cfg.stiffness_n_m * (q1 + q0) * dx;
            let physical_loss = cfg.damping_n_s_m * dt * step.speed_m_s.powi(2);
            let numerical_loss = 0.5 * total_mass * dv * dv + 0.5 * cfg.stiffness_n_m * dx * dx;
            near(
                storage + physical_loss + numerical_loss,
                applied * dx,
                3e-10,
            );
            near(
                step.fluid_kinetic_change_j,
                0.5 * fluid_mass * (step.speed_m_s + v0) * dv,
                3e-10,
            );
            near(step.mechanical_work_residual_j, 0.0, 2e-10);
            v0 = step.speed_m_s;
            q0 = q1;
        }
    }
}

#[test]
fn open_pressure_drives_body_with_the_original_traction_datum() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let fluid = material(3.0, 1.1);
    let cfg = GuidedStokesConfig {
        stiffness_n_m: 0.0,
        damping_n_s_m: 0.0,
        ..config()
    };
    let edges: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|&edge| mesh.edges[edge].normal[0] < -0.5)
        .collect();
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -9.0 * n);
    let exterior: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|edge| !edges.contains(edge))
        .map(|edge| MaterialStokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let mut motion = GuidedStokesMotion::new(
        &system,
        &fluid,
        MaterialFluxNumerics::default(),
        mass(),
        body(0.0),
        vec![0.0; system.n_u],
        &edges,
        cfg,
    )
    .unwrap();
    let step = motion
        .step(
            0.001,
            &zero,
            &exterior,
            &[None],
            0.0,
            GuidedStokesControls::default(),
            || false,
        )
        .unwrap();
    assert!(step.speed_m_s < 0.0 && step.fluid_force_n < 0.0);
    near(
        mass().mass() * step.speed_m_s / 0.001,
        step.fluid_force_n,
        2e-8,
    );
    assert!(step.interface_work_j > 0.0);
    assert_eq!(step.boundary_load.pressure_offsets_pa, vec![(0, None)]);
    near(
        step.body_kinetic_change_j + step.mechanical_time_discretization_loss_j,
        step.interface_work_j,
        2e-10,
    );
    let before = motion.state().clone();
    assert!(
        motion
            .step(
                0.001,
                &zero,
                &exterior,
                &[Some(0.0)],
                0.0,
                GuidedStokesControls::default(),
                || false
            )
            .is_err()
    );
    assert_eq!(
        motion.state(),
        &before,
        "traction datum must not acquire a second pressure offset"
    );
}

struct CountedFluid {
    state: FluidTransportState,
    queries: Cell<usize>,
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn viscous_lid_motion_reacts_back_and_late_cancellation_preserves_continuation() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let fluid = CountedFluid {
        state: material(2.0, 0.5),
        queries: Cell::new(0),
    };
    let cfg = GuidedStokesConfig {
        stiffness_n_m: 0.0,
        damping_n_s_m: 0.0,
        ..config()
    };
    let edges: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|&edge| mesh.edges[edge].normal[1] > 0.5)
        .collect();
    let exterior: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|edge| !edges.contains(edge))
        .map(|edge| MaterialStokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let make = || {
        GuidedStokesMotion::new(
            &system,
            &fluid,
            MaterialFluxNumerics::default(),
            mass(),
            body(0.1),
            vec![0.0; system.n_u],
            &edges,
            cfg,
        )
        .unwrap()
    };
    let mut motion = make();
    let mut reference = make();
    let queries = fluid.queries.get();
    let first = motion
        .step(
            0.002,
            &zero,
            &exterior,
            &[Some(4.0)],
            0.0,
            GuidedStokesControls::default(),
            || false,
        )
        .unwrap();
    reference
        .step(
            0.002,
            &zero,
            &exterior,
            &[Some(4.0)],
            0.0,
            GuidedStokesControls::default(),
            || false,
        )
        .unwrap();
    assert!(first.speed_m_s > 0.0 && first.speed_m_s < 0.1);
    assert!(first.fluid_force_n < 0.0 && first.flow.kinetic_energy_j_per_m() > 0.0);
    near(first.guide_force_n[0], 0.0, 1e-13);
    near(
        first.guide_force_n[1],
        -cfg.depth_m * first.boundary_load.stress.force_n_per_m[1],
        1e-13,
    );
    let before = motion.state().clone();
    let mut polls = 0;
    let refused = motion.step(
        0.002,
        &zero,
        &exterior,
        &[Some(4.0)],
        0.0,
        GuidedStokesControls::default(),
        || {
            polls += 1;
            polls == 5
        },
    );
    assert!(matches!(refused, Err(GuidedStokesError::Cancelled)));
    assert_eq!(
        polls, 5,
        "cancellation must reach the last prepublication boundary"
    );
    assert_eq!(motion.state(), &before);
    let next = motion
        .step(
            0.002,
            &zero,
            &exterior,
            &[Some(4.0)],
            0.0,
            GuidedStokesControls::default(),
            || false,
        )
        .unwrap();
    let expected = reference
        .step(
            0.002,
            &zero,
            &exterior,
            &[Some(4.0)],
            0.0,
            GuidedStokesControls::default(),
            || false,
        )
        .unwrap();
    assert_eq!(motion.state(), reference.state());
    assert_eq!(next.flow.discrete().x, expected.flow.discrete().x);
    assert_eq!(
        fluid.queries.get(),
        queries,
        "source transport is bound once, not re-queried during Schur trials"
    );
}

#[test]
fn geometry_bound_and_guide_admission_refuse_without_partial_motion() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let fluid = material(1.0, 0.4);
    let edges = boundary_edges(&mesh);
    let cfg = GuidedStokesConfig {
        max_displacement_m: 1e-5,
        ..config()
    };
    let make = || {
        GuidedStokesMotion::new(
            &system,
            &fluid,
            MaterialFluxNumerics::default(),
            mass(),
            body(0.0),
            vec![0.0; system.n_u],
            &edges,
            cfg,
        )
        .unwrap()
    };
    let mut motion = make();
    let before = motion.state().clone();
    assert!(matches!(
        motion.step(
            0.01,
            &zero,
            &[],
            &[Some(0.0)],
            100.0,
            GuidedStokesControls::default(),
            || false
        ),
        Err(GuidedStokesError::DisplacementLimit { .. })
    ));
    assert_eq!(motion.state(), &before);
    let mut reference = make();
    motion
        .step(
            0.001,
            &zero,
            &[],
            &[Some(100.0)],
            0.1,
            GuidedStokesControls::default(),
            || false,
        )
        .unwrap();
    reference
        .step(
            0.001,
            &zero,
            &[],
            &[Some(0.0)],
            0.1,
            GuidedStokesControls::default(),
            || false,
        )
        .unwrap();
    near(
        motion.state().displacement_m(),
        reference.state().displacement_m(),
        1e-12,
    );
    near(
        motion.state().body().linear_momentum_world().x,
        reference.state().body().linear_momentum_world().x,
        1e-10,
    );
    let bad = RigidBodyState::new(Pose::identity(), Vec3::new(0.0, 1.0, 0.0), Vec3::ZERO).unwrap();
    assert!(
        GuidedStokesMotion::new(
            &system,
            &fluid,
            MaterialFluxNumerics::default(),
            mass(),
            bad,
            vec![0.0; system.n_u],
            &edges,
            cfg
        )
        .is_err()
    );
    let rotating =
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)).unwrap();
    assert!(
        GuidedStokesMotion::new(
            &system,
            &fluid,
            MaterialFluxNumerics::default(),
            mass(),
            rotating,
            vec![0.0; system.n_u],
            &edges,
            cfg
        )
        .is_err()
    );
    assert!(
        GuidedStokesMotion::new(
            &system,
            &fluid,
            MaterialFluxNumerics::default(),
            mass(),
            body(0.1),
            vec![0.0; system.n_u],
            &edges,
            cfg
        )
        .is_err()
    );
}
