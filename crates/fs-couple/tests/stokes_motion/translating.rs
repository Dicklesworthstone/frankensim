//! Original nonlinear flow and native mechanics over finite world translations.

use super::{CountedFluid, body, boundary_edges, mass, material, near, uniform_history, zero};
use fs_couple::stokes_motion::translating::{
    GuidedTranslationConfig, GuidedTranslationControls, GuidedTranslationError,
    GuidedTranslationMotion,
};
use fs_flux::{
    FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary, NavierStokesControls, TriMesh,
};
use fs_solid::Mesh2;
use std::cell::Cell;

fn config() -> GuidedTranslationConfig {
    GuidedTranslationConfig {
        direction: [1.0, 0.0],
        depth_m: 0.3,
        stiffness_n_m: 0.0,
        damping_n_s_m: 0.0,
        velocity_tolerance_m_s: 1e-8,
    }
}
fn controls() -> GuidedTranslationControls {
    GuidedTranslationControls {
        speed_bracket_m_s: [-3.0, 3.0],
        maximum_flow_solves: 32,
        force_absolute_n: 2e-8,
        force_relative: 1e-8,
        work_absolute_j: 2e-9,
        work_relative: 1e-8,
    }
}
fn fluid_controls() -> NavierStokesControls {
    NavierStokesControls {
        relative_tolerance: 1e-10,
        ..NavierStokesControls::default()
    }
}

#[test]
fn finite_container_translation_resolves_fluid_inertia_and_total_momentum() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.2, 0.8, 3, 2));
    let system = FluxSystem::new(&mesh);
    let edges = boundary_edges(&mesh);
    let area: f64 = mesh.areas.iter().sum();
    let dt = 0.2;
    for (density, depth, sign, direction) in
        [(2.0, 0.3, 1.0, [1.0, 0.0]), (5.0, 0.6, -1.0, [0.6, 0.8])]
    {
        let applied = 0.9 * sign;
        let source = CountedFluid {
            state: material(density, 0.8),
            queries: Cell::new(0),
        };
        let cfg = GuidedTranslationConfig {
            direction,
            depth_m: depth,
            ..config()
        };
        let initial = 0.7 * sign;
        let velocity = direction.map(|d| d * initial);
        let initial_body = fs_mbd::RigidBodyState::new(
            fs_mbd::Pose::identity(),
            fs_mbd::Vec3::new(
                mass().mass() * velocity[0],
                mass().mass() * velocity[1],
                0.0,
            ),
            fs_mbd::Vec3::ZERO,
        )
        .unwrap();
        let mut history = vec![0.0; system.n_u];
        for (edge, geometry) in mesh.edges.iter().enumerate() {
            history[2 * edge] = velocity[0] * geometry.normal[0] + velocity[1] * geometry.normal[1];
        }
        let mut owner = GuidedTranslationMotion::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            mass(),
            initial_body,
            &history,
            &edges,
            cfg,
        )
        .unwrap();
        let queries = source.queries.get();
        let fluid_mass = density * area * depth;
        let total_mass = mass().mass() + fluid_mass;
        let mut speed = initial;
        let mut position = 0.0;
        let mut previous_total_momentum = direction.map(|d| d * total_mass * initial);
        for _ in 0..10 {
            let expected = speed + dt * applied / total_mass;
            let step = owner
                .step(
                    dt,
                    &zero,
                    &[],
                    &[Some(91.0)],
                    applied,
                    fluid_controls(),
                    controls(),
                    || false,
                )
                .unwrap();
            position += dt * expected;
            near(step.speed_m_s, expected, 4e-8);
            let center = owner.body().pose().position_world();
            for i in 0..2 {
                near(
                    owner.flow_state().translation_m()[i],
                    direction[i] * position,
                    5e-8,
                );
                near([center.x, center.y][i], direction[i] * position, 5e-8);
            }
            near(
                step.fluid_force_n,
                -fluid_mass * (expected - speed) / dt,
                4e-7,
            );
            let body_momentum = owner.body().linear_momentum_world();
            let momentum: [f64; 2] = core::array::from_fn(|i| {
                [body_momentum.x, body_momentum.y][i]
                    + depth * step.flow.world_momentum_n_s_per_m[i]
            });
            for i in 0..2 {
                near(
                    momentum[i] - previous_total_momentum[i],
                    direction[i] * applied * dt,
                    3e-8,
                );
                near(step.guide_force_n[i], 0.0, 4e-7);
            }
            near(
                step.guide_torque_n_m,
                -(0.6 * direction[1] - 0.4 * direction[0]) * step.fluid_force_n,
                4e-7,
            );
            let (error, divergence) = system
                .velocity_error(&step.flow.world_velocity_dofs, &|_| {
                    direction.map(|d| d * expected)
                });
            assert!(error < 4e-8 && divergence < 4e-8);
            // Independent complete-box energy identity: uniform fluid has no
            // viscous heat; both finite masses retain their own BE inertia loss.
            let dv = step.speed_m_s - speed;
            let storage = 0.5 * total_mass * (step.speed_m_s + speed) * dv;
            let time_loss = 0.5 * total_mass * dv * dv;
            near(storage + time_loss, applied * dt * step.speed_m_s, 4e-8);
            near(step.mechanical_work_residual_j, 0.0, 3e-9);
            assert!(
                step.flow_solves <= 4,
                "linear inertia must not require many fluid trials"
            );
            speed = step.speed_m_s;
            previous_total_momentum = momentum;
        }
        assert!(
            position.abs() > 1.2,
            "the mesh must actually move farther than its width"
        );
        assert_eq!(
            source.queries.get(),
            queries,
            "one immutable material owner throughout all trials"
        );
    }
}

#[test]
fn endpoint_world_forcing_and_spring_work_follow_the_moving_domain() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.2, 0.8, 3, 2));
    let system = FluxSystem::new(&mesh);
    let edges = boundary_edges(&mesh);
    let source = material(2.0, 0.8);
    let cfg = GuidedTranslationConfig {
        stiffness_n_m: 0.4,
        damping_n_s_m: 0.2,
        ..config()
    };
    let volume = 1.2 * 0.8 * cfg.depth_m;
    let total_mass = mass().mass() + 2.0 * volume;
    let mut owner = GuidedTranslationMotion::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        mass(),
        body(0.3),
        &uniform_history(&mesh, 0.3),
        &edges,
        cfg,
    )
    .unwrap();
    let dt = 0.15;
    let applied = 0.4;
    let intercept = 0.7;
    let gradient = 0.3;
    let force = |world: [f64; 2]| [intercept + gradient * world[0], 0.0];
    let mut speed = 0.3;
    let mut position = 0.0;
    for _ in 0..5 {
        // Integrate the original linear volumetric force over the ACTUAL
        // translated rectangle. The implicit force depends on X1=X0+dt*v1.
        let expected = (total_mass * speed
            + dt * (applied + volume * (intercept + gradient * (0.6 + position))
                - cfg.stiffness_n_m * position))
            / (total_mass
                + cfg.damping_n_s_m * dt
                + (cfg.stiffness_n_m - volume * gradient) * dt * dt);
        let step = owner
            .step(
                dt,
                &force,
                &[],
                &[Some(0.0)],
                applied,
                fluid_controls(),
                controls(),
                || false,
            )
            .unwrap();
        let next = position + dt * expected;
        near(step.speed_m_s, expected, 5e-8);
        near(step.displacement_m, next, 5e-8);
        near(
            step.spring_energy_change_j,
            0.5 * cfg.stiffness_n_m * (next * next - position * position),
            2e-8,
        );
        near(
            step.damping_loss_j,
            cfg.damping_n_s_m * dt * expected * expected,
            2e-8,
        );
        near(
            step.fluid_force_n,
            volume * (intercept + gradient * (0.6 + next)) - 2.0 * volume * (expected - speed) / dt,
            5e-7,
        );
        near(step.mechanical_work_residual_j, 0.0, 3e-9);
        position = step.displacement_m;
        speed = step.speed_m_s;
    }
}

#[test]
fn nonlinear_viscous_flow_reacts_on_the_continuing_body() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let source = material(2.0, 0.3);
    let edges: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|&e| mesh.edges[e].normal[1] > 0.5)
        .collect();
    // The outer translating control surface has zero world velocity data;
    // left/right crossing fluxes are physical and use relative advection.
    let exterior: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|e| !edges.contains(e))
        .map(|edge| MaterialNavierStokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let mut owner = GuidedTranslationMotion::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        mass(),
        body(0.1),
        &vec![0.0; system.n_u],
        &edges,
        config(),
    )
    .unwrap();
    let mut previous = 0.1;
    let mut limits = controls();
    limits.speed_bracket_m_s = [0.0, 0.15];
    for _ in 0..3 {
        let step = owner
            .step(
                0.002,
                &zero,
                &exterior,
                &[Some(0.0)],
                0.0,
                fluid_controls(),
                limits,
                || false,
            )
            .unwrap();
        assert!(step.speed_m_s > 0.0 && step.speed_m_s < previous);
        assert!(step.fluid_force_n < 0.0 && step.flow.world_kinetic_energy_j_per_m > 0.0);
        assert!(step.flow.relative_flow.nonlinear_residual <= 1e-9);
        assert!(step.flow.relative_flow.iterations > 1);
        near(
            step.momentum.linear_increment_world.x / 0.002,
            step.fluid_force_n,
            3e-8,
        );
        near(
            step.guide_force_n[1],
            -config().depth_m * step.boundary_load.stress.force_n_per_m[1],
            1e-12,
        );
        near(
            step.guide_torque_n_m,
            -config().depth_m * step.boundary_load.stress.moment_n,
            1e-12,
        );
        near(step.mechanical_work_residual_j, 0.0, 3e-9);
        previous = step.speed_m_s;
    }
}

#[test]
fn rejected_root_and_late_cancellation_keep_both_original_histories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let source = material(2.0, 0.5);
    let edges = boundary_edges(&mesh);
    // A tiny RELATIVE transverse momentum can still violate the authored
    // absolute guide-speed tolerance by many orders at a large axial speed.
    let off_guide = fs_mbd::RigidBodyState::new(
        fs_mbd::Pose::identity(),
        fs_mbd::Vec3::new(mass().mass() * 1e12, mass().mass() * 1e-3, 0.0),
        fs_mbd::Vec3::ZERO,
    )
    .unwrap();
    assert!(
        GuidedTranslationMotion::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            mass(),
            off_guide,
            &uniform_history(&mesh, 1e12),
            &edges,
            config(),
        )
        .is_err()
    );
    let make = || {
        GuidedTranslationMotion::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            mass(),
            body(0.0),
            &vec![0.0; system.n_u],
            &edges,
            config(),
        )
        .unwrap()
    };
    let mut owner = make();
    let mut reference = make();
    for motion in [&mut owner, &mut reference] {
        motion
            .step(
                0.02,
                &zero,
                &[],
                &[Some(0.0)],
                0.7,
                fluid_controls(),
                controls(),
                || false,
            )
            .unwrap();
    }
    let before = (owner.body(), owner.flow_state().clone());
    let mut limited = controls();
    limited.maximum_flow_solves = 1;
    assert!(matches!(
        owner.step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            0.7,
            fluid_controls(),
            limited,
            || false
        ),
        Err(GuidedTranslationError::SolveLimit { .. })
    ));
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    let mut unbracketed = controls();
    unbracketed.speed_bracket_m_s = [1.0, 2.0];
    assert!(matches!(
        owner.step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            0.7,
            fluid_controls(),
            unbracketed,
            || false
        ),
        Err(GuidedTranslationError::Unbracketed { .. })
    ));
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    assert!(
        owner
            .step(
                0.02,
                &zero,
                &[],
                &[None],
                0.7,
                fluid_controls(),
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    let mut polls = 0;
    let expected = reference
        .step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            0.7,
            fluid_controls(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        owner.step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            0.7,
            fluid_controls(),
            controls(),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(GuidedTranslationError::Cancelled)
    ));
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    let actual = owner
        .step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            0.7,
            fluid_controls(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(owner.body(), reference.body());
    assert_eq!(owner.flow_state(), reference.flow_state());
    assert_eq!(
        actual.flow.world_velocity_dofs,
        expected.flow.world_velocity_dofs
    );
}
