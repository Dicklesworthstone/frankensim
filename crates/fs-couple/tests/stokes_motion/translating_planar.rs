//! Finite two-component ALE/body dynamics; authored numerical material fixtures.

use super::{CountedFluid, boundary_edges, mass, material, near, zero};
use fs_couple::stokes_motion::translating::planar::{
    PlanarTranslationConfig, PlanarTranslationControls, PlanarTranslationError,
    PlanarTranslationMotion,
};
use fs_flux::{
    FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary, NavierStokesControls,
    OpenBoundaryPolicy, TriMesh,
};
use fs_mbd::{Pose, RigidBodyState, Vec3};
use fs_solid::Mesh2;
use std::cell::Cell;

fn config() -> PlanarTranslationConfig {
    PlanarTranslationConfig {
        depth_m: 0.3,
        stiffness_n_m: [0.0; 2],
        damping_n_s_m: [0.0; 2],
        velocity_tolerance_m_s: 1e-8,
    }
}
fn controls() -> PlanarTranslationControls {
    PlanarTranslationControls {
        velocity_bounds_m_s: [[-32.0, 32.0]; 2],
        jacobian_step_m_s: [0.2, 0.2],
        maximum_flow_solves: 32,
        force_absolute_n: 3e-8,
        force_relative: 1e-8,
        work_absolute_j: 3e-9,
        work_relative: 1e-8,
    }
}
fn fluid_controls() -> NavierStokesControls {
    NavierStokesControls {
        relative_tolerance: 1e-10,
        ..NavierStokesControls::default()
    }
}
fn native_body(velocity: [f64; 2]) -> RigidBodyState {
    RigidBodyState::new(
        Pose::identity(),
        Vec3::new(
            mass().mass() * velocity[0],
            mass().mass() * velocity[1],
            0.0,
        ),
        Vec3::ZERO,
    )
    .unwrap()
}
fn history(mesh: &TriMesh, velocity: [f64; 2]) -> Vec<f64> {
    let mut values = vec![0.0; 2 * mesh.edges.len()];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        values[2 * edge] = velocity[0] * geometry.normal[0] + velocity[1] * geometry.normal[1];
    }
    values
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[test]
fn freely_turning_translation_resolves_both_fluid_inertias_and_combined_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.2, 0.8, 3, 2));
    let system = FluxSystem::new(&mesh);
    let edges = boundary_edges(&mesh);
    let area: f64 = mesh.areas.iter().sum();
    let dt = 0.2;
    for (density, depth, boost) in [(2.0, 0.3, [0.0, 0.0]), (5.0, 0.6, [8.0, -4.0])] {
        let source = CountedFluid {
            state: material(density, 0.8),
            queries: Cell::new(0),
        };
        let cfg = PlanarTranslationConfig {
            depth_m: depth,
            ..config()
        };
        let initial = [0.6 + boost[0], -0.7 + boost[1]];
        let mut owner = PlanarTranslationMotion::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            mass(),
            native_body(initial),
            &history(&mesh, initial),
            &edges,
            cfg,
        )
        .unwrap();
        let queries = source.queries.get();
        let fluid_mass = density * area * depth;
        let total_mass = mass().mass() + fluid_mass;
        let mut velocity = initial;
        let mut translation = [0.0; 2];
        let mut total_momentum = velocity.map(|v| total_mass * v);
        for n in 0..10 {
            // The applied direction changes independently of the current
            // velocity. A single translation guide cannot follow this path.
            let applied = [0.4 + 0.1 * n as f64, -0.3 + 0.07 * n as f64];
            let expected: [f64; 2] =
                core::array::from_fn(|i| velocity[i] + dt * applied[i] / total_mass);
            let old_body = owner.body();
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
            let native_increment = step
                .body_after
                .linear_momentum_world()
                .sub(old_body.linear_momentum_world());
            let p = owner.body().linear_momentum_world();
            let combined = [
                p.x + depth * step.flow.world_momentum_n_s_per_m[0],
                p.y + depth * step.flow.world_momentum_n_s_per_m[1],
            ];
            for i in 0..2 {
                translation[i] += dt * expected[i];
                near(step.velocity_m_s[i], expected[i], 8e-8);
                near(step.translation_m[i], translation[i], 2e-7);
                near(owner.flow_state().translation_m()[i], translation[i], 2e-7);
                near(
                    step.fluid_force_n[i],
                    -fluid_mass * (expected[i] - velocity[i]) / dt,
                    6e-7,
                );
                near(combined[i] - total_momentum[i], applied[i] * dt, 8e-8);
                near(
                    [native_increment.x, native_increment.y][i] / dt,
                    step.fluid_force_n[i] + applied[i],
                    5e-8,
                );
            }
            let center = owner.body().pose().position_world();
            near(center.x, translation[0], 2e-7);
            near(center.y, translation[1], 2e-7);
            assert_eq!(
                owner.body().pose().orientation(),
                old_body.pose().orientation()
            );
            assert_eq!(owner.body().angular_momentum_body(), Vec3::ZERO);
            near(
                step.restraint_torque_n_m,
                -(0.6 * step.fluid_force_n[1] - 0.4 * step.fluid_force_n[0]),
                6e-7,
            );
            let (error, divergence) =
                system.velocity_error(&step.flow.world_velocity_dofs, &|_| expected);
            assert!(error < 2e-7 && divergence < 2e-7);
            let dv = core::array::from_fn(|i| step.velocity_m_s[i] - velocity[i]);
            let mean_twice = core::array::from_fn(|i| step.velocity_m_s[i] + velocity[i]);
            let storage = 0.5 * total_mass * dot(mean_twice, dv);
            let numerical = 0.5 * total_mass * dot(dv, dv);
            near(
                storage + numerical,
                dt * dot(applied, step.velocity_m_s),
                2e-7,
            );
            near(step.mechanical_work_residual_j, 0.0, 2e-8);
            near(owner.flow_state().time_s(), (n + 1) as f64 * dt, 1e-14);
            assert!(
                step.flow_solves <= 4,
                "an affine coupled endpoint needs at most four native trials"
            );
            velocity = step.velocity_m_s;
            total_momentum = combined;
        }
        assert!(translation[0].hypot(translation[1]) > 1.2);
        assert_eq!(
            source.queries.get(),
            queries,
            "trials must retain the original material binding"
        );
    }
}

#[test]
fn off_diagonal_world_force_closes_one_coupled_endpoint_with_springs_and_dashpots() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.2, 0.8, 3, 2));
    let system = FluxSystem::new(&mesh);
    let edges = boundary_edges(&mesh);
    let source = material(2.0, 0.8);
    let cfg = PlanarTranslationConfig {
        stiffness_n_m: [0.4, 0.8],
        damping_n_s_m: [0.2, 0.1],
        ..config()
    };
    // A symmetric gradient is conservative, so an impermeable uniformly
    // translating container has no shear flow. Its integrated WORLD force
    // nevertheless couples x and y endpoint positions in both equations.
    let gradient = [[0.3, 0.65], [0.65, -0.2]];
    let intercept = [0.7, -0.4];
    let force = |x: [f64; 2]| core::array::from_fn(|i| intercept[i] + dot(gradient[i], x));
    let initial = [0.3, -0.15];
    let mut owner = PlanarTranslationMotion::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        mass(),
        native_body(initial),
        &history(&mesh, initial),
        &edges,
        cfg,
    )
    .unwrap();
    let dt = 0.2;
    let applied = [0.4, -0.2];
    let volume = 1.2 * 0.8 * cfg.depth_m;
    let fluid_mass = 2.0 * volume;
    let total_mass = mass().mass() + fluid_mass;
    let mut v0 = initial;
    let mut q0 = [0.0; 2];
    for _ in 0..6 {
        let center0 = [0.6 + q0[0], 0.4 + q0[1]];
        let mut matrix = [[0.0; 2]; 2];
        let rhs: [f64; 2] = core::array::from_fn(|i| {
            total_mass * v0[i]
                + dt * (applied[i] + volume * (intercept[i] + dot(gradient[i], center0))
                    - cfg.stiffness_n_m[i] * q0[i])
        });
        for i in 0..2 {
            for j in 0..2 {
                matrix[i][j] = -dt * dt * volume * gradient[i][j];
            }
            matrix[i][i] += total_mass + dt * cfg.damping_n_s_m[i] + dt * dt * cfg.stiffness_n_m[i];
        }
        // Independent physical 2x2 oracle, not the production finite-difference
        // response or LU workspace. Both off-diagonal entries are nonzero.
        let determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        let expected = [
            (rhs[0] * matrix[1][1] - matrix[0][1] * rhs[1]) / determinant,
            (matrix[0][0] * rhs[1] - rhs[0] * matrix[1][0]) / determinant,
        ];
        let q1: [f64; 2] = core::array::from_fn(|i| q0[i] + dt * expected[i]);
        let endpoint_force = force([0.6 + q1[0], 0.4 + q1[1]]).map(|f| volume * f);
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
        let mut spring_change = 0.0;
        let mut spring_time_loss = 0.0;
        let mut damping = 0.0;
        for i in 0..2 {
            near(step.velocity_m_s[i], expected[i], 6e-8);
            near(step.translation_m[i], q1[i], 1e-7);
            near(
                step.fluid_force_n[i],
                endpoint_force[i] - fluid_mass * (expected[i] - v0[i]) / dt,
                6e-7,
            );
            let dx = q1[i] - q0[i];
            spring_change += 0.5 * cfg.stiffness_n_m[i] * (q1[i] + q0[i]) * dx;
            spring_time_loss += 0.5 * cfg.stiffness_n_m[i] * dx * dx;
            damping += dt * cfg.damping_n_s_m[i] * expected[i] * expected[i];
        }
        near(step.spring_energy_change_j, spring_change, 1e-8);
        near(step.damping_loss_j, damping, 1e-8);
        let dv = core::array::from_fn(|i| expected[i] - v0[i]);
        let sums = core::array::from_fn(|i| expected[i] + v0[i]);
        let storage = 0.5 * total_mass * dot(sums, dv) + spring_change;
        let time_loss = 0.5 * total_mass * dot(dv, dv) + spring_time_loss;
        let input = core::array::from_fn(|i| applied[i] + endpoint_force[i]);
        near(
            storage + damping + time_loss,
            dt * dot(input, expected),
            2e-10,
        );
        near(step.mechanical_work_residual_j, 0.0, 5e-9);
        assert!(step.flow_solves <= 4);
        v0 = step.velocity_m_s;
        q0 = step.translation_m;
    }
}

#[test]
fn open_pressure_and_tangential_viscosity_act_on_the_same_free_body() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let source = material(2.0, 0.5);
    let edges: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|&e| mesh.edges[e].normal[0] < -0.5)
        .collect();
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -2.0 * n);
    let exterior: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .filter(|e| !edges.contains(e))
        .map(|edge| MaterialNavierStokesBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    // Tangential motion has zero initial normal trace on the left face; the
    // initially resting liquid resolves its physical tangential drag.
    let initial = [0.0, 0.15];
    let mut owner = PlanarTranslationMotion::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        mass(),
        native_body(initial),
        &vec![0.0; system.n_u],
        &edges,
        config(),
    )
    .unwrap();
    let mut previous = initial;
    for _ in 0..3 {
        let step = owner
            .step(
                0.002,
                &zero,
                &exterior,
                &[None],
                [0.0; 2],
                fluid_controls(),
                controls(),
                || false,
            )
            .unwrap();
        assert!(
            step.velocity_m_s[0] < previous[0],
            "pressure must accelerate normal momentum"
        );
        assert!(step.velocity_m_s[1] > 0.0 && step.velocity_m_s[1] < previous[1]);
        assert!(step.fluid_force_n[0] < 0.0 && step.fluid_force_n[1] < 0.0);
        assert!(step.flow.relative_flow.iterations > 1);
        near(
            step.momentum.linear_increment_world.x / 0.002,
            step.fluid_force_n[0],
            6e-8,
        );
        near(
            step.momentum.linear_increment_world.y / 0.002,
            step.fluid_force_n[1],
            6e-8,
        );
        near(
            step.restraint_torque_n_m,
            -config().depth_m * step.boundary_load.stress.moment_n,
            1e-12,
        );
        near(step.mechanical_work_residual_j, 0.0, 5e-9);
        assert_eq!(step.boundary_load.pressure_offsets_pa, vec![(0, None)]);
        assert!(step.flow.world_kinetic_energy_j_per_m > 0.0);
        previous = step.velocity_m_s;
    }
}

#[test]
fn native_refusal_bounds_budget_and_final_cancellation_preserve_both_owners() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let source = material(2.0, 0.5);
    let edges = boundary_edges(&mesh);
    let make = || {
        PlanarTranslationMotion::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            mass(),
            native_body([0.0; 2]),
            &vec![0.0; system.n_u],
            &edges,
            config(),
        )
        .unwrap()
    };
    let mut owner = make();
    let mut reference = make();
    let applied = [0.7, -0.4];
    for motion in [&mut owner, &mut reference] {
        motion
            .step(
                0.02,
                &zero,
                &[],
                &[Some(0.0)],
                applied,
                fluid_controls(),
                controls(),
                || false,
            )
            .unwrap();
    }
    let before = (owner.body(), owner.flow_state().clone());
    let limited = PlanarTranslationControls {
        maximum_flow_solves: 1,
        ..controls()
    };
    assert!(matches!(
        owner.step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            applied,
            fluid_controls(),
            limited,
            || false
        ),
        Err(PlanarTranslationError::SolveLimit { flow_solves: 1, .. })
    ));
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    let bounded = PlanarTranslationControls {
        velocity_bounds_m_s: [[1.0, 2.0], [-2.0, -1.0]],
        ..controls()
    };
    assert!(
        owner
            .step(
                0.02,
                &zero,
                &[],
                &[Some(0.0)],
                applied,
                fluid_controls(),
                bounded,
                || false,
            )
            .is_err(),
        "the physical root lies outside the authored velocity box"
    );
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    assert!(
        owner
            .step(
                0.02,
                &zero,
                &[],
                &[None],
                applied,
                fluid_controls(),
                controls(),
                || false,
            )
            .is_err(),
        "native mean-zero pressure requires the original datum"
    );
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    let mut polls = 0;
    let expected = reference
        .step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            applied,
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
            applied,
            fluid_controls(),
            controls(),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(PlanarTranslationError::Cancelled)
    ));
    assert_eq!((owner.body(), owner.flow_state().clone()), before);
    let actual = owner
        .step(
            0.02,
            &zero,
            &[],
            &[Some(0.0)],
            applied,
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
    for invalid_body in [
        RigidBodyState::new(Pose::identity(), Vec3::new(0.0, 0.0, 0.1), Vec3::ZERO).unwrap(),
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, 0.1)).unwrap(),
    ] {
        assert!(
            PlanarTranslationMotion::new(
                &system,
                &source,
                MaterialFluxNumerics::default(),
                mass(),
                invalid_body,
                &vec![0.0; system.n_u],
                &edges,
                config(),
            )
            .is_err(),
            "unsupported momenta must never be projected away"
        );
    }
}
