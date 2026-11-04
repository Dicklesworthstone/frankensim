//! G1/G3/G4: genuine native rotating fluid/body momentum and atomic acceptance.

use fs_couple::stokes_motion::rigid::{
    PlanarRigidConfig, PlanarRigidControls, PlanarRigidError, PlanarRigidMotion,
};
use fs_couple::stokes_motion::translating::planar::{
    PlanarTranslationConfig, PlanarTranslationControls, PlanarTranslationMotion,
};
use fs_flux::{
    FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary, NavierStokesControls,
    OpenBoundaryPolicy, TriMesh,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solid::Mesh2;

#[path = "rigid/constraints.rs"]
mod constraints;
#[path = "rigid/rheology.rs"]
mod rheology;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} versus {b:.16e}");
}
fn rotate(q: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    [q[0] * v[0] - q[1] * v[1], q[1] * v[0] + q[0] * v[1]]
}
fn inverse(q: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    rotate([q[0], -q[1]], v)
}
fn spin(omega: f64, r: [f64; 2]) -> [f64; 2] {
    [-omega * r[1], omega * r[0]]
}
fn zero(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn mass() -> MassProperties {
    MassProperties::new(0.7, Vec3::ZERO, Vec3::new(0.2, 0.25, 0.3)).unwrap()
}
fn body(pivot: [f64; 2], v: [f64; 2], omega: f64, angle: f64) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            Vec3::new(pivot[0], pivot[1], 0.0),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), angle).unwrap(),
        )
        .unwrap(),
        Vec3::new(v[0] * mass().mass(), v[1] * mass().mass(), 0.0),
        Vec3::new(0.0, 0.0, 0.3 * omega),
    )
    .unwrap()
}
fn config() -> PlanarRigidConfig {
    PlanarRigidConfig {
        depth_m: 0.4,
        translation_stiffness_n_m: [0.0; 2],
        rotation_stiffness_n_m_rad: 0.0,
        translation_damping_n_s_m: [0.0; 2],
        rotation_damping_n_m_s_rad: 0.0,
        velocity_tolerance_m_s: 1e-9,
        angular_velocity_tolerance_rad_s: 1e-9,
    }
}
fn controls() -> PlanarRigidControls {
    PlanarRigidControls {
        velocity_bounds_m_s: [[-4.0, 4.0]; 2],
        angular_velocity_bounds_rad_s: [-2.0, 2.0],
        jacobian_step_m_s: [1e-4; 2],
        jacobian_step_rad_s: 1e-4,
        maximum_flow_solves: 64,
        force_absolute_n: 2e-8,
        torque_absolute_n_m: 2e-8,
        momentum_relative: 1e-9,
        work_absolute_j: 2e-8,
        work_relative: 1e-9,
    }
}
fn fluid_controls() -> NavierStokesControls {
    NavierStokesControls {
        relative_tolerance: 1e-11,
        ..NavierStokesControls::default()
    }
}
fn history(mesh: &TriMesh, field: impl Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut values = vec![0.0; 2 * mesh.edges.len()];
    for (edge, g) in mesh.edges.iter().enumerate() {
        let a = field(mesh.verts[g.verts.0]);
        let b = field(mesh.verts[g.verts.1]);
        for axis in 0..2 {
            values[2 * edge] += 0.5 * (a[axis] + b[axis]) * g.normal[axis];
            values[2 * edge + 1] += (b[axis] - a[axis]) * g.normal[axis] / 12.0;
        }
    }
    values
}
fn edges(mesh: &TriMesh) -> Vec<usize> {
    (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect()
}

/// An independent affine solution has rigid rotation plus Couette shear.
/// The bottom wall transfers nonzero force AND torque about an offset COM.
/// The authored force is the exact covariant BDF1 manufactured forcing, with
/// nonlinear relative advection included; no native response supplies the oracle.
#[test]
fn rotating_shear_closes_three_actual_momenta_and_physical_slip_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (rho, mu, dt) = (2.0, 0.3, 0.03);
    let material = FluidTransportState::try_new(rho, 1000.0, 0.5, mu).unwrap();
    let pivot = [0.35, 0.4];
    let selected: Vec<_> = edges(&mesh)
        .into_iter()
        .filter(|&e| mesh.edges[e].normal[1] < -0.5)
        .collect();
    let mut cfg = config();
    cfg.translation_stiffness_n_m = [0.8, 1.3];
    cfg.rotation_stiffness_n_m_rad = 0.4;
    cfg.translation_damping_n_s_m = [0.07, 0.09];
    cfg.rotation_damping_n_m_s_rad = 0.05;
    let (omega0, omega1) = (0.1, 0.4);
    let angle: f64 = dt * omega1;
    let rotation = [angle.cos(), angle.sin()];
    for beta in [None, Some(1.2), Some(0.0)] {
        let gamma = if beta == Some(0.0) { 0.0 } else { 0.2 };
        let tau = mu * gamma;
        let slip = match beta {
            Some(b) if b > 0.0 => tau / b,
            Some(_) => 0.08,
            None => 0.0,
        };
        for boost in [[0.0; 2], [0.7, -0.4]] {
            let v0 = [0.05 + boost[0], -0.03 + boost[1]];
            let v1 = [-0.02 + boost[0], 0.04 + boost[1]];
            let end = [pivot[0] + dt * v1[0], pivot[1] + dt * v1[1]];
            let initial_field = |y: [f64; 2]| {
                let w = spin(omega0, [y[0] - pivot[0], y[1] - pivot[1]]);
                [v0[0] + w[0] + gamma * y[1] + slip, v0[1] + w[1]]
            };
            let initial = history(&mesh, initial_field);
            let v1_ref = inverse(rotation, v1);
            let exact = |y: [f64; 2]| {
                let w = spin(omega1, [y[0] - pivot[0], y[1] - pivot[1]]);
                [v1_ref[0] + w[0] + gamma * y[1] + slip, v1_ref[1] + w[1]]
            };
            let force = |x: [f64; 2]| {
                let r = inverse(rotation, [x[0] - end[0], x[1] - end[1]]);
                let y = [pivot[0] + r[0], pivot[1] + r[1]];
                let current = exact(y);
                let old = inverse(rotation, initial_field(y));
                let advective = [0.0, omega1 * (gamma * y[1] + slip)];
                rotate(
                    rotation,
                    [
                        rho * ((current[0] - old[0]) / dt + advective[0]),
                        rho * ((current[1] - old[1]) / dt + advective[1]),
                    ],
                )
            };
            let traction = |_: [f64; 2], normal: [f64; 2]| {
                let n = inverse(rotation, normal);
                rotate(rotation, [tau * n[1], tau * n[0]])
            };
            let exterior: Vec<_> = edges(&mesh)
                .into_iter()
                .filter(|e| !selected.contains(e))
                .map(|edge| MaterialNavierStokesBoundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                })
                .collect();
            let native = body(pivot, v0, omega0, 0.37);
            let mut owner = match beta {
                None => PlanarRigidMotion::new(
                    &system,
                    &material,
                    MaterialFluxNumerics::default(),
                    mass(),
                    native,
                    &initial,
                    &selected,
                    cfg,
                )
                .unwrap(),
                Some(b) => PlanarRigidMotion::new_with_slip(
                    &system,
                    &material,
                    MaterialFluxNumerics::default(),
                    mass(),
                    native,
                    &initial,
                    &selected,
                    cfg,
                    b,
                )
                .unwrap(),
            };
            let physical_force = rotate(rotation, [cfg.depth_m * tau, 0.0]);
            let physical_torque = cfg.depth_m * pivot[1] * tau;
            let external: [f64; 2] = core::array::from_fn(|i| {
                mass().mass() * (v1[i] - v0[i]) / dt - physical_force[i]
                    + cfg.translation_stiffness_n_m[i] * dt * v1[i]
                    + cfg.translation_damping_n_s_m[i] * v1[i]
            });
            let external_torque = 0.3 * (omega1 - omega0) / dt - physical_torque
                + cfg.rotation_stiffness_n_m_rad * angle
                + cfg.rotation_damping_n_m_s_rad * omega1;
            let step = owner
                .step(
                    dt,
                    &force,
                    &exterior,
                    &[None],
                    external,
                    external_torque,
                    fluid_controls(),
                    controls(),
                    || false,
                )
                .unwrap();
            for i in 0..2 {
                near(step.velocity_m_s[i], v1[i], 3e-7);
                near(step.fluid_force_n[i], physical_force[i], 3e-7);
                near(step.flow.pivot_world_m[i], end[i], 2e-8);
            }
            near(step.angular_velocity_rad_s, omega1, 3e-7);
            near(step.fluid_torque_n_m, physical_torque, 3e-7);
            near(step.rotation_rad, angle, 2e-8);
            near(
                step.body_after.angular_momentum_body().z - native.angular_momentum_body().z,
                0.3 * (omega1 - omega0),
                1e-7,
            );
            let (error, divergence) =
                system.velocity_error(step.flow.flow.solution.velocity_dofs(), &exact);
            assert!(error < 3e-7 && divergence < 3e-7, "{error}, {divergence}");
            let expected_work = dt
                * (physical_force[0] * v1[0]
                    + physical_force[1] * v1[1]
                    + physical_torque * omega1);
            near(step.interface_work_j, expected_work, 3e-8);
            let expected_slip = match beta {
                Some(b) if b > 0.0 => dt * cfg.depth_m * tau * tau / b,
                _ => 0.0,
            };
            near(step.interface_slip_dissipation_j, expected_slip, 2e-8);
            // Rigid spin and a Galilean boost contribute no viscous heat.
            // The unit-area affine shear dissipates mu*gamma^2 per depth.
            near(
                step.bulk_viscous_dissipation_j,
                dt * cfg.depth_m * mu * gamma * gamma,
                2e-9,
            );
            near(step.mechanical_work_residual_j, 0.0, 3e-8);
            let flux = step.flow.boundary_transport_world(&selected, end).unwrap();
            near(flux.outward_mass_flux_kg_per_m_s, 0.0, 3e-8);
            // The body material point and the moving native domain share pose.
            let original = Vec3::new(0.2, 0.6, 0.0);
            let arm = native
                .pose()
                .orientation()
                .rotate_world_to_body(original.sub(native.pose().position_world()));
            let actual = step.body_after.pose().point_world_from_body(arm).unwrap();
            let r = rotate(rotation, [original.x - pivot[0], original.y - pivot[1]]);
            near(actual.x, end[0] + r[0], 3e-8);
            near(actual.y, end[1] + r[1], 3e-8);
        }
    }
}

/// A closed uniformly translating fluid supplies its real added inertia.
/// The zero-spin limit must agree with the existing two-coordinate ALE owner,
/// which now shares the same bounded response search.
#[test]
fn zero_spin_limit_matches_native_translation_and_independent_added_mass() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let rho = 2.0;
    let material = FluidTransportState::try_new(rho, 1000.0, 0.5, 0.3).unwrap();
    let selected = edges(&mesh);
    let v0 = [0.04, -0.02];
    let initial = history(&mesh, |_| v0);
    let native = body([0.5, 0.5], v0, 0.0, 0.0);
    let cfg = config();
    let mut owner = PlanarRigidMotion::new(
        &system,
        &material,
        MaterialFluxNumerics::default(),
        mass(),
        native,
        &initial,
        &selected,
        cfg,
    )
    .unwrap();
    let mut baseline = PlanarTranslationMotion::new(
        &system,
        &material,
        MaterialFluxNumerics::default(),
        mass(),
        native,
        &initial,
        &selected,
        PlanarTranslationConfig {
            depth_m: cfg.depth_m,
            stiffness_n_m: [0.0; 2],
            damping_n_s_m: [0.0; 2],
            velocity_tolerance_m_s: cfg.velocity_tolerance_m_s,
        },
    )
    .unwrap();
    let (dt, applied) = (0.03, [0.5, -0.2]);
    let c = controls();
    let step = owner
        .step(
            dt,
            &zero,
            &[],
            &[Some(0.0)],
            applied,
            0.0,
            fluid_controls(),
            c,
            || false,
        )
        .unwrap();
    let other = baseline
        .step(
            dt,
            &zero,
            &[],
            &[Some(0.0)],
            applied,
            fluid_controls(),
            PlanarTranslationControls {
                velocity_bounds_m_s: c.velocity_bounds_m_s,
                jacobian_step_m_s: c.jacobian_step_m_s,
                maximum_flow_solves: c.maximum_flow_solves,
                force_absolute_n: c.force_absolute_n,
                force_relative: c.momentum_relative,
                work_absolute_j: c.work_absolute_j,
                work_relative: c.work_relative,
            },
            || false,
        )
        .unwrap();
    let fluid_mass = rho * cfg.depth_m;
    for i in 0..2 {
        let exact = v0[i] + dt * applied[i] / (mass().mass() + fluid_mass);
        near(step.velocity_m_s[i], exact, 3e-7);
        near(step.velocity_m_s[i], other.velocity_m_s[i], 3e-7);
        near(
            step.fluid_force_n[i],
            -fluid_mass * applied[i] / (mass().mass() + fluid_mass),
            3e-7,
        );
        near(step.translation_m[i], other.translation_m[i], 2e-8);
    }
    near(step.angular_velocity_rad_s, 0.0, 3e-7);
    near(step.fluid_torque_n_m, 0.0, 3e-7);
    near(step.mechanical_work_residual_j, 0.0, 3e-8);
}

/// A discrete manufactured rigidly rotating fluid continues beyond a half turn.
/// The forcing owns the covariant BDF1 temporal defect; this is a pose/history
/// test, not a claim of unforced exact angular conservation for the whole fluid.
#[test]
fn finite_rotation_continuation_retains_unwrapped_angle_and_original_body_axes() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let rho = 2.0;
    let material = FluidTransportState::try_new(rho, 1000.0, 0.5, 0.3).unwrap();
    let pivot = [0.5, 0.5];
    let (omega, dt) = (0.6, 0.7);
    let field = |y: [f64; 2]| spin(omega, [y[0] - pivot[0], y[1] - pivot[1]]);
    let initial = history(&mesh, field);
    let mut owner = PlanarRigidMotion::new(
        &system,
        &material,
        MaterialFluxNumerics::default(),
        mass(),
        body(pivot, [0.0; 2], omega, 2.8),
        &initial,
        &edges(&mesh),
        config(),
    )
    .unwrap();
    let delta: f64 = dt * omega;
    let turn = [delta.cos(), delta.sin()];
    for index in 1..=9 {
        let angle = f64::from(index) * delta;
        let q = [angle.cos(), angle.sin()];
        let force = |x: [f64; 2]| {
            let r = inverse(q, [x[0] - pivot[0], x[1] - pivot[1]]);
            let u = spin(omega, r);
            let old = inverse(turn, u);
            rotate(q, [rho * (u[0] - old[0]) / dt, rho * (u[1] - old[1]) / dt])
        };
        let step = owner
            .step(
                dt,
                &force,
                &[],
                &[Some(0.0)],
                [0.0; 2],
                0.0,
                fluid_controls(),
                controls(),
                || false,
            )
            .unwrap();
        near(step.rotation_rad, angle, 3e-7);
        near(step.angular_velocity_rad_s, omega, 3e-7);
        let axis = step
            .body_after
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        near(axis.x, (2.8 + angle).cos(), 3e-7);
        near(axis.y, (2.8 + angle).sin(), 3e-7);
        let (error, divergence) =
            system.velocity_error(step.flow.flow.solution.velocity_dofs(), &field);
        assert!(error < 3e-7 && divergence < 3e-7);
        near(step.mechanical_work_residual_j, 0.0, 3e-8);
    }
    assert!(owner.rotation_rad() > core::f64::consts::PI);
    near(owner.flow_state().time_s(), 9.0 * dt, 1e-12);
}

#[test]
fn malformed_planar_data_budget_and_final_cancellation_preserve_the_joint_owner() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let material = FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.3).unwrap();
    let selected = edges(&mesh);
    let initial = vec![0.0; system.n_u];
    let native = body([0.5, 0.5], [0.0; 2], 0.0, 0.0);
    let make = || {
        PlanarRigidMotion::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            mass(),
            native,
            &initial,
            &selected,
            config(),
        )
        .unwrap()
    };
    let mut owner = make();
    let before = owner.flow_state().clone();
    let mut limited = controls();
    limited.maximum_flow_solves = 1;
    assert!(matches!(
        owner.step(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.0,
            fluid_controls(),
            limited,
            || false
        ),
        Err(PlanarRigidError::SolveLimit { .. })
    ));
    assert_eq!(owner.body(), native);
    assert_eq!(owner.flow_state(), &before);
    assert_eq!(owner.rotation_rad(), 0.0);
    // A closed pressure component requires an explicit absolute offset.
    assert!(
        owner
            .step(
                0.03,
                &zero,
                &[],
                &[None],
                [0.5, -0.2],
                0.0,
                fluid_controls(),
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(owner.flow_state(), &before);
    let mut reference = make();
    let mut polls = 0;
    reference
        .step(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.0,
            fluid_controls(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        owner.step(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.0,
            fluid_controls(),
            controls(),
            || {
                observed += 1;
                observed == polls
            }
        ),
        Err(PlanarRigidError::Cancelled)
    ));
    assert_eq!(observed, polls);
    assert_eq!(owner.body(), native);
    assert_eq!(owner.flow_state(), &before);
    assert_eq!(owner.rotation_rad(), 0.0);
    owner
        .step(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.0,
            fluid_controls(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(owner.body(), reference.body());
    assert_eq!(owner.flow_state(), reference.flow_state());
    assert_eq!(owner.rotation_rad(), reference.rotation_rad());

    let tipped = RigidBodyState::new(
        Pose::new(
            native.pose().position_world(),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.1).unwrap(),
        )
        .unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap();
    assert!(
        PlanarRigidMotion::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            mass(),
            tipped,
            &initial,
            &selected,
            config()
        )
        .is_err()
    );
    let spinning = body([0.5, 0.5], [0.0; 2], 0.4, 0.0);
    let mut wrong = history(&mesh, |y| spin(0.4, [y[0] - 0.5, y[1] - 0.5]));
    wrong[2 * selected[0] + 1] += 1e-3;
    assert!(matches!(
        PlanarRigidMotion::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            mass(),
            spinning,
            &wrong,
            &selected,
            config()
        ),
        Err(PlanarRigidError::InterfaceVelocity { .. })
    ));
    for beta in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PlanarRigidMotion::new_with_slip(
                &system,
                &material,
                MaterialFluxNumerics::default(),
                mass(),
                native,
                &initial,
                &selected,
                config(),
                beta
            )
            .is_err()
        );
    }
}
