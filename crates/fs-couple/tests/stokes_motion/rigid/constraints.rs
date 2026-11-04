//! G1/G3/G4: ideal supports close actual native fluid/body momentum.

use super::*;
use fs_couple::stokes_motion::rigid::{PlanarRigidConstraints, PlanarRigidSupport};

fn stationary_controls() -> PlanarRigidControls {
    PlanarRigidControls {
        velocity_bounds_m_s: [[0.0; 2]; 2],
        angular_velocity_bounds_rad_s: [0.0; 2],
        jacobian_step_m_s: [0.0; 2],
        jacobian_step_rad_s: 0.0,
        maximum_flow_solves: 1,
        ..controls()
    }
}

/// A guided body accelerates the closed unit-area fluid. The independent
/// effective inertia is M + rho*area*depth; transverse force and torque go
/// to the guide. Swap x/y to exercise both reduced search embeddings.
#[test]
fn g1_g3_guided_piston_retains_added_mass_and_skips_locked_axis_flow_probes() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let rho = 2.0;
    let material = FluidTransportState::try_new(rho, 1000.0, 0.5, 0.3).unwrap();
    let selected = edges(&mesh);
    let cfg = config();
    let dt = 0.03;
    for axis in 0..2 {
        let other = 1 - axis;
        let mut v0 = [0.0; 2];
        v0[axis] = 0.04;
        let initial = history(&mesh, |_| v0);
        let native = body([0.5; 2], v0, 0.0, 0.0);
        let make = || {
            PlanarRigidMotion::new(
                &system,
                &material,
                MaterialFluxNumerics::default(),
                mass(),
                native,
                &initial,
                &selected,
                cfg,
            )
            .unwrap()
        };
        let mut guided = make();
        let mut constraints = PlanarRigidConstraints::FIXED;
        constraints.translation_locked[axis] = false;
        guided.set_constraints(constraints).unwrap();
        let mut limits = stationary_controls();
        limits.velocity_bounds_m_s[axis] = [-1.0, 1.0];
        limits.jacobian_step_m_s[axis] = 1e-4;
        limits.maximum_flow_solves = 16;
        let mut external = [0.0; 2];
        external[axis] = 0.5;
        external[other] = -0.2;
        let torque = 0.07;
        let step = guided
            .step(
                dt,
                &zero,
                &[],
                &[Some(0.0)],
                external,
                torque,
                fluid_controls(),
                limits,
                || false,
            )
            .unwrap();
        let expected = v0[axis] + dt * external[axis] / (mass().mass() + rho * cfg.depth_m);
        near(step.velocity_m_s[axis], expected, 3e-8);
        assert_eq!(step.velocity_m_s[other], 0.0);
        assert_eq!(step.angular_velocity_rad_s, 0.0);
        assert_eq!(step.flow.pivot_world_m[other], 0.5);
        assert_eq!(step.rotation_rad, 0.0);
        near(
            step.fluid_force_n[axis],
            -rho * cfg.depth_m * (expected - v0[axis]) / dt,
            3e-8,
        );
        assert_eq!(step.support.force_n[axis], 0.0);
        near(step.support.force_n[other], -external[other], 3e-8);
        near(step.support.torque_n_m, -torque, 3e-8);
        near(step.support.impulse_n_s[other], -dt * external[other], 3e-9);
        near(step.support.angular_impulse_n_m_s, -dt * torque, 3e-9);
        assert_eq!(step.support.work_j, 0.0);
        near(step.mechanical_work_residual_j, 0.0, 3e-9);
        near(
            step.flow.world_momentum_n_s_per_m[axis],
            rho * expected,
            3e-8,
        );
        // All-free parity under the loads that actually do work on the guide.
        let mut free = make();
        external[other] = 0.0;
        let unrestricted = free
            .step(
                dt,
                &zero,
                &[],
                &[Some(0.0)],
                external,
                0.0,
                fluid_controls(),
                controls(),
                || false,
            )
            .unwrap();
        near(
            unrestricted.velocity_m_s[axis],
            step.velocity_m_s[axis],
            3e-8,
        );
        assert_eq!(unrestricted.support, PlanarRigidSupport::default());
        assert!(
            step.flow_solves + 2 <= unrestricted.flow_solves,
            "guide {}, free {}",
            step.flow_solves,
            unrestricted.flow_solves,
        );
    }
}

/// Independent changing shear and changing spin on a shaft. The actual
/// nonlinear fluid stress supplies torque to the rotor and force to its
/// support. The retained law is evaluated at endpoint strain and position.
#[test]
fn g1_shaft_reaction_and_nonlinear_shear_share_original_momentum_and_slip_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (rho, dt) = (2.0, 0.03);
    let material = FluidTransportState::try_new(rho, 1000.0, 0.5, 0.8).unwrap();
    let pivot = [0.35, 0.4];
    let selected: Vec<_> = edges(&mesh)
        .into_iter()
        .filter(|&e| mesh.edges[e].normal[1] < -0.5)
        .collect();
    let (gamma0, gamma1) = (0.1, 0.65);
    let (omega0, omega1) = (0.1, 0.4);
    let angle: f64 = dt * omega1;
    let rotation = [angle.cos(), angle.sin()];
    let mut cfg = config();
    cfg.rotation_stiffness_n_m_rad = 0.4;
    cfg.rotation_damping_n_m_s_rad = 0.05;
    for (exponent, beta) in [(0.55, Some(0.9)), (1.3, None)] {
        // Explicit numerical constitutive fixture, not a sourced material card.
        let viscosity =
            |g: f64| 0.06 + 0.74 * (1.0 + (1.8 * g).powi(2)).powf(0.5 * (exponent - 1.0));
        let mu = viscosity(gamma1);
        let tau = mu * gamma1;
        let slip0 = beta.map_or(0.0, |b| viscosity(gamma0) * gamma0 / b);
        let slip1 = beta.map_or(0.0, |b| tau / b);
        let old = |y: [f64; 2]| {
            let w = spin(omega0, [y[0] - pivot[0], y[1] - pivot[1]]);
            [w[0] + gamma0 * y[1] + slip0, w[1]]
        };
        let exact = |y: [f64; 2]| {
            let w = spin(omega1, [y[0] - pivot[0], y[1] - pivot[1]]);
            [w[0] + gamma1 * y[1] + slip1, w[1]]
        };
        let initial = history(&mesh, old);
        let force = |x: [f64; 2]| {
            let arm = inverse(rotation, [x[0] - pivot[0], x[1] - pivot[1]]);
            let y = [pivot[0] + arm[0], pivot[1] + arm[1]];
            let now = exact(y);
            let prior = inverse(rotation, old(y));
            rotate(
                rotation,
                [
                    rho * (now[0] - prior[0]) / dt,
                    rho * ((now[1] - prior[1]) / dt + omega1 * (gamma1 * y[1] + slip1)),
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
        let native = body(pivot, [0.0; 2], omega0, 0.37);
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
        owner
            .set_constraints(PlanarRigidConstraints::SHAFT)
            .unwrap();
        let physical_force = rotate(rotation, [cfg.depth_m * tau, 0.0]);
        let physical_torque = cfg.depth_m * pivot[1] * tau;
        let external = [0.13, -0.07];
        let external_torque = 0.3 * (omega1 - omega0) / dt - physical_torque
            + cfg.rotation_stiffness_n_m_rad * angle
            + cfg.rotation_damping_n_m_s_rad * omega1;
        let mut limits = controls();
        limits.velocity_bounds_m_s = [[0.0; 2]; 2];
        limits.jacobian_step_m_s = [0.0; 2];
        let step = owner
            .step_with_rheology(
                dt,
                &force,
                &exterior,
                &[None],
                external,
                external_torque,
                fluid_controls(),
                limits,
                |point, reference| {
                    FluidTransportState::try_new(
                        reference.density_kg_m3(),
                        reference.specific_heat_j_kg_k(),
                        reference.thermal_conductivity_w_m_k(),
                        viscosity(point.shear_rate_per_s),
                    )
                },
                || false,
            )
            .unwrap();
        assert_eq!(step.velocity_m_s, [0.0; 2]);
        assert_eq!(step.flow.pivot_world_m, pivot);
        near(step.angular_velocity_rad_s, omega1, 4e-7);
        near(step.rotation_rad, angle, 3e-8);
        for i in 0..2 {
            near(step.fluid_force_n[i], physical_force[i], 4e-7);
            near(
                step.support.force_n[i],
                -physical_force[i] - external[i],
                4e-7,
            );
            near(
                step.support.impulse_n_s[i],
                dt * step.support.force_n[i],
                1e-14,
            );
        }
        near(step.fluid_torque_n_m, physical_torque, 4e-7);
        assert_eq!(step.support.torque_n_m, 0.0);
        assert_eq!(step.support.work_j, 0.0);
        near(step.interface_work_j, dt * physical_torque * omega1, 4e-8);
        near(
            step.interface_slip_dissipation_j,
            beta.map_or(0.0, |b| dt * cfg.depth_m * tau * tau / b),
            4e-8,
        );
        near(step.mechanical_work_residual_j, 0.0, 4e-8);
        let (error, div) = system.velocity_error(step.flow.flow.solution.velocity_dofs(), &exact);
        assert!(error < 4e-7 && div < 4e-7, "{error} {div}");
        for sample in step.flow.rheology_samples() {
            near(sample.point().shear_rate_per_s, gamma1, 4e-6);
            near(sample.transport().dynamic_viscosity_pa_s(), mu, 4e-6);
        }
        assert_eq!(owner.constraints(), PlanarRigidConstraints::SHAFT);
    }
}

/// With no free body coordinates there is one real fluid solve and no Schur
/// solve. The moving liquid still exerts its physical shear force and torque.
#[test]
fn g1_fixed_body_returns_support_wrench_in_one_fluid_solve() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (mu, gamma, dt) = (0.4, 0.3, 0.2);
    let tau = mu * gamma;
    let material = FluidTransportState::try_new(2.0, 1000.0, 0.5, mu).unwrap();
    let pivot = [0.35, 0.4];
    let selected: Vec<_> = edges(&mesh)
        .into_iter()
        .filter(|&e| mesh.edges[e].normal[1] < -0.5)
        .collect();
    let initial = history(&mesh, |y| [gamma * y[1], 0.0]);
    let traction = |_: [f64; 2], n: [f64; 2]| [tau * n[1], tau * n[0]];
    let exterior: Vec<_> = edges(&mesh)
        .into_iter()
        .filter(|e| !selected.contains(e))
        .map(|edge| MaterialNavierStokesBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    let native = body(pivot, [0.0; 2], 0.0, 0.0);
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
    owner
        .set_constraints(PlanarRigidConstraints::FIXED)
        .unwrap();
    let external = [1.2, -0.3];
    let torque = 0.7;
    let step = owner
        .step(
            dt,
            &zero,
            &exterior,
            &[None],
            external,
            torque,
            fluid_controls(),
            stationary_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(step.flow_solves, 1);
    assert_eq!(step.body_after, native);
    near(step.fluid_force_n[0], cfg.depth_m * tau, 2e-8);
    near(step.fluid_torque_n_m, cfg.depth_m * pivot[1] * tau, 2e-8);
    near(
        step.support.force_n[0],
        -external[0] - cfg.depth_m * tau,
        2e-8,
    );
    near(step.support.force_n[1], -external[1], 2e-8);
    near(
        step.support.torque_n_m,
        -torque - cfg.depth_m * pivot[1] * tau,
        2e-8,
    );
    assert_eq!(step.support.work_j, 0.0);
    assert_eq!(step.body_kinetic_change_j, 0.0);
    assert_eq!(step.interface_work_j, 0.0);
    assert_eq!(step.mechanical_work_residual_j, 0.0);
    let (error, div) = system.velocity_error(step.flow.flow.solution.velocity_dofs(), &|y| {
        [gamma * y[1], 0.0]
    });
    assert!(error < 2e-8 && div < 2e-8, "{error} {div}");
    // Removing an ideal support does not kick or reset the accepted state.
    let accepted = owner.flow_state().clone();
    owner.set_constraints(PlanarRigidConstraints::FREE).unwrap();
    assert_eq!(owner.body(), native);
    assert_eq!(owner.flow_state(), &accepted);
}

/// Support activation cannot erase momentum; failed flow/support admission and
/// final cancellation preserve the original joint state and allow exact retry.
#[test]
fn g4_stationary_activation_reaction_overflow_and_final_cancel_are_atomic() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let material = FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.3).unwrap();
    let selected = edges(&mesh);
    let initial = history(&mesh, |_| [0.0; 2]);
    let native = body([0.5; 2], [0.0; 2], 0.0, 0.0);
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
    let v = [1e-14, 0.0];
    let moving_history = history(&mesh, |_| v);
    let mut moving = PlanarRigidMotion::new(
        &system,
        &material,
        MaterialFluxNumerics::default(),
        mass(),
        body([0.5; 2], v, 0.0, 0.0),
        &moving_history,
        &selected,
        config(),
    )
    .unwrap();
    let before_moving = moving.body();
    let before_flow = moving.flow_state().clone();
    assert!(
        moving
            .set_constraints(PlanarRigidConstraints::FIXED)
            .is_err()
    );
    assert_eq!(moving.constraints(), PlanarRigidConstraints::FREE);
    assert_eq!(moving.body(), before_moving);
    assert_eq!(moving.flow_state(), &before_flow);

    let mut owner = make();
    owner
        .set_constraints(PlanarRigidConstraints::FIXED)
        .unwrap();
    let before = owner.flow_state().clone();
    assert!(
        owner
            .step(
                2.0,
                &zero,
                &[],
                &[Some(0.0)],
                [f64::MAX, 0.0],
                0.0,
                fluid_controls(),
                stationary_controls(),
                || false,
            )
            .is_err(),
        "finite force with overflowing support impulse must refuse"
    );
    assert_eq!(owner.body(), native);
    assert_eq!(owner.flow_state(), &before);
    let mut bad = stationary_controls();
    bad.velocity_bounds_m_s[0] = [1.0, 2.0];
    assert!(
        owner
            .step(
                0.03,
                &zero,
                &[],
                &[Some(0.0)],
                [0.5, -0.2],
                0.07,
                fluid_controls(),
                bad,
                || false,
            )
            .is_err()
    );
    assert_eq!(owner.flow_state(), &before);

    let mut fresh = make();
    fresh
        .set_constraints(PlanarRigidConstraints::FIXED)
        .unwrap();
    let mut polls = 0;
    let expected = fresh
        .step(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.07,
            fluid_controls(),
            stationary_controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        owner.step(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.07,
            fluid_controls(),
            stationary_controls(),
            || {
                count += 1;
                count == polls
            },
        ),
        Err(PlanarRigidError::Cancelled)
    ));
    assert_eq!(count, polls);
    assert_eq!(owner.constraints(), PlanarRigidConstraints::FIXED);
    assert_eq!(owner.body(), native);
    assert_eq!(owner.flow_state(), &before);
    let retry = owner
        .step(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.07,
            fluid_controls(),
            stationary_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(owner.body(), fresh.body());
    assert_eq!(owner.flow_state(), fresh.flow_state());
    assert_eq!(retry.support, expected.support);
}
