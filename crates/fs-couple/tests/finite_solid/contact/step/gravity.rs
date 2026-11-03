//! Original spatial fields through actual mixed mechanical motion and impacts.
use super::*;
use fs_couple::finite_solid::contact::scene::{FiniteSolidSphereScene, SpherePairLaw};
use fs_couple::finite_solid::contact::step::free::{FreeFlightControls, FreeFlightError};
use fs_couple::finite_solid::gravity::StationaryGravity;
use fs_mbd::gravity::{
    CentralGravity, GravityError, GravityLoad, RigidGravityField, UniformGravity, step_in_field,
    step_in_field_with_loads,
};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use std::cell::Cell;
use std::sync::Arc;

const FRAME: u64 = 713;

// Instrument the original physical owner; never replace its force or potential.
struct ObservedField<'a> {
    field: &'a dyn RigidGravityField,
    samples: Cell<usize>,
    cancel: Option<(&'a CancelGate, usize)>,
}
impl RigidGravityField for ObservedField<'_> {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        let load = self.field.load(frame, pose, properties);
        self.samples.set(self.samples.get() + 1);
        if let Some((gate, at)) = self.cancel {
            if self.samples.get() == at {
                gate.request();
            }
        }
        load
    }
}
fn flight_settings(energy: f64) -> FreeFlightControls {
    FreeFlightControls {
        solid: controls(1e-9),
        energy_tolerance_j: energy,
        ..Default::default()
    }
}
fn central_force(x: Vec3, mass: f64, mu: f64) -> Vec3 {
    let radius = norm(x);
    x.scale(-mu * mass / radius.powi(3))
}
fn kinetic(body: RigidBodyState, property: MassProperties) -> f64 {
    body.linear_momentum_world()
        .dot(body.linear_momentum_world())
        / (2.0 * property.mass())
        + 0.5
            * body
                .angular_momentum_body()
                .dot(property.angular_velocity_body(body.angular_momentum_body()))
}
fn world_momentum(body: RigidBodyState) -> Vec3 {
    body.pose()
        .position_world()
        .cross(body.linear_momentum_world())
        .add(
            body.pose()
                .orientation()
                .rotate_body_to_world(body.angular_momentum_body()),
        )
}

#[test]
fn uniform_field_matches_legacy_map_with_bound_solid_gravity_and_pressure() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let gravity = Gravity::new(Vec3::new(0.4, -0.7, 0.2)).unwrap();
        let field = UniformGravity::new(FRAME, gravity, Vec3::new(1.0, -2.0, 3.0)).unwrap();
        let boundary = Arc::new(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap(),
        );
        let solid = FiniteSolidBody::new_with_gravity_and_sources(
            problem(&material, &[], &[]),
            None,
            Some((boundary, TetPressure::Uniform(0.3))),
            StationaryGravity::new(&field, FRAME),
            &[[0.0; 3]; 4],
            &[[0.1, 0.2, -0.1]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let property = properties();
        let body = body_at([1.5, -0.4, 0.7], Vec3::ZERO, Vec3::new(0.2, 0.1, -0.3));
        let wrench = Wrench {
            force_world: Vec3::new(0.2, 0.3, -0.1),
            torque_body: Vec3::new(0.0, 0.0, 0.4),
        };
        let observed = ObservedField {
            field: &field,
            samples: Cell::new(0),
            cancel: None,
        };
        let old = solid.native_state().clone();
        let expected = solid
            .rigid_free_flight_candidate_with_pressure(
                &[body],
                &[property],
                &[],
                [0.0; 3],
                0.3,
                &[wrench],
                gravity,
                0.01,
                [0.2, -0.3, 0.1],
                flight_settings(1e-8),
                cx,
            )
            .unwrap();
        let (next, rigid, report) = solid
            .rigid_free_flight_candidate_in_field(
                &[body],
                &[property],
                &[],
                [0.0; 3],
                Some(0.3),
                &[wrench],
                &observed,
                FRAME,
                0.01,
                [0.2, -0.3, 0.1],
                flight_settings(1e-8),
                cx,
            )
            .unwrap();
        assert_eq!(observed.samples.get(), 2);
        assert_eq!(rigid, expected.1);
        assert_eq!(next.native_state(), expected.0.native_state());
        near(report.source_work_j, expected.2.source_work_j, 3e-14);
        let dx = rigid[0]
            .pose()
            .position_world()
            .sub(body.pose().position_world());
        near(
            report.rigid_gravity_work_j,
            gravity.acceleration_world().scale(2.0).dot(dx),
            2e-14,
        );
        assert!(report.solid.gravity.is_some());
        assert!(report.solid.pressure_work_j.abs() > 1e-9);
        near(
            report.source_work_j,
            report.solid.source_work_j + report.rigid_held_work_j + report.rigid_gravity_work_j,
            2e-14,
        );
        near(report.rigid_gravity_work_defect_j, 0.0, 2e-14);
        near(report.rigid_held_force_work_defect_j, 0.0, 2e-14);
        same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 1e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 1e-10);
        assert_eq!(solid.native_state(), &old);
    });
}

#[test]
fn nonuniform_field_resamples_native_endpoints_and_opposing_work_errors_cannot_cancel() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let property = properties();
        let mass = property.mass();
        let mu = 8.0;
        let h = 0.2;
        let x0 = Vec3::new(2.0, 0.3, 0.0);
        let v0 = Vec3::new(0.2, 0.5, 0.0);
        let force = Vec3::new(1.0, 0.0, 0.0);
        let body = body_at([x0.x, x0.y, x0.z], Vec3::ZERO, v0);
        let field = CentralGravity::new(FRAME, Vec3::ZERO, mu).unwrap();
        // Independent scalar Newtonian force and velocity-Verlet oracle. All
        // expected errors and both gates are computed before the candidate.
        let g0 = central_force(x0, mass, mu);
        let vm = v0.add(force.add(g0).scale(0.5 * h / mass));
        let x1 = x0.add(vm.scale(h));
        let g1 = central_force(x1, mass, mu);
        let v1 = vm.add(force.add(g1).scale(0.5 * h / mass));
        let potential_work = mu * mass * (1.0 / norm(x1) - 1.0 / norm(x0));
        let gravity_kick = 0.25 * h * (g0.dot(v0.add(vm)) + g1.dot(vm.add(v1)));
        let held_work = force.dot(x1.sub(x0));
        let held_kick = 0.25 * h * force.dot(v0.add(vm.scale(2.0)).add(v1));
        let gravity_defect = gravity_kick - potential_work;
        let held_defect = held_kick - held_work;
        let absolute_defect = gravity_defect.abs() + held_defect.abs();
        assert!(gravity_defect * held_defect < 0.0);
        let tight_gate = 0.5 * (absolute_defect + (gravity_defect + held_defect).abs());
        let original = solid.native_state().clone();
        let wrench = Wrench {
            force_world: force,
            torque_body: Vec3::ZERO,
        };
        let run = |energy| {
            solid.rigid_free_flight_candidate_in_field(
                &[body],
                &[property],
                &[],
                [0.0; 3],
                None,
                &[wrench],
                &field,
                FRAME,
                h,
                [0.0; 3],
                flight_settings(energy),
                cx,
            )
        };
        let (_, next, report) = run(absolute_defect + 1e-9).unwrap();
        same(next[0].pose().position_world(), x1, 2e-14);
        same(next[0].linear_momentum_world(), v1.scale(mass), 2e-14);
        same(
            report.external_impulse_kg_m_s,
            g0.add(g1).scale(0.5 * h).add(force.scale(h)),
            2e-14,
        );
        same(
            report.external_angular_impulse_kg_m2_s,
            x0.cross(force.add(g0).scale(0.5 * h))
                .add(x1.cross(force.add(g1).scale(0.5 * h))),
            2e-14,
        );
        near(report.rigid_gravity_work_j, potential_work, 3e-14);
        near(report.rigid_gravity_work_defect_j, gravity_defect, 3e-14);
        near(report.rigid_held_force_work_defect_j, held_defect, 3e-14);
        near(report.absolute_integration_defect_j, absolute_defect, 4e-14);
        near(report.work_energy_residual_j, 0.0, 4e-14);
        let Err(FreeFlightError::Audit(rejected)) = run(tight_gate) else {
            panic!("opposing physical quadratures must not cancel through admission");
        };
        assert!(rejected.integration_defect_j.abs() < tight_gate);
        assert!(rejected.absolute_integration_defect_j > tight_gate);
        assert_eq!(solid.native_state(), &original);
        assert_eq!(body.pose().position_world(), x0);
    });
}

#[test]
fn explicitly_declared_quadrupole_retains_native_tidal_spin_and_orbital_exchange() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let property = MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.3, 0.4, 0.5),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.35).unwrap(),
        )
        .unwrap();
        let body = RigidBodyState::new(
            Pose::new(
                Vec3::new(5.0, 1.0, 0.3),
                UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 1.0), 0.4).unwrap(),
            )
            .unwrap(),
            Vec3::new(0.1, 0.7, -0.2),
            Vec3::new(0.2, -0.1, 0.3),
        )
        .unwrap();
        let field = CentralGravity::new(FRAME, Vec3::ZERO, 12.0)
            .unwrap()
            .quadrupole(1.0, 0.25)
            .unwrap();
        let observed = ObservedField {
            field: &field,
            samples: Cell::new(0),
            cancel: None,
        };
        let native =
            step_in_field_with_loads(&observed, FRAME, body, property, Wrench::ZERO, 0.01, || {
                false
            })
            .unwrap();
        assert_eq!(observed.samples.get(), 2);
        assert!(norm(native.load_before.wrench.torque_body) > 1e-4);
        assert_eq!(
            native.step,
            step_in_field(&field, FRAME, body, property, Wrench::ZERO, 0.01, || false).unwrap()
        );
        let (_, next, report) = solid
            .rigid_free_flight_candidate_in_field(
                &[body],
                &[property],
                &[],
                [0.0; 3],
                None,
                &[Wrench::ZERO],
                &field,
                FRAME,
                0.01,
                [0.0; 3],
                flight_settings(1e-5),
                cx,
            )
            .unwrap();
        assert_eq!(next[0], native.step.state_after);
        same(world_momentum(next[0]), world_momentum(body), 2e-12);
        same(report.external_angular_impulse_kg_m2_s, Vec3::ZERO, 2e-13);
        assert!(
            norm(
                next[0]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(next[0].angular_momentum_body())
                    .sub(
                        body.pose()
                            .orientation()
                            .rotate_body_to_world(body.angular_momentum_body())
                    )
            ) > 1e-6
        );
        near(
            report.rigid_gravity_work_j,
            native.load_before.potential_energy_j - native.load_after.potential_energy_j,
            0.0,
        );
        near(
            kinetic(next[0], property) - kinetic(body, property),
            report.rigid_gravity_work_j + report.integration_defect_j,
            2e-13,
        );
        near(report.work_energy_residual_j, 0.0, 2e-13);
    });
}

#[test]
fn actual_mixed_collision_retains_nonuniform_field_across_retries_and_continuation() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = [
            body_at([-0.625, -2.0, -2.0], Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)),
            body_at([0.625, -2.0, -2.0], Vec3::ZERO, Vec3::new(-1.0, 0.0, 0.0)),
        ];
        let scene = FiniteSolidSphereScene::new(
            &[sphere(3, 0, 0.0, 0.0), sphere(9, 1, 0.0, 0.0)],
            &[SpherePairLaw {
                first: 3,
                second: 9,
                friction_coefficient: 0.0,
                restitution: 1.0,
            }],
            1,
        )
        .unwrap();
        let field = CentralGravity::new(FRAME, Vec3::ZERO, 0.01).unwrap();
        let observed = ObservedField {
            field: &field,
            samples: Cell::new(0),
            cancel: None,
        };
        let mut c = settings();
        c.maximum_step_s = 0.05;
        c.free.energy_tolerance_j = 1e-6;
        c.integration_tolerance_j = 1e-5;
        let (next, rigid, report) = solid
            .mixed_sphere_trajectory_candidate_in_field(
                &surface,
                &scene,
                &initial,
                &[properties(); 2],
                &[],
                [0.0; 3],
                None,
                &[Wrench::ZERO; 2],
                &observed,
                FRAME,
                1.0,
                [0.0; 3],
                c,
                cx,
            )
            .unwrap();
        assert!(report.rejected_steps > 0);
        assert_eq!(
            observed.samples.get(),
            4 * (report.steps.len() + report.rejected_steps)
        );
        assert_eq!(
            report
                .steps
                .iter()
                .filter(|s| s.impact.as_ref().is_some_and(|r| r
                    .contact
                    .impulses
                    .iter()
                    .any(|j| j[0] > 1.0)))
                .count(),
            1
        );
        assert!(rigid[0].linear_momentum_world().x < -1.99);
        assert!(rigid[1].linear_momentum_world().x > 1.99);
        assert!(rigid.iter().all(|b| b.linear_momentum_world().y > 0.0));
        let potential = |b: RigidBodyState| -0.02 / norm(b.pose().position_world());
        let work: f64 = initial
            .into_iter()
            .zip(&rigid)
            .map(|(a, &b)| potential(a) - potential(b))
            .sum();
        near(report.source_work_j, work, 3e-13);
        near(report.unallocated_contact_loss_j, 0.0, 1e-10);
        near(report.work_energy_residual_j, 0.0, 1e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 1e-10);
        near(next.native_state().t, 1.0, 2e-15);
        let continued = next
            .mixed_sphere_trajectory_candidate_in_field(
                &surface,
                &scene,
                &rigid,
                &[properties(); 2],
                &[],
                [0.0; 3],
                None,
                &[Wrench::ZERO; 2],
                &field,
                FRAME,
                0.1,
                [0.0; 3],
                c,
                cx,
            )
            .unwrap();
        near(continued.0.native_state().t, 1.1, 2e-15);
        assert!(continued.1[0].pose().position_world().x < rigid[0].pose().position_world().x);
        assert_eq!(solid.native_state().steps, 0);
    });
}

#[test]
fn native_frame_and_late_field_domain_refusals_publish_no_mechanical_candidate() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let old = solid.native_state().clone();
        let center = Vec3::new(0.0, -2.0, -2.0);
        let field = CentralGravity::new(FRAME, center, 0.001)
            .unwrap()
            .quadrupole(0.3, 0.5)
            .unwrap();
        let observed = ObservedField {
            field: &field,
            samples: Cell::new(0),
            cancel: None,
        };
        let property = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.01, 0.01, 0.01)).unwrap();
        let body = body_at([1.0, -2.0, -2.0], Vec3::ZERO, Vec3::new(-1.0, 0.0, 0.0));
        assert!(matches!(
            solid.rigid_free_flight_candidate_in_field(
                &[body],
                &[property],
                &[],
                [0.0; 3],
                None,
                &[Wrench::ZERO],
                &field,
                FRAME + 1,
                0.1,
                [0.0; 3],
                flight_settings(1e-4),
                cx,
            ),
            Err(FreeFlightError::Gravity(GravityError::FrameMismatch))
        ));
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(1, 0, 0.0, 0.0)], &[], 0).unwrap();
        let mut c = settings();
        c.maximum_step_s = 0.15;
        c.free.energy_tolerance_j = 1e-4;
        c.integration_tolerance_j = 1e-3;
        let result = solid.mixed_sphere_trajectory_candidate_in_field(
            &surface,
            &scene,
            &[body],
            &[property],
            &[],
            [0.0; 3],
            None,
            &[Wrench::ZERO],
            &observed,
            FRAME,
            0.6,
            [0.0; 3],
            c,
            cx,
        );
        assert!(matches!(
            result,
            Err(FiniteSolidSphereTrajectoryError::Free(
                FreeFlightError::Gravity(GravityError::OutsideExpansionDomain)
            ))
        ));
        assert!(
            observed.samples.get() >= 6,
            "refusal follows completed internal intervals"
        );
        assert_eq!(solid.native_state(), &old);
        assert_eq!(body.pose().position_world(), Vec3::new(1.0, -2.0, -2.0));
    });
}

#[test]
fn cancellation_after_the_last_native_field_sample_keeps_both_original_owners() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let old = solid.native_state().clone();
        let field = CentralGravity::new(FRAME, Vec3::ZERO, 1.0).unwrap();
        let observed = ObservedField {
            field: &field,
            samples: Cell::new(0),
            cancel: Some((gate, 4)),
        };
        let bodies = [
            body_at([2.0, 0.0, 0.0], Vec3::ZERO, Vec3::ZERO),
            body_at([0.0, 2.0, 0.0], Vec3::ZERO, Vec3::ZERO),
        ];
        assert!(matches!(
            solid.rigid_free_flight_candidate_in_field(
                &bodies,
                &[properties(); 2],
                &[],
                [0.0; 3],
                None,
                &[Wrench::ZERO; 2],
                &observed,
                FRAME,
                0.01,
                [0.0; 3],
                flight_settings(1e-7),
                cx,
            ),
            Err(FreeFlightError::Cancelled)
        ));
        assert_eq!(observed.samples.get(), 4);
        assert_eq!(solid.native_state(), &old);
        assert!(
            bodies
                .iter()
                .all(|b| b.linear_momentum_world() == Vec3::ZERO)
        );
    });
}
