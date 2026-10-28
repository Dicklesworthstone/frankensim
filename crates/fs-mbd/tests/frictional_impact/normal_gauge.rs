//! G1/G3: a four-foot body can have unique motion without unique foot loads.
use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactDeterminacy, PsdVelocityCriterion};
use fs_mbd::frictional_impact::step::{RigidContactStepControls, step_rigid_frictional_contacts};
use fs_mbd::{Gravity, Wrench};

fn square() -> [RigidFrictionalContact; 4] {
    [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)].map(|(x, y)| {
        let mut c = plane(
            RigidImpactPoint {
                body: 0,
                arm_body: Vec3::new(x, y, -1.0),
            },
            0.1,
            0.0,
        );
        c.id = match (x > 0.0, y > 0.0) {
            (false, false) => 10,
            (true, false) => 20,
            (true, true) => 30,
            (false, true) => 40,
        };
        c.activation_tolerance_m = 1e-7;
        c
    })
}

#[test]
fn g1_four_foot_contact_determines_motion_without_selecting_a_pressure_distribution() {
    // An actual uniform 2 m cube of mass 3 kg has all principal inertias 2 kg m².
    // Four corner witnesses retain the one-dimensional normal reaction gauge;
    // no foot is deleted and no centre of pressure is supplied to the solver.
    let properties = [property(3.0, Vec3::new(2.0, 2.0, 2.0))];
    let contacts = square();
    for speed in [0.02, 0.1, 3.0] {
        let initial = [body(
            Vec3::new(0.0, 0.0, 1.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(3.0 * speed, 0.0, -3.0),
        )];
        let mut bodies = initial;
        let report = apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &contacts,
            controls(),
            &CancelGate::new(),
        )
        .unwrap_or_else(|error| panic!("four-foot speed={speed}: {error:?}"));
        let FrictionalContactDeterminacy::UniqueVelocity(witness) = &report.contact.determinacy
        else {
            panic!("the four foot loads are indeterminate even though body motion is unique");
        };
        let PsdVelocityCriterion::CandidateFrictionDirection {
            energy_coercivity_lower_bound,
            ..
        } = witness.criterion
        else {
            panic!("frictionful redundant normals require the stronger global witness");
        };
        assert_eq!(witness.dimension, 12);
        assert_eq!(witness.rank, 6);
        assert!(energy_coercivity_lower_bound > 0.8);
        let impulse_x = -(3.0 * speed).min(0.3);
        close_vector(
            bodies[0].linear_momentum_world(),
            Vec3::new(3.0 * speed + impulse_x, 0.0, 0.0),
        );
        close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
        let impulses = &report.contact.impulses;
        close(impulses.iter().map(|p| p[0]).sum(), 3.0);
        close(impulses.iter().map(|p| p[1]).sum(), impulse_x);
        close(
            impulses[1][0] + impulses[2][0] - impulses[0][0] - impulses[3][0],
            -impulse_x,
        );
        let before = 1.5 * (1.0 + speed * speed);
        let after = (3.0 * speed + impulse_x).powi(2) / 6.0;
        close(report.contact.dissipated_energy, before - after);
        close(report.total_energy_residual_j, 0.0);

        let mut permuted_bodies = initial;
        let repeated = apply_rigid_frictional_impacts(
            &mut permuted_bodies,
            &properties,
            &[contacts[2], contacts[0], contacts[3], contacts[1]],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(bodies, permuted_bodies);
        assert_eq!(report.contact, repeated.contact);
    }
}

#[test]
fn g3_four_foot_oblique_friction_retains_the_circular_cone_and_both_load_transfers() {
    let properties = [property(3.0, Vec3::new(2.0, 2.0, 2.0))];
    let contacts = square();
    // A genuine two-axis tangent direction exercises the circular Coulomb cone.
    // Clipping its coordinates independently would give the wrong motion.
    let direction = Vec3::new(0.6, 0.8, 0.0);
    for speed in [0.02, 0.1, 3.0] {
        let momentum = direction.scale(3.0 * speed).add(Vec3::new(0.0, 0.0, -3.0));
        let mut bodies = [body(
            Vec3::new(0.0, 0.0, 1.0),
            UnitQuaternion::IDENTITY,
            momentum,
        )];
        let report = apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &contacts,
            controls(),
            &CancelGate::new(),
        )
        .unwrap_or_else(|error| panic!("oblique four-foot speed={speed}: {error:?}"));
        let magnitude = (3.0 * speed).min(0.3);
        close_vector(
            bodies[0].linear_momentum_world(),
            direction.scale(3.0 * speed - magnitude),
        );
        close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
        let p = &report.contact.impulses;
        close(
            p[1][0] + p[2][0] - p[0][0] - p[3][0],
            magnitude * direction.x,
        );
        close(
            p[2][0] + p[3][0] - p[0][0] - p[1][0],
            magnitude * direction.y,
        );
        close(report.support_work_j, 0.0);
        close(report.total_energy_residual_j, 0.0);
    }
}

#[test]
fn g1_loaded_four_foot_body_crosses_friction_threshold_and_refines_its_real_trajectory() {
    let properties = [property(3.0, Vec3::new(2.0, 2.0, 2.0))];
    let contacts = square();
    let config = RigidContactStepControls {
        impact: controls(),
        drift_energy_tolerance_j: 1e-10,
        drift_angular_momentum_tolerance_nms: 1e-10,
        work_energy_tolerance_j: 3e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
    };
    let mut errors = [0.0; 3];
    for (rung, count) in [8_usize, 16, 32].into_iter().enumerate() {
        let dt = 1.0 / count as f64;
        for gx in [0.05_f64, 0.1, 0.2] {
            let initial_position = Vec3::new(0.0, 0.0, 1.0 + 1e-8);
            let mut bodies = [body(initial_position, UnitQuaternion::IDENTITY, Vec3::ZERO)];
            for step in 0..count {
                let report = step_rigid_frictional_contacts(
                    &mut bodies,
                    &properties,
                    &[Wrench::ZERO],
                    Gravity::new(Vec3::new(gx, 0.0, -1.0)).unwrap(),
                    &contacts,
                    dt,
                    config,
                    |path| -> Result<(), &'static str> {
                        // This contains the ENTIRE uniform cube. At every point
                        // of the original constant-angular-speed drift, a corner
                        // moves by at most sqrt(3)*|omega|*dt. COM height is affine.
                        let q = path.post_contact[0].pose().orientation();
                        let mut lower = f64::INFINITY;
                        for x in [-1.0, 1.0] {
                            for y in [-1.0, 1.0] {
                                for z in [-1.0, 1.0] {
                                    lower = lower.min(q.rotate_body_to_world(Vec3::new(x, y, z)).z);
                                }
                            }
                        }
                        let z0 = path.post_contact[0].pose().position_world().z;
                        let z1 = path.endpoint[0].pose().position_world().z;
                        let omega = path.post_contact[0]
                            .angular_momentum_body()
                            .norm_squared()
                            .sqrt()
                            / 2.0;
                        if z0.min(z1) + lower - 3.0_f64.sqrt() * omega * dt < 0.0 {
                            return Err("whole cube can intersect the support during drift");
                        }
                        Ok(())
                    },
                    &CancelGate::new(),
                )
                .unwrap_or_else(|error| panic!("gx={gx}, count={count}, step={step}: {error:?}"));
                close(
                    report.impact.contact.impulses.iter().map(|p| p[0]).sum(),
                    3.0 * dt,
                );
                close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
                close(report.total_energy_residual_j, 0.0);
            }
            let acceleration = (gx - 0.1).max(0.0);
            close_vector(
                bodies[0].linear_momentum_world(),
                Vec3::new(3.0 * acceleration, 0.0, 0.0),
            );
            // Independent first-order force/contact/drift recurrence, over 1 s.
            let x_discrete = 0.5 * acceleration * (1.0 + dt);
            close(bodies[0].pose().position_world().x, x_discrete);
            close(bodies[0].pose().position_world().z, initial_position.z);
            if gx > 0.1 {
                errors[rung] = (bodies[0].pose().position_world().x - 0.5 * acceleration).abs();
            }
        }
    }
    close(errors[0], 0.00625);
    close(errors[1], 0.003125);
    close(errors[2], 0.0015625);
}
