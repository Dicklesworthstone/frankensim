//! G1/G3/G4: sustained contact with original articulated inertia and real drift.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};
use fs_mbd::articulated_contacts::friction::ArticulatedFrictionalContact;
use fs_mbd::articulated_contacts::step::{
    ArticulatedContactPath, ArticulatedContactStepConfig, ArticulatedContactStepError,
    ArticulatedContactStepResult, step_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::trajectory::ArticulatedTrajectorySystem;
use fs_mbd::articulated_contacts::{ContactCounterpart, LinkContactPoint};
use std::cell::Cell;

const CLEARANCE: f64 = 1e-8;
fn inertia(mass: f64, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::new(0.0, 0.5, 0.0),
        Mat3 {
            m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
        },
    )
    .unwrap()
}
fn carrier(limits: Option<JointLimits>) -> ArticulatedModel {
    // An exposed housing with a horizontally translating internal carriage.
    // Its two feet are on y=0 and every link COM is at y=1/2 in its own frame.
    ArticulatedModel::new(vec![
        Link::new(
            "housing",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, 0.5),
        ),
        Link::new(
            "carriage",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), limits).unwrap(),
            inertia(1.0, 0.02),
        ),
    ])
    .unwrap()
}
fn block(mass: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "block",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(mass, mass / 3.0),
    )])
    .unwrap()
}
fn state(velocity: f64) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), Vec3::new(0.0, CLEARANCE, 0.0)).unwrap(),
            Twist::new(Vec3::default(), Vec3::new(velocity, 0.0, 0.0)),
        ),
        joint_position: vec![0.0],
        joint_velocity: vec![0.0],
    }
}
fn config(duration: f64) -> ArticulatedContactStepConfig {
    ArticulatedContactStepConfig {
        duration,
        angular_origin_world: Vec3::new(3.0, -2.0, 1.0),
        contact: FrictionalContactConfig::default(),
        drift_energy_tolerance: 1e-8,
        energy_tolerance: 1e-8,
        linear_momentum_tolerance: 1e-8,
        angular_momentum_tolerance: 1e-8,
    }
}
fn feet(mu: f64) -> [ArticulatedFrictionalContact; 2] {
    [-1.0, 1.0].map(|x| ArticulatedFrictionalContact {
        id: if x < 0.0 { 10 } else { 20 },
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::new(x, 0.0, 0.0),
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::default(),
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-7,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::Target { speed: 0.0 },
    })
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 3e-8,
        "{actual:?} != {expected:?}"
    );
}
fn audit<P>(result: &ArticulatedContactStepResult<P>) {
    close(result.kinetic_work_residual, result.drift_energy_defect);
    vector(result.linear_momentum_residual, Vec3::default());
    vector(result.angular_momentum_residual, Vec3::default());
    vector(result.drift_angular_momentum_defect, Vec3::default());
    assert!(
        result
            .initial_gaps
            .iter()
            .chain(&result.endpoint_gaps)
            .all(|gap| gap.signed_gap >= 0.0)
    );
}
fn translated_housing(path: &ArticulatedContactPath<'_>) -> Result<(), &'static str> {
    // This manufactured family has only a horizontal prismatic coordinate,
    // zero angular motion and constant post-contact velocities. Native drift
    // is therefore straight translation for every link and the complete
    // housing lower face has constant clearance. Sampling checks the native
    // realization of that analytic family; it is not a general CCD argument.
    for state in path.post_contact.iter().chain(path.endpoint) {
        if state.base.twist_body.angular.norm() > 1e-10
            || state.base.twist_body.linear.y.abs() > 1e-10
            || state.base.world_from_base.translation().y < 0.5 * CLEARANCE
        {
            return Err("path left the admitted horizontal translation family");
        }
    }
    let midpoint = path
        .sample(0.5, || false)
        .map_err(|_| "native sample refused")?;
    for ((first, middle), last) in path.post_contact.iter().zip(midpoint).zip(path.endpoint) {
        vector(
            middle.base.world_from_base.translation(),
            (first.base.world_from_base.translation() + last.base.world_from_base.translation())
                .scale(0.5),
        );
        vector(middle.base.twist_body.angular, Vec3::default());
    }
    Ok(())
}

#[test]
fn gravity_and_actuation_support_the_carrier_while_its_internal_carriage_accelerates() {
    let model = carrier(None);
    let zero = [Wrench::default(); 2];
    let efforts = [2.0];
    let contacts = feet(0.5);
    let dt = 1.0 / 32.0;
    let mut current = state(0.0);
    let mut work = 0.0;
    let mut displacement_work = 0.0;
    for step in 0..8 {
        let q_before = current.joint_position[0];
        let result = step_articulated_frictional_contacts(
            &[ArticulatedTrajectorySystem {
                model: &model,
                state: &current,
                loads: FreeFloatingMotionLoads {
                    generalized_force: &efforts,
                    external_body_wrench: &zero,
                    gravity_world: Vec3::new(0.0, -8.0, 0.0),
                },
            }],
            &contacts,
            config(dt),
            translated_housing,
            || false,
        )
        .unwrap();
        audit(&result);
        let impulses = &result.contact.contact.impulses;
        close(impulses[0][0] + impulses[1][0], 24.0 * dt);
        close(impulses[0][1] + impulses[1][1], 2.0 * dt);
        // Global angular balance independently fixes the normal load transfer.
        close(impulses[1][0] - impulses[0][0], (8.0 * q_before - 1.0) * dt);
        vector(
            result.states_after[0].base.twist_body.linear,
            Vec3::default(),
        );
        vector(
            result.states_after[0].base.twist_body.angular,
            Vec3::default(),
        );
        close(
            result.states_after[0].joint_velocity[0],
            2.0 * dt * (step + 1) as f64,
        );
        close(result.gravity_displacement_work, 0.0);
        close(result.gravity_impulse_work, 96.0 * dt * dt);
        // The free motor kick also gives the 2 kg housing speed -F*dt/2.
        // Support removes that 1/2*m*v^2 = F^2*dt^2/4 along with the
        // gravity kick. This first-order reaction loss is not friction heat.
        close(
            result.contact.contact.dissipated_energy,
            result.gravity_impulse_work + dt * dt,
        );
        work += result.actuator_impulse_work;
        displacement_work += result.actuator_displacement_work;
        current = result.states_after[0].clone();
    }
    close(current.joint_position[0], 2.0 * dt * dt * 36.0);
    close(work, 0.5 * 0.5_f64.powi(2) + 8.0 * dt * dt);
    close(displacement_work - work, 8.0 * dt * dt);
}

#[test]
fn sliding_keeps_the_free_carriage_momentum_and_converges_at_first_order() {
    let model = carrier(None);
    let wrenches = [Wrench::default(); 2];
    let effort = [0.0];
    let contacts = feet(0.2);
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let dt = 0.25 / count as f64;
        let mut current = state(3.0);
        for _ in 0..count {
            let result = step_articulated_frictional_contacts(
                &[ArticulatedTrajectorySystem {
                    model: &model,
                    state: &current,
                    loads: FreeFloatingMotionLoads {
                        generalized_force: &effort,
                        external_body_wrench: &wrenches,
                        gravity_world: Vec3::new(0.0, -8.0, 0.0),
                    },
                }],
                &contacts,
                config(dt),
                translated_housing,
                || false,
            )
            .unwrap();
            audit(&result);
            close(
                result.contact.contact.impulses.iter().map(|p| p[0]).sum(),
                24.0 * dt,
            );
            close(
                result.contact.contact.impulses.iter().map(|p| p[1]).sum(),
                -4.8 * dt,
            );
            current = result.states_after[0].clone();
        }
        let root_x = current.base.world_from_base.translation().x;
        close(current.base.twist_body.linear.x, 3.0 - 2.4 * 0.25);
        close(current.joint_velocity[0], 2.4 * 0.25);
        close(root_x + current.joint_position[0], 3.0 * 0.25);
        close(root_x, 3.0 * 0.25 - 0.5 * 2.4 * 0.25 * (0.25 + dt));
        errors.push((root_x - (3.0 * 0.25 - 0.5 * 2.4 * 0.25_f64.powi(2))).abs());
    }
    eprintln!("sliding carrier position errors: {errors:?}");
    close(errors[0] / errors[1], 2.0);
    close(errors[1] / errors[2], 2.0);
}

#[test]
fn separated_trees_transmit_static_load_through_dynamic_contact_before_drifting() {
    let lower = block(2.0);
    let upper = block(1.0);
    let mut a = state(0.0);
    a.joint_position.clear();
    a.joint_velocity.clear();
    let mut b = a.clone();
    b.base.world_from_base =
        Se3::from_parts(So3::identity(), Vec3::new(0.0, 1.0 + 2.0 * CLEARANCE, 0.0)).unwrap();
    let zero = [Wrench::default()];
    let pressure = [Wrench::new(Vec3::default(), Vec3::new(0.0, -3.0, 0.0))];
    let mut contacts = feet(0.0).to_vec();
    contacts.push(ArticulatedFrictionalContact {
        id: 30,
        point: LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::new(0.0, 1.0, 0.0),
        }),
        ..contacts[0]
    });
    let dt = 1.0 / 16.0;
    for _ in 0..16 {
        let loads = |wrenches| FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: wrenches,
            gravity_world: Vec3::new(0.0, -8.0, 0.0),
        };
        let result = step_articulated_frictional_contacts(
            &[
                ArticulatedTrajectorySystem {
                    model: &lower,
                    state: &a,
                    loads: loads(&zero),
                },
                ArticulatedTrajectorySystem {
                    model: &upper,
                    state: &b,
                    loads: loads(&pressure),
                },
            ],
            &contacts,
            config(dt),
            |path| {
                // Every accepted spatial velocity is zero: the entire two-box
                // numerical drift is stationary and all exposed faces stay apart.
                for state in path.post_contact.iter().chain(path.endpoint) {
                    vector(state.base.twist_body.linear, Vec3::default());
                    vector(state.base.twist_body.angular, Vec3::default());
                }
                Ok::<_, &'static str>(())
            },
            || false,
        )
        .unwrap();
        audit(&result);
        close(result.contact.contact.impulses[2][0], 11.0 * dt);
        close(
            result.contact.contact.impulses[0][0] + result.contact.contact.impulses[1][0],
            27.0 * dt,
        );
        vector(result.external_linear_impulse_world, Vec3::default());
        a = result.states_after[0].clone();
        b = result.states_after[1].clone();
    }
    close(a.base.world_from_base.translation().y, CLEARANCE);
    close(
        b.base.world_from_base.translation().y,
        1.0 + 2.0 * CLEARANCE,
    );
}

#[test]
fn translating_support_work_and_contact_permutation_preserve_the_original_solution() {
    let model = carrier(None);
    let current = state(1.0);
    let wrenches = [Wrench::default(); 2];
    let effort = [2.0];
    let mut contacts = feet(0.5);
    for c in &mut contacts {
        c.other = ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(1.0, 0.0, 0.0),
        };
    }
    let systems = [ArticulatedTrajectorySystem {
        model: &model,
        state: &current,
        loads: FreeFloatingMotionLoads {
            generalized_force: &effort,
            external_body_wrench: &wrenches,
            gravity_world: Vec3::new(0.0, -8.0, 0.0),
        },
    }];
    let first = step_articulated_frictional_contacts(
        &systems,
        &contacts,
        config(0.125),
        translated_housing,
        || false,
    )
    .unwrap();
    contacts.reverse();
    let second = step_articulated_frictional_contacts(
        &systems,
        &contacts,
        config(0.125),
        translated_housing,
        || false,
    )
    .unwrap();
    assert_eq!(first.states_after, second.states_after);
    assert_eq!(first.contact.contact, second.contact.contact);
    assert_eq!(first.initial_gaps, second.initial_gaps);
    close(first.contact.plane_work, 0.25);
    close(
        first.states_after[0].base.world_from_base.translation().x,
        0.125,
    );
    audit(&first);
}

#[test]
fn whole_scene_rotation_retains_spatial_momentum_and_all_work_terms() {
    let model = carrier(None);
    let current = state(1.0);
    let wrenches = [Wrench::default(); 2];
    let effort = [0.0];
    let contacts = feet(0.2);
    let base = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &current,
            loads: FreeFloatingMotionLoads {
                generalized_force: &effort,
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, -8.0, 0.0),
            },
        }],
        &contacts,
        config(0.125),
        translated_housing,
        || false,
    )
    .unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.1, 0.3))).unwrap();
    let rotate = |v| rotation.rotate(v).unwrap();
    let mut rotated = current.clone();
    rotated.base.world_from_base =
        Se3::from_parts(rotation, rotate(current.base.world_from_base.translation())).unwrap();
    let contacts = contacts.map(|c| ArticulatedFrictionalContact {
        normal_world: rotate(c.normal_world),
        tangent_world: rotate(c.tangent_world),
        ..c
    });
    let settings = ArticulatedContactStepConfig {
        angular_origin_world: rotate(config(0.125).angular_origin_world),
        ..config(0.125)
    };
    let result = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &rotated,
            loads: FreeFloatingMotionLoads {
                generalized_force: &effort,
                external_body_wrench: &wrenches,
                gravity_world: rotate(Vec3::new(0.0, -8.0, 0.0)),
            },
        }],
        &contacts,
        settings,
        |path| {
            // A fixed rotation of the same admitted horizontal translation family.
            let endpoint = &path.endpoint[0];
            vector(
                endpoint.base.world_from_base.translation(),
                rotate(base.states_after[0].base.world_from_base.translation()),
            );
            Ok::<_, &'static str>(())
        },
        || false,
    )
    .unwrap();
    audit(&result);
    vector(
        result.momentum_after.linear_world,
        rotate(base.momentum_after.linear_world),
    );
    vector(
        result.momentum_after.angular_about_origin_world,
        rotate(base.momentum_after.angular_about_origin_world),
    );
    vector(
        result.gravity_angular_splitting_impulse_world,
        rotate(base.gravity_angular_splitting_impulse_world),
    );
    close(result.gravity_impulse_work, base.gravity_impulse_work);
    close(
        result.contact.contact.dissipated_energy,
        base.contact.contact.dissipated_energy,
    );
    close(
        result.states_after[0].joint_position[0],
        base.states_after[0].joint_position[0],
    );
}

#[test]
fn late_geometry_cancellation_and_penetration_never_publish_a_partial_step() {
    let model = carrier(None);
    let current = state(0.0);
    let saved = current.clone();
    let wrenches = [Wrench::default(); 2];
    let effort = [0.0];
    let contacts = feet(0.5);
    let systems = [ArticulatedTrajectorySystem {
        model: &model,
        state: &current,
        loads: FreeFloatingMotionLoads {
            generalized_force: &effort,
            external_body_wrench: &wrenches,
            gravity_world: Vec3::new(0.0, -8.0, 0.0),
        },
    }];
    let refusal = step_articulated_frictional_contacts(
        &systems,
        &contacts,
        config(0.1),
        |path| {
            assert_eq!(path.endpoint.len(), 1);
            Err::<(), _>("omitted feature intersects")
        },
        || false,
    );
    assert!(matches!(
        refusal,
        Err(ArticulatedContactStepError::Path(
            "omitted feature intersects"
        ))
    ));
    let cancel = Cell::new(false);
    let refusal = step_articulated_frictional_contacts(
        &systems,
        &contacts,
        config(0.1),
        |_| {
            cancel.set(true);
            Ok::<_, &'static str>(())
        },
        || cancel.get(),
    );
    assert!(matches!(
        refusal,
        Err(ArticulatedContactStepError::Cancelled)
    ));
    assert_eq!(current, saved);
    let mut elevated = current.clone();
    elevated.base.world_from_base =
        Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.01, 0.0)).unwrap();
    let refusal = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            state: &elevated,
            ..systems[0]
        }],
        &contacts,
        config(0.1),
        |_| Ok::<_, &'static str>(()),
        || false,
    );
    assert!(matches!(
        refusal,
        Err(ArticulatedContactStepError::Penetration { endpoint: true, .. })
    ));
    let mut penetrating = current.clone();
    penetrating.base.world_from_base =
        Se3::from_parts(So3::identity(), Vec3::new(0.0, -1e-12, 0.0)).unwrap();
    let refusal = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            state: &penetrating,
            ..systems[0]
        }],
        &contacts,
        config(0.1),
        |_| Ok::<_, &'static str>(()),
        || false,
    );
    assert!(matches!(
        refusal,
        Err(ArticulatedContactStepError::Penetration {
            endpoint: false,
            ..
        })
    ));
    let retry = step_articulated_frictional_contacts(
        &systems,
        &contacts,
        config(0.1),
        translated_housing,
        || false,
    )
    .unwrap();
    audit(&retry);

    // A freely rotating/sliding tree has nonlinear internal drift. A budget
    // below its actual energy defect refuses without invoking geometry.
    let mut curved = current.clone();
    curved.base.twist_body.angular = Vec3::new(0.3, 0.2, 0.4);
    curved.joint_position[0] = 0.2;
    curved.joint_velocity[0] = 0.7;
    let consulted = Cell::new(false);
    let narrow = ArticulatedContactStepConfig {
        drift_energy_tolerance: 1e-14,
        ..config(0.1)
    };
    let refusal = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            state: &curved,
            loads: FreeFloatingMotionLoads {
                gravity_world: Vec3::default(),
                ..systems[0].loads
            },
            ..systems[0]
        }],
        &[],
        narrow,
        |_| {
            consulted.set(true);
            Ok::<_, &'static str>(())
        },
        || false,
    );
    assert!(matches!(
        refusal,
        Err(ArticulatedContactStepError::Balance(_))
    ));
    assert!(!consulted.get());
}

#[test]
fn force_limits_are_validated_before_impulse_conversion() {
    let model = carrier(Some(JointLimits::new(-1.0, 1.0, 100.0, 0.5).unwrap()));
    let current = state(0.0);
    let wrenches = [Wrench::default(); 2];
    let effort = [1.0];
    let result = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &current,
            loads: FreeFloatingMotionLoads {
                generalized_force: &effort,
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, -8.0, 0.0),
            },
        }],
        &feet(0.5),
        config(1e-6),
        translated_housing,
        || false,
    );
    assert!(matches!(
        result,
        Err(ArticulatedContactStepError::Dynamics {
            error: ArticulatedError::JointEffortLimitViolation { .. },
            ..
        })
    ));
}

#[test]
fn off_origin_body_load_keeps_force_moment_and_spin_impulse_once() {
    let model = block(2.0);
    let mut current = state(0.6);
    current.joint_position.clear();
    current.joint_velocity.clear();
    // Force at the actual COM plus a pure 0.4 N m body torque. About the
    // link origin, c cross F + tau = (0,0,-0.1), not merely tau.
    let wrenches = [Wrench::new(
        Vec3::new(0.0, 0.0, -0.1),
        Vec3::new(1.0, 0.0, 0.0),
    )];
    let result = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &current,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::default(),
            },
        }],
        &[],
        config(0.2),
        |_| Ok::<_, &'static str>(()), // Declared empty collision scene.
        || false,
    )
    .unwrap();
    audit(&result);
    vector(
        result.external_linear_impulse_world,
        Vec3::new(0.2, 0.0, 0.0),
    );
    let arm = Vec3::new(0.0, 0.5 + CLEARANCE, 0.0) - config(0.2).angular_origin_world;
    vector(
        result.external_angular_impulse_about_origin_world,
        arm.cross(Vec3::new(0.2, 0.0, 0.0)) + Vec3::new(0.0, 0.0, 0.08),
    );
    vector(
        result.states_after[0].base.twist_body.angular,
        Vec3::new(0.0, 0.0, 0.12),
    );
    close(result.body_wrench_impulse_work, 0.1348);
    close(result.gravity_impulse_work, 0.0);
    close(result.actuator_impulse_work, 0.0);
}
