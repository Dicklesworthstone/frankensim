//! Automatic contact on translating native articulated links.
//!
//! A nonrotating prescribed base and fixed/prismatic joints leave every link
//! orientation and prismatic mass matrix constant. Held native efforts, body
//! wrenches and gravity therefore give constant joint accelerations. The
//! original ABA and material-point map supply each complete quadratic centre
//! path; no independent slider mass or sampled trajectory is substituted.

use super::planes::{AcceleratingPlane, MixedSphereSceneResult, SpherePlaneLaw};
use super::{
    CentreCommand, ContactSphere, MixedSphereControls, MixedSphereError, Pair, SpherePairLaw,
    State, System, admit_carrier, admit_pairs, admit_window, base, centre_velocity, finite_vector,
    poll,
};
use crate::articulated::{
    ArticulatedError, JointKind, Kinematics, forward_dynamics, forward_kinematics,
};
use crate::articulated_contacts::trajectory::ArticulatedTrajectoryError;
use crate::articulated_contacts::{ArticulatedContactError, LinkContactPoint};
use fs_ga::Vec3;
use fs_math::dd::Dd;

/// One explicitly configured spherical collision surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SceneSphere {
    /// The existing whole-tree balanced carrier, with its original enclosure,
    /// inertia and free/prescribed motion admission. Only one per system.
    Carrier {
        system_index: usize,
        sphere: ContactSphere,
    },
    /// A spherical surface fixed to one link of a nonrotating prescribed
    /// fixed/prismatic tree. point_body is the sphere CENTRE in link axes.
    /// The radius is collision geometry, not an inferred enclosure of the link
    /// inertia or any omitted link shape. Native COM offsets remain unchanged.
    TranslatingLink {
        centre: LinkContactPoint,
        sphere: ContactSphere,
    },
}

/// Advance declared native link spheres and original balanced carriers.
///
/// Every system must have configured geometry and must use exactly one mode:
/// one original Carrier, or one or more TranslatingLink surfaces. The latter
/// admits arbitrary fixed link offsets/rotations, native COM offsets, full
/// inertia tensors and branched or serial fixed/prismatic trees. Its prescribed
/// base must have exactly zero angular velocity and angular acceleration;
/// translation has constant world acceleration. Held original joint efforts,
/// body-frame wrenches and world gravity retain their native meanings.
///
/// Complete pair coverage is required over the declared spheres, including
/// distinct moving links of the SAME native tree. To avoid inventing internal
/// collision between welded portions of one rigid body, at most one sphere
/// may be attached to a connected component of fixed joints. Carrier and link
/// geometry cannot overlap on one system. At most seven spheres, 21 pairs and
/// 21 planes use the existing bounded shared event loop.
///
/// All configured sphere pairs and sphere/plane pairs have explicit laws.
/// Omitted link shapes are not checked or inferred to fit a declared sphere.
/// Native ABA supplies coupled link-point accelerations; the full quadratic
/// centre paths and actual native endpoints/material witnesses must pass the
/// original geometry owners. Native initial, midpoint-predictor and endpoint
/// joint admission also bounds every coordinate extremum: an interior turning
/// point at t* has q(t*)=q0+v0*t*/2, between q0 and the native predictor
/// q0+v0*h/2. No additional configuration is integrated or projected.
/// A compensated fixed-axis reconstruction must agree with native Cartesian
/// acceleration at the physical acceleration scale; unresolved subtraction
/// of large gravity or support terms refuses before collision discovery.
///
/// Impacts retain the original simultaneous native mass coupling, joint rates,
/// full prescribed-support impulses and signed plane/base/actuator work.
/// Every successful result can continue from its actual states, clock and
/// rebased planes. A late geometry, joint-limit, solver, budget or cancellation
/// refusal publishes no partial state.
///
/// This is floating analytic coverage of the CONFIGURED surfaces. It is not
/// interval-certified CCD or support for rotating offset links, revolute or
/// helical link paths, sustained contact, arbitrary support commands or heat
/// partition. Those require their own motion/contact owners.
#[allow(clippy::too_many_arguments)]
pub fn advance_translating_sphere_scene(
    systems: &[System<'_>],
    geometry: &[SceneSphere],
    pair_laws: &[SpherePairLaw],
    planes: &[AcceleratingPlane],
    plane_laws: &[SpherePlaneLaw],
    start_time_s: f64,
    controls: MixedSphereControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixedSphereSceneResult, MixedSphereError> {
    poll(&mut cancelled)?;
    let prepared = PreparedGeometry::translating(
        systems,
        geometry,
        pair_laws,
        start_time_s,
        controls,
        !planes.is_empty(),
    )?;
    super::advance_prepared(
        systems,
        prepared,
        planes,
        plane_laws,
        start_time_s,
        controls,
        cancelled,
    )
}

#[derive(Clone, Copy)]
pub(super) struct SphereAttachment {
    pub(super) centre: LinkContactPoint,
}

#[derive(Clone, Copy, Default)]
pub(super) struct CentrePath {
    pub(super) position_world: Vec3,
    pub(super) velocity_world: Vec3,
    pub(super) acceleration_world: Vec3,
}

pub(super) struct NativeGeometry {
    pub(super) centres: Vec<CentrePath>,
    /// Complete original system order during an event; empty otherwise.
    pub(super) kinematics: Vec<Kinematics>,
}

pub(super) struct PreparedGeometry {
    pub(super) spheres: Vec<ContactSphere>,
    pub(super) pairs: Vec<Pair>,
    pub(super) commands: Vec<CentreCommand>,
    pub(super) attachments: Vec<SphereAttachment>,
    translating: Vec<bool>,
}

impl PreparedGeometry {
    pub(super) fn carriers(
        systems: &[System<'_>],
        spheres: &[ContactSphere],
        laws: &[SpherePairLaw],
        clock: f64,
        controls: MixedSphereControls,
        allow_single: bool,
    ) -> Result<Self, MixedSphereError> {
        let (pairs, commands) =
            super::admit(systems, spheres, laws, clock, controls, allow_single)?;
        Ok(Self {
            spheres: spheres.to_vec(),
            pairs,
            commands,
            attachments: (0..systems.len())
                .map(|system_index| SphereAttachment {
                    centre: LinkContactPoint {
                        system_index,
                        link_index: 0,
                        point_body: Vec3::default(),
                    },
                })
                .collect(),
            translating: vec![false; systems.len()],
        })
    }

    fn translating(
        systems: &[System<'_>],
        geometry: &[SceneSphere],
        laws: &[SpherePairLaw],
        clock: f64,
        controls: MixedSphereControls,
        allow_single: bool,
    ) -> Result<Self, MixedSphereError> {
        admit_window(systems, clock, controls)?;
        if geometry.is_empty() || geometry.len() > 7 || (geometry.len() < 2 && !allow_single) {
            return Err(MixedSphereError::Invalid(
                "bounded pair or plane geometry required",
            ));
        }
        let mut spheres = Vec::with_capacity(geometry.len());
        let mut attachments = Vec::with_capacity(geometry.len());
        let mut modes = vec![None; systems.len()];
        let mut carriers = vec![None; systems.len()];
        for entry in geometry {
            let (sphere, centre, translating) = match *entry {
                SceneSphere::Carrier {
                    system_index,
                    sphere,
                } => (
                    sphere,
                    LinkContactPoint {
                        system_index,
                        link_index: 0,
                        point_body: Vec3::default(),
                    },
                    false,
                ),
                SceneSphere::TranslatingLink { centre, sphere } => (sphere, centre, true),
            };
            if centre.system_index >= systems.len()
                || centre.link_index >= systems[centre.system_index].model.link_count()
                || !finite_vector(centre.point_body)
            {
                return Err(MixedSphereError::Invalid("invalid sphere link attachment"));
            }
            let index = centre.system_index;
            if modes[index].is_some_and(|mode| mode != translating)
                || (!translating && carriers[index].is_some())
            {
                return Err(MixedSphereError::Invalid(
                    "one geometry mode per native system required",
                ));
            }
            modes[index] = Some(translating);
            if !translating {
                carriers[index] = Some(sphere);
            }
            spheres.push(sphere);
            attachments.push(SphereAttachment { centre });
        }
        let pairs = admit_pairs(&spheres, laws, controls)?;
        let mut commands = Vec::with_capacity(systems.len());
        for (i, system) in systems.iter().enumerate() {
            match modes[i] {
                None => {
                    return Err(MixedSphereError::Invalid(
                        "every native system needs declared geometry",
                    ));
                }
                Some(false) => commands.push(admit_carrier(
                    system,
                    carriers[i].ok_or(MixedSphereError::Invalid("missing carrier geometry"))?,
                )?),
                Some(true) => {
                    let State::Prescribed(state) = system.state else {
                        return Err(MixedSphereError::Invalid(
                            "link spheres require a prescribed translating tree",
                        ));
                    };
                    if state.base.twist_body.angular != Vec3::default()
                        || state.base.acceleration_body.angular != Vec3::default()
                        || system.model.links().iter().any(|link| {
                            !matches!(link.joint().kind(), JointKind::Fixed | JointKind::Prismatic)
                        })
                    {
                        return Err(MixedSphereError::Invalid(
                            "nonrotating fixed/prismatic link paths required",
                        ));
                    }
                    if system.loads.external_body_wrench.len() != system.model.link_count()
                        || system.loads.generalized_force.len() != system.model.dof_count()
                    {
                        return Err(MixedSphereError::Invalid(
                            "native load dimensions must match model",
                        ));
                    }
                    let acceleration_world = state
                        .base
                        .world_from_base
                        .rotation()
                        .rotate(state.base.acceleration_body.linear)?;
                    if !finite_vector(acceleration_world) {
                        return Err(MixedSphereError::Invalid(
                            "nonfinite translating support acceleration",
                        ));
                    }
                    commands.push(CentreCommand {
                        acceleration_world,
                        angular_velocity_world: Vec3::default(),
                    });
                }
            }
        }
        // Fixed-connected links are one rigid component, not independently
        // colliding bodies. Every remaining declared pair is actually checked.
        for (i, first) in attachments.iter().enumerate() {
            for second in &attachments[i + 1..] {
                if first.centre.system_index == second.centre.system_index {
                    let model = systems[first.centre.system_index].model;
                    let component = |mut link: usize| {
                        while model.links()[link].joint().kind() == JointKind::Fixed {
                            let Some(parent) = model.links()[link].parent() else {
                                break;
                            };
                            link = parent;
                        }
                        link
                    };
                    if component(first.centre.link_index) == component(second.centre.link_index) {
                        return Err(MixedSphereError::Invalid(
                            "one sphere per fixed-joint rigid component required",
                        ));
                    }
                }
            }
        }
        Ok(Self {
            spheres,
            pairs,
            commands,
            attachments,
            translating: modes.into_iter().map(|mode| mode == Some(true)).collect(),
        })
    }

    /// Native physical accelerations are evaluated once per accepted segment.
    /// Constant-mass, fixed-orientation admission makes them complete paths,
    /// rather than a Taylor guess about general articulated motion.
    pub(super) fn paths(
        &self,
        systems: &[System<'_>],
        states: &[State],
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<CentrePath>, MixedSphereError> {
        let mut centres = vec![CentrePath::default(); self.spheres.len()];
        for (i, (system, state)) in systems.iter().zip(states).enumerate() {
            poll(cancelled)?;
            if self.translating[i] {
                let State::Prescribed(state) = state else {
                    return Err(MixedSphereError::Invalid("native boundary kind changed"));
                };
                let kinematics = forward_kinematics(
                    system.model,
                    state.base,
                    &state.joint_position,
                    &state.joint_velocity,
                )
                .map_err(|error| native_error(i, time, error))?;
                poll(cancelled)?;
                let dynamics = forward_dynamics(
                    system.model,
                    state.base,
                    &state.joint_position,
                    &state.joint_velocity,
                    system.loads.generalized_force,
                    system.loads.gravity_world,
                    system.loads.external_body_wrench,
                )
                .map_err(|error| native_error(i, time, error))?;
                poll(cancelled)?;
                for (slot, attachment) in self.attachments.iter().enumerate() {
                    let point = attachment.centre;
                    if point.system_index != i {
                        continue;
                    }
                    let motion = dynamics
                        .prepare_point_motion(
                            &kinematics,
                            point.link_index,
                            system.loads.gravity_world,
                        )
                        .and_then(|prepared| prepared.point_motion(point.point_body))
                        .map_err(|error| native_error(i, time, error))?;
                    if motion.angular_velocity_world != Vec3::default()
                        || motion.angular_acceleration_world != Vec3::default()
                    {
                        return Err(MixedSphereError::Invalid(
                            "native link left translating motion class",
                        ));
                    }
                    self.admit_acceleration(
                        system,
                        &kinematics,
                        &dynamics.generalized_acceleration,
                        point,
                        motion.point_acceleration_world,
                        time,
                        cancelled,
                    )?;
                    centres[slot] = CentrePath {
                        position_world: motion.point_world,
                        velocity_world: motion.point_velocity_world,
                        acceleration_world: motion.point_acceleration_world,
                    };
                }
            } else {
                self.carrier_path(i, state, &mut centres)?;
            }
        }
        poll(cancelled)?;
        Ok(centres)
    }

    /// The admitted FK map is affine in q. Pair native qdd with those original
    /// fixed world axes independently of ABA's gravity-shifted body vector.
    /// This is an admission check, never a replacement acceleration or state.
    #[allow(clippy::too_many_arguments)]
    fn admit_acceleration(
        &self,
        system: &System<'_>,
        kinematics: &Kinematics,
        acceleration: &[f64],
        point: LinkContactPoint,
        native: Vec3,
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MixedSphereError> {
        let base = self.commands[point.system_index].acceleration_world;
        let mut expected = [base.x, base.y, base.z].map(Dd::from_f64);
        let mut absolute = [base.x.abs(), base.y.abs(), base.z.abs()].map(Dd::from_f64);
        let mut link = point.link_index;
        let mut count = 1_usize;
        loop {
            if count % 64 == 1 {
                poll(cancelled)?;
            }
            if let Some(dof) = system.model.dof_index(link) {
                let axis = kinematics.world_from_link[link]
                    .rotation()
                    .rotate(system.model.links()[link].joint().motion_subspace().linear)?;
                for (i, coefficient) in [axis.x, axis.y, axis.z].into_iter().enumerate() {
                    let term = Dd::from_f64(coefficient) * Dd::from_f64(acceleration[dof]);
                    expected[i] = expected[i] + term;
                    absolute[i] = absolute[i]
                        + Dd::from_f64(coefficient.abs()) * Dd::from_f64(acceleration[dof].abs());
                }
            }
            count += 1;
            let Some(parent) = system.model.links()[link].parent() else {
                break;
            };
            link = parent;
        }
        let native = [native.x, native.y, native.z];
        let differences: [f64; 3] =
            core::array::from_fn(|i| (Dd::from_f64(native[i]) - expected[i]).to_f64().abs());
        let physical = expected.map(|value| value.to_f64().abs());
        let original_terms = absolute.map(Dd::to_f64);
        if !differences
            .into_iter()
            .chain(physical)
            .chain(original_terms)
            .all(f64::is_finite)
        {
            return Err(MixedSphereError::UnresolvedLinkPath {
                system_index: point.system_index,
                link_index: point.link_index,
                time,
                acceleration_residual: f64::INFINITY,
                allowed_acceleration_residual: 0.0,
            });
        }
        let residual = differences.into_iter().fold(0.0, f64::max);
        let scale = physical.into_iter().fold(0.0, f64::max);
        let terms = original_terms.into_iter().fold(0.0, f64::max);
        // The ordinary rounding allowance follows the resulting physical
        // acceleration, not a huge cancelled gravity/background. Original-term
        // DD arithmetic has its separate second-order roundoff allowance.
        let allowed = (256.0 * f64::EPSILON * count as f64) * scale
            + (256.0 * f64::EPSILON * f64::EPSILON * count as f64) * terms;
        if !residual.is_finite() || !allowed.is_finite() || residual > allowed {
            return Err(MixedSphereError::UnresolvedLinkPath {
                system_index: point.system_index,
                link_index: point.link_index,
                time,
                acceleration_residual: residual,
                allowed_acceleration_residual: allowed,
            });
        }
        Ok(())
    }

    fn carrier_path(
        &self,
        system: usize,
        state: &State,
        centres: &mut [CentrePath],
    ) -> Result<(), MixedSphereError> {
        for (slot, attachment) in self.attachments.iter().enumerate() {
            if attachment.centre.system_index == system {
                centres[slot] = CentrePath {
                    position_world: base(state).world_from_base.translation(),
                    velocity_world: centre_velocity(state)?,
                    acceleration_world: self.commands[system].acceleration_world,
                };
            }
        }
        Ok(())
    }

    /// Reconstruct real endpoints once; event FK is reused by every witness.
    pub(super) fn snapshot(
        &self,
        systems: &[System<'_>],
        states: &[State],
        event: bool,
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NativeGeometry, MixedSphereError> {
        let mut centres = vec![CentrePath::default(); self.spheres.len()];
        let mut kinematics = Vec::with_capacity(if event { systems.len() } else { 0 });
        for (i, (system, state)) in systems.iter().zip(states).enumerate() {
            poll(cancelled)?;
            if event || self.translating[i] {
                let input = state.contact_system(system.model);
                let native = forward_kinematics(
                    system.model,
                    base(state),
                    input.joint_position,
                    input.joint_velocity,
                )
                .map_err(|error| native_error(i, time, error))?;
                if self.translating[i] {
                    for (slot, attachment) in self.attachments.iter().enumerate() {
                        let point = attachment.centre;
                        if point.system_index != i {
                            continue;
                        }
                        let pose = native.world_from_link[point.link_index];
                        let twist = native.body_twist[point.link_index];
                        if twist.angular != Vec3::default() {
                            return Err(MixedSphereError::Invalid(
                                "native endpoint left translating motion class",
                            ));
                        }
                        centres[slot] = CentrePath {
                            position_world: pose.transform_point(point.point_body)?,
                            velocity_world: pose.rotation().rotate(twist.linear)?,
                            // Only the original zero-duration endpoint owner
                            // sees this; each next segment recomputes native ABA.
                            acceleration_world: Vec3::default(),
                        };
                    }
                }
                if event {
                    kinematics.push(native);
                }
            }
            if !self.translating[i] {
                self.carrier_path(i, state, &mut centres)?;
            }
        }
        if centres
            .iter()
            .any(|path| !finite_vector(path.position_world) || !finite_vector(path.velocity_world))
        {
            return Err(MixedSphereError::Invalid(
                "nonfinite native sphere endpoint",
            ));
        }
        poll(cancelled)?;
        Ok(NativeGeometry {
            centres,
            kinematics,
        })
    }
}

fn native_error(system_index: usize, time: f64, error: ArticulatedError) -> MixedSphereError {
    ArticulatedTrajectoryError::Contact {
        event_index: 0,
        time,
        error: ArticulatedContactError::Dynamics {
            system_index,
            error,
        },
    }
    .into()
}
