//! Complete mixed spherical scenes against independently accelerating planes.
//!
//! Every plane keeps its normal fixed and translates along a quadratic world
//! path. Both normal approach and tangential conveyor speed use the actual
//! instantaneous plane motion. The parent supplies the original native smooth
//! dynamics, one simultaneous body/plane impact and DD window accounting.

use super::{
    CentreCommand, ContactSphere, MixedArticulatedTrajectoryResult, MixedSphereControls,
    MixedSphereError, Observation, SpherePairLaw, State, System, Work, base, centre_velocity,
    finite_vector, poll,
};
use crate::articulated::Kinematics;
use crate::articulated::contact_geometry::prepare_contact_frame;
use crate::articulated_contacts::friction::ArticulatedFrictionalContact;
use crate::articulated_contacts::trajectory::ArticulatedTrajectoryError;
use crate::articulated_contacts::{ContactCounterpart, LinkContactPoint, unit_normal};
use fs_contact::ballistic::{BallisticContactInput, first_ballistic_contact};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_ga::Vec3;

/// An infinite plane with constant world translation acceleration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcceleratingPlane {
    /// Unique plane identity, independent of the sphere ID namespace.
    pub id: u64,
    /// A point on the plane at the start of this call, in metres.
    pub point_world: Vec3,
    /// World material velocity at the start, including conveyor motion, m/s.
    pub velocity_world: Vec3,
    /// Constant world material acceleration, including tangent components, m/s².
    pub acceleration_world: Vec3,
    /// Nonzero fixed normal into the permitted half-space.
    pub normal_world: Vec3,
    /// Nonparallel tangent hint for the native circular Coulomb law.
    pub tangent_world: Vec3,
}
impl AcceleratingPlane {
    /// Rebase the same plane command to a later window without advancing twice.
    /// No rotation or change to the original normal/tangent hints is introduced.
    pub fn at_elapsed(self, elapsed_s: f64) -> Result<Self, MixedSphereError> {
        let sample = self.sample(elapsed_s)?;
        Ok(Self {
            point_world: sample.point,
            velocity_world: sample.velocity,
            ..self
        })
    }

    fn sample(self, time: f64) -> Result<PlaneSample, MixedSphereError> {
        if !time.is_finite() || time < 0.0 {
            return Err(MixedSphereError::Invalid(
                "finite nonnegative plane time required",
            ));
        }
        let point = self.point_world
            + self.velocity_world.scale(time)
            + self.acceleration_world.scale(time).scale(0.5 * time);
        let velocity = self.velocity_world + self.acceleration_world.scale(time);
        if !finite_vector(point) || !finite_vector(velocity) {
            return Err(MixedSphereError::Invalid(
                "unrepresentable accelerating plane state",
            ));
        }
        Ok(PlaneSample { point, velocity })
    }
}

/// Explicit material law for one sphere against one plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpherePlaneLaw {
    pub sphere: u64,
    pub plane: u64,
    pub friction_coefficient: f64,
    pub restitution: f64,
}

/// Stable geometry identity for one simultaneous native contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum MixedSphereSceneContactKey {
    SpherePair([u64; 2]),
    SpherePlane { sphere: u64, plane: u64 },
}

/// Material witnesses admitted together at one automatically found event.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedSphereSceneEventGeometry {
    /// Canonical body pairs first, then sphere/plane ID pairs; same order as the
    /// native event's contact_ids and contact-space vectors.
    pub contacts: Vec<MixedSphereSceneContactKey>,
    pub maximum_gap_m: f64,
}

/// Complete native window with reusable external plane commands.
#[derive(Debug, Clone)]
pub struct MixedSphereSceneResult {
    pub trajectory: MixedArticulatedTrajectoryResult,
    pub event_geometry: Vec<MixedSphereSceneEventGeometry>,
    /// Commands rebased to the accepted endpoint, in the original plane order.
    /// Use these with states_after and time_after_s for the next window.
    pub planes_after: Vec<AcceleratingPlane>,
    pub smooth_steps_attempted: usize,
    pub geometry_queries: usize,
    /// Minimum complete accepted-path, actual endpoint and material-witness gap.
    /// This is a floating analytic observation, not an interval enclosure.
    pub minimum_gap_m: f64,
}

/// Discover and advance simultaneous body-pair and accelerating-plane impacts.
///
/// All original mixed-carrier admission and native motion/impulse gates remain.
/// One sphere is sufficient when planes are supplied. Every unordered sphere
/// pair and every sphere/plane combination requires exactly one explicit law.
/// At most 21 planes and seven spheres are admitted; the original active-contact
/// cap remains 21 and refuses rather than dropping an oversized contact set.
/// The shared geometry-query budget charges every plane and body-pair query.
///
/// A plane's complete normal gap is quadratic. The existing ballistic owner
/// finds first closing contact; complete path checks include an interior
/// minimum and isolated grazing. Actual native endpoints and root material
/// witnesses must also remain nonpenetrating. Nearby body and plane features
/// enter ONE native simultaneous solve, using the plane's event-time material
/// velocity in signed support work. A plane's acceleration is never inserted
/// into a body impulse or confused with its instantaneous velocity.
///
/// Input planes belong to this window's start, independent of the native clock.
/// planes_after explicitly rebases them for continuation. Repeated impacts,
/// bounded event-time retries, actual native momentum/energy increments, and
/// cancellation remain with the parent's shared event loop. Sustained contact,
/// rotating planes, arbitrary articulated offsets, interval-certified CCD and
/// thermal allocation require other owners. All failure paths publish no state.
#[allow(clippy::too_many_arguments)]
pub fn advance_mixed_sphere_scene(
    systems: &[System<'_>],
    spheres: &[ContactSphere],
    sphere_pair_laws: &[SpherePairLaw],
    planes: &[AcceleratingPlane],
    plane_laws: &[SpherePlaneLaw],
    start_time_s: f64,
    controls: MixedSphereControls,
    cancelled: impl FnMut() -> bool,
) -> Result<MixedSphereSceneResult, MixedSphereError> {
    super::advance_scene(
        systems,
        spheres,
        sphere_pair_laws,
        planes,
        plane_laws,
        start_time_s,
        controls,
        cancelled,
    )
}

#[derive(Clone, Copy)]
struct PlaneSample {
    point: Vec3,
    velocity: Vec3,
}

#[derive(Clone, Copy)]
struct PlanePair {
    sphere: usize,
    plane: usize,
    id: u64,
    law: SpherePlaneLaw,
}

pub(super) struct PreparedPlanes {
    planes: Vec<AcceleratingPlane>,
    normals: Vec<Vec3>,
    pairs: Vec<PlanePair>,
}
impl PreparedPlanes {
    pub(super) fn admit(
        spheres: &[ContactSphere],
        planes: &[AcceleratingPlane],
        laws: &[SpherePlaneLaw],
        body_pair_count: usize,
    ) -> Result<Self, MixedSphereError> {
        if planes.len() > 21 || laws.len() != spheres.len() * planes.len() {
            return Err(MixedSphereError::Invalid(
                "complete bounded sphere/plane laws required",
            ));
        }
        let mut ids: Vec<_> = planes.iter().map(|plane| plane.id).collect();
        ids.sort_unstable();
        if ids.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(MixedSphereError::Invalid("duplicate plane identity"));
        }
        let mut normals = Vec::with_capacity(planes.len());
        for plane in planes {
            if !finite_vector(plane.point_world)
                || !finite_vector(plane.velocity_world)
                || !finite_vector(plane.acceleration_world)
            {
                return Err(MixedSphereError::Invalid("finite plane motion required"));
            }
            let normal = unit_normal(plane.normal_world)
                .map_err(|_| MixedSphereError::Invalid("invalid plane normal"))?;
            let tangent = unit_normal(plane.tangent_world)
                .map_err(|_| MixedSphereError::Invalid("invalid plane tangent"))?;
            unit_normal(normal.cross(tangent))
                .map_err(|_| MixedSphereError::Invalid("parallel plane tangent"))?;
            normals.push(normal);
        }
        let mut ordered = laws.to_vec();
        ordered.sort_unstable_by_key(|law| (law.sphere, law.plane));
        let mut pairs = Vec::with_capacity(ordered.len());
        for (i, law) in ordered.into_iter().enumerate() {
            if !law.friction_coefficient.is_finite()
                || law.friction_coefficient < 0.0
                || !law.restitution.is_finite()
                || !(0.0..=1.0).contains(&law.restitution)
                || pairs.last().is_some_and(|pair: &PlanePair| {
                    (pair.law.sphere, pair.law.plane) == (law.sphere, law.plane)
                })
            {
                return Err(MixedSphereError::Invalid(
                    "invalid or duplicate sphere/plane law",
                ));
            }
            let sphere = spheres
                .iter()
                .position(|sphere| sphere.id == law.sphere)
                .ok_or(MixedSphereError::Invalid(
                    "plane law references unknown sphere",
                ))?;
            let plane = planes
                .iter()
                .position(|plane| plane.id == law.plane)
                .ok_or(MixedSphereError::Invalid(
                    "plane law references unknown plane",
                ))?;
            pairs.push(PlanePair {
                sphere,
                plane,
                id: (body_pair_count + i) as u64,
                law,
            });
        }
        Ok(Self {
            planes: planes.to_vec(),
            normals,
            pairs,
        })
    }

    fn samples(&self, time: f64) -> Result<Vec<PlaneSample>, MixedSphereError> {
        self.planes.iter().map(|plane| plane.sample(time)).collect()
    }

    pub(super) fn endpoints(
        &self,
        elapsed: f64,
    ) -> Result<Vec<AcceleratingPlane>, MixedSphereError> {
        self.planes
            .iter()
            .map(|plane| plane.at_elapsed(elapsed))
            .collect()
    }

    fn input(
        &self,
        pair: PlanePair,
        sphere: ContactSphere,
        state: &State,
        command: CentreCommand,
        sample: PlaneSample,
        duration: f64,
    ) -> Result<BallisticContactInput, MixedSphereError> {
        let normal = self.normals[pair.plane];
        Ok(BallisticContactInput {
            gap: normal.dot(base(state).world_from_base.translation() - sample.point)
                - sphere.radius_m,
            relative_normal_speed: normal.dot(centre_velocity(state)? - sample.velocity),
            relative_normal_acceleration: normal
                .dot(command.acceleration_world - self.planes[pair.plane].acceleration_world),
            duration_s: duration,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn encounter(
        &self,
        spheres: &[ContactSphere],
        states: &[State],
        commands: &[CentreCommand],
        elapsed: f64,
        horizon: f64,
        controls: MixedSphereControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<(f64, u64)>, MixedSphereError> {
        let samples = self.samples(elapsed)?;
        let mut earliest: Option<(f64, u64)> = None;
        for &pair in &self.pairs {
            poll(cancelled)?;
            work.query(controls)?;
            let input = self.input(
                pair,
                spheres[pair.sphere],
                &states[pair.sphere],
                commands[pair.sphere],
                samples[pair.plane],
                horizon,
            )?;
            let hit =
                first_ballistic_contact(input).map_err(|error| MixedSphereError::PlanePath {
                    sphere: pair.law.sphere,
                    plane: pair.law.plane,
                    time: elapsed,
                    error,
                })?;
            if let Some(hit) = hit {
                if hit.relative_normal_speed < 0.0 || hit.time_s == 0.0 {
                    if hit.relative_normal_speed.abs()
                        <= controls.trajectory.contact.velocity_tolerance
                    {
                        return Err(MixedSphereError::UnresolvedPlaneContact {
                            sphere: pair.law.sphere,
                            plane: pair.law.plane,
                            time: elapsed + hit.time_s,
                        });
                    }
                    if earliest.is_none_or(|(time, _)| hit.time_s < time) {
                        earliest = Some((hit.time_s, pair.id));
                    }
                }
                // An isolated future tangent touch has zero impulse. Its
                // complete nonpenetrating path is still checked below.
            }
        }
        Ok(earliest)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn observe(
        &self,
        spheres: &[ContactSphere],
        before: &[State],
        after: &[State],
        commands: &[CentreCommand],
        kinematics: &[Kinematics],
        duration: f64,
        elapsed: f64,
        event: bool,
        controls: MixedSphereControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
        observation: &mut Observation,
    ) -> Result<bool, MixedSphereError> {
        let first = self.samples(elapsed)?;
        let last = self.samples(elapsed + duration)?;
        for &pair in &self.pairs {
            poll(cancelled)?;
            work.query(controls)?;
            let input = self.input(
                pair,
                spheres[pair.sphere],
                &before[pair.sphere],
                commands[pair.sphere],
                first[pair.plane],
                duration,
            )?;
            let at = |time: f64| {
                (0.5 * input.relative_normal_acceleration)
                    .mul_add(time, input.relative_normal_speed)
                    .mul_add(time, input.gap)
            };
            let mut path_gap = input.gap.min(at(duration));
            if input.relative_normal_acceleration > 0.0 {
                let t = -input.relative_normal_speed / input.relative_normal_acceleration;
                if t > 0.0 && t < duration {
                    path_gap = path_gap.min(at(t));
                }
            }
            let normal = self.normals[pair.plane];
            let actual_gap = normal.dot(
                base(&after[pair.sphere]).world_from_base.translation() - last[pair.plane].point,
            ) - spheres[pair.sphere].radius_m;
            if !path_gap.is_finite() || !actual_gap.is_finite() {
                return Err(MixedSphereError::Invalid(
                    "nonfinite complete sphere/plane path gap",
                ));
            }
            observation.minimum_gap = observation.minimum_gap.min(path_gap).min(actual_gap);
            if observation.minimum_gap < 0.0 {
                return Ok(false);
            }
            if event && actual_gap <= controls.event_gap_tolerance_m {
                work.query(controls)?;
                let root = kinematics[pair.sphere].world_from_link[0];
                let arm = root
                    .rotation()
                    .inverse()
                    .rotate(normal.scale(-spheres[pair.sphere].radius_m))?;
                let contact = ArticulatedFrictionalContact {
                    id: pair.id,
                    point: LinkContactPoint {
                        system_index: pair.sphere,
                        link_index: 0,
                        point_body: arm,
                    },
                    other: ContactCounterpart::PrescribedPlane {
                        point_world: last[pair.plane].point,
                        velocity_world: last[pair.plane].velocity,
                    },
                    normal_world: self.planes[pair.plane].normal_world,
                    tangent_world: self.planes[pair.plane].tangent_world,
                    activation_tolerance: controls.event_gap_tolerance_m,
                    friction_coefficient: pair.law.friction_coefficient,
                    normal_law: NormalVelocityLaw::NewtonRestitution {
                        coefficient: pair.law.restitution,
                    },
                };
                let frame = prepare_contact_frame(
                    contact.point,
                    contact.other,
                    contact.normal_world,
                    contact.tangent_world,
                    pair.id as usize,
                    kinematics,
                )
                .map_err(|error| ArticulatedTrajectoryError::Contact {
                    event_index: 0,
                    time: elapsed + duration,
                    error: error.into(),
                })?;
                let gap = frame.rows[0].gap;
                observation.minimum_gap = observation.minimum_gap.min(gap);
                if gap < 0.0 {
                    return Ok(false);
                }
                if gap <= controls.event_gap_tolerance_m {
                    observation.contacts.push(contact);
                    observation
                        .keys
                        .push(MixedSphereSceneContactKey::SpherePlane {
                            sphere: pair.law.sphere,
                            plane: pair.law.plane,
                        });
                    observation.maximum_gap = observation.maximum_gap.max(actual_gap).max(gap);
                }
            }
        }
        poll(cancelled)?;
        Ok(true)
    }
}
