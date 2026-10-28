//! Automatic trajectories of finite native spheres under constant centre loads.
//!
//! Common acceleration cancels from each relative centre path; body-specific
//! forces give a quadratic relative curve. The existing contact owners find
//! its first entry and all distance extrema. All touching pairs then enter one
//! simultaneous normal/Coulomb event. Native accelerated drift, free rotation
//! and point impulses keep their original inertia and audits.
//! This floating analytic lane is not interval-certified general-shape CCD.

use super::{
    RigidFrictionalContact, RigidFrictionalImpactError, RigidFrictionalImpactReport,
    RigidImpactControls, RigidImpactCounterpart, RigidImpactPoint, apply_rigid_frictional_impacts,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_contact::swept_sphere::accelerated::{AcceleratedSphereInput, sphere_path};
use fs_contact::swept_sphere::{SweptSphereError, SweptSphereInput, first_sphere_contact};
use fs_exec::CancelGate;

/// A complete spherical boundary centred at the native centre of mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactSphere {
    pub id: u64,
    pub radius_m: f64,
}

/// Explicit, source-resolved isotropic law for one ordered sphere pair.
/// Every unordered pair must appear exactly once, with first < second.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePairLaw {
    pub first: u64,
    pub second: u64,
    pub friction_coefficient: f64,
    pub restitution: f64,
}

/// Bounded event work and independent native-physics gates.
#[derive(Clone, Copy, Debug)]
pub struct SphereContactControls {
    pub impact: RigidImpactControls,
    pub maximum_pairs: usize,
    pub maximum_events: usize,
    /// Maximum outward event bracket in metres. No negative gap is accepted.
    /// Nearby contacts enter the same event at this declared spatial resolution.
    pub event_gap_tolerance_m: f64,
    pub drift_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
}

/// One simultaneous event, with canonical stable sphere identities.
#[derive(Clone, Debug)]
pub struct SphereContactEvent {
    pub elapsed_s: f64,
    pub pairs: Vec<[u64; 2]>,
    pub maximum_gap_m: f64,
    pub impact: RigidFrictionalImpactReport,
}

/// Whole trajectory receipt; impact loss is deliberately unallocated.
#[derive(Clone, Debug)]
pub struct SphereContactReport {
    pub duration_s: f64,
    pub events: Vec<SphereContactEvent>,
    pub pair_queries: usize,
    pub minimum_path_gap_m: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Signed work of the original constant gravitational forces along the
    /// actual accepted centre displacements. No impact loss is added to it.
    pub gravity_work_j: f64,
    pub gravity_impulse_ns: Vec3,
    /// Integrated gravitational torque about the world origin, evaluated from
    /// the original total mass moment, momentum and external centre forces.
    pub gravity_angular_impulse_nms: Vec3,
    /// Signed work of the additional body-specific constant centre forces.
    pub centre_force_work_j: f64,
    pub centre_force_impulse_ns: Vec3,
    /// Exact piecewise quadratic orbital torque integral of the centre forces.
    /// Individual collision impulses change this integral between events.
    pub centre_force_angular_impulse_nms: Vec3,
    pub unallocated_impact_loss_j: f64,
    /// Measured native drift kinetic change minus all applied force work.
    pub drift_energy_defect_j: f64,
    pub total_energy_residual_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// Every refusal leaves the supplied bodies unchanged.
#[derive(Debug)]
pub enum SphereContactError {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Detection(SweptSphereError),
    Impact(RigidFrictionalImpactError),
    Penetration {
        pair: [u64; 2],
        gap_m: f64,
    },
    EventBudget,
    /// Root, native witness, or floating clock cannot resolve the next event.
    Unresolved(&'static str),
    Cancelled,
    Audit(&'static str),
}
impl From<DynamicsError> for SphereContactError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl From<SweptSphereError> for SphereContactError {
    fn from(e: SweptSphereError) -> Self {
        Self::Detection(e)
    }
}
impl From<RigidFrictionalImpactError> for SphereContactError {
    fn from(e: RigidFrictionalImpactError) -> Self {
        Self::Impact(e)
    }
}
impl core::fmt::Display for SphereContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native sphere contact: {self:?}")
    }
}
impl std::error::Error for SphereContactError {}

#[derive(Clone, Copy)]
struct Pair {
    a: usize,
    b: usize,
    ids: [u64; 2],
    radius: f64,
    law: SpherePairLaw,
}

/// Advance all native bodies through automatically discovered sphere collisions.
///
/// The declared spheres are the COMPLETE collision geometry. Their centres are
/// the original COMs; anisotropic principal inertia is allowed when it fits
/// inside the sphere. Every pair has an explicit immutable law. No caller
/// supplies witnesses, impact times, contact normals, or path-admission code.
///
/// Pair enumeration and impulse accumulation use stable sphere IDs. Each free
/// segment checks the analytic interior minimum separation as well as actual
/// endpoints. The event time may move outward by at most 64 adjacent floating
/// values to resolve native nonnegative witness gaps; no body is projected.
/// The actual accepted gaps are recorded and bounded by event_gap_tolerance_m.
/// Initial penetration, missed or unresolved events and repeated zero-time
/// impulses refuse atomically. Events closer than the declared spatial
/// resolution are treated simultaneously, never by arbitrary pair priority.
///
/// Every body receives the same declared constant gravitational acceleration;
/// its force remains its original mass times gravity. The complete relative
/// path is linear between events, while each actual centre follows its native
/// accelerated trajectory. Gravity::ZERO explicitly selects free flight.
///
/// Gravity work and its total linear/angular impulse are independently checked
/// across the complete trajectory. Internal impulses cancel in the total COM
/// equation, so its gravitational torque integral is analytic even across
/// multiple frictional events. No prescribed support, differential force,
/// rolling/spin resistance, finite-patch history, deformation, or thermal
/// allocation is inferred. Native free attitude remains approximate for
/// asymmetric bodies, with its original checks and refinement requirement.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_sphere_contacts(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    gravity: Gravity,
    duration_s: f64,
    controls: SphereContactControls,
    gate: &CancelGate,
) -> Result<SphereContactReport, SphereContactError> {
    if bodies.len() > controls.impact.maximum_bodies {
        return Err(SphereContactError::Invalid("body work bound"));
    }
    advance_forced_sphere_contacts(
        bodies,
        properties,
        spheres,
        laws,
        gravity,
        &vec![Vec3::ZERO; bodies.len()],
        duration_s,
        controls,
        gate,
    )
}

/// Advance original native spheres under uniform gravity and an explicit
/// constant world force at EACH original centre of mass, in body-slice order.
///
/// Unlike common gravity, different force/mass ratios curve the relative path.
/// The original contact owner isolates all distance extrema on each quadratic
/// path and brackets the first entering contact. Equal accelerations retain the
/// cheaper straight relative-path query. No timesteps sample over collisions.
/// All touching pairs still enter one original simultaneous Coulomb event.
///
/// Force work uses actual accepted centre displacements; orbital angular
/// impulse integrates each accepted quadratic segment, including changes in
/// velocity at collisions. No force acts at an invented arm and no spin torque
/// is supplied. Persistent loaded contact, arbitrary shapes, force fields,
/// time-varying forcing, rolling resistance and impact heat allocation remain
/// outside this event-driven lane. An unresolved zero-speed loaded contact
/// refuses instead of tunneling or manufacturing a support impulse.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_forced_sphere_contacts(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    gravity: Gravity,
    centre_forces_world_n: &[Vec3],
    duration_s: f64,
    controls: SphereContactControls,
    gate: &CancelGate,
) -> Result<SphereContactReport, SphereContactError> {
    poll(gate)?;
    gravity.validate()?;
    let pairs = prepare(bodies, properties, spheres, laws, duration_s, controls)?;
    if centre_forces_world_n.len() != bodies.len()
        || centre_forces_world_n.iter().any(|f| !f.is_finite())
    {
        return Err(SphereContactError::Invalid(
            "one finite centre force per original body required",
        ));
    }
    let accelerations: Vec<_> = centre_forces_world_n
        .iter()
        .zip(properties)
        .map(|(&force, p)| force.scale(p.mass().recip()))
        .collect();
    if accelerations.iter().any(|a| !a.is_finite()) {
        return Err(SphereContactError::Invalid("unresolvable force/mass ratio"));
    }
    let (before, linear_before, angular_before) = super::totals(bodies, properties, gate)?;
    let mut staged = bodies.to_vec();
    let mut elapsed = 0.0;
    let mut events = Vec::new();
    let mut queries = 0usize;
    let mut minimum_gap = gaps(&staged, &pairs, gate)?;
    let mut loss = 0.0;
    let mut drift_defect = 0.0;
    let mut force_angular = Vec3::ZERO;
    loop {
        poll(gate)?;
        let remaining = duration_s - elapsed;
        let velocities = staged
            .iter()
            .zip(properties)
            .map(|(&b, &p)| b.center_of_mass_velocity_world(p))
            .collect::<Result<Vec<_>, _>>()?;
        let mut hit: Option<f64> = None;
        for pair in &pairs {
            poll(gate)?;
            queries = queries
                .checked_add(1)
                .ok_or(SphereContactError::EventBudget)?;
            if let Some(time) = first_hit(&staged, &velocities, &accelerations, *pair, remaining)? {
                hit = Some(hit.map_or(time, |old| old.min(time)));
            }
        }
        let Some(mut h) = hit else {
            if remaining > 0.0 {
                minimum_gap = minimum_gap.min(path_gap(
                    &staged,
                    &velocities,
                    &accelerations,
                    &pairs,
                    remaining,
                    gate,
                )?);
                let next = drift(
                    &staged,
                    properties,
                    gravity,
                    centre_forces_world_n,
                    remaining,
                    gate,
                )?;
                minimum_gap = minimum_gap.min(gaps(&next, &pairs, gate)?);
                drift_defect += drift_audit(
                    &staged,
                    &next,
                    properties,
                    gravity,
                    centre_forces_world_n,
                    remaining,
                    controls,
                    gate,
                )?;
                force_angular = force_angular.add(
                    centre_force_exchange(
                        &staged,
                        &next,
                        properties,
                        gravity,
                        centre_forces_world_n,
                        remaining,
                        gate,
                    )?
                    .angular,
                );
                staged = next;
            }
            break;
        };
        if events.len() >= controls.maximum_events {
            return Err(SphereContactError::EventBudget);
        }
        if h > 0.0 && elapsed + h == elapsed {
            return Err(SphereContactError::Unresolved("event clock cannot advance"));
        }
        // Resolve the actual native surface witnesses from the nonpenetrating
        // side of the analytic root. Root adjustment changes time, never pose.
        let mut candidate = None;
        for _ in 0..=64 {
            poll(gate)?;
            let next = if h == 0.0 {
                staged.clone()
            } else {
                drift(&staged, properties, gravity, centre_forces_world_n, h, gate)?
            };
            if let (Ok(gap), Ok(path_minimum)) = (
                gaps(&next, &pairs, gate),
                path_gap(&staged, &velocities, &accelerations, &pairs, h, gate),
            ) {
                let (contacts, ids, maximum_gap, witnesses_nonnegative) =
                    contacts(&next, spheres, &pairs, controls, gate)?;
                if !contacts.is_empty() && witnesses_nonnegative {
                    candidate = Some((next, contacts, ids, maximum_gap, gap, path_minimum));
                    break;
                }
            }
            if h == 0.0 {
                break;
            }
            h = h.next_down();
        }
        let Some((mut next, contacts, ids, maximum_gap, gap, path_minimum)) = candidate else {
            return Err(SphereContactError::Unresolved(
                "event geometry or outward root bracket",
            ));
        };
        minimum_gap = minimum_gap.min(path_minimum).min(gap);
        drift_defect += drift_audit(
            &staged,
            &next,
            properties,
            gravity,
            centre_forces_world_n,
            h,
            controls,
            gate,
        )?;
        force_angular = force_angular.add(
            centre_force_exchange(
                &staged,
                &next,
                properties,
                gravity,
                centre_forces_world_n,
                h,
                gate,
            )?
            .angular,
        );
        let impact = apply_rigid_frictional_impacts(
            &mut next,
            properties,
            &contacts,
            controls.impact,
            gate,
        )?;
        if impact
            .contact
            .impulses
            .iter()
            .all(|p| p.iter().all(|x| *x == 0.0))
        {
            return Err(SphereContactError::Unresolved(
                "closing geometry produced no resolving impulse",
            ));
        }
        loss += impact.contact.dissipated_energy;
        elapsed += h;
        if !elapsed.is_finite() || elapsed > duration_s {
            return Err(SphereContactError::Unresolved(
                "event time outside duration",
            ));
        }
        events.push(SphereContactEvent {
            elapsed_s: elapsed,
            pairs: ids,
            maximum_gap_m: maximum_gap,
            impact,
        });
        staged = next;
        if elapsed == duration_s {
            break;
        }
    }
    let (after, linear_after, angular_after) = super::totals(&staged, properties, gate)?;
    let supplied = gravity_exchange(
        bodies,
        &staged,
        properties,
        gravity,
        centre_forces_world_n,
        duration_s,
        gate,
    )?;
    let mut forced = centre_force_exchange(
        bodies,
        &staged,
        properties,
        gravity,
        centre_forces_world_n,
        duration_s,
        gate,
    )?;
    // Work and linear impulse are independent of event timing, but a body's
    // orbital force torque depends on its actual piecewise trajectory.
    forced.angular = force_angular;
    let energy_residual = after - before - supplied.work - forced.work + loss;
    let linear = linear_after
        .sub(linear_before)
        .sub(supplied.impulse)
        .sub(forced.impulse);
    let angular = angular_after
        .sub(angular_before)
        .sub(supplied.angular)
        .sub(forced.angular);
    if !energy_residual.is_finite()
        || !loss.is_finite()
        || !drift_defect.is_finite()
        || !linear.is_finite()
        || !angular.is_finite()
        || energy_residual.abs() > controls.total_energy_tolerance_j
        || super::magnitude(linear) > controls.linear_momentum_tolerance_ns
        || super::magnitude(angular) > controls.angular_momentum_tolerance_nms
    {
        return Err(SphereContactError::Audit(
            "whole trajectory energy or momentum",
        ));
    }
    poll(gate)?;
    bodies.copy_from_slice(&staged);
    Ok(SphereContactReport {
        duration_s,
        events,
        pair_queries: queries,
        minimum_path_gap_m: minimum_gap,
        kinetic_before_j: before,
        kinetic_after_j: after,
        gravity_work_j: supplied.work,
        gravity_impulse_ns: supplied.impulse,
        gravity_angular_impulse_nms: supplied.angular,
        centre_force_work_j: forced.work,
        centre_force_impulse_ns: forced.impulse,
        centre_force_angular_impulse_nms: forced.angular,
        unallocated_impact_loss_j: loss,
        drift_energy_defect_j: drift_defect,
        total_energy_residual_j: energy_residual,
        linear_momentum_residual_ns: linear,
        angular_momentum_residual_nms: angular,
    })
}

fn prepare(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    duration: f64,
    c: SphereContactControls,
) -> Result<Vec<Pair>, SphereContactError> {
    let n = bodies.len();
    let count = n
        .checked_mul(n.saturating_sub(1))
        .and_then(|v| v.checked_div(2))
        .ok_or(SphereContactError::Invalid("pair count overflow"))?;
    if n == 0
        || n != properties.len()
        || n != spheres.len()
        || n > c.impact.maximum_bodies
        || count != laws.len()
        || count > c.maximum_pairs
        || c.maximum_pairs > 65_536
        || c.maximum_events == 0
        || c.maximum_events > 4096
        || !duration.is_finite()
        || duration <= 0.0
        || [
            c.event_gap_tolerance_m,
            c.drift_energy_tolerance_j,
            c.total_energy_tolerance_j,
            c.linear_momentum_tolerance_ns,
            c.angular_momentum_tolerance_nms,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || c.event_gap_tolerance_m > c.impact.geometry_tolerance_m
    {
        return Err(SphereContactError::Invalid(
            "bounded complete spheres, laws, duration and gates required",
        ));
    }
    let mut order: Vec<_> = (0..n).collect();
    order.sort_by_key(|&i| spheres[i].id);
    for (k, &i) in order.iter().enumerate() {
        properties[i].validate()?;
        let s = spheres[i];
        let moments = properties[i].principal_inertia_body();
        let bound = 2.0 * properties[i].mass() * s.radius_m * s.radius_m;
        if !s.radius_m.is_finite()
            || s.radius_m <= 0.0
            || !bound.is_finite()
            || moments.x + moments.y + moments.z > bound
            || (k > 0 && spheres[order[k - 1]].id == s.id)
        {
            return Err(SphereContactError::Invalid(
                "distinct spheres and native inertia within radius required",
            ));
        }
    }
    let mut sorted = laws.to_vec();
    sorted.sort_by_key(|l| (l.first, l.second));
    let mut pairs = Vec::with_capacity(count);
    for (j, &a) in order.iter().enumerate() {
        for &b in &order[j + 1..] {
            let law = sorted[pairs.len()];
            let ids = [spheres[a].id, spheres[b].id];
            let radius = spheres[a].radius_m + spheres[b].radius_m;
            if [law.first, law.second] != ids
                || !law.friction_coefficient.is_finite()
                || law.friction_coefficient < 0.0
                || !law.restitution.is_finite()
                || !(0.0..=1.0).contains(&law.restitution)
                || !radius.is_finite()
            {
                return Err(SphereContactError::Invalid(
                    "one valid explicit law per canonical sphere pair required",
                ));
            }
            pairs.push(Pair {
                a,
                b,
                ids,
                radius,
                law,
            });
        }
    }
    Ok(pairs)
}

fn relative(bodies: &[RigidBodyState], p: Pair) -> Vec3 {
    bodies[p.a]
        .pose()
        .position_world()
        .sub(bodies[p.b].pose().position_world())
}
fn checked_gap(r: Vec3, pair: Pair) -> Result<f64, SphereContactError> {
    let gap = r.stable_norm("sphere separation")? - pair.radius;
    if !gap.is_finite() {
        return Err(SphereContactError::Unresolved("separation range"));
    }
    if gap < 0.0 {
        return Err(SphereContactError::Penetration {
            pair: pair.ids,
            gap_m: gap,
        });
    }
    Ok(gap)
}
fn gaps(
    bodies: &[RigidBodyState],
    pairs: &[Pair],
    gate: &CancelGate,
) -> Result<f64, SphereContactError> {
    let mut minimum = f64::INFINITY;
    for &pair in pairs {
        poll(gate)?;
        minimum = minimum.min(checked_gap(relative(bodies, pair), pair)?);
    }
    Ok(minimum)
}

fn first_hit(
    bodies: &[RigidBodyState],
    velocities: &[Vec3],
    accelerations: &[Vec3],
    pair: Pair,
    horizon: f64,
) -> Result<Option<f64>, SphereContactError> {
    let r = relative(bodies, pair);
    checked_gap(r, pair)?;
    let v = velocities[pair.a].sub(velocities[pair.b]);
    let motion = SweptSphereInput {
        relative_position_world: super::components(r),
        relative_velocity_world: super::components(v),
        radius_sum: pair.radius,
        duration_s: horizon,
    };
    let acceleration = accelerations[pair.a].sub(accelerations[pair.b]);
    if acceleration != Vec3::ZERO {
        return Ok(sphere_path(AcceleratedSphereInput {
            motion,
            relative_acceleration_world: super::components(acceleration),
        })?
        .first_entering_contact
        .map(|c| c.time_s));
    }
    let result = first_sphere_contact(motion)?;
    match result {
        // A genuine tangent touch never enters the interior and has no impulse.
        Some(hit) if hit.relative_normal_speed < 0.0 => Ok(Some(hit.time_s)),
        _ => Ok(None),
    }
}

fn path_gap(
    bodies: &[RigidBodyState],
    velocities: &[Vec3],
    accelerations: &[Vec3],
    pairs: &[Pair],
    h: f64,
    gate: &CancelGate,
) -> Result<f64, SphereContactError> {
    let mut minimum = f64::INFINITY;
    for &p in pairs {
        poll(gate)?;
        let r = relative(bodies, p);
        let acceleration = accelerations[p.a].sub(accelerations[p.b]);
        if acceleration != Vec3::ZERO {
            let path = sphere_path(AcceleratedSphereInput {
                motion: SweptSphereInput {
                    relative_position_world: super::components(r),
                    relative_velocity_world: super::components(
                        velocities[p.a].sub(velocities[p.b]),
                    ),
                    radius_sum: p.radius,
                    duration_s: h,
                },
                relative_acceleration_world: super::components(acceleration),
            })?;
            if path.minimum_gap < 0.0 {
                return Err(SphereContactError::Penetration {
                    pair: p.ids,
                    gap_m: path.minimum_gap,
                });
            }
            minimum = minimum.min(path.minimum_gap);
            continue;
        }
        let d = velocities[p.a].sub(velocities[p.b]).scale(h);
        let scale = r
            .stable_norm("path separation")?
            .max(d.stable_norm("path displacement")?);
        if scale == 0.0 {
            return Err(SphereContactError::Unresolved("zero path scale"));
        }
        let rs = r.scale(scale.recip());
        let ds = d.scale(scale.recip());
        let square = ds.norm_squared();
        let fraction = if square == 0.0 {
            0.0
        } else {
            (-rs.dot(ds) / square).clamp(0.0, 1.0)
        };
        minimum = minimum.min(checked_gap(r.add(d.scale(fraction)), p)?);
    }
    Ok(minimum)
}

fn contacts(
    bodies: &[RigidBodyState],
    spheres: &[ContactSphere],
    pairs: &[Pair],
    controls: SphereContactControls,
    gate: &CancelGate,
) -> Result<(Vec<RigidFrictionalContact>, Vec<[u64; 2]>, f64, bool), SphereContactError> {
    let mut result = Vec::new();
    let mut ids = Vec::new();
    let mut maximum_gap = 0.0_f64;
    let mut nonnegative = true;
    for (index, &p) in pairs.iter().enumerate() {
        poll(gate)?;
        let r = relative(bodies, p);
        let gap = checked_gap(r, p)?;
        if gap > controls.event_gap_tolerance_m {
            continue;
        }
        let n = r.normalized("sphere contact normal")?;
        let a = RigidImpactPoint {
            body: p.a,
            arm_body: bodies[p.a]
                .pose()
                .orientation()
                .rotate_world_to_body(n.scale(-spheres[p.a].radius_m)),
        };
        let b = RigidImpactPoint {
            body: p.b,
            arm_body: bodies[p.b]
                .pose()
                .orientation()
                .rotate_world_to_body(n.scale(spheres[p.b].radius_m)),
        };
        let pa = bodies[p.a].pose().point_world_from_body(a.arm_body)?;
        let pb = bodies[p.b].pose().point_world_from_body(b.arm_body)?;
        nonnegative &= pa.sub(pb).dot(n) >= 0.0;
        let axis = if n.x.abs() < 0.8 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        result.push(RigidFrictionalContact {
            id: index as u64,
            point: a,
            other: RigidImpactCounterpart::Body(b),
            normal_world: n,
            tangent_world: n.cross(axis),
            activation_tolerance_m: controls.event_gap_tolerance_m,
            friction_coefficient: p.law.friction_coefficient,
            normal_law: NormalVelocityLaw::NewtonRestitution {
                coefficient: p.law.restitution,
            },
        });
        ids.push(p.ids);
        maximum_gap = maximum_gap.max(gap);
    }
    Ok((result, ids, maximum_gap, nonnegative))
}

fn drift(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    gravity: Gravity,
    centre_forces: &[Vec3],
    h: f64,
    gate: &CancelGate,
) -> Result<Vec<RigidBodyState>, SphereContactError> {
    let integrator = RigidBodyIntegrator::new(gravity);
    bodies
        .iter()
        .zip(properties)
        .zip(centre_forces)
        .map(|((&b, &p), &force)| {
            poll(gate)?;
            Ok(integrator
                .step(
                    b,
                    p,
                    Wrench {
                        force_world: force,
                        torque_body: Vec3::ZERO,
                    },
                    h,
                )?
                .state_after)
        })
        .collect()
}
fn drift_audit(
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    properties: &[MassProperties],
    gravity: Gravity,
    centre_forces: &[Vec3],
    duration_s: f64,
    c: SphereContactControls,
    gate: &CancelGate,
) -> Result<f64, SphereContactError> {
    let a = super::totals(before, properties, gate)?;
    let b = super::totals(after, properties, gate)?;
    let supplied = gravity_exchange(
        before,
        after,
        properties,
        gravity,
        centre_forces,
        duration_s,
        gate,
    )?;
    let forced = centre_force_exchange(
        before,
        after,
        properties,
        gravity,
        centre_forces,
        duration_s,
        gate,
    )?;
    let energy = b.0 - a.0 - supplied.work - forced.work;
    if !energy.is_finite()
        || energy.abs() > c.drift_energy_tolerance_j
        || super::magnitude(b.1.sub(a.1).sub(supplied.impulse).sub(forced.impulse))
            > c.linear_momentum_tolerance_ns
        || super::magnitude(b.2.sub(a.2).sub(supplied.angular).sub(forced.angular))
            > c.angular_momentum_tolerance_nms
    {
        return Err(SphereContactError::Audit("native constant-force drift"));
    }
    Ok(energy)
}

struct GravityExchange {
    work: f64,
    impulse: Vec3,
    angular: Vec3,
}

fn gravity_exchange(
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    properties: &[MassProperties],
    gravity: Gravity,
    centre_forces: &[Vec3],
    duration_s: f64,
    gate: &CancelGate,
) -> Result<GravityExchange, SphereContactError> {
    let mut exchange = GravityExchange {
        work: 0.0,
        impulse: Vec3::ZERO,
        angular: Vec3::ZERO,
    };
    let g = gravity.acceleration_world();
    if g == Vec3::ZERO {
        return Ok(exchange);
    }
    for (((start, end), property), &centre_force) in
        before.iter().zip(after).zip(properties).zip(centre_forces)
    {
        poll(gate)?;
        let force = g.scale(property.mass());
        let impulse = force.scale(duration_s);
        exchange.work += force.dot(
            end.pose()
                .position_world()
                .sub(start.pose().position_world()),
        );
        exchange.impulse = exchange.impulse.add(impulse);
        // Across internal impacts, sum(m*x) still follows the closed system's
        // uniformly accelerated COM equation. The quadratic g term has zero
        // cross product with g. Thus the TOTAL torque integral depends only
        // on the original positions and total momentum, not contact allocation.
        exchange.angular = exchange
            .angular
            .add(start.pose().position_world().cross(impulse))
            .add(
                start
                    .linear_momentum_world()
                    .cross(g.scale(duration_s))
                    .scale(0.5 * duration_s),
            )
            // The original closed-system mass moment also accelerates under
            // the constant sum of centre forces. Internal impulses cancel.
            .add(centre_force.cross(g).scale((duration_s / 6.0) * duration_s * duration_s));
    }
    if !exchange.work.is_finite() || !exchange.impulse.is_finite() || !exchange.angular.is_finite()
    {
        return Err(SphereContactError::Audit(
            "gravitational work or impulse overflow",
        ));
    }
    Ok(exchange)
}

#[allow(clippy::too_many_arguments)]
fn centre_force_exchange(
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    properties: &[MassProperties],
    gravity: Gravity,
    forces: &[Vec3],
    h: f64,
    gate: &CancelGate,
) -> Result<GravityExchange, SphereContactError> {
    let mut exchange = GravityExchange {
        work: 0.0,
        impulse: Vec3::ZERO,
        angular: Vec3::ZERO,
    };
    for (((start, end), &p), &force) in before.iter().zip(after).zip(properties).zip(forces) {
        poll(gate)?;
        if force == Vec3::ZERO {
            continue;
        }
        let impulse = force.scale(h);
        exchange.work += force.dot(
            end.pose()
                .position_world()
                .sub(start.pose().position_world()),
        );
        exchange.impulse = exchange.impulse.add(impulse);
        // Integral x(t)×F dt = x0×F h + v0×F h²/2 + g×F h³/6.
        // The body's own centre-force acceleration is parallel to F.
        exchange.angular = exchange
            .angular
            .add(start.pose().position_world().cross(impulse))
            .add(
                start
                    .center_of_mass_velocity_world(p)?
                    .cross(impulse)
                    .scale(0.5 * h),
            )
            .add(
                gravity
                    .acceleration_world()
                    .cross(impulse)
                    .scale((h / 6.0) * h),
            );
    }
    if !exchange.work.is_finite() || !exchange.impulse.is_finite() || !exchange.angular.is_finite()
    {
        return Err(SphereContactError::Audit(
            "centre force work or impulse overflow",
        ));
    }
    Ok(exchange)
}
fn poll(gate: &CancelGate) -> Result<(), SphereContactError> {
    if gate.is_requested() {
        Err(SphereContactError::Cancelled)
    } else {
        Ok(())
    }
}
