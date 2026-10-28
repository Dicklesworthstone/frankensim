//! Automatic force-free trajectories of finite native spheres.
//!
//! Every pair is searched over its complete linear centre path. The existing
//! scalar ballistic owner isolates zeros of the squared-distance polynomial;
//! all touching pairs then enter one simultaneous normal/Coulomb event. Native
//! free rotation and point impulses keep their original inertia and audits.
//! This floating analytic lane is not interval-certified general-shape CCD.

use super::{
    RigidFrictionalContact, RigidFrictionalImpactError, RigidFrictionalImpactReport,
    RigidImpactControls, RigidImpactCounterpart, RigidImpactPoint, apply_rigid_frictional_impacts,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};
use fs_contact::ballistic::{
    BallisticContactError, BallisticContactInput, first_ballistic_contact,
};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
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
    pub unallocated_impact_loss_j: f64,
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
    Detection(BallisticContactError),
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
impl From<BallisticContactError> for SphereContactError {
    fn from(e: BallisticContactError) -> Self {
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
/// This is force-free impulsive dynamics, with no sustained external load,
/// rolling/spin resistance, finite-patch history, shape deformation, or thermal
/// allocation. The original contact solver's admissibility/determinacy limits
/// apply. Native free attitude is approximate for general asymmetric bodies;
/// its original energy/momentum checks and time refinement remain necessary.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_sphere_contacts(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    duration_s: f64,
    controls: SphereContactControls,
    gate: &CancelGate,
) -> Result<SphereContactReport, SphereContactError> {
    poll(gate)?;
    let pairs = prepare(bodies, properties, spheres, laws, duration_s, controls)?;
    let (before, linear_before, angular_before) = super::totals(bodies, properties, gate)?;
    let mut staged = bodies.to_vec();
    let mut elapsed = 0.0;
    let mut events = Vec::new();
    let mut queries = 0usize;
    let mut minimum_gap = gaps(&staged, &pairs, gate)?;
    let mut loss = 0.0;
    let mut drift_defect = 0.0;
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
            if let Some(time) = first_hit(&staged, &velocities, *pair, remaining)? {
                hit = Some(hit.map_or(time, |old| old.min(time)));
            }
        }
        let Some(mut h) = hit else {
            if remaining > 0.0 {
                minimum_gap =
                    minimum_gap.min(path_gap(&staged, &velocities, &pairs, remaining, gate)?);
                let next = drift(&staged, properties, remaining, gate)?;
                minimum_gap = minimum_gap.min(gaps(&next, &pairs, gate)?);
                drift_defect += drift_audit(&staged, &next, properties, controls, gate)?;
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
                drift(&staged, properties, h, gate)?
            };
            if let Ok(gap) = gaps(&next, &pairs, gate) {
                let (contacts, ids, maximum_gap, witnesses_nonnegative) =
                    contacts(&next, spheres, &pairs, controls, gate)?;
                if !contacts.is_empty() && witnesses_nonnegative {
                    candidate = Some((next, contacts, ids, maximum_gap, gap));
                    break;
                }
            }
            if h == 0.0 {
                break;
            }
            h = h.next_down();
        }
        let Some((mut next, contacts, ids, maximum_gap, gap)) = candidate else {
            return Err(SphereContactError::Unresolved(
                "event geometry or outward root bracket",
            ));
        };
        minimum_gap = minimum_gap
            .min(path_gap(&staged, &velocities, &pairs, h, gate)?)
            .min(gap);
        drift_defect += drift_audit(&staged, &next, properties, controls, gate)?;
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
    let energy_residual = after - before + loss;
    let linear = linear_after.sub(linear_before);
    let angular = angular_after.sub(angular_before);
    if !energy_residual.is_finite()
        || !loss.is_finite()
        || !drift_defect.is_finite()
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
    pair: Pair,
    horizon: f64,
) -> Result<Option<f64>, SphereContactError> {
    let r = relative(bodies, pair);
    let gap = checked_gap(r, pair)?;
    let v = velocities[pair.a].sub(velocities[pair.b]);
    let distance = r.stable_norm("sphere separation")?;
    let speed = v.stable_norm("relative sphere velocity")?;
    if horizon == 0.0 || speed == 0.0 {
        return Ok(if gap == 0.0 && r.dot(v) < 0.0 {
            Some(0.0)
        } else {
            None
        });
    }
    // Use dimensionless time t/horizon and common length scaling. Form the
    // constant as (d-R)(d+R) to retain small resolved initial separations.
    let scale = distance.max(pair.radius).max(speed * horizon);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(SphereContactError::Unresolved("sweep scaling"));
    }
    let rs = r.scale(scale.recip());
    let vs = v.scale(horizon / scale);
    let constant = (gap / scale) * (distance / scale + pair.radius / scale);
    if gap != 0.0 && constant == 0.0 {
        return Err(SphereContactError::Unresolved("sweep loses nonzero gap"));
    }
    let result = first_ballistic_contact(BallisticContactInput {
        gap: constant,
        relative_normal_speed: 2.0 * rs.dot(vs),
        relative_normal_acceleration: 2.0 * vs.norm_squared(),
        duration_s: 1.0,
    })?;
    match result {
        // A genuine tangent touch never enters the interior and has no impulse.
        Some(hit) if hit.relative_normal_speed < 0.0 => {
            let time = hit.time_s * horizon;
            if !time.is_finite() || (hit.time_s > 0.0 && time == 0.0) {
                return Err(SphereContactError::Unresolved("event time range"));
            }
            Ok(Some(time))
        }
        _ => Ok(None),
    }
}

fn path_gap(
    bodies: &[RigidBodyState],
    velocities: &[Vec3],
    pairs: &[Pair],
    h: f64,
    gate: &CancelGate,
) -> Result<f64, SphereContactError> {
    let mut minimum = f64::INFINITY;
    for &p in pairs {
        poll(gate)?;
        let r = relative(bodies, p);
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
    h: f64,
    gate: &CancelGate,
) -> Result<Vec<RigidBodyState>, SphereContactError> {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    bodies
        .iter()
        .zip(properties)
        .map(|(&b, &p)| {
            poll(gate)?;
            Ok(integrator.step(b, p, Wrench::ZERO, h)?.state_after)
        })
        .collect()
}
fn drift_audit(
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    properties: &[MassProperties],
    c: SphereContactControls,
    gate: &CancelGate,
) -> Result<f64, SphereContactError> {
    let a = super::totals(before, properties, gate)?;
    let b = super::totals(after, properties, gate)?;
    let energy = b.0 - a.0;
    if !energy.is_finite()
        || energy.abs() > c.drift_energy_tolerance_j
        || super::magnitude(b.1.sub(a.1)) > c.linear_momentum_tolerance_ns
        || super::magnitude(b.2.sub(a.2)) > c.angular_momentum_tolerance_nms
    {
        return Err(SphereContactError::Audit("native free drift"));
    }
    Ok(energy)
}
fn poll(gate: &CancelGate) -> Result<(), SphereContactError> {
    if gate.is_requested() {
        Err(SphereContactError::Cancelled)
    } else {
        Ok(())
    }
}
