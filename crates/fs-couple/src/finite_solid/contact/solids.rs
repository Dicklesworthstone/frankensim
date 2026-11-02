//! Simultaneous contact among original finite-deformation solids.
//!
//! Every material-point trace uses its body's supported consistent mass. Signed
//! incidence assembles every cross-contact block before ONE original Coulomb
//! solve; one summed impulse then updates each native body. Distinct points on
//! the same body retain both signed traces and their cross terms. No lumped
//! point masses, pairwise restitution sweep or configuration repair is used.
//!
//! Witnesses are caller-supplied current material points. This is an atomic
//! fixed-configuration event, not collision detection or a trajectory step.
//! Actual native endpoint velocities, quadratic work and total linear/angular
//! momentum are checked, retaining stationary-support reactions. Finite-gap
//! friction couples cannot be hidden by relocating the supplied witnesses.

use super::{array, contact_frame, magnitude, vector, world};
use crate::finite_solid::impulse::{
    FiniteSolidImpulse, FiniteSolidImpulseControls, FiniteSolidMomentum, FiniteSolidPoint,
};
use crate::finite_solid::{FiniteSolidBody, FiniteSolidError};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResult,
    NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::Cx;
use fs_mbd::Vec3;
use fs_solver::SolveReport;

/// One original material witness in the input body array.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidContactPoint {
    pub body: usize,
    pub point: FiniteSolidPoint,
}

/// A pair of current material witnesses on the same contact normal line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidPairContact {
    pub id: u64,
    pub first: SolidContactPoint,
    pub second: SolidContactPoint,
    /// From the second witness toward the first; positive velocity separates.
    pub normal_world: Vec3,
    pub tangent_world: Vec3,
    pub activation_tolerance_m: f64,
    /// Explicit interface coefficient; no bulk-to-friction inference is made.
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// Bounded assembly and native contact/mass-solve admission controls.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidContactAssemblyControls {
    /// At most eight original prepared bodies in this bounded event owner.
    pub maximum_solids: usize,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
    pub solid_impulse: FiniteSolidImpulseControls,
    pub contact: FrictionalContactConfig,
}
impl Default for FiniteSolidContactAssemblyControls {
    fn default() -> Self {
        Self {
            maximum_solids: 8,
            geometry_tolerance_m: 1e-10,
            velocity_tolerance_m_s: 1e-8,
            linear_tolerance_kg_m_s: 1e-9,
            angular_tolerance_kg_m2_s: 1e-9,
            energy_tolerance_j: 1e-8,
            solid_impulse: FiniteSolidImpulseControls::default(),
            contact: FrictionalContactConfig::default(),
        }
    }
}

/// Contact arrays use stable ID order; body arrays retain input order.
#[derive(Clone, Debug)]
pub struct FiniteSolidContactAssemblyReport {
    pub contact_ids: Vec<u64>,
    /// Complete signed mobility, row-major 3N by 3N [1/kg].
    pub contact_mobility_per_kg: Vec<f64>,
    pub solid_mass_solves: Vec<Vec<SolveReport>>,
    pub contact: FrictionalContactResult,
    pub solid_impulses: Vec<FiniteSolidImpulse>,
    /// External stationary-support impulse about the caller's common origin.
    pub support_impulse: FiniteSolidMomentum,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Sum of actual native quadratic increments, not total-energy subtraction.
    pub kinetic_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub dissipation_residual_j: f64,
    /// Actual change minus external stationary-support reaction.
    pub linear_momentum_residual_kg_m_s: Vec3,
    pub angular_momentum_residual_kg_m2_s: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
    /// Opposing native-body impulse errors cannot cancel these absolute gates.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    /// Sum of each supplied pair's |(x_first - x_second) cross impulse|.
    /// Opposing finite-gap friction couples cannot cancel this angular gate.
    pub absolute_contact_couple_kg_m2_s: f64,
    pub absolute_work_residual_j: f64,
}

#[derive(Debug)]
pub enum FiniteSolidContactAssemblyError {
    Invalid(&'static str),
    Solid {
        body: usize,
        error: FiniteSolidError,
    },
    Contact(FrictionalContactError),
    Cancelled,
    /// No body candidate is returned after a failed physical audit.
    Audit(Box<FiniteSolidContactAssemblyReport>),
}
impl core::fmt::Display for FiniteSolidContactAssemblyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-solid assembly contact: {self:?}")
    }
}
impl std::error::Error for FiniteSolidContactAssemblyError {}

/// Resolve all supplied deformable/deformable contacts together.
///
/// Bodies may have different source materials, reference densities and fixed
/// supports, but must share one physical time. Sources, configuration, native
/// acceleration and accepted-step history are unchanged by this instantaneous
/// event. Distinct witnesses on one body are allowed without dropping their
/// signed cross-point response. The original contact solver retains its range,
/// nonuniqueness, restitution and passivity refusals; representative redundant
/// impulses are not promoted to a unique per-contact heat allocation.
#[allow(clippy::too_many_lines)]
pub fn contact_candidate<'a>(
    bodies: &[FiniteSolidBody<'a>],
    contacts: &[FiniteSolidPairContact],
    origin_m: [f64; 3],
    controls: FiniteSolidContactAssemblyControls,
    cx: &Cx<'_>,
) -> Result<
    (Vec<FiniteSolidBody<'a>>, FiniteSolidContactAssemblyReport),
    FiniteSolidContactAssemblyError,
> {
    use FiniteSolidContactAssemblyError as E;
    poll(cx)?;
    if bodies.is_empty()
        || controls.maximum_solids == 0
        || controls.maximum_solids > 8
        || bodies.len() > controls.maximum_solids
        || controls.contact.maximum_contacts == 0
        || controls.contact.maximum_contacts > 21
        || contacts.len() > controls.contact.maximum_contacts
        || !vector(origin_m).is_finite()
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.linear_tolerance_kg_m_s,
            controls.angular_tolerance_kg_m2_s,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(E::Invalid(
            "bounded bodies/contacts and positive finite physical controls required",
        ));
    }
    let time = bodies[0].native_state().t;
    if !time.is_finite() || bodies.iter().any(|b| b.native_state().t != time) {
        return Err(E::Invalid(
            "contact bodies must share one finite physical time",
        ));
    }
    // Bounds above cap global rows at 63 and total signed point rows at 126.
    // Native mass solves separately enforce their original workspace budgets.
    let dimension = 3 * contacts.len();
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_by_key(|&i| contacts[i].id);
    let mut rows = vec![Vec::new(); bodies.len()];
    let mut global_rows = vec![Vec::new(); bodies.len()];
    let mut axes = Vec::with_capacity(contacts.len());
    let mut separations = Vec::with_capacity(contacts.len());
    let mut laws = Vec::with_capacity(contacts.len());
    for (slot, &index) in order.iter().enumerate() {
        poll(cx)?;
        let c = contacts[index];
        if c.first.body >= bodies.len()
            || c.second.body >= bodies.len()
            || c.first == c.second
            || (slot > 0 && contacts[order[slot - 1]].id == c.id)
        {
            return Err(E::Invalid(
                "distinct valid material witnesses and unique contact IDs required",
            ));
        }
        let frame = contact_frame(c.normal_world, c.tangent_world).map_err(E::Invalid)?;
        let mut positions = [Vec3::ZERO; 2];
        let mut velocities = [Vec3::ZERO; 2];
        for (endpoint, point) in [c.first, c.second].into_iter().enumerate() {
            let body = &bodies[point.body];
            let source_error = |error| E::Solid {
                body: point.body,
                error,
            };
            positions[endpoint] = vector(
                body.point_position_m(point.point, cx)
                    .map_err(source_error)?,
            );
            velocities[endpoint] = vector(
                body.point_velocity_m_s(point.point, cx)
                    .map_err(source_error)?,
            );
            let sign = if endpoint == 0 { 1.0 } else { -1.0 };
            for (axis, direction) in frame.into_iter().enumerate() {
                rows[point.body].push((point.point, array(direction.scale(sign))));
                global_rows[point.body].push(3 * slot + axis);
            }
        }
        let separation = positions[0].sub(positions[1]);
        let gap = frame[0].dot(separation);
        if magnitude(separation.sub(frame[0].scale(gap))) > controls.geometry_tolerance_m {
            return Err(E::Invalid(
                "material witnesses do not share a contact normal line",
            ));
        }
        let relative = velocities[0].sub(velocities[1]);
        laws.push(FrictionalContact {
            id: c.id,
            gap,
            activation_tolerance: c.activation_tolerance_m,
            relative_velocity: frame.map(|axis| relative.dot(axis)),
            friction_coefficient: c.friction_coefficient,
            normal_law: c.normal_law,
        });
        axes.push(frame);
        separations.push(separation);
    }
    let mut matrix = vec![0.0; dimension * dimension];
    let mut solid_mass_solves = Vec::with_capacity(bodies.len());
    for (index, body) in bodies.iter().enumerate() {
        poll(cx)?;
        let mobility = body
            .point_mobility_matrix_per_kg(
                &rows[index],
                6 * controls.contact.maximum_contacts,
                controls.solid_impulse.mass_solve,
                cx,
            )
            .map_err(|error| E::Solid { body: index, error })?;
        scatter_mobility(
            &mut matrix,
            dimension,
            &global_rows[index],
            &mobility.matrix_per_kg,
        )
        .map_err(E::Invalid)?;
        solid_mass_solves.push(mobility.solves);
    }
    let contact = solve_frictional_contacts(&laws, &matrix, controls.contact, || {
        cx.checkpoint().is_err()
    })
    .map_err(E::Contact)?;
    let mut impulses = vec![Vec::new(); bodies.len()];
    let mut absolute_contact_couple = 0.0;
    for (slot, &index) in order.iter().enumerate() {
        poll(cx)?;
        let c = contacts[index];
        let impulse = world(axes[slot], contact.impulses[slot]);
        absolute_contact_couple += magnitude(separations[slot].cross(impulse));
        impulses[c.first.body].push((c.first.point, array(impulse)));
        impulses[c.second.body].push((c.second.point, array(impulse.scale(-1.0))));
    }
    let mut next = Vec::with_capacity(bodies.len());
    let mut solid_impulses = Vec::with_capacity(bodies.len());
    let mut before = 0.0;
    let mut after = 0.0;
    let mut change = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut support_linear = Vec3::ZERO;
    let mut support_angular = Vec3::ZERO;
    let mut absolute_linear = 0.0;
    let mut absolute_angular = 0.0;
    let mut absolute_work = 0.0;
    for (index, body) in bodies.iter().enumerate() {
        poll(cx)?;
        // One summed impulse preserves all within-body cross-contact work.
        let (candidate, receipt) = body
            .point_impulses_candidate(&impulses[index], origin_m, controls.solid_impulse, cx)
            .map_err(|error| E::Solid { body: index, error })?;
        before += receipt.kinetic_before_j;
        after += receipt.kinetic_after_j;
        change += receipt.kinetic_change_j;
        // This is actual delta momentum minus native support reaction, without
        // subtracting two large accepted momenta under a common boost.
        linear = linear
            .add(vector(receipt.applied.linear_kg_m_s))
            .add(vector(receipt.linear_residual_kg_m_s));
        angular = angular
            .add(vector(receipt.applied.angular_kg_m2_s))
            .add(vector(receipt.angular_residual_kg_m2_s));
        support_linear = support_linear.add(vector(receipt.support.linear_kg_m_s));
        support_angular = support_angular.add(vector(receipt.support.angular_kg_m2_s));
        absolute_linear += magnitude(vector(receipt.linear_residual_kg_m_s));
        absolute_angular += magnitude(vector(receipt.angular_residual_kg_m2_s));
        absolute_work += receipt.work_residual_j.abs();
        next.push(candidate);
        solid_impulses.push(receipt);
    }
    let mut relative_work = 0.0;
    let mut endpoint_error = 0.0_f64;
    for (slot, &index) in order.iter().enumerate() {
        poll(cx)?;
        let c = contacts[index];
        let first = next[c.first.body]
            .point_velocity_m_s(c.first.point, cx)
            .map_err(|error| E::Solid {
                body: c.first.body,
                error,
            })?;
        let second = next[c.second.body]
            .point_velocity_m_s(c.second.point, cx)
            .map_err(|error| E::Solid {
                body: c.second.body,
                error,
            })?;
        let relative = vector(first).sub(vector(second));
        for axis in 0..3 {
            let velocity = relative.dot(axes[slot][axis]);
            endpoint_error = endpoint_error
                .max((velocity - contact.relative_velocities_after[slot][axis]).abs());
            relative_work += (0.5 * laws[slot].relative_velocity[axis] + 0.5 * velocity)
                * contact.impulses[slot][axis];
        }
    }
    let solver_work = contact.relative_work;
    let report = FiniteSolidContactAssemblyReport {
        contact_ids: laws.iter().map(|c| c.id).collect(),
        contact_mobility_per_kg: matrix,
        solid_mass_solves,
        solid_impulses,
        support_impulse: FiniteSolidMomentum {
            origin_m,
            linear_kg_m_s: array(support_linear),
            angular_kg_m2_s: array(support_angular),
        },
        kinetic_before_j: before,
        kinetic_after_j: after,
        kinetic_change_j: change,
        relative_impulse_work_j: relative_work,
        work_energy_residual_j: change - relative_work,
        dissipation_residual_j: change + contact.dissipated_energy,
        linear_momentum_residual_kg_m_s: linear,
        angular_momentum_residual_kg_m2_s: angular,
        endpoint_velocity_residual_m_s: endpoint_error,
        absolute_linear_residual_kg_m_s: absolute_linear,
        absolute_angular_residual_kg_m2_s: absolute_angular,
        absolute_contact_couple_kg_m2_s: absolute_contact_couple,
        absolute_work_residual_j: absolute_work,
        contact,
    };
    if [
        before,
        after,
        change,
        relative_work,
        endpoint_error,
        report.work_energy_residual_j,
        report.dissipation_residual_j,
        magnitude(linear),
        magnitude(angular),
        magnitude(support_linear),
        magnitude(support_angular),
        absolute_linear,
        absolute_angular,
        absolute_contact_couple,
        absolute_work,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || endpoint_error > controls.velocity_tolerance_m_s
        || magnitude(linear) > controls.linear_tolerance_kg_m_s
        || magnitude(angular) > controls.angular_tolerance_kg_m2_s
        || absolute_linear > controls.linear_tolerance_kg_m_s
        || absolute_angular > controls.angular_tolerance_kg_m2_s
        || absolute_contact_couple > controls.angular_tolerance_kg_m2_s
        || absolute_work > controls.energy_tolerance_j
        || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
        || report.dissipation_residual_j.abs() > controls.energy_tolerance_j
        || (relative_work - solver_work).abs() > controls.energy_tolerance_j
        || change > controls.energy_tolerance_j
    {
        return Err(E::Audit(Box::new(report)));
    }
    poll(cx)?;
    Ok((next, report))
}

fn poll(cx: &Cx<'_>) -> Result<(), FiniteSolidContactAssemblyError> {
    cx.checkpoint()
        .map_err(|_| FiniteSolidContactAssemblyError::Cancelled)
}

// Signed directions are already in the original point-mobility queries. A
// same-body pair repeats each global row; BOTH cross terms must accumulate.
// Preserve every native directed entry: no symmetry or diagonal repair.
fn scatter_mobility(
    global: &mut [f64],
    dimension: usize,
    rows: &[usize],
    local: &[f64],
) -> Result<(), &'static str> {
    if global.len() != dimension * dimension
        || local.len() != rows.len() * rows.len()
        || rows.iter().any(|&row| row >= dimension)
    {
        return Err("native point mobility dimensions do not match signed incidence");
    }
    for (i, &row) in rows.iter().enumerate() {
        for (j, &column) in rows.iter().enumerate() {
            let entry = &mut global[row * dimension + column];
            *entry += local[i * rows.len() + j];
            if !entry.is_finite() {
                return Err("assembled contact mobility is not representable");
            }
        }
    }
    Ok(())
}
