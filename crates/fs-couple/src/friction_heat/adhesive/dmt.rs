//! Driven native DMT normal/shear motion and original finite-body heat.
//!
//! The original pair solves its complete normal/shear mobility once. Only its
//! actual junction-friction loss enters the existing massless-face transport
//! solve; reversible normal adhesion/elastic storage is not caloric heat.
//! Native momenta, original contact history and all finite-cell temperatures
//! are admitted together, including actual work from prescribed body wrenches.
//!
//! Mechanical sources and cell geometry are frozen over this increment. The
//! native affine contact-discretization path supplies its exact mean real area
//! to both contact conduction and ambient surface coverage. The true endpoint
//! area is installed with the
//! caloric endpoint. This remains backward-Euler transport with endpoint
//! temperatures/conductivities, not an implicit thermomechanical accuracy
//! guarantee. The returned endpoint route continues on these original cells.
//! The pose/topology owner still consumes native travel and rebinds moved face
//! centres. Original source temperature continuation belongs at that boundary.

use super::super::FrictionHeatRoute;
use super::super::cohesive::{CohesiveHeatError, solve_frozen_cohesive_heat};
use super::super::network::{FrictionHeatNetworkReport, FrictionHeatNetworkStep};
use super::area::{ContactAreaEvolution, ContactAreaPath};
use fs_contact::normal_patch::dmt::coupled::{DmtCoupledError, DmtCoupledState};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_math::{dd::Dd, det};
use fs_mbd::tangential_friction::adhesive::normal_motion::{
    PreparedRigidDmtPair, RigidDmtContact, RigidDmtControls, RigidDmtError, RigidDmtReport,
};
use fs_mbd::{MassProperties, RigidBodyState, Vec3, Wrench};
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, MaterialId};

pub use super::super::cohesive::{
    CohesiveHeatConductivity as RigidDmtHeatConductivity,
    CohesiveHeatControls as RigidDmtThermalControls, CohesiveHeatSource as RigidDmtHeatSource,
};

/// Explicit original native-surface to thermal-material assignment, in A/B order.
#[derive(Clone, Copy)]
pub struct RigidDmtHeatBody<'a> {
    pub body_id: u64,
    pub thermal_material: MaterialId,
    pub mechanical_material: &'a MaterialStateId,
}

/// One actual DMT pair and its coincident incoming solid-contact thermal face.
#[derive(Clone, Copy)]
pub struct RigidDmtHeatContact<'a> {
    pub contact: RigidDmtContact<'a>,
    pub bodies: [RigidDmtHeatBody<'a>; 2],
    pub route: &'a FrictionHeatRoute,
    /// Positive numerical seed only; the face temperature is solved.
    pub temperature_guess_k: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidDmtHeatControls {
    pub mechanical: RigidDmtControls,
    pub thermal: RigidDmtThermalControls,
    /// Maximum |A(t)-A(0)|/A(0) over the original affine normal path.
    /// Finite and nonnegative; never silently clips the actual mechanical area.
    pub maximum_fractional_area_change: f64,
}

/// Original true solid areas, excluding the noncontact attractive annulus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDmtHeatAreas {
    /// Original physical footprint before this increment.
    pub incoming_m2: f64,
    pub endpoint_m2: f64,
    pub path_mean_m2: f64,
    /// Exact extrema for the admitted monotone affine normal coordinate.
    pub path_extrema_m2: [f64; 2],
    pub maximum_fractional_change: f64,
}

#[derive(Debug)]
pub enum RigidDmtHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidDmtError),
    Heat(CohesiveHeatError<E>),
    AreaChange { observed: f64, maximum: f64 },
}
impl<E> From<RigidDmtError> for RigidDmtHeatError<E> {
    fn from(error: RigidDmtError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<CohesiveHeatError<E>> for RigidDmtHeatError<E> {
    fn from(error: CohesiveHeatError<E>) -> Self {
        Self::Heat(error)
    }
}
impl<E> From<DepositionError> for RigidDmtHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Heat(error.into())
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidDmtHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native DMT heat refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidDmtHeatError<E> {}

pub struct HeatedRigidDmtReport<P> {
    /// Original native normal storage, applied work, friction loss and travel.
    pub mechanical: RigidDmtReport,
    pub heat: FrictionHeatNetworkReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    pub target_time_s: f64,
    pub areas: RigidDmtHeatAreas,
    /// Actual endpoint footprint on the unchanged original finite cells.
    /// Retains reference conductance per area; the next callback resolves k.
    pub route_after: FrictionHeatRoute,
    /// Delta K + Delta U_normal + Delta H - external body work - other heat.
    pub combined_energy_residual_j: f64,
    /// Sum of the independent native mechanical and transport residuals.
    pub summed_energy_residual_j: f64,
}

/// Complete staged native and finite-cell endpoint; no accepted owner changed.
pub struct HeatedRigidDmtCandidate<'t, P> {
    /// Original contact A/B order, independent of native array indexing.
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: DmtCoupledState,
    pub thermal: DepositionSimulation<'t>,
    pub report: HeatedRigidDmtReport<P>,
}

/// Solve one original driven pair and its genuine finite thermal receivers.
///
/// The full straight normal path must have positive real area. Touch, gap,
/// detachment and thermal-face/source/body mismatches refuse before publication.
/// Original k*A/d and ambient exposed area use the exact time mean of the
/// native affine contact-discretization path. Accepted coverage and
/// `route_after` retain its true endpoint area, original cells, material and
/// mass. Endpoint coverage
/// must fit every original face and the complete shared-cell surface budget.
/// Scheduled topology events at that endpoint refuse a changed patch until
/// their original owner processes them. There is no heat fraction, added
/// capacity or invented traction distribution.
/// Unequal body temperatures can make one body's signed heat negative.
///
/// `resolve` supplies endpoint bulk conductivities and OTHER thermal sources,
/// retaining staged source receipts in its payload. It must not publish history
/// or reinsert friction loss. Thermal trials never repeat mechanical evolution.
/// The callback sees the mean-area thermal geometry used by this residual.
/// Source coefficients remain frozen even when the solved temperature changes;
/// use the native DMT temperature refresh at the next accepted geometry boundary
/// and refine both time and thermal space for feedback/transport accuracy.
///
/// Poses remain held as in `PreparedRigidDmtPair`. The returned actual body,
/// marker and reaction-point travel belongs to the surrounding pose owner.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_dmt_heat<'t, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDmtHeatContact<'_>,
    history: DmtCoupledState,
    wrenches: [Wrench; 2],
    thermal: &DepositionSimulation<'t>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidDmtHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<RigidDmtHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidDmtCandidate<'t, P>, RigidDmtHeatError<E>> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    let target = accepted_time_s + dt_s;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= accepted_time_s
        || !controls.maximum_fractional_area_change.is_finite()
        || controls.maximum_fractional_area_change < 0.0
        || !contact.temperature_guess_k.is_finite()
        || contact.temperature_guess_k <= 0.0
    {
        return Err(RigidDmtHeatError::Invalid(
            "synchronized positive duration, area bound and thermal seed required",
        ));
    }
    contact.route.validate(thermal)?;
    // Use the thermal owner's represented interval for the actual native step.
    let mut prepared = PreparedRigidDmtPair::new(
        bodies,
        properties,
        contact.contact,
        target - accepted_time_s,
        controls.mechanical,
        gate,
    )?;
    let before = contact
        .contact
        .normal
        .evaluate(history.position_m()[0])
        .map_err(|e| RigidDmtError::from(DmtCoupledError::from(e)))?;
    if before.indentation_m <= 0.0 || before.contact_area_m2 <= 0.0 {
        return Err(RigidDmtHeatError::Invalid(
            "a finite solid DMT junction is required for this conducting face",
        ));
    }
    validate_binding(
        bodies,
        properties,
        &prepared,
        contact,
        thermal,
        before.contact_area_m2,
        controls.mechanical.geometry_tolerance_m,
    )?;
    let mechanical = prepared.advance_with_wrenches(history, wrenches, gate)?;
    let path = mechanical.report.contact.path;
    let a0 = before.contact_area_m2;
    let a1 = path.contact.contact_area_m2;
    if path.contact.indentation_m <= 0.0 || a1 <= 0.0 {
        return Err(RigidDmtHeatError::Invalid(
            "touch, opening and gap transport need their original topology owner",
        ));
    }
    // On attached DMT contact A = pi*R*delta. Its affine normal path is
    // monotone, so these endpoints are the true extrema of the entire path.
    let fractional = (a1 - a0).abs() / a0;
    if !fractional.is_finite() || fractional > controls.maximum_fractional_area_change {
        return Err(RigidDmtHeatError::AreaChange {
            observed: fractional,
            maximum: controls.maximum_fractional_area_change,
        });
    }
    let areas = RigidDmtHeatAreas {
        incoming_m2: a0,
        endpoint_m2: a1,
        path_mean_m2: path.mean_real_area_m2,
        path_extrema_m2: [a0.min(a1), a0.max(a1)],
        maximum_fractional_change: fractional,
    };
    let loss = mechanical.report.contact.friction_dissipation_j;
    if !loss.is_finite() || loss < 0.0 || loss != mechanical.report.dissipated_energy_j {
        return Err(RigidDmtHeatError::Invalid(
            "native source and actual body friction loss differ",
        ));
    }
    let area_path = ContactAreaEvolution::new(
        thermal,
        &[ContactAreaPath {
            route: contact.route,
            mean_area_m2: path.mean_real_area_m2,
            endpoint_area_m2: a1,
        }],
        controls.thermal.maximum_interfaces,
        check,
    )?;
    let mean_thermal = area_path.thermal();
    let routes: Vec<_> = area_path.routes().iter().collect();
    let network = FrictionHeatNetworkStep::new(
        mean_thermal,
        target,
        &routes,
        controls.thermal.maximum_interfaces,
        check,
    )?;
    let solved = solve_frozen_cohesive_heat(
        mean_thermal,
        &network,
        &[contact.temperature_guess_k],
        &[loss],
        controls.thermal,
        energy_guess_j,
        resolve,
        gate,
    )?;
    let heat = &solved.heat;
    let native = &mechanical.report;
    let combined = (Dd::from_f64(native.kinetic_energy_change_j)
        + Dd::from_f64(native.normal_storage_change_j)
        + Dd::from_f64(heat.thermal.stored_energy_change_j)
        - Dd::from_f64(native.external_work_j)
        - Dd::from_f64(heat.external_heat_j))
    .to_f64();
    let summed = (Dd::from_f64(native.audit.total_energy_residual_j)
        + Dd::from_f64(heat.total_energy_residual_j))
    .to_f64();
    if heat.dissipated_j != loss
        || [combined, summed]
            .into_iter()
            .any(|r| !r.is_finite() || r.abs() > controls.thermal.combined_energy_tolerance_j)
    {
        return Err(RigidDmtHeatError::Invalid(
            "original DMT normal/shear motion and finite heat do not balance",
        ));
    }
    let (thermal, mut routes_after) = area_path.finish(solved.thermal, check)?;
    // This adapter admits exactly one original contact throughout.
    let route_after = routes_after.remove(0);
    check()?;
    Ok(HeatedRigidDmtCandidate {
        bodies_after: mechanical.bodies_after,
        state_after: mechanical.state_after,
        thermal,
        report: HeatedRigidDmtReport {
            mechanical: mechanical.report,
            heat: solved.heat,
            source: solved.source,
            nonlinear: solved.nonlinear,
            thermal_trials: solved.trials,
            target_time_s: target,
            areas,
            route_after,
            combined_energy_residual_j: combined,
            summed_energy_residual_j: summed,
        },
    })
}

/// Publish native momenta, original DMT history and finite temperatures once.
/// Every source, geometry, caloric, budget or cancellation failure rolls back
/// all three accepted inputs. Callback payloads remain caller-owned until return.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_dmt_heat<'t, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDmtHeatContact<'_>,
    history: &mut DmtCoupledState,
    wrenches: [Wrench; 2],
    thermal: &mut DepositionSimulation<'t>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidDmtHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<RigidDmtHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidDmtReport<P>, RigidDmtHeatError<E>> {
    let candidate = propose_rigid_dmt_heat(
        bodies,
        properties,
        contact,
        *history,
        wrenches,
        thermal,
        accepted_time_s,
        dt_s,
        controls,
        energy_guess_j,
        resolve,
        gate,
    )?;
    if gate.is_requested() {
        return Err(DepositionError::Cancelled.into());
    }
    bodies[contact.contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.contact.binding.body_b] = candidate.bodies_after[1];
    *history = candidate.state_after;
    *thermal = candidate.thermal;
    Ok(candidate.report)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn validate_binding<E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    prepared: &PreparedRigidDmtPair<'_>,
    contact: RigidDmtHeatContact<'_>,
    thermal: &DepositionSimulation<'_>,
    area: f64,
    tolerance: f64,
) -> Result<(), RigidDmtHeatError<E>> {
    let source = contact.contact.normal.normal_state().elastic().interface();
    if contact.bodies[0].body_id == contact.bodies[1].body_id
        || contact
            .bodies
            .iter()
            .zip(source.surface_materials())
            .any(|(body, original)| body.mechanical_material != original)
        || (area - contact.route.contact_area_m2()).abs()
            > 64.0 * f64::EPSILON * area.max(contact.route.contact_area_m2())
    {
        return Err(RigidDmtHeatError::Invalid(
            "original ordered materials and matching incoming solid area required",
        ));
    }
    let geometry = contact.route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or(RigidDmtHeatError::Invalid(
            "thermal route has no shared face",
        ))?;
    let normal = contact.contact.frame.normal();
    let normal = Vec3::new(normal[0], normal[1], normal[2]);
    let ids = [
        contact.contact.binding.body_a,
        contact.contact.binding.body_b,
    ];
    for (side, id) in contact.route.cells().into_iter().enumerate() {
        let cell = thermal.cells().get(&id).ok_or(RigidDmtHeatError::Invalid(
            "missing original thermal receiver",
        ))?;
        if cell.body_id != contact.bodies[side].body_id
            || cell.material_id != contact.bodies[side].thermal_material
        {
            return Err(RigidDmtHeatError::Invalid(
                "thermal receiver changed original A/B body or material assignment",
            ));
        }
        let point = bodies[ids[side]]
            .point_kinematics(properties[ids[side]], prepared.reaction_arms_body()[side])
            .map_err(RigidDmtError::from)?
            .point_world;
        let g = geometry[side];
        let axis = Vec3::new(
            g.end_m[0] - g.start_m[0],
            g.end_m[1] - g.start_m[1],
            g.end_m[2] - g.start_m[2],
        );
        if magnitude(point.sub(Vec3::new(centre[0], centre[1], centre[2]))) > tolerance
            || magnitude(axis.cross(normal)) > tolerance
        {
            return Err(RigidDmtHeatError::Invalid(
                "native reaction plane and thermal face centre/normal must coincide",
            ));
        }
        // The native normal is ordered A -> B. Opposing thermal faces alone
        // also admit the reversed physical assignment. Use the actual interior
        // half-vector, independent of each cell's start/end parameterization.
        let interior = if g.start_m == centre {
            g.end_m
        } else {
            g.start_m
        };
        let inward = Vec3::new(
            0.5 * (interior[0] - centre[0]),
            0.5 * (interior[1] - centre[1]),
            0.5 * (interior[2] - centre[2]),
        );
        let signed_depth = inward.dot(normal) * if side == 0 { -1.0 } else { 1.0 };
        if !signed_depth.is_finite() || signed_depth <= tolerance {
            return Err(RigidDmtHeatError::Invalid(
                "thermal A/B interiors must retain the native A-to-B normal ordering",
            ));
        }
    }
    Ok(())
}

fn magnitude(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
