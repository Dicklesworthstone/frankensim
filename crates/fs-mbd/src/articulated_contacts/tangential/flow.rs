//! Continuous dry friction on the original coupled articulated inertia.
//!
//! The shared contact owner resolves every exact arrest, holding reaction and
//! breakaway under held Coulomb capacities and drives. Its full native mobility
//! is `J M^-1 J^T`, including joint and shared-body recoil. The integrated
//! impulses update the original ABA states once; continuous travel and work
//! retain all intervening events instead of pairing only the outer velocities.
//!
//! These are physical scalar guides/contact ports, not independent components
//! of a two-dimensional Coulomb disk. Normal contact, evolving geometry and
//! material conditions remain with the enclosing mechanical/thermal owners.

use super::{
    ArticulatedContactError, ArticulatedContactSystem, ArticulatedTangentialContact, EventResponse,
    EventWorkspace, finite, poll, prepare_row,
};
use crate::articulated::FreeFloatingImpulseResponse;
use fs_contact::tangential::impact::multiple::relaxation::stiction::{
    CoupledStictionFlowControls, CoupledStictionFlowError, CoupledStictionFlowInput,
    CoupledStictionFlowReport, MAX_STICTION_FLOW_PORTS, advance_coupled_stiction,
};
use fs_ga::Vec3;
use fs_math::dd::Dd;
use fs_tribo::FrictionLaw;

/// A source-declared dry-friction port and an optional held actuator force.
#[derive(Clone, Copy, Debug)]
pub struct ArticulatedDryFrictionSite<'a> {
    /// Unique stable identity; reports and matrix rows use ascending ID order.
    pub id: u64,
    pub contact: ArticulatedTangentialContact,
    /// Prescribed compressive load, N; this does not enforce normal kinematics.
    pub normal_force_n: f64,
    /// The original constant Coulomb law, including distinct breakaway/sliding
    /// coefficients. Rate-dependent or history-bearing laws are not frozen here.
    pub law: &'a FrictionLaw,
    /// Independent actuator force along this port, N. Positive acts on `point`
    /// and negatively on a finite counterpart. It is not extra friction capacity.
    pub drive_force_n: f64,
}

/// Original flow budgets plus explicit native configuration/response gates.
#[derive(Clone, Copy, Debug)]
pub struct ArticulatedDryFrictionControls {
    pub flow: CoupledStictionFlowControls,
    /// Largest interval for which the caller admits the frozen projection, s.
    pub maximum_duration_s: f64,
    /// Bounds each actual port's complete piecewise-affine path, m/s.
    pub maximum_slip_speed_m_s: f64,
    /// Common-point admission; no witness is moved or projected, m.
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    /// Actual native kinetic change minus source/boundary work plus loss, J.
    pub energy_tolerance_j: f64,
}

/// Atomic native and continuous-path candidate for one physical phase.
#[derive(Clone, Debug)]
#[must_use]
pub struct ArticulatedDryFrictionResult {
    /// Native body/joint velocity changes in input-system order. Poses are fixed.
    pub responses: Vec<FreeFloatingImpulseResponse>,
    pub contact_ids: Vec<u64>,
    /// Complete row-major native mobility, in canonical contact-ID order, 1/kg.
    pub mobility_per_kg: Vec<f64>,
    /// Positive relative motion is receiving point minus counterpart along the
    /// authored tangent. This is the opposite of `TangentialLaw`'s slip convention.
    /// Per-port continuous work is retained; only its loss may feed a heat owner.
    pub path: CoupledStictionFlowReport,
    /// Independently measured native endpoint, including roundoff residuals, m/s.
    pub relative_speeds_after_m_s: Vec<f64>,
    pub kinetic_energy_change_j: f64,
    /// Signed work supplied by prescribed plane material, including drive, J.
    pub support_work_j: f64,
    pub energy_balance_residual_j: f64,
    pub linear_momentum_residual: Vec3,
    pub angular_momentum_residual: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
}

/// No partial native update or force-law state publishes on refusal.
#[derive(Debug)]
pub enum ArticulatedDryFrictionError {
    Mechanics(ArticulatedContactError),
    Flow(CoupledStictionFlowError),
    UnsupportedLaw { contact_id: u64 },
}
impl From<ArticulatedContactError> for ArticulatedDryFrictionError {
    fn from(error: ArticulatedContactError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<CoupledStictionFlowError> for ArticulatedDryFrictionError {
    fn from(error: CoupledStictionFlowError) -> Self {
        Self::Flow(error)
    }
}
impl core::fmt::Display for ArticulatedDryFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous articulated dry friction: {self:?}")
    }
}
impl std::error::Error for ArticulatedDryFrictionError {}

/// Advance up to six coupled scalar friction ports through exact stop events.
///
/// Capacities come from the original Coulomb law and supplied normal load. All
/// resting ports resolve together; incomparable static/kinetic continuations
/// retain the original owner's ambiguity refusal. There is no capture speed,
/// diagonal-mass approximation, impulse clipping or velocity reversal at arrest.
/// Applied force work, moving-support work and physical friction loss remain
/// distinct, and the original native momentum/energy audit must close.
///
/// Bodies, joint positions, contact arms and source data are held during this
/// phase. Returned travel is the frozen-projection trajectory, not an SE(3)
/// pose update. Geometry/load changes require refinement by the enclosing
/// integrator. Numerical resting residuals are reported, never silently reset.
/// Elastic contacts continue to use the retained-history tangential owner.
pub fn advance_articulated_dry_friction(
    systems: &[ArticulatedContactSystem<'_>],
    sites: &[ArticulatedDryFrictionSite<'_>],
    duration_s: f64,
    controls: ArticulatedDryFrictionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedDryFrictionResult, ArticulatedDryFrictionError> {
    poll(&mut cancelled)?;
    if systems.is_empty()
        || systems.len() > 128
        || sites.len() > MAX_STICTION_FLOW_PORTS
        || !duration_s.is_finite()
        || duration_s < 0.0
        || duration_s > controls.maximum_duration_s
        || [
            controls.maximum_duration_s,
            controls.maximum_slip_speed_m_s,
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(ArticulatedContactError::Invalid("bounded dry-friction phase").into());
    }
    controls
        .flow
        .flow
        .validate()
        .map_err(CoupledStictionFlowError::from)?;
    let mut order: Vec<_> = (0..sites.len()).collect();
    order.sort_by_key(|&i| sites[i].id);
    if order
        .windows(2)
        .any(|pair| sites[pair[0]].id == sites[pair[1]].id)
    {
        return Err(ArticulatedContactError::Invalid("duplicate dry-friction contact ID").into());
    }
    let mut workspace = EventWorkspace::new(systems, &mut cancelled)?;
    let mut rows = Vec::with_capacity(sites.len());
    let mut speeds = Vec::with_capacity(sites.len());
    let mut static_capacity = Vec::with_capacity(sites.len());
    let mut kinetic_capacity = Vec::with_capacity(sites.len());
    let mut drives = Vec::with_capacity(sites.len());
    for &index in &order {
        poll(&mut cancelled)?;
        let site = sites[index];
        site.law.validate().map_err(ArticulatedContactError::from)?;
        if !matches!(site.law, FrictionLaw::Coulomb { .. }) {
            return Err(ArticulatedDryFrictionError::UnsupportedLaw {
                contact_id: site.id,
            });
        }
        if !site.normal_force_n.is_finite() || site.normal_force_n < 0.0 {
            return Err(
                ArticulatedContactError::Invalid("nonnegative normal load required").into(),
            );
        }
        let mu_s = site
            .law
            .static_coefficient()
            .map_err(ArticulatedContactError::from)?;
        let mu_k = site
            .law
            .zero_speed_kinetic_coefficient()
            .map_err(ArticulatedContactError::from)?;
        let capacity = |mu: f64| -> Result<f64, ArticulatedContactError> {
            let force = finite(mu * site.normal_force_n)?;
            if force == 0.0 && mu != 0.0 && site.normal_force_n != 0.0 {
                return Err(ArticulatedContactError::Unresolved(
                    "representable friction capacity",
                ));
            }
            Ok(force)
        };
        static_capacity.push(capacity(mu_s)?);
        kinetic_capacity.push(capacity(mu_k)?);
        drives.push(finite(site.drive_force_n)?);
        let row = prepare_row(site.contact, &workspace, controls.geometry_tolerance_m)?;
        speeds.push(finite(row.velocity(&workspace.before))?);
        rows.push(row);
    }
    let mobility = workspace.delassus(systems, &rows, &mut cancelled)?;
    let path = advance_coupled_stiction(
        &CoupledStictionFlowInput {
            velocity_before_m_s: &speeds,
            mobility_per_kg: &mobility,
            static_capacity_n: &static_capacity,
            kinetic_capacity_n: &kinetic_capacity,
            drive_n: &drives,
            duration_s,
        },
        controls.flow,
        &mut cancelled,
    )?;
    let count = sites.len();
    if path.flow.maximum_speed_m_s[..count]
        .iter()
        .any(|&v| v > controls.maximum_slip_speed_m_s)
    {
        return Err(ArticulatedContactError::Unresolved("dry-friction path speed limit").into());
    }
    let impulses: Vec<_> = (0..count)
        .map(|i| {
            finite(
                (Dd::from_f64(path.flow.friction_impulse_ns[i])
                    + Dd::from_f64(path.flow.drive_impulse_ns[i]))
                .to_f64(),
            )
        })
        .collect::<Result<_, _>>()?;
    let mechanics = workspace.apply_with_source_work(
        systems,
        &rows,
        EventResponse {
            impulses: &impulses,
            predicted_speeds: &path.flow.velocity_after_m_s[..count],
            stored_energy_change: 0.0,
            dissipated_energy: path.flow.dissipated_j,
            velocity_tolerance: controls.velocity_tolerance_m_s,
            energy_tolerance: controls.energy_tolerance_j,
        },
        path.flow.total_drive_work_j,
        &mut cancelled,
    )?;
    if mechanics
        .speeds_after
        .iter()
        .any(|v| v.abs() > controls.maximum_slip_speed_m_s)
    {
        return Err(ArticulatedContactError::Unresolved("native dry-friction speed limit").into());
    }
    poll(&mut cancelled)?;
    Ok(ArticulatedDryFrictionResult {
        responses: mechanics.responses,
        contact_ids: order.iter().map(|&i| sites[i].id).collect(),
        mobility_per_kg: mobility,
        path,
        relative_speeds_after_m_s: mechanics.speeds_after,
        kinetic_energy_change_j: mechanics.kinetic_energy_change,
        support_work_j: mechanics.plane_work,
        energy_balance_residual_j: mechanics.energy_balance_residual,
        linear_momentum_residual: mechanics.linear_momentum_residual,
        angular_momentum_residual: mechanics.angular_momentum_residual,
        endpoint_velocity_residual_m_s: mechanics.endpoint_velocity_residual,
    })
}
