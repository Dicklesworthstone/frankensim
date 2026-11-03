//! Coupled normal and tangential impulses on complete solid trajectories.
//!
//! The original inverse mass or complete eigenbasis supplies every entry of
//! J M^-1 J^T. The existing fs-contact owner solves the non-associated Coulomb
//! law, including normal/shear and cross-contact coupling. Its candidate is
//! applied through the original impulse trial; actual contact velocities,
//! final-slip traction, original momentum and work must pass before acceptance.

use super::spectral::CompleteSpectralSolid;
use super::{ExactLinearSolid, LinearSolidImpulseError, LinearSolidImpulseReport, impulse};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResiduals,
    FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_math::det;

/// One declared contact in an orthonormal normal-first point frame.
///
/// The three dimensionless covectors follow `free_dofs()` for complete solids,
/// or all represented coordinates for a surrounding supported owner. They map
/// native velocity to [normal, tangent 1, tangent 2] m/s; their transpose maps
/// contact impulse to nodal kg m/s. Geometry owns the unit frame, shape
/// functions, event time and completeness of the contact set. These properties
/// cannot be inferred from arbitrary algebraic rows.
#[derive(Clone, Copy, Debug)]
pub struct SolidFrictionalContact<'a> {
    /// Stable geometry-owned identity; duplicate IDs refuse.
    pub id: u64,
    pub velocity_weights: [&'a [f64]; 3],
    /// Nonnegative separation [m]; penetration is not repaired by an impulse.
    pub gap_m: f64,
    pub activation_tolerance_m: f64,
    /// Opposing-surface velocity in the SAME contact frame, held throughout
    /// the instantaneous event. Relative velocity is J*v minus this value.
    pub prescribed_surface_velocity_mps: [f64; 3],
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// Native law controls plus independent actual-state and adapter work limits.
#[derive(Clone, Copy, Debug)]
pub struct SolidFrictionalImpactControls {
    pub contact: FrictionalContactConfig,
    pub maximum_impulse_error_kg_m_s: f64,
    /// Maximum component defect between actual and predicted contact velocity.
    /// This does not relax the contact law's physical residual tolerances.
    pub maximum_contact_speed_error_mps: f64,
    pub maximum_work_error_j: f64,
    /// Conservative adapter scratch; original impulse audits and the contact
    /// solver's separately bounded, at-most-63-row scratch retain their owners.
    pub maximum_workspace_bytes: usize,
    /// Scalar-product bound for full contact mobility, velocity and transpose
    /// assembly. Excludes the native impulse audit and solver iterations.
    pub maximum_response_multiply_terms: u64,
}

/// Refusal publishes no physical, modal, source or cumulative-account changes.
#[derive(Debug)]
pub enum SolidFrictionalImpactError {
    Invalid(&'static str),
    Budget {
        resource: &'static str,
        required: u128,
        maximum: u128,
    },
    Native(LinearSolidImpulseError),
    Contact(FrictionalContactError),
    Speed {
        contact: usize,
        component: usize,
        defect_mps: f64,
        allowance_mps: f64,
    },
    Work {
        defect_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for SolidFrictionalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native solid Coulomb impact refused: {self:?}")
    }
}
impl core::error::Error for SolidFrictionalImpactError {}
impl From<LinearSolidImpulseError> for SolidFrictionalImpactError {
    fn from(error: LinearSolidImpulseError) -> Self {
        if matches!(error, LinearSolidImpulseError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Native(error)
        }
    }
}
impl From<FrictionalContactError> for SolidFrictionalImpactError {
    fn from(error: FrictionalContactError) -> Self {
        if matches!(error, FrictionalContactError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Contact(error)
        }
    }
}

/// Shared actual-state receipt for complete and supported solid consumers.
///
/// Retain `law.determinacy` and its total-work witness. A velocity-only result
/// supplies one checked impulse representative, not unique individual loads
/// or local source/heat allocations. Prescribed solid-support reaction work,
/// when present, remains separate in that owner's original impulse receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct SolidFrictionalContactImpulseReport {
    pub law: FrictionalContactResult,
    pub relative_velocities_after_mps: Vec<[f64; 3]>,
    pub maximum_contact_speed_error_mps: f64,
    /// Coulomb residuals reconstructed from the actual native endpoint.
    pub measured_residuals: FrictionalContactResiduals,
    /// Sum p dot opposing-surface velocity for the selected representative [J].
    pub prescribed_surface_work_j: f64,
    /// Sum p dot (relative velocity before + actual relative velocity after)/2.
    pub relative_work_j: f64,
    /// Nonnegative total impact loss, never a friction/normal heat partition.
    pub dissipated_energy_j: f64,
    pub passivity_residual_j: f64,
    /// Native energy tolerance plus 64 eps times absolute midpoint work terms.
    pub allowed_relative_work_j: f64,
    pub relative_work_error_j: f64,
    pub work_balance_error_j: f64,
}

/// One accepted event; contact loss is not double-counted as native viscosity.
#[derive(Clone, Debug)]
pub struct SolidFrictionalImpactReport {
    pub contact: SolidFrictionalContactImpulseReport,
    pub impulse: LinearSolidImpulseReport,
}

type Result<T> = core::result::Result<T, SolidFrictionalImpactError>;

fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| SolidFrictionalImpactError::Cancelled)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolidFrictionalImpactError::Invalid(
            "unrepresentable contact mobility, velocity, impulse or work",
        ))
    }
}
fn rounded(value: Dd) -> Result<f64> {
    finite(value.hi)?;
    finite(value.lo)?;
    finite(value.to_f64())
}
fn dot(a: &[f64], b: &[f64], cx: &Cx<'_>) -> Result<Dd> {
    let mut sum = Dd::ZERO;
    for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
        if i % 256 == 0 {
            poll(cx)?;
        }
        sum = sum + Dd::from_f64(a) * Dd::from_f64(b);
    }
    rounded(sum)?;
    Ok(sum)
}
fn relative_velocity(
    contact: &SolidFrictionalContact<'_>,
    velocity: &[f64],
    cx: &Cx<'_>,
) -> Result<[f64; 3]> {
    let mut result = [0.0; 3];
    for (a, value) in result.iter_mut().enumerate() {
        *value = rounded(
            dot(contact.velocity_weights[a], velocity, cx)?
                - Dd::from_f64(contact.prescribed_surface_velocity_mps[a]),
        )?;
    }
    Ok(result)
}
fn bound(resource: &'static str, required: u128, maximum: u128) -> Result<()> {
    if required > maximum {
        Err(SolidFrictionalImpactError::Budget {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}

pub(super) struct FrictionalImpulseProposal {
    pub(super) impulse_kg_m_s: Vec<f64>,
    law: FrictionalContactResult,
    relative_velocities_before_mps: Vec<[f64; 3]>,
    order: Vec<usize>,
}

/// Native velocity and mobility action use the SAME original coordinate order.
/// Supported owners retain real prescribed velocities with zero incremental
/// mobility there, and recover cross-mass reactions in their original trial.
pub(super) fn propose(
    velocity_mps: &[f64],
    contacts: &[SolidFrictionalContact<'_>],
    controls: SolidFrictionalImpactControls,
    mut inverse_mass_action: impl FnMut(
        &[f64],
        &Cx<'_>,
    ) -> core::result::Result<Vec<f64>, LinearSolidImpulseError>,
    cx: &Cx<'_>,
) -> Result<FrictionalImpulseProposal> {
    poll(cx)?;
    let (n, m) = (velocity_mps.len(), contacts.len());
    if n == 0 || !(1..=21).contains(&m) || m > controls.contact.maximum_contacts {
        return Err(SolidFrictionalImpactError::Invalid(
            "nonempty original state and 1..=21 admitted contact frames required",
        ));
    }
    if [
        controls.maximum_impulse_error_kg_m_s,
        controls.maximum_contact_speed_error_mps,
        controls.maximum_work_error_j,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(SolidFrictionalImpactError::Invalid(
            "finite nonnegative impulse, speed and work allowances required",
        ));
    }
    let (nn, dd) = (n as u128, (3 * m) as u128);
    bound(
        "Coulomb adapter workspace bytes",
        256 * nn + 1024 * dd + 16 * dd * dd,
        controls.maximum_workspace_bytes as u128,
    )?;
    bound(
        "Coulomb response scalar-product terms",
        (2 * dd)
            .saturating_mul(nn.saturating_mul(nn))
            .saturating_add((dd * dd + 10 * dd).saturating_mul(nn)),
        controls.maximum_response_multiply_terms as u128,
    )?;
    let mut order: Vec<_> = (0..m).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    if order
        .windows(2)
        .any(|pair| contacts[pair[0]].id == contacts[pair[1]].id)
    {
        return Err(SolidFrictionalImpactError::Invalid("duplicate contact ID"));
    }
    let mut native = Vec::with_capacity(m);
    let mut relative_velocities_before_mps = Vec::with_capacity(m);
    for contact in contacts {
        poll(cx)?;
        if contact
            .velocity_weights
            .iter()
            .any(|row| row.len() != n || row.iter().any(|v| !v.is_finite()))
            || contact
                .prescribed_surface_velocity_mps
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(SolidFrictionalImpactError::Invalid(
                "three complete finite covectors and finite opposing-surface velocities required",
            ));
        }
        let relative_velocity = relative_velocity(contact, velocity_mps, cx)?;
        relative_velocities_before_mps.push(relative_velocity);
        native.push(FrictionalContact {
            id: contact.id,
            gap: contact.gap_m,
            activation_tolerance: contact.activation_tolerance_m,
            relative_velocity,
            friction_coefficient: contact.friction_coefficient,
            normal_law: contact.normal_law,
        });
    }
    let d = 3 * m;
    let mut delassus = vec![0.0; d * d];
    for &j in &order {
        for b in 0..3 {
            poll(cx)?;
            let response = inverse_mass_action(contacts[j].velocity_weights[b], cx)?;
            if response.len() != n {
                return Err(SolidFrictionalImpactError::Invalid(
                    "native mobility shape changed",
                ));
            }
            for &i in &order {
                for a in 0..3 {
                    delassus[(3 * i + a) * d + 3 * j + b] =
                        rounded(dot(contacts[i].velocity_weights[a], &response, cx)?)?;
                }
            }
        }
    }
    let law = solve_frictional_contacts(&native, &delassus, controls.contact, || {
        cx.checkpoint().is_err()
    })?;
    let mut impulse_kg_m_s = Vec::with_capacity(n);
    for coordinate in 0..n {
        if coordinate % 256 == 0 {
            poll(cx)?;
        }
        let mut value = Dd::ZERO;
        for &i in &order {
            for a in 0..3 {
                value = value
                    + Dd::from_f64(contacts[i].velocity_weights[a][coordinate])
                        * Dd::from_f64(law.impulses[i][a]);
            }
        }
        impulse_kg_m_s.push(rounded(value)?);
    }
    poll(cx)?;
    Ok(FrictionalImpulseProposal {
        impulse_kg_m_s,
        law,
        relative_velocities_before_mps,
        order,
    })
}

/// Recheck the original Coulomb law at the actual native endpoint. The shared
/// solver remains the sole law/branch/determinacy owner; this independently
/// prevents finite-basis or native-state roundoff from bypassing its tolerances.
fn measured_law(
    proposal: &FrictionalImpulseProposal,
    contacts: &[SolidFrictionalContact<'_>],
    velocities: &[[f64; 3]],
    config: FrictionalContactConfig,
    cx: &Cx<'_>,
) -> Result<FrictionalContactResiduals> {
    let mut r = FrictionalContactResiduals::default();
    for &i in &proposal.law.active_indices {
        poll(cx)?;
        let (p, v) = (proposal.law.impulses[i], velocities[i]);
        let g = finite(v[0] - proposal.law.normal_targets[i])?;
        let radius = finite(contacts[i].friction_coefficient * p[0])?;
        let speed = finite(det::hypot(v[1], v[2]))?;
        r.normal_impulse_violation = r.normal_impulse_violation.max((-p[0]).max(0.0));
        r.normal_velocity_violation = r.normal_velocity_violation.max((-g).max(0.0));
        r.normal_complementarity = r.normal_complementarity.max(finite(p[0] * g)?.abs());
        r.cone_violation = r
            .cone_violation
            .max(finite(det::hypot(p[1], p[2]) - radius)?.max(0.0));
        let duality = rounded(
            Dd::from_f64(radius) * Dd::from_f64(speed)
                + Dd::from_f64(p[1]) * Dd::from_f64(v[1])
                + Dd::from_f64(p[2]) * Dd::from_f64(v[2]),
        )?;
        r.tangential_duality_gap = r.tangential_duality_gap.max(duality.abs());
        if speed > config.velocity_tolerance {
            let residual = finite(det::hypot(
                finite(p[1] + radius * (v[1] / speed))?,
                finite(p[2] + radius * (v[2] / speed))?,
            ))?;
            r.sliding_impulse_residual = r.sliding_impulse_residual.max(residual);
        }
    }
    if r.normal_impulse_violation > config.impulse_tolerance
        || r.normal_velocity_violation > config.velocity_tolerance
        || r.normal_complementarity > config.complementarity_tolerance
        || r.cone_violation > config.impulse_tolerance
        || r.sliding_impulse_residual > config.impulse_tolerance
        || r.tangential_duality_gap > config.complementarity_tolerance
    {
        return Err(FrictionalContactError::ResidualNotSatisfied(r).into());
    }
    Ok(r)
}

pub(super) fn audit(
    proposal: FrictionalImpulseProposal,
    contacts: &[SolidFrictionalContact<'_>],
    velocity_after_mps: &[f64],
    applied_work_j: f64,
    controls: SolidFrictionalImpactControls,
    cx: &Cx<'_>,
) -> Result<SolidFrictionalContactImpulseReport> {
    poll(cx)?;
    if velocity_after_mps.len() != proposal.impulse_kg_m_s.len()
        || contacts.len() != proposal.law.impulses.len()
    {
        return Err(SolidFrictionalImpactError::Invalid(
            "native impact candidate shape changed",
        ));
    }
    let mut relative_velocities_after_mps = vec![[0.0; 3]; contacts.len()];
    let mut maximum_contact_speed_error_mps = 0.0_f64;
    let (mut prescribed_work, mut relative_work, mut work_scale) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
    for &i in &proposal.order {
        poll(cx)?;
        let after = relative_velocity(&contacts[i], velocity_after_mps, cx)?;
        relative_velocities_after_mps[i] = after;
        for a in 0..3 {
            let defect = finite(after[a] - proposal.law.relative_velocities_after[i][a])?;
            maximum_contact_speed_error_mps = maximum_contact_speed_error_mps.max(defect.abs());
            if defect.abs() > controls.maximum_contact_speed_error_mps {
                return Err(SolidFrictionalImpactError::Speed {
                    contact: i,
                    component: a,
                    defect_mps: defect,
                    allowance_mps: controls.maximum_contact_speed_error_mps,
                });
            }
            let p = Dd::from_f64(proposal.law.impulses[i][a]);
            prescribed_work =
                prescribed_work + p * Dd::from_f64(contacts[i].prescribed_surface_velocity_mps[a]);
            let work = p
                * (Dd::from_f64(0.5 * proposal.relative_velocities_before_mps[i][a])
                    + Dd::from_f64(0.5 * after[a]));
            relative_work = relative_work + work;
            work_scale = work_scale + Dd::from_f64(rounded(work)?.abs());
        }
    }
    let measured_residuals = measured_law(
        &proposal,
        contacts,
        &relative_velocities_after_mps,
        controls.contact,
        cx,
    )?;
    let prescribed_surface_work_j = rounded(prescribed_work)?;
    let relative_work_j = rounded(relative_work)?;
    let relative_work_error_j = rounded(relative_work - Dd::from_f64(proposal.law.relative_work))?;
    let work_balance_error_j =
        rounded(Dd::from_f64(applied_work_j) - prescribed_work - relative_work)?;
    for defect_j in [relative_work_error_j, work_balance_error_j] {
        if defect_j.abs() > controls.maximum_work_error_j {
            return Err(SolidFrictionalImpactError::Work {
                defect_j,
                allowance_j: controls.maximum_work_error_j,
            });
        }
    }
    let allowed_relative_work_j =
        finite(controls.contact.energy_tolerance + 64.0 * f64::EPSILON * rounded(work_scale)?)?;
    if relative_work_j > allowed_relative_work_j {
        return Err(SolidFrictionalImpactError::Work {
            defect_j: relative_work_j,
            allowance_j: allowed_relative_work_j,
        });
    }
    poll(cx)?;
    Ok(SolidFrictionalContactImpulseReport {
        law: proposal.law,
        relative_velocities_after_mps,
        maximum_contact_speed_error_mps,
        measured_residuals,
        prescribed_surface_work_j,
        relative_work_j,
        dissipated_energy_j: (-relative_work_j).max(0.0),
        passivity_residual_j: relative_work_j.max(0.0),
        allowed_relative_work_j,
        relative_work_error_j,
        work_balance_error_j,
    })
}

impl ExactLinearSolid {
    /// Resolve simultaneous normal and two-axis Coulomb impulses on the
    /// original prepared trajectory. All original cross-mass, normal/shear and
    /// intercontact mobility remains in one fs-contact solve. There is no
    /// sequential normal/friction clipping or added contact compliance.
    ///
    /// Native momentum/work, measured contact velocity and the physical Coulomb
    /// residuals gate a private trial. No positions, clocks, source/material
    /// preparation or earlier viscosity losses change. Prescribed surface work
    /// is distinct from passive relative loss; no heat destination is invented.
    /// Geometry, event discovery and sustained contact are separate owners.
    pub fn apply_frictional_impacts(
        &mut self,
        contacts: &[SolidFrictionalContact<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidFrictionalImpactReport> {
        let proposal = propose(
            self.velocity_mps(),
            contacts,
            controls,
            |row, cx| impulse::action(&self.inverse_mass, row, cx),
            cx,
        )?;
        let candidate = self.trial_momentum_impulse(
            &proposal.impulse_kg_m_s,
            controls.maximum_impulse_error_kg_m_s,
            cx,
        )?;
        let contact = audit(
            proposal,
            contacts,
            &candidate.velocity_mps,
            candidate.report.energy.supplied_work_j,
            controls,
            cx,
        )?;
        Ok(SolidFrictionalImpactReport {
            contact,
            impulse: self.accept_momentum_impulse(candidate),
        })
    }
}

impl CompleteSpectralSolid {
    /// The same simultaneous Coulomb event on every admitted rigid and elastic
    /// mode, without truncation, refactorization or cloning prepared operators.
    /// Actual original-coordinate readback still gates a rounded input basis.
    /// See the dense owner's method for work, geometry and no-claim boundaries.
    pub fn apply_frictional_impacts(
        &mut self,
        contacts: &[SolidFrictionalContact<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidFrictionalImpactReport> {
        let proposal = propose(
            self.velocity_mps(),
            contacts,
            controls,
            |row, cx| self.momentum_response(row, cx),
            cx,
        )?;
        let candidate = self.trial_momentum_impulse(
            &proposal.impulse_kg_m_s,
            controls.maximum_impulse_error_kg_m_s,
            cx,
        )?;
        let contact = audit(
            proposal,
            contacts,
            &candidate.velocity_mps,
            candidate.report.energy.supplied_work_j,
            controls,
            cx,
        )?;
        Ok(SolidFrictionalImpactReport {
            contact,
            impulse: self.accept_momentum_impulse(candidate),
        })
    }
}
