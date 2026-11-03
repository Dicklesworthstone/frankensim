//! Native complete-solid response to a declared simultaneous normal impact.
//!
//! Geometry supplies contact covectors and an already admissible event. The
//! original prepared inverse/full eigenbasis supplies J M^-1 J^T; fs-contact
//! owns complementarity, restitution, rank admission and relative-work limits.
//! The original impulse owner then checks M delta-v. Actual contact velocities
//! and work must pass before either physical or modal state is published.

use super::spectral::CompleteSpectralSolid;
use super::{ExactLinearSolid, LinearSolidImpulseError, LinearSolidImpulseReport, impulse};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_exec::Cx;
use fs_math::dd::Dd;

/// One contact normal covector in the native owner's original coordinate order.
/// For complete solids that is `free_dofs()`; a surrounding supported owner
/// may use its full represented order. The dimensionless row maps nodal m/s to
/// normal m/s, and its transpose maps a normal impulse to nodal kg m/s.
/// Geometry owns its unit-normal/shape-function construction and active-set
/// completeness. This adapter does not infer geometry from an arbitrary row.
#[derive(Clone, Copy, Debug)]
pub struct SolidNormalContact<'a> {
    pub normal_velocity_weights: &'a [f64],
    /// Nonnegative separation [m]; penetration always refuses.
    pub gap_m: f64,
    pub activation_tolerance_m: f64,
    /// Prescribed opposing-surface speed in the same separating direction.
    /// Relative approach is J*v minus this speed, with a negative sign.
    /// It is held over the instantaneous event; no position is advanced.
    pub prescribed_surface_speed_mps: f64,
    pub restitution: f64,
}

/// Explicit numerical and capacity limits, independent of physical restitution.
#[derive(Clone, Copy, Debug)]
pub struct SolidNormalImpactControls {
    pub contact: MultiImpactConfig,
    pub maximum_impulse_error_kg_m_s: f64,
    pub maximum_contact_speed_error_mps: f64,
    pub maximum_work_error_j: f64,
    /// Conservative adapter scratch. Existing native impulse audits and the
    /// shared solver's bounded (at most 64-contact) workspace remain separate.
    pub maximum_workspace_bytes: usize,
    /// Conservative scalar-product term bound for contact mobility assembly
    /// and its force transpose; excludes the existing native impulse audit
    /// and the contact solver's separately bounded iterations.
    pub maximum_response_multiply_terms: u64,
}

/// A refused event leaves the original trajectory and all accounts unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum SolidNormalImpactError {
    Invalid(&'static str),
    Budget {
        resource: &'static str,
        required: u128,
        maximum: u128,
    },
    Native(LinearSolidImpulseError),
    Contact(MultiImpactError),
    Speed {
        contact: usize,
        defect_mps: f64,
        allowance_mps: f64,
    },
    Work {
        defect_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for SolidNormalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native solid normal impact refused: {self:?}")
    }
}
impl core::error::Error for SolidNormalImpactError {}
impl From<LinearSolidImpulseError> for SolidNormalImpactError {
    fn from(error: LinearSolidImpulseError) -> Self {
        if matches!(error, LinearSolidImpulseError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Native(error)
        }
    }
}
impl From<MultiImpactError> for SolidNormalImpactError {
    fn from(error: MultiImpactError) -> Self {
        if matches!(error, MultiImpactError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Contact(error)
        }
    }
}

/// Contact-law and measured-state account shared with supported solid adapters.
/// A redundant law's impulse vector is only a checked representative; inspect
/// `law.determinacy` before attributing individual contact loads or source work.
#[derive(Clone, Debug, PartialEq)]
pub struct SolidContactImpulseReport {
    pub law: MultiImpactResult,
    /// Measured from the actual accepted native velocity, including inactive rows.
    pub relative_normal_speeds_after_mps: Vec<f64>,
    pub maximum_contact_speed_error_mps: f64,
    pub maximum_complementarity_speed_residual_mps: f64,
    /// Sum p*u_surface for the selected impulse representative [J]. A supported
    /// adapter retains constrained-reaction work in its own native receipt.
    /// Redundant individual source allocations need additional range admission.
    pub prescribed_surface_work_j: f64,
    /// Sum p*(relative_speed_before + measured_after)/2 [J].
    pub relative_work_j: f64,
    /// Nonnegative lost relative mechanical energy; never assigned to heat here.
    pub dissipated_energy_j: f64,
    pub passivity_residual_j: f64,
    /// Actual relative work minus the contact solver's predicted relative work.
    pub relative_work_error_j: f64,
    /// Native applied impulse work minus prescribed-surface and relative work.
    pub work_balance_error_j: f64,
}

/// One accepted event. Contact loss is not added again to native dissipation:
/// the original trajectory already accounts for the supplied impulse work.
#[derive(Clone, Debug)]
pub struct SolidNormalImpactReport {
    pub contact: SolidContactImpulseReport,
    pub impulse: LinearSolidImpulseReport,
}

type Result<T> = core::result::Result<T, SolidNormalImpactError>;

fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| SolidNormalImpactError::Cancelled)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolidNormalImpactError::Invalid(
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
    let mut value = Dd::ZERO;
    for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
        if i % 256 == 0 {
            poll(cx)?;
        }
        value = value + Dd::from_f64(a) * Dd::from_f64(b);
    }
    rounded(value)?;
    Ok(value)
}
fn speed(contact: &SolidNormalContact<'_>, velocity: &[f64], cx: &Cx<'_>) -> Result<f64> {
    rounded(
        dot(contact.normal_velocity_weights, velocity, cx)?
            - Dd::from_f64(contact.prescribed_surface_speed_mps),
    )
}
fn bound(resource: &'static str, required: u128, maximum: u128) -> Result<()> {
    if required > maximum {
        Err(SolidNormalImpactError::Budget {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}

pub(super) struct NormalImpulseProposal {
    pub(super) impulse_kg_m_s: Vec<f64>,
    law: MultiImpactResult,
    relative_speeds_before_mps: Vec<f64>,
}

/// The velocity and returned mobility live in the SAME original coordinates.
/// Supported adapters return zero response on prescribed rows, retaining all
/// cross-mass coupling in their native trial/reaction owner.
pub(super) fn propose(
    velocity_mps: &[f64],
    contacts: &[SolidNormalContact<'_>],
    controls: SolidNormalImpactControls,
    mut inverse_mass_action: impl FnMut(
        &[f64],
        &Cx<'_>,
    ) -> core::result::Result<Vec<f64>, LinearSolidImpulseError>,
    cx: &Cx<'_>,
) -> Result<NormalImpulseProposal> {
    poll(cx)?;
    let (n, m) = (velocity_mps.len(), contacts.len());
    if n == 0 || m == 0 || m > 64 || m > controls.contact.solver.max_variables {
        return Err(SolidNormalImpactError::Invalid(
            "nonempty original state and 1..=64 admitted contact rows required",
        ));
    }
    if [
        controls.maximum_impulse_error_kg_m_s,
        controls.maximum_contact_speed_error_mps,
        controls.maximum_work_error_j,
        controls.contact.energy_tolerance,
        controls.contact.solver.gradient_tolerance,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(SolidNormalImpactError::Invalid(
            "finite nonnegative impulse, speed and work allowances required",
        ));
    }
    let (nn, mm) = (n as u128, m as u128);
    bound(
        "contact adapter workspace bytes",
        256 * nn + 512 * mm + 16 * mm * mm,
        controls.maximum_workspace_bytes as u128,
    )?;
    bound(
        "contact response scalar-product terms",
        (2 * mm)
            .saturating_mul(nn.saturating_mul(nn))
            .saturating_add((mm * mm + 4 * mm).saturating_mul(nn)),
        controls.maximum_response_multiply_terms as u128,
    )?;
    let mut native = Vec::with_capacity(m);
    let mut relative_speeds_before_mps = Vec::with_capacity(m);
    for contact in contacts {
        poll(cx)?;
        if contact.normal_velocity_weights.len() != n
            || contact
                .normal_velocity_weights
                .iter()
                .any(|x| !x.is_finite())
            || !contact.prescribed_surface_speed_mps.is_finite()
        {
            return Err(SolidNormalImpactError::Invalid(
                "complete finite contact covectors and prescribed surface speeds required",
            ));
        }
        let relative_normal_speed = speed(contact, velocity_mps, cx)?;
        relative_speeds_before_mps.push(relative_normal_speed);
        native.push(NormalContactImpact {
            gap: contact.gap_m,
            activation_tolerance: contact.activation_tolerance_m,
            relative_normal_speed,
            restitution: contact.restitution,
        });
    }
    let mut delassus = vec![0.0; m * m];
    for (j, contact) in contacts.iter().enumerate() {
        poll(cx)?;
        let response = inverse_mass_action(contact.normal_velocity_weights, cx)?;
        if response.len() != n {
            return Err(SolidNormalImpactError::Invalid(
                "native mobility shape changed",
            ));
        }
        for (i, row) in contacts.iter().enumerate() {
            delassus[i * m + j] = rounded(dot(row.normal_velocity_weights, &response, cx)?)?;
        }
    }
    let law = solve_normal_impacts(&native, &delassus, controls.contact, || {
        cx.checkpoint().is_err()
    })?;
    let mut impulse_kg_m_s = Vec::with_capacity(n);
    for coordinate in 0..n {
        if coordinate % 256 == 0 {
            poll(cx)?;
        }
        let value =
            contacts
                .iter()
                .zip(&law.normal_impulses)
                .fold(Dd::ZERO, |sum, (contact, &p)| {
                    sum + Dd::from_f64(contact.normal_velocity_weights[coordinate])
                        * Dd::from_f64(p)
                });
        impulse_kg_m_s.push(rounded(value)?);
    }
    poll(cx)?;
    Ok(NormalImpulseProposal {
        impulse_kg_m_s,
        law,
        relative_speeds_before_mps,
    })
}

pub(super) fn audit(
    proposal: NormalImpulseProposal,
    contacts: &[SolidNormalContact<'_>],
    velocity_after_mps: &[f64],
    applied_work_j: f64,
    controls: SolidNormalImpactControls,
    cx: &Cx<'_>,
) -> Result<SolidContactImpulseReport> {
    poll(cx)?;
    if velocity_after_mps.len() != proposal.impulse_kg_m_s.len()
        || contacts.len() != proposal.law.normal_impulses.len()
    {
        return Err(SolidNormalImpactError::Invalid(
            "native impact candidate shape changed",
        ));
    }
    let mut relative_normal_speeds_after_mps = Vec::with_capacity(contacts.len());
    let mut maximum_contact_speed_error_mps = 0.0_f64;
    let (mut prescribed_work, mut relative_work) = (Dd::ZERO, Dd::ZERO);
    for (i, contact) in contacts.iter().enumerate() {
        poll(cx)?;
        let after = speed(contact, velocity_after_mps, cx)?;
        let defect = finite(after - proposal.law.relative_normal_speeds_after[i])?;
        maximum_contact_speed_error_mps = maximum_contact_speed_error_mps.max(defect.abs());
        if defect.abs() > controls.maximum_contact_speed_error_mps {
            return Err(SolidNormalImpactError::Speed {
                contact: i,
                defect_mps: defect,
                allowance_mps: controls.maximum_contact_speed_error_mps,
            });
        }
        relative_normal_speeds_after_mps.push(after);
        let p = Dd::from_f64(proposal.law.normal_impulses[i]);
        prescribed_work = prescribed_work + p * Dd::from_f64(contact.prescribed_surface_speed_mps);
        relative_work = relative_work
            + p * (Dd::from_f64(0.5 * proposal.relative_speeds_before_mps[i])
                + Dd::from_f64(0.5 * after));
    }
    let mut maximum_complementarity_speed_residual_mps = 0.0_f64;
    let contact_allowance = finite(
        controls.maximum_contact_speed_error_mps + controls.contact.solver.gradient_tolerance,
    )?;
    for &i in &proposal.law.active_indices {
        poll(cx)?;
        let gradient =
            finite(relative_normal_speeds_after_mps[i] - proposal.law.restitution_targets[i])?;
        let residual = if proposal.law.normal_impulses[i] > 0.0 {
            gradient.abs()
        } else {
            (-gradient).max(0.0)
        };
        maximum_complementarity_speed_residual_mps =
            maximum_complementarity_speed_residual_mps.max(residual);
        if residual > contact_allowance {
            return Err(SolidNormalImpactError::Speed {
                contact: i,
                defect_mps: residual,
                allowance_mps: contact_allowance,
            });
        }
    }
    let prescribed_surface_work_j = rounded(prescribed_work)?;
    let relative_work_j = rounded(relative_work)?;
    let relative_work_error_j = rounded(relative_work - Dd::from_f64(proposal.law.relative_work))?;
    let work_balance_error_j =
        rounded(Dd::from_f64(applied_work_j) - prescribed_work - relative_work)?;
    for defect_j in [relative_work_error_j, work_balance_error_j] {
        if defect_j.abs() > controls.maximum_work_error_j {
            return Err(SolidNormalImpactError::Work {
                defect_j,
                allowance_j: controls.maximum_work_error_j,
            });
        }
    }
    if relative_work_j > controls.contact.energy_tolerance {
        return Err(SolidNormalImpactError::Work {
            defect_j: relative_work_j,
            allowance_j: controls.contact.energy_tolerance,
        });
    }
    poll(cx)?;
    Ok(SolidContactImpulseReport {
        law: proposal.law,
        relative_normal_speeds_after_mps,
        maximum_contact_speed_error_mps,
        maximum_complementarity_speed_residual_mps,
        prescribed_surface_work_j,
        relative_work_j,
        dissipated_energy_j: (-relative_work_j).max(0.0),
        passivity_residual_j: relative_work_j.max(0.0),
        relative_work_error_j,
        work_balance_error_j,
    })
}

impl ExactLinearSolid {
    /// Resolve one declared simultaneous frictionless event against prescribed
    /// surfaces. The prepared consistent inverse supplies all cross-contact
    /// mobility; no sequential kicks, diagonal mass or propagator rebuild is
    /// introduced. Covectors follow `free_dofs()` and fixed supports remain held.
    ///
    /// The shared contact law owns restitution/complementarity. Actual native
    /// momentum, normal speeds and work must all pass before publication.
    /// Clock, positions and material/source preparation stay fixed. Geometry,
    /// event timing, friction, sustained contact and loss-to-heat partition are
    /// separate owners. Contact loss is not added again to native dissipation.
    pub fn apply_normal_impacts(
        &mut self,
        contacts: &[SolidNormalContact<'_>],
        controls: SolidNormalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidNormalImpactReport> {
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
        Ok(SolidNormalImpactReport {
            contact,
            impulse: self.accept_momentum_impulse(candidate),
        })
    }
}

impl CompleteSpectralSolid {
    /// The same simultaneous normal law on every admitted elastic and rigid
    /// mode. Complete-basis mobility is borrowed; no basis truncation, new
    /// factorization or prepared-operator clone is introduced. The original M
    /// and measured velocity still gate the accepted event, so basis roundoff
    /// cannot silently change the contact response. See the dense owner's
    /// `apply_normal_impacts` for geometry, work and no-claim boundaries.
    pub fn apply_normal_impacts(
        &mut self,
        contacts: &[SolidNormalContact<'_>],
        controls: SolidNormalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidNormalImpactReport> {
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
        Ok(SolidNormalImpactReport {
            contact,
            impulse: self.accept_momentum_impulse(candidate),
        })
    }
}
