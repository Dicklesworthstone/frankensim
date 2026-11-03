//! One simultaneous Coulomb event across independent complete solid owners.
//!
//! Point maps supply reciprocal impulses to finite participants. The existing
//! frictional-impact owner assembles the full contact Delassus matrix and solves
//! the original non-associated Coulomb law. Each original trajectory supplies
//! its own mass action and unpublished impulse trial; no prepared matrix is
//! copied or refactorized. All actual endpoints and work accounts pass before
//! any trajectory, modal state, reaction history or source account is accepted.

use super::frictional_impact::{
    self, SolidFrictionalContact, SolidFrictionalContactImpulseReport,
    SolidFrictionalImpactControls, SolidFrictionalImpactError,
};
use super::impulse::{self, ExactImpulseCandidate};
use super::spectral::{CompleteSpectralSolid, SpectralImpulseCandidate};
use super::supported::{
    SupportedImpulseCandidate, SupportedSolidError, SupportedSolidImpulse, SupportedSpectralSolid,
};
use super::{ExactLinearSolid, LinearSolidImpulseError, LinearSolidImpulseReport};
use core::ops::Range;
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::Cx;
use fs_math::dd::Dd;

/// Original accepted trajectories, exclusively borrowed for this one event.
/// Exact and spectral point maps use `free_dofs()` order; supported maps use
/// all `represented_dofs()`, including actual prescribed support velocities.
pub enum SolidContactBody<'a> {
    Exact(&'a mut ExactLinearSolid),
    Spectral(&'a mut CompleteSpectralSolid),
    Supported(&'a mut SupportedSpectralSolid),
}

/// A geometry-owned point velocity in a shared orthonormal contact frame.
/// The three covectors are normal, tangent 1 and tangent 2, respectively.
/// Geometry owns frame orientation, interpolation and coincident event points.
#[derive(Clone, Copy, Debug)]
pub struct SolidContactPoint<'a> {
    pub body: usize,
    pub velocity_weights: [&'a [f64]; 3],
}

/// Either another finite participant or an external velocity source.
#[derive(Clone, Copy, Debug)]
pub enum SolidContactOther<'a> {
    Body(SolidContactPoint<'a>),
    /// Opposing point velocity [m/s] in the first point's contact frame.
    PrescribedVelocity([f64; 3]),
}

/// One declared pair. Relative velocity is first point minus second point.
/// A positive normal impulse acts on the first point and its negative on the
/// second finite point. Both maps must express the SAME normal/tangent frame.
#[derive(Clone, Copy, Debug)]
pub struct SolidContactPair<'a> {
    pub id: u64,
    pub first: SolidContactPoint<'a>,
    pub second: SolidContactOther<'a>,
    pub gap_m: f64,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// Network packing and candidate limits in addition to the shared law limits.
#[derive(Clone, Copy, Debug)]
pub struct SolidContactNetworkControls {
    pub impact: SolidFrictionalImpactControls,
    pub maximum_bodies: usize,
    pub maximum_total_dofs: usize,
    /// Additional conservative scratch for packed point maps, states and body
    /// trials. Shared impact/solver scratch retains its separate budget.
    pub maximum_workspace_bytes: usize,
}

/// Original per-body receipt, in caller body order. Supported reaction work
/// remains distinguishable from the applied finite-pair contact work.
#[derive(Clone, Debug)]
pub enum SolidContactBodyImpulse {
    Free(LinearSolidImpulseReport),
    Supported(SupportedSolidImpulse),
}

/// Accepted instantaneous network event at the common original clock.
#[derive(Clone, Debug)]
pub struct SolidContactNetworkReport {
    pub time_s: f64,
    pub bodies: Vec<SolidContactBodyImpulse>,
    /// Retains native determinacy and total-work witnesses. A velocity-only
    /// law result does not establish unique per-contact loads or heat shares.
    pub contact: SolidFrictionalContactImpulseReport,
    pub applied_work_j: f64,
    pub support_work_j: f64,
    pub kinetic_change_j: f64,
    /// Applied plus support work minus actual native kinetic increments [J].
    pub work_balance_error_j: f64,
}

/// Every refusal leaves every borrowed owner at its previously accepted state.
#[derive(Debug)]
pub enum SolidContactNetworkError {
    Invalid(&'static str),
    Budget {
        resource: &'static str,
        required: u128,
        maximum: u128,
    },
    UnsynchronizedTime {
        body: usize,
        expected_s: f64,
        actual_s: f64,
    },
    Contact(SolidFrictionalImpactError),
    Native {
        body: usize,
        error: LinearSolidImpulseError,
    },
    Supported {
        body: usize,
        error: SupportedSolidError,
    },
    Work {
        defect_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for SolidContactNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete solid contact network refused: {self:?}")
    }
}
impl core::error::Error for SolidContactNetworkError {}

type Result<T> = core::result::Result<T, SolidContactNetworkError>;

fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| SolidContactNetworkError::Cancelled)
}

fn reserve<T>(count: usize) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| SolidContactNetworkError::Invalid("contact network allocation refused"))?;
    Ok(values)
}

fn zeros(count: usize) -> Result<Vec<f64>> {
    let mut values = reserve(count)?;
    values.resize(count, 0.0);
    Ok(values)
}

fn bound(resource: &'static str, required: u128, maximum: u128) -> Result<()> {
    if required > maximum {
        Err(SolidContactNetworkError::Budget {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}

fn represented(value: Dd) -> Result<f64> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(SolidContactNetworkError::Invalid(
            "unrepresentable network work",
        ))
    }
}

fn contact_error(error: SolidFrictionalImpactError) -> SolidContactNetworkError {
    if matches!(error, SolidFrictionalImpactError::Cancelled) {
        SolidContactNetworkError::Cancelled
    } else {
        SolidContactNetworkError::Contact(error)
    }
}

impl SolidContactBody<'_> {
    fn velocity(&self) -> &[f64] {
        match self {
            Self::Exact(body) => body.velocity_mps(),
            Self::Spectral(body) => body.velocity_mps(),
            Self::Supported(body) => body.velocity_mps(),
        }
    }

    fn time(&self) -> f64 {
        match self {
            Self::Exact(body) => body.time_s(),
            Self::Spectral(body) => body.time_s(),
            Self::Supported(body) => body.time_s(),
        }
    }

    fn response(
        &self,
        impulse: &[f64],
        cx: &Cx<'_>,
    ) -> core::result::Result<Vec<f64>, LinearSolidImpulseError> {
        match self {
            Self::Exact(body) => impulse::action(&body.inverse_mass, impulse, cx),
            Self::Spectral(body) => body.momentum_response(impulse, cx),
            Self::Supported(body) => body.momentum_response(impulse, cx),
        }
    }

    fn trial<'b>(
        &'b mut self,
        index: usize,
        impulse: &[f64],
        allowance: f64,
        cx: &Cx<'_>,
    ) -> Result<Candidate<'b>> {
        let native_error = |error| {
            if matches!(error, LinearSolidImpulseError::Cancelled) {
                SolidContactNetworkError::Cancelled
            } else {
                SolidContactNetworkError::Native { body: index, error }
            }
        };
        match self {
            Self::Exact(body) => {
                let trial = body
                    .trial_momentum_impulse(impulse, allowance, cx)
                    .map_err(native_error)?;
                Ok(Candidate::Exact { body, trial })
            }
            Self::Spectral(body) => {
                let trial = body
                    .trial_momentum_impulse(impulse, allowance, cx)
                    .map_err(native_error)?;
                Ok(Candidate::Spectral { body, trial })
            }
            Self::Supported(body) => {
                let trial = body
                    .trial_fixed_support_momentum_impulse(impulse, allowance, cx)
                    .map_err(|error| {
                        if matches!(error, SupportedSolidError::Cancelled) {
                            SolidContactNetworkError::Cancelled
                        } else {
                            SolidContactNetworkError::Supported { body: index, error }
                        }
                    })?;
                Ok(Candidate::Supported { body, trial })
            }
        }
    }
}

// Keeping each native owner and its own trial in one variant makes the final
// publication infallible; an index cannot accidentally apply another body's
// candidate or a different representation's modal state.
enum Candidate<'a> {
    Exact {
        body: &'a mut ExactLinearSolid,
        trial: ExactImpulseCandidate,
    },
    Spectral {
        body: &'a mut CompleteSpectralSolid,
        trial: SpectralImpulseCandidate,
    },
    Supported {
        body: &'a mut SupportedSpectralSolid,
        trial: SupportedImpulseCandidate,
    },
}

impl Candidate<'_> {
    fn velocity(&self) -> &[f64] {
        match self {
            Self::Exact { trial, .. } => &trial.velocity_mps,
            Self::Spectral { trial, .. } => &trial.velocity_mps,
            Self::Supported { trial, .. } => &trial.velocity_mps,
        }
    }

    fn work(&self) -> [f64; 3] {
        match self {
            Self::Exact { trial, .. } => [
                trial.report.energy.supplied_work_j,
                0.0,
                trial.report.kinetic_change_j,
            ],
            Self::Spectral { trial, .. } => [
                trial.report.energy.supplied_work_j,
                0.0,
                trial.report.kinetic_change_j,
            ],
            Self::Supported { trial, .. } => [
                trial.energy.applied_work_j,
                trial.energy.support_work_j,
                trial.kinetic_change_j,
            ],
        }
    }

    fn accept(self) -> SolidContactBodyImpulse {
        match self {
            Self::Exact { body, trial } => {
                SolidContactBodyImpulse::Free(body.accept_momentum_impulse(trial))
            }
            Self::Spectral { body, trial } => {
                SolidContactBodyImpulse::Free(body.accept_momentum_impulse(trial))
            }
            Self::Supported { body, trial } => {
                SolidContactBodyImpulse::Supported(body.accept_momentum_impulse(trial))
            }
        }
    }
}

fn add_point(
    rows: &mut [Vec<f64>; 3],
    point: SolidContactPoint<'_>,
    ranges: &[Range<usize>],
    sign: f64,
    cx: &Cx<'_>,
) -> Result<()> {
    let range = ranges
        .get(point.body)
        .ok_or(SolidContactNetworkError::Invalid(
            "contact point body index out of range",
        ))?;
    for axis in 0..3 {
        if point.velocity_weights[axis].len() != range.len() {
            return Err(SolidContactNetworkError::Invalid(
                "contact point requires every original native coordinate",
            ));
        }
        for (local, &weight) in point.velocity_weights[axis].iter().enumerate() {
            if local % 256 == 0 {
                poll(cx)?;
            }
            let target = &mut rows[axis][range.start + local];
            *target += sign * weight;
            if !weight.is_finite() || !target.is_finite() {
                return Err(SolidContactNetworkError::Invalid(
                    "finite representable contact point maps required",
                ));
            }
        }
    }
    Ok(())
}

/// Resolve 1..=21 declared pairs at the identical accepted time of every body.
/// Finite pairs can share bodies or connect two points of the same body; their
/// cross-contact and normal/shear coupling remains in one original law solve.
/// Each mass action visits only nonzero participating blocks. Complete rigid
/// and elastic modes remain active, and prescribed support rates remain fixed.
///
/// This supplies an instantaneous event, not a detector or event-time finder.
/// Caller geometry must provide all simultaneous contacts and physically
/// consistent point maps. For arbitrary algebraic maps, no world linear or
/// angular momentum claim can be inferred beyond the original DOF balances.
/// No impact loss is silently assigned to heat or native viscous dissipation.
pub fn apply_frictional_impacts(
    bodies: &mut [SolidContactBody<'_>],
    contacts: &[SolidContactPair<'_>],
    controls: SolidContactNetworkControls,
    cx: &Cx<'_>,
) -> Result<SolidContactNetworkReport> {
    poll(cx)?;
    if bodies.is_empty() || !(1..=21).contains(&contacts.len()) {
        return Err(SolidContactNetworkError::Invalid(
            "nonempty native body network and 1..=21 contact pairs required",
        ));
    }
    bound(
        "network bodies",
        bodies.len() as u128,
        controls.maximum_bodies as u128,
    )?;
    let mut n = 0usize;
    let time_s = bodies[0].time();
    for (index, body) in bodies.iter().enumerate() {
        poll(cx)?;
        if !body.time().is_finite() || body.time() != time_s {
            return Err(SolidContactNetworkError::UnsynchronizedTime {
                body: index,
                expected_s: time_s,
                actual_s: body.time(),
            });
        }
        n = n
            .checked_add(body.velocity().len())
            .ok_or(SolidContactNetworkError::Invalid(
                "network dimension overflow",
            ))?;
    }
    bound(
        "network native DOFs",
        n as u128,
        controls.maximum_total_dofs as u128,
    )?;
    let (nn, mm, bb) = (n as u128, contacts.len() as u128, bodies.len() as u128);
    bound(
        "network workspace bytes",
        256 * nn + 24 * mm * nn + 2048 * bb + 1024 * mm,
        controls.maximum_workspace_bytes as u128,
    )?;

    let mut ranges = reserve(bodies.len())?;
    let mut before = reserve(n)?;
    for body in bodies.iter() {
        let start = before.len();
        for block in body.velocity().chunks(256) {
            poll(cx)?;
            before.extend_from_slice(block);
        }
        ranges.push(start..before.len());
    }
    let mut point_rows = reserve(contacts.len())?;
    for contact in contacts {
        poll(cx)?;
        let mut rows = [zeros(n)?, zeros(n)?, zeros(n)?];
        add_point(&mut rows, contact.first, &ranges, 1.0, cx)?;
        if let SolidContactOther::Body(second) = contact.second {
            add_point(&mut rows, second, &ranges, -1.0, cx)?;
        }
        point_rows.push(rows);
    }
    let mut native_contacts = reserve(contacts.len())?;
    for (contact, rows) in contacts.iter().zip(&point_rows) {
        native_contacts.push(SolidFrictionalContact {
            id: contact.id,
            velocity_weights: [&rows[0], &rows[1], &rows[2]],
            gap_m: contact.gap_m,
            activation_tolerance_m: contact.activation_tolerance_m,
            prescribed_surface_velocity_mps: match contact.second {
                SolidContactOther::Body(_) => [0.0; 3],
                SolidContactOther::PrescribedVelocity(velocity) => velocity,
            },
            friction_coefficient: contact.friction_coefficient,
            normal_law: contact.normal_law,
        });
    }
    let proposal = frictional_impact::propose(
        &before,
        &native_contacts,
        controls.impact,
        |impulse, cx| {
            let mut response = Vec::new();
            response.try_reserve_exact(n).map_err(|_| {
                LinearSolidImpulseError::Invalid("contact network response allocation refused")
            })?;
            response.resize(n, 0.0);
            for (body, range) in bodies.iter().zip(&ranges) {
                impulse::poll(cx)?;
                let input = &impulse[range.clone()];
                let mut nonzero = false;
                for block in input.chunks(256) {
                    impulse::poll(cx)?;
                    if block.iter().any(|value| *value != 0.0) {
                        nonzero = true;
                        break;
                    }
                }
                if !nonzero {
                    continue;
                }
                let block = body.response(input, cx)?;
                if block.len() != range.len() {
                    return Err(LinearSolidImpulseError::Invalid(
                        "native response shape changed",
                    ));
                }
                response[range.clone()].copy_from_slice(&block);
            }
            Ok(response)
        },
        cx,
    )
    .map_err(contact_error)?;

    // Allocate all publication storage before borrowing owners into candidates.
    let mut receipts = reserve(bodies.len())?;
    let mut candidates = reserve(bodies.len())?;
    let mut actual_after = reserve(n)?;
    let mut work = [Dd::ZERO; 3];
    for (index, (body, range)) in bodies.iter_mut().zip(&ranges).enumerate() {
        poll(cx)?;
        let trial = body.trial(
            index,
            &proposal.impulse_kg_m_s[range.clone()],
            controls.impact.maximum_impulse_error_kg_m_s,
            cx,
        )?;
        for block in trial.velocity().chunks(256) {
            poll(cx)?;
            actual_after.extend_from_slice(block);
        }
        for (sum, value) in work.iter_mut().zip(trial.work()) {
            *sum = *sum + Dd::from_f64(value);
        }
        candidates.push(trial);
    }
    let [applied_work_j, support_work_j, kinetic_change_j] = [
        represented(work[0])?,
        represented(work[1])?,
        represented(work[2])?,
    ];
    let work_balance_error_j = represented(work[0] + work[1] - work[2])?;
    if work_balance_error_j.abs() > controls.impact.maximum_work_error_j {
        return Err(SolidContactNetworkError::Work {
            defect_j: work_balance_error_j,
            allowance_j: controls.impact.maximum_work_error_j,
        });
    }
    let contact = frictional_impact::audit(
        proposal,
        &native_contacts,
        &actual_after,
        applied_work_j,
        controls.impact,
        cx,
    )
    .map_err(contact_error)?;
    poll(cx)?;
    // No fallible work, allocation or cancellation point after this boundary.
    for candidate in candidates {
        receipts.push(candidate.accept());
    }
    Ok(SolidContactNetworkReport {
        time_s,
        bodies: receipts,
        contact,
        applied_work_j,
        support_work_j,
        kinetic_change_j,
        work_balance_error_j,
    })
}
