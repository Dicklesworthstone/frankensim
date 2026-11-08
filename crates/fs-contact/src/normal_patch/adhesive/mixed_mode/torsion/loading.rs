//! Quasistatic prescribed twist on the original attached I/III JKR branch.
//!
//! Normal approach is held by the apparatus. The first twist boundary is
//! retained exactly: a radial fold at negative approach can release to a
//! separated state; the opening edge at positive approach cannot. No dynamic
//! impulse, heat fraction, partial slip or rebonding law is inferred.
//!
//! Radius decreases with |phi|. Endpoints and, on reversal, zero twist enclose
//! every radius on the accepted segment. The original radius screens are
//! monotone or convex on that interval, so those extrema admit the whole path.
//! Work uses the original factored potential difference, including at a fold;
//! the singular equilibrium tangent is never inverted.

use super::{
    JkrTorsionalError, JkrTorsionalLimit, JkrTorsionalLimitKind, JkrTorsionalResponse,
    SourceJkrMixedMode, path::JkrTorsionalPathError,
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalLoadingControls {
    /// Original equilibrium queries per operation; a reversal needs at most two.
    pub maximum_evaluations: usize,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrTorsionalLoadingError {
    Source(JkrTorsionalError),
    Path(JkrTorsionalPathError),
    Invalid(&'static str),
    StateMismatch,
    Budget,
    Cancelled,
    /// Opening-edge loss of applicability is not a detachment instability.
    NotAtRadialFold,
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<JkrTorsionalError> for JkrTorsionalLoadingError {
    fn from(value: JkrTorsionalError) -> Self {
        Self::Source(value)
    }
}
impl From<JkrTorsionalPathError> for JkrTorsionalLoadingError {
    fn from(value: JkrTorsionalPathError) -> Self {
        Self::Path(value)
    }
}
impl core::fmt::Display for JkrTorsionalLoadingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prescribed torsional JKR loading: {self:?}")
    }
}
impl std::error::Error for JkrTorsionalLoadingError {}

/// A source-owned attached or separated history, never a caller-edited report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalLoadingState {
    identity: ContentHash,
    indentation_m: f64,
    twist_angle_rad: f64,
    contact: Option<JkrTorsionalResponse>,
    released_energy_j: f64,
    steps: u64,
}
impl JkrTorsionalLoadingState {
    pub const fn model_identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn indentation_m(&self) -> f64 {
        self.indentation_m
    }
    pub const fn twist_angle_rad(&self) -> f64 {
        self.twist_angle_rad
    }
    pub const fn contact(&self) -> Option<JkrTorsionalResponse> {
        self.contact
    }
    /// Accumulated mechanical release, awaiting an explicit downstream destination.
    pub const fn released_energy_j(&self) -> f64 {
        self.released_energy_j
    }
    pub const fn steps(&self) -> u64 {
        self.steps
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.contact.map_or(0.0, |point| point.reversible_energy_j)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalLoadingReport {
    pub before: Option<JkrTorsionalResponse>,
    pub contact: Option<JkrTorsionalResponse>,
    pub boundary: Option<JkrTorsionalLimitKind>,
    pub requested_twist_rad: f64,
    pub accepted_twist_rad: f64,
    /// Caller travel beyond the first boundary; never silently consumed.
    pub untraversed_twist_rad: f64,
    /// Integral M dphi on the smooth attached part; held approach does no work.
    pub prescribed_work_j: f64,
    /// Mechanical potential lost at a held-coordinate snap, not assigned heat.
    pub released_energy_j: f64,
    pub storage_change_j: f64,
    /// Delta U + release - prescribed work.
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub minimum_contact_radius_m: f64,
    pub maximum_contact_radius_m: f64,
    pub evaluations: usize,
}

/// Accept the sealed successor only after all consumers accept this proposal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalLoadingStep {
    previous: JkrTorsionalLoadingState,
    state: JkrTorsionalLoadingState,
    report: JkrTorsionalLoadingReport,
}
impl JkrTorsionalLoadingStep {
    pub const fn previous_state(&self) -> JkrTorsionalLoadingState {
        self.previous
    }
    pub const fn state(&self) -> JkrTorsionalLoadingState {
        self.state
    }
    pub const fn report(&self) -> &JkrTorsionalLoadingReport {
        &self.report
    }
}

#[derive(Clone, Copy, Debug)]
pub struct JkrTorsionalLoading<'a> {
    source: &'a SourceJkrMixedMode,
    indentation_m: f64,
    boundary: Option<JkrTorsionalLimit>,
    controls: JkrTorsionalLoadingControls,
}
impl<'a> JkrTorsionalLoading<'a> {
    /// Prepare the original fixed-approach apparatus and its ideal twist limit.
    /// Source applicability may bind before that ideal limit.
    pub fn new(
        source: &'a SourceJkrMixedMode,
        indentation_m: f64,
        controls: JkrTorsionalLoadingControls,
        gate: &CancelGate,
    ) -> Result<Self, JkrTorsionalLoadingError> {
        poll(gate)?;
        if controls.maximum_evaluations == 0
            || [
                controls.energy_absolute_tolerance_j,
                controls.energy_relative_tolerance,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(JkrTorsionalLoadingError::Invalid("loading controls"));
        }
        let boundary = source.torsional_displacement_limit(indentation_m)?;
        poll(gate)?;
        Ok(Self {
            source,
            indentation_m,
            boundary,
            controls,
        })
    }
    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    pub const fn limit(&self) -> Option<JkrTorsionalLimit> {
        self.boundary
    }
    pub const fn controls(&self) -> JkrTorsionalLoadingControls {
        self.controls
    }

    /// Start from an actual equilibrium, optionally exactly at its boundary.
    pub fn initial_state(
        &self,
        twist_angle_rad: f64,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalLoadingState, JkrTorsionalLoadingError> {
        let contact = self.evaluate(twist_angle_rad, &mut 0, gate)?;
        poll(gate)?;
        Ok(JkrTorsionalLoadingState {
            identity: self.source.identity(),
            indentation_m: self.indentation_m,
            twist_angle_rad,
            contact: Some(contact),
            released_energy_j: 0.0,
            steps: 0,
        })
    }

    /// Follow a straight prescribed twist segment, stopping at its first
    /// original attached boundary. Reversals preserve and return elastic work.
    /// A detached history remains detached; twisting alone cannot rebond it.
    pub fn advance_to_twist(
        &self,
        previous: JkrTorsionalLoadingState,
        requested_twist_rad: f64,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalLoadingStep, JkrTorsionalLoadingError> {
        self.admit(previous, gate)?;
        if !requested_twist_rad.is_finite()
            || !(requested_twist_rad - previous.twist_angle_rad).is_finite()
        {
            return Err(JkrTorsionalLoadingError::Invalid("finite prescribed twist"));
        }
        let mut accepted = requested_twist_rad;
        let mut evaluations = 0;
        let mut radius = [0.0; 2];
        let mut work = 0.0;
        let contact = if let Some(before) = previous.contact {
            if let Some(limit) = self.boundary {
                accepted = accepted.clamp(-limit.twist_angle_rad, limit.twist_angle_rad);
            }
            let after = if accepted == previous.twist_angle_rad {
                before
            } else {
                self.evaluate(accepted, &mut evaluations, gate)?
            };
            radius = [
                before.contact_radius_m.min(after.contact_radius_m),
                before.contact_radius_m.max(after.contact_radius_m),
            ];
            if before.twist_angle_rad.is_sign_negative() != after.twist_angle_rad.is_sign_negative()
                && before.twist_angle_rad != 0.0
                && after.twist_angle_rad != 0.0
            {
                let zero = self.evaluate(0.0, &mut evaluations, gate)?;
                radius[1] = radius[1].max(zero.contact_radius_m);
            }
            // Nominal strain a/R + sqrt(C/a) is convex; the other
            // original radius screens are monotone. Endpoint extrema suffice.
            radius = [radius[0].next_down(), radius[1].next_up()];
            for a in radius {
                self.source
                    .normal_source()
                    .check_applicability(a)
                    .map_err(JkrTorsionalError::from)?;
            }
            work = self.source.torsional_endpoint_work(before, after)?;
            Some(after)
        } else {
            None
        };
        let state = JkrTorsionalLoadingState {
            twist_angle_rad: accepted,
            contact,
            steps: next_step(previous.steps)?,
            ..previous
        };
        self.finish(
            previous,
            state,
            requested_twist_rad,
            work,
            0.0,
            radius,
            evaluations,
            gate,
        )
    }

    /// Release only the authenticated radial fold at negative approach.
    ///
    /// The original equilibrium potential supplies the release. At this fold,
    /// U=E*a*(delta^2/3+a^4/(5R^2))>0. This follows from the original
    /// stationarity and zero radial curvature; it is not a fitted loss fraction.
    /// The apparatus holds both coordinates, so no mechanical work is done
    /// during the jump. The receiver must explicitly route released energy.
    /// Releasing twice, releasing at the opening edge, or rebonding is refused.
    pub fn detach_at_fold(
        &self,
        previous: JkrTorsionalLoadingState,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalLoadingStep, JkrTorsionalLoadingError> {
        self.admit(previous, gate)?;
        let contact = previous
            .contact
            .filter(|point| {
                point.boundary == Some(JkrTorsionalLimitKind::RadialFold)
                    && point.indentation_m < 0.0
            })
            .ok_or(JkrTorsionalLoadingError::NotAtRadialFold)?;
        let release = contact.reversible_energy_j;
        if !release.is_finite() || release < 0.0 {
            return Err(JkrTorsionalLoadingError::Invalid(
                "nonnegative fold release",
            ));
        }
        let total = (Dd::from_f64(previous.released_energy_j) + Dd::from_f64(release)).to_f64();
        if !total.is_finite() {
            return Err(JkrTorsionalLoadingError::Invalid("cumulative fold release"));
        }
        let state = JkrTorsionalLoadingState {
            contact: None,
            released_energy_j: total,
            steps: next_step(previous.steps)?,
            ..previous
        };
        self.finish(
            previous,
            state,
            previous.twist_angle_rad,
            0.0,
            release,
            [0.0, contact.contact_radius_m],
            0,
            gate,
        )
    }

    fn admit(
        &self,
        state: JkrTorsionalLoadingState,
        gate: &CancelGate,
    ) -> Result<(), JkrTorsionalLoadingError> {
        poll(gate)?;
        if state.identity != self.source.identity()
            || state.indentation_m.to_bits() != self.indentation_m.to_bits()
        {
            return Err(JkrTorsionalLoadingError::StateMismatch);
        }
        Ok(())
    }

    fn evaluate(
        &self,
        twist: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalResponse, JkrTorsionalLoadingError> {
        poll(gate)?;
        if *evaluations >= self.controls.maximum_evaluations {
            return Err(JkrTorsionalLoadingError::Budget);
        }
        *evaluations += 1;
        Ok(self
            .source
            .evaluate_torsion_displacement_controlled(self.indentation_m, twist)?)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        previous: JkrTorsionalLoadingState,
        state: JkrTorsionalLoadingState,
        requested: f64,
        work: f64,
        release: f64,
        radius: [f64; 2],
        evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalLoadingStep, JkrTorsionalLoadingError> {
        let change = (Dd::from_f64(state.stored_energy_j())
            - Dd::from_f64(previous.stored_energy_j()))
        .to_f64();
        let residual = (Dd::from_f64(change) + Dd::from_f64(release) - Dd::from_f64(work)).to_f64();
        let mut scale = work.abs() + release;
        for point in [previous.contact, state.contact].into_iter().flatten() {
            scale += point.normal_elastic_energy_j + point.torsional_elastic_energy_j
                - point.surface_energy_j;
        }
        let tolerance = self.controls.energy_absolute_tolerance_j
            + self.controls.energy_relative_tolerance * scale;
        if [change, residual, scale, tolerance]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(JkrTorsionalLoadingError::Invalid("finite work balance"));
        }
        if residual.abs() > tolerance {
            return Err(JkrTorsionalLoadingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        poll(gate)?;
        Ok(JkrTorsionalLoadingStep {
            previous,
            state,
            report: JkrTorsionalLoadingReport {
                before: previous.contact,
                contact: state.contact,
                boundary: state.contact.and_then(|point| point.boundary),
                requested_twist_rad: requested,
                accepted_twist_rad: state.twist_angle_rad,
                untraversed_twist_rad: requested - state.twist_angle_rad,
                prescribed_work_j: work,
                released_energy_j: release,
                storage_change_j: change,
                energy_residual_j: residual,
                energy_tolerance_j: tolerance,
                minimum_contact_radius_m: radius[0],
                maximum_contact_radius_m: radius[1],
                evaluations,
            },
        })
    }
}

fn next_step(steps: u64) -> Result<u64, JkrTorsionalLoadingError> {
    steps
        .checked_add(1)
        .ok_or(JkrTorsionalLoadingError::Invalid("loading step overflow"))
}
fn poll(gate: &CancelGate) -> Result<(), JkrTorsionalLoadingError> {
    if gate.is_requested() {
        Err(JkrTorsionalLoadingError::Cancelled)
    } else {
        Ok(())
    }
}
