//! Reversible vector shear-stage motion through the original finite apparatus.
//!
//! The normal stage stays fixed while the shear stage follows one straight
//! segment. The actual tip, contact area and both spring stores re-equilibrate.
//! An accepted segment remains on the connected full-stick branch throughout;
//! crossing its first slip/fold or an applicability bound refuses atomically.
//! This is quasistatic mechanical loading, with no parameter change, heat,
//! reattachment, or invented post-instability dynamics.

/// First-event continuation with retained boundary storage and unused travel.
pub mod events;

use super::{
    JkrMixedModeError, JkrPullOffBranch, JkrPullOffError, JkrPullOffFixture, JkrPullOffResponse,
    JkrPullOffState, Parameters, Work, finite,
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;

/// Immutable successor after an admitted complete shear-stage segment.
/// Accept `response().state()` together with `source()`; neither trial mutates
/// the original owner or its accepted accumulated energy release.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrShearStageTransition {
    original_fixture_identity: ContentHash,
    source: JkrPullOffFixture,
    before: JkrPullOffResponse,
    response: JkrPullOffResponse,
    mechanical_work_j: f64,
    evaluations: usize,
    event: Option<events::JkrShearStageEvent>,
    requested_stage_m: [f64; 2],
    untraversed_stage_m: [f64; 2],
}
impl JkrShearStageTransition {
    pub const fn original_fixture_identity(&self) -> ContentHash {
        self.original_fixture_identity
    }
    pub const fn source(&self) -> &JkrPullOffFixture {
        &self.source
    }
    pub fn into_source(self) -> JkrPullOffFixture {
        self.source
    }
    pub const fn before(&self) -> &JkrPullOffResponse {
        &self.before
    }
    pub const fn response(&self) -> &JkrPullOffResponse {
        &self.response
    }
    /// Integral T dot dD_t, equal to the complete potential change on this
    /// reversible segment. The held normal stage does no external work.
    pub const fn mechanical_work_j(&self) -> f64 {
        self.mechanical_work_j
    }
    /// All original-state, successor-construction, path and endpoint probes
    /// share the original fixture's single maximum_evaluations budget.
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    /// First physical boundary located by the event-continuation entry point.
    pub const fn event(&self) -> Option<events::JkrShearStageEvent> {
        self.event
    }
    pub const fn requested_stage_m(&self) -> [f64; 2] {
        self.requested_stage_m
    }
    /// Motion left for a separately selected slip/dynamic owner, in metres.
    pub const fn untraversed_stage_m(&self) -> [f64; 2] {
        self.untraversed_stage_m
    }
}

impl JkrPullOffFixture {
    /// Move the two shear stages along a straight segment at held normal stage.
    ///
    /// The original source claims, geometry, spring stiffnesses and bounds stay
    /// fixed. Stage direction may turn or reverse; actual elastic shifts are
    /// solved instead of resetting contact history to the prepared preload.
    /// Detached states stay detached and their previous release is retained.
    /// No-op motion preserves the complete fixture and canonical response.
    ///
    /// Continuous branch check: with c=2E*/Kn, eliminating the normal tip gives
    /// delta(a)=(Dn+c*a^3/(3R))/(1+c*a) and opening g=a^2/R-delta>0.
    /// On that tensile domain g''>0, so f=2*pi*W*a-E* g^2 is strictly concave.
    /// Required |D_t|^2=2*(1+k*a/Kt)^2*f/k is strictly log-concave wherever
    /// positive: it has only one maximum/fold. The stable branch has decreasing
    /// radius as |D_t| increases. The slip ratio and both shear screens then
    /// increase; delta'(a)=c*g/(1+c*a)>0 and all original radius screens are
    /// monotone. Thus admissibility is one interval of shear magnitudes.
    ///
    /// The endpoints check its upper end. The segment's minimum magnitude
    /// checks contact expansion against the original maximum radius, including
    /// reversals whose endpoints alone would miss excessive interior contact.
    /// The lower magnitude is conservatively rounded; the existing numerical
    /// equilibrium/fold brackets remain floating-point estimates, not interval
    /// certificates. No stage scanning or unbounded adaptive retry is used.
    pub fn advance_tangential_stage(
        &self,
        previous: JkrPullOffState,
        requested_stage_m: [f64; 2],
        gate: &CancelGate,
    ) -> Result<JkrShearStageTransition, JkrPullOffError> {
        if previous.identity != self.identity {
            return Err(JkrPullOffError::StateMismatch);
        }
        for value in requested_stage_m {
            finite(value)?;
        }
        let mut budget = Work::new(self.limits.maximum_evaluations, gate);
        let before = self.shear_stage_before(previous, &mut budget)?;
        self.move_shear_from_before(before, requested_stage_m, &mut budget)
    }

    fn shear_stage_before(
        &self,
        previous: JkrPullOffState,
        budget: &mut Work<'_>,
    ) -> Result<JkrPullOffResponse, JkrPullOffError> {
        budget.poll()?;
        let mut before = match previous.branch {
            JkrPullOffBranch::Attached => self.attached(previous.normal_stage_m, budget)?,
            JkrPullOffBranch::Detached => {
                self.detached(previous.normal_stage_m, previous.released_energy_j)?
            }
        };
        before.state.released_energy_j = previous.released_energy_j;
        Ok(before)
    }

    fn move_shear_from_before(
        &self,
        before: JkrPullOffResponse,
        requested_stage_m: [f64; 2],
        budget: &mut Work<'_>,
    ) -> Result<JkrShearStageTransition, JkrPullOffError> {
        let previous = before.state;
        let (source, mut response) = if requested_stage_m == self.tangential_stage {
            (self.clone(), before)
        } else {
            if previous.branch == JkrPullOffBranch::Attached
                && let Some(minimum) =
                    minimum_segment_magnitude(self.tangential_stage, requested_stage_m)
            {
                let parameters = Parameters {
                    stage: [minimum, 0.0],
                    ..self.parameters()
                };
                let maximum =
                    parameters.kinematics(self.limits.maximum_contact_radius_m, budget)?;
                if maximum.ratio <= 1.0 || previous.normal_stage_m > maximum.stage {
                    return Err(JkrPullOffError::OutsideStageEnvelope);
                }
                self.contact
                    .normal_source()
                    .validate_indentation(maximum.delta)
                    .map_err(JkrMixedModeError::from)?;
            }
            let source = Self::new_with_budget(
                self.contact.clone(),
                self.shear.clone(),
                self.normal_stiffness,
                self.tangential_stiffness,
                requested_stage_m,
                self.limits,
                budget,
            )?;
            let response = match previous.branch {
                JkrPullOffBranch::Attached => source.attached(previous.normal_stage_m, budget)?,
                JkrPullOffBranch::Detached => {
                    source.detached(previous.normal_stage_m, previous.released_energy_j)?
                }
            };
            (source, response)
        };
        response.state.released_energy_j = previous.released_energy_j;
        let mechanical_work_j = response.total_stored_energy_j - before.total_stored_energy_j;
        finite(mechanical_work_j)?;
        budget.poll()?;
        Ok(JkrShearStageTransition {
            original_fixture_identity: self.identity,
            source,
            before,
            response,
            mechanical_work_j,
            evaluations: budget.evaluations,
            event: None,
            requested_stage_m,
            untraversed_stage_m: [0.0; 2],
        })
    }
}

/// An endpoint minimum is already checked by its equilibrium. Otherwise give
/// a lower bound on the interior minimum, in metres. Scaling avoids squared
/// overflow and subtracting opposite large stage coordinates. The fixed margin
/// covers normalization, the 2-D projection and norm rounding; it deliberately
/// refuses paths too close to the radius cap to establish their interior side.
fn minimum_segment_magnitude(from: [f64; 2], to: [f64; 2]) -> Option<f64> {
    let scale = from.into_iter().chain(to).map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 {
        return None;
    }
    let p = from.map(|x| x / scale);
    let q = to.map(|x| x / scale);
    let d = [q[0] - p[0], q[1] - p[1]];
    let dot_p = p[0] * d[0] + p[1] * d[1];
    let dot_q = q[0] * d[0] + q[1] * d[1];
    if dot_p > 32.0 * f64::EPSILON || dot_q < -32.0 * f64::EPSILON {
        return None;
    }
    let squared = d[0] * d[0] + d[1] * d[1];
    let t = if squared > 0.0 {
        (-dot_p / squared).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let closest = [(1.0 - t) * p[0] + t * q[0], (1.0 - t) * p[1] + t * q[1]];
    let minimum = (closest[0].hypot(closest[1]) - 64.0 * f64::EPSILON).max(0.0) * scale;
    Some(if minimum > 0.0 {
        minimum.next_down()
    } else {
        0.0
    })
}
