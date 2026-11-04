//! Transfer the original normal trajectory's sealed jump history.
use super::{
    JkrSnapHeat, JkrSnapHeatCandidate, JkrSnapHeatControls, JkrSnapHeatError, JkrSnapHeatStamp,
};
use fs_contact::normal_patch::adhesive::motion::{JkrNormalMotionStep, MAX_JKR_NORMAL_EVENTS};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

impl<'s> JkrSnapHeat<'s> {
    /// Stage snap heat from a complete accepted normal-motion candidate.
    ///
    /// The motion owner seals its initial/final source states and every branch
    /// event. Their ordering and completeness cannot be changed through the
    /// step's public observations. This adapter regenerates original release
    /// energies; it never trusts caller-editable work or force reports.
    /// The first accepted interval binds the full motion continuation, including
    /// source, mass, policy, time and step count. Every later interval must start
    /// at that exact token. A saved closed cycle cannot be replayed with a fresh
    /// heat stamp. Once bound, prescribed-path `propose`/`apply` refuse instead
    /// of silently discarding the motion chronology.
    ///
    /// Only jump thresholds and the final approach need replay: all omitted
    /// motion stays on one reversible branch and creates no snap release. This
    /// is an equivalent heat history, not a reconstruction of smooth force work
    /// or the elapsed trajectory. In particular, a within-branch reversal needs
    /// no caller-authored endpoint. Every actual attach/detach event remains.
    ///
    /// All heat enters at the supplied thermal state's current time. Advance
    /// transport separately when it matters; this method does not invent heat
    /// evolution between event times or update the frozen normal source.
    ///
    /// Use the returned candidate with native `step_with_admission`, installing
    /// its owners only after mechanical acceptance. Any failure leaves the
    /// accepted contact/heat history and native thermal state unchanged.
    pub fn propose_normal_motion<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        motion: &JkrNormalMotionStep,
        controls: JkrSnapHeatControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrSnapHeatCandidate<'s, 'm>, JkrSnapHeatError> {
        if motion.initial_state().contact_state() != self.state()
            || self
                .normal_motion_state
                .is_some_and(|previous| previous != motion.initial_state())
        {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        let mut endpoints = [0.0; MAX_JKR_NORMAL_EVENTS + 1];
        for (index, event) in motion.events().enumerate() {
            endpoints[index] = event.before().indentation_m();
        }
        endpoints[motion.event_count()] = motion.state().contact_state().indentation_m();
        let mut candidate = self.propose_segments(
            thermal,
            expected,
            &endpoints[..=motion.event_count()],
            motion.state().contact_state(),
            controls,
            checkpoint,
        )?;
        // Bind the source regeneration to every retained original event, not
        // merely to matching branch names at the two interval endpoints.
        if candidate.report().jumps.len() != motion.event_count()
            || candidate
                .report()
                .jumps
                .iter()
                .zip(motion.events())
                .any(|(regenerated, event)| *regenerated != event.jump)
        {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        candidate.history.normal_motion_state = Some(motion.state());
        Ok(candidate)
    }
}
