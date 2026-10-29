//! Original lip/cup/valve mechanics observed as modal mouth volume flow.
//!
//! This is not a microphone: mode zero is total mouth flow, and the other
//! channels retain the bore's wall-normalized radial basis. A physical exterior
//! observer must use that basis, the declared aperture and its own load model.
use super::{BrassControl, BrassEvent, BrassSequence, BrassSequenceConfig, LipIsland,
    MmBankLoad, MmLineConfig, PressureRenderer, RenderError, invalid, physical, validate_events};
use crate::mm_line::MmLoad;
use fs_duct::{Duct, LossModel, Segment};
use fs_material::gas::GasState;
use fs_vfit::discretize::DelayedFilter;

mod realize;

/// Explicit output-realization limits, separate from the original input bank.
#[derive(Clone, Copy, Debug)]
pub struct MouthFlowConfig {
    /// Last-quarter energy divided by complete energy, for EACH mouth mode.
    pub tail_tolerance: f64,
    /// Relative change of the common uncontaminated tap prefix, for EACH mode.
    pub tap_change_tolerance: f64,
    /// Power-of-two FFT ceiling, in 512..=65536; at least two rungs are required.
    pub max_taps: usize,
    /// Sum of path x mode x ceiling taps admitted before allocation. Each
    /// runtime FIR also stores one history sample per tap. Realization scratch
    /// is bounded separately by one path's mode count and FFT ceiling.
    pub max_total_taps: usize,
}
impl Default for MouthFlowConfig {
    fn default() -> Self {
        Self { tail_tolerance: 1e-10, tap_change_tolerance: 1e-5,
            max_taps: 65_536, max_total_taps: 1_048_576 }
    }
}

/// Observed time-image convergence; not a continuum or source-validity bound.
#[derive(Clone, Debug)]
pub struct MouthFlowRung {
    /// Common tap count for this path's modal outputs.
    pub n_fft: usize,
    /// Tail-energy fractions in original modal order.
    pub tail_energy: Vec<f64>,
    /// Common-prefix changes in original modal order; infinity on first rung.
    pub tap_change: Vec<f64>,
}

/// One original valve path's output realization.
#[derive(Clone, Debug)]
pub struct MouthFlowReceipt {
    /// Original path label.
    pub label: String,
    /// Actual final staircase radius used by the modal owner [m].
    pub mouth_radius_m: f64,
    /// FIR length of every mouth-mode output for this path.
    pub n_fft: usize,
    /// Full declared output convergence history, without clipping or renormalizing.
    pub rungs: Vec<MouthFlowRung>,
}

/// A finite performance with its actual outgoing wave and all mouth modal flows.
///
/// Every output filter receives the same outgoing pressure history from sample
/// zero, including inactive paths. Valve changes select/crossfade those histories
/// at the original source sample, with the original bank's ramp. This is the
/// existing LTI path-switch approximation, NOT resolved physical valve travel.
/// All paths must have the same actual discrete aperture radius and mode basis.
///
/// The original lip, cup and input reflection remain the only mechanics. A
/// failed physical or output callback poisons continuation, since a prefix may
/// have advanced. Size refusals occur before any state or buffer is changed.
/// The source's diagnostics retain their existing callback dependence. There
/// is no allocation-free/hard-real-time or full exterior-radiation claim.
pub struct BrassMouthSequence {
    source: BrassSequence,
    filters: Vec<Vec<DelayedFilter>>,
    receipts: Vec<MouthFlowReceipt>,
    modes: usize,
    active: usize,
    previous: usize,
    fade_remaining: usize,
    fade_total: usize,
    next_event: usize,
    completed: u64,
    poisoned: bool,
}
impl BrassMouthSequence {
    /// Bind output transfers to exactly the same original sources as the voice.
    /// Analytic or strict tabulated loads are admitted; diagnostic splices and
    /// matched continuations refuse. The original input bank and this output
    /// each keep their own convergence test. No scalar mouth-flow replacement.
    ///
    /// # Errors
    /// Original voice admission, invalid source/schedule/budget, mismatched
    /// aperture bases, incomplete tables, or unconverged output time images.
    #[allow(clippy::too_many_arguments)]
    pub fn new<'a>(ducts: &[Duct], labels: &[&str], gas: &GasState,
        loads: impl Into<MmBankLoad<'a>>, line: &MmLineConfig, lip: LipIsland,
        cup_volume_m3: f64, sequence: BrassSequenceConfig, events: Vec<BrassEvent>,
        output: MouthFlowConfig) -> Result<Self, RenderError>
    {
        validate_events(sequence, line.sample_rate_hz, ducts.len(), &lip, &events)?;
        let loads = loads.into();
        realize::admit(ducts, labels, gas, loads, line, output)?;
        let source = BrassSequence::new(ducts, labels, gas, loads, line, lip,
            cup_volume_m3, sequence, events)?;
        let mut filters = Vec::with_capacity(ducts.len());
        let mut receipts: Vec<MouthFlowReceipt> = Vec::with_capacity(ducts.len());
        for (combo, duct) in ducts.iter().enumerate() {
            let (path, receipt) = realize::path(duct, labels[combo], gas,
                realize::load_at(loads, combo), line, source.voice().port_impedance(), output)?;
            if receipts.first().is_some_and(|r| r.mouth_radius_m != receipt.mouth_radius_m) {
                return Err(invalid("mouth-mode crossfades require the same discrete aperture basis"));
            }
            filters.push(path);
            receipts.push(receipt);
        }
        Ok(Self { source, filters, receipts, modes: line.n_modes,
            active: sequence.initial_combo, previous: sequence.initial_combo,
            fade_remaining: 0, fade_total: 0, next_event: 0, completed: 0, poisoned: false })
    }

    /// Original mechanics, immutable so output history cannot be bypassed.
    #[must_use]
    pub const fn source(&self) -> &BrassSequence { &self.source }
    /// Original modal order, including the plane mode at index zero.
    #[must_use]
    pub const fn mode_count(&self) -> usize { self.modes }
    /// Successfully completed source samples (not scalar output-channel values).
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.completed }
    /// Path-local convergence and aperture metadata.
    #[must_use]
    pub fn receipts(&self) -> &[MouthFlowReceipt] { &self.receipts }

    /// Render aligned mouthpiece pressure [Pa], outgoing inlet pressure [Pa],
    /// and sample-major mouth modal flow [m³/s]. For M samples, the first two
    /// buffers have M entries and `mouth_flow` has M * `mode_count()` entries.
    /// None of these streams is an exterior microphone pressure.
    ///
    /// # Errors
    /// Empty, oversized, mismatched or over-horizon buffers before mutation;
    /// original physical/output failure poisons continuation after any prefix.
    pub fn block(&mut self, mouthpiece: &mut [f64], outgoing: &mut [f64],
        mouth_flow: &mut [f64]) -> Result<(), RenderError>
    {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if mouthpiece.is_empty() { return Err(RenderError::EmptyBlock); }
        if mouthpiece.len() != outgoing.len()
            || mouthpiece.len().checked_mul(self.modes) != Some(mouth_flow.len())
            || mouthpiece.len() > self.source.max_block_len()
        {
            return Err(invalid("mouth-flow callback needs aligned buffers within source capacity"));
        }
        self.source.validate_sample_count(mouthpiece.len() as u64)?;
        // No observer can alter the physical source. Its exact event splits
        // happen first; the observer then walks that same source sample range.
        if let Err(error) = self.source.block_outgoing(mouthpiece, outgoing) {
            self.poisoned = true;
            return Err(error);
        }
        if let Err(error) = self.observe(outgoing, mouth_flow) {
            self.poisoned = true;
            return Err(error);
        }
        self.completed += mouthpiece.len() as u64;
        Ok(())
    }

    fn observe(&mut self, outgoing: &[f64], output: &mut [f64]) -> Result<(), RenderError> {
        for (offset, (&wave, frame)) in outgoing.iter().zip(output.chunks_mut(self.modes)).enumerate() {
            let sample = self.completed + offset as u64;
            while let Some(event) = self.source.events().get(self.next_event).copied() {
                if event.sample != sample { break; }
                if let BrassControl::SetValve { combo, fade_samples } = event.control {
                    if combo != self.active {
                        self.previous = self.active;
                        self.active = combo;
                        self.fade_total = fade_samples;
                        self.fade_remaining = fade_samples;
                    }
                }
                self.next_event += 1;
            }
            let mut current = [0.0; fs_duct::modal::MAX_MODES];
            let mut old = current;
            for (combo, filters) in self.filters.iter_mut().enumerate() {
                for (mode, filter) in filters.iter_mut().enumerate() {
                    let value = filter.push(wave).map_err(|e|
                        physical(format!("modal mouth-flow filter refused: {e:?}")))?;
                    if !value.is_finite() { return Err(invalid("modal mouth flow left the finite set")); }
                    if combo == self.active { current[mode] = value; }
                    if combo == self.previous { old[mode] = value; }
                }
            }
            if self.fade_remaining > 0 {
                self.fade_remaining -= 1;
                let weight = 1.0 - self.fade_remaining as f64 / self.fade_total as f64;
                for mode in 0..self.modes {
                    frame[mode] = (1.0 - weight) * old[mode] + weight * current[mode];
                }
            } else {
                frame.copy_from_slice(&current[..self.modes]);
            }
            if frame.iter().any(|x| !x.is_finite()) {
                return Err(invalid("crossfaded modal mouth flow left the finite set"));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
