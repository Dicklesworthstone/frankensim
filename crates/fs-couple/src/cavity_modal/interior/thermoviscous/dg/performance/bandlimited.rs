//! One shared filter for the original selected port pressure and terminal flows.
//!
//! This adapter is independent of resolved NSF geometry and controls. It never
//! filters mechanical input or changes an accepted field interval. All source
//! channels are filtered before the existing physical Rayleigh receiver, whose
//! geometry, material propagation, differentiation and delay remain its own.
//! The finite window consumes exactly the original covering source intervals;
//! filter startup is zero history and there is no future-data or tail padding.
//! This is observation interpolation, not energy-conserving coupling across
//! different mechanical clocks or a certificate of the source's bandwidth.

use super::{
    AcousticRealizeError, BandlimitedObservationPolicy, BandlimitedResampler, CancelGate,
    RenderError, sizing,
};

pub(super) struct BandlimitedDgObservation {
    resampler: BandlimitedResampler,
    frame: Vec<f64>,
}

impl BandlimitedDgObservation {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        physical_dt_s: f64,
        output_rate_hz: u32,
        output_samples: u64,
        source_intervals: usize,
        terminal_channels: usize,
        policy: BandlimitedObservationPolicy,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        let channels = terminal_channels
            .checked_add(1)
            .ok_or_else(|| sizing("resolved field observation channel count overflows"))?;
        let resampler = BandlimitedResampler::new_cancellable(
            physical_dt_s,
            output_rate_hz,
            channels,
            output_samples,
            policy,
            &mut || gate.is_requested(),
        )
        .map_err(|what| {
            if gate.is_requested() {
                RenderError::Voice(AcousticRealizeError::Cancelled)
            } else {
                sizing(what)
            }
        })?;
        if source_intervals != resampler.source_intervals() {
            return Err(sizing(
                "bandlimited field controls must cover exactly the original physical output window",
            ));
        }
        let mut frame = Vec::new();
        frame
            .try_reserve_exact(channels)
            .map_err(|_| sizing("resolved field observation frame allocation failed"))?;
        frame.resize(channels, 0.0);
        Ok(Self { resampler, frame })
    }

    pub(super) fn resampler(&self) -> &BandlimitedResampler {
        &self.resampler
    }

    pub(super) fn required_source_intervals(&self) -> Result<usize, RenderError> {
        self.resampler.required_source_intervals().map_err(sizing)
    }

    pub(super) fn push_source(
        &mut self,
        pressure: f64,
        terminal_flows: &[f64],
    ) -> Result<(), RenderError> {
        if terminal_flows.len() != self.frame.len() - 1 {
            return Err(sizing(
                "bandlimited field observation lost its original terminal channels",
            ));
        }
        self.frame[0] = pressure;
        self.frame[1..].copy_from_slice(terminal_flows);
        self.resampler.push_frame(&self.frame).map_err(sizing)
    }

    pub(super) fn render_next(&mut self, terminal_flows: &mut [f64]) -> Result<f64, RenderError> {
        if terminal_flows.len() != self.frame.len() - 1 {
            return Err(sizing(
                "bandlimited field output lost its original terminal channels",
            ));
        }
        self.resampler
            .render_next(&mut self.frame)
            .map_err(sizing)?;
        terminal_flows.copy_from_slice(&self.frame[1..]);
        Ok(self.frame[0])
    }
}

#[cfg(test)]
#[path = "bandlimited/tests.rs"]
mod tests;
