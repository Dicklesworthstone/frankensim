//! One common temporal filter for every original source channel, before the
//! Rayleigh receivers differentiate, delay and sum their distinct surfaces.

use super::{RenderError, sizing, source_intervals_for};
use crate::pcm_wav::observation::resample::{BandlimitedObservationPolicy, BandlimitedResampler};

/// Cold channel order: inlet pressure, original terminal flows, then each
/// selected reservoir's original tile flows. Aperture/tile identities remain
/// owned and checked by their existing physical observers.
pub(in crate::nsf_dg) struct BandlimitedObservation {
    resampler: BandlimitedResampler,
    physical_dt_s: f64,
    output_rate_hz: u32,
    output_samples: u64,
    terminal_channels: usize,
    reservoir_channels: Vec<usize>,
    source: Vec<f64>,
}

impl BandlimitedObservation {
    pub(in crate::nsf_dg) fn new(
        physical_dt_s: f64,
        output_rate_hz: u32,
        output_samples: u64,
        terminal_channels: usize,
        reservoir_channels: impl IntoIterator<Item = usize>,
        policy: BandlimitedObservationPolicy,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        let reservoir_channels: Vec<_> = reservoir_channels.into_iter().collect();
        let channels = reservoir_channels.iter().try_fold(
            terminal_channels
                .checked_add(1)
                .ok_or_else(|| sizing("resolved observation source channel count overflow"))?,
            |count, &next| {
                count
                    .checked_add(next)
                    .ok_or_else(|| sizing("resolved observation source channel count overflow"))
            },
        )?;
        let resampler = BandlimitedResampler::new_cancellable(
            physical_dt_s,
            output_rate_hz,
            channels,
            output_samples,
            policy,
            cancelled,
        )
        .map_err(sizing)?;
        if resampler.source_intervals()
            != source_intervals_for(physical_dt_s, output_rate_hz, output_samples)?
        {
            return Err(sizing(
                "bandlimited observation must retain the original complete physical source window",
            ));
        }
        let mut source = Vec::new();
        source
            .try_reserve_exact(channels)
            .map_err(|_| sizing("bandlimited observation channel scratch allocation failed"))?;
        source.resize(channels, 0.0);
        Ok(Self {
            resampler,
            physical_dt_s,
            output_rate_hz,
            output_samples,
            terminal_channels,
            reservoir_channels,
            source,
        })
    }

    pub(in crate::nsf_dg) fn resampler(&self) -> &BandlimitedResampler {
        &self.resampler
    }

    /// Transfer the actual unadvanced filter, without rebuilding its kernel or
    /// changing its source clock, channel basis, finite window or zero history.
    pub(in crate::nsf_dg) fn validate_cold_basis(
        &self,
        physical_dt_s: f64,
        output_rate_hz: u32,
        output_samples: u64,
        terminal_channels: usize,
        reservoir_channels: impl IntoIterator<Item = usize>,
    ) -> Result<(), RenderError> {
        if self.physical_dt_s.to_bits() != physical_dt_s.to_bits()
            || self.output_rate_hz != output_rate_hz
            || self.output_samples != output_samples
            || self.terminal_channels != terminal_channels
            || !self
                .reservoir_channels
                .iter()
                .copied()
                .eq(reservoir_channels)
            || self.resampler.accepted_source_intervals() != 0
            || self.resampler.rendered_samples() != 0
        {
            return Err(sizing(
                "cold bandlimited observation requires identical original clocks, channels and history",
            ));
        }
        Ok(())
    }

    /// Physical coverage of the next output interval, independent of FIR delay.
    /// The mechanics still consume every authored interval through that end.
    pub(in crate::nsf_dg) fn required_source_intervals(&self) -> Result<usize, RenderError> {
        self.resampler.required_source_intervals().map_err(sizing)
    }

    pub(in crate::nsf_dg) fn push_source<'s>(
        &mut self,
        pressure: f64,
        terminals: &[f64],
        reservoirs: impl IntoIterator<Item = &'s [f64]>,
    ) -> Result<(), RenderError> {
        if terminals.len() != self.terminal_channels {
            return Err(sizing(
                "bandlimited observation lost its original terminal channels",
            ));
        }
        self.source[0] = pressure;
        let mut next = 1 + terminals.len();
        self.source[1..next].copy_from_slice(terminals);
        let mut groups = 0;
        for (group, values) in reservoirs.into_iter().enumerate() {
            if self.reservoir_channels.get(group) != Some(&values.len()) {
                return Err(sizing(
                    "bandlimited observation lost its original reservoir tile channels",
                ));
            }
            self.source[next..next + values.len()].copy_from_slice(values);
            next += values.len();
            groups += 1;
        }
        if groups != self.reservoir_channels.len() {
            return Err(sizing(
                "bandlimited observation lost an original reservoir surface",
            ));
        }
        self.resampler.push_frame(&self.source).map_err(sizing)
    }

    pub(in crate::nsf_dg) fn render_next(&mut self) -> Result<&[f64], RenderError> {
        self.resampler
            .render_next(&mut self.source)
            .map_err(sizing)?;
        Ok(&self.source)
    }
}

#[cfg(test)]
#[path = "bandlimited/tests.rs"]
mod tests;
