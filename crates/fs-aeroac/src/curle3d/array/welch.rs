//! Measured source factors use the original prepared pressure/velocity paths.

use super::*;
use crate::curle3d::quadrupole::PreparedQuadrupoleObserver;
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

const MIXED_QUADRUPOLE_SCOPE: &str = "Coherently recorded 3D compact volume/force sources and \
symmetric Lighthill volume stresses; homogeneous stationary lossless free-space point quadrature; \
caller supplies disjoint physical mechanisms with no permeable-surface double counting; \
bin-integrated Welch pressure/velocity moments, not deterministic cross-frequency phases; \
no source validation, boundaries, moving sources, mean-flow propagation, absorption, \
unresolved SGS sound or measured-SPL qualification.";

const MIXED_REFLECTED_QUADRUPOLE_SCOPE: &str = "Coherently recorded 3D compact volume/force \
sources and symmetric Lighthill volume stresses in one homogeneous stationary lossless \
ideal perpendicular-wall or pi/m-wedge domain; all physical source and image paths retain \
their measured coherence. Caller supplies disjoint mechanisms without double counting. \
Bin-integrated Welch pressure/velocity moments, not deterministic cross-frequency phases; \
no source validation, finite wall impedance, general diffraction, enclosed-room reverberation, \
motion, absorption, source feedback, unresolved SGS sound or measured-SPL qualification.";

impl PreparedCompactSourceObserver {
    /// Propagate one measured interior Welch bin through the original medium.
    ///
    /// Channels are source-major `[Q, Fx, Fy, Fz]`: physical volume injection
    /// in m³/s followed by force ON the fluid in N for each prepared position.
    /// All channels must share the original sampling clock. A missing mechanism
    /// is an explicitly zero channel; no span conversion or force is inferred.
    ///
    /// Each recorded factor is propagated coherently through the same cached
    /// pressure/velocity coefficients, including any molecular dispersion or
    /// ideal-plane image paths. Source, force-axis and volume/force cross terms
    /// survive. Images never become additional independent source channels.
    ///
    /// The one-sided /Hz factors are scaled by sqrt(2 df) exactly once before
    /// the existing peak-phasor statistics. Results are bin-integrated second
    /// moments, not PSD/Hz or coherent-gain tone amplitudes. Window leakage and
    /// the original source/propagation scope remain. Adding bins assumes a
    /// stationary resolved spectrum; this supplies no cross-frequency beats,
    /// confidence interval or statistical independence of overlapping segments.
    ///
    /// Allocation-free O(sources * segments), with source-boundary cancellation
    /// and no FFT, Green-function evaluation or dense source covariance.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        plane::poll(&mut cancelled)?;
        let measured = admit(
            spectrum,
            bin,
            self.frequency_hz,
            self.source_count(),
            4,
            1,
            budget,
        )?;
        let mut statistics = StatisticsSum::default();
        for factor in measured.factors() {
            plane::poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            for (transfer, channels) in self.transfers.iter().zip(factor.chunks_exact(4)) {
                plane::poll(&mut cancelled)?;
                transfer.accumulate(source(channels, measured)?, &mut sums)?;
            }
            statistics.add(sums.map(ComplexSum::value))?;
        }
        plane::poll(&mut cancelled)?;
        statistics.finish(self.scope())
    }

    /// Coherently observe recorded compact sources AND volume-stress sources.
    ///
    /// One synchronized record contains source-major compact `[Q, Fx, Fy, Fz]`
    /// blocks (m³/s, N), followed by quadrupole `[xx, yy, zz, xy, xz, yz]`
    /// blocks (N m). Each group follows its own prepared source order; shear
    /// entries are not pre-doubled. Every channel uses the same clock, window
    /// and segment boundaries. Supply disjoint physical mechanisms: adding a
    /// volume source already represented by a permeable acoustic surface or
    /// another source model would double count its radiation.
    ///
    /// For each original Welch factor, BOTH mechanisms enter one pressure and
    /// full vector-velocity sum before second moments are formed. Cross-cell,
    /// cross-component and cross-mechanism interference therefore survive.
    /// The original peak-equivalent scaling is applied once per coefficient.
    /// Outputs are bin-integrated moments, not PSD/Hz or deterministic waveform
    /// phases across bins. Unknown source correlation is never independence.
    ///
    /// Both observers must have exactly the same acoustic state, receiver and
    /// frequency and identical ideal boundary descriptions (including wedge
    /// order and plane ordering). Free space, perpendicular walls and pi/m
    /// wedges are admitted. Molecular dispersion/absorption remains unavailable
    /// to the quadrupole owner. All source-model, compact-cell, finite-record,
    /// window and aliasing limits remain; no SGS sound is inferred.
    ///
    /// Allocation-free O((compact sources + cells) * segments), with bounded
    /// cancellation and a checked COMBINED source/factor work budget. Existing
    /// cached transfers and double-double reductions are reused; no second
    /// spectrum, source-factor vectors or transform is constructed.
    pub fn evaluate_welch_bin_with_quadrupoles(
        &self,
        quadrupoles: &PreparedQuadrupoleObserver,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        plane::poll(&mut cancelled)?;
        if self.acoustic() != quadrupoles.acoustic()
            || self.observer_position_m() != quadrupoles.observer_position_m()
            || self.frequency_hz().to_bits() != quadrupoles.frequency_hz().to_bits()
            || self.molecular_state().is_some()
            || self.reflection_planes() != quadrupoles.reflection_planes()
            || self.reflection_wedge_order() != quadrupoles.reflection_wedge_order()
        {
            return Err(invalid(
                "mixed recorded sources require one lossless medium, receiver, frequency and boundary domain",
            ));
        }
        let measured = spectrum
            .peak_bin_hz(bin, self.frequency_hz())
            .map_err(AeroacError::Welch)?;
        let compact_channels = self
            .source_count()
            .checked_mul(4)
            .ok_or_else(|| invalid("recorded compact-source channel count overflow"))?;
        if quadrupoles
            .source_count()
            .checked_mul(6)
            .and_then(|n| n.checked_add(compact_channels))
            != Some(measured.channel_count())
        {
            return Err(invalid(
                "mixed recorded sources require ordered compact then quadrupole channels",
            ));
        }
        if measured.segment_count() > budget.max_modes
            || self
                .source_count()
                .checked_add(quadrupoles.source_count())
                .and_then(|n| n.checked_mul(measured.segment_count()))
                .is_none_or(|work| work > budget.max_source_mode_pairs)
        {
            return Err(invalid(
                "combined recorded source coherence exceeds work budget",
            ));
        }
        let mut statistics = StatisticsSum::default();
        for factor in measured.factors() {
            plane::poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            let (compact, stress) = factor.split_at(compact_channels);
            for (transfer, channels) in self.transfers.iter().zip(compact.chunks_exact(4)) {
                plane::poll(&mut cancelled)?;
                transfer.accumulate(source(channels, measured)?, &mut sums)?;
            }
            quadrupoles.accumulate_welch_factor(stress, measured, &mut sums, &mut cancelled)?;
            statistics.add(sums.map(ComplexSum::value))?;
        }
        plane::poll(&mut cancelled)?;
        statistics.finish(if self.reflection_planes().is_empty() {
            MIXED_QUADRUPOLE_SCOPE
        } else {
            MIXED_REFLECTED_QUADRUPOLE_SCOPE
        })
    }
}

/// Same measured-bin admission for compact, half-space and convected owners.
/// `paths_per_source` counts separately applied paths, not folded image paths.
pub(super) fn admit(
    spectrum: &WelchCrossSpectrum,
    bin: usize,
    frequency_hz: f64,
    sources: usize,
    channels_per_source: usize,
    paths_per_source: usize,
    budget: CoherenceBudget,
) -> Result<WelchPeakBin<'_>, AeroacError> {
    let measured = spectrum
        .peak_bin_hz(bin, frequency_hz)
        .map_err(|_| invalid("recorded compact source bin/frequency"))?;
    if sources.checked_mul(channels_per_source) != Some(measured.channel_count()) {
        return Err(invalid(
            "Welch source channels must match the prepared compact observer",
        ));
    }
    if measured.segment_count() > budget.max_modes
        || sources
            .checked_mul(measured.segment_count())
            .and_then(|work| work.checked_mul(paths_per_source))
            .is_none_or(|work| work > budget.max_source_mode_pairs)
    {
        return Err(invalid(
            "recorded compact source coherence exceeds work budget",
        ));
    }
    Ok(measured)
}

pub(super) fn source(
    channels: &[fs_fft::C64],
    bin: WelchPeakBin<'_>,
) -> Result<CompactSourcePhasor, AeroacError> {
    let value = |i: usize| -> Result<C64, AeroacError> {
        let peak = bin
            .peak_factor(channels[i])
            .map_err(|_| invalid("recorded compact source peak factor is unrepresentable"))?;
        Ok(C64::new(peak.re, peak.im))
    };
    Ok(CompactSourcePhasor {
        volume_flow_m3_s: value(0)?,
        force_on_fluid_n: [value(1)?, value(2)?, value(3)?],
    })
}
