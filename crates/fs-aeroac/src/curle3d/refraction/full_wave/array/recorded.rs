//! Original measured source factors through prepared refracting p/u transfers.

use super::*;
use crate::curle3d::refraction::full_wave::quadrupole::array::{
    FullWaveQuadrupoleField, PreparedFullWaveQuadrupoleObserver,
};
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

const MIXED_TRANSMITTED_SCOPE: &str = "Coherently recorded stationary compact volume \
flows [m³/s], forces on fluid [N] and symmetric Lighthill stress moments [N m]; full-wave \
transmission across one infinite plane between original stationary homogeneous lossless \
fluids, retaining near, critical and evanescent fields. Disjoint physical mechanisms \
without permeable-surface double counting; original clocks and channel order. Bin-integrated \
Welch pressure/velocity moments, not deterministic cross-frequency phases. Per-source \
quadrature estimates and tail bounds give no relative guarantee at interference nulls. \
No source-size or SGS validation, absorption, mean flow, finite boundaries, source loading, \
interval certificate or measured-SPL qualification.";

const MIXED_REFLECTED_SCOPE: &str = "Coherently recorded stationary compact volume flows \
[m³/s], forces on fluid [N] and symmetric Lighthill stress moments [N m]; reflected-only \
full-wave field in the incident fluid at one infinite planar lossless fluid interface, \
retaining near, critical and evanescent fields. Disjoint physical mechanisms without \
permeable-surface double counting; original clocks and channel order. Bin-integrated \
Welch pressure/velocity moments, not deterministic cross-frequency phases. Per-source \
quadrature estimates and tail bounds give no relative guarantee at interference nulls. \
No direct field, source-size or SGS validation, absorption, mean flow, finite boundaries, \
source loading, interval certificate or measured-SPL qualification.";

const MIXED_TOTAL_SCOPE: &str = "Coherently recorded stationary compact volume flows \
[m³/s], forces on fluid [N] and symmetric Lighthill stress moments [N m]; coherent \
direct-plus-reflected full-wave field in the incident fluid at one infinite planar lossless \
fluid interface. Disjoint physical mechanisms without permeable-surface double counting; \
original clocks and channel order. Bin-integrated Welch pressure/velocity moments, not \
deterministic cross-frequency phases. Reflected-transfer estimates give no relative \
guarantee at direct/reflected or cross-mechanism interference nulls. No singular direct \
self-field, source-size or SGS validation, absorption, mean flow, finite boundaries, \
source loading, interval certificate or measured-SPL qualification.";

impl PreparedFullWaveCompactSourceObserver {
    /// Propagate one original measured interior Welch bin through this field.
    ///
    /// Source-major channels are `[Q, Fx, Fy, Fz]`: volume injection [m³/s] and
    /// force ON fluid [N], in exactly the prepared source order. All channels
    /// share one sample clock, window and segment boundaries. The caller binds
    /// the recorded source geometry; absent mechanisms are zero channels.
    ///
    /// Original source and mechanism cross terms survive in pressure and full
    /// vector velocity for transmitted, reflected-only and incident-total fields.
    /// The original peak-bin owner applies sqrt(2 df) once per coefficient. The
    /// result is a bin-integrated second moment, not PSD/Hz, a tone-amplitude fit
    /// or deterministic phases across frequencies. All finite-record, window,
    /// aliasing and prepared propagation limitations remain.
    ///
    /// Allocation-free O(sources * segments), with checked work and bounded
    /// cancellation. Cached transfers are reused without new FFTs, quadrature,
    /// source vectors, propagation delays or dense source covariance matrices.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        let channels = self
            .source_count()
            .checked_mul(4)
            .ok_or_else(|| invalid("refracting recorded compact channel count overflow"))?;
        let measured = admit_record(
            spectrum,
            bin,
            self.frequency_hz(),
            self.source_count(),
            channels,
            budget,
        )?;
        let mut statistics = StatisticsSum::default();
        for factor in measured.factors() {
            poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            self.accumulate_welch_factor(factor, measured, &mut sums, &mut cancelled)?;
            statistics.add(sums.map(ComplexSum::value))?;
        }
        poll(&mut cancelled)?;
        statistics.finish(self.scope())
    }

    /// Coherently observe original recorded compact AND volume-stress sources.
    ///
    /// One synchronized Welch record contains source-major compact `[Q,Fx,Fy,Fz]`
    /// blocks [m³/s,N], followed by stress `[xx,yy,zz,xy,xz,yz]` blocks [N m].
    /// Each group follows its own prepared positions in order; physical shear
    /// entries are not pre-doubled and cell volumes are already included. Supply
    /// disjoint mechanisms: a volume already represented by a permeable surface
    /// or another source model must not be counted again.
    ///
    /// Both prepared observers must retain identical original fluid states,
    /// plane origin and oriented normal, receiver, exact frequency and field
    /// kind. Transmitted, reflected-only and incident-total fields are admitted;
    /// mixing their different physical observations refuses.
    ///
    /// For each original factor, BOTH banks enter one double-double pressure
    /// and vector-velocity reduction BEFORE moments are formed. Cross-cell,
    /// cross-component, Q/F/stress and direct/reflected interference survive.
    /// Peak scaling is applied once, on the common original clock. Unknown
    /// correlation is never treated as independence; separate spectra cannot
    /// reconstruct the missing cross-mechanism phase.
    ///
    /// Bin-integrated moments retain the original window and finite-record
    /// limitations, not deterministic cross-frequency beats or inferred SGS
    /// radiation. Allocation-free O((compact sources + cells) * segments), with
    /// a checked COMBINED source/factor budget and bounded cancellation. Neither
    /// source solution, FFT, propagation quadrature nor covariance is rebuilt.
    pub fn evaluate_welch_bin_with_quadrupoles(
        &self,
        quadrupoles: &PreparedFullWaveQuadrupoleObserver,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        let compact_interface = self.interface();
        let stress_interface = quadrupoles.interface();
        let matching_field = matches!(
            (self.field_kind(), quadrupoles.field_kind()),
            (
                FullWaveCompactFieldKind::Transmitted,
                FullWaveQuadrupoleField::Transmitted
            ) | (
                FullWaveCompactFieldKind::Reflected,
                FullWaveQuadrupoleField::Reflected
            ) | (
                FullWaveCompactFieldKind::IncidentTotal,
                FullWaveQuadrupoleField::IncidentTotal
            )
        );
        if compact_interface.incident_state() != stress_interface.incident_state()
            || compact_interface.receiving_state() != stress_interface.receiving_state()
            || compact_interface.point_m() != stress_interface.point_m()
            || compact_interface.receiving_normal() != stress_interface.receiving_normal()
            || self.observer_position_m() != quadrupoles.observer_position_m()
            || self.frequency_hz().to_bits() != quadrupoles.frequency_hz().to_bits()
            || !matching_field
        {
            return Err(invalid(
                "mixed recorded sources require one interface, receiver, frequency and field kind",
            ));
        }
        let compact_channels = self
            .source_count()
            .checked_mul(4)
            .ok_or_else(|| invalid("refracting recorded compact channel count overflow"))?;
        let channels = quadrupoles
            .source_count()
            .checked_mul(6)
            .and_then(|n| n.checked_add(compact_channels))
            .ok_or_else(|| invalid("refracting recorded mixed channel count overflow"))?;
        let sources = self
            .source_count()
            .checked_add(quadrupoles.source_count())
            .ok_or_else(|| invalid("refracting recorded mixed source count overflow"))?;
        let measured = admit_record(
            spectrum,
            bin,
            self.frequency_hz(),
            sources,
            channels,
            budget,
        )?;
        let mut statistics = StatisticsSum::default();
        for factor in measured.factors() {
            poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            let (compact, stress) = factor.split_at(compact_channels);
            self.accumulate_welch_factor(compact, measured, &mut sums, &mut cancelled)?;
            quadrupoles.accumulate_welch_factor(stress, measured, &mut sums, &mut cancelled)?;
            statistics.add(sums.map(ComplexSum::value))?;
        }
        poll(&mut cancelled)?;
        statistics.finish(match self.field_kind() {
            FullWaveCompactFieldKind::Transmitted => MIXED_TRANSMITTED_SCOPE,
            FullWaveCompactFieldKind::Reflected => MIXED_REFLECTED_SCOPE,
            FullWaveCompactFieldKind::IncidentTotal => MIXED_TOTAL_SCOPE,
        })
    }

    fn accumulate_welch_factor(
        &self,
        factor: &[fs_fft::C64],
        measured: WelchPeakBin<'_>,
        sums: &mut [ComplexSum; 4],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        for (index, (transfer, channels)) in self
            .transfers
            .iter()
            .zip(factor.chunks_exact(4))
            .enumerate()
        {
            if index.is_multiple_of(256) {
                poll(cancelled)?;
            }
            let mut values = [C64::ZERO; 4];
            for (coefficient, value) in channels.iter().zip(&mut values) {
                let peak = measured
                    .peak_factor(*coefficient)
                    .map_err(AeroacError::Welch)?;
                *value = C64::new(peak.re, peak.im);
            }
            transfer.accumulate(
                CompactSourcePhasor {
                    volume_flow_m3_s: values[0],
                    force_on_fluid_n: [values[1], values[2], values[3]],
                },
                sums,
            )?;
        }
        Ok(())
    }
}

fn admit_record(
    spectrum: &WelchCrossSpectrum,
    bin: usize,
    frequency_hz: f64,
    sources: usize,
    channels: usize,
    budget: CoherenceBudget,
) -> Result<WelchPeakBin<'_>, AeroacError> {
    let measured = spectrum
        .peak_bin_hz(bin, frequency_hz)
        .map_err(AeroacError::Welch)?;
    if channels != measured.channel_count() {
        return Err(invalid(
            "refracting recorded source channels must match both prepared source orders",
        ));
    }
    if measured.segment_count() > budget.max_modes
        || sources
            .checked_mul(measured.segment_count())
            .is_none_or(|work| work > budget.max_source_mode_pairs)
    {
        return Err(invalid(
            "combined refracting recorded source coherence exceeds work budget",
        ));
    }
    Ok(measured)
}
