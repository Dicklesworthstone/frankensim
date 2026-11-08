//! Original joint compact/stress recordings through a finite fluid layer.

use super::super::super::array::admit_record;
use super::super::{LayerQuadrupoleField, PreparedLayeredFullWaveQuadrupoleObserver};
use super::*;
use crate::curle3d::statistics::{ComplexSum, StatisticsSum};
use fs_fft::welch::WelchCrossSpectrum;

const MIXED_TRANSMITTED_SCOPE: &str = "Jointly recorded stationary compact volume \
flows [m³/s], forces on fluid [N] and symmetric Lighthill stress moments [N m]; \
full spherical-wave transmission through an infinite planar homogeneous lossless \
nontrapping fluid layer between original stationary fluids. Internal returns, \
near/evanescent fields and all source/mechanism phases remain coherent. Disjoint \
physical mechanisms without permeable-surface double counting; original clock \
and channel order. Bin-integrated Welch pressure/velocity moments retain finite-record \
and window limits, not deterministic cross-frequency phases. Per-source estimates \
give no relative guarantee at interference nulls. No source-size or SGS validation, \
absorption, mean flow, finite edges, source loading, interval certificate or \
measured-SPL qualification.";

const MIXED_REFLECTED_SCOPE: &str = "Jointly recorded stationary compact volume \
flows [m³/s], forces on fluid [N] and symmetric Lighthill stress moments [N m]; \
reflected-only spherical-wave field before an infinite planar homogeneous lossless \
nontrapping fluid layer. Every internal return retains its phase within the original \
joint factor; disjoint physical mechanisms without permeable-surface double counting. \
Bin-integrated Welch pressure/velocity moments retain original clock, channel order \
and finite-record limits. No direct field, deterministic cross-frequency phases, \
source-size or SGS validation, absorption, mean flow, finite edges, source loading, \
relative guarantee at interference nulls, interval certificate or measured-SPL qualification.";

const MIXED_TOTAL_SCOPE: &str = "Jointly recorded stationary compact volume flows \
[m³/s], forces on fluid [N] and symmetric Lighthill stress moments [N m]; coherent \
direct-plus-reflected spherical-wave field before an infinite planar homogeneous \
lossless nontrapping fluid layer. Original clock, channel order and correlation retain \
all tensor, mechanism and internal/direct-path cross terms. Disjoint physical sources \
without permeable-surface double counting. Bin-integrated Welch moments retain \
finite-record limits, not deterministic cross-frequency phases. Reflected-transfer \
estimates give no relative total-field guarantee at nulls. No singular total self-field, \
source-size or SGS validation, absorption, mean flow, finite edges, source loading, \
interval certificate or measured-SPL qualification.";

impl PreparedLayeredFullWaveCompactSourceObserver {
    /// Propagate one original interior Welch bin through the prepared layer field.
    ///
    /// Source-major channels are [Q,Fx,Fy,Fz] in m³/s and N, in the prepared
    /// original position order. The caller binds that geometry; all channels
    /// share one clock, window and segment history. Missing mechanisms are zero.
    ///
    /// The shared native factor reducer keeps Q/F, source and internal-path
    /// interference before forming full pressure/velocity moments. sqrt(2 df)
    /// is applied once, and the recorded bin must exactly match this frequency.
    /// Results retain bin-integrated power and finite-record limitations;
    /// no new FFT, quadrature, propagation delay or independence is introduced.
    /// Application is allocation-free with bounded work and cancellation.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view()
            .evaluate_welch_bin(spectrum, bin, budget, cancelled)
    }

    /// Observe one synchronized compact AND volume-stress record through this layer.
    ///
    /// Compact [Q,Fx,Fy,Fz] blocks [m³/s,N] precede physical stress blocks
    /// [xx,yy,zz,xy,xz,yz] [N m]. Each bank keeps its own original source order;
    /// stress moments already include cell volume and have no engineering-shear
    /// multiplier. Supply disjoint mechanisms without permeable-surface recounting.
    ///
    /// Both observers must share all three original fluids, first face and
    /// oriented normal, layer thickness, receiver, exact frequency and the
    /// transmitted, reflected-only or incident-total field. The original native
    /// maps retain every internal return. Both banks enter one DD pressure and
    /// vector-velocity sum within each factor BEFORE quadratic moments.
    ///
    /// The combined source/segment budget is admitted before propagation and
    /// peak scaling occurs once. Original window/aliasing limits remain;
    /// separate spectra cannot reconstruct missing cross-mechanism phase.
    /// Application allocates nothing and introduces no waveform continuation,
    /// cross-frequency phase, extra FFT or propagation quadrature.
    pub fn evaluate_welch_bin_with_quadrupoles(
        &self,
        quadrupoles: &PreparedLayeredFullWaveQuadrupoleObserver,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        let compact_layer = self.layer();
        let stress_layer = quadrupoles.layer();
        let matching_field = matches!(
            (self.field_kind(), quadrupoles.field_kind()),
            (FullWaveCompactFieldKind::Transmitted, LayerQuadrupoleField::Transmitted)
                | (FullWaveCompactFieldKind::Reflected, LayerQuadrupoleField::Reflected)
                | (FullWaveCompactFieldKind::IncidentTotal, LayerQuadrupoleField::IncidentTotal)
        );
        if compact_layer.incident_state() != stress_layer.incident_state()
            || compact_layer.layer_state() != stress_layer.layer_state()
            || compact_layer.receiving_state() != stress_layer.receiving_state()
            || compact_layer.first_face_point_m() != stress_layer.first_face_point_m()
            || compact_layer.receiving_normal() != stress_layer.receiving_normal()
            || compact_layer.thickness_m() != stress_layer.thickness_m()
            || self.observer_position_m() != quadrupoles.observer_position_m()
            || self.frequency_hz().to_bits() != quadrupoles.frequency_hz().to_bits()
            || !matching_field
        {
            return Err(invalid(
                "mixed recorded sources require one physical layer, receiver, frequency and field",
            ));
        }
        let compact_channels = self.source_count().checked_mul(4)
            .ok_or_else(|| invalid("layer recorded compact channel count overflow"))?;
        let channels = quadrupoles.source_count().checked_mul(6)
            .and_then(|count| count.checked_add(compact_channels))
            .ok_or_else(|| invalid("layer recorded mixed channel count overflow"))?;
        let sources = self.source_count().checked_add(quadrupoles.source_count())
            .ok_or_else(|| invalid("layer recorded mixed source count overflow"))?;
        let measured = admit_record(
            spectrum, bin, self.frequency_hz(), sources, channels, budget,
        )?;
        let mut statistics = StatisticsSum::default();
        for factor in measured.factors() {
            poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            let (compact, stress) = factor.split_at(compact_channels);
            self.view()
                .accumulate_welch_factor(compact, measured, &mut sums, &mut cancelled)?;
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
}
