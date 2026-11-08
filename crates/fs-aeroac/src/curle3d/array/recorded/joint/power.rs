//! Source-bound sound power of one original compact/stress Welch recording.
//!
//! Every original factor enters the existing analytic Q/F/T sphere integral
//! coherently, including between the two physical source layouts. The native
//! Welch peak conversion supplies sqrt(2 df) exactly once; the existing power
//! operator supplies the peak-phasor time average. Summing its factor powers
//! gives watts in one bin, and division by the original df gives W/Hz.
//! No angular pressure-squared approximation or dense covariance is introduced.

use super::{RecordedJointSourceSpectrum, invalid, plane};
use crate::AeroacError;
use crate::curle3d::power::{
    CompactMultipolePhasor, PreparedMultipoleRadiatedPower, RadiatedPowerBudget,
    RadiatedPowerModeBudget,
};
use core::ops::Range;
use fs_fft::welch::WelchPeakBin;
use fs_math::{c64::C64, dd::Dd};

/// Authority of this lossless, source-bound Welch power observation.
pub const JOINT_RADIATED_POWER_SCOPE: &str = "Full-sphere sound power of one original joint compact volume/force and Lighthill-stress recording in a homogeneous stationary lossless exterior. All measured within-bin source/mechanism cross terms are retained. Complete native bins integrate W/Hz to W; no cross-frequency phase, confidence interval, source validation, unresolved turbulence, enclosed-volume double counting, motion, absorption, scattering, mutual loading or measured-SWL qualification is supplied.";

/// One native positive-frequency Welch bin after coherent spatial integration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointRadiatedPowerBin {
    /// Original Welch index; DC and the even-length Nyquist bin are excluded.
    pub bin: usize,
    pub frequency_hz: f64,
    /// Original uniform spectral spacing, not an independently chosen bandwidth.
    pub bin_width_hz: f64,
    /// Full-sphere radiated power assigned to this complete native bin [W].
    pub power_w: f64,
    /// One-sided sound-power spectral density, `power_w / bin_width_hz` [W/Hz].
    pub power_density_w_hz: f64,
    pub scope: &'static str,
}

/// Bounds the complete frequency sweep before any result allocation.
#[derive(Clone, Copy, Debug)]
pub struct JointRadiatedPowerBandBudget {
    /// Original combined compact-plus-stress source/pair preparation bounds.
    pub spatial: RadiatedPowerBudget,
    /// Original factor count and source-pair/factor bounds for EACH bin.
    pub modes: RadiatedPowerModeBudget,
    pub max_bins: usize,
    /// Combined source count squared times original segments times requested bins.
    /// Counts silent/coincident sources and zero factors as real admitted work.
    pub max_source_pair_mode_bin_products: usize,
}

/// An explicitly requested contiguous range of complete native spectral bins.
#[derive(Clone, Debug, PartialEq)]
pub struct JointRadiatedPowerBand {
    /// Ordered original bins; the same record/window/normalization supplies all.
    pub bins: Vec<JointRadiatedPowerBin>,
    /// Sum of the bin powers, reduced in double-double precision [W].
    pub power_w: f64,
    /// Edges of the selected complete DFT cells [Hz], not interpolated limits.
    /// For bins `a..b`, these are `(a-0.5)*df` and `(b-0.5)*df`.
    pub frequency_edges_hz: [f64; 2],
    pub scope: &'static str,
}

/// One prepared frequency bound by borrow to its original physical record.
///
/// Positions, source order, material, sample clock and measured factors cannot
/// be replaced independently. Preparation uses the existing analytic multipole
/// owner; repeated evaluation fills one retained O(N) source workspace, with no
/// allocation, trigonometry, FFT, covariance matrix or receiver quadrature.
/// Original factors remain coherent within each frequency; their powers add.
/// A cancelled evaluation publishes no partial result and can be retried.
pub struct PreparedJointRadiatedPower<'a> {
    record: &'a RecordedJointSourceSpectrum,
    bin: usize,
    operator: PreparedMultipoleRadiatedPower,
    sources: Vec<CompactMultipolePhasor>,
}

impl RecordedJointSourceSpectrum {
    /// Prepare full-sphere power from this record's actual two source layouts.
    ///
    /// Applies equally to manually acquired records and the original native
    /// FW-H surface/exterior-volume recorder's result. The compact sources are
    /// followed by the volume stresses without altering their geometry, physical
    /// units or phases. The caller remains responsible for disjoint source
    /// mechanisms; enclosed stresses must not be counted again outside FW-H.
    ///
    /// This is the homogeneous stationary LOSSLESS free-space power model.
    /// Molecular absorption or reflecting boundaries are not inferred from a
    /// separate microphone operator. Both even-length endpoints refuse; the
    /// highest bin of an odd record is an ordinary admitted positive frequency.
    pub fn prepare_radiated_power(
        &self,
        bin: usize,
        budget: RadiatedPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedJointRadiatedPower<'_>, AeroacError> {
        plane::poll(&mut cancelled)?;
        let frequency = self.bin_frequency(bin)?;
        let (count, _) = source_work(self, budget)?;
        let mut positions = Vec::new();
        positions
            .try_reserve_exact(count)
            .map_err(|_| invalid("joint power geometry allocation"))?;
        for &position in self
            .compact_positions_m
            .iter()
            .chain(&self.quadrupole_positions_m)
        {
            plane::poll(&mut cancelled)?;
            positions.push(position);
        }
        let operator = PreparedMultipoleRadiatedPower::new_cancellable(
            self.acoustic,
            frequency,
            &positions,
            budget,
            &mut cancelled,
        )?;
        let mut sources = Vec::new();
        sources
            .try_reserve_exact(count)
            .map_err(|_| invalid("joint power source workspace allocation"))?;
        for _ in 0..count {
            plane::poll(&mut cancelled)?;
            sources.push(CompactMultipolePhasor::default());
        }
        plane::poll(&mut cancelled)?;
        Ok(PreparedJointRadiatedPower {
            record: self,
            bin,
            operator,
            sources,
        })
    }

    /// Integrate an explicit nonempty range of complete interior Welch bins.
    ///
    /// Preparation and factor work are bounded over the WHOLE requested band,
    /// before allocation. One spatial operator/workspace is retained at a time:
    /// memory is O(N² + N + bins), rather than O(bins*N²) or channel covariance.
    /// No partial spectrum is returned after cancellation or numerical refusal.
    ///
    /// Powers add as stationary spectral second moments. This does not recover
    /// cross-frequency beats, infer independent overlapping records, interpolate
    /// partial bins, compensate leakage or supply missing DC/Nyquist power.
    pub fn radiated_power_band(
        &self,
        bins: Range<usize>,
        budget: JointRadiatedPowerBandBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<JointRadiatedPowerBand, AeroacError> {
        plane::poll(&mut cancelled)?;
        let count = bins
            .end
            .checked_sub(bins.start)
            .filter(|&n| n != 0 && n <= budget.max_bins)
            .ok_or_else(|| invalid("joint power requires a bounded nonempty bin range"))?;
        self.bin_frequency(bins.start)?;
        self.bin_frequency(bins.end - 1)?;
        let (_, pairs) = source_work(self, budget.spatial)?;
        let mode_work = admit_modes(pairs, self.spectrum.segment_count(), budget.modes)?;
        if mode_work
            .checked_mul(count)
            .is_none_or(|work| work > budget.max_source_pair_mode_bin_products)
        {
            return Err(invalid(
                "joint power band exceeds complete pair/factor/bin budget",
            ));
        }
        let df = self.spectrum.bin_width_hz();
        let frequency_edges_hz = [(bins.start as f64 - 0.5) * df, (bins.end as f64 - 0.5) * df];
        if frequency_edges_hz.iter().any(|edge| !edge.is_finite()) {
            return Err(invalid("joint power frequency edges are not representable"));
        }
        let mut rows = Vec::new();
        rows.try_reserve_exact(count)
            .map_err(|_| invalid("joint power band allocation"))?;
        let mut power = Dd::ZERO;
        for bin in bins {
            plane::poll(&mut cancelled)?;
            let row = self
                .prepare_radiated_power(bin, budget.spatial, &mut cancelled)?
                .evaluate(budget.modes, &mut cancelled)?;
            power = power + Dd::from_f64(row.power_w);
            rows.push(row);
        }
        plane::poll(&mut cancelled)?;
        let power_w = checked_power(power)?;
        Ok(JointRadiatedPowerBand {
            bins: rows,
            power_w,
            frequency_edges_hz,
            scope: JOINT_RADIATED_POWER_SCOPE,
        })
    }
}

impl PreparedJointRadiatedPower<'_> {
    /// Borrow the original record, including its source identity and sample clock.
    pub const fn record(&self) -> &RecordedJointSourceSpectrum {
        self.record
    }

    /// Native frequency index retained at preparation.
    pub const fn bin(&self) -> usize {
        self.bin
    }

    /// Borrow the exact analytic spatial operator without rebuilding its geometry.
    pub const fn power_operator(&self) -> &PreparedMultipoleRadiatedPower {
        &self.operator
    }

    /// Evaluate every original factor coherently, then add their physical powers.
    ///
    /// The one-sided `/sqrt(Hz)` factors use the native `sqrt(2*df)` peak
    /// conversion once per component. The original operator retains every
    /// compact/stress cross term and its peak time average; there is no extra
    /// factor of two, density, source volume, derivative or propagation delay.
    pub fn evaluate(
        &mut self,
        budget: RadiatedPowerModeBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<JointRadiatedPowerBin, AeroacError> {
        plane::poll(&mut cancelled)?;
        let count = self.sources.len();
        let pairs = count
            .checked_mul(count)
            .ok_or_else(|| invalid("joint power source-pair count overflow"))?;
        admit_modes(pairs, self.record.spectrum.segment_count(), budget)?;
        let frequency = self.operator.compact_power().frequency_hz();
        let measured = self
            .record
            .spectrum
            .peak_bin_hz(self.bin, frequency)
            .map_err(AeroacError::Welch)?;
        let compact_count = self.record.compact_positions_m.len();
        let mut power = Dd::ZERO;
        for factor in measured.factors() {
            plane::poll(&mut cancelled)?;
            for (source, channels) in self.sources[..compact_count]
                .iter_mut()
                .zip(factor[..4 * compact_count].chunks_exact(4))
            {
                plane::poll(&mut cancelled)?;
                source.compact.volume_flow_m3_s = physical(measured, channels[0])?;
                for (force, &coefficient) in source
                    .compact
                    .force_on_fluid_n
                    .iter_mut()
                    .zip(&channels[1..])
                {
                    *force = physical(measured, coefficient)?;
                }
            }
            for (source, channels) in self.sources[compact_count..]
                .iter_mut()
                .zip(factor[4 * compact_count..].chunks_exact(6))
            {
                plane::poll(&mut cancelled)?;
                for (stress, &coefficient) in
                    source.quadrupole.stress_moment_n_m.iter_mut().zip(channels)
                {
                    *stress = physical(measured, coefficient)?;
                }
            }
            let contribution = self
                .operator
                .evaluate_cancellable(&self.sources, &mut cancelled)?;
            power = power + Dd::from_f64(contribution.power_w);
        }
        plane::poll(&mut cancelled)?;
        let power_w = checked_power(power)?;
        let bin_width_hz = self.record.spectrum.bin_width_hz();
        let power_density_w_hz = (power / Dd::from_f64(bin_width_hz)).to_f64();
        if !power_density_w_hz.is_finite() || (power_w != 0.0 && power_density_w_hz == 0.0) {
            return Err(invalid("joint power spectral density is not representable"));
        }
        Ok(JointRadiatedPowerBin {
            bin: self.bin,
            frequency_hz: frequency,
            bin_width_hz,
            power_w,
            power_density_w_hz,
            scope: JOINT_RADIATED_POWER_SCOPE,
        })
    }
}

fn source_work(
    record: &RecordedJointSourceSpectrum,
    budget: RadiatedPowerBudget,
) -> Result<(usize, usize), AeroacError> {
    let count = record
        .compact_positions_m
        .len()
        .checked_add(record.quadrupole_positions_m.len())
        .ok_or_else(|| invalid("joint power source count overflow"))?;
    let pairs = count
        .checked_mul(count)
        .ok_or_else(|| invalid("joint power source-pair count overflow"))?;
    if count > budget.max_sources || pairs > budget.max_source_pairs {
        return Err(invalid("joint power exceeds combined source/pair budget"));
    }
    Ok((count, pairs))
}

fn admit_modes(
    pairs: usize,
    modes: usize,
    budget: RadiatedPowerModeBudget,
) -> Result<usize, AeroacError> {
    let work = pairs
        .checked_mul(modes)
        .filter(|&work| work <= budget.max_source_pair_mode_products)
        .ok_or_else(|| invalid("joint power factors exceed pair/factor budget"))?;
    if modes > budget.max_modes {
        return Err(invalid("joint power factors exceed mode count budget"));
    }
    Ok(work)
}

fn physical(bin: WelchPeakBin<'_>, value: fs_fft::C64) -> Result<C64, AeroacError> {
    let value = bin.peak_factor(value).map_err(AeroacError::Welch)?;
    Ok(C64::new(value.re, value.im))
}

fn checked_power(power: Dd) -> Result<f64, AeroacError> {
    let value = power.to_f64();
    if !power.hi.is_finite() || !power.lo.is_finite() || !value.is_finite() || value < 0.0 {
        return Err(invalid(
            "joint radiated power is not representable and nonnegative",
        ));
    }
    Ok(value)
}
