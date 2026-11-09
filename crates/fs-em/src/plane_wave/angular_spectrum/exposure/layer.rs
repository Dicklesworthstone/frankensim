//! Finite-window local film heat and reversible field/polarization storage.
//!
//! These are settled prescribed harmonics through the original film operator,
//! not a switch-on transient. Ohmic and each original polarization reservoir dissipate
//! after all spatial and temporal interference. Boundary power loss also
//! contains the signed change of electromagnetic/polarization storage.
//! Admitted uniform and triangular bands evaluate field populations and reservoir
//! heat through shared orthogonal temporal modes. Native angular frequencies
//! and the distinction between real fields and proper factors are retained;
//! a real band admits its full carrier, not just its slow beat. Wider
//! windows keep the original pair law. Signed flux and reversible storage
//! use cross products of the same temporal modes, without a positivity clamp.
//! The factor is not a null-relative error certificate and adds no material
//! or switch-on approximation.

use super::*;
use crate::plane_wave::layers::FIELD_PREPARATION_WORKSPACE_BYTES;
use fs_fft::harmonic_window::{
    HarmonicBandRef, HarmonicBandWorkspace, MAX_UNIFORM_HARMONIC_BAND_MODES,
    MAX_UNIFORM_HARMONIC_BAND_TONES,
};

mod reservoir;
use reservoir::{LossLaw, PoleTransfer};

mod temporal_moments;

mod energy;
pub use energy::{
    AngularSpectrumEnergyVolumeLimits, AngularSpectrumEnergyVolumeSample,
    AngularSpectrumExposureLayerEnergyDensity, AngularSpectrumExposureLayerEnergyLimits,
    AngularSpectrumExposureLayerVolumeEnergy, PreparedAngularSpectrumExposureLayerEnergy,
    PreparedAngularSpectrumExposureLayerVolumeEnergy,
};

mod integrated_heat;
pub use integrated_heat::{
    AngularSpectrumExposureLayerHeatLoadSpectrum, AngularSpectrumHeatAreaQuadrature,
    AngularSpectrumHeatDepthSample, AngularSpectrumHeatVolumeQuadrature,
    AngularSpectrumHeatWindowLimits, AngularSpectrumIntegratedHeatLimits,
};

mod heat_spectrum;
pub use heat_spectrum::{
    AngularSpectrumExposureHeatChannel, AngularSpectrumExposureHeatSpectrumLimits,
    AngularSpectrumExposureLayerHeatSpectrum,
};

/// Whole-bank preparation and depth-query limits. A frequency group belongs
/// to one mean/factor bank; equal frequencies across banks share one film cache.
/// At most 4096 authored groups are admitted, bounding the preparation sort.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureLayerLimits {
    pub maximum_frequency_groups: usize,
    pub maximum_frequency_layer_mode_pairs: usize,
    pub maximum_point_line_pairs: usize,
    /// N * sum(F_r²) * (3 + D + 2*L), with D Debye and L Lorentz poles, plus
    /// sum_real F_r*N*(width+height) for a nonzero common carrier. The three
    /// base channels cover the vector field and its instantaneous E/H stores.
    /// Eligible short-band banks reserve N*(6+3*D+5*L)*(F_r+1)*17 for
    /// vector value/derivative projections and signed cross contractions.
    /// At most 64 proper or 32 real tones use this bounded factor. This bounds
    /// temporal projections before the query's shutter is known.
    pub maximum_point_reservoir_frequency_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// Mean local irreversible heat and signed reversible storage rate [W/m³].
/// Storage includes epsilon_infinity |E|²/2, mu |H|²/2 and each
/// |P_p|²/(2 Delta_epsilon_p), plus each Lorentz pole's
/// |P_dot_p|²/(2 omega0_p² Delta_epsilon_p). It excludes thermal energy
/// and parameter work. Real fields are contracted before window averaging.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureLayerDissipation {
    /// Window/ensemble mean of sigma |E_real|², always nonnegative.
    pub ohmic_heat_w_m3: f64,
    /// Window/ensemble mean of the original poles' irreversible friction heat.
    pub polarization_heat_w_m3: f64,
    /// Mean irreversible heating rate; excludes reversible storage changes.
    pub total_heat_w_m3: f64,
    /// Signed detector-weighted mean dU/dt. Only a uniform window gives
    /// (U(end)-U(start))/duration; a shaped window retains its sensitivity.
    /// This is reversible storage rate, not deposited heat.
    pub storage_rate_w_m3: f64,
}

/// Finite-window field/heat on one depth plane of the selected original film.
/// Transverse heat/storage integrals are original-grid quadratures in W/m of
/// depth, not total film watts; their products can require spatial refinement.
/// The field's normal watts use continuous modal cell integration separately.
#[derive(Debug)]
pub struct AngularSpectrumExposureLayerField {
    field: AngularSpectrumExposureField,
    dissipation: Vec<AngularSpectrumExposureLayerDissipation>,
    absorbed_power_per_depth_w_m: f64,
    storage_rate_per_depth_w_m: f64,
}
impl AngularSpectrumExposureLayerField {
    pub fn field(&self) -> &AngularSpectrumExposureField {
        &self.field
    }
    pub fn dissipation(&self) -> &[AngularSpectrumExposureLayerDissipation] {
        &self.dissipation
    }
    pub fn absorbed_power_per_depth_w_m(&self) -> f64 {
        self.absorbed_power_per_depth_w_m
    }
    pub fn storage_rate_per_depth_w_m(&self) -> f64 {
        self.storage_rate_per_depth_w_m
    }
}

struct FrequencyLayer<'b, 's, 'm> {
    omega: f64,
    layer: Option<PreparedAngularSpectrumLayer<'b, 's, 'm>>,
    poles: Vec<PoleTransfer>,
}

/// One selected layer with one original TE/TM cache per distinct frequency
/// needing a one-sided response, shared across the mean and covariance factors.
/// Two-sided lines reuse their already prepared native two-port columns.
/// Original source states and frequency groups remain borrowed. Repeated depths
/// do not solve the films or reevaluate constitutive laws again.
pub struct PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
    anchor: &'b PlanarAngularSpectrum<'s, 'm>,
    mean: Option<&'e AngularSpectrumExposure<'b, 's, 'm>>,
    factors: &'e [AngularSpectrumExposure<'b, 's, 'm>],
    statistics: Option<&'e [AngularSpectrumFactorStatistics]>,
    layer_index: usize,
    thickness_m: f64,
    law: LossLaw<'m>,
    frequencies: Vec<FrequencyLayer<'b, 's, 'm>>,
    workspace_bytes: usize,
    preparation_layer_mode_pairs: usize,
    point_line_pairs: usize,
    point_reservoir_frequency_pairs: usize,
    field_fft_work: u64,
}

impl<'b, 's, 'm> AngularSpectrumExposure<'b, 's, 'm> {
    pub fn prepare_layer(
        &self,
        layer_index: usize,
        limits: AngularSpectrumExposureLayerLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumExposureLayer<'_, 'b, 's, 'm>, AngularSpectrumError> {
        PreparedAngularSpectrumExposureLayer::new(
            self.anchor,
            Some(self),
            &[],
            None,
            layer_index,
            limits,
            &mut cancelled,
        )
    }
}

impl<'e, 'b, 's, 'm> AngularSpectrumExposureEnsemble<'e, 'b, 's, 'm> {
    pub fn prepare_layer(
        &self,
        layer_index: usize,
        limits: AngularSpectrumExposureLayerLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm>, AngularSpectrumError> {
        PreparedAngularSpectrumExposureLayer::new(
            self.anchor,
            self.mean,
            self.factors,
            self.statistics,
            layer_index,
            limits,
            &mut cancelled,
        )
    }
}

#[derive(Clone, Copy)]
struct MomentKernel<P> {
    product: P,
    // For proper fields this is half the exact frequency difference, so a
    // common large carrier never enters a slow covariance-storage derivative.
    storage_frequency: Dd,
}

// Real peak phasors use BOTH signed frequencies. Their full carrier span must
// fit the shared factor; a short beat alone never licenses cycle averaging.
fn positive_band_tones(frequencies: usize, proper: bool) -> Option<usize> {
    if frequencies < 2 {
        return None;
    }
    let maximum = MAX_UNIFORM_HARMONIC_BAND_TONES / if proper { 1 } else { 2 };
    (frequencies <= maximum).then_some(frequencies)
}

struct PositiveLayerMoments<'w> {
    factor: HarmonicBandRef<'w>,
    frequencies: usize,
    proper: bool,
    derivative_frequencies: [Dd; MAX_UNIFORM_HARMONIC_BAND_TONES],
}

impl<'w> PositiveLayerMoments<'w> {
    fn prepare(
        window: HarmonicWindow,
        frequencies: &[f64],
        proper: bool,
        workspace: &'w mut HarmonicBandWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<Self>, AngularSpectrumError> {
        let Some(_) = positive_band_tones(frequencies.len(), proper) else {
            return Ok(None);
        };
        poll(cancelled)?;
        let factor = if proper {
            window.prepare_band_rad_s_in(frequencies, workspace)
        } else {
            window.prepare_real_band_rad_s_in(frequencies, workspace)
        }
        .map_err(|_| bad("film exposure positive temporal factor is unrepresentable"))?;
        poll(cancelled)?;
        Ok(factor.map(|factor| {
            // A common proper carrier contributes only an imaginary multiple
            // of each self product, hence exactly zero real storage rate.
            // Remove it before differentiation, never after large products.
            let reference = if proper {
                Dd::from_f64(frequencies[0])
            } else {
                Dd::ZERO
            };
            let mut derivative_frequencies = [Dd::ZERO; MAX_UNIFORM_HARMONIC_BAND_TONES];
            for (target, &frequency) in derivative_frequencies.iter_mut().zip(frequencies) {
                *target = Dd::from_f64(frequency) - reference;
            }
            Self {
                factor,
                frequencies: frequencies.len(),
                proper,
                derivative_frequencies,
            }
        }))
    }

    // Inputs keep the native exp(+i*omega*t) EM convention. Proper random
    // factors have <Re(E)^2>=<|E|^2>/2. The real factor uses both signed
    // frequencies algebraically, retaining every physical sum-frequency term.
    fn square(
        &self,
        field: impl Fn(usize) -> [Dd; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, MomentError> {
        self.square_scaled(field, Dd::ONE, cancelled)
    }

    // Positive material normalization joins the real/proper temporal weight
    // before products, preserving both extreme range and authored sigma.
    fn square_scaled(
        &self,
        field: impl Fn(usize) -> [Dd; 2],
        material_scale: Dd,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, MomentError> {
        if cancelled() {
            return Err(MomentError::Cancelled);
        }
        let mut amplitudes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
        let half = Dd::from_f64(0.5);
        let mut nonzero = false;
        for i in 0..self.frequencies {
            let [real, imaginary] = field(i);
            amplitudes[i] = [real, -imaginary];
            nonzero |= real != Dd::ZERO || imaginary != Dd::ZERO;
        }
        if !nonzero {
            return Ok(Dd::ZERO);
        }
        let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        let count = self.factor.mode_count();
        self.factor
            .project_parts_dd(&amplitudes[..self.frequencies], &mut modes[..count])
            .map_err(MomentError::Product)?;
        let scale = material_scale * if self.proper { half } else { Dd::ONE };
        let mut power = Dd::ZERO;
        for (mode, &field) in modes[..count].iter().enumerate() {
            power = power
                + self
                    .factor
                    .mode_cross_parts_dd_scaled(mode, field, field, scale)
                    .map_err(MomentError::Product)?[0];
        }
        if cancelled() {
            return Err(MomentError::Cancelled);
        }
        if !power.to_f64().is_finite() {
            return Err(MomentError::Unrepresentable);
        }
        Ok(power)
    }
}

// Filtering in DD retains the original product residual through a reservoir
// null. Rounding each complex field/filter product first would lose it before
// the positive temporal factor has a chance to resolve the surviving heat.
fn filtered_parts(field: C64, transfer: C64) -> [Dd; 2] {
    let (a, b) = (Dd::from_f64(field.re), Dd::from_f64(field.im));
    let (c, d) = (Dd::from_f64(transfer.re), Dd::from_f64(transfer.im));
    [a * c - b * d, a * d + b * c]
}

impl<'e, 'b, 's, 'm> PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
    fn new(
        anchor: &'b PlanarAngularSpectrum<'s, 'm>,
        mean: Option<&'e AngularSpectrumExposure<'b, 's, 'm>>,
        factors: &'e [AngularSpectrumExposure<'b, 's, 'm>],
        statistics: Option<&'e [AngularSpectrumFactorStatistics]>,
        layer_index: usize,
        limits: AngularSpectrumExposureLayerLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(cancelled)?;
        let selected = anchor
            .source()
            .layer(layer_index)
            .ok_or_else(|| bad("selected field layer is outside the original stack"))?;
        let law = LossLaw::new(selected.response)?;
        let inertial_scratch_bytes = if law.has_inertial_storage() {
            core::mem::size_of::<[C64; 3]>()
        } else {
            0
        };
        let n = anchor.modes.len();
        let mut groups = 0_usize;
        let mut point_line_pairs = 0_usize;
        let mut pairs = 0_usize;
        let mut carrier_pairs = 0_usize;
        let mut positive_pairs = 0_usize;
        let mut field_fft_work = 0_u64;
        let mut query_workspace = anchor.workspace_bytes();
        let extra_kernel_bytes = core::mem::size_of::<MomentKernel<HarmonicProduct>>()
            .max(core::mem::size_of::<MomentKernel<ProperHarmonicProduct>>())
            .saturating_sub(core::mem::size_of::<HarmonicProduct>());
        for (bank, proper) in ensemble::banks(mean, factors, statistics) {
            poll(cancelled)?;
            if bank.anchor.grid() != anchor.grid()
                || bank.anchor.carrier_rad_m != anchor.carrier_rad_m
                || !bank
                    .anchor
                    .source()
                    .same_sources_and_geometry(anchor.source())
            {
                return Err(bad(
                    "film exposure requires the same original materials, layers, grid and carrier",
                ));
            }
            groups = groups
                .checked_add(bank.frequency_count())
                .ok_or_else(|| bad("film exposure frequency count overflow"))?;
            point_line_pairs = n
                .checked_mul(bank.line_count())
                .and_then(|v| v.checked_add(point_line_pairs))
                .ok_or_else(|| bad("film exposure line work overflow"))?;
            let f = bank.frequency_count();
            let band_tones = positive_band_tones(f, proper);
            if let Some(tones) = band_tones {
                positive_pairs = n
                    .checked_mul(6 + law.pole_count() + 2 * law.reservoir_pair_channels())
                    .and_then(|v| v.checked_mul(tones + 1))
                    .and_then(|v| v.checked_mul(MAX_UNIFORM_HARMONIC_BAND_MODES))
                    .and_then(|v| v.checked_add(positive_pairs))
                    .ok_or_else(|| bad("film exposure positive moment work overflow"))?;
            }
            let bank_pairs = n
                .checked_mul(f)
                .and_then(|v| v.checked_mul(f))
                .ok_or_else(|| bad("film exposure pair work overflow"))?;
            pairs = pairs
                .checked_add(bank_pairs)
                .ok_or_else(|| bad("film exposure pair work overflow"))?;
            if !proper {
                carrier_pairs = bank
                    .point_frequency_pairs()
                    .checked_sub(bank_pairs)
                    .and_then(|v| carrier_pairs.checked_add(v))
                    .ok_or_else(|| bad("film exposure carrier work overflow"))?;
            }
            field_fft_work = field_fft_work
                .checked_add(bank.field_fft_work())
                .ok_or_else(|| bad("film exposure FFT work overflow"))?;
            let scratch = f
                .checked_mul(f)
                .and_then(|v| v.checked_mul(extra_kernel_bytes))
                .and_then(|v| {
                    v.checked_add(f.checked_mul(
                        core::mem::size_of::<usize>()
                            + core::mem::size_of::<[C64; 6]>()
                            + inertial_scratch_bytes,
                    )?)
                })
                .and_then(|v| v.checked_add(bank.workspace_bytes()))
                .ok_or_else(|| bad("film exposure workspace overflow"))?;
            query_workspace = query_workspace.max(scratch);
        }
        let point_reservoir_frequency_pairs = pairs
            .checked_mul(law.reservoir_pair_channels() + 3)
            .and_then(|v| v.checked_add(carrier_pairs))
            .and_then(|v| v.checked_add(positive_pairs))
            .ok_or_else(|| bad("film exposure reservoir work overflow"))?;
        let index_bytes = groups
            .checked_mul(core::mem::size_of::<ExposureLine<'_, '_, '_>>())
            .ok_or_else(|| bad("film exposure workspace overflow"))?;
        if groups > 4096
            || groups > limits.maximum_frequency_groups
            || point_line_pairs > limits.maximum_point_line_pairs
            || point_reservoir_frequency_pairs > limits.maximum_point_reservoir_frequency_pairs
            || field_fft_work > limits.maximum_fft_work
            || index_bytes
                .checked_add(query_workspace)
                .is_none_or(|v| v > limits.maximum_workspace_bytes)
        {
            return Err(bad("film exposure exceeds its work or workspace budget"));
        }
        let mut representatives: Vec<ExposureLine<'b, 's, 'm>> = Vec::new();
        representatives
            .try_reserve_exact(groups)
            .map_err(|_| bad("film exposure frequency allocation failed"))?;
        for bank in mean.into_iter().chain(factors.iter()) {
            for group in 0..bank.frequency_count() {
                poll(cancelled)?;
                representatives.push(
                    bank.lines
                        .get(bank.source_indices[bank.group_offsets[group]]),
                );
            }
        }
        // If both kinds share a frequency, keep a one-sided representative:
        // it needs a shared cache, whereas two-sided lines already own theirs.
        representatives.sort_unstable_by(|a, b| {
            a.beam()
                .source()
                .angular_frequency_rad_s()
                .total_cmp(&b.beam().source().angular_frequency_rad_s())
                .then_with(|| a.is_two_sided().cmp(&b.is_two_sided()))
        });
        representatives.dedup_by(|a, b| {
            a.beam().source().angular_frequency_rad_s()
                == b.beam().source().angular_frequency_rad_s()
        });
        poll(cancelled)?;
        let count = representatives.len();
        let one_sided_count = representatives
            .iter()
            .filter(|line| !line.is_two_sided())
            .count();
        let per_frequency_pairs = n
            .checked_mul(2)
            .and_then(|v| v.checked_mul(anchor.source().layer_count()))
            .ok_or_else(|| bad("film exposure preparation work overflow"))?;
        let preparation_layer_mode_pairs = one_sided_count
            .checked_mul(per_frequency_pairs)
            .ok_or_else(|| bad("film exposure preparation work overflow"))?;
        let shared_layer_bytes = n
            .checked_mul(core::mem::size_of::<
                super::super::layers::PreparedMode<'_, '_>,
            >())
            .and_then(|v| v.checked_mul(one_sided_count))
            .ok_or_else(|| bad("film exposure workspace overflow"))?;
        let cache_bytes = law
            .pole_count()
            .checked_mul(core::mem::size_of::<PoleTransfer>())
            .and_then(|v| v.checked_add(core::mem::size_of::<FrequencyLayer<'_, '_, '_>>()))
            .and_then(|v| v.checked_mul(count))
            .and_then(|v| v.checked_add(shared_layer_bytes))
            .ok_or_else(|| bad("film exposure workspace overflow"))?;
        // Conservative peak: shared cache, one original bank query, source
        // indices, fixed preparation scratch, all sums and both output planes.
        // Other borrowed bank owners and allocator overhead are excluded.
        let workspace_bytes = n
            .checked_mul(
                core::mem::size_of::<[Dd; 12]>()
                    + core::mem::size_of::<AngularSpectrumExposurePoint>()
                    + core::mem::size_of::<AngularSpectrumExposureLayerDissipation>(),
            )
            .and_then(|v| v.checked_add(cache_bytes))
            .and_then(|v| v.checked_add(query_workspace))
            .and_then(|v| v.checked_add(index_bytes))
            .and_then(|v| v.checked_add(FIELD_PREPARATION_WORKSPACE_BYTES))
            .ok_or_else(|| bad("film exposure workspace overflow"))?;
        if preparation_layer_mode_pairs > limits.maximum_frequency_layer_mode_pairs
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "film exposure exceeds its preparation or workspace budget",
            ));
        }
        let mut frequencies = Vec::new();
        frequencies
            .try_reserve_exact(count)
            .map_err(|_| bad("film exposure cache allocation failed"))?;
        for line in representatives {
            poll(cancelled)?;
            let beam = line.beam();
            let omega = beam.source().angular_frequency_rad_s();
            let mut poles = Vec::new();
            poles
                .try_reserve_exact(law.pole_count())
                .map_err(|_| bad("film exposure reservoir allocation failed"))?;
            for &pole in law.debye_poles {
                poles.push(PoleTransfer::new(pole, omega)?);
            }
            for &pole in law.lorentz_poles {
                poles.push(PoleTransfer::lorentz(pole, omega)?);
            }
            let layer = if line.is_two_sided() {
                None
            } else {
                Some(beam.prepare_layer(
                    layer_index,
                    AngularSpectrumLayerLimits {
                        maximum_layer_mode_pairs: per_frequency_pairs,
                        maximum_workspace_bytes: limits.maximum_workspace_bytes,
                        maximum_fft_work: limits.maximum_fft_work,
                    },
                    &mut *cancelled,
                )?)
            };
            frequencies.push(FrequencyLayer {
                omega,
                layer,
                poles,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            anchor,
            mean,
            factors,
            statistics,
            layer_index,
            thickness_m: selected.thickness_m,
            law,
            frequencies,
            workspace_bytes,
            preparation_layer_mode_pairs,
            point_line_pairs,
            point_reservoir_frequency_pairs,
            field_fft_work,
        })
    }

    /// Original anchor stack. Every admitted frequency retains these same
    /// material owners, film order and geometry; no detached effective medium
    /// is substituted when a downstream volume owner checks its binding.
    pub fn source(&self) -> &'s PlanarDielectricLayers<'m> {
        self.anchor.source()
    }

    pub fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub fn frequency_count(&self) -> usize {
        self.frequencies.len()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn preparation_layer_mode_pairs(&self) -> usize {
        self.preparation_layer_mode_pairs
    }
    pub fn point_line_pairs(&self) -> usize {
        self.point_line_pairs
    }
    pub fn point_reservoir_frequency_pairs(&self) -> usize {
        self.point_reservoir_frequency_pairs
    }
    pub fn field_fft_work(&self) -> u64 {
        self.field_fft_work
    }

    /// Mean irreversible heat and signed storage rate on an interior plane.
    /// Integrating heat PLUS storage rate through depth gives the net boundary
    /// power entering that film. Heat alone generally differs for a finite
    /// shutter. This does not assign switch-on state or thermal parameter work.
    pub fn field(
        &self,
        depth_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureLayerField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.thickness_m {
            return Err(bad("field depth is outside the selected layer"));
        }
        let n = self.anchor.modes.len();
        let mut sums = Vec::new();
        sums.try_reserve_exact(n)
            .map_err(|_| bad("film exposure sum allocation failed"))?;
        for point in 0..n {
            if point % 64 == 0 {
                poll(&mut cancelled)?;
            }
            sums.push([Dd::ZERO; 12]);
        }
        let mut normal_power = Dd::ZERO;
        for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
            if proper {
                normal_power = normal_power
                    + self.accumulate_bank(
                        bank,
                        depth_m,
                        true,
                        window,
                        |a, b| window.proper_product_kernel(a, b),
                        proper_product,
                        ProperHarmonicProduct::mean_real_product_parts_dd,
                        &mut sums,
                        &mut cancelled,
                    )?;
            } else {
                normal_power = normal_power
                    + self.accumulate_bank(
                        bank,
                        depth_m,
                        false,
                        window,
                        |a, b| window.product_kernel(a, b),
                        product,
                        HarmonicProduct::mean_real_product_parts_dd,
                        &mut sums,
                        &mut cancelled,
                    )?;
            }
        }
        let mut points = Vec::new();
        let mut dissipation = Vec::new();
        points
            .try_reserve_exact(n)
            .map_err(|_| bad("film exposure field allocation failed"))?;
        dissipation
            .try_reserve_exact(n)
            .map_err(|_| bad("film exposure heat allocation failed"))?;
        let mut totals = [Dd::ZERO; 2];
        for (point, sum) in sums.into_iter().enumerate() {
            if point % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let value = sum.map(Dd::to_f64);
            let heat = (sum[9] + sum[10]).to_f64();
            if value.iter().any(|v| !v.is_finite())
                || !heat.is_finite()
                || value[..6].iter().any(|v| *v < 0.0)
                || value[9] < 0.0
                || value[10] < 0.0
            {
                return Err(bad("film exposure heat or storage rate is unrepresentable"));
            }
            points.push(AngularSpectrumExposurePoint {
                electric_mean_squared_v2_m2: [value[0], value[1], value[2]],
                magnetic_mean_squared_a2_m2: [value[3], value[4], value[5]],
                poynting_w_m2: [value[6], value[7], value[8]],
            });
            dissipation.push(AngularSpectrumExposureLayerDissipation {
                ohmic_heat_w_m3: value[9],
                polarization_heat_w_m3: value[10],
                total_heat_w_m3: heat,
                storage_rate_w_m3: value[11],
            });
            totals[0] = totals[0] + sum[9] + sum[10];
            totals[1] = totals[1] + sum[11];
        }
        let grid = self.anchor.grid();
        let totals = totals.map(|v| v.to_f64() * grid.spacing_m[0] * grid.spacing_m[1]);
        let normal_power_w = normal_power.to_f64();
        if !normal_power_w.is_finite() || totals.iter().any(|v| !v.is_finite()) {
            return Err(bad("integrated film exposure is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(AngularSpectrumExposureLayerField {
            field: AngularSpectrumExposureField {
                grid,
                points,
                normal_power_w,
            },
            dissipation,
            absorbed_power_per_depth_w_m: totals[0],
            storage_rate_per_depth_w_m: totals[1],
        })
    }

    fn accumulate_bank<P: Copy>(
        &self,
        bank: &AngularSpectrumExposure<'b, 's, 'm>,
        depth: f64,
        proper: bool,
        window: HarmonicWindow,
        make_kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
        product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
        product_parts: impl Fn(P, [Dd; 2], [Dd; 2]) -> Result<Dd, fs_fft::harmonic_window::WindowError>,
        sums: &mut [[Dd; 12]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, AngularSpectrumError> {
        poll(cancelled)?;
        let f = bank.frequency_count();
        if f == 0 {
            return Ok(Dd::ZERO);
        }
        // The factor is reused across every spatial point and original loss
        // reservoir. Its bounded stack workspace adds no query allocation.
        let mut positive_workspace =
            positive_band_tones(f, proper).map(|_| HarmonicBandWorkspace::new());
        let positive = if let Some(workspace) = positive_workspace.as_mut() {
            PositiveLayerMoments::prepare(window, &bank.frequencies, proper, workspace, cancelled)?
        } else {
            None
        };
        let mut kernels = Vec::new();
        kernels
            .try_reserve_exact(f * f)
            .map_err(|_| bad("film exposure kernel allocation failed"))?;
        for &a in &bank.frequencies {
            for &b in &bank.frequencies {
                poll(cancelled)?;
                let storage_frequency = if proper {
                    (Dd::from_f64(a) - Dd::from_f64(b)) * Dd::from_f64(0.5)
                } else {
                    Dd::from_f64(a)
                };
                kernels.push(MomentKernel {
                    product: make_kernel(a, b)
                        .map_err(|_| bad("film exposure temporal kernel is unrepresentable"))?,
                    storage_frequency,
                });
            }
        }
        let mut mapping = Vec::new();
        let mut fields = Vec::new();
        let mut filtered = Vec::new();
        let mut kinetic = Vec::new();
        mapping
            .try_reserve_exact(f)
            .map_err(|_| bad("film exposure index allocation failed"))?;
        fields
            .try_reserve_exact(f)
            .map_err(|_| bad("film exposure field allocation failed"))?;
        if positive.is_none() {
            filtered
                .try_reserve_exact(f)
                .map_err(|_| bad("film exposure reservoir field allocation failed"))?;
            if self.law.has_inertial_storage() {
                kinetic
                    .try_reserve_exact(f)
                    .map_err(|_| bad("film exposure inertial field allocation failed"))?;
            }
        }
        for group in 0..f {
            poll(cancelled)?;
            let index = self
                .frequencies
                .binary_search_by(|cache| cache.omega.total_cmp(&bank.frequencies[group]))
                .map_err(|_| bad("original film frequency cache is missing"))?;
            let layer = &self.frequencies[index].layer;
            fields.push(
                bank.frequency_modal_field_with(group, cancelled, |line, cancelled| {
                    line.modal_layer_field(layer.as_ref(), self.layer_index, depth, cancelled)
                })?,
            );
            mapping.push(index);
            if positive.is_none() {
                filtered.push([C64::ZERO; 6]);
                if self.law.has_inertial_storage() {
                    kinetic.push([C64::ZERO; 3]);
                }
            }
        }
        let normal_power = {
            let carrier = bank.carrier_fields(&fields, proper, cancelled)?;
            integrate_cell_power(
                self.anchor.grid(),
                &fields,
                carrier.as_ref(),
                &kernels,
                |kernel, a, b| product_parts(kernel.product, a, b),
                cancelled,
            )?
        };
        for (group, field) in fields.iter_mut().enumerate() {
            bank.synthesize_frequency_field(group, field, cancelled)?;
        }
        for (point, sum) in sums.iter_mut().enumerate() {
            poll(cancelled)?;
            let projected = positive
                .as_ref()
                .map(|positive| positive.project_field(|i| fields[i][point], cancelled))
                .transpose()
                .map_err(moment_error)?;
            let moments = if let (Some(positive), Some(projected)) = (&positive, &projected) {
                projected
                    .moments(positive, cancelled)
                    .map_err(moment_error)?
            } else {
                layer_point_moments(
                    f,
                    &kernels,
                    |i| fields[i][point],
                    |k, a, b| product(k.product, a, b),
                    None,
                    cancelled,
                )
                .map_err(moment_error)?
            };
            let mut storage = if let (Some(positive), Some(projected)) = (&positive, &projected) {
                projected
                    .storage_rate(
                        positive,
                        self.law.epsilon_infinity,
                        self.law.permeability,
                        cancelled,
                    )
                    .map_err(moment_error)?
            } else {
                instantaneous_storage_rate(
                    &kernels,
                    f,
                    |i| fields[i][point],
                    &product,
                    self.law.epsilon_infinity,
                    self.law.permeability,
                    cancelled,
                )
                .map_err(moment_error)?
            };
            let mut polarization_heat = Dd::ZERO;
            for pole in 0..self.law.pole_count() {
                poll(cancelled)?;
                let (heat, rate) = if let Some(positive) = &positive {
                    positive
                        .reservoir_moments(
                            |i| {
                                [
                                    fields[i][point][0],
                                    fields[i][point][1],
                                    fields[i][point][2],
                                ]
                            },
                            |i| self.frequencies[mapping[i]].poles[pole],
                            pole >= self.law.debye_poles.len(),
                            cancelled,
                        )
                        .map_err(moment_error)?
                } else {
                    for i in 0..f {
                        if i % 64 == 0 {
                            poll(cancelled)?;
                        }
                        let transfer = self.frequencies[mapping[i]].poles[pole];
                        for c in 0..3 {
                            filtered[i][c] = fields[i][point][c] * transfer.heat;
                            filtered[i][c + 3] = fields[i][point][c] * transfer.storage;
                            if !kinetic.is_empty() {
                                kinetic[i][c] = fields[i][point][c] * transfer.kinetic_storage;
                            }
                        }
                    }
                    pole_moments_with_heat(
                        &kernels,
                        &filtered,
                        (!kinetic.is_empty()).then_some(kinetic.as_slice()),
                        None,
                        &product,
                        cancelled,
                    )
                    .map_err(moment_error)?
                };
                polarization_heat = polarization_heat + heat;
                storage = storage + rate;
            }
            for c in 0..9 {
                sum[c] = sum[c] + moments[c];
            }
            sum[9] = sum[9]
                + Dd::from_f64(self.law.conductivity) * (moments[0] + moments[1] + moments[2]);
            sum[10] = sum[10] + polarization_heat;
            sum[11] = sum[11] + storage;
        }
        Ok(normal_power)
    }
}

fn moment_error(error: MomentError) -> AngularSpectrumError {
    match error {
        MomentError::Cancelled => AngularSpectrumError::Cancelled,
        MomentError::Product(_) | MomentError::Unrepresentable => {
            bad("film reservoir product is unrepresentable")
        }
    }
}

fn rotate_i(z: C64) -> C64 {
    C64::new(-z.im, z.re)
}

// The same projected field supplies both its diagonal populations and signed
// Poynting components. A cross moment can be negative, including real backflow.
fn layer_point_moments<P: Copy>(
    frequencies: usize,
    kernels: &[P],
    field_at: impl Fn(usize) -> [C64; 6],
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    positive: Option<&PositiveLayerMoments<'_>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 9], MomentError> {
    let Some(positive) = positive else {
        return integrate_point_dd(frequencies, kernels, field_at, product, cancelled);
    };
    positive
        .project_field(field_at, cancelled)?
        .moments(positive, cancelled)
}

fn instantaneous_storage_rate<P: Copy>(
    kernels: &[MomentKernel<P>],
    f: usize,
    field: impl Fn(usize) -> [C64; 6],
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    epsilon: f64,
    mu: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, MomentError> {
    let mut electric = Dd::ZERO;
    let mut magnetic = Dd::ZERO;
    for i in 0..f {
        let a = field(i);
        for j in 0..f {
            if j % 64 == 0 && cancelled() {
                return Err(MomentError::Cancelled);
            }
            let b = field(j);
            let kernel = kernels[i * f + j];
            for c in 0..3 {
                electric = electric
                    + product(kernel.product, rotate_i(a[c]), b[c])? * kernel.storage_frequency;
                magnetic = magnetic
                    + product(kernel.product, rotate_i(a[c + 3]), b[c + 3])?
                        * kernel.storage_frequency;
            }
        }
    }
    Ok(electric * Dd::from_f64(epsilon) + magnetic * Dd::from_f64(mu))
}

#[cfg(test)]
fn pole_moments<P: Copy>(
    kernels: &[MomentKernel<P>],
    fields: &[[C64; 6]],
    kinetic: Option<&[[C64; 3]]>,
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Dd, Dd), MomentError> {
    pole_moments_with_heat(kernels, fields, kinetic, None, product, cancelled)
}

fn pole_moments_with_heat<P: Copy>(
    kernels: &[MomentKernel<P>],
    fields: &[[C64; 6]],
    kinetic: Option<&[[C64; 3]]>,
    positive_heat: Option<Dd>,
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Dd, Dd), MomentError> {
    let f = fields.len();
    if kinetic.is_some_and(|values| values.len() != f) {
        return Err(MomentError::Unrepresentable);
    }
    let mut heat = positive_heat.unwrap_or(Dd::ZERO);
    let mut storage = Dd::ZERO;
    for (i, a) in fields.iter().enumerate() {
        for (j, b) in fields.iter().enumerate() {
            if j % 64 == 0 && cancelled() {
                return Err(MomentError::Cancelled);
            }
            let kernel = kernels[i * f + j];
            for c in 0..3 {
                if positive_heat.is_none() {
                    heat = heat + product(kernel.product, a[c], b[c])?;
                }
                storage = storage
                    + product(kernel.product, rotate_i(a[c + 3]), b[c + 3])?
                        * kernel.storage_frequency;
                if let Some(values) = kinetic {
                    storage = storage
                        + product(kernel.product, rotate_i(values[i][c]), values[j][c])?
                            * kernel.storage_frequency;
                }
            }
        }
    }
    if !heat.to_f64().is_finite() || heat.to_f64() < 0.0 || !storage.to_f64().is_finite() {
        return Err(MomentError::Unrepresentable);
    }
    Ok((heat, storage))
}

#[cfg(test)]
mod lorentz_tests;

#[cfg(test)]
mod short_window_tests;

#[cfg(test)]
mod temporal_moments_tests;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn proper_slow_beat_storage_ignores_a_large_common_carrier() {
        let window = HarmonicWindow::new(0.17, 0.12).unwrap();
        let delta = 0.25;
        let half_width = 0.5 * delta * window.duration_s();
        // E_x = H_y = exp(i wa t) + i exp(i wb t), epsilon=mu=1.
        // Proper ensemble energy is 1 - sin(delta*t), so its rate is
        // -delta*cos(delta*t), without any common-carrier contribution.
        let expected = -delta
            * (delta * (window.start_s() + 0.5 * window.duration_s())).cos()
            * half_width.sin()
            / half_width;
        let fields = [
            [
                C64::ONE,
                C64::ZERO,
                C64::ZERO,
                C64::ZERO,
                C64::ONE,
                C64::ZERO,
            ],
            [
                C64::new(0.0, 1.0),
                C64::ZERO,
                C64::ZERO,
                C64::ZERO,
                C64::new(0.0, 1.0),
                C64::ZERO,
            ],
        ];
        let mut observed = Vec::new();
        for carrier in [0.0, 1_125_899_906_842_624.0] {
            let frequencies = [carrier + 2.0, carrier + 2.25];
            let mut kernels = Vec::new();
            for a in frequencies {
                for b in frequencies {
                    kernels.push(MomentKernel {
                        product: window.proper_product_kernel(a, b).unwrap(),
                        storage_frequency: (Dd::from_f64(a) - Dd::from_f64(b)) * Dd::from_f64(0.5),
                    });
                }
            }
            let actual = instantaneous_storage_rate(
                &kernels,
                2,
                |i| fields[i],
                proper_product,
                1.0,
                1.0,
                &mut || false,
            )
            .unwrap()
            .to_f64();
            assert!((actual - expected).abs() < 5e-15, "{actual} != {expected}");
            observed.push(actual);
        }
        assert_eq!(observed[0], observed[1]);
    }

    #[test]
    fn real_polarization_current_heat_and_storage_satisfy_the_reservoir_ode() {
        use fs_material::field_storage::dielectric::DebyePole;
        let pole = DebyePole {
            permittivity_increment_f_m: 0.7,
            relaxation_time_s: 0.4,
            capacitance_f: 1.0,
            resistance_ohm: 0.4,
        };
        let frequencies = [2.0, 2.3];
        let electric = [C64::new(0.8, 0.2), C64::new(-0.3, 0.6)];
        let window = HarmonicWindow::new(0.17, 0.12).unwrap();
        let mut kernels = Vec::new();
        for a in frequencies {
            for b in frequencies {
                kernels.push(MomentKernel {
                    product: window.product_kernel(a, b).unwrap(),
                    storage_frequency: Dd::from_f64(a),
                });
            }
        }
        let filtered = frequencies
            .iter()
            .zip(electric)
            .map(|(&omega, e)| {
                let transfer = PoleTransfer::new(pole, omega).unwrap();
                [
                    e * transfer.heat,
                    C64::ZERO,
                    C64::ZERO,
                    e * transfer.storage,
                    C64::ZERO,
                    C64::ZERO,
                ]
            })
            .collect::<Vec<_>>();
        let (heat, storage) =
            pole_moments(&kernels, &filtered, None, product, &mut || false).unwrap();
        let mut electrical_work = Dd::ZERO;
        for i in 0..2 {
            for j in 0..2 {
                // Independent original ODE: P = Delta E/(1+i omega tau).
                let polarization = electric[j].scale(0.7) / C64::new(1.0, frequencies[j] * 0.4);
                let current = rotate_i(polarization).scale(frequencies[j]);
                electrical_work = electrical_work
                    + product(kernels[i * 2 + j].product, electric[i], current).unwrap();
            }
        }
        assert!(heat.to_f64() > 0.0);
        assert!((heat + storage - electrical_work).to_f64().abs() < 2e-14);
    }
}

