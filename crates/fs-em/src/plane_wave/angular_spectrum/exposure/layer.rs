//! Finite-window local film heat and reversible field/polarization storage.
//!
//! These are settled prescribed harmonics through the original film operator,
//! not a switch-on transient. Ohmic and each original Debye reservoir dissipate
//! after all spatial and temporal interference. Boundary power loss also
//! contains the signed change of electromagnetic/polarization storage.

use super::*;
use crate::plane_wave::layers::{FIELD_PREPARATION_WORKSPACE_BYTES, PreparedLayerField};

mod reservoir;
use reservoir::{LossLaw, PoleTransfer};

mod integrated_heat;
pub use integrated_heat::{
    AngularSpectrumExposureLayerHeatLoadSpectrum, AngularSpectrumHeatAreaQuadrature,
    AngularSpectrumHeatDepthSample, AngularSpectrumHeatVolumeQuadrature,
    AngularSpectrumIntegratedHeatLimits,
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
    /// N * sum(F_r²) * (3 + P), with P original Debye poles. The three base
    /// channels cover the vector field and its instantaneous E/H stores.
    pub maximum_point_reservoir_frequency_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// Mean local irreversible heat and signed reversible storage rate [W/m³].
/// Storage includes epsilon_infinity |E|²/2, mu |H|²/2 and each
/// |P_p|²/(2 Delta_epsilon_p). It excludes thermal energy and parameter work.
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
        let n = anchor.modes.len();
        let mut groups = 0_usize;
        let mut point_line_pairs = 0_usize;
        let mut pairs = 0_usize;
        let mut field_fft_work = 0_u64;
        let mut query_workspace = anchor.workspace_bytes();
        let extra_kernel_bytes = core::mem::size_of::<MomentKernel<HarmonicProduct>>()
            .max(core::mem::size_of::<MomentKernel<ProperHarmonicProduct>>())
            .saturating_sub(core::mem::size_of::<HarmonicProduct>());
        for bank in mean.into_iter().chain(factors.iter()) {
            poll(cancelled)?;
            if bank.anchor.grid() != anchor.grid()
                || !bank
                    .anchor
                    .source()
                    .same_sources_and_geometry(anchor.source())
            {
                return Err(bad(
                    "film exposure requires the same original materials, layers and grid",
                ));
            }
            groups = groups
                .checked_add(bank.frequency_count())
                .ok_or_else(|| bad("film exposure frequency count overflow"))?;
            point_line_pairs = n
                .checked_mul(bank.line_count())
                .and_then(|v| v.checked_add(point_line_pairs))
                .ok_or_else(|| bad("film exposure line work overflow"))?;
            pairs = pairs
                .checked_add(bank.point_frequency_pairs())
                .ok_or_else(|| bad("film exposure pair work overflow"))?;
            field_fft_work = field_fft_work
                .checked_add(bank.field_fft_work())
                .ok_or_else(|| bad("film exposure FFT work overflow"))?;
            let f = bank.frequency_count();
            let scratch = f
                .checked_mul(f)
                .and_then(|v| v.checked_mul(extra_kernel_bytes))
                .and_then(|v| {
                    v.checked_add(f.checked_mul(
                        core::mem::size_of::<usize>() + core::mem::size_of::<[C64; 6]>(),
                    )?)
                })
                .and_then(|v| v.checked_add(bank.workspace_bytes()))
                .ok_or_else(|| bad("film exposure workspace overflow"))?;
            query_workspace = query_workspace.max(scratch);
        }
        let point_reservoir_frequency_pairs = pairs
            .checked_mul(law.poles.len() + 3)
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
            .checked_mul(core::mem::size_of::<[PreparedLayerField<'_, '_>; 2]>())
            .and_then(|v| v.checked_mul(one_sided_count))
            .ok_or_else(|| bad("film exposure workspace overflow"))?;
        let cache_bytes = law
            .poles
            .len()
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
                .try_reserve_exact(law.poles.len())
                .map_err(|_| bad("film exposure reservoir allocation failed"))?;
            for &pole in law.poles {
                poles.push(PoleTransfer::new(pole, omega)?);
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
        mapping
            .try_reserve_exact(f)
            .map_err(|_| bad("film exposure index allocation failed"))?;
        fields
            .try_reserve_exact(f)
            .map_err(|_| bad("film exposure field allocation failed"))?;
        filtered
            .try_reserve_exact(f)
            .map_err(|_| bad("film exposure reservoir field allocation failed"))?;
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
            filtered.push([C64::ZERO; 6]);
        }
        let normal_power = integrate_cell_power(
            self.anchor.grid(),
            &fields,
            &kernels,
            |kernel, a, b| product_parts(kernel.product, a, b),
            cancelled,
        )?;
        for (group, field) in fields.iter_mut().enumerate() {
            bank.synthesize_frequency_field(group, field, cancelled)?;
        }
        for (point, sum) in sums.iter_mut().enumerate() {
            poll(cancelled)?;
            let moments = integrate_point_dd(
                f,
                &kernels,
                |i| fields[i][point],
                |k, a, b| product(k.product, a, b),
                cancelled,
            )
            .map_err(moment_error)?;
            let mut storage = instantaneous_storage_rate(
                &kernels,
                f,
                |i| fields[i][point],
                &product,
                self.law.epsilon_infinity,
                self.law.permeability,
                cancelled,
            )
            .map_err(moment_error)?;
            let mut polarization_heat = Dd::ZERO;
            for pole in 0..self.law.poles.len() {
                poll(cancelled)?;
                for i in 0..f {
                    if i % 64 == 0 {
                        poll(cancelled)?;
                    }
                    let transfer = self.frequencies[mapping[i]].poles[pole];
                    for c in 0..3 {
                        filtered[i][c] = fields[i][point][c] * transfer.heat;
                        filtered[i][c + 3] = fields[i][point][c] * transfer.storage;
                    }
                }
                let (heat, rate) =
                    pole_moments(&kernels, &filtered, &product, cancelled).map_err(moment_error)?;
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

fn pole_moments<P: Copy>(
    kernels: &[MomentKernel<P>],
    fields: &[[C64; 6]],
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Dd, Dd), MomentError> {
    let f = fields.len();
    let mut heat = Dd::ZERO;
    let mut storage = Dd::ZERO;
    for (i, a) in fields.iter().enumerate() {
        for (j, b) in fields.iter().enumerate() {
            if j % 64 == 0 && cancelled() {
                return Err(MomentError::Cancelled);
            }
            let kernel = kernels[i * f + j];
            for c in 0..3 {
                heat = heat + product(kernel.product, a[c], b[c])?;
                storage = storage
                    + product(kernel.product, rotate_i(a[c + 3]), b[c + 3])?
                        * kernel.storage_frequency;
            }
        }
    }
    if !heat.to_f64().is_finite() || heat.to_f64() < 0.0 || !storage.to_f64().is_finite() {
        return Err(MomentError::Unrepresentable);
    }
    Ok((heat, storage))
}

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
        let (heat, storage) = pole_moments(&kernels, &filtered, product, &mut || false).unwrap();
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
