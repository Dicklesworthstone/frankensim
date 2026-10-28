//! Material propagation on each original Rayleigh surface-to-receiver path.
//! Geometry and signed source weights are never collapsed to a mean distance.
use super::{BaffledPressure, MAX_HISTORY, RayleighMedium, SurfaceSample, path};
use crate::air_path::{GasAbsorptionModel, ObserverFilter, ObserverFilterReceipt};
use fs_material::{gas::GasState, relaxation::RelaxationState};
use fs_math::det;
use std::{collections::BTreeMap, f64::consts::TAU};

/// Complete cold-work and retained-response bounds, in addition to the
/// original receiver's geometry, input-count and causal-history limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GasPropagationLimits {
    /// Distinct exact path lengths requiring a source-law filter design.
    pub max_unique_paths: usize,
    /// Sum of retained source-filter taps across distinct paths.
    pub max_filter_taps: usize,
    /// Scalar products in path convolution and signed modal accumulation.
    pub max_kernel_terms: usize,
}
impl GasPropagationLimits {
    /// Explicit cold-preparation budget for a finite surface quadrature.
    /// This bounds construction and retained responses; it does not promise
    /// real-time construction or playback. Inspect the admitted receiver's
    /// actual `multiply_adds_per_sample` for the selected geometry and gas.
    pub const fn surface_reference() -> Self {
        Self {
            max_unique_paths: 16_384,
            max_filter_taps: 4_000_000,
            max_kernel_terms: 500_000_000,
        }
    }
}

/// Explicit constitutive choice. Generic classical transport never acquires
/// terrestrial-air or molecular-relaxation data from a material name.
#[derive(Clone, Debug)]
pub struct GasPropagation {
    pub law: GasAbsorptionModel,
    /// Required only by the molecular law; bound to the exact same gas state.
    pub relaxation: Option<RelaxationState>,
    pub limits: GasPropagationLimits,
}
impl GasPropagation {
    /// Classical viscothermal propagation from rho, c, mu, k and cp.
    /// The original filter receipt explicitly retains omitted relaxation.
    pub const fn classical(limits: GasPropagationLimits) -> Self {
        Self {
            law: GasAbsorptionModel::StokesKirchhoff,
            relaxation: None,
            limits,
        }
    }

    /// Shared native record tail: LAW MAX_PATHS MAX_FILTER_TAPS MAX_KERNEL_TERMS.
    /// Classical transport uses the complete bound gas. Lossless is an explicit
    /// comparison; species relaxation data and terrestrial-air eligibility are
    /// never inferred from a name or from temperature alone.
    pub fn read_record(fields: &mut std::str::SplitAsciiWhitespace<'_>)
        -> Result<Self, &'static str>
    {
        let law = match fields.next() {
            Some("stokes-kirchhoff") => GasAbsorptionModel::StokesKirchhoff,
            Some("lossless") => GasAbsorptionModel::Lossless,
            _ => return Err("observer_gas requires stokes-kirchhoff or lossless and three work budgets"),
        };
        let ceiling = GasPropagationLimits::surface_reference();
        let mut count = |maximum| fields.next().and_then(|s| s.parse::<usize>().ok())
            .filter(|&n| n > 0 && n <= maximum)
            .ok_or("observer_gas work budget is outside its positive shared ceiling");
        let limits = GasPropagationLimits {
            max_unique_paths: count(ceiling.max_unique_paths)?,
            max_filter_taps: count(ceiling.max_filter_taps)?,
            max_kernel_terms: count(ceiling.max_kernel_terms)?,
        };
        if fields.next().is_some() { return Err("unexpected observer_gas fields"); }
        Ok(Self { law, relaxation: None, limits })
    }
}

/// Original physical law and each source-derived path filter, before the
/// common numerical sinc/padding. Ordering follows first surface occurrence.
#[derive(Clone, Debug)]
pub struct GasPropagationReport {
    pub gas: GasState,
    pub law: GasAbsorptionModel,
    /// One index per original surface sample into `filters`; only BIT-EQUAL
    /// computed distances share a filter. No distance binning or tolerance.
    pub surface_filter_indices: Vec<usize>,
    pub filters: Vec<ObserverFilterReceipt>,
    /// Additional whole-sample padding needed to make ALL residual delays
    /// causal. It is common to every path, preserving their relative phase.
    /// The receiver includes it in `filter_delay_samples`, exactly once.
    pub additional_latency_samples: usize,
    pub kernel_terms: usize,
}

struct PreparedPath {
    response: Vec<f64>,
    receipt: ObserverFilterReceipt,
    residual_delay_samples: f64,
    first: usize,
    kernel: Vec<f64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), String> {
    if cancelled() {
        Err("material Rayleigh propagation cancelled".into())
    } else {
        Ok(())
    }
}

impl BaffledPressure {
    /// Align separately prepared physical apertures to one numerical padding.
    /// Called only on cold material-path receivers before publication. Relative
    /// physical travel and gas dispersion stay intact; the same reported latency
    /// then accompanies the mouth and terminals when their pressure is summed.
    pub(crate) fn align_gas_latency(
        &mut self, report: &mut GasPropagationReport, latency: usize,
    ) -> Result<(), String> {
        let extra = latency.checked_sub(report.additional_latency_samples)
            .ok_or("material receiver alignment cannot remove causal padding")?;
        if extra == 0 { return Ok(()); }
        let last = self.delay_samples.1.checked_add(extra)
            .filter(|&last| last < MAX_HISTORY)
            .ok_or("aligned material receiver exceeds its causal-history budget")?;
        let frames = last + 1;
        let samples = frames.checked_mul(self.modes).ok_or("aligned receiver history overflow")?;
        self.history = vec![0.0; samples];
        self.history_frames = frames;
        for (delay, _) in &mut self.kernels { *delay += extra; }
        self.delay_samples = (self.delay_samples.0 + extra, last);
        self.material_path_delay += extra;
        report.additional_latency_samples = latency;
        self.prepare_taps();
        Ok(())
    }

    /// Shared bandlimited Rayleigh observer with source-gas propagation on
    /// EVERY physical quadrature path. The original ObserverFilter owns all
    /// absorption, dispersion, source admission and response-tail certificates.
    ///
    /// For path r, its filter is convolved with the residual travel delay
    /// `r/c - filter.low_frequency_delay_s`. If a path needs negative residual
    /// delay, ALL paths receive the same explicitly reported causal padding.
    /// Thus each path's low-frequency delay is r/c plus the common numerical
    /// latency, up to the existing sinc interpolation error. No atmospheric
    /// filter is applied to the already summed pressure or counted twice.
    ///
    /// Radiation still uses rho*dA/(2*pi*r), signed original mode shapes and
    /// the fixed infinite-baffle model. This adds no mechanical reaction, room
    /// scattering, inferred molecular data or continuum-accuracy certificate.
    /// Complete preparation is cold and atomic; sample stepping uses the SAME
    /// original allocation-free receiver. Cancellation polls between bounded
    /// filter designs, path convolutions and surface tiles.
    #[allow(clippy::too_many_arguments)]
    pub fn new_multirate_gas_differentiating(
        surface: &[SurfaceSample],
        inputs: usize,
        mechanics_rate: u32,
        output_rate: u32,
        position_m: [f64; 3],
        gas: &GasState,
        propagation: &GasPropagation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, GasPropagationReport), String> {
        poll(cancelled)?;
        gas.validate_coherence()
            .map_err(|error| format!("Rayleigh propagation gas: {error}"))?;
        let limits = propagation.limits;
        if limits.max_unique_paths == 0
            || limits.max_filter_taps == 0
            || limits.max_kernel_terms == 0
        {
            return Err("material Rayleigh propagation needs positive finite work budgets".into());
        }
        if propagation.relaxation.is_some()
            && propagation.law != GasAbsorptionModel::MolecularRelaxation
        {
            return Err(
                "a supplied gas relaxation state requires the molecular propagation law".into(),
            );
        }
        // Reuse complete original geometry, clocks and allocation admission.
        let mut receiver = Self::new_multirate(
            surface,
            inputs,
            mechanics_rate,
            output_rate,
            position_m,
            RayleighMedium::from(gas),
        )?;
        let rate = f64::from(output_rate);
        let mut by_distance = BTreeMap::<u64, usize>::new();
        let mut indices = Vec::with_capacity(surface.len());
        let mut distances = Vec::new();
        for sample in surface {
            poll(cancelled)?;
            let distance = det::sqrt(
                (0..3)
                    .map(|i| (position_m[i] - sample.position_m[i]).powi(2))
                    .sum(),
            );
            let index = if let Some(&index) = by_distance.get(&distance.to_bits()) {
                index
            } else {
                if distances.len() == limits.max_unique_paths {
                    return Err("Rayleigh surface exceeds its declared distinct-path budget".into());
                }
                let index = distances.len();
                distances.push(distance);
                by_distance.insert(distance.to_bits(), index);
                index
            };
            indices.push(index);
        }
        let mut prepared = Vec::with_capacity(distances.len());
        let mut total_filter_taps = 0usize;
        let mut minimum_delay = 0.0_f64;
        for distance in distances {
            poll(cancelled)?;
            let filter = ObserverFilter::for_law(
                propagation.law,
                gas,
                propagation.relaxation.as_ref(),
                distance,
                output_rate,
            )
            .map_err(|error| format!("Rayleigh gas path at {distance} m: {error}"))?;
            total_filter_taps = total_filter_taps
                .checked_add(filter.taps().len())
                .ok_or("Rayleigh source response extent overflow")?;
            if total_filter_taps > limits.max_filter_taps {
                return Err(
                    "Rayleigh gas paths exceed their declared source-response budget".into(),
                );
            }
            let receipt = *filter.receipt();
            let delay = rate * receipt.receiver_delay_s;
            if !delay.is_finite() || !receipt.low_frequency_delay_s.is_finite() {
                return Err("Rayleigh gas path delay is unrepresentable".into());
            }
            minimum_delay = minimum_delay.min(delay);
            prepared.push(PreparedPath {
                response: filter.taps().to_vec(),
                receipt,
                residual_delay_samples: delay,
                first: 0,
                kernel: Vec::new(),
            });
        }
        let padding = (-minimum_delay).ceil();
        if !padding.is_finite() || padding >= MAX_HISTORY as f64 {
            return Err("Rayleigh gas paths exceed the shared causal-padding budget".into());
        }
        let padding = padding as usize;
        let mut remaining_terms = limits.max_kernel_terms;
        for prepared in &mut prepared {
            poll(cancelled)?;
            (prepared.first, prepared.kernel) = path::compose(
                &prepared.response,
                prepared.residual_delay_samples + padding as f64,
                MAX_HISTORY,
                &mut remaining_terms,
            )?;
            // The convolved path now owns this response. Release its source
            // taps while the other exact-distance kernels are prepared.
            prepared.response = Vec::new();
        }
        let mut bins = BTreeMap::<usize, Vec<f64>>::new();
        for (sample, &index) in surface.iter().zip(&indices) {
            poll(cancelled)?;
            let prepared = &prepared[index];
            let terms = prepared
                .kernel
                .len()
                .checked_mul(inputs)
                .ok_or("Rayleigh material projection work overflow")?;
            if terms > remaining_terms {
                return Err("Rayleigh material projection exceeds its declared kernel work".into());
            }
            remaining_terms -= terms;
            let gain = gas.density * sample.area_m2 / (TAU * prepared.receipt.range_m);
            for (tap, &weight) in prepared.kernel.iter().enumerate() {
                let row = bins
                    .entry(prepared.first + tap)
                    .or_insert_with(|| vec![0.0; inputs]);
                for (value, &shape) in row.iter_mut().zip(&sample.mode_shape) {
                    *value += gain * weight * shape;
                }
            }
        }
        if bins.values().flatten().any(|value| !value.is_finite()) {
            return Err("material Rayleigh kernel is unrepresentable".into());
        }
        let first = *bins
            .first_key_value()
            .ok_or("empty material Rayleigh kernel")?
            .0;
        let last = *bins
            .last_key_value()
            .ok_or("empty material Rayleigh kernel")?
            .0;
        let frames = last + 1;
        let samples = frames
            .checked_mul(inputs)
            .ok_or("material Rayleigh history overflow")?;
        receiver.kernels = bins.into_iter().collect();
        receiver.history = vec![0.0; samples];
        receiver.history_frames = frames;
        receiver.delay_samples = (first, last);
        receiver.interpolation_delay = path::HALF;
        receiver.material_path_delay = padding;
        receiver.differentiate_velocity = true;
        receiver.set_input_ratio((mechanics_rate / output_rate) as usize)?;
        receiver.prepare_taps();
        poll(cancelled)?;
        Ok((
            receiver,
            GasPropagationReport {
                gas: *gas,
                law: propagation.law,
                surface_filter_indices: indices,
                filters: prepared.into_iter().map(|path| path.receipt).collect(),
                additional_latency_samples: padding,
                kernel_terms: limits.max_kernel_terms - remaining_terms,
            },
        ))
    }
}

#[cfg(test)]
#[path = "baffled_pressure_gas_tests.rs"]
mod tests;
