//! Exact per-tile volume flows through the original bounded Rayleigh receiver.
//! Partitioning the input basis changes only summation grouping. Every tile
//! retains its area, position, flow sign, travel delay and selected gas law.
//! This is the same stationary infinite-baffle observation, not a new exterior.

use super::{BaffledPressure, MAX_HISTORY, MAX_RAYLEIGH_INPUTS, MAX_TAP_RING, RayleighMedium, SurfaceSample};
use super::gas::{GasPropagation, GasPropagationReport};
use fs_material::gas::GasState;

/// The existing shared surface quadrature ceiling, independent of basis size.
pub const MAX_RAYLEIGH_FLOW_TILES: usize = 120_000;

/// One actual physical tile; its input is outward volume flow [m³/s].
/// The original surface integrator receives its normal velocity Q/area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowSurfaceSample {
    pub position_m: [f64; 3],
    pub area_m2: f64,
}

/// Aggregate retained storage bounds for all original pieces together.
/// Each piece also keeps the original receiver's cold-construction bounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TilePressureLimits {
    pub maximum_tiles: usize,
    /// Floating-point propagation, filter, work-buffer and tap-history values.
    pub maximum_history_values: usize,
    /// Floating-point kernel coefficients, including the transposed layout.
    pub maximum_kernel_values: usize,
}

impl TilePressureLimits {
    /// Explicit bounded campaign allowance; this is not a performance claim.
    /// History and kernel ceilings are 256 MiB each, excluding object metadata.
    pub const fn surface_reference() -> Self {
        Self {
            maximum_tiles: MAX_RAYLEIGH_FLOW_TILES,
            maximum_history_values: 1 << 25,
            maximum_kernel_values: 1 << 25,
        }
    }
}

/// One physical microphone over arbitrary independently resolved surface tiles.
/// Original receivers keep at most 128 inputs each. No tiles are merged into
/// a piston or resampled spatially. A one-piece surface is bit-identical to
/// its original receiver; larger surfaces retain fixed source-order sums.
/// Large banks reuse the original one-channel pressure differentiator within
/// each piece, commuting the same linear filter without one FIR per tile.
/// Failed sample previews advance no piece; all hot buffers are preallocated.
#[derive(Clone)]
pub struct BaffledTilePressure {
    receivers: Vec<BaffledPressure>,
    inputs: usize,
    limits: TilePressureLimits,
    history_values: usize,
    kernel_values: usize,
    material_latency: Option<usize>,
    started: bool,
    pub position_m: [f64; 3],
    /// Complete original FIR support, including physical path differences.
    pub delay_samples: (usize, usize),
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), String> {
    if cancelled() { Err("Rayleigh tile preparation or observation cancelled".into()) }
    else { Ok(()) }
}
fn add(a: usize, b: usize) -> Result<usize, String> {
    a.checked_add(b).ok_or_else(|| "Rayleigh tile storage count overflow".into())
}

fn surface(
    samples: &[FlowSurfaceSample], cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<SurfaceSample>, String> {
    let mut result = Vec::with_capacity(samples.len());
    for (index, sample) in samples.iter().enumerate() {
        poll(cancelled)?;
        if !sample.area_m2.is_finite() || sample.area_m2 <= 0.0 {
            return Err("Rayleigh flow tile requires its finite positive physical area".into());
        }
        let mut mode_shape = vec![0.0; samples.len()];
        mode_shape[index] = 1.0 / sample.area_m2;
        result.push(SurfaceSample {
            position_m: sample.position_m, area_m2: sample.area_m2, mode_shape,
        });
    }
    Ok(result)
}

// Count the exact retained floating-point layouts after a cold padding shift.
// No receiver changes until the complete bank's target layout is admitted.
fn values_after(receiver: &BaffledPressure, extra: usize) -> Result<(usize, usize), String> {
    let frames = receiver.history_frames.checked_add(extra)
        .filter(|&frames| frames <= MAX_HISTORY).ok_or("Rayleigh tile causal history overflow")?;
    let history = frames.checked_mul(receiver.modes).ok_or("Rayleigh tile history size overflow")?;
    let ring = frames.checked_mul(receiver.kernels.len()).filter(|&n| n <= MAX_TAP_RING);
    let mut histories = add(history, ring.unwrap_or(0))?;
    for count in [
        receiver.previous_velocity.len(), receiver.acceleration.len(),
        receiver.decimator.retained_values(),
        receiver.pressure_derivative.as_ref().map_or(0, |filter| filter.retained_values()),
    ] { histories = add(histories, count)?; }
    let mut coefficients = 0usize;
    for (_, row) in &receiver.kernels { coefficients = add(coefficients, row.len())?; }
    if ring.is_some() {
        coefficients = add(coefficients, receiver.kernels.len().checked_mul(receiver.modes)
            .ok_or("Rayleigh tile packed kernel overflow")?)?;
    }
    Ok((histories, coefficients))
}

impl BaffledTilePressure {
    fn empty(
        inputs: usize, position_m: [f64; 3], limits: TilePressureLimits,
    ) -> Result<Self, String> {
        if inputs == 0 || inputs > limits.maximum_tiles
            || limits.maximum_tiles > MAX_RAYLEIGH_FLOW_TILES
            || limits.maximum_history_values == 0 || limits.maximum_kernel_values == 0
        {
            return Err("Rayleigh tile count and aggregate retained-storage bounds are invalid".into());
        }
        Ok(Self {
            receivers: Vec::with_capacity(inputs.div_ceil(MAX_RAYLEIGH_INPUTS)),
            inputs, limits, history_values: 0, kernel_values: 0,
            material_latency: None, started: false, position_m, delay_samples: (0, 0),
        })
    }

    fn admit_values(&self, histories: usize, kernels: usize) -> Result<(), String> {
        if histories > self.limits.maximum_history_values || kernels > self.limits.maximum_kernel_values {
            return Err("Rayleigh tiles exceed their aggregate retained history or kernel allowance".into());
        }
        Ok(())
    }

    fn push(&mut self, receiver: BaffledPressure) -> Result<(), String> {
        let (history, kernels) = values_after(&receiver, 0)?;
        let history = add(self.history_values, history)?;
        let kernels = add(self.kernel_values, kernels)?;
        self.admit_values(history, kernels)?;
        self.history_values = history;
        self.kernel_values = kernels;
        self.delay_samples = if self.receivers.is_empty() { receiver.delay_samples }
            else { (self.delay_samples.0.min(receiver.delay_samples.0),
                self.delay_samples.1.max(receiver.delay_samples.1)) };
        self.receivers.push(receiver);
        Ok(())
    }

    /// Observe signed original tile flows on one regular output clock.
    /// The caller owns any conservative physical-to-output time projection.
    pub fn new_differentiating(
        samples: &[FlowSurfaceSample], rate: u32, position_m: [f64; 3],
        medium: RayleighMedium, limits: TilePressureLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, String> {
        poll(cancelled)?;
        let mut bank = Self::empty(samples.len(), position_m, limits)?;
        for chunk in samples.chunks(MAX_RAYLEIGH_INPUTS) {
            let basis = surface(chunk, cancelled)?;
            let mut receiver = BaffledPressure::new_differentiating(
                &basis, chunk.len(), rate, position_m, medium)?;
            if samples.len() > MAX_RAYLEIGH_INPUTS { receiver.differentiate_pressure()?; }
            bank.push(receiver)?;
        }
        poll(cancelled)?;
        Ok(bank)
    }

    /// Apply the SAME original material law to every original physical path.
    /// Path/filter/kernel-work allowances are shared across the entire bank,
    /// including repeated actual path preparations in different pieces.
    /// Reports retain exact per-tile source order; equal distances are never
    /// rounded or placed in approximate bins. Only common numerical padding
    /// is aligned, preserving physical travel times and dispersion.
    #[allow(clippy::too_many_arguments)]
    pub fn new_gas_differentiating(
        samples: &[FlowSurfaceSample], rate: u32, position_m: [f64; 3],
        gas: &GasState, propagation: &GasPropagation, limits: TilePressureLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, GasPropagationReport), String> {
        poll(cancelled)?;
        let mut bank = Self::empty(samples.len(), position_m, limits)?;
        let mut remaining = propagation.limits;
        let mut combined: Option<GasPropagationReport> = None;
        for chunk in samples.chunks(MAX_RAYLEIGH_INPUTS) {
            let basis = surface(chunk, cancelled)?;
            let mut bounded = propagation.clone();
            bounded.limits = remaining;
            let (mut receiver, report) = BaffledPressure::new_multirate_gas_differentiating(
                &basis, chunk.len(), rate, rate, position_m, gas, &bounded, cancelled)?;
            remaining.max_unique_paths = remaining.max_unique_paths.checked_sub(report.filters.len())
                .ok_or("Rayleigh tile bank exceeded its shared path allowance")?;
            let filter_taps = report.filters.iter().try_fold(0usize,
                |sum, filter| add(sum, filter.taps))?;
            remaining.max_filter_taps = remaining.max_filter_taps.checked_sub(filter_taps)
                .ok_or("Rayleigh tile bank exceeded its shared filter allowance")?;
            remaining.max_kernel_terms = remaining.max_kernel_terms.checked_sub(report.kernel_terms)
                .ok_or("Rayleigh tile bank exceeded its shared kernel-work allowance")?;
            if samples.len() > MAX_RAYLEIGH_INPUTS { receiver.differentiate_pressure()?; }
            bank.push(receiver)?;
            if let Some(total) = &mut combined {
                let offset = total.filters.len();
                for index in report.surface_filter_indices {
                    total.surface_filter_indices.push(add(offset, index)?);
                }
                total.filters.extend(report.filters);
                total.additional_latency_samples =
                    total.additional_latency_samples.max(report.additional_latency_samples);
                total.kernel_terms = add(total.kernel_terms, report.kernel_terms)?;
            } else { combined = Some(report); }
        }
        let report = combined.expect("admitted nonempty tile surface");
        bank.pad_to(report.additional_latency_samples)?;
        poll(cancelled)?;
        Ok((bank, report))
    }

    fn pad_to(&mut self, latency: usize) -> Result<(), String> {
        if self.started { return Err("Rayleigh tile latency alignment requires a cold receiver".into()); }
        let (mut histories, mut kernels) = (0usize, 0usize);
        for receiver in &self.receivers {
            let extra = latency.checked_sub(receiver.material_path_delay)
                .ok_or("Rayleigh tile alignment cannot remove material padding")?;
            let values = values_after(receiver, extra)?;
            histories = add(histories, values.0)?;
            kernels = add(kernels, values.1)?;
        }
        self.admit_values(histories, kernels)?;
        // Every fallible shift extent was checked above for every piece.
        for receiver in &mut self.receivers {
            let extra = latency - receiver.material_path_delay;
            receiver.add_gas_padding(extra)?;
        }
        self.history_values = histories;
        self.kernel_values = kernels;
        self.material_latency = Some(latency);
        self.delay_samples = (
            self.receivers.iter().map(|r| r.delay_samples.0).min().expect("nonempty bank"),
            self.receivers.iter().map(|r| r.delay_samples.1).max().expect("nonempty bank"),
        );
        Ok(())
    }

    /// Align this complete material surface with another original aperture.
    /// No physical path delay or source phase is removed or reinterpreted.
    pub(crate) fn align_gas_latency(
        &mut self, report: &mut GasPropagationReport, latency: usize,
    ) -> Result<(), String> {
        if self.material_latency != Some(report.additional_latency_samples) {
            return Err("Rayleigh tile material report does not match its receiver".into());
        }
        self.pad_to(latency)?;
        report.additional_latency_samples = latency;
        Ok(())
    }

    pub const fn input_count(&self) -> usize { self.inputs }
    pub fn receiver_count(&self) -> usize { self.receivers.len() }
    pub fn filter_delay_samples(&self) -> f64 {
        self.receivers.first().expect("nonempty bank").filter_delay_samples()
    }
    pub fn multiply_adds_per_sample(&self) -> usize {
        self.receivers.iter().map(BaffledPressure::multiply_adds_per_sample).sum()
    }
    pub const fn retained_values(&self) -> (usize, usize) {
        (self.history_values, self.kernel_values)
    }

    pub fn step_trace(&mut self, flows: &[f64]) -> Result<f64, String> {
        self.step_trace_with_cancel(flows, &mut || false)
    }

    /// Preview every piece and the complete pressure sum before committing
    /// any original observer. Invalid inputs, overflow and cancellation retain
    /// the last accepted history in every piece.
    pub fn step_trace_with_cancel(
        &mut self, flows: &[f64], cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, String> {
        if flows.len() != self.inputs || flows.iter().any(|q| !q.is_finite()) {
            return Err("Rayleigh tile observation requires every finite original flow".into());
        }
        let mut offset = 0usize;
        let mut sum: Option<f64> = None;
        for receiver in &mut self.receivers {
            poll(cancelled)?;
            let end = offset + receiver.modes;
            let value = receiver.preview_trace(&flows[offset..end])?;
            let value = sum.map_or(value, |prior| prior + value);
            if !value.is_finite() { return Err("Rayleigh tile pressure sum overflow".into()); }
            sum = Some(value);
            offset = end;
        }
        poll(cancelled)?;
        offset = 0;
        for receiver in &mut self.receivers {
            let end = offset + receiver.modes;
            receiver.accept_trace(&flows[offset..end]);
            offset = end;
        }
        self.started = true;
        Ok(sum.expect("nonempty bank"))
    }
}

#[cfg(test)]
#[path = "baffled_pressure_tiles_tests.rs"]
mod tests;
