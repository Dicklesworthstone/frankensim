//! Finite-window p/u observations of the original prepared BEM transfers.
//!
//! Source phasors retain the common exp(-i omega t) clock. Deterministic
//! fields include sum- and difference-frequency terms; explicitly proper
//! fluctuation factors include their cross-frequency coherence. No new
//! propagation solve, plane-wave approximation or source statistics are inferred.

use super::super::{ProducerUncertainty, ReceiverSpec, SourceId, TransferError};
use super::{MAX_PRODUCTS, PreparedHarmonicObservation};
use fs_blake3::ContentHash;
use fs_fft::harmonic_window::{
    HarmonicProduct, HarmonicWindow, ProperHarmonicProduct, WindowError,
};
use fs_math::{c64::C64, dd::Dd};

mod array;
pub use array::{
    HarmonicPressureArrayLimits, HarmonicPressureChannel, HarmonicPressureChannelMoments,
    HarmonicPressureTap, MAX_PRESSURE_ARRAY_CHANNELS, MAX_PRESSURE_ARRAY_TAPS,
    PreparedHarmonicPressureArray, WindowedHarmonicPressureArray,
};

mod readout;
pub use readout::{
    HarmonicAveragedField, HarmonicReadoutCovariance, HarmonicReadoutLimits,
    HarmonicReceiverPair, MAX_READOUT_RECEIVER_PAIRS, WindowedHarmonicReadout,
};

/// Bounds the deterministic in-place sort used for frequency grouping.
pub const MAX_EXPOSURE_LINES: usize = 4096;
/// Empty factors count too: an authored factor bank must be bounded.
pub const MAX_EXPOSURE_FACTORS: usize = 256;

/// Source coefficients through one existing frequency-specific physical owner.
/// Coordinates retain that owner's units and order. All lines share the
/// declared emission clock; propagation phase is already in the transfer.
#[derive(Clone, Copy)]
pub struct HarmonicExposureLine<'a> {
    pub observation: &'a PreparedHarmonicObservation,
    pub coefficients: &'a [C64],
}

/// Whole-preparation and per-window admission. Counts include supplied lines
/// before equal-frequency coalescing, including silent lines and factors.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicExposureLimits {
    pub maximum_lines: usize,
    pub maximum_factors: usize,
    /// Four complex p/u source-transfer products per source, receiver and line.
    pub maximum_transfer_products: usize,
    /// Receiver count times sum of squared line counts over mean/factor banks.
    /// Each pair has a fixed four-component (ten unique moments) contraction.
    pub maximum_receiver_frequency_pairs: usize,
    /// Logical heap payload: retained fields, line/factor vector headers,
    /// grouping indices, query accumulators and one returned point array.
    /// Borrowed owners/inputs, stack values and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Uncentered real-field moments at one physical receiver. Pressure/velocity
/// covariance is not inferred from a plane-wave impedance, even near a body.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicExposureMoments {
    pub mean_pressure_pa: f64,
    pub mean_particle_velocity_m_s: [f64; 3],
    pub mean_square_pressure_pa2: f64,
    /// E[<u_i(t) u_j(t)>], real symmetric, including off-diagonal terms.
    pub velocity_second_m2_s2: [[f64; 3]; 3],
    /// E[<p(t) u(t)>], signed instantaneous energy flux averaged over the window.
    /// This is not a cycle-averaged complex/reactive-intensity surrogate.
    pub mean_intensity_w_m2: [f64; 3],
}

struct Line {
    omega: f64,
    // p, ux, uy, uz, each [real, imaginary]. Low parts survive all products.
    fields: Vec<[Dd; 8]>,
}

struct Bank {
    lines: Vec<Line>,
}

/// Prepared fields for a deterministic mean and explicitly joint proper
/// fluctuation factors: E[xi_r]=0, E[xi_r conj(xi_s)]=delta_rs,
/// E[xi_r xi_s]=0. A factor spans ALL its frequencies and source coordinates.
/// Factors add second moments; mean/fluctuation cross terms vanish.
/// Weights belong in the peak amplitudes, not in an unscaled PSD/Hz array.
/// No Gaussian distribution or stationarity is required or inferred.
///
/// Every line must match the canonical owner's producer hash, ordered sources,
/// exact receiver geometry and original medium. These checks prevent mixing
/// different geometries/operators under the same receiver labels. Their
/// original spatial/solver uncertainty remains unresolved.
///
/// Preparation contracts existing transfers once in O(receivers * sources *
/// lines), without a BEM solve. Equal frequencies field-sum in source order
/// before products. Repeated queries share one temporal kernel per frequency
/// pair across all receivers and components. Double-double sums and products
/// improve cancellation, but are not a relative error certificate near nulls.
pub struct PreparedHarmonicExposure<'a> {
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    mean: Bank,
    factors: Vec<Bank>,
    workspace_bytes: usize,
    receiver_frequency_pairs: usize,
}

/// Complete observation; identities and source clock borrow the immutable owner.
pub struct WindowedHarmonicExposure<'a> {
    window: HarmonicWindow,
    receivers: &'a [ReceiverSpec],
    sources: &'a [SourceId],
    source_clock: &'a str,
    producer_hash: ContentHash,
    points: Vec<HarmonicExposureMoments>,
}

impl WindowedHarmonicExposure<'_> {
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        self.receivers
    }
    pub fn sources(&self) -> &[SourceId] {
        self.sources
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        ProducerUncertainty::Unresolved {
            source_hash: self.producer_hash,
        }
    }
    pub fn points(&self) -> &[HarmonicExposureMoments] {
        &self.points
    }
}

impl<'a> PreparedHarmonicExposure<'a> {
    pub fn new(
        reference: &'a PreparedHarmonicObservation,
        source_clock: &'a str,
        mean: &[HarmonicExposureLine<'_>],
        proper_factors: &[&[HarmonicExposureLine<'_>]],
        limits: HarmonicExposureLimits,
    ) -> Result<Self, TransferError> {
        Self::new_cancellable(
            reference,
            source_clock,
            mean,
            proper_factors,
            limits,
            || false,
        )
    }

    /// No partial prepared object is published after cancellation or refusal.
    /// The canonical owner also supplies physical scope for an explicitly dark
    /// field, so no geometry or medium needs to be invented for an empty bank.
    #[allow(clippy::too_many_lines)] // Atomic admission, identity checks and preparation.
    pub fn new_cancellable(
        reference: &'a PreparedHarmonicObservation,
        source_clock: &'a str,
        mean: &[HarmonicExposureLine<'_>],
        proper_factors: &[&[HarmonicExposureLine<'_>]],
        limits: HarmonicExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        if source_clock.trim().is_empty() {
            return Err(bad("BEM exposure needs an explicit common source clock"));
        }
        admit(proper_factors.len(), MAX_EXPOSURE_FACTORS)?;
        admit(proper_factors.len(), limits.maximum_factors)?;
        let mut lines = 0usize;
        let mut pairs = 0usize;
        let mut largest = 0usize;
        for bank in core::iter::once(mean).chain(proper_factors.iter().copied()) {
            lines = lines.checked_add(bank.len()).ok_or_else(overflow)?;
            pairs = bank
                .len()
                .checked_mul(bank.len())
                .and_then(|count| pairs.checked_add(count))
                .ok_or_else(overflow)?;
            largest = largest.max(bank.len());
        }
        let receivers = reference.receivers.len();
        let transfer_products = 4usize
            .checked_mul(receivers)
            .and_then(|count| count.checked_mul(reference.sources.len()))
            .and_then(|count| count.checked_mul(lines))
            .ok_or_else(overflow)?;
        let receiver_pairs = receivers.checked_mul(pairs).ok_or_else(overflow)?;
        admit(lines, MAX_EXPOSURE_LINES)?;
        admit(lines, limits.maximum_lines)?;
        admit(transfer_products, MAX_PRODUCTS)?;
        admit(transfer_products, limits.maximum_transfer_products)?;
        admit(receiver_pairs, MAX_PRODUCTS)?;
        admit(receiver_pairs, limits.maximum_receiver_frequency_pairs)?;
        let per_line = receivers
            .checked_mul(core::mem::size_of::<[Dd; 8]>())
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<Line>()))
            .ok_or_else(overflow)?;
        let workspace_bytes = lines
            .checked_mul(per_line)
            .and_then(|bytes| {
                bytes.checked_add(
                    proper_factors
                        .len()
                        .checked_mul(core::mem::size_of::<Bank>())?,
                )
            })
            .and_then(|bytes| {
                bytes.checked_add(largest.checked_mul(core::mem::size_of::<usize>())?)
            })
            .and_then(|bytes| {
                bytes.checked_add(receivers.checked_mul(
                    core::mem::size_of::<[Dd; 14]>()
                        + core::mem::size_of::<HarmonicExposureMoments>(),
                )?)
            })
            .ok_or_else(overflow)?;
        admit(workspace_bytes, limits.maximum_workspace_bytes)?;

        for bank in core::iter::once(mean).chain(proper_factors.iter().copied()) {
            for line in bank {
                poll(&mut cancelled)?;
                let owner = line.observation;
                if owner.producer_hash != reference.producer_hash
                    || owner.sources != reference.sources
                    || owner.receivers != reference.receivers
                    || owner.medium.density.to_bits() != reference.medium.density.to_bits()
                    || owner.medium.sound_speed.to_bits() != reference.medium.sound_speed.to_bits()
                {
                    return Err(bad(
                        "BEM exposure line belongs to a different physical owner",
                    ));
                }
                owner.factors(line.coefficients, 1)?;
            }
        }
        let mean = prepare_bank(mean, receivers, &mut cancelled)?;
        let mut factors = storage(proper_factors.len())?;
        for bank in proper_factors {
            factors.push(prepare_bank(bank, receivers, &mut cancelled)?);
        }
        let actual_pairs = core::iter::once(&mean)
            .chain(factors.iter())
            .map(|bank| bank.lines.len() * bank.lines.len())
            .sum::<usize>();
        poll(&mut cancelled)?;
        Ok(Self {
            reference,
            source_clock,
            mean,
            factors,
            workspace_bytes,
            receiver_frequency_pairs: receivers * actual_pairs,
        })
    }

    pub fn reference(&self) -> &PreparedHarmonicObservation {
        self.reference
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn factor_count(&self) -> usize {
        self.factors.len()
    }
    /// Actual receiver/frequency work after exact-frequency grouping.
    pub fn receiver_frequency_pairs(&self) -> usize {
        self.receiver_frequency_pairs
    }

    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedHarmonicExposure<'_>, TransferError> {
        self.observe_window_cancellable(window, || false)
    }

    /// No solve, sampling or implicit carrier-cycle average occurs. Cancellation
    /// discards the local result and leaves this immutable owner reusable.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedHarmonicExposure<'_>, TransferError> {
        poll(&mut cancelled)?;
        let receivers = self.reference.receivers.len();
        let mut sums = storage(receivers)?;
        sums.resize(receivers, [Dd::ZERO; 14]);
        for line in &self.mean.lines {
            poll(&mut cancelled)?;
            // Integrate against an exact unit DC field: do not round the mean
            // phase to C64 before cancellation between large spectral lines.
            let kernel = window
                .product_kernel(line.omega, 0.0)
                .map_err(window_error)?;
            for (sum, field) in sums.iter_mut().zip(&line.fields) {
                for c in 0..4 {
                    sum[c] = sum[c]
                        + kernel
                            .mean_real_product_parts_dd(component(field, c), [Dd::ONE, Dd::ZERO])
                            .map_err(window_error)?;
                }
            }
        }
        for (bank_index, bank) in core::iter::once(&self.mean)
            .chain(self.factors.iter())
            .enumerate()
        {
            poll(&mut cancelled)?;
            for (i, a) in bank.lines.iter().enumerate() {
                for (j, b) in bank.lines[..=i].iter().enumerate() {
                    poll(&mut cancelled)?;
                    let kernel = if bank_index == 0 {
                        Kernel::Deterministic(
                            window
                                .product_kernel(a.omega, b.omega)
                                .map_err(window_error)?,
                        )
                    } else {
                        Kernel::Proper(
                            window
                                .proper_product_kernel(a.omega, b.omega)
                                .map_err(window_error)?,
                        )
                    };
                    for ((sum, af), bf) in sums.iter_mut().zip(&a.fields).zip(&b.fields) {
                        let mut index = 4;
                        for row in 0..4 {
                            for col in row..4 {
                                let mut value =
                                    kernel.product(component(af, row), component(bf, col))?;
                                if i != j {
                                    value = if row == col {
                                        value * Dd::from_f64(2.0)
                                    } else {
                                        value
                                            + kernel
                                                .product(component(af, col), component(bf, row))?
                                    };
                                }
                                sum[index] = sum[index] + value;
                                index += 1;
                            }
                        }
                    }
                }
            }
        }
        let mut points = storage(receivers)?;
        for sum in sums {
            poll(&mut cancelled)?;
            if sum
                .iter()
                .any(|value| !finite(*value) || !value.to_f64().is_finite())
            {
                return Err(bad("BEM exposure moment is unrepresentable"));
            }
            let mut second = [[0.0; 4]; 4];
            let mut index = 4;
            for row in 0..4 {
                for col in row..4 {
                    let value = sum[index].to_f64();
                    if row == col && value < 0.0 {
                        return Err(bad("BEM exposure self moment became negative"));
                    }
                    second[row][col] = value;
                    second[col][row] = value;
                    index += 1;
                }
            }
            points.push(HarmonicExposureMoments {
                mean_pressure_pa: sum[0].to_f64(),
                mean_particle_velocity_m_s: core::array::from_fn(|c| sum[c + 1].to_f64()),
                mean_square_pressure_pa2: second[0][0],
                velocity_second_m2_s2: core::array::from_fn(|row| {
                    core::array::from_fn(|col| second[row + 1][col + 1])
                }),
                mean_intensity_w_m2: core::array::from_fn(|c| second[0][c + 1]),
            });
        }
        poll(&mut cancelled)?;
        Ok(WindowedHarmonicExposure {
            window,
            receivers: &self.reference.receivers,
            sources: &self.reference.sources,
            source_clock: self.source_clock,
            producer_hash: self.reference.producer_hash,
            points,
        })
    }
}

fn prepare_bank(
    input: &[HarmonicExposureLine<'_>],
    receivers: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Bank, TransferError> {
    poll(cancelled)?;
    let mut indices = storage(input.len())?;
    indices.extend(0..input.len());
    // The explicit index tie-break preserves authored order at equal frequency.
    indices.sort_unstable_by(|a, b| {
        input[*a]
            .observation
            .omega
            .total_cmp(&input[*b].observation.omega)
            .then_with(|| a.cmp(b))
    });
    let mut lines: Vec<Line> = storage(input.len())?;
    for index in indices {
        poll(cancelled)?;
        let line = input[index];
        let owner = line.observation;
        if lines
            .last()
            .is_none_or(|last| last.omega.to_bits() != owner.omega.to_bits())
        {
            let mut fields = storage(receivers)?;
            fields.resize(receivers, [Dd::ZERO; 8]);
            lines.push(Line {
                omega: owner.omega,
                fields,
            });
        }
        let fields = &mut lines.last_mut().expect("frequency group exists").fields;
        for (receiver, field) in fields.iter_mut().enumerate() {
            poll(cancelled)?;
            for (source, coefficient) in line.coefficients.iter().enumerate() {
                if source % 16 == 0 {
                    poll(cancelled)?;
                }
                let transfer = owner.transfers[source * receivers + receiver];
                let cr = Dd::from_f64(coefficient.re);
                let ci = Dd::from_f64(coefficient.im);
                for c in 0..4 {
                    let tr = Dd::from_f64(transfer[c].re);
                    let ti = Dd::from_f64(transfer[c].im);
                    field[2 * c] = field[2 * c] + tr * cr - ti * ci;
                    field[2 * c + 1] = field[2 * c + 1] + tr * ci + ti * cr;
                }
            }
        }
    }
    if lines
        .iter()
        .flat_map(|line| &line.fields)
        .flatten()
        .any(|value| !finite(*value))
    {
        return Err(bad("BEM exposure coherent field is unrepresentable"));
    }
    Ok(Bank { lines })
}

#[derive(Clone, Copy)]
enum Kernel {
    Deterministic(HarmonicProduct),
    Proper(ProperHarmonicProduct),
}

impl Kernel {
    fn product(self, a: [Dd; 2], b: [Dd; 2]) -> Result<Dd, TransferError> {
        match self {
            Self::Deterministic(kernel) => kernel.mean_real_product_parts_dd(a, b),
            Self::Proper(kernel) => kernel.mean_real_product_parts_dd(a, b),
        }
        .map_err(window_error)
    }
}

fn component(field: &[Dd; 8], c: usize) -> [Dd; 2] {
    [field[2 * c], field[2 * c + 1]]
}

fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

fn storage<T>(count: usize) -> Result<Vec<T>, TransferError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| bad("BEM exposure allocation failed"))?;
    Ok(values)
}

fn admit(requested: usize, cap: usize) -> Result<(), TransferError> {
    if requested > cap {
        Err(TransferError::Cap { requested, cap })
    } else {
        Ok(())
    }
}

fn bad(what: &'static str) -> TransferError {
    TransferError::Invalid(what)
}

fn overflow() -> TransferError {
    bad("BEM exposure work/storage count overflow")
}

fn window_error(_: WindowError) -> TransferError {
    bad("BEM exposure temporal moment is unrepresentable")
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransferError> {
    if cancelled() {
        Err(TransferError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
