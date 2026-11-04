//! Explicit stationary sums of selected original loaded Welch bins.

use super::*;
use fs_fft::welch::WelchCrossSpectrum;

mod microphones;
mod radiation;
mod statistics;
pub use microphones::PreparedRadiatingWelchBandMicrophones;
pub use radiation::{PreparedRadiatingWelchBandFarField, RadiatingWelchBandPowerMoments};
pub use statistics::RadiatingWelchBandStatistics;

const MAX_BINS: usize = 4096;

/// Bounds the complete borrowed bank, including the opening caches and cold
/// mappings already admitted by its original bin preparations. No factor bank
/// is copied or recomputed. Observation workspaces have separate byte limits.
#[derive(Clone, Copy, Debug)]
pub struct RadiatingWelchBandLimits {
    pub maximum_bins: usize,
    pub maximum_total_factor_bytes: usize,
    pub maximum_total_mapping_products: usize,
}

/// Exact discrete coverage. The selected indices remain available through
/// `bins()`. Missing interior bins contribute no claimed result; they are not
/// inferred silent or filled by interpolation. DC and Nyquist are always absent.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RadiatingWelchBandCoverage {
    pub selected_bins: usize,
    pub available_interior_bins: usize,
}
impl RadiatingWelchBandCoverage {
    pub fn omitted_interior_bins(self) -> usize {
        self.available_interior_bins - self.selected_bins
    }
    /// Even when true, this does not include DC or Nyquist or certify measured
    /// bandwidth, window leakage, stationarity, or unrecorded frequencies.
    pub fn covers_all_interior_bins(self) -> bool {
        self.selected_bins == self.available_interior_bins
    }
}

/// A declared stationary sum over explicitly selected, strictly increasing bins
/// of ONE exact recording and ONE original physical network/exterior model.
/// Each bin retains its own frequency, loaded flow transfer, original source
/// laws and observations. Conventional propagation must use the same loss law;
/// molecular propagation must borrow the exact same original reservoir owner.
/// Canonical labels/units and physical p/u receiver and
/// opening identities must match. Matching labels does not fit or interpolate
/// any source impedance; those original per-frequency laws remain inspectable.
///
/// This operation assumes that cross-frequency products average to zero in the
/// requested stationary statistic. It never derives that assumption from Welch
/// segments, constructs cross-frequency coherence, or synthesizes a waveform.
/// The selected bins can have gaps. Their ALREADY integrated moments are summed
/// once without another df, one-sided, or peak-to-mean-square multiplier.
#[derive(Clone, Copy)]
pub struct PreparedRadiatingWelchBand<'b, 'p, 'a> {
    bins: &'b [PreparedRadiatingWelchBin<'p, 'a>],
    factor_bytes: usize,
    mapping_products: usize,
}

impl<'b, 'p, 'a> PreparedRadiatingWelchBand<'b, 'p, 'a> {
    /// Explicitly declare the stationary-bin interpretation. All identities
    /// and limits are checked even when every selected record bin is silent.
    pub fn stationary(
        bins: &'b [PreparedRadiatingWelchBin<'p, 'a>],
        limits: RadiatingWelchBandLimits,
        gate: &CancelGate,
    ) -> Result<Self, ExteriorNetworkError> {
        poll(gate)?;
        if bins.is_empty() || bins.len() > MAX_BINS || bins.len() > limits.maximum_bins {
            return Err(bad("stationary exterior band needs bounded explicit bins"));
        }
        let first = &bins[0];
        let mut factor_bytes = 0;
        let mut mapping_products = 0;
        for (index, bin) in bins.iter().enumerate() {
            poll(gate)?;
            if !core::ptr::eq(bin.spectrum(), first.spectrum())
                || !core::ptr::eq(bin.owner.model, first.owner.model)
                || !same_propagation(bin.owner.propagation, first.owner.propagation)
                || bin.source_coordinates() != first.source_coordinates()
                || bin.owner.ports() != first.owner.ports()
                || bin.owner.observer.sources() != first.owner.observer.sources()
                || bin.owner.observer.receivers() != first.owner.observer.receivers()
            {
                return Err(bad(
                    "band bins require the same recording, physical model, propagation and source/receiver identities",
                ));
            }
            if index != 0 && bins[index - 1].bin() >= bin.bin() {
                return Err(bad("band indices must be explicit, unique and increasing"));
            }
            factor_bytes = add(factor_bytes, bin.factor_bytes())?;
            mapping_products = add(mapping_products, bin.mapping_products())?;
        }
        if factor_bytes > limits.maximum_total_factor_bytes
            || mapping_products > limits.maximum_total_mapping_products
        {
            return Err(bad(
                "complete band exceeds its original factor/mapping allowance",
            ));
        }
        poll(gate)?;
        Ok(Self {
            bins,
            factor_bytes,
            mapping_products,
        })
    }

    pub fn bins(&self) -> &'b [PreparedRadiatingWelchBin<'p, 'a>] {
        self.bins
    }
    pub fn spectrum(&self) -> &WelchCrossSpectrum {
        self.bins[0].spectrum()
    }
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        self.bins[0].source_coordinates()
    }
    pub fn ports(&self) -> &[String] {
        self.bins[0].owner.ports()
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        self.bins[0].owner.observer.receivers()
    }
    pub fn coverage(&self) -> RadiatingWelchBandCoverage {
        RadiatingWelchBandCoverage {
            selected_bins: self.bins.len(),
            available_interior_bins: self.spectrum().bin_count() - 2,
        }
    }
    pub fn factor_bytes(&self) -> usize {
        self.factor_bytes
    }
    pub fn mapping_products(&self) -> usize {
        self.mapping_products
    }
}

fn same_propagation(a: Propagation<'_>, b: Propagation<'_>) -> bool {
    match (a, b) {
        (Propagation::Duct(a), Propagation::Duct(b)) => a == b,
        (Propagation::Molecular(a), Propagation::Molecular(b)) => core::ptr::eq(a, b),
        _ => false,
    }
}

fn add(a: usize, b: usize) -> Result<usize, ExteriorNetworkError> {
    a.checked_add(b)
        .ok_or_else(|| bad("stationary band count overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize, ExteriorNetworkError> {
    a.checked_mul(b)
        .ok_or_else(|| bad("stationary band count overflow"))
}
fn bytes<T>(n: usize) -> Result<usize, ExteriorNetworkError> {
    mul(n, core::mem::size_of::<T>())
}
fn storage(required: usize, maximum: usize) -> Result<(), ExteriorNetworkError> {
    if required > maximum {
        Err(bad("complete band observation workspace exceeds allowance"))
    } else {
        Ok(())
    }
}
fn reserve<T>(n: usize) -> Result<Vec<T>, ExteriorNetworkError> {
    let mut out = Vec::new();
    out.try_reserve_exact(n)
        .map_err(|_| bad("band observation allocation refused"))?;
    Ok(out)
}
fn zeros<T: Default + Clone>(n: usize) -> Result<Vec<T>, ExteriorNetworkError> {
    let mut out = reserve(n)?;
    out.resize(n, T::default());
    Ok(out)
}

// Each bin's original budget remains authoritative. The whole-band allowance
// conservatively counts the SUM OF THOSE ALLOWANCES plus one product-equivalent
// per scalar/complex reduction entry. Unspent per-bin allowances are not pooled
// or used to weaken another bin's original gate. Check before any observation.
fn work(
    bins: usize,
    per_bin: usize,
    reduced_entries: usize,
    budget: HarmonicWorkBudget,
) -> Result<(), ExteriorNetworkError> {
    let required = mul(bins, add(per_bin, reduced_entries)?)?;
    if required > budget.maximum_products {
        Err(bad("complete band observation exceeds work allowance"))
    } else {
        Ok(())
    }
}

// Reuse the original complex Neumaier reduction. At most MAX_BINS values lie
// between these polls; per-bin physical contractions retain their own polling.
fn reduce(
    values: impl Iterator<Item = C64>,
    gate: &CancelGate,
) -> Result<C64, ExteriorNetworkError> {
    poll(gate)?;
    let value = sum(values)?;
    poll(gate)?;
    Ok(value)
}

// Streaming version of the same componentwise Neumaier reduction used above.
// It avoids retaining bins*microphones² or bins*directions extra observations.
#[derive(Clone, Copy, Default)]
struct Accumulator {
    value: C64,
    correction: C64,
}
impl Accumulator {
    fn add(&mut self, value: C64) {
        fn component(sum: &mut f64, correction: &mut f64, value: f64) {
            let next = *sum + value;
            *correction += if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            };
            *sum = next;
        }
        component(&mut self.value.re, &mut self.correction.re, value.re);
        component(&mut self.value.im, &mut self.correction.im, value.im);
    }
    fn finish(self) -> Result<C64, ExteriorNetworkError> {
        let result = self.value + self.correction;
        if finite(result) {
            Ok(result)
        } else {
            Err(bad("stationary band moment overflow"))
        }
    }
}
