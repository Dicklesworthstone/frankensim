//! Recorded generator coherence through the original loaded opening flows.

use super::*;
use fs_fft::welch::WelchCrossSpectrum;

const MAX_SEGMENTS: usize = 4096;
const MAX_FACTOR_ENTRIES: usize = MAX_SEGMENTS * MAX_NETWORK_PORTS;
const MAX_MAPPING_PRODUCTS: usize = MAX_SEGMENTS * MAX_NETWORK_PORTS * (MAX_NETWORK_PORTS + 2);

/// Whole-bin bounds for cold source-to-opening composition. The original duct,
/// receiver, microphone and angular owners retain their own observation limits.
/// No dense source covariance or source-by-receiver factor bank is allocated.
#[derive(Clone, Copy, Debug)]
pub struct RadiatingWelchLimits {
    /// Recorded segment factors; also subject to the fixed ceiling of 4096.
    /// These are Welch averaging factors, not independent statistical trials.
    pub maximum_segments: usize,
    /// Product allowance: segments * sources * (openings + 2), including the
    /// shared peak conversion's two scaling products per source and segment.
    pub maximum_mapping_products: usize,
    /// Logical bytes for segments * openings peak complex flow factors.
    /// At most 4096 * 32 complex values (2 MiB), independent of receiver count.
    pub maximum_factor_bytes: usize,
}

/// One measured interior bin bound to the exact loaded network and recording.
/// The cache contains actual outward opening-flow peak factors [m³/s], including
/// all original source impedance, interior propagation and mutual radiation
/// loading. It cannot be attached to another physical preparation.
///
/// Cold composition is performed once. Repeated field, microphone and far-field
/// observations allocate nothing and perform no source, graph or BEM solves.
/// Results are bin-integrated moments/powers, not densities per Hz. No temporal
/// realization, cross-frequency coherence or stationarity is inferred.
pub struct PreparedRadiatingWelchBin<'p, 'a> {
    pub(super) owner: &'p PreparedRadiatingNetwork<'a>,
    spectrum: &'p WelchCrossSpectrum,
    bin: usize,
    pub(super) segments: usize,
    pub(super) opening_factors: Vec<C64>,
    mapping_products: usize,
}

impl<'a> PreparedRadiatingNetwork<'a> {
    /// Bind measured source channels in exact canonical label AND SI quantity
    /// order. Pressure channels are open-circuit source pressures before their
    /// complete source impedance; prescribed flow channels retain graph-inflow
    /// m³/s. No channel permutation or internal-pressure proxy is inferred.
    ///
    /// Physical frequency, order, units, whole-bin work and storage are checked
    /// before allocation, including for silent records. DC/Nyquist refuse. The
    /// shared Welch owner applies sqrt(2*df) once; the original observation owners
    /// subsequently apply the peak time-average half once. All segment factors
    /// remain available together for actual-covariance angular convergence.
    pub fn prepare_welch_bin<'p>(
        &'p self,
        spectrum: &'p WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[NetworkSourceCoordinate],
        limits: RadiatingWelchLimits,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingWelchBin<'p, 'a>, ExteriorNetworkError> {
        poll(gate)?;
        let view = spectrum
            .peak_bin_rad_s(bin, self.omega())
            .map_err(|_| bad("recorded exterior bin must match the positive physical frequency"))?;
        let sources = self.source_coordinates().len();
        if recorded_order != self.source_coordinates() || view.channel_count() != sources {
            return Err(bad(
                "recorded exterior labels, SI source quantities or order differ",
            ));
        }
        let segments = view.segment_count();
        let entries = segments
            .checked_mul(self.opening_indices.len())
            .ok_or_else(|| bad("recorded opening factor count overflow"))?;
        let mapping_products = segments
            .checked_mul(sources)
            .and_then(|count| count.checked_mul(self.opening_indices.len() + 2))
            .ok_or_else(|| bad("recorded opening mapping work overflow"))?;
        let bytes = entries
            .checked_mul(core::mem::size_of::<C64>())
            .ok_or_else(|| bad("recorded opening factor storage overflow"))?;
        if segments > MAX_SEGMENTS
            || segments > limits.maximum_segments
            || entries > MAX_FACTOR_ENTRIES
            || mapping_products > MAX_MAPPING_PRODUCTS
            || mapping_products > limits.maximum_mapping_products
            || bytes > limits.maximum_factor_bytes
        {
            return Err(bad(
                "recorded opening factors exceed whole-bin work or storage bounds",
            ));
        }
        let mut opening_factors = Vec::new();
        opening_factors
            .try_reserve_exact(entries)
            .map_err(|_| bad("recorded opening factor allocation refused"))?;
        let mut excitation = [C64::ZERO; MAX_NETWORK_PORTS];
        for factor in view.factors() {
            poll(gate)?;
            for (out, &value) in excitation[..sources].iter_mut().zip(factor) {
                let value = view
                    .peak_factor(value)
                    .map_err(|_| bad("recorded exterior source peak is unrepresentable"))?;
                *out = C64::new(value.re, value.im);
            }
            for &port in &self.opening_indices {
                // Original graph inflow is the negative of outward radiation.
                // Sum complex flows before applying any quadratic observable.
                opening_factors.push(-sum(excitation[..sources].iter().enumerate().map(
                    |(source, &value)| {
                        self.network.flow_transfer()[port * sources + source] * value
                    },
                ))?);
            }
        }
        poll(gate)?;
        Ok(PreparedRadiatingWelchBin {
            owner: self,
            spectrum,
            bin,
            segments,
            opening_factors,
            mapping_products,
        })
    }
}

impl PreparedRadiatingWelchBin<'_, '_> {
    /// Exact immutable recording; its Welch/window/segment provenance is kept.
    pub fn spectrum(&self) -> &WelchCrossSpectrum {
        self.spectrum
    }
    pub fn bin(&self) -> usize {
        self.bin
    }
    pub fn omega(&self) -> f64 {
        self.owner.omega()
    }
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        self.owner.source_coordinates()
    }
    pub fn segment_count(&self) -> usize {
        self.segments
    }
    /// Actual logical cache storage; no additional receiver-sized factor bank.
    pub fn factor_bytes(&self) -> usize {
        self.opening_factors.len() * core::mem::size_of::<C64>()
    }
    /// Source-to-opening and peak-conversion products admitted once.
    pub fn mapping_products(&self) -> usize {
        self.mapping_products
    }

    /// Original signed generator/source/interface/network work plus physical
    /// receiver p/u statistics of the SAME measured excitation. The receiver
    /// fields use cached outward flows, never a surrogate internal pressure.
    /// Exterior surface work already belongs to the original external-load
    /// ledger and is reported without debiting it a second time.
    ///
    /// The network and field budgets independently bound their original
    /// contractions; neither repeats cold opening composition. Cancellation or
    /// any refusal invalidates all output getters, even after a prior success.
    pub fn evaluate_into(
        &self,
        network_budget: NetworkCoherenceBudget,
        field_budget: HarmonicWorkBudget,
        output: &mut RadiatingSourceStatistics<'_, '_>,
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        output.valid = false;
        if !std::ptr::eq(self.owner, output.owner) {
            return Err(bad(
                "measured statistics buffer belongs to another preparation",
            ));
        }
        poll(gate)?;
        self.owner.network.observe_welch_bin_into(
            self.spectrum,
            self.bin,
            self.owner.source_coordinates(),
            network_budget,
            &mut output.network,
            || gate.is_requested(),
        )?;
        self.owner.observer.evaluate_coherence_modes(
            &self.opening_factors,
            self.segments,
            field_budget,
            &mut output.fields,
            || gate.is_requested(),
        )?;
        self.owner.finish_statistics(output, gate)
    }
}

impl PreparedRadiatingMicrophones<'_, '_> {
    /// Full Hermitian microphone bin CSD [Pa²] from actual loaded opening
    /// flows. Original near-field p/u directionality, delays and receiver
    /// cross terms are retained; mapped openings are not independent powers.
    pub fn welch_bin(
        &mut self,
        measured: &PreparedRadiatingWelchBin<'_, '_>,
        budget: HarmonicWorkBudget,
        gate: &CancelGate,
    ) -> Result<&[C64], ExteriorNetworkError> {
        if !std::ptr::eq(self.owner, measured.owner) {
            return Err(bad(
                "measured microphone bin belongs to another preparation",
            ));
        }
        poll(gate)?;
        Ok(self.array.evaluate_coherence_modes(
            &measured.opening_factors,
            measured.segments,
            budget,
            || gate.is_requested(),
        )?)
    }
}
