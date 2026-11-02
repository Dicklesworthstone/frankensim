//! Broadband interference of the original native network interior fields.
//!
//! Frequencies are prepared together from one network, medium binding, probe
//! list and source clock. A spectral line names an index in that immutable
//! bank, so a foreign geometry cannot enter merely through matching labels.
//! No frequency-bin bandwidth, random phase or stationary spectrum is inferred.

use super::*;

/// Maximum frequency-specific native loaded solves in one prepared spectrum.
pub const MAX_FIELD_FREQUENCIES: usize = 256;
/// Maximum supplied mean/factor lines, including repeated or silent lines.
pub const MAX_FIELD_SPECTRAL_LINES: usize = 4096;
const MAX_SPECTRAL_TRANSFERS: usize = 16_777_216;
const MAX_SPECTRAL_PRODUCTS: usize = 64_000_000;

/// Original source/load laws at one strictly positive angular frequency.
/// Frequency-dependent impedances remain native inputs, not interpolated data.
#[derive(Clone, Copy)]
pub struct NetworkFieldFrequency<'a> {
    pub omega: f64,
    pub loss: LossModel,
    pub sources: &'a [PortSource],
    pub source_bank: Option<&'a CoupledSourceBank>,
    pub loads: &'a [PortTermination],
    pub coupled_load: Option<&'a CoupledTermination>,
}

/// Whole-spectrum admission, in addition to the original graph/probe budgets.
#[derive(Clone, Copy, Debug)]
pub struct NetworkFieldSpectrumBudget {
    pub maximum_frequencies: usize,
    /// Complex p/U transfer entries summed over all retained frequency banks.
    /// Original per-bank source/port storage and preparation scratch retain
    /// their native hard caps; this is not an allocator-overhead byte account.
    pub maximum_retained_transfer_entries: usize,
}

/// All frequencies share these original physical geometry/clock inputs.
pub struct NetworkFieldSpectrumRequest<'a> {
    pub source_clock: &'a str,
    pub probes: &'a [NetworkFieldProbe],
    /// Strictly increasing positive omega; indices retain this authored order.
    pub frequencies: &'a [NetworkFieldFrequency<'a>],
    pub field_budget: NetworkFieldBudget,
    pub spectrum_budget: NetworkFieldSpectrumBudget,
}

/// One prepared native source-to-field response per frequency.
#[derive(Clone, Debug)]
pub struct PreparedNetworkFieldSpectrum {
    fields: Vec<PreparedNetworkFields>,
    retained_transfer_entries: usize,
}

impl PreparedNetworkFieldSpectrum {
    pub fn new(
        net: &DuctNetwork,
        gas: &GasState,
        request: &NetworkFieldSpectrumRequest<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        Self::prepare(net, gas, None, request, cancelled)
    }

    pub fn new_regional(
        net: &DuctNetwork,
        media: &NetworkMedia,
        request: &NetworkFieldSpectrumRequest<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        media.validate_for(net)?;
        Self::prepare(net, media.exterior_gas(), Some(media), request, cancelled)
    }

    fn prepare(
        net: &DuctNetwork,
        gas: &GasState,
        media: Option<&NetworkMedia>,
        request: &NetworkFieldSpectrumRequest<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(&mut cancelled, 0)?;
        let count = request.frequencies.len();
        if count == 0
            || count > MAX_FIELD_FREQUENCIES
            || count > request.spectrum_budget.maximum_frequencies
            || request.source_clock.trim().is_empty()
        {
            return Err(bad(
                "network field spectrum needs a clock and bounded nonempty frequencies",
            ));
        }
        let mut previous = 0.0;
        let mut entries = 0usize;
        for (index, frequency) in request.frequencies.iter().enumerate() {
            poll(&mut cancelled, index)?;
            if !frequency.omega.is_finite() || frequency.omega <= previous {
                return Err(bad(
                    "network field frequencies must be positive and strictly increasing",
                ));
            }
            previous = frequency.omega;
            let sources = frequency
                .sources
                .len()
                .checked_add(
                    frequency
                        .source_bank
                        .map_or(0, |b| b.impedance.labels.len()),
                )
                .ok_or_else(count_overflow)?;
            let additional = request
                .probes
                .len()
                .checked_mul(sources)
                .and_then(|n| n.checked_mul(2))
                .ok_or_else(count_overflow)?;
            entries = entries.checked_add(additional).ok_or_else(count_overflow)?;
        }
        if entries > MAX_SPECTRAL_TRANSFERS
            || entries > request.spectrum_budget.maximum_retained_transfer_entries
        {
            return Err(bad(
                "network field spectrum exceeds its retained transfer budget",
            ));
        }
        let mut fields: Vec<PreparedNetworkFields> = storage(count)?;
        for (index, frequency) in request.frequencies.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let local = NetworkFieldRequest {
                omega: frequency.omega,
                loss: frequency.loss,
                source_clock: request.source_clock,
                sources: frequency.sources,
                source_bank: frequency.source_bank,
                loads: frequency.loads,
                coupled_load: frequency.coupled_load,
                probes: request.probes,
                budget: request.field_budget,
            };
            let field = if let Some(media) = media {
                PreparedNetworkFields::new_regional(net, media, &local, &mut cancelled)?
            } else {
                PreparedNetworkFields::new(net, gas, &local, &mut cancelled)?
            };
            if let Some(first) = fields.first() {
                if field.port_response().source_coordinates()
                    != first.port_response().source_coordinates()
                    || field.locations() != first.locations()
                {
                    return Err(bad(
                        "network spectrum source coordinates or physical probes differ",
                    ));
                }
            }
            fields.push(field);
        }
        poll(&mut cancelled, count)?;
        Ok(Self {
            fields,
            retained_transfer_entries: entries,
        })
    }

    pub fn frequency_count(&self) -> usize {
        self.fields.len()
    }
    pub fn frequency_fields(&self, index: usize) -> Option<&PreparedNetworkFields> {
        self.fields.get(index)
    }
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        self.fields[0].port_response().source_coordinates()
    }
    pub fn locations(&self) -> &[NetworkFieldLocation] {
        self.fields[0].locations()
    }
    pub fn source_clock(&self) -> &str {
        self.fields[0].source_clock()
    }
    pub fn retained_transfer_entries(&self) -> usize {
        self.retained_transfer_entries
    }

    /// Contract source vectors once and retain coherent DD field sums for
    /// repeated windows. Every spectral line names this owner's frequency bank.
    pub fn prepare_observation(
        &self,
        mean: &[NetworkFieldSpectralLine<'_>],
        factors: &[NetworkFieldSpectralFactor<'_>],
        limits: NetworkFieldExposureLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedNetworkFieldExposure<'_>, DuctError> {
        PreparedNetworkFieldExposure::new(self, mean, factors, limits, cancelled)
    }
}

/// Complete source-coordinate excitation at one frequency of the same owner.
#[derive(Clone, Copy, Debug)]
pub struct NetworkFieldSpectralLine<'a> {
    pub frequency_index: usize,
    /// Peak SI coordinates in PreparedNetworkFieldSpectrum::source_coordinates order.
    pub excitation: &'a [C64],
}

/// Joint second-moment model of one zero-mean latent source coefficient.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NetworkFieldFactorStatistics {
    /// Real E[x]=0, E[x²]=1; retains phase-locked amplitude noise and sum terms.
    Real,
    /// Complex E[z]=E[z²]=0, E[|z|²]=1; retains difference-frequency beats.
    ProperComplex,
}

/// One latent coefficient multiplies ALL lines/source coordinates in this bank.
/// Real factors have joint E[x_r*x_s]=delta_rs. Proper factors have joint
/// E[z_r*conj(z_s)]=delta_rs AND E[z_r*z_s]=0. Mixed real/proper second moments
/// vanish. These are authored moment assumptions, not a distribution model.
/// Source amplitudes include statistical square-root weights, never PSD/Hz.
#[derive(Clone, Copy, Debug)]
pub struct NetworkFieldSpectralFactor<'a> {
    pub lines: &'a [NetworkFieldSpectralLine<'a>],
    pub statistics: NetworkFieldFactorStatistics,
}

/// Preparation and per-window bounds for contracted field observations.
#[derive(Clone, Copy, Debug)]
pub struct NetworkFieldExposureLimits {
    pub maximum_lines: usize,
    pub maximum_factors: usize,
    /// Two complex source contractions per probe/source/supplied line.
    pub maximum_transfer_products: usize,
    /// Probes times the sum of squared line counts over all mean/factor banks.
    pub maximum_probe_frequency_pairs: usize,
    /// Logical numerical heap payload: retained DD fields/vector headers,
    /// grouping indices, query accumulators and complete returned moments.
    /// The borrowed native spectrum, allocator overhead and earlier outputs
    /// are excluded.
    pub maximum_workspace_bytes: usize,
}

struct SpectralLine {
    omega: f64,
    // p.re, p.im, U.re, U.im at every original spatial probe.
    fields: Vec<[Dd; 4]>,
}
struct SpectralBank {
    lines: Vec<SpectralLine>,
    statistics: NetworkFieldFactorStatistics,
}

/// Original spectral owner plus one deterministic mean and real/proper factors.
/// A general PSD covariance of stacked real/imaginary coordinates can be
/// authored as real factors. Gaussianity and stationarity are not inferred.
pub struct PreparedNetworkFieldExposure<'a> {
    spectrum: &'a PreparedNetworkFieldSpectrum,
    mean: SpectralBank,
    factors: Vec<SpectralBank>,
    workspace_bytes: usize,
    probe_frequency_pairs: usize,
}

/// Real-field moments over the ensemble and normalized physical window.
pub struct WindowedNetworkFieldExposure<'a> {
    window: HarmonicWindow,
    spectrum: &'a PreparedNetworkFieldSpectrum,
    moments: Vec<NetworkFieldMoments>,
}

impl WindowedNetworkFieldExposure<'_> {
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }
    pub fn locations(&self) -> &[NetworkFieldLocation] {
        self.spectrum.locations()
    }
    pub fn source_clock(&self) -> &str {
        self.spectrum.source_clock()
    }
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        self.spectrum.source_coordinates()
    }
    pub fn moments(&self) -> &[NetworkFieldMoments] {
        &self.moments
    }
}

impl<'a> PreparedNetworkFieldExposure<'a> {
    fn new(
        spectrum: &'a PreparedNetworkFieldSpectrum,
        mean: &[NetworkFieldSpectralLine<'_>],
        factors: &[NetworkFieldSpectralFactor<'_>],
        limits: NetworkFieldExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(&mut cancelled, 0)?;
        if factors.len() > MAX_FIELD_FACTORS || factors.len() > limits.maximum_factors {
            return Err(bad("network spectral factor count exceeds its budget"));
        }
        let mut lines = 0usize;
        let mut pairs = 0usize;
        let mut largest = 0usize;
        for bank in core::iter::once(mean).chain(factors.iter().map(|f| f.lines)) {
            lines = lines.checked_add(bank.len()).ok_or_else(count_overflow)?;
            pairs = bank
                .len()
                .checked_mul(bank.len())
                .and_then(|n| pairs.checked_add(n))
                .ok_or_else(count_overflow)?;
            largest = largest.max(bank.len());
        }
        let probes = spectrum.locations().len();
        let source_count = spectrum.source_coordinates().len();
        let transfer_products = probes
            .checked_mul(source_count)
            .and_then(|n| n.checked_mul(lines))
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(count_overflow)?;
        let probe_pairs = probes.checked_mul(pairs).ok_or_else(count_overflow)?;
        if lines > MAX_FIELD_SPECTRAL_LINES
            || lines > limits.maximum_lines
            || transfer_products > MAX_SPECTRAL_PRODUCTS
            || transfer_products > limits.maximum_transfer_products
            || probe_pairs > MAX_SPECTRAL_PRODUCTS
            || probe_pairs > limits.maximum_probe_frequency_pairs
        {
            return Err(bad("network spectral observation work budget exceeded"));
        }
        let per_line = probes
            .checked_mul(core::mem::size_of::<[Dd; 4]>())
            .and_then(|n| n.checked_add(core::mem::size_of::<SpectralLine>()))
            .ok_or_else(count_overflow)?;
        let workspace_bytes = lines
            .checked_mul(per_line)
            .and_then(|n| {
                n.checked_add(
                    factors
                        .len()
                        .checked_mul(core::mem::size_of::<SpectralBank>())?,
                )
            })
            .and_then(|n| n.checked_add(largest.checked_mul(core::mem::size_of::<usize>())?))
            .and_then(|n| {
                n.checked_add(probes.checked_mul(
                    core::mem::size_of::<[Dd; 5]>() + core::mem::size_of::<NetworkFieldMoments>(),
                )?)
            })
            .ok_or_else(count_overflow)?;
        if workspace_bytes > limits.maximum_workspace_bytes {
            return Err(bad(
                "network spectral observation workspace budget exceeded",
            ));
        }
        for bank in core::iter::once(mean).chain(factors.iter().map(|f| f.lines)) {
            for (index, line) in bank.iter().enumerate() {
                poll(&mut cancelled, index)?;
                let field = spectrum
                    .frequency_fields(line.frequency_index)
                    .ok_or_else(|| bad("network spectral line frequency index is out of range"))?;
                field.validate_excitation(line.excitation)?;
            }
        }
        let mean = prepare_bank(
            spectrum,
            mean,
            NetworkFieldFactorStatistics::Real,
            &mut cancelled,
        )?;
        let mut prepared = storage(factors.len())?;
        for factor in factors {
            prepared.push(prepare_bank(
                spectrum,
                factor.lines,
                factor.statistics,
                &mut cancelled,
            )?);
        }
        let actual_pairs = core::iter::once(&mean)
            .chain(prepared.iter())
            .map(|b| b.lines.len() * b.lines.len())
            .sum::<usize>()
            * probes;
        poll(&mut cancelled, lines)?;
        Ok(Self {
            spectrum,
            mean,
            factors: prepared,
            workspace_bytes,
            probe_frequency_pairs: actual_pairs,
        })
    }

    pub fn spectrum(&self) -> &PreparedNetworkFieldSpectrum {
        self.spectrum
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn factor_count(&self) -> usize {
        self.factors.len()
    }
    pub fn probe_frequency_pairs(&self) -> usize {
        self.probe_frequency_pairs
    }

    /// Means, uncentered squares and signed p*U over the actual window.
    /// Both beat and carrier-sum terms survive for the deterministic/real banks.
    /// Properness alone removes a proper factor's sum terms. No power or RMS
    /// values are combined before coherent fields.
    ///
    /// The output power includes reversible storage flow; it is not local
    /// heat. Nonuniform profiles give detector-weighted averages, not an
    /// unweighted work integral when multiplied by duration.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedNetworkFieldExposure<'_>, DuctError> {
        poll(&mut cancelled, 0)?;
        let probes = self.spectrum.locations().len();
        let mut sums = storage(probes)?;
        sums.resize(probes, [Dd::ZERO; 5]);
        for line in &self.mean.lines {
            poll(&mut cancelled, 0)?;
            let dc = window
                .product_kernel(line.omega, 0.0)
                .map_err(window_error)?;
            for (index, (sum, f)) in sums.iter_mut().zip(&line.fields).enumerate() {
                if index % 64 == 0 {
                    poll(&mut cancelled, index)?;
                }
                sum[0] = sum[0]
                    + dc.mean_real_product_parts_dd([f[0], f[1]], [Dd::ONE, Dd::ZERO])
                        .map_err(window_error)?;
                sum[1] = sum[1]
                    + dc.mean_real_product_parts_dd([f[2], f[3]], [Dd::ONE, Dd::ZERO])
                        .map_err(window_error)?;
            }
        }
        for bank in core::iter::once(&self.mean).chain(self.factors.iter()) {
            for (i, a) in bank.lines.iter().enumerate() {
                for (j, b) in bank.lines[..=i].iter().enumerate() {
                    poll(&mut cancelled, i)?;
                    let kernel = match bank.statistics {
                        NetworkFieldFactorStatistics::Real => PairKernel::Real(
                            window
                                .product_kernel(a.omega, b.omega)
                                .map_err(window_error)?,
                        ),
                        NetworkFieldFactorStatistics::ProperComplex => PairKernel::Proper(
                            window
                                .proper_product_kernel(a.omega, b.omega)
                                .map_err(window_error)?,
                        ),
                    };
                    for (index, ((sum, a), b)) in
                        sums.iter_mut().zip(&a.fields).zip(&b.fields).enumerate()
                    {
                        if index % 64 == 0 {
                            poll(&mut cancelled, index)?;
                        }
                        add_pair(sum, a, b, kernel, i != j)?;
                    }
                }
            }
        }
        let mut moments = storage(probes)?;
        for (index, (sum, location)) in sums.into_iter().zip(self.spectrum.locations()).enumerate()
        {
            poll(&mut cancelled, index)?;
            if sum.iter().any(|v| !v.hi.is_finite() || !v.lo.is_finite()) {
                return Err(DuctError::Singular);
            }
            let [mp, mq, pp, qq, pq] = sum.map(Dd::to_f64);
            if pp < 0.0 || qq < 0.0 {
                return Err(bad("network spectral self moment became negative"));
            }
            let area = location.area_m2;
            let out = NetworkFieldMoments {
                mean_pressure_pa: mp,
                mean_volume_flow_m3_s: mq,
                mean_particle_velocity_m_s: mq / area,
                mean_square_pressure_pa2: pp,
                mean_square_volume_flow_m6_s2: qq,
                mean_square_particle_velocity_m2_s2: (qq / area) / area,
                mean_power_w: pq,
                mean_intensity_w_m2: pq / area,
            };
            if [
                mp,
                mq,
                pp,
                qq,
                pq,
                out.mean_particle_velocity_m_s,
                out.mean_square_particle_velocity_m2_s2,
                out.mean_intensity_w_m2,
            ]
            .iter()
            .any(|x| !x.is_finite())
            {
                return Err(DuctError::Singular);
            }
            moments.push(out);
        }
        poll(&mut cancelled, probes)?;
        Ok(WindowedNetworkFieldExposure {
            window,
            spectrum: self.spectrum,
            moments,
        })
    }
}

fn prepare_bank(
    spectrum: &PreparedNetworkFieldSpectrum,
    input: &[NetworkFieldSpectralLine<'_>],
    statistics: NetworkFieldFactorStatistics,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpectralBank, DuctError> {
    let mut indices = storage(input.len())?;
    indices.extend(0..input.len());
    indices.sort_unstable_by_key(|&i| (input[i].frequency_index, i));
    let mut lines: Vec<SpectralLine> = storage(input.len())?;
    for index in indices {
        poll(cancelled, index)?;
        let line = input[index];
        let owner = &spectrum.fields[line.frequency_index];
        let omega = owner.port_response().omega();
        if lines
            .last()
            .is_none_or(|last| last.omega.to_bits() != omega.to_bits())
        {
            let mut fields = storage(spectrum.locations().len())?;
            fields.resize(spectrum.locations().len(), [Dd::ZERO; 4]);
            lines.push(SpectralLine { omega, fields });
        }
        for (probe, target) in lines
            .last_mut()
            .expect("frequency group exists")
            .fields
            .iter_mut()
            .enumerate()
        {
            poll(cancelled, probe)?;
            let [p, q] = owner.parts(probe, line.excitation)?;
            for (sum, value) in target.iter_mut().zip([p[0], p[1], q[0], q[1]]) {
                *sum = *sum + value;
                if !sum.hi.is_finite() || !sum.lo.is_finite() {
                    return Err(DuctError::Singular);
                }
            }
        }
    }
    Ok(SpectralBank { lines, statistics })
}

#[derive(Clone, Copy)]
enum PairKernel {
    Real(HarmonicProduct),
    Proper(ProperHarmonicProduct),
}
impl PairKernel {
    fn product(self, a: [Dd; 2], b: [Dd; 2]) -> Result<Dd, DuctError> {
        match self {
            Self::Real(k) => k.mean_real_product_parts_dd(a, b),
            Self::Proper(k) => k.mean_real_product_parts_dd(a, b),
        }
        .map_err(window_error)
    }
}

fn add_pair(
    sum: &mut [Dd; 5],
    a: &[Dd; 4],
    b: &[Dd; 4],
    kernel: PairKernel,
    distinct: bool,
) -> Result<(), DuctError> {
    let (pa, qa, pb, qb) = ([a[0], a[1]], [a[2], a[3]], [b[0], b[1]], [b[2], b[3]]);
    let weight = Dd::from_f64(if distinct { 2.0 } else { 1.0 });
    sum[2] = sum[2] + weight * kernel.product(pa, pb)?;
    sum[3] = sum[3] + weight * kernel.product(qa, qb)?;
    let mut power = kernel.product(pa, qb)?;
    if distinct {
        power = power + kernel.product(qa, pb)?;
    }
    sum[4] = sum[4] + power;
    Ok(())
}

fn storage<T>(n: usize) -> Result<Vec<T>, DuctError> {
    let mut out = Vec::new();
    out.try_reserve_exact(n)
        .map_err(|_| bad("network spectral allocation failed"))?;
    Ok(out)
}
fn count_overflow() -> DuctError {
    bad("network spectral work or storage count overflow")
}
