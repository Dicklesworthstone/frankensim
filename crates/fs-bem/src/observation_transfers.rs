//! Exterior observation realization, never an acoustic load. The producer
//! supplies exp(-iwt) BEM fields; conjugation is performed exactly once at the
//! DFT boundary. A declared C4 spectral taper and explicit reconstruction
//! latency make the sampled observation finite-band. Neither is hidden in a
//! fitted allpass or in a passivity repair. Accuracy is measured, not inferred
//! from FFT size. Spatial/producer uncertainty remains separate.
use crate::{
    helmholtz::{self, HelmholtzError, Medium, laplace},
    near_field::Geometry,
    panel3d::SpherePanels,
    radiation_policy::GeometryPolicy,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_fft::{C64 as FftComplex, Fft};
use fs_math::{c64::C64, det};
use fs_vfit::observation::{
    FirError, FirKernel, MAX_OBSERVATION_ENTRIES, MAX_OBSERVATION_INPUTS, MAX_OBSERVATION_OUTPUTS,
    ObservationFir,
};
use std::collections::BTreeMap;

mod reciprocal;
pub use reciprocal::ReciprocalBemObservation;
mod scattering;
pub use scattering::{PointSourceSpec, ScatteringObservation};
mod wall;
pub use wall::{
    BoundaryObservationBake, LocalBoundary, LocalBoundaryObservation, LocalWall, LowerBandWindow,
};

/// Largest admitted DFT/CQ size. Work budgets may be stricter.
pub const MAX_OBSERVATION_FFT: usize = 65_536;
/// At most eight fixed, doubling ladder rungs; at least four are required.
pub const MAX_OBSERVATION_RUNGS: usize = 8;

/// Physical source coordinate. Velocity traces are supplied separately, in
/// m/s per unit coordinate rate. Aperture traces integrate to exactly one m³/s
/// per unit flow, checked against the physical panel areas.
/// Point-volume-flow coordinates instead belong to [`ScatteringObservation`]
/// and are never admitted as panel traces.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SourceId {
    StructuralMode(String),
    ApertureVolumeFlow(String),
    /// Prescribed point-source outward volume flow [m³/s], not a panel trace.
    PointVolumeFlow(String),
}
impl SourceId {
    fn name(&self) -> &str {
        match self {
            Self::StructuralMode(s) | Self::ApertureVolumeFlow(s) | Self::PointVolumeFlow(s) => s,
        }
    }
}
/// All four components are retained, in order p, dp/dx, dp/dy, dp/dz.
#[derive(Debug, Clone, PartialEq)]
pub struct ReceiverSpec {
    pub id: String,
    pub position_m: [f64; 3],
}

/// Producer error is separate from FIR convergence. An unresolved spatial
/// error is preserved explicitly, never replaced with a guessed percentage.
/// Such a bake qualifies only realization of its declared discrete operator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProducerUncertainty {
    ExactAnalytic,
    RelativeEstimate {
        relative: f64,
        source_hash: ContentHash,
    },
    Unresolved {
        source_hash: ContentHash,
    },
}

/// Required declaration. No atmosphere, source, receiver, uncertainty, band,
/// latency or convergence control is defaulted. The geometric delays are
/// source-major, receiver-minor, conservative earliest arrivals [s].
#[derive(Debug, Clone)]
pub struct TransferDeclaration {
    pub producer_hash: ContentHash,
    pub sources: Vec<SourceId>,
    pub receivers: Vec<ReceiverSpec>,
    pub fs_hz: f64,
    pub geometric_delays_s: Vec<f64>,
    pub pass_band_upper_hz: f64,
    pub taper_stop_hz: f64,
    pub filter_latency_samples: usize,
    pub fft_lengths: Vec<usize>,
    pub tail_energy_tolerance: f64,
    pub dtft_relative_tolerance: f64,
    pub producer_uncertainty: ProducerUncertainty,
    pub maximum_frequency_entries: usize,
}

/// No partial operator is returned by a failed bake or CQ run.
#[derive(Debug, Clone, PartialEq)]
pub enum TransferError {
    Invalid(&'static str),
    Helmholtz(HelmholtzError),
    Material(fs_material::MaterialError),
    Mechanical(fs_phs::PhsError),
    Fir(FirError),
    Cap { requested: usize, cap: usize },
    NotConverged { tail: f64, wrap: f64, dtft: f64 },
    Cancelled,
}
impl From<HelmholtzError> for TransferError {
    fn from(e: HelmholtzError) -> Self {
        Self::Helmholtz(e)
    }
}
impl From<fs_material::MaterialError> for TransferError {
    fn from(error: fs_material::MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<fs_phs::PhsError> for TransferError {
    fn from(error: fs_phs::PhsError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<FirError> for TransferError {
    fn from(e: FirError) -> Self {
        Self::Fir(e)
    }
}
impl core::fmt::Display for TransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "observation transfer: {self:?}")
    }
}
impl std::error::Error for TransferError {}

/// One computed DFT rung; fractions are the worst physical quantity, never a
/// mixture of Pa and Pa/m. Tail is energy after N/4; wrap is energy in the
/// second half. The Cartesian gradient is ONE vector quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct LadderRung {
    pub n_fft: usize,
    pub tail_energy_fraction: f64,
    pub wrap_energy_fraction: f64,
    pub dtft_relative_change: Option<f64>,
}
/// Portable data for fs-vfit/receiver consumers. Every source and every p/grad
/// component is present. This carries Estimate-only numerical information; no
/// passivity or full spatial-convergence authority is minted.
#[derive(Debug, Clone)]
pub struct ObservationTransferBake {
    pub declaration_hash: ContentHash,
    pub sources: Vec<SourceId>,
    pub receivers: Vec<ReceiverSpec>,
    pub fs_hz: f64,
    pub geometric_delays_s: Vec<f64>,
    pub filter_latency_samples: usize,
    pub pass_band_upper_hz: f64,
    pub taper_stop_hz: f64,
    pub producer_uncertainty: ProducerUncertainty,
    pub kernels: Vec<FirKernel>,
    pub ladder: Vec<LadderRung>,
}
impl ObservationTransferBake {
    /// Real causal replay outside the mechanics/oscillation loop.
    pub fn observer(&self) -> Result<ObservationFir, TransferError> {
        Ok(ObservationFir::new(
            self.sources.len(),
            4 * self.receivers.len(),
            self.kernels.clone(),
        )?)
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn phase(x: f64) -> C64 {
    C64::new(det::cos(x), det::sin(x))
}
fn cap(requested: usize, limit: usize) -> Result<(), TransferError> {
    if requested > limit {
        Err(TransferError::Cap {
            requested,
            cap: limit,
        })
    } else {
        Ok(())
    }
}
/// Cold FFT-size cap admission; rejects cap+1 before any allocation/solve.
pub fn admit_observation_fft(n: usize) -> Result<(), TransferError> {
    cap(n, MAX_OBSERVATION_FFT)?;
    if n < 32 || !n.is_power_of_two() {
        return Err(TransferError::Invalid(
            "DFT size must be a power of two >=32",
        ));
    }
    Ok(())
}
impl TransferDeclaration {
    /// Allocation-free declaration/cap admission; returns the complete number
    /// of source-to-receiver-component kernels, without solving or truncating.
    pub fn admit_work(&self) -> Result<usize, TransferError> {
        cap(self.sources.len(), MAX_OBSERVATION_INPUTS)?;
        cap(
            self.receivers
                .len()
                .checked_mul(4)
                .ok_or(TransferError::Invalid("receiver shape overflow"))?,
            MAX_OBSERVATION_OUTPUTS,
        )?;
        cap(self.fft_lengths.len(), MAX_OBSERVATION_RUNGS)?;
        let pairs = self
            .sources
            .len()
            .checked_mul(self.receivers.len())
            .ok_or(TransferError::Invalid("source shape overflow"))?;
        if pairs == 0
            || self.geometric_delays_s.len() != pairs
            || self.producer_hash == ContentHash([0; 32])
            || self
                .sources
                .iter()
                .any(|s| s.name().is_empty() || s.name().len() > 1024)
            || self.receivers.iter().any(|r| {
                r.id.is_empty() || r.id.len() > 1024 || r.position_m.iter().any(|x| !x.is_finite())
            })
            || !self.fs_hz.is_finite()
            || self.fs_hz <= 0.0
            || !self.pass_band_upper_hz.is_finite()
            || self.pass_band_upper_hz <= 0.0
            || !self.taper_stop_hz.is_finite()
            || self.taper_stop_hz <= self.pass_band_upper_hz
            || self.taper_stop_hz >= 0.5 * self.fs_hz
            || self.geometric_delays_s.iter().any(|x| {
                !x.is_finite() || *x < 0.0 || *x * self.fs_hz > MAX_OBSERVATION_ENTRIES as f64
            })
            || self.filter_latency_samples == 0
            || self.filter_latency_samples > MAX_OBSERVATION_FFT / 4
            || self.fft_lengths.len() < 4
            || self
                .fft_lengths
                .windows(2)
                .any(|n| n[0].checked_mul(2) != Some(n[1]))
            || !self.tail_energy_tolerance.is_finite()
            || !(0.0..=1e-10).contains(&self.tail_energy_tolerance)
            || self.tail_energy_tolerance == 0.0
            || !self.dtft_relative_tolerance.is_finite()
            || !(0.0..=1e-4).contains(&self.dtft_relative_tolerance)
            || self.dtft_relative_tolerance == 0.0
            || match self.producer_uncertainty {
                ProducerUncertainty::ExactAnalytic => false,
                ProducerUncertainty::RelativeEstimate {
                    relative,
                    source_hash,
                } => !relative.is_finite() || relative < 0.0 || source_hash == ContentHash([0; 32]),
                ProducerUncertainty::Unresolved { source_hash } => {
                    source_hash == ContentHash([0; 32])
                }
            }
        {
            return Err(TransferError::Invalid(
                "incomplete or invalid transfer declaration",
            ));
        }
        for &n in &self.fft_lengths {
            admit_observation_fft(n)?;
        }
        let n = *self
            .fft_lengths
            .last()
            .ok_or(TransferError::Invalid("empty DFT ladder"))?;
        let entries = pairs
            .checked_mul(4)
            .and_then(|v| v.checked_mul(n / 2 + 1))
            .ok_or(TransferError::Invalid("frequency shape overflow"))?;
        cap(
            entries,
            self.maximum_frequency_entries.min(MAX_OBSERVATION_ENTRIES),
        )?;
        let history = n / 2
            + self
                .geometric_delays_s
                .iter()
                .map(|d| (d * self.fs_hz).floor() as usize)
                .max()
                .unwrap_or(0);
        cap(
            history
                .checked_mul(self.sources.len())
                .ok_or(TransferError::Invalid("history shape overflow"))?,
            MAX_OBSERVATION_ENTRIES,
        )?;
        Ok(pairs * 4)
    }
    /// Hash every explicit numerical/physical declaration, including source
    /// coordinate kind, receiver components and the producer's identity.
    pub fn declaration_hash(&self) -> ContentHash {
        let mut h = DomainHasher::new("fs-bem:observation-transfer:v1");
        h.update(self.producer_hash.as_bytes());
        for n in [
            self.sources.len(),
            self.receivers.len(),
            self.fft_lengths.len(),
        ] {
            h.update(&(n as u64).to_le_bytes());
        }
        for s in &self.sources {
            h.update(&[match s {
                SourceId::StructuralMode(_) => 0,
                SourceId::ApertureVolumeFlow(_) => 1,
                SourceId::PointVolumeFlow(_) => 2,
            }]);
            h.update(&(s.name().len() as u64).to_le_bytes());
            h.update(s.name().as_bytes());
        }
        for r in &self.receivers {
            h.update(&(r.id.len() as u64).to_le_bytes());
            h.update(r.id.as_bytes());
            for x in r.position_m {
                h.update(&x.to_bits().to_le_bytes());
            }
        }
        for x in [
            self.fs_hz,
            self.pass_band_upper_hz,
            self.taper_stop_hz,
            self.tail_energy_tolerance,
            self.dtft_relative_tolerance,
        ]
        .into_iter()
        .chain(self.geometric_delays_s.iter().copied())
        {
            h.update(&x.to_bits().to_le_bytes());
        }
        match self.producer_uncertainty {
            ProducerUncertainty::ExactAnalytic => h.update(&[0]),
            ProducerUncertainty::RelativeEstimate {
                relative,
                source_hash,
            } => {
                h.update(&[1]);
                h.update(&relative.to_bits().to_le_bytes());
                h.update(source_hash.as_bytes());
            }
            ProducerUncertainty::Unresolved { source_hash } => {
                h.update(&[2]);
                h.update(source_hash.as_bytes());
            }
        }
        for n in [self.filter_latency_samples, self.maximum_frequency_entries]
            .into_iter()
            .chain(self.fft_lengths.iter().copied())
        {
            h.update(&(n as u64).to_le_bytes());
        }
        h.finalize()
    }
}

/// C4 transition, exactly unity throughout the declared pass band and exactly
/// zero before Nyquist. This is a declared reconstruction window, not a claim
/// about the radiator. Its explicit latency is retained in the output data.
fn observation_taper(f_hz: f64, pass_hz: f64, stop_hz: f64) -> f64 {
    if f_hz <= pass_hz {
        return 1.0;
    }
    if f_hz >= stop_hz {
        return 0.0;
    }
    let x = (f_hz - pass_hz) / (stop_hz - pass_hz);
    // Evaluate the upper half by symmetry to avoid cancellation near x=1.
    if x <= 0.5 {
        1.0 - observation_transition_low(x)
    } else {
        observation_transition_low(1.0 - x)
    }
}

fn observation_transition_low(t: f64) -> f64 {
    t.powi(5) * (126.0 + t * (-420.0 + t * (540.0 + t * (-315.0 + 70.0 * t))))
}

/// Channel ranges of the physical quantities, in kernel order p, dp/dx, dp/dy,
/// dp/dz per source/receiver pair: pressure [Pa] alone, and the Cartesian
/// gradient [Pa/m] as one vector. A component's fraction of its own energy is
/// frame-dependent, and an on-axis receiver's zero-by-symmetry components hold
/// only roundoff whose self-relative tail is ~3/4 at every resolution, so gates
/// measure each component against its whole vector.
fn quantities(channels: usize) -> impl Iterator<Item = core::ops::Range<usize>> {
    (0..channels / 4).flat_map(|pair| [4 * pair..4 * pair + 1, 4 * pair + 1..4 * pair + 4])
}

fn dtft_change(
    a: &[FirKernel],
    b: &[FirKernel],
    decl: &TransferDeclaration,
    lower_pass_hz: f64,
) -> Result<f64, TransferError> {
    let mut largest = 0.0_f64;
    for range in quantities(a.len().min(b.len())) {
        let mut samples = Vec::with_capacity(127 * range.len());
        let mut scale = 0.0_f64;
        for (a, b) in a[range.clone()].iter().zip(&b[range]) {
            for i in 0..127 {
                let w = if lower_pass_hz == 0.0 {
                    // Preserve the historical full-band arithmetic exactly.
                    core::f64::consts::TAU * decl.pass_band_upper_hz * (i as f64 + 0.371)
                        / (127.0 * decl.fs_hz)
                } else {
                    core::f64::consts::TAU
                        * (lower_pass_hz
                            + (decl.pass_band_upper_hz - lower_pass_hz) * (i as f64 + 0.371)
                                / 127.0)
                        / decl.fs_hz
                };
                let av = a.transfer_at(w)?;
                let bv = b.transfer_at(w)?;
                scale = scale
                    .max(av.re.abs())
                    .max(av.im.abs())
                    .max(bv.re.abs())
                    .max(bv.im.abs());
                samples.push((av, bv));
            }
        }
        // Normalize before subtracting/squaring: finite large/small physical
        // units must not turn a convergence discrepancy into inf/inf or 0/0.
        let (mut error, mut norm) = (0.0, 0.0);
        if scale > 0.0 {
            for (av, bv) in samples {
                let av = C64::new(av.re / scale, av.im / scale);
                let bv = C64::new(bv.re / scale, bv.im / scale);
                error += (av - bv).norm_sq();
                norm += bv.norm_sq();
            }
        }
        let relative = if norm == 0.0 {
            if error == 0.0 { 0.0 } else { f64::INFINITY }
        } else {
            det::sqrt(error / norm)
        };
        largest = largest.max(relative);
    }
    Ok(largest)
}

/// Bake complete source-major, receiver-major p/grad fields from the existing
/// frequency owner on nested dense DFT grids. Callback convention is exp(-iwt).
/// Every rung is retained; there is no early success before four rungs.
/// # Errors
/// Declaration/cap, producer/shape/nonfinite, cancellation, or a final tail,
/// wrap-around or inter-rung DTFT discrepancy above the declared gates.
pub fn bake_observation_transfers(
    mut solve_at: impl FnMut(f64) -> Result<Vec<C64>, HelmholtzError>,
    decl: &TransferDeclaration,
    cancelled: impl FnMut() -> bool,
) -> Result<ObservationTransferBake, TransferError> {
    bake_observation_transfers_with_errors(
        |omega| solve_at(omega).map_err(Into::into),
        decl,
        0.0,
        cancelled,
    )
}

// The same numerical owner, allowing native constitutive refusals to retain
// their types. The wall adapter admits its explicit nonzero passband first.
fn bake_observation_transfers_with_errors(
    mut solve_at: impl FnMut(f64) -> Result<Vec<C64>, TransferError>,
    decl: &TransferDeclaration,
    lower_pass_hz: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ObservationTransferBake, TransferError> {
    let channels = decl.admit_work()?;
    let finest = *decl
        .fft_lengths
        .last()
        .ok_or(TransferError::Invalid("empty DFT ladder"))?;
    let mut cache = BTreeMap::<usize, Vec<C64>>::new();
    let mut previous: Option<Vec<FirKernel>> = None;
    let mut ladder = Vec::new();
    for &n in &decl.fft_lengths {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        let mut spectra = vec![vec![FftComplex::default(); n]; channels];
        for bin in 0..=n / 2 {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let f = decl.fs_hz * bin as f64 / n as f64;
            let window = observation_taper(f, decl.pass_band_upper_hz, decl.taper_stop_hz);
            if window == 0.0 {
                continue;
            }
            let key = bin * (finest / n);
            if let std::collections::btree_map::Entry::Vacant(entry) = cache.entry(key) {
                // Static velocity produces no exterior acoustic pressure or
                // pressure gradient: q=-rho*s*v vanishes at s=0, p(infty)=0.
                let row = if bin == 0 {
                    vec![C64::ZERO; channels]
                } else {
                    solve_at(core::f64::consts::TAU * f)?
                };
                if row.len() != channels || row.iter().any(|&z| !finite(z)) {
                    return Err(TransferError::Invalid("producer field shape/nonfinite"));
                }
                entry.insert(row);
            }
            let row = &cache[&key];
            for c in 0..channels {
                let pair = c / 4;
                let integer = (decl.geometric_delays_s[pair] * decl.fs_hz).floor();
                let shift = integer - decl.filter_latency_samples as f64;
                let h = row[c].conj() * phase(core::f64::consts::TAU * f * shift / decl.fs_hz);
                spectra[c][bin] = FftComplex::new(h.re * window, h.im * window);
                if bin > 0 && bin < n / 2 {
                    spectra[c][n - bin] = spectra[c][bin].conj();
                }
            }
        }
        let fft = Fft::new(n);
        let mut scratch = vec![FftComplex::default(); n];
        let mut kernels = Vec::with_capacity(channels);
        for (c, spectrum) in spectra.iter_mut().enumerate() {
            fft.inverse(spectrum, &mut scratch);
            if spectrum
                .iter()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
            {
                return Err(TransferError::Invalid("nonfinite IFFT kernel"));
            }
            let scale = spectrum
                .iter()
                .fold(0.0_f64, |s, v| s.max(v.re.abs()).max(v.im.abs()));
            let energy = |v: &FftComplex| {
                if scale == 0.0 {
                    0.0
                } else {
                    (v.re / scale).powi(2) + (v.im / scale).powi(2)
                }
            };
            let total = spectrum.iter().map(energy).sum::<f64>();
            let imaginary = if scale == 0.0 {
                0.0
            } else {
                spectrum.iter().map(|v| (v.im / scale).powi(2)).sum::<f64>()
            };
            if imaginary > 1e-24 * total {
                return Err(TransferError::Invalid("nonreal observation kernel"));
            }
            let pair = c / 4;
            kernels.push(FirKernel {
                input: pair / decl.receivers.len(),
                output: (pair % decl.receivers.len()) * 4 + c % 4,
                integer_delay_samples: (decl.geometric_delays_s[pair] * decl.fs_hz).floor()
                    as usize,
                taps: spectrum[..n / 2].iter().map(|v| v.re).collect(),
            });
        }
        let (mut tail, mut wrap) = (0.0_f64, 0.0_f64);
        for range in quantities(channels) {
            let group = &spectra[range];
            let scale = group
                .iter()
                .flatten()
                .fold(0.0_f64, |s, v| s.max(v.re.abs()).max(v.im.abs()));
            if scale == 0.0 {
                continue;
            }
            let energy = |v: &FftComplex| (v.re / scale).powi(2) + (v.im / scale).powi(2);
            let part = |from: usize| {
                group
                    .iter()
                    .flat_map(|s| &s[from..])
                    .map(energy)
                    .sum::<f64>()
            };
            let total = part(0);
            tail = tail.max(part(n / 4) / total);
            wrap = wrap.max(part(n / 2) / total);
        }
        let change = previous
            .as_ref()
            .map(|old| dtft_change(old, &kernels, decl, lower_pass_hz))
            .transpose()?;
        ladder.push(LadderRung {
            n_fft: n,
            tail_energy_fraction: tail,
            wrap_energy_fraction: wrap,
            dtft_relative_change: change,
        });
        previous = Some(kernels);
    }
    let last = ladder
        .last()
        .ok_or(TransferError::Invalid("empty ladder"))?;
    let change = last.dtft_relative_change.unwrap_or(f64::INFINITY);
    if last.tail_energy_fraction > decl.tail_energy_tolerance
        || last.wrap_energy_fraction >= 1e-10
        || change >= decl.dtft_relative_tolerance
    {
        return Err(TransferError::NotConverged {
            tail: last.tail_energy_fraction,
            wrap: last.wrap_energy_fraction,
            dtft: change,
        });
    }
    if cancelled() {
        return Err(TransferError::Cancelled);
    }
    Ok(ObservationTransferBake {
        declaration_hash: decl.declaration_hash(),
        sources: decl.sources.clone(),
        receivers: decl.receivers.clone(),
        fs_hz: decl.fs_hz,
        geometric_delays_s: decl.geometric_delays_s.clone(),
        filter_latency_samples: decl.filter_latency_samples,
        pass_band_upper_hz: decl.pass_band_upper_hz,
        taper_stop_hz: decl.taper_stop_hz,
        producer_uncertainty: decl.producer_uncertainty,
        kernels: previous.ok_or(TransferError::Invalid("empty kernel ladder"))?,
        ladder,
    })
}

/// A single source/geometry adapter for both dense Fourier samples and CQ.
/// Real frequency samples can also feed the existing broadband rational owner;
/// this adapter does not contain a second fit or a special instrument model.
pub struct BemObservation<'a> {
    surface: &'a SpherePanels,
    fields: Vec<Vec<C64>>,
    points: Vec<[f64; 3]>,
    medium: Medium,
    minimum_clearance_m: f64,
    clearances_s: Vec<f64>,
    sources: Vec<SourceId>,
    receivers: Vec<ReceiverSpec>,
    policy: Option<GeometryPolicy<'a>>,
    panel_rule: helmholtz::PanelRule,
}
impl<'a> BemObservation<'a> {
    /// Bind all physical traces to admitted exterior receivers. Unit aperture
    /// volume flow is integrated with the SAME panel areas used by the solver.
    pub fn new(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[SourceId],
        traces: &[Vec<f64>],
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
    ) -> Result<Self, TransferError> {
        Self::new_admitted(
            surface,
            medium,
            sources,
            traces,
            receivers,
            minimum_clearance_m,
            false,
        )
    }

    fn new_admitted(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[SourceId],
        traces: &[Vec<f64>],
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        allow_zero_structural_trace: bool,
    ) -> Result<Self, TransferError> {
        if sources.is_empty()
            || sources.len() != traces.len()
            || sources.len() > MAX_OBSERVATION_INPUTS
            || receivers.iter().any(|r| r.id.is_empty())
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
        {
            return Err(TransferError::Invalid(
                "invalid complete BEM source declaration",
            ));
        }
        let points: Vec<_> = receivers.iter().map(|r| r.position_m).collect();
        let geometry = Geometry::new(surface, &points, minimum_clearance_m)?;
        for (source, trace) in sources.iter().zip(traces) {
            if matches!(source, SourceId::PointVolumeFlow(_))
                || source.name().is_empty()
                || trace.len() != surface.areas().len()
                || trace.iter().any(|v| !v.is_finite())
                || (trace.iter().all(|&v| v == 0.0)
                    && !(allow_zero_structural_trace
                        && matches!(source, SourceId::StructuralMode(_))))
            {
                return Err(TransferError::Invalid("missing/nonfinite source trace"));
            }
            if matches!(source, SourceId::ApertureVolumeFlow(_)) {
                let flow = trace
                    .iter()
                    .zip(surface.areas())
                    .map(|(v, a)| v * a)
                    .sum::<f64>();
                if (flow - 1.0).abs() > 1e-12 {
                    return Err(TransferError::Invalid(
                        "aperture trace must integrate to unit SI volume flow",
                    ));
                }
            }
        }
        let clearances_s = geometry
            .clearances_m()
            .iter()
            .map(|d| d / medium.sound_speed)
            .collect();
        Ok(Self {
            surface,
            fields: traces
                .iter()
                .map(|f| f.iter().copied().map(C64::from_re).collect())
                .collect(),
            points,
            medium,
            minimum_clearance_m,
            clearances_s,
            sources: sources.to_vec(),
            receivers: receivers.to_vec(),
            policy: None,
            panel_rule: helmholtz::PanelRule::Centroid,
        })
    }
    /// Bind the existing geometry-selected boundary owner and an explicit panel
    /// rule. This is the same route used by closed-body modal loading: triangle
    /// CBIE below the geometric exclusion bound, then ruled Burton--Miller.
    /// Receiver pressure/gradient retain the existing centroid Green rows.
    /// Their spatial uncertainty is not removed by boundary quadrature. A
    /// structural mode with an identically zero normal trace is retained as a
    /// zero observation column; aperture flow must still integrate to one.
    ///
    /// The legacy constructor, its operator hash and CQ behavior are unchanged.
    /// This real-frequency route has no matching Laplace continuation; CQ is
    /// refused instead of silently substituting the legacy boundary operator.
    pub fn new_routed(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[SourceId],
        traces: &[Vec<f64>],
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        panel_rule: helmholtz::PanelRule,
    ) -> Result<Self, TransferError> {
        if matches!(
            panel_rule,
            helmholtz::PanelRule::TriangleNear { max_pairs: 0 }
                | helmholtz::PanelRule::Triangle { max_pairs: 0 }
        ) {
            return Err(TransferError::Invalid(
                "empty observation boundary quadrature allowance",
            ));
        }
        let policy = GeometryPolicy::new(surface)?;
        let mut owner = Self::new_admitted(
            surface,
            medium,
            sources,
            traces,
            receivers,
            minimum_clearance_m,
            true,
        )?;
        owner.policy = Some(policy);
        owner.panel_rule = panel_rule;
        Ok(owner)
    }
    /// Bind ONE ruled Burton--Miller boundary operator at every frequency. A
    /// transfer bake needs one discrete operator across its band: the routed
    /// owner switches from triangle CBIE to Burton--Miller at the geometric
    /// bound, and the two discretizations differ there by their own errors, a
    /// step in every transfer that no DFT ladder converges. The panel rule must
    /// integrate triangles. Otherwise as [`Self::new_routed`]: the same closed
    /// outward admission, zero structural traces stay columns, receiver rows
    /// keep their centroid Green rows, and CQ is refused.
    pub fn new_ruled_burton_miller(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[SourceId],
        traces: &[Vec<f64>],
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        panel_rule: helmholtz::PanelRule,
    ) -> Result<Self, TransferError> {
        if matches!(
            panel_rule,
            helmholtz::PanelRule::Centroid
                | helmholtz::PanelRule::TriangleNear { max_pairs: 0 }
                | helmholtz::PanelRule::Triangle { max_pairs: 0 }
        ) {
            return Err(TransferError::Invalid(
                "ruled Burton-Miller observation needs a nonempty triangle rule",
            ));
        }
        GeometryPolicy::new(surface)?;
        let mut owner = Self::new_admitted(
            surface,
            medium,
            sources,
            traces,
            receivers,
            minimum_clearance_m,
            true,
        )?;
        owner.panel_rule = panel_rule;
        Ok(owner)
    }
    /// Only [`Self::new_ruled_burton_miller`] binds a triangle rule without a policy.
    fn ruled_burton_miller(&self) -> bool {
        self.policy.is_none() && !matches!(self.panel_rule, helmholtz::PanelRule::Centroid)
    }
    /// Bake this bound BEM owner. A cross-wired source, receiver, geometric
    /// delay or operator hash is refused before any frequency solve.
    pub fn bake(
        &self,
        decl: &TransferDeclaration,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ObservationTransferBake, TransferError> {
        if decl.sources != self.sources
            || decl.receivers != self.receivers
            || decl.geometric_delays_s != self.geometric_delays_s()
            || decl.producer_hash != self.operator_hash()
            || matches!(
                decl.producer_uncertainty,
                ProducerUncertainty::ExactAnalytic
            )
        {
            return Err(TransferError::Invalid(
                "cross-wired BEM transfer declaration",
            ));
        }
        bake_observation_transfers(|w| self.frequency_sample(w), decl, cancelled)
    }
    /// Conservative geometric delay for every complete source/receiver pair.
    pub fn geometric_delays_s(&self) -> Vec<f64> {
        (0..self.fields.len())
            .flat_map(|_| self.clearances_s.iter().copied())
            .collect()
    }
    /// Hash physical geometry, medium and every signed source trace. Compose
    /// this with the caller's material/source provenance in the declaration.
    pub fn operator_hash(&self) -> ContentHash {
        let mut h = DomainHasher::new(if self.policy.is_some() {
            "fs-bem:observation-boundary:routed:v1"
        } else if self.ruled_burton_miller() {
            "fs-bem:observation-boundary:ruled-burton-miller:v1"
        } else {
            "fs-bem:observation-boundary:v1"
        });
        h.update(&(self.surface.areas().len() as u64).to_le_bytes());
        h.update(&(self.fields.len() as u64).to_le_bytes());
        for value in [
            self.medium.density,
            self.medium.sound_speed,
            self.minimum_clearance_m,
        ]
        .into_iter()
        .chain(self.surface.centroids().iter().flatten().copied())
        .chain(self.surface.normals().iter().flatten().copied())
        .chain(self.surface.areas().iter().copied())
        .chain(self.fields.iter().flat_map(|row| row.iter().map(|z| z.re)))
        {
            h.update(&value.to_bits().to_le_bytes());
        }
        for triangle in self.surface.triangles().into_iter().flatten() {
            for &value in triangle.iter().flatten() {
                h.update(&value.to_bits().to_le_bytes());
            }
        }
        if self.policy.is_some() || self.ruled_burton_miller() {
            let (tag, pairs) = match self.panel_rule {
                helmholtz::PanelRule::Centroid => (0_u64, 0),
                helmholtz::PanelRule::TriangleNear { max_pairs } => (1, max_pairs),
                helmholtz::PanelRule::Triangle { max_pairs } => (2, max_pairs),
            };
            h.update(&tag.to_le_bytes());
            h.update(&(pairs as u64).to_le_bytes());
        }
        h.finalize()
    }
    fn flatten(fields: Vec<laplace::PointField>, output: &mut Vec<C64>) {
        for p in fields {
            output.push(p.pressure_pa);
            output.extend(p.gradient_pa_per_m);
        }
    }
    /// exp(-iwt) complete p/grad rows through the production radiation owner.
    pub fn frequency_sample(&self, omega: f64) -> Result<Vec<C64>, HelmholtzError> {
        if !omega.is_finite() || omega < 0.0 {
            return Err(HelmholtzError::BadParameter {
                what: "invalid observation frequency",
            });
        }
        if omega == 0.0 {
            return Ok(vec![C64::ZERO; 4 * self.points.len() * self.fields.len()]);
        }
        let k = omega / self.medium.sound_speed;
        let formulation = match &self.policy {
            Some(policy) => policy.formulation(k)?,
            None => helmholtz::Formulation::BurtonMiller,
        };
        let mut out = Vec::new();
        for batch in self
            .fields
            .chunks(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH)
        {
            let fields: Vec<_> = batch.iter().map(Vec::as_slice).collect();
            for solution in helmholtz::solve_radiation_batch_with_rule(
                self.surface,
                k,
                self.medium,
                &fields,
                formulation,
                self.panel_rule,
            )? {
                Self::flatten(
                    laplace::observe_frequency(
                        self.surface,
                        &solution,
                        &self.points,
                        self.minimum_clearance_m,
                    )?,
                    &mut out,
                );
            }
        }
        Ok(out)
    }
    /// CQ reference samples, with the same full source basis and receiver rows.
    /// A per-receiver advance may extract a known propagation prefix exactly.
    pub fn laplace_sample(&self, s: C64, advances_s: &[f64]) -> Result<Vec<C64>, HelmholtzError> {
        if self.policy.is_some() {
            return Err(HelmholtzError::BadParameter {
                what: "geometry-routed observation has no matching Laplace boundary operator",
            });
        }
        if self.ruled_burton_miller() {
            return Err(HelmholtzError::BadParameter {
                what: "ruled Burton-Miller observation has no matching Laplace boundary operator",
            });
        }
        let mut out = Vec::new();
        for batch in self
            .fields
            .chunks(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH)
        {
            let fields: Vec<_> = batch.iter().map(Vec::as_slice).collect();
            for solution in
                laplace::solve_radiation_laplace_batch(self.surface, s, self.medium, &fields)?
            {
                Self::flatten(
                    laplace::observe_laplace(
                        self.surface,
                        &solution,
                        &self.points,
                        self.minimum_clearance_m,
                        advances_s,
                    )?,
                    &mut out,
                );
            }
        }
        Ok(out)
    }
}

/// BDF2 CQ reference (Lubich, Numer. Math. 52, 1988). The generating function
/// is delta(zeta)=3/2-2*zeta+zeta²/2. Evaluate the original Laplace-domain BEM
/// operator at every contour point; never substitute a fitted FIR/rational
/// model. Input histories begin at rest. Smooth compatible onset is required
/// for the uncorrected method's formal second order.
///
/// `damping_exponent` means rho=exp(-exponent/N). The output is the first
/// input.len() time samples; N>=2*len prevents circular input wrap. CQ contour
/// alias/roundoff and time error remain separate refinement responsibilities.
pub fn cq_time_response(
    mut solve_at: impl FnMut(C64) -> Result<Vec<C64>, HelmholtzError>,
    input: &[Vec<f64>],
    outputs: usize,
    dt_s: f64,
    n_fft: usize,
    damping_exponent: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<Vec<f64>>, TransferError> {
    admit_observation_fft(n_fft)?;
    let steps = input.first().map_or(0, Vec::len);
    let entries = n_fft
        .checked_mul(input.len().max(outputs))
        .ok_or(TransferError::Invalid("CQ shape overflow"))?;
    cap(entries, MAX_OBSERVATION_ENTRIES)?;
    cap(input.len(), MAX_OBSERVATION_INPUTS)?;
    cap(outputs, MAX_OBSERVATION_OUTPUTS)?;
    if input.is_empty()
        || outputs == 0
        || steps == 0
        || steps > n_fft / 2
        || input
            .iter()
            .any(|r| r.len() != steps || r.iter().any(|v| !v.is_finite()))
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !damping_exponent.is_finite()
        || !(8.0..=24.0).contains(&damping_exponent)
    {
        return Err(TransferError::Invalid(
            "invalid CQ input, horizon or contour",
        ));
    }
    let rho = det::exp(-damping_exponent / n_fft as f64);
    let fft = Fft::new(n_fft);
    let mut scratch = vec![FftComplex::default(); n_fft];
    let mut source = Vec::with_capacity(input.len());
    for row in input {
        let mut data = vec![FftComplex::default(); n_fft];
        let mut weight = 1.0;
        for (x, &v) in data.iter_mut().zip(row) {
            x.re = v * weight;
            weight *= rho;
        }
        fft.forward(&mut data, &mut scratch);
        source.push(data);
    }
    let mut target = vec![vec![FftComplex::default(); n_fft]; outputs];
    for bin in 0..=n_fft / 2 {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        let z = phase(-core::f64::consts::TAU * bin as f64 / n_fft as f64).scale(rho);
        let s = (C64::from_re(1.5) - z.scale(2.0) + (z * z).scale(0.5)).scale(1.0 / dt_s);
        let h = solve_at(s)?;
        if h.len() != input.len() * outputs || h.iter().any(|&v| !finite(v)) {
            return Err(TransferError::Invalid("CQ operator shape/nonfinite"));
        }
        for output in 0..outputs {
            let mut y = C64::ZERO;
            for (i, x) in source.iter().enumerate() {
                y = y + h[i * outputs + output] * C64::new(x[bin].re, x[bin].im);
            }
            // Real Laplace operators obey H(conj(s))=conj(H(s)). At the two
            // real contour nodes, roundoff imaginary parts must be negligible.
            if (bin == 0 || bin == n_fft / 2) && y.im.abs() > 1e-10 * (1.0 + y.re.abs()) {
                return Err(TransferError::Invalid(
                    "CQ operator is not a real causal field",
                ));
            }
            target[output][bin] = FftComplex::new(
                y.re,
                if bin == 0 || bin == n_fft / 2 {
                    0.0
                } else {
                    y.im
                },
            );
            if bin > 0 && bin < n_fft / 2 {
                target[output][n_fft - bin] = target[output][bin].conj();
            }
        }
    }
    let mut result = Vec::with_capacity(outputs);
    for mut data in target {
        fft.inverse(&mut data, &mut scratch);
        let mut weight = 1.0;
        let mut row = Vec::with_capacity(steps);
        for v in &data[..steps] {
            row.push(v.re / weight);
            weight *= rho;
        }
        if row.iter().any(|v| !v.is_finite()) {
            return Err(TransferError::Invalid("nonfinite CQ response"));
        }
        result.push(row);
    }
    Ok(result)
}
