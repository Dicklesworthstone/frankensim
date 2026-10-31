//! Linear propagation of source coherence through the actual BEM field.
//!
//! At one common positive frequency, source peak phasors have second moment
//! C = B B*. Each column of B is propagated coherently; independent columns
//! are combined only after forming field products. Thus direct/scattered and
//! cross-source interference survive, without assembling a dense source C.
//! All physical averages contain the peak-to-mean-square factor 1/2. These
//! are spectral-line quantities, not a silently normalized power density.

use super::*;
use crate::near_field::{FirstOrder, MAX_RECEIVERS};

const MAX_MODES: usize = 4096;
const MAX_PRODUCTS: usize = 64_000_000;

/// Explicit ceiling on complex transfer/outer-product terms per evaluation.
/// Source counts, receiver counts and coherence rank also have fixed caps.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicWorkBudget {
    pub maximum_products: usize,
}

/// A peak phasor in the common exp(-iωt) convention, at a physical receiver.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicField {
    pub pressure_pa: C64,
    pub particle_velocity_m_s: [C64; 3],
}

/// Time/ensemble averages at one frequency. No plane-wave approximation is
/// made: a standing wave can have pressure power but zero active intensity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicFieldStatistics {
    pub mean_square_pressure_pa2: f64,
    /// 1/2 E[u_i conj(u_j)], including complex off-diagonal correlations.
    pub velocity_cross_m2_s2: [[C64; 3]; 3],
    /// 1/2 E[p conj(u)]; real = active, imaginary = reactive intensity under
    /// exp(-iωt). Components retain their signs; no positive clipping occurs.
    pub complex_intensity_w_m2: [C64; 3],
}

/// Immutable source-to-receiver p/u transfer bank of one existing BEM owner.
/// Preparation performs all boundary solves and receiver quadrature. Applying
/// sources/factors is allocation-free and never reassembles the BEM operator.
/// Spatial/model uncertainty is still unresolved by this algebraic observer.
pub struct PreparedHarmonicObservation {
    omega: f64,
    medium: Medium,
    sources: Vec<SourceId>,
    receivers: Vec<ReceiverSpec>,
    producer_hash: ContentHash,
    // Source-major, receiver-minor; p, ux, uy, uz per declared source unit.
    transfers: Vec<[C64; 4]>,
}

pub(super) fn admit_frequency(omega: f64, medium: Medium) -> Result<(), TransferError> {
    let momentum = omega * medium.density;
    let impedance = medium.density * medium.sound_speed;
    if !omega.is_finite()
        || omega <= 0.0
        || !medium.density.is_finite()
        || medium.density <= 0.0
        || !medium.sound_speed.is_finite()
        || medium.sound_speed <= 0.0
        || !momentum.is_finite()
        || momentum <= 0.0
        || !(1.0 / momentum).is_finite()
        || !impedance.is_finite()
        || impedance <= 0.0
    {
        return Err(TransferError::Invalid(
            "invalid harmonic frequency or medium",
        ));
    }
    Ok(())
}

// Neumaier summation retains weak residuals after cancellation, separately
// for each physical component. Do not reduce absolute amplitudes instead.
#[derive(Clone, Copy, Debug, Default)]
struct Sum {
    value: C64,
    correction: C64,
}
impl Sum {
    fn add(&mut self, z: C64) {
        fn component(sum: &mut f64, correction: &mut f64, x: f64) {
            let next = *sum + x;
            *correction += if sum.abs() >= x.abs() {
                (*sum - next) + x
            } else {
                (x - next) + *sum
            };
            *sum = next;
        }
        component(&mut self.value.re, &mut self.correction.re, z.re);
        component(&mut self.value.im, &mut self.correction.im, z.im);
    }
    fn finish(self) -> C64 {
        self.value + self.correction
    }
}

#[derive(Default)]
struct StatisticsSum {
    pressure: Sum,
    velocity: [[Sum; 3]; 3],
    intensity: [Sum; 3],
}
impl StatisticsSum {
    // The field has already been scaled by sqrt(peak second moment / 2).
    fn add(&mut self, field: [C64; 4]) {
        self.pressure.add(C64::from_re(field[0].norm_sq()));
        for i in 0..3 {
            self.intensity[i].add(field[0] * field[i + 1].conj());
            for j in i..3 {
                let product = if i == j {
                    C64::from_re(field[i + 1].norm_sq())
                } else {
                    field[i + 1] * field[j + 1].conj()
                };
                self.velocity[i][j].add(product);
            }
        }
    }
    fn finish(self) -> Result<HarmonicFieldStatistics, TransferError> {
        let mut out = HarmonicFieldStatistics {
            mean_square_pressure_pa2: self.pressure.finish().re,
            complex_intensity_w_m2: self.intensity.map(Sum::finish),
            ..HarmonicFieldStatistics::default()
        };
        for i in 0..3 {
            for j in i..3 {
                let value = self.velocity[i][j].finish();
                out.velocity_cross_m2_s2[i][j] = value;
                out.velocity_cross_m2_s2[j][i] = value.conj();
            }
        }
        if !out.mean_square_pressure_pa2.is_finite()
            || out
                .velocity_cross_m2_s2
                .iter()
                .flatten()
                .any(|&z| !finite(z))
            || out.complex_intensity_w_m2.iter().any(|&z| !finite(z))
        {
            return Err(TransferError::Invalid("harmonic field statistics overflow"));
        }
        Ok(out)
    }
}

impl PreparedHarmonicObservation {
    pub(super) fn from_samples(
        omega: f64,
        medium: Medium,
        sources: &[SourceId],
        receivers: &[ReceiverSpec],
        producer_hash: ContentHash,
        samples: Vec<C64>,
    ) -> Result<Self, TransferError> {
        admit_frequency(omega, medium)?;
        cap(sources.len(), MAX_OBSERVATION_INPUTS)?;
        cap(receivers.len(), MAX_RECEIVERS)?;
        if sources.is_empty()
            || receivers.is_empty()
            || samples.len() != 4 * sources.len() * receivers.len()
            || samples.iter().any(|&z| !finite(z))
        {
            return Err(TransferError::Invalid(
                "incomplete harmonic field source basis",
            ));
        }
        // Euler under exp(-iωt): u = grad(p)/(iωρ), not its conjugate.
        let euler = C64::new(0.0, -1.0 / (omega * medium.density));
        let transfers: Vec<_> = samples
            .chunks_exact(4)
            .map(|v| [v[0], euler * v[1], euler * v[2], euler * v[3]])
            .collect();
        if transfers.iter().flatten().any(|&z| !finite(z)) {
            return Err(TransferError::Invalid(
                "harmonic velocity conversion overflow",
            ));
        }
        Ok(Self {
            omega,
            medium,
            sources: sources.to_vec(),
            receivers: receivers.to_vec(),
            producer_hash,
            transfers,
        })
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Ordered physical input coordinates. Each coefficient/factor entry uses
    /// the SI unit of this source coordinate; mixed source kinds stay distinct.
    pub fn sources(&self) -> &[SourceId] {
        &self.sources
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        &self.receivers
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        ProducerUncertainty::Unresolved {
            source_hash: self.producer_hash,
        }
    }
    fn shape(&self, output: usize) -> Result<(), TransferError> {
        if output != self.receivers.len() {
            return Err(TransferError::Invalid("harmonic receiver output shape"));
        }
        Ok(())
    }
    fn work(
        &self,
        modes: usize,
        components: usize,
        cross: bool,
        budget: HarmonicWorkBudget,
    ) -> Result<(), TransferError> {
        cap(modes, MAX_MODES)?;
        cap(budget.maximum_products, MAX_PRODUCTS)?;
        if budget.maximum_products == 0 {
            return Err(TransferError::Invalid("empty harmonic work allowance"));
        }
        let r = self.receivers.len();
        let per_mode =
            components * self.sources.len() * r + if cross { r * (r + 1) / 2 } else { 10 * r };
        let requested = per_mode
            .checked_mul(modes)
            .ok_or(TransferError::Invalid("harmonic work overflow"))?;
        cap(requested, budget.maximum_products)
    }
    fn factors(&self, b: &[C64], modes: usize) -> Result<(), TransferError> {
        cap(modes, MAX_MODES)?;
        if b.len() != modes * self.sources.len() || b.iter().any(|&z| !finite(z)) {
            return Err(TransferError::Invalid(
                "harmonic coherence-factor shape or value",
            ));
        }
        Ok(())
    }
    fn powers(&self, powers: &[f64]) -> Result<(), TransferError> {
        if powers.len() != self.sources.len() || powers.iter().any(|&x| !x.is_finite() || x < 0.0) {
            return Err(TransferError::Invalid(
                "independent source peak second moments",
            ));
        }
        Ok(())
    }
    fn field(&self, coefficients: &[C64], receiver: usize) -> Result<[C64; 4], TransferError> {
        let mut sums = [Sum::default(); 4];
        for (source, &coefficient) in coefficients.iter().enumerate() {
            let transfer = self.transfers[source * self.receivers.len() + receiver];
            for c in 0..4 {
                sums[c].add(transfer[c] * coefficient);
            }
        }
        let field = sums.map(Sum::finish);
        if field.iter().any(|&z| !finite(z)) {
            return Err(TransferError::Invalid("coherent harmonic field overflow"));
        }
        Ok(field)
    }

    /// O(sources × receivers), allocation-free. Sources share this observer's
    /// frequency and phase reference. On failure, discard the output buffer.
    pub fn evaluate_coherent(
        &self,
        coefficients: &[C64],
        output: &mut [HarmonicField],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        self.shape(output.len())?;
        self.factors(coefficients, 1)?;
        for (receiver, out) in output.iter_mut().enumerate() {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let f = self.field(coefficients, receiver)?;
            *out = HarmonicField {
                pressure_pa: f[0],
                particle_velocity_m_s: [f[1], f[2], f[3]],
            };
        }
        Ok(())
    }

    /// B is mode-major, source-minor: B[mode * sources + source]. Source peak
    /// second moment is B B*. A deterministic coherent drive is one mode; an
    /// uncorrelated source bank can use `evaluate_independent_sources` directly.
    /// Zero modes explicitly mean silence. O(sources × receivers × modes), no
    /// dense source covariance or allocation. Discard output after any error.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[C64],
        modes: usize,
        budget: HarmonicWorkBudget,
        output: &mut [HarmonicFieldStatistics],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        self.shape(output.len())?;
        self.work(modes, 4, false, budget)?;
        self.factors(factors, modes)?;
        for (receiver, out) in output.iter_mut().enumerate() {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let mut sum = StatisticsSum::default();
            for factor in factors.chunks_exact(self.sources.len()) {
                if cancelled() {
                    return Err(TransferError::Cancelled);
                }
                sum.add(
                    self.field(factor, receiver)?
                        .map(|v| v.scale(core::f64::consts::FRAC_1_SQRT_2)),
                );
            }
            *out = sum.finish()?;
        }
        Ok(())
    }

    /// Diagonal source covariance in O(sources × receivers), without forming
    /// a diagonal factor matrix. Inputs are E[|source peak phasor|²], in each
    /// source coordinate's squared SI unit, not RMS amplitudes or dB values.
    pub fn evaluate_independent_sources(
        &self,
        peak_second_moments: &[f64],
        budget: HarmonicWorkBudget,
        output: &mut [HarmonicFieldStatistics],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        self.shape(output.len())?;
        self.powers(peak_second_moments)?;
        self.work(1, 14, false, budget)?;
        for (receiver, out) in output.iter_mut().enumerate() {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let mut sum = StatisticsSum::default();
            for (source, &power) in peak_second_moments.iter().enumerate() {
                let scale = det::sqrt(power) * core::f64::consts::FRAC_1_SQRT_2;
                sum.add(
                    self.transfers[source * self.receivers.len() + receiver]
                        .map(|v| v.scale(scale)),
                );
            }
            *out = sum.finish()?;
        }
        Ok(())
    }

    /// Prepare one existing ideal first-order microphone per receiver. This
    /// projects the full near-field p/u before taking cross spectra, preserving
    /// pressure/velocity cross terms and microphone-array phase differences.
    pub fn microphones(
        &self,
        patterns: &[FirstOrder],
    ) -> Result<PreparedMicrophoneArray<'_>, TransferError> {
        self.shape(patterns.len())?;
        let mut transfers = Vec::with_capacity(self.transfers.len());
        for (pair, &f) in self.transfers.iter().enumerate() {
            transfers.push(patterns[pair % patterns.len()].observe(
                f[0],
                [f[1], f[2], f[3]],
                self.medium,
            )?);
        }
        let entries = patterns.len() * patterns.len();
        Ok(PreparedMicrophoneArray {
            field: self,
            transfers,
            sums: vec![Sum::default(); entries],
            spectra: vec![C64::ZERO; entries],
        })
    }
}

/// Prepared microphone transfers and reusable cross-spectrum workspace.
/// Returned matrices are receiver-major C_ij = 1/2 E[m_i conj(m_j)] [Pa²].
/// A beamformer y = w* m therefore has power w* C w. Spatial source/BEM
/// uncertainty remains that of `field`; this is not a coherence certificate.
pub struct PreparedMicrophoneArray<'a> {
    field: &'a PreparedHarmonicObservation,
    transfers: Vec<C64>,
    sums: Vec<Sum>,
    spectra: Vec<C64>,
}
impl PreparedMicrophoneArray<'_> {
    fn add_mode(&mut self, amplitudes: &[C64]) {
        let r = amplitudes.len();
        for i in 0..r {
            for j in i..r {
                let product = if i == j {
                    C64::from_re(amplitudes[i].norm_sq())
                } else {
                    amplitudes[i] * amplitudes[j].conj()
                };
                self.sums[i * r + j].add(product);
            }
        }
    }
    fn finish(&mut self) -> Result<&[C64], TransferError> {
        let r = self.field.receivers.len();
        for i in 0..r {
            for j in i..r {
                let value = self.sums[i * r + j].finish();
                if !finite(value) {
                    return Err(TransferError::Invalid("microphone cross spectrum overflow"));
                }
                self.spectra[i * r + j] = value;
                self.spectra[j * r + i] = value.conj();
            }
        }
        Ok(&self.spectra)
    }
    /// Factorized source coherence, with no runtime allocation. Work is
    /// O(modes × (sources × receivers + receivers²)); each source/receiver
    /// amplitude is computed once per mode, never once per receiver pair.
    pub fn evaluate_coherence_modes(
        &mut self,
        factors: &[C64],
        modes: usize,
        budget: HarmonicWorkBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[C64], TransferError> {
        self.field.work(modes, 1, true, budget)?;
        self.field.factors(factors, modes)?;
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        self.sums.fill(Sum::default());
        let r = self.field.receivers.len();
        for factor in factors.chunks_exact(self.field.sources.len()) {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let mut amplitudes = [C64::ZERO; MAX_RECEIVERS];
            for (receiver, amplitude) in amplitudes[..r].iter_mut().enumerate() {
                let mut sum = Sum::default();
                for (source, &coefficient) in factor.iter().enumerate() {
                    sum.add(self.transfers[source * r + receiver] * coefficient);
                }
                *amplitude = sum.finish().scale(core::f64::consts::FRAC_1_SQRT_2);
            }
            self.add_mode(&amplitudes[..r]);
        }
        self.finish()
    }
    /// Independent sources can still produce correlated microphones. Retains
    /// that full spatial cross spectrum in O(sources × receivers²) work, without
    /// constructing or multiplying a dense source covariance.
    pub fn evaluate_independent_sources(
        &mut self,
        peak_second_moments: &[f64],
        budget: HarmonicWorkBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[C64], TransferError> {
        self.field.powers(peak_second_moments)?;
        // Each source supplies one projected vector and one upper outer product.
        let r = self.field.receivers.len();
        cap(budget.maximum_products, MAX_PRODUCTS)?;
        if budget.maximum_products == 0 {
            return Err(TransferError::Invalid("empty harmonic work allowance"));
        }
        cap(
            self.field.sources.len() * (r + r * (r + 1) / 2),
            budget.maximum_products,
        )?;
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        self.sums.fill(Sum::default());
        for (source, &power) in peak_second_moments.iter().enumerate() {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let scale = det::sqrt(power) * core::f64::consts::FRAC_1_SQRT_2;
            let mut amplitudes = [C64::ZERO; MAX_RECEIVERS];
            for (receiver, amplitude) in amplitudes[..r].iter_mut().enumerate() {
                *amplitude = self.transfers[source * r + receiver].scale(scale);
            }
            self.add_mode(&amplitudes[..r]);
        }
        self.finish()
    }
}

#[cfg(test)]
mod tests;
