//! Reference-assisted interference prediction from the original Welch factors.
//!
//! At each bin fit target = sum(reference * gain) by pivoted, reorthogonalized
//! QR in segment-factor space. Do not form or invert a channel CSD and do not
//! obtain residual power by subtracting two large powers. The learned map can
//! be applied to a separate compatible record, preserving all residual cross
//! spectra and their positive-semidefinite factorization.
//!
//! This is a finite-record linear least-squares map, not source separation,
//! a causal time-domain filter or a claim that reference-correlated content is
//! unwanted. Same-record residuals can overfit; one segment has rank at most one.
//! Stationarity, physical channel identities/units and holdout choice belong to
//! the caller. No cross-frequency covariance or statistical confidence is added.

use super::{
    C64, Dd, WelchChannelTransferBudget, WelchChannelTransferRequirements, WelchConfig,
    WelchCrossSpectrum, WelchError, allocate, poll,
};

/// Bounds for learning a complete frequency-dependent reference transfer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WelchReferenceFitBudget {
    pub maximum_reference_channels: usize,
    pub maximum_target_channels: usize,
    /// Conservative segment-vector visits and small triangular-algebra work.
    /// Includes both reorthogonalization passes and pivot norm scans.
    pub maximum_work_units: usize,
    /// Retained coefficients and peak reusable QR workspace, including headers.
    /// Borrowed input factors and allocator metadata are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Admission computed without allocation or reading the recorded factors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WelchReferenceFitRequirements {
    pub coefficient_values: usize,
    pub work_units: usize,
    pub workspace_bytes: usize,
}

/// Learned target <- reference map under the original exp(-i omega t) law.
///
/// Coefficients and source order remain owned by this immutable model. Numerical
/// rank is declared per bin. Dependent references receive a deterministic basic
/// solution, not a minimum-norm or uniquely identifiable source decomposition.
#[derive(Clone, Debug)]
pub struct WelchReferenceTransfer {
    config: WelchConfig,
    timestep_s: f64,
    input_channels: usize,
    training_segments: usize,
    targets: Vec<usize>,
    references: Vec<usize>,
    ranks: Vec<usize>,
    gains: Vec<Z>,
}

impl WelchCrossSpectrum {
    /// Conservative bounds for M references, T targets and S retained segments.
    /// QR stores O(M*S + M²) scratch; coefficients use O(bins*T*M) storage.
    /// There is no O(input_channels²) covariance matrix or new Fourier transform.
    pub fn reference_fit_requirements(
        &self,
        target_channels: usize,
        reference_channels: usize,
    ) -> Result<WelchReferenceFitRequirements, WelchError> {
        let (t, m, s, b) = (
            target_channels,
            reference_channels,
            self.segments,
            self.bin_count(),
        );
        if t == 0 || m == 0 || t > self.channels || m > self.channels {
            return Err(WelchError::Invalid(
                "nonempty target/reference channel sets required",
            ));
        }
        let mm = size(m.checked_mul(m))?;
        let coefficients = size(b.checked_mul(t).and_then(|n| n.checked_mul(m)))?;
        // Full-rank upper bound. Early rank termination only reduces work.
        // Reference loading/normalization and QR: S*(4*M²+8*M), plus M².
        // Each target: S*(M+2) + M² + 3*M. Index validation is counted too.
        let qr_visits = size(
            mm.checked_mul(4)
                .and_then(|n| m.checked_mul(8).and_then(|k| n.checked_add(k))),
        )?;
        let target_visits = size(t.checked_mul(size(m.checked_add(2))?))?;
        let segment_work = size(
            qr_visits
                .checked_add(target_visits)
                .and_then(|n| n.checked_mul(s)),
        )?;
        let target_algebra = size(t.checked_mul(size(mm.checked_add(size(m.checked_mul(3))?))?))?;
        let per_bin = size(
            segment_work
                .checked_add(mm)
                .and_then(|n| n.checked_add(target_algebra)),
        )?;
        let index_work = size(
            mm.checked_add(size(t.checked_mul(t))?)
                .and_then(|n| t.checked_mul(m).and_then(|k| n.checked_add(k)))
                .and_then(|n| n.checked_add(m))
                .and_then(|n| n.checked_add(t)),
        )?;
        let work_units = size(
            b.checked_mul(per_bin)
                .and_then(|n| n.checked_add(index_work)),
        )?;
        let complex_values = size(
            coefficients
                .checked_add(size(m.checked_mul(s))?)
                .and_then(|n| n.checked_add(mm))
                .and_then(|n| n.checked_add(s))
                .and_then(|n| n.checked_add(m)),
        )?;
        let complex_bytes = size(complex_values.checked_mul(core::mem::size_of::<Z>()))?;
        let indices = size(
            m.checked_mul(2)
                .and_then(|n| n.checked_add(t))
                .and_then(|n| n.checked_add(b)),
        )?;
        let bytes = size(
            complex_bytes
                .checked_add(size(indices.checked_mul(core::mem::size_of::<usize>()))?)
                .and_then(|n| {
                    m.checked_mul(core::mem::size_of::<Dd>() + core::mem::size_of::<f64>())
                        .and_then(|k| n.checked_add(k))
                })
                .and_then(|n| {
                    n.checked_add(
                        core::mem::size_of::<WelchReferenceTransfer>() + core::mem::size_of::<Qr>(),
                    )
                }),
        )?;
        Ok(WelchReferenceFitRequirements {
            coefficient_values: coefficients,
            work_units,
            workspace_bytes: bytes,
        })
    }

    /// Learn complex gains predicting each target from the declared references.
    ///
    /// Channel sets are nonempty, unique and disjoint. Each reference is scaled
    /// to unit segment norm before rank selection, so changes of physical units
    /// do not silently discard a weak but independent reference. A zero channel
    /// contributes no direction. relative_rank_tolerance in (0,1) applies to the
    /// residual norm of these normalized columns; it is a numerical selection
    /// rule, not a noise floor or an uncertainty estimate.
    ///
    /// Pivoted modified Gram-Schmidt repeats every projection before accepting
    /// its residual. The same DD factors supply gains and later cancellation.
    /// All budgets are admitted before allocation. Cancellation is polled at
    /// most every 256 vector/algebra entries and before publishing the model.
    pub fn fit_reference_transfer(
        &self,
        targets: &[usize],
        references: &[usize],
        relative_rank_tolerance: f64,
        budget: WelchReferenceFitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WelchReferenceTransfer, WelchError> {
        poll(&mut cancelled)?;
        if !relative_rank_tolerance.is_finite()
            || relative_rank_tolerance <= 0.0
            || relative_rank_tolerance >= 1.0
        {
            return Err(WelchError::Invalid(
                "reference rank tolerance must be finite and in (0,1)",
            ));
        }
        let needs = self.reference_fit_requirements(targets.len(), references.len())?;
        if references.len() > budget.maximum_reference_channels
            || targets.len() > budget.maximum_target_channels
        {
            return Err(WelchError::BudgetExceeded("reference fit channels"));
        }
        if needs.work_units > budget.maximum_work_units {
            return Err(WelchError::BudgetExceeded("reference fit algebra"));
        }
        if needs.workspace_bytes > budget.maximum_workspace_bytes {
            return Err(WelchError::BudgetExceeded("reference fit workspace"));
        }
        let mut work = Poller {
            cancelled: &mut cancelled,
            since_poll: 0,
        };
        for (i, &channel) in references.iter().enumerate() {
            work.tick()?;
            if channel >= self.channels {
                return Err(WelchError::Invalid("reference channel outside record"));
            }
            for &other in &references[..i] {
                work.tick()?;
                if channel == other {
                    return Err(WelchError::Invalid("duplicate reference channel"));
                }
            }
        }
        for (i, &channel) in targets.iter().enumerate() {
            work.tick()?;
            if channel >= self.channels {
                return Err(WelchError::Invalid("target channel outside record"));
            }
            for &other in targets[..i].iter().chain(references) {
                work.tick()?;
                if channel == other {
                    return Err(WelchError::Invalid(
                        "targets must be unique and separate from references",
                    ));
                }
            }
        }
        let (m, t, s) = (references.len(), targets.len(), self.segments);
        let mut gains = allocate(needs.coefficient_values, Z::ZERO)?;
        let mut ranks = allocate(self.bin_count(), 0usize)?;
        let mut owned_targets = allocate(t, 0usize)?;
        owned_targets.copy_from_slice(targets);
        let mut owned_references = allocate(m, 0usize)?;
        owned_references.copy_from_slice(references);
        let mut qr = Qr::new(m, s)?;
        for (bin, rank) in ranks.iter_mut().enumerate() {
            work.poll()?;
            let start = bin * s * self.channels;
            let input = &self.factors[start..start + s * self.channels];
            qr.load(input, self.channels, references, &mut work)?;
            *rank = qr.factor(relative_rank_tolerance, &mut work)?;
            for (out, &target) in targets.iter().enumerate() {
                let row = &mut gains[(bin * t + out) * m..(bin * t + out + 1) * m];
                qr.fit_target(input, self.channels, target, *rank, row, &mut work)?;
            }
        }
        work.poll()?;
        Ok(WelchReferenceTransfer {
            config: self.config,
            timestep_s: self.timestep_s,
            input_channels: self.channels,
            training_segments: self.segments,
            targets: owned_targets,
            references: owned_references,
            ranks,
            gains,
        })
    }
}

impl WelchReferenceTransfer {
    /// Original target channel indices, in residual output order.
    pub fn target_channels(&self) -> &[usize] {
        &self.targets
    }
    /// Original reference channel indices, in coefficient order.
    pub fn reference_channels(&self) -> &[usize] {
        &self.references
    }
    /// Number of retained training segments, without an independence claim.
    pub const fn training_segment_count(&self) -> usize {
        self.training_segments
    }
    /// Number of physical Welch bins represented by this transfer.
    pub fn bin_count(&self) -> usize {
        self.ranks.len()
    }
    /// Selected numerical reference rank at a bin (at most training segments).
    pub fn rank(&self, bin: usize) -> Result<usize, WelchError> {
        self.ranks
            .get(bin)
            .copied()
            .ok_or(WelchError::Invalid("bin outside reference transfer"))
    }
    /// Predicted target = sum(gain * reference), with NO implicit conjugation.
    /// Indices address the declared target/reference lists, not input channels.
    /// The application path retains the unrounded DD coefficient internally.
    pub fn gain(&self, bin: usize, target: usize, reference: usize) -> Result<C64, WelchError> {
        if bin >= self.bin_count()
            || target >= self.targets.len()
            || reference >= self.references.len()
        {
            return Err(WelchError::Invalid(
                "coefficient outside reference transfer",
            ));
        }
        self.gains[(bin * self.targets.len() + target) * self.references.len() + reference]
            .rounded()
    }

    fn compatible(&self, record: &WelchCrossSpectrum) -> Result<(), WelchError> {
        if record.config != self.config
            || record.timestep_s.to_bits() != self.timestep_s.to_bits()
            || record.channels != self.input_channels
        {
            return Err(WelchError::Invalid(
                "reference transfer needs the same channel layout, Welch configuration and clock",
            ));
        }
        Ok(())
    }

    /// Output storage and products for applying this map to a compatible record.
    /// coefficient_values is zero: this map owns its coefficients; no callback
    /// or frequency matrix is allocated. Training and application may have
    /// different segment counts. Borrowed model/record storage is excluded.
    pub fn residual_requirements(
        &self,
        record: &WelchCrossSpectrum,
    ) -> Result<WelchChannelTransferRequirements, WelchError> {
        self.compatible(record)?;
        let factors = size(
            record
                .bin_count()
                .checked_mul(record.segments)
                .and_then(|n| n.checked_mul(self.targets.len())),
        )?;
        Ok(WelchChannelTransferRequirements {
            coefficient_values: 0,
            complex_products: size(factors.checked_mul(self.references.len()))?,
            output_factors: factors,
            workspace_bytes: size(factors.checked_mul(core::mem::size_of::<C64>()))?,
        })
    }

    /// Form target - predicted_reference on EVERY original segment factor.
    ///
    /// This can use the training record or an independent compatible holdout.
    /// DD gain products and sums are rounded once per residual factor. Thus
    /// residual cross spectra remain PSD by construction, including correlations
    /// among multiple targets, without subtracting almost equal power matrices.
    /// DC and even-length Nyquist stay real; odd-length top bins remain complex.
    ///
    /// No sample clocks, bin widths, original Welch normalization or segment
    /// boundaries change. This does not implement a causal finite-record filter
    /// or reconstruct edge transients. Invalid input, cancellation or arithmetic
    /// overflow returns no partial spectrum and leaves both inputs unchanged.
    pub fn residual_spectrum(
        &self,
        record: &WelchCrossSpectrum,
        budget: WelchChannelTransferBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WelchCrossSpectrum, WelchError> {
        poll(&mut cancelled)?;
        let needs = self.residual_requirements(record)?;
        if self.targets.len() > budget.maximum_output_channels
            || needs.complex_products > budget.maximum_complex_products
            || needs.workspace_bytes > budget.maximum_workspace_bytes
        {
            return Err(WelchError::BudgetExceeded("reference residual transfer"));
        }
        let mut factors = allocate(needs.output_factors, C64::default())?;
        let mut work = Poller {
            cancelled: &mut cancelled,
            since_poll: 0,
        };
        let (t, m) = (self.targets.len(), self.references.len());
        for bin in 0..record.bin_count() {
            work.poll()?;
            for segment in 0..record.segments {
                let start = (bin * record.segments + segment) * record.channels;
                let input = &record.factors[start..start + record.channels];
                for (target, &channel) in self.targets.iter().enumerate() {
                    let mut value = Z::from(input[channel]);
                    let row = &self.gains[(bin * t + target) * m..(bin * t + target + 1) * m];
                    for (&gain, &reference) in row.iter().zip(&self.references) {
                        work.tick()?;
                        value = value.sub(gain.mul(Z::from(input[reference])));
                    }
                    factors[(bin * record.segments + segment) * t + target] = value.rounded()?;
                }
            }
        }
        work.poll()?;
        Ok(WelchCrossSpectrum {
            config: record.config,
            timestep_s: record.timestep_s,
            bin_width_hz: record.bin_width_hz,
            channels: t,
            segments: record.segments,
            input_samples: record.input_samples,
            used_samples: record.used_samples,
            factors,
        })
    }
}

struct Qr {
    m: usize,
    segments: usize,
    columns: Vec<Z>,
    upper: Vec<Z>,
    target: Vec<Z>,
    rhs: Vec<Z>,
    scales: Vec<f64>,
    norms: Vec<Dd>,
    permutation: Vec<usize>,
}

impl Qr {
    fn new(m: usize, segments: usize) -> Result<Self, WelchError> {
        Ok(Self {
            m,
            segments,
            columns: allocate(size(m.checked_mul(segments))?, Z::ZERO)?,
            upper: allocate(size(m.checked_mul(m))?, Z::ZERO)?,
            target: allocate(segments, Z::ZERO)?,
            rhs: allocate(m, Z::ZERO)?,
            scales: allocate(m, 0.0)?,
            norms: allocate(m, Dd::ZERO)?,
            permutation: allocate(m, 0)?,
        })
    }

    fn load(
        &mut self,
        input: &[C64],
        channels: usize,
        references: &[usize],
        work: &mut Poller<'_, impl FnMut() -> bool>,
    ) -> Result<(), WelchError> {
        for value in &mut self.upper {
            work.tick()?;
            *value = Z::ZERO;
        }
        for (j, &channel) in references.iter().enumerate() {
            self.permutation[j] = j;
            let scale = channel_scale(input, channels, channel, work)?;
            self.scales[j] = scale;
            let column = &mut self.columns[j * self.segments..(j + 1) * self.segments];
            for (r, value) in column.iter_mut().enumerate() {
                work.tick()?;
                *value = if scale == 0.0 {
                    Z::ZERO
                } else {
                    Z::from(input[r * channels + channel]).divide(Dd::from_f64(scale))?
                };
            }
            let norm = vector_norm(column, work)?;
            self.norms[j] = norm;
            for value in column {
                work.tick()?;
                if positive(norm) {
                    *value = value.divide(norm)?;
                }
            }
        }
        Ok(())
    }

    fn factor(
        &mut self,
        tolerance: f64,
        work: &mut Poller<'_, impl FnMut() -> bool>,
    ) -> Result<usize, WelchError> {
        let (m, s) = (self.m, self.segments);
        let mut rank = 0;
        for k in 0..m.min(s) {
            let (mut pivot, mut largest) = (k, Dd::ZERO);
            for j in k..m {
                let norm = vector_norm(&self.columns[j * s..(j + 1) * s], work)?;
                if positive(norm - largest) {
                    largest = norm;
                    pivot = j;
                }
            }
            if !positive(largest - Dd::from_f64(tolerance)) {
                break;
            }
            if pivot != k {
                for r in 0..s {
                    work.tick()?;
                    self.columns.swap(k * s + r, pivot * s + r);
                }
                for row in 0..k {
                    work.tick()?;
                    self.upper.swap(row * m + k, row * m + pivot);
                }
                self.permutation.swap(k, pivot);
            }
            self.upper[k * m + k] = Z {
                re: largest,
                im: Dd::ZERO,
            };
            for value in &mut self.columns[k * s..(k + 1) * s] {
                work.tick()?;
                *value = value.divide(largest)?;
            }
            for j in k + 1..m {
                // Reorthogonalize this residual against the same accepted
                // direction; cancellation is retained in DD between passes.
                for _ in 0..2 {
                    let mut projection = Z::ZERO;
                    for r in 0..s {
                        work.tick()?;
                        projection = projection
                            .add(self.columns[k * s + r].conj().mul(self.columns[j * s + r]));
                    }
                    self.upper[k * m + j] = self.upper[k * m + j].add(projection).checked()?;
                    for r in 0..s {
                        work.tick()?;
                        self.columns[j * s + r] = self.columns[j * s + r]
                            .sub(self.columns[k * s + r].mul(projection))
                            .checked()?;
                    }
                }
            }
            rank += 1;
        }
        Ok(rank)
    }

    fn fit_target(
        &mut self,
        input: &[C64],
        channels: usize,
        channel: usize,
        rank: usize,
        gains: &mut [Z],
        work: &mut Poller<'_, impl FnMut() -> bool>,
    ) -> Result<(), WelchError> {
        if rank == 0 {
            return Ok(());
        }
        let scale = channel_scale(input, channels, channel, work)?;
        if scale == 0.0 {
            return Ok(());
        }
        for (r, value) in self.target.iter_mut().enumerate() {
            work.tick()?;
            *value = Z::from(input[r * channels + channel]).divide(Dd::from_f64(scale))?;
        }
        for k in 0..rank {
            let mut projection = Z::ZERO;
            for r in 0..self.segments {
                work.tick()?;
                projection = projection.add(
                    self.columns[k * self.segments + r]
                        .conj()
                        .mul(self.target[r]),
                );
            }
            self.rhs[k] = projection.checked()?;
        }
        for k in (0..rank).rev() {
            let mut value = self.rhs[k];
            for j in k + 1..rank {
                work.tick()?;
                value = value.sub(self.upper[k * self.m + j].mul(self.rhs[j]));
            }
            work.tick()?;
            self.rhs[k] = value.divide(self.upper[k * self.m + k].re)?;
        }
        for k in 0..rank {
            work.tick()?;
            let original = self.permutation[k];
            let normalized = self.rhs[k].divide(self.norms[original])?;
            let gain = Z {
                re: rescale(normalized.re, scale, self.scales[original])?,
                im: rescale(normalized.im, scale, self.scales[original])?,
            };
            gain.rounded()?;
            gains[original] = gain;
        }
        Ok(())
    }
}

fn channel_scale(
    input: &[C64],
    channels: usize,
    channel: usize,
    work: &mut Poller<'_, impl FnMut() -> bool>,
) -> Result<f64, WelchError> {
    let mut scale = 0.0f64;
    for row in input.chunks_exact(channels) {
        work.tick()?;
        scale = scale.max(row[channel].re.abs()).max(row[channel].im.abs());
    }
    Ok(scale)
}

fn vector_norm(
    values: &[Z],
    work: &mut Poller<'_, impl FnMut() -> bool>,
) -> Result<Dd, WelchError> {
    let mut scale = 0.0f64;
    for value in values {
        work.tick()?;
        scale = scale
            .max(value.re.hi.abs())
            .max(value.re.lo.abs())
            .max(value.im.hi.abs())
            .max(value.im.lo.abs());
    }
    if scale == 0.0 {
        return Ok(Dd::ZERO);
    }
    let mut square = Dd::ZERO;
    for &value in values {
        work.tick()?;
        let z = value.divide(Dd::from_f64(scale))?;
        square = square + z.re * z.re + z.im * z.im;
    }
    let norm = square.sqrt() * Dd::from_f64(scale);
    if !finite(norm) || !positive(norm) {
        Err(WelchError::Unrepresentable)
    } else {
        Ok(norm)
    }
}

// Try bounded evaluation orders rather than overflowing an intermediate ratio
// when the final physical gain is still representable. DD extends precision,
// not exponent range; a genuinely unrepresentable nonzero gain refuses.
fn rescale(value: Dd, numerator: f64, denominator: f64) -> Result<Dd, WelchError> {
    if value == Dd::ZERO {
        return Ok(Dd::ZERO);
    }
    let (n, d) = (Dd::from_f64(numerator), Dd::from_f64(denominator));
    for result in [(value * n) / d, (value / d) * n, value * (n / d)] {
        if finite(result) && result.to_f64().is_finite() && result.to_f64() != 0.0 {
            return Ok(result);
        }
    }
    Err(WelchError::Unrepresentable)
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Z {
    re: Dd,
    im: Dd,
}
impl Z {
    const ZERO: Self = Self {
        re: Dd::ZERO,
        im: Dd::ZERO,
    };
    fn from(z: C64) -> Self {
        Self {
            re: Dd::from_f64(z.re),
            im: Dd::from_f64(z.im),
        }
    }
    fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
    fn add(self, b: Self) -> Self {
        Self {
            re: self.re + b.re,
            im: self.im + b.im,
        }
    }
    fn sub(self, b: Self) -> Self {
        Self {
            re: self.re - b.re,
            im: self.im - b.im,
        }
    }
    fn mul(self, b: Self) -> Self {
        Self {
            re: self.re * b.re - self.im * b.im,
            im: self.re * b.im + self.im * b.re,
        }
    }
    fn divide(self, b: Dd) -> Result<Self, WelchError> {
        let result = Self {
            re: self.re / b,
            im: self.im / b,
        }
        .checked()?;
        if (self.re != Dd::ZERO && result.re == Dd::ZERO)
            || (self.im != Dd::ZERO && result.im == Dd::ZERO)
        {
            Err(WelchError::Unrepresentable)
        } else {
            Ok(result)
        }
    }
    fn checked(self) -> Result<Self, WelchError> {
        if finite(self.re) && finite(self.im) {
            Ok(self)
        } else {
            Err(WelchError::Unrepresentable)
        }
    }
    fn rounded(self) -> Result<C64, WelchError> {
        self.checked()?;
        let value = C64::new(self.re.to_f64(), self.im.to_f64());
        if value.re.is_finite() && value.im.is_finite() {
            Ok(value)
        } else {
            Err(WelchError::Unrepresentable)
        }
    }
}
fn finite(v: Dd) -> bool {
    v.hi.is_finite() && v.lo.is_finite()
}
fn positive(v: Dd) -> bool {
    v.hi > 0.0 || (v.hi == 0.0 && v.lo > 0.0)
}
fn size(v: Option<usize>) -> Result<usize, WelchError> {
    v.ok_or(WelchError::BudgetExceeded(
        "reference transfer shape overflow",
    ))
}
struct Poller<'a, F: FnMut() -> bool> {
    cancelled: &'a mut F,
    since_poll: usize,
}
impl<F: FnMut() -> bool> Poller<'_, F> {
    fn tick(&mut self) -> Result<(), WelchError> {
        self.since_poll += 1;
        if self.since_poll == 256 {
            self.poll()?;
        }
        Ok(())
    }
    fn poll(&mut self) -> Result<(), WelchError> {
        self.since_poll = 0;
        poll(self.cancelled)
    }
}

#[cfg(test)]
#[path = "reference/tests.rs"]
mod tests;
