//! Original-initial-value moments of a harmonically driven relaxing state.

use super::*;

pub const MAX_HARMONIC_DECAY_TONES: usize = 64;
const DEGREE: usize = 32;
type Pair = [Dd; 2];

/// Original angular source clock, with a separately retained interval width.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPhaseInterval {
    pub start: Pair,
    pub end: Pair,
    /// end-start, evaluated directly even when that subtraction would vanish.
    pub increment: Pair,
}

/// exp(-i omega t) at both endpoints, without forming an f64 start+duration.
pub fn phase_angular_interval_from_clock_dd(
    omega: f64,
    start: Dd,
    duration: f64,
) -> Result<HarmonicPhaseInterval, WindowError> {
    if !omega.is_finite() {
        return Err(WindowError::InvalidFrequency);
    }
    if !start.hi.is_finite() || !start.lo.is_finite()
        || !duration.is_finite() || duration < 0.0
    {
        return Err(WindowError::InvalidWindow);
    }
    let start = checked(Dd::from_pair(start.hi, start.lo))?;
    let angle = checked_source_product(Dd::from_f64(omega), start)?;
    let beta = checked_source_product(Dd::from_f64(omega), Dd::from_f64(duration))?;
    let (a, b) = proper_phase(angle);
    let (c, d) = proper_phase(beta);
    let (_, half_sine) = proper_phase(beta * Dd::from_f64(0.5));
    let change = [-Dd::from_f64(2.0) * checked_source_product(half_sine, half_sine)?, d];
    Ok(HarmonicPhaseInterval {
        start: [a, b],
        end: multiply([a, b], [c, d])?,
        increment: multiply([a, b], change)?,
    })
}

/// Explicit admission before allocation or temporal kernel preparation.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicDecayIntervalLimits {
    pub maximum_tones: usize,
    /// Retained struct and vector capacity; excludes caller inputs and stack.
    pub maximum_workspace_bytes: usize,
    /// Conservative scalar-kernel/series allowance for preparation plus one query.
    /// Transcendental primitives count as one; this is not a wall-clock limit.
    pub maximum_work: usize,
}

/// Means and original endpoint of q'+q/tau=Re(sum drive_n exp(-i omega_n t)).
#[derive(Clone, Copy, Debug)]
pub struct HarmonicDecayMoments {
    pub end: Dd,
    /// Direct endpoint change, retaining changes below endpoint resolution.
    pub increment: Dd,
    pub mean_square: Dd,
    /// Mean q(t)*Re(sum work_n exp(-i omega_n t)), independently integrated.
    pub mean_mixed: Dd,
    duration: Dd,
}

impl HarmonicDecayMoments {
    pub fn square_integral(self) -> Result<Dd, WindowError> {
        checked_source_product(self.mean_square, self.duration)
    }

    pub fn mixed_integral(self) -> Result<Dd, WindowError> {
        checked_source_product(self.mean_mixed, self.duration)
    }
}

#[derive(Clone, Copy, Debug)]
struct Tone {
    omega: f64,
    beta: Dd,
    phase: HarmonicPhaseInterval,
    causal: Pair,
    front: Pair,
    particular: Pair,
}

/// Reusable, allocation-free-query moments for one original relaxation time.
///
/// Frequencies are signed angular frequencies, amplitudes are real-field peak
/// phasors, and drive is the ORIGINAL forcing rate, not a reconstructed settled
/// state. A Maxwell spring therefore supplies drive=-i*omega*strain and its
/// retained original initial spring strain. No artificial equilibration occurs.
///
/// Short intervals use the ODE Taylor recurrence, then an orthogonal polynomial
/// square of the COMPLETE state. Its admitted analytic remainder is below
/// 1e-28 of that state's L2 norm. Longer intervals use DD closed pair/decay
/// moments and refuse cancellation with less than 1e-24 of the absolute term
/// sum remaining. Neither path clamps negative energy, samples a time grid, or
/// provides an interval/null-relative roundoff certificate. Nonfinite results
/// and unresolved positive states refuse. Callers may scale their physical
/// projections before querying to stay within the DD exponent range.
///
/// A finite positive tau and representable h/tau and omega*h are required.
/// Zero duration returns unchanged state and zero means/integrals by convention.
/// Checked integral methods are optional: a representable normalized mean need
/// not be lost merely because multiplication by h would underflow. Polynomial
/// terms below exponent range may vanish; an active state's entire square may
/// not. Mixed moments are signed and have no relative guarantee near a null.
#[derive(Debug)]
pub struct PreparedHarmonicDecayInterval {
    duration: Dd,
    alpha: Dd,
    decay: Dd,
    decrement: Dd,
    short: bool,
    tones: Vec<Tone>,
    products: Vec<HarmonicProduct>,
}

impl PreparedHarmonicDecayInterval {
    pub fn required_workspace_bytes(n: usize) -> Result<usize, ProductSpectrumError> {
        if n > MAX_HARMONIC_DECAY_TONES { return Err(ProductSpectrumError::Budget); }
        Ok(core::mem::size_of::<Self>() + n * core::mem::size_of::<Tone>()
            + n * (n + 1) / 2 * core::mem::size_of::<HarmonicProduct>())
    }

    pub fn required_work(n: usize) -> Result<usize, ProductSpectrumError> {
        if n > MAX_HARMONIC_DECAY_TONES { return Err(ProductSpectrumError::Budget); }
        Ok(128 * (n + 1) * (n + 1) + 512 * (n + 1) + 8192)
    }

    pub fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>() + self.tones.capacity() * core::mem::size_of::<Tone>()
            + self.products.capacity() * core::mem::size_of::<HarmonicProduct>()
    }

    pub fn phases(&self, index: usize) -> Option<HarmonicPhaseInterval> {
        self.tones.get(index).map(|tone| tone.phase)
    }

    pub fn new(
        angular_frequencies: &[f64], start: Dd, duration: f64, relaxation_time: f64,
        limits: HarmonicDecayIntervalLimits,
    ) -> Result<Self, ProductSpectrumError> {
        Self::new_cancellable(angular_frequencies, start, duration, relaxation_time, limits, &mut || false)
    }

    pub fn new_cancellable(
        angular_frequencies: &[f64], start: Dd, duration: f64, relaxation_time: f64,
        limits: HarmonicDecayIntervalLimits, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ProductSpectrumError> {
        spectrum_poll(cancelled)?;
        let n = angular_frequencies.len();
        if n > limits.maximum_tones
            || Self::required_workspace_bytes(n)? > limits.maximum_workspace_bytes
            || Self::required_work(n)? > limits.maximum_work
        { return Err(ProductSpectrumError::Budget); }
        if !relaxation_time.is_finite() || relaxation_time <= 0.0
            || !duration.is_finite() || duration < 0.0
            || !start.hi.is_finite() || !start.lo.is_finite()
        { return Err(WindowError::InvalidWindow.into()); }
        let width = Dd::from_f64(duration);
        let alpha = checked(width / Dd::from_f64(relaxation_time))?;
        if duration > 0.0 && is_zero(alpha) { return Err(WindowError::Unrepresentable.into()); }
        let (decay, decrement) = decay_parts(alpha)?;
        let mut tones = Vec::new();
        tones.try_reserve_exact(n).map_err(|_| ProductSpectrumError::Allocation)?;
        let mut short = alpha.hi <= 0.5;
        for &omega in angular_frequencies {
            spectrum_poll(cancelled)?;
            let phase = phase_angular_interval_from_clock_dd(omega, start, duration)?;
            let beta = checked_source_product(Dd::from_f64(omega), width)?;
            short &= beta.hi.abs() <= 0.5;
            tones.push(Tone { omega, beta, phase, causal: [Dd::ZERO; 2],
                front: [Dd::ZERO; 2], particular: [Dd::ZERO; 2] });
        }
        let mut products = Vec::new();
        if duration > 0.0 {
            for tone in &mut tones {
                spectrum_poll(cancelled)?;
                tone.causal = integration_kernels(alpha, tone.beta, decrement)?.1;
                if !short {
                    tone.particular = divide([width, Dd::ZERO], alpha, -tone.beta)?;
                    let (_, sine) = proper_phase(tone.beta);
                    let (_, half_sine) = proper_phase(tone.beta * Dd::from_f64(0.5));
                    let numerator = [checked(decrement + decay * Dd::from_f64(2.0)
                        * half_sine * half_sine)?, checked(-decay * sine)?];
                    tone.front = divide(numerator, alpha, tone.beta)?;
                }
            }
            if !short {
                products.try_reserve_exact(n * (n + 1) / 2)
                    .map_err(|_| ProductSpectrumError::Allocation)?;
                let window = HarmonicWindow::new(0.0, duration)?;
                for i in 0..n {
                    spectrum_poll(cancelled)?;
                    for j in 0..=i {
                        products.push(window.product_kernel(tones[i].omega, tones[j].omega)?);
                    }
                }
            }
        }
        spectrum_poll(cancelled)?;
        Ok(Self { duration: width, alpha, decay, decrement, short, tones, products })
    }

    pub fn evaluate(&self, initial: Dd, drive: &[Pair], work: &[Pair])
        -> Result<HarmonicDecayMoments, ProductSpectrumError>
    {
        self.evaluate_cancellable(initial, drive, work, &mut || false)
    }

    pub fn evaluate_cancellable(
        &self, initial: Dd, drive: &[Pair], work: &[Pair],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicDecayMoments, ProductSpectrumError> {
        spectrum_poll(cancelled)?;
        check_amplitude(initial)?;
        let initial = checked(Dd::from_pair(initial.hi, initial.lo))?;
        let drives = self.local_fields(drive)?;
        let works = self.local_fields(work)?;
        if is_zero(self.duration) {
            return Ok(HarmonicDecayMoments { end: initial, increment: Dd::ZERO,
                mean_square: Dd::ZERO, mean_mixed: Dd::ZERO, duration: self.duration });
        }
        let mut convolution = Dd::ZERO;
        for (i, tone) in self.tones.iter().enumerate() {
            spectrum_poll(cancelled)?;
            convolution = checked(convolution + integrate_source_tone(self.duration, drives[i], tone.causal)?)?;
        }
        let end = checked(decayed(initial, self.alpha, self.decay)? + convolution)?;
        let increment = checked(-checked_source_product(initial, self.decrement)? + convolution)?;
        let (mean_square, mean_mixed) = if self.short {
            self.polynomial(initial, &drives, &works, cancelled)?
        } else { self.closed(initial, &drives, &works, cancelled)? };
        let active = !is_zero(initial) || self.tones.iter().enumerate().any(|(i, tone)|
            !is_zero(drives[i][0]) || (tone.omega != 0.0 && !is_zero(drives[i][1])));
        if mean_square.lt(Dd::ZERO) || (active && is_zero(mean_square)) {
            return Err(WindowError::Unrepresentable.into());
        }
        spectrum_poll(cancelled)?;
        Ok(HarmonicDecayMoments { end, increment, mean_square, mean_mixed, duration: self.duration })
    }

    // Coalesce equal real-waveform clocks BEFORE rotation or products.
    // Opposite frequencies describe the same real carrier with conjugated phasors.
    fn local_fields(&self, fields: &[Pair]) -> Result<[Pair; MAX_HARMONIC_DECAY_TONES], WindowError> {
        if fields.len() != self.tones.len() { return Err(WindowError::InvalidBandShape); }
        let mut local = [[Dd::ZERO; 2]; MAX_HARMONIC_DECAY_TONES];
        for (i, values) in fields.iter().enumerate() {
            let first = (0..i).find(|&j| self.tones[j].omega.abs() == self.tones[i].omega.abs()).unwrap_or(i);
            for k in 0..2 {
                check_amplitude(values[k])?;
                let value = Dd::from_pair(values[k].hi, values[k].lo);
                let value = if k == 1 && self.tones[first].omega != self.tones[i].omega { -value } else { value };
                local[first][k] = checked(local[first][k] + value)?;
            }
        }
        for (tone, values) in self.tones.iter().zip(&mut local) {
            *values = multiply(*values, tone.phase.start)?;
        }
        Ok(local)
    }

    fn polynomial(&self, initial: Dd, drive: &[Pair], work: &[Pair],
        cancelled: &mut impl FnMut() -> bool) -> Result<(Dd, Dd), ProductSpectrumError>
    {
        let mut forcing = [Dd::ZERO; DEGREE + 1];
        let mut polynomial = [Dd::ZERO; DEGREE + 1];
        let mut other = [Dd::ZERO; DEGREE + 1];
        let mut envelope = Dd::ZERO;
        let mut largest_beta = Dd::ZERO;
        for (i, tone) in self.tones.iter().enumerate() {
            spectrum_poll(cancelled)?;
            let mut d = [checked_source_product(drive[i][0], self.duration)?,
                checked_source_product(drive[i][1], self.duration)?];
            let mut w = work[i];
            envelope = checked(envelope + d[0].abs() + d[1].abs())?;
            if largest_beta.lt(tone.beta.abs()) { largest_beta = tone.beta.abs(); }
            for k in 0..=DEGREE {
                forcing[k] = checked(forcing[k] + d[0])?;
                other[k] = checked(other[k] + w[0])?;
                let r = tone.beta / Dd::from_f64((k + 1) as f64);
                d = [checked(d[1] * r)?, checked(-d[0] * r)?];
                w = [checked(w[1] * r)?, checked(-w[0] * r)?];
            }
        }
        polynomial[0] = initial;
        for k in 0..DEGREE {
            polynomial[k + 1] = checked((forcing[k] - self.alpha * polynomial[k])
                / Dd::from_f64((k + 1) as f64))?;
        }
        let a = legendre(&polynomial)?;
        let b = legendre(&other)?;
        let mut square = Dd::ZERO;
        let mut mixed = Dd::ZERO;
        for k in 0..=DEGREE {
            let divisor = Dd::from_f64((2 * k + 1) as f64);
            square = checked(square + checked(a[k] * a[k])? / divisor)?;
            mixed = checked(mixed + checked(a[k] * b[k])? / divisor)?;
        }
        // Coefficient majorant for the original ODE, not a settled/transient
        // subtraction: |c_k| <= |q0| alpha^k/k! + S(alpha+B)^(k-1)/k!.
        let rate = checked(self.alpha + largest_beta)?;
        let mut initial_tail = initial.abs();
        let mut driven_tail = envelope;
        for k in 1..=DEGREE + 1 {
            initial_tail = checked(initial_tail * self.alpha / Dd::from_f64(k as f64))?;
            if k > 1 { driven_tail = checked(driven_tail * rate / Dd::from_f64(k as f64))?; }
        }
        let tail = checked((initial_tail + driven_tail)
            / (Dd::ONE - rate / Dd::from_f64((DEGREE + 2) as f64)))?;
        if checked(square.sqrt() * Dd::from_f64(1e-28))?.lt(tail) {
            return Err(WindowError::Unrepresentable.into());
        }
        Ok((square, mixed))
    }

    fn closed(&self, initial: Dd, drive: &[Pair], work: &[Pair],
        cancelled: &mut impl FnMut() -> bool) -> Result<(Dd, Dd), ProductSpectrumError>
    {
        let mut amplitudes = [[Dd::ZERO; 2]; MAX_HARMONIC_DECAY_TONES];
        let mut transient = initial;
        for (i, tone) in self.tones.iter().enumerate() {
            amplitudes[i] = multiply(drive[i], tone.particular)?;
            transient = checked(transient - amplitudes[i][0])?;
        }
        let gain = if self.alpha.hi > 373.0 {
            checked(Dd::from_f64(0.5) / self.alpha)?
        } else {
            let twice = self.alpha * Dd::from_f64(2.0);
            checked(decay_parts(twice)?.1 / twice)?
        };
        let mut square = checked_source_product(checked_source_product(transient, gain)?, transient)?;
        let mut envelope = square.abs();
        let mut mixed = Dd::ZERO;
        for (i, tone) in self.tones.iter().enumerate() {
            spectrum_poll(cancelled)?;
            let cross = checked(Dd::from_f64(2.0) * transient * multiply(amplitudes[i], tone.front)?[0])?;
            square = checked(square + cross)?;
            envelope = checked(envelope + cross.abs())?;
            mixed = checked(mixed + transient * multiply(work[i], tone.front)?[0])?;
            for j in 0..=i {
                let kernel = self.products[i * (i + 1) / 2 + j];
                let factor = Dd::from_f64(if i == j { 1.0 } else { 2.0 });
                let term = checked(factor * kernel.mean_real_product_parts_dd(amplitudes[i], amplitudes[j])?)?;
                square = checked(square + term)?;
                envelope = checked(envelope + term.abs())?;
                mixed = checked(mixed + kernel.mean_real_product_parts_dd(amplitudes[i], work[j])?)?;
                if i != j {
                    mixed = checked(mixed + kernel.mean_real_product_parts_dd(work[i], amplitudes[j])?)?;
                }
            }
        }
        if square.lt(checked(envelope * Dd::from_f64(1e-24))?) {
            return Err(WindowError::Unrepresentable.into());
        }
        Ok((square, mixed))
    }
}

fn multiply(a: Pair, b: Pair) -> Result<Pair, WindowError> {
    let (real, imaginary) = rotate(a, b)?;
    Ok([real, imaginary])
}

// Smith-scaled complex quotient; no alpha²+beta² or reciprocal tau.
fn divide(a: Pair, real: Dd, imaginary: Dd) -> Result<Pair, WindowError> {
    let (r, ar, ai, swap) = if !real.abs().lt(imaginary.abs()) {
        (checked(imaginary / real)?, checked(a[0] / real)?, checked(a[1] / real)?, false)
    } else {
        (checked(real / imaginary)?, checked(a[0] / imaginary)?, checked(a[1] / imaginary)?, true)
    };
    let denominator = checked(Dd::ONE + r * r)?;
    if swap {
        Ok([checked((ar * r + ai) / denominator)?, checked((ai * r - ar) / denominator)?])
    } else {
        Ok([checked((ar + ai * r) / denominator)?, checked((ai - ar * r) / denominator)?])
    }
}

fn decayed(value: Dd, alpha: Dd, decay: Dd) -> Result<Dd, WindowError> {
    if alpha.hi <= 350.0 || is_zero(value) { return checked(value * decay); }
    if alpha.hi > 1492.0 { return Ok(Dd::ZERO); }
    // A subnormal exp(-alpha) can lose a representable large initial state.
    // Each third-decay stays normal throughout this finite output range.
    let third = decay_parts(alpha / Dd::from_f64(3.0))?.0;
    checked(checked(checked(value * third)? * third)? * third)
}

// Power -> shifted Legendre: integral u^k P_l(2u-1) du
// = (k!)²/((k-l)!(k+l+1)!). Weights are positive; no monomial Gram subtraction.
fn legendre(power: &[Dd; DEGREE + 1]) -> Result<[Dd; DEGREE + 1], WindowError> {
    let mut modes = [Dd::ZERO; DEGREE + 1];
    for (k, &coefficient) in power.iter().enumerate() {
        let mut weight = Dd::ONE / Dd::from_f64((k + 1) as f64);
        for (l, mode) in modes.iter_mut().enumerate().take(k + 1) {
            *mode = checked(*mode + weight * coefficient)?;
            weight = checked(weight * Dd::from_f64(((2 * l + 3) * (k - l)) as f64)
                / Dd::from_f64(((2 * l + 1) * (k + l + 2)) as f64))?;
        }
    }
    Ok(modes)
}

#[cfg(test)]
mod tests;
