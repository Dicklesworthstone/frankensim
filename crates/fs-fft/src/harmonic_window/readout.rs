//! Original-carrier linear readouts of a narrow, rapidly oscillating band.
//!
//! Write omega_i=omega_c+delta_i and t=t_mid+h*x, h=T/2. Expand only the
//! envelope: B_n=sum_i A_i*exp(-i*delta_i*t_mid)*(delta_i*h)^n/n!.
//! The readout is exp(-i*omega_c*t_mid)*sum_n (-i)^n*B_n*M_n(C), where
//! C=omega_c*h and M_n is the normalized profile integral of x^n*exp(-i*C*x).
//! Coherent fields therefore cancel before their small surviving envelope
//! modes are multiplied by fast-carrier integrals. Independently rounded
//! sinc values at each original frequency are never subtracted.
//!
//! At |delta_i*h|<=1/2 the omitted Taylor tail is at most
//! 2*b^N/N! * sum_i |A_i|. Admission chooses N<=36 with this exact-arithmetic
//! bound below 1e-51. This does not bound floating arithmetic, phase uncertainty
//! or relative error at an arbitrary coherent null. |C|>=64 keeps every upward
//! moment recurrence multiplier n/|C| at most 35/64. |C|<=2^500 keeps inverse
//! carrier squares normal. Unsupported counts/bands/carriers return None for
//! the caller's original linear kernels; no observation is silently truncated.

use super::{Dd, HarmonicWindow, HarmonicWindowProfile, WindowError, check_amplitude, checked};

/// Maximum original angular-frequency sources in a coherent linear readout.
pub const MAX_HARMONIC_READOUT_TONES: usize = 64;
/// Maximum retained envelope degree plus one.
pub const MAX_HARMONIC_READOUT_TERMS: usize = 36;
/// Fixed temporary centered-field payload of one observation, excluding scalars.
pub const HARMONIC_READOUT_QUERY_SCRATCH_BYTES: usize =
    core::mem::size_of::<[[Dd; 2]; MAX_HARMONIC_READOUT_TONES]>();

/// Fixed reusable coefficients for one original-clock, normalized readout.
///
/// Preparation and observations allocate nothing. Scalar locals, reference
/// metadata, compiler temporaries and call-stack overhead are excluded from
/// the numerical payload. Refusal may alter scratch; the next successful
/// preparation overwrites every admitted coefficient before returning.
pub struct HarmonicReadoutWorkspace {
    rotations: [[Dd; 2]; MAX_HARMONIC_READOUT_TONES],
    coefficients: [Dd; MAX_HARMONIC_READOUT_TONES * MAX_HARMONIC_READOUT_TERMS],
    moments: [Dd; MAX_HARMONIC_READOUT_TERMS],
}

impl HarmonicReadoutWorkspace {
    #[must_use]
    pub const fn new() -> Self {
        Self {
            rotations: [[Dd::ZERO; 2]; MAX_HARMONIC_READOUT_TONES],
            coefficients: [Dd::ZERO; MAX_HARMONIC_READOUT_TONES * MAX_HARMONIC_READOUT_TERMS],
            moments: [Dd::ZERO; MAX_HARMONIC_READOUT_TERMS],
        }
    }
}

impl Default for HarmonicReadoutWorkspace {
    fn default() -> Self {
        Self::new()
    }
}

/// Borrowed linear factor retaining the physical carrier and both DD clock parts.
///
/// It integrates an analytic field, not its power or covariance. Take the real
/// component for a real peak-phasor field. Proper factors use both components
/// of this actual readout when forming covariance of averaged measurements.
#[derive(Clone, Copy)]
pub struct HarmonicReadoutRef<'a> {
    rotations: &'a [[Dd; 2]],
    coefficients: &'a [Dd],
    moments: &'a [Dd],
    carrier_rotation: [Dd; 2],
}

struct Spec {
    carrier: Dd,
    carrier_phase: Dd,
    terms: usize,
}

fn admit(window: HarmonicWindow, frequencies: &[f64]) -> Result<Option<Spec>, WindowError> {
    if frequencies.iter().any(|frequency| !frequency.is_finite()) {
        return Err(WindowError::InvalidFrequency);
    }
    if !(2..=MAX_HARMONIC_READOUT_TONES).contains(&frequencies.len()) {
        return Ok(None);
    }
    let minimum = frequencies.iter().copied().fold(f64::INFINITY, f64::min);
    let maximum = frequencies
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let half = Dd::from_f64(0.5);
    let carrier = checked(Dd::from_f64(minimum) * half + Dd::from_f64(maximum) * half)?;
    let carrier_phase = window.half_width_extended(carrier)?;
    let absolute_phase = carrier_phase.abs();
    let maximum_phase = Dd::from_f64(f64::from_bits(1523_u64 << 52)); // 2^500.
    if absolute_phase.lt(Dd::from_f64(64.0)) || maximum_phase.lt(absolute_phase) {
        return Ok(None);
    }
    let width = window
        .half_width_extended(checked(Dd::from_f64(maximum) - carrier)?)?
        .abs();
    if half.lt(width) {
        return Ok(None);
    }
    let mut tail = Dd::ONE;
    for terms in 1..=MAX_HARMONIC_READOUT_TERMS {
        tail = tail * width / Dd::from_f64(terms as f64);
        if tail.lt(Dd::from_f64(5e-52)) {
            return Ok(Some(Spec {
                carrier,
                carrier_phase,
                terms,
            }));
        }
    }
    // The maximum admitted width needs fewer than 36 terms. Keep refusal
    // explicit if a future bound or arithmetic policy changes that invariant.
    Ok(None)
}

impl HarmonicWindow {
    /// Admit an angular-frequency linear readout without allocating or preparing.
    ///
    /// `None` selects the original linear kernels: fewer than two or more than
    /// 64 tones, envelope half-width above 1/2, or absolute carrier half-phase
    /// outside [64,2^500]. Uniform and triangular sensitivities are normalized.
    /// Nonfinite inputs or unrepresentable phase products remain explicit errors.
    pub fn high_carrier_readout_term_count_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<usize>, WindowError> {
        Ok(admit(self, frequencies)?.map(|spec| spec.terms))
    }

    /// Prepare a coherent original-carrier readout in caller-owned fixed scratch.
    pub fn prepare_high_carrier_readout_rad_s_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut HarmonicReadoutWorkspace,
    ) -> Result<Option<HarmonicReadoutRef<'a>>, WindowError> {
        let Some(spec) = admit(self, frequencies)? else {
            return Ok(None);
        };
        let count = frequencies.len();
        for (source, &frequency) in frequencies.iter().enumerate() {
            let offset = checked(Dd::from_f64(frequency) - spec.carrier)?;
            let width = self.half_width_extended(offset)?;
            // Keep both frequency components, epoch and half-duration
            // rotations separate. Forming one DD frequency*clock product
            // can discard a resolved low-frequency*low-clock phase.
            workspace.rotations[source] = self.rotation_parts(offset)?;
            let mut coefficient = Dd::ONE;
            for term in 0..spec.terms {
                workspace.coefficients[term * count + source] = coefficient;
                coefficient = coefficient * width / Dd::from_f64((term + 1) as f64);
            }
        }
        moments(
            self.profile,
            spec.carrier_phase,
            &mut workspace.moments[..spec.terms],
        )?;
        Ok(Some(HarmonicReadoutRef {
            rotations: &workspace.rotations[..count],
            coefficients: &workspace.coefficients[..count * spec.terms],
            moments: &workspace.moments[..spec.terms],
            carrier_rotation: self.rotation_parts(spec.carrier)?,
        }))
    }
}

// For positive-half integrals J_n=A_n-i*B_n, integration by parts gives
// A_n=(sin(C)-n*B_(n-1))/C and B_n=(n*A_(n-1)-cos(C))/C.
// Uniform full-interval moments use A_even and -i*B_odd. For triangular
// sensitivity K_n=J_n-J_(n+1), avoid subtracting its two bright endpoint
// terms: Aw_n=(B_n-n*Bw_(n-1))/C, Bw_n=(n*Aw_(n-1)-A_n)/C.
// Its full moments are 2*Aw_even and -2i*Bw_odd. Multiplication by (-i)^n
// leaves a real coefficient with signs +,-,-,+ repeated every four terms.
fn moments(
    profile: HarmonicWindowProfile,
    phase: Dd,
    output: &mut [Dd],
) -> Result<(), WindowError> {
    let [cosine, negative_sine] = super::phase_rad_s_dd_from_clock_dd(phase, Dd::ONE)?;
    let sine = -negative_sine;
    let inverse = Dd::ONE / phase;
    let sine_over = sine * inverse;
    let cosine_over = cosine * inverse;
    let mut real = sine_over;
    let mut sine_part = (Dd::ONE - cosine) * inverse;
    let mut weighted_real = sine_part * inverse;
    let mut weighted_sine = (Dd::ONE - real) * inverse;
    // Uniform odd moments need the boundary at -1 as well as +1, whereas
    // J_0's imaginary part contains the positive-half lower endpoint.
    let mut uniform = sine_over;
    for (term, out) in output.iter_mut().enumerate() {
        if term != 0 {
            let n = Dd::from_f64(term as f64);
            let next_real = sine_over - (n * inverse) * sine_part;
            let next_sine = (n * inverse) * real - cosine_over;
            real = next_real;
            sine_part = next_sine;
            let next_weighted_real = (sine_part - n * weighted_sine) * inverse;
            let next_weighted_sine = (n * weighted_real - real) * inverse;
            weighted_real = next_weighted_real;
            weighted_sine = next_weighted_sine;
            uniform = if term % 2 == 0 {
                sine_over - (n * inverse) * uniform
            } else {
                (n * inverse) * uniform - cosine_over
            };
        }
        let value = match profile {
            HarmonicWindowProfile::Uniform => uniform,
            HarmonicWindowProfile::Triangular => {
                Dd::from_f64(2.0)
                    * if term % 2 == 0 {
                        weighted_real
                    } else {
                        weighted_sine
                    }
            }
        };
        *out = checked(if term % 4 == 1 || term % 4 == 2 {
            -value
        } else {
            value
        })?;
    }
    Ok(())
}

impl HarmonicReadoutRef<'_> {
    #[must_use]
    pub fn source_count(self) -> usize {
        self.rotations.len()
    }

    #[must_use]
    pub fn term_count(self) -> usize {
        self.moments.len()
    }

    /// Four real products per source/final rotation and two per projection term.
    #[must_use]
    pub fn mean_real_products(self) -> usize {
        4 * (self.source_count() + 1) + 2 * (self.source_count() + 1) * self.term_count()
    }

    /// Actual complex mean on the original source clock, with no proper half factor.
    pub fn mean_parts_dd(self, fields: &[[Dd; 2]]) -> Result<[Dd; 2], WindowError> {
        if fields.len() != self.source_count() {
            return Err(WindowError::InvalidBandShape);
        }
        let mut centered = [[Dd::ZERO; 2]; MAX_HARMONIC_READOUT_TONES];
        for (out, (&[ar, ai], &[cr, ci])) in
            centered.iter_mut().zip(fields.iter().zip(self.rotations))
        {
            check_amplitude(ar)?;
            check_amplitude(ai)?;
            *out = [checked(ar * cr - ai * ci)?, checked(ar * ci + ai * cr)?];
        }
        let mut result = [Dd::ZERO; 2];
        for (term, &moment) in self.moments.iter().enumerate() {
            let mut mode = [Dd::ZERO; 2];
            for (source, field) in centered[..fields.len()].iter().enumerate() {
                let coefficient = self.coefficients[term * fields.len() + source];
                mode[0] = mode[0] + coefficient * field[0];
                mode[1] = mode[1] + coefficient * field[1];
            }
            result[0] = result[0] + checked(mode[0])? * moment;
            result[1] = result[1] + checked(mode[1])? * moment;
        }
        let [cr, ci] = self.carrier_rotation;
        Ok([
            checked(result[0] * cr - result[1] * ci)?,
            checked(result[0] * ci + result[1] * cr)?,
        ])
    }
}

#[cfg(test)]
mod tests;
