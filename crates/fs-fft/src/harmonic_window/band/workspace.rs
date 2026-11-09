//! Bounded preparation without heap allocation, sharing all cached contractions.

use super::basis::{self, Units};
use super::*;

/// Reusable fixed storage for one positive temporal factor, with no heap allocation.
///
/// Its 19,456-byte numerical payload holds at most 64 original source rotations
/// and 64*17 coefficients. A prepared borrow prevents reuse until observation
/// is complete. Refusal can alter scratch; the next preparation overwrites every
/// admitted entry before returning a new borrow.
/// Profile-aware triangular preparation uses an additional bounded 6,864-byte
/// numerical scratch payload for its polynomial recurrence and Taylor rows;
/// it changes neither this retained payload nor per-observation storage.
pub struct UniformHarmonicBandWorkspace {
    rotations: [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES],
    coefficients: [Dd; MAX_UNIFORM_HARMONIC_BAND_TONES * MAX_UNIFORM_HARMONIC_BAND_MODES],
}

impl UniformHarmonicBandWorkspace {
    /// Create zeroed fixed scratch; preparation and observations allocate nothing.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            rotations: [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES],
            coefficients: [Dd::ZERO;
                MAX_UNIFORM_HARMONIC_BAND_TONES * MAX_UNIFORM_HARMONIC_BAND_MODES],
        }
    }
}

impl Default for UniformHarmonicBandWorkspace {
    fn default() -> Self {
        Self::new()
    }
}

/// Immutable temporal factor borrowing a caller-owned fixed workspace.
///
/// Analytic preparation removes a common carrier only from second moments.
/// Real preparation retains the absolute carrier and uses conjugate parity:
/// even modes retain the rotated real amplitude, odd modes the imaginary one.
/// This equals the mirrored +/-frequency real-field expansion without halving
/// tiny amplitudes or doubling source storage. Its mode-zero real coefficient
/// is the actual window mean. Other modes omit their common (-i)^mode phase.
#[derive(Clone, Copy)]
pub struct UniformHarmonicBandRef<'a> {
    pub(super) rotations: &'a [[Dd; 2]],
    pub(super) coefficients: &'a [Dd],
    pub(super) modes: usize,
    pub(super) real_field: bool,
}

fn prepare_in<'a>(
    window: HarmonicWindow,
    frequencies: &[f64],
    units: Units,
    real_field: bool,
    workspace: &'a mut UniformHarmonicBandWorkspace,
) -> Result<Option<UniformHarmonicBandRef<'a>>, WindowError> {
    let Some(spec) = basis::admit(window, frequencies, units, real_field)? else {
        return Ok(None);
    };
    prepare_spec_in(window, frequencies, units, spec, workspace)
}

pub(super) fn prepare_profile_in<'a>(
    window: HarmonicWindow,
    frequencies: &[f64],
    units: Units,
    real_field: bool,
    workspace: &'a mut UniformHarmonicBandWorkspace,
) -> Result<Option<UniformHarmonicBandRef<'a>>, WindowError> {
    let Some(spec) = basis::admit_profile(window, frequencies, units, real_field)? else {
        return Ok(None);
    };
    prepare_spec_in(window, frequencies, units, spec, workspace)
}

fn prepare_spec_in<'a>(
    window: HarmonicWindow,
    frequencies: &[f64],
    units: Units,
    spec: basis::Spec,
    workspace: &'a mut UniformHarmonicBandWorkspace,
) -> Result<Option<UniformHarmonicBandRef<'a>>, WindowError> {
    let rotations = &mut workspace.rotations[..frequencies.len()];
    let coefficients = &mut workspace.coefficients[..frequencies.len() * spec.modes];
    basis::fill(window, frequencies, units, &spec, rotations, coefficients)?;
    Ok(Some(UniformHarmonicBandRef {
        rotations,
        coefficients,
        modes: spec.modes,
        real_field: spec.real_field,
    }))
}

impl HarmonicWindow {
    /// Prepare an analytic Hz band in fixed scratch, without allocating.
    pub fn prepare_uniform_band_hz_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut UniformHarmonicBandWorkspace,
    ) -> Result<Option<UniformHarmonicBandRef<'a>>, WindowError> {
        prepare_in(self, frequencies, Units::Hertz, false, workspace)
    }

    /// Prepare an analytic original-rad/s band in fixed scratch, without allocating.
    pub fn prepare_uniform_band_rad_s_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut UniformHarmonicBandWorkspace,
    ) -> Result<Option<UniformHarmonicBandRef<'a>>, WindowError> {
        prepare_in(self, frequencies, Units::Radians, false, workspace)
    }

    /// Prepare actual real-peak Hz products in fixed scratch; no half factor.
    /// The same 3..32-tone full-frequency admission as the cached real factor applies.
    pub fn prepare_uniform_real_band_hz_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut UniformHarmonicBandWorkspace,
    ) -> Result<Option<UniformHarmonicBandRef<'a>>, WindowError> {
        prepare_in(self, frequencies, Units::Hertz, true, workspace)
    }

    /// Prepare actual real-peak original-rad/s products, without allocating.
    /// Sum-frequency contributions remain; no carrier-period averaging is applied.
    pub fn prepare_uniform_real_band_rad_s_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut UniformHarmonicBandWorkspace,
    ) -> Result<Option<UniformHarmonicBandRef<'a>>, WindowError> {
        prepare_in(self, frequencies, Units::Radians, true, workspace)
    }
}

impl UniformHarmonicBandRef<'_> {
    /// Original source count; real parity does not double this count.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.rotations.len()
    }

    /// Number of retained orthogonal temporal modes.
    #[must_use]
    pub const fn mode_count(&self) -> usize {
        self.modes
    }

    /// Logical source projection terms plus one squared norm per mode.
    #[must_use]
    pub fn power_work(&self) -> usize {
        (self.source_count() + 1) * self.modes
    }

    /// Two projections plus one complex cross moment per mode.
    #[must_use]
    pub fn cross_work(&self) -> usize {
        (2 * self.source_count() + 1) * self.modes
    }

    /// Two projections plus the four real entries of a coherency per mode.
    #[must_use]
    pub fn coherency_work(&self) -> usize {
        (2 * self.source_count() + 4) * self.modes
    }

    /// Project original field components into unweighted temporal modes.
    ///
    /// Fields must have `source_count()` entries and output `mode_count()`.
    /// Real preparation accepts original peak phasors directly. Its even
    /// outputs are real, odd outputs imaginary, and mode zero is the actual
    /// field mean. Scratch can change on refusal.
    pub fn project_parts_dd(
        &self,
        fields: &[[Dd; 2]],
        output: &mut [[Dd; 2]],
    ) -> Result<(), WindowError> {
        if fields.len() != self.source_count() || output.len() != self.modes {
            return Err(WindowError::InvalidBandShape);
        }
        let mut centered = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
        for (i, (&[ar, ai], &[cr, ci])) in fields.iter().zip(self.rotations).enumerate() {
            check_amplitude(ar)?;
            check_amplitude(ai)?;
            centered[i] = [checked(ar * cr - ai * ci)?, checked(ar * ci + ai * cr)?];
        }
        for (mode, out) in output.iter_mut().enumerate() {
            *out = [Dd::ZERO; 2];
            for (i, field) in centered[..fields.len()].iter().enumerate() {
                let coefficient = self.coefficients[mode * fields.len() + i];
                if self.real_field {
                    let component = mode % 2;
                    out[component] = out[component] + coefficient * field[component];
                } else {
                    out[0] = out[0] + coefficient * field[0];
                    out[1] = out[1] + coefficient * field[1];
                }
            }
            checked(out[0])?;
            checked(out[1])?;
        }
        Ok(())
    }

    /// `(2*mode+1)*left*conj(right)` for one retained projected mode.
    pub fn mode_cross_parts_dd(
        &self,
        mode: usize,
        left: [Dd; 2],
        right: [Dd; 2],
    ) -> Result<[Dd; 2], WindowError> {
        self.mode_cross_parts_dd_scaled(mode, left, right, Dd::ONE)
    }

    /// Apply a finite positive physical scale inside the weighted modal product.
    ///
    /// Analytic proper-complex peak factors use scale 0.5. Actual real factors
    /// use scale 1. Applying the scale before multiplication preserves a finite
    /// half-power even when its unscaled analytic squared norm would overflow.
    pub fn mode_cross_parts_dd_scaled(
        &self,
        mode: usize,
        left: [Dd; 2],
        right: [Dd; 2],
        positive_scale: Dd,
    ) -> Result<[Dd; 2], WindowError> {
        if mode >= self.modes {
            return Err(WindowError::InvalidBandShape);
        }
        check_scale(positive_scale)?;
        for amplitude in left.into_iter().chain(right) {
            check_amplitude(amplitude)?;
        }
        weighted_cross(
            left,
            right,
            checked(Dd::from_f64((2 * mode + 1) as f64) * positive_scale)?,
        )
    }

    /// Nonnegative mean squared norm of the prepared analytic or actual real field.
    pub fn mean_square_parts_dd(&self, fields: &[[Dd; 2]]) -> Result<Dd, WindowError> {
        self.mean_square_parts_dd_scaled(fields, Dd::ONE)
    }

    /// Mean squared norm with a finite positive scale inside every modal product.
    pub fn mean_square_parts_dd_scaled(
        &self,
        fields: &[[Dd; 2]],
        positive_scale: Dd,
    ) -> Result<Dd, WindowError> {
        check_scale(positive_scale)?;
        let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        self.project_parts_dd(fields, &mut modes[..self.modes])?;
        let mut result = Dd::ZERO;
        for (l, &field) in modes[..self.modes].iter().enumerate() {
            let weight = checked(Dd::from_f64((2 * l + 1) as f64) * positive_scale)?;
            result = checked(result + weighted_cross(field, field, weight)?[0])?;
        }
        Ok(result)
    }

    /// Mean `left*conj(right)`; actual real preparation retains all sum tones.
    pub fn mean_cross_parts_dd(
        &self,
        left: &[[Dd; 2]],
        right: &[[Dd; 2]],
    ) -> Result<[Dd; 2], WindowError> {
        self.mean_cross_parts_dd_scaled(left, right, Dd::ONE)
    }

    /// Complex cross moment with a finite positive scale inside every product.
    pub fn mean_cross_parts_dd_scaled(
        &self,
        left: &[[Dd; 2]],
        right: &[[Dd; 2]],
        positive_scale: Dd,
    ) -> Result<[Dd; 2], WindowError> {
        check_scale(positive_scale)?;
        let mut a = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        let mut b = a;
        self.project_parts_dd(left, &mut a[..self.modes])?;
        self.project_parts_dd(right, &mut b[..self.modes])?;
        let mut result = [Dd::ZERO; 2];
        for l in 0..self.modes {
            let weight = checked(Dd::from_f64((2 * l + 1) as f64) * positive_scale)?;
            let cross = weighted_cross(a[l], b[l], weight)?;
            result[0] = checked(result[0] + cross[0])?;
            result[1] = checked(result[1] + cross[1])?;
        }
        Ok(result)
    }
}

fn check_scale(scale: Dd) -> Result<(), WindowError> {
    check_amplitude(scale)?;
    if Dd::ZERO.lt(scale) {
        Ok(())
    } else {
        Err(WindowError::Unrepresentable)
    }
}
