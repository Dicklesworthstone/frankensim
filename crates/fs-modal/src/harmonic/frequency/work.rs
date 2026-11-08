//! Second derivatives of signed physical work, with products scaled before summation.

use super::super::{HarmonicResponseError, finite, invalid, parameter::PowerSum};
use fs_math::c64::C64;

#[derive(Default)]
pub(super) struct ComplexSum {
    real: PowerSum,
    imaginary: PowerSum,
}

impl ComplexSum {
    pub(super) fn product(
        &mut self,
        a: C64,
        b: C64,
        multiplicity: f64,
    ) -> Result<(), HarmonicResponseError> {
        self.real
            .half_inner(C64::new(a.re, -a.im), b, 2.0 * multiplicity)?;
        self.imaginary
            .half_inner(C64::new(a.im, a.re), b, 2.0 * multiplicity)
    }

    pub(super) fn value(self) -> Result<C64, HarmonicResponseError> {
        Ok(C64::new(self.real.value()?, self.imaginary.value()?))
    }
}

/// -i(2*q_omega + omega*q_omegaomega), with cancellation before range recovery.
pub(super) fn velocity_curvature(
    omega: f64,
    first: C64,
    second: C64,
) -> Result<C64, HarmonicResponseError> {
    let mut sum = ComplexSum::default();
    sum.product(C64::from_re(omega), second, 1.0)?;
    sum.product(C64::from_re(2.0), first, 1.0)?;
    let value = sum.value()?;
    Ok(C64::new(value.im, -value.re))
}

/// Second derivative of `Re(flow^H effort)/2` in the caller's parameter units.
///
/// All arrays are original power-conjugate phasors or their first/second
/// derivatives in the same ordering. Includes both cross terms; signed
/// boundary work is not an inferred constitutive heat or storage split.
/// Binary-normalized double-double products retain finite totals when the
/// individual differentiated products overflow. No temporary arrays are used.
/// Cancellation and invalid/unrepresentable input return no partial result.
#[allow(clippy::too_many_arguments)]
pub fn harmonic_mean_power_curvature_with_cancel(
    flow: &[C64],
    flow_derivative: &[C64],
    flow_second_derivative: &[C64],
    effort: &[C64],
    effort_derivative: &[C64],
    effort_second_derivative: &[C64],
    mut cancelled: impl FnMut() -> bool,
) -> Result<f64, HarmonicResponseError> {
    if cancelled() {
        return Err(HarmonicResponseError::Cancelled);
    }
    let n = flow.len();
    if n == 0
        || [
            flow_derivative,
            flow_second_derivative,
            effort,
            effort_derivative,
            effort_second_derivative,
        ]
        .iter()
        .any(|v| v.len() != n)
    {
        return Err(invalid(
            "matching nonempty power-conjugate frequency jets required",
        ));
    }
    let mut total = PowerSum::default();
    for i in 0..n {
        if cancelled() {
            return Err(HarmonicResponseError::Cancelled);
        }
        if [
            flow[i],
            flow_derivative[i],
            flow_second_derivative[i],
            effort[i],
            effort_derivative[i],
            effort_second_derivative[i],
        ]
        .into_iter()
        .any(|z| !finite(z))
        {
            return Err(invalid(
                "finite physical power phasors and derivatives required",
            ));
        }
        total.half_inner(flow_second_derivative[i], effort[i], 1.0)?;
        total.half_inner(flow_derivative[i], effort_derivative[i], 2.0)?;
        total.half_inner(flow[i], effort_second_derivative[i], 1.0)?;
    }
    if cancelled() {
        return Err(HarmonicResponseError::Cancelled);
    }
    total.value()
}

#[cfg(test)]
#[path = "work_tests.rs"]
mod tests;
