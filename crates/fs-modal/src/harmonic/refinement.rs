//! Bounded residual correction with the already admitted physical factor.
//!
//! A weak observation near an antiresonance can lose relative accuracy during
//! back substitution even when the complete pencil is well conditioned. Form
//! its original equation residual in DD, solve that correction with the SAME
//! factor, and retain the caller's componentwise residual threshold.

use super::parameter::PowerSum;
use super::{
    HarmonicResponseError, PreparedHarmonicResponse, finite, invalid, poll, scale_complex,
};
use crate::quadratic_fraction;
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd};

/// One correction only; the caller bounds passes and re-audits the original
/// equation before publication. Neither the factor nor physical inputs change.
pub(super) fn correct(
    prepared: &PreparedHarmonicResponse,
    solution: &mut [C64],
    rhs: &[C64],
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    poll(cx)?;
    let exponent = solution
        .iter()
        .chain(rhs)
        .flat_map(|value| [value.re, value.im])
        .filter(|value| *value != 0.0)
        .map(|value| quadratic_fraction(value).1)
        .max()
        .unwrap_or(0);
    // Exact binary normalization avoids squaring a dimensional response or
    // forming dimensional matrix products before cancellation. Two O(n)
    // scratch vectors; no second factor or dense operator is constructed.
    let mut normalized = Vec::with_capacity(prepared.n);
    for &value in solution.iter() {
        poll(cx)?;
        normalized.push(scale_complex(value, -exponent)?);
    }
    let mut correction = Vec::with_capacity(prepared.n);
    for i in 0..prepared.n {
        poll(cx)?;
        let b = scale_complex(rhs[i], -exponent)?;
        let mut real = PowerSum::default();
        let mut imaginary = PowerSum::default();
        real.add_scalar(b.re)?;
        imaginary.add_scalar(b.im)?;
        for j in prepared.operator_columns(i) {
            let a = prepared.operator_entry(i, j)?;
            let x = normalized[j];
            // b - A*x, keeping all four signed quadrature products until
            // cancellation. The factor 2 exactly removes half_inner's 1/2.
            real.half_inner(C64::new(-a.re, a.im), x, 2.0)?;
            imaginary.half_inner(C64::new(-a.im, -a.re), x, 2.0)?;
        }
        correction.push(C64::new(real.value()?, imaginary.value()?));
    }
    poll(cx)?;
    prepared.solve_factored(&mut correction, cx)?;
    for i in 0..prepared.n {
        poll(cx)?;
        let corrected = C64::new(
            (Dd::from_f64(normalized[i].re) + Dd::from_f64(correction[i].re)).to_f64(),
            (Dd::from_f64(normalized[i].im) + Dd::from_f64(correction[i].im)).to_f64(),
        );
        if !finite(corrected) {
            return Err(invalid("harmonic residual correction is not representable"));
        }
        solution[i] = scale_complex(corrected, exponent)?;
    }
    poll(cx)
}
