//! Normalized first-order coherent observations in original physical units.
//!
//! Match the complete parameter-jet owner: each real factor stays normalized
//! until all signed terms have been accumulated. Receiver weights form a
//! linear transfer row, without complex conjugation. Scratch is constant.

use super::super::{HarmonicResponseError, finite, invalid, parameter::PowerSum};
use super::HarmonicTransferSample;
use fs_math::c64::C64;

pub(super) fn real(
    receiver: &[f64],
    values: &[C64],
    derivatives: &[C64],
) -> Result<HarmonicTransferSample, HarmonicResponseError> {
    if receiver.len() != values.len() {
        return Err(invalid(
            "one finite fixed receiver weight per coordinate required",
        ));
    }
    project(values, derivatives, |i| {
        (C64::from_re(receiver[i]), C64::ZERO)
    })
}

pub(super) fn complex(
    receiver: &[C64],
    receiver_angular_derivative: Option<&[C64]>,
    values: &[C64],
    derivatives: &[C64],
) -> Result<HarmonicTransferSample, HarmonicResponseError> {
    if receiver.len() != values.len()
        || receiver_angular_derivative.is_some_and(|slope| slope.len() != values.len())
    {
        return Err(invalid(
            "one finite complex receiver weight and slope per coordinate required",
        ));
    }
    project(values, derivatives, |i| {
        (
            receiver[i],
            receiver_angular_derivative.map_or(C64::ZERO, |slope| slope[i]),
        )
    })
}

fn project(
    values: &[C64],
    derivatives: &[C64],
    mut receiver_at: impl FnMut(usize) -> (C64, C64),
) -> Result<HarmonicTransferSample, HarmonicResponseError> {
    let n = values.len();
    if n == 0 || derivatives.len() != n {
        return Err(invalid(
            "matching nonempty physical field and frequency derivative slices required",
        ));
    }
    let mut sums: [PowerSum; 4] = core::array::from_fn(|_| PowerSum::default());
    for i in 0..n {
        let (receiver, receiver_derivative) = receiver_at(i);
        let value = values[i];
        let derivative = derivatives[i];
        if [receiver, receiver_derivative, value, derivative]
            .into_iter()
            .any(|z| !finite(z))
        {
            return Err(invalid(
                "finite physical receiver and field frequency jets required",
            ));
        }
        // H = r*x; H_omega = r*x_omega + r_omega*x.
        // half_inner(..., 2) normalizes every factor before multiplying;
        // it never forms either the raw complex product or a premature half.
        // Keep the same term order as the complete parameter-jet projection.
        for (slot, weight, field) in [
            (0, receiver, value),
            (2, receiver, derivative),
            (2, receiver_derivative, value),
        ] {
            sums[slot].half_inner(C64::new(weight.re, -weight.im), field, 2.0)?;
            sums[slot + 1].half_inner(C64::new(weight.im, weight.re), field, 2.0)?;
        }
    }
    let [real, imaginary, slope_real, slope_imaginary] = sums;
    let amplitude = C64::new(real.value()?, imaginary.value()?);
    let angular_derivative = C64::new(slope_real.value()?, slope_imaginary.value()?);
    if !finite(amplitude) || !finite(angular_derivative) {
        return Err(invalid(
            "physical transfer or its derivative is not representable",
        ));
    }
    Ok(HarmonicTransferSample {
        amplitude,
        angular_derivative,
    })
}
