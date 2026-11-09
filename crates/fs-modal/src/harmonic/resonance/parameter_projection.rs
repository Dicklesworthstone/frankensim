//! Complete coherent receiver jets for peak and cancellation-notch sensitivities.

use super::super::{HarmonicResponseError, finite, invalid, parameter::PowerSum};
use super::{HarmonicTransferParameterJet, HarmonicTransferSample};
use fs_math::c64::C64;

/// Borrowed field derivatives in one fixed physical coordinate ordering.
///
/// Derivatives are partial derivatives at fixed angular frequency or fixed
/// parameter, as indicated. In particular, the mixed derivative is not the
/// total derivative along a moving resonance. Every slice has the same
/// nonzero length. This carrier can describe displacement, velocity, pressure
/// or another field whose physical owner supplies the complete analytic jet.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicTransferFieldParameterJet<'a> {
    /// Original complex field.
    pub values: &'a [C64],
    /// Partial derivative per rad/s.
    pub angular_derivative: &'a [C64],
    /// Second partial derivative per (rad/s)^2.
    pub angular_second_derivative: &'a [C64],
    /// Partial derivative per caller-declared unit of the physical parameter.
    pub parameter_derivative: &'a [C64],
    /// Mixed partial derivative per (rad/s) per parameter unit.
    pub angular_parameter_derivative: &'a [C64],
}

/// A coherent complex receiver and its authored frequency/parameter law.
///
/// Weights act as a transfer row, without conjugation. Each optional derivative
/// is either a complete slice in the same ordering as weights, or None to
/// explicitly declare that derivative zero at this frequency and parameter.
/// A changing sensor, phase delay or radiation observation must supply its
/// complete derivatives; they are never inferred from the field or samples.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReceiverParameterJet<'a> {
    /// Complex linear transfer weights.
    pub weights: &'a [C64],
    /// Partial receiver derivative per rad/s.
    pub angular_derivative: Option<&'a [C64]>,
    /// Second partial receiver derivative per (rad/s)^2.
    pub angular_second_derivative: Option<&'a [C64]>,
    /// Partial receiver derivative per parameter unit.
    pub parameter_derivative: Option<&'a [C64]>,
    /// Mixed receiver derivative per (rad/s) per parameter unit.
    pub angular_parameter_derivative: Option<&'a [C64]>,
}

impl HarmonicTransferParameterJet {
    /// Observe a complete field jet through a fixed real linear receiver.
    ///
    /// Signed cancellation is retained before rounding. The work is linear in
    /// the number of coordinates, with constant scratch and no allocation.
    /// The receiver and coordinate basis are held fixed in both variables.
    pub fn project_real(
        receiver: &[f64],
        field: HarmonicTransferFieldParameterJet<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        if receiver.len() != field.values.len() {
            return Err(invalid(
                "one real receiver weight per field coordinate required",
            ));
        }
        project(field, |i| {
            [
                C64::from_re(receiver[i]),
                C64::ZERO,
                C64::ZERO,
                C64::ZERO,
                C64::ZERO,
            ]
        })
    }

    /// Observe a complete field jet through fixed complex transfer weights.
    ///
    /// This is a linear transfer row, not a Hermitian inner product.
    pub fn project_complex(
        receiver: &[C64],
        field: HarmonicTransferFieldParameterJet<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        if receiver.len() != field.values.len() {
            return Err(invalid(
                "one complex receiver weight per field coordinate required",
            ));
        }
        project(field, |i| {
            [receiver[i], C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO]
        })
    }

    /// Observe a field and a receiver that both depend on frequency and parameter.
    ///
    /// For H = sum(r*q), this retains the complete product rules:
    /// H_w = sum(r_w*q + r*q_w),
    /// H_ww = sum(r_ww*q + 2*r_w*q_w + r*q_ww),
    /// H_t = sum(r_t*q + r*q_t), and
    /// H_wt = sum(r_wt*q + r_w*q_t + r_t*q_w + r*q_wt).
    ///
    /// Every signed complex product uses the harmonic owner's binary-normalized
    /// double-double accumulator. Finite cancellation can survive individually
    /// overflowing products; only the final physical observation is recovered.
    /// No source, load, dissipation or physical work is changed by observation.
    pub fn project_complex_with_derivatives(
        receiver: HarmonicReceiverParameterJet<'_>,
        field: HarmonicTransferFieldParameterJet<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        let n = field.values.len();
        if receiver.weights.len() != n
            || [
                receiver.angular_derivative,
                receiver.angular_second_derivative,
                receiver.parameter_derivative,
                receiver.angular_parameter_derivative,
            ]
            .into_iter()
            .flatten()
            .any(|values| values.len() != n)
        {
            return Err(invalid(
                "complete receiver derivative slices in field order required",
            ));
        }
        project(field, |i| {
            [
                receiver.weights[i],
                receiver.angular_derivative.map_or(C64::ZERO, |v| v[i]),
                receiver
                    .angular_second_derivative
                    .map_or(C64::ZERO, |v| v[i]),
                receiver.parameter_derivative.map_or(C64::ZERO, |v| v[i]),
                receiver
                    .angular_parameter_derivative
                    .map_or(C64::ZERO, |v| v[i]),
            ]
        })
    }
}

#[derive(Default)]
struct ComplexSum {
    real: PowerSum,
    imaginary: PowerSum,
}

impl ComplexSum {
    fn product(
        &mut self,
        receiver: C64,
        field: C64,
        multiplicity: f64,
    ) -> Result<(), HarmonicResponseError> {
        // half_inner normalizes all four real factors before multiplying;
        // neither a large coherent product nor a subnormal input is rounded
        // prematurely by first materializing receiver*field or dividing by 2.
        self.real.half_inner(
            C64::new(receiver.re, -receiver.im),
            field,
            2.0 * multiplicity,
        )?;
        self.imaginary.half_inner(
            C64::new(receiver.im, receiver.re),
            field,
            2.0 * multiplicity,
        )
    }

    fn value(self) -> Result<C64, HarmonicResponseError> {
        let value = C64::new(self.real.value()?, self.imaginary.value()?);
        if !finite(value) {
            return Err(invalid(
                "physical transfer parameter jet is not representable",
            ));
        }
        Ok(value)
    }
}

fn project(
    field: HarmonicTransferFieldParameterJet<'_>,
    mut receiver_at: impl FnMut(usize) -> [C64; 5],
) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
    let n = field.values.len();
    if n == 0
        || [
            field.angular_derivative,
            field.angular_second_derivative,
            field.parameter_derivative,
            field.angular_parameter_derivative,
        ]
        .iter()
        .any(|v| v.len() != n)
    {
        return Err(invalid(
            "matching nonempty physical field derivative slices required",
        ));
    }
    let mut sums: [ComplexSum; 5] = core::array::from_fn(|_| ComplexSum::default());
    for i in 0..n {
        let x = [
            field.values[i],
            field.angular_derivative[i],
            field.angular_second_derivative[i],
            field.parameter_derivative[i],
            field.angular_parameter_derivative[i],
        ];
        let r = receiver_at(i);
        if x.into_iter().chain(r).any(|z| !finite(z)) {
            return Err(invalid(
                "finite physical receiver and field parameter jets required",
            ));
        }
        for axis in 0..5 {
            sums[axis].product(r[0], x[axis], 1.0)?;
        }
        sums[1].product(r[1], x[0], 1.0)?;
        sums[2].product(r[2], x[0], 1.0)?;
        sums[2].product(r[1], x[1], 2.0)?;
        sums[3].product(r[3], x[0], 1.0)?;
        sums[4].product(r[4], x[0], 1.0)?;
        sums[4].product(r[1], x[3], 1.0)?;
        sums[4].product(r[3], x[1], 1.0)?;
    }
    let [value, angular, angular_second, parameter, mixed] = sums;
    Ok(HarmonicTransferParameterJet {
        sample: HarmonicTransferSample {
            amplitude: value.value()?,
            angular_derivative: angular.value()?,
        },
        angular_second_derivative: angular_second.value()?,
        parameter_derivative: parameter.value()?,
        angular_parameter_derivative: mixed.value()?,
    })
}
