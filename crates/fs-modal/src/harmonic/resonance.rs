//! Resolve peaks and notches of an actual complex forced transfer.
//!
//! The stationary condition is `Re(conj(H) * dH/domega) = 0`. It retains
//! interference between modes, source phase, nonproportional damping and any
//! dispersive load already present in the evaluator. It is not a pole search:
//! displacement, velocity and pressure can peak at different frequencies.
//!
//! The caller supplies a continuous differentiable transfer and either brackets
//! ONE isolated extremum/crossing or supplies a grid resolving multiple peaks
//! and notches. Opposite endpoint signs are checked; uniqueness and absence of
//! unresolved poles are not certified. Refinement is local, without a spectral
//! completeness claim, damping floor or assumed Lorentzian line shape. Each evaluation
//! keeps its original solve/domain/work limits, and cancellation publishes no
//! partial result. The returned bracket is numerical sign evidence, not an
//! interval enclosure accounting for the evaluator's forward error.

use super::{HarmonicFrequencyResponse, HarmonicResponseError, finite, invalid};
use fs_exec::Cx;
use fs_math::c64::C64;

mod query;
pub use query::{
    HalfPowerBandwidth, HarmonicTransferSample, ResonanceControls, ResonanceError,
    TransferExtremum, TransferExtremumKind, TransferSurvey,
    refine_half_power_bandwidth_with_cancel, refine_transfer_extremum_with_cancel,
    survey_transfer_extrema_with_cancel,
};

mod sensitivity;
pub use sensitivity::{
    HalfPowerBandwidthParameterResponse, HarmonicTransferParameterJet,
    TransferExtremumParameterResponse,
    refine_half_power_bandwidth_with_parameter_derivative_with_cancel,
    refine_transfer_extremum_with_parameter_derivative_with_cancel,
};
mod native_parameter_projection;
mod parameter_projection;
pub use parameter_projection::{HarmonicReceiverParameterJet, HarmonicTransferFieldParameterJet};

impl HarmonicTransferSample {
    /// Observe any fixed real linear receiver from a complete physical field
    /// and its frequency tangent. Values and weights retain their caller-owned
    /// units and ordering; source/receiver frequency dependence must already be
    /// included in the supplied tangent. Compensated accumulation preserves
    /// weak surviving terms beside large coherent cancellations.
    pub fn project_real(
        receiver: &[f64],
        values: &[C64],
        derivatives: &[C64],
    ) -> Result<Self, HarmonicResponseError> {
        project(receiver, values, derivatives)
    }

    /// Observe a fixed complex linear receiver, such as a coherent phased
    /// sensor combination. The receiver is applied as a transfer row,
    /// `sum receiver[i] * values[i]` (not a Hermitian inner product).
    /// Frequency dependence belongs in `project_complex_with_derivative`.
    pub fn project_complex(
        receiver: &[C64],
        values: &[C64],
        derivatives: &[C64],
    ) -> Result<Self, HarmonicResponseError> {
        project_complex(receiver, None, values, derivatives)
    }

    /// Observe a complex receiver whose transfer weights vary with angular
    /// frequency. `receiver_angular_derivative` is d(receiver)/domega in the
    /// same ordering, and the returned tangent retains both chain-rule terms:
    /// `sum(receiver' * value + receiver * value')`.
    ///
    /// This is observation only: it changes no mechanical solve, load,
    /// dissipation or source work. The caller owns the receiver's physical
    /// phase/delay/filter law and its frequency support.
    pub fn project_complex_with_derivative(
        receiver: &[C64],
        receiver_angular_derivative: &[C64],
        values: &[C64],
        derivatives: &[C64],
    ) -> Result<Self, HarmonicResponseError> {
        project_complex(
            receiver,
            Some(receiver_angular_derivative),
            values,
            derivatives,
        )
    }
}

impl HarmonicFrequencyResponse {
    /// Observe a fixed real linear displacement receiver, including signed
    /// cancellation between modes. This uses the complete solved frequency
    /// tangent, not natural frequencies or modal amplitude magnitudes.
    pub fn displacement_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_real(
            receiver,
            &self.response.displacement,
            &self.displacement_derivative,
        )
    }

    /// Observe a fixed real linear velocity receiver. The derivative includes
    /// both `-i q` and `-i omega q_omega`, so its peak need not equal the
    /// displacement peak of the same structure and forcing.
    pub fn velocity_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_real(
            receiver,
            &self.response.velocity,
            &self.velocity_derivative,
        )
    }

    /// Observe displacement through fixed complex transfer weights.
    pub fn displacement_transfer_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_complex(
            receiver,
            &self.response.displacement,
            &self.displacement_derivative,
        )
    }

    /// Observe displacement through a frequency-dependent complex transfer row.
    pub fn displacement_transfer_complex_with_derivative(
        &self,
        receiver: &[C64],
        receiver_angular_derivative: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_complex_with_derivative(
            receiver,
            receiver_angular_derivative,
            &self.response.displacement,
            &self.displacement_derivative,
        )
    }

    /// Observe velocity through fixed complex transfer weights.
    pub fn velocity_transfer_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_complex(
            receiver,
            &self.response.velocity,
            &self.velocity_derivative,
        )
    }

    /// Observe velocity through a frequency-dependent complex transfer row.
    pub fn velocity_transfer_complex_with_derivative(
        &self,
        receiver: &[C64],
        receiver_angular_derivative: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_complex_with_derivative(
            receiver,
            receiver_angular_derivative,
            &self.response.velocity,
            &self.velocity_derivative,
        )
    }
}

fn project(
    receiver: &[f64],
    values: &[C64],
    derivatives: &[C64],
) -> Result<HarmonicTransferSample, HarmonicResponseError> {
    if receiver.len() != values.len()
        || derivatives.len() != values.len()
        || receiver.is_empty()
        || receiver.iter().any(|x| !x.is_finite())
    {
        return Err(invalid(
            "one finite fixed receiver weight per coordinate required",
        ));
    }
    let mut sums = [0.0; 4];
    let mut correction = [0.0; 4];
    for ((&weight, &value), &derivative) in receiver.iter().zip(values).zip(derivatives) {
        let terms = [
            weight * value.re,
            weight * value.im,
            weight * derivative.re,
            weight * derivative.im,
        ];
        for ((sum, error), term) in sums.iter_mut().zip(&mut correction).zip(terms) {
            let next = *sum + term;
            *error += if sum.abs() >= term.abs() {
                (*sum - next) + term
            } else {
                (term - next) + *sum
            };
            *sum = next;
        }
    }
    let amplitude = C64::new(sums[0] + correction[0], sums[1] + correction[1]);
    let angular_derivative = C64::new(sums[2] + correction[2], sums[3] + correction[3]);
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

fn project_complex(
    receiver: &[C64],
    receiver_angular_derivative: Option<&[C64]>,
    values: &[C64],
    derivatives: &[C64],
) -> Result<HarmonicTransferSample, HarmonicResponseError> {
    if receiver.len() != values.len()
        || receiver_angular_derivative.is_some_and(|slope| slope.len() != values.len())
        || derivatives.len() != values.len()
        || receiver.is_empty()
        || receiver.iter().any(|z| !finite(*z))
        || receiver_angular_derivative.is_some_and(|slope| slope.iter().any(|z| !finite(*z)))
    {
        return Err(invalid(
            "one finite complex receiver weight and slope per coordinate required",
        ));
    }
    let mut sums = [0.0; 4];
    let mut correction = [0.0; 4];
    let mut add = |slot: usize, term: f64| {
        let next = sums[slot] + term;
        correction[slot] += if sums[slot].abs() >= term.abs() {
            (sums[slot] - next) + term
        } else {
            (term - next) + sums[slot]
        };
        sums[slot] = next;
    };
    for (i, ((&weight, &value), &derivative)) in
        receiver.iter().zip(values).zip(derivatives).enumerate()
    {
        let weight_derivative = receiver_angular_derivative.map_or(C64::ZERO, |slope| slope[i]);
        add(0, weight.re * value.re);
        add(0, -weight.im * value.im);
        add(1, weight.re * value.im);
        add(1, weight.im * value.re);

        add(2, weight_derivative.re * value.re);
        add(2, -weight_derivative.im * value.im);
        add(2, weight.re * derivative.re);
        add(2, -weight.im * derivative.im);
        add(3, weight_derivative.re * value.im);
        add(3, weight_derivative.im * value.re);
        add(3, weight.re * derivative.im);
        add(3, weight.im * derivative.re);
    }
    let amplitude = C64::new(sums[0] + correction[0], sums[1] + correction[1]);
    let angular_derivative = C64::new(sums[2] + correction[2], sums[3] + correction[3]);
    if !finite(amplitude) || !finite(angular_derivative) {
        return Err(invalid(
            "complex physical transfer or its derivative is not representable",
        ));
    }
    Ok(HarmonicTransferSample {
        amplitude,
        angular_derivative,
    })
}

/// Refine ONE isolated extremum of the supplied actual complex response.
/// Strict endpoint slopes are required: a zero-slope endpoint has only
/// one-sided evidence and refuses. A sampled lossless pole propagates its
/// evaluator's singularity error; an unresolved divergent slope refuses at
/// floating-point spacing. Neither becomes a finite response peak.
///
/// This executor-context entry uses the same owner as
/// [`refine_transfer_extremum_with_cancel`].
pub fn refine_transfer_extremum<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<TransferExtremum, ResonanceError<E>> {
    refine_transfer_extremum_with_cancel(
        bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
    )
}

/// Discover and refine the response peaks and notches resolved by a supplied
/// frequency grid, retaining the actual complex transfer and its tangent.
///
/// This executor-context entry shares the grid cache, global evaluation and
/// extremum budgets, and resolution limits of
/// [`survey_transfer_extrema_with_cancel`]. The grid is not a completeness
/// certificate; unresolved peak/notch pairs can remain between its samples.
pub fn survey_transfer_extrema<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<TransferSurvey, ResonanceError<E>> {
    survey_transfer_extrema_with_cancel(
        frequencies_rad_s,
        controls,
        maximum_extrema,
        || cx.checkpoint().is_err(),
        evaluate,
    )
}

/// Refine a separately bracketed crossing on each side of this response peak.
/// Use the SAME immutable physical transfer/source/receiver as for `peak`.
/// The lower bracket must rise through half `|H_peak|^2`; the upper must fall
/// through it. The caller isolates each crossing from neighbouring resonances.
/// Missing crossings, a DC/zero peak, a notch and a singular evaluator refuse.
///
/// This executor-context entry uses the same owner as
/// [`refine_half_power_bandwidth_with_cancel`].
pub fn refine_half_power_bandwidth<E>(
    peak: &TransferExtremum,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<HalfPowerBandwidth, ResonanceError<E>> {
    refine_half_power_bandwidth_with_cancel(
        peak,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
    )
}

/// Refine and differentiate the actual driven peak or cancellation notch.
///
/// This executor-context entry uses the same bounded first-derivative search
/// and one reserved final complete jet as its callback counterpart. Source,
/// receiver and parameter must describe the same immutable physical model.
pub fn refine_transfer_extremum_with_parameter_derivative<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
    evaluate_jet: impl FnMut(f64) -> Result<HarmonicTransferParameterJet, E>,
) -> Result<TransferExtremumParameterResponse, ResonanceError<E>> {
    refine_transfer_extremum_with_parameter_derivative_with_cancel(
        bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_jet,
    )
}

/// Differentiate both actual half-power crossings, width and observed Q.
///
/// Two final edge jets are reserved within this query's evaluation cap. The
/// previously resolved peak retains its own independent evaluation count.
pub fn refine_half_power_bandwidth_with_parameter_derivative<E>(
    peak: &TransferExtremumParameterResponse,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
    evaluate_jet: impl FnMut(f64) -> Result<HarmonicTransferParameterJet, E>,
) -> Result<HalfPowerBandwidthParameterResponse, ResonanceError<E>> {
    refine_half_power_bandwidth_with_parameter_derivative_with_cancel(
        peak,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_jet,
    )
}
