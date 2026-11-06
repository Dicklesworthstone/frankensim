//! Existing peak/linewidth refinement applied to one physical output sideband.

use super::{C64, PhsPortSidebandError, PreparedPhsPortSidebands, admit_forcing};
use core::cell::RefCell;
use fs_modal::harmonic::resonance::{
    HalfPowerBandwidth, HarmonicTransferSample, ResonanceControls, ResonanceError,
    TransferExtremum, refine_half_power_bandwidth_with_cancel,
    refine_transfer_extremum_with_cancel,
};
use fs_orbit::sideband::sensitivity::SidebandCarrierDerivativeControls;

/// One fixed-source, fixed-port harmonic-transfer query on the prepared orbit.
#[derive(Clone, Debug)]
pub struct PhsPortResonanceRequest<'a> {
    /// Complete original-port peak efforts at the preparation's signed input
    /// sideband indices. These complex phasors stay FIXED throughout the query;
    /// no source filter, amplitude normalization or frequency slope is inferred.
    pub forcing: &'a [C64],
    /// Original physical port to observe, without an additional receiver filter.
    pub output_port: usize,
    /// Fixed signed output index k in the preparation's retained output range.
    pub output_sideband: isize,
    /// Bracket of one isolated extremum in ACTUAL positive output frequency
    /// w=carrier+k*omega [rad/s]. The carrier itself may be negative.
    pub peak_bracket_rad_s: (f64, f64),
    /// Optional rising/falling brackets for half of the SAME peak's |H|^2.
    /// These also use actual positive output frequency, not carrier frequency.
    pub half_power_brackets_rad_s: Option<((f64, f64), (f64, f64))>,
    /// One TOTAL evaluation allowance shared by the peak and both flanks.
    pub controls: ResonanceControls,
    /// Existing additional per-evaluation carrier-derivative limits.
    pub derivative_controls: SidebandCarrierDerivativeControls,
    /// Original physical port output units per angular frequency, for every
    /// port. Each response/derivative pair retains its existing output checks.
    pub derivative_envelope_tolerances: &'a [f64],
}

/// Actual physical output-frequency peak and optional response linewidth.
#[derive(Clone, Debug)]
pub struct PhsPortResonanceResponse {
    /// The sample uses the modal owner's exp(-i*w*t) convention. Both complex
    /// amplitude and derivative were converted from the source's exp(+i*w*t).
    pub peak: TransferExtremum,
    /// Response quality factor is w_peak/linewidth for THIS original port and
    /// forcing. It is not a stored-energy, modal damping or Floquet decay Q.
    pub bandwidth: Option<HalfPowerBandwidth>,
    /// Peak evaluations plus both flank evaluations, bounded by one request cap.
    pub evaluations: usize,
}

type Result<T> = core::result::Result<T, ResonanceError<PhsPortSidebandError>>;

impl PreparedPhsPortSidebands<'_, '_> {
    /// Refine an original physical output resonance using the shared modal owner.
    ///
    /// Every evaluator call factors carrier=w-k*omega using this SAME prepared
    /// complete tangent, then obtains the physical response and analytic carrier
    /// derivative from that one factor. Since the base clock and k stay fixed,
    /// d(carrier)/dw=1. The selected original-port coefficient and its slope are
    /// both conjugated into the modal owner's exp(-i*w*t) convention.
    ///
    /// The entire requested frequency span must lie strictly inside one open
    /// interval (j*omega/2,(j+1)*omega/2), j>=0, with positive lower endpoint.
    /// This conservatively excludes commensurate opposite signed-line
    /// coincidences, where the real-signal receiver would combine different
    /// coefficients at an isolated frequency. No unrelated sidebands are summed,
    /// folded or used as a changing receiver. Sources and receivers with their
    /// own frequency dependence require their complete additional chain rule;
    /// this request explicitly fixes the probe phasors and one original port.
    ///
    /// The one evaluation cap includes the peak and BOTH half-power flanks.
    /// Each evaluation retains the existing core, added derivative and physical
    /// observation work/payload limits. Thus the campaign work is bounded by the
    /// evaluation cap times those per-evaluation allowances; tangent preparation
    /// precedes this call. No physical response arrays are retained between
    /// evaluations, and at most one carrier factor is live at once.
    ///
    /// This is local numerical evidence for the caller's fixed retained order
    /// and resolved base. Brackets must isolate continuous extrema/crossings.
    /// It adds no harmonic-tail/forward-error enclosure, global peak search,
    /// Floquet stability verdict, finite-amplitude prediction or damping fit.
    /// The caller establishes base and sideband-order adequacy across the whole
    /// span; no adaptive state-order refinement runs inside this peak query.
    ///
    /// # Errors
    /// Invalid physical receiver or frequency domain, original physical/numerical
    /// refusal, insufficient shared evaluations, invalid local sign brackets or
    /// cancellation. No partial peak/linewidth result is published on failure.
    pub fn refine_port_resonance(
        &mut self,
        request: &PhsPortResonanceRequest<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsPortResonanceResponse> {
        if cancelled() {
            return Err(ResonanceError::Cancelled);
        }
        let system = self.source.problem.original().system;
        let controls = &self.source.controls;
        let ports = system.port_dim();
        if request.output_port >= ports
            || request.output_sideband.unsigned_abs() > controls.output_harmonics
        {
            return Err(ResonanceError::Invalid(
                "one original physical port and retained output sideband required",
            ));
        }
        admit_forcing(system, request.forcing, controls).map_err(ResonanceError::Evaluation)?;
        if request.derivative_envelope_tolerances.len() != ports
            || request
                .derivative_envelope_tolerances
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(ResonanceError::Invalid(
                "positive output derivative tolerances for every original port required",
            ));
        }
        let omega = self.source.orbit.omega;
        admit_frequency_span(request, omega)?;
        let coefficient = (request.output_sideband + controls.output_harmonics as isize) as usize
            * ports
            + request.output_port;
        let maximum = request.controls.maximum_evaluations;
        // Both the shared query and long physical evaluation poll the SAME
        // callback. Each borrow ends at one poll; no nested borrow is retained.
        let gate = RefCell::new(cancelled);
        let mut evaluate = |frequency: f64| {
            let carrier = (-(request.output_sideband as f64)).mul_add(omega, frequency);
            let mut factor = self.factor(carrier, || (gate.borrow_mut())())?;
            let result = factor.solve_with_carrier_derivative(
                request.forcing,
                None,
                request.derivative_controls,
                request.derivative_envelope_tolerances,
                || (gate.borrow_mut())(),
            )?;
            Ok(HarmonicTransferSample {
                amplitude: result.response.port_coefficients[coefficient].conj(),
                angular_derivative: result.port_angular_derivative[coefficient].conj(),
            })
        };
        refine_query(
            request,
            maximum,
            &mut || (gate.borrow_mut())(),
            &mut evaluate,
        )
        .map_err(|error| campaign_error(error, maximum))
    }
}

// Keep orchestration and the shared peak/flank allowance identical for fixed
// and adaptively assessed physical evaluators. Search mathematics stays modal.
pub(super) fn refine_query<E>(
    request: &PhsPortResonanceRequest<'_>,
    maximum: usize,
    cancelled: &mut impl FnMut() -> bool,
    evaluate: &mut impl FnMut(f64) -> core::result::Result<HarmonicTransferSample, E>,
) -> core::result::Result<PhsPortResonanceResponse, ResonanceError<E>> {
    let mut controls = request.controls;
    controls.maximum_evaluations = maximum;
    let peak = refine_transfer_extremum_with_cancel(
        request.peak_bracket_rad_s,
        controls,
        &mut *cancelled,
        &mut *evaluate,
    )?;
    let bandwidth = match request.half_power_brackets_rad_s {
        None => None,
        Some((lower, upper)) => {
            let remaining = maximum - peak.evaluations;
            if remaining < 4 {
                return Err(ResonanceError::EvaluationBudget {
                    maximum_evaluations: maximum,
                });
            }
            controls.maximum_evaluations = remaining;
            Some(refine_half_power_bandwidth_with_cancel(
                &peak,
                lower,
                upper,
                controls,
                &mut *cancelled,
                &mut *evaluate,
            )?)
        }
    };
    let evaluations = peak.evaluations + bandwidth.as_ref().map_or(0, |b| b.evaluations);
    if cancelled() {
        return Err(ResonanceError::Cancelled);
    }
    Ok(PhsPortResonanceResponse {
        peak,
        bandwidth,
        evaluations,
    })
}

pub(super) fn campaign_error(
    error: ResonanceError<PhsPortSidebandError>,
    maximum: usize,
) -> ResonanceError<PhsPortSidebandError> {
    match error {
        ResonanceError::Evaluation(PhsPortSidebandError::Cancelled) => ResonanceError::Cancelled,
        ResonanceError::EvaluationBudget { .. } => ResonanceError::EvaluationBudget {
            maximum_evaluations: maximum,
        },
        other => other,
    }
}

pub(super) fn admit_frequency_span(
    request: &PhsPortResonanceRequest<'_>,
    omega: f64,
) -> Result<()> {
    let valid = |(lower, upper): (f64, f64)| {
        lower.is_finite() && upper.is_finite() && lower > 0.0 && upper > lower
    };
    if !valid(request.peak_bracket_rad_s)
        || request
            .half_power_brackets_rad_s
            .is_some_and(|(lower, upper)| !valid(lower) || !valid(upper))
    {
        return Err(ResonanceError::Invalid(
            "finite increasing strictly positive physical frequency brackets required",
        ));
    }
    let (mut lower, mut upper) = request.peak_bracket_rad_s;
    if let Some((lo, hi)) = request.half_power_brackets_rad_s {
        lower = lower.min(lo.0).min(hi.0);
        upper = upper.max(lo.1).max(hi.1);
    }
    let spacing = 0.5 * omega;
    let cell = (lower / spacing).floor();
    let cell_lower = cell * spacing;
    let cell_upper = (cell + 1.0) * spacing;
    if !spacing.is_finite()
        || spacing <= 0.0
        || !cell.is_finite()
        || cell + 1.0 == cell
        || !cell_lower.is_finite()
        || !cell_upper.is_finite()
        || lower <= cell_lower
        || upper >= cell_upper
    {
        return Err(ResonanceError::Invalid(
            "entire physical frequency span must lie strictly between consecutive omega/2 boundaries without signed-line coincidences",
        ));
    }
    Ok(())
}
