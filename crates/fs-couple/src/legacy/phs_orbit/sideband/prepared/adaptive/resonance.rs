//! Local physical peak searches with a common resolved order at every query.

use super::super::super::{ObservationBudget, admit_storage_with_retained, input_work};
use super::super::resonance::{admit_frequency_span, campaign_error, refine_query};
use super::super::sensitivity::{
    admit_port_tolerances, observe_pair, ode_carrier_derivative, trace_payload,
};
use super::super::{
    C64, PhsPortResonanceRequest, PhsPortResonanceResponse,
    PhsPortSidebandCarrierDerivativeResponse, PhsPortSidebandError, PhsPortSidebandSource,
    PreparedPeriodicLinearization, SidebandScales, admit_forcing, input_forcing,
};
use super::{admit_observations, observer_controls};
use core::cell::RefCell;
use fs_math::{dd::Dd, det};
use fs_modal::harmonic::resonance::{HarmonicTransferSample, ResonanceError};
use fs_orbit::sideband::{
    SidebandError,
    sensitivity::refinement::{SidebandCarrierRefinementCheck, SidebandCarrierRefinementControls},
};

/// Accuracy goals for a local physical peak and optional half-power campaign.
#[derive(Clone, Debug)]
pub struct PhsPortAdaptiveResonanceRequest<'a> {
    /// Original fixed probe/receiver, physical-frequency brackets, derivative
    /// limits and one TOTAL frequency-evaluation cap for all attempted rounds.
    pub resonance: PhsPortResonanceRequest<'a>,
    pub refinement: &'a SidebandCarrierRefinementControls,
    /// Actual low/high output envelope changes, for EVERY original port, in
    /// original port units. All retained signed output coefficients count once.
    pub port_change_tolerances: &'a [f64],
    /// Corresponding changes in port units per angular frequency.
    pub derivative_port_change_tolerances: &'a [f64],
    /// Independent inter-order change allowance for the peak, both half-power
    /// crossings and linewidth. It includes the full numerical root-bracket
    /// widths, so coarse root estimates cannot imply finer location agreement.
    pub absolute_angular_change_tolerance_rad_s: f64,
    pub relative_angular_change_tolerance: f64,
}

/// Worst measured ratios across the frequencies actually queried in one round.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhsPortResonanceRefinementCheck {
    pub lower_harmonics: usize,
    pub harmonics: usize,
    pub evaluations: usize,
    /// Maximum across full equation, grid and finite omitted-band checks.
    pub maximum_equation_ratio: f64,
    pub maximum_derivative_equation_ratio: f64,
    pub maximum_response_change_ratio: f64,
    pub maximum_derivative_change_ratio: f64,
    pub maximum_port_change_ratio: f64,
    pub maximum_port_derivative_change_ratio: f64,
    /// Present only after this and an earlier order completed the SAME query.
    pub peak_change_ratio: Option<f64>,
    /// Maximum of lower crossing, upper crossing and linewidth changes.
    pub bandwidth_change_ratio: Option<f64>,
    pub completed: bool,
}

/// A complete local query with physical pair and actual frequency-location
/// convergence evidence, without an infinite-tail or global-peak enclosure.
#[derive(Clone, Debug)]
pub struct AdaptivePhsPortResonanceResponse {
    /// The accepted final round. Its evaluation count is that round alone.
    pub response: PhsPortResonanceResponse,
    pub state_harmonics: usize,
    pub refinement_trace: Vec<PhsPortResonanceRefinementCheck>,
    /// ALL attempted frequency queries, including rejected rounds and flanks.
    pub evaluations: usize,
    /// Actual low/high pair-solve attempts; at most twice `evaluations`.
    pub pair_solves: usize,
    /// TOTAL original OUTPUT observation Hessian actions over the campaign.
    /// Generic preparation/assessment tangent callbacks are accounted separately.
    pub total_hessian_actions: usize,
}

type Result<T> = core::result::Result<T, ResonanceError<PhsPortSidebandError>>;

impl PhsPortSidebandSource<'_> {
    /// Refine one peak/linewidth while resolving the complete physical pair at
    /// EVERY queried frequency and preserving a common order within each search.
    ///
    /// A round prepares ONE high-order original tangent and uses its common
    /// verified fine quadrature to actually factor/solve both the lower and
    /// higher retained systems at identical carriers. Only one LU is live at a
    /// time. The lower result is not a clipped/padded copy of the higher result.
    /// Original response AND derivative equations are sampled independently;
    /// all state changes and all original G^T Hxx output coefficient changes
    /// have their own physical tolerances. A failure of these resolution checks
    /// restarts the entire query at the next common order. No root-finder step
    /// sees a transfer assembled from frequency-dependent accepted orders.
    ///
    /// The peak and optional BOTH half-power flanks must also agree between two
    /// completed, distinct-order searches. The maximum displacement between
    /// their returned numerical brackets must fit the separate location goals;
    /// equal coarse midpoint estimates cannot hide unresolved root widths.
    /// The first high order is min(2*initial+1,max), so
    /// a constant transfer family starting at one normally completes searches
    /// at three and seven. An intervening unresolved round does not erase the
    /// earlier completed comparison. Maximum order without all checks refuses.
    ///
    /// The original probe support, base, clock, phase and selected port/sideband
    /// stay fixed. Observed positive frequency is w=carrier+k*omega; the existing
    /// signed-line-coincidence exclusion and modal convention conversion apply.
    /// The authored physical output bandwidth remains independent and never
    /// grows silently. Its sampling/truncation refusals remain physical errors.
    ///
    /// One evaluation means one requested frequency with up to TWO actual pair
    /// solves. The request cap includes all restarted rounds and both flanks.
    /// Source OUTPUT observation actions and physical arithmetic are cumulative
    /// across the ENTIRE call. Each generic call retains its core and additional
    /// derivative allowance: total numerical work is bounded conservatively by
    /// (2*evaluations+8)*core.maximum_work plus 2*evaluations*derivative.maximum_work
    /// plus the adapter work cap, with representability checked before physics.
    /// These exclude source-callback internal work as in the underlying owners.
    /// Generic payload caps add to the adapter cap. One tangent, one live LU,
    /// two transient physical pairs, fixed mapped input and bounded history are
    /// admitted; no unbounded frequency-to-response cache is retained.
    ///
    /// This measures finite sampled local convergence; it does not certify an
    /// infinite harmonic tail, a global maximum, Floquet stability or a damping Q.
    ///
    /// # Errors
    /// Original physical/numerical refusal, invalid fixed receiver/domain or
    /// local sign brackets, cumulative budget exhaustion, unresolved harmonic
    /// or frequency-location changes, or cancellation. No partial peak escapes.
    pub fn refine_port_resonance_adaptive(
        &mut self,
        request: &PhsPortAdaptiveResonanceRequest<'_>,
        scales: &SidebandScales,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptivePhsPortResonanceResponse> {
        if cancelled() {
            return Err(ResonanceError::Cancelled);
        }
        admit_request(self, request)?;
        let query = &request.resonance;
        let system = self.problem.original().system;
        let (n, m) = (system.state_dim(), system.port_dim());
        let initial = self.controls.sideband.harmonics;
        let maximum = query.controls.maximum_evaluations;
        let omega = self.orbit.omega;
        let coefficient =
            (query.output_sideband + self.controls.output_harmonics as isize) as usize * m
                + query.output_port;
        let static_retained = 8 * m as u128
            + 8 * core::mem::size_of::<PhsPortResonanceRefinementCheck>() as u128
            + core::mem::size_of::<PhsPortResonanceResponse>() as u128;
        // Reserve the fixed mapped source and owned round controls/history before
        // the first allocation; the original borrowed forcing is caller storage.
        admit_storage_with_retained(
            n,
            m,
            &self.controls,
            self.controls.initial_samples,
            static_retained + 16 * n as u128 * (2 * initial + 1) as u128,
        )
        .map_err(ResonanceError::Evaluation)?;
        let rhs = input_forcing(system, query.forcing, initial, &mut cancelled)
            .map_err(ResonanceError::Evaluation)?;
        let mut budget = ObservationBudget {
            actions: 0,
            work: input_work(n, m, initial),
            retained_bytes: static_retained,
        };
        let mut trace = Vec::with_capacity(8);
        let mut previous_completed: Option<PhsPortResonanceResponse> = None;
        let mut evaluations = 0_usize;
        let mut pair_solves = 0_usize;
        let mut lower = initial;
        let mut high = (2 * lower + 1).min(request.refinement.maximum_harmonics);
        let gate = RefCell::new(cancelled);
        loop {
            if (gate.borrow_mut())() {
                return Err(ResonanceError::Cancelled);
            }
            if maximum - evaluations < 3 {
                return Err(ResonanceError::EvaluationBudget {
                    maximum_evaluations: maximum,
                });
            }
            let count = n * (2 * high + 1);
            let retained = static_retained + 16 * rhs.capacity() as u128 + 32 * count as u128;
            budget.retained_bytes = retained;
            let mut observer =
                observer_controls(self, high, retained).map_err(ResonanceError::Evaluation)?;
            observer.sideband.samples_per_period = observer
                .sideband
                .samples_per_period
                .max(2 * (3 * high + 2).max(self.orbit.coeffs[0].len() - 1) + 1);
            let padding_work = 64 * (count + rhs.len()) as u128;
            admit_observations(n, m, &observer, &budget, 4, padding_work)
                .map_err(ResonanceError::Evaluation)?;
            budget.work += padding_work;
            let mut padded = vec![C64::ZERO; count];
            let offset = (high - initial) * n;
            padded[offset..offset + rhs.len()].copy_from_slice(&rhs);
            let slope = vec![C64::ZERO; count];
            let low_offset = (high - lower) * n;
            let low_end = low_offset + n * (2 * lower + 1);
            let prepared = PreparedPeriodicLinearization::new(
                &self.problem,
                self.orbit,
                observer.sideband,
                scales,
                || (gate.borrow_mut())(),
            )
            .map_err(PhsPortSidebandError::from)
            .map_err(ResonanceError::Evaluation)?;
            let mut check = PhsPortResonanceRefinementCheck {
                lower_harmonics: lower,
                harmonics: high,
                evaluations: 0,
                maximum_equation_ratio: 0.0,
                maximum_derivative_equation_ratio: 0.0,
                maximum_response_change_ratio: 0.0,
                maximum_derivative_change_ratio: 0.0,
                maximum_port_change_ratio: 0.0,
                maximum_port_derivative_change_ratio: 0.0,
                peak_change_ratio: None,
                bandwidth_change_ratio: None,
                completed: false,
            };
            let remaining = maximum - evaluations;
            let round = {
                let mut evaluate = |frequency: f64| {
                    evaluations += 1;
                    check.evaluations += 1;
                    let derivative_work =
                        16 * n as u128 * n as u128 * (2 * lower + 2 * high + 2) as u128;
                    let comparison_work =
                        256 * m as u128 * (2 * observer.output_harmonics + 1) as u128;
                    admit_observations(
                        n,
                        m,
                        &observer,
                        &budget,
                        4,
                        derivative_work + comparison_work,
                    )?;
                    budget.work += derivative_work + comparison_work;
                    let carrier = (-(query.output_sideband as f64)).mul_add(omega, frequency);
                    pair_solves += 1;
                    let low = prepared.solve_carrier_derivative_assessed_at_order(
                        lower,
                        carrier,
                        &padded[low_offset..low_end],
                        &slope[low_offset..low_end],
                        query.derivative_controls,
                        &request.refinement.tolerances,
                        |_, out| ode_carrier_derivative(n, out),
                        || (gate.borrow_mut())(),
                    )?;
                    pair_solves += 1;
                    let mut current = prepared.solve_carrier_derivative_assessed_at_order(
                        high,
                        carrier,
                        &padded,
                        &slope,
                        query.derivative_controls,
                        &request.refinement.tolerances,
                        |_, out| ode_carrier_derivative(n, out),
                        || (gate.borrow_mut())(),
                    )?;
                    current.compare_with(&low, &request.refinement.tolerances)?;
                    accumulate(&mut check, current.check);
                    if !current.check.is_resolved() {
                        return Err(Evaluation::Refine);
                    }
                    // High state arrays remain live while the original low pair
                    // is observed; then low's ENTIRE physical pair remains live
                    // during high observation. Admit both lifetimes explicitly.
                    budget.retained_bytes = retained
                        + 16 * (current.response.response.coefficients.capacity() as u128
                            + current.response.angular_derivative.coefficients.capacity() as u128);
                    observer.sideband.harmonics = lower;
                    let low = observe_pair(
                        self,
                        carrier,
                        low.response,
                        &observer,
                        query.derivative_envelope_tolerances,
                        low.base_residual,
                        low.tangent_sampling_defect,
                        &mut budget,
                        &mut || (gate.borrow_mut())(),
                    )?;
                    budget.retained_bytes = retained + pair_payload(&low);
                    observer.sideband.harmonics = high;
                    let current = observe_pair(
                        self,
                        carrier,
                        current.response,
                        &observer,
                        query.derivative_envelope_tolerances,
                        current.base_residual,
                        current.tangent_sampling_defect,
                        &mut budget,
                        &mut || (gate.borrow_mut())(),
                    )?;
                    let port_change = output_change(
                        &low.response.port_coefficients,
                        &current.response.port_coefficients,
                        request.port_change_tolerances,
                    );
                    let derivative_change = output_change(
                        &low.port_angular_derivative,
                        &current.port_angular_derivative,
                        request.derivative_port_change_tolerances,
                    );
                    check.maximum_port_change_ratio =
                        check.maximum_port_change_ratio.max(port_change);
                    check.maximum_port_derivative_change_ratio = check
                        .maximum_port_derivative_change_ratio
                        .max(derivative_change);
                    budget.retained_bytes = retained;
                    if port_change > 1.0 || derivative_change > 1.0 {
                        return Err(Evaluation::Refine);
                    }
                    Ok(HarmonicTransferSample {
                        amplitude: current.response.port_coefficients[coefficient].conj(),
                        angular_derivative: current.port_angular_derivative[coefficient].conj(),
                    })
                };
                refine_query(
                    query,
                    remaining,
                    &mut || (gate.borrow_mut())(),
                    &mut evaluate,
                )
            };
            // Keep no cached tangent, large forcing arrays or response pairs
            // across rounds. The earlier COMPLETED scalar query alone survives.
            drop(prepared);
            let complete = match round {
                Ok(current) => Some(current),
                Err(ResonanceError::Evaluation(Evaluation::Refine)) => None,
                Err(error) => return Err(convert_error(error, maximum)),
            };
            if let Some(current) = complete {
                check.completed = true;
                if let Some(previous) = &previous_completed {
                    check.peak_change_ratio = Some(if current.peak.kind == previous.peak.kind {
                        frequency_change(
                            current.peak.bracket_rad_s,
                            previous.peak.bracket_rad_s,
                            request,
                        )
                    } else {
                        f64::INFINITY
                    });
                    if let (Some(a), Some(b)) = (&current.bandwidth, &previous.bandwidth) {
                        check.bandwidth_change_ratio = Some(
                            frequency_change(a.lower_bracket_rad_s, b.lower_bracket_rad_s, request)
                                .max(frequency_change(
                                    a.upper_bracket_rad_s,
                                    b.upper_bracket_rad_s,
                                    request,
                                ))
                                .max(frequency_change(
                                    (
                                        a.upper_bracket_rad_s.0 - a.lower_bracket_rad_s.1,
                                        a.upper_bracket_rad_s.1 - a.lower_bracket_rad_s.0,
                                    ),
                                    (
                                        b.upper_bracket_rad_s.0 - b.lower_bracket_rad_s.1,
                                        b.upper_bracket_rad_s.1 - b.lower_bracket_rad_s.0,
                                    ),
                                    request,
                                )),
                        );
                    }
                }
                let resolved = check.peak_change_ratio.is_some_and(|r| r <= 1.0)
                    && (query.half_power_brackets_rad_s.is_none()
                        || check.bandwidth_change_ratio.is_some_and(|r| r <= 1.0));
                trace.push(check);
                if resolved {
                    if (gate.borrow_mut())() {
                        return Err(ResonanceError::Cancelled);
                    }
                    return Ok(AdaptivePhsPortResonanceResponse {
                        response: current,
                        state_harmonics: high,
                        refinement_trace: trace,
                        evaluations,
                        pair_solves,
                        total_hessian_actions: budget.actions,
                    });
                }
                previous_completed = Some(current);
            } else {
                trace.push(check);
            }
            if high == request.refinement.maximum_harmonics {
                return Err(ResonanceError::Evaluation(
                    PhsPortSidebandError::ResonanceUnresolved { last: check },
                ));
            }
            lower = high;
            high = (2 * high + 1).min(request.refinement.maximum_harmonics);
        }
    }
}

#[derive(Debug)]
enum Evaluation {
    Physical(PhsPortSidebandError),
    Refine,
}
impl From<PhsPortSidebandError> for Evaluation {
    fn from(error: PhsPortSidebandError) -> Self {
        Self::Physical(error)
    }
}
impl From<SidebandError> for Evaluation {
    fn from(error: SidebandError) -> Self {
        Self::Physical(error.into())
    }
}

fn convert_error(
    error: ResonanceError<Evaluation>,
    maximum: usize,
) -> ResonanceError<PhsPortSidebandError> {
    let mapped = match error {
        ResonanceError::Invalid(message) => ResonanceError::Invalid(message),
        ResonanceError::InvalidBracket(message) => ResonanceError::InvalidBracket(message),
        ResonanceError::UnresolvedStationaryInterval => {
            ResonanceError::UnresolvedStationaryInterval
        }
        ResonanceError::Evaluation(Evaluation::Physical(error)) => {
            ResonanceError::Evaluation(error)
        }
        // Resolution restarts are handled by the caller before conversion.
        ResonanceError::Evaluation(Evaluation::Refine) => {
            ResonanceError::Invalid("unhandled common-order refinement")
        }
        ResonanceError::EvaluationBudget { .. } => ResonanceError::EvaluationBudget {
            maximum_evaluations: maximum,
        },
        ResonanceError::ExtremumBudget { maximum_extrema } => {
            ResonanceError::ExtremumBudget { maximum_extrema }
        }
        ResonanceError::Cancelled => ResonanceError::Cancelled,
    };
    campaign_error(mapped, maximum)
}

fn admit_request(
    source: &PhsPortSidebandSource<'_>,
    request: &PhsPortAdaptiveResonanceRequest<'_>,
) -> Result<()> {
    let query = &request.resonance;
    let system = source.problem.original().system;
    let (n, m) = (system.state_dim(), system.port_dim());
    if query.output_port >= m
        || query.output_sideband.unsigned_abs() > source.controls.output_harmonics
    {
        return Err(ResonanceError::Invalid(
            "one original physical port and retained output sideband required",
        ));
    }
    admit_forcing(system, query.forcing, &source.controls).map_err(ResonanceError::Evaluation)?;
    for tolerances in [
        query.derivative_envelope_tolerances,
        request.port_change_tolerances,
        request.derivative_port_change_tolerances,
    ] {
        admit_port_tolerances(m, tolerances).map_err(ResonanceError::Evaluation)?;
    }
    let tolerances = &request.refinement.tolerances;
    if request.refinement.maximum_harmonics <= source.controls.sideband.harmonics
        || request.refinement.maximum_harmonics > 64
        || [
            &tolerances.equation_tolerances,
            &tolerances.state_tolerances,
            &tolerances.derivative_equation_tolerances,
            &tolerances.derivative_state_tolerances,
        ]
        .iter()
        .any(|v| v.len() != n || v.iter().any(|x| !x.is_finite() || *x <= 0.0))
    {
        return Err(ResonanceError::Invalid(
            "increasing bounded orders and independent full physical pair tolerances required",
        ));
    }
    let valid_tolerances = |absolute: f64, relative: f64| {
        absolute.is_finite()
            && absolute >= 0.0
            && relative.is_finite()
            && (0.0..1.0).contains(&relative)
            && (absolute > 0.0 || relative > 0.0)
    };
    if !valid_tolerances(
        request.absolute_angular_change_tolerance_rad_s,
        request.relative_angular_change_tolerance,
    ) || !valid_tolerances(
        query.controls.absolute_angular_tolerance_rad_s,
        query.controls.relative_angular_tolerance,
    ) || query.controls.maximum_evaluations < 3
        || query.controls.maximum_evaluations > usize::MAX / 2
    {
        return Err(ResonanceError::Invalid(
            "independent finite root/location tolerances and bounded shared evaluation count required",
        ));
    }
    let extra = query.derivative_controls;
    if !extra.maximum_relative_residual.is_finite()
        || extra.maximum_relative_residual <= 0.0
        || extra.maximum_relative_residual >= 1.0
        || !extra.maximum_omitted_residual.is_finite()
        || extra.maximum_omitted_residual <= 0.0
    {
        return Err(ResonanceError::Invalid(
            "positive derivative limits with backward error below one required",
        ));
    }
    let calls = 2 * query.controls.maximum_evaluations as u128;
    if (calls + 8)
        .checked_mul(source.controls.sideband.maximum_work as u128)
        .and_then(|core| {
            calls
                .checked_mul(extra.maximum_work as u128)
                .and_then(|added| core.checked_add(added))
        })
        .and_then(|work| work.checked_add(source.controls.max_observation_work as u128))
        .is_none()
    {
        return Err(ResonanceError::Evaluation(PhsPortSidebandError::Budget(
            "total bounded resonance arithmetic is not representable",
        )));
    }
    admit_frequency_span(query, source.orbit.omega)
}

fn accumulate(
    check: &mut PhsPortResonanceRefinementCheck,
    measured: SidebandCarrierRefinementCheck,
) {
    let ratio = |c: &fs_orbit::sideband::refinement::SidebandRefinementCheck| {
        c.equation_residual_ratio
            .max(c.sampling_defect_ratio)
            .max(c.omitted_residual_ratio)
    };
    check.maximum_equation_ratio = check.maximum_equation_ratio.max(ratio(&measured.response));
    check.maximum_derivative_equation_ratio = check
        .maximum_derivative_equation_ratio
        .max(ratio(&measured.angular_derivative));
    check.maximum_response_change_ratio = check.maximum_response_change_ratio.max(
        measured
            .response
            .state_change_ratio
            .unwrap_or(f64::INFINITY),
    );
    check.maximum_derivative_change_ratio = check.maximum_derivative_change_ratio.max(
        measured
            .angular_derivative
            .state_change_ratio
            .unwrap_or(f64::INFINITY),
    );
}

fn output_change(previous: &[C64], current: &[C64], tolerances: &[f64]) -> f64 {
    let ports = tolerances.len();
    let mut maximum = 0.0_f64;
    for (port, tolerance) in tolerances.iter().enumerate() {
        let mut norm = 0.0;
        for (a, b) in previous
            .iter()
            .skip(port)
            .step_by(ports)
            .zip(current.iter().skip(port).step_by(ports))
        {
            let re = (Dd::from_f64(a.re) - Dd::from_f64(b.re)).to_f64() / tolerance;
            let im = (Dd::from_f64(a.im) - Dd::from_f64(b.im)).to_f64() / tolerance;
            norm = det::hypot(det::hypot(norm, re), im);
        }
        if !norm.is_finite() {
            return f64::INFINITY;
        }
        maximum = maximum.max(norm);
    }
    maximum
}

fn frequency_change(
    a: (f64, f64),
    b: (f64, f64),
    request: &PhsPortAdaptiveResonanceRequest<'_>,
) -> f64 {
    let allowance = request.absolute_angular_change_tolerance_rad_s.max(
        request.relative_angular_change_tolerance
            * a.0.abs().max(a.1.abs()).max(b.0.abs()).max(b.1.abs()),
    );
    let lower = (Dd::from_f64(a.0) - Dd::from_f64(b.1)).to_f64() / allowance;
    let upper = (Dd::from_f64(a.1) - Dd::from_f64(b.0)).to_f64() / allowance;
    if !lower.is_finite() || !upper.is_finite() {
        return f64::INFINITY;
    }
    lower.abs().max(upper.abs())
}

fn pair_payload(pair: &PhsPortSidebandCarrierDerivativeResponse) -> u128 {
    16 * (pair.response.state_response.coefficients.capacity() as u128
        + pair.state_angular_derivative.coefficients.capacity() as u128
        + pair.response.port_coefficients.capacity() as u128
        + pair.port_angular_derivative.capacity() as u128)
        + trace_payload(
            &pair.response.observation_trace,
            pair.response.observation_trace.capacity(),
        )
        + trace_payload(
            &pair.derivative_observation_trace,
            pair.derivative_observation_trace.capacity(),
        )
}
