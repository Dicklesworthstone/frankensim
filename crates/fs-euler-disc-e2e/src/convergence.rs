//! Typed outcome, refinement, and calibration-readiness analysis for Euler runs.
//!
//! This module deliberately analyses already-executed [`CoupledRun`] values; it
//! neither extends a run nor fits physical parameters.  In particular, a time
//! horizon is a right-censoring boundary, not a terminal spin time, and a
//! small energy defect is not evidence of trajectory convergence or physical
//! validity.

use core::fmt;

use crate::coupled_runner::{CoupledNumericalRefusalReason, CoupledRun, CoupledTerminal};

/// Relative tolerance used only to bind a retained floating checkpoint time
/// back to the caller's declared horizon. It is deliberately far tighter than
/// any scientific timing acceptance band, while allowing accumulated binary
/// timestep roundoff in a long fixed-step run.
pub const HORIZON_MATCH_RELATIVE_TOLERANCE: f64 = 1.0e-9;

/// A physically defined terminal condition observed by the current runner.
///
/// Future eventful mechanics may add distinct physical terminal kinds without
/// changing the meaning of right-censoring or numerical refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhysicalTerminal {
    /// The declared inclination event was reached.
    InclinationThreshold,
}

/// The mutually exclusive disposition of one run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RunOutcome {
    /// A physical terminal event was observed at the stated simulation time.
    PhysicalTerminal {
        kind: PhysicalTerminal,
        event_time_s: f64,
    },
    /// The run was still live when its declared horizon elapsed.
    ///
    /// This is a lower bound on a future terminal time, not an observed
    /// duration and not a failed run.
    RightCensored { censor_time_s: f64 },
    /// The numerical lane refused before a physical terminal event.
    NumericalRefusal {
        last_valid_time_s: f64,
        reason: CoupledNumericalRefusalReason,
    },
}

impl RunOutcome {
    /// Coarse outcome class suitable for determining whether refinement rows
    /// are even comparable.
    #[must_use]
    pub const fn class(self) -> OutcomeClass {
        match self {
            Self::PhysicalTerminal { .. } => OutcomeClass::PhysicalTerminal,
            Self::RightCensored { .. } => OutcomeClass::RightCensored,
            Self::NumericalRefusal { .. } => OutcomeClass::NumericalRefusal,
        }
    }

    /// Returns an observed physical terminal time, never a censoring time.
    #[must_use]
    pub const fn observed_terminal_time_s(self) -> Option<f64> {
        match self {
            Self::PhysicalTerminal { event_time_s, .. } => Some(event_time_s),
            Self::RightCensored { .. } | Self::NumericalRefusal { .. } => None,
        }
    }
}

/// Coarse outcome class used only for comparability diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutcomeClass {
    PhysicalTerminal,
    RightCensored,
    NumericalRefusal,
}

/// Refusal emitted when a retained run cannot be classified safely.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutcomeError {
    /// The runner checkpoint did not retain a finite non-negative time.
    InvalidCheckpointTime,
}

impl fmt::Display for OutcomeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for OutcomeError {}

/// Classifies a coupled run without treating horizon completion as a terminal
/// event.
pub fn classify_outcome(run: &CoupledRun) -> Result<RunOutcome, OutcomeError> {
    let time_s = run.checkpoint.time_s;
    if !time_s.is_finite() || time_s < 0.0 {
        return Err(OutcomeError::InvalidCheckpointTime);
    }
    Ok(match run.terminal {
        CoupledTerminal::TerminalInclination => RunOutcome::PhysicalTerminal {
            kind: PhysicalTerminal::InclinationThreshold,
            event_time_s: time_s,
        },
        CoupledTerminal::HorizonReached => RunOutcome::RightCensored {
            censor_time_s: time_s,
        },
        CoupledTerminal::NumericalRefusal { reason } => RunOutcome::NumericalRefusal {
            last_valid_time_s: time_s,
            reason,
        },
    })
}

/// A duration ordering that is permitted only for two observed, like-kind
/// physical terminal events.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObservedDurationOrdering {
    LeftShorter,
    Equal,
    LeftLonger,
}

/// Why an apparent duration ranking is not scientifically comparable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RankingRefusal {
    LeftNotObservedPhysicalTerminal { class: OutcomeClass },
    RightNotObservedPhysicalTerminal { class: OutcomeClass },
    InvalidLeftObservedTerminalTime,
    InvalidRightObservedTerminalTime,
    DifferentPhysicalTerminalKinds,
}

/// Conservative result of comparing terminal-time information with
/// right-censoring retained explicitly.
///
/// `Indeterminate` is a successful, information-preserving outcome: it means
/// the available event and lower-bound information cannot establish a strict
/// ordering. It must not be converted to a tie or a ranking.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CensorAwareDurationOrdering {
    ProvenLeftShorter,
    EqualObserved,
    ProvenLeftLonger,
    Indeterminate,
}

/// A comparison refusal distinct from an indeterminate censoring result.
/// Numerical failures and malformed manually constructed outcomes do not carry
/// admissible lower-bound information.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CensorAwareRankingRefusal {
    LeftNumericalRefusal,
    RightNumericalRefusal,
    InvalidLeftPhysicalTerminalTime,
    InvalidRightPhysicalTerminalTime,
    InvalidLeftCensorTime,
    InvalidRightCensorTime,
    DifferentPhysicalTerminalKinds,
}

/// Compare two observed physical terminal times.
///
/// The function intentionally refuses a comparison involving right-censoring
/// or a numerical refusal.  A later preregistered censor-aware score is a
/// separate analysis layer; it must not be improvised here.
pub fn compare_observed_durations(
    left: RunOutcome,
    right: RunOutcome,
) -> Result<ObservedDurationOrdering, RankingRefusal> {
    let RunOutcome::PhysicalTerminal {
        kind: left_kind,
        event_time_s: left_time,
    } = left
    else {
        return Err(RankingRefusal::LeftNotObservedPhysicalTerminal {
            class: left.class(),
        });
    };
    let RunOutcome::PhysicalTerminal {
        kind: right_kind,
        event_time_s: right_time,
    } = right
    else {
        return Err(RankingRefusal::RightNotObservedPhysicalTerminal {
            class: right.class(),
        });
    };
    if !left_time.is_finite() || left_time < 0.0 {
        return Err(RankingRefusal::InvalidLeftObservedTerminalTime);
    }
    if !right_time.is_finite() || right_time < 0.0 {
        return Err(RankingRefusal::InvalidRightObservedTerminalTime);
    }
    if left_kind != right_kind {
        return Err(RankingRefusal::DifferentPhysicalTerminalKinds);
    }
    Ok(if left_time < right_time {
        ObservedDurationOrdering::LeftShorter
    } else if left_time > right_time {
        ObservedDurationOrdering::LeftLonger
    } else {
        ObservedDurationOrdering::Equal
    })
}

/// Compare two terminal-time observations while respecting right-censoring.
///
/// An observed event at `t` proves that it is shorter than a different run
/// only when `t` is strictly below that run's observed event time or retained
/// censoring lower bound. Two censored runs, or an observed event overlapping
/// a censoring boundary, remain indeterminate. A numerical refusal is not a
/// censoring record and therefore refuses the comparison.
pub fn compare_censor_aware_durations(
    left: RunOutcome,
    right: RunOutcome,
) -> Result<CensorAwareDurationOrdering, CensorAwareRankingRefusal> {
    validate_duration_outcome(left, true)?;
    validate_duration_outcome(right, false)?;
    if let (
        RunOutcome::PhysicalTerminal {
            kind: left_kind, ..
        },
        RunOutcome::PhysicalTerminal {
            kind: right_kind, ..
        },
    ) = (left, right)
    {
        if left_kind != right_kind {
            return Err(CensorAwareRankingRefusal::DifferentPhysicalTerminalKinds);
        }
    }
    match (left, right) {
        (
            RunOutcome::PhysicalTerminal {
                event_time_s: left_time,
                ..
            },
            RunOutcome::PhysicalTerminal {
                event_time_s: right_time,
                ..
            },
        ) => Ok(if left_time < right_time {
            CensorAwareDurationOrdering::ProvenLeftShorter
        } else if left_time > right_time {
            CensorAwareDurationOrdering::ProvenLeftLonger
        } else {
            CensorAwareDurationOrdering::EqualObserved
        }),
        (
            RunOutcome::PhysicalTerminal {
                event_time_s: left_time,
                ..
            },
            RunOutcome::RightCensored {
                censor_time_s: right_lower_bound,
            },
        ) if left_time < right_lower_bound => Ok(CensorAwareDurationOrdering::ProvenLeftShorter),
        (
            RunOutcome::RightCensored {
                censor_time_s: left_lower_bound,
            },
            RunOutcome::PhysicalTerminal {
                event_time_s: right_time,
                ..
            },
        ) if left_lower_bound > right_time => Ok(CensorAwareDurationOrdering::ProvenLeftLonger),
        (
            RunOutcome::PhysicalTerminal { .. } | RunOutcome::RightCensored { .. },
            RunOutcome::PhysicalTerminal { .. } | RunOutcome::RightCensored { .. },
        ) => Ok(CensorAwareDurationOrdering::Indeterminate),
        (RunOutcome::NumericalRefusal { .. }, _) | (_, RunOutcome::NumericalRefusal { .. }) => {
            unreachable!("numerical refusals are rejected before comparison")
        }
    }
}

fn validate_duration_outcome(
    outcome: RunOutcome,
    left: bool,
) -> Result<(), CensorAwareRankingRefusal> {
    match outcome {
        RunOutcome::NumericalRefusal { .. } => Err(if left {
            CensorAwareRankingRefusal::LeftNumericalRefusal
        } else {
            CensorAwareRankingRefusal::RightNumericalRefusal
        }),
        RunOutcome::PhysicalTerminal { event_time_s, .. }
            if !event_time_s.is_finite() || event_time_s < 0.0 =>
        {
            Err(if left {
                CensorAwareRankingRefusal::InvalidLeftPhysicalTerminalTime
            } else {
                CensorAwareRankingRefusal::InvalidRightPhysicalTerminalTime
            })
        }
        RunOutcome::RightCensored { censor_time_s }
            if !censor_time_s.is_finite() || censor_time_s < 0.0 =>
        {
            Err(if left {
                CensorAwareRankingRefusal::InvalidLeftCensorTime
            } else {
                CensorAwareRankingRefusal::InvalidRightCensorTime
            })
        }
        RunOutcome::PhysicalTerminal { .. } | RunOutcome::RightCensored { .. } => Ok(()),
    }
}

/// Declared rule for extending a live trajectory without converting a censor
/// boundary into an unbounded computation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HorizonContinuationPolicy {
    /// First declared horizon in seconds.
    pub initial_horizon_s: f64,
    /// Hard final horizon in seconds. A live run at this point is censored.
    pub maximum_horizon_s: f64,
    /// Strictly greater than one multiplier applied between continuation rungs.
    pub multiplier: f64,
    /// Maximum number of continuation requests after the first horizon.
    pub maximum_extensions: u32,
}

/// Invalid continuation-policy field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HorizonPolicyError {
    InvalidInitialHorizon,
    InvalidMaximumHorizon,
    InvalidMultiplier,
    InvalidCurrentHorizon,
    CurrentHorizonExceedsMaximum,
    InvalidCensorTime,
    CensorTimeDoesNotMatchCurrentHorizon,
}

impl fmt::Display for HorizonPolicyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for HorizonPolicyError {}

impl HorizonContinuationPolicy {
    /// Validates the finite, bounded continuation schedule.
    pub fn validate(self) -> Result<(), HorizonPolicyError> {
        if !self.initial_horizon_s.is_finite() || self.initial_horizon_s <= 0.0 {
            return Err(HorizonPolicyError::InvalidInitialHorizon);
        }
        if !self.maximum_horizon_s.is_finite() || self.maximum_horizon_s < self.initial_horizon_s {
            return Err(HorizonPolicyError::InvalidMaximumHorizon);
        }
        if !self.multiplier.is_finite() || self.multiplier <= 1.0 {
            return Err(HorizonPolicyError::InvalidMultiplier);
        }
        Ok(())
    }

    /// Computes the next declared horizon for a right-censored run.
    ///
    /// `None` means the retained censoring horizon is final. A physical event
    /// or numerical refusal is never continued by this policy.
    pub fn next_horizon_s(
        self,
        outcome: RunOutcome,
        current_horizon_s: f64,
        extensions_already_requested: u32,
    ) -> Result<Option<f64>, HorizonPolicyError> {
        self.validate()?;
        if !current_horizon_s.is_finite() || current_horizon_s <= 0.0 {
            return Err(HorizonPolicyError::InvalidCurrentHorizon);
        }
        if current_horizon_s > self.maximum_horizon_s {
            return Err(HorizonPolicyError::CurrentHorizonExceedsMaximum);
        }
        let RunOutcome::RightCensored { censor_time_s } = outcome else {
            return Ok(None);
        };
        if !censor_time_s.is_finite() || censor_time_s < 0.0 {
            return Err(HorizonPolicyError::InvalidCensorTime);
        }
        if !same_declared_horizon(censor_time_s, current_horizon_s) {
            return Err(HorizonPolicyError::CensorTimeDoesNotMatchCurrentHorizon);
        }
        if extensions_already_requested >= self.maximum_extensions
            || current_horizon_s >= self.maximum_horizon_s
        {
            return Ok(None);
        }
        let expanded = current_horizon_s * self.multiplier;
        Ok(Some(expanded.min(self.maximum_horizon_s)))
    }
}

/// Caller-declared evidence that all three numerical rungs stayed in one
/// smooth mode. It is intentionally separate from terminal classification.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RefinementMode {
    /// A named smooth branch, identical at every rung before order is reported.
    Smooth { branch_id: String },
    /// An eventful or switching branch: compare differences, but do not infer
    /// a smooth observed order from them.
    Eventful { reason: String },
    /// The caller could not establish a common branch.
    Unresolved { reason: String },
}

/// Non-dimensionalization scales for reported final-QoI and energy deltas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConvergenceScales {
    pub inclination_rad: f64,
    pub precession_rad_per_s: f64,
    pub spin_rad_per_s: f64,
    pub work_j: f64,
    pub energy_j: f64,
}

/// Three resolution rungs for a single fixed physical case.
#[derive(Clone, Debug)]
pub struct ThreeRungConvergence<'a> {
    pub coarse: &'a CoupledRun,
    pub fine: &'a CoupledRun,
    pub reference: &'a CoupledRun,
    pub coarse_timestep_s: f64,
    pub fine_timestep_s: f64,
    pub reference_timestep_s: f64,
    pub mode: RefinementMode,
    pub scales: ConvergenceScales,
}

/// One normalized final-QoI difference between two rungs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinalQoiDelta {
    pub inclination: f64,
    pub precession: f64,
    pub spin: f64,
}

/// Work and energy differences between two rungs, normalized by declared
/// scales. Channel order is gravity, contact, rolling, base, gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WorkEnergyDelta {
    pub channel_work: [f64; 5],
    pub energy_defect: f64,
}

/// Observed order is only meaningful for a smooth common branch and positive
/// finite successive differences.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ObservedOrder {
    Available {
        inclination: f64,
        precession: f64,
        spin: f64,
        channel_work_linf: f64,
        energy_defect: f64,
    },
    NotApplicable {
        reason: OrderUnavailableReason,
    },
}

/// Why the module deliberately withheld an observed order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrderUnavailableReason {
    TerminalClassDisagreement,
    IncompatibleTerminalDisposition,
    NonSmoothOrUnresolvedMode,
    NonHalvingTimesteps,
    ExactOrNonMonotoneDifferences,
}

/// Convergence report for `h`, `h/2`, and `h/4` runs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConvergenceReceipt {
    pub coarse_outcome: RunOutcome,
    pub fine_outcome: RunOutcome,
    pub reference_outcome: RunOutcome,
    pub terminal_class_agreement: bool,
    /// Present only when all three runs reached the same physical terminal kind.
    pub coarse_fine_event_time_delta_s: Option<f64>,
    /// Present only when all three runs reached the same physical terminal kind.
    pub fine_reference_event_time_delta_s: Option<f64>,
    pub coarse_fine_qoi: FinalQoiDelta,
    pub fine_reference_qoi: FinalQoiDelta,
    pub coarse_fine_work_energy: WorkEnergyDelta,
    pub fine_reference_work_energy: WorkEnergyDelta,
    pub observed_order: ObservedOrder,
}

/// Refusal from refinement analysis. This is software/input validation, not a
/// physical model failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConvergenceError {
    Outcome(OutcomeError),
    MissingRetainedSample { rung: &'static str },
    InvalidTimestep { rung: &'static str },
    InvalidScale { field: &'static str },
    NonFiniteRetainedValue { rung: &'static str },
    NonFiniteNormalizedDelta { field: &'static str },
}

impl fmt::Display for ConvergenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ConvergenceError {}

/// Analyses three fixed-case resolution rungs without applying an arbitrary
/// pass/fail threshold. Consumers must compare these declared deltas against
/// a preregistered acceptance band appropriate to the QoI and regime.
pub fn analyse_three_rung_convergence(
    input: ThreeRungConvergence<'_>,
) -> Result<ConvergenceReceipt, ConvergenceError> {
    validate_timestep(input.coarse_timestep_s, "coarse")?;
    validate_timestep(input.fine_timestep_s, "fine")?;
    validate_timestep(input.reference_timestep_s, "reference")?;
    validate_scales(input.scales)?;

    let coarse_outcome = classify_outcome(input.coarse).map_err(ConvergenceError::Outcome)?;
    let fine_outcome = classify_outcome(input.fine).map_err(ConvergenceError::Outcome)?;
    let reference_outcome = classify_outcome(input.reference).map_err(ConvergenceError::Outcome)?;
    let coarse = retained_final(input.coarse, "coarse")?;
    let fine = retained_final(input.fine, "fine")?;
    let reference = retained_final(input.reference, "reference")?;

    let terminal_class_agreement = coarse_outcome.class() == fine_outcome.class()
        && fine_outcome.class() == reference_outcome.class();
    let terminal_disposition_orderable =
        orderable_terminal_disposition(coarse_outcome, fine_outcome, reference_outcome);
    let event_times = same_physical_terminal_times(coarse_outcome, fine_outcome, reference_outcome);
    let coarse_fine_qoi = qoi_delta(coarse, fine, input.scales)?;
    let fine_reference_qoi = qoi_delta(fine, reference, input.scales)?;
    let coarse_fine_work_energy = work_energy_delta(input.coarse, input.fine, input.scales)?;
    let fine_reference_work_energy = work_energy_delta(input.fine, input.reference, input.scales)?;
    let observed_order = observed_order(
        terminal_class_agreement,
        terminal_disposition_orderable,
        &input.mode,
        input.coarse_timestep_s,
        input.fine_timestep_s,
        input.reference_timestep_s,
        coarse_fine_qoi,
        fine_reference_qoi,
        coarse_fine_work_energy,
        fine_reference_work_energy,
    );

    Ok(ConvergenceReceipt {
        coarse_outcome,
        fine_outcome,
        reference_outcome,
        terminal_class_agreement,
        coarse_fine_event_time_delta_s: event_times.map(|times| (times.0 - times.1).abs()),
        fine_reference_event_time_delta_s: event_times.map(|times| (times.1 - times.2).abs()),
        coarse_fine_qoi,
        fine_reference_qoi,
        coarse_fine_work_energy,
        fine_reference_work_energy,
        observed_order,
    })
}

fn validate_timestep(value: f64, rung: &'static str) -> Result<(), ConvergenceError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(ConvergenceError::InvalidTimestep { rung })
    }
}

fn validate_scales(scales: ConvergenceScales) -> Result<(), ConvergenceError> {
    for (field, value) in [
        ("inclination_rad", scales.inclination_rad),
        ("precession_rad_per_s", scales.precession_rad_per_s),
        ("spin_rad_per_s", scales.spin_rad_per_s),
        ("work_j", scales.work_j),
        ("energy_j", scales.energy_j),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(ConvergenceError::InvalidScale { field });
        }
    }
    Ok(())
}

fn retained_final<'run>(
    run: &'run CoupledRun,
    rung: &'static str,
) -> Result<&'run crate::coupled_runner::CoupledSample, ConvergenceError> {
    let sample = run
        .samples
        .last()
        .ok_or(ConvergenceError::MissingRetainedSample { rung })?;
    let values = [
        sample.inclination_rad,
        sample.precession_rad_per_s,
        sample.spin_rad_per_s,
        sample.energy_defect_j,
    ];
    if values.into_iter().all(f64::is_finite)
        && run
            .checkpoint
            .accumulated_channel_work_j
            .into_iter()
            .all(f64::is_finite)
    {
        Ok(sample)
    } else {
        Err(ConvergenceError::NonFiniteRetainedValue { rung })
    }
}

fn qoi_delta(
    left: &crate::coupled_runner::CoupledSample,
    right: &crate::coupled_runner::CoupledSample,
    scales: ConvergenceScales,
) -> Result<FinalQoiDelta, ConvergenceError> {
    Ok(FinalQoiDelta {
        inclination: normalized_abs_delta(
            left.inclination_rad,
            right.inclination_rad,
            scales.inclination_rad,
            "inclination",
        )?,
        precession: normalized_abs_delta(
            left.precession_rad_per_s,
            right.precession_rad_per_s,
            scales.precession_rad_per_s,
            "precession",
        )?,
        spin: normalized_abs_delta(
            left.spin_rad_per_s,
            right.spin_rad_per_s,
            scales.spin_rad_per_s,
            "spin",
        )?,
    })
}

fn work_energy_delta(
    left: &CoupledRun,
    right: &CoupledRun,
    scales: ConvergenceScales,
) -> Result<WorkEnergyDelta, ConvergenceError> {
    let mut channel_work = [0.0; 5];
    for (index, value) in channel_work.iter_mut().enumerate() {
        *value = normalized_abs_delta(
            left.checkpoint.accumulated_channel_work_j[index],
            right.checkpoint.accumulated_channel_work_j[index],
            scales.work_j,
            "channel_work",
        )?;
    }
    Ok(WorkEnergyDelta {
        channel_work,
        energy_defect: normalized_abs_delta(
            left.checkpoint.accumulated_energy_defect_j,
            right.checkpoint.accumulated_energy_defect_j,
            scales.energy_j,
            "energy_defect",
        )?,
    })
}

fn normalized_abs_delta(
    left: f64,
    right: f64,
    scale: f64,
    field: &'static str,
) -> Result<f64, ConvergenceError> {
    let value = (left - right).abs() / scale;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ConvergenceError::NonFiniteNormalizedDelta { field })
    }
}

fn same_physical_terminal_times(
    coarse: RunOutcome,
    fine: RunOutcome,
    reference: RunOutcome,
) -> Option<(f64, f64, f64)> {
    match (coarse, fine, reference) {
        (
            RunOutcome::PhysicalTerminal {
                kind: coarse_kind,
                event_time_s: coarse_time,
            },
            RunOutcome::PhysicalTerminal {
                kind: fine_kind,
                event_time_s: fine_time,
            },
            RunOutcome::PhysicalTerminal {
                kind: reference_kind,
                event_time_s: reference_time,
            },
        ) if coarse_kind == fine_kind && fine_kind == reference_kind => {
            Some((coarse_time, fine_time, reference_time))
        }
        _ => None,
    }
}

/// A physical event needs the same event kind; a censored trajectory needs the
/// same declared censoring time. Numerical refusals never yield an order.
fn orderable_terminal_disposition(
    coarse: RunOutcome,
    fine: RunOutcome,
    reference: RunOutcome,
) -> bool {
    match (coarse, fine, reference) {
        (
            RunOutcome::PhysicalTerminal {
                kind: coarse_kind, ..
            },
            RunOutcome::PhysicalTerminal {
                kind: fine_kind, ..
            },
            RunOutcome::PhysicalTerminal {
                kind: reference_kind,
                ..
            },
        ) => coarse_kind == fine_kind && fine_kind == reference_kind,
        (
            RunOutcome::RightCensored {
                censor_time_s: coarse_time,
            },
            RunOutcome::RightCensored {
                censor_time_s: fine_time,
            },
            RunOutcome::RightCensored {
                censor_time_s: reference_time,
            },
        ) => {
            same_declared_horizon(coarse_time, fine_time)
                && same_declared_horizon(fine_time, reference_time)
        }
        _ => false,
    }
}

#[allow(clippy::too_many_arguments)]
fn observed_order(
    terminal_class_agreement: bool,
    terminal_disposition_orderable: bool,
    mode: &RefinementMode,
    coarse_timestep_s: f64,
    fine_timestep_s: f64,
    reference_timestep_s: f64,
    coarse_fine_qoi: FinalQoiDelta,
    fine_reference_qoi: FinalQoiDelta,
    coarse_fine_work_energy: WorkEnergyDelta,
    fine_reference_work_energy: WorkEnergyDelta,
) -> ObservedOrder {
    if !terminal_class_agreement {
        return ObservedOrder::NotApplicable {
            reason: OrderUnavailableReason::TerminalClassDisagreement,
        };
    }
    if !terminal_disposition_orderable {
        return ObservedOrder::NotApplicable {
            reason: OrderUnavailableReason::IncompatibleTerminalDisposition,
        };
    }
    if !matches!(mode, RefinementMode::Smooth { .. }) {
        return ObservedOrder::NotApplicable {
            reason: OrderUnavailableReason::NonSmoothOrUnresolvedMode,
        };
    }
    if !approximately_halves(coarse_timestep_s, fine_timestep_s)
        || !approximately_halves(fine_timestep_s, reference_timestep_s)
    {
        return ObservedOrder::NotApplicable {
            reason: OrderUnavailableReason::NonHalvingTimesteps,
        };
    }
    let coarse_work = linf(coarse_fine_work_energy.channel_work);
    let fine_work = linf(fine_reference_work_energy.channel_work);
    let orders = [
        observed_order_scalar(coarse_fine_qoi.inclination, fine_reference_qoi.inclination),
        observed_order_scalar(coarse_fine_qoi.precession, fine_reference_qoi.precession),
        observed_order_scalar(coarse_fine_qoi.spin, fine_reference_qoi.spin),
        observed_order_scalar(coarse_work, fine_work),
        observed_order_scalar(
            coarse_fine_work_energy.energy_defect,
            fine_reference_work_energy.energy_defect,
        ),
    ];
    let [
        Some(inclination),
        Some(precession),
        Some(spin),
        Some(channel_work_linf),
        Some(energy_defect),
    ] = orders
    else {
        return ObservedOrder::NotApplicable {
            reason: OrderUnavailableReason::ExactOrNonMonotoneDifferences,
        };
    };
    ObservedOrder::Available {
        inclination,
        precession,
        spin,
        channel_work_linf,
        energy_defect,
    }
}

fn approximately_halves(coarse: f64, fine: f64) -> bool {
    (coarse - 2.0 * fine).abs() <= 16.0 * f64::EPSILON * coarse.abs().max(fine.abs())
}

fn same_declared_horizon(censor_time_s: f64, current_horizon_s: f64) -> bool {
    (censor_time_s - current_horizon_s).abs()
        <= HORIZON_MATCH_RELATIVE_TOLERANCE * censor_time_s.abs().max(current_horizon_s.abs())
}

fn linf(values: [f64; 5]) -> f64 {
    values.into_iter().fold(0.0_f64, f64::max)
}

fn observed_order_scalar(coarse_difference: f64, fine_difference: f64) -> Option<f64> {
    if coarse_difference.is_finite()
        && fine_difference.is_finite()
        && coarse_difference > fine_difference
        && fine_difference > 0.0
    {
        Some((coarse_difference / fine_difference).log2())
    } else {
        None
    }
}

/// Named evidence categories that must be retained before calibration or blind
/// scoring is even structurally ready.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CalibrationEvidenceKind {
    Specimen,
    Rig,
    Instrument,
    RawObservations,
    ObservationCovariance,
    CalibrationPartition,
    BlindHoldout,
}

impl CalibrationEvidenceKind {
    /// Stable diagnostic name for a required evidence category.
    ///
    /// These names are for an explicit no-data record only.  They do not
    /// promote a declared identity into retained data or physical authority.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::Specimen => "specimen",
            Self::Rig => "rig",
            Self::Instrument => "instrument",
            Self::RawObservations => "raw-observations",
            Self::ObservationCovariance => "observation-covariance",
            Self::CalibrationPartition => "calibration-partition",
            Self::BlindHoldout => "blind-holdout",
        }
    }
}

/// An evidence identity is only a declared binding. It is not an authenticated
/// artifact, a calibrated measurement, or a physical-validation promotion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DeclaredEvidence {
    Present { identity: String },
    Missing,
}

/// The minimum artifact bindings required before a later calibration pipeline
/// may attempt fitting or a blind-score join.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CalibrationReadinessInput {
    pub specimen: DeclaredEvidence,
    pub rig: DeclaredEvidence,
    pub instrument: DeclaredEvidence,
    pub raw_observations: DeclaredEvidence,
    pub observation_covariance: DeclaredEvidence,
    pub calibration_partition: DeclaredEvidence,
    pub blind_holdout: DeclaredEvidence,
}

/// Report every absent prerequisite in stable admission order.
///
/// [`admit_calibration_readiness`] returns the first invalid field so callers
/// can fix one malformed binding directly.  A no-data report instead needs the
/// complete absence set: otherwise a human-facing diagnostic can silently
/// become stale as its typed input evolves.  This function performs no fitting
/// and treats invalid-but-present identities separately from absence.
#[must_use]
pub fn missing_calibration_evidence(
    input: &CalibrationReadinessInput,
) -> Vec<CalibrationEvidenceKind> {
    [
        (CalibrationEvidenceKind::Specimen, &input.specimen),
        (CalibrationEvidenceKind::Rig, &input.rig),
        (CalibrationEvidenceKind::Instrument, &input.instrument),
        (
            CalibrationEvidenceKind::RawObservations,
            &input.raw_observations,
        ),
        (
            CalibrationEvidenceKind::ObservationCovariance,
            &input.observation_covariance,
        ),
        (
            CalibrationEvidenceKind::CalibrationPartition,
            &input.calibration_partition,
        ),
        (CalibrationEvidenceKind::BlindHoldout, &input.blind_holdout),
    ]
    .into_iter()
    .filter_map(|(kind, evidence)| matches!(evidence, DeclaredEvidence::Missing).then_some(kind))
    .collect()
}

/// Structurally complete calibration-data bindings. This deliberately carries
/// no numeric observation or parameter data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CalibrationReadinessReceipt {
    pub specimen_identity: String,
    pub rig_identity: String,
    pub instrument_identity: String,
    pub raw_observations_identity: String,
    pub observation_covariance_identity: String,
    pub calibration_partition_identity: String,
    pub blind_holdout_identity: String,
}

/// Readiness failures are explicit rather than silently falling back to
/// synthetic or target-fitted data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CalibrationReadinessError {
    MissingEvidence { kind: CalibrationEvidenceKind },
    InvalidEvidenceIdentity { kind: CalibrationEvidenceKind },
    PartitionAlias,
}

impl fmt::Display for CalibrationReadinessError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for CalibrationReadinessError {}

/// Admit the complete list of declared evidence bindings.
///
/// This function intentionally stops at structural readiness: independently
/// authenticated custody, calibrated metrology, likelihood construction, and
/// physical validation are separate required stages.
pub fn admit_calibration_readiness(
    input: CalibrationReadinessInput,
) -> Result<CalibrationReadinessReceipt, CalibrationReadinessError> {
    let specimen_identity = evidence_identity(input.specimen, CalibrationEvidenceKind::Specimen)?;
    let rig_identity = evidence_identity(input.rig, CalibrationEvidenceKind::Rig)?;
    let instrument_identity =
        evidence_identity(input.instrument, CalibrationEvidenceKind::Instrument)?;
    let raw_observations_identity = evidence_identity(
        input.raw_observations,
        CalibrationEvidenceKind::RawObservations,
    )?;
    let observation_covariance_identity = evidence_identity(
        input.observation_covariance,
        CalibrationEvidenceKind::ObservationCovariance,
    )?;
    let calibration_partition_identity = evidence_identity(
        input.calibration_partition,
        CalibrationEvidenceKind::CalibrationPartition,
    )?;
    let blind_holdout_identity =
        evidence_identity(input.blind_holdout, CalibrationEvidenceKind::BlindHoldout)?;
    if calibration_partition_identity == blind_holdout_identity {
        return Err(CalibrationReadinessError::PartitionAlias);
    }
    Ok(CalibrationReadinessReceipt {
        specimen_identity,
        rig_identity,
        instrument_identity,
        raw_observations_identity,
        observation_covariance_identity,
        calibration_partition_identity,
        blind_holdout_identity,
    })
}

fn evidence_identity(
    evidence: DeclaredEvidence,
    kind: CalibrationEvidenceKind,
) -> Result<String, CalibrationReadinessError> {
    match evidence {
        DeclaredEvidence::Missing => Err(CalibrationReadinessError::MissingEvidence { kind }),
        DeclaredEvidence::Present { identity }
            if identity.trim().is_empty() || identity.chars().any(char::is_control) =>
        {
            Err(CalibrationReadinessError::InvalidEvidenceIdentity { kind })
        }
        DeclaredEvidence::Present { identity } => Ok(identity),
    }
}
