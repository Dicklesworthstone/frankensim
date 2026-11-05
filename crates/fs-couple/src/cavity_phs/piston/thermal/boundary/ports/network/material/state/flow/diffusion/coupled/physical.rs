//! One complete-assembly time controller for closed and reservoir transport.
//!
//! The physical state comparison is mandatory. A transport may add directional
//! stream errors, but cannot replace the original storage-coordinate checks.
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::{
    GasTransportError, GasTransportStep,
};

/// Control and accounting failures contain no constitutive callback error.
/// Keeping this set separate lets every adapter retain its original typed
/// physical refusal without an erased source or an unreachable conversion arm.
#[derive(Debug)]
pub(crate) enum PhysicalControlError {
    Invalid(&'static str),
    Balance(&'static str),
    Circuit(crate::material_circuit::MaterialCircuitError),
    Gas {
        chamber: u64,
        error: fs_material::MaterialError,
    },
    Cancelled,
}
impl<E, S> From<PhysicalControlError> for CoupledPhysicalError<E, S> {
    fn from(error: PhysicalControlError) -> Self {
        match error {
            PhysicalControlError::Invalid(reason) => Self::Invalid(reason),
            PhysicalControlError::Balance(reason) => Self::Balance(reason),
            PhysicalControlError::Circuit(error) => Self::Circuit(error),
            PhysicalControlError::Gas { chamber, error } => Self::Gas { chamber, error },
            PhysicalControlError::Cancelled => Self::Cancelled,
        }
    }
}
impl From<PhysicalControlError> for CoupledDiffusionError {
    fn from(error: PhysicalControlError) -> Self {
        match error {
            PhysicalControlError::Invalid(reason) => Self::Invalid(reason),
            PhysicalControlError::Balance(reason) => Self::Balance(reason),
            PhysicalControlError::Circuit(error) => Self::Circuit(error),
            PhysicalControlError::Gas { chamber, error } => Self::Gas { chamber, error },
            PhysicalControlError::Cancelled => Self::Cancelled,
        }
    }
}

/// Classify the original failure by type without inspecting or cloning a
/// caller's source payload. A smaller step cannot cure a constitutive refusal.
pub(crate) trait PhysicalStepError {
    fn is_cancelled(&self) -> bool;
    fn refinable(&self) -> bool;
}
impl PhysicalStepError for MaterialFlowError {
    fn is_cancelled(&self) -> bool {
        matches!(self, Self::Cancelled)
    }
    fn refinable(&self) -> bool {
        super::refinable(self)
    }
}
impl<E> PhysicalStepError for MaterialGasLawFlowError<E> {
    fn is_cancelled(&self) -> bool {
        matches!(self, Self::Cancelled)
    }
    fn refinable(&self) -> bool {
        match self {
            Self::Exchange(error) => feedback_refinable(error),
            Self::Transport(error) => native_transport_refinable(error),
            #[cfg(feature = "molecular-transport")]
            Self::Molecular(
                fs_gas::material::transport::parcel::ParcelMolecularError::Transport(error),
            ) => native_transport_refinable(error),
            _ => false,
        }
    }
}
fn native_transport_refinable(error: &GasTransportError) -> bool {
    // Delegate the native finite-inventory classification unchanged. This clone
    // contains only the native transport error, never a caller-selected source.
    <GasTransportStep as TransportLedger>::refinable(&MaterialFlowError::Transport(error.clone()))
}

/// One actual accepted physical half step and its original source receipts.
/// Coarse probes and rejected trials never enter the accepted trajectory.
pub struct PhysicalSubstep<A> {
    pub start_s: f64,
    pub end_s: f64,
    pub advance: A,
}

pub(crate) trait PhysicalEstimate: Copy + core::fmt::Debug {
    fn from_state(state: CoupledDiffusionEstimate) -> Self;
    fn maximum(self) -> f64;
    fn accumulate(&mut self, other: Self);
}
impl PhysicalEstimate for CoupledDiffusionEstimate {
    fn from_state(state: Self) -> Self {
        state
    }
    fn maximum(self) -> f64 {
        CoupledDiffusionEstimate::maximum(self)
    }
    fn accumulate(&mut self, other: Self) {
        CoupledDiffusionEstimate::accumulate(self, other);
    }
}

pub(crate) struct PhysicalWindow<A, E> {
    pub start_s: f64,
    pub end_s: f64,
    pub evaluations: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    pub estimated_error: E,
    pub steps: Vec<PhysicalSubstep<A>>,
}

impl<'m> MaterialExchangeState<'m> {
    /// Run the original first-order physical split with complete-state local
    /// error control and optional source-owned integrated-stream comparisons.
    /// The returned owners remain private until the caller checks its window
    /// balances and publishes them. Native residual budgets never scale with
    /// the number of accepted steps.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn coupled_physical_window_candidate<
        A,
        E: PhysicalEstimate,
        S: PhysicalStepError,
    >(
        &self,
        cx: &Cx<'_>,
        target_s: f64,
        drives: &[ChamberDrive],
        guess_j: &[f64],
        time: CoupledDiffusionControls,
        mut evolve: impl FnMut(&Self, f64, &[ChamberDrive], &[f64]) -> Result<(Self, A), S>,
        mut additional_error: impl FnMut(CoupledDiffusionEstimate, &A, &A, &A, f64) -> Result<E, S>,
        transport_refinable: impl Fn(&S) -> bool,
    ) -> Result<(Self, PhysicalWindow<A, E>), CoupledPhysicalError<E, S>> {
        use CoupledPhysicalError as Failure;
        check(cx)?;
        time.validate()?;
        let start_s = self.time_s();
        let duration = finite(target_s - start_s)?;
        if duration <= 0.0
            || guess_j.len() != self.thermal.cells().len()
            || guess_j.iter().any(|x| !x.is_finite())
            || drives.iter().any(|d| !d.heat_j.is_finite())
            || self
                .thermal
                .pending_events()
                .first()
                .is_some_and(|e| e.time_s <= target_s)
        {
            return Err(Failure::Invalid(
                "one forward fixed-assembly window, finite inputs and no pending events through target required",
            ));
        }
        let cap = time
            .maximum_evaluations
            .min(self.thermal.work_budget().max_substeps);
        let mut next = Self::from_parts(
            self.thermal.clone(),
            self.circuit.clone(),
            self.chambers.clone(),
        );
        let mut steps = Vec::new();
        let (mut evaluations, mut accepted_intervals, mut rejected_intervals) = (0, 0, 0);
        let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
        let mut estimated_error = E::from_state(CoupledDiffusionEstimate::default());
        let mut step_s = time.maximum_step_s.min(duration);
        while next.time_s() < target_s {
            check(cx)?;
            let before = next.time_s();
            let remaining = finite(target_s - before)?;
            let dt = step_s.min(remaining);
            let end = if dt == remaining {
                target_s
            } else {
                finite(before + dt)?
            };
            let dt = finite(end - before)?;
            let midpoint = finite(before + 0.5 * (end - before))?;
            if midpoint <= before || midpoint >= end {
                return Err(Failure::Invalid(
                    "coupled diffusion half-step clock is unrepresentable",
                ));
            }
            let fine_error_factor = error::fine_error_factor(midpoint - before, end - midpoint)
                .ok_or(Failure::Invalid(
                    "coupled diffusion half-step ratio is unrepresentable",
                ))?;
            let accepted_substeps = steps.len();
            let mut evaluate = |state: &Self,
                                end: f64,
                                evaluations: &mut usize|
             -> Result<(Self, A), CoupledPhysicalError<E, S>> {
                check(cx)?;
                if *evaluations >= cap {
                    return Err(Failure::WorkLimit {
                        time_s: before,
                        evaluations: *evaluations,
                        accepted_substeps,
                    });
                }
                let left = (state.time_s() - start_s) / duration;
                let right = if end == target_s {
                    1.0
                } else {
                    (end - start_s) / duration
                };
                let subdrives: Vec<_> = drives
                    .iter()
                    .map(|drive| ChamberDrive {
                        heat_j: drive.heat_j * right - drive.heat_j * left,
                        ..*drive
                    })
                    .collect();
                let fraction = (end - state.time_s()) / duration;
                let guess: Vec<_> = guess_j.iter().map(|value| value * fraction).collect();
                *evaluations += 1;
                let (candidate, report) =
                    evolve(state, end, &subdrives, &guess).map_err(|error| {
                        if error.is_cancelled() {
                            Failure::Cancelled
                        } else {
                            Failure::Step {
                                start_s: state.time_s(),
                                end_s: end,
                                error,
                            }
                        }
                    })?;
                fixed(self, &candidate)?;
                Ok((candidate, report))
            };
            let comparison = (|| {
                let (coarse, coarse_report) = evaluate(&next, end, &mut evaluations)?;
                let (half, first) = evaluate(&next, midpoint, &mut evaluations)?;
                let (fine, second) = evaluate(&half, end, &mut evaluations)?;
                let state_error = compare(&next, &coarse, &fine, time, fine_error_factor, cx)?;
                let estimate = additional_error(
                    state_error,
                    &coarse_report,
                    &first,
                    &second,
                    fine_error_factor,
                )
                .map_err(|error| {
                    if error.is_cancelled() {
                        Failure::Cancelled
                    } else {
                        Failure::Step {
                            start_s: before,
                            end_s: end,
                            error,
                        }
                    }
                })?;
                // Additional stream coordinates cannot suppress the mandatory
                // complete physical-state estimate, even for a future adapter.
                let ratio = state_error.maximum().max(estimate.maximum());
                Ok((fine, first, second, estimate, ratio))
            })();
            let reduction = match comparison {
                Ok((fine, first, second, estimate, ratio)) if ratio <= 1.0 => {
                    estimated_error.accumulate(estimate);
                    accepted_intervals += 1;
                    minimum = minimum.min(midpoint - before).min(end - midpoint);
                    maximum = maximum.max(midpoint - before).max(end - midpoint);
                    steps.push(PhysicalSubstep {
                        start_s: before,
                        end_s: midpoint,
                        advance: first,
                    });
                    steps.push(PhysicalSubstep {
                        start_s: midpoint,
                        end_s: end,
                        advance: second,
                    });
                    next = fine;
                    let factor = error::step_factor(ratio, true);
                    step_s =
                        error::bounded_step_s(dt, factor, time.minimum_step_s, time.maximum_step_s);
                    continue;
                }
                Ok((_, _, _, estimate, ratio)) => {
                    let failure = Failure::Accuracy {
                        start_s: before,
                        step_s: dt,
                        estimate,
                    };
                    (error::step_factor(ratio, false), failure)
                }
                Err(Failure::Step {
                    start_s,
                    end_s,
                    error,
                }) if error.refinable() || transport_refinable(&error) => (
                    0.5,
                    Failure::Step {
                        start_s,
                        end_s,
                        error,
                    },
                ),
                Err(error) => return Err(error),
            };
            rejected_intervals += 1;
            if dt <= time.minimum_step_s {
                return Err(reduction.1);
            }
            let smaller = error::bounded_step_s(dt, reduction.0, time.minimum_step_s, dt);
            if !smaller.is_finite() || smaller >= dt || before + smaller >= end {
                return Err(reduction.1);
            }
            step_s = smaller;
        }
        check(cx)?;
        Ok((
            next,
            PhysicalWindow {
                start_s,
                end_s: target_s,
                evaluations,
                accepted_intervals,
                rejected_intervals,
                minimum_accepted_substep_s: minimum,
                maximum_accepted_substep_s: maximum,
                estimated_error,
                steps,
            },
        ))
    }
}
