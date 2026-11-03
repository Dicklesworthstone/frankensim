//! Simultaneous bilinear batch reactions competing for one shared reagent.
//!
//! Each original [`BilinearBatchReaction`] consumes the shared reagent and one
//! distinct other species. No product may be a reactant anywhere in this batch.
//! All rates use the same staged inventory at the same absolute elapsed time:
//! `dξ_i/dt = k_i n_Ai n_Bi / V`. There is no independent-recovery clipping,
//! reaction priority, substep product funding, or inferred chemical mechanism.
//!
//! Classical RK4 with step doubling advances extents. Every stage is rebuilt
//! from the original stocks by [`ReactionNetwork`]; accepted changes therefore
//! retain its exact equation-conservation owner and shared-demand checks.
//! Positive initial reactants must remain positive at finite time. Unsupported
//! trial stages halve the step; work exhaustion returns no partial result.
//!
//! The caller supplies a maximum step, attempt budget and extent-error budget.
//! The difference of one full and two half steps, divided by 15, estimates the
//! local error of the accepted half steps. Time-weighted acceptance and a check
//! on accumulated estimates enforce that *estimated* budget. This is not a
//! rigorous global-error bound, stiffness solver, or calibration certificate.
//! All original source conditions remain required at every derivative stage.
//! Reactants decrease monotonically, so their admitted source endpoints cover
//! the intervening path under the original rectangular applicability domains.

use super::{BatchReactionError, BilinearBatchReaction, OWNED_AXES};
use fs_matdb::QueryPoint;
use fs_qty::{Amount, Dimensionless, Temperature, Time, Volume};
use fs_thermochem::{
    ReactionId, SpeciesId,
    stoichiometry::network::{NetworkAdvance, ReactionNetwork, ReactionNetworkError},
};
use std::collections::{BTreeMap, BTreeSet};

/// Bounded mixed-feed scope; the conserved inventory owner also bounds species.
pub const MAX_BATCH_REACTIONS: usize = 32;
/// Hard limit on caller-authorized full/half-step trials, including rejections.
pub const MAX_STEP_ATTEMPTS: usize = 65_536;

/// Explicit numerical controls, retained in the returned result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatchIntegrationBudget {
    /// Maximum physical duration of any proposed full step.
    pub maximum_step: Time,
    /// Absolute accumulated local-error estimate budget per reaction [mol extent].
    pub absolute_extent_error: Amount,
    /// Additional fraction of that reaction's initial stoichiometric capacity.
    pub relative_extent_error: Dimensionless,
    /// Maximum full/half-step attempts, accepted and rejected together.
    pub maximum_attempts: usize,
}

/// Atomic refusal; the original inventory is never modified.
#[derive(Clone, Debug, PartialEq)]
pub enum BatchNetworkError {
    /// Unsupported topology, incomplete conditions, or invalid numerical controls.
    Invalid(&'static str),
    /// Original source-card admission failed for the named reaction.
    Source {
        /// Original reaction whose declared support refused the state.
        reaction: ReactionId,
        /// Unchanged refusal from the original bilinear source owner.
        refusal: BatchReactionError,
    },
    /// Original conserved-network admission failed.
    Inventory(ReactionNetworkError),
    /// A positive initial reactant vanished in a finite-time trial.
    FiniteTimeDepletion {
        /// Consumed species with a positive initial stock.
        species: SpeciesId,
    },
    /// Positive rates, extents, clocks or error scales lost numerical range.
    Arithmetic(&'static str),
    /// No complete result was produced within the authorized work.
    BudgetExhausted {
        /// Trials performed, including rejections.
        attempts: usize,
        /// Successful steps kept only in the unreturned candidate.
        accepted_steps: usize,
        /// Why the last trial rejected or why another step was needed.
        last_trial: &'static str,
    },
}
impl core::fmt::Display for BatchNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shared-reagent batch: {self:?}")
    }
}
impl std::error::Error for BatchNetworkError {}
impl From<ReactionNetworkError> for BatchNetworkError {
    fn from(error: ReactionNetworkError) -> Self {
        Self::Inventory(error)
    }
}
type Result<T> = core::result::Result<T, BatchNetworkError>;

/// Original per-reaction source coordinates, in canonical reaction order.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionSourcePath {
    /// Original reaction identity.
    pub reaction: ReactionId,
    /// Actual initial concentrations and time zero.
    pub initial: QueryPoint,
    /// Actual final concentrations and total requested duration.
    pub final_point: QueryPoint,
}

/// Completed numerical batch, with original source and conservation owners.
#[derive(Clone, Debug)]
pub struct BatchNetworkAdvance {
    model: BilinearBatchNetwork,
    /// Complete canonical inventory before any reaction progressed.
    pub initial_inventory: Vec<(SpeciesId, Amount)>,
    /// Original network accounting: extents, gross flows and final inventory.
    pub stage: NetworkAdvance,
    /// Original source support for each actual concentration/time path.
    pub source_paths: Vec<ReactionSourcePath>,
    /// Fixed reference volume, shared by all empirical availability variables.
    pub volume: Volume,
    /// Total requested physical duration, never reset at numerical substeps.
    pub elapsed: Time,
    /// Fixed admitted absolute temperature in kelvin.
    pub temperature: Temperature,
    /// Caller-authorized numerical controls.
    pub budget: BatchIntegrationBudget,
    /// Accepted full time intervals, each evaluated using two half steps.
    pub accepted_steps: usize,
    /// Rejected full/half-step trials.
    pub rejected_steps: usize,
    /// Sum of accepted local estimates in mol extent; not certified global error.
    pub estimated_extent_errors: Vec<(ReactionId, Amount)>,
}
impl BatchNetworkAdvance {
    /// Exact original source cards, reactant mappings and conserved equations.
    pub fn model(&self) -> &BilinearBatchNetwork {
        &self.model
    }
}

/// Bounded, canonically ordered set of original bilinear kinetic declarations.
#[derive(Clone, Debug)]
pub struct BilinearBatchNetwork {
    models: Vec<BilinearBatchReaction>,
    network: ReactionNetwork,
    shared_reagent: SpeciesId,
    reactant_indices: Vec<[usize; 2]>,
}
impl BilinearBatchNetwork {
    /// Join distinct reaction IDs sharing exactly one consumed reagent. Each
    /// other reactant must be unique; products cannot feed any reaction here.
    /// More general feedback and equilibrium networks require another model.
    pub fn new(mut models: Vec<BilinearBatchReaction>, shared_reagent: SpeciesId) -> Result<Self> {
        if models.is_empty() || models.len() > MAX_BATCH_REACTIONS {
            return Err(BatchNetworkError::Invalid(
                "requires 1..=32 original bilinear models",
            ));
        }
        models.sort_by(|a, b| a.reaction.id().cmp(b.reaction.id()));
        let network = ReactionNetwork::new(models.iter().map(|m| m.reaction.clone()).collect())?;
        let mut other_reactants = BTreeSet::new();
        let mut all_reactants = BTreeSet::from([shared_reagent.clone()]);
        for model in &models {
            let shared = model
                .reactants
                .iter()
                .position(|id| id == &shared_reagent)
                .ok_or(BatchNetworkError::Invalid(
                    "every reaction must consume the shared reagent",
                ))?;
            let other = model.reactants[1 - shared].clone();
            if !other_reactants.insert(other.clone()) {
                return Err(BatchNetworkError::Invalid(
                    "each other reactant must be unique",
                ));
            }
            all_reactants.insert(other);
        }
        if models.iter().any(|model| {
            model
                .reaction
                .terms()
                .iter()
                .any(|term| term.coefficient > 0 && all_reactants.contains(&term.species))
        }) {
            return Err(BatchNetworkError::Invalid(
                "products cannot feed reactants in this batch",
            ));
        }
        let reactant_indices = models
            .iter()
            .map(|model| {
                core::array::from_fn(|i| {
                    network
                        .species()
                        .binary_search(&model.reactants[i])
                        .expect("original reaction participant belongs to its network")
                })
            })
            .collect();
        Ok(Self {
            models,
            network,
            shared_reagent,
            reactant_indices,
        })
    }

    /// Original models, sorted by original reaction identity.
    pub fn models(&self) -> &[BilinearBatchReaction] {
        &self.models
    }
    /// Existing exact conservation and source-bound heat-accounting owner.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.network
    }
    /// Explicit finite reagent inventory shared by all rates.
    pub fn shared_reagent(&self) -> &SpeciesId {
        &self.shared_reagent
    }

    /// Integrate all reaction rates simultaneously from one complete inventory.
    /// Supply exactly one additional-condition query per original reaction ID;
    /// owned time/volume/temperature/concentration coordinates are supplied here.
    /// Source validity is checked at absolute elapsed stage times and the full
    /// requested horizon. No partial result survives any late refusal.
    pub fn advance(
        &self,
        inventory: &[(SpeciesId, Amount)],
        volume: Volume,
        elapsed: Time,
        temperature: Temperature,
        conditions: &[(ReactionId, QueryPoint)],
        budget: BatchIntegrationBudget,
    ) -> Result<BatchNetworkAdvance> {
        validate_inputs(volume, elapsed, temperature, budget)?;
        if conditions.len() != self.models.len() {
            return Err(BatchNetworkError::Invalid(
                "supply one condition query per reaction",
            ));
        }
        let mut by_id = BTreeMap::new();
        for (id, point) in conditions {
            if by_id.insert(id, point).is_some() {
                return Err(BatchNetworkError::Invalid("duplicate reaction conditions"));
            }
        }
        let ordered = self
            .models
            .iter()
            .map(|model| {
                let point = *by_id
                    .get(model.reaction.id())
                    .ok_or(BatchNetworkError::Invalid(
                        "missing original reaction conditions",
                    ))?;
                if point.axes().keys().any(|axis| {
                    OWNED_AXES.contains(&axis.as_str())
                        || model.source.validity.bound(axis).is_none()
                }) {
                    return Err(BatchNetworkError::Invalid(
                        "conditions must contain only additional source-declared axes",
                    ));
                }
                Ok(point)
            })
            .collect::<Result<Vec<_>>>()?;
        let zeros = vec![0.0; self.models.len()];
        let initial = self
            .network
            .advance_inventory(inventory, &self.requests(&zeros))?
            .final_inventory;
        let problem = Problem {
            owner: self,
            initial: &initial,
            conditions: &ordered,
            volume,
            temperature,
        };
        let (_, initial_points) = problem.state(&zeros, 0.0)?;
        // Admit the requested total time before attempting numerical work.
        // These concentrations are only used to check the unchanged axes and
        // horizon; the actual final concentrations are checked again below.
        problem.state(&zeros, elapsed.value())?;
        let limits = self
            .models
            .iter()
            .zip(&self.reactant_indices)
            .map(|(model, indices)| {
                let capacity = (initial[indices[0]].1.value() / model.coefficients[0] as f64)
                    .min(initial[indices[1]].1.value() / model.coefficients[1] as f64);
                let limit = budget
                    .relative_extent_error
                    .value()
                    .mul_add(capacity, budget.absolute_extent_error.value());
                if !limit.is_finite() || limit <= 0.0 {
                    return Err(BatchNetworkError::Arithmetic(
                        "extent error budget range loss",
                    ));
                }
                Ok(limit)
            })
            .collect::<Result<Vec<_>>>()?;
        let mut extents = zeros;
        let mut errors = vec![0.0; self.models.len()];
        let inactive = elapsed.value() == 0.0
            || problem
                .rates(&extents, 0.0)?
                .iter()
                .all(|rate| *rate == 0.0);
        let mut clock = if inactive { elapsed.value() } else { 0.0 };
        let mut step = budget.maximum_step.value().min(elapsed.value());
        let mut attempts = 0;
        let mut accepted = 0;
        let mut rejected = 0;
        let mut last_trial = "another physical step is required";
        while clock < elapsed.value() {
            if attempts == budget.maximum_attempts {
                return Err(BatchNetworkError::BudgetExhausted {
                    attempts,
                    accepted_steps: accepted,
                    last_trial,
                });
            }
            let end = (clock + step).min(elapsed.value());
            let h = end - clock;
            if !end.is_finite()
                || end <= clock
                || clock + h * 0.5 == clock
                || clock + h * 0.5 == end
            {
                return Err(BatchNetworkError::Arithmetic(
                    "physical integration clock lost resolution",
                ));
            }
            attempts += 1;
            let trial = (|| {
                let coarse = problem.rk4(&extents, clock, end)?;
                let mid = clock + h * 0.5;
                let first = problem.rk4(&extents, clock, mid)?;
                let fine = problem.rk4(&first, mid, end)?;
                Ok((coarse, fine))
            })();
            let (coarse, fine) = match trial {
                Ok(values) => values,
                Err(error) => {
                    let Some(reason) = retry_reason(&error) else {
                        return Err(error);
                    };
                    rejected += 1;
                    last_trial = reason;
                    step = h * 0.5;
                    continue;
                }
            };
            let estimates = coarse
                .iter()
                .zip(&fine)
                .map(|(a, b)| local_error_estimate(*a, *b))
                .collect::<Result<Vec<_>>>()?;
            let fraction = h / elapsed.value();
            let fits = estimates.iter().enumerate().all(|(i, error)| {
                error.is_finite()
                    && *error <= limits[i] * fraction
                    && errors[i] + error <= limits[i]
            });
            if !fits {
                rejected += 1;
                last_trial = "local or accumulated extent-error estimate exceeded its budget";
                step = h * 0.5;
                continue;
            }
            let grow = estimates
                .iter()
                .enumerate()
                .all(|(i, error)| *error <= limits[i] * fraction / 32.0);
            for (total, local) in errors.iter_mut().zip(estimates) {
                *total += local;
            }
            extents = fine;
            clock = end;
            accepted += 1;
            last_trial = "another physical step is required";
            step = if grow {
                (h * 2.0).min(budget.maximum_step.value())
            } else {
                h
            };
        }
        let (stage, final_points) = problem.state(&extents, elapsed.value())?;
        let source_paths = self
            .models
            .iter()
            .zip(initial_points)
            .zip(final_points)
            .map(|((model, initial), final_point)| ReactionSourcePath {
                reaction: model.reaction.id().clone(),
                initial,
                final_point,
            })
            .collect();
        let estimated_extent_errors = self
            .models
            .iter()
            .zip(errors)
            .map(|(model, error)| (model.reaction.id().clone(), Amount::new(error)))
            .collect();
        Ok(BatchNetworkAdvance {
            model: self.clone(),
            initial_inventory: initial,
            stage,
            source_paths,
            volume,
            elapsed,
            temperature,
            budget,
            accepted_steps: accepted,
            rejected_steps: rejected,
            estimated_extent_errors,
        })
    }

    fn requests(&self, extents: &[f64]) -> Vec<(ReactionId, Amount)> {
        self.models
            .iter()
            .zip(extents)
            .map(|(model, value)| (model.reaction.id().clone(), Amount::new(*value)))
            .collect()
    }
}

fn validate_inputs(
    volume: Volume,
    elapsed: Time,
    temperature: Temperature,
    budget: BatchIntegrationBudget,
) -> Result<()> {
    if !volume.value().is_finite()
        || volume.value() <= 0.0
        || !temperature.value().is_finite()
        || temperature.value() <= 0.0
        || !elapsed.value().is_finite()
        || elapsed.value() < 0.0
        || !budget.maximum_step.value().is_finite()
        || budget.maximum_step.value() <= 0.0
        || !budget.absolute_extent_error.value().is_finite()
        || budget.absolute_extent_error.value() <= 0.0
        || !budget.relative_extent_error.value().is_finite()
        || !(0.0..1.0).contains(&budget.relative_extent_error.value())
        || !(1..=MAX_STEP_ATTEMPTS).contains(&budget.maximum_attempts)
    {
        return Err(BatchNetworkError::Invalid(
            "requires finite positive volume/kelvin/step/absolute error, nonnegative time, relative error in [0,1), and a bounded positive attempt count",
        ));
    }
    Ok(())
}

fn retry_reason(error: &BatchNetworkError) -> Option<&'static str> {
    match error {
        BatchNetworkError::Inventory(ReactionNetworkError::InsufficientInventory { .. }) => {
            Some("trial exceeded shared inventory")
        }
        BatchNetworkError::Source {
            refusal: BatchReactionError::OutsideValidity { .. },
            ..
        } => Some("trial left original source support"),
        BatchNetworkError::FiniteTimeDepletion { .. } => {
            Some("trial erased positive finite-time reactant stock")
        }
        _ => None,
    }
}

struct Problem<'a> {
    owner: &'a BilinearBatchNetwork,
    initial: &'a [(SpeciesId, Amount)],
    conditions: &'a [&'a QueryPoint],
    volume: Volume,
    temperature: Temperature,
}
impl Problem<'_> {
    fn state(&self, extents: &[f64], time: f64) -> Result<(NetworkAdvance, Vec<QueryPoint>)> {
        let stage = self
            .owner
            .network
            .advance_inventory(self.initial, &self.owner.requests(extents))?;
        let mut points = Vec::with_capacity(self.owner.models.len());
        for (i, model) in self.owner.models.iter().enumerate() {
            let indices = self.owner.reactant_indices[i];
            let n = core::array::from_fn(|j| stage.final_inventory[indices[j]].1.value());
            for j in 0..2 {
                if self.initial[indices[j]].1.value() > 0.0 && n[j] <= 0.0 {
                    return Err(BatchNetworkError::FiniteTimeDepletion {
                        species: model.reactants[j].clone(),
                    });
                }
            }
            let query = (|| {
                let point = model.point(
                    n,
                    self.volume,
                    Time::new(time),
                    self.temperature,
                    self.conditions[i],
                )?;
                model.admit(&point, "network-stage")?;
                Ok(point)
            })()
            .map_err(|refusal| BatchNetworkError::Source {
                reaction: model.reaction.id().clone(),
                refusal,
            })?;
            points.push(query);
        }
        Ok((stage, points))
    }

    fn rates(&self, extents: &[f64], time: f64) -> Result<Vec<f64>> {
        let (stage, _) = self.state(extents, time)?;
        self.owner
            .models
            .iter()
            .enumerate()
            .map(|(i, model)| {
                let indices = self.owner.reactant_indices[i];
                let a = stage.final_inventory[indices[0]].1.value();
                let b = stage.final_inventory[indices[1]].1.value();
                if a == 0.0 || b == 0.0 || model.rate.value() == 0.0 {
                    return Ok(0.0);
                }
                let rate = (model.rate.value() * (a / self.volume.value())) * b;
                if !rate.is_finite() || rate <= 0.0 {
                    return Err(BatchNetworkError::Arithmetic(
                        "positive simultaneous rate range loss",
                    ));
                }
                Ok(rate)
            })
            .collect()
    }

    fn rk4(&self, start: &[f64], begin: f64, end: f64) -> Result<Vec<f64>> {
        let h = end - begin;
        let mid = begin + h * 0.5;
        if h <= 0.0 || mid <= begin || mid >= end {
            return Err(BatchNetworkError::Arithmetic(
                "RK stage clock lost resolution",
            ));
        }
        let first = self.rates(start, begin)?;
        let second = self.rates(&increment(start, &first, h * 0.5)?, mid)?;
        let third = self.rates(&increment(start, &second, h * 0.5)?, mid)?;
        let fourth = self.rates(&increment(start, &third, h)?, end)?;
        let weighted = first
            .iter()
            .enumerate()
            .map(|(i, first)| weighted_rate([*first, second[i], third[i], fourth[i]]))
            .collect::<Result<Vec<_>>>()?;
        let result = increment(start, &weighted, h)?;
        self.state(&result, end)?;
        Ok(result)
    }
}

fn weighted_rate(rates: [f64; 4]) -> Result<f64> {
    // Scale before averaging: dividing each positive subnormal rate by 3 or 6
    // can erase every contribution even when their weighted mean is representable.
    let scale = rates.into_iter().fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok(0.0);
    }
    let normalized =
        (rates[0] / scale + 2.0 * (rates[1] / scale) + 2.0 * (rates[2] / scale) + rates[3] / scale)
            / 6.0;
    let rate = scale * normalized;
    if !rate.is_finite() || rate <= 0.0 {
        return Err(BatchNetworkError::Arithmetic(
            "positive weighted RK rate range loss",
        ));
    }
    Ok(rate)
}

fn local_error_estimate(coarse: f64, fine: f64) -> Result<f64> {
    let difference = (coarse - fine).abs();
    let estimate = difference / 15.0;
    if !estimate.is_finite() || (difference > 0.0 && estimate == 0.0) {
        return Err(BatchNetworkError::Arithmetic(
            "positive local-error estimate range loss",
        ));
    }
    Ok(estimate)
}

fn increment(start: &[f64], rates: &[f64], h: f64) -> Result<Vec<f64>> {
    start
        .iter()
        .zip(rates)
        .map(|(before, rate)| {
            let after = rate.mul_add(h, *before);
            if !after.is_finite() || after < *before || (*rate > 0.0 && after == *before) {
                return Err(BatchNetworkError::Arithmetic(
                    "positive extent increment lost numerical range",
                ));
            }
            Ok(after)
        })
        .collect()
}

#[cfg(test)]
#[path = "network/tests.rs"]
mod tests;
