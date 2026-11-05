//! Grind-to-spec by bounded inversion of the original batch population balance.
//!
//! A target includes a named source class and every no-larger source class.
//! Its mass cannot escape to coarser classes in the admitted daughter matrix,
//! so cumulative target mass is nondecreasing for this constant-rate model.
//! Each trial calls the original `prepare` and `transform` owners; no new
//! kinetic law, arbitrary sieve interpolation or plant scale-up is introduced.
//!
//! A successful nonzero solve retains an evaluated below-target/meeting-target
//! time bracket no wider than the requested tolerance. Comparisons include
//! twice the original Poisson tail plus a caller-declared rounding allowance.
//! These are nominal binary64 bounds, not an interval certificate or allowance
//! for experimental/model uncertainty. An ambiguous comparison refuses rather
//! than declaring a threshold crossed. No exact earliest time is claimed.

use super::{
    GrindingKineticsError, KineticGrindingPass, SelectionBreakageModel, UniformizationBudget,
};
use crate::comminution::SpecificGrindingEnergy;
use crate::separation::ParticleStream;
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, Energy, Mass, Power, QuantitySpec, Time};

/// A cumulative target on the original source discretization.
#[derive(Clone, Debug, PartialEq)]
pub struct GrindTarget {
    /// Include this original class and every class of equal or smaller diameter.
    /// A representative/sink label is not an independently resolved sieve cut.
    pub coarsest_class: String,
    /// Required fraction of total product dry solids, in (0,1].
    pub minimum_fraction: Dimensionless,
}

/// Explicit search window, numerical allowances and total forward-call budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrindTargetSettings {
    /// Positive upper contact time; the search starts at zero and never expands.
    pub maximum_contact_time: Time,
    /// Positive maximum width of a successful nonzero time bracket.
    pub time_tolerance: Time,
    /// Nonnegative absolute fraction allowance for source/arithmetic rounding.
    /// This caller assumption is not a computed error certificate or model error.
    pub fraction_roundoff_allowance: Dimensionless,
    /// Total original `prepare`/`transform` trials, including endpoints; at least 2.
    pub maximum_forward_evaluations: usize,
    /// Unmodified numerical-work and tail limits for every original preparation.
    pub forward_budget: UniformizationBudget,
}

/// An actual evaluated point, with its nominal numerical comparison allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrindTargetPoint {
    /// Time passed to the original solver.
    pub contact_time: Time,
    /// Original product mass in the target classes divided by total product mass.
    pub mass_fraction: Dimensionless,
    /// Original preparation's nominal Poisson-tail probability bound.
    pub poisson_tail_bound: Dimensionless,
    /// Twice the tail bound plus the declared fraction rounding allowance.
    pub nominal_fraction_error: Dimensionless,
}
impl GrindTargetPoint {
    /// Nominal lower comparison value; deliberately not clipped to [0,1].
    pub fn lower_fraction(self) -> f64 {
        self.mass_fraction.value() - self.nominal_fraction_error.value()
    }
    /// Nominal upper comparison value; deliberately not clipped to [0,1].
    pub fn upper_fraction(self) -> f64 {
        self.mass_fraction.value() + self.nominal_fraction_error.value()
    }
}

/// Retained original evaluations on opposite sides of the requested target.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrindTargetBracket {
    /// Its upper fraction comparison is strictly below the requested fraction.
    pub below: GrindTargetPoint,
    /// Its lower fraction comparison meets the requested fraction.
    pub meeting: GrindTargetPoint,
}

/// Invalid input, original forward refusal or unresolved bounded inverse.
#[derive(Clone, Debug, PartialEq)]
pub enum GrindTargetError {
    /// Invalid target, settings, energy inputs or missing source class.
    Invalid(&'static str),
    /// The original preparation or conservative stream owner refused a trial.
    Forward {
        /// Actual attempted contact time.
        contact_time: Time,
        /// Unmodified original refusal.
        error: GrindingKineticsError,
    },
    /// The target remains above the upper comparison at the admitted time limit.
    Unreachable {
        /// Requested mass fraction.
        target: Dimensionless,
        /// Actual result at the caller's maximum contact time.
        at_limit: GrindTargetPoint,
    },
    /// The total forward-call budget ended before the time tolerance was met.
    EvaluationBudget {
        /// Retained below/meeting bracket; no successful time is returned.
        bracket: GrindTargetBracket,
        /// Number of original forward evaluations performed.
        evaluations: usize,
    },
    /// A threshold or time refinement cannot be resolved at the declared accuracy.
    NumericalResolution {
        /// Concrete comparison or arithmetic failure.
        reason: &'static str,
        /// Last actual evaluated point.
        point: GrindTargetPoint,
        /// Retained enclosing bracket, if both endpoints have been established.
        bracket: Option<GrindTargetBracket>,
    },
    /// Caller cancellation observed between or within original forward calls.
    Cancelled,
}
impl core::fmt::Display for GrindTargetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "grinding target: {self:?}")
    }
}
impl std::error::Error for GrindTargetError {}

/// Accepted original product and pass, with the actual bounded inverse result.
#[derive(Clone, Debug, PartialEq)]
pub struct GrindTargetSolution {
    /// Original prepared pass at zero or the meeting endpoint; retains its sources.
    pub pass: KineticGrindingPass,
    /// Original transformed feed, including every element and carrier water.
    pub product: ParticleStream,
    /// Requested source-bin target.
    pub target: GrindTarget,
    /// Exact included source class identities in their canonical ID order.
    pub included_classes: Vec<String>,
    /// Accepted original point, including its nominal tail/rounding comparison.
    pub point: GrindTargetPoint,
    /// None only when the initial feed already meets the nominal target at zero.
    pub bracket: Option<GrindTargetBracket>,
    /// Total original forward calls performed, including the initial zero trial.
    pub evaluations: usize,
    /// Original search and forward budgets, unchanged.
    pub settings: GrindTargetSettings,
}

/// Energy accounting from explicitly supplied constant net grinding power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatchGrindingEnergy {
    /// Caller-supplied net power; neither inferred from kinetics nor motor rating.
    pub net_power: Power,
    /// Caller-supplied batch dry mass; the PSD flow scale is not a batch inventory.
    pub dry_batch_mass: Mass,
    /// Energy at the retained below endpoint, or zero for an already-met target.
    pub below_energy: Energy,
    /// Energy at the accepted meeting time, `net_power * contact_time`.
    pub meeting_energy: Energy,
    /// Meeting-time energy divided by the explicit dry batch mass, in J/kg.
    pub specific_energy: SpecificGrindingEnergy,
}

impl SelectionBreakageModel {
    /// Find a target-crossing time bracket on the admitted source bins.
    ///
    /// Both source conditions and feed layout pass through the original owners
    /// at every trial. A target already met at zero needs only that trial; otherwise
    /// the maximum time must also be admitted. Unreachable or numerically ambiguous
    /// targets, exhausted work, cancellation and any source refusal return errors.
    pub fn time_to_target(
        &self,
        feed: &ParticleStream,
        target: GrindTarget,
        conditions: &QueryPoint,
        settings: GrindTargetSettings,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GrindTargetSolution, GrindTargetError> {
        settings.validate(&target)?;
        if conditions.axes().contains_key(super::CONTACT_TIME_AXIS) {
            return Err(GrindTargetError::Invalid(
                "the inverse owns grinding-contact-time",
            ));
        }
        let classes = self.daughters().classes();
        let limit = classes
            .iter()
            .find(|c| c.id == target.coarsest_class)
            .ok_or(GrindTargetError::Invalid("unknown target source class"))?
            .diameter;
        let included: Vec<_> = classes.iter().map(|c| c.diameter <= limit).collect();
        // Check the no-escape property on the admitted matrix, not ID ordering.
        for (j, column) in self.daughters().columns().iter().enumerate() {
            if !included[j] || self.rates()[j].value() == 0.0 {
                continue;
            }
            for transfer in &column.transfers {
                let i = classes
                    .binary_search_by(|c| c.id.cmp(&transfer.output))
                    .map_err(|_| GrindTargetError::Invalid("unknown daughter class"))?;
                if transfer.mass_fraction.value() > 0.0 && !included[i] {
                    return Err(GrindTargetError::Invalid(
                        "target mass can escape to a coarser class",
                    ));
                }
            }
        }
        let mut evaluate = |time: Time| {
            if cancelled() {
                return Err(GrindTargetError::Cancelled);
            }
            let map_error = |error| match error {
                GrindingKineticsError::Cancelled => GrindTargetError::Cancelled,
                error => GrindTargetError::Forward {
                    contact_time: time,
                    error,
                },
            };
            let pass = self
                .prepare(time, conditions, settings.forward_budget, &mut cancelled)
                .map_err(map_error)?;
            if cancelled() {
                return Err(GrindTargetError::Cancelled);
            }
            let product = pass.transform(feed).map_err(map_error)?;
            if cancelled() {
                return Err(GrindTargetError::Cancelled);
            }
            let mut total = super::Sum::default();
            let mut selected = super::Sum::default();
            for (class, &keep) in product.classes().iter().zip(&included) {
                total.add(class.dry_solids().value()).map_err(map_error)?;
                if keep {
                    selected
                        .add(class.dry_solids().value())
                        .map_err(map_error)?;
                }
            }
            let mass = total.value().map_err(map_error)?;
            if mass <= 0.0 {
                return Err(GrindTargetError::Invalid("target needs positive dry feed"));
            }
            let selected_mass = selected.value().map_err(map_error)?;
            let fraction = selected_mass / mass;
            if !fraction.is_finite()
                || !(0.0..=1.0).contains(&fraction)
                || (selected_mass > 0.0 && fraction == 0.0)
            {
                return Err(GrindTargetError::Invalid(
                    "unrepresentable cumulative product fraction",
                ));
            }
            let point = GrindTargetPoint {
                contact_time: time,
                mass_fraction: Dimensionless::new(fraction),
                poisson_tail_bound: pass.poisson_tail_bound,
                nominal_fraction_error: Dimensionless::new(
                    2.0 * pass.poisson_tail_bound.value()
                        + settings.fraction_roundoff_allowance.value(),
                ),
            };
            Ok((pass, product, point))
        };
        let finish = |(pass, product, point), bracket, evaluations| GrindTargetSolution {
            pass,
            product,
            point,
            bracket,
            evaluations,
            settings,
            target: target.clone(),
            included_classes: classes
                .iter()
                .zip(&included)
                .filter(|(_, keep)| **keep)
                .map(|(c, _)| c.id.clone())
                .collect(),
        };
        let requested = target.minimum_fraction.value();
        let unresolved = |point, bracket, reason| GrindTargetError::NumericalResolution {
            reason,
            point,
            bracket,
        };
        let initial = evaluate(Time::new(0.0))?;
        if initial.2.lower_fraction() >= requested {
            return Ok(finish(initial, None, 1));
        }
        if initial.2.upper_fraction() >= requested {
            return Err(unresolved(
                initial.2,
                None,
                "initial target comparison is numerically ambiguous",
            ));
        }
        let mut upper = evaluate(settings.maximum_contact_time)?;
        let mut evaluations = 2;
        if upper.2.upper_fraction() < initial.2.lower_fraction() {
            return Err(unresolved(
                upper.2,
                None,
                "maximum-time fraction violates the initial monotone comparison",
            ));
        }
        if upper.2.upper_fraction() < requested {
            return Err(GrindTargetError::Unreachable {
                target: target.minimum_fraction,
                at_limit: upper.2,
            });
        }
        if upper.2.lower_fraction() < requested {
            return Err(unresolved(
                upper.2,
                None,
                "maximum-time target comparison is numerically ambiguous",
            ));
        }
        let mut bracket = GrindTargetBracket {
            below: initial.2,
            meeting: upper.2,
        };
        loop {
            let lo = bracket.below.contact_time.value();
            let hi = bracket.meeting.contact_time.value();
            if hi - lo <= settings.time_tolerance.value() {
                return Ok(finish(upper, Some(bracket), evaluations));
            }
            if evaluations >= settings.maximum_forward_evaluations {
                return Err(GrindTargetError::EvaluationBudget {
                    bracket,
                    evaluations,
                });
            }
            let middle = lo + (hi - lo) * 0.5;
            if middle <= lo || middle >= hi {
                return Err(unresolved(
                    bracket.meeting,
                    Some(bracket),
                    "no representable interior time remains",
                ));
            }
            let trial = evaluate(Time::new(middle))?;
            evaluations += 1;
            let point = trial.2;
            if point.upper_fraction() < bracket.below.lower_fraction()
                || point.lower_fraction() > bracket.meeting.upper_fraction()
            {
                return Err(unresolved(
                    point,
                    Some(bracket),
                    "trial violates monotone endpoint comparisons",
                ));
            }
            if point.upper_fraction() < requested {
                bracket.below = point;
            } else if point.lower_fraction() >= requested {
                bracket.meeting = point;
                upper = trial;
            } else {
                return Err(unresolved(
                    point,
                    Some(bracket),
                    "trial target comparison is numerically ambiguous",
                ));
            }
        }
    }
}

impl GrindTargetSolution {
    /// Quote energy from declared constant net power and batch mass at this time.
    /// If the admitted point declares `dry-charge-mass`, the supplied mass must
    /// match it exactly. This accounting neither predicts power nor changes the
    /// source calibration; it excludes idle time, auxiliaries and drive losses.
    pub fn energy_at_constant_net_power(
        &self,
        net_power: Power,
        dry_batch_mass: Mass,
    ) -> Result<BatchGrindingEnergy, GrindTargetError> {
        let p = net_power.value();
        let m = dry_batch_mass.value();
        if !p.is_finite() || p <= 0.0 || !m.is_finite() || m <= 0.0 {
            return Err(GrindTargetError::Invalid(
                "net power and batch mass must be finite and positive",
            ));
        }
        if let Some(&source_mass) = self.pass.point().axes().get("dry-charge-mass") {
            if source_mass != m
                || self.pass.point().axis_quantities().get("dry-charge-mass")
                    != Some(&QuantitySpec::dimensional(Mass::DIMS))
            {
                return Err(GrindTargetError::Invalid(
                    "batch mass conflicts with admitted dry-charge-mass",
                ));
            }
        }
        let energy = |time: Time| {
            let value = p * time.value();
            if !value.is_finite() || (time.value() > 0.0 && value == 0.0) {
                Err(GrindTargetError::Invalid(
                    "net energy loses numerical range",
                ))
            } else {
                Ok(Energy::new(value))
            }
        };
        let meeting_energy = energy(self.pass.elapsed())?;
        let below_energy = energy(
            self.bracket
                .map_or(Time::new(0.0), |b| b.below.contact_time),
        )?;
        let specific = meeting_energy.value() / m;
        if !specific.is_finite() || (meeting_energy.value() > 0.0 && specific == 0.0) {
            return Err(GrindTargetError::Invalid(
                "specific energy loses numerical range",
            ));
        }
        Ok(BatchGrindingEnergy {
            net_power,
            dry_batch_mass,
            below_energy,
            meeting_energy,
            specific_energy: SpecificGrindingEnergy::new(specific),
        })
    }
}

impl GrindTargetSettings {
    fn validate(self, target: &GrindTarget) -> Result<(), GrindTargetError> {
        let fraction = target.minimum_fraction.value();
        let max_time = self.maximum_contact_time.value();
        let tolerance = self.time_tolerance.value();
        let allowance = self.fraction_roundoff_allowance.value();
        if !fraction.is_finite()
            || !(0.0 < fraction && fraction <= 1.0)
            || !max_time.is_finite()
            || max_time <= 0.0
            || !tolerance.is_finite()
            || tolerance <= 0.0
            || tolerance > max_time
            || !allowance.is_finite()
            || !(0.0..1.0).contains(&allowance)
            || self.maximum_forward_evaluations < 2
        {
            return Err(GrindTargetError::Invalid(
                "invalid fraction, time window, allowance or forward-call budget",
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "target/tests.rs"]
mod tests;
