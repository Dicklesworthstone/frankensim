//! One-reaction equilibrium with finite stocks and declared ideal activities.
//!
//! The caller explicitly chooses a fixed dilute-solution amount scale or a
//! pure condensed phase for EVERY participant. The formation sources supply
//! standard chemical potentials, not evidence validating those activity laws.
//! This is an isothermal, fixed-reference-pressure, single reaction closure;
//! it is not a multireaction aqueous speciation or kinetic model.
//!
//! For n_i = n_i0 + nu_i xi, dG/dxi = RT(ln Q - ln K). Within the admitted
//! interval, d(ln Q)/dxi = sum_aqueous(nu_i^2/n_i) >= 0. Bisection therefore
//! finds the constrained minimum without a guessed conversion. Pure-phase
//! inventory limits require an inequality, not a fabricated Q = K equality.
//! References: IUPAC <https://doi.org/10.1351/goldbook.C01023> and USGS
//! <https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-13.htm>.

use std::collections::BTreeMap;

use fs_math::det;
use fs_qty::{Amount, Dimensionless, Dims, QuantitySpec};

use super::{FormationGibbsSource, StandardReactionGibbs, Sum};
use crate::SpeciesId;
use crate::stoichiometry::{ConservedReaction, ReactionBalanceError, Result, representable};

/// Explicit caller assumption, separate from the source's standard-state data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum IdealActivityModel {
    /// a = n / n_ref with a constant solvent basis and activity coefficient one.
    /// `amount_at_unit_activity` is positive [mol] and must match the selected
    /// aqueous standard state. It is NOT inferred from tank volume or molarity.
    /// Different scales are allowed only as explicit caller declarations;
    /// the solver does not authenticate the physical consistency of a mixture.
    IdealDilute {
        /// Amount that would have unit relative activity on the declared basis.
        amount_at_unit_activity: Amount,
    },
    /// Unit activity for the selected pure solid or liquid, including its
    /// one-sided appearance/disappearance limit. Stocks are still finite.
    PureCondensed,
}

/// Complete explicit state and activity assumption for one active participant.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumInventory {
    /// Same identity as the conserved reaction and formation source.
    pub species: SpeciesId,
    /// Finite nonnegative initial amount [mol]; absent phases use zero.
    pub amount: Amount,
    /// Caller-selected law, checked against the selected source's phase.
    pub activity: IdealActivityModel,
}

/// Caller-owned numerical tolerances and bounded work, not model uncertainty.
#[derive(Clone, Copy, Debug)]
pub struct EquilibriumSolveOptions {
    /// Maximum final bracket width [mol of the equation as written].
    pub extent_tolerance: Amount,
    /// Maximum absolute ln(Q/K) for an interior solution.
    pub log_residual_tolerance: Dimensionless,
    /// At most 4096 bisections; exhausted work or floating resolution refuses.
    pub max_iterations: usize,
}

/// Outcome of minimization along ONLY the admitted reaction coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EquilibriumExtentStatus {
    /// Finite activity residual and extent bracket both meet caller tolerances.
    Interior,
    /// Forward motion remains favored at the last admitted pure-phase stock.
    ForwardInventoryLimit,
    /// Reverse motion remains favored at the last admitted pure-phase stock.
    ReverseInventoryLimit,
    /// Neither direction has the required stocks; no activity equality claimed.
    InventoryBlocked,
    /// All activities are constant and Delta_r G is zero over the full interval.
    NonUnique,
}

/// Numerical equilibrium result, retaining source evidence and assumptions.
#[derive(Clone, Debug)]
pub struct SingleReactionEquilibrium {
    /// Original standard Gibbs energy, common reference state and all receipts.
    pub standard: StandardReactionGibbs,
    /// Canonically ordered initial stocks and explicit activity models.
    pub initial: Vec<EquilibriumInventory>,
    /// Canonically ordered final stocks, evaluated by FMA without clipping.
    pub final_inventory: Vec<(SpeciesId, Amount)>,
    /// Signed extent [mol]; negative reverses the written equation.
    pub extent: Amount,
    /// Final numerical bracket; full feasible interval for a nonunique result.
    /// This is not an interval-arithmetic or thermodynamic uncertainty certificate.
    pub extent_bracket: [Amount; 2],
    /// Interior root, one-sided inventory constraint, or degenerate state.
    pub status: EquilibriumExtentStatus,
    /// Finite ln(Q/K) at the returned state, or None for an inventory-blocked
    /// state with undefined zero activities. Boundary values need not be zero.
    pub log_quotient_over_equilibrium: Option<Dimensionless>,
    /// Number of bisection iterations used, excluding boundary checks.
    pub iterations: usize,
}

impl ConservedReaction {
    /// Minimize ideal-model free energy along this one conserved reaction.
    ///
    /// All active species and source packs are required, including zero stocks
    /// for absent products. Source-owned phase flags must match the explicit
    /// law: aqueous for `IdealDilute`, pure solid/liquid for `PureCondensed`.
    /// Gases, solid solutions and concentration-dependent activity coefficients
    /// need other models. A fixed solvent basis and unit solvent activity are
    /// approximations the caller must justify; neither is plant calibration.
    ///
    /// Signed extent bounds are rounded inward using FMA stock checks. A
    /// boundary result may therefore retain a sub-ulp-scale remainder; it never
    /// changes a negative stock to zero. A dilute-species limiting endpoint
    /// cannot masquerade as equilibrium when an interior root is unresolved.
    pub fn equilibrate_ideal(
        &self,
        sources: &[FormationGibbsSource<'_>],
        inventory: &[EquilibriumInventory],
        options: EquilibriumSolveOptions,
    ) -> Result<SingleReactionEquilibrium> {
        if !options.extent_tolerance.value().is_finite()
            || options.extent_tolerance.value() <= 0.0
            || !options.log_residual_tolerance.value().is_finite()
            || options.log_residual_tolerance.value() <= 0.0
            || options.max_iterations == 0
            || options.max_iterations > 4096
        {
            return Err(ReactionBalanceError::Invalid(
                "invalid equilibrium solve budget",
            ));
        }
        let standard = self.standard_gibbs_change(sources)?;
        if inventory.len() != self.terms().len() {
            return Err(ReactionBalanceError::Invalid(
                "supply every equilibrium inventory",
            ));
        }
        let mut by_species = BTreeMap::new();
        for entry in inventory {
            if !entry.amount.value().is_finite() || entry.amount.value() < 0.0 {
                return Err(ReactionBalanceError::Invalid(
                    "invalid equilibrium inventory",
                ));
            }
            if by_species.insert(&entry.species, entry).is_some() {
                return Err(ReactionBalanceError::Invalid(
                    "duplicate equilibrium inventory",
                ));
            }
        }
        let mut initial = Vec::with_capacity(inventory.len());
        let mut limits = Vec::with_capacity(inventory.len());
        let mut lower = f64::NEG_INFINITY;
        let mut upper = f64::INFINITY;
        for evidence in standard.terms() {
            let term = &evidence.term;
            let entry =
                by_species
                    .get(&term.species)
                    .ok_or_else(|| ReactionBalanceError::Species {
                        species: term.species.clone(),
                        reason: "missing equilibrium inventory",
                    })?;
            let source = sources
                .iter()
                .find(|source| *source.species == term.species)
                .ok_or(ReactionBalanceError::Invalid("missing equilibrium source"))?;
            let claim = source
                .pack
                .claims()
                .claim(evidence.answer.receipt.selected)
                .ok_or(ReactionBalanceError::Invalid(
                    "missing selected equilibrium claim",
                ))?;
            let owns = |axis: &str| {
                claim.validity.bound(axis) == Some((1.0, 1.0))
                    && claim.validity.axis_quantities().get(axis)
                        == Some(&QuantitySpec::dimensional(Dims::NONE))
            };
            match entry.activity {
                IdealActivityModel::IdealDilute {
                    amount_at_unit_activity,
                } => {
                    if !owns("source-phase-aqueous")
                        || !amount_at_unit_activity.value().is_finite()
                        || amount_at_unit_activity.value() <= 0.0
                    {
                        return Err(ReactionBalanceError::Invalid(
                            "ideal dilute law needs an aqueous source and positive amount scale",
                        ));
                    }
                }
                IdealActivityModel::PureCondensed => {
                    if !owns("source-pure-species")
                        || !(owns("source-phase-solid") || owns("source-phase-liquid"))
                    {
                        return Err(ReactionBalanceError::Invalid(
                            "unit condensed activity needs a pure solid or liquid source",
                        ));
                    }
                }
            }
            let coefficient = term.coefficient as f64;
            let mut limit = entry.amount.value() / coefficient.abs();
            representable(
                limit,
                entry.amount.value() > 0.0,
                "equilibrium extent bound",
            )?;
            if limit.mul_add(coefficient.abs(), -entry.amount.value()) > 0.0 {
                limit = limit.next_down();
            }
            representable(
                limit,
                entry.amount.value() > 0.0,
                "equilibrium extent bound",
            )?;
            if coefficient < 0.0 {
                upper = upper.min(limit);
            } else {
                lower = lower.max(-limit);
            }
            limits.push(limit);
            initial.push((**entry).clone());
        }
        let state = State {
            standard: &standard,
            initial: &initial,
            limits: &limits,
        };
        let mut bracket = [lower, upper];
        let (extent, status, residual, iterations) = if lower == upper {
            (0.0, EquilibriumExtentStatus::InventoryBlocked, None, 0)
        } else if initial
            .iter()
            .all(|entry| entry.activity == IdealActivityModel::PureCondensed)
            && standard.log_equilibrium_constant().value() == 0.0
        {
            (0.0, EquilibriumExtentStatus::NonUnique, Some(0.0), 0)
        } else {
            let at_lower = state.endpoint(lower, false)?;
            let at_upper = state.endpoint(upper, true)?;
            if let Residual::Finite(value) = at_lower
                && value >= 0.0
            {
                bracket = [lower, lower];
                (
                    lower,
                    EquilibriumExtentStatus::ReverseInventoryLimit,
                    Some(value),
                    0,
                )
            } else if let Residual::Finite(value) = at_upper
                && value <= 0.0
            {
                bracket = [upper, upper];
                (
                    upper,
                    EquilibriumExtentStatus::ForwardInventoryLimit,
                    Some(value),
                    0,
                )
            } else {
                let (extent, value, count) = state.bisect(&mut bracket, options)?;
                (
                    extent,
                    EquilibriumExtentStatus::Interior,
                    Some(value),
                    count,
                )
            }
        };
        let final_inventory = state.inventory(extent, true)?;
        Ok(SingleReactionEquilibrium {
            standard,
            initial,
            final_inventory,
            extent: Amount::new(extent),
            extent_bracket: bracket.map(Amount::new),
            status,
            log_quotient_over_equilibrium: residual.map(Dimensionless::new),
            iterations,
        })
    }
}

#[derive(Clone, Copy)]
enum Residual {
    NegativeInfinity,
    Finite(f64),
    PositiveInfinity,
    Undefined,
}
impl Residual {
    fn finite(self) -> Option<f64> {
        if let Self::Finite(value) = self {
            Some(value)
        } else {
            None
        }
    }
}

struct State<'a> {
    standard: &'a StandardReactionGibbs,
    initial: &'a [EquilibriumInventory],
    limits: &'a [f64],
}

impl State<'_> {
    fn inventory(&self, extent: f64, check_resolution: bool) -> Result<Vec<(SpeciesId, Amount)>> {
        self.standard
            .terms()
            .iter()
            .zip(self.initial)
            .map(|(source, entry)| {
                let amount = (source.term.coefficient as f64).mul_add(extent, entry.amount.value());
                representable(amount, false, "equilibrium final inventory")?;
                if amount < 0.0 {
                    return Err(ReactionBalanceError::Arithmetic(
                        "negative equilibrium stock",
                    ));
                }
                if check_resolution && extent != 0.0 && amount == entry.amount.value() {
                    return Err(ReactionBalanceError::Arithmetic(
                        "equilibrium change below stock resolution",
                    ));
                }
                Ok((entry.species.clone(), Amount::new(amount)))
            })
            .collect()
    }

    fn residual(&self, extent: f64) -> Result<Residual> {
        let inventory = self.inventory(extent, false)?;
        let mut negative_infinity = false;
        let mut positive_infinity = false;
        let mut sum = Sum::default();
        sum.add(-self.standard.log_equilibrium_constant().value())?;
        for ((source, entry), (_, amount)) in self
            .standard
            .terms()
            .iter()
            .zip(self.initial)
            .zip(inventory)
        {
            if let IdealActivityModel::IdealDilute {
                amount_at_unit_activity,
            } = entry.activity
            {
                let coefficient = source.term.coefficient as f64;
                if amount.value() == 0.0 {
                    if coefficient > 0.0 {
                        negative_infinity = true;
                    } else {
                        positive_infinity = true;
                    }
                } else {
                    // Separate logs avoid overflow/underflow of n / n_ref.
                    let logarithm =
                        det::ln(amount.value()) - det::ln(amount_at_unit_activity.value());
                    sum.add(coefficient * logarithm)?;
                }
            }
        }
        Ok(match (negative_infinity, positive_infinity) {
            (true, true) => Residual::Undefined,
            (true, false) => Residual::NegativeInfinity,
            (false, true) => Residual::PositiveInfinity,
            (false, false) => Residual::Finite(sum.finish()?),
        })
    }

    fn endpoint(&self, extent: f64, forward: bool) -> Result<Residual> {
        // A finite binary64 endpoint rounds inward. Do not mistake a tiny
        // remaining aqueous stock for a genuine phase-exhaustion equilibrium.
        for ((source, entry), limit) in self
            .standard
            .terms()
            .iter()
            .zip(self.initial)
            .zip(self.limits)
        {
            if matches!(entry.activity, IdealActivityModel::IdealDilute { .. })
                && (source.term.coefficient < 0) == forward
                && *limit == extent.abs()
            {
                return Ok(if forward {
                    Residual::PositiveInfinity
                } else {
                    Residual::NegativeInfinity
                });
            }
        }
        self.residual(extent)
    }

    fn bisect(
        &self,
        bracket: &mut [f64; 2],
        options: EquilibriumSolveOptions,
    ) -> Result<(f64, f64, usize)> {
        for iteration in 1..=options.max_iterations {
            // Overflow-safe, including opposite signs and subnormal stocks.
            let midpoint = bracket[0].midpoint(bracket[1]);
            if midpoint <= bracket[0] || midpoint >= bracket[1] {
                return Err(ReactionBalanceError::Arithmetic(
                    "equilibrium root below extent resolution",
                ));
            }
            let residual =
                self.residual(midpoint)?
                    .finite()
                    .ok_or(ReactionBalanceError::Arithmetic(
                        "undefined interior equilibrium activity",
                    ))?;
            if residual == 0.0 {
                *bracket = [midpoint, midpoint];
                return Ok((midpoint, residual, iteration));
            }
            if residual < 0.0 {
                bracket[0] = midpoint;
            } else {
                bracket[1] = midpoint;
            }
            if residual.abs() <= options.log_residual_tolerance.value()
                && bracket[1] - bracket[0] <= options.extent_tolerance.value()
            {
                return Ok((midpoint, residual, iteration));
            }
        }
        Err(ReactionBalanceError::Invalid(
            "equilibrium iteration budget exhausted",
        ))
    }
}
