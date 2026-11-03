//! Source-bound, finite-batch progress for an explicit two-reactant rate law.
//!
//! For the conserved equation `a A + b B -> products`, this owner solves
//! `dξ/dt = k n_A n_B / V`, with `n_A = n_A0 - a ξ` and
//! `n_B = n_B0 - b ξ`, at fixed volume and temperature. Both empirical orders
//! are one; they are declared by this law, not inferred from `a` and `b`.
//! The equal- and unequal-capacity analytic solutions deplete both stocks
//! together. A recovery curve clipped to a reagent ceiling is not this law.
//!
//! For a mineral, inventory divided by reference batch volume is an empirical
//! availability variable. It is not pure-solid activity, exposed surface area,
//! an elementary reaction mechanism, or a prediction of mineral kinetics.
//! Original source support is required; admission does not establish calibration.
//! The empirical rate-law convention is described by IUPAC:
//! <https://doi.org/10.1351/goldbook.R05141> and
//! <https://doi.org/10.1351/goldbook.O04322>. Neither supplies an ore-specific k.
//!
//! Scope: one irreversible equation, exactly two consumed species, constant
//! k/volume/temperature and other source conditions, no inflow or outflow.
//! Equilibria, speciation, surface evolution and coupled reactions remain
//! separate models. Element/charge conservation comes from the original
//! [`ConservedReaction`]; binary64 progression is not an interval certificate.

/// Simultaneous finite-batch competition for one explicitly shared reagent.
#[path = "reaction_batch/network.rs"]
pub mod network;

/// Spherical mineral dissolution with film, surface and product-layer resistance.
#[path = "reaction_batch/shrinking_core.rs"]
pub mod shrinking_core;

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::det;
use fs_qty::{Amount, AmountConcentration, Qty, QuantitySpec, Temperature, Time, Volume};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{ConservedReaction, ReactionAmount, ReactionBalanceError},
};

/// Empirical bilinear rate coefficient, in m³/(mol s).
pub type BilinearRateConstant = Qty<3, 0, -1, 0, 0, -1>;
/// Law identifier; version 1 fixes both concentration orders at one.
pub const BILINEAR_BATCH_LAW: &str = "bilinear-batch-reaction";
/// The sole dimensioned model parameter.
pub const RATE_PARAMETER: &str = "rate-constant";
/// Seconds since this batch advance started; the whole interval is admitted.
pub const ELAPSED_TIME_AXIS: &str = "elapsed-time";
/// Fixed reference volume in cubic metres.
pub const VOLUME_AXIS: &str = "batch-volume";
/// Fixed absolute temperature in kelvin.
pub const TEMPERATURE_AXIS: &str = "temperature";
/// Actual A inventory divided by reference volume, in mol/m³.
pub const REACTANT_A_CONCENTRATION_AXIS: &str = "reactant-a-concentration";
/// Actual B inventory divided by reference volume, in mol/m³.
pub const REACTANT_B_CONCENTRATION_AXIS: &str = "reactant-b-concentration";

const OWNED_AXES: [&str; 5] = [
    ELAPSED_TIME_AXIS,
    VOLUME_AXIS,
    TEMPERATURE_AXIS,
    REACTANT_A_CONCENTRATION_AXIS,
    REACTANT_B_CONCENTRATION_AXIS,
];

/// Invalid model/state, unsupported source conditions, or numerical range loss.
#[derive(Clone, Debug, PartialEq)]
pub enum BatchReactionError {
    /// Original source-card validation or typed query admission failed.
    Source(MatDbError),
    /// The existing conserved-inventory owner refused the supplied state.
    Reaction(ReactionBalanceError),
    /// Invalid law, reactant mapping, dimensions, or physical input.
    Invalid(&'static str),
    /// A required initial/final coordinate is missing, mistyped or unsupported.
    OutsideValidity {
        /// Which endpoint refused; both must be admitted before returning.
        endpoint: &'static str,
    },
    /// Positive progress or remaining stock cannot be represented faithfully.
    Arithmetic(&'static str),
}
impl core::fmt::Display for BatchReactionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bilinear batch reaction: {self:?}")
    }
}
impl std::error::Error for BatchReactionError {}
impl From<MatDbError> for BatchReactionError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<ReactionBalanceError> for BatchReactionError {
    fn from(error: ReactionBalanceError) -> Self {
        Self::Reaction(error)
    }
}
type Result<T> = core::result::Result<T, BatchReactionError>;

/// Immutable original reaction, explicit A/B mapping and original source card.
#[derive(Clone, Debug)]
pub struct BilinearBatchReaction {
    reaction: ConservedReaction,
    source: ConstitutiveModelCard,
    reactants: [SpeciesId; 2],
    coefficients: [i128; 2],
    rate: BilinearRateConstant,
}

/// Actual stock account for one explicitly mapped reactant.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactantAccount {
    /// Exact species identity from the original equation.
    pub species: SpeciesId,
    /// Positive consumption coefficient in the equation as written.
    pub coefficient: i128,
    /// Admitted initial inventory, in mol of this species.
    pub initial: Amount,
    /// Gross amount consumed, from original stoichiometry and returned extent.
    pub consumed: Amount,
    /// Actual final inventory returned by the conserved-inventory owner.
    pub remaining: Amount,
}

/// A completed, immutable batch calculation with its original source evidence.
#[derive(Clone, Debug)]
pub struct BatchReactionAdvance {
    model: BilinearBatchReaction,
    /// Complete canonical inventory before the advance, including products.
    pub initial_inventory: Vec<(SpeciesId, Amount)>,
    /// Complete canonical inventory after the advance, including products.
    pub final_inventory: Vec<(SpeciesId, Amount)>,
    /// Mol of the original integer equation as written.
    pub extent: Amount,
    /// Conservatively admitted stoichiometric availability bound, not a rate.
    pub limiting_extent: Amount,
    /// Gross consumed/formed amounts; coefficient sign identifies direction.
    pub amounts: Vec<ReactionAmount>,
    /// Explicit A then B stock accounts.
    pub reactants: [ReactantAccount; 2],
    /// Constant reference volume actually used.
    pub volume: Volume,
    /// Actual finite contact time.
    pub elapsed: Time,
    /// Constant absolute temperature actually admitted.
    pub temperature: Temperature,
    /// Time zero and actual initial reactant concentrations.
    pub initial_point: QueryPoint,
    /// Contact time and actual final reactant concentrations.
    pub final_point: QueryPoint,
}
impl BatchReactionAdvance {
    /// Original source hashes/provenance, explicit mapping, and exact equation.
    /// The equation's conservation certificate remains at `model().reaction()`.
    pub fn model(&self) -> &BilinearBatchReaction {
        &self.model
    }
}

impl BilinearBatchReaction {
    /// Admit one conserved equation with exactly two distinct consumed species.
    /// The mapping binds source concentration axes to actual species; no role
    /// or kinetic order is inferred from names or stoichiometric coefficients.
    /// Requires version 1, zero hidden state, one nonnegative typed k, original
    /// source hashes/provenance, and finite typed bounds for every source axis.
    pub fn from_card(
        reaction: ConservedReaction,
        reactant_a: SpeciesId,
        reactant_b: SpeciesId,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        source.validate()?;
        if source.law.0 != BILINEAR_BATCH_LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.parameters.len() != 1
            || source.sources.is_empty()
            || source.provenance.source.trim().is_empty()
            || source.provenance.license.trim().is_empty()
        {
            return Err(BatchReactionError::Invalid(
                "requires version 1 bilinear law, zero state, one parameter and original sources",
            ));
        }
        let rate = source
            .parameters
            .get(RATE_PARAMETER)
            .ok_or(BatchReactionError::Invalid(
                "missing bilinear rate-constant",
            ))?;
        if rate.dims != BilinearRateConstant::DIMS || !rate.value.is_finite() || rate.value < 0.0 {
            return Err(BatchReactionError::Invalid(
                "rate must be finite, nonnegative and in m3/(mol s)",
            ));
        }
        if source.validity.bounds().iter().any(|(axis, (lo, hi))| {
            !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
                || !source.validity.axis_quantities().contains_key(axis)
        }) {
            return Err(BatchReactionError::Invalid(
                "every source axis needs finite bounds and an explicit quantity",
            ));
        }
        for (axis, dims) in [
            (ELAPSED_TIME_AXIS, Time::DIMS),
            (VOLUME_AXIS, Volume::DIMS),
            (TEMPERATURE_AXIS, Temperature::DIMS),
            (REACTANT_A_CONCENTRATION_AXIS, AmountConcentration::DIMS),
            (REACTANT_B_CONCENTRATION_AXIS, AmountConcentration::DIMS),
        ] {
            let (lo, hi) = source
                .validity
                .bound(axis)
                .ok_or(BatchReactionError::Invalid(
                    "declare all five owned applicability axes",
                ))?;
            if source.validity.axis_quantities().get(axis) != Some(&QuantitySpec::dimensional(dims))
                || lo < 0.0
                || (axis == ELAPSED_TIME_AXIS && (lo != 0.0 || hi <= 0.0))
                || ((axis == VOLUME_AXIS || axis == TEMPERATURE_AXIS) && lo <= 0.0)
            {
                return Err(BatchReactionError::Invalid(
                    "owned axes require canonical quantities, time from zero, positive volume/kelvin and nonnegative concentrations",
                ));
            }
        }
        if reactant_a == reactant_b
            || reaction
                .terms()
                .iter()
                .filter(|t| t.coefficient < 0)
                .count()
                != 2
        {
            return Err(BatchReactionError::Invalid(
                "requires exactly two distinct consumed species",
            ));
        }
        let reactants = [reactant_a, reactant_b];
        let mut coefficients = [0; 2];
        for (i, species) in reactants.iter().enumerate() {
            coefficients[i] = -reaction
                .terms()
                .iter()
                .find(|term| &term.species == species && term.coefficient < 0)
                .ok_or(BatchReactionError::Invalid(
                    "each mapped species must be a reactant in the original equation",
                ))?
                .coefficient;
        }
        let rate = BilinearRateConstant::new(rate.value);
        Ok(Self {
            reaction,
            source,
            reactants,
            coefficients,
            rate,
        })
    }

    /// Exact original conserved equation, including its conservation evidence.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Entire original source declaration; no provenance is synthesized.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Explicit species identities bound to source concentration axes A and B.
    pub fn reactants(&self) -> &[SpeciesId; 2] {
        &self.reactants
    }
    /// Declared empirical coefficient, without any inferred Arrhenius behavior.
    pub const fn rate(&self) -> BilinearRateConstant {
        self.rate
    }

    /// Advance a complete molar inventory without mutating it.
    ///
    /// All active species, including initially absent products, must be supplied
    /// exactly once. Supply other source coordinates through `conditions`; the
    /// five owned axes must not be supplied there. Temperature and all other
    /// conditions are fixed; no pH, activity or surface area is inferred.
    ///
    /// Both endpoints are checked against the original rectangular source
    /// domain. Time increases and both reactant inventories decrease
    /// monotonically, so this admits the entire path under the declared law.
    /// Even zero-progress cases require valid source conditions. An extent that
    /// erases positive finite-time remaining stock refuses instead of rounding
    /// to complete exhaustion or clipping to a stoichiometric ceiling.
    pub fn advance(
        &self,
        inventory: &[(SpeciesId, Amount)],
        volume: Volume,
        elapsed: Time,
        temperature: Temperature,
        conditions: &QueryPoint,
    ) -> Result<BatchReactionAdvance> {
        if !volume.value().is_finite()
            || volume.value() <= 0.0
            || !elapsed.value().is_finite()
            || elapsed.value() < 0.0
            || !temperature.value().is_finite()
            || temperature.value() <= 0.0
        {
            return Err(BatchReactionError::Invalid(
                "volume and absolute kelvin must be finite positive; time finite nonnegative",
            ));
        }
        if conditions.axes().keys().any(|axis| {
            OWNED_AXES.contains(&axis.as_str()) || self.source.validity.bound(axis).is_none()
        }) {
            return Err(BatchReactionError::Invalid(
                "conditions must contain only additional source-declared axes",
            ));
        }
        let initial = self
            .reaction
            .advance_inventory(inventory, Amount::new(0.0))?;
        let n = [
            stock(&initial, &self.reactants[0]),
            stock(&initial, &self.reactants[1]),
        ];
        let initial_point = self.point(n, volume, Time::new(0.0), temperature, conditions)?;
        self.admit(&initial_point, "initial")?;
        // Time support is known before evaluating progress, so an unsupported
        // duration cannot masquerade as numerical exhaustion of the reaction.
        if elapsed.value()
            > self
                .source
                .validity
                .bound(ELAPSED_TIME_AXIS)
                .expect("admitted card contains elapsed-time support")
                .1
        {
            return Err(BatchReactionError::OutsideValidity { endpoint: "final" });
        }
        let limiting_extent = self.reaction.limiting_extent(&[
            (self.reactants[0].clone(), Amount::new(n[0])),
            (self.reactants[1].clone(), Amount::new(n[1])),
        ])?;
        let coefficients = [self.coefficients[0] as f64, self.coefficients[1] as f64];
        let progresses =
            elapsed.value() > 0.0 && self.rate.value() > 0.0 && n[0] > 0.0 && n[1] > 0.0;
        let extent = if progresses {
            analytic_extent(
                n,
                coefficients,
                self.rate.value(),
                volume.value(),
                elapsed.value(),
            )?
        } else {
            0.0
        };
        let extent = Amount::new(extent);
        let final_inventory = self.reaction.advance_inventory(&initial, extent)?;
        let remaining = [
            stock(&final_inventory, &self.reactants[0]),
            stock(&final_inventory, &self.reactants[1]),
        ];
        if progresses && remaining.iter().any(|n| *n <= 0.0) {
            return Err(BatchReactionError::Arithmetic(
                "scalar extent erased positive finite-time reactant inventory",
            ));
        }
        let final_point = self.point(remaining, volume, elapsed, temperature, conditions)?;
        self.admit(&final_point, "final")?;
        let amounts = self.reaction.amounts(extent)?;
        let reactants = core::array::from_fn(|i| ReactantAccount {
            species: self.reactants[i].clone(),
            coefficient: self.coefficients[i],
            initial: Amount::new(n[i]),
            consumed: Amount::new(coefficients[i] * extent.value()),
            remaining: Amount::new(remaining[i]),
        });
        Ok(BatchReactionAdvance {
            model: self.clone(),
            initial_inventory: initial,
            final_inventory,
            extent,
            limiting_extent,
            amounts,
            reactants,
            volume,
            elapsed,
            temperature,
            initial_point,
            final_point,
        })
    }

    fn point(
        &self,
        n: [f64; 2],
        volume: Volume,
        elapsed: Time,
        temperature: Temperature,
        conditions: &QueryPoint,
    ) -> Result<QueryPoint> {
        let concentrations = [n[0] / volume.value(), n[1] / volume.value()];
        if concentrations
            .iter()
            .zip(n)
            .any(|(c, n)| !c.is_finite() || (n > 0.0 && *c == 0.0))
        {
            return Err(BatchReactionError::Arithmetic(
                "reactant concentration range loss",
            ));
        }
        let mut point = conditions.clone();
        for (axis, dims, value) in [
            (ELAPSED_TIME_AXIS, Time::DIMS, elapsed.value()),
            (VOLUME_AXIS, Volume::DIMS, volume.value()),
            (TEMPERATURE_AXIS, Temperature::DIMS, temperature.value()),
            (
                REACTANT_A_CONCENTRATION_AXIS,
                AmountConcentration::DIMS,
                concentrations[0],
            ),
            (
                REACTANT_B_CONCENTRATION_AXIS,
                AmountConcentration::DIMS,
                concentrations[1],
            ),
        ] {
            point = point.with_quantity(axis, QuantitySpec::dimensional(dims), value)?;
        }
        Ok(point)
    }

    fn admit(&self, point: &QueryPoint, endpoint: &'static str) -> Result<()> {
        if self
            .source
            .validity
            .contains_typed(point.axes(), point.axis_quantities())
        {
            Ok(())
        } else {
            Err(BatchReactionError::OutsideValidity { endpoint })
        }
    }
}

// Only called with the complete inventory admitted by the original reaction.
fn stock(inventory: &[(SpeciesId, Amount)], species: &SpeciesId) -> f64 {
    inventory
        .iter()
        .find(|(id, _)| id == species)
        .expect("admitted reaction inventory contains each mapped reactant")
        .1
        .value()
}

fn analytic_extent(n: [f64; 2], coefficients: [f64; 2], k: f64, v: f64, t: f64) -> Result<f64> {
    // A=nA/a, B=nB/b; dξ/dt=(kab/V)(A-ξ)(B-ξ).
    let capacities = [n[0] / coefficients[0], n[1] / coefficients[1]];
    let m = capacities[0].min(capacities[1]);
    let big = capacities[0].max(capacities[1]);
    let z = ((((k / v) * coefficients[0]) * coefficients[1]) * big) * t;
    if !m.is_finite() || m <= 0.0 || !big.is_finite() || !z.is_finite() || z <= 0.0 {
        return Err(BatchReactionError::Arithmetic(
            "positive capacity or kinetic scale range loss",
        ));
    }
    let (fraction, remainder) = if m == big {
        (z / (1.0 + z), m / (1.0 + z))
    } else {
        let delta = (big - m) / big;
        let exponent = delta * z;
        let u = -det::expm1(-exponent);
        let denominator = delta + (m / big) * u;
        (
            u / denominator,
            m * (delta / denominator) * det::exp(-exponent),
        )
    };
    let extent = m * fraction;
    if !extent.is_finite()
        || extent <= 0.0
        || extent >= m
        || !remainder.is_finite()
        || remainder <= 0.0
    {
        return Err(BatchReactionError::Arithmetic(
            "positive progress and positive finite-time remainder must both be representable",
        ));
    }
    Ok(extent)
}

#[cfg(test)]
#[path = "reaction_batch/tests.rs"]
mod tests;
