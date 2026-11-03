//! Finite-reagent chemistry in one isothermal, full, perfectly mixed tank.
//!
//! The molar balance is `dn_s/dt = Q C_in,s - (Q/V)n_s + Σ ν_sj r_j`,
//! with explicitly declared local rates `r_j = k_j n_Aj n_Bj / V` in mol/s.
//! Every species, including undissolved minerals and dissolved products, shares
//! the same ideal residence. Mineral inventory/volume is empirical availability,
//! not solid activity, particle surface, or an inferred mineral mechanism.
//!
//! Symmetric exact-transport / simultaneous-RK4-reaction / exact-transport
//! splitting is second order in general. One full step versus two half steps,
//! divided by three, estimates errors in stocks, integrated outlets and extents.
//! Time-weighted local budgets and their accumulated estimates bound accepted
//! work; these estimates are not certified global errors. Source support is
//! checked at every split endpoint and reaction stage. Rectangular support
//! covers the monotone transport and reaction paths of each split subproblem,
//! not a rigorous enclosure of the unsplit continuous trajectory.
//!
//! The original [`ReactionNetwork`] owns exact atom/charge conservation and
//! gross-stock admission at every reaction stage. Products cannot feed any
//! reactant in this scope. Flow and chemistry never clip or repair inventories;
//! independently integrated outlets expose balance residuals. Late failure
//! returns no partial result and cannot mutate the caller's inventory.
//!
//! The general mole balance and ideal mixing assumptions follow Fogler,
//! *Elements of Chemical Reaction Engineering*, chapter 1, equations 1-4/1-8:
//! <https://www.pearson.com/content/dam/one-dot-com/one-dot-com/hei-samples/science/9780135337554.pdf>.
//! Empirical orders follow <https://doi.org/10.1351/goldbook.R05141> and
//! <https://doi.org/10.1351/goldbook.O04322>; neither supplies mineral rates.
//! No heat feedback, speciation, gas transfer, solids slip, variable volume or
//! steady-card/batch-age promotion is included. Call again for a changed feed.

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::{det, eft::two_sum};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Qty, QuantitySpec, Temperature, Time, Volume,
    VolumetricFlowRate,
};
use fs_thermochem::{
    ReactionId, SpeciesId,
    stoichiometry::{
        ConservedReaction,
        network::{NetworkAdvance, ReactionNetwork, ReactionNetworkError},
    },
};
use std::collections::{BTreeMap, BTreeSet};

/// Local rate coefficient in m³/(mol s), with both empirical orders equal to one.
pub type BilinearRateConstant = Qty<3, 0, -1, 0, 0, -1>;
/// Distinct local-law identity; batch elapsed-time laws are never converted.
pub const LOCAL_BILINEAR_LAW: &str = "local-bilinear-reaction";
/// Sole model parameter, with [`BilinearRateConstant`] dimensions.
pub const RATE_PARAMETER: &str = "rate-constant";
/// Fixed active full-tank volume, in m³.
pub const VOLUME_AXIS: &str = "tank-volume";
/// Fixed admitted absolute temperature, in kelvin.
pub const TEMPERATURE_AXIS: &str = "temperature";
/// Actual A stock divided by active tank volume, in mol/m³.
pub const REACTANT_A_CONCENTRATION_AXIS: &str = "reactant-a-concentration";
/// Actual B stock divided by active tank volume, in mol/m³.
pub const REACTANT_B_CONCENTRATION_AXIS: &str = "reactant-b-concentration";
const OWNED_AXES: [&str; 4] = [
    VOLUME_AXIS,
    TEMPERATURE_AXIS,
    REACTANT_A_CONCENTRATION_AXIS,
    REACTANT_B_CONCENTRATION_AXIS,
];
/// Maximum independent mineral pathways sharing one reagent.
pub const MAX_FLOW_REACTIONS: usize = 32;
/// Hard bound on accepted plus rejected full/half-step trials.
pub const MAX_STEP_ATTEMPTS: usize = 65_536;

/// Source, physical input, numerical, or original conserved-inventory refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ReactionFlowError {
    /// Malformed model, incomplete axes, or invalid physical/numerical controls.
    Invalid(&'static str),
    /// Original material-card validation failed.
    Source(MatDbError),
    /// Original local law does not admit the actual staged conditions.
    OutsideValidity { reaction: ReactionId },
    /// Original conserved-network refusal, without changing its meaning.
    Inventory(ReactionNetworkError),
    /// An explicit trial erased a positive finite-time reactant inventory.
    TrialDepletion,
    /// Unrepresentable positive flow, progress, clock, or numerical estimate.
    Arithmetic(&'static str),
    /// No complete result could be produced in the authorized work.
    BudgetExhausted {
        attempts: usize,
        accepted_steps: usize,
        last_trial: &'static str,
    },
}
impl core::fmt::Display for ReactionFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reactive stirred tank: {self:?}")
    }
}
impl std::error::Error for ReactionFlowError {}
impl From<MatDbError> for ReactionFlowError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<ReactionNetworkError> for ReactionFlowError {
    fn from(error: ReactionNetworkError) -> Self {
        Self::Inventory(error)
    }
}
type Result<T> = core::result::Result<T, ReactionFlowError>;

/// Caller-authorized numerical work and accumulated local-estimate budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowIntegrationBudget {
    /// Maximum duration of one proposed complete splitting step.
    pub maximum_step: Time,
    /// Absolute estimated-error allowance for each inventory/outlet/extent [mol].
    pub absolute_amount_error: Amount,
    /// Additional fraction of each component's total available/possible amount.
    pub relative_amount_error: Dimensionless,
    /// Maximum full/half trials, including rejected trials.
    pub maximum_attempts: usize,
}

/// Original local rate declaration and explicit species-to-axis mapping.
#[derive(Clone, Debug)]
pub struct LocalBilinearReaction {
    reaction: ConservedReaction,
    source: ConstitutiveModelCard,
    reactants: [SpeciesId; 2],
    rate: BilinearRateConstant,
}
impl LocalBilinearReaction {
    /// Admit a version-1 local law, zero hidden state, two consumed species and
    /// one finite nonnegative typed k. All source axes need finite typed bounds.
    pub fn from_card(
        reaction: ConservedReaction,
        reactant_a: SpeciesId,
        reactant_b: SpeciesId,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        source.validate()?;
        if source.law.0 != LOCAL_BILINEAR_LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.parameters.len() != 1
            || source.sources.is_empty()
            || source.provenance.source.trim().is_empty()
            || source.provenance.license.trim().is_empty()
        {
            return Err(ReactionFlowError::Invalid(
                "requires original version-1 local law, zero state and one parameter",
            ));
        }
        let rate = source
            .parameters
            .get(RATE_PARAMETER)
            .ok_or(ReactionFlowError::Invalid("missing rate-constant"))?;
        if rate.dims != BilinearRateConstant::DIMS || !rate.value.is_finite() || rate.value < 0.0 {
            return Err(ReactionFlowError::Invalid(
                "k must be finite nonnegative m3/(mol s)",
            ));
        }
        for (axis, (lo, hi)) in source.validity.bounds() {
            if !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
                || !source.validity.axis_quantities().contains_key(axis)
                || [
                    "elapsed-time",
                    "batch-volume",
                    "residence-time",
                    "mean-residence-time",
                    "tank-count",
                ]
                .contains(&axis.as_str())
            {
                return Err(ReactionFlowError::Invalid(
                    "local source axes must be finite typed conditions, not batch age or residence",
                ));
            }
        }
        for (axis, dims) in [
            (VOLUME_AXIS, Volume::DIMS),
            (TEMPERATURE_AXIS, Temperature::DIMS),
            (REACTANT_A_CONCENTRATION_AXIS, AmountConcentration::DIMS),
            (REACTANT_B_CONCENTRATION_AXIS, AmountConcentration::DIMS),
        ] {
            let (lo, _) = source
                .validity
                .bound(axis)
                .ok_or(ReactionFlowError::Invalid(
                    "declare all four local applicability axes",
                ))?;
            if source.validity.axis_quantities().get(axis) != Some(&QuantitySpec::dimensional(dims))
                || lo < 0.0
                || ((axis == VOLUME_AXIS || axis == TEMPERATURE_AXIS) && lo <= 0.0)
            {
                return Err(ReactionFlowError::Invalid(
                    "local axes require positive volume/kelvin and nonnegative concentrations",
                ));
            }
        }
        let reactants = [reactant_a, reactant_b];
        if reactants[0] == reactants[1]
            || reaction
                .terms()
                .iter()
                .filter(|t| t.coefficient < 0)
                .count()
                != 2
            || reactants.iter().any(|id| {
                !reaction
                    .terms()
                    .iter()
                    .any(|t| &t.species == id && t.coefficient < 0)
            })
        {
            return Err(ReactionFlowError::Invalid(
                "map exactly the two distinct consumed species",
            ));
        }
        let rate = BilinearRateConstant::new(rate.value);
        Ok(Self {
            reaction,
            source,
            reactants,
            rate,
        })
    }
    /// Original equation and its exact conservation certificate.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original source hashes, provenance, parameters and applicability domain.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// A then B, as explicitly bound to the source concentration coordinates.
    pub fn reactants(&self) -> &[SpeciesId; 2] {
        &self.reactants
    }
    /// Explicit constant k; no temperature extrapolation is inferred.
    pub const fn rate(&self) -> BilinearRateConstant {
        self.rate
    }
}

/// Fixed-volume CSTR; runtime inventories are explicit call inputs/results.
#[derive(Clone, Debug)]
pub struct ReactiveCstr {
    volume: Volume,
    models: Vec<LocalBilinearReaction>,
    network: ReactionNetwork,
    indices: Vec<[usize; 2]>,
    shared_reagent: SpeciesId,
}

/// Completed interval with independently integrated transport and chemistry.
#[derive(Clone, Debug)]
pub struct ReactionFlowAdvance {
    model: ReactiveCstr,
    /// Complete canonical initial tank inventory [mol].
    pub initial_inventory: Vec<(SpeciesId, Amount)>,
    /// Complete canonical final tank inventory; input for the next interval.
    pub final_inventory: Vec<(SpeciesId, Amount)>,
    /// Actual total incoming amount `Q * elapsed * C_in`, including zeroes.
    pub inlet_amounts: Vec<(SpeciesId, Amount)>,
    /// Independently integrated outlet amounts; never a balance residual plug.
    pub outlet_amounts: Vec<(SpeciesId, Amount)>,
    /// Cumulative reaction accounting funded by initial stocks plus all inlet.
    /// Its `final_inventory` is BEFORE subtracting outlet, not the final tank.
    pub reactions: NetworkAdvance,
    /// `initial + inlet + formed - consumed - outlet - final`, in mol/species.
    pub balance_residuals: Vec<(SpeciesId, Amount)>,
    /// Actual source coordinates at the physical interval boundaries.
    pub source_endpoints: Vec<(ReactionId, QueryPoint, QueryPoint)>,
    /// Piecewise constant actual inlet/outlet volumetric flow.
    pub flow: VolumetricFlowRate,
    /// Completed requested duration.
    pub elapsed: Time,
    /// Fixed absolute temperature for this interval.
    pub temperature: Temperature,
    /// Original numerical authorization.
    pub budget: FlowIntegrationBudget,
    /// Accepted two-half-step candidates.
    pub accepted_steps: usize,
    /// Rejected candidates; none changed caller inventory.
    pub rejected_steps: usize,
    /// Accumulated estimated local inventory errors, not global certificates.
    pub estimated_inventory_errors: Vec<(SpeciesId, Amount)>,
    /// Accumulated estimated local outlet errors.
    pub estimated_outlet_errors: Vec<(SpeciesId, Amount)>,
    /// Accumulated estimated local extent errors on original equation bases.
    pub estimated_extent_errors: Vec<(ReactionId, Amount)>,
}
impl ReactionFlowAdvance {
    /// Original local laws, full source cards and conserved network.
    pub fn model(&self) -> &ReactiveCstr {
        &self.model
    }
}

impl ReactiveCstr {
    /// Join 1..=32 pathways consuming one explicit shared reagent and distinct
    /// other reactants. Products cannot be reactants anywhere in this scope.
    pub fn new(
        volume: Volume,
        mut models: Vec<LocalBilinearReaction>,
        shared_reagent: SpeciesId,
    ) -> Result<Self> {
        if !volume.value().is_finite()
            || volume.value() <= 0.0
            || models.is_empty()
            || models.len() > MAX_FLOW_REACTIONS
        {
            return Err(ReactionFlowError::Invalid(
                "positive finite volume and 1..=32 local models required",
            ));
        }
        models.sort_by(|a, b| a.reaction.id().cmp(b.reaction.id()));
        let network = ReactionNetwork::new(models.iter().map(|m| m.reaction.clone()).collect())?;
        let mut reactants = BTreeSet::from([shared_reagent.clone()]);
        for model in &models {
            let shared = model
                .reactants
                .iter()
                .position(|s| s == &shared_reagent)
                .ok_or(ReactionFlowError::Invalid(
                    "every pathway must consume the shared reagent",
                ))?;
            if !reactants.insert(model.reactants[1 - shared].clone()) {
                return Err(ReactionFlowError::Invalid(
                    "each other reactant must be unique",
                ));
            }
        }
        if models.iter().any(|m| {
            m.reaction
                .terms()
                .iter()
                .any(|t| t.coefficient > 0 && reactants.contains(&t.species))
        }) {
            return Err(ReactionFlowError::Invalid(
                "products cannot feed reactants in this scope",
            ));
        }
        let indices = models
            .iter()
            .map(|m| {
                core::array::from_fn(|i| {
                    network
                        .species()
                        .binary_search(&m.reactants[i])
                        .expect("original reaction participant")
                })
            })
            .collect();
        Ok(Self {
            volume,
            models,
            network,
            indices,
            shared_reagent,
        })
    }
    /// Constant active full volume; chemistry does not infer volume changes.
    pub const fn volume(&self) -> Volume {
        self.volume
    }
    /// Original local models in canonical reaction-ID order.
    pub fn models(&self) -> &[LocalBilinearReaction] {
        &self.models
    }
    /// Original exact conservation and shared gross-demand owner.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.network
    }
    /// Explicit shared finite reagent identity.
    pub fn shared_reagent(&self) -> &SpeciesId {
        &self.shared_reagent
    }

    /// Advance one constant-feed/flow/temperature interval atomically. Supply all
    /// network species exactly once in both inventory and inlet concentration,
    /// including zero products, and one additional-condition query per reaction.
    /// At Q=0 chemistry continues; at k=0 transport uses its exact exponential.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        inventory: &[(SpeciesId, Amount)],
        inlet: &[(SpeciesId, AmountConcentration)],
        flow: VolumetricFlowRate,
        elapsed: Time,
        temperature: Temperature,
        conditions: &[(ReactionId, QueryPoint)],
        budget: FlowIntegrationBudget,
    ) -> Result<ReactionFlowAdvance> {
        validate(flow, elapsed, temperature, budget)?;
        let zeros = vec![0.0; self.models.len()];
        let initial = self
            .network
            .advance_inventory(inventory, &self.requests(&zeros))?
            .final_inventory;
        // Reuse the original complete-species admission before using concentrations.
        let inlet_axis = inlet
            .iter()
            .map(|(id, c)| (id.clone(), Amount::new(c.value())))
            .collect::<Vec<_>>();
        let inlet_axis = self
            .network
            .advance_inventory(&inlet_axis, &self.requests(&zeros))?
            .final_inventory;
        let equilibrium = inlet_axis
            .iter()
            .map(|(_, c)| product(c.value(), self.volume.value()))
            .collect::<Result<Vec<_>>>()?;
        let incoming_volume = product(flow.value(), elapsed.value())?;
        let incoming = inlet_axis
            .iter()
            .map(|(_, c)| product(c.value(), incoming_volume))
            .collect::<Result<Vec<_>>>()?;
        let available = initial
            .iter()
            .zip(&incoming)
            .map(|((id, n), fed)| Ok((id.clone(), Amount::new(add(n.value(), *fed)?))))
            .collect::<Result<Vec<_>>>()?;
        if conditions.len() != self.models.len() {
            return Err(ReactionFlowError::Invalid(
                "one additional-condition query per reaction required",
            ));
        }
        let mut by_id = BTreeMap::new();
        for (id, point) in conditions {
            if by_id.insert(id, point).is_some() {
                return Err(ReactionFlowError::Invalid("duplicate reaction conditions"));
            }
        }
        let ordered = self
            .models
            .iter()
            .map(|m| {
                let point = *by_id
                    .get(m.reaction.id())
                    .ok_or(ReactionFlowError::Invalid("missing reaction conditions"))?;
                if point.axes().keys().any(|a| {
                    OWNED_AXES.contains(&a.as_str()) || m.source.validity.bound(a).is_none()
                }) {
                    return Err(ReactionFlowError::Invalid(
                        "only additional source-declared axes belong in conditions",
                    ));
                }
                Ok(point)
            })
            .collect::<Result<Vec<_>>>()?;
        let problem = Problem {
            owner: self,
            equilibrium: &equilibrium,
            flow,
            temperature,
            conditions: &ordered,
        };
        let initial_points = problem.points(&initial)?;
        let limits = self.error_limits(&available, budget)?;
        let ns = initial.len();
        let mut state = Trial {
            inventory: initial.clone(),
            outlet: vec![0.0; ns],
            extents: zeros,
        };
        let mut errors = vec![0.0; 2 * ns + self.models.len()];
        let mut clock = 0.0;
        let mut step = budget.maximum_step.value().min(elapsed.value());
        let mut accepted = 0;
        let mut rejected = 0;
        let mut last_trial = "another physical step is required";
        if self.models.iter().all(|m| m.rate.value() == 0.0) && elapsed.value() > 0.0 {
            problem.transport(&mut state, elapsed.value())?;
            clock = elapsed.value();
            accepted = 1;
        }
        while clock < elapsed.value() {
            if accepted + rejected == budget.maximum_attempts {
                return Err(ReactionFlowError::BudgetExhausted {
                    attempts: accepted + rejected,
                    accepted_steps: accepted,
                    last_trial,
                });
            }
            // Absorb clock-rounding residue into the final physical step.
            let end = if clock + step >= elapsed.value() * (1.0 - 8.0 * f64::EPSILON) {
                elapsed.value()
            } else {
                clock + step
            };
            let h = end - clock;
            if !end.is_finite() || h <= 0.0 || clock + h * 0.5 == clock || h * 0.5 == 0.0 {
                return Err(ReactionFlowError::Arithmetic(
                    "integration clock lost resolution",
                ));
            }
            let start = Trial {
                inventory: state.inventory.clone(),
                outlet: vec![0.0; ns],
                extents: vec![0.0; self.models.len()],
            };
            let trial = (|| {
                let coarse = problem.split(start.clone(), h)?;
                let first = problem.split(start, h * 0.5)?;
                Ok((coarse, problem.split(first, h * 0.5)?))
            })();
            let (coarse, fine) = match trial {
                Ok(pair) => pair,
                Err(error) if retryable(&error) => {
                    rejected += 1;
                    last_trial = match error {
                        ReactionFlowError::OutsideValidity { .. } => {
                            "trial left original local-source support"
                        }
                        ReactionFlowError::TrialDepletion => {
                            "trial erased positive finite-time reagent stock"
                        }
                        _ => "trial exceeded shared gross inventory",
                    };
                    step = h * 0.5;
                    continue;
                }
                Err(error) => return Err(error),
            };
            let estimates = coarse
                .values()
                .zip(fine.values())
                .map(|(a, b)| {
                    let difference = (a - b).abs();
                    let error = difference / 3.0;
                    if !error.is_finite() || (difference > 0.0 && error == 0.0) {
                        return Err(ReactionFlowError::Arithmetic(
                            "local error estimate lost range",
                        ));
                    }
                    Ok(error)
                })
                .collect::<Result<Vec<_>>>()?;
            let fraction = h / elapsed.value();
            if estimates
                .iter()
                .enumerate()
                .any(|(i, e)| *e > limits[i] * fraction || errors[i] + e > limits[i])
            {
                rejected += 1;
                last_trial = "local or accumulated amount-error estimate exceeded its allowance";
                step = h * 0.5;
                continue;
            }
            let grow = estimates
                .iter()
                .enumerate()
                .all(|(i, e)| *e <= limits[i] * fraction / 8.0);
            for (total, e) in errors.iter_mut().zip(estimates) {
                *total += e;
            }
            for (total, n) in state.outlet.iter_mut().zip(fine.outlet) {
                *total = add(*total, n)?;
            }
            for (total, n) in state.extents.iter_mut().zip(fine.extents) {
                *total = add(*total, n)?;
            }
            state.inventory = fine.inventory;
            clock = end;
            accepted += 1;
            last_trial = "another physical step is required";
            step = if grow {
                (h * 2.0).min(budget.maximum_step.value())
            } else {
                h
            };
        }
        let final_points = problem.points(&state.inventory)?;
        let reactions = self
            .network
            .advance_inventory(&available, &self.requests(&state.extents))?;
        let residuals = (0..ns)
            .map(|i| {
                residual(&[
                    initial[i].1.value(),
                    incoming[i],
                    reactions.formed[i].1.value(),
                    -reactions.consumed[i].1.value(),
                    -state.outlet[i],
                    -state.inventory[i].1.value(),
                ])
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(ReactionFlowAdvance {
            model: self.clone(),
            initial_inventory: initial,
            final_inventory: state.inventory,
            inlet_amounts: self.species_amounts(&incoming),
            outlet_amounts: self.species_amounts(&state.outlet),
            reactions,
            balance_residuals: self.species_amounts(&residuals),
            source_endpoints: self
                .models
                .iter()
                .zip(initial_points)
                .zip(final_points)
                .map(|((m, a), b)| (m.reaction.id().clone(), a, b))
                .collect(),
            flow,
            elapsed,
            temperature,
            budget,
            accepted_steps: accepted,
            rejected_steps: rejected,
            estimated_inventory_errors: self.species_amounts(&errors[..ns]),
            estimated_outlet_errors: self.species_amounts(&errors[ns..2 * ns]),
            estimated_extent_errors: self.requests(&errors[2 * ns..]),
        })
    }
    fn requests(&self, values: &[f64]) -> Vec<(ReactionId, Amount)> {
        self.models
            .iter()
            .zip(values)
            .map(|(m, x)| (m.reaction.id().clone(), Amount::new(*x)))
            .collect()
    }
    fn species_amounts(&self, values: &[f64]) -> Vec<(SpeciesId, Amount)> {
        self.network
            .species()
            .iter()
            .zip(values)
            .map(|(id, x)| (id.clone(), Amount::new(*x)))
            .collect()
    }
    fn error_limits(
        &self,
        available: &[(SpeciesId, Amount)],
        budget: FlowIntegrationBudget,
    ) -> Result<Vec<f64>> {
        let mut stocks = available.iter().map(|(_, n)| n.value()).collect::<Vec<_>>();
        let mut capacities = Vec::new();
        for model in &self.models {
            let capacity = model
                .reaction
                .terms()
                .iter()
                .filter(|t| t.coefficient < 0)
                .map(|t| {
                    let i = self
                        .network
                        .species()
                        .binary_search(&t.species)
                        .expect("original participant");
                    available[i].1.value() / -(t.coefficient as f64)
                })
                .fold(f64::INFINITY, f64::min);
            capacities.push(capacity);
            for term in model.reaction.terms().iter().filter(|t| t.coefficient > 0) {
                let i = self
                    .network
                    .species()
                    .binary_search(&term.species)
                    .expect("original participant");
                stocks[i] = add(stocks[i], product(term.coefficient as f64, capacity)?)?;
            }
        }
        stocks
            .iter()
            .chain(&stocks)
            .chain(&capacities)
            .map(|scale| {
                let limit = budget
                    .relative_amount_error
                    .value()
                    .mul_add(*scale, budget.absolute_amount_error.value());
                if !limit.is_finite() || limit <= 0.0 {
                    return Err(ReactionFlowError::Arithmetic("error allowance lost range"));
                }
                Ok(limit)
            })
            .collect()
    }
}

#[derive(Clone)]
struct Trial {
    inventory: Vec<(SpeciesId, Amount)>,
    outlet: Vec<f64>,
    extents: Vec<f64>,
}
impl Trial {
    fn values(&self) -> impl Iterator<Item = f64> + '_ {
        self.inventory
            .iter()
            .map(|(_, n)| n.value())
            .chain(self.outlet.iter().copied())
            .chain(self.extents.iter().copied())
    }
}
struct Problem<'a> {
    owner: &'a ReactiveCstr,
    equilibrium: &'a [f64],
    flow: VolumetricFlowRate,
    temperature: Temperature,
    conditions: &'a [&'a QueryPoint],
}
impl Problem<'_> {
    fn points(&self, inventory: &[(SpeciesId, Amount)]) -> Result<Vec<QueryPoint>> {
        self.owner
            .models
            .iter()
            .enumerate()
            .map(|(j, m)| {
                let indices = self.owner.indices[j];
                let concentrations =
                    indices.map(|i| inventory[i].1.value() / self.owner.volume.value());
                let mut point = self.conditions[j].clone();
                for (axis, dims, value) in [
                    (VOLUME_AXIS, Volume::DIMS, self.owner.volume.value()),
                    (
                        TEMPERATURE_AXIS,
                        Temperature::DIMS,
                        self.temperature.value(),
                    ),
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
                    if !value.is_finite() {
                        return Err(ReactionFlowError::Arithmetic(
                            "source concentration overflow",
                        ));
                    }
                    point = point.with_quantity(axis, QuantitySpec::dimensional(dims), value)?;
                }
                if concentrations
                    .iter()
                    .zip(indices)
                    .any(|(c, i)| *c == 0.0 && inventory[i].1.value() > 0.0)
                {
                    return Err(ReactionFlowError::Arithmetic(
                        "positive concentration underflow",
                    ));
                }
                if !m
                    .source
                    .validity
                    .contains_typed(point.axes(), point.axis_quantities())
                {
                    return Err(ReactionFlowError::OutsideValidity {
                        reaction: m.reaction.id().clone(),
                    });
                }
                Ok(point)
            })
            .collect()
    }
    fn transport(&self, state: &mut Trial, h: f64) -> Result<()> {
        if self.flow.value() == 0.0 || h == 0.0 {
            return Ok(());
        }
        let x = product(self.flow.value(), h)? / self.owner.volume.value();
        if !x.is_finite() || x <= 0.0 {
            return Err(ReactionFlowError::Arithmetic("flow exposure lost range"));
        }
        let e = det::exp(-x);
        let u = -det::expm1(-x);
        // Integral of inlet-fed stock: x-(1-exp(-x)), evaluated without small-x cancellation.
        let w = if x < 0.125 {
            let mut term = product(x, x)? * 0.5;
            let mut sum = term;
            for n in 3..=16 {
                term *= -x / n as f64;
                sum += term;
            }
            sum
        } else {
            x - u
        };
        for (i, (_, amount)) in state.inventory.iter_mut().enumerate() {
            let before = amount.value();
            let target = self.equilibrium[i];
            let after = if before == target {
                before
            } else {
                add(product(before, e)?, product(target, u)?)?
            };
            if (before > 0.0 && e == 0.0) || (after == before && before != target) {
                return Err(ReactionFlowError::Arithmetic(
                    "finite transport change lost resolution",
                ));
            }
            state.outlet[i] = add(
                state.outlet[i],
                add(product(before, u)?, product(target, w)?)?,
            )?;
            *amount = Amount::new(after);
        }
        self.points(&state.inventory)?;
        Ok(())
    }
    fn reaction_state(
        &self,
        initial: &[(SpeciesId, Amount)],
        extents: &[f64],
    ) -> Result<NetworkAdvance> {
        let stage = self
            .owner
            .network
            .advance_inventory(initial, &self.owner.requests(extents))?;
        for indices in &self.owner.indices {
            if indices
                .iter()
                .any(|i| initial[*i].1.value() > 0.0 && stage.final_inventory[*i].1.value() <= 0.0)
            {
                return Err(ReactionFlowError::TrialDepletion);
            }
        }
        self.points(&stage.final_inventory)?;
        Ok(stage)
    }
    fn rates(&self, initial: &[(SpeciesId, Amount)], extents: &[f64]) -> Result<Vec<f64>> {
        let stage = self.reaction_state(initial, extents)?;
        self.owner
            .models
            .iter()
            .zip(&self.owner.indices)
            .map(|(m, indices)| {
                let a = stage.final_inventory[indices[0]].1.value();
                let b = stage.final_inventory[indices[1]].1.value();
                if m.rate.value() == 0.0 || a == 0.0 || b == 0.0 {
                    return Ok(0.0);
                }
                product(product(m.rate.value(), a / self.owner.volume.value())?, b)
            })
            .collect()
    }
    fn split(&self, mut state: Trial, h: f64) -> Result<Trial> {
        self.transport(&mut state, h * 0.5)?;
        let zero = vec![0.0; self.owner.models.len()];
        let k1 = self.rates(&state.inventory, &zero)?;
        let k2 = self.rates(&state.inventory, &increments(&k1, h * 0.5)?)?;
        let k3 = self.rates(&state.inventory, &increments(&k2, h * 0.5)?)?;
        let k4 = self.rates(&state.inventory, &increments(&k3, h)?)?;
        let extents = (0..zero.len())
            .map(|i| {
                let scale = k1[i].max(k2[i]).max(k3[i]).max(k4[i]);
                let rate = if scale == 0.0 {
                    0.0
                } else {
                    product(
                        scale,
                        (k1[i] / scale
                            + 2.0 * (k2[i] / scale)
                            + 2.0 * (k3[i] / scale)
                            + k4[i] / scale)
                            / 6.0,
                    )?
                };
                product(rate, h)
            })
            .collect::<Result<Vec<_>>>()?;
        state.inventory = self
            .reaction_state(&state.inventory, &extents)?
            .final_inventory;
        for (total, n) in state.extents.iter_mut().zip(extents) {
            *total = add(*total, n)?;
        }
        self.transport(&mut state, h * 0.5)?;
        Ok(state)
    }
}
fn increments(rates: &[f64], h: f64) -> Result<Vec<f64>> {
    rates.iter().map(|r| product(*r, h)).collect()
}
fn product(a: f64, b: f64) -> Result<f64> {
    let value = a * b;
    if !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        return Err(ReactionFlowError::Arithmetic(
            "positive product lost numerical range",
        ));
    }
    Ok(value)
}
fn add(a: f64, b: f64) -> Result<f64> {
    let value = a + b;
    if !value.is_finite() || (b > 0.0 && value == a) {
        return Err(ReactionFlowError::Arithmetic(
            "positive amount increment lost numerical range",
        ));
    }
    Ok(value)
}
fn residual(values: &[f64]) -> Result<f64> {
    let (mut sum, mut correction) = (0.0, 0.0);
    for value in values {
        let (next, error) = two_sum(sum, *value);
        sum = next;
        correction += error;
    }
    let result = sum + correction;
    if !result.is_finite() {
        return Err(ReactionFlowError::Arithmetic("balance residual overflow"));
    }
    Ok(result)
}
fn retryable(error: &ReactionFlowError) -> bool {
    matches!(
        error,
        ReactionFlowError::OutsideValidity { .. }
            | ReactionFlowError::TrialDepletion
            | ReactionFlowError::Inventory(ReactionNetworkError::InsufficientInventory { .. })
    )
}
fn validate(
    flow: VolumetricFlowRate,
    elapsed: Time,
    temperature: Temperature,
    budget: FlowIntegrationBudget,
) -> Result<()> {
    if !flow.value().is_finite()
        || flow.value() < 0.0
        || !elapsed.value().is_finite()
        || elapsed.value() < 0.0
        || !temperature.value().is_finite()
        || temperature.value() <= 0.0
        || !budget.maximum_step.value().is_finite()
        || budget.maximum_step.value() <= 0.0
        || !budget.absolute_amount_error.value().is_finite()
        || budget.absolute_amount_error.value() <= 0.0
        || !budget.relative_amount_error.value().is_finite()
        || !(0.0..1.0).contains(&budget.relative_amount_error.value())
        || !(1..=MAX_STEP_ATTEMPTS).contains(&budget.maximum_attempts)
    {
        return Err(ReactionFlowError::Invalid(
            "finite nonnegative flow/time, positive kelvin/step/absolute error, relative error in [0,1), bounded attempts required",
        ));
    }
    Ok(())
}

#[cfg(test)]
#[path = "reaction_flow/tests.rs"]
mod tests;
