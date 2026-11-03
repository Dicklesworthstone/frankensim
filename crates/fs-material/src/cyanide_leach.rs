//! Finite isothermal gold cyanidation coupled to the original ideal acid/base owner.
//!
//! The explicitly supplied empirical law is `r_Au = k n_Au c_CN_free c_O2`.
//! Its three first-order factors are model choices, not stoichiometric orders or
//! coefficients inferred from a mineral name. Terminology follows IUPAC:
//! <https://doi.org/10.1351/goldbook.R05141> and
//! <https://doi.org/10.1351/goldbook.O04322>. No measured ore rate is provided.
//! Accessible gold means the caller's reactive pool, not total head-grade gold.
//!
//! Cumulative progress x is mol Au dissolved. The original conserved equation
//! owns consumption/formation at extent x/4. Cyanide available outside the gold
//! complex decreases by 2x and dissolved O2 by x/4. Optional competing metal
//! pools redistribute that finite cyanide between HCN, CN and metal complexes.
//! Free CN at a rate stage comes from simultaneous equilibrium, never from
//! initial free CN or a fixed cyanide-consumption factor.
//! The new Au(CN)2- stock enters that charge balance as an explicit spectator.
//! Other analytical families and declared base spectators remain present.
//!
//! The solvent activity basis n_ref is an independent amount, not solution
//! volume. Solvent activity stays one. Stoichiometric OH formation is reported
//! separately from actual equilibrium OH. The gas-transfer adapter supplies
//! finite oxygen exchange and dosing. No closed solvent H/O balance, cyanide
//! destruction, metal redox/dissolution, passivation, heat feedback,
//! concentrated-solution accuracy, or ore calibration is asserted.
//!
//! Bounded RK4 full/two-half trials retain the positive two-half solution.
//! The difference divided by 15 estimates its local Au error; time-weighted
//! allowances control the accumulated estimates. These are not certified error
//! bounds and exclude speciation/source/model uncertainty. Every rate stage and
//! accepted endpoint is source admitted. No inventory is clipped or normalized.

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Qty, QuantitySpec, Temperature, Time, Volume,
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction, ReactionAmount, ReactionBalanceError,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::acid_base::{
            AcidBaseError, AcidBaseSpeciation, IdealAcidBaseSystem, SpeciationOptions,
            SpectatorIon,
            complexation::{LigandSpeciation, LigandSpeciationError},
        },
        thermodynamics::gas_liquid::GasLiquidError,
    },
};
use std::{collections::BTreeSet, sync::Arc};

/// Finite competing metal pools sharing cyanide and the original pH owner.
pub mod competition;
/// Finite gas/liquid oxygen exchange coupled to the original cyanidation owner.
pub mod gas_transfer;

/// Prescribed-pressure oxygen uptake and stripping with a maintained gas reservoir.
pub mod aeration;

/// Metered dissolved NaCN/NaOH additions with retained original charge and clock.
pub mod dosing;

/// Explicit empirical coefficient in m6/(mol2 s).
pub type CyanidationRateConstant = Qty<6, 0, -1, 0, 0, -2>;
/// Accessible gold dissolution rate in mol Au/s.
pub type GoldDissolutionRate = Qty<0, 0, -1, 0, 0, 1>;
/// Version-one source identity for this particular empirical law.
pub const CYANIDATION_LAW: &str = "speciated-gold-cyanidation";
/// Sole, nonnegative coefficient; no built-in ore constant exists.
pub const RATE_PARAMETER: &str = "rate-constant";
/// Actual ideal-model pH, dimensionless.
pub const PH_AXIS: &str = "pH";
/// Actual uncomplexed cyanide ion concentration, mol/m3.
pub const FREE_CYANIDE_AXIS: &str = "free-cyanide-concentration";
/// Actual dissolved molecular oxygen concentration, mol/m3.
pub const OXYGEN_AXIS: &str = "dissolved-oxygen-concentration";
/// Actual remaining accessible gold, mol Au.
pub const GOLD_AXIS: &str = "accessible-gold-inventory";
/// Absolute batch clock, seconds since original initialization.
pub const ELAPSED_TIME_AXIS: &str = "elapsed-time";
/// Actual fixed temperature, using dimensional kelvin coordinates.
pub const TEMPERATURE_AXIS: &str = "temperature";
/// Actual fixed solution volume, m3.
pub const VOLUME_AXIS: &str = "solution-volume";
/// Optional source axis for actual uncomplexed HCN+CN analytical amount, mol.
pub const TOTAL_CYANIDE_AXIS: &str = "uncomplexed-cyanide-inventory";
/// Optional source axis for actual gold-complex amount, mol.
pub const GOLD_COMPLEX_AXIS: &str = "gold-cyanide-complex-inventory";
/// Optional source axis for the explicit solvent/activity-reference amount, mol.
pub const ACTIVITY_REFERENCE_AXIS: &str = "amount-at-unit-activity";

const AU: &str = "Au(s)";
const CN: &str = "CN[-](aq)";
const HCN: &str = "HCN(aq)";
const OXYGEN: &str = "O2(aq)";
const COMPLEX: &str = "Au(CN)2[-](aq)";
const HYDROXIDE: &str = "OH[-](aq)";

/// A refusal leaves the accepted state and clock unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum CyanidationError {
    /// Original kinetic source-card or typed-coordinate admission failure.
    Source(MatDbError),
    /// Original coupled acid/base solver refusal.
    Speciation(AcidBaseError),
    /// Original finite metal/ligand/pH equilibrium refused the staged state.
    Complexation(LigandSpeciationError),
    /// Original gas/aqueous equilibrium or source admission refused.
    GasEquilibrium(GasLiquidError),
    /// Original exact chemical equation/amount refusal.
    Reaction(ReactionBalanceError),
    /// Original shared-species definition admission failed.
    Network(ReactionNetworkError),
    /// Missing, inconsistent or nonphysical declarations.
    Invalid(&'static str),
    /// Actual staged conditions lack original kinetic support.
    OutsideValidity {
        axis: String,
        value: f64,
        minimum: f64,
        maximum: f64,
    },
    /// An explicit trial would exhaust a positive finite-time reactant; reduce its step.
    TrialDepletion,
    /// Required positive stock, progress, clock or estimate cannot be represented.
    Arithmetic(&'static str),
    /// The explicit integration/work limit was exhausted without a complete advance.
    Budget,
    /// Cancellation was observed before publication.
    Cancelled,
}
impl core::fmt::Display for CyanidationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "speciated gold cyanidation: {self:?}")
    }
}
impl std::error::Error for CyanidationError {}
impl From<MatDbError> for CyanidationError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e)
    }
}
impl From<AcidBaseError> for CyanidationError {
    fn from(e: AcidBaseError) -> Self {
        Self::Speciation(e)
    }
}
impl From<LigandSpeciationError> for CyanidationError {
    fn from(e: LigandSpeciationError) -> Self {
        Self::Complexation(e)
    }
}
impl From<GasLiquidError> for CyanidationError {
    fn from(e: GasLiquidError) -> Self {
        Self::GasEquilibrium(e)
    }
}
impl From<ReactionBalanceError> for CyanidationError {
    fn from(e: ReactionBalanceError) -> Self {
        Self::Reaction(e)
    }
}
impl From<ReactionNetworkError> for CyanidationError {
    fn from(e: ReactionNetworkError) -> Self {
        Self::Network(e)
    }
}
type Result<T> = core::result::Result<T, CyanidationError>;

/// Immutable empirical law, original exact reaction, and source-backed acid/base system.
#[derive(Clone, Debug)]
pub struct GoldCyanidation {
    chemistry: IdealAcidBaseSystem,
    cyanide_family: String,
    source: ConstitutiveModelCard,
    reaction: ConservedReaction,
    rate: CyanidationRateConstant,
}
impl GoldCyanidation {
    /// Bind the explicit HCN/CN family and the dissolved-oxygen Elsner equation.
    /// Shared-species admission verifies the original CN/water/OH definitions.
    pub fn from_card(
        chemistry: IdealAcidBaseSystem,
        cyanide_family: impl Into<String>,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        source.validate()?;
        require(
            source.law.0 == CYANIDATION_LAW
                && source.law_version == 1
                && source.state_schema_version == 0
                && source.initial_state == InitialStatePolicy::ZeroInternalState
                && source.parameters.len() == 1
                && !source.sources.is_empty()
                && !source.provenance.source.trim().is_empty()
                && !source.provenance.license.trim().is_empty(),
            "require a sourced version-one empirical law with one rate parameter",
        )?;
        let parameter = source
            .parameters
            .get(RATE_PARAMETER)
            .ok_or(CyanidationError::Invalid("missing rate constant"))?;
        require(
            parameter.dims == CyanidationRateConstant::DIMS,
            "rate constant must be m6/(mol2 s)",
        )?;
        nonnegative(
            parameter.value,
            "rate constant must be finite and nonnegative",
        )?;
        let rate = CyanidationRateConstant::new(parameter.value);
        for (axis, &(low, high)) in source.validity.bounds() {
            require(
                low.is_finite()
                    && high.is_finite()
                    && low <= high
                    && source.validity.axis_quantities().contains_key(axis),
                "all source axes need finite ordered typed bounds",
            )?;
        }
        for (index, &(axis, dims)) in owned_axes().iter().enumerate() {
            let Some((low, high)) = source.validity.bound(axis) else {
                require(index >= 7, "missing required kinetic source axis")?;
                continue;
            };
            require(
                source.validity.axis_quantities().get(axis)
                    == Some(&QuantitySpec::dimensional(dims))
                    && (axis == PH_AXIS || low >= 0.0)
                    && (!matches!(
                        axis,
                        TEMPERATURE_AXIS | VOLUME_AXIS | ACTIVITY_REFERENCE_AXIS
                    ) || low > 0.0)
                    && (axis != ELAPSED_TIME_AXIS || (low == 0.0 && high > 0.0)),
                "invalid owned kinetic source axis",
            )?;
        }
        let cyanide_family = cyanide_family.into();
        let family = chemistry
            .families()
            .iter()
            .find(|f| f.id() == cyanide_family)
            .ok_or(CyanidationError::Invalid(
                "unknown uncomplexed-cyanide family",
            ))?;
        require(
            family.steps().len() == 1
                && family.species().len() == 2
                && family.species()[0].as_str() == HCN
                && family.species()[1].as_str() == CN
                && chemistry.water().proton().as_str() == "H+(aq)"
                && chemistry.water().water().as_str() == "H2O(l)"
                && chemistry.water().hydroxide().as_str() == HYDROXIDE,
            "requires the explicit HCN/CN family and original water/proton/hydroxide roles",
        )?;
        let reaction = MineralReaction::GoldCyanidationDissolvedOxygen.build()?;
        let participants: BTreeSet<_> = reaction.terms().iter().map(|term| &term.species).collect();
        require(
            chemistry
                .families()
                .iter()
                .filter(|family| family.id() != cyanide_family)
                .all(|family| {
                    family
                        .species()
                        .iter()
                        .all(|species| !participants.contains(species))
                }),
            "other acid families cannot own a cyanidation reactant or product",
        )?;
        ReactionNetwork::new(
            chemistry
                .reaction_network()
                .reactions()
                .iter()
                .cloned()
                .chain(std::iter::once(reaction.clone()))
                .collect(),
        )?;
        let temperature = chemistry.water().standard().temperature().value();
        let (low, high) = source
            .validity
            .bound(TEMPERATURE_AXIS)
            .expect("required source axis");
        require(
            (low..=high).contains(&temperature),
            "kinetic source does not cover the acid/base reference temperature",
        )?;
        Ok(Self {
            chemistry,
            cyanide_family,
            source,
            reaction,
            rate,
        })
    }
    /// Original source-backed coupled acid/base owner and all Gibbs receipts.
    pub fn chemistry(&self) -> &IdealAcidBaseSystem {
        &self.chemistry
    }
    /// Explicit uncomplexed HCN+CN family identity.
    pub fn cyanide_family(&self) -> &str {
        &self.cyanide_family
    }
    /// Original independent empirical rate source and applicability.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Original exact dissolved-oxygen gold cyanidation equation.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original nonnegative empirical coefficient.
    pub const fn rate_constant(&self) -> CyanidationRateConstant {
        self.rate
    }
}

/// Original accessible reactants and explicitly declared fixed solution basis.
#[derive(Clone, Debug, PartialEq)]
pub struct CyanidationCharge {
    /// Reactive/accessible gold pool, mol Au; permanently locked gold is excluded.
    pub accessible_gold: Amount,
    /// Supplied HCN+CN, mol, before optional competing-metal equilibration.
    /// Excludes ligands in the separately supplied initial gold complex.
    pub uncomplexed_cyanide: Amount,
    /// Dissolved molecular oxygen, mol O2; no headspace or transfer is inferred.
    pub dissolved_oxygen: Amount,
    /// Initial Au(CN)2- inventory, mol complex.
    pub gold_complex: Amount,
    /// Common ideal activity basis; for unit molality, solvent kg times 1 mol/kg.
    pub amount_at_unit_activity: Amount,
    /// Fixed volume used ONLY to turn dissolved amounts into rate concentrations.
    pub solution_volume: Volume,
    /// Must equal the retained acid/base Gibbs-source temperature.
    pub temperature: Temperature,
    /// Explicit fixed counterions, including those accompanying added metal.
    /// Omit the internally owned gold complex and optional competing species.
    pub spectators: Vec<SpectatorIon>,
    /// Every other original acid-family total, exactly once, including zero totals.
    pub other_families: Vec<(String, Amount)>,
}

/// One actual, source-admitted simultaneous state, including the complete speciation.
#[derive(Clone, Debug)]
pub struct CyanidationState {
    /// Absolute time from original batch initialization.
    pub elapsed: Time,
    /// Cumulative mol Au dissolved; the written-equation extent is this divided by four.
    pub gold_dissolved: Amount,
    /// Remaining original accessible gold, mol.
    pub accessible_gold: Amount,
    /// Remaining uncomplexed HCN+CN family, mol.
    pub uncomplexed_cyanide: Amount,
    /// Remaining dissolved O2 molecules, mol.
    pub dissolved_oxygen: Amount,
    /// Initial plus newly formed Au(CN)2- complex, mol.
    pub gold_complex: Amount,
    /// Actual uncomplexed cyanide-ion amount, mol, from the shared equilibrium.
    pub free_cyanide: Amount,
    /// Original coupled speciation with actual OH, pH, families and spectator stocks.
    pub speciation: AcidBaseSpeciation,
    /// Optional simultaneous competing-metal distributions and cyanide balance.
    /// Its aqueous report is the same complete solution as `speciation`.
    pub competition: Option<LigandSpeciation>,
    /// Complete actual kinetic source coordinates.
    pub conditions: QueryPoint,
    /// Actual empirical mol Au/s at this state.
    pub rate: GoldDissolutionRate,
}
/// Recomputed residuals in explicitly tracked analytical pools, not solvent H/O.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CyanidationConservation {
    /// Remaining accessible Au plus complex Au minus original total Au.
    pub gold: Amount,
    /// Uncomplexed family plus metal-bound ligands and twice gold-complex stock,
    /// minus original cyanide. Includes the reported equilibrium solve residual.
    pub cyanide: Amount,
    /// Fixed canonical Na+ spectator minus original and delivered sodium.
    pub sodium: Amount,
    /// Remaining O2 plus original-equation O2 consumption minus original DO.
    pub dissolved_oxygen: Amount,
}
/// Accuracy estimates and bounded runtime work, independent of source uncertainty.
#[derive(Clone, Copy, Debug)]
pub struct CyanidationNumerics {
    /// Maximum trial duration.
    pub maximum_step: Time,
    /// Absolute cumulative Au local-error-estimate allowance for this call.
    pub absolute_gold_error: Amount,
    /// Additional fraction of the current finite Au progress ceiling.
    pub relative_gold_error: Dimensionless,
    /// Maximum accepted plus rejected full/two-half trials, at most 65,536.
    pub maximum_attempts: usize,
    /// Shared rate/speciation cancellation-work units, at most ten million.
    pub maximum_work: usize,
}
impl Default for CyanidationNumerics {
    fn default() -> Self {
        Self {
            maximum_step: Time::new(5.0),
            absolute_gold_error: Amount::new(1e-12),
            relative_gold_error: Dimensionless::new(1e-7),
            maximum_attempts: 8192,
            maximum_work: 2_000_000,
        }
    }
}
#[derive(Clone, Debug)]
struct Layout {
    model: GoldCyanidation,
    charge: CyanidationCharge,
    conditions: QueryPoint,
    speciation: SpeciationOptions,
    competition: Option<competition::BoundCompetition>,
    limit: f64,
    reagents: dosing::ReagentInputs,
}
/// Persistent closed-bath progress and clock; no partial trial state is public.
#[derive(Clone, Debug)]
pub struct GoldCyanidationBatch {
    layout: Arc<Layout>,
    state: CyanidationState,
    accumulated_error: f64,
}
/// One atomic accepted advance with original exact-equation accounting.
#[derive(Clone, Debug)]
pub struct CyanidationAdvance {
    layout: Arc<Layout>,
    /// Clock at the start of this call.
    pub previous_elapsed: Time,
    /// Accepted absolute clock.
    pub elapsed: Time,
    /// Mol Au dissolved during this call, distinct from cumulative state progress.
    pub interval_gold_dissolved: Amount,
    /// Accepted final stocks, source point and simultaneous acid/base solution.
    pub state: CyanidationState,
    /// Cumulative written-equation extent, mol; Au dissolved / 4.
    pub reaction_extent: Amount,
    /// Original equation's cumulative consumed/formed amounts. Stoichiometric OH
    /// is distinct from the actual equilibrium OH in `state.speciation`.
    pub stoichiometric_amounts: Vec<ReactionAmount>,
    /// Recomputed analytical-pool residuals; no finite solvent balance is implied.
    pub conservation: CyanidationConservation,
    /// Largest accepted local full/two-half difference divided by 15, mol Au.
    pub maximum_local_error: Amount,
    /// Sum of accepted local estimates across the retained history, not a certified bound.
    pub accumulated_error_estimate: Amount,
    /// Accepted full intervals (each retaining its two-half approximation).
    pub accepted_steps: usize,
    /// Rejected over-depleting or inaccurate trials.
    pub rejected_steps: usize,
    /// Work charged, including nested speciation and rejected trials.
    pub work: usize,
}
impl CyanidationAdvance {
    /// Cumulative accepted reagent inputs, separate from the original charge.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.reagents
    }
    /// Original immutable kinetic and equilibrium source owners.
    pub fn model(&self) -> &GoldCyanidation {
        &self.layout.model
    }
}
impl GoldCyanidationBatch {
    /// Validate the complete original charge and solve its initial acid/base state.
    pub fn new(
        model: GoldCyanidation,
        charge: CyanidationCharge,
        conditions: QueryPoint,
        speciation: SpeciationOptions,
    ) -> Result<Self> {
        Self::new_inner(model, charge, conditions, speciation, None)
    }

    fn new_inner(
        model: GoldCyanidation,
        charge: CyanidationCharge,
        conditions: QueryPoint,
        speciation: SpeciationOptions,
        competition: Option<competition::BoundCompetition>,
    ) -> Result<Self> {
        for value in [
            charge.accessible_gold,
            charge.uncomplexed_cyanide,
            charge.dissolved_oxygen,
            charge.gold_complex,
        ] {
            nonnegative(
                value.value(),
                "original reactant/complex stocks must be finite and nonnegative",
            )?;
        }
        positive(
            charge.amount_at_unit_activity.value(),
            "explicit solvent/activity amount must be positive",
        )?;
        positive(
            charge.solution_volume.value(),
            "solution volume must be positive",
        )?;
        positive(
            charge.temperature.value(),
            "absolute temperature must be positive",
        )?;
        require(
            charge.temperature == model.chemistry.water().standard().temperature(),
            "bath temperature must equal the retained equilibrium reference temperature",
        )?;
        require(
            owned_axes()
                .iter()
                .all(|(axis, _)| !conditions.axes().contains_key(*axis)),
            "caller cannot override actual kinetic coordinates",
        )?;
        let forbidden: BTreeSet<_> = model.reaction.terms().iter().map(|t| &t.species).collect();
        require(
            charge.spectators.len() < 256
                && charge
                    .spectators
                    .iter()
                    .all(|s| !forbidden.contains(&s.species)),
            "reserve reactive species and the gold-complex spectator for the original owners",
        )?;
        require(
            charge
                .other_families
                .iter()
                .all(|(id, _)| id != &model.cyanide_family),
            "uncomplexed cyanide total is separately owned",
        )?;
        // These finite totals also admit the analytical residual scale.
        sum([charge.accessible_gold.value(), charge.gold_complex.value()].into_iter())?;
        sum([
            charge.uncomplexed_cyanide.value(),
            product(2.0, charge.gold_complex.value())?,
        ]
        .into_iter())?;
        let mut limit = charge
            .accessible_gold
            .value()
            .min(product(0.5, charge.uncomplexed_cyanide.value())?);
        if charge.dissolved_oxygen.value() < product(0.25, limit)? {
            limit = product(4.0, charge.dissolved_oxygen.value())?;
        }
        let layout = Arc::new(Layout {
            model,
            charge,
            conditions,
            speciation,
            competition,
            limit,
            reagents: dosing::ReagentInputs::default(),
        });
        let mut never = || false;
        let mut work = Work {
            used: 0,
            maximum: usize::MAX,
            cancelled: &mut never,
            stopped: None,
        };
        let state = evaluate(&layout, 0.0, 0.0, &mut work)?;
        Ok(Self {
            layout,
            state,
            accumulated_error: 0.0,
        })
    }
    /// Original kinetic and source-backed equilibrium owners.
    pub fn model(&self) -> &GoldCyanidation {
        &self.layout.model
    }
    /// Original accessible stocks, fixed solvent basis and explicit counterions.
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.charge
    }
    /// Cumulative delivered NaCN/NaOH and sodium, not a replacement initial charge.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.reagents
    }

    /// Add finite already-dissolved equivalents and re-speciate at the retained
    /// clock/progress, including any original competing-metal equilibrium.
    pub fn dose_reagent(
        &mut self,
        reagent: dosing::DissolvedSodiumReagent,
        amount: Amount,
        maximum_work: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<dosing::ReagentDoseReceipt> {
        let mut work = dosing::work(maximum_work, &mut cancelled)?;
        let mut layout = dosing::prepare(&self.layout, reagent, amount, &mut work)?;
        if layout.charge.dissolved_oxygen.value() < product(0.25, layout.limit)? {
            layout.limit = product(4.0, layout.charge.dissolved_oxygen.value())?;
        }
        let layout = Arc::new(layout);
        let state = evaluate(
            &layout,
            self.elapsed().value(),
            self.state.gold_dissolved.value(),
            &mut work,
        )?;
        let receipt = dosing::receipt(layout.clone(), reagent, amount, &state, &mut work)?;
        self.layout = layout;
        self.state = state;
        Ok(receipt)
    }
    /// Current private-history snapshot through an immutable reference.
    pub fn state(&self) -> &CyanidationState {
        &self.state
    }
    /// Accepted absolute clock.
    pub const fn elapsed(&self) -> Time {
        self.state.elapsed
    }
    /// Current Au/cyanide/original-DO ceiling, not a finite-time conversion promise.
    pub fn maximum_gold_dissolution(&self) -> Amount {
        Amount::new(self.layout.limit)
    }

    /// Advance with simultaneous stage speciation; any failure leaves the bath intact.
    pub fn advance(
        &mut self,
        duration: Time,
        numerics: CyanidationNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CyanidationAdvance> {
        nonnegative(
            duration.value(),
            "advance duration must be finite and nonnegative",
        )?;
        let tolerance = numerics.tolerance(self.layout.limit)?;
        let end = self.state.elapsed.value() + duration.value();
        positive_or_zero(end, "absolute clock overflows")?;
        if duration.value() > 0.0 && end == self.state.elapsed.value() {
            return Err(CyanidationError::Arithmetic(
                "positive duration is below clock resolution",
            ));
        }
        let mut work = Work {
            used: 0,
            maximum: numerics.maximum_work,
            cancelled: &mut cancelled,
            stopped: None,
        };
        work.tick()?;
        // Refuse an unsupported requested absolute clock before doing rate work.
        let (low, high) = self
            .layout
            .model
            .source
            .validity
            .bound(ELAPSED_TIME_AXIS)
            .expect("required time axis");
        if end < low || end > high {
            return Err(CyanidationError::OutsideValidity {
                axis: ELAPSED_TIME_AXIS.into(),
                value: end,
                minimum: low,
                maximum: high,
            });
        }
        let mut state = self.state.clone();
        let mut step = numerics.maximum_step.value().min(duration.value());
        let (mut accepted, mut rejected, mut largest, mut accumulated) =
            (0, 0, 0.0_f64, self.accumulated_error);
        while state.elapsed.value() < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CyanidationError::Budget);
            }
            let remaining = end - state.elapsed.value();
            step = step.min(remaining);
            let finish = if step == remaining {
                end
            } else {
                state.elapsed.value() + step
            };
            let middle = state.elapsed.value() + 0.5 * (finish - state.elapsed.value());
            if middle <= state.elapsed.value() || finish <= middle {
                return Err(CyanidationError::Arithmetic(
                    "positive substep is below clock resolution",
                ));
            }
            let trial = (|| -> Result<_> {
                let full = rk4(&self.layout, &state, finish, &mut work)?;
                let first = rk4(&self.layout, &state, middle, &mut work)?;
                let refined = rk4(&self.layout, &first, finish, &mut work)?;
                Ok((full, refined))
            })();
            let (full, refined) = match trial {
                Err(CyanidationError::TrialDepletion) => {
                    rejected += 1;
                    step = product(step, 0.5)?;
                    continue;
                }
                other => other?,
            };
            let error = quotient(
                (refined.gold_dissolved.value() - full.gold_dissolved.value()).abs(),
                15.0,
            )?;
            let local = product(
                tolerance,
                quotient(finish - state.elapsed.value(), duration.value())?,
            )?;
            if error > local {
                rejected += 1;
                step = product(step, 0.5)?;
                continue;
            }
            state = refined;
            accepted += 1;
            largest = largest.max(error);
            accumulated = sum([accumulated, error].into_iter())?;
            if error < local / 32.0 {
                step = (step * 2.0).min(numerics.maximum_step.value());
            }
        }
        let extent = Amount::new(quotient(state.gold_dissolved.value(), 4.0)?);
        let amounts = self.layout.model.reaction.amounts(extent)?;
        let conservation = conservation(&self.layout, &state, &amounts)?;
        let mut answer = CyanidationAdvance {
            layout: self.layout.clone(),
            previous_elapsed: self.state.elapsed,
            elapsed: state.elapsed,
            interval_gold_dissolved: Amount::new(
                state.gold_dissolved.value() - self.state.gold_dissolved.value(),
            ),
            state: state.clone(),
            reaction_extent: extent,
            stoichiometric_amounts: amounts,
            conservation,
            maximum_local_error: Amount::new(largest),
            accumulated_error_estimate: Amount::new(accumulated),
            accepted_steps: accepted,
            rejected_steps: rejected,
            work: work.used,
        };
        work.tick()?;
        answer.work = work.used;
        self.state = state;
        self.accumulated_error = accumulated;
        Ok(answer)
    }
}

fn rk4(
    layout: &Layout,
    initial: &CyanidationState,
    end: f64,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CyanidationState> {
    let duration = end - initial.elapsed.value();
    let half = product(duration, 0.5)?;
    let middle = initial.elapsed.value() + half;
    let x = initial.gold_dissolved.value();
    let r1 = initial.rate.value();
    let second = evaluate(layout, middle, progress(x, half, r1, layout.limit)?, work)?;
    let third = evaluate(
        layout,
        middle,
        progress(x, half, second.rate.value(), layout.limit)?,
        work,
    )?;
    let fourth = evaluate(
        layout,
        end,
        progress(x, duration, third.rate.value(), layout.limit)?,
        work,
    )?;
    let average = sum([
        quotient(r1, 6.0)?,
        quotient(second.rate.value(), 3.0)?,
        quotient(third.rate.value(), 3.0)?,
        quotient(fourth.rate.value(), 6.0)?,
    ]
    .into_iter())?;
    evaluate(
        layout,
        end,
        progress(x, duration, average, layout.limit)?,
        work,
    )
}
fn progress(x: f64, time: f64, rate: f64, limit: f64) -> Result<f64> {
    if rate == 0.0 {
        return Ok(x);
    }
    let increment = time * rate;
    let next = x + increment;
    if !next.is_finite() || next >= limit {
        return Err(CyanidationError::TrialDepletion);
    }
    if increment == 0.0 || next == x {
        return Err(CyanidationError::Arithmetic(
            "positive Au progress is below numeric resolution",
        ));
    }
    Ok(next)
}
fn evaluate(
    layout: &Layout,
    elapsed: f64,
    x: f64,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CyanidationState> {
    evaluate_with_oxygen(layout, elapsed, x, None, work)
}

// The gas/liquid owner supplies an independently integrated dissolved stock;
// the closed-batch path retains its original stoichiometric oxygen remainder.
fn evaluate_with_oxygen(
    layout: &Layout,
    elapsed: f64,
    x: f64,
    oxygen_stock: Option<f64>,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CyanidationState> {
    work.tick()?;
    let amounts = layout
        .model
        .reaction
        .amounts(Amount::new(quotient(x, 4.0)?))?;
    let c = &layout.charge;
    let gold = remainder(c.accessible_gold.value(), amount(&amounts, AU), x)?;
    let cyanide = remainder(dosing::cyanide_total(layout)?, amount(&amounts, CN), x)?;
    let oxygen = match oxygen_stock {
        Some(value) => positive_or_zero(value, "integrated dissolved oxygen loses range")?,
        None => remainder(c.dissolved_oxygen.value(), amount(&amounts, OXYGEN), x)?,
    };
    let complex = sum([c.gold_complex.value(), amount(&amounts, COMPLEX)].into_iter())?;
    let mut spectators = dosing::spectators(layout)?;
    spectators.push(SpectatorIon {
        species: SpeciesId::new(COMPLEX).expect("fixed species identity"),
        charge: -1,
        amount: Amount::new(complex),
    });
    let equilibrium = match &layout.competition {
        Some(metals) => metals
            .model
            .speciate(
                Amount::new(cyanide),
                &metals.totals,
                &c.other_families,
                &spectators,
                c.amount_at_unit_activity,
                metals.numerics,
                || work.interrupted(),
            )
            .map(|result| {
                (
                    result.aqueous.clone(),
                    result.uncomplexed_ligand.value(),
                    Some(result),
                )
            })
            .map_err(CyanidationError::from),
        None => {
            let mut totals = c.other_families.clone();
            totals.push((layout.model.cyanide_family.clone(), Amount::new(cyanide)));
            layout
                .model
                .chemistry
                .speciate(
                    &totals,
                    &spectators,
                    c.amount_at_unit_activity,
                    layout.speciation,
                    || work.interrupted(),
                )
                .map(|result| (result, cyanide, None))
                .map_err(CyanidationError::from)
        }
    };
    if let Some(error) = work.stopped.take() {
        return Err(error);
    }
    let (speciation, cyanide, competition) = equilibrium?;
    let family = speciation
        .families
        .iter()
        .find(|f| f.id == layout.model.cyanide_family)
        .expect("original family retained");
    let free_cyanide = family
        .forms
        .iter()
        .find(|f| f.species.as_str() == CN)
        .expect("original CN member retained")
        .amount;
    let free_concentration = quotient(free_cyanide.value(), c.solution_volume.value())?;
    let oxygen_concentration = quotient(oxygen, c.solution_volume.value())?;
    let values = [
        speciation.ph.value(),
        free_concentration,
        oxygen_concentration,
        gold,
        elapsed,
        c.temperature.value(),
        c.solution_volume.value(),
        cyanide,
        complex,
        c.amount_at_unit_activity.value(),
    ];
    let mut point = layout.conditions.clone();
    for ((axis, dims), value) in owned_axes().into_iter().zip(values) {
        point = point.with_quantity(axis, QuantitySpec::dimensional(dims), value)?;
    }
    for (axis, &(minimum, maximum)) in layout.model.source.validity.bounds() {
        let value = *point.axes().get(axis).ok_or(CyanidationError::Invalid(
            "missing original kinetic source condition",
        ))?;
        require(
            point.axis_quantities().get(axis)
                == layout.model.source.validity.axis_quantities().get(axis),
            "kinetic source condition has wrong quantity",
        )?;
        if value < minimum || value > maximum {
            return Err(CyanidationError::OutsideValidity {
                axis: axis.clone(),
                value,
                minimum,
                maximum,
            });
        }
    }
    let rate = if layout.model.rate.value() == 0.0
        || gold == 0.0
        || free_concentration == 0.0
        || oxygen_concentration == 0.0
    {
        0.0
    } else {
        product(
            product(
                product(layout.model.rate.value(), gold)?,
                free_concentration,
            )?,
            oxygen_concentration,
        )?
    };
    Ok(CyanidationState {
        elapsed: Time::new(elapsed),
        gold_dissolved: Amount::new(x),
        accessible_gold: Amount::new(gold),
        uncomplexed_cyanide: Amount::new(cyanide),
        dissolved_oxygen: Amount::new(oxygen),
        gold_complex: Amount::new(complex),
        free_cyanide,
        speciation,
        competition,
        conditions: point,
        rate: GoldDissolutionRate::new(rate),
    })
}
fn conservation(
    layout: &Layout,
    state: &CyanidationState,
    amounts: &[ReactionAmount],
) -> Result<CyanidationConservation> {
    let c = &layout.charge;
    let gold = sum([
        state.accessible_gold.value(),
        state.gold_complex.value(),
        -c.accessible_gold.value(),
        -c.gold_complex.value(),
    ]
    .into_iter())?;
    let (cyanide, sodium) = dosing::balances(layout, state)?;
    let oxygen = sum([
        state.dissolved_oxygen.value(),
        amount(amounts, OXYGEN),
        -c.dissolved_oxygen.value(),
    ]
    .into_iter())?;
    for (residual, scale) in [
        (gold, c.accessible_gold.value().max(c.gold_complex.value())),
        (
            cyanide - ligand_balance_residual(state),
            dosing::cyanide_total(layout)?.max(product(2.0, c.gold_complex.value())?),
        ),
        (oxygen, c.dissolved_oxygen.value()),
    ] {
        if residual.abs() > 128.0 * f64::EPSILON * scale {
            return Err(CyanidationError::Arithmetic(
                "analytical inventory balance does not close",
            ));
        }
    }
    Ok(CyanidationConservation {
        gold: Amount::new(gold),
        cyanide: Amount::new(cyanide),
        sodium: Amount::new(sodium),
        dissolved_oxygen: Amount::new(oxygen),
    })
}
fn bound_cyanide(state: &CyanidationState) -> f64 {
    state
        .competition
        .as_ref()
        .map_or(0.0, |result| result.bound_ligand.value())
}
// The independently recomputed analytical balance includes the equilibrium
// owner's explicitly bounded ligand residual. Verify the remaining accounting
// at roundoff scale without erasing or normalizing the published residual.
fn ligand_balance_residual(state: &CyanidationState) -> f64 {
    state
        .competition
        .as_ref()
        .map_or(0.0, |result| result.ligand_balance_residual.value())
}
fn amount(amounts: &[ReactionAmount], species: &str) -> f64 {
    amounts
        .iter()
        .find(|a| a.term.species.as_str() == species)
        .expect("original gold equation participant")
        .amount
        .value()
}
fn remainder(initial: f64, consumed: f64, x: f64) -> Result<f64> {
    let value = initial - consumed;
    if value < 0.0 || (x > 0.0 && initial > 0.0 && value == 0.0) {
        Err(CyanidationError::TrialDepletion)
    } else {
        Ok(value)
    }
}
fn owned_axes() -> [(&'static str, Dims); 10] {
    [
        (PH_AXIS, Dims::NONE),
        (FREE_CYANIDE_AXIS, AmountConcentration::DIMS),
        (OXYGEN_AXIS, AmountConcentration::DIMS),
        (GOLD_AXIS, Amount::DIMS),
        (ELAPSED_TIME_AXIS, Time::DIMS),
        (TEMPERATURE_AXIS, Temperature::DIMS),
        (VOLUME_AXIS, Volume::DIMS),
        (TOTAL_CYANIDE_AXIS, Amount::DIMS),
        (GOLD_COMPLEX_AXIS, Amount::DIMS),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS),
    ]
}
impl CyanidationNumerics {
    fn tolerance(self, limit: f64) -> Result<f64> {
        positive(self.maximum_step.value(), "maximum step must be positive")?;
        positive(
            self.absolute_gold_error.value(),
            "absolute Au error allowance must be positive",
        )?;
        require(
            self.relative_gold_error.value().is_finite()
                && (0.0..=0.01).contains(&self.relative_gold_error.value())
                && (1..=65_536).contains(&self.maximum_attempts)
                && (1..=10_000_000).contains(&self.maximum_work),
            "invalid relative accuracy or work limit",
        )?;
        sum([
            self.absolute_gold_error.value(),
            product(self.relative_gold_error.value(), limit)?,
        ]
        .into_iter())
    }
}
struct Work<'a, F: FnMut() -> bool> {
    used: usize,
    maximum: usize,
    cancelled: &'a mut F,
    stopped: Option<CyanidationError>,
}
impl<F: FnMut() -> bool> Work<'_, F> {
    fn interrupted(&mut self) -> bool {
        if self.stopped.is_some() {
            return true;
        }
        if (self.cancelled)() {
            self.stopped = Some(CyanidationError::Cancelled);
        } else if self.used >= self.maximum {
            self.stopped = Some(CyanidationError::Budget);
        } else {
            self.used += 1;
        }
        self.stopped.is_some()
    }
    fn tick(&mut self) -> Result<()> {
        if self.interrupted() {
            Err(self.stopped.take().expect("interrupted work cause"))
        } else {
            Ok(())
        }
    }
}
fn require(condition: bool, what: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(CyanidationError::Invalid(what))
    }
}
fn nonnegative(value: f64, what: &'static str) -> Result<()> {
    require(value.is_finite() && value >= 0.0, what)
}
fn positive(value: f64, what: &'static str) -> Result<()> {
    require(value.is_finite() && value > 0.0, what)
}
fn positive_or_zero(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(CyanidationError::Arithmetic(what))
    }
}
fn product(a: f64, b: f64) -> Result<f64> {
    let result = positive_or_zero(a * b, "positive product overflows")?;
    if a > 0.0 && b > 0.0 && result == 0.0 {
        Err(CyanidationError::Arithmetic("positive product underflows"))
    } else {
        Ok(result)
    }
}
fn quotient(a: f64, b: f64) -> Result<f64> {
    let result = positive_or_zero(a / b, "positive quotient loses range")?;
    if a > 0.0 && result == 0.0 {
        Err(CyanidationError::Arithmetic("positive quotient underflows"))
    } else {
        Ok(result)
    }
}
fn sum(values: impl Iterator<Item = f64>) -> Result<f64> {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = total + value;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    let result = total + correction;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(CyanidationError::Arithmetic("inventory sum overflows"))
    }
}

#[cfg(test)]
#[path = "cyanide_leach/tests.rs"]
mod tests;
