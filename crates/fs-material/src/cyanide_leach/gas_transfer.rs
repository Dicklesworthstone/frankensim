//! Finite oxygen headspace exchange during source-bound gold cyanidation.
//!
//! The original chemical owner supplies r [mol Au/s] at every stage. With
//! D [mol O2] in the liquid and G [mol O2] in a fixed ideal headspace,
//! `J=kLa*(D_equilibrium(G)-D)`, `dD/dt=J-r/4`, `dG/dt=-J`, `dx/dt=r`.
//! A fourth coordinate independently integrates signed gas-to-liquid transfer.
//! Original Gibbs evidence or a directly measured Henry solubility supplies the
//! partition; a distinct local source supplies kLa [1/s]. The liquid-film form follows the EPA oxygen-transfer
//! proceedings, Section II, "Fundamental Oxygen Transfer Model":
//! <https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=9101KJ0Z.TXT>.
//! No clean-water coefficient is silently assigned to an ore slurry.
//!
//! Source-admitted RK4 stages jointly update gas, liquid, reaction and transfer
//! accounts. Full/two-half difference divided by 15 estimates local error in
//! `|dx|+4*(|dD|+|dG|+|dTransfer|)` Au-equivalent moles. This is not a rigorous
//! global error, source uncertainty or continuous source-domain enclosure.
//! Fixed temperature, solvent activity basis and volumes are explicit. Oxygen
//! can be dosed into the finite headspace between advances; no unmetered
//! reservoir, air composition, bubble hydrodynamics or thermal balance exists.

use super::*;
use fs_blake3::ContentHash;
use fs_matdb::{ClaimSelection, MaterialAnswer, NormalizedPack, PropertyKey};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Frequency, Pressure};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use fs_thermochem::stoichiometry::thermodynamics::gas_liquid::{
    GasLiquidError, IdealGasLiquidPartition,
};

/// Distinct, local liquid-film oxygen transfer source identity.
pub const OXYGEN_TRANSFER_LAW: &str = "local-oxygen-gas-liquid-transfer";
/// Source parameter in inverse seconds; not a gas-flow rate.
pub const KLA_PARAMETER: &str = "volumetric-transfer-rate";
/// Actual oxygen partial pressure in the finite headspace, Pa.
pub const PARTIAL_PRESSURE_AXIS: &str = "oxygen-partial-pressure";
/// Actual fixed headspace volume, m3.
pub const GAS_VOLUME_AXIS: &str = "gas-volume";
/// Signed oxygen transfer, positive from gas to liquid, mol O2/s.
pub type OxygenTransferRate = Qty<0, 0, -1, 0, 0, 1>;
/// IUPAC H_s^cp = dissolved amount concentration / gas partial pressure,
/// mol/(m3 Pa), not its reciprocal or an activity-based equilibrium constant.
/// <https://doi.org/10.1351/goldbook.09024>.
pub type HenrySolubilityCp = Qty<-2, -1, 2, 0, 0, 1>;
/// Original material property; shared convention with the dissolved-gas owner.
pub const HENRY_SOLUBILITY_PROPERTY: &str = "henry_molar_solubility";
const HENRY_FLAGS: [(&str, f64); 5] = [
    ("source-solvent-pure-water", 1.0),
    ("source-henry-infinite-dilution", 1.0),
    ("source-solute-molecular-oxygen", 1.0),
    ("engineering-reference-model", 1.0),
    ("source-specimen-qualified", 0.0),
];

/// Any refused operation leaves the entire accepted state and accounts intact.
#[derive(Clone, Debug, PartialEq)]
pub enum GasTransferError {
    /// Original reaction, speciation, kinetic source, numerical or cancellation refusal.
    Cyanidation(CyanidationError),
    /// Original gas/aqueous equilibrium-source or ideal-partition refusal.
    Partition(GasLiquidError),
}
impl core::fmt::Display for GasTransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "oxygen-fed cyanidation: {self:?}")
    }
}
impl std::error::Error for GasTransferError {}
impl From<CyanidationError> for GasTransferError {
    fn from(error: CyanidationError) -> Self {
        Self::Cyanidation(error)
    }
}
impl From<GasLiquidError> for GasTransferError {
    fn from(error: GasLiquidError) -> Self {
        Self::Partition(error)
    }
}
impl From<MatDbError> for GasTransferError {
    fn from(error: MatDbError) -> Self {
        Self::Cyanidation(error.into())
    }
}
impl From<ReactionBalanceError> for GasTransferError {
    fn from(error: ReactionBalanceError) -> Self {
        Self::Cyanidation(error.into())
    }
}
type Result<T> = core::result::Result<T, GasTransferError>;

/// Direct source evidence for molecular O2 solubility in pure water at infinite
/// dilution. The source does not invent a finite concentration/pressure range,
/// salt correction, heat balance or formation Gibbs energy. Applying it to an
/// ore liquor remains an explicitly declared engineering reference assumption.
#[derive(Clone, Debug)]
pub struct OxygenHenrySource {
    pack: Arc<NormalizedPack>,
    pack_hash: ContentHash,
    point: QueryPoint,
    answer: MaterialAnswer,
}
impl OxygenHenrySource {
    /// Resolve the exact Hcp quantity from original normalized source evidence.
    /// The selected claim itself must declare the solute, solvent, infinite-
    /// dilution convention and typed absolute-temperature support. An initial
    /// selection is explicit; later evaluations pin that exact original claim.
    pub fn from_pack(
        pack: NormalizedPack,
        point: QueryPoint,
        selection: ClaimSelection,
    ) -> Result<Self> {
        let key = PropertyKey::new(HENRY_SOLUBILITY_PROPERTY, HenrySolubilityCp::DIMS);
        let answer = match selection {
            ClaimSelection::Policy(policy) => pack.claims().query_typed(&key, &point, policy)?,
            ClaimSelection::Pinned(id) => pack.claims().query_pinned_typed(&key, &point, id)?,
        };
        pack.claims().verify_receipt(&answer.receipt)?;
        require(
            answer.receipt.observation_backed && !answer.receipt.source_hashes.is_empty(),
            "direct Henry solubility needs original observation-backed evidence",
        )?;
        let claim = pack
            .claims()
            .claim(answer.receipt.selected)
            .expect("resolved original claim");
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        require(
            point.axis_quantities().get(TEMPERATURE_AXIS) == Some(&temperature)
                && claim.validity.axis_quantities().get(TEMPERATURE_AXIS) == Some(&temperature)
                && claim.validity.bound(TEMPERATURE_AXIS).is_some(),
            "Henry source needs its own semantic absolute-temperature support",
        )?;
        for (flag, expected) in HENRY_FLAGS {
            require(
                point.axes().get(flag) == Some(&expected)
                    && point.axis_quantities().get(flag)
                        == Some(&QuantitySpec::dimensional(Dims::NONE))
                    && claim.validity.bound(flag) == Some((expected, expected))
                    && claim.validity.axis_quantities().get(flag)
                        == Some(&QuantitySpec::dimensional(Dims::NONE)),
                "Henry claim must identify molecular O2, pure water, infinite dilution and unqualified engineering reference use",
            )?;
        }
        require(
            point.axes().len() == HENRY_FLAGS.len() + 1
                && claim.validity.bounds().len() == HENRY_FLAGS.len() + 1
                && point
                    .axes()
                    .keys()
                    .all(|axis| point.axis_quantities().contains_key(axis))
                && claim.validity.bounds().iter().all(|(axis, &(lo, hi))| {
                    lo.is_finite()
                        && hi.is_finite()
                        && lo <= hi
                        && claim.validity.axis_quantities().contains_key(axis)
                }),
            "Henry reference supports only typed temperature and its five declared source flags",
        )?;
        positive(
            answer.evidence.value.value,
            "Hcp must be finite and positive",
        )?;
        Ok(Self {
            pack_hash: pack.content_hash(),
            pack: Arc::new(pack),
            point,
            answer,
        })
    }
    /// Immutable original pack, including normalization and provenance records.
    pub fn pack(&self) -> &NormalizedPack {
        &self.pack
    }
    /// Whole-pack identity at source admission.
    pub const fn pack_hash(&self) -> ContentHash {
        self.pack_hash
    }
    /// Original typed source query, including all source qualification flags.
    pub fn point(&self) -> &QueryPoint {
        &self.point
    }
    /// Original source answer and replayable receipt; uncertainty stays unchanged.
    pub fn answer(&self) -> &MaterialAnswer {
        &self.answer
    }

    pub(super) fn at_runtime(&self, actual: &QueryPoint) -> super::Result<MaterialAnswer> {
        let mut point = QueryPoint::new();
        for (axis, quantity) in self.point.axis_quantities() {
            let value = *actual.axes().get(axis).ok_or(CyanidationError::Invalid(
                "missing actual Henry source condition",
            ))?;
            // The physical owner uses dimensional kelvin. This explicit bridge
            // preserves the source's semantic absolute-temperature descriptor.
            let runtime_quantity = if axis == TEMPERATURE_AXIS {
                QuantitySpec::dimensional(Temperature::DIMS)
            } else {
                *quantity
            };
            require(
                actual.axis_quantities().get(axis) == Some(&runtime_quantity),
                "actual Henry source condition has wrong quantity",
            )?;
            point = point.with_quantity(axis, *quantity, value)?;
        }
        let answer = self.pack.claims().query_pinned_typed(
            &PropertyKey::new(HENRY_SOLUBILITY_PROPERTY, HenrySolubilityCp::DIMS),
            &point,
            self.answer.receipt.selected,
        )?;
        positive(
            answer.evidence.value.value,
            "actual Hcp must be finite and positive",
        )?;
        Ok(answer)
    }
}

#[derive(Clone, Debug)]
pub(super) enum OxygenPartition {
    Thermodynamic(IdealGasLiquidPartition),
    Henry(OxygenHenrySource),
}

/// Original oxygen partition evidence and an independent liquid-film transfer source.
#[derive(Clone, Debug)]
pub struct OxygenTransfer {
    partition: OxygenPartition,
    source: ConstitutiveModelCard,
    kla: Frequency,
}
impl OxygenTransfer {
    /// Bind a constant nonnegative kLa to O2(g) -> O2(aq), without inferring a
    /// kinetic coefficient from thermodynamic data or an oxygen/carbon name.
    pub fn from_card(
        equilibrium: IdealGasLiquidPartition,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        require(
            equilibrium.gas().as_str() == "O2(g)" && equilibrium.aqueous().as_str() == OXYGEN,
            "require original molecular oxygen gas/aqueous partition",
        )?;
        Self::with_partition(OxygenPartition::Thermodynamic(equilibrium), source)
    }
    /// Use directly measured Hcp with an independently declared kLa. No gas or
    /// aqueous formation Gibbs values are synthesized from this measurement.
    pub fn from_henry_source(
        henry: OxygenHenrySource,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        Self::with_partition(OxygenPartition::Henry(henry), source)
    }
    fn with_partition(partition: OxygenPartition, source: ConstitutiveModelCard) -> Result<Self> {
        source.validate()?;
        require(
            source.law.0 == OXYGEN_TRANSFER_LAW
                && source.law_version == 1
                && source.state_schema_version == 0
                && source.initial_state == InitialStatePolicy::ZeroInternalState
                && source.parameters.len() == 1
                && !source.sources.is_empty()
                && !source.provenance.source.trim().is_empty()
                && !source.provenance.license.trim().is_empty(),
            "require original oxygen partition and distinct version-one local transfer source",
        )?;
        let parameter = source
            .parameters
            .get(KLA_PARAMETER)
            .ok_or(CyanidationError::Invalid(
                "missing volumetric transfer rate",
            ))?;
        require(
            parameter.dims == Frequency::DIMS,
            "kLa requires inverse-second dimensions",
        )?;
        nonnegative(parameter.value, "kLa must be finite and nonnegative")?;
        let kla = Frequency::new(parameter.value);
        for (axis, &(low, high)) in source.validity.bounds() {
            require(
                low.is_finite()
                    && high.is_finite()
                    && low <= high
                    && source.validity.axis_quantities().contains_key(axis)
                    && ![ELAPSED_TIME_AXIS, "residence-time", "batch-volume"]
                        .contains(&axis.as_str()),
                "local transfer source needs finite typed conditions, not batch age or residence",
            )?;
        }
        for (axis, dims, strictly_positive) in transfer_axes() {
            let (low, high) = source
                .validity
                .bound(axis)
                .ok_or(CyanidationError::Invalid(
                    "missing local oxygen transfer source axis",
                ))?;
            require(
                source.validity.axis_quantities().get(axis)
                    == Some(&QuantitySpec::dimensional(dims))
                    && low >= 0.0
                    && high > 0.0
                    && (!strictly_positive || low > 0.0),
                "invalid local transfer source support",
            )?;
        }
        Ok(Self {
            partition,
            source,
            kla,
        })
    }
    /// Original gas/aqueous Gibbs owner when that route was explicitly selected.
    /// Direct Henry measurements do not carry a fabricated Gibbs model.
    pub fn equilibrium(&self) -> Option<&IdealGasLiquidPartition> {
        match &self.partition {
            OxygenPartition::Thermodynamic(value) => Some(value),
            OxygenPartition::Henry(_) => None,
        }
    }
    /// Original measured Henry source, if that independent route was selected.
    pub fn henry_source(&self) -> Option<&OxygenHenrySource> {
        match &self.partition {
            OxygenPartition::Henry(value) => Some(value),
            OxygenPartition::Thermodynamic(_) => None,
        }
    }
    /// Original transfer source and its complete applicability box.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Declared kLa, inverse seconds; zero disables exchange without changing chemistry.
    pub const fn kla(&self) -> Frequency {
        self.kla
    }
}

/// Numerical controls for the joint Au/O2 stocks and signed transfer account.
#[derive(Clone, Copy, Debug)]
pub struct AerationNumerics {
    /// Maximum RK4 trial duration.
    pub maximum_step: Time,
    /// Absolute time-weighted local-error sum allowance in Au-equivalent moles.
    pub absolute_amount_error: Amount,
    /// Additional fraction of initial Au plus four times all supplied O2.
    pub relative_amount_error: Dimensionless,
    /// Accepted and rejected trials, at most 65,536.
    pub maximum_attempts: usize,
    /// Shared stage/speciation work, at most ten million.
    pub maximum_work: usize,
}
impl Default for AerationNumerics {
    fn default() -> Self {
        Self {
            maximum_step: Time::new(1.0),
            absolute_amount_error: Amount::new(1e-12),
            relative_amount_error: Dimensionless::new(1e-7),
            maximum_attempts: 8192,
            maximum_work: 2_000_000,
        }
    }
}
impl AerationNumerics {
    fn tolerance(self, scale: f64) -> Result<f64> {
        // Same hard numerical bounds as the original solver; error semantics
        // above deliberately include oxygen, not only the original Au progress.
        Ok(CyanidationNumerics {
            maximum_step: self.maximum_step,
            absolute_gold_error: self.absolute_amount_error,
            relative_gold_error: self.relative_amount_error,
            maximum_attempts: self.maximum_attempts,
            maximum_work: self.maximum_work,
        }
        .tolerance(scale)?)
    }
}

/// Simultaneous chemical, gas and liquid state at one accepted absolute clock.
#[derive(Clone, Debug)]
pub struct AeratedCyanidationState {
    /// Original complete solution state, with actual free CN and pH.
    pub solution: CyanidationState,
    /// Oxygen remaining in the finite gas volume, mol O2.
    pub gas_oxygen: Amount,
    /// Ideal oxygen partial pressure; total gas pressure is not inferred.
    pub oxygen_partial_pressure: Pressure,
    /// Liquid oxygen amount in equilibrium with the CURRENT gas stock.
    pub equilibrium_dissolved_oxygen: Amount,
    /// Actual signed exchange rate, mol O2/s; zero when contact is disabled.
    pub oxygen_transfer_rate: OxygenTransferRate,
    /// Signed independently integrated net gas-to-liquid transfer since attachment.
    pub cumulative_oxygen_transfer: Amount,
    /// Original transfer-source conditions at this state.
    pub transfer_conditions: QueryPoint,
    /// Actual temperature/context query and original pinned Hcp receipt for the
    /// direct measured route. `None` means the retained Gibbs owner was used.
    pub henry_answer: Option<MaterialAnswer>,
    /// Explicit contact switch for this state; it changes no source coefficient.
    pub gas_contact_active: bool,
}
/// Recomputed independent balances; no finite solvent H/O or heat claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AerationConservation {
    /// Remaining solid Au plus dissolved complex Au minus original total.
    pub gold: Amount,
    /// Uncomplexed cyanide plus gold/competing-metal ligands minus original total.
    /// Retains any reported, tolerance-checked ligand-equilibrium residual.
    pub cyanide: Amount,
    /// Fixed canonical Na+ spectator minus original and delivered sodium.
    pub sodium: Amount,
    /// Liquid oxygen + reaction consumption - initial liquid oxygen - net transfer.
    pub liquid_oxygen: Amount,
    /// Gas oxygen + net transfer - initial gas oxygen - all headspace doses.
    pub gas_oxygen: Amount,
    /// Gas + liquid + original-equation consumption - all initial/supplied O2.
    pub total_oxygen: Amount,
}
/// A completed interval retaining the original chemical and transfer models.
#[derive(Clone, Debug)]
pub struct AerationAdvance {
    layout: Arc<GasLayout>,
    /// Accepted clock before this interval.
    pub previous_elapsed: Time,
    /// Accepted final simultaneous state and source points.
    pub state: AeratedCyanidationState,
    /// Gold dissolved during THIS interval, mol Au.
    pub interval_gold_dissolved: Amount,
    /// Signed gas-to-liquid exchange during THIS interval, mol O2.
    pub interval_oxygen_transfer: Amount,
    /// Actual metered O2 additions to the headspace since attachment.
    pub cumulative_oxygen_dose: Amount,
    /// Cumulative extent of the original written Elsner equation, mol.
    pub reaction_extent: Amount,
    /// Original equation's consumed/formed amounts, distinct from equilibrium OH.
    pub stoichiometric_amounts: Vec<ReactionAmount>,
    /// Independently recomputed Au, cyanide and three oxygen residuals.
    pub conservation: AerationConservation,
    /// Largest accepted local full/two-half estimate, Au-equivalent mol.
    pub maximum_local_error: Amount,
    /// Sum of accepted local estimates since attachment; not a global certificate.
    pub accumulated_error_estimate: Amount,
    /// Accepted two-half-step updates.
    pub accepted_steps: usize,
    /// Rejected inaccurate or depleting trials.
    pub rejected_steps: usize,
    /// Charged stage/speciation work, including unsuccessful trials.
    pub work: usize,
}
impl AerationAdvance {
    /// Cumulative accepted dissolved-reagent inputs at this accepted endpoint.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.chemical.reagents
    }
    /// Original reaction and source-backed chemical rate/acid-base owner.
    pub fn chemical_model(&self) -> &GoldCyanidation {
        &self.layout.chemical.model
    }
    /// Original phase-partition evidence and separate local transfer source.
    pub fn transfer_model(&self) -> &OxygenTransfer {
        &self.layout.transfer
    }
    /// Original liquid charge and independent activity basis.
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.chemical.charge
    }
    /// Fixed volume used to calculate the oxygen partial pressure.
    pub fn gas_volume(&self) -> Volume {
        self.layout.gas_volume
    }
    /// Oxygen initially attached in the headspace, before metered additions.
    pub fn initial_gas_oxygen(&self) -> Amount {
        self.layout.initial_gas
    }
}

#[derive(Clone, Debug)]
struct GasLayout {
    chemical: Layout,
    transfer: OxygenTransfer,
    gas_volume: Volume,
    initial_gas: Amount,
    conditions: QueryPoint,
}
/// Persistent coupled leach and finite oxygen headspace, with atomic updates.
#[derive(Clone, Debug)]
pub struct AeratedGoldCyanidation {
    layout: Arc<GasLayout>,
    state: AeratedCyanidationState,
    cumulative_dose: Amount,
    accumulated_error: f64,
}
impl AeratedGoldCyanidation {
    /// Attach a finite oxygen headspace to an admitted, possibly already advanced
    /// bath. Original charge, chemistry, absolute time and prior reaction progress
    /// are retained. The new transfer account starts at zero at attachment.
    pub fn from_batch(
        batch: GoldCyanidationBatch,
        transfer: OxygenTransfer,
        gas_oxygen: Amount,
        gas_volume: Volume,
        transfer_conditions: QueryPoint,
    ) -> Result<Self> {
        nonnegative(
            gas_oxygen.value(),
            "gas oxygen must be finite and nonnegative",
        )?;
        positive(gas_volume.value(), "headspace volume must be positive")?;
        if let Some(equilibrium) = transfer.equilibrium() {
            let aqueous_reference = batch.model().chemistry().water().standard();
            let gas_reference = equilibrium.standard();
            require(
                gas_reference.temperature() == aqueous_reference.temperature()
                    && gas_reference.reference_pressure() == aqueous_reference.reference_pressure()
                    && gas_reference.reference_axis() == aqueous_reference.reference_axis(),
                "oxygen partition and aqueous chemistry need a common T, pressure and reference",
            )?;
            ReactionNetwork::new(vec![
                batch.model().reaction().clone(),
                equilibrium.reaction().clone(),
            ])
            .map_err(CyanidationError::from)?;
        }
        require(
            transfer_axes()
                .iter()
                .all(|(axis, _, _)| !transfer_conditions.axes().contains_key(*axis))
                && owned_axes()
                    .iter()
                    .all(|(axis, _)| !transfer_conditions.axes().contains_key(*axis)),
            "caller cannot override actual gas/liquid source coordinates",
        )?;
        let mut chemical = (*batch.layout).clone();
        chemical.limit = chemical
            .charge
            .accessible_gold
            .value()
            .min(product(0.5, dosing::cyanide_total(&chemical)?)?);
        sum([chemical.charge.dissolved_oxygen.value(), gas_oxygen.value()].into_iter())?;
        let layout = Arc::new(GasLayout {
            chemical,
            transfer,
            gas_volume,
            initial_gas: gas_oxygen,
            conditions: transfer_conditions,
        });
        let mut never = || false;
        let mut work = Work {
            used: 0,
            maximum: 1_000_000,
            cancelled: &mut never,
            stopped: None,
        };
        let raw = [
            batch.state.gold_dissolved.value(),
            batch.state.dissolved_oxygen.value(),
            gas_oxygen.value(),
            0.0,
        ];
        let state = evaluate(&layout, batch.elapsed().value(), raw, true, &mut work)?;
        Ok(Self {
            layout,
            state,
            cumulative_dose: Amount::new(0.0),
            accumulated_error: 0.0,
        })
    }
    /// Original source-bound chemical rate and acid/base system.
    pub fn chemical_model(&self) -> &GoldCyanidation {
        &self.layout.chemical.model
    }
    /// Original partition and separate local transfer source.
    pub fn transfer_model(&self) -> &OxygenTransfer {
        &self.layout.transfer
    }
    /// Original liquid charge, including its independent activity basis.
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.chemical.charge
    }
    /// Cumulative reagent inputs remain separate from the original liquid charge.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.chemical.reagents
    }

    /// Add finite dissolved NaCN/NaOH equivalents while retaining gas inventory,
    /// oxygen-dose/transfer accounts, contact switch, clock and reaction progress.
    /// All original chemical/metal and transfer sources admit the post-dose state.
    pub fn dose_reagent(
        &mut self,
        reagent: dosing::DissolvedSodiumReagent,
        amount: Amount,
        maximum_work: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<dosing::ReagentDoseReceipt> {
        let mut work = dosing::work(maximum_work, &mut cancelled)?;
        let mut layout = (*self.layout).clone();
        layout.chemical = dosing::prepare(&layout.chemical, reagent, amount, &mut work)?;
        let state = evaluate(
            &layout,
            self.state.solution.elapsed.value(),
            coordinates(&self.state),
            self.state.gas_contact_active,
            &mut work,
        )?;
        let candidate = Self {
            layout: Arc::new(layout),
            state,
            cumulative_dose: self.cumulative_dose,
            accumulated_error: self.accumulated_error,
        };
        candidate.balances(&candidate.state, candidate.cumulative_dose.value())?;
        let receipt = dosing::receipt(
            Arc::new(candidate.layout.chemical.clone()),
            reagent,
            amount,
            &candidate.state.solution,
            &mut work,
        )?;
        *self = candidate;
        Ok(receipt)
    }
    /// Accepted chemical/gas/transfer state, immutable to callers.
    pub fn state(&self) -> &AeratedCyanidationState {
        &self.state
    }
    /// Cumulative, actually accepted oxygen additions to the gas phase.
    pub const fn cumulative_oxygen_dose(&self) -> Amount {
        self.cumulative_dose
    }
    /// Current overall stoichiometric ceiling, not a finite-time conversion promise.
    pub fn maximum_gold_dissolution(&self) -> Result<Amount> {
        Ok(Amount::new(
            self.layout
                .chemical
                .limit
                .min(product(4.0, self.oxygen_supplied()?)?),
        ))
    }
    fn oxygen_supplied(&self) -> Result<f64> {
        Ok(sum([
            self.initial_charge().dissolved_oxygen.value(),
            self.layout.initial_gas.value(),
            self.cumulative_dose.value(),
        ]
        .into_iter())?)
    }
    /// Meter an explicit amount of molecular O2 into the fixed gas volume at
    /// the admitted temperature. No gas purity, air flow or supply is inferred.
    /// The liquid/chemical state and clock remain continuous; source refusal
    /// or cancellation leaves both inventory and dose account untouched.
    pub fn dose_oxygen(
        &mut self,
        oxygen: Amount,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<()> {
        nonnegative(oxygen.value(), "oxygen dose must be finite and nonnegative")?;
        let mut work = Work {
            used: 0,
            maximum: 1_000_000,
            cancelled: &mut cancelled,
            stopped: None,
        };
        work.tick()?;
        let gas = sum([self.state.gas_oxygen.value(), oxygen.value()].into_iter())?;
        let dose = sum([self.cumulative_dose.value(), oxygen.value()].into_iter())?;
        require(
            oxygen.value() == 0.0
                || (gas > self.state.gas_oxygen.value() && dose > self.cumulative_dose.value()),
            "positive dose is below inventory resolution",
        )?;
        let mut raw = coordinates(&self.state);
        raw[2] = gas;
        let state = evaluate(
            &self.layout,
            self.state.solution.elapsed.value(),
            raw,
            self.state.gas_contact_active,
            &mut work,
        )?;
        self.balances(&state, dose)?;
        work.tick()?;
        self.state = state;
        self.cumulative_dose = Amount::new(dose);
        Ok(())
    }
    /// Advance all reaction, gas and liquid coordinates together. Disabling
    /// contact sets exchange to zero while dissolution can consume existing DO.
    pub fn advance(
        &mut self,
        duration: Time,
        gas_contact_active: bool,
        numerics: AerationNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AerationAdvance> {
        nonnegative(
            duration.value(),
            "advance duration must be finite and nonnegative",
        )?;
        let scale = sum([
            self.initial_charge().accessible_gold.value(),
            product(4.0, self.oxygen_supplied()?)?,
        ]
        .into_iter())?;
        let tolerance = numerics.tolerance(scale)?;
        let start = self.state.solution.elapsed.value();
        let end = positive_or_zero(start + duration.value(), "absolute clock overflows")?;
        require(
            duration.value() == 0.0 || end > start,
            "positive duration is below clock resolution",
        )?;
        let mut work = Work {
            used: 0,
            maximum: numerics.maximum_work,
            cancelled: &mut cancelled,
            stopped: None,
        };
        work.tick()?;
        let mut state = evaluate(
            &self.layout,
            start,
            coordinates(&self.state),
            gas_contact_active,
            &mut work,
        )?;
        let (mut accepted, mut rejected, mut largest, mut accumulated) =
            (0, 0, 0.0_f64, self.accumulated_error);
        let mut step = numerics.maximum_step.value().min(duration.value());
        while state.solution.elapsed.value() < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CyanidationError::Budget.into());
            }
            let now = state.solution.elapsed.value();
            let remaining = end - now;
            step = step.min(remaining);
            // Avoid leaving a final remainder too short to represent two half
            // steps. Split the real remaining interval before that tail arises.
            if step < remaining && remaining - step < 0.5 * step {
                step = 0.5 * remaining;
            }
            let finish = if step == remaining { end } else { now + step };
            let middle = now.midpoint(finish);
            if middle <= now || finish <= middle {
                return Err(CyanidationError::Arithmetic(
                    "positive half-step loses clock resolution",
                )
                .into());
            }
            let trial = (|| -> Result<_> {
                let full = rk4(&self.layout, &state, finish, gas_contact_active, &mut work)?;
                let first = rk4(&self.layout, &state, middle, gas_contact_active, &mut work)?;
                let refined = rk4(&self.layout, &first, finish, gas_contact_active, &mut work)?;
                Ok((full, refined))
            })();
            let (full, refined) = match trial {
                Err(GasTransferError::Cyanidation(CyanidationError::TrialDepletion)) => {
                    rejected += 1;
                    step = product(step, 0.5)?;
                    continue;
                }
                other => other?,
            };
            let a = coordinates(&full);
            let b = coordinates(&refined);
            let error = quotient(
                sum((0..4).map(|i| (a[i] - b[i]).abs() * if i == 0 { 1.0 } else { 4.0 }))?,
                15.0,
            )?;
            let allowance = product(tolerance, quotient(finish - now, duration.value())?)?;
            if error > allowance {
                rejected += 1;
                step = product(step, 0.5)?;
                continue;
            }
            state = refined;
            accepted += 1;
            largest = largest.max(error);
            accumulated = sum([accumulated, error].into_iter())?;
            if error < allowance / 32.0 {
                step = (step * 2.0).min(numerics.maximum_step.value());
            }
        }
        let extent = Amount::new(quotient(state.solution.gold_dissolved.value(), 4.0)?);
        let amounts = self.chemical_model().reaction.amounts(extent)?;
        let conservation = self.balances(&state, self.cumulative_dose.value())?;
        work.tick()?;
        let answer = AerationAdvance {
            layout: self.layout.clone(),
            previous_elapsed: self.state.solution.elapsed,
            interval_gold_dissolved: Amount::new(
                state.solution.gold_dissolved.value() - self.state.solution.gold_dissolved.value(),
            ),
            interval_oxygen_transfer: Amount::new(
                state.cumulative_oxygen_transfer.value()
                    - self.state.cumulative_oxygen_transfer.value(),
            ),
            cumulative_oxygen_dose: self.cumulative_dose,
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
        self.state = state;
        self.accumulated_error = accumulated;
        Ok(answer)
    }
    fn balances(&self, state: &AeratedCyanidationState, dose: f64) -> Result<AerationConservation> {
        let c = self.initial_charge();
        let s = &state.solution;
        let reaction = self
            .chemical_model()
            .reaction
            .amounts(Amount::new(quotient(s.gold_dissolved.value(), 4.0)?))?;
        let consumed = amount(&reaction, OXYGEN);
        let gold = sum([
            s.accessible_gold.value(),
            s.gold_complex.value(),
            -c.accessible_gold.value(),
            -c.gold_complex.value(),
        ]
        .into_iter())?;
        let (cyanide, sodium) = dosing::balances(&self.layout.chemical, s)?;
        let liquid = sum([
            s.dissolved_oxygen.value(),
            consumed,
            -c.dissolved_oxygen.value(),
            -state.cumulative_oxygen_transfer.value(),
        ]
        .into_iter())?;
        let gas = sum([
            state.gas_oxygen.value(),
            state.cumulative_oxygen_transfer.value(),
            -self.layout.initial_gas.value(),
            -dose,
        ]
        .into_iter())?;
        let total = sum([
            s.dissolved_oxygen.value(),
            state.gas_oxygen.value(),
            consumed,
            -c.dissolved_oxygen.value(),
            -self.layout.initial_gas.value(),
            -dose,
        ]
        .into_iter())?;
        let oxygen_scale = sum([
            c.dissolved_oxygen.value(),
            self.layout.initial_gas.value(),
            dose,
        ]
        .into_iter())?;
        for (residual, scale) in [
            (gold, c.accessible_gold.value().max(c.gold_complex.value())),
            (
                cyanide - ligand_balance_residual(s),
                dosing::cyanide_total(&self.layout.chemical)?
                    .max(product(2.0, c.gold_complex.value())?),
            ),
            (liquid, oxygen_scale),
            (gas, oxygen_scale),
            (total, oxygen_scale),
        ] {
            if residual.abs() > 4096.0 * f64::EPSILON * scale {
                return Err(CyanidationError::Arithmetic(
                    "coupled Au/cyanide/oxygen balance loses closure",
                )
                .into());
            }
        }
        Ok(AerationConservation {
            gold: Amount::new(gold),
            cyanide: Amount::new(cyanide),
            sodium: Amount::new(sodium),
            liquid_oxygen: Amount::new(liquid),
            gas_oxygen: Amount::new(gas),
            total_oxygen: Amount::new(total),
        })
    }
}

fn transfer_axes() -> [(&'static str, Dims, bool); 6] {
    [
        (OXYGEN_AXIS, AmountConcentration::DIMS, false),
        (PARTIAL_PRESSURE_AXIS, Pressure::DIMS, false),
        (TEMPERATURE_AXIS, Temperature::DIMS, true),
        (VOLUME_AXIS, Volume::DIMS, true),
        (GAS_VOLUME_AXIS, Volume::DIMS, true),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, true),
    ]
}
fn coordinates(state: &AeratedCyanidationState) -> [f64; 4] {
    [
        state.solution.gold_dissolved.value(),
        state.solution.dissolved_oxygen.value(),
        state.gas_oxygen.value(),
        state.cumulative_oxygen_transfer.value(),
    ]
}
fn evaluate(
    layout: &GasLayout,
    elapsed: f64,
    raw: [f64; 4],
    active: bool,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<AeratedCyanidationState> {
    work.tick()?;
    for &value in &raw[..3] {
        nonnegative(value, "negative or nonfinite coupled inventory")?;
    }
    require(raw[3].is_finite(), "nonfinite signed oxygen transfer")?;
    let solution = evaluate_with_oxygen(&layout.chemical, elapsed, raw[0], Some(raw[1]), work)?;
    let c = &layout.chemical.charge;
    let pressure = match &layout.transfer.partition {
        OxygenPartition::Thermodynamic(equilibrium) => {
            equilibrium.ideal_partial_pressure(Amount::new(raw[2]), layout.gas_volume)?
        }
        OxygenPartition::Henry(_) => Pressure::new(product(
            quotient(raw[2], layout.gas_volume.value())?,
            product(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, c.temperature.value())?,
        )?),
    };
    let values = [
        quotient(raw[1], c.solution_volume.value())?,
        pressure.value(),
        c.temperature.value(),
        c.solution_volume.value(),
        layout.gas_volume.value(),
        c.amount_at_unit_activity.value(),
    ];
    let mut point = layout.conditions.clone();
    // Known evolving chemistry coordinates retain the original solution
    // owner's actual values. A caller's fixed extra condition cannot stand in
    // for changing pH, free cyanide, gold loading or absolute elapsed time.
    for (axis, dims) in owned_axes() {
        point = point.with_quantity(
            axis,
            QuantitySpec::dimensional(dims),
            solution.conditions.axes()[axis],
        )?;
    }
    for ((axis, dims, _), value) in transfer_axes().into_iter().zip(values) {
        point = point.with_quantity(axis, QuantitySpec::dimensional(dims), value)?;
    }
    let (target, henry_answer) = match &layout.transfer.partition {
        OxygenPartition::Thermodynamic(equilibrium) => (
            equilibrium.equilibrium_aqueous_amount(pressure, c.amount_at_unit_activity)?,
            None,
        ),
        OxygenPartition::Henry(henry) => {
            let answer = henry.at_runtime(&point)?;
            // Hcp is amount concentration per pressure, so physical solution
            // volume converts to mol. The independent activity basis is NOT used.
            let dissolved = product(
                product(answer.evidence.value.value, pressure.value())?,
                c.solution_volume.value(),
            )?;
            (Amount::new(dissolved), Some(answer))
        }
    };
    for (axis, &(minimum, maximum)) in layout.transfer.source.validity.bounds() {
        let value = *point.axes().get(axis).ok_or(CyanidationError::Invalid(
            "missing transfer source condition",
        ))?;
        require(
            point.axis_quantities().get(axis)
                == layout.transfer.source.validity.axis_quantities().get(axis),
            "transfer source condition has wrong quantity",
        )?;
        if value < minimum || value > maximum {
            return Err(CyanidationError::OutsideValidity {
                axis: axis.clone(),
                value,
                minimum,
                maximum,
            }
            .into());
        }
    }
    let rate = if !active || layout.transfer.kla.value() == 0.0 || target.value() == raw[1] {
        0.0
    } else {
        let difference = target.value() - raw[1];
        product(layout.transfer.kla.value(), difference.abs())?.copysign(difference)
    };
    Ok(AeratedCyanidationState {
        solution,
        gas_oxygen: Amount::new(raw[2]),
        oxygen_partial_pressure: pressure,
        equilibrium_dissolved_oxygen: target,
        oxygen_transfer_rate: OxygenTransferRate::new(rate),
        cumulative_oxygen_transfer: Amount::new(raw[3]),
        transfer_conditions: point,
        henry_answer,
        gas_contact_active: active,
    })
}
fn derivatives(state: &AeratedCyanidationState) -> Result<[f64; 4]> {
    let r = state.solution.rate.value();
    let j = state.oxygen_transfer_rate.value();
    let consumption = quotient(r, 4.0)?;
    Ok([r, sum([j, -consumption].into_iter())?, -j, j])
}
fn update(initial: [f64; 4], h: f64, rate: [f64; 4], limit: f64) -> Result<[f64; 4]> {
    let mut next = initial;
    next[0] = progress(initial[0], h, rate[0], limit)?;
    for i in 1..4 {
        if rate[i] == 0.0 {
            continue;
        }
        let delta = h * rate[i];
        let value = initial[i] + delta;
        if !delta.is_finite() || !value.is_finite() {
            return Err(CyanidationError::Arithmetic("coupled oxygen update overflows").into());
        }
        if i < 3 && (value < 0.0 || (initial[i] > 0.0 && rate[i] < 0.0 && value == 0.0)) {
            return Err(CyanidationError::TrialDepletion.into());
        }
        if delta == 0.0 || value == initial[i] {
            return Err(CyanidationError::Arithmetic(
                "nonzero oxygen update is below inventory resolution",
            )
            .into());
        }
        next[i] = if value == 0.0 { 0.0 } else { value };
    }
    Ok(next)
}
fn rk4(
    layout: &GasLayout,
    initial: &AeratedCyanidationState,
    end: f64,
    active: bool,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<AeratedCyanidationState> {
    let start = initial.solution.elapsed.value();
    let h = end - start;
    let middle = start.midpoint(end);
    let y = coordinates(initial);
    let k1 = derivatives(initial)?;
    let stage2 = evaluate(
        layout,
        middle,
        update(y, 0.5 * h, k1, layout.chemical.limit)?,
        active,
        work,
    )?;
    let k2 = derivatives(&stage2)?;
    let stage3 = evaluate(
        layout,
        middle,
        update(y, 0.5 * h, k2, layout.chemical.limit)?,
        active,
        work,
    )?;
    let k3 = derivatives(&stage3)?;
    let stage4 = evaluate(
        layout,
        end,
        update(y, h, k3, layout.chemical.limit)?,
        active,
        work,
    )?;
    let k4 = derivatives(&stage4)?;
    let mut average = [0.0; 4];
    for i in 0..4 {
        average[i] = sum([k1[i] / 6.0, k2[i] / 3.0, k3[i] / 3.0, k4[i] / 6.0].into_iter())?;
    }
    evaluate(
        layout,
        end,
        update(y, h, average, layout.chemical.limit)?,
        active,
        work,
    )
}
