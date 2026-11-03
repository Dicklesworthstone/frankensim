//! Cyanidation in one full, isothermal, ideally mixed slurry tank.
//!
//! Every declared analytical component has the same residence: `dn/dt =
//! Q*C_in - Q*n/V + nu*r`. The local empirical rate is explicitly declared as
//! `r_Au = k*n_Au*c_CN_free*c_O2`; it is not inferred from a batch-age card.
//! Original acid/base and finite-metal equilibrium determine actual free CN
//! and pH at every reaction stage. Metal families, counterions and background
//! acid families travel with the slurry, including during feed changes.
//!
//! The ideal mole balance follows Fogler, chapter 1, equations 1-4/1-8:
//! <https://www.pearson.com/content/dam/one-dot-com/one-dot-com/hei-samples/science/9780135337554.pdf>.
//! The source supplies no ore rate, accessible-gold fraction or plant calibration.
//! Fixed liquid volume, temperature and solvent activity basis are explicit.
//! Flow and every feed concentration use that same declared solution-volume
//! basis; entrained solids share its residence. Solid displacement is not inferred.
//! Oxygen comes only from initial liquid and the feed. No gas exchange, carbon
//! adsorption, solids slip, dissolution of competing metals, variable volume,
//! solvent H/O balance or heat feedback is inferred.
//!
//! Exact transport / RK4 reaction / exact transport splitting is second order.
//! Full-step/two-half differences divided by three estimate errors in all
//! stocks, independently integrated exports and reaction progress. Source
//! admission covers every numerical split endpoint and reaction stage, not a
//! rigorous enclosure of the unsplit path. Local numerical estimates exclude
//! equilibrium/source/model uncertainty. No inventories or balances are repaired.

use super::*;
use fs_math::det;
use fs_qty::VolumetricFlowRate;
use fs_thermochem::stoichiometry::thermodynamics::acid_base::complexation::IdealLigandSystem;
use std::collections::BTreeMap;

/// Continuous leaching with a source-declared maintained oxygen reservoir.
pub mod aeration;

/// Continuous aerated leaching with explicitly retained carbon populations.
pub mod cil;

/// Independent local authority, distinct from the original batch-age law.
pub const LOCAL_CYANIDATION_LAW: &str = "local-speciated-gold-cyanidation";
const FORBIDDEN_AXES: [&str; 5] = [
    ELAPSED_TIME_AXIS,
    "batch-age",
    "residence-time",
    "mean-residence-time",
    "batch-volume",
];

/// Original chemical sources plus an independently authored local kinetic card.
#[derive(Clone, Debug)]
pub struct LocalGoldCyanidation {
    chemistry: IdealAcidBaseSystem,
    cyanide_family: String,
    source: ConstitutiveModelCard,
    reaction: ConservedReaction,
    rate: CyanidationRateConstant,
}
impl LocalGoldCyanidation {
    /// Require the local law, typed current-state support and one explicit k.
    /// The supplied card is retained unchanged; batch age is not reactor residence.
    pub fn from_card(
        chemistry: IdealAcidBaseSystem,
        cyanide_family: impl Into<String>,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        source.validate()?;
        require(
            source.law.0 == LOCAL_CYANIDATION_LAW
                && source.law_version == 1
                && source.state_schema_version == 0
                && source.initial_state == InitialStatePolicy::ZeroInternalState
                && source.parameters.len() == 1
                && !source.sources.is_empty()
                && !source.provenance.source.trim().is_empty()
                && !source.provenance.license.trim().is_empty(),
            "require an original version-one local cyanidation law with one coefficient",
        )?;
        let parameter = source
            .parameters
            .get(RATE_PARAMETER)
            .ok_or(CyanidationError::Invalid(
                "missing local cyanidation coefficient",
            ))?;
        require(
            parameter.dims == CyanidationRateConstant::DIMS,
            "local cyanidation coefficient must have m6/(mol2 s) dimensions",
        )?;
        nonnegative(
            parameter.value,
            "local cyanidation coefficient must be finite and nonnegative",
        )?;
        let rate = CyanidationRateConstant::new(parameter.value);
        for (axis, &(lo, hi)) in source.validity.bounds() {
            require(
                lo.is_finite()
                    && hi.is_finite()
                    && lo <= hi
                    && source.validity.axis_quantities().contains_key(axis)
                    && !FORBIDDEN_AXES.contains(&axis.as_str()),
                "local kinetic support must be finite and typed, not a batch or residence clock",
            )?;
        }
        for (i, (axis, dims)) in local_axes().into_iter().enumerate() {
            let Some((lo, hi)) = source.validity.bound(axis) else {
                require(i >= 6, "missing required local kinetic coordinate")?;
                continue;
            };
            require(
                source.validity.axis_quantities().get(axis)
                    == Some(&QuantitySpec::dimensional(dims))
                    && (axis == PH_AXIS || lo >= 0.0)
                    && (!matches!(
                        axis,
                        TEMPERATURE_AXIS | VOLUME_AXIS | ACTIVITY_REFERENCE_AXIS
                    ) || lo > 0.0)
                    && hi >= lo,
                "invalid local kinetic coordinate support",
            )?;
        }
        let cyanide_family = cyanide_family.into();
        let reaction = cyanidation_reaction(&chemistry, &cyanide_family)?;
        let (lo, hi) = source
            .validity
            .bound(TEMPERATURE_AXIS)
            .expect("required T support");
        require(
            (lo..=hi).contains(&chemistry.water().standard().temperature().value()),
            "local kinetics must cover the retained equilibrium temperature",
        )?;
        Ok(Self {
            chemistry,
            cyanide_family,
            source,
            reaction,
            rate,
        })
    }
    pub fn chemistry(&self) -> &IdealAcidBaseSystem {
        &self.chemistry
    }
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    pub fn cyanide_family(&self) -> &str {
        &self.cyanide_family
    }
    pub const fn rate_constant(&self) -> CyanidationRateConstant {
        self.rate
    }
}

/// One original fixed-charge species in an explicitly declared incoming stream.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectatorFeed {
    pub species: SpeciesId,
    pub charge: i32,
    pub concentration: AmountConcentration,
}
/// Constant incoming analytical composition during one accepted flow interval.
/// Every original family, spectator and competing-metal ID is required exactly
/// once, including zeros. Concentrations and flow use the declared solution-volume basis.
#[derive(Clone, Debug, PartialEq)]
pub struct CyanidationFeed {
    pub accessible_gold: AmountConcentration,
    /// All ligand outside the separately declared Au(CN)2 stock, before binding.
    pub cyanide_ligand: AmountConcentration,
    pub dissolved_oxygen: AmountConcentration,
    pub gold_complex: AmountConcentration,
    pub other_families: Vec<(String, AmountConcentration)>,
    pub spectators: Vec<SpectatorFeed>,
    pub competing_metals: Vec<(String, AmountConcentration)>,
}
/// Extensive analytical stocks or stream amounts, in canonical component order.
/// `cyanide_ligand` includes HCN, CN and competing-metal ligands, excluding gold ligands.
#[derive(Clone, Debug, PartialEq)]
pub struct CyanidationInventory {
    pub accessible_gold: Amount,
    pub cyanide_ligand: Amount,
    pub dissolved_oxygen: Amount,
    pub gold_complex: Amount,
    pub other_families: Vec<(String, Amount)>,
    pub spectators: Vec<SpectatorIon>,
    pub competing_metals: Vec<(String, Amount)>,
}
/// Actual tank solution; cumulative reaction is separate from remaining gold complex.
#[derive(Clone, Debug)]
pub struct CyanidationFlowState {
    pub elapsed: Time,
    pub inventory: CyanidationInventory,
    pub uncomplexed_cyanide: Amount,
    pub free_cyanide: Amount,
    pub speciation: AcidBaseSpeciation,
    pub competition: Option<LigandSpeciation>,
    pub conditions: QueryPoint,
    pub rate: GoldDissolutionRate,
}
/// Independent material balances. Residual fields may be signed.
#[derive(Clone, Debug, PartialEq)]
pub struct FlowConservation {
    pub gold: Amount,
    /// Includes the original equilibrium owner's reported ligand residual.
    pub cyanide: Amount,
    pub dissolved_oxygen: Amount,
    /// Fixed canonical Na+ spectator only; no new elemental speciation is inferred.
    pub sodium: Amount,
    /// Final + outlet - initial - inlet - reaction contribution for every component.
    pub components: CyanidationInventory,
}
/// Shared bounded integration controls, distinct from empirical uncertainty.
#[derive(Clone, Copy, Debug)]
pub struct FlowNumerics {
    pub maximum_step: Time,
    /// L1 amount allowance over stocks, exports and gold reaction progress.
    pub absolute_amount_error: Amount,
    pub relative_amount_error: Dimensionless,
    pub maximum_attempts: usize,
    pub maximum_work: usize,
}
impl Default for FlowNumerics {
    fn default() -> Self {
        Self {
            maximum_step: Time::new(1.0),
            absolute_amount_error: Amount::new(1e-10),
            relative_amount_error: Dimensionless::new(1e-6),
            maximum_attempts: 8192,
            maximum_work: 2_000_000,
        }
    }
}
#[derive(Clone, Debug)]
struct FlowLayout {
    model: LocalGoldCyanidation,
    original_charge: CyanidationCharge,
    charge: CyanidationCharge,
    conditions: QueryPoint,
    speciation: SpeciationOptions,
    competition: Option<competition::BoundCompetition>,
    initial: Vec<f64>,
}
/// Persistent mixed-tank inventory with independently accumulated incoming and outgoing streams.
#[derive(Clone, Debug)]
pub struct CyanidationCstr {
    layout: Arc<FlowLayout>,
    raw: Vec<f64>,
    state: CyanidationFlowState,
    inlet: Vec<f64>,
    outlet: Vec<f64>,
    gold_dissolved: f64,
    accumulated_error: f64,
}
// Internal source callbacks keep their original error type across the liquid
// owner. Public standalone calls still expose exactly CyanidationError.
enum FlowAdvanceError<E> {
    Cyanidation(CyanidationError),
    Admission(E),
}
impl<E> From<CyanidationError> for FlowAdvanceError<E> {
    fn from(error: CyanidationError) -> Self {
        Self::Cyanidation(error)
    }
}

/// One accepted interval. No partially accepted trial or failed dose is published.
#[derive(Clone, Debug)]
pub struct CyanidationFlowAdvance {
    layout: Arc<FlowLayout>,
    pub previous_elapsed: Time,
    pub elapsed: Time,
    pub flow: VolumetricFlowRate,
    pub feed: CyanidationFeed,
    pub state: CyanidationFlowState,
    pub inlet: CyanidationInventory,
    pub outlet: CyanidationInventory,
    pub cumulative_inlet: CyanidationInventory,
    pub cumulative_outlet: CyanidationInventory,
    pub interval_gold_dissolved: Amount,
    pub cumulative_gold_dissolved: Amount,
    /// Original written equation at cumulative Au progress divided by four.
    pub stoichiometric_amounts: Vec<ReactionAmount>,
    pub interval_conservation: FlowConservation,
    pub cumulative_conservation: FlowConservation,
    pub maximum_local_error: Amount,
    pub accumulated_error_estimate: Amount,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub work: usize,
}
impl CyanidationFlowAdvance {
    pub fn model(&self) -> &LocalGoldCyanidation {
        &self.layout.model
    }
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.original_charge
    }
    /// Original finite-metal equilibrium families and their retained source receipts.
    pub fn competing_chemistry(&self) -> Option<&IdealLigandSystem> {
        self.layout.competition.as_ref().map(|metals| &metals.model)
    }
}
impl CyanidationCstr {
    pub fn new(
        model: LocalGoldCyanidation,
        charge: CyanidationCharge,
        conditions: QueryPoint,
        speciation: SpeciationOptions,
    ) -> Result<Self> {
        Self::new_inner(model, charge, conditions, speciation, None)
    }
    pub fn new_with_competing_metals(
        model: LocalGoldCyanidation,
        charge: CyanidationCharge,
        conditions: QueryPoint,
        metals: competition::CompetingMetals,
    ) -> Result<Self> {
        require(
            !metals.families.is_empty(),
            "competing metals require declared association families",
        )?;
        require(
            metals.totals.len() < 256,
            "flow metal component count is bounded",
        )?;
        require(
            metals.families.iter().all(|family| {
                family.species().iter().all(|species| {
                    !model
                        .reaction
                        .terms()
                        .iter()
                        .any(|term| &term.species == species)
                })
            }),
            "competing metals cannot own original leach stocks",
        )?;
        let chemistry = IdealLigandSystem::new(
            model.chemistry.clone(),
            model.cyanide_family.clone(),
            metals.families,
        )?;
        let mut totals = metals.totals;
        totals.sort_by(|a, b| a.0.cmp(&b.0));
        let speciation = metals.numerics.aqueous;
        Self::new_inner(
            model,
            charge,
            conditions,
            speciation,
            Some(competition::BoundCompetition {
                model: chemistry,
                totals,
                numerics: metals.numerics,
            }),
        )
    }
    fn new_inner(
        model: LocalGoldCyanidation,
        mut charge: CyanidationCharge,
        conditions: QueryPoint,
        speciation: SpeciationOptions,
        competition: Option<competition::BoundCompetition>,
    ) -> Result<Self> {
        positive(
            charge.solution_volume.value(),
            "flowing solution volume must be positive",
        )?;
        positive(
            charge.amount_at_unit_activity.value(),
            "flowing activity basis must be positive",
        )?;
        require(
            charge.temperature == model.chemistry.water().standard().temperature(),
            "flowing bath temperature must equal the original equilibrium temperature",
        )?;
        require(
            conditions.axes().keys().all(|axis| {
                !FORBIDDEN_AXES.contains(&axis.as_str())
                    && local_axes().iter().all(|(owned, _)| axis != owned)
                    && model.source.validity.bound(axis).is_some()
            }),
            "only additional original local source coordinates belong in flow conditions",
        )?;
        require(
            charge.spectators.len() < 256 && charge.other_families.len() < 256,
            "flow background component count is bounded",
        )?;
        require(
            charge.spectators.iter().all(|ion| {
                !model
                    .reaction
                    .terms()
                    .iter()
                    .any(|term| term.species == ion.species)
            }),
            "flow spectators cannot duplicate original leach stocks",
        )?;
        require(
            charge
                .other_families
                .iter()
                .all(|(id, _)| id != &model.cyanide_family),
            "flowing cyanide analytical total is separately owned",
        )?;
        let original_charge = charge.clone();
        charge.other_families.sort_by(|a, b| a.0.cmp(&b.0));
        charge.spectators.sort_by(|a, b| a.species.cmp(&b.species));
        let mut raw = vec![
            charge.accessible_gold.value(),
            charge.uncomplexed_cyanide.value(),
            charge.dissolved_oxygen.value(),
            charge.gold_complex.value(),
        ];
        raw.extend(charge.other_families.iter().map(|(_, n)| n.value()));
        raw.extend(charge.spectators.iter().map(|ion| ion.amount.value()));
        if let Some(metals) = &competition {
            raw.extend(metals.totals.iter().map(|(_, n)| n.value()));
        }
        for value in &raw {
            nonnegative(*value, "flow inventories must be finite and nonnegative")?;
        }
        sum(raw.iter().copied())?;
        let layout = Arc::new(FlowLayout {
            model,
            original_charge,
            charge,
            conditions,
            speciation,
            competition,
            initial: raw.clone(),
        });
        let mut never = || false;
        let mut work = Work {
            used: 0,
            maximum: 1_000_000,
            cancelled: &mut never,
            stopped: None,
        };
        let state = evaluate(&layout, &raw, 0.0, &mut work)?;
        Ok(Self {
            inlet: vec![0.0; raw.len()],
            outlet: vec![0.0; raw.len()],
            layout,
            raw,
            state,
            gold_dissolved: 0.0,
            accumulated_error: 0.0,
        })
    }
    pub fn model(&self) -> &LocalGoldCyanidation {
        &self.layout.model
    }
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.original_charge
    }
    pub fn state(&self) -> &CyanidationFlowState {
        &self.state
    }
    pub fn elapsed(&self) -> Time {
        self.state.elapsed
    }
    pub fn cumulative_inlet(&self) -> CyanidationInventory {
        inventory(&self.layout, &self.inlet)
    }
    pub fn cumulative_outlet(&self) -> CyanidationInventory {
        inventory(&self.layout, &self.outlet)
    }
    pub fn cumulative_gold_dissolved(&self) -> Amount {
        Amount::new(self.gold_dissolved)
    }
    pub fn competing_chemistry(&self) -> Option<&IdealLigandSystem> {
        self.layout.competition.as_ref().map(|metals| &metals.model)
    }

    /// Advance one constant feed/flow interval atomically. A stopped flow leaves
    /// chemical dissolution active; no reactor clock is promoted to source age.
    pub fn advance(
        &mut self,
        duration: Time,
        flow: VolumetricFlowRate,
        feed: &CyanidationFeed,
        numerics: FlowNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CyanidationFlowAdvance> {
        self.advance_admitted(
            duration,
            flow,
            feed,
            numerics,
            &mut cancelled,
            &mut |_| Ok::<(), CyanidationError>(()),
        )
        .map_err(|error| match error {
            FlowAdvanceError::Cyanidation(error) | FlowAdvanceError::Admission(error) => error,
        })
    }

    // The liquid owner remains responsible for its own interval and cumulative
    // accounts. Coupled owners may add source admission without replacing the
    // accepted transport/reaction algorithm or erasing typed source errors.
    fn advance_admitted<E>(
        &mut self,
        duration: Time,
        flow: VolumetricFlowRate,
        feed: &CyanidationFeed,
        numerics: FlowNumerics,
        mut cancelled: &mut dyn FnMut() -> bool,
        admit: &mut dyn FnMut(&CyanidationFlowState) -> core::result::Result<(), E>,
    ) -> core::result::Result<CyanidationFlowAdvance, FlowAdvanceError<E>> {
        nonnegative(
            duration.value(),
            "flow duration must be finite and nonnegative",
        )?;
        nonnegative(flow.value(), "slurry flow must be finite and nonnegative")?;
        let start = self.elapsed().value();
        let end = add(start, duration.value())?;
        let mut work = dosing::work(numerics.maximum_work, &mut cancelled)?;
        let concentration = feed_values(&self.layout, feed, &mut work)?;
        let mut source = |state: &CyanidationFlowState| {
            admit(state).map_err(FlowAdvanceError::Admission)
        };
        source(&self.state)?;
        let volume = self.layout.charge.solution_volume.value();
        let exposure = quotient(product(flow.value(), duration.value())?, volume)?;
        let target = concentration
            .iter()
            .map(|c| product(*c, volume))
            .collect::<Result<Vec<_>>>()?;
        let incoming = target
            .iter()
            .map(|n| product(*n, exposure))
            .collect::<Result<Vec<_>>>()?;
        let scale = sum(self.raw.iter().chain(&incoming).copied())?;
        let tolerance = CyanidationNumerics {
            maximum_step: numerics.maximum_step,
            absolute_gold_error: numerics.absolute_amount_error,
            relative_gold_error: numerics.relative_amount_error,
            maximum_attempts: numerics.maximum_attempts,
            maximum_work: numerics.maximum_work,
        }
        .tolerance(scale)?;
        let problem = Problem {
            layout: &self.layout,
            target: &target,
            flow,
        };
        let mut current = Trial {
            raw: self.raw.clone(),
            outlet: vec![0.0; self.raw.len()],
            gold: 0.0,
        };
        let (mut clock, mut step) = (start, numerics.maximum_step.value().min(duration.value()));
        let (mut accepted, mut rejected) = (0, 0);
        let (mut error_sum, mut largest) = (0.0_f64, 0.0_f64);
        // The transport-only case is exact; it needs no adaptive reaction trials.
        if self.model().rate.value() == 0.0 && duration.value() > 0.0 {
            problem.transport_admitted(&mut current, duration.value(), end, &mut work, &mut source)?;
            clock = end;
            accepted = 1;
        }
        while clock < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CyanidationError::Budget.into());
            }
            let remaining = end - clock;
            // Integrate the last two steps as a pair, without leaving a
            // floating-point remainder too short for chemistry or transport.
            let final_pair = remaining > step && remaining - step < step;
            let finish = if remaining <= step {
                end
            } else if final_pair {
                clock.midpoint(end)
            } else {
                clock + step
            };
            let middle = clock.midpoint(finish);
            let h = finish - clock;
            if h <= 0.0 || middle <= clock || middle >= finish {
                return Err(CyanidationError::Arithmetic(
                    "flow integration clock loses resolution",
                ).into());
            }
            let initial = Trial {
                raw: current.raw.clone(),
                outlet: vec![0.0; current.raw.len()],
                gold: 0.0,
            };
            let trial = (|| {
                let coarse = problem.split_admitted(initial.clone(), clock, finish, &mut work, &mut source)?;
                let first = problem.split_admitted(initial, clock, middle, &mut work, &mut source)?;
                let fine = problem.split_admitted(first, middle, finish, &mut work, &mut source)?;
                Ok((coarse, fine))
            })();
            let (coarse, fine) = match trial {
                Ok(pair) => pair,
                Err(FlowAdvanceError::Cyanidation(
                    CyanidationError::TrialDepletion | CyanidationError::OutsideValidity { .. },
                )) => {
                    rejected += 1;
                    step = h * 0.5;
                    continue;
                }
                Err(error) => return Err(error),
            };
            let difference = sum(coarse
                .values()
                .zip(fine.values())
                .map(|(a, b)| (a - b).abs()))?;
            let error = quotient(difference, 3.0)?;
            let allowance = product(tolerance, quotient(h, duration.value())?)?;
            let next_error = sum([error_sum, error].into_iter())?;
            if error > allowance || next_error > tolerance {
                rejected += 1;
                step = h * 0.5;
                continue;
            }
            current.raw = fine.raw;
            for (total, increment) in current.outlet.iter_mut().zip(fine.outlet) {
                *total = add(*total, increment)?;
            }
            current.gold = add(current.gold, fine.gold)?;
            largest = largest.max(error);
            error_sum = next_error;
            accepted += 1;
            clock = finish;
            step = if final_pair {
                (end - finish).min(numerics.maximum_step.value())
            } else if error <= allowance / 8.0 {
                (h * 2.0).min(numerics.maximum_step.value())
            } else {
                h
            };
        }
        let state = evaluate(&self.layout, &current.raw, end, &mut work)?;
        source(&state)?;
        let cumulative_inlet = self
            .inlet
            .iter()
            .zip(&incoming)
            .map(|(a, b)| add(*a, *b))
            .collect::<Result<Vec<_>>>()?;
        let cumulative_outlet = self
            .outlet
            .iter()
            .zip(&current.outlet)
            .map(|(a, b)| add(*a, *b))
            .collect::<Result<Vec<_>>>()?;
        let gold = add(self.gold_dissolved, current.gold)?;
        let accumulated_error = sum([self.accumulated_error, error_sum].into_iter())?;
        let interval_conservation = conservation(
            &self.layout,
            &state,
            &self.raw,
            &incoming,
            &current.outlet,
            current.gold,
        )?;
        let cumulative_conservation = conservation(
            &self.layout,
            &state,
            &self.layout.initial,
            &cumulative_inlet,
            &cumulative_outlet,
            gold,
        )?;
        let stoichiometric_amounts = self
            .model()
            .reaction
            .amounts(Amount::new(quotient(gold, 4.0)?))
            .map_err(CyanidationError::from)?;
        let mut answer = CyanidationFlowAdvance {
            layout: self.layout.clone(),
            previous_elapsed: Time::new(start),
            elapsed: Time::new(end),
            flow,
            feed: feed.clone(),
            state: state.clone(),
            inlet: inventory(&self.layout, &incoming),
            outlet: inventory(&self.layout, &current.outlet),
            cumulative_inlet: inventory(&self.layout, &cumulative_inlet),
            cumulative_outlet: inventory(&self.layout, &cumulative_outlet),
            interval_gold_dissolved: Amount::new(current.gold),
            cumulative_gold_dissolved: Amount::new(gold),
            stoichiometric_amounts,
            interval_conservation,
            cumulative_conservation,
            maximum_local_error: Amount::new(largest),
            accumulated_error_estimate: Amount::new(accumulated_error),
            accepted_steps: accepted,
            rejected_steps: rejected,
            work: 0,
        };
        work.tick()?;
        answer.work = work.used;
        self.raw = current.raw;
        self.state = state;
        self.inlet = cumulative_inlet;
        self.outlet = cumulative_outlet;
        self.gold_dissolved = gold;
        self.accumulated_error = accumulated_error;
        Ok(answer)
    }
}

fn local_axes() -> Vec<(&'static str, Dims)> {
    owned_axes()
        .into_iter()
        .filter(|(axis, _)| *axis != ELAPSED_TIME_AXIS)
        .collect()
}
fn inventory(layout: &FlowLayout, raw: &[f64]) -> CyanidationInventory {
    let mut i = 4;
    let other_families = layout
        .charge
        .other_families
        .iter()
        .map(|(id, _)| {
            let value = (id.clone(), Amount::new(raw[i]));
            i += 1;
            value
        })
        .collect();
    let spectators = layout
        .charge
        .spectators
        .iter()
        .map(|ion| {
            let value = SpectatorIon {
                species: ion.species.clone(),
                charge: ion.charge,
                amount: Amount::new(raw[i]),
            };
            i += 1;
            value
        })
        .collect();
    let competing_metals = layout.competition.as_ref().map_or_else(Vec::new, |metals| {
        metals
            .totals
            .iter()
            .map(|(id, _)| {
                let value = (id.clone(), Amount::new(raw[i]));
                i += 1;
                value
            })
            .collect()
    });
    CyanidationInventory {
        accessible_gold: Amount::new(raw[0]),
        cyanide_ligand: Amount::new(raw[1]),
        dissolved_oxygen: Amount::new(raw[2]),
        gold_complex: Amount::new(raw[3]),
        other_families,
        spectators,
        competing_metals,
    }
}
fn feed_values(
    layout: &FlowLayout,
    feed: &CyanidationFeed,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<Vec<f64>> {
    work.tick()?;
    let ordered = |values: &[(String, AmountConcentration)],
                   required: &[(String, Amount)]|
     -> Result<Vec<f64>> {
        require(
            values.len() == required.len(),
            "feed analytical family count must match the original layout",
        )?;
        let map: BTreeMap<_, _> = values.iter().map(|(id, c)| (id, c)).collect();
        require(
            map.len() == values.len() && values.len() == required.len(),
            "flow feed needs every original analytical family exactly once",
        )?;
        required
            .iter()
            .map(|(id, _)| {
                map.get(id)
                    .map(|c| c.value())
                    .ok_or(CyanidationError::Invalid(
                        "missing original feed family identity",
                    ))
            })
            .collect()
    };
    let mut values = vec![
        feed.accessible_gold.value(),
        feed.cyanide_ligand.value(),
        feed.dissolved_oxygen.value(),
        feed.gold_complex.value(),
    ];
    values.extend(ordered(
        &feed.other_families,
        &layout.charge.other_families,
    )?);
    require(
        feed.spectators.len() == layout.charge.spectators.len(),
        "feed spectator count must match the original layout",
    )?;
    let spectators: BTreeMap<_, _> = feed
        .spectators
        .iter()
        .map(|ion| (&ion.species, ion))
        .collect();
    require(
        spectators.len() == feed.spectators.len()
            && feed.spectators.len() == layout.charge.spectators.len(),
        "flow feed needs every original spectator exactly once",
    )?;
    for original in &layout.charge.spectators {
        let ion = spectators
            .get(&original.species)
            .ok_or(CyanidationError::Invalid("missing original feed spectator"))?;
        require(
            ion.charge == original.charge,
            "feed spectator charge cannot change its original identity",
        )?;
        values.push(ion.concentration.value());
    }
    values.extend(ordered(
        &feed.competing_metals,
        layout
            .competition
            .as_ref()
            .map_or(&[], |m| m.totals.as_slice()),
    )?);
    for value in &values {
        work.tick()?;
        nonnegative(*value, "feed concentrations must be finite and nonnegative")?;
    }
    Ok(values)
}
fn evaluate(
    layout: &FlowLayout,
    raw: &[f64],
    elapsed: f64,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CyanidationFlowState> {
    work.tick()?;
    let inventory = inventory(layout, raw);
    let eq = speciate_current(
        &layout.model.chemistry,
        &layout.model.cyanide_family,
        layout.competition.as_ref(),
        CurrentSolution {
            cyanide: inventory.cyanide_ligand,
            gold_complex: inventory.gold_complex,
            other_families: &inventory.other_families,
            spectators: &inventory.spectators,
            metal_totals: &inventory.competing_metals,
            amount_at_unit_activity: layout.charge.amount_at_unit_activity,
        },
        layout.speciation,
        work,
    )?;
    let cn = quotient(
        eq.free_cyanide.value(),
        layout.charge.solution_volume.value(),
    )?;
    let oxygen = quotient(raw[2], layout.charge.solution_volume.value())?;
    let values = [
        eq.speciation.ph.value(),
        cn,
        oxygen,
        raw[0],
        layout.charge.temperature.value(),
        layout.charge.solution_volume.value(),
        eq.uncomplexed_cyanide.value(),
        raw[3],
        layout.charge.amount_at_unit_activity.value(),
    ];
    let mut point = layout.conditions.clone();
    for ((axis, dims), value) in local_axes().into_iter().zip(values) {
        point = point.with_quantity(axis, QuantitySpec::dimensional(dims), value)?;
    }
    for (axis, &(minimum, maximum)) in layout.model.source.validity.bounds() {
        let value = *point.axes().get(axis).ok_or(CyanidationError::Invalid(
            "missing original local source coordinate",
        ))?;
        require(
            point.axis_quantities().get(axis)
                == layout.model.source.validity.axis_quantities().get(axis),
            "local source coordinate has the wrong quantity",
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
    let rate = if layout.model.rate.value() == 0.0 || raw[0] == 0.0 || cn == 0.0 || oxygen == 0.0 {
        0.0
    } else {
        product(
            product(product(layout.model.rate.value(), raw[0])?, cn)?,
            oxygen,
        )?
    };
    Ok(CyanidationFlowState {
        elapsed: Time::new(elapsed),
        inventory,
        uncomplexed_cyanide: eq.uncomplexed_cyanide,
        free_cyanide: eq.free_cyanide,
        speciation: eq.speciation,
        competition: eq.competition,
        conditions: point,
        rate: GoldDissolutionRate::new(rate),
    })
}
#[derive(Clone)]
struct Trial {
    raw: Vec<f64>,
    outlet: Vec<f64>,
    gold: f64,
}
impl Trial {
    fn values(&self) -> impl Iterator<Item = f64> + '_ {
        self.raw
            .iter()
            .chain(&self.outlet)
            .copied()
            .chain(std::iter::once(self.gold))
    }
}
struct Problem<'a> {
    layout: &'a FlowLayout,
    target: &'a [f64],
    flow: VolumetricFlowRate,
}
impl Problem<'_> {
    fn transport(
        &self,
        state: &mut Trial,
        h: f64,
        time: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<()> {
        self.transport_admitted(state, h, time, work, &mut |_| Ok(()))
    }
    fn transport_admitted<E: From<CyanidationError>>(
        &self,
        state: &mut Trial,
        h: f64,
        time: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
        admit: &mut impl FnMut(&CyanidationFlowState) -> core::result::Result<(), E>,
    ) -> core::result::Result<(), E> {
        if self.flow.value() == 0.0 || h == 0.0 {
            return Ok(());
        }
        let x = quotient(
            product(self.flow.value(), h)?,
            self.layout.charge.solution_volume.value(),
        )?;
        let e = det::exp(-x);
        let u = -det::expm1(-x);
        // Integral of freshly supplied stock, without cancellation for small exposure.
        let w = if x < 0.125 {
            let mut term = product(x, x)? * 0.5;
            let mut value = term;
            for n in 3..=16 {
                term *= -x / n as f64;
                value += term;
            }
            positive_or_zero(value, "transport integral loses range")?
        } else {
            positive_or_zero(x - u, "transport integral loses range")?
        };
        for i in 0..state.raw.len() {
            work.tick()?;
            let before = state.raw[i];
            let target = self.target[i];
            let (after, exported) = if before == target {
                (before, product(before, x)?)
            } else {
                if before > 0.0 && e == 0.0 {
                    return Err(CyanidationError::Arithmetic(
                        "finite transport survival is below inventory resolution",
                    )
                    .into());
                }
                // Form the convex update around its nearer endpoint. Separately
                // weighting two almost equal stocks can create artificial drift.
                // A sub-ulp net change near equilibrium may round to the old
                // stock while gross inlet and outlet transport remain resolved.
                let after = if u <= 0.5 {
                    let change = product((target - before).abs(), u)?;
                    sum([before, if target > before { change } else { -change }].into_iter())?
                } else {
                    let change = product((before - target).abs(), e)?;
                    sum([target, if before > target { change } else { -change }].into_iter())?
                };
                (
                    after,
                    sum([product(before, u)?, product(target, w)?].into_iter())?,
                )
            };
            state.raw[i] = after;
            state.outlet[i] = add(state.outlet[i], exported)?;
        }
        admit(&evaluate(self.layout, &state.raw, time, work)?)?;
        Ok(())
    }
    fn reacted(&self, initial: &[f64], gold: f64) -> Result<Vec<f64>> {
        let amounts = self
            .layout
            .model
            .reaction
            .amounts(Amount::new(quotient(gold, 4.0)?))?;
        let mut raw = initial.to_vec();
        for (i, species) in [(0, AU), (1, CN), (2, OXYGEN)] {
            let consumed = amount(&amounts, species);
            raw[i] = remainder(initial[i], consumed, gold)?;
            if consumed > 0.0 && raw[i] == initial[i] {
                return Err(CyanidationError::Arithmetic(
                    "positive leach consumption is below stock resolution",
                ));
            }
        }
        raw[3] = add(initial[3], amount(&amounts, COMPLEX))?;
        Ok(raw)
    }
    fn reaction(
        &self,
        state: &mut Trial,
        start: f64,
        end: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<()> {
        self.reaction_admitted(state, start, end, work, &mut |_| Ok(()))
    }
    fn reaction_admitted<E: From<CyanidationError>>(
        &self,
        state: &mut Trial,
        start: f64,
        end: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
        admit: &mut impl FnMut(&CyanidationFlowState) -> core::result::Result<(), E>,
    ) -> core::result::Result<(), E> {
        let h = end - start;
        let middle = start.midpoint(end);
        let mut rate =
            |gold: f64, time: f64, work: &mut Work<'_, _>| -> core::result::Result<f64, E> {
                let solution = evaluate(self.layout, &self.reacted(&state.raw, gold)?, time, work)?;
                admit(&solution)?;
                Ok(solution.rate.value())
            };
        let a = rate(0.0, start, work)?;
        let b = rate(product(product(h, 0.5)?, a)?, middle, work)?;
        let c = rate(product(product(h, 0.5)?, b)?, middle, work)?;
        let d = rate(product(h, c)?, end, work)?;
        let average = sum([
            quotient(a, 6.0)?,
            quotient(b, 3.0)?,
            quotient(c, 3.0)?,
            quotient(d, 6.0)?,
        ]
        .into_iter())?;
        let gold = product(h, average)?;
        state.raw = self.reacted(&state.raw, gold)?;
        admit(&evaluate(self.layout, &state.raw, end, work)?)?;
        state.gold = add(state.gold, gold)?;
        Ok(())
    }
    fn split_admitted<E: From<CyanidationError>>(
        &self,
        mut state: Trial,
        start: f64,
        end: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
        admit: &mut impl FnMut(&CyanidationFlowState) -> core::result::Result<(), E>,
    ) -> core::result::Result<Trial, E> {
        let h = end - start;
        self.transport_admitted(&mut state, product(h, 0.5)?, start, work, admit)?;
        self.reaction_admitted(&mut state, start, end, work, admit)?;
        self.transport_admitted(&mut state, product(h, 0.5)?, end, work, admit)?;
        Ok(state)
    }
}
fn conservation(
    layout: &FlowLayout,
    state: &CyanidationFlowState,
    initial: &[f64],
    incoming: &[f64],
    outgoing: &[f64],
    gold: f64,
) -> Result<FlowConservation> {
    let current = &state.inventory;
    let mut final_raw = vec![
        current.accessible_gold.value(),
        current.cyanide_ligand.value(),
        current.dissolved_oxygen.value(),
        current.gold_complex.value(),
    ];
    final_raw.extend(current.other_families.iter().map(|(_, n)| n.value()));
    final_raw.extend(current.spectators.iter().map(|ion| ion.amount.value()));
    final_raw.extend(current.competing_metals.iter().map(|(_, n)| n.value()));
    let amounts = layout
        .model
        .reaction
        .amounts(Amount::new(quotient(gold, 4.0)?))?;
    let mut delta = vec![0.0; initial.len()];
    for (i, species) in [(0, AU), (1, CN), (2, OXYGEN)] {
        delta[i] = -amount(&amounts, species);
    }
    delta[3] = amount(&amounts, COMPLEX);
    let mut residuals = Vec::with_capacity(initial.len());
    for i in 0..initial.len() {
        let residual = sum([
            final_raw[i],
            outgoing[i],
            -initial[i],
            -incoming[i],
            -delta[i],
        ]
        .into_iter())?;
        let scale = sum([initial[i], incoming[i], delta[i].abs()].into_iter())?;
        require_balance(residual, scale)?;
        residuals.push(residual);
    }
    let ligand_residual = state
        .competition
        .as_ref()
        .map_or(0.0, |s| s.ligand_balance_residual.value());
    let bound = state
        .competition
        .as_ref()
        .map_or(0.0, |s| s.bound_ligand.value());
    let actual_ligand_residual =
        sum([state.uncomplexed_cyanide.value(), bound, -final_raw[1]].into_iter())?;
    require_balance(actual_ligand_residual - ligand_residual, final_raw[1])?;
    let total_gold = sum([residuals[0], residuals[3]].into_iter())?;
    let cyanide = sum([
        residuals[1],
        product(2.0, residuals[3].abs())? * residuals[3].signum(),
        actual_ligand_residual,
    ]
    .into_iter())?;
    let components = inventory(layout, &residuals);
    let sodium = components
        .spectators
        .iter()
        .find(|s| s.species.as_str() == dosing::SODIUM)
        .map_or(Amount::new(0.0), |s| s.amount);
    Ok(FlowConservation {
        gold: Amount::new(total_gold),
        cyanide: Amount::new(cyanide),
        dissolved_oxygen: Amount::new(residuals[2]),
        sodium,
        components,
    })
}
fn require_balance(residual: f64, scale: f64) -> Result<()> {
    if residual.abs() > 1024.0 * f64::EPSILON * scale {
        Err(CyanidationError::Arithmetic(
            "flow analytical balance does not close",
        ))
    } else {
        Ok(())
    }
}
fn add(original: f64, increment: f64) -> Result<f64> {
    let value = positive_or_zero(
        original + increment,
        "flow inventory accumulation overflows",
    )?;
    if increment > 0.0 && value <= original {
        Err(CyanidationError::Arithmetic(
            "positive flow increment is below inventory resolution",
        ))
    } else {
        Ok(value)
    }
}

#[cfg(test)]
#[path = "flow/tests.rs"]
mod tests;
