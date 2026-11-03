//! Prescribed-pressure oxygen uptake and stripping during gold cyanidation.
//!
//! At fixed liquid volume the explicitly sourced closure is
//! d(n_O2)/dt = kLa (n_sat - n_O2) - r_Au/4. The original gas-solubility owner
//! supplies n_sat = n_ref K p_O2/p_ref from retained formation-Gibbs sources.
//! The original cyanidation owner supplies actual speciation, r_Au and its exact
//! chemical equation. Gold and total cyanide ligand are closed finite pools;
//! a maintained gas reservoir supplies or removes oxygen.
//!
//! Model structure: US EPA, Aeration Equipment Evaluation, EPA/600/2-88/022,
//! section 4, Basic Theory, equation 16:
//! <https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P100GGR5.TXT>.
//! That clean-water report does not calibrate a slurry, ore, tank or aerator.
//! kLa is an independent source declaration in s^-1, never inferred from
//! equilibrium thermochemistry. No finite headspace, gas-feed limit, gas power,
//! changing volume, heat feedback or finite solvent balance is asserted.
//!
//! Bounded full/two-half RK4 trials independently evolve Au, DO and gross oxygen
//! in/out. Raw differences control all four components: division by 15 would
//! presume smooth error scaling across the gross-flux reversal kink. Retained
//! two-half states remain positive without clipping. Local estimates and their
//! accumulated sums are not certified global bounds or model uncertainties.

use super::*;
use fs_qty::{Frequency, Pressure};
use fs_thermochem::stoichiometry::thermodynamics::gas_liquid::IdealGasLiquidPartition;

/// Version-one completely mixed, prescribed-pressure oxygen-transfer closure.
pub const OXYGEN_TRANSFER_LAW: &str = "prescribed-pressure-oxygen-transfer";
/// Independent nonnegative volumetric transfer coefficient, s^-1.
pub const KLA_PARAMETER: &str = "volumetric-transfer-coefficient";
/// Actual oxygen partial pressure, Pa; distinct from source reference pressure.
pub const OXYGEN_PARTIAL_PRESSURE_AXIS: &str = "oxygen-partial-pressure";
/// Optional actual equilibrium dissolved-O2 concentration, mol/m3.
pub const SATURATION_CONCENTRATION_AXIS: &str = "oxygen-saturation-concentration";

/// Immutable original gas equilibrium and independent transfer-rate source.
#[derive(Clone, Debug)]
pub struct OxygenTransfer {
    solubility: IdealGasLiquidPartition,
    source: ConstitutiveModelCard,
    coefficient: Frequency,
}
impl OxygenTransfer {
    /// Require one explicitly sourced kLa and finite, typed support for actual
    /// pO2, DO concentration, absolute clock, temperature and liquid volume.
    pub fn from_card(
        solubility: IdealGasLiquidPartition,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        source.validate()?;
        require(
            solubility.gas().as_str() == "O2(g)" && solubility.aqueous().as_str() == OXYGEN,
            "oxygen transfer requires the original O2(g)/O2(aq) identities",
        )?;
        require(
            source.law.0 == OXYGEN_TRANSFER_LAW
                && source.law_version == 1
                && source.state_schema_version == 0
                && source.initial_state == InitialStatePolicy::ZeroInternalState
                && source.parameters.len() == 1
                && !source.sources.is_empty()
                && !source.provenance.source.trim().is_empty()
                && !source.provenance.license.trim().is_empty(),
            "require a sourced version-one transfer law with one coefficient",
        )?;
        let parameter = source
            .parameters
            .get(KLA_PARAMETER)
            .ok_or(CyanidationError::Invalid("missing transfer coefficient"))?;
        require(
            parameter.dims == Frequency::DIMS,
            "kLa must have inverse-time units",
        )?;
        nonnegative(parameter.value, "kLa must be finite and nonnegative")?;
        let coefficient = Frequency::new(parameter.value);
        for (axis, &(low, high)) in source.validity.bounds() {
            require(
                low.is_finite()
                    && high.is_finite()
                    && low <= high
                    && source.validity.axis_quantities().contains_key(axis),
                "all transfer-source axes need finite ordered typed bounds",
            )?;
        }
        for (index, (axis, dims)) in transfer_axes().into_iter().enumerate() {
            let Some((low, high)) = source.validity.bound(axis) else {
                require(index >= 5, "missing required oxygen-transfer source axis")?;
                continue;
            };
            require(
                source.validity.axis_quantities().get(axis)
                    == Some(&QuantitySpec::dimensional(dims))
                    && low >= 0.0
                    && (!matches!(
                        axis,
                        TEMPERATURE_AXIS | VOLUME_AXIS | ACTIVITY_REFERENCE_AXIS
                    ) || low > 0.0)
                    && (axis != ELAPSED_TIME_AXIS || (low == 0.0 && high > 0.0)),
                "invalid owned oxygen-transfer source axis",
            )?;
        }
        let (low, high) = source
            .validity
            .bound(TEMPERATURE_AXIS)
            .expect("required transfer temperature axis");
        require(
            (low..=high).contains(&solubility.standard().temperature().value()),
            "transfer source must cover the gas equilibrium temperature",
        )?;
        Ok(Self {
            solubility,
            source,
            coefficient,
        })
    }
    /// Original gas equilibrium, exact equation and all source-selection receipts.
    pub fn solubility(&self) -> &IdealGasLiquidPartition {
        &self.solubility
    }
    /// Original independent transfer coefficient and its applicability.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Explicit kLa; zero means no oxygen crosses the boundary.
    pub const fn coefficient(&self) -> Frequency {
        self.coefficient
    }
}

/// Componentwise numerical estimates, all in mol; no physical uncertainty claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AerationErrors {
    pub gold: Amount,
    pub dissolved_oxygen: Amount,
    pub oxygen_in: Amount,
    pub oxygen_out: Amount,
}
impl AerationErrors {
    fn from_values(v: [f64; 4]) -> Self {
        Self {
            gold: Amount::new(v[0]),
            dissolved_oxygen: Amount::new(v[1]),
            oxygen_in: Amount::new(v[2]),
            oxygen_out: Amount::new(v[3]),
        }
    }
}

/// Explicit local-estimate allowances and bounded work for one atomic call.
#[derive(Clone, Copy, Debug)]
pub struct AerationNumerics {
    pub maximum_step: Time,
    pub absolute_gold_error: Amount,
    pub absolute_oxygen_error: Amount,
    /// Applies independently to gross inward and gross outward O2.
    pub absolute_transfer_error: Amount,
    /// Relative to original Au/CN ceiling, current DO/saturation, and transfer
    /// stock scale, respectively; never relative to a large historical total.
    pub relative_error: Dimensionless,
    pub maximum_attempts: usize,
    pub maximum_work: usize,
}
impl Default for AerationNumerics {
    fn default() -> Self {
        Self {
            maximum_step: Time::new(5.0),
            absolute_gold_error: Amount::new(1e-12),
            absolute_oxygen_error: Amount::new(1e-12),
            absolute_transfer_error: Amount::new(1e-12),
            relative_error: Dimensionless::new(1e-7),
            maximum_attempts: 8192,
            maximum_work: 2_000_000,
        }
    }
}
impl AerationNumerics {
    fn tolerances(self, limit: f64, oxygen_scale: f64) -> Result<[f64; 4]> {
        positive(
            self.maximum_step.value(),
            "maximum aeration step must be positive",
        )?;
        require(
            self.relative_error.value().is_finite()
                && (0.0..=0.01).contains(&self.relative_error.value())
                && (1..=65_536).contains(&self.maximum_attempts)
                && (1..=10_000_000).contains(&self.maximum_work),
            "invalid aeration relative accuracy or work limit",
        )?;
        let transfer_scale = sum([oxygen_scale, quotient(limit, 4.0)?].into_iter())?;
        let mut result = [0.0; 4];
        for (i, (absolute, scale)) in [
            (self.absolute_gold_error, limit),
            (self.absolute_oxygen_error, oxygen_scale),
            (self.absolute_transfer_error, transfer_scale),
            (self.absolute_transfer_error, transfer_scale),
        ]
        .into_iter()
        .enumerate()
        {
            positive(
                absolute.value(),
                "absolute aeration error allowances must be positive",
            )?;
            result[i] = sum([
                absolute.value(),
                product(self.relative_error.value(), scale)?,
            ]
            .into_iter())?;
        }
        Ok(result)
    }
}

/// Independently recomputed Au, cyanide-ligand and open-boundary O2 residuals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AerationConservation {
    pub gold: Amount,
    pub cyanide: Amount,
    /// Fixed canonical Na+ spectator minus original and delivered sodium.
    pub sodium: Amount,
    /// DO + original-equation O2 consumption + gross out - initial DO - gross in.
    pub dissolved_oxygen: Amount,
}

/// Persistent chemistry, actual DO, cumulative gas exchange and absolute clock.
#[derive(Clone, Debug)]
pub struct AeratedGoldCyanidationBatch {
    layout: Arc<Layout>,
    transfer: Arc<OxygenTransfer>,
    state: CyanidationState,
    oxygen_in: f64,
    oxygen_out: f64,
    accumulated_errors: [f64; 4],
    balance_allowance: f64,
}

/// One accepted, atomic interval. Original chemistry and transfer sources are retained.
#[derive(Clone, Debug)]
pub struct AerationAdvance {
    layout: Arc<Layout>,
    transfer: Arc<OxygenTransfer>,
    pub previous_elapsed: Time,
    pub elapsed: Time,
    pub interval_gold_dissolved: Amount,
    pub state: CyanidationState,
    pub reaction_extent: Amount,
    pub stoichiometric_amounts: Vec<ReactionAmount>,
    /// Prescribed partial pressure used throughout this call.
    pub oxygen_partial_pressure: Pressure,
    /// Fixed source-derived equilibrium target for this call.
    pub oxygen_saturation: Amount,
    pub oxygen_in: Amount,
    pub oxygen_out: Amount,
    /// Signed cumulative oxygen crossing the gas boundary, inward positive.
    pub net_oxygen_transfer: Amount,
    pub interval_oxygen_in: Amount,
    pub interval_oxygen_out: Amount,
    pub conservation: AerationConservation,
    /// Accumulated floating-point roundoff allowance for the independent O2
    /// balance, separate from integration error estimates; no stock adjustment.
    pub oxygen_balance_allowance: Amount,
    pub maximum_local_errors: AerationErrors,
    pub accumulated_error_estimates: AerationErrors,
    /// Actual admitted endpoint coordinates for the independent kLa source.
    pub transfer_conditions: QueryPoint,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub work: usize,
}
impl AerationAdvance {
    pub fn model(&self) -> &GoldCyanidation {
        &self.layout.model
    }
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
    /// Cumulative accepted dissolved-reagent inputs at this endpoint.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.reagents
    }
}

impl AeratedGoldCyanidationBatch {
    /// Reuse original charge/speciation admission, then retain an open oxygen
    /// boundary. Gold's ceiling is min(original accessible Au, original CN/2);
    /// original DO no longer imposes an exhaustion ceiling.
    pub fn new(
        model: GoldCyanidation,
        transfer: OxygenTransfer,
        charge: CyanidationCharge,
        conditions: QueryPoint,
        speciation: SpeciationOptions,
    ) -> Result<Self> {
        Self::from_batch(
            GoldCyanidationBatch::new(model, charge, conditions, speciation)?,
            transfer,
        )
    }

    /// Attach a maintained-pressure reservoir to a freshly admitted bath,
    /// retaining any simultaneous finite metal/ligand competition. The original
    /// charge and source owners remain unchanged; gas exchange starts at zero.
    pub fn from_batch(closed: GoldCyanidationBatch, transfer: OxygenTransfer) -> Result<Self> {
        require(
            closed.elapsed().value() == 0.0 && closed.state.gold_dissolved.value() == 0.0,
            "maintained-pressure attachment requires a freshly admitted bath",
        )?;
        let model = &closed.layout.model;
        let charge = &closed.layout.charge;
        let conditions = &closed.layout.conditions;
        let external_gas = transfer.solubility.gas();
        require(
            model.chemistry.families().iter().all(|family| {
                family
                    .species()
                    .iter()
                    .all(|species| species != external_gas)
            }) && charge
                .spectators
                .iter()
                .all(|ion| &ion.species != external_gas)
                && closed
                    .layout
                    .competition
                    .as_ref()
                    .is_none_or(|competition| {
                        competition.model.families().iter().all(|family| {
                            family
                                .species()
                                .iter()
                                .all(|species| species != external_gas)
                        })
                    }),
            "the prescribed external gas cannot also belong to a liquid acid family or spectator pool",
        )?;
        let gas = transfer.solubility.standard();
        let water = model.chemistry.water().standard();
        require(
            [OXYGEN_PARTIAL_PRESSURE_AXIS, SATURATION_CONCENTRATION_AXIS]
                .iter()
                .all(|axis| {
                    !conditions.axes().contains_key(*axis)
                        && model.source.validity.bound(axis).is_none()
                }),
            "gas pressure and saturation are owned by the transfer source, not fixed chemical conditions",
        )?;
        require(
            gas.temperature() == water.temperature()
                && gas.reference_pressure() == water.reference_pressure()
                && gas.reference_axis() == water.reference_axis(),
            "gas and acid/base sources need the same temperature, reference pressure and elemental convention",
        )?;
        let aqueous_network = closed
            .layout
            .competition
            .as_ref()
            .map_or(model.chemistry.reaction_network(), |competition| {
                competition.model.reaction_network()
            });
        ReactionNetwork::new(
            aqueous_network
                .reactions()
                .iter()
                .cloned()
                .chain([
                    model.reaction.clone(),
                    transfer.solubility.reaction().clone(),
                ])
                .collect(),
        )?;
        let mut layout = (*closed.layout).clone();
        layout.limit = layout
            .charge
            .accessible_gold
            .value()
            .min(product(0.5, dosing::cyanide_total(&layout)?)?);
        Ok(Self {
            layout: Arc::new(layout),
            transfer: Arc::new(transfer),
            state: closed.state,
            oxygen_in: 0.0,
            oxygen_out: 0.0,
            accumulated_errors: [0.0; 4],
            balance_allowance: 0.0,
        })
    }
    pub fn model(&self) -> &GoldCyanidation {
        &self.layout.model
    }
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.charge
    }
    pub fn state(&self) -> &CyanidationState {
        &self.state
    }
    pub fn elapsed(&self) -> Time {
        self.state.elapsed
    }
    pub fn maximum_gold_dissolution(&self) -> Amount {
        Amount::new(self.layout.limit)
    }
    pub fn oxygen_in(&self) -> Amount {
        Amount::new(self.oxygen_in)
    }
    pub fn oxygen_out(&self) -> Amount {
        Amount::new(self.oxygen_out)
    }
    /// Cumulative reagent inputs, separate from the original initial charge.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.reagents
    }

    /// Re-speciate finite dissolved additions without advancing time or resetting
    /// gross oxygen exchange. The chemical source is checked now; no gas boundary
    /// is retained between calls, so the next `advance` admits its supplied
    /// pressure and transfer context against the new actual solution.
    pub fn dose_reagent(
        &mut self,
        reagent: dosing::DissolvedSodiumReagent,
        amount: Amount,
        maximum_work: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<dosing::ReagentDoseReceipt> {
        let mut work = dosing::work(maximum_work, &mut cancelled)?;
        let layout = Arc::new(dosing::prepare(&self.layout, reagent, amount, &mut work)?);
        let state = evaluate_with_oxygen(
            &layout,
            self.elapsed().value(),
            self.state.gold_dissolved.value(),
            Some(self.state.dissolved_oxygen.value()),
            &mut work,
        )?;
        let amounts = layout
            .model
            .reaction
            .amounts(Amount::new(quotient(state.gold_dissolved.value(), 4.0)?))?;
        open_conservation(
            &layout,
            &state,
            &amounts,
            self.oxygen_in,
            self.oxygen_out,
            self.balance_allowance,
        )?;
        let receipt = dosing::receipt(layout.clone(), reagent, amount, &state, &mut work)?;
        self.layout = layout;
        self.state = state;
        Ok(receipt)
    }

    /// Advance at prescribed fixed pO2. Its value may change on a later call.
    /// Additional kLa source coordinates cannot override actual chemical/gas
    /// coordinates. Every stage is admitted by both original kinetic sources.
    /// Cancellation, unsupported conditions and exhausted budgets publish nothing.
    pub fn advance(
        &mut self,
        duration: Time,
        oxygen_partial_pressure: Pressure,
        transfer_conditions: &QueryPoint,
        numerics: AerationNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AerationAdvance> {
        nonnegative(
            duration.value(),
            "aeration duration must be finite and nonnegative",
        )?;
        let end = sum([self.state.elapsed.value(), duration.value()].into_iter())?;
        require(
            duration.value() == 0.0 || end > self.state.elapsed.value(),
            "positive aeration duration must advance the absolute clock",
        )?;
        let saturation = self.transfer.solubility.equilibrium_aqueous_amount(
            oxygen_partial_pressure,
            self.layout.charge.amount_at_unit_activity,
        )?;
        let tolerances = numerics.tolerances(
            self.layout.limit,
            self.state.dissolved_oxygen.value().max(saturation.value()),
        )?;
        let boundary = Boundary {
            transfer: &self.transfer,
            pressure: oxygen_partial_pressure,
            saturation: saturation.value(),
            extra: transfer_conditions,
        };
        boundary.validate_extras(&self.layout.conditions)?;
        let mut work = Work {
            used: 0,
            maximum: numerics.maximum_work,
            cancelled: &mut cancelled,
            stopped: None,
        };
        work.tick()?;
        let mut current = Sample {
            chemistry: self.state.clone(),
            inward: 0.0,
            outward: 0.0,
        };
        // Zero-duration calls still validate actual current sources and pressure.
        boundary.rate(&current.chemistry)?;
        check_end_clock(&self.layout.model.source, end)?;
        check_end_clock(&self.transfer.source, end)?;
        let mut step = numerics.maximum_step.value().min(duration.value());
        let (mut accepted, mut rejected) = (0, 0);
        let mut largest = [0.0_f64; 4];
        let mut accumulated = self.accumulated_errors;
        let mut balance_allowance = self.balance_allowance;
        while current.chemistry.elapsed.value() < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CyanidationError::Budget);
            }
            let start = current.chemistry.elapsed.value();
            let remaining = end - start;
            step = step.min(remaining);
            let finish = if step == remaining { end } else { start + step };
            let middle = start + 0.5 * (finish - start);
            if middle <= start || finish <= middle {
                return Err(CyanidationError::Arithmetic(
                    "aeration substep is below clock resolution",
                ));
            }
            let trial = (|| -> Result<_> {
                let full = aerated_rk4(&self.layout, &boundary, &current, finish, &mut work)?;
                let first = aerated_rk4(&self.layout, &boundary, &current, middle, &mut work)?;
                let refined = aerated_rk4(&self.layout, &boundary, &first, finish, &mut work)?;
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
            let mut errors = [0.0; 4];
            let mut local = [0.0; 4];
            for (i, (a, b)) in full.values().into_iter().zip(refined.values()).enumerate() {
                errors[i] = (a - b).abs();
                local[i] = product(tolerances[i], quotient(finish - start, duration.value())?)?;
            }
            if errors
                .iter()
                .zip(local)
                .any(|(error, tolerance)| *error > tolerance)
            {
                rejected += 1;
                step = product(step, 0.5)?;
                continue;
            }
            // Linear O2 invariants are preserved by RK weights up to roundoff.
            // Accumulate a per-accepted-step arithmetic allowance; never repair DO.
            let scale = [
                self.layout.charge.dissolved_oxygen.value(),
                current.chemistry.dissolved_oxygen.value(),
                refined.chemistry.dissolved_oxygen.value(),
                quotient(refined.chemistry.gold_dissolved.value(), 4.0)?,
                self.oxygen_in,
                self.oxygen_out,
                refined.inward,
                refined.outward,
            ]
            .into_iter()
            .fold(0.0_f64, f64::max);
            balance_allowance =
                sum([balance_allowance, product(256.0 * f64::EPSILON, scale)?].into_iter())?;
            current = refined;
            accepted += 1;
            for i in 0..4 {
                largest[i] = largest[i].max(errors[i]);
                accumulated[i] = sum([accumulated[i], errors[i]].into_iter())?;
            }
            if errors
                .iter()
                .zip(local)
                .all(|(error, tolerance)| *error < tolerance / 32.0)
            {
                step = (step * 2.0).min(numerics.maximum_step.value());
            }
        }
        let (_, final_point) = boundary.rate(&current.chemistry)?;
        let oxygen_in = add_positive(self.oxygen_in, current.inward)?;
        let oxygen_out = add_positive(self.oxygen_out, current.outward)?;
        let extent = Amount::new(quotient(current.chemistry.gold_dissolved.value(), 4.0)?);
        let amounts = self.layout.model.reaction.amounts(extent)?;
        let conservation = open_conservation(
            &self.layout,
            &current.chemistry,
            &amounts,
            oxygen_in,
            oxygen_out,
            balance_allowance,
        )?;
        let mut answer = AerationAdvance {
            layout: self.layout.clone(),
            transfer: self.transfer.clone(),
            previous_elapsed: self.state.elapsed,
            elapsed: current.chemistry.elapsed,
            interval_gold_dissolved: Amount::new(
                current.chemistry.gold_dissolved.value() - self.state.gold_dissolved.value(),
            ),
            state: current.chemistry.clone(),
            reaction_extent: extent,
            stoichiometric_amounts: amounts,
            oxygen_partial_pressure,
            oxygen_saturation: saturation,
            oxygen_in: Amount::new(oxygen_in),
            oxygen_out: Amount::new(oxygen_out),
            net_oxygen_transfer: Amount::new(sum([oxygen_in, -oxygen_out].into_iter())?),
            interval_oxygen_in: Amount::new(current.inward),
            interval_oxygen_out: Amount::new(current.outward),
            conservation,
            oxygen_balance_allowance: Amount::new(balance_allowance),
            maximum_local_errors: AerationErrors::from_values(largest),
            accumulated_error_estimates: AerationErrors::from_values(accumulated),
            transfer_conditions: final_point,
            accepted_steps: accepted,
            rejected_steps: rejected,
            work: work.used,
        };
        work.tick()?;
        answer.work = work.used;
        self.state = current.chemistry;
        self.oxygen_in = oxygen_in;
        self.oxygen_out = oxygen_out;
        self.accumulated_errors = accumulated;
        self.balance_allowance = balance_allowance;
        Ok(answer)
    }
}

struct Boundary<'a> {
    transfer: &'a OxygenTransfer,
    pressure: Pressure,
    saturation: f64,
    extra: &'a QueryPoint,
}
impl Boundary<'_> {
    fn validate_extras(&self, chemical_extra: &QueryPoint) -> Result<()> {
        for (axis, value) in self.extra.axes() {
            require(
                self.extra.axis_quantities().contains_key(axis)
                    && !owned_axes()
                        .iter()
                        .chain(transfer_axes().iter())
                        .any(|(owned, _)| axis == owned),
                "caller cannot override an owned chemical or transfer coordinate",
            )?;
            if let Some(existing) = chemical_extra.axes().get(axis) {
                require(
                    existing == value
                        && chemical_extra.axis_quantities().get(axis)
                            == self.extra.axis_quantities().get(axis),
                    "shared additional source coordinates must agree",
                )?;
            }
        }
        Ok(())
    }
    fn rate(&self, state: &CyanidationState) -> Result<(f64, QueryPoint)> {
        let mut point = state.conditions.clone();
        for (axis, value) in self.extra.axes() {
            point = point.with_quantity(axis, self.extra.axis_quantities()[axis], *value)?;
        }
        point = point
            .with_quantity(
                OXYGEN_PARTIAL_PRESSURE_AXIS,
                QuantitySpec::dimensional(Pressure::DIMS),
                self.pressure.value(),
            )?
            .with_quantity(
                SATURATION_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                quotient(self.saturation, state.conditions.axes()[VOLUME_AXIS])?,
            )?;
        for (axis, &(minimum, maximum)) in self.transfer.source.validity.bounds() {
            let value = *point.axes().get(axis).ok_or(CyanidationError::Invalid(
                "missing original transfer source condition",
            ))?;
            require(
                point.axis_quantities().get(axis)
                    == self.transfer.source.validity.axis_quantities().get(axis),
                "transfer source condition has wrong quantity",
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
        let difference = self.saturation - state.dissolved_oxygen.value();
        let magnitude = product(self.transfer.coefficient.value(), difference.abs())?;
        Ok((
            if difference < 0.0 {
                -magnitude
            } else {
                magnitude
            },
            point,
        ))
    }
}

#[derive(Clone)]
struct Sample {
    chemistry: CyanidationState,
    // Integration of new transfer starts at zero on each advance, preventing
    // large historical totals from hiding local transfer errors or small rates.
    inward: f64,
    outward: f64,
}
impl Sample {
    fn values(&self) -> [f64; 4] {
        [
            self.chemistry.gold_dissolved.value(),
            self.chemistry.dissolved_oxygen.value(),
            self.inward,
            self.outward,
        ]
    }
}

fn derivative(boundary: &Boundary<'_>, sample: &Sample) -> Result<[f64; 4]> {
    let (flux, _) = boundary.rate(&sample.chemistry)?;
    Ok([
        sample.chemistry.rate.value(),
        sum([flux, -quotient(sample.chemistry.rate.value(), 4.0)?].into_iter())?,
        flux.max(0.0),
        (-flux).max(0.0),
    ])
}
fn displaced(
    layout: &Layout,
    boundary: &Boundary<'_>,
    initial: &Sample,
    elapsed: f64,
    rate: [f64; 4],
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<Sample> {
    let dt = elapsed - initial.chemistry.elapsed.value();
    let x = progress(
        initial.chemistry.gold_dissolved.value(),
        dt,
        rate[0],
        layout.limit,
    )?;
    let change = product(dt, rate[1].abs())?;
    let oxygen = sum([
        initial.chemistry.dissolved_oxygen.value(),
        if rate[1] < 0.0 { -change } else { change },
    ]
    .into_iter())?;
    let old = initial.chemistry.dissolved_oxygen.value();
    if oxygen < 0.0 || (old > 0.0 && oxygen == 0.0) || oxygen > old.max(boundary.saturation) {
        return Err(CyanidationError::TrialDepletion);
    }
    if change > 0.0 && oxygen == old {
        return Err(CyanidationError::Arithmetic(
            "nonzero DO change is below numeric resolution",
        ));
    }
    let chemistry = evaluate_with_oxygen(layout, elapsed, x, Some(oxygen), work)?;
    let sample = Sample {
        chemistry,
        inward: add_positive(initial.inward, product(dt, rate[2])?)?,
        outward: add_positive(initial.outward, product(dt, rate[3])?)?,
    };
    // Admit endpoints even if they do not become a subsequent rate stage.
    boundary.rate(&sample.chemistry)?;
    Ok(sample)
}
fn aerated_rk4(
    layout: &Layout,
    boundary: &Boundary<'_>,
    initial: &Sample,
    end: f64,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<Sample> {
    let middle = stage_midpoint(initial.chemistry.elapsed.value(), end)?;
    let first = derivative(boundary, initial)?;
    let second = displaced(layout, boundary, initial, middle, first, work)?;
    let second = derivative(boundary, &second)?;
    let third = displaced(layout, boundary, initial, middle, second, work)?;
    let third = derivative(boundary, &third)?;
    let fourth = displaced(layout, boundary, initial, end, third, work)?;
    let fourth = derivative(boundary, &fourth)?;
    let mut average = [0.0; 4];
    for i in 0..4 {
        average[i] = sum([
            signed_fraction(first[i], 6.0)?,
            signed_fraction(second[i], 3.0)?,
            signed_fraction(third[i], 3.0)?,
            signed_fraction(fourth[i], 6.0)?,
        ]
        .into_iter())?;
    }
    displaced(layout, boundary, initial, end, average, work)
}
fn stage_midpoint(start: f64, end: f64) -> Result<f64> {
    let middle = start + product(end - start, 0.5)?;
    if middle <= start || middle >= end {
        Err(CyanidationError::Arithmetic(
            "RK midpoint is below absolute clock resolution",
        ))
    } else {
        Ok(middle)
    }
}
fn signed_fraction(value: f64, divisor: f64) -> Result<f64> {
    let magnitude = quotient(value.abs(), divisor)?;
    Ok(if value < 0.0 { -magnitude } else { magnitude })
}
fn add_positive(initial: f64, increment: f64) -> Result<f64> {
    let next = sum([initial, increment].into_iter())?;
    if increment > 0.0 && next == initial {
        Err(CyanidationError::Arithmetic(
            "positive cumulative transfer is below numeric resolution",
        ))
    } else {
        Ok(next)
    }
}
fn check_end_clock(source: &ConstitutiveModelCard, end: f64) -> Result<()> {
    let (minimum, maximum) = source
        .validity
        .bound(ELAPSED_TIME_AXIS)
        .expect("required absolute source clock");
    if end < minimum || end > maximum {
        Err(CyanidationError::OutsideValidity {
            axis: ELAPSED_TIME_AXIS.into(),
            value: end,
            minimum,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn transfer_axes() -> [(&'static str, Dims); 7] {
    [
        (OXYGEN_PARTIAL_PRESSURE_AXIS, Pressure::DIMS),
        (OXYGEN_AXIS, AmountConcentration::DIMS),
        (ELAPSED_TIME_AXIS, Time::DIMS),
        (TEMPERATURE_AXIS, Temperature::DIMS),
        (VOLUME_AXIS, Volume::DIMS),
        (SATURATION_CONCENTRATION_AXIS, AmountConcentration::DIMS),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS),
    ]
}
fn open_conservation(
    layout: &Layout,
    state: &CyanidationState,
    amounts: &[ReactionAmount],
    oxygen_in: f64,
    oxygen_out: f64,
    oxygen_allowance: f64,
) -> Result<AerationConservation> {
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
        oxygen_out,
        -c.dissolved_oxygen.value(),
        -oxygen_in,
    ]
    .into_iter())?;
    for (residual, allowance) in [
        (
            gold,
            128.0 * f64::EPSILON * c.accessible_gold.value().max(c.gold_complex.value()),
        ),
        (
            cyanide - ligand_balance_residual(state),
            128.0
                * f64::EPSILON
                * dosing::cyanide_total(layout)?.max(product(2.0, c.gold_complex.value())?),
        ),
        (oxygen, oxygen_allowance),
    ] {
        if residual.abs() > allowance {
            return Err(CyanidationError::Arithmetic(
                "open-boundary analytical balance does not close",
            ));
        }
    }
    Ok(AerationConservation {
        gold: Amount::new(gold),
        cyanide: Amount::new(cyanide),
        sodium: Amount::new(sodium),
        dissolved_oxygen: Amount::new(oxygen),
    })
}

#[cfg(test)]
mod tests;
