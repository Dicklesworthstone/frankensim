//! A declared reaction drives the original finite gas and inertial wall.
//!
//! The existing mass-action owner supplies `r(T, p_i/p0)` in mol/(m3 s).
//! At one accepted state, `dξ = dt V r` and `n_i' = n_i + ν_i dξ`.
//! The original gas owner reconstructs temperature at conserved mass, volume
//! and TOTAL internal energy, including species formation references. Chemical
//! energy is therefore not added again as a heat input. The new composition,
//! EOS, caloric, transport and acoustic state then drives the original piston
//! pressure/work solve. One physical interval is used by both stages.
//!
//! This is a first-order reaction-then-motion split for one homogeneous,
//! caller-declared stoichiometric mass-action law. Its rate scale is declared
//! data, not an inferred Arrhenius law or a validated engine mechanism. Donor
//! consumption, source-domain, entropy and extensive-balance checks refuse
//! an unresolved step instead of clipping species or renormalizing energy.
//! Refine the interval to assess temporal error. No flame/shock propagation,
//! heterogeneous chemistry, stiff mechanism integration or equilibrium sound
//! speed is claimed. The gas owner still enforces its dilute gas/T/p envelope.
//!
//! Control-volume species/energy convention:
//! <https://cantera.org/stable/reference/reactors/ideal-gas-reactor.html>.
//! Rate activities and units are those of `fs_thermochem::equilibrium`, not
//! concentration-based rate constants copied without their basis conversion.

/// Chemical composition coupled to finite material-wall heat and piston motion.
pub mod thermal;

use super::{GasPiston, PistonError, PistonInput, PistonSolve, PistonStep};
use fs_material::MaterialError;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_qty::{Dimensionless, Temperature};
use fs_thermochem::equilibrium::{
    DeclaredStoichiometricMassActionEvaluationV1, DeclaredStoichiometricMassActionV1,
    IdealGasSpeciesActivityV1, StoichiometricMassActionErrorV1,
};
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

/// Reaction authority and a complete source catalog, including initially
/// absent products. This owns no additional temperature or species inventory.
/// The catalog's fractions identify its source set; the parcel owns composition.
#[derive(Clone, Debug, PartialEq)]
pub struct GasReaction {
    catalog: ThermochemicalGasModel,
    kinetics: DeclaredStoichiometricMassActionV1,
}

/// Physical step admission. Consumption is measured from actual reacting
/// inventory; neither this bound nor entropy admission is an accuracy estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasReactionSolve {
    /// Original fixed-density internal-energy reconstruction controls.
    pub gas: GasEnergySolve,
    /// Maximum fraction of any consumed species in `(0,1)` per explicit step.
    pub maximum_consumed_fraction: f64,
    /// Maximum discrepancy of any final `n_i` from its stoichiometric value [mol].
    pub species_tolerance_mol: f64,
    /// Maximum tolerated roundoff decrease in closed reaction entropy [J/K].
    pub entropy_tolerance_j_k: f64,
}

/// No refusal mutates accepted gas, material inertia, geometry or momentum.
#[derive(Clone, Debug, PartialEq)]
pub enum GasReactionError {
    Invalid(&'static str),
    Gas(GasEnergyError),
    Rates(StoichiometricMassActionErrorV1),
    Piston(PistonError),
    /// Shorten the physical interval; the consumed species is never clipped.
    ConsumptionLimit {
        species: SpeciesId,
        fraction: f64,
        maximum: f64,
    },
    Cancelled,
}
impl core::fmt::Display for GasReactionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reactive gas coupling: {self:?}")
    }
}
impl std::error::Error for GasReactionError {}
impl From<GasEnergyError> for GasReactionError {
    fn from(error: GasEnergyError) -> Self {
        if error == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(error)
        }
    }
}
impl From<MaterialError> for GasReactionError {
    fn from(error: MaterialError) -> Self {
        GasEnergyError::from(error).into()
    }
}
impl From<StoichiometricMassActionErrorV1> for GasReactionError {
    fn from(error: StoichiometricMassActionErrorV1) -> Self {
        Self::Rates(error)
    }
}
impl From<PistonError> for GasReactionError {
    fn from(error: PistonError) -> Self {
        if error == PistonError::Cancelled {
            Self::Cancelled
        } else {
            Self::Piston(error)
        }
    }
}
fn invalid(message: &'static str) -> GasReactionError {
    GasReactionError::Invalid(message)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasReactionError> {
    if cancelled() {
        Err(GasReactionError::Cancelled)
    } else {
        Ok(())
    }
}

/// Actual initial/final inventory, including inert and initially absent species.
#[derive(Clone, Debug, PartialEq)]
pub struct GasReactionSpecies {
    pub species: SpeciesId,
    pub initial_mol: f64,
    pub final_mol: f64,
    /// Exact declared coefficient times the accepted extent [mol]; zero for inert gas.
    pub stoichiometric_change_mol: f64,
}

/// Original kinetic evaluation and measured closed reaction balances.
#[derive(Clone, Debug, PartialEq)]
pub struct GasReactionStep {
    pub rate: DeclaredStoichiometricMassActionEvaluationV1,
    /// Signed forward reaction extent [mol] over this one physical interval.
    pub extent_mol: f64,
    pub species: Vec<GasReactionSpecies>,
    pub maximum_consumed_fraction: f64,
    pub maximum_species_residual_mol: f64,
    /// Actual total mass discrepancy (the original parcel preserves mass exactly) [kg].
    pub mass_residual_kg: f64,
    /// Actual internal-energy discrepancy at fixed geometry, bounded about zero [J].
    pub stored_energy_change_j: f64,
    /// Measured entropy change of this isolated chemical stage [J/K].
    pub entropy_change_j_k: f64,
}

/// A proposed gas state usable by the existing flow, thermal and acoustic paths.
#[derive(Clone, Debug, PartialEq)]
pub struct GasReactionCandidate {
    gas: ClosedGasParcel,
    step: GasReactionStep,
}
impl GasReactionCandidate {
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn step(&self) -> &GasReactionStep {
        &self.step
    }
    pub fn into_parts(self) -> (ClosedGasParcel, GasReactionStep) {
        (self.gas, self.step)
    }
}

impl GasReaction {
    /// Bind the existing exact reaction certificate to identical NASA-9
    /// species models in the gas catalog. Additional catalog species are inert.
    /// Element/charge conservation belongs to the existing certificate. This
    /// additionally checks supplied molar masses balance the selected column;
    /// a conserved formula alone does not authenticate caller-supplied masses.
    pub fn new(
        catalog: ThermochemicalGasModel,
        kinetics: DeclaredStoichiometricMassActionV1,
    ) -> Result<Self, GasReactionError> {
        let components = catalog.thermodynamics().components();
        let mut balance = 0.0;
        let mut scale = 0.0;
        for term in kinetics.reverse_rate_closure().equilibrium().terms() {
            let index = components
                .binary_search_by(|s| s.species().cmp(term.species()))
                .map_err(|_| {
                    invalid("reaction species absent from the original gas source catalog")
                })?;
            if &components[index] != term.model() {
                return Err(invalid(
                    "reaction and gas must share identical species caloric sources and reference conventions",
                ));
            }
            let contribution = term.coefficient() as f64 * term.model().molar_mass().value();
            balance += contribution;
            scale += contribution.abs();
        }
        if !positive(scale) || !balance.is_finite() || balance.abs() > 128.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "declared species molar masses do not conserve reaction mass within roundoff",
            ));
        }
        Ok(Self { catalog, kinetics })
    }
    pub const fn catalog(&self) -> &ThermochemicalGasModel {
        &self.catalog
    }
    pub const fn kinetics(&self) -> &DeclaredStoichiometricMassActionV1 {
        &self.kinetics
    }

    // Validate the live subset against the full catalog. An absent product's
    // sources remain in the immutable catalog, and are restored when it appears.
    fn inventory(&self, gas: &ClosedGasParcel) -> Result<Vec<f64>, GasReactionError> {
        let components = self.catalog.thermodynamics().components();
        let model = gas.model();
        let active = model.thermodynamics().components();
        let total = gas.mass_kg()
            / gas
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value();
        if !positive(total) {
            return Err(invalid("positive represented total gas amount required"));
        }
        let mut amounts = vec![0.0; components.len()];
        for (i, species) in active.iter().enumerate() {
            let index = components
                .binary_search_by(|s| s.species().cmp(species.species()))
                .map_err(|_| {
                    invalid("live gas has a species outside the reaction source catalog")
                })?;
            if &components[index] != species
                || self.catalog.transport_materials()[index] != model.transport_materials()[i]
                || self
                    .catalog
                    .source_packs()
                    .iter()
                    .find(|p| &p.species == species.species())
                    != model
                        .source_packs()
                        .iter()
                        .find(|p| &p.species == species.species())
            {
                return Err(invalid(
                    "live gas and reaction catalog disagree on original species sources",
                ));
            }
            let amount = total * model.thermodynamics().mole_composition().fractions()[i];
            if !positive(amount) {
                return Err(invalid("positive gas species inventory is unrepresentable"));
            }
            amounts[index] = amount;
        }
        let expected_water = self
            .catalog
            .water_species()
            .filter(|id| active.iter().any(|s| s.species() == *id));
        if model.water_species() != expected_water {
            return Err(invalid(
                "live gas and reaction catalog disagree on explicit water identity",
            ));
        }
        Ok(amounts)
    }

    /// Evolve one reaction at fixed volume and total internal energy. Inputs
    /// are immutable and publication belongs to the surrounding coupled solve.
    /// Actual partial pressures come from this parcel, not a caller override.
    pub fn parcel_candidate(
        &self,
        gas: &ClosedGasParcel,
        duration_s: f64,
        controls: GasReactionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasReactionCandidate, GasReactionError> {
        poll(&mut cancelled)?;
        controls.gas.validate()?;
        if !positive(duration_s)
            || !positive(controls.species_tolerance_mol)
            || !positive(controls.entropy_tolerance_j_k)
            || !positive(controls.maximum_consumed_fraction)
            || controls.maximum_consumed_fraction >= 1.0
        {
            return Err(invalid(
                "positive reaction time/tolerances and a consumed-fraction limit in (0,1) required",
            ));
        }
        let initial = self.inventory(gas)?;
        let components = self.catalog.thermodynamics().components();
        let terms = self.kinetics.reverse_rate_closure().equilibrium().terms();
        let total: f64 = initial.iter().sum();
        let pressure_ratio = gas.state().state().pressure
            / self
                .kinetics
                .reverse_rate_closure()
                .equilibrium()
                .convention()
                .reference_pressure()
                .value();
        let mut activities = Vec::with_capacity(terms.len());
        for term in terms {
            let index = components
                .binary_search_by(|s| s.species().cmp(term.species()))
                .map_err(|_| invalid("bound reaction species disappeared"))?;
            let value = (initial[index] / total) * pressure_ratio;
            if initial[index] > 0.0 && !positive(value) {
                return Err(invalid(
                    "positive species partial pressure is unrepresentable",
                ));
            }
            activities.push(IdealGasSpeciesActivityV1::new(
                term.species().clone(),
                Dimensionless::new(value),
            )?);
        }
        poll(&mut cancelled)?;
        let rate = self.kinetics.evaluate(
            Temperature::new(gas.state().state().temperature),
            activities,
        )?;
        let progress = rate
            .net_progress_rate()
            .ok_or_else(|| {
                invalid("reaction rate is log-only; no finite extent can be integrated")
            })?
            .value();
        let extent = (progress * gas.volume_m3()) * duration_s;
        if !extent.is_finite() || (progress != 0.0 && extent == 0.0) {
            return Err(invalid("reaction extent overflowed or lost a nonzero rate"));
        }
        let mut changes = vec![0.0; components.len()];
        let mut expected = initial.clone();
        let mut consumed: f64 = 0.0;
        for term in terms {
            let index = components
                .binary_search_by(|s| s.species().cmp(term.species()))
                .map_err(|_| invalid("bound reaction species disappeared"))?;
            let change = term.coefficient() as f64 * extent;
            if !change.is_finite() {
                return Err(invalid(
                    "stoichiometric species increment is unrepresentable",
                ));
            }
            if change < 0.0 {
                let fraction = -change / initial[index];
                if !fraction.is_finite() || fraction > controls.maximum_consumed_fraction {
                    return Err(GasReactionError::ConsumptionLimit {
                        species: term.species().clone(),
                        fraction,
                        maximum: controls.maximum_consumed_fraction,
                    });
                }
                consumed = consumed.max(fraction);
            }
            changes[index] = change;
            expected[index] += change;
            if !expected[index].is_finite()
                || expected[index] < 0.0
                || (change != 0.0 && expected[index] == initial[index])
            {
                return Err(invalid(
                    "reaction lost a species increment or made an invalid inventory",
                ));
            }
        }
        let mut next = gas.clone();
        if extent != 0.0 {
            let total: f64 = expected.iter().sum();
            if !positive(total) {
                return Err(invalid("reaction produced no represented total amount"));
            }
            let mut fractions: Vec<_> = expected.iter().map(|n| n / total).collect();
            let sum: f64 = fractions.iter().sum();
            if (sum - 1.0).abs() > 256.0 * f64::EPSILON {
                return Err(invalid(
                    "reaction mole fractions exceed normalization roundoff",
                ));
            }
            if sum != 1.0 {
                let last = expected
                    .iter()
                    .rposition(|n| *n > 0.0)
                    .ok_or_else(|| invalid("reaction has no active species"))?;
                fractions[last] = 1.0 - fractions[..last].iter().sum::<f64>();
            }
            if fractions
                .iter()
                .zip(&expected)
                .any(|(x, n)| *n > 0.0 && !positive(*x))
            {
                return Err(invalid(
                    "reaction would silently discard a positive trace species",
                ));
            }
            let composition =
                Composition::new(CompositionBasis::MoleFraction, fractions).map_err(|_| {
                    invalid("reaction composition does not satisfy the original quantity owner")
                })?;
            let model = self.catalog.with_composition(composition)?;
            next.replace_model_conserving_mass_energy(model, controls.gas, &mut cancelled)?;
        }
        let actual = self.inventory(&next)?;
        let mut maximum_residual: f64 = 0.0;
        let species = components
            .iter()
            .enumerate()
            .map(|(i, s)| {
                maximum_residual = maximum_residual.max((actual[i] - expected[i]).abs());
                GasReactionSpecies {
                    species: s.species().clone(),
                    initial_mol: initial[i],
                    final_mol: actual[i],
                    stoichiometric_change_mol: changes[i],
                }
            })
            .collect();
        let energy = gas.mass_kg()
            * (next.state().specific_internal_energy_j_kg()
                - gas.state().specific_internal_energy_j_kg());
        let entropy = gas.mass_kg()
            * (next.state().specific_entropy_j_kg_k()? - gas.state().specific_entropy_j_kg_k()?);
        if maximum_residual > controls.species_tolerance_mol
            || !energy.is_finite()
            || energy.abs() > gas.mass_kg() * controls.gas.energy_tolerance_j_kg
            || !entropy.is_finite()
            || entropy < -controls.entropy_tolerance_j_k
        {
            return Err(invalid(
                "closed reaction species, total internal energy or entropy did not close; refine the interval",
            ));
        }
        poll(&mut cancelled)?;
        Ok(GasReactionCandidate {
            step: GasReactionStep {
                rate,
                extent_mol: extent,
                species,
                maximum_consumed_fraction: consumed,
                maximum_species_residual_mol: maximum_residual,
                mass_residual_kg: next.mass_kg() - gas.mass_kg(),
                stored_energy_change_j: energy,
                entropy_change_j_k: entropy,
            },
            gas: next,
        })
    }
}

/// Reaction and original pressure/work reports enclose the same finite gas.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactivePistonStep {
    pub reaction: GasReactionStep,
    /// Original mechanical step starting after the fixed-volume reaction stage.
    pub piston: PistonStep,
    /// Whole-interval internal-energy change, including the chemical stage [J].
    pub gas_energy_change_j: f64,
    /// Actual gas + kinetic change + exterior work - force work - external heat [J].
    pub total_energy_residual_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ReactivePistonCandidate<P = GasPiston> {
    piston: P,
    step: ReactivePistonStep,
}
impl<P> ReactivePistonCandidate<P> {
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    pub const fn step(&self) -> &ReactivePistonStep {
        &self.step
    }
    pub fn into_parts(self) -> (P, ReactivePistonStep) {
        (self.piston, self.step)
    }
}
impl ReactivePistonCandidate {
    pub(super) fn map_piston<P>(
        self,
        bind: impl FnOnce(GasPiston) -> P,
    ) -> ReactivePistonCandidate<P> {
        ReactivePistonCandidate {
            piston: bind(self.piston),
            step: self.step,
        }
    }
}
impl GasPiston {
    // The chemical stage has no separate physical clock or port tick. Retain
    // its entropy and actual energy roundoff in the enclosing accepted account.
    fn record_reactive_account(
        &self,
        next: &mut Self,
        input: PistonInput,
        piston: &PistonStep,
        reaction: &GasReactionStep,
        gas_energy_j: f64,
        tolerance_j: f64,
    ) -> Result<(), GasReactionError> {
        let storage = gas_energy_j + piston.kinetic_energy_change_j;
        let entropy = piston.gas.entropy_change_j_k + reaction.entropy_change_j_k;
        if !storage.is_finite() || !entropy.is_finite() {
            return Err(invalid("reactive gas storage or entropy account is unrepresentable"));
        }
        next.ports = self.ports.clone();
        next.ports.accept(
            crate::port_adoption::StepAccount {
                tick: next.ports.next_tick(),
                duration_s: input.duration_s,
                storage_change_j: Some(storage),
                loss_j: Some(0.0),
                included_work_j: 0.0,
                entropy_change_j_k: entropy,
                tolerance_j: Some(tolerance_j),
                unavailable: None,
            },
            [
                crate::port_adoption::PortTransfer::work(
                    piston.external_force_work_j - piston.exterior_pressure_work_j,
                    true,
                ),
                crate::port_adoption::PortTransfer::heat(
                    input.heat_j,
                    piston.midpoint.state().temperature,
                ),
            ],
        );
        Ok(())
    }

    /// React once, then solve the original simultaneous gas pressure/work and
    /// wall momentum equations. Every mechanical trial sees the SAME reacted
    /// gas; Newton iterations never burn another fraction of the inventory.
    /// `input.heat_j` is external heat only, never a reaction enthalpy credit.
    pub fn reactive_candidate(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        kinetics: GasReactionSolve,
        controls: PistonSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactivePistonCandidate, GasReactionError> {
        let (gas, reaction) = reaction
            .parcel_candidate(&self.gas, input.duration_s, kinetics, &mut cancelled)?
            .into_parts();
        let mut next = self.clone();
        next.gas = gas;
        let piston = next.step(input, controls, &mut cancelled)?;
        let gas_energy = self.gas.mass_kg()
            * (next.gas.state().specific_internal_energy_j_kg()
                - self.gas.state().specific_internal_energy_j_kg());
        let total = gas_energy + piston.kinetic_energy_change_j + piston.exterior_pressure_work_j
            - piston.external_force_work_j
            - input.heat_j;
        if !total.is_finite() || total.abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "reactive gas and original mechanical work do not close whole-interval energy",
            ));
        }
        self.record_reactive_account(&mut next, input, &piston, &reaction,
            gas_energy, controls.energy_tolerance_j)?;
        poll(&mut cancelled)?;
        Ok(ReactivePistonCandidate {
            piston: next,
            step: ReactivePistonStep {
                reaction,
                piston,
                gas_energy_change_j: gas_energy,
                total_energy_residual_j: total,
            },
        })
    }
}

#[cfg(test)]
mod tests;
