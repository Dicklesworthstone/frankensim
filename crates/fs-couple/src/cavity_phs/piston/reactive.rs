//! Declared reaction channels drive the original finite gas and inertial wall.
//!
//! The existing mass-action owner supplies `r(T, p_i/p0)` in mol/(m3 s).
//! By default `dξ_r = dt V r_r` uses the accepted initial rate; the optional
//! [`implicit`] paths solve the single or coupled endpoint rates instead. Both
//! update `n_i' = n_i + sum_r ν_ir dξ_r` with the same conservative owner.
//! The original gas owner reconstructs temperature at conserved mass, volume
//! and TOTAL internal energy, including species formation references. Chemical
//! energy is therefore not added again as a heat input. The new composition,
//! EOS, caloric, transport and acoustic state then drives the original piston
//! pressure/work solve. One physical interval is used by both stages.
//! The [`heat`] and [`boundary`] entries instead expose actual heated or
//! moving-volume endpoints to enclosing simultaneous thermal/mechanical solves.
//!
//! This is a first-order reaction-then-motion split for homogeneous,
//! caller-declared stoichiometric mass-action laws. Each rate scale is declared
//! data, optionally resolved from a material card at the actual gas state;
//! neither an Arrhenius law nor an engine mechanism is inferred. Donor
//! consumption, source-domain, entropy and extensive-balance checks refuse
//! an unresolved step instead of clipping species or renormalizing energy.
//! Refine the interval to assess temporal error. No flame/shock propagation,
//! heterogeneous chemistry or equilibrium sound
//! speed is claimed. The gas owner still enforces its dilute gas/T/p envelope.
//!
//! Control-volume species/energy convention:
//! <https://cantera.org/stable/reference/reactors/ideal-gas-reactor.html>.
//! Rate activities and units are those of `fs_thermochem::equilibrium`, not
//! concentration-based rate constants copied without their basis conversion.

/// Chemical composition coupled to finite material-wall heat and piston motion.
pub mod thermal;
/// Chemical composition coupled to the shared circuit and finite thermal network.
pub mod network;
/// Source-backed kinetic scales evaluated at the actual reacting gas state.
pub mod material_rate;
/// Competing reactions sharing one original inventory and energy reconstruction.
pub mod mechanism;
/// Bounded temporal refinement of the existing conservative chemistry update.
pub mod subcycling;
/// Backward-Euler extent with the same conservative gas and source-rate owners.
pub mod implicit;
/// Endpoint chemistry with explicitly balanced heat on the original gas owner.
pub mod heat;
/// Same-endpoint chemistry with explicit heat, volume and mechanical work.
pub mod boundary;
/// Endpoint chemistry and native piston pressure work solved on one interval.
pub mod coupled;

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
    rate_material: Option<material_rate::MaterialRateScale>,
    additional_channels: Vec<mechanism::Channel>,
    subcycling: Option<subcycling::ReactionSubcycling>,
    implicit: Option<implicit::ImplicitReactionSolve>,
    implicit_mechanism: Option<implicit::coupled::ImplicitMechanismSolve>,
}

/// Physical step admission. Consumption is measured from actual reacting
/// inventory; neither this bound nor entropy admission is an accuracy estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasReactionSolve {
    /// Original fixed-density internal-energy reconstruction controls.
    pub gas: GasEnergySolve,
    /// Maximum fraction withdrawn across ALL channels per substep, in `(0,1)`.
    /// Explicit/scalar steps use the original inventory. Coupled implicit steps
    /// also count simultaneous production as available stock, so an intermediate
    /// can be produced and consumed within the same endpoint solve.
    pub maximum_consumed_fraction: f64,
    /// Maximum discrepancy of any final `n_i` from its stoichiometric value [mol].
    pub species_tolerance_mol: f64,
    /// Maximum tolerated roundoff decrease in closed reaction entropy [J/K].
    pub entropy_tolerance_j_k: f64,
}
impl GasReactionSolve {
    fn validate(self, duration_s: f64) -> Result<(), GasReactionError> {
        self.gas.validate()?;
        if !positive(duration_s)
            || !positive(self.species_tolerance_mol)
            || !positive(self.entropy_tolerance_j_k)
            || !positive(self.maximum_consumed_fraction)
            || self.maximum_consumed_fraction >= 1.0
        {
            return Err(invalid(
                "positive reaction time/tolerances and a consumed-fraction limit in (0,1) required",
            ));
        }
        Ok(())
    }
}

/// No refusal mutates accepted gas, material inertia, geometry or momentum.
#[derive(Clone, Debug, PartialEq)]
pub enum GasReactionError {
    Invalid(&'static str),
    Gas(GasEnergyError),
    Rates(StoichiometricMassActionErrorV1),
    /// Original material source, dimensional or state-domain refusal.
    RateMaterial(fs_material::state_point::MaterialStatePointError),
    Piston(PistonError),
    /// Shorten the physical interval; the consumed species is never clipped.
    ConsumptionLimit {
        species: SpeciesId,
        fraction: f64,
        maximum: f64,
    },
    /// A finite chemical step decreases entropy beyond its rounding allowance.
    EntropyDecrease { change_j_k: f64, tolerance_j_k: f64 },
    /// No complete error-admitted trajectory fits the explicit work limits.
    SubcyclingLimit { evaluations: usize, substeps: usize, last_error_ratio: Option<f64> },
    /// The admitted molar-extent bracket does not contain an implicit root.
    ImplicitOutsideBracket { lower_residual_mol: f64, upper_residual_mol: f64 },
    /// Bounded nonlinear work ended before the amount residual was admitted.
    ImplicitNotConverged { iterations: usize, residual_mol: f64 },
    /// Total coupled endpoint/Jacobian work, including rejected trials, is spent.
    ImplicitEvaluationLimit { evaluations: usize, maximum: usize },
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
    /// Sum of all declared coefficients times channel extents [mol]; zero for inert gas.
    pub stoichiometric_change_mol: f64,
}

/// Original kinetic evaluation and measured closed reaction balances.
#[derive(Clone, Debug, PartialEq)]
pub struct GasReactionStep {
    /// Identity of the first canonical channel; legacy scalar fields describe it.
    pub primary_reaction: fs_thermochem::ReactionId,
    /// Explicit initial or implicit endpoint rate; with subcycling, the first
    /// accepted substep's rate. The corresponding substep records its method.
    pub rate: DeclaredStoichiometricMassActionEvaluationV1,
    /// Actual rate-scale material evaluation; None on the original fixed-scale path.
    pub rate_material: Option<fs_material::state_point::ResolvedMaterialStatePoint>,
    /// Signed extent [mol] of the PRIMARY channel, not a sum across reactions.
    pub extent_mol: f64,
    /// Remaining canonical channel rates/extents; empty for the original API.
    pub additional_reactions: Vec<mechanism::ReactionChannelStep>,
    pub species: Vec<GasReactionSpecies>,
    /// Maximum gross consumption fraction of an accepted chemical substep.
    pub maximum_consumed_fraction: f64,
    pub maximum_species_residual_mol: f64,
    /// Actual total mass discrepancy (the original parcel preserves mass exactly) [kg].
    pub mass_residual_kg: f64,
    /// Actual internal-energy change [J]. The ordinary adiabatic entry bounds
    /// it about zero; heated and moving-boundary entries balance it against
    /// their separately reported heat supply and delivered mechanical work.
    pub stored_energy_change_j: f64,
    /// Measured gas entropy change [J/K]. For a heated endpoint this includes
    /// heat transport; its separate report subtracts the endpoint heat entropy.
    pub entropy_change_j_k: f64,
    /// Actual accepted substeps and estimate; None on the original one-step path.
    pub subcycling: Option<subcycling::ReactionSubcyclingReport>,
    /// Endpoint equation for one implicit step. A subcycled summary leaves this
    /// None and retains each original equation in its actual substeps instead.
    pub implicit: Option<implicit::ImplicitReactionReport>,
    /// Coupled endpoint equations, distinct from the legacy scalar bracket.
    /// Subcycling retains these only in each actual accepted substep.
    pub implicit_mechanism: Option<implicit::coupled::ImplicitMechanismReport>,
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
        Ok(Self { catalog, kinetics, rate_material: None, additional_channels: Vec::new(), subcycling: None, implicit: None, implicit_mechanism: None })
    }
    pub const fn catalog(&self) -> &ThermochemicalGasModel {
        &self.catalog
    }
    /// Original PRIMARY-channel kinetics. With a material source, the effective scale
    /// and its source evaluation are retained in each GasReactionStep instead.
    pub const fn kinetics(&self) -> &DeclaredStoichiometricMassActionV1 {
        &self.kinetics
    }

    // Validate the live subset against the full catalog. An absent product's
    // sources remain in the immutable catalog, and are restored when it appears.
    pub(crate) fn inventory(&self, gas: &ClosedGasParcel) -> Result<Vec<f64>, GasReactionError> {
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

    /// Evolve every reaction on one original inventory at fixed total energy. Inputs
    /// are immutable and publication belongs to the surrounding coupled solve.
    /// Actual partial pressures come from this parcel, not a caller override.
    /// An explicitly bound subcycling policy refines chemistry here, before
    /// any existing mechanical or thermal consumer advances.
    pub fn parcel_candidate(
        &self,
        gas: &ClosedGasParcel,
        duration_s: f64,
        controls: GasReactionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasReactionCandidate, GasReactionError> {
        if let Some(policy) = self.subcycling {
            return subcycling::candidate(self, gas, duration_s, controls, policy, &mut cancelled);
        }
        self.single_step_candidate(gas, duration_s, controls, cancelled)
    }

    // One numerical/constitutive primitive for direct and refined chemistry.
    fn single_step_candidate(
        &self, gas: &ClosedGasParcel, duration_s: f64, controls: GasReactionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasReactionCandidate, GasReactionError> {
        poll(&mut cancelled)?;
        controls.validate(duration_s)?;
        if let Some(policy) = self.implicit {
            return implicit::candidate(self, gas, duration_s, controls, policy, &mut cancelled);
        }
        if let Some(policy) = &self.implicit_mechanism {
            return implicit::coupled::candidate(self, gas, duration_s, controls, policy, &mut cancelled);
        }
        let initial = self.inventory(gas)?;
        let progress = self.progress(gas, &initial, duration_s, controls, &mut cancelled)?;
        let next = self.reconstruct_inventory(gas, &progress.changes, &progress.expected,
            controls, &mut cancelled)?;
        self.admit_progress(gas, &initial, next, progress, controls, &mut cancelled)
    }

    // Shared caloric/source reconstruction for explicit steps and implicit
    // trials. Only the final admitted candidate is checked for entropy growth;
    // nonlinear bracket samples are not physical accepted trajectories.
    fn reconstruct_inventory(
        &self, gas: &ClosedGasParcel, changes: &[f64], expected: &[f64],
        controls: GasReactionSolve, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ClosedGasParcel, GasReactionError> {
        self.reconstruct_inventory_with_heat(gas, changes, expected, 0.0, controls, cancelled)
            .map(|(gas, _)| gas)
    }

    fn reconstruct_inventory_with_heat(
        &self, gas: &ClosedGasParcel, changes: &[f64], expected: &[f64], heat_j: f64,
        controls: GasReactionSolve, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(ClosedGasParcel, usize), GasReactionError> {
        self.reconstruct_inventory_at_boundary(
            gas, changes, expected,
            boundary::GasReactionBoundary {
                final_volume_m3: gas.volume_m3(), heat_j, work_by_gas_j: 0.0,
            },
            controls, cancelled,
        )
    }

    fn reconstruct_inventory_at_boundary(
        &self, gas: &ClosedGasParcel, changes: &[f64], expected: &[f64],
        boundary: boundary::GasReactionBoundary,
        controls: GasReactionSolve, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(ClosedGasParcel, usize), GasReactionError> {
        boundary.validate()?;
        let mut next = gas.clone();
        let mut iterations = 0;
        if changes.iter().any(|change| *change != 0.0) {
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
                .zip(expected)
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
            iterations = next.replace_model_with_heat_work_volume(
                model, boundary.final_volume_m3, boundary.heat_j, boundary.work_by_gas_j,
                controls.gas, cancelled,
            )?.iterations;
        } else if boundary.final_volume_m3 != gas.volume_m3() || boundary.work_by_gas_j != 0.0 {
            let model = next.model().clone();
            iterations = next.replace_model_with_heat_work_volume(
                model, boundary.final_volume_m3, boundary.heat_j, boundary.work_by_gas_j,
                controls.gas, cancelled,
            )?.iterations;
        } else if boundary.heat_j != 0.0 {
            iterations = next.add_heat(boundary.heat_j, controls.gas, cancelled)?.iterations;
        }
        Ok((next, iterations))
    }

    fn admit_progress(
        &self, gas: &ClosedGasParcel, initial: &[f64], next: ClosedGasParcel,
        progress: mechanism::Progress, controls: GasReactionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasReactionCandidate, GasReactionError> {
        self.admit_progress_with_heat(gas, initial, next, progress, 0.0, controls, cancelled)
    }

    fn admit_progress_with_heat(
        &self, gas: &ClosedGasParcel, initial: &[f64], next: ClosedGasParcel,
        progress: mechanism::Progress, heat_j: f64, controls: GasReactionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasReactionCandidate, GasReactionError> {
        self.admit_progress_at_boundary(
            gas, initial, next, progress,
            boundary::GasReactionBoundary {
                final_volume_m3: gas.volume_m3(), heat_j, work_by_gas_j: 0.0,
            },
            true, controls, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_progress_at_boundary(
        &self, gas: &ClosedGasParcel, initial: &[f64], next: ClosedGasParcel,
        progress: mechanism::Progress, boundary: boundary::GasReactionBoundary,
        admit_entropy: bool, controls: GasReactionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasReactionCandidate, GasReactionError> {
        if next.volume_m3().to_bits() != boundary.final_volume_m3.to_bits()
            || next.mass_kg().to_bits() != gas.mass_kg().to_bits()
        {
            return Err(invalid("reaction endpoint lost its original mass or prescribed volume"));
        }
        let mechanism::Progress { primary, additional, changes, expected, consumed } = progress;
        let components = self.catalog.thermodynamics().components();
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
        let energy_residual = (energy + boundary.work_by_gas_j) - boundary.heat_j;
        // Backward Euler uses the actual endpoint temperature for heat entropy
        // transport. This residual includes its nonnegative time-discretization
        // entropy; it is not a measurement of chemical entropy production alone.
        let heat_entropy = boundary.heat_j / next.state().state().temperature;
        let entropy_residual = entropy - heat_entropy;
        if maximum_residual > controls.species_tolerance_mol
            || !energy.is_finite()
            || !energy_residual.is_finite()
            || energy_residual.abs() > gas.mass_kg() * controls.gas.energy_tolerance_j_kg
            || !entropy.is_finite()
            || !heat_entropy.is_finite()
            || !entropy_residual.is_finite()
        {
            return Err(invalid(
                "reaction species, internal energy/heat/work or entropy did not close; refine the interval",
            ));
        }
        if admit_entropy && entropy_residual < -controls.entropy_tolerance_j_k {
            return Err(GasReactionError::EntropyDecrease {
                change_j_k: entropy_residual, tolerance_j_k: controls.entropy_tolerance_j_k,
            });
        }
        poll(cancelled)?;
        Ok(GasReactionCandidate {
            step: GasReactionStep {
                primary_reaction: primary.reaction_id,
                rate: primary.rate,
                rate_material: primary.rate_material,
                extent_mol: primary.extent_mol,
                additional_reactions: additional,
                species,
                maximum_consumed_fraction: consumed,
                maximum_species_residual_mol: maximum_residual,
                mass_residual_kg: next.mass_kg() - gas.mass_kg(),
                stored_energy_change_j: energy,
                entropy_change_j_k: entropy,
                subcycling: None,
                implicit: None,
                implicit_mechanism: None,
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


