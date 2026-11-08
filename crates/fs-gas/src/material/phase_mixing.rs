//! Adiabatic mixing of two original phase-equilibrium streams at common pressure.
//!
//! The native material mixer retains all species source cards and PR data,
//! cancels conserved formation/entropy references, and resolves phase change.
//! This owner supplies steady flow conservation and physical entropy acceptance.

use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::equilibrium::{
        PhaseEnthalpyControls,
        mixing::{
            PhaseMixingControl, PhaseMixingError, PhaseMixingFeed, PhaseMixingSolution,
            mix_phase_equilibrium_at_pressure,
        },
    },
};
use fs_math::eft::two_sum;
use fs_qty::{MassFlowRate, Temperature};
use fs_thermochem::{SemanticError, SpeciesId};
use std::collections::BTreeMap;

/// A steady original feed: pressure, temperature, composition and mass flow.
pub use super::phase_heat_exchanger::MixtureHeatExchangerInlet as MixtureMixingInlet;

/// Original source policy and datum-invariant conservation allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureMixingControl {
    /// Shared original phase policy, output bracket, numerical reference and
    /// whole-call original PT allowance.
    pub phase: PhaseMixingControl,
    /// Absolute outgoing-minus-incoming enthalpy-flow allowance [W].
    pub absolute_power_tolerance_w: f64,
    /// Relative allowance scaled by native reference-free caloric increments,
    /// never by arbitrary absolute formation enthalpies.
    pub relative_power_tolerance: f64,
    /// Absolute total and individual-species flow allowance [kg/s].
    pub absolute_mass_flow_tolerance_kg_s: f64,
    /// Relative mass-flow allowance for each conserved species and total mass.
    pub relative_mass_flow_tolerance: f64,
    /// Numerical allowance for negative total original entropy production [W/K].
    pub entropy_tolerance_w_k: f64,
}
impl MixtureMixingControl {
    /// Declare original phase calorics and a source-supported numerical reference.
    #[must_use]
    pub fn new(reference_temperature: Temperature, equilibrium: PhaseEnthalpyControls) -> Self {
        Self {
            phase: PhaseMixingControl::new(reference_temperature, equilibrium),
            absolute_power_tolerance_w: 1e-5,
            relative_power_tolerance: 1e-9,
            absolute_mass_flow_tolerance_kg_s: 1e-12,
            relative_mass_flow_tolerance: 1e-12,
            entropy_tolerance_w_k: 1e-7,
        }
    }
    /// Validate finite physical/source allowances before selecting any phase.
    pub fn validate(self) -> Result<(), MixtureMixingError> {
        self.phase
            .validate()
            .map_err(MixtureMixingError::Material)?;
        if !positive(self.absolute_power_tolerance_w)
            || !positive(self.absolute_mass_flow_tolerance_kg_s)
            || !valid_relative(self.relative_power_tolerance)
            || !valid_relative(self.relative_mass_flow_tolerance)
            || !self.entropy_tolerance_w_k.is_finite()
            || self.entropy_tolerance_w_k < 0.0
        {
            return Err(MixtureMixingError::Invalid(
                "steady mixture-mixing tolerances",
            ));
        }
        Ok(())
    }
}

/// No partial outlet or replacement original source is published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureMixingError {
    /// Invalid input, pressure mismatch or unrepresentable positive inventory.
    Invalid(&'static str),
    /// Native source join, original phase, caloric or bounded-work refusal.
    Material(PhaseMixingError),
    /// Original quantity owner's mass/mole composition conversion refused.
    Composition(SemanticError),
    /// Actual outgoing mass or a conserved species failed its flow allowance.
    MassBalance {
        species: Option<SpeciesId>,
        residual_kg_s: f64,
        tolerance_kg_s: f64,
    },
    /// Original outgoing-minus-incoming enthalpy flow failed closure.
    Energy { residual_w: f64, tolerance_w: f64 },
    /// Original outgoing-minus-incoming entropy violated its allowance.
    Entropy {
        production_w_k: f64,
        tolerance_w_k: f64,
    },
    /// Cancellation before the full original outlet was published.
    Cancelled,
}
impl core::fmt::Display for MixtureMixingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adiabatic original mixture streams: {self:?}")
    }
}
impl std::error::Error for MixtureMixingError {}
type Error = MixtureMixingError;

/// Canonical original species ledger across the two incoming streams.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureSpeciesFlow {
    /// Original conserved species identity.
    pub species: SpeciesId,
    /// Actual inlet mass flows in feed order [kg/s].
    pub inlet_mass_flow_kg_s: [f64; 2],
    /// Actual mass flow from the accepted mixed composition [kg/s].
    pub outlet_mass_flow_kg_s: f64,
    /// Outgoing minus both incoming species flows [kg/s].
    pub residual_kg_s: f64,
}

/// A reusable mixed original model, phase outlet and physical flow balances.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureMixingStep {
    /// Owned original mixture, inlet states, accepted outlet and native work.
    /// This model can directly feed downstream PH, machinery and nozzle owners.
    pub mixture: PhaseMixingSolution,
    /// Sum of both positive supplied inlet mass flows.
    pub mass_flow: MassFlowRate,
    /// Canonical per-species incoming/outgoing flow ledger.
    pub species: Vec<MixtureSpeciesFlow>,
    /// Sum of accepted species mass flows minus supplied total flow [kg/s].
    pub mass_residual_kg_s: f64,
    /// Datum-free outgoing-minus-incoming enthalpy flow [W].
    pub energy_residual_w: f64,
    /// Ideal composition-mixing part of entropy production [W/K].
    pub ideal_mixing_entropy_flow_w_k: f64,
    /// Complete original outgoing-minus-incoming entropy [W/K], unclipped.
    pub entropy_production_w_k: f64,
}

/// Mix two steady streams adiabatically at their common fixed pressure.
///
/// Supplied flows, species and total original enthalpy are conserved. The
/// original PR departures retain real-gas heat of mixing and the original
/// phase owner may convert either inlet into a different outlet phase state.
/// The output temperature need not be an inlet-temperature average.
///
/// This model has negligible chamber holdup and velocity/potential changes,
/// no shaft work, external heat, reaction or pressure drop. Stream enthalpy
/// contains pressure-flow work exactly once. Distinct inlet pressures require
/// their own upstream valve/pressure-change owner and are rejected here.
/// An optional original PR catalog supplies known cross-species data for
/// disjoint feeds; missing or conflicting source data are native refusals.
/// Original formation and entropy datums cancel before rounding. Species
/// rounding is audited separately with its declared physical flow allowance.
pub fn mix_adiabatic_mixture_streams(
    inlets: [MixtureMixingInlet<'_>; 2],
    interaction_catalog: Option<&ThermochemicalPengRobinsonModel>,
    mut controls: MixtureMixingControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureMixingStep, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if inlets.iter().any(|s| {
        !positive(s.mass_flow.value())
            || !positive(s.pressure.value())
            || !positive(s.temperature.value())
    }) || inlets[0].pressure != inlets[1].pressure
    {
        return Err(Error::Invalid(
            "positive flows and temperatures at exactly common pressure required",
        ));
    }
    let flow = inlets[0].mass_flow.value() + inlets[1].mass_flow.value();
    if !positive(flow) {
        return Err(Error::Invalid("unrepresentable total mixed mass flow"));
    }
    controls.phase.equilibrium.enthalpy.enthalpy_tolerance_j_kg = controls
        .phase
        .equilibrium
        .enthalpy
        .enthalpy_tolerance_j_kg
        .min(0.25 * controls.absolute_power_tolerance_w / flow);
    controls.validate()?;
    let mixture = mix_phase_equilibrium_at_pressure(
        inlets.map(|s| PhaseMixingFeed {
            model: s.model,
            temperature: s.temperature,
            weight: s.mass_flow.value(),
        }),
        inlets[0].pressure,
        interaction_catalog,
        controls.phase,
        &mut cancelled,
    )
    .map_err(|e| {
        if e == PhaseMixingError::Cancelled {
            Error::Cancelled
        } else {
            Error::Material(e)
        }
    })?;
    let mut incoming: BTreeMap<SpeciesId, [f64; 2]> = BTreeMap::new();
    for (side, inlet) in inlets.iter().enumerate() {
        for (species, value) in species_flows(inlet.model, inlet.mass_flow.value())? {
            incoming.entry(species).or_insert([0.0; 2])[side] = value;
        }
    }
    let mut species = Vec::with_capacity(incoming.len());
    let mut outgoing = Vec::with_capacity(incoming.len());
    for (id, value) in species_flows(&mixture.model, flow)? {
        poll(&mut cancelled)?;
        let before = incoming
            .remove(&id)
            .ok_or(Error::Invalid("mixed source invented a species"))?;
        let residual = sum([value, -before[0], -before[1]]);
        let tolerance = controls.absolute_mass_flow_tolerance_kg_s
            + controls.relative_mass_flow_tolerance * value.abs().max(before[0] + before[1]);
        if !positive(tolerance) || !residual.is_finite() || residual.abs() > tolerance {
            return Err(Error::MassBalance {
                species: Some(id),
                residual_kg_s: residual,
                tolerance_kg_s: tolerance,
            });
        }
        outgoing.push(value);
        species.push(MixtureSpeciesFlow {
            species: id,
            inlet_mass_flow_kg_s: before,
            outlet_mass_flow_kg_s: value,
            residual_kg_s: residual,
        });
    }
    if !incoming.is_empty() {
        return Err(Error::Invalid("mixed source lost a conserved species"));
    }
    let mass_residual = sum(outgoing.into_iter().chain([-flow]));
    let mass_tolerance =
        controls.absolute_mass_flow_tolerance_kg_s + controls.relative_mass_flow_tolerance * flow;
    if !positive(mass_tolerance)
        || !mass_residual.is_finite()
        || mass_residual.abs() > mass_tolerance
    {
        return Err(Error::MassBalance {
            species: None,
            residual_kg_s: mass_residual,
            tolerance_kg_s: mass_tolerance,
        });
    }
    let energy = flow * mixture.enthalpy_residual_j_kg;
    let accepted = mixture
        .outlet
        .accepted_increment_j_kg()
        .ok_or(Error::Invalid("missing original mixed enthalpy increment"))?;
    let scale = flow
        * accepted
            .abs()
            .max(mixture.target_enthalpy_increment_j_kg.abs());
    let energy_tolerance =
        controls.absolute_power_tolerance_w + controls.relative_power_tolerance * scale;
    if !positive(energy_tolerance) || !energy.is_finite() || energy.abs() > energy_tolerance {
        return Err(Error::Energy {
            residual_w: energy,
            tolerance_w: energy_tolerance,
        });
    }
    let entropy = flow * mixture.entropy_change_j_kg_k;
    let mixing_entropy = flow * mixture.ideal_mixing_entropy_j_kg_k;
    if !entropy.is_finite()
        || !mixing_entropy.is_finite()
        || entropy < -controls.entropy_tolerance_w_k
    {
        return Err(Error::Entropy {
            production_w_k: entropy,
            tolerance_w_k: controls.entropy_tolerance_w_k,
        });
    }
    poll(&mut cancelled)?;
    Ok(MixtureMixingStep {
        mixture,
        mass_flow: MassFlowRate::new(flow),
        species,
        mass_residual_kg_s: mass_residual,
        energy_residual_w: energy,
        ideal_mixing_entropy_flow_w_k: mixing_entropy,
        entropy_production_w_k: entropy,
    })
}

fn species_flows(
    model: &ThermochemicalPengRobinsonModel,
    flow: f64,
) -> Result<Vec<(SpeciesId, f64)>, Error> {
    let original = model.source().thermodynamics();
    let masses: Vec<_> = original
        .components()
        .iter()
        .map(|s| s.molar_mass())
        .collect();
    let fractions = original
        .mole_composition()
        .to_mass_fractions(&masses)
        .map_err(Error::Composition)?;
    original
        .components()
        .iter()
        .zip(fractions.fractions())
        .map(|(s, fraction)| {
            let value = flow * fraction;
            if !positive(value) {
                return Err(Error::Invalid(
                    "positive original species flow was lost or overflowed",
                ));
            }
            Ok((s.species().clone(), value))
        })
        .collect()
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut total = 0.0;
    let mut correction = 0.0;
    for value in values {
        let (next, error) = two_sum(total, value);
        total = next;
        correction += error;
    }
    total + correction
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn valid_relative(x: f64) -> bool {
    x.is_finite() && (0.0..1.0).contains(&x)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
