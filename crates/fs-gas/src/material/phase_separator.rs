//! Steady adiabatic routing of an original equilibrium liquid/vapor split.
//!
//! Both outlets retain their actual compositions and molar masses. Flow
//! enthalpy includes pressure work once; routing adds no latent heat, shaft
//! work or second ideal-mixing entropy term. No separator geometry is inferred.

use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError,
        equilibrium::{
            ResolvedThermochemicalPhaseState,
            separation::{PhaseSeparationError, PhaseSeparationSolution},
        },
    },
};
use fs_math::eft::two_sum;
use fs_qty::MassFlowRate;
use fs_thermochem::eos::peng_robinson::phase_screen::PengRobinsonPhaseScreenOptionsV1;
use fs_thermochem::{SemanticError, SpeciesId};
use std::collections::BTreeMap;

/// Original source, temperature, pressure and supplied positive inlet mass flow.
pub use super::phase_heat_exchanger::MixtureHeatExchangerInlet as MixtureSeparatorInlet;

/// Original phase policy and physical total/species/enthalpy/entropy allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureSeparatorControl {
    /// Original bounded phase policy. A reused inlet must retain this same policy.
    pub phase: PengRobinsonPhaseScreenOptionsV1,
    /// Maximum absolute outgoing-minus-incoming enthalpy flow [W].
    pub power_tolerance_w: f64,
    /// Absolute allowance on total and each original species mass flow [kg/s].
    pub absolute_mass_flow_tolerance_kg_s: f64,
    /// Relative mass-flow allowance, scaled separately for every conserved species.
    pub relative_mass_flow_tolerance: f64,
    /// Maximum absolute entropy-flow defect [W/K]. Routing an already admitted
    /// equilibrium split is ideally reversible; neither sign is clipped.
    pub entropy_tolerance_w_k: f64,
}
impl MixtureSeparatorControl {
    /// Use an explicitly supplied original phase screen and physical tolerances.
    pub fn new(phase: PengRobinsonPhaseScreenOptionsV1) -> Self {
        Self {
            phase,
            power_tolerance_w: 1e-5,
            absolute_mass_flow_tolerance_kg_s: 1e-12,
            relative_mass_flow_tolerance: 1e-12,
            entropy_tolerance_w_k: 1e-7,
        }
    }
    /// Admit the original source policy and finite physical tolerances.
    pub fn validate(self) -> Result<(), MixtureSeparatorError> {
        self.phase
            .validate()
            .map_err(|e| MixtureSeparatorError::Equilibrium(ThermochemicalFlashError::Phase(e)))?;
        if !positive(self.power_tolerance_w)
            || !positive(self.absolute_mass_flow_tolerance_kg_s)
            || !self.relative_mass_flow_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_mass_flow_tolerance)
            || !self.entropy_tolerance_w_k.is_finite()
            || self.entropy_tolerance_w_k < 0.0
        {
            return Err(MixtureSeparatorError::Invalid(
                "separator conservation allowances",
            ));
        }
        Ok(())
    }
}

/// Refusal before either physical outlet is published.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureSeparatorError {
    Invalid(&'static str),
    Equilibrium(ThermochemicalFlashError),
    Separation(PhaseSeparationError),
    Composition(SemanticError),
    MassBalance {
        species: Option<SpeciesId>,
        residual_kg_s: f64,
        tolerance_kg_s: f64,
    },
    Energy {
        residual_w: f64,
        tolerance_w: f64,
    },
    Entropy {
        residual_w_k: f64,
        tolerance_w_k: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MixtureSeparatorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adiabatic mixture separator: {self:?}")
    }
}
impl std::error::Error for MixtureSeparatorError {}
type Error = MixtureSeparatorError;

/// One original species before and after physical phase routing.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureSeparatedSpeciesFlow {
    pub species: SpeciesId,
    pub inlet_mass_flow_kg_s: f64,
    /// Actual liquid and vapor species flows [kg/s], in that order.
    pub outlet_mass_flow_kg_s: [f64; 2],
    pub residual_kg_s: f64,
}

/// Both reusable source/caloric outlets and their uncorrected physical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureSeparatorStep {
    pub separation: PhaseSeparationSolution,
    pub inlet_mass_flow: MassFlowRate,
    /// Original mole fractions converted using each phase's actual molar mass.
    pub outlet_mass_flows: [MassFlowRate; 2],
    pub species: Vec<MixtureSeparatedSpeciesFlow>,
    pub mass_residual_kg_s: f64,
    pub energy_residual_w: f64,
    pub entropy_residual_w_k: f64,
    /// New whole-inlet selections: one normally, zero for a reused source state.
    pub inlet_pt_evaluations: usize,
}
impl MixtureSeparatorStep {
    /// Borrow the liquid or vapor as an existing mixer/throttle/heat-exchanger
    /// feed. Its source is owned by this report; no caller source is modified.
    pub fn outlet(&self, index: usize) -> Option<MixtureSeparatorInlet<'_>> {
        self.separation
            .outlets
            .get(index)
            .map(|phase| MixtureSeparatorInlet {
                model: phase.model(),
                temperature: phase.reference().temperature(),
                pressure: phase.reference().pressure(),
                mass_flow: self.outlet_mass_flows[index],
            })
    }
}

/// Separate an original steady equilibrium feed at unchanged T,p.
///
/// The source owner determines a nontrivial split; homogeneous and unresolved
/// phase-boundary feeds refuse rather than acquire a guessed vapor quality.
/// Both outlets conserve the supplied total and species mass flows and original
/// flow enthalpy. Their ideal demixing entropy was already in the inlet flash.
/// No external heat, work, pressure loss or kinetic/potential change is modeled.
/// Residence time, entrainment, droplet separation and geometry need their own
/// physical models. The retained finite phase screen is not a global certificate.
pub fn separate_adiabatic_mixture_stream(
    inlet: MixtureSeparatorInlet<'_>,
    controls: MixtureSeparatorControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureSeparatorStep, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if !positive(inlet.mass_flow.value()) {
        return Err(Error::Invalid(
            "positive finite separator mass flow required",
        ));
    }
    let state = inlet
        .model
        .evaluate_phase_equilibrium(
            inlet.temperature,
            inlet.pressure,
            controls.phase,
            &mut cancelled,
        )
        .map_err(|e| match e {
            ThermochemicalFlashError::Cancelled => Error::Cancelled,
            other => Error::Equilibrium(other),
        })?;
    let mut result = separate_resolved_mixture_stream(
        inlet.model,
        state,
        inlet.mass_flow,
        controls,
        &mut cancelled,
    )?;
    result.inlet_pt_evaluations = 1;
    poll(&mut cancelled)?;
    Ok(result)
}

/// Route an already accepted throttle/PH outlet without repeating its flash.
/// Source identity, retained phase policy and all flow balances are still checked.
pub fn separate_resolved_mixture_stream(
    model: &ThermochemicalPengRobinsonModel,
    state: ResolvedThermochemicalPhaseState,
    mass_flow: MassFlowRate,
    controls: MixtureSeparatorControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureSeparatorStep, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let flow = mass_flow.value();
    if !positive(flow) || state.phase().options() != controls.phase {
        return Err(Error::Invalid(
            "positive flow and unchanged inlet phase policy required",
        ));
    }
    let separation = model
        .separate_phase_equilibrium(state, &mut cancelled)
        .map_err(|e| {
            if e == PhaseSeparationError::Cancelled {
                Error::Cancelled
            } else {
                Error::Separation(e)
            }
        })?;
    let flows = separation.mass_fractions.map(|weight| flow * weight);
    if flows.iter().any(|x| !positive(*x)) {
        return Err(Error::Invalid("positive phase flow was lost or overflowed"));
    }
    let mut inventory: BTreeMap<SpeciesId, [f64; 3]> = species_flows(model, flow)?
        .into_iter()
        .map(|(id, value)| (id, [value, 0.0, 0.0]))
        .collect();
    for (side, phase) in separation.outlets.iter().enumerate() {
        for (id, value) in species_flows(phase.model(), flows[side])? {
            inventory
                .get_mut(&id)
                .ok_or(Error::Invalid("separator source invented a species"))?[side + 1] = value;
        }
    }
    let mut species = Vec::with_capacity(inventory.len());
    for (id, [before, liquid, vapor]) in inventory {
        poll(&mut cancelled)?;
        let residual = sum([liquid, vapor, -before]);
        admit_mass(
            Some(id.clone()),
            residual,
            before.max(liquid + vapor),
            controls,
        )?;
        species.push(MixtureSeparatedSpeciesFlow {
            species: id,
            inlet_mass_flow_kg_s: before,
            outlet_mass_flow_kg_s: [liquid, vapor],
            residual_kg_s: residual,
        });
    }
    let mass_residual = sum([flows[0], flows[1], -flow]);
    admit_mass(None, mass_residual, flow, controls)?;
    let energy = flow * separation.enthalpy_residual_j_kg;
    let entropy = flow * separation.entropy_residual_j_kg_k;
    if !energy.is_finite() || energy.abs() > controls.power_tolerance_w {
        return Err(Error::Energy {
            residual_w: energy,
            tolerance_w: controls.power_tolerance_w,
        });
    }
    if !entropy.is_finite() || entropy.abs() > controls.entropy_tolerance_w_k {
        return Err(Error::Entropy {
            residual_w_k: entropy,
            tolerance_w_k: controls.entropy_tolerance_w_k,
        });
    }
    poll(&mut cancelled)?;
    Ok(MixtureSeparatorStep {
        separation,
        inlet_mass_flow: mass_flow,
        outlet_mass_flows: flows.map(MassFlowRate::new),
        species,
        mass_residual_kg_s: mass_residual,
        energy_residual_w: energy,
        entropy_residual_w_k: entropy,
        inlet_pt_evaluations: 0,
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
                    "original positive species flow is unrepresentable",
                ));
            }
            Ok((s.species().clone(), value))
        })
        .collect()
}
fn admit_mass(
    species: Option<SpeciesId>,
    residual: f64,
    scale: f64,
    controls: MixtureSeparatorControl,
) -> Result<(), Error> {
    let tolerance =
        controls.absolute_mass_flow_tolerance_kg_s + controls.relative_mass_flow_tolerance * scale;
    if !positive(tolerance) || !residual.is_finite() || residual.abs() > tolerance {
        return Err(Error::MassBalance {
            species,
            residual_kg_s: residual,
            tolerance_kg_s: tolerance,
        });
    }
    Ok(())
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let (next, error) = two_sum(total, value);
        total = next;
        correction += error;
    }
    total + correction
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
