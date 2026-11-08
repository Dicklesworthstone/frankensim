//! Conservative phase calorics for two source-compatible material portions.
//!
//! The original catalog join owns species, source and PR pair compatibility.
//! Conserved species reference enthalpies and entropies cancel before rounding;
//! the original PH owner resolves the mixed phase state without a Cp surrogate.

use super::super::{
    ThermochemicalFlashError, ThermochemicalPengRobinsonModel,
    enthalpy::{EnthalpyWork, MAX_ENTHALPY_FLASH_EVALUATIONS, Target},
    failure,
};
use super::{PhaseEnthalpyControls, PhaseEnthalpySolution, ResolvedThermochemicalPhaseState};
use crate::thermochemical_gas::energy::{
    GasEnergyError, parcel::transfer::mixed_model_from_portions,
};
use crate::thermochemical_gas::nonideal::energy::parcel::transfer::mixing::model_for_source;
use fs_math::{det, eft::two_sum};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Two inlet selections, one mixed reference, and a bounded original PH solve.
pub const MAX_PHASE_MIXING_PT_EVALUATIONS: usize = 3 + MAX_ENTHALPY_FLASH_EVALUATIONS;

/// An original feed whose positive weight is proportional to its supplied mass.
///
/// Only weight ratios affect the caloric closure. A steady-flow consumer can
/// supply its two mass-flow values; no time interval or vessel is inferred here.
#[derive(Clone, Copy, Debug)]
pub struct PhaseMixingFeed<'m> {
    /// Immutable original species calorics, active composition and PR data.
    pub model: &'m ThermochemicalPengRobinsonModel,
    /// Source-supported inlet temperature at the common mixing pressure.
    pub temperature: Temperature,
    /// Positive finite proportional mass weight.
    pub weight: f64,
}

/// Original phase policy, caloric reference and whole-call PT allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseMixingControl {
    /// Common source-supported temperature used to cancel species references.
    /// It must lie in the inclusive outlet temperature bracket. This is a
    /// numerical reference, not a prescribed physical mixing temperature.
    pub reference_temperature: Temperature,
    /// Original outlet PH bracket, tolerances and phase policy. The same phase
    /// policy selects both inlet states and the mixed reference.
    pub equilibrium: PhaseEnthalpyControls,
    /// Whole-call PT selections, including both inlets and the mixed reference.
    pub maximum_pt_evaluations: usize,
}
impl PhaseMixingControl {
    /// Declare the numerical reference and original source controls.
    #[must_use]
    pub fn new(reference_temperature: Temperature, equilibrium: PhaseEnthalpyControls) -> Self {
        Self {
            reference_temperature,
            equilibrium,
            maximum_pt_evaluations: MAX_PHASE_MIXING_PT_EVALUATIONS,
        }
    }
    /// Admit controls before any source work.
    pub fn validate(self) -> Result<(), PhaseMixingError> {
        self.equilibrium
            .validate()
            .map_err(|e| phase_error(2, 0, e))?;
        let t = self.reference_temperature.value();
        let bracket = self
            .equilibrium
            .enthalpy
            .temperature_bracket
            .map(|t| t.value());
        if !positive(t)
            || t < bracket[0]
            || t > bracket[1]
            || !(1..=MAX_PHASE_MIXING_PT_EVALUATIONS).contains(&self.maximum_pt_evaluations)
        {
            return Err(PhaseMixingError::Invalid(
                "phase-mixing reference or work controls",
            ));
        }
        Ok(())
    }
}

/// Original data, phase or finite-work refusal; no partial mixed state is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseMixingError {
    /// Invalid input or a positive amount lost to floating-point representation.
    Invalid(&'static str),
    /// The existing species/NASA/transport/pack/PR catalog join refused.
    Catalog(GasEnergyError),
    /// Original phase/caloric refusal. Stages 0/1 are inlets, 2 is the mixed
    /// reference and 3 is the native PH solve. Nested work is retained.
    Phase {
        stage: usize,
        completed_pt_evaluations: usize,
        source: ThermochemicalFlashError,
    },
    /// Whole-call PT limit reached before another original selection.
    WorkLimit { pt_evaluations: usize },
    /// Cancellation before the complete result was published.
    Cancelled,
}
impl core::fmt::Display for PhaseMixingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original phase mixing: {self:?}")
    }
}
impl std::error::Error for PhaseMixingError {}
type Error = PhaseMixingError;

/// An owned mixed original model, original phase states and datum-free balances.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseMixingSolution {
    /// Reusable original source, now bound to the conserved mixed composition.
    pub model: ThermochemicalPengRobinsonModel,
    /// Both original inlet states, in supplied order.
    pub inlets: [ResolvedThermochemicalPhaseState; 2],
    /// Native mixed PH endpoint. Its work counts only the final caloric solve;
    /// whole-call counts below include the retained inlet/reference selections.
    pub outlet: PhaseEnthalpySolution,
    /// Proportional weights divided by their total.
    pub mass_weights: [f64; 2],
    /// Numerical reference used for the original species increments.
    pub reference_temperature: Temperature,
    /// Required mixed-reference enthalpy increment [J/kg].
    pub target_enthalpy_increment_j_kg: f64,
    /// Native mixed enthalpy minus the conserved incoming enthalpy [J/kg],
    /// evaluated after conserved species formation references cancel.
    pub enthalpy_residual_j_kg: f64,
    /// Ideal composition-mixing contribution [J/(kg K)]. Original phase
    /// redistribution and thermal equilibration remain in the total below.
    pub ideal_mixing_entropy_j_kg_k: f64,
    /// Mixed-reference entropy minus incoming entropy, excluding the ideal
    /// composition contribution [J/(kg K)]; not itself a production bound.
    pub thermal_phase_reference_entropy_j_kg_k: f64,
    /// Complete mixed-minus-incoming original entropy [J/(kg K)], unclipped.
    /// A composing process owner supplies its physical acceptance allowance.
    pub entropy_change_j_kg_k: f64,
    /// Original PT selections across every inlet, reference and outlet query.
    pub pt_evaluations: usize,
    /// Sum of original EOS work ceilings, not an observed EOS count.
    pub eos_evaluation_upper_bound: usize,
}

/// Resolve adiabatic mixing at a supplied common pressure using original laws.
///
/// The species/transport/source-pack join is shared with existing parcel mass
/// transfer. Identical species must have identical original source records;
/// unsupported combinations and conflicting formation/entropy conventions
/// refuse. An optional ORIGINAL PR catalog can provide cross-species critical
/// and interaction entries absent from disjoint feeds; its NASA/transport data
/// are not mixed in. Missing pair interactions are never invented as zero.
///
/// At one reference temperature, conserved ideal species h/s datums cancel
/// analytically. Only source-integrated temperature increments, original phase
/// departures and ideal composition entropy enter the target. The already
/// selected mixed reference feeds the existing PH scalar owner, avoiding a
/// repeated anchor flash. Original phase evidence is retained without claiming
/// global phase stability. This supplies calorics, not flow selection, pressure
/// drop, momentum balance, kinetics or a vessel time step.
pub fn mix_phase_equilibrium_at_pressure(
    feeds: [PhaseMixingFeed<'_>; 2],
    pressure: Pressure,
    interaction_catalog: Option<&ThermochemicalPengRobinsonModel>,
    controls: PhaseMixingControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhaseMixingSolution, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if !positive(pressure.value())
        || feeds
            .iter()
            .any(|f| !positive(f.weight) || !positive(f.temperature.value()))
    {
        return Err(Error::Invalid(
            "positive original mixing pressure, weights and temperatures required",
        ));
    }
    let total = feeds[0].weight + feeds[1].weight;
    let weights = feeds.map(|f| f.weight / total);
    if !positive(total) || weights.iter().any(|w| !positive(*w)) {
        return Err(Error::Invalid("unrepresentable positive mixing inventory"));
    }
    // Identical immutable sources already share one exact canonical inventory.
    // Preserve it directly: a redundant mass/mole round trip can perturb an
    // exact no-mixing state. PR catalogs are still checked below.
    let joined = if feeds[0].model.source() == feeds[1].model.source() {
        feeds[0].model.source().clone()
    } else {
        mixed_model_from_portions(
            &[
                (feeds[0].model.source(), feeds[0].weight),
                (feeds[1].model.source(), feeds[1].weight),
            ],
            total,
        )
        .map_err(Error::Catalog)?
    };
    let mut catalogs = vec![feeds[0].model, feeds[1].model];
    if let Some(catalog) = interaction_catalog {
        catalogs.push(catalog);
    }
    let model = model_for_source(&catalogs, joined).map_err(Error::Catalog)?;
    model
        .admit_enthalpy_flash(
            pressure,
            controls.equilibrium.enthalpy,
            Some(controls.reference_temperature),
            &mut cancelled,
        )
        .map_err(|e| phase_error(2, 0, e))?;
    let mut pt = 0;
    let mut eos = 0;
    let mut select = |source: &ThermochemicalPengRobinsonModel, t: Temperature, stage| {
        poll(&mut cancelled)?;
        if pt == controls.maximum_pt_evaluations {
            return Err(Error::WorkLimit { pt_evaluations: pt });
        }
        let sample = source
            .caloric_phase_sample(
                t,
                pressure,
                controls.equilibrium.phase_options(),
                &mut cancelled,
            )
            .map_err(|e| phase_error(stage, pt, e))?;
        pt += 1;
        eos += sample.eos_evaluation_upper_bound;
        Ok(sample)
    };
    let samples = [
        select(feeds[0].model, feeds[0].temperature, 0)?,
        select(feeds[1].model, feeds[1].temperature, 1)?,
    ];
    let anchor = select(&model, controls.reference_temperature, 2)?;
    let anchor_mass = model.mechanics().molar_mass().value();
    let anchor_h = anchor.residual_h / anchor_mass;
    let anchor_s = anchor.entropy_correction / anchor_mass;
    let mut enthalpy = [0.0; 2];
    let mut entropy = [0.0; 2];
    for i in 0..2 {
        poll(&mut cancelled)?;
        let source = feeds[i].model.source().thermodynamics();
        let dh = source
            .specific_enthalpy_change(controls.reference_temperature, feeds[i].temperature)
            .map_err(|e| {
                phase_error(
                    i,
                    pt,
                    failure("mixing original enthalpy increment", e).into(),
                )
            })?
            .value();
        let ds = source
            .specific_entropy_change(
                controls.reference_temperature,
                pressure,
                feeds[i].temperature,
                pressure,
            )
            .map_err(|e| {
                phase_error(
                    i,
                    pt,
                    failure("mixing original entropy increment", e).into(),
                )
            })?
            .value();
        let mass = feeds[i].model.mechanics().molar_mass().value();
        enthalpy[i] = weights[i] * (dh + (samples[i].residual_h / mass - anchor_h));
        entropy[i] = weights[i] * ((anchor_s - samples[i].entropy_correction / mass) - ds);
    }
    let target = sum(enthalpy);
    let reference_entropy = sum(entropy);
    let mixing_entropy = ideal_composition_entropy(&feeds, weights, &model);
    if [target, reference_entropy, mixing_entropy]
        .iter()
        .any(|x| !x.is_finite())
    {
        return Err(Error::Invalid(
            "unrepresentable reference-free mixing target",
        ));
    }
    let remaining = controls.maximum_pt_evaluations - pt;
    if remaining == 0 && target != 0.0 {
        return Err(Error::WorkLimit { pt_evaluations: pt });
    }
    let mut caloric = controls.equilibrium.enthalpy;
    caloric.maximum_flash_evaluations = caloric.maximum_flash_evaluations.min(remaining.max(1));
    let outcome = model
        .solve_flash_enthalpy(
            pressure,
            Target::Increment {
                initial: anchor,
                requested: target,
            },
            caloric,
            EnthalpyWork::default(),
            &mut |t, p, c| {
                model.caloric_phase_sample(t, p, controls.equilibrium.phase_options(), c)
            },
            &mut cancelled,
        )
        .map_err(|e| phase_error(3, pt, e))?;
    let outlet = model
        .publish_phase_enthalpy(outcome, &mut cancelled)
        .map_err(|e| phase_error(3, pt, e))?;
    pt += outlet.pt_evaluations();
    eos += outlet.eos_evaluation_upper_bound();
    let entropy_change = sum([
        outlet
            .entropy_change_j_kg_k()
            .ok_or(Error::Invalid("missing mixed source entropy increment"))?,
        reference_entropy,
        mixing_entropy,
    ]);
    if !entropy_change.is_finite() {
        return Err(Error::Invalid("nonfinite original mixing entropy"));
    }
    let [first, second] = samples;
    let inlets = [
        feeds[0]
            .model
            .publish_caloric_phase(first, &mut cancelled)
            .map_err(|e| phase_error(0, pt, e))?,
        feeds[1]
            .model
            .publish_caloric_phase(second, &mut cancelled)
            .map_err(|e| phase_error(1, pt, e))?,
    ];
    let residual = outlet.residual_j_kg();
    poll(&mut cancelled)?;
    Ok(PhaseMixingSolution {
        model,
        inlets,
        outlet,
        mass_weights: weights,
        reference_temperature: controls.reference_temperature,
        target_enthalpy_increment_j_kg: target,
        enthalpy_residual_j_kg: residual,
        ideal_mixing_entropy_j_kg_k: mixing_entropy,
        thermal_phase_reference_entropy_j_kg_k: reference_entropy,
        entropy_change_j_kg_k: entropy_change,
        pt_evaluations: pt,
        eos_evaluation_upper_bound: eos,
    })
}

fn ideal_composition_entropy(
    feeds: &[PhaseMixingFeed<'_>; 2],
    weights: [f64; 2],
    mixed: &ThermochemicalPengRobinsonModel,
) -> f64 {
    sum(feeds.iter().enumerate().map(|(i, feed)| {
        let components = feed.model.source().thermodynamics().components();
        let fractions = feed.model.mechanics().mole_fractions();
        let divergence = sum(mixed
            .source()
            .thermodynamics()
            .components()
            .iter()
            .zip(mixed.mechanics().mole_fractions())
            .map(|(species, &z)| {
                let x = components
                    .binary_search_by(|s| s.species().cmp(species.species()))
                    .map(|index| fractions[index])
                    .unwrap_or(0.0);
                relative_entropy_term(x, z)
            }));
        R * (weights[i] / feed.model.mechanics().molar_mass().value()) * divergence
    }))
}

// x ln(x/y) - x + y. Linear terms sum to zero in each canonical composition.
// A bounded small-difference series retains quadratic mixing entropy when the
// compositions are nearly equal; positive trace amounts are never thresholded.
fn relative_entropy_term(x: f64, y: f64) -> f64 {
    if x == 0.0 {
        return y;
    }
    if x == y {
        return 0.0;
    }
    let r = (x - y) / y;
    if r.abs() <= 0.25 {
        let mut power = r * r;
        let mut terms = [0.0; 31];
        for (index, term) in terms.iter_mut().enumerate() {
            let n = index + 2;
            *term = power / ((n * (n - 1)) as f64);
            power *= -r;
        }
        y * sum(terms)
    } else {
        x * (det::ln(x) - det::ln(y)) - x + y
    }
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
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn phase_error(
    stage: usize,
    completed_pt_evaluations: usize,
    source: ThermochemicalFlashError,
) -> Error {
    if source == ThermochemicalFlashError::Cancelled {
        Error::Cancelled
    } else {
        Error::Phase {
            stage,
            completed_pt_evaluations,
            source,
        }
    }
}
