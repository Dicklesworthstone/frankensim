//! Evaluated prompt spontaneous-fission neutrons on the original finite stock.
//!
//! MF1/MT456 gives neutrons per spontaneous fission, while the constant MF5/MT18
//! component weights and outgoing densities describe their energies. Native
//! decay events include all supplied parent branches, so the effective source
//! uses B_SF/B_all without changing either original branch sum. Ordinary MF8
//! neutron spectra are already per parent transformation and receive no second
//! branching factor. Explicit radioactive products provide delayed emission;
//! aggregate delayed-family yields are not added to those products again.
//!
//! All moments use the original analytical emission integrator and conservative
//! representative-energy projection. Neutron energy must be funded by the
//! parent's explicitly routed heavy-particle radiation. Source probabilities,
//! densities and multiplicities are retained without normalization repair.

use std::{collections::BTreeSet, sync::Arc};

use super::{BoundEvaluatedBurnupDecay, DecayAssemblyError, EvaluatedBurnupDecayAssembly};
use crate::burnup::BurnupModel;
use crate::diffusion::burnup::DiffusionBurnupError;
use crate::diffusion::burnup::decay::BurnupDecayNeutronEmission;
use crate::endf::EndfError;
use crate::endf::decay::DecayEvaluation;
use crate::endf::decay::model::emissions::neutrons::{
    DecayNeutronBindingError, DecayNeutronProjection, project_group_moments,
};
use crate::endf::decay::model::emissions::{
    DecayEmissionCoverage, DecayEmissionLimits, integrate_continuum, integrate_spectrum, visit,
};
use crate::endf::decay::spectra::{DecayRadiation, DecaySpectrum};
use crate::endf::neutron_yield::FissionNeutronKind;
use crate::endf::spectrum::PositiveIntegral;
use crate::endf::spontaneous_neutrons::{
    SpontaneousFissionNeutronSpectrum, SpontaneousFissionNeutronYield,
};
use crate::fusion_neutrons::NeutronKineticModel;
use crate::irradiation::NuclideId;
use crate::{Accumulator, FissionError, nonnegative, sum};

/// Matched original prompt multiplicity and outgoing spectrum for one SF parent.
#[derive(Clone, Debug)]
pub struct SpontaneousFissionPromptNeutronBinding {
    pub parent: NuclideId,
    pub multiplicity: Arc<SpontaneousFissionNeutronYield>,
    pub spectrum: Arc<SpontaneousFissionNeutronSpectrum>,
}

/// Preparation bounds and an explicit absolute probability-normalization
/// allowance. Admitted discrepancies remain observable and are never repaired.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaluatedBurnupNeutronBudget {
    pub limits: DecayEmissionLimits,
    pub normalization_tolerance: f64,
}

/// Original branching and retained source normalization, in binding order.
#[derive(Clone, Debug, PartialEq)]
pub struct SpontaneousFissionNeutronNormalization {
    pub spontaneous_branch_sum: f64,
    pub total_decay_branch_sum: f64,
    /// Conversion from all native decay events to native SF events.
    pub fission_events_per_decay: f64,
    pub component_probability_sum: f64,
    pub component_integrals: Arc<[f64]>,
    pub weighted_spectrum_integral: f64,
    pub declared_particles_per_fission: f64,
    pub integrated_particles_per_fission: f64,
    pub multiplicity_residual: f64,
}

/// One immutable evaluated source on an original native inventory/group grid.
/// It prepares coefficients only; evolution and all energy ledgers stay native.
#[derive(Clone, Debug)]
pub struct EvaluatedBurnupNeutrons {
    data_id: Arc<str>,
    assembly: Arc<EvaluatedBurnupDecayAssembly>,
    kinetic_model: Arc<NeutronKineticModel>,
    spontaneous_sources: Arc<[SpontaneousFissionPromptNeutronBinding]>,
    normalization: Arc<[SpontaneousFissionNeutronNormalization]>,
    coverage: Arc<[DecayEmissionCoverage]>,
    confirmed_nonemitting_parents: Arc<[usize]>,
    projection: Arc<[DecayNeutronProjection]>,
    neutron_yields: Arc<[Vec<f64>]>,
}

impl EvaluatedBurnupNeutrons {
    pub fn new(
        data_id: String,
        assembly: Arc<EvaluatedBurnupDecayAssembly>,
        kinetic_model: Arc<NeutronKineticModel>,
        spontaneous_sources: Vec<SpontaneousFissionPromptNeutronBinding>,
        confirmed_nonemitting_parents: &[usize],
        budget: EvaluatedBurnupNeutronBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EvaluatedBurnupNeutronError> {
        prepare(
            data_id,
            assembly,
            kinetic_model,
            spontaneous_sources,
            confirmed_nonemitting_parents,
            budget,
            &mut checkpoint,
        )
    }

    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    pub fn assembly(&self) -> &Arc<EvaluatedBurnupDecayAssembly> {
        &self.assembly
    }
    pub fn kinetic_model(&self) -> &Arc<NeutronKineticModel> {
        &self.kinetic_model
    }
    pub fn spontaneous_sources(&self) -> &[SpontaneousFissionPromptNeutronBinding] {
        &self.spontaneous_sources
    }
    pub fn normalization(&self) -> &[SpontaneousFissionNeutronNormalization] {
        &self.normalization
    }
    pub fn coverage(&self) -> &[DecayEmissionCoverage] {
        &self.coverage
    }
    pub fn confirmed_nonemitting_parents(&self) -> &[usize] {
        &self.confirmed_nonemitting_parents
    }
    pub fn projection(&self) -> &[DecayNeutronProjection] {
        &self.projection
    }
    pub fn neutron_yields(&self) -> &[Vec<f64>] {
        &self.neutron_yields
    }

    /// Retain the exact evaluated decay mapping while admitting reordered native
    /// channels and appended neutron reactions. The resulting emission targets
    /// that same finite inventory; no stock or history is recreated.
    pub fn bind_burnup(
        self: &Arc<Self>,
        data_id: String,
        model: Arc<BurnupModel>,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BoundEvaluatedBurnupNeutrons, EvaluatedBurnupNeutronError> {
        checkpoint()?;
        let binding = self
            .assembly
            .bind_burnup(Arc::clone(&model), &mut checkpoint)?;
        let emission = Arc::new(BurnupDecayNeutronEmission::new(
            data_id,
            model,
            self.neutron_yields.to_vec(),
        )?);
        checkpoint()?;
        Ok(BoundEvaluatedBurnupNeutrons {
            source: Arc::clone(self),
            binding,
            emission,
        })
    }
}

#[derive(Clone, Debug)]
pub struct BoundEvaluatedBurnupNeutrons {
    source: Arc<EvaluatedBurnupNeutrons>,
    binding: BoundEvaluatedBurnupDecay,
    emission: Arc<BurnupDecayNeutronEmission>,
}

impl BoundEvaluatedBurnupNeutrons {
    pub fn source(&self) -> &Arc<EvaluatedBurnupNeutrons> {
        &self.source
    }
    pub fn burnup_binding(&self) -> &BoundEvaluatedBurnupDecay {
        &self.binding
    }
    pub fn emission(&self) -> &Arc<BurnupDecayNeutronEmission> {
        &self.emission
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum EvaluatedBurnupNeutronError {
    Invalid(&'static str),
    Normalization {
        parent: usize,
        /// None denotes the component probability sum or weighted spectrum.
        component: Option<usize>,
        value: f64,
        tolerance: f64,
    },
    Assembly(DecayAssemblyError),
    Source(EndfError),
    Projection(DecayNeutronBindingError),
    Coupling(DiffusionBurnupError),
}

impl From<DecayAssemblyError> for EvaluatedBurnupNeutronError {
    fn from(error: DecayAssemblyError) -> Self {
        Self::Assembly(error)
    }
}
impl From<EndfError> for EvaluatedBurnupNeutronError {
    fn from(error: EndfError) -> Self {
        Self::Source(error)
    }
}
impl From<FissionError> for EvaluatedBurnupNeutronError {
    fn from(error: FissionError) -> Self {
        Self::Source(error.into())
    }
}
impl From<DecayNeutronBindingError> for EvaluatedBurnupNeutronError {
    fn from(error: DecayNeutronBindingError) -> Self {
        Self::Projection(error)
    }
}
impl From<DiffusionBurnupError> for EvaluatedBurnupNeutronError {
    fn from(error: DiffusionBurnupError) -> Self {
        Self::Coupling(error)
    }
}
impl core::fmt::Display for EvaluatedBurnupNeutronError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "evaluated burnup neutron source refused: {self:?}")
    }
}
impl std::error::Error for EvaluatedBurnupNeutronError {}

fn prepare(
    data_id: String,
    assembly: Arc<EvaluatedBurnupDecayAssembly>,
    kinetic_model: Arc<NeutronKineticModel>,
    spontaneous_sources: Vec<SpontaneousFissionPromptNeutronBinding>,
    confirmed_nonemitting_parents: &[usize],
    budget: EvaluatedBurnupNeutronBudget,
    checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
) -> Result<EvaluatedBurnupNeutrons, EvaluatedBurnupNeutronError> {
    let mut poll = || checkpoint();
    let checkpoint = &mut poll;
    checkpoint()?;
    if data_id.trim().is_empty()
        || !nonnegative(budget.normalization_tolerance)
        || budget.normalization_tolerance >= 1.0
        || assembly.model().groups().as_ref() != kinetic_model.diffusion().groups().as_ref()
    {
        return Err(invalid(
            "identified source, matching neutron grid and normalization allowance in [0,1) required",
        ));
    }
    let limits = budget.limits;
    let parents = assembly.nuclides().len();
    let groups = kinetic_model.group_energies_j().len();
    let components = spontaneous_sources.iter().try_fold(0_usize, |n, source| {
        n.checked_add(source.spectrum.components().len())
    });
    let values = groups
        .checked_add(4)
        .and_then(|n| n.checked_mul(parents))
        .and_then(|n| n.checked_add(groups.checked_mul(4)?))
        .and_then(|n| n.checked_add(components?))
        .and_then(|n| n.checked_add(spontaneous_sources.len().checked_mul(9)?));
    if groups > limits.max_groups
        || spontaneous_sources.len() > parents
        || confirmed_nonemitting_parents.len() > parents
        || values.is_none_or(|n| n > limits.max_group_values)
    {
        return Err(FissionError::Budget("evaluated burnup neutron storage").into());
    }
    let edges = kinetic_model.diffusion().groups().edges_j();
    let mut work = 0;
    let mut sf_by_parent = vec![None; parents];
    for (index, source) in spontaneous_sources.iter().enumerate() {
        visit(&mut work, limits, checkpoint)?;
        let parent = assembly.index(source.parent).ok_or_else(|| {
            invalid("SF neutron source parent is outside its evaluated inventory")
        })?;
        if sf_by_parent[parent].replace(index).is_some()
            || source.multiplicity.kind() != FissionNeutronKind::Prompt
            || source.spectrum.kind() != FissionNeutronKind::Prompt
            || source.multiplicity.metadata() != source.spectrum.metadata()
        {
            return Err(invalid(
                "unique matched prompt SF yield and spectrum required; total/delayed sources cannot substitute",
            ));
        }
        let original = assembly.evaluations()[parent].metadata();
        let metadata = source.multiplicity.metadata();
        if metadata.nuclide != original.nuclide
            || metadata.nuclide != source.parent
            || metadata.material_number != original.material_number
            || metadata.atomic_weight_ratio != original.atomic_weight_ratio
        {
            return Err(invalid(
                "SF neutron data and original MF8 parent identity or mass disagree",
            ));
        }
    }
    let mut confirmed = vec![false; parents];
    for &parent in confirmed_nonemitting_parents {
        visit(&mut work, limits, checkpoint)?;
        if parent >= parents || confirmed[parent] {
            return Err(invalid(
                "nonemitting declarations require unique original parent indices",
            ));
        }
        confirmed[parent] = true;
    }
    let mut normalization = vec![None; spontaneous_sources.len()];
    let mut coverage = Vec::with_capacity(parents);
    let mut projection = Vec::with_capacity(parents);
    let mut neutron_yields = Vec::with_capacity(parents);
    for (parent, evaluation) in assembly.evaluations().iter().enumerate() {
        visit(&mut work, limits, checkpoint)?;
        let mut all_branches = Accumulator::new(0.0);
        let mut sf_branches = Accumulator::new(0.0);
        let mut ordinary_neutron_mode = false;
        for mode in evaluation.modes() {
            visit(&mut work, limits, checkpoint)?;
            all_branches.add(mode.probability)?;
            if mode.rtyp == 6.0 {
                sf_branches.add(mode.probability)?;
            } else if mode.probability > 0.0 && mode.rtyp.to_string().contains('5') {
                ordinary_neutron_mode = true;
            }
        }
        let stable = evaluation.law().decay_constant().value() == 0.0;
        let sf = sf_branches.value();
        let all = all_branches.value();
        let spectrum = evaluation
            .spectra()
            .iter()
            .find(|spectrum| spectrum.radiation == DecayRadiation::Neutron);
        if confirmed[parent] && (stable || sf > 0.0 || ordinary_neutron_mode || spectrum.is_some())
        {
            return Err(invalid(
                "nonemitting declarations apply only to unreported parents without neutron-producing modes",
            ));
        }
        if (sf > 0.0) != sf_by_parent[parent].is_some() {
            return Err(invalid(
                "every positive SF parent needs exactly one prompt source, and no unused source is admitted",
            ));
        }
        let mut counts = vec![PositiveIntegral::new(); groups];
        let mut energies = counts.clone();
        let mut ordinary_reported = false;
        if !stable {
            if let Some(spectrum) = spectrum {
                if sf > 0.0 {
                    // Empty SF STYP5 records contain an average only. They do
                    // not replace MF5 or prove an ordinary beta-neutron zero.
                    ordinary_reported =
                        admit_mixed_spectrum(evaluation, spectrum, &mut work, limits, checkpoint)?;
                } else {
                    ordinary_reported = true;
                }
                if ordinary_reported {
                    integrate_spectrum(
                        spectrum,
                        edges,
                        &mut counts,
                        &mut energies,
                        &mut work,
                        limits,
                        checkpoint,
                    )?;
                }
            }
            if ordinary_neutron_mode && !ordinary_reported {
                return Err(DecayNeutronBindingError::Unreported {
                    parent,
                    nuclide: assembly.nuclides()[parent],
                }
                .into());
            }
            if sf == 0.0 && !ordinary_reported && !confirmed[parent] {
                return Err(DecayNeutronBindingError::Unreported {
                    parent,
                    nuclide: assembly.nuclides()[parent],
                }
                .into());
            }
            if ordinary_neutron_mode {
                let mut particles = PositiveIntegral::new();
                for count in &counts {
                    particles.add_product(1.0, count.value()?, 1.0, 1.0)?;
                }
                if particles.value()? == 0.0 {
                    return Err(invalid(
                        "a positive ordinary neutron-emitting mode cannot have empty neutron coverage",
                    ));
                }
            }
            if let Some(source) = sf_by_parent[parent] {
                normalization[source] = Some(integrate_sf(
                    parent,
                    &spontaneous_sources[source],
                    sf,
                    all,
                    edges,
                    &mut counts,
                    &mut energies,
                    &mut work,
                    budget,
                    checkpoint,
                )?);
            }
        }
        coverage.push(if stable {
            DecayEmissionCoverage::Stable
        } else if ordinary_reported || sf > 0.0 {
            DecayEmissionCoverage::Reported
        } else {
            DecayEmissionCoverage::Unreported
        });
        let counts = counts
            .into_iter()
            .map(PositiveIntegral::value)
            .collect::<Result<Vec<_>, _>>()?;
        let energies = energies
            .into_iter()
            .map(PositiveIntegral::value)
            .collect::<Result<Vec<_>, _>>()?;
        let mut available = PositiveIntegral::new();
        if !stable {
            let heavy = evaluation.average_energies()[2]
                .value_j
                .ok_or_else(|| invalid("unknown original heavy-particle energy"))?;
            available.add_product(
                1.0,
                heavy,
                assembly.routing()[parent].heavy_particles.radiation,
                1.0,
            )?;
        }
        let (row, moments) = project_group_moments(
            parent,
            &counts,
            &energies,
            kinetic_model.group_energies_j(),
            available.value()?,
            limits,
            &mut work,
            checkpoint,
        )?;
        neutron_yields.push(row);
        projection.push(moments);
    }
    checkpoint()?;
    Ok(EvaluatedBurnupNeutrons {
        data_id: data_id.into(),
        assembly,
        kinetic_model,
        spontaneous_sources: spontaneous_sources.into(),
        normalization: normalization
            .into_iter()
            .map(|item| item.expect("every SF source was used"))
            .collect::<Vec<_>>()
            .into(),
        coverage: coverage.into(),
        confirmed_nonemitting_parents: confirmed_nonemitting_parents.into(),
        projection: projection.into(),
        neutron_yields: neutron_yields.into(),
    })
}

fn admit_mixed_spectrum(
    evaluation: &DecayEvaluation,
    spectrum: &DecaySpectrum,
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EvaluatedBurnupNeutronError> {
    let mut reported = false;
    let mut ordinary_modes = BTreeSet::new();
    for mode in evaluation.modes() {
        visit(work, limits, checkpoint)?;
        if mode.probability > 0.0 && !mode.rtyp.to_string().contains('6') {
            ordinary_modes.insert(mode.rtyp.to_bits());
        }
    }
    let known = |rtyp: f64| ordinary_modes.contains(&rtyp.to_bits());
    if spectrum.discrete_normalization.value > 0.0 {
        for line in spectrum.discrete.iter() {
            visit(work, limits, checkpoint)?;
            if line.intensity.value > 0.0 {
                if !known(line.decay_mode) {
                    return Err(invalid(
                        "MF8 neutron line on an SF parent is not an identified disjoint ordinary decay mode",
                    ));
                }
                reported = true;
            }
        }
    }
    if spectrum.continuum_normalization.value > 0.0 {
        let continuum = spectrum.continuum.as_ref().ok_or_else(|| {
            invalid("positive ordinary neutron normalization requires a continuum")
        })?;
        let mut positive = false;
        for point in continuum.points.iter() {
            visit(work, limits, checkpoint)?;
            positive |= point.density_per_j > 0.0;
        }
        if positive {
            if !known(continuum.decay_mode) {
                return Err(invalid(
                    "MF8 neutron continuum on an SF parent is not an identified disjoint ordinary decay mode",
                ));
            }
            reported = true;
        }
    }
    Ok(reported)
}

fn integrate_sf(
    parent: usize,
    source: &SpontaneousFissionPromptNeutronBinding,
    sf_branch_sum: f64,
    all_branch_sum: f64,
    edges: &[f64],
    counts: &mut [PositiveIntegral],
    energies: &mut [PositiveIntegral],
    work: &mut usize,
    budget: EvaluatedBurnupNeutronBudget,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<SpontaneousFissionNeutronNormalization, EvaluatedBurnupNeutronError> {
    let fraction = sf_branch_sum / all_branch_sum;
    if !fraction.is_finite() || fraction <= 0.0 {
        return Err(FissionError::Arithmetic.into());
    }
    let mut probability_sum = Accumulator::new(0.0);
    let mut weighted_integral = PositiveIntegral::new();
    let mut integrals = Vec::with_capacity(source.spectrum.components().len());
    let nu = source.multiplicity.yield_per_fission();
    for (component_index, component) in source.spectrum.components().iter().enumerate() {
        visit(work, budget.limits, checkpoint)?;
        probability_sum.add(component.probability)?;
        let mut component_counts = vec![PositiveIntegral::new(); counts.len()];
        let mut component_energies = component_counts.clone();
        integrate_continuum(
            DecayRadiation::Neutron,
            &component.points,
            &component.interpolation,
            1.0,
            edges,
            &mut component_counts,
            &mut component_energies,
            work,
            budget.limits,
            checkpoint,
        )?;
        let mut integral = PositiveIntegral::new();
        for group in 0..counts.len() {
            visit(work, budget.limits, checkpoint)?;
            let count = component_counts[group].value()?;
            let energy = component_energies[group].value()?;
            integral.add_product(1.0, count, 1.0, 1.0)?;
            counts[group].add_triple_product(fraction, nu, component.probability, count, 1.0)?;
            energies[group].add_triple_product(fraction, nu, component.probability, energy, 1.0)?;
        }
        let integral = integral.value()?;
        admit_normalization(
            parent,
            Some(component_index),
            integral,
            budget.normalization_tolerance,
        )?;
        weighted_integral.add_product(1.0, component.probability, integral, 1.0)?;
        integrals.push(integral);
    }
    let probability_sum = probability_sum.value();
    let weighted_integral = weighted_integral.value()?;
    admit_normalization(
        parent,
        None,
        probability_sum,
        budget.normalization_tolerance,
    )?;
    admit_normalization(
        parent,
        None,
        weighted_integral,
        budget.normalization_tolerance,
    )?;
    let mut emitted = PositiveIntegral::new();
    emitted.add_product(1.0, nu, weighted_integral, 1.0)?;
    let emitted = emitted.value()?;
    Ok(SpontaneousFissionNeutronNormalization {
        spontaneous_branch_sum: sf_branch_sum,
        total_decay_branch_sum: all_branch_sum,
        fission_events_per_decay: fraction,
        component_probability_sum: probability_sum,
        component_integrals: integrals.into(),
        weighted_spectrum_integral: weighted_integral,
        declared_particles_per_fission: nu,
        integrated_particles_per_fission: emitted,
        multiplicity_residual: sum(&[emitted, -nu]),
    })
}

fn admit_normalization(
    parent: usize,
    component: Option<usize>,
    value: f64,
    tolerance: f64,
) -> Result<(), EvaluatedBurnupNeutronError> {
    if !nonnegative(value) || (value - 1.0).abs() > tolerance {
        return Err(EvaluatedBurnupNeutronError::Normalization {
            parent,
            component,
            value,
            tolerance,
        });
    }
    Ok(())
}

fn invalid(message: &'static str) -> EvaluatedBurnupNeutronError {
    EvaluatedBurnupNeutronError::Invalid(message)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
pub(crate) use tests::material_fixture;
