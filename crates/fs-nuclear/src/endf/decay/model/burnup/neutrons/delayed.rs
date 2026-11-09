//! Evaluated SF delayed families on the original finite radioactive inventory.
//!
//! Preparation retains the NSUB=4 MF1/MT455 family order and MF5/MT455 laws.
//! A component with probability p and integrated density A creates nu*p*A
//! precursors per SF event. Only its conditional native spectrum divides by A;
//! source probabilities, density integrals and total multiplicity are not repaired.
//!
//! Attachment explicitly reserves existing outward heavy-particle energy.
//! Every positive independent SF fragment must lie outside the retained nuclide
//! cohort, preventing a second delayed-emission history for the same products.
//! The original native precursor storage, movement and decay own all evolution.

use super::*;
use crate::diffusion::DelayedGroup;
use crate::irradiation::NeutronGroups;

/// Prepared source families, independent of their final native family indices.
/// Build these against the intended grid/representative energies, then use
/// families() to construct the final explicit-energy neutron model.
#[derive(Clone, Debug)]
pub struct EvaluatedSpontaneousFissionDelayedNeutrons {
    data_id: Arc<str>,
    multiplicity: Arc<SpontaneousFissionNeutronYield>,
    spectrum: Arc<SpontaneousFissionNeutronSpectrum>,
    groups: Arc<NeutronGroups>,
    representative_energies_j: Arc<[f64]>,
    families: Arc<[DelayedGroup]>,
    precursor_yields_per_fission: Arc<[f64]>,
    component_integrals: Arc<[f64]>,
    component_projections: Arc<[DecayNeutronProjection]>,
    component_probability_sum: f64,
    weighted_spectrum_integral: f64,
    source_energy_j_per_fission: f64,
}

impl EvaluatedSpontaneousFissionDelayedNeutrons {
    pub fn new(
        data_id: String,
        multiplicity: Arc<SpontaneousFissionNeutronYield>,
        spectrum: Arc<SpontaneousFissionNeutronSpectrum>,
        kinetic_grid: &NeutronKineticModel,
        budget: EvaluatedBurnupNeutronBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EvaluatedBurnupNeutronError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || multiplicity.kind() != FissionNeutronKind::Delayed
            || spectrum.kind() != FissionNeutronKind::Delayed
            || multiplicity.metadata() != spectrum.metadata()
            || !nonnegative(budget.normalization_tolerance)
            || budget.normalization_tolerance >= 1.0
        {
            return Err(invalid(
                "identified matching delayed MF1/MF5 sources and a finite allowance required",
            ));
        }
        let groups = kinetic_grid.diffusion().groups();
        let energies = kinetic_grid.group_energies_j();
        let components = spectrum.components();
        let lambdas = multiplicity.decay_constants_s_inv();
        let values = groups
            .len()
            .checked_add(7)
            .and_then(|width| width.checked_mul(components.len()))
            .and_then(|retained| retained.checked_add(groups.len().checked_mul(7)?))
            .and_then(|values| values.checked_add(3))
            .ok_or(FissionError::Budget("evaluated delayed family values"))?;
        if components.is_empty()
            || components.len() != lambdas.len()
            || components.len() > fs_la::decay::MAX_DECAY_NODES
            || groups.len() > budget.limits.max_groups
            || values > budget.limits.max_group_values
        {
            return Err(FissionError::Budget("evaluated delayed family shape").into());
        }
        let mut work = 0;
        let mut probability = Accumulator::new(0.0);
        let mut weighted_integral = PositiveIntegral::new();
        let mut source_energy = PositiveIntegral::new();
        let mut families = Vec::with_capacity(components.len());
        let mut yields = Vec::with_capacity(components.len());
        let mut integrals = Vec::with_capacity(components.len());
        let mut projections = Vec::with_capacity(components.len());
        for (index, (component, &lambda)) in components.iter().zip(lambdas).enumerate() {
            visit(&mut work, budget.limits, &mut checkpoint)?;
            if !lambda.is_finite() || lambda <= 0.0 || !nonnegative(component.probability) {
                return Err(invalid(
                    "positive finite native family rates and nonnegative probabilities required",
                ));
            }
            probability.add(component.probability)?;
            let mut counts = (0..groups.len())
                .map(|_| PositiveIntegral::new())
                .collect::<Vec<_>>();
            let mut energy = (0..groups.len())
                .map(|_| PositiveIntegral::new())
                .collect::<Vec<_>>();
            spectra::integrate_component(
                component,
                groups.edges_j(),
                &mut counts,
                &mut energy,
                &mut work,
                budget.limits,
                &mut checkpoint,
            )?;
            let particles = counts
                .into_iter()
                .map(|value| value.value())
                .collect::<Result<Vec<_>, _>>()?;
            let moments = energy
                .into_iter()
                .map(|value| value.value())
                .collect::<Result<Vec<_>, _>>()?;
            let area = sum(&particles);
            if !area.is_finite() || area <= 0.0 {
                return Err(invalid(
                    "each delayed family requires a positive finite conditional spectrum",
                ));
            }
            admit_normalization(0, Some(index), area, budget.normalization_tolerance)?;
            weighted_integral.add_product(1.0, component.probability, area, 1.0)?;
            source_energy.add_triple_product(
                1.0,
                multiplicity.yield_per_fission(),
                component.probability,
                sum(&moments),
                1.0,
            )?;
            let (projected, projection) = project_group_moments(
                0,
                &particles,
                &moments,
                energies,
                f64::MAX,
                budget.limits,
                &mut work,
                &mut checkpoint,
            )?;
            let mut conditional = Vec::with_capacity(groups.len());
            for &count in &projected {
                visit(&mut work, budget.limits, &mut checkpoint)?;
                conditional.push(scale(count, 1.0, area)?);
            }
            if (sum(&conditional) - 1.0).abs() > 64.0 * f64::EPSILON {
                return Err(invalid(
                    "projected conditional family does not satisfy native normalization",
                ));
            }
            let mut yield_ = PositiveIntegral::new();
            yield_.add_triple_product(
                1.0,
                multiplicity.yield_per_fission(),
                component.probability,
                area,
                1.0,
            )?;
            yields.push(yield_.value()?);
            families.push(DelayedGroup {
                data_id: format!(
                    "{} | MAT{} MF1/MT455 MF5/MT455 family{}",
                    data_id,
                    spectrum.metadata().material_number,
                    index,
                ),
                decay_constant_s_inv: lambda,
                spectrum: conditional,
            });
            integrals.push(area);
            projections.push(projection);
        }
        let probability_sum = probability.value();
        let weighted = weighted_integral.value()?;
        admit_normalization(0, None, probability_sum, budget.normalization_tolerance)?;
        admit_normalization(0, None, weighted, budget.normalization_tolerance)?;
        checkpoint()?;
        Ok(Self {
            data_id: data_id.into(),
            multiplicity,
            spectrum,
            groups: Arc::clone(groups),
            representative_energies_j: energies.to_vec().into(),
            families: families.into(),
            precursor_yields_per_fission: yields.into(),
            component_integrals: integrals.into(),
            component_projections: projections.into(),
            component_probability_sum: probability_sum,
            weighted_spectrum_integral: weighted,
            source_energy_j_per_fission: source_energy.value()?,
        })
    }

    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    pub fn multiplicity(&self) -> &Arc<SpontaneousFissionNeutronYield> {
        &self.multiplicity
    }
    pub fn spectrum(&self) -> &Arc<SpontaneousFissionNeutronSpectrum> {
        &self.spectrum
    }
    pub fn families(&self) -> &[DelayedGroup] {
        &self.families
    }
    pub fn precursor_yields_per_fission(&self) -> &[f64] {
        &self.precursor_yields_per_fission
    }
    pub fn component_integrals(&self) -> &[f64] {
        &self.component_integrals
    }
    /// Original unweighted component moments and conservative projection defects.
    pub fn component_projections(&self) -> &[DecayNeutronProjection] {
        &self.component_projections
    }
    pub fn component_probability_sum(&self) -> f64 {
        self.component_probability_sum
    }
    pub fn weighted_spectrum_integral(&self) -> f64 {
        self.weighted_spectrum_integral
    }
    /// Retained original spectrum energy, before the native family projection.
    pub fn source_energy_j_per_fission(&self) -> f64 {
        self.source_energy_j_per_fission
    }
}

/// Caller declaration that a bounded part of the original outward EHP account
/// funds delayed neutron storage. The source records do not establish this
/// partition themselves. No local heat, electromagnetic energy, or future
/// exported-fragment decay supplies additional funding.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DelayedSpontaneousFissionFunding {
    ParentHeavyRadiation { maximum_energy_j_per_fission: f64 },
}

/// Original prepared family order mapped uniquely into the final native model.
#[derive(Clone, Debug)]
pub struct DelayedSpontaneousFissionBinding {
    pub parent: NuclideId,
    pub source: Arc<EvaluatedSpontaneousFissionDelayedNeutrons>,
    pub native_family_indices: Vec<usize>,
    pub energy_funding: DelayedSpontaneousFissionFunding,
}

pub(super) fn attach(
    mut source: EvaluatedBurnupNeutrons,
    bindings: Vec<DelayedSpontaneousFissionBinding>,
    budget: EvaluatedBurnupNeutronBudget,
    checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
) -> Result<EvaluatedBurnupNeutrons, EvaluatedBurnupNeutronError> {
    let mut poll = || checkpoint();
    let checkpoint = &mut poll;
    checkpoint()?;
    if !source.delayed_sources.is_empty()
        || !nonnegative(budget.normalization_tolerance)
        || budget.normalization_tolerance >= 1.0
    {
        return Err(invalid(
            "delayed SF bindings are declared once with a finite normalization allowance",
        ));
    }
    let parents = source.assembly.nuclides().len();
    let native = source.kinetic_model.diffusion().delayed_groups();
    let values = native
        .len()
        .checked_add(1)
        .and_then(|width| parents.checked_mul(width))
        .ok_or(FissionError::Budget(
            "finite delayed precursor source values",
        ))?;
    if bindings.len() > parents
        || source.kinetic_model.group_energies_j().len() > budget.limits.max_groups
        || values > budget.limits.max_group_values
    {
        return Err(FissionError::Budget("finite delayed precursor source shape").into());
    }
    let mut work = 0;
    let mut precursor_yields = vec![vec![0.0; native.len()]; parents];
    let mut reserved = vec![0.0; parents];
    let mut used_parents = BTreeSet::new();
    let mut used_families = BTreeSet::new();
    for binding in &bindings {
        visit(&mut work, budget.limits, checkpoint)?;
        let prepared = &binding.source;
        let parent = source.assembly.index(binding.parent).ok_or_else(|| {
            invalid("delayed SF parent must be represented in the original inventory")
        })?;
        let evaluation = &source.assembly.evaluations()[parent];
        let metadata = prepared.multiplicity.metadata();
        if !used_parents.insert(parent)
            || binding.parent != metadata.nuclide
            || metadata.nuclide != evaluation.metadata().nuclide
            || metadata.material_number != evaluation.metadata().material_number
            || metadata.atomic_weight_ratio != evaluation.metadata().atomic_weight_ratio
            || prepared.groups.as_ref() != source.kinetic_model.diffusion().groups().as_ref()
            || prepared.representative_energies_j.as_ref()
                != source.kinetic_model.group_energies_j()
            || binding.native_family_indices.len() != prepared.families.len()
        {
            return Err(invalid(
                "delayed SF source identity, unique parent, original grid, energies and family map must match",
            ));
        }
        admit_normalization(
            parent,
            None,
            prepared.component_probability_sum,
            budget.normalization_tolerance,
        )?;
        admit_normalization(
            parent,
            None,
            prepared.weighted_spectrum_integral,
            budget.normalization_tolerance,
        )?;
        for (index, &area) in prepared.component_integrals.iter().enumerate() {
            visit(&mut work, budget.limits, checkpoint)?;
            admit_normalization(parent, Some(index), area, budget.normalization_tolerance)?;
        }
        let prompt = source
            .spontaneous_sources
            .iter()
            .position(|original| original.parent == binding.parent)
            .ok_or_else(|| {
                invalid("delayed SF retains the original matched prompt source and branch fraction")
            })?;
        let fraction = source.normalization[prompt].fission_events_per_decay;
        if !fraction.is_finite() || fraction <= 0.0 {
            return Err(invalid(
                "a delayed SF binding requires a positive original pure SF branch",
            ));
        }
        if prepared
            .precursor_yields_per_fission
            .iter()
            .any(|&yield_| yield_ > 0.0)
        {
            let fragments = source
                .assembly
                .spontaneous_fission_products()
                .iter()
                .find(|products| products.parent == binding.parent)
                .ok_or_else(|| {
                    invalid("delayed SF requires the original independent product boundary")
                })?;
            for (&id, value) in fragments
                .yields
                .products()
                .iter()
                .zip(fragments.yields.values())
            {
                visit(&mut work, budget.limits, checkpoint)?;
                // The immutable assembly already proves every positive
                // unrepresented source product is explicitly exported.
                if value.yield_per_fission > 0.0 && source.assembly.index(id).is_some() {
                    return Err(invalid(
                        "delayed SF families require every positive independent fragment to be explicitly exported",
                    ));
                }
            }
        }
        let mut native_energy = PositiveIntegral::new();
        for (family, &index) in binding.native_family_indices.iter().enumerate() {
            visit(&mut work, budget.limits, checkpoint)?;
            if !used_families.insert(index) || native.get(index) != Some(&prepared.families[family])
            {
                return Err(invalid(
                    "each delayed source family requires its unique exact native family",
                ));
            }
            let yield_ = prepared.precursor_yields_per_fission[family];
            native_energy.add_product(
                1.0,
                yield_,
                source.kinetic_model.precursor_energies_j()[index],
                1.0,
            )?;
            precursor_yields[parent][index] = scale(yield_, fraction, 1.0)?;
        }
        let native_energy = native_energy.value()?;
        let DelayedSpontaneousFissionFunding::ParentHeavyRadiation {
            maximum_energy_j_per_fission,
        } = binding.energy_funding;
        if !nonnegative(maximum_energy_j_per_fission)
            || native_energy > maximum_energy_j_per_fission
            || prepared.source_energy_j_per_fission > maximum_energy_j_per_fission
        {
            return Err(invalid(
                "both original and native delayed energy must fit the explicit per-fission funding cap",
            ));
        }
        let delayed_native = scale(native_energy, fraction, 1.0)?;
        let delayed_original = scale(prepared.source_energy_j_per_fission, fraction, 1.0)?;
        let mut immediate_native = PositiveIntegral::new();
        for (&count, &energy) in source.neutron_yields[parent]
            .iter()
            .zip(source.kinetic_model.group_energies_j())
        {
            visit(&mut work, budget.limits, checkpoint)?;
            immediate_native.add_product(1.0, count, energy, 1.0)?;
        }
        let immediate = immediate_native
            .value()?
            .max(source.projection[parent].kinetic_energy_j_per_decay);
        let heavy = evaluation.average_energies()[2]
            .value_j
            .ok_or_else(|| invalid("unknown original heavy-particle energy"))?;
        let available = scale(
            heavy,
            source.assembly.routing()[parent].heavy_particles.radiation,
            1.0,
        )?;
        let mut combined = Accumulator::new(immediate);
        combined.add(delayed_native.max(delayed_original))?;
        if combined.value() > available {
            return Err(invalid(
                "immediate and delayed neutron energy exceed the original outward heavy-particle account",
            ));
        }
        reserved[parent] = delayed_native;
    }
    checkpoint()?;
    source.delayed_sources = bindings.into();
    source.precursor_yields = precursor_yields.into();
    source.delayed_energy_j_per_decay = reserved.into();
    Ok(source)
}

fn scale(a: f64, b: f64, divisor: f64) -> Result<f64, EvaluatedBurnupNeutronError> {
    let mut result = PositiveIntegral::new();
    result.add_product(1.0, a, b, divisor)?;
    Ok(result.value()?)
}

#[cfg(test)]
pub(crate) use super::tests::delayed::material_delayed_sf_fixture;
