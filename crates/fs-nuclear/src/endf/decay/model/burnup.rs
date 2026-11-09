//! Evaluated decay channels with multiple immediate nuclear products.
//!
//! Spontaneous fission is an original BurnupRate::Decay channel with supplied
//! independent fragment yields. Its branching belongs to the parent hazard;
//! product multiplicities are per fission and are not weighted a second time.
//! The single-residual DecayModel is neither fabricated nor broadened here.
//! All evolution, irradiation, movement and energy remain the original burnup,
//! neutron and material owners' work.
//!
//! ENDF MF8/MT457 average energy is per parent transformation, summed over all
//! modes. The existing explicit parent-average routing is retained for each
//! branch; a mode's Q does not supply extra heat. MF8/MT454 independent products
//! precede their subsequent decays, unlike cumulative MT459 products. Fragment
//! yields establish neither neutron multiplicity nor neutron spectrum; these
//! require separately supplied data, including the MF1/MF5 SF data where used.
//! See ENDF-102 (2023), sections 8.3–8.4:
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use super::{
    DecayAssemblyError, DecayEnergyRouting, DecayModelEntry, DecayProductRetention,
    products::{RetentionBinding, resolve_mode},
    routed_energy,
};
use crate::burnup::transport::BurnupTransportStep;
use crate::burnup::{
    BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct, BurnupRate,
    BurnupReaction, BurnupStep,
};
use crate::endf::decay::DecayEvaluation;
use crate::endf::fission_yields::spontaneous::SpontaneousFissionYieldEvaluation;
use crate::irradiation::{NeutronGroups, NuclideId};
use crate::{Accumulator, FissionError, nonnegative, sum};

/// Evaluated prompt-fission and ordinary decay neutrons on the native inventory.
pub mod neutrons;

/// Independent products for one original spontaneously fissioning parent.
/// Every positive source yield is represented or explicitly exported. Export
/// means outside this inventory model; it does not make a product stable or
/// automatically deposit/escape its future decay energy.
#[derive(Clone, Debug)]
pub struct SpontaneousFissionProductBinding {
    pub parent: NuclideId,
    pub yields: Arc<SpontaneousFissionYieldEvaluation>,
    /// Unique, unrepresented source products in caller declaration order.
    pub exported_products: Vec<NuclideId>,
}

/// Exact source mode of an original native spontaneous-reaction channel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EvaluatedDecayModeBinding {
    pub evaluation_index: usize,
    pub mode_index: usize,
    pub reaction_index: usize,
    pub spontaneous_fission_binding: Option<usize>,
}

/// Immutable source data and the original admitted multi-product reaction model.
/// It creates no new inventory, propagator or material energy account.
#[derive(Clone, Debug)]
pub struct EvaluatedBurnupDecayAssembly {
    model: Arc<BurnupModel>,
    evaluations: Arc<[Arc<DecayEvaluation>]>,
    routing: Arc<[DecayEnergyRouting]>,
    nuclides: Arc<[NuclideId]>,
    indices: BTreeMap<NuclideId, usize>,
    fission_products: Arc<[SpontaneousFissionProductBinding]>,
    exported_indices: Arc<[Vec<usize>]>,
    mode_sources: Arc<[EvaluatedDecayModeBinding]>,
    product_retention: RetentionBinding,
    limits: BurnupLimits,
}

impl EvaluatedBurnupDecayAssembly {
    /// Admit complete evaluated decay channels directly through BurnupModel.
    /// Ordinary modes retain their original residual resolver. Pure RTYP=6 uses
    /// explicitly supplied independent SF products; compound SF modes and
    /// nonzero SF RFS (without one interpretable residual) refuse.
    /// No isotope, branch, spectrum, product, temperature domain or stable
    /// endpoint is synthesized. Supplied zero-probability modes create no event.
    pub fn new(
        data_id: String,
        entries: Vec<DecayModelEntry>,
        fission_products: Vec<SpontaneousFissionProductBinding>,
        groups: Arc<NeutronGroups>,
        temperature_range_k: [f64; 2],
        limits: BurnupLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DecayAssemblyError> {
        checkpoint()?;
        if data_id.trim().is_empty() || entries.is_empty() {
            return Err(invalid(
                "an identified nonempty evaluated burnup assembly is required",
            ));
        }
        if entries.len() > limits.propagation.maximum_nodes
            || entries.len() > fs_la::decay::MAX_DECAY_NODES
            || groups.len() > limits.max_neutron_groups
            || fission_products.len() > entries.len()
        {
            return Err(DecayAssemblyError::Budget("evaluated burnup source shape"));
        }
        let mut indices = BTreeMap::new();
        let mut modes = 0_usize;
        for (index, entry) in entries.iter().enumerate() {
            checkpoint()?;
            let id = entry.evaluation.metadata().nuclide;
            if indices.insert(id, index).is_some() {
                return Err(DecayAssemblyError::DuplicateNuclide(id));
            }
            modes = modes
                .checked_add(entry.evaluation.modes().len())
                .ok_or(DecayAssemblyError::Budget("evaluated burnup decay modes"))?;
            if modes > limits.max_reactions {
                return Err(DecayAssemblyError::Budget("evaluated burnup decay modes"));
            }
            entry.routing.light_particles.validate()?;
            entry.routing.electromagnetic.validate()?;
            entry.routing.heavy_particles.validate()?;
        }
        let mut source_by_parent = BTreeMap::new();
        let mut retained = Vec::with_capacity(fission_products.len());
        let mut exported_indices = Vec::with_capacity(fission_products.len());
        for (binding_index, binding) in fission_products.iter().enumerate() {
            checkpoint()?;
            if binding.yields.metadata().nuclide != binding.parent
                || !indices.contains_key(&binding.parent)
                || source_by_parent.insert(binding.parent, binding_index).is_some()
            {
                return Err(invalid(
                    "SF yields require a unique original represented parent and isomer",
                ));
            }
            if binding.exported_products.len() > binding.yields.products().len() {
                return Err(invalid("SF exports must be original evaluated products"));
            }
            let mut exports = BTreeSet::new();
            let mut export_indices = Vec::with_capacity(binding.exported_products.len());
            for &id in &binding.exported_products {
                checkpoint()?;
                let index = binding
                    .yields
                    .products()
                    .binary_search(&id)
                    .map_err(|_| invalid("an SF export is absent from its original yield source"))?;
                if indices.contains_key(&id) || !exports.insert(id) {
                    return Err(invalid(
                        "SF exports must be unique and absent from the represented inventory",
                    ));
                }
                export_indices.push(index);
            }
            let mut products = Vec::new();
            for (&id, value) in binding
                .yields
                .products()
                .iter()
                .zip(binding.yields.values())
            {
                checkpoint()?;
                if value.yield_per_fission == 0.0 {
                    continue;
                }
                if let Some(&nuclide) = indices.get(&id) {
                    if products.len() >= limits.max_products {
                        return Err(DecayAssemblyError::Budget(
                            "retained spontaneous-fission products",
                        ));
                    }
                    products.push(BurnupProduct {
                        nuclide,
                        yield_per_event: value.yield_per_fission,
                    });
                } else if !exports.contains(&id) {
                    return Err(DecayAssemblyError::MissingDaughter {
                        parent: binding.parent,
                        daughter: id,
                    });
                }
            }
            retained.push(products);
            exported_indices.push(export_indices);
        }
        let mut used = vec![false; fission_products.len()];
        let mut nuclides = Vec::with_capacity(entries.len());
        let mut reactions = Vec::with_capacity(modes);
        let mut mode_sources = Vec::with_capacity(modes);
        let mut products_used = 0_usize;
        for (parent, entry) in entries.iter().enumerate() {
            checkpoint()?;
            let source = &entry.evaluation;
            let id = source.metadata().nuclide;
            nuclides.push(BurnupNuclide { id, law: source.law() });
            let energy_per_event = routed_energy(source, entry.routing)?;
            for (mode_index, mode) in source.modes().iter().enumerate() {
                checkpoint()?;
                if mode.probability == 0.0 {
                    continue;
                }
                let (products, spontaneous_fission_binding) = if mode.rtyp == 6.0 {
                    if mode.daughter_isomer != 0 {
                        return Err(invalid(
                            "nonzero spontaneous-fission RFS has no supported single-residual interpretation",
                        ));
                    }
                    let &binding = source_by_parent.get(&id).ok_or_else(|| {
                        invalid("a positive spontaneous-fission mode needs independent product yields")
                    })?;
                    used[binding] = true;
                    (retained[binding].clone(), Some(binding))
                } else {
                    let resolved = resolve_mode(id, mode.rtyp, mode.daughter_isomer)?;
                    let &daughter = indices
                        .get(&resolved.residual)
                        .ok_or(DecayAssemblyError::MissingDaughter {
                            parent: id,
                            daughter: resolved.residual,
                        })?;
                    (
                        RetentionBinding::default()
                            .products(daughter, resolved.charged)
                            .as_slice()
                            .to_vec(),
                        None,
                    )
                };
                products_used = products_used
                    .checked_add(products.len())
                    .ok_or(DecayAssemblyError::Budget("evaluated burnup products"))?;
                if products_used > limits.max_products
                    || products_used > limits.propagation.maximum_edges
                {
                    return Err(DecayAssemblyError::Budget("evaluated burnup products"));
                }
                let mut label = format!(
                    "{} | {} | MAT{} MF8/MT457 mode{}",
                    data_id,
                    source.data_id(),
                    source.metadata().material_number,
                    mode_index,
                );
                if let Some(binding) = spontaneous_fission_binding {
                    let yields = &fission_products[binding].yields;
                    label.push_str(&format!(
                        " | {} | MAT{} MF8/MT454 independent SF",
                        yields.data_id(),
                        yields.metadata().material_number,
                    ));
                }
                let reaction_index = reactions.len();
                reactions.push(BurnupReaction {
                    data_id: label,
                    parent,
                    rate: BurnupRate::Decay { branch_fraction: mode.probability },
                    products,
                    energy_per_event,
                });
                mode_sources.push(EvaluatedDecayModeBinding {
                    evaluation_index: parent,
                    mode_index,
                    reaction_index,
                    spontaneous_fission_binding,
                });
            }
        }
        if used.iter().any(|used| !used) {
            return Err(invalid(
                "an SF product binding has no positive pure spontaneous-fission mode",
            ));
        }
        let model = Arc::new(BurnupModel::new(
            BurnupData {
                data_id,
                groups,
                temperature_min_k: temperature_range_k[0],
                temperature_max_k: temperature_range_k[1],
                nuclides,
                reactions,
            },
            limits,
        )?);
        let (evaluations, routing): (Vec<_>, Vec<_>) = entries
            .into_iter()
            .map(|entry| (entry.evaluation, entry.routing))
            .unzip();
        let nuclides = model
            .nuclides()
            .iter()
            .map(|nuclide| nuclide.id)
            .collect::<Vec<_>>();
        checkpoint()?;
        Ok(Self {
            model,
            evaluations: evaluations.into(),
            routing: routing.into(),
            nuclides: nuclides.into(),
            indices,
            fission_products: fission_products.into(),
            exported_indices: exported_indices.into(),
            mode_sources: mode_sources.into(),
            product_retention: RetentionBinding::default(),
            limits,
        })
    }

    pub fn model(&self) -> &Arc<BurnupModel> { &self.model }
    pub fn evaluations(&self) -> &[Arc<DecayEvaluation>] { &self.evaluations }
    pub fn routing(&self) -> &[DecayEnergyRouting] { &self.routing }
    pub fn nuclides(&self) -> &[NuclideId] { &self.nuclides }
    pub fn index(&self, nuclide: NuclideId) -> Option<usize> {
        self.indices.get(&nuclide).copied()
    }
    pub fn spontaneous_fission_products(&self) -> &[SpontaneousFissionProductBinding] {
        &self.fission_products
    }
    pub fn mode_sources(&self) -> &[EvaluatedDecayModeBinding] { &self.mode_sources }
    pub fn burnup_product_retention(&self) -> DecayProductRetention {
        self.product_retention.policy()
    }

    /// Copy the original native data for explicitly appended neutron reactions.
    /// No source, channel identity, range or energy is regenerated from labels.
    pub fn burnup_data(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupData, DecayAssemblyError> {
        clone_data(&self.model, &mut checkpoint)
    }

    /// Retain ordinary emitted He4/H1 through the same admission and product
    /// merger as the residual-lineage assembly. SF fragments remain exactly the
    /// independent yield source; this policy does not infer extra SF products.
    pub fn with_burnup_product_retention(
        mut self,
        policy: DecayProductRetention,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DecayAssemblyError> {
        checkpoint()?;
        let retention = RetentionBinding::new(
            policy, &self.indices, &self.evaluations, &mut checkpoint,
        )?;
        let mut data = clone_data(&self.model, &mut checkpoint)?;
        for source in self.mode_sources.iter() {
            checkpoint()?;
            if source.spontaneous_fission_binding.is_none() {
                let evaluation = &self.evaluations[source.evaluation_index];
                let mode = evaluation.modes()[source.mode_index];
                let resolved = resolve_mode(
                    evaluation.metadata().nuclide, mode.rtyp, mode.daughter_isomer,
                )?;
                let daughter = self.indices[&resolved.residual];
                data.reactions[source.reaction_index].products = retention
                    .products(daughter, resolved.charged)
                    .as_slice()
                    .to_vec();
            }
        }
        let model = Arc::new(BurnupModel::new(data, self.limits)?);
        checkpoint()?;
        self.model = model;
        self.product_retention = retention;
        Ok(self)
    }

    /// Bind every original decay channel to an admitted reaction model. Added
    /// neutron reactions and reordered channels are permitted; isotope order,
    /// half-lives, group/domain declarations, products, labels, probabilities
    /// and energy routes must match exactly. Unknown/duplicate decays refuse.
    pub fn bind_burnup(
        self: &Arc<Self>,
        model: Arc<BurnupModel>,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BoundEvaluatedBurnupDecay, DecayAssemblyError> {
        checkpoint()?;
        if model.nuclides() != self.model.nuclides()
            || model.groups() != self.model.groups()
            || model.temperature_min_k() != self.model.temperature_min_k()
            || model.temperature_max_k() != self.model.temperature_max_k()
        {
            return Err(invalid(
                "bound burnup identities, laws, grid or temperature domain differ from the evaluated source",
            ));
        }
        if model.reactions().len() > self.limits.max_reactions {
            return Err(DecayAssemblyError::Budget("bound evaluated decay channels"));
        }
        let mut expected = BTreeMap::new();
        for (index, reaction) in self.model.reactions().iter().enumerate() {
            checkpoint()?;
            if expected
                .insert((reaction.parent, reaction.data_id.as_str()), index)
                .is_some()
            {
                return Err(invalid("ambiguous original evaluated decay channel"));
            }
        }
        let mut reaction_indices = vec![0; expected.len()];
        for (index, reaction) in model.reactions().iter().enumerate() {
            checkpoint()?;
            if let BurnupRate::Decay { .. } = reaction.rate {
                let original_index = expected
                    .remove(&(reaction.parent, reaction.data_id.as_str()))
                    .ok_or_else(|| invalid("bound burnup has an unknown or duplicate decay channel"))?;
                let original = &self.model.reactions()[original_index];
                if reaction.rate != original.rate
                    || reaction.energy_per_event != original.energy_per_event
                    || !same_products(&original.products, &reaction.products, &mut checkpoint)?
                {
                    return Err(invalid(
                        "bound decay branch, products or energy differ from the original evaluation",
                    ));
                }
                reaction_indices[original_index] = index;
            }
        }
        if !expected.is_empty() {
            return Err(invalid("bound burnup omits an evaluated decay channel"));
        }
        checkpoint()?;
        Ok(BoundEvaluatedBurnupDecay {
            source: Arc::clone(self),
            model,
            reaction_indices: reaction_indices.into(),
        })
    }
}

/// Read-only source interpretation of the original accepted native trajectory.
#[derive(Clone, Debug)]
pub struct BoundEvaluatedBurnupDecay {
    source: Arc<EvaluatedBurnupDecayAssembly>,
    model: Arc<BurnupModel>,
    reaction_indices: Arc<[usize]>,
}

/// Explicit bounds for source-observable arrays and noncancelling count errors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpontaneousFissionReportBudget {
    /// SF-event and exported-product scalar values together, across all cells.
    pub max_values: usize,
    pub max_work: usize,
    /// Parent-event and fragment-count arithmetic defects [expected nuclei].
    /// Applies only to this source projection: native trajectory defects remain
    /// in the original BurnupStep and are not certified a second time here.
    /// Evaluated yield uncertainties remain source metadata, not this allowance.
    pub nuclei_tolerance: f64,
}

/// Actual spontaneous fissions and independent exported births. All arrays use
/// source binding order, then each binding's declared exported-product order.
/// Subsequent exterior decay and its energy are outside this source report.
#[derive(Clone, Debug, PartialEq)]
pub struct SpontaneousFissionStep {
    pub duration_s: f64,
    pub fissions: Vec<f64>,
    pub exported_product_births: Vec<Vec<f64>>,
    /// Additional arithmetic defects of this source projection only [nuclei].
    pub nuclei_roundoff: f64,
}

impl BoundEvaluatedBurnupDecay {
    pub fn source(&self) -> &Arc<EvaluatedBurnupDecayAssembly> { &self.source }
    pub fn model(&self) -> &Arc<BurnupModel> { &self.model }
    /// Native channel index for every original source reaction, in source order.
    pub fn reaction_indices(&self) -> &[usize] { &self.reaction_indices }

    pub fn from_burnup_step(
        &self,
        model: &Arc<BurnupModel>,
        step: &BurnupStep,
        budget: SpontaneousFissionReportBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<SpontaneousFissionStep, DecayAssemblyError> {
        checkpoint()?;
        self.admit_report(model, 1, budget)?;
        self.apply_step(step, budget, &mut 0, &mut checkpoint)
    }

    /// Attribute events and exterior fragment births to the original event cell;
    /// all compartments share the same work, output and numerical allowances.
    pub fn from_transport_step(
        &self,
        model: &Arc<BurnupModel>,
        step: &BurnupTransportStep,
        budget: SpontaneousFissionReportBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<SpontaneousFissionStep>, DecayAssemblyError> {
        checkpoint()?;
        self.admit_report(model, step.burnup.len(), budget)?;
        if step.burnup.is_empty()
            || !nonnegative(step.duration_s)
            || step.burnup.iter().any(|local| local.duration_s != step.duration_s)
        {
            return Err(invalid("matching nonempty native transport report intervals required"));
        }
        let mut work = 0;
        let mut error = Accumulator::new(0.0);
        let mut reports = Vec::with_capacity(step.burnup.len());
        for local in &step.burnup {
            let report = self.apply_step(local, budget, &mut work, &mut checkpoint)?;
            error.add(report.nuclei_roundoff)?;
            reports.push(report);
        }
        if upper(error)? > budget.nuclei_tolerance {
            return Err(DecayAssemblyError::Budget("collective SF source count roundoff"));
        }
        checkpoint()?;
        Ok(reports)
    }

    fn admit_report(
        &self,
        model: &Arc<BurnupModel>,
        cells: usize,
        budget: SpontaneousFissionReportBudget,
    ) -> Result<(), DecayAssemblyError> {
        if !Arc::ptr_eq(model, &self.model) || !nonnegative(budget.nuclei_tolerance) {
            return Err(invalid(
                "SF reports require their exact native model and a finite nonnegative count allowance",
            ));
        }
        let values = self.source.exported_indices.iter()
            .try_fold(self.source.fission_products.len(), |total, row| total.checked_add(row.len()))
            .and_then(|values| values.checked_mul(cells))
            .ok_or(DecayAssemblyError::Budget("SF report values"))?;
        if values > budget.max_values {
            return Err(DecayAssemblyError::Budget("SF report values"));
        }
        // Bound per-cell container storage even for an ordinary-only source.
        if cells > budget.max_work {
            return Err(DecayAssemblyError::Budget("SF report compartments"));
        }
        Ok(())
    }

    fn apply_step(
        &self,
        step: &BurnupStep,
        budget: SpontaneousFissionReportBudget,
        work: &mut usize,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<SpontaneousFissionStep, DecayAssemblyError> {
        visit(work, budget, checkpoint)?;
        if !nonnegative(step.duration_s)
            || step.reaction_counts.len() != self.model.reactions().len()
        {
            return Err(invalid("matching native reaction counts and finite duration required"));
        }
        for &count in &step.reaction_counts {
            visit(work, budget, checkpoint)?;
            if !nonnegative(count) || (step.duration_s == 0.0 && count != 0.0) {
                return Err(invalid("finite physical event counts required"));
            }
        }
        let n = self.source.fission_products.len();
        let mut events = vec![Accumulator::new(0.0); n];
        let mut event_errors = vec![Accumulator::new(0.0); n];
        for source in self.source.mode_sources.iter() {
            visit(work, budget, checkpoint)?;
            if let Some(binding) = source.spontaneous_fission_binding {
                let count = step.reaction_counts[self.reaction_indices[source.reaction_index]];
                event_errors[binding].add(events[binding].add(count)?.abs())?;
            }
        }
        let mut fissions = Vec::with_capacity(n);
        let mut exported_product_births = Vec::with_capacity(n);
        let mut error = Accumulator::new(0.0);
        for (binding, (events, mut local_error)) in events.into_iter().zip(event_errors).enumerate() {
            visit(work, budget, checkpoint)?;
            let count = checked(events.value())?;
            local_error.add(checked(sum(&[events.hi, events.lo, -count]))?.abs())?;
            let count_error = upper(local_error)?;
            error.add(count_error)?;
            let mut products = Vec::with_capacity(self.source.exported_indices[binding].len());
            for &product in &self.source.exported_indices[binding] {
                visit(work, budget, checkpoint)?;
                let yield_ = self.source.fission_products[binding]
                    .yields.values()[product].yield_per_fission;
                let count = product_count(count, yield_, &mut error)?;
                let propagated = product_count(count_error, yield_, &mut error)?;
                error.add(propagated)?;
                products.push(count);
            }
            fissions.push(count);
            exported_product_births.push(products);
        }
        let nuclei_roundoff = upper(error)?;
        if nuclei_roundoff > budget.nuclei_tolerance {
            return Err(DecayAssemblyError::Budget("SF source count roundoff"));
        }
        checkpoint()?;
        Ok(SpontaneousFissionStep {
            duration_s: step.duration_s,
            fissions,
            exported_product_births,
            nuclei_roundoff,
        })
    }
}

fn same_products(
    original: &[BurnupProduct],
    bound: &[BurnupProduct],
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, DecayAssemblyError> {
    if original.len() != bound.len() {
        return Ok(false);
    }
    let mut expected = BTreeMap::new();
    for product in original {
        checkpoint()?;
        expected.insert(product.nuclide, product.yield_per_event);
    }
    for product in bound {
        checkpoint()?;
        if expected.remove(&product.nuclide) != Some(product.yield_per_event) {
            return Ok(false);
        }
    }
    Ok(expected.is_empty())
}

fn clone_data(
    model: &BurnupModel,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<BurnupData, DecayAssemblyError> {
    checkpoint()?;
    let mut nuclides = Vec::with_capacity(model.nuclides().len());
    for nuclide in model.nuclides() {
        checkpoint()?;
        nuclides.push(nuclide.clone());
    }
    let mut reactions = Vec::with_capacity(model.reactions().len());
    for reaction in model.reactions() {
        checkpoint()?;
        reactions.push(reaction.clone());
    }
    checkpoint()?;
    Ok(BurnupData {
        data_id: model.data_id().to_owned(),
        groups: Arc::clone(model.groups()),
        temperature_min_k: model.temperature_min_k(),
        temperature_max_k: model.temperature_max_k(),
        nuclides,
        reactions,
    })
}

fn product_count(a: f64, b: f64, error: &mut Accumulator) -> Result<f64, DecayAssemblyError> {
    let value = checked(a * b)?;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        return Err(FissionError::Arithmetic.into());
    }
    error.add(a.mul_add(b, -value).abs())?;
    Ok(value)
}
fn checked(value: f64) -> Result<f64, DecayAssemblyError> {
    if value.is_finite() { Ok(value) } else { Err(FissionError::Arithmetic.into()) }
}
fn upper(value: Accumulator) -> Result<f64, DecayAssemblyError> {
    let result = checked(value.value())?;
    Ok(if sum(&[value.hi, value.lo, -result]) > 0.0 {
        checked(result.next_up())?
    } else {
        result
    })
}
fn visit(
    work: &mut usize,
    budget: SpontaneousFissionReportBudget,
    checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
) -> Result<(), DecayAssemblyError> {
    checkpoint()?;
    *work = work.checked_add(1).ok_or(DecayAssemblyError::Budget("SF report work"))?;
    if *work > budget.max_work {
        return Err(DecayAssemblyError::Budget("SF report work"));
    }
    Ok(())
}
fn invalid(message: &'static str) -> DecayAssemblyError {
    DecayAssemblyError::Invalid(message)
}

#[cfg(test)]
mod tests;
