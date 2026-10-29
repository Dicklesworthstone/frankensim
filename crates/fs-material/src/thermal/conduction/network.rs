//! Lower prepared layer fluxes into the original conservative enthalpy residual.
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use super::{LayerConductionConfig, LayerConductionResponse, TwoLayerConduction};
use super::source::{LayerSourceError, ResolvedConductionLayer};
use fs_blake3::ContentHash;
use fs_matdb::PropertyLineSlice;
use fs_thermal::enthalpy::{CellGeometry, CellId, ContactPatch, DepositionError,
    DepositionSimulation, MaterialId, ThermalLaw, backward_euler::StepCoefficients};

/// Explicit original link and homogeneous distances, in the link's a/b order.
/// The contact resistance is caller-declared here; it is not inferred from the
/// existing conductance or confused with either bulk layer resistance.
pub struct LayerLinkBinding<'a, 'm> {
    pub link: u64,
    pub layers: [&'a ResolvedConductionLayer<'m>; 2],
    pub distances_m: [f64; 2],
    pub contact_area_resistance_m2_k_w: f64,
}
struct PreparedLink<'a, 'm> {
    binding: LayerLinkBinding<'a, 'm>,
    patch: ContactPatch,
    materials: [MaterialId; 2],
    geometry: [CellGeometry; 2],
}

/// Original heat graph with an opt-in subset of source-integrated links.
/// No state, source, time-stepper or second conservation owner is introduced.
/// Untargeted links retain their authored coefficients. Reprepare after contact,
/// material or geometry changes; h and the original thermal clock may evolve.
pub struct PreparedLayerConduction<'a, 'm> {
    links: Vec<PreparedLink<'a, 'm>>,
}

/// Physical interface response and the source charts used by that exact trial.
#[derive(Clone, Debug, PartialEq)]
pub struct LayerLinkSample {
    pub link: u64,
    pub cells: [CellId; 2],
    pub response: LayerConductionResponse,
    /// Actual G*(Ta-Tb) supplied to the original thermal kernel. It agrees with
    /// the integrated power to the declared heat tolerance, including roundoff.
    pub applied_power_a_to_b_w: f64,
    pub source_lines: [Arc<PropertyLineSlice>; 2],
    pub material_card_identities: [ContentHash; 2],
    pub caloric_curve_identities: [ContentHash; 2],
}

impl<'a, 'm> PreparedLayerConduction<'a, 'm> {
    /// Original link, source layers and geometric distances. Coupled consumers
    /// use these immutable bindings to admit the same physical cells and source
    /// choices before replacing an existing bulk-conduction approximation.
    pub fn bindings(&self) -> impl ExactSizeIterator<Item = &LayerLinkBinding<'a, 'm>> {
        self.links.iter().map(|link| &link.binding)
    }

    /// Bind bounded unique links to their actual caloric objects and coverage.
    pub fn new(thermal: &DepositionSimulation<'_>, mut bindings: Vec<LayerLinkBinding<'a, 'm>>,
        max_links: usize, mut checkpoint: impl FnMut() -> Result<(), DepositionError>)
        -> Result<Self, LayerSourceError> {
        checkpoint()?;
        if bindings.len() > max_links || bindings.len() > thermal.links().len() {
            return Err(DepositionError::Budget("integrated conduction links").into());
        }
        bindings.sort_unstable_by_key(|b| b.link);
        let mut ids = BTreeSet::new();
        let mut links = Vec::with_capacity(bindings.len());
        for binding in bindings {
            checkpoint()?;
            if !ids.insert(binding.link) { return Err(LayerSourceError::Invalid("duplicate integrated conduction link")); }
            let patch = thermal.contact_patch(binding.link)
                .ok_or(LayerSourceError::Invalid("inactive integrated conduction link"))?;
            let a = &thermal.cells()[&patch.link.a];
            let b = &thermal.cells()[&patch.link.b];
            let prepared = PreparedLink { binding, patch,
                materials: [a.material_id, b.material_id], geometry: [a.geometry, b.geometry] };
            prepared.validate(thermal)?;
            TwoLayerConduction::new(prepared.binding.layers.map(|layer| layer.conductivity_curve()),
                prepared.binding.distances_m, patch.link.area_m2, prepared.binding.contact_area_resistance_m2_k_w)?;
            links.push(prepared);
        }
        checkpoint()?;
        Ok(Self { links })
    }

    /// Endpoint coefficient callback for BackwardEulerStep::evaluate, or any
    /// coupled owner consuming the same coefficient maps. Recompute the physical
    /// interface at EVERY trial; do not freeze the initial sample for an interval.
    /// The returned secants are flux lowerings, not analytic Jacobian entries.
    /// The response separately supplies the full endpoint derivatives.
    pub fn resolve(&self, thermal: &DepositionSimulation<'_>, config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>)
        -> Result<StepCoefficients<Vec<LayerLinkSample>>, LayerSourceError> {
        checkpoint()?;
        self.links.len().checked_mul(config.max_iterations)
            .ok_or(DepositionError::Budget("integrated conduction work overflow"))?;
        let mut conductances_w_k = BTreeMap::new();
        let mut payload = Vec::with_capacity(self.links.len());
        for link in &self.links {
            checkpoint()?;
            let t = link.validate(thermal)?;
            let binding = &link.binding;
            let model = TwoLayerConduction::new(binding.layers.map(|layer| layer.conductivity_curve()),
                binding.distances_m, link.patch.link.area_m2, binding.contact_area_resistance_m2_k_w)?;
            let response = model.evaluate(t, config, &mut checkpoint)?;
            let applied = response.conductance_w_k * (t[0]-t[1]);
            let allowance = config.absolute_heat_tolerance_w
                + config.relative_heat_tolerance*applied.abs().max(response.power_a_to_b_w.abs());
            if !applied.is_finite() || !allowance.is_finite()
                || (applied-response.power_a_to_b_w).abs() > allowance {
                return Err(DepositionError::Arithmetic.into());
            }
            conductances_w_k.insert(binding.link, response.conductance_w_k);
            payload.push(LayerLinkSample { link: binding.link, cells: [link.patch.link.a, link.patch.link.b],
                response, applied_power_a_to_b_w: applied,
                source_lines: binding.layers.map(|layer| Arc::clone(layer.source_line())),
                material_card_identities: binding.layers.map(|layer| layer.caloric_curve().material_card_identity()),
                caloric_curve_identities: binding.layers.map(|layer| layer.caloric_curve().identity()) });
        }
        checkpoint()?;
        Ok(StepCoefficients { sources_w: BTreeMap::new(), conductances_w_k, payload })
    }
}
impl PreparedLink<'_, '_> {
    fn validate(&self, thermal: &DepositionSimulation<'_>) -> Result<[f64; 2], LayerSourceError> {
        if thermal.contact_patch(self.binding.link) != Some(self.patch) {
            return Err(LayerSourceError::Invalid("integrated layer contact changed; reprepare the binding"));
        }
        let mut temperatures = [0.0; 2];
        for (side, id) in [self.patch.link.a, self.patch.link.b].into_iter().enumerate() {
            let cell = thermal.cells().get(&id).ok_or(LayerSourceError::Invalid("inactive integrated layer cell"))?;
            let layer = self.binding.layers[side];
            if cell.geometry != self.geometry[side] || cell.material_id != self.materials[side]
                || !thermal.uses_material_law(cell.material_id, layer.caloric_curve() as &dyn ThermalLaw) {
                return Err(LayerSourceError::Invalid("integrated layers must retain original geometry and caloric objects"));
            }
            temperatures[side] = layer.temperature_at_enthalpy(cell.specific_h_j_kg)?;
        }
        Ok(temperatures)
    }
}
