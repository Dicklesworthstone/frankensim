//! Prescribed constitutive temperatures on retained longitudinal material memory.
//! Only source relaxation times may change. Density and every storage modulus
//! require validity-wide constant claims, not equality at two sampled endpoints.
//! The original pinned evaluator owns curves, units, context and refusal.
use super::*;
use crate::render::schedule::force::coupled::coefficients::rates::{
    MaterialRelaxationTime, MaterialRelaxationTimeReport,
};
use fs_matdb::{EvaluationDecision, MaterialCard, QueryPoint};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, ResolvedMaterialStatePoint};
use fs_qty::{QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// One original layer's prescribed temperature at an accepted sample boundary.
/// Supply every bound layer in `members()` order, then spatial layer order.
/// This is a constitutive condition, not a request to create or advance heat.
pub struct FrameRelaxationLayerTemperature<'a> {
    pub element: usize,
    pub layer: usize,
    pub card: &'a MaterialCard,
    /// Existing absolute-temperature coordinate; all other coordinates stay fixed.
    pub temperature_axis: &'a str,
    pub temperature_k: f64,
}

impl<F> FrameRelaxationBinding<'_, '_, F> {
    /// Re-evaluate the SAME layer claims at prescribed temperatures and update
    /// their original native Maxwell rates without erasing any accepted state.
    ///
    /// Every card, property key, source selector, coordinate convention and
    /// frequency support is retained. Missing/foreign/out-of-domain sources
    /// refuse; no automatic material swap or nominal coefficient fallback.
    /// Every density, Young/equilibrium modulus and branch increment must be
    /// a validity-wide scalar constant. Only time claims may have T dependence.
    /// Even coincident values from a temperature-dependent storage law refuse:
    /// that law needs its own internal-energy/entropy and thermal-work owner.
    ///
    /// All source queries and native preparation complete before either owner
    /// changes. Failure leaves the binding, every fluid momentum, material internal strain
    /// and the mechanical sample clock untouched. Foreign frames' branches remain
    /// unchanged. Source snapshots advance even when the numeric rates do not.
    /// The source model remains immutable reference geometry, not a second state.
    ///
    /// Hold conditions over each mechanical interval and split imposed events
    /// at boundaries. Smooth thermal histories require time refinement. This
    /// rate-only operation does NOT model thermal expansion, caloric storage,
    /// self-heating, temperature gradients inside a layer, or phase changes.
    pub fn set_layer_temperatures(
        &mut self,
        mechanical: &mut CoupledModalSystem,
        layers: &[FrameRelaxationLayerTemperature<'_>],
        gate: &CancelGate,
    ) -> Result<MaterialRelaxationTimeReport> {
        poll(gate)?;
        let count = self.members.iter().try_fold(0usize, |n, member|
            n.checked_add(member.source.elastic().layers().len())
                .ok_or_else(|| invalid("source temperature layer count overflow")))?;
        if layers.len() != count || mechanical.mode_count() != self.native_modes {
            return Err(invalid("temperature update needs every original layer and unchanged complete mechanical layout"));
        }
        let native = mechanical.material_relaxation_branches().get(self.memory.clone())
            .ok_or_else(|| invalid("temperature update lost original native material block"))?;
        if native.len() != self.branches.len() || native.len() != self.coordinates.len()
            || mechanical.material_memory().get(self.memory.clone()).is_none() {
            return Err(invalid("temperature update lost complete material histories"));
        }
        for (actual, expected) in native.iter().zip(&self.branches) {
            poll(gate)?;
            if actual.stiffness.to_bits() != expected.stiffness.to_bits()
                || actual.relaxation_time_s.to_bits() != expected.relaxation_time_s.to_bits()
                || actual.projection != expected.projection {
                return Err(invalid("native material coefficients changed outside their original source owner"));
            }
        }
        let nyquist = positive(core::f64::consts::PI / mechanical.sample_period_s())?;
        let mut refreshed = Vec::with_capacity(self.members.len());
        let mut at = 0;
        for member in &self.members {
            poll(gate)?;
            let section = member.source.elastic();
            let mut next = Vec::with_capacity(section.layers().len());
            for (layer_index, (layer, roles)) in section.layers().iter()
                .zip(member.source.properties()).enumerate() {
                poll(gate)?;
                let input = &layers[at];
                at += 1;
                if input.element != member.element || input.layer != layer_index {
                    return Err(invalid("temperature inputs must preserve original member and spatial layer order"));
                }
                require_constant_storage(&layer.material, roles)?;
                let point = temperature_query(&layer.material, input.temperature_axis, input.temperature_k)?;
                let material = layer.material.at_query_point(input.card, &point)?;
                require_constant_storage(&material, roles)?;
                next.push(ConcentricMaterialLayer { outer_radius_m: layer.outer_radius_m, material });
            }
            let source = ConcentricElasticSection::try_new(section.inner_radius_m(), next)?
                .with_prony_bending(member.source.properties())?;
            let band = source.omega_band_rad_s();
            if band.0 != 0.0 || band.1 < nyquist {
                return Err(invalid("refreshed relaxation sources must support DC through the original sampled Nyquist band"));
            }
            // Re-querying original constant claims must preserve every stored
            // coefficient. No tolerance can silently admit new kinetic/storage law.
            if source.elastic().layer_integrals() != section.layer_integrals() {
                return Err(invalid("rate-only refresh changed original reference mass or elastic storage"));
            }
            for (old, new) in member.source.spectra().iter().zip(source.spectra()) {
                if old.e_inf.to_bits() != new.e_inf.to_bits() || old.terms.len() != new.terms.len()
                    || old.terms.iter().zip(&new.terms).any(|(a,b)| a.0.to_bits() != b.0.to_bits()) {
                    return Err(invalid("rate-only refresh changed source modulus increments or memory topology"));
                }
            }
            refreshed.push(FrameRelaxationMember { element: member.element, source });
        }
        let mut updates = Vec::with_capacity(self.coordinates.len());
        for member in &refreshed {
            for (layer, law) in member.source.spectra().iter().enumerate() {
                for (arm, &(increment, tau)) in law.terms.iter().enumerate() {
                    poll(gate)?;
                    if increment == 0.0 { continue; }
                    for channel in [FrameStrainChannel::Axial, FrameStrainChannel::SymmetricBending,
                        FrameStrainChannel::AntisymmetricBending] {
                        let index = updates.len();
                        let coordinate = self.coordinates.get(index)
                            .ok_or_else(|| invalid("refreshed source added an unowned material coordinate"))?;
                        if coordinate.element != member.element || coordinate.layer != layer
                            || coordinate.source_arm != arm || coordinate.channel != channel
                            || coordinate.native_memory != self.memory.start + index {
                            return Err(invalid("refreshed source changed original material coordinate order"));
                        }
                        updates.push(MaterialRelaxationTime { branch: coordinate.native_memory,
                            relaxation_time_s: tau });
                    }
                }
            }
        }
        if updates.len() != self.coordinates.len() {
            return Err(invalid("refreshed source omitted an original material coordinate"));
        }
        poll(gate)?;
        let report = mechanical.set_material_relaxation_times(&updates, gate)?;
        // The native owner checked cancellation immediately before publication.
        // No fallible operation or new cancellation point may split these moves.
        for (branch, update) in self.branches.iter_mut().zip(updates) {
            branch.relaxation_time_s = update.relaxation_time_s;
        }
        self.members = refreshed;
        Ok(report)
    }
}

fn require_constant_storage(material: &ResolvedMaterialStatePoint,
    roles: &LayerRelaxationProperties) -> Result<()> {
    for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, roles.equilibrium_modulus.as_str()]
        .into_iter().chain(roles.branches.iter().map(|pair| pair.0.as_str())) {
        let property = material.property(name)
            .ok_or_else(|| invalid("temperature refresh lacks an original storage coefficient"))?;
        if !property.source_axes().is_empty()
            || property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity {
            return Err(invalid("rate-only temperature feedback requires source-declared constant density and storage moduli"));
        }
    }
    Ok(())
}

fn temperature_query(material: &ResolvedMaterialStatePoint, axis: &str, t: f64) -> Result<QueryPoint> {
    let absolute = QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static));
    let quantity = material.axis_quantities().get(axis).copied();
    if !t.is_finite() || t <= 0.0 || !material.query_point().iter().any(|(name,_)| name == axis)
        || !(quantity == Some(absolute) || axis == "T"
            && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS))) {
        return Err(invalid("existing absolute-temperature coordinate and finite positive Kelvin required"));
    }
    if material.query_point().iter().any(|(name,_)| name != axis && (name == "T"
        || material.axis_quantities().get(name).is_some_and(|q| q.dims() == Temperature::DIMS))) {
        return Err(invalid("a layer rate update cannot leave a second temperature coordinate unowned"));
    }
    material.query().with(axis, t).map_err(|_| invalid("original temperature query could not be retained"))
}

#[cfg(test)]
mod tests;
