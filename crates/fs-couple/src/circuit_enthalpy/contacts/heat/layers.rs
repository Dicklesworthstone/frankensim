//! Source-integrated k(T) fluxes on the original composed heat network.
//! The material owner retains the interface solve and source charts; this
//! adapter only checks correspondence and replaces the selected native links.

use super::*;
use fs_material::thermal::conduction::{
    LayerConductionConfig,
    network::{LayerLinkSample, PreparedLayerConduction},
};
use std::collections::BTreeMap;

pub(super) struct IntegratedConduction<'n, 'a> {
    model: &'n PreparedLayerConduction<'n, 'a>,
    controls: LayerConductionConfig,
}

impl<'n, 'a> PreparedHeatNetwork<'n, 'a> {
    /// Replace selected endpoint-k approximations by the existing source-curve
    /// layer flux solve, in EVERY thermal trial and accepted observer epoch.
    ///
    /// The layer bindings must retain this assembly's original cells, caloric
    /// objects, conductivity keys/claims, full fixed conditions and distances.
    /// Their explicit constant contact resistance belongs to the layer solve;
    /// the same link cannot also carry a separate material thermal contact.
    /// Untargeted conduction, electrical interfaces and radiation are unchanged.
    ///
    /// The original owner admits homogeneous isotropic source-free layers in
    /// one pure phase with fixed geometry. It integrates k(T) through the actual
    /// interface temperatures. This does not add spatial storage within those
    /// layers, arbitrary tensor heat flow, or phase-boundary propagation.
    /// Source-domain exits and cancellation refuse the whole coupled candidate.
    pub fn with_integrated_conduction(
        mut self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        model: &'n PreparedLayerConduction<'n, 'a>,
        controls: LayerConductionConfig,
    ) -> Result<Self, HeatNetworkError> {
        check(cx)?;
        if self.integrated_conduction.is_some() || model.bindings().len() == 0 {
            return Err(HeatNetworkError::Invalid(
                "one nonempty integrated conduction binding is required",
            ));
        }
        let initial = self.resolve(cx, thermal)?;
        self.admit_integrated_conduction(cx, thermal, model, &initial.bulk.source)?;
        self.integrated_conduction = Some(IntegratedConduction { model, controls });
        self.resolve(cx, thermal)?;
        check(cx)?;
        Ok(self)
    }

    // Immutable material/geometry correspondence is checked once. The original
    // owners retain source pins, and the layer network checks live caloric and
    // contact geometry at each trial without rebuilding this correspondence.
    fn admit_integrated_conduction(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        model: &PreparedLayerConduction<'_, '_>,
        sample: &HeatNetworkSample,
    ) -> Result<(), HeatNetworkError> {
        check(cx)?;
        let invalid = HeatNetworkError::Invalid;
        let bulk = self.network.electrical.bulk;
        for binding in model.bindings() {
            check(cx)?;
            if self
                .thermal_contacts
                .iter()
                .any(|contact| contact.link == binding.link)
            {
                return Err(invalid(
                    "integrated layer contact resistance cannot share a separate thermal-contact owner",
                ));
            }
            let original = bulk
                .links
                .iter()
                .find(|link| link.link == binding.link)
                .ok_or(invalid(
                    "integrated conduction link is absent from the material assembly",
                ))?;
            if binding.distances_m != [original.distance_a_m, original.distance_b_m] {
                return Err(invalid(
                    "integrated conduction must retain original layer distances",
                ));
            }
            let link = thermal.links().get(&binding.link).ok_or(invalid(
                "integrated conduction needs the original active thermal link",
            ))?;
            for (side, id) in [link.a, link.b].into_iter().enumerate() {
                let cell = bulk
                    .cells
                    .iter()
                    .find(|cell| cell.cell == id)
                    .ok_or(invalid("integrated layer has no original material cell"))?;
                let source = sample
                    .material
                    .cells
                    .iter()
                    .find(|cell| cell.cell == id)
                    .ok_or(invalid("integrated layer has no original material sample"))?;
                let layer = binding.layers[side];
                let line = layer.source_line();
                if !std::ptr::eq(layer.caloric_curve(), cell.curve)
                    || line.axis != cell.point.temperature_axis
                {
                    return Err(invalid(
                        "integrated layer must retain the same caloric object and temperature axis",
                    ));
                }
                let (key, _) = cell.conductivity.as_ref().ok_or(invalid(
                    "integrated layer requires the original conductivity binding",
                ))?;
                let property = source
                    .conductivity
                    .as_ref()
                    .and_then(|material| material.property(key.name()))
                    .ok_or(invalid(
                        "integrated layer lost its original conductivity source",
                    ))?;
                let (t, answer) = line
                    .knots
                    .first()
                    .ok_or(invalid("integrated layer has no retained source interval"))?;
                let point = cell
                    .point
                    .at(source.phase)
                    .map_err(ElectricalContactNetworkError::State)?
                    .with(cell.point.temperature_axis, *t)
                    .map_err(|_| invalid("integrated layer source context is invalid"))?;
                if answer.receipt.property != key.name()
                    || answer.receipt.selected != property.answer().receipt.selected
                    || answer.evidence.value.quantity != key.quantity()
                    || answer.receipt.axis_quantities != *point.axis_quantities()
                    || answer.receipt.query_point.len() != point.axes().len()
                    || answer.receipt.query_point.iter().any(|(axis, value)| {
                        point
                            .axes()
                            .get(axis)
                            .is_none_or(|x| x.to_bits() != value.to_bits())
                    })
                {
                    return Err(invalid(
                        "integrated and bulk conductivity must share the same property, source and full fixed conditions",
                    ));
                }
            }
        }
        check(cx)?;
        Ok(())
    }

    pub(super) fn resolve_integrated_conduction(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        conductances: &mut BTreeMap<u64, f64>,
        sample: &mut MaterialNetworkSample,
    ) -> Result<Vec<LayerLinkSample>, HeatNetworkError> {
        let Some(integrated) = &self.integrated_conduction else {
            return Ok(Vec::new());
        };
        check(cx)?;
        let invalid = HeatNetworkError::Invalid;
        let evaluated = integrated.model.resolve(thermal, integrated.controls, || {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)
        })?;
        for (&id, &value) in &evaluated.conductances_w_k {
            check(cx)?;
            *conductances.get_mut(&id).ok_or(invalid(
                "integrated conductivity lost its native link coefficient",
            ))? = value;
            // Endpoint k samples remain truthful, while the reported applied
            // secant now describes the original source-integrated interface.
            let index = sample
                .links
                .binary_search_by_key(&id, |link| link.link)
                .map_err(|_| invalid("integrated conductivity lost its native link sample"))?;
            sample.links[index].conductance_w_k = value;
        }
        check(cx)?;
        Ok(evaluated.payload)
    }
}
