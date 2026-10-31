//! Finite-body receiving cells for original Maxwell resistance channels.
//!
//! The source geometry supplies the receiving material volume. All heat of
//! that body enters its explicit uniform cell; there is no fitted partition.
//! Existing conduction can carry it into a contact film in the SAME nonlinear
//! thermal/contact trial. Spatially nonuniform body temperature needs a finer
//! material projection; a modal energy fraction is not such a projection.
//! The existing source integral admits held-effort modal/Maxwell trajectories;
//! external radiation states and varying-force integrals remain unsupported.

use super::*;
use crate::render::schedule::force::coupled::coefficients::heat::{MaterialHeatIntegral, MaterialHeatBudget};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::Temperature;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// One original material body and its independently integrated branch losses.
/// The geometry producer must bind the complete branch range for this volume.
/// Coefficients remain source-declared constants over the receiving cell's
/// temperature domain; a temperature-dependent stored modulus requires its
/// own thermoelastic history continuation and is refused here.
pub struct MaterialRelaxationHeatDestination<'a> {
    pub cell: CellId,
    pub material: &'a MaterialCard,
    pub caloric: &'a EquilibriumEnthalpyPhaseCurve,
    pub temperature_axis: &'a str,
    pub density_property: &'a str,
    pub volume_m3: f64,
    pub source: ResolvedMaterialStatePoint,
    pub channel: MaterialHeatIntegral,
}

impl MaterialRelaxationHeatDestination<'_> {
    /// Check this source body's caloric law, material domain, physical volume
    /// and receiving mass without advancing either owner. Description compilers
    /// may use this before output; execution repeats admission at every trial.
    pub fn admit<E>(&self, thermal: &DepositionSimulation<'_>) -> Result<(), LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        if self.channel.groups().len() != 1 || self.channel.groups()[0].is_empty() {
            return Err(bad("one complete nonempty material group must own each finite body volume"));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static));
        let temperature_quantity = self.source.axis_quantities().get(self.temperature_axis);
        let owns_temperature = match temperature_quantity {
            Some(q) => *q == absolute || self.temperature_axis == "T"
                && *q == QuantitySpec::dimensional(Temperature::DIMS),
            None => self.temperature_axis == "T",
        };
        if !owns_temperature || self.source.state_coordinate(self.temperature_axis).is_none() {
            return Err(bad("material heat requires its original absolute-temperature coordinate"));
        }
        let cell = thermal.cells().get(&self.cell).ok_or(bad("material heat receiving cell is inactive"))?;
        if !thermal.uses_material_law(cell.material_id, self.caloric)
            || self.source.card_identity() != self.material.content_hash() {
            return Err(bad("material heat must retain the original source card and caloric object"));
        }
        let phase = self.caloric.state_at_specific_enthalpy(cell.specific_h_j_kg)?;
        if !phase.is_fully_solid() || phase.material_card_identity() != self.material.content_hash() {
            return Err(bad("hereditary material heat requires the same fully solid source material"));
        }
        let density = self.source.property(self.density_property)
            .ok_or(bad("material heat requires its original source density"))?;
        if density.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS)
            || !same_size(cell.geometry.volume_m3, self.volume_m3)
            || !same_size(cell.mass_kg, self.volume_m3 * density.value_si())
            || !same_size(phase.bulk_density_kg_m3(), density.value_si())
            || self.caloric.knots().iter().any(|k| !same_size(k.bulk_density_kg_m3, density.value_si())) {
            return Err(bad("material heat cell must contain its actual source volume and fixed reference mass"));
        }
        for property in self.source.properties() {
            let claim = self.material.claims().claim(property.answer().receipt.selected)
                .ok_or(bad("material heat source claim is absent from its original card"))?;
            if claim.value.axes().any(|(_, dims)| dims == Temperature::DIMS) {
                return Err(bad("temperature-dependent body coefficients need their original thermoelastic continuation"));
            }
        }
        let knots = self.caloric.knots();
        for temperature in [phase.temperature_k(), knots.first().expect("caloric knots").temperature_k,
            knots.last().expect("caloric knots").temperature_k] {
            let query = self.source.query().with(self.temperature_axis, temperature)
                .map_err(|_| bad("material heat temperature must retain its original source schema"))?;
            let refreshed = self.source.at_query_point(self.material, &query)?;
            if refreshed.properties().iter().zip(self.source.properties()).any(|(a, b)|
                a.value_si().to_bits() != b.value_si().to_bits()) {
                return Err(bad("material heat changed a fixed original constitutive coefficient"));
            }
        }
        Ok(())
    }
}

impl MultiFrictionModalSystem {
    /// Reuse the original material-group heat owner and its aggregate resource
    /// budget. Prepare one complete group per finite body destination and bind
    /// it to the original source geometry/material. External radiation states
    /// retain that owner's explicit refusal; none is replaced with an endpoint
    /// force or silently excluded from the underlying trajectory.
    pub fn prepare_material_heat(
        &self, branches: &[core::ops::Range<usize>], budget: MaterialHeatBudget, gate: &CancelGate,
    ) -> Result<MaterialHeatIntegral, ModalFrictionError> {
        Ok(self.core.network.prepare_material_heat(branches, budget, gate)?)
    }
}
