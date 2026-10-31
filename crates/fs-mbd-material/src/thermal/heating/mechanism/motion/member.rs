//! One finite material inventory shared by axial and torsional joint reductions.
use super::{Error, positive};
use crate::thermal::{
    ExpandedArticulated,
    heating::{HeatingError, ThermalRegion},
};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::section::ConcentricElasticSection;
use fs_material::state_point::{
    DENSITY_PROPERTY, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
};
use fs_mbd::articulated::JointKind;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeSet;

#[derive(Clone)]
pub(super) struct MemberInventory<'a> {
    pub dof: usize,
    pub regions: Vec<ThermalRegion<'a>>,
    pub locations: Vec<(usize, usize)>,
    masses: Vec<f64>,
}
impl<'a> MemberInventory<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        link_index: usize,
        kind: JointKind,
        cells: &[CellId],
        section: &ConcentricElasticSection,
        length_m: f64,
        configuration: &ExpandedArticulated<'_>,
        regions: &[Vec<ThermalRegion<'a>>],
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if !positive(length_m)
            || cells.len() != section.layers().len()
            || cells.iter().collect::<BTreeSet<_>>().len() != cells.len()
        {
            return Err(Error::Invalid(
                "member needs positive reference length and one distinct finite cell per layer",
            ));
        }
        let model = configuration.model();
        let link = model
            .links()
            .get(link_index)
            .ok_or(Error::Invalid("member child link does not exist"))?;
        let parent = link
            .parent()
            .ok_or(Error::Invalid("member needs a parent-child joint"))?;
        if link.joint().kind() != kind {
            return Err(Error::Invalid(
                "member joint kind differs from its physical coordinate",
            ));
        }
        let dof = model
            .dof_index(link_index)
            .ok_or(Error::Invalid("member joint has no coordinate"))?;
        let mut selected = Vec::with_capacity(cells.len());
        let mut locations = Vec::with_capacity(cells.len());
        let mut masses = Vec::with_capacity(cells.len());
        for (i, cell_id) in cells.iter().enumerate() {
            checkpoint()?;
            let (body, region) = [link_index, parent]
                .into_iter()
                .find_map(|body| {
                    regions[body]
                        .iter()
                        .find(|r| r.cell_id == *cell_id)
                        .map(|r| (body, r))
                })
                .ok_or(Error::Invalid(
                    "member layers must be actual regions of their child or direct parent",
                ))?;
            let original = configuration.bodies()[body].reference();
            let index = original
                .regions()
                .iter()
                .position(|r| r.key() == region.key)
                .ok_or(Error::Invalid("member lost its original region"))?;
            let material = original.regions()[index].material();
            let layer = &section.layers()[i].material;
            if layer.card_identity() != region.card.content_hash()
                || layer.query_point() != material.query_point()
                || layer.axis_quantities() != material.axis_quantities()
                || layer
                    .property(DENSITY_PROPERTY)
                    .map(|p| p.answer().receipt.selected)
                    != material
                        .property(DENSITY_PROPERTY)
                        .map(|p| p.answer().receipt.selected)
            {
                return Err(Error::Invalid(
                    "member layer must retain its original card, context and density source",
                ));
            }
            let expanded = configuration.bodies()[body]
                .regions()
                .iter()
                .find(|r| r.reference_index() == index)
                .ok_or(Error::Invalid("member lost its expanded region"))?;
            if layer.property(ALPHA).map(|p| p.answer().receipt.selected)
                != expanded
                    .expansion()
                    .reference()
                    .resolved()
                    .property(ALPHA)
                    .map(|p| p.answer().receipt.selected)
            {
                return Err(Error::Invalid(
                    "member and carrier geometry must share their original expansion source",
                ));
            }
            let integral = &section.layer_integrals()[i];
            let volume = integral.area_m2() * length_m;
            let mass = integral.line_density_kg_m() * length_m;
            let cell = thermal
                .cells()
                .get(cell_id)
                .ok_or(Error::Invalid("member thermal cell is inactive"))?;
            if !same_inventory(volume, region.geometry.volume_m3)
                || !same_inventory(mass, cell.mass_kg)
            {
                return Err(Error::Invalid(
                    "member section must represent its existing reference volume and mass exactly once",
                ));
            }
            masses.push(cell.mass_kg);
            locations.push((body, index));
            selected.push(region.clone());
        }
        Ok(Self {
            dof,
            regions: selected,
            locations,
            masses,
        })
    }

    pub fn phases(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<EquilibriumPhaseState>, Error> {
        let mut phases = Vec::with_capacity(self.regions.len());
        for (region, &mass) in self.regions.iter().zip(&self.masses) {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&region.cell_id)
                .ok_or(Error::Invalid("member thermal cell is inactive"))?;
            if cell.material_id != region.material_id
                || cell.body_id != region.body_id
                || cell.segment_id != region.segment_id
                || cell.mass_kg != mass
                || cell.geometry != region.geometry
                || !thermal.uses_material_law(cell.material_id, region.caloric.curve())
            {
                return Err(Error::Invalid(
                    "member must retain its exact finite cell and caloric inventory",
                ));
            }
            let phase = region
                .caloric
                .curve()
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(HeatingError::from)?;
            if !phase.is_fully_solid()
                || phase.material_card_identity() != region.card.content_hash()
            {
                return Err(Error::Invalid(
                    "member needs its original solid caloric state",
                ));
            }
            phases.push(phase);
        }
        Ok(phases)
    }
}

fn same_inventory(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 128.0 * f64::EPSILON * a.max(b)
}
