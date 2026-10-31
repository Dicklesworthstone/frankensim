//! Source-bound axial storage on the original articulated/thermal timeline.
//!
//! A perfectly bonded section is condensed onto one prismatic coordinate.
//! The original section owner supplies both force and potential; no modulus,
//! expansion integral, mass or time integrator is reconstructed here. Source
//! coefficients freeze during a mechanical substep. Their endpoint change is
//! explicit parameter work, not reciprocal thermoelastic heat.
use super::{Error, finite, positive};
use crate::thermal::{
    ExpandedArticulated,
    heating::{HeatingError, ThermalRegion},
};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section_temperature::thermoelastic::ConcentricThermoelasticLaw;
use fs_material::state_point::{
    DENSITY_PROPERTY, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
};
use fs_mbd::articulated::JointKind;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeSet;

/// Finite axial compliance between a prismatic child and its direct parent.
///
/// The member's reference material inventory is ALREADY carried by these links.
/// One existing complete thermal region represents each concentric layer.
/// Its reference volume must be area * reference length and its actual mass
/// must match the original section density. No extra member mass is added.
/// This is lumped-mass compliance on the existing joint. Rigid-link inertia
/// and prescribed placement retain their own meaning; distributed deformation,
/// longitudinal rod modes and FEM are outside this small-strain reduction.
#[derive(Clone, Debug)]
pub struct AxialJointMember<'a> {
    pub link: usize,
    /// Existing finite cells, in spatial section-layer order, without repeats.
    pub cells: Vec<CellId>,
    /// Original source-pinned section. Scalar reference rho/alpha and scalar or
    /// globally affine E(T) follow this owner's existing constitutive scope.
    pub law: ConcentricThermoelasticLaw<'a>,
    pub reference_length_m: f64,
    /// q at zero total engineering strain, at the source's stress-free T.
    /// Positive q increases member length with unit displacement Jacobian.
    pub reference_joint_position_m: f64,
    /// Bound on total, thermal and each layer's elastic strain, strictly (0,1).
    pub maximum_absolute_strain: f64,
}

/// Actual accepted endpoint, including original layer source and phase answers.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialMemberState {
    pub link: usize,
    pub cells: Vec<CellId>,
    pub phases: Vec<EquilibriumPhaseState>,
    pub section: ConcentricThermoelasticSection,
    pub total_axial_strain: f64,
    /// Generalized force on the child; the native joint owns parent reaction.
    pub restoring_effort_n: f64,
    /// Frozen-temperature mechanical potential, including self-equilibrated
    /// layer stress. Not the full thermoelastic internal energy psi + T*s.
    pub elastic_energy_j: f64,
}

#[derive(Clone)]
pub(super) struct BoundAxial<'a> {
    definition: AxialJointMember<'a>,
    pub dof: usize,
    pub state: AxialMemberState,
    regions: Vec<ThermalRegion<'a>>,
    /// Original body and reference-region indices, in section-layer order.
    locations: Vec<(usize, usize)>,
    masses: Vec<f64>,
}
impl<'a> BoundAxial<'a> {
    pub fn new(
        definition: AxialJointMember<'a>,
        configuration: &ExpandedArticulated<'_>,
        regions: &[Vec<ThermalRegion<'a>>],
        thermal: &DepositionSimulation<'_>,
        position: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if !positive(definition.reference_length_m)
            || !definition.reference_joint_position_m.is_finite()
            || !positive(definition.maximum_absolute_strain)
            || definition.maximum_absolute_strain >= 1.0
            || definition.cells.len() != definition.law.reference().layers().len()
            || definition.cells.iter().collect::<BTreeSet<_>>().len() != definition.cells.len()
        {
            return Err(Error::Invalid(
                "axial member needs finite reference geometry, bounded strain and one distinct cell per layer",
            ));
        }
        let model = configuration.model();
        let link = model
            .links()
            .get(definition.link)
            .ok_or(Error::Invalid("axial member child link does not exist"))?;
        let parent = link
            .parent()
            .ok_or(Error::Invalid("axial member needs a parent-child joint"))?;
        if link.joint().kind() != JointKind::Prismatic {
            return Err(Error::Invalid("axial member requires a prismatic joint"));
        }
        let dof = model
            .dof_index(definition.link)
            .ok_or(Error::Invalid("axial member joint has no coordinate"))?;
        let mut selected = Vec::with_capacity(definition.cells.len());
        let mut locations = Vec::with_capacity(definition.cells.len());
        let mut masses = Vec::with_capacity(definition.cells.len());
        for (i, cell_id) in definition.cells.iter().enumerate() {
            checkpoint()?;
            let (body, region) = [definition.link, parent]
                .into_iter()
                .find_map(|body| {
                    regions[body]
                        .iter()
                        .find(|r| r.cell_id == *cell_id)
                        .map(|r| (body, r))
                })
                .ok_or(Error::Invalid(
                    "axial layers must be actual regions of their child or direct parent",
                ))?;
            let original = configuration.bodies()[body].reference();
            let index = original
                .regions()
                .iter()
                .position(|r| r.key() == region.key)
                .ok_or(Error::Invalid("axial member lost its original region"))?;
            let material = original.regions()[index].material();
            let layer = &definition.law.reference().layers()[i].material;
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
                    "axial layer must retain its original material, state context and density source",
                ));
            }
            let expanded = configuration.bodies()[body]
                .regions()
                .iter()
                .find(|r| r.reference_index() == index)
                .ok_or(Error::Invalid("axial member lost its expanded region"))?;
            if layer.property(ALPHA).map(|p| p.answer().receipt.selected)
                != expanded
                    .expansion()
                    .reference()
                    .resolved()
                    .property(ALPHA)
                    .map(|p| p.answer().receipt.selected)
            {
                return Err(Error::Invalid(
                    "axial eigenstrain and link geometry must share their original expansion source",
                ));
            }
            let integral = &definition.law.reference().layer_integrals()[i];
            let volume = finite(integral.area_m2() * definition.reference_length_m)?;
            let mass = finite(integral.line_density_kg_m() * definition.reference_length_m)?;
            let cell = thermal
                .cells()
                .get(cell_id)
                .ok_or(Error::Invalid("axial member thermal cell is inactive"))?;
            if !same_inventory(volume, region.geometry.volume_m3)
                || !same_inventory(mass, cell.mass_kg)
            {
                return Err(Error::Invalid(
                    "axial section must represent the existing reference volume and mass exactly once",
                ));
            }
            masses.push(cell.mass_kg);
            locations.push((body, index));
            selected.push(region.clone());
        }
        let phases = phases(&selected, &masses, thermal, checkpoint)?;
        let strain = strain(&definition, position[dof])?;
        let temperatures = phases
            .iter()
            .map(|phase| phase.temperature_k())
            .collect::<Vec<_>>();
        let section = definition
            .law
            .at_temperatures_k(&temperatures, strain, definition.maximum_absolute_strain)?
            .section;
        let state = AxialMemberState {
            link: definition.link,
            cells: definition.cells.clone(),
            phases,
            section,
            total_axial_strain: strain,
            restoring_effort_n: 0.0,
            elastic_energy_j: 0.0,
        };
        let mut bound = Self {
            definition,
            dof,
            state,
            regions: selected,
            locations,
            masses,
        };
        bound.validate_expansion(configuration)?;
        bound.update_position(position[dof])?;
        checkpoint()?;
        Ok(bound)
    }

    /// The original force and its exact same frozen-temperature potential.
    /// No source query, hashing, allocation or separate stiffness fit occurs.
    pub fn response(&self, position: f64) -> Result<(f64, f64), Error> {
        let strain = strain(&self.definition, position)?;
        let (force, energy_per_length) = self
            .state
            .section
            .axial_force_and_energy_at_strain(strain, self.definition.maximum_absolute_strain)?;
        let energy = finite(energy_per_length * self.definition.reference_length_m)?;
        if energy_per_length != 0.0 && energy == 0.0 {
            return Err(Error::Invalid("axial stored energy underflow"));
        }
        Ok((-force, energy))
    }

    pub fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        configuration: &ExpandedArticulated<'_>,
        position: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let phases = phases(&self.regions, &self.masses, thermal, checkpoint)?;
        let total_strain = strain(&self.definition, position[self.dof])?;
        let temperatures = phases
            .iter()
            .map(|phase| phase.temperature_k())
            .collect::<Vec<_>>();
        let section = self
            .definition
            .law
            .at_temperatures_k(
                &temperatures,
                total_strain,
                self.definition.maximum_absolute_strain,
            )?
            .section;
        checkpoint()?;
        self.state.section = section;
        self.state.phases = phases;
        self.validate_expansion(configuration)?;
        self.update_position(position[self.dof])?;
        Ok(())
    }

    fn update_position(&mut self, position: f64) -> Result<(), Error> {
        let (force, energy) = self.response(position)?;
        self.state.total_axial_strain = strain(&self.definition, position)?;
        self.state.restoring_effort_n = force;
        self.state.elastic_energy_j = energy;
        Ok(())
    }

    fn validate_expansion(&self, configuration: &ExpandedArticulated<'_>) -> Result<(), Error> {
        for (&(body, index), path) in self
            .locations
            .iter()
            .zip(self.state.section.expansion_paths())
        {
            let actual = configuration.bodies()[body]
                .regions()
                .iter()
                .find(|r| r.reference_index() == index)
                .ok_or(Error::Invalid("axial member lost its source geometry"))?;
            if path.free_linear_strain() != actual.expansion().free_linear_strain() {
                return Err(Error::Invalid(
                    "axial and geometric expansion must use the same actual temperature path",
                ));
            }
        }
        Ok(())
    }
}

fn strain(definition: &AxialJointMember<'_>, position: f64) -> Result<f64, Error> {
    let displacement = finite(position - definition.reference_joint_position_m)?;
    let strain = finite(displacement / definition.reference_length_m)?;
    if displacement != 0.0 && strain == 0.0 {
        return Err(Error::Invalid("axial strain underflow"));
    }
    Ok(strain)
}

fn same_inventory(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 128.0 * f64::EPSILON * a.max(b)
}

fn phases(
    regions: &[ThermalRegion<'_>],
    masses: &[f64],
    thermal: &DepositionSimulation<'_>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<EquilibriumPhaseState>, Error> {
    let mut phases = Vec::with_capacity(regions.len());
    for (region, &mass) in regions.iter().zip(masses) {
        checkpoint()?;
        let cell = thermal
            .cells()
            .get(&region.cell_id)
            .ok_or(Error::Invalid("axial member thermal cell is inactive"))?;
        if cell.material_id != region.material_id
            || cell.body_id != region.body_id
            || cell.segment_id != region.segment_id
            || cell.mass_kg != mass
            || cell.geometry != region.geometry
            || !thermal.uses_material_law(cell.material_id, region.caloric.curve())
        {
            return Err(Error::Invalid(
                "axial member must retain its exact finite cell and caloric inventory",
            ));
        }
        let phase = region
            .caloric
            .curve()
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(HeatingError::from)?;
        if !phase.is_fully_solid() || phase.material_card_identity() != region.card.content_hash() {
            return Err(Error::Invalid(
                "axial member needs its original solid caloric state",
            ));
        }
        phases.push(phase);
    }
    Ok(phases)
}
