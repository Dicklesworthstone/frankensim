//! Source-bound axial storage on the original articulated/thermal timeline.
//!
//! A perfectly bonded section is condensed onto one prismatic coordinate.
//! The original section owner supplies both force and potential; no modulus,
//! expansion integral, mass or time integrator is reconstructed here. Source
//! coefficients freeze during a mechanical substep. Their endpoint change is
//! explicit parameter work, not reciprocal thermoelastic heat.
use super::{Error, finite, member::MemberInventory, positive};
use crate::thermal::{ExpandedArticulated, heating::ThermalRegion};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section_temperature::thermoelastic::ConcentricThermoelasticLaw;
use fs_mbd::articulated::JointKind;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

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
    inventory: MemberInventory<'a>,
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
        {
            return Err(Error::Invalid(
                "axial member needs finite reference geometry, bounded strain and one distinct cell per layer",
            ));
        }
        let inventory = MemberInventory::new(
            definition.link,
            JointKind::Prismatic,
            &definition.cells,
            definition.law.reference(),
            definition.reference_length_m,
            configuration,
            regions,
            thermal,
            checkpoint,
        )?;
        let dof = inventory.dof;
        let phases = inventory.phases(thermal, checkpoint)?;
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
            inventory,
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
        let phases = self.inventory.phases(thermal, checkpoint)?;
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
            .inventory
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

