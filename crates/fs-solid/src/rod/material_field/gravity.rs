//! Reference-mass gravity for the existing spatial rod equilibrium.
//!
//! For a segment with constant reference line density, integrating the linear
//! centreline shape functions gives half its mass to each incident node. This
//! preserves both total weight and its first moment, including at the clamped
//! root. The load is `m_node * g`, not current density times deformed length.
//! It is a consistent static body-load vector, NOT a dynamic mass matrix.

use fs_exec::CancelGate;
use fs_material::section::ConcentricShearSection;

use super::{SpatialMaterialRod, SpatialRodError, checkpoint};
use super::supports::RodSupportReaction;
use crate::rod::Rod;
use crate::rod::derivatives::RodEnergy;

mod buoyancy;
pub use buoyancy::{ReferenceDisplacementVolume, SubmergedFluidLoad};

/// Immutable reference inventory shared by thermal snapshots of one specimen.
/// Original sections retain the actual density receipts and reference geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct RodReferenceMass {
    sections: Vec<ConcentricShearSection>,
    segment_length_m: f64,
    segment_mass_kg: Vec<f64>,
    nodal_mass_kg: Vec<f64>,
    total_mass_kg: f64,
}

impl RodReferenceMass {
    pub(super) fn try_new(
        segment_length_m: f64, sections: &[ConcentricShearSection],
    ) -> Result<Self, SpatialRodError> {
        let mut segment_mass_kg = Vec::with_capacity(sections.len());
        let mut nodal_mass_kg = vec![0.0; sections.len() + 1];
        let mut total_mass_kg = 0.0;
        for (segment, section) in sections.iter().enumerate() {
            let mass = segment_length_m * section.elastic().line_density_kg_m();
            let half = 0.5 * mass;
            if !mass.is_finite() || mass <= 0.0 || half == 0.0 {
                return Err(SpatialRodError::Invalid("reference segment mass is unrepresentable"));
            }
            segment_mass_kg.push(mass);
            nodal_mass_kg[segment] += half;
            nodal_mass_kg[segment + 1] += half;
            total_mass_kg += mass;
        }
        if !total_mass_kg.is_finite() || !nodal_mass_kg.iter().all(|m| m.is_finite()) {
            return Err(SpatialRodError::Invalid("reference rod mass sum is unrepresentable"));
        }
        Ok(Self { sections: sections.to_vec(), segment_length_m, segment_mass_kg,
            nodal_mass_kg, total_mass_kg })
    }
    /// Original sections at the mass-reference state, NOT current heated data.
    #[must_use]
    pub fn sections(&self) -> &[ConcentricShearSection] { &self.sections }
    /// Reference segment length [m].
    #[must_use]
    pub const fn segment_length_m(&self) -> f64 { self.segment_length_m }
    /// Material mass per segment [kg], unchanged by temperature or elastic motion.
    #[must_use]
    pub fn segment_mass_kg(&self) -> &[f64] { &self.segment_mass_kg }
    /// Integrated centreline shape-function weights [kg], including root zero.
    #[must_use]
    pub fn nodal_mass_kg(&self) -> &[f64] { &self.nodal_mass_kg }
    /// Total reference mass [kg].
    #[must_use]
    pub const fn total_mass_kg(&self) -> f64 { self.total_mass_kg }
}

#[derive(Clone, Debug)]
pub(super) struct GravityLoad {
    acceleration_m_s2: [f64; 3],
    force_n: Vec<[f64; 3]>,
    buoyancy: Option<SubmergedFluidLoad>,
}

impl SpatialMaterialRod {
    /// Reference inventory, retained even when current material density changes.
    #[must_use]
    pub fn reference_mass(&self) -> &RodReferenceMass { &self.reference_mass }

    /// Set an explicit uniform world-frame gravity vector on a new static
    /// specimen. Zero removes its forces; there is no implicit terrestrial preset.
    /// This replaces the old gravity field, never accumulates another weight.
    /// Reference mass and every support are preserved. An explicitly bound
    /// submerged fluid is retained and its buoyancy follows the new gravity.
    /// Without that binding no buoyancy is inferred.
    ///
    /// Gravity is a permanent dead load in every equilibrium stage; `steps` in
    /// `solve_static` continues to ramp the separately supplied TIP load only.
    /// This is not dynamic activation or a work-conserving change of environment.
    ///
    /// # Errors
    /// Nonfinite acceleration, force underflow/overflow, potential overflow, or
    /// cancellation refuses without changing the original specimen.
    pub fn with_gravity(
        &self, acceleration_m_s2: [f64; 3], gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        if !acceleration_m_s2.iter().all(|v| v.is_finite()) {
            return Err(SpatialRodError::Invalid("gravity must be finite world-frame acceleration in m/s2"));
        }
        let fluid = self.submerged_fluid_load()
            .map(|load| (load.transport_state(), load.volume_model()));
        let mut candidate = self.clone();
        candidate.gravity = if acceleration_m_s2.iter().all(|v| *v == 0.0) && fluid.is_none() {
            None
        } else {
            let mut masses = self.reference_mass.nodal_mass_kg().to_vec();
            for link in &self.tension_links {
                checkpoint(gate)?;
                for attachment in link.attachments() {
                    masses[attachment.node] += 0.5 * link.reference_mass_kg();
                }
            }
            let mut force_n = Vec::with_capacity(self.rod.positions.len());
            for mass in masses {
                checkpoint(gate)?;
                let force = acceleration_m_s2.map(|g| mass * g);
                if force.iter().zip(acceleration_m_s2).any(|(&f, g)|
                    !f.is_finite() || (g != 0.0 && f == 0.0))
                {
                    return Err(SpatialRodError::Invalid("reference-mass gravity force is unrepresentable"));
                }
                force_n.push(force);
            }
            Some(GravityLoad { acceleration_m_s2, force_n, buoyancy: None })
        };
        if let Some((transport, volume_model)) = fluid {
            let load = SubmergedFluidLoad::try_new(&candidate, transport, volume_model, gate)?;
            candidate.gravity.as_mut().expect("retained fluid keeps its gravity field").buoyancy = Some(load);
        }
        candidate.equilibrium_energy_for(&candidate.rod)?;
        checkpoint(gate)?;
        Ok(candidate)
    }

    /// Explicit gravity vector [m/s2]; absence means no gravitational loading.
    #[must_use]
    pub fn gravity_m_s2(&self) -> Option<[f64; 3]> {
        self.gravity.as_ref().map(|load| load.acceleration_m_s2)
    }
    /// World-frame body forces [N], including root and attached-link inventory.
    /// Rotating link offsets also contribute `gravity_nodal_moments_body_nm`.
    #[must_use]
    pub fn gravity_nodal_forces_n(&self) -> Option<&[[f64; 3]]> {
        self.gravity.as_ref().map(|load| load.force_n.as_slice())
    }
    /// Gravity potential change [J] relative to the actual unloaded reference pose.
    /// Kept separate from `energy_j()`, which remains stored elastic energy.
    pub fn gravity_potential_j(&self) -> Result<f64, SpatialRodError> {
        self.gravity_potential_for(&self.rod)
    }
    pub(super) fn gravity_potential_for(&self, rod: &Rod) -> Result<f64, SpatialRodError> {
        let Some(load) = &self.gravity else { return Ok(0.0); };
        let mut potential = 0.0;
        for (node, (position, force)) in rod.positions.iter().zip(&load.force_n).enumerate() {
            let reference = self.reference_geometry.as_ref().map_or(
                [node as f64 * rod.l0, 0.0, 0.0], |geometry| geometry.positions_m()[node]);
            for c in 0..3 { potential -= force[c] * (position[c] - reference[c]); }
        }
        // Nodal forces already include link translations; only rotating offsets
        // remain. Adding complete endpoint potentials would count weight twice.
        if !self.tension_links.is_empty() {
            potential += self.link_gravity_offset_potential(rod, load.acceleration_m_s2)?;
        }
        if !potential.is_finite() {
            return Err(SpatialRodError::Invalid("gravity potential is unrepresentable"));
        }
        Ok(potential)
    }

    // Same objective in every Newton probe, line search and reaction recovery.
    pub(super) fn equilibrium_energy_for(&self, rod: &Rod) -> Result<f64, SpatialRodError> {
        let internal = self.energy_for(rod)?;
        if self.gravity.is_none() { return Ok(internal); }
        let value = internal + self.gravity_potential_for(rod)? + self.buoyancy_potential_for(rod)?;
        if !value.is_finite() {
            return Err(SpatialRodError::Invalid("elastic plus gravity potential is unrepresentable"));
        }
        Ok(value)
    }

    // One endpoint force for the native rotational gradient/tangent and the
    // independent derivative oracle. The translational terms are assembled
    // separately from their public nodal force fields, exactly once.
    pub(super) fn link_body_load_force_n(&self, index: usize) -> Option<[f64; 3]> {
        self.gravity.as_ref().map(|load| {
            let weight = load.acceleration_m_s2.map(|g|
                0.5 * self.tension_links[index].reference_mass_kg() * g);
            load.buoyancy.as_ref().map_or(weight, |fluid| {
                let buoyancy = fluid.link_endpoint_forces_n()[index];
                std::array::from_fn(|axis| weight[axis] + buoyancy[axis])
            })
        })
    }

    /// Recover root-clamp force and body-frame moment through the SAME
    /// analytical assembly as Newton and additional-support recovery. Segment,
    /// root-link and gravity contributions enter once each, without any strain
    /// or coordinate perturbation outside the admitted physical state.
    /// Tip loads have zero root derivative. This is a numerical reaction, not a
    /// forced whole-body balance closure. The accepted pose is never modified.
    ///
    /// # Errors
    /// Source-domain, canonical-rotation, arithmetic or cancellation failures.
    pub fn root_support_reaction(
        &self, gate: &CancelGate,
    ) -> Result<RodSupportReaction, SpatialRodError> {
        checkpoint(gate)?;
        let energy = super::derivatives::SpatialEnergy::new(self);
        let assembled = energy.linearize(&self.rod, &[], false, &|| checkpoint(gate))?
            .ok_or(SpatialRodError::Invalid("native rod omitted its analytical reaction assembly"))?;
        let reaction = assembled.gradient[0];
        if !reaction.iter().all(|v| v.is_finite()) {
            return Err(SpatialRodError::Invalid("root reaction is unrepresentable"));
        }
        checkpoint(gate)?;
        Ok(RodSupportReaction { node: 0, force_n: [Some(reaction[0]), Some(reaction[1]), Some(reaction[2])],
            moment_body_nm: Some([reaction[3], reaction[4], reaction[5]]) })
    }
}
