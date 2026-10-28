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
use crate::rod::derivatives::{RodEnergy, perturb};

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
}

impl SpatialMaterialRod {
    /// Reference inventory, retained even when current material density changes.
    #[must_use]
    pub fn reference_mass(&self) -> &RodReferenceMass { &self.reference_mass }

    /// Set an explicit uniform world-frame gravity vector on a new static
    /// specimen. Zero removes gravity; there is no implicit terrestrial preset.
    /// This replaces the old gravity field, never accumulates another weight.
    /// Reference mass and every support are preserved. No buoyancy is inferred.
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
        let mut candidate = self.clone();
        candidate.gravity = if acceleration_m_s2.iter().all(|v| *v == 0.0) {
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
            Some(GravityLoad { acceleration_m_s2, force_n })
        };
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
    fn gravity_potential_for(&self, rod: &Rod) -> Result<f64, SpatialRodError> {
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
        let value = internal + self.gravity_potential_for(rod)?;
        if !value.is_finite() {
            return Err(SpatialRodError::Invalid("elastic plus gravity potential is unrepresentable"));
        }
        Ok(value)
    }

    /// Recover root-clamp force and body-frame moment through the SAME local
    /// derivative owner as Newton and additional-support recovery. Only incident
    /// segment energies are differenced; root link forces and gravitational
    /// force/offset moments are included analytically, once each.
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
        energy.value(&self.rod)?;
        let mut plus = self.rod.clone();
        let mut minus = self.rod.clone();
        let mut reaction = [0.0; 6];
        for (coordinate, value) in reaction.iter_mut().enumerate() {
            checkpoint(gate)?;
            let eps = super::super::fd_step(coordinate, self.probe_translation_scale_m(), 3e-6)?;
            perturb(&mut plus, 0, coordinate, eps)?;
            perturb(&mut minus, 0, coordinate, -eps)?;
            *value = energy.derivative(&self.rod, &plus, &minus, 0, coordinate, eps)?;
            plus.positions[0] = self.rod.positions[0];
            minus.positions[0] = self.rod.positions[0];
            plus.frames[0] = self.rod.frames[0];
            minus.frames[0] = self.rod.frames[0];
        }
        if !reaction.iter().all(|v| v.is_finite()) {
            return Err(SpatialRodError::Invalid("root reaction is unrepresentable"));
        }
        checkpoint(gate)?;
        Ok(RodSupportReaction { node: 0, force_n: [Some(reaction[0]), Some(reaction[1]), Some(reaction[2])],
            moment_body_nm: Some([reaction[3], reaction[4], reaction[5]]) })
    }
}
