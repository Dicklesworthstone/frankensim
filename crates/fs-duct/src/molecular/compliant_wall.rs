//! Molecular gas and source-owned locally reacting structural walls.
//!
//! The gas retains its original molecular/shear/thermal constitutive owner.
//! WallModel adds the mechanical volume admittance in parallel with that gas
//! compression: Y_wall = 2 pi r sqrt(1 + r'^2) / Z_surface. This is wall
//! displacement, not a second thermal boundary layer or a fitted gas loss.
//! Every bore, profile quadrature point and finite chimney uses this same map.
//!
//! A wall law may include structural dissipation or radiation; the acoustic
//! work delivered to it must not be described as gas heat. The boundary remains
//! locally reacting: shell propagation, nonlocal wall coupling, large wall
//! motion and thermoelastic coupling require their respective source owners.

use super::*;
use crate::profile::RiccatiSteps;
use crate::{WallLaw, WallModel};
use fs_material::gas_field::AxialGasField;
use fs_material::relaxation::RelaxationAbsorption;
use std::sync::Arc;

/// Shared immutable material-wall response. Its original owner supplies the
/// impedance, source data and frequency/geometric applicability admission.
/// Thread-safe ownership lets cloned duct/runtime bindings retain one law.
pub type MolecularCompliantWall = Arc<dyn WallLaw + Send + Sync>;

impl core::fmt::Debug for MolecularDuct {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("MolecularDuct")
            .field("state", &self.state)
            .field("field_binding", &self.field_binding)
            .field("points", &self.points)
            .field("pieces", &self.pieces)
            .field("reference_impedance", &self.reference_impedance)
            .field("compliant_wall", &self.compliant_wall.is_some())
            .finish()
    }
}

impl MolecularDuct {
    /// Bind a source-owned mechanical wall to the original molecular gas.
    /// `thermal_wall` controls gas shear/heat diffusion independently; None
    /// still permits wall displacement. ConeSlices is live for either wall
    /// effect, and profile geometry supplies its actual surface-area factor.
    ///
    /// All existing response, chimney-flow and coupled-exterior APIs use this
    /// same immutable binding. The wall is evaluated at each actual radius and
    /// frequency, so a source-domain refusal propagates without substitution.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_compliant_wall(
        duct: &Duct,
        state: RelaxationState,
        thermal_wall: Option<MolecularWallResolution>,
        wall: MolecularCompliantWall,
        cones: ConeSlices,
        maximum_blocks: usize,
    ) -> Result<Self, MolecularDuctError> {
        Self::new_with_compliant_wall_controlled(
            duct,
            state,
            thermal_wall,
            wall,
            cones,
            maximum_blocks,
            &mut || false,
        )
    }

    /// Cancellable cold binding; frequency-dependent wall admission occurs in
    /// the response call, before any completed response can be returned.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_compliant_wall_controlled(
        duct: &Duct,
        state: RelaxationState,
        thermal_wall: Option<MolecularWallResolution>,
        wall: MolecularCompliantWall,
        cones: ConeSlices,
        maximum_blocks: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularDuctError> {
        Self::new_bound_controlled(
            duct,
            state,
            thermal_wall,
            Some(wall),
            cones,
            maximum_blocks,
            cancelled,
        )
    }

    /// Primitive temperature/composition fields and a mechanical material wall
    /// on the same original field quadrature. Gas-state and radial-spectrum
    /// caches remain shared; surface slopes belong to geometry, not gas keys.
    #[allow(clippy::too_many_arguments)]
    pub fn from_gas_field_with_compliant_wall(
        duct: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        thermal_wall: Option<MolecularWallResolution>,
        wall: MolecularCompliantWall,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        maximum_blocks: usize,
    ) -> Result<Self, MolecularDuctError> {
        Self::from_gas_field_with_compliant_wall_controlled(
            duct,
            field,
            law,
            thermal_wall,
            wall,
            steps_per_span,
            uniform_cones,
            maximum_blocks,
            &mut || false,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn from_gas_field_with_compliant_wall_controlled(
        duct: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        thermal_wall: Option<MolecularWallResolution>,
        wall: MolecularCompliantWall,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        maximum_blocks: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularDuctError> {
        Self::from_gas_field_bound_controlled(
            duct,
            field,
            law,
            thermal_wall,
            Some(wall),
            steps_per_span,
            uniform_cones,
            maximum_blocks,
            cancelled,
        )
    }

    /// The exact source-owned wall binding; None means a mechanically rigid
    /// surface, independently of any retained gas shear/thermal diffusion.
    pub fn compliant_wall(&self) -> Option<&(dyn WallLaw + Send + Sync)> {
        self.compliant_wall.as_deref()
    }

    /// Cache the *composed* response for identical material/radius/slope
    /// samples within one frequency preparation. A radius shared by two
    /// different surface slopes never aliases the physical wall area.
    pub(super) fn wall_wave(
        &self,
        point: usize,
        slant: f64,
        omega: f64,
        waves: &[SegmentWave],
        cache: &mut BTreeMap<(usize, u64), SegmentWave>,
    ) -> Result<SegmentWave, MolecularDuctError> {
        let Some(wall) = &self.compliant_wall else {
            return Ok(waves[point]);
        };
        let key = (point, slant.to_bits());
        if let Some(&wave) = cache.get(&key) {
            return Ok(wave);
        }
        let wave = WallModel::Law(wall.as_ref()).apply(
            waves[point],
            self.points[point].radius,
            omega,
            slant,
        )?;
        cache.insert(key, wave);
        Ok(wave)
    }
}

#[cfg(test)]
#[path = "compliant_wall_tests.rs"]
mod tests;
