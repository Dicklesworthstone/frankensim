//! Explicit specimen direction and source-bounded reversible recoil material.
//!
//! These are per-cell constitutive choices on the existing planar owner. The
//! original FEEC stiffness, current loading, boundary admission and checked CG
//! remain shared. No separate permanent-magnet field solver is introduced.

use super::{MagneticError, finite};
use fs_material::{
    field_storage::permanent::ResolvedRecoilMaterial, state_point::ResolvedMaterialStatePoint,
};

/// Constitutive model for one original triangle. Linear induced material and
/// permanent magnets may share a mesh, with their distinct source roles kept.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarMagneticMaterial {
    /// Induced B=mu H using the canonical source-resolved magnetic permeability.
    Linear(ResolvedMaterialStatePoint),
    /// Reversible B=mu_recoil H+Br n, within all original recoil limits.
    Recoil(PlanarRecoilMaterial),
}

/// A complete source law and its independently supplied specimen orientation.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanarRecoilMaterial {
    material: ResolvedRecoilMaterial,
    direction: [f64; 2],
    remanence_t: [f64; 2],
}
impl PlanarRecoilMaterial {
    /// Require an explicit unit vector in the original planar mesh frame. A
    /// norm error at most 1e-12 is treated as roundoff and normalized once.
    /// Direction is never inferred from a material name or geometry ordering.
    pub fn new(
        material: ResolvedRecoilMaterial,
        direction: [f64; 2],
    ) -> Result<Self, MagneticError> {
        let norm = direction[0].hypot(direction[1]);
        if !direction.iter().all(|v| v.is_finite())
            || !norm.is_finite()
            || (norm - 1.0).abs() > 1e-12
        {
            return Err(MagneticError::Invalid(
                "permanent-magnet direction must be a finite unit vector in the mesh frame",
            ));
        }
        let direction = direction.map(|v| v / norm);
        let remanence_t = direction.map(|v| material.remanence_t() * v);
        for value in remanence_t {
            finite(value)?;
            finite(value / material.permeability_h_m())?;
        }
        if remanence_t == [0.0; 2] {
            return Err(MagneticError::Invalid(
                "remanence vector is not representable",
            ));
        }
        Ok(Self {
            material,
            direction,
            remanence_t,
        })
    }
    /// Complete original material law, query coordinates and property receipts.
    pub fn material(&self) -> &ResolvedRecoilMaterial {
        &self.material
    }
    /// Admitted unit remanence direction in original mesh coordinates.
    pub const fn direction(&self) -> [f64; 2] {
        self.direction
    }
    /// Source magnitude times admitted specimen direction [T].
    pub const fn remanence_t(&self) -> [f64; 2] {
        self.remanence_t
    }
}

/// Reference used by the reported fixed-state constitutive energy. These are
/// field-storage potentials, not a magnetic body's total formation energy,
/// circuit source work, irreversible heat or a force/coenergy calculation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MagneticEnergyReference {
    /// Induced cells use their source-integrated energy, with zero at B=0.
    ZeroFlux,
    /// Magnet cells use (B-Br).H/2=mu_recoil*|H|²/2, with zero at B=Br;
    /// induced cells retain their zero-flux energy. Br and material conditions
    /// are held fixed.
    FixedRemanence,
}
