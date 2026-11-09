//! Explicit reference strain and support motion for native axial elasticity.
//!
//! These inputs describe a frozen small-strain state. A thermal caller supplies
//! the stress-free strain obtained from its own expansion law and temperature
//! reference; no coefficient, stress or altered Young's modulus is inferred.

use super::{Member, checkpoint};
use crate::{GroundStructure, LayoutCase, TrussConstructionError};
use fs_exec::Cx;

/// Borrowed imposed state in the complete original ground ordering.
///
/// Omitted fields mean zero everywhere without allocating zero-filled arrays.
/// A supplied field must include every original member/node, including removed
/// members. Values on removed members are admitted but contribute no force.
/// Prescribed motion is allowed only on coordinates marked supported by the
/// supplied [LayoutCase]; free-coordinate entries must be exactly zero.
#[derive(Debug, Clone, Copy, Default)]
pub struct ElasticImposedState<'a> {
    /// Axial stress-free engineering strain, positive for natural expansion.
    ///
    /// This can represent thermal expansion, fabrication misfit or another
    /// explicitly supplied eigenstrain. It changes the stress-free extension
    /// by original_length * strain, without changing geometry or source moduli.
    pub member_stress_free_strain: Option<&'a [f64]>,
    /// Prescribed global x/y displacement at each supported coordinate [m].
    ///
    /// Original coordinates remain the reference geometry. This is boundary
    /// motion in the infinitesimal model, not a geometry update or rotation.
    pub prescribed_displacements_m: Option<&'a [[f64; 2]]>,
}

impl ElasticImposedState<'_> {
    pub(super) fn admit(
        self,
        ground: &GroundStructure,
        case: &LayoutCase,
        cx: &Cx<'_>,
    ) -> Result<(), TrussConstructionError> {
        checkpoint(cx)?;
        if let Some(strains) = self.member_stress_free_strain {
            if strains.len() != ground.members().len() {
                return Err(TrussConstructionError::VectorLength {
                    field: "elastic stress-free strains",
                    expected: ground.members().len(),
                    actual: strains.len(),
                });
            }
            for (index, strain) in strains.iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                if !strain.is_finite() {
                    return Err(super::invalid(
                        "elastic stress-free strains",
                        "must be finite",
                    ));
                }
            }
        }
        if let Some(displacements) = self.prescribed_displacements_m {
            if displacements.len() != ground.nodes().len() {
                return Err(TrussConstructionError::VectorLength {
                    field: "elastic prescribed displacements",
                    expected: ground.nodes().len(),
                    actual: displacements.len(),
                });
            }
            for (node, displacement) in displacements.iter().enumerate() {
                if node % 256 == 0 {
                    checkpoint(cx)?;
                }
                for (axis, &value) in displacement.iter().enumerate() {
                    if !value.is_finite() || (!case.supported()[node][axis] && value != 0.0) {
                        return Err(super::invalid(
                            "elastic prescribed displacements",
                            "require finite values and zero on free coordinates",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Equivalent axial load before projecting onto free nodal coordinates.
    pub(super) fn equivalent_force(
        self,
        index: usize,
        length: f64,
        member: &Member,
    ) -> Result<f64, TrussConstructionError> {
        let reference_extension = length * self.stress_free_strain(index);
        let prescribed_extension = member.extension(|slot| self.displacement(slot / 2, slot % 2));
        let force = member.stiffness * (reference_extension - prescribed_extension);
        if [reference_extension, prescribed_extension, force]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(super::invalid(
                "elastic imposed loads",
                "reference extension, prescribed extension and equivalent force must be finite",
            ));
        }
        Ok(force)
    }

    pub(super) fn is_supplied(self) -> bool {
        self.member_stress_free_strain.is_some() || self.prescribed_displacements_m.is_some()
    }

    pub(super) fn stress_free_strain(self, member: usize) -> f64 {
        self.member_stress_free_strain
            .map_or(0.0, |values| values[member])
    }

    pub(super) fn displacement(self, node: usize, axis: usize) -> f64 {
        self.prescribed_displacements_m
            .map_or(0.0, |values| values[node][axis])
    }
}
