//! Reversible multiwinding saturation derived from the original material path.
//!
//! Each winding has total linkage `lambda_i = N_i Phi + L_leak_i I_i`.
//! All signed ampere-turns act on ONE physical core and gap:
//! `sum N_i I_i = length H(Phi/area) + R_gap Phi`.
//! The native coupled storage integrates this potential and retains the
//! original Kirchhoff flux states, exact finite-step work and analytic tangent.
//! Core energy is never duplicated per winding. These fixed-geometry,
//! field-only laws do not infer hysteresis, iron loss, fringing or motion.

use super::{
    MAX_MATERIAL_CIRCUIT_BRANCHES, MaterialCircuitError, MaterialReversibleInductor,
    ReversibleInductorGeometry, finite, invalid,
};
use fs_material::field_storage::nonlinear::ReversibleMagneticMaterial;
use fs_phs::circuit::inductance::{MagneticStorageBlock, SharedCoreInductance};
use std::sync::Arc;

/// One fixed common magnetic path; air-gap reluctance is explicitly supplied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReversibleCoreGeometry {
    /// Uniform core cross-section [m²].
    pub core_area_m2: f64,
    /// Uniform material path length [m].
    pub core_length_m: f64,
    /// Additional series gap reluctance [1/H]; zero neglects the gap.
    pub gap_reluctance_per_h: f64,
}

/// One original circuit branch linking the common flux with signed turns.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedCoreWinding {
    /// Index in the ORIGINAL CircuitGraph, not compressed inductor order.
    pub branch: usize,
    /// Nonzero effective turns, signed relative to the branch orientation.
    pub turns: f64,
    /// Positive explicitly authored leakage [H]. Zero-leakage singular
    /// constraints need an ideal-transformer representation, not guessed loss.
    pub leakage_inductance_h: f64,
}

/// Immutable source, uniform geometry, windings and the native shared potential.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialSharedCoreInductor {
    // The existing single-turn, zero-leakage source mapper already owns exact
    // signed Phi/MMF knots, field-only admission and source B/H observations.
    // This private reference is not an additional circuit branch or store.
    reference: MaterialReversibleInductor,
    windings: Arc<[SharedCoreWinding]>,
    law: SharedCoreInductance,
}

/// Endpoint observations in declared winding order. All energy components
/// partition the native magnetic store; they must not be added to it again.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedCoreOperatingPoint {
    pub flux_linkages_wb: Vec<f64>,
    pub currents_a: Vec<f64>,
    pub common_flux_wb: f64,
    pub core_flux_density_t: f64,
    pub core_field_strength_a_m: f64,
    /// Reciprocal row-major matrix d(lambda_i)/d(I_j) [H]. At source kinks,
    /// uses the native one-sided convention rather than claiming smoothness.
    pub differential_inductance_h: Vec<f64>,
    pub core_energy_j: f64,
    pub gap_energy_j: f64,
    pub leakage_energy_j: Vec<f64>,
    pub stored_energy_j: f64,
}

impl MaterialSharedCoreInductor {
    /// Retain every original H(B) breakpoint and signed winding orientation.
    /// The same existing material admission refuses unsupported thermal
    /// free-energy dependence, extrapolation and nonrepresentable geometry.
    /// Explicit positive leakage gives independent nonsingular flux states.
    pub fn new(
        windings: Vec<SharedCoreWinding>,
        material: ReversibleMagneticMaterial,
        geometry: ReversibleCoreGeometry,
    ) -> Result<Self, MaterialCircuitError> {
        if windings.is_empty()
            || windings.len() > MAX_MATERIAL_CIRCUIT_BRANCHES
            || windings.iter().enumerate().any(|(i, winding)| {
                winding.branch >= MAX_MATERIAL_CIRCUIT_BRANCHES
                    || windings[..i]
                        .iter()
                        .any(|other| other.branch == winding.branch)
            })
        {
            return Err(invalid(
                "bounded distinct original shared-core winding branches required",
            ));
        }
        let reference = MaterialReversibleInductor::new(
            windings[0].branch,
            material,
            ReversibleInductorGeometry {
                turns: 1.0,
                core_area_m2: geometry.core_area_m2,
                core_length_m: geometry.core_length_m,
                gap_reluctance_per_h: geometry.gap_reluctance_per_h,
                leakage_inductance_h: 0.0,
            },
        )?;
        let law = SharedCoreInductance::try_new(
            windings.iter().map(|winding| winding.turns).collect(),
            windings
                .iter()
                .map(|winding| winding.leakage_inductance_h)
                .collect(),
            reference.curve().knots().to_vec(),
        )?;
        Ok(Self {
            reference,
            windings: windings.into(),
            law,
        })
    }

    pub fn material(&self) -> &ReversibleMagneticMaterial {
        self.reference.material()
    }

    pub fn geometry(&self) -> ReversibleCoreGeometry {
        let geometry = self.reference.geometry();
        ReversibleCoreGeometry {
            core_area_m2: geometry.core_area_m2,
            core_length_m: geometry.core_length_m,
            gap_reluctance_per_h: geometry.gap_reluctance_per_h,
        }
    }

    pub fn windings(&self) -> &[SharedCoreWinding] {
        &self.windings
    }

    /// Complete common-core-plus-leakage potential, shared with native assembly.
    pub fn law(&self) -> &SharedCoreInductance {
        &self.law
    }

    pub fn binding(&self) -> MagneticStorageBlock {
        MagneticStorageBlock {
            branches: self.windings.iter().map(|winding| winding.branch).collect(),
            law: self.law.clone(),
        }
    }

    /// Initial total linkages [Wb] from simultaneous winding currents [A].
    /// Saturation is determined by their summed signed ampere-turns.
    pub fn flux_linkages_wb(&self, currents_a: &[f64]) -> Result<Vec<f64>, MaterialCircuitError> {
        Ok(self.law.flux_linkages_wb(currents_a)?)
    }

    /// Observe the common material and all leakage stores without a field
    /// solve, new fit, loss model or state transition.
    pub fn operating_point(
        &self,
        flux_linkages_wb: &[f64],
    ) -> Result<SharedCoreOperatingPoint, MaterialCircuitError> {
        let common_flux_wb = self.law.common_flux_wb(flux_linkages_wb)?;
        let reference = self.reference.operating_point(common_flux_wb)?;
        let currents_a = self.law.currents_a(flux_linkages_wb)?;
        let leakage_energy_j = self
            .windings
            .iter()
            .zip(&currents_a)
            .map(|(winding, &current)| {
                finite((0.5 * winding.leakage_inductance_h * current) * current)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(SharedCoreOperatingPoint {
            flux_linkages_wb: flux_linkages_wb.to_vec(),
            currents_a,
            common_flux_wb,
            core_flux_density_t: reference.core_flux_density_t,
            core_field_strength_a_m: reference.core_field_strength_a_m,
            differential_inductance_h: self.law.differential_inductance_h(flux_linkages_wb)?,
            core_energy_j: reference.core_energy_j,
            gap_energy_j: reference.gap_energy_j,
            leakage_energy_j,
            stored_energy_j: self.law.stored_energy_j(flux_linkages_wb)?,
        })
    }
}
