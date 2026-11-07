//! Source-resolved Kelvin--Voigt loss on the original CST/DKT shell.
//!
//! Each homogeneous facet retains its measured reduced plane-stress tensor,
//! fixed material state, thickness and in-plane orientation. Membrane and
//! bending resultants use the same physical thickness; no Poisson ratio or
//! modal decay time is used to invent a dissipative coefficient. The complete
//! modal matrix and passive channels come from the existing shell projection.

use super::{ModalPlateViscosity, PlateViscosityLimits, ShellViscosity, invalid, owner, poll};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::plane_stress_viscosity::PlaneStressViscosity;
use fs_modal::ModePair;
use fs_plate::shell::{ShellMesh, ShellModel};
use std::sync::Arc;

/// One homogeneous facet's source-bound physical loss, prepared once.
/// Clones share the immutable material record rather than copying its claims.
#[derive(Clone, Debug)]
pub struct ShellMaterialViscosity {
    source: Arc<PlaneStressViscosity>,
    thickness_m: f64,
    material_angle_rad: f64,
    section: ShellViscosity,
}

impl ShellMaterialViscosity {
    /// Integrate an explicitly characterized reduced tensor through thickness.
    /// `material_angle_rad` is counterclockwise from the original facet's local
    /// x axis toward its local y axis, about the oriented facet normal. Material
    /// axes are those already declared by `source.plane()`.
    ///
    /// This is a homogeneous plane-stress law about its midsurface: A = h V,
    /// B = 0 and D = h³ V / 12. It adds neither elastic stiffness nor inertia.
    /// The supplied mesh/model must already carry the corresponding physical
    /// geometry and thickness. Arbitrary 3-D viscoelastic plane-stress reduction
    /// and eccentric layered sections have separate constitutive owners.
    ///
    /// # Errors
    /// Invalid thickness/orientation, unrepresentable resultants, or the source
    /// owner's passivity refusal. Source frequency support is checked against
    /// actual retained modes by [`ModalShellMaterialViscosity::project`].
    pub fn new(
        source: Arc<PlaneStressViscosity>,
        thickness_m: f64,
        material_angle_rad: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let bending = source
            .bending_matrix_n_m_s(thickness_m, material_angle_rad)
            .map_err(owner)?;
        let reduced = source
            .rotated_matrix_pa_s(material_angle_rad)
            .map_err(owner)?;
        let membrane = reduced.map(|value| thickness_m * value);
        if membrane
            .iter()
            .zip(reduced)
            .any(|(a, v)| !a.is_finite() || (v != 0.0 && *a == 0.0))
        {
            return Err(invalid(
                "source shell membrane viscosity is unrepresentable",
            ));
        }
        let section = ShellViscosity::new(membrane, bending).map_err(owner)?;
        Ok(Self {
            source,
            thickness_m,
            material_angle_rad,
            section,
        })
    }

    /// Original material state, selected properties, plane and source support.
    pub fn source(&self) -> &PlaneStressViscosity {
        &self.source
    }
    /// Actual homogeneous through-thickness extent [m].
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }
    /// Material orientation in the original facet frame [rad].
    pub const fn material_angle_rad(&self) -> f64 {
        self.material_angle_rad
    }
    /// Prepared physical A/D resultants, with exact zero coupling B.
    pub const fn section(&self) -> &ShellViscosity {
        &self.section
    }
}

/// Complete modal material loss with its original facet sources and frequencies.
/// The adapter is cold preparation; propagation uses the existing passive
/// channels and heat accounting, with no material lookup or mesh work per step.
#[derive(Clone, Debug)]
pub struct ModalShellMaterialViscosity {
    facets: Vec<ShellMaterialViscosity>,
    angular_frequencies_rad_s: Vec<f64>,
    omega_band_rad_s: (f64, f64),
    modal: ModalPlateViscosity,
}

impl ModalShellMaterialViscosity {
    /// Bind one source-resolved viscosity per original facet to its retained
    /// eigenpairs. Actual modal frequencies must lie in every selected source
    /// band, even when a source tensor is zero. The query point alone is not an
    /// applicability test. All source states remain immutable in the result.
    ///
    /// Modes retain their order and mass normalization. The existing spatial
    /// and modal owners check the original mass Gram and retain every cross
    /// term, including directional and zero loss. Numerical drilling stiffness
    /// contributes no material damping. The same bare CST/DKT mesh, model and
    /// eigenpairs are required; a public `ShellModel` does not encode mesh or
    /// material identity, so those correspondences remain caller obligations.
    ///
    /// # Errors
    /// Missing facets, invalid/unsupported frequencies, invalid original basis,
    /// exhausted projection work, cancellation or constitutive refusal.
    pub fn project(
        mesh: &ShellMesh,
        facets: &[ShellMaterialViscosity],
        model: &ShellModel,
        modes: &[ModePair],
        limits: PlateViscosityLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if facets.len() != mesh.tris.len()
            || facets.is_empty()
            || modes.is_empty()
            || modes.len() > limits.max_modes
            || !(1..=64).contains(&limits.max_modes)
            || modes.iter().any(|mode| mode.phi.len() != model.free)
            || model
                .free
                .checked_mul(modes.len())
                .is_none_or(|entries| entries > limits.max_projection_terms)
        {
            return Err(invalid(
                "source shell loss needs bounded original facets and modes",
            ));
        }
        let mut band = (0.0_f64, f64::INFINITY);
        for facet in facets {
            poll(gate)?;
            let support = facet.source.omega_band_rad_s();
            band = (band.0.max(support.0), band.1.min(support.1));
        }
        if band.1 < band.0 {
            return Err(invalid(
                "shell viscosity sources have no common frequency support",
            ));
        }
        let mut frequencies = Vec::with_capacity(modes.len());
        for mode in modes {
            poll(gate)?;
            let omega = fs_math::det::sqrt(mode.lambda);
            if !omega.is_finite() || omega < band.0 || omega > band.1 {
                return Err(invalid(
                    "retained shell frequency lies outside viscosity source support",
                ));
            }
            frequencies.push(omega);
        }
        let sections: Vec<_> = facets.iter().map(|facet| facet.section.clone()).collect();
        let basis: Vec<_> = modes.iter().map(|mode| mode.phi.clone()).collect();
        let modal =
            ModalPlateViscosity::from_cst_dkt_shell(mesh, &sections, model, &basis, limits, gate)?;
        let facets = facets.to_vec();
        poll(gate)?;
        Ok(Self {
            facets,
            angular_frequencies_rad_s: frequencies,
            omega_band_rad_s: band,
            modal,
        })
    }

    /// Original facet order, source states, thicknesses and orientations.
    pub fn facets(&self) -> &[ShellMaterialViscosity] {
        &self.facets
    }
    /// Frequencies of the exact retained eigenpairs used for applicability [rad/s].
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.angular_frequencies_rad_s
    }
    /// Intersection of all selected viscosity source bands [rad/s].
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }
    /// Full modal matrix and existing physical dashpot/heat connections.
    pub const fn modal(&self) -> &ModalPlateViscosity {
        &self.modal
    }
}

#[cfg(test)]
mod tests;
