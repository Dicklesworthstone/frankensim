//! Prescribed stress-free strain changes on an installed membrane.
//!
//! The original plane-stress reference transaction changes the natural
//! configuration at fixed in-plane deformation and fixed material-frame viscous
//! metrics, then restores transverse traction equilibrium. No physical time
//! passes and no relaxation heat is generated. The same admitted endpoint is
//! condensed into the existing Green-strain sheet coordinates, without another
//! transverse solve or a new constitutive law.
//!
//! Fixed-full-F reference work, exact-volume constraint actuation, native
//! transverse loading and its signed endpoint defect remain separate. A spatial
//! owner must integrate and gate these accounts before accepting the event.
//! This is a prescribed split event, not a continuous thermal path, a growth
//! evolution law, a change of material coordinates or global force equilibrium.

use super::ViscoMembraneState;
use crate::hyper::cards::FiniteStrainPotential;
use crate::hyper::incompressible_plane_stress::{
    IncompressiblePlaneStressPotential, generic_response,
};
use crate::hyper::plane_stress::{
    PlaneStressResponse, PlaneStressSettings, plane_stress_response_from,
};
use crate::sheet::hyperelastic::{admitted_representative, response_from_plane_stress};
use crate::sheet::{SheetError, SheetResponse};
use crate::visco_hyper::ViscoHyperelastic;

/// Complete immutable reference-change proposal per original reference area.
/// The native cumulative mechanical/reference/loss ledgers remain in `state`.
#[derive(Clone, Debug)]
pub struct ViscoMembraneReferenceChange {
    pub state: ViscoMembraneState,
    /// Zero-duration endpoint potential, second-Piola resultants and exact
    /// Green-strain tangent. Physical stored energy is reported separately.
    pub sheet: SheetResponse,
    /// Physical Piola stress in the canonical Cholesky spatial frame [Pa].
    /// Exact-volume states include their compression-positive pressure term.
    pub first_piola_pa: [f64; 9],
    pub pressure_pa: Option<f64>,
    pub stored_before_j_m2: f64,
    pub stored_after_j_m2: f64,
    /// Signed native reference input at fixed FULL deformation. Adding
    /// `native_transverse_loading_work_j_m2` gives the stored-energy change.
    pub reference_work_j_m2: f64,
    /// Signed pressure-volume reference input for exact incompressibility;
    /// zero for finite compliance. It is not relaxation heat.
    pub constraint_reference_work_j_m2: f64,
    /// Native frozen-metric loading increment during transverse equilibration.
    /// For exact volume this includes constraint actuation and is not solely
    /// work of physical transverse tractions.
    pub native_transverse_loading_work_j_m2: f64,
    /// Physical endpoint transverse Piola dotted with its deformation change.
    pub transverse_endpoint_work_j_m2: f64,
    /// Transverse endpoint plus constraint-reference work minus native loading.
    /// The spatial owner must gate this signed finite-event integration defect.
    pub endpoint_work_defect_j_m2: f64,
    pub relative_transverse_residual: f64,
}

impl ViscoHyperelastic {
    /// Change prescribed natural strain while retaining installed membrane
    /// strain, original thickness and every accepted material-frame metric.
    ///
    /// Call on the TARGET law, supplying the previous law and its unmodified
    /// membrane history. Equilibrium energy, elastic domain and ordered branch
    /// energies must match; prescribed reference and relaxation times may
    /// change under the original reference owner's admission rules. All input
    /// states remain unchanged on refusal, including a late tangent refusal.
    ///
    /// `strain = [E11,E22,2E12]` must reproduce the accepted canonical in-plane
    /// deformation exactly through the shared membrane representative. An event
    /// cannot silently add in-plane loading or rename the material frame. The
    /// native transverse solver owns the seed, stable branch, domain and work
    /// limits. Its exact accepted deformation is reused for tangent condensation;
    /// no second transverse solve can change the proposed physical endpoint.
    ///
    /// Work is scaled by the ORIGINAL positive reference thickness. A signed
    /// reference input or endpoint defect is never made into physical heat.
    /// Large split events can leave the native intermediate elastic domain even
    /// when a reduced endpoint exists; reduce the event rather than widen it.
    pub fn change_membrane_reference_from(
        &self,
        previous: &Self,
        accepted: &ViscoMembraneState,
        strain: [f64; 3],
        reference_thickness_m: f64,
        settings: PlaneStressSettings,
    ) -> Result<ViscoMembraneReferenceChange, SheetError> {
        if !reference_thickness_m.is_finite() || reference_thickness_m <= 0.0 {
            return Err(SheetError::Invalid(
                "finite positive original membrane reference thickness required",
            ));
        }
        let in_plane = admitted_representative(strain)?;
        let old_f = accepted.deformation_gradient();
        if [old_f[0], old_f[1], old_f[3], old_f[4]] != in_plane {
            return Err(SheetError::Invalid(
                "reference event requires unchanged canonical in-plane membrane deformation",
            ));
        }
        let (
            state,
            reduced,
            first_piola_pa,
            pressure_pa,
            reference_work_pa,
            constraint_reference_work_pa,
            native_transverse_loading_work_pa,
            transverse_endpoint_work_pa,
            endpoint_work_defect_pa,
        ) = match accepted {
            ViscoMembraneState::Full(old) => {
                let change = self.change_plane_stress_reference_from(
                    previous,
                    old,
                    Some([old_f[2], old_f[5], old_f[8]]),
                    settings,
                )?;
                let tangent = self
                    .freeze(&change.material.state, 0.0)?
                    .consistent_tangent(&change.plane_stress.deformation_gradient)?;
                let reduced = plane_stress_response_from(change.plane_stress, tangent)?;
                (
                    ViscoMembraneState::Full(change.material.state),
                    reduced,
                    change.material.first_piola_pa,
                    None,
                    change.reference_work_pa,
                    change.constraint_reference_work_pa,
                    change.material.loading_work_pa,
                    change.transverse_endpoint_work_pa,
                    change.endpoint_work_defect_pa,
                )
            }
            ViscoMembraneState::Incompressible(old) => {
                let change = self.change_incompressible_plane_stress_reference_from(
                    previous,
                    old,
                    Some([old_f[2], old_f[5]]),
                    settings,
                )?;
                let remainder = self
                    .freeze_remainder(&change.material.state, 0.0)?
                    .incompressible_remainder_response(&change.plane_stress.deformation_gradient)?;
                let reduced = generic_response(change.plane_stress, Some(remainder), &in_plane)?;
                (
                    ViscoMembraneState::Incompressible(change.material.state),
                    PlaneStressResponse {
                        deformation_gradient: reduced.deformation_gradient,
                        energy_density_pa: reduced.energy_density_pa,
                        piola: reduced.piola,
                        tangent: reduced.tangent,
                        relative_transverse_residual: reduced.relative_transverse_residual,
                        iterations: reduced.iterations,
                    },
                    reduced.full_piola,
                    Some(reduced.pressure_pa),
                    change.reference_work_pa,
                    change.constraint_reference_work_pa,
                    change.material.loading_work_pa,
                    change.transverse_endpoint_work_pa,
                    change.endpoint_work_defect_pa,
                )
            }
        };
        let sheet = response_from_plane_stress(strain, reduced, reference_thickness_m)?;
        let [
            stored_before_j_m2,
            stored_after_j_m2,
            reference_work_j_m2,
            constraint_reference_work_j_m2,
            native_transverse_loading_work_j_m2,
            transverse_endpoint_work_j_m2,
            endpoint_work_defect_j_m2,
        ] = [
            accepted.stored_pa(),
            state.stored_pa(),
            reference_work_pa,
            constraint_reference_work_pa,
            native_transverse_loading_work_pa,
            transverse_endpoint_work_pa,
            endpoint_work_defect_pa,
        ]
        .map(|work| reference_thickness_m * work);
        if ![
            stored_before_j_m2,
            stored_after_j_m2,
            reference_work_j_m2,
            constraint_reference_work_j_m2,
            native_transverse_loading_work_j_m2,
            transverse_endpoint_work_j_m2,
            endpoint_work_defect_j_m2,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(SheetError::Invalid(
                "unrepresentable membrane reference storage or work",
            ));
        }
        Ok(ViscoMembraneReferenceChange {
            state,
            sheet,
            first_piola_pa,
            pressure_pa,
            stored_before_j_m2,
            stored_after_j_m2,
            reference_work_j_m2,
            constraint_reference_work_j_m2,
            native_transverse_loading_work_j_m2,
            transverse_endpoint_work_j_m2,
            endpoint_work_defect_j_m2,
            relative_transverse_residual: reduced.relative_transverse_residual,
        })
    }
}
