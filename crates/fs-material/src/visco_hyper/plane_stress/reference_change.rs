//! Prescribed natural-strain events at fixed in-plane deformation and memory.
//!
//! First change the reference at fixed full F through the original material
//! transaction, then re-equilibrate the transverse column with its native
//! zero-duration plane-stress update. Both stages retain the same viscous
//! metrics: no physical time passes and no branch relaxation heat is generated.
//! This is an explicitly split prescribed event, not a resolved continuous
//! thermal path, growth law, caloric model or coordinate-frame remapping.
//!
//! Exact incompressibility adds a physical reference-actuation term. For the
//! native compression-positive pressure, `P = P_remainder - p cof(F)` and
//! `dPsi = P:dF + (dPsi/dG)_F:dG + p d(det G)`. The fixed-F remainder transaction
//! owns the middle term; `constraint_reference_work_pa` separately reports the
//! endpoint approximation `p_new (det G_new - det G_old)`. It is signed reference
//! input, not heat or a bulk penalty. The private off-constraint intermediate
//! is only an algebraic evaluation of the declared native remainder extension.
//! The final endpoint satisfies the new exact volume and transverse tractions.
//!
//! The native state retains its established fixed-F reference-work account and
//! frozen-metric loading ledger. For exact volume, that loading ledger includes
//! constraint actuation and is not physical transverse traction work alone.
//! The difference between endpoint mechanical plus constraint-reference work
//! and native loading is returned separately as a signed numerical defect.
//! A spatial owner must integrate and gate it before accepting the proposal.
//! Smaller events are needed when the defect is too large; this local update
//! does not certify a finite event's thermomechanical path or global equilibrium.
//! Original material domains also apply to the fixed-F target evaluation. A
//! large split event can therefore be refused even if a final reduced endpoint
//! exists; this method never widens the elastic domain to bridge that event.

use crate::MaterialError;
use crate::hyper::incompressible_plane_stress::IncompressiblePlaneStressState;
use crate::hyper::plane_stress::{PlaneStressSettings, PlaneStressState};
use crate::visco_hyper::{
    ViscoHyperelastic, ViscoRemainderState, ViscoRemainderStepResponse, ViscoState,
    ViscoStepResponse,
};

/// Immutable finite-compliance reference-change and transverse-equilibrium proposal.
#[derive(Clone, Debug)]
pub struct ViscoPlaneStressReferenceChange {
    /// Final native zero-duration update. `loading_work_pa` is the transverse
    /// frozen-metric energy increment; no branch relaxation loss is added.
    pub material: ViscoStepResponse,
    /// Converged target-reference reduction at the accepted in-plane F.
    pub plane_stress: PlaneStressState,
    /// Signed native reference input at fixed full F, per original volume (Pa).
    pub reference_work_pa: f64,
    /// Zero for finite compliance: there is no exact-volume constraint multiplier.
    pub constraint_reference_work_pa: f64,
    /// Physical endpoint `P_new[:,2] . (F_new[:,2]-F_old[:,2])`, in Pa.
    pub transverse_endpoint_work_pa: f64,
    /// Endpoint transverse work minus native loading, in Pa. This signed
    /// integration defect is separate from physical relaxation loss.
    pub endpoint_work_defect_pa: f64,
}

/// Immutable exact-volume reference change with its physical constraint work.
#[derive(Clone, Debug)]
pub struct ViscoIncompressiblePlaneStressReferenceChange {
    /// Native bulk-free update and final history. Its loading increment includes
    /// constraint actuation; its Piola stress excludes pressure. No heat is added.
    pub material: ViscoRemainderStepResponse,
    /// Exact target volume, reduced potential and full pressure-bearing stress.
    pub plane_stress: IncompressiblePlaneStressState,
    /// Fixed-full-F remainder reference input (Pa), retaining the native ledger.
    pub reference_work_pa: f64,
    /// Signed `p_new * (det G_new - det G_old)` (Pa), using the compression-positive
    /// native pressure. Total reference input is this plus `reference_work_pa`.
    /// This endpoint approximation is reported once, never as viscous heat.
    pub constraint_reference_work_pa: f64,
    /// Physical full-Piola transverse endpoint work (Pa), including pressure.
    pub transverse_endpoint_work_pa: f64,
    /// `transverse_endpoint_work_pa + constraint_reference_work_pa
    /// - material.loading_work_pa` (Pa). It is an endpoint quadrature/splitting
    /// defect, not the first-order physical pressure-volume reference input.
    pub endpoint_work_defect_pa: f64,
}

impl ViscoHyperelastic {
    /// Change prescribed G while retaining accepted in-plane deformation and
    /// viscous memory, then restore finite-compliance transverse equilibrium.
    ///
    /// Call on the target law with an unmodified thin-history output owned by
    /// `previous`. The original reference transaction admits identical physical
    /// equilibrium/ordered branch energies and domains; relaxation times may
    /// change. The same fixed natural-coordinate metrics and accumulated loss
    /// survive both zero-duration stages. An optional transverse seed may reuse
    /// the accepted column. Otherwise the native step uses its target-volume
    /// thickness guess and accepted shears; the reducer owns fallback and convergence.
    ///
    /// `reference_work_pa + material.loading_work_pa` is the storage increment.
    /// Physical endpoint work uses the returned full Piola. Their reported
    /// signed endpoint defect must be assessed separately by the caller, never
    /// added to material heat. Failure at either stage leaves `accepted` intact.
    pub fn change_plane_stress_reference_from(
        &self,
        previous: &Self,
        accepted: &ViscoState,
        initial_transverse: Option<[f64; 3]>,
        settings: PlaneStressSettings,
    ) -> Result<ViscoPlaneStressReferenceChange, MaterialError> {
        let in_plane = accepted_plane(&accepted.deformation_gradient)?;
        let reference = self.change_reference_from(previous, accepted)?;
        let endpoint = self.step_plane_stress_response(
            &reference.state,
            &in_plane,
            0.0,
            initial_transverse,
            settings,
        )?;
        let (transverse_endpoint_work_pa, endpoint_work_defect_pa) = work_accounts(
            &accepted.deformation_gradient,
            &endpoint.material.state.deformation_gradient,
            &endpoint.material.first_piola_pa,
            reference.reference_work_pa,
            0.0,
            endpoint.material.loading_work_pa,
        )?;
        Ok(ViscoPlaneStressReferenceChange {
            material: endpoint.material,
            plane_stress: endpoint.plane_stress,
            reference_work_pa: reference.reference_work_pa,
            constraint_reference_work_pa: 0.0,
            transverse_endpoint_work_pa,
            endpoint_work_defect_pa,
        })
    }

    /// Change prescribed G and restore exact-volume plane stress without a bulk
    /// penalty, branch relaxation, or invented pressure heat.
    ///
    /// The accepted history must already have thin kinematics and the previous
    /// exact volume. The final state has `det(F)=det(G_target)`. All original
    /// reference-change and remainder-split admission rules remain in force.
    /// Only the two transverse shears are seeded; thickness is imposed by the
    /// target natural volume. The private fixed-F intermediate may be off that
    /// constraint, and must never be published as a physical mixed endpoint.
    ///
    /// The unchanged native loading ledger includes pressure-volume reference
    /// actuation. Its physical separation is explicit in the response: total
    /// reference input is `reference_work_pa + constraint_reference_work_pa`,
    /// and adding physical transverse endpoint work minus `endpoint_work_defect_pa`
    /// recovers the storage increment. The caller must gate the signed defect;
    /// neither it nor the constraint input is nonnegative relaxation heat.
    pub fn change_incompressible_plane_stress_reference_from(
        &self,
        previous: &Self,
        accepted: &ViscoRemainderState,
        initial_shear: Option<[f64; 2]>,
        settings: PlaneStressSettings,
    ) -> Result<ViscoIncompressiblePlaneStressReferenceChange, MaterialError> {
        let in_plane = accepted_plane(&accepted.deformation_gradient)?;
        let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
        let old_volume = previous.reference_deformation().volume_ratio();
        let old_thickness = old_volume / determinant;
        if !(old_thickness.is_finite() && old_thickness > 0.0)
            || (accepted.deformation_gradient[8] / old_thickness - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(state_error(
                "accepted history violates its previous exact volume",
            ));
        }
        let reference = self.change_remainder_reference_from(previous, accepted)?;
        let endpoint = self.step_incompressible_plane_stress_response(
            &reference.state,
            &in_plane,
            0.0,
            initial_shear,
            settings,
        )?;
        let constraint_reference_work_pa = endpoint.plane_stress.pressure_pa
            * (self.reference_deformation().volume_ratio() - old_volume);
        let (transverse_endpoint_work_pa, endpoint_work_defect_pa) = work_accounts(
            &accepted.deformation_gradient,
            &endpoint.material.state.deformation_gradient,
            &endpoint.plane_stress.full_piola,
            reference.reference_work_pa,
            constraint_reference_work_pa,
            endpoint.material.loading_work_pa,
        )?;
        Ok(ViscoIncompressiblePlaneStressReferenceChange {
            material: endpoint.material,
            plane_stress: endpoint.plane_stress,
            reference_work_pa: reference.reference_work_pa,
            constraint_reference_work_pa,
            transverse_endpoint_work_pa,
            endpoint_work_defect_pa,
        })
    }
}

fn state_error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("viscoelastic plane-stress reference change: {what}"),
    }
}

fn accepted_plane(f: &[f64; 9]) -> Result<[f64; 4], MaterialError> {
    let determinant = f[0] * f[4] - f[1] * f[3];
    if !f.iter().all(|v| v.is_finite())
        || f[6] != 0.0
        || f[7] != 0.0
        || !(f[8] > 0.0 && determinant.is_finite() && determinant > 0.0)
    {
        return Err(state_error(
            "accepted history does not have admitted thin kinematics",
        ));
    }
    Ok([f[0], f[1], f[3], f[4]])
}

fn work_accounts(
    old: &[f64; 9],
    new: &[f64; 9],
    piola: &[f64; 9],
    reference_work: f64,
    constraint_work: f64,
    loading_work: f64,
) -> Result<(f64, f64), MaterialError> {
    let endpoint: f64 = [2, 5, 8]
        .iter()
        .map(|&i| piola[i] * (new[i] - old[i]))
        .sum();
    let defect = endpoint + constraint_work - loading_work;
    if ![
        endpoint,
        constraint_work,
        reference_work + constraint_work,
        defect,
    ]
    .iter()
    .all(|v| v.is_finite())
    {
        return Err(state_error(
            "unrepresentable reference or endpoint work account",
        ));
    }
    Ok((endpoint, defect))
}
