//! Original ordered shear claims applied over the normal owner's real disk.

/// Pinned thermal junction strengths on explicitly frozen normal geometry.
pub mod temperature;

/// Original junction claims rebound to a thermally updated normal-contact disk.
pub mod normal_temperature;

use core::f64::consts::PI;

use fs_tribo::adhesive::spinning::{
    AdhesiveSlidingPatch, AdhesiveSlidingPatchParameters, AdhesiveSlidingWrench,
};

use super::{AdhesiveShearError, BoundNormalInterfaceState, ContentHash, SourceAdhesiveShear};
use crate::normal_patch::{adhesive::JkrResponse, dmt::DmtResponse, maugis::MaugisResponse};

/// Original immutable source owners plus one actual circular solid-contact area.
/// Both constitutive leaves use the same local uniform-shear law. The static
/// leaf supplies support points of its full force/moment capacity set; evaluating
/// it alone does not select a static reaction. Geometry changes require a new
/// normal response, and no temperature, load or provenance is silently refreshed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceAdhesiveSlidingPatch<'a> {
    normal: &'a BoundNormalInterfaceState,
    shear: &'a SourceAdhesiveShear,
    patch: AdhesiveSlidingPatch,
    static_patch: AdhesiveSlidingPatch,
    contact_area_m2: f64,
}

impl<'a> SourceAdhesiveSlidingPatch<'a> {
    /// Original ordered normal source and its unchanged physical state query.
    #[must_use]
    pub const fn normal_state(&self) -> &'a BoundNormalInterfaceState {
        self.normal
    }

    /// Original static/kinetic junction claims with all resolution receipts.
    #[must_use]
    pub const fn shear_source(&self) -> &'a SourceAdhesiveShear {
        self.shear
    }

    /// Actual selected kinetic strength and true circular radius.
    #[must_use]
    pub const fn patch(&self) -> &AdhesiveSlidingPatch {
        &self.patch
    }

    /// Static support law at the independently selected static junction strength.
    #[must_use]
    pub const fn static_patch(&self) -> &AdhesiveSlidingPatch {
        &self.static_patch
    }

    /// Original solid area, excluding noncontact attraction or cohesive annuli.
    #[must_use]
    pub const fn contact_area_m2(&self) -> f64 {
        self.contact_area_m2
    }

    /// Physical two-axis force, normal-axis spin moment and unallocated heat rate.
    pub fn evaluate(
        &self,
        slip_velocity_m_s: [f64; 2],
        spin_rad_s: f64,
    ) -> Result<AdhesiveSlidingWrench, AdhesiveShearError> {
        Ok(self.patch.evaluate(slip_velocity_m_s, spin_rad_s)?)
    }
}

impl SourceAdhesiveShear {
    /// JKR attached area carries spin resistance at zero or tensile net load.
    pub fn sliding_patch_jkr<'a>(
        &'a self,
        normal: &'a BoundNormalInterfaceState,
        response: &JkrResponse,
    ) -> Result<SourceAdhesiveSlidingPatch<'a>, AdhesiveShearError> {
        self.sliding_patch(
            normal,
            response.normal_state_identity(),
            response.contact_radius_m,
            response.contact_area_m2,
        )
    }

    /// DMT finite-gap attraction has no solid disk and therefore no spin torque.
    pub fn sliding_patch_dmt<'a>(
        &'a self,
        normal: &'a BoundNormalInterfaceState,
        response: &DmtResponse,
    ) -> Result<SourceAdhesiveSlidingPatch<'a>, AdhesiveShearError> {
        self.sliding_patch(
            normal,
            response.normal_state_identity(),
            response.contact_radius_m,
            response.contact_area_m2,
        )
    }

    /// The Maugis cohesive annulus is excluded from the junction-shear disk.
    pub fn sliding_patch_maugis<'a>(
        &'a self,
        normal: &'a BoundNormalInterfaceState,
        response: &MaugisResponse,
    ) -> Result<SourceAdhesiveSlidingPatch<'a>, AdhesiveShearError> {
        self.sliding_patch(
            normal,
            response.normal_state_identity(),
            response.contact_radius_m,
            response.contact_area_m2,
        )
    }

    /// Preserve the complete Maugis phase, including an attractive noncontact
    /// state with no solid disk. The original model token is required because
    /// noncontacting equilibria deliberately carry no solid-contact response.
    /// Their cohesive annulus never generates shear, spin moment or heat.
    pub fn sliding_patch_maugis_equilibrium<'a>(
        &'a self,
        model: &'a crate::normal_patch::maugis::SourceMaugisSphereContact,
        response: &crate::normal_patch::maugis::continuation::MaugisEquilibrium,
    ) -> Result<SourceAdhesiveSlidingPatch<'a>, AdhesiveShearError> {
        if response.model_identity() != model.identity() {
            return Err(AdhesiveShearError::Invalid {
                reason: "Maugis equilibrium must retain its original normal model",
            });
        }
        match response.contact_response() {
            Some(contact) => self.sliding_patch_maugis(model.normal_state(), contact),
            None => self.sliding_patch(
                model.normal_state(),
                model.normal_state().identity(),
                0.0,
                0.0,
            ),
        }
    }

    fn sliding_patch<'a>(
        &'a self,
        normal: &'a BoundNormalInterfaceState,
        response_identity: ContentHash,
        radius: f64,
        area: f64,
    ) -> Result<SourceAdhesiveSlidingPatch<'a>, AdhesiveShearError> {
        if normal.identity() != self.normal_identity || response_identity != self.normal_identity {
            return Err(AdhesiveShearError::Invalid {
                reason: "spinning adhesive patch must retain its original normal and shear source state",
            });
        }
        let circular_area = PI * radius * radius;
        if !radius.is_finite()
            || radius < 0.0
            || !area.is_finite()
            || area < 0.0
            || !circular_area.is_finite()
            || (radius == 0.0) != (area == 0.0)
            || (area - circular_area).abs() > 32.0 * f64::EPSILON * area.max(circular_area)
        {
            return Err(AdhesiveShearError::Invalid {
                reason: "spinning adhesive junction requires the actual circular solid area",
            });
        }
        // SourceAdhesiveShear selects only the measured/authored constant
        // junction strengths. It never adds a pressure-dependent term here.
        let patch = |strength| {
            AdhesiveSlidingPatch::new(AdhesiveSlidingPatchParameters {
                contact_radius_m: radius,
                shear_strength_pa: strength,
            })
        };
        Ok(SourceAdhesiveSlidingPatch {
            normal,
            shear: self,
            patch: patch(self.law.kinetic_shear_strength_pa())?,
            static_patch: patch(self.law.static_shear_strength_pa())?,
            contact_area_m2: area,
        })
    }
}
