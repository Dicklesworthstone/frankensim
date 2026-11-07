//! Passive anisotropic absorption from the original two complex index axes.
//!
//! This extends the same planar interface and global optic axis used by the
//! transparent owner. It solves complex Maxwell dispersion, not transparent
//! refraction followed by an empirical attenuation multiplier. For
//! exp(+i omega t-i k.r), admitted transmitted modes have Im(kz)<=0.
//! Polarizations interfere before Poynting flux and dielectric heating are
//! evaluated. Scalar real permeability carries no magnetic loss here.
//!
//! This is one homogeneous, reciprocal, uniaxial half-space with a lossless
//! incident port. There is no absorbing incident port, finite lossy stack,
//! gain, biaxial/Jordan model, causal pulse reservoir or thermal time step.
//! Spatial energy-flow reference: Passler, Jeannin and Paarmann (2020),
//! <https://arxiv.org/abs/2002.03832>. The time/phase signs here remain the
//! native exp(+i omega t-i k.r) convention.

use super::{C0, Incidence, OrientedUniaxialInterface, ReflectionError, Result, positive};
use fs_math::{c64::C64, det};

mod kernel;
pub use kernel::{AbsorbingCrystalField, AbsorbingCrystalMode, AbsorbingCrystalResponse};

/// One original material/frequency/illumination snapshot. Point queries have
/// fixed work, allocate nothing and never resolve source laws or solve the
/// boundary again. The borrowed interface retains both original spectra.
pub struct PreparedAbsorbingCrystalField<'s, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    boundary: AbsorbingCrystalResponse,
}

impl<'s, 'm> PreparedAbsorbingCrystalField<'s, 'm> {
    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.boundary.omega_rad_s()
    }

    pub const fn boundary_response(&self) -> &AbsorbingCrystalResponse {
        &self.boundary
    }

    /// Physical Cartesian point [m], with the interface at z=0. The original
    /// half-space selection and peak-phasor convention are shared with the
    /// transparent field owner. Positive heat is evaluated from the original
    /// loss tensor and the selected coherent electric field, in W/m³.
    /// Choosing one eigenmode returns that mode's heat; the heat of the total
    /// field generally differs from the sum of the two separate mode heats.
    pub fn field(
        &self,
        side: super::OrientedCrystalFieldSide,
        position_m: [f64; 3],
    ) -> Result<AbsorbingCrystalField> {
        self.boundary
            .field(side, position_m)
            .map_err(ReflectionError::Invalid)
    }
}

impl<'m> OrientedUniaxialInterface<'m> {
    /// Prepare a passive absorbing crystal from both original complex
    /// principal-index laws. The original transparent query still refuses
    /// loss; this query explicitly admits it. A positive cosine retains tiny
    /// grazing inputs without reconstructing cosine from the tangential index.
    pub fn prepare_absorbing_field(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingCrystalField<'_, 'm>> {
        self.prepare_absorbing(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            incident_tangential_e_v_m,
        )
    }

    /// Signed conserved kx/k0 counterpart. Its magnitude must lie strictly
    /// inside the original lossless incident light cone. Source support and
    /// interpolation-gap limits are checked before the fixed Maxwell solve.
    pub fn prepare_absorbing_field_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingCrystalField<'_, 'm>> {
        self.prepare_absorbing(
            omega_rad_s,
            Incidence::TangentialIndex(tangential_index),
            incident_tangential_e_v_m,
        )
    }

    fn prepare_absorbing(
        &self,
        omega_rad_s: f64,
        incidence: Incidence,
        electric: [C64; 2],
    ) -> Result<PreparedAbsorbingCrystalField<'_, 'm>> {
        let omega = positive(omega_rad_s)?;
        let epsilon = self.incident.lossless_permittivity(omega)?;
        let ni = positive((C0 * det::sqrt(self.incident_mu_h_m)) * det::sqrt(epsilon))?;
        let indices = self.crystal.complex_principal_indices(omega)?;
        let permittivities = self.crystal.principal_permittivities_f_m(omega)?;
        let boundary = kernel::solve(kernel::Input {
            omega,
            ni,
            indices,
            mui: self.incident_mu_h_m,
            mut_: self.crystal.permeability_h_m(),
            axis: self.optic_axis,
            incidence,
            electric,
            loss_f_m: permittivities.map(|value| -value.im),
        })
        .map_err(ReflectionError::Invalid)?;
        Ok(PreparedAbsorbingCrystalField {
            source: self,
            boundary,
        })
    }
}
