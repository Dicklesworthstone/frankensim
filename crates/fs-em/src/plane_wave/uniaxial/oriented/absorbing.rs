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
pub use kernel::{
    AbsorbingCrystalField, AbsorbingCrystalIlluminationPower, AbsorbingCrystalMode,
    AbsorbingCrystalResponse, AbsorbingCrystalStatistics,
};

mod illumination;
pub use illumination::{
    AbsorbingCrystalIllumination, AbsorbingCrystalIlluminationSlabAbsorption,
    CrystalIlluminationError,
};

mod coherent;
pub use coherent::AbsorbingCrystalCoherentIllumination;

/// One original material/frequency/illumination snapshot. Point queries have
/// fixed work, allocate nothing and never resolve source laws or solve the
/// boundary again. The borrowed interface retains both original spectra.
pub struct PreparedAbsorbingCrystalField<'s, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    boundary: AbsorbingCrystalResponse,
}

/// Cycle-mean dielectric heat in a finite prism of the original crystal
/// half-space. The face area is measured in the interface plane and depths
/// are the original Cartesian z coordinates. This control volume introduces
/// no exit interface or reflected wave. Its material and illumination remain
/// frozen; a thermal consumer must declare the time over which it uses this
/// power and bind the volume to the same original material snapshot.
pub struct AbsorbingCrystalSlabAbsorption<'w, 's, 'm> {
    prepared_field: &'w PreparedAbsorbingCrystalField<'s, 'm>,
    area_m2: f64,
    start_depth_m: f64,
    end_depth_m: f64,
    volume_m3: f64,
    absorbed_power_w: f64,
}

impl<'w, 's, 'm> AbsorbingCrystalSlabAbsorption<'w, 's, 'm> {
    pub const fn prepared_field(&self) -> &'w PreparedAbsorbingCrystalField<'s, 'm> {
        self.prepared_field
    }

    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.prepared_field.source()
    }

    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    pub const fn start_depth_m(&self) -> f64 {
        self.start_depth_m
    }

    pub const fn end_depth_m(&self) -> f64 {
        self.end_depth_m
    }

    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }

    /// Integral of omega/2 E*·[-Im(epsilon)]·E over the declared prism [W].
    /// Ordinary and extraordinary fields interfere before heat is integrated.
    pub const fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
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

    /// Analytically integrate the coherent transmitted dielectric heat over
    /// a finite positive area and `0 <= start_depth_m < end_depth_m`.
    /// The two modal self terms and their interference use the original loss
    /// tensor. Fixed, allocation-free exponential integrals retain weak loss,
    /// thin intervals and opaque tails without subtracting nearly equal port
    /// powers, depth quadrature, material queries or another boundary solve.
    pub fn slab_absorption(
        &self,
        area_m2: f64,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<AbsorbingCrystalSlabAbsorption<'_, 's, 'm>> {
        let volume_m3 = kernel::slab_volume_m3(area_m2, start_depth_m, end_depth_m)
            .map_err(ReflectionError::Invalid)?;
        let absorbed_power_w = self
            .boundary
            .slab_absorbed_power_w(area_m2, start_depth_m, end_depth_m)
            .map_err(ReflectionError::Invalid)?;
        Ok(AbsorbingCrystalSlabAbsorption {
            prepared_field: self,
            area_m2,
            start_depth_m,
            end_depth_m,
            volume_m3,
            absorbed_power_w,
        })
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
            [1.0, 0.0],
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
            [1.0, 0.0],
            incident_tangential_e_v_m,
        )
    }

    /// Prepare any transverse incidence azimuth on this same original
    /// interface. The finite nonzero [x,y] incidence-plane axis is normalized;
    /// the positive tangential wavevector points along it. Input amplitudes
    /// are GLOBAL Cartesian [Ex,Ey] peak fields, not local p/s amplitudes.
    /// The original optic axis remains fixed in the global material frame.
    /// At normal incidence the plane axis has no physical effect.
    pub fn prepare_absorbing_field_in_plane(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incidence_plane_axis: [f64; 2],
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingCrystalField<'_, 'm>> {
        self.prepare_absorbing(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            incidence_plane_axis,
            incident_tangential_e_v_m,
        )
    }

    /// Signed k_t/k0 along an authored transverse incidence-plane axis.
    /// Global [Ex,Ey], original source support and incident light-cone
    /// admission are identical to `prepare_absorbing_field_in_plane`.
    /// Reversing both the axis and signed index preserves the same wave.
    pub fn prepare_absorbing_field_at_tangential_index_in_plane(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        incidence_plane_axis: [f64; 2],
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingCrystalField<'_, 'm>> {
        self.prepare_absorbing(
            omega_rad_s,
            Incidence::TangentialIndex(tangential_index),
            incidence_plane_axis,
            incident_tangential_e_v_m,
        )
    }

    fn prepare_absorbing(
        &self,
        omega_rad_s: f64,
        incidence: Incidence,
        incidence_plane_axis: [f64; 2],
        electric: [C64; 2],
    ) -> Result<PreparedAbsorbingCrystalField<'_, 'm>> {
        let omega = positive(omega_rad_s)?;
        let epsilon = self.incident.lossless_permittivity(omega)?;
        let ni = positive((C0 * det::sqrt(self.incident_mu_h_m)) * det::sqrt(epsilon))?;
        let indices = self.crystal.complex_principal_indices(omega)?;
        let permittivities = self.crystal.principal_permittivities_f_m(omega)?;
        let boundary = kernel::solve_in_plane(
            kernel::Input {
                omega,
                ni,
                indices,
                mui: self.incident_mu_h_m,
                mut_: self.crystal.permeability_h_m(),
                axis: self.optic_axis,
                incidence,
                electric,
                loss_f_m: permittivities.map(|value| -value.im),
            },
            incidence_plane_axis,
        )
        .map_err(ReflectionError::Invalid)?;
        Ok(PreparedAbsorbingCrystalField {
            source: self,
            boundary,
        })
    }
}
