//! Vector Maxwell refraction into an arbitrarily oriented transparent crystal.
//!
//! Coordinates and time convention match the principal-plane owner: z>0 is
//! crystal, the incident wavevector lies in x-z, and exp(+i omega t-i k.r).
//! The source supplies both principal indices and scalar permeability. An
//! out-of-plane optic axis couples the two incident polarizations; solving
//! both tangential E and H conditions retains this coupling and both outgoing
//! eigenmodes. Eigenmode phase and energy directions need not coincide.
//! [`PreparedOrientedCrystalField`] propagates the solved fields coherently
//! through their actual half-spaces, retaining both modal phase and decay.
//!
//! This is one planar, lossless, homogeneous half-space. It does not supply
//! absorbing/biaxial crystals, multilayers, transient fields or ray traversal.
//! [`spectrum`] prepares coherent incident directions against this same owner
//! and propagates both eigenmodes to physical points in either half-space.
//! Maxwell's anisotropic dispersion and eigenvectors: MIT 6.974 §2.8,
//! <https://ocw.mit.edu/courses/6-974-fundamentals-of-photonics-quantum-electronics-spring-2006/8e54d1625e9e71eb5b4e4998e6967e4e_chapter2.pdf>.

use super::kernel::Incidence;
use super::{DielectricResponse, ReflectionError, UniaxialOpticalResponse, positive};
use fs_material::field_storage::linear_magnetic_permeability;
use fs_math::{c64::C64, det};

mod kernel;
use kernel::Input;
pub use kernel::{
    OrientedCrystalField, OrientedCrystalFieldSide, OrientedCrystalMode, OrientedCrystalResponse,
};

/// Source-bound coherent fields for sparse, authored incident directions.
pub mod spectrum;

/// Coherent scattering and internal fields of one finite transparent crystal.
pub mod slab;

/// Coherent stacks of finite crystals with independent original global axes.
pub mod stack;

const C0: f64 = 299_792_458.0;
type Result<T> = core::result::Result<T, ReflectionError>;

/// One evaluation of the original spectral laws, shared by all prepared modes.
#[derive(Clone, Copy)]
struct ResolvedFrequency {
    omega: f64,
    ni: f64,
    no: f64,
    ne: f64,
    mui: f64,
    mut_: f64,
}

impl ResolvedFrequency {
    fn solve(
        self,
        axis: [f64; 3],
        incidence: Incidence,
        electric: [C64; 2],
    ) -> Result<OrientedCrystalResponse> {
        kernel::solve(Input {
            omega: self.omega,
            ni: self.ni,
            no: self.no,
            ne: self.ne,
            mui: self.mui,
            mut_: self.mut_,
            axis,
            incidence,
            electric,
        })
        .map_err(ReflectionError::Invalid)
    }
}

/// Original source laws and an explicitly authored Cartesian optic axis.
pub struct OrientedUniaxialInterface<'a> {
    incident: DielectricResponse<'a>,
    crystal: UniaxialOpticalResponse<'a>,
    incident_mu_h_m: f64,
    optic_axis: [f64; 3],
}

/// One source-bound frequency and illumination, solved once at z=0.
/// Repeated spatial queries retain both transmitted modes and allocate nothing;
/// they do not resolve materials or solve the Maxwell boundary system again.
pub struct PreparedOrientedCrystalField<'s, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    omega_rad_s: f64,
    boundary: OrientedCrystalResponse,
}

impl<'s, 'm> PreparedOrientedCrystalField<'s, 'm> {
    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    pub const fn boundary_response(&self) -> &OrientedCrystalResponse {
        &self.boundary
    }

    /// Sample authored Cartesian [x,y,z] in metres, with z=0 at the interface.
    /// Incident-side queries require z<=0; crystal queries require z>=0.
    /// Either boundary trace may be requested at z=0. The plane wave has ky=0,
    /// so fields are independent of y even when its energy walks off in y.
    /// The authored incident complex amplitude is anchored at [0,0,0].
    ///
    /// These are complex peak phasors with the original exp(+i omega t)
    /// convention. Components interfere before their signed cycle-mean
    /// Poynting vector is formed; evanescent decay does not represent heat.
    pub fn field(
        &self,
        side: OrientedCrystalFieldSide,
        position_m: [f64; 3],
    ) -> Result<OrientedCrystalField> {
        kernel::field(&self.boundary, side, position_m).map_err(ReflectionError::Invalid)
    }
}

impl<'a> OrientedUniaxialInterface<'a> {
    pub fn try_new(
        incident: DielectricResponse<'a>,
        crystal: UniaxialOpticalResponse<'a>,
        optic_axis: [f64; 3],
    ) -> Result<Self> {
        let incident_mu_h_m = linear_magnetic_permeability(incident.material())?;
        Ok(Self {
            incident,
            crystal,
            incident_mu_h_m,
            optic_axis: kernel::unit(optic_axis).map_err(ReflectionError::Invalid)?,
        })
    }

    pub const fn incident_response(&self) -> DielectricResponse<'a> {
        self.incident
    }
    pub const fn crystal_response(&self) -> UniaxialOpticalResponse<'a> {
        self.crystal
    }
    pub const fn optic_axis(&self) -> [f64; 3] {
        self.optic_axis
    }

    /// Positive incident cosine, kx>=0, and complex tangential [Ex,Ey] [V/m].
    /// The input retains arbitrary linear/elliptical polarization and phase.
    /// Zero illumination, exact grazing and unsupported spectral data refuse.
    /// Evaluation is fixed-size and allocation-free.
    pub fn response(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<OrientedCrystalResponse> {
        self.evaluate(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            incident_tangential_e_v_m,
        )
    }

    /// Signed conserved kx/k0 and complex tangential [Ex,Ey] [V/m].
    /// This admits exact critical indices without an angle round-trip. Its
    /// magnitude must remain below the original incident phase index.
    pub fn response_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<OrientedCrystalResponse> {
        self.evaluate(
            omega_rad_s,
            Incidence::TangentialIndex(tangential_index),
            incident_tangential_e_v_m,
        )
    }

    /// Resolve the original spectral support and prepare one coherent field
    /// for repeated physical point queries. Positive cosine implies kx>=0.
    pub fn prepare_field(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<PreparedOrientedCrystalField<'_, 'a>> {
        Ok(PreparedOrientedCrystalField {
            source: self,
            omega_rad_s,
            boundary: self.response(omega_rad_s, incidence_cosine, incident_tangential_e_v_m)?,
        })
    }

    /// Conserved signed kx/k0 counterpart, including exact critical inputs.
    /// The preparation retains the same original source-law refusals.
    pub fn prepare_field_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<PreparedOrientedCrystalField<'_, 'a>> {
        Ok(PreparedOrientedCrystalField {
            source: self,
            omega_rad_s,
            boundary: self.response_at_tangential_index(
                omega_rad_s,
                tangential_index,
                incident_tangential_e_v_m,
            )?,
        })
    }

    fn evaluate(
        &self,
        omega_rad_s: f64,
        incidence: Incidence,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<OrientedCrystalResponse> {
        self.resolve_frequency(omega_rad_s)?.solve(
            self.optic_axis,
            incidence,
            incident_tangential_e_v_m,
        )
    }

    fn resolve_frequency(&self, omega_rad_s: f64) -> Result<ResolvedFrequency> {
        let omega = positive(omega_rad_s)?;
        let epsilon = self.incident.lossless_permittivity(omega)?;
        let ni = positive((C0 * det::sqrt(self.incident_mu_h_m)) * det::sqrt(epsilon))?;
        let [no, ne] = self.crystal.phase_indices(omega)?;
        Ok(ResolvedFrequency {
            omega,
            ni,
            no,
            ne,
            mui: self.incident_mu_h_m,
            mut_: self.crystal.permeability_h_m(),
        })
    }
}
