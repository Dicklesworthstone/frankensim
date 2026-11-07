//! A finite transparent uniaxial slab with coherent internal reflections.
//!
//! The original arbitrarily oriented crystal occupies 0<z<d between two
//! isotropic lossless half-spaces. Tangential E and H are matched at BOTH
//! faces, retaining polarization coupling, etalon resonances and evanescent
//! tunneling. Forward modes are anchored at the entrance, backward modes at
//! the exit; no exponentially growing propagation factor is constructed.
//!
//! Exact coalesced crystal critical modes need a generalized (Jordan) basis
//! and explicitly refuse here. Numerically unresolved boundary systems also
//! refuse; positive propagating output flux must remain representable.
//! The existing half-space owner retains its original critical-mode
//! support. This owner supplies no absorbing anisotropy, multilayer crystal
//! stack, spatially varying axis or causal time-domain reservoir.

use super::{DielectricResponse, Incidence, OrientedUniaxialInterface, ReflectionError, positive};
use fs_material::field_storage::linear_magnetic_permeability;
use fs_math::{c64::C64, det};

mod kernel;
pub use kernel::{SlabRegion, SlabWave, UniaxialSlabResponse};

mod spectrum;
pub use spectrum::PreparedUniaxialSlabSpectrum;

const C0: f64 = 299_792_458.0;
type Result<T> = core::result::Result<T, ReflectionError>;

/// One original entrance/crystal owner, original exit law and physical depth.
pub struct OrientedUniaxialSlab<'s, 'm> {
    entrance: &'s OrientedUniaxialInterface<'m>,
    exit: DielectricResponse<'m>,
    exit_mu_h_m: f64,
    thickness_m: f64,
}

impl<'s, 'm> OrientedUniaxialSlab<'s, 'm> {
    /// Retain the exact spectral laws, permeability and fixed global optic
    /// axis. Thickness must be positive, finite and representable.
    pub fn try_new(
        entrance: &'s OrientedUniaxialInterface<'m>,
        exit: DielectricResponse<'m>,
        thickness_m: f64,
    ) -> Result<Self> {
        Ok(Self {
            entrance,
            exit,
            exit_mu_h_m: linear_magnetic_permeability(exit.material())?,
            thickness_m: positive(thickness_m)?,
        })
    }

    pub const fn entrance_interface(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.entrance
    }
    pub const fn exit_response(&self) -> DielectricResponse<'m> {
        self.exit
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// One x-z incident plane with positive cosine and original complex Ex/Ey
    /// peak fields [V/m]. The returned field snapshot requires no further
    /// material evaluation or boundary solve at different observation points.
    pub fn response(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<UniaxialSlabResponse> {
        kernel::solve(
            self.resolve_frequency(omega_rad_s)?,
            Incidence::Cosine(incidence_cosine),
            incident_tangential_e_v_m,
        )
        .map_err(ReflectionError::Invalid)
    }

    /// Signed conserved kx/k0 without an angle round-trip. The incident port
    /// must propagate; the original exit may propagate, be critical or decay.
    pub fn response_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_phase_index: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<UniaxialSlabResponse> {
        kernel::solve(
            self.resolve_frequency(omega_rad_s)?,
            Incidence::TangentialIndex(tangential_phase_index),
            incident_tangential_e_v_m,
        )
        .map_err(ReflectionError::Invalid)
    }

    fn resolve_frequency(&self, omega_rad_s: f64) -> Result<kernel::Input> {
        let crystal = self.entrance.resolve_frequency(omega_rad_s)?;
        let epsilon = self.exit.lossless_permittivity(crystal.omega)?;
        let exit_index = positive((C0 * det::sqrt(self.exit_mu_h_m)) * det::sqrt(epsilon))?;
        Ok(kernel::Input {
            crystal,
            exit_index,
            exit_mu_h_m: self.exit_mu_h_m,
            thickness_m: self.thickness_m,
            axis: self.entrance.optic_axis(),
        })
    }
}
