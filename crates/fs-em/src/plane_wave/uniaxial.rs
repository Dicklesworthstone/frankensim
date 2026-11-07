//! Source-bound refraction into a transparent uniaxial crystal.
//!
//! The incident medium is isotropic and lossless at the queried frequency.
//! The crystal's explicitly authored optic axis lies in the incidence plane:
//! TE excites its ordinary mode and TM its extraordinary mode. Each response
//! retains signed complex boundary fields, normal power and the distinct
//! phase-normal and Poynting directions. A tilted extraordinary wave can carry
//! positive normal energy with a negative normal phase component.
//!
//! The interface is z=0, crystal z>0; [x,z] is an explicit right-handed
//! principal section with y perpendicular to it. Tangential electric axes
//! remain fixed for reflected and transmitted fields, as in the original
//! isotropic plane-wave owner. A negative tangential query changes incidence,
//! not the authored crystal axis. At normal incidence x is caller-selected.
//!
//! This is one homogeneous half-space interface. The [`oriented`] owner
//! additionally resolves out-of-plane axes and coupled vector polarization.
//! Biaxial crystals, absorbing anisotropy, mode-coupled multilayers, spatial
//! ray traversal and dispersive group-delay magnitudes require other owners.
//! The existing isotropic interface/film/beam APIs retain their own domain.

use super::optical::UniaxialOpticalResponse;
use super::{DielectricResponse, Polarization, ReflectionError, positive};
use fs_material::{
    field_storage::linear_magnetic_permeability, state_point::ResolvedMaterialStatePoint,
};
use fs_math::det;

mod kernel;
pub use kernel::{UniaxialInterfaceResponse, UniaxialMode, UniaxialWaveKind};

/// Coupled ordinary/extraordinary transmission for an arbitrary crystal axis.
pub mod oriented;

/// One original isotropic incident law, one two-axis crystal source, and one
/// fixed geometry-authored optic axis in the principal incidence section.
pub struct PrincipalPlaneUniaxialInterface<'a> {
    incident: DielectricResponse<'a>,
    crystal: UniaxialOpticalResponse<'a>,
    incident_mu_h_m: f64,
    optic_axis: [f64; 2],
}

impl<'a> PrincipalPlaneUniaxialInterface<'a> {
    /// Keep the original complete source owners and normalize a finite nonzero
    /// [tangential, normal] optic axis. Reversing the axis describes the same
    /// uniaxial crystal; no fast/slow ordering is inferred from its sign.
    /// Spectral losslessness and support are checked at every query.
    pub fn try_new(
        incident: DielectricResponse<'a>,
        crystal: UniaxialOpticalResponse<'a>,
        optic_axis: [f64; 2],
    ) -> Result<Self, ReflectionError> {
        let incident_mu_h_m = linear_magnetic_permeability(incident.material())?;
        let optic_axis = kernel::unit_axis(optic_axis).map_err(ReflectionError::Invalid)?;
        Ok(Self {
            incident,
            crystal,
            incident_mu_h_m,
            optic_axis,
        })
    }

    pub const fn incident_response(&self) -> DielectricResponse<'a> {
        self.incident
    }

    pub fn incident_material(&self) -> &'a ResolvedMaterialStatePoint {
        self.incident.material()
    }

    pub const fn crystal_response(&self) -> UniaxialOpticalResponse<'a> {
        self.crystal
    }

    /// Normalized [x,z] optic axis in the caller's actual principal section.
    pub const fn optic_axis(&self) -> [f64; 2] {
        self.optic_axis
    }

    /// Response for positive incident cosine in (0,1], with kx>=0.
    /// Exact grazing refuses because the incident normal power is zero.
    /// This path retains tiny positive grazing cosines without reconstructing
    /// them by subtracting a rounded sin² from one.
    pub fn response(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<UniaxialInterfaceResponse, ReflectionError> {
        self.evaluate(
            omega_rad_s,
            kernel::Incidence::Cosine(incidence_cosine),
            polarization,
        )
    }

    /// Response for an explicitly conserved signed tangential phase index
    /// kx/k0, where k0=omega/c0. Its magnitude must be strictly below the
    /// original incident phase index. This coordinate admits exact critical
    /// data and signed incidence without an angle round-trip.
    pub fn response_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        polarization: Polarization,
    ) -> Result<UniaxialInterfaceResponse, ReflectionError> {
        self.evaluate(
            omega_rad_s,
            kernel::Incidence::TangentialIndex(tangential_index),
            polarization,
        )
    }

    fn evaluate(
        &self,
        omega_rad_s: f64,
        incidence: kernel::Incidence,
        polarization: Polarization,
    ) -> Result<UniaxialInterfaceResponse, ReflectionError> {
        let omega = positive(omega_rad_s)?;
        let incident_epsilon = self.incident.lossless_permittivity(omega)?;
        // The original incident law supplies epsilon and sourced mu; no
        // nonmagnetic assumption or constant-index replacement is inserted.
        let incident_index = positive(
            (kernel::LIGHT_SPEED_M_S * det::sqrt(self.incident_mu_h_m))
                * det::sqrt(incident_epsilon),
        )?;
        let [ordinary_index, extraordinary_index] = self.crystal.phase_indices(omega)?;
        kernel::solve(kernel::Input {
            omega_rad_s: omega,
            incident_index,
            incident_mu_h_m: self.incident_mu_h_m,
            ordinary_index,
            extraordinary_index,
            crystal_mu_h_m: self.crystal.permeability_h_m(),
            optic_axis: self.optic_axis,
            incidence,
            mode: match polarization {
                Polarization::Te => UniaxialMode::Ordinary,
                Polarization::Tm => UniaxialMode::Extraordinary,
            },
        })
        .map_err(ReflectionError::Invalid)
    }
}
