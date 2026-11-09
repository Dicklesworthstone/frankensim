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
//! Local dispersive group transport reuses the native carrier and both original
//! principal spectral slopes. Biaxial crystals, absorbing anisotropy,
//! mode-coupled multilayers, spatial ray traversal and broadband pulse evolution
//! require other owners.
//! The existing isotropic interface/film/beam APIs retain their own domain.

use super::optical::UniaxialOpticalResponse;
use super::{DielectricResponse, Polarization, ReflectionError, positive};
use fs_material::{
    field_storage::linear_magnetic_permeability, state_point::ResolvedMaterialStatePoint,
};
use fs_math::det;

mod kernel;
pub use kernel::{UniaxialInterfaceResponse, UniaxialMode, UniaxialWaveKind};

mod dispersion;
pub use dispersion::UniaxialFrequencyDispersion;

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

/// Original transparent principal-plane carrier with its local bulk group
/// transport. The material/optic-axis owner is borrowed unchanged.
#[derive(Clone, Copy)]
pub struct UniaxialGroupResponse<'s, 'm> {
    source: &'s PrincipalPlaneUniaxialInterface<'m>,
    omega_rad_s: f64,
    carrier: UniaxialInterfaceResponse,
    principal_group_indices: [f64; 2],
    group_velocity_m_s: Option<[f64; 2]>,
}

impl<'s, 'm> UniaxialGroupResponse<'s, 'm> {
    pub const fn source(self) -> &'s PrincipalPlaneUniaxialInterface<'m> {
        self.source
    }
    pub const fn angular_frequency_rad_s(self) -> f64 {
        self.omega_rad_s
    }
    /// Unchanged native Snell/Fresnel solution, including phase and Poynting directions.
    pub const fn carrier(self) -> UniaxialInterfaceResponse {
        self.carrier
    }
    /// Original [ordinary, extraordinary] n+omega*dn/domega.
    pub const fn principal_group_indices(self) -> [f64; 2] {
        self.principal_group_indices
    }
    /// Signed local grad_k(omega) in the original [x,z] frame [m/s].
    /// None denotes an evanescent wave or a vanishing dispersive denominator,
    /// which has no finite group velocity. A critical bulk mode can have a
    /// finite tangent velocity without carrying normal transmitted power.
    ///
    /// Both principal spectral slopes enter an extraordinary wave. Phase and
    /// group directions differ through crystal walk-off. Alignment with the
    /// Poynting direction requires a positive dispersive denominator in this
    /// transparent reciprocal model; anomalous signed slopes remain signed.
    /// Finite-band index samples do not establish a causal energy reservoir,
    /// signal-front velocity, broadband packet evolution or an interface delay.
    pub const fn group_velocity_m_s(self) -> Option<[f64; 2]> {
        self.group_velocity_m_s
    }
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

    /// Add the original principal spectral slopes to the native carrier.
    /// The derivative is grad_k(omega) of the homogeneous crystal dispersion
    /// at frozen material/optic-axis state, not an incident-angle frequency
    /// sweep or a normal-transfer/Fresnel phase delay. Preparation allocates
    /// nothing and never substitutes a phase index for its group index.
    /// Original source endpoints, interpolation cusps and positive loss refuse.
    pub fn group_velocity_response(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<UniaxialGroupResponse<'_, 'a>, ReflectionError> {
        let (carrier, indices) = self.evaluate_with_indices(
            omega_rad_s,
            kernel::Incidence::Cosine(incidence_cosine),
            polarization,
        )?;
        self.prepare_group_response(omega_rad_s, carrier, [indices[1], indices[2]])
    }

    /// The same group transport with explicit signed kx/k0. Exact critical
    /// coordinates and negative incidence retain the native carrier convention.
    pub fn group_velocity_response_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        polarization: Polarization,
    ) -> Result<UniaxialGroupResponse<'_, 'a>, ReflectionError> {
        let (carrier, indices) = self.evaluate_with_indices(
            omega_rad_s,
            kernel::Incidence::TangentialIndex(tangential_index),
            polarization,
        )?;
        self.prepare_group_response(omega_rad_s, carrier, [indices[1], indices[2]])
    }

    fn prepare_group_response(
        &self,
        omega_rad_s: f64,
        carrier: UniaxialInterfaceResponse,
        indices: [f64; 2],
    ) -> Result<UniaxialGroupResponse<'_, 'a>, ReflectionError> {
        let principal_group_indices = self.crystal.group_indices(omega_rad_s)?;
        let group_velocity_m_s = kernel::group_velocity(
            carrier,
            omega_rad_s,
            indices,
            principal_group_indices,
            self.optic_axis,
        )
        .map_err(ReflectionError::Invalid)?;
        Ok(UniaxialGroupResponse {
            source: self,
            omega_rad_s,
            carrier,
            principal_group_indices,
            group_velocity_m_s,
        })
    }

    fn evaluate(
        &self,
        omega_rad_s: f64,
        incidence: kernel::Incidence,
        polarization: Polarization,
    ) -> Result<UniaxialInterfaceResponse, ReflectionError> {
        self.evaluate_with_indices(omega_rad_s, incidence, polarization)
            .map(|(carrier, _)| carrier)
    }

    // Share the native carrier and its original phase samples with group
    // transport; do not repeat the interface root solve.
    fn evaluate_with_indices(
        &self,
        omega_rad_s: f64,
        incidence: kernel::Incidence,
        polarization: Polarization,
    ) -> Result<(UniaxialInterfaceResponse, [f64; 3]), ReflectionError> {
        let omega = positive(omega_rad_s)?;
        let incident_epsilon = self.incident.lossless_permittivity(omega)?;
        // The original incident law supplies epsilon and sourced mu; no
        // nonmagnetic assumption or constant-index replacement is inserted.
        let incident_index = positive(
            (kernel::LIGHT_SPEED_M_S * det::sqrt(self.incident_mu_h_m))
                * det::sqrt(incident_epsilon),
        )?;
        let [ordinary_index, extraordinary_index] = self.crystal.phase_indices(omega)?;
        let carrier = kernel::solve(kernel::Input {
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
        .map_err(ReflectionError::Invalid)?;
        Ok((carrier, [incident_index, ordinary_index, extraordinary_index]))
    }
}
