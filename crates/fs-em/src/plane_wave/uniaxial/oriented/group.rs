//! Original-source bulk group transport for the native vector carrier.
//! This differentiates each transparent spectral eigenbranch at fixed material
//! and optic-axis state, without resolving the Maxwell boundary system again.

use super::{
    Incidence, OrientedCrystalResponse, OrientedUniaxialInterface,
    PreparedOrientedCrystalField, ReflectionError, Result, kernel,
};
use fs_math::c64::C64;

/// The unchanged native boundary carrier and separate local group velocities.
/// Both original principal spectra and the authored global optic axis remain
/// bound through the original interface owner. Exact equal phase indices with
/// unequal group slopes refuse for an out-of-plane optic axis: the native
/// degenerate TE/TM carrier is then not a spectral eigenbasis. Principal-plane
/// crossings remain admitted.
#[derive(Clone, Copy)]
pub struct OrientedCrystalGroupResponse<'s, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    omega_rad_s: f64,
    boundary: OrientedCrystalResponse,
    principal_group_indices: [f64; 2],
    group_velocities_m_s: [Option<[f64; 3]>; 2],
}

impl<'s, 'm> OrientedCrystalGroupResponse<'s, 'm> {
    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Original complex fields, native eigenmodes and normal power fractions.
    pub const fn boundary_response(&self) -> &OrientedCrystalResponse {
        &self.boundary
    }

    /// Original [ordinary, extraordinary] n+omega*dn/domega.
    pub const fn principal_group_indices(&self) -> [f64; 2] {
        self.principal_group_indices
    }

    /// Separate ordinary/extraordinary grad_k(omega) in the original [x,y,z]
    /// frame [m/s], including out-of-plane walk-off despite conserved ky=0.
    /// Unexcited and evanescent modes return None, as do exactly zero spectral
    /// denominators. Critical modes retain exact zero normal group velocity.
    ///
    /// Signed anomalous slopes are retained. Alignment with Poynting flux
    /// requires positive dispersive energy in this transparent reciprocal
    /// model. These are modal bulk derivatives, not a coherent-superposition
    /// average, interface delay, broadband solution or signal-front speed.
    pub const fn group_velocities_m_s(&self) -> [Option<[f64; 3]>; 2] {
        self.group_velocities_m_s
    }
}

impl<'a> OrientedUniaxialInterface<'a> {
    /// Solve one native carrier and retain its original modal group transport.
    /// The original two-sided spectral domain and loss refusals remain in force.
    pub fn group_velocity_response(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<OrientedCrystalGroupResponse<'_, 'a>> {
        self.prepare_group_velocity(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            incident_tangential_e_v_m,
        )
    }

    /// Signed conserved kx/k0 counterpart, including exact critical coordinates.
    pub fn group_velocity_response_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<OrientedCrystalGroupResponse<'_, 'a>> {
        self.prepare_group_velocity(
            omega_rad_s,
            Incidence::TangentialIndex(tangential_index),
            incident_tangential_e_v_m,
        )
    }

    fn prepare_group_velocity(
        &self,
        omega_rad_s: f64,
        incidence: Incidence,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<OrientedCrystalGroupResponse<'_, 'a>> {
        let original = self.resolve_frequency(omega_rad_s)?;
        let boundary = original.solve(self.optic_axis, incidence, incident_tangential_e_v_m)?;
        prepare(
            self,
            original.omega,
            boundary,
            [original.no, original.ne],
        )
    }
}

impl<'s, 'm> PreparedOrientedCrystalField<'s, 'm> {
    /// Add original spectral derivatives to this exact prepared boundary.
    /// Materials are borrowed unchanged; no second boundary or root solve runs.
    pub fn group_velocity_response(&self) -> Result<OrientedCrystalGroupResponse<'s, 'm>> {
        let indices = self.source.crystal.phase_indices(self.omega_rad_s)?;
        prepare(self.source, self.omega_rad_s, self.boundary, indices)
    }
}

fn prepare<'s, 'm>(
    source: &'s OrientedUniaxialInterface<'m>,
    omega_rad_s: f64,
    boundary: OrientedCrystalResponse,
    indices: [f64; 2],
) -> Result<OrientedCrystalGroupResponse<'s, 'm>> {
    let principal_group_indices = source.crystal.group_indices(omega_rad_s)?;
    let group_velocities_m_s = kernel::mode_velocities(
        &boundary,
        omega_rad_s,
        indices,
        principal_group_indices,
        source.optic_axis,
    )
    .map_err(ReflectionError::Invalid)?;
    Ok(OrientedCrystalGroupResponse {
        source,
        omega_rad_s,
        boundary,
        principal_group_indices,
        group_velocities_m_s,
    })
}

