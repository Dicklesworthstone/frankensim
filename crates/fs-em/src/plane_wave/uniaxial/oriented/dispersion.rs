//! Separate normal-frequency propagation jets of the original vector modes.
//!
//! The native Maxwell boundary fields remain unchanged. Each outgoing branch
//! retains its own phase, delay, GDD and evanescent attenuation derivatives;
//! their coherent sum is not assigned a scalar pulse or group delay.

use super::super::UniaxialMode;
use super::super::dispersion::{
    MaterialDerivatives, kernel as dispersion_kernel, material_derivatives,
};
use super::{
    Incidence, OrientedCrystalResponse, OrientedUniaxialInterface, PreparedOrientedCrystalField,
    ReflectionError, ResolvedFrequency, Result,
};
use crate::plane_wave::dispersion::{
    BulkPropagation, IncidenceFrequencyConstraint, propagation_at_distance,
};
use fs_math::c64::C64;

#[cfg(test)]
mod tests;

/// Original arbitrary-axis carrier with independent local modal frequency jets.
/// The original incident law, both principal crystal laws and Cartesian optic
/// axis remain bound through the unchanged interface. Preparation solves the
/// native Maxwell boundary once and allocates no storage.
#[derive(Clone, Copy)]
pub struct OrientedCrystalFrequencyDispersion<'s, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    omega_rad_s: f64,
    constraint: IncidenceFrequencyConstraint,
    boundary: OrientedCrystalResponse,
    normal_derivatives: [[C64; 2]; 2],
}

impl<'s, 'm> OrientedCrystalFrequencyDispersion<'s, 'm> {
    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    pub const fn frequency_constraint(&self) -> IncidenceFrequencyConstraint {
        self.constraint
    }

    /// Unchanged native complex fields, phase roots and normal powers.
    pub const fn boundary_response(&self) -> &OrientedCrystalResponse {
        &self.boundary
    }

    /// Complex d(kz)/domega [s/m], ordinary then extraordinary. These are
    /// derivatives of the normal phase at a fixed transverse observation
    /// coordinate, not slowness along the generally different energy rays.
    /// Eigenmode derivatives remain defined when its excitation is zero.
    pub const fn transmitted_normal_wavenumber_derivatives_s_m(&self) -> [C64; 2] {
        [self.normal_derivatives[0][0], self.normal_derivatives[1][0]]
    }

    /// Complex d²(kz)/domega² [s²/m], ordinary then extraordinary. Fixed
    /// transverse wavenumber can introduce geometric normal-transfer GVD even
    /// for constant principal indices. Evanescent derivatives remain complex.
    pub const fn transmitted_normal_wavenumber_second_derivatives_s2_m(&self) -> [C64; 2] {
        [self.normal_derivatives[0][1], self.normal_derivatives[1][1]]
    }

    /// Propagation exp(-i*kz*z) of one original eigenmode through a nonnegative
    /// normal distance [m], including its local delay, GDD and attenuation jets.
    /// This does not differentiate the boundary amplitude or combine modal
    /// delays. In the admitted exactly isotropic limit the ordinary and
    /// extraordinary slots are the native deterministic TE and TM modes.
    ///
    /// Signed slopes are retained. They do not specify a signal-front speed,
    /// finite-crystal residence time or broadband pulse trajectory.
    pub fn transmitted_propagation(
        &self,
        mode: UniaxialMode,
        normal_distance_m: f64,
    ) -> Result<BulkPropagation> {
        let index = match mode {
            UniaxialMode::Ordinary => 0,
            UniaxialMode::Extraordinary => 1,
        };
        propagation_at_distance(
            self.boundary.transmitted[index].normal_wavenumber_rad_m,
            self.omega_rad_s,
            self.normal_derivatives[index],
            normal_distance_m,
        )
    }
}

impl<'a> OrientedUniaxialInterface<'a> {
    /// Solve one original vector carrier and differentiate both outgoing
    /// normal branches analytically under the explicit frequency constraint.
    /// FixedAngle requires the incident law's original lossless frequency jet;
    /// FixedTransverseWavenumber requires no incident frequency derivative.
    /// Both principal source domains, endpoints, cusps and loss refusals remain
    /// authoritative. No frequency differencing or second boundary solve runs.
    ///
    /// Either exactly critical branch refuses the complete two-mode jet query,
    /// even if that mode is unexcited; the original carrier remains available.
    /// At equal principal phase indices with an out-of-plane optic axis, both
    /// first and second local principal jets must coincide. Otherwise the
    /// native degenerate TE/TM basis is not an admitted spectral eigenbasis.
    /// Equal-phase principal-plane crossings remain supported.
    pub fn frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incident_tangential_e_v_m: [C64; 2],
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<OrientedCrystalFrequencyDispersion<'_, 'a>> {
        self.prepare_frequency_dispersion(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            incident_tangential_e_v_m,
            constraint,
        )
    }

    /// Signed carrier kx/k0 counterpart without an angle round-trip.
    /// FixedAngle follows the incident direction of this original carrier;
    /// FixedTransverseWavenumber holds its dimensional transverse wavenumber.
    pub fn frequency_dispersion_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        incident_tangential_e_v_m: [C64; 2],
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<OrientedCrystalFrequencyDispersion<'_, 'a>> {
        self.prepare_frequency_dispersion(
            omega_rad_s,
            Incidence::TangentialIndex(tangential_index),
            incident_tangential_e_v_m,
            constraint,
        )
    }

    fn prepare_frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence: Incidence,
        incident_tangential_e_v_m: [C64; 2],
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<OrientedCrystalFrequencyDispersion<'_, 'a>> {
        let original = self.resolve_frequency(omega_rad_s)?;
        let boundary = original.solve(self.optic_axis, incidence, incident_tangential_e_v_m)?;
        prepare(self, original, incidence, boundary, constraint)
    }
}

impl<'s, 'm> PreparedOrientedCrystalField<'s, 'm> {
    /// Add local frequency jets to this exact prepared vector boundary.
    /// The original cosine or signed tangent is retained during field
    /// preparation, including positive grazing cosines that cannot be recovered
    /// from rounded transverse wavenumbers. No Maxwell boundary or carrier-root
    /// solve is repeated; original material jets are evaluated at this frequency.
    ///
    /// The same critical-branch and spectral-eigenbasis refusals as the direct
    /// interface query apply. Coherent field sampling is unchanged.
    pub fn frequency_dispersion(
        &self,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<OrientedCrystalFrequencyDispersion<'s, 'm>> {
        let original = self.source.resolve_frequency(self.omega_rad_s)?;
        prepare(
            self.source,
            original,
            self.incidence,
            self.boundary,
            constraint,
        )
    }
}

fn prepare<'s, 'm>(
    source: &'s OrientedUniaxialInterface<'m>,
    original: ResolvedFrequency,
    incidence: Incidence,
    boundary: OrientedCrystalResponse,
    constraint: IncidenceFrequencyConstraint,
) -> Result<OrientedCrystalFrequencyDispersion<'s, 'm>> {
    let derivatives =
        material_derivatives(source.incident, source.crystal, original.omega, constraint)?;
    prepare_with_derivatives(
        source,
        original,
        incidence,
        boundary,
        constraint,
        &derivatives,
    )
}

/// Share the original material evaluation with the full boundary-jet owner.
pub(super) fn prepare_with_derivatives<'s, 'm>(
    source: &'s OrientedUniaxialInterface<'m>,
    original: ResolvedFrequency,
    incidence: Incidence,
    boundary: OrientedCrystalResponse,
    constraint: IncidenceFrequencyConstraint,
    derivatives: &MaterialDerivatives,
) -> Result<OrientedCrystalFrequencyDispersion<'s, 'm>> {
    if original.no == original.ne
        && source.optic_axis[1] != 0.0
        && (derivatives.relative_permittivity[1] != derivatives.relative_permittivity[2]
            || derivatives.principal_wavenumber[0] != derivatives.principal_wavenumber[1])
    {
        return Err(ReflectionError::Invalid(
            "equal principal phase indices with unequal frequency jets need a spectral eigenbasis for an out-of-plane optic axis",
        ));
    }
    let mut normal_derivatives = [[C64::ZERO; 2]; 2];
    for (index, mode) in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary]
        .into_iter()
        .enumerate()
    {
        let carrier = boundary.transmitted[index];
        normal_derivatives[index] =
            dispersion_kernel::normal_derivatives(dispersion_kernel::Input {
                omega: original.omega,
                indices: [original.ni, original.no, original.ne],
                relative_permittivity_derivatives: derivatives.relative_permittivity,
                principal_wavenumber_derivatives: derivatives.principal_wavenumber,
                incidence,
                axis: source.optic_axis,
                mode,
                kind: carrier.wave_kind,
                normal_wavenumber: carrier.normal_wavenumber_rad_m,
                constraint,
            })
            .map_err(ReflectionError::Invalid)?;
    }
    Ok(OrientedCrystalFrequencyDispersion {
        source,
        omega_rad_s: original.omega,
        constraint,
        boundary,
        normal_derivatives,
    })
}
