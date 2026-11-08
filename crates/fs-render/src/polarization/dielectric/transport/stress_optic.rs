//! Mechanical stress birefringence on the original homogeneous polarized beam.
//!
//! This is the isotropic linear extension of the sourced 7980 differential-index
//! coefficient. The caller supplies the original uniform Cauchy stress for this
//! segment; successive physical pieces may have different principal axes and
//! must be applied in traversal order. Stress is not averaged across pieces.
//!
//! Scalar CHARMS phase-index provenance and stress-optic calibration are
//! distinct. CHARMS supports 400--2600 nm; the stress source reports only 193,
//! 248 and 633 nm. Their present native intersection is the 633 nm source line.
//! A CHARMS temperature is not a certification of the stress coefficient at
//! that temperature: the stress experiment did not state a temperature law.
//!
//! The original ray, frame, medium and Beer advance remain authoritative.
//! Only relative polarization phase is added; common stress-induced index/group
//! shifts, double-ray splitting, walk-off and interface changes are not supplied.

use super::{BeamError, BeamMedium, PolarizedBeam, physical_roundoff};
use crate::dielectric::GlassProvenance;
use crate::polarization::Stokes;
use fs_material::Voigt;
use fs_material::silica_stress_optic::{Corning7980StressOptics, SilicaStressOpticError};
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StressOpticBeamError {
    /// The source coefficient belongs to 7980, not an arbitrary scalar glass.
    WrongMaterial,
    Source(SilicaStressOpticError),
    Beam(BeamError),
}
impl From<SilicaStressOpticError> for StressOpticBeamError {
    fn from(error: SilicaStressOpticError) -> Self {
        Self::Source(error)
    }
}
impl From<BeamError> for StressOpticBeamError {
    fn from(error: BeamError) -> Self {
        Self::Beam(error)
    }
}
impl core::fmt::Display for StressOpticBeamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native Corning 7980 stress optics: {self:?}")
    }
}
impl core::error::Error for StressOpticBeamError {}

impl PolarizedBeam {
    /// Advance through one uniformly stressed original Corning 7980 segment.
    ///
    /// Stress is the accepted mechanical owner's tension-positive tensor-shear
    /// Voigt vector in world coordinates, in Pa. This application adopts the
    /// isotropic linear constitutive extension; it does not certify a stress or
    /// temperature envelope beyond the source's stated laboratory conditions.
    ///
    /// Positive R makes the larger transverse principal stress the higher-n,
    /// slow axis. With exp(-i*omega*t) and V=2 Im(E_h conj(E_v)), horizontal
    /// tension turns +45-degree linear input toward POSITIVE V. The same rule
    /// holds in every right-handed native beam frame, including reverse rays.
    pub fn advance_corning7980_stress_optic(
        self,
        stress_pa: &Voigt,
        distance_m: f64,
    ) -> Result<Self, StressOpticBeamError> {
        let BeamMedium::Homogeneous(glass) = self.medium() else {
            return Err(BeamError::GradientFieldRequired.into());
        };
        if glass.provenance() != GlassProvenance::CharmsCorning7980V1 {
            return Err(StressOpticBeamError::WrongMaterial);
        }
        let law = Corning7980StressOptics::at_source_wavelength(self.wavelength_nm())?;
        let h = self.frame().horizontal();
        let v = self.frame().vertical();
        let retardance = law
            .project(stress_pa, [h.x, h.y, h.z], [v.x, v.y, v.z])?
            .retardance(distance_m)?;
        let advanced = self.advance(distance_m)?;
        if retardance.phase_rad == 0.0 {
            return Ok(advanced);
        }
        let phase = retardance.phase_rad.rem_euclid(core::f64::consts::TAU);
        // One half-angle evaluation also preserves 1-cos(phase) for weak
        // retardance where a direct subtraction would erase the mixing term.
        let (half_sine, half_cosine) = det::sin_cos(0.5 * phase);
        let one_minus_cosine = 2.0 * half_sine * half_sine;
        let sine = 2.0 * half_sine * half_cosine;
        let cosine = 1.0 - one_minus_cosine;
        let a = retardance.cos_double_slow_axis;
        let b = retardance.sin_double_slow_axis;
        let input = advanced.power();
        let parallel = a * input.q + b * input.u;
        // Rodrigues rotation of (Q,U,V), with no atan2, eigenvectors, matrix
        // allocation, common phase or per-mode power split. The sign follows
        // E_slow *= exp(+i*phase) under the native Jones/Stokes convention.
        let output = Stokes {
            i: input.i,
            q: cosine * input.q + sine * b * input.v + one_minus_cosine * a * parallel,
            u: cosine * input.u - sine * a * input.v + one_minus_cosine * b * parallel,
            v: cosine * input.v + sine * (a * input.u - b * input.q),
        };
        let power = physical_roundoff(output, input.i)?;
        Ok(Self { power, ..advanced })
    }
}

#[cfg(test)]
mod tests;
