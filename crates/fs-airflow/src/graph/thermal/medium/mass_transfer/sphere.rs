//! Explicit dilute mass-transfer correlations for an isolated spherical body.
//!
//! Frössling, *Über die Verdunstung fallender Tropfen*, Gerlands Beiträge zur
//! Geophysik 52, 170–216 (1938), Eq. (4.4), Table 4.2 and section 3.2.9:
//! <https://www.pdf-archive.com/2019/11/10/nils-frossling-ber-die-verdunstung-fallender-tropfen/nils-frossling-ber-die-verdunstung-fallender-tropfen.pdf>.
//! Its wind factor is `1+0.276 sqrt(Re) Sc^(1/3)`, hence diameter-based
//! `Sh=2+0.552 sqrt(Re) Sc^(1/3)`. The experiments span `Re=2..800` and
//! `D/nu=0.373..1.687`; this evaluator conservatively admits `Sc=0.6..2.6`.
//! The original derivation assumes weak vapor loading, negligible surface
//! blowing, uniform spherical surface, constant properties and quasi-steady
//! low-compressibility flow. Continuum, negligible buoyancy/deformation and
//! those physical assumptions remain caller authority.
//!
//! This is a direct mass-transfer model, not a heat-conductivity analogy. It
//! adds no transport source evidence. There is no automatic correlation,
//! zero-speed fallback, turbulent extension, or strong-blowing correction.

/// An explicitly selected model; there is deliberately no `Default`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ForcedSphereMassTransfer {
    Frossling1938,
}

/// Physical diameter and externally prescribed gas-to-body relative speed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForcedSphereGeometry {
    pub diameter_m: f64,
    pub relative_speed_m_s: f64,
}
impl ForcedSphereGeometry {
    pub fn validate(self) -> Result<(), ForcedSphereFilmError> {
        if !positive(self.diameter_m) || !positive(self.relative_speed_m_s) {
            return Err(ForcedSphereFilmError::Invalid(
                "forced sphere requires positive finite diameter and relative speed",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ForcedSphereFilmError {
    Invalid(&'static str),
    OutsideCorrelation { reynolds: f64, schmidt: f64 },
}
impl core::fmt::Display for ForcedSphereFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "forced spherical mass transfer: {self:?}")
    }
}
impl std::error::Error for ForcedSphereFilmError {}

/// Admitted transport groups and the exact inputs used in this evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForcedSphereFilm {
    model: ForcedSphereMassTransfer,
    geometry: ForcedSphereGeometry,
    kinematic_viscosity_m2_s: f64,
    diffusivity_m2_s: f64,
    reynolds: f64,
    schmidt: f64,
    sherwood_number: f64,
    coefficient_m_s: f64,
}
impl ForcedSphereFilm {
    pub const fn model(&self) -> ForcedSphereMassTransfer {
        self.model
    }
    pub const fn geometry(&self) -> ForcedSphereGeometry {
        self.geometry
    }
    pub const fn kinematic_viscosity_m2_s(&self) -> f64 {
        self.kinematic_viscosity_m2_s
    }
    pub const fn diffusivity_m2_s(&self) -> f64 {
        self.diffusivity_m2_s
    }
    pub const fn reynolds(&self) -> f64 {
        self.reynolds
    }
    pub const fn schmidt(&self) -> f64 {
        self.schmidt
    }
    pub const fn sherwood_number(&self) -> f64 {
        self.sherwood_number
    }
    pub const fn mass_transfer_coefficient_m_s(&self) -> f64 {
        self.coefficient_m_s
    }
}

/// Evaluate `Re=u*d/nu`, `Sc=nu/D_AB` and the selected diameter-based Sherwood.
///
/// Molecular diffusivity remains independent of momentum and heat transport.
/// Positive transport and prescribed speed are required; values outside the
/// selected correlation's admitted groups refuse without clamping or switching.
/// The numerical closure performs bounded scalar work and no allocation.
pub fn evaluate(
    model: ForcedSphereMassTransfer,
    geometry: ForcedSphereGeometry,
    kinematic_viscosity_m2_s: f64,
    diffusivity_m2_s: f64,
) -> Result<ForcedSphereFilm, ForcedSphereFilmError> {
    geometry.validate()?;
    if !positive(kinematic_viscosity_m2_s) || !positive(diffusivity_m2_s) {
        return Err(ForcedSphereFilmError::Invalid(
            "positive finite viscosity and molecular diffusivity required",
        ));
    }
    let reynolds =
        (geometry.relative_speed_m_s / kinematic_viscosity_m2_s) * geometry.diameter_m;
    let schmidt = kinematic_viscosity_m2_s / diffusivity_m2_s;
    if !positive(reynolds) || !positive(schmidt) {
        return Err(ForcedSphereFilmError::Invalid(
            "forced sphere transport groups are unrepresentable",
        ));
    }
    let sherwood_number = match model {
        ForcedSphereMassTransfer::Frossling1938 => {
            if !(2.0..=800.0).contains(&reynolds) || !(0.6..=2.6).contains(&schmidt) {
                return Err(ForcedSphereFilmError::OutsideCorrelation { reynolds, schmidt });
            }
            2.0 + 0.552 * fs_math::det::sqrt(reynolds)
                * fs_math::det::pow(schmidt, 1.0 / 3.0)
        }
    };
    let coefficient_m_s = (diffusivity_m2_s / geometry.diameter_m) * sherwood_number;
    if !positive(coefficient_m_s) {
        return Err(ForcedSphereFilmError::Invalid(
            "forced sphere mass-transfer coefficient is unrepresentable",
        ));
    }
    Ok(ForcedSphereFilm {
        model,
        geometry,
        kinematic_viscosity_m2_s,
        diffusivity_m2_s,
        reynolds,
        schmidt,
        sherwood_number,
        coefficient_m_s,
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
mod tests;
