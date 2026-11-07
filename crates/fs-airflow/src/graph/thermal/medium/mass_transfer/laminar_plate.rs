//! Admitted laminar binary-film transfer over a physical flat-plate strip.
//!
//! Seader and Henley, *Separation Process Principles*, 2nd ed., Wiley,
//! section 3.4, Eqs. (3-132), (3-137)--(3-138):
//! <https://catalogimages.wiley.com/images/db/pdf/0471464805.ch3.pdf>.
//! `Sh_L = 0.664 sqrt(Re_L) Sc^(1/3)`, `k_m = Sh_L D_AB/L`.
//! This implementation admits `0 < Re_L <= 5e5` and conservatively `Sc >= 0.6`
//! within the source's `Sc > 0.5` analogy. The source assumes an attached thin
//! laminar boundary layer, zero pressure gradient, uniform wall concentration,
//! constant properties and transfer weak enough to leave the velocity profile
//! undisturbed. Continuum validity and those physical assumptions are caller
//! authority. This numerical evaluator adds no material source evidence.

use fs_thermal::sublimation::SurfaceKinetics;

/// Streamwise positions retain the actual flow-facing physical leading edge.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarPlateGeometry {
    pub start_m: f64,
    pub end_m: f64,
    /// Positive relative free-stream speed parallel to the plate.
    pub speed_m_s: f64,
}
impl LaminarPlateGeometry {
    /// Reject zero flow: stagnant diffusion needs a separately authored geometry.
    pub fn validate(self) -> Result<(), LaminarPlateFilmError> {
        if !self.start_m.is_finite()
            || self.start_m < 0.0
            || !positive(self.end_m)
            || self.end_m <= self.start_m
            || !positive(self.speed_m_s)
        {
            return Err(LaminarPlateFilmError::Invalid(
                "ordered nonnegative plate positions and positive finite speed required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LaminarPlateFilmError {
    Invalid(&'static str),
    OutsideCorrelation { reynolds: f64, schmidt: f64 },
}
impl core::fmt::Display for LaminarPlateFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "laminar plate gas film: {self:?}")
    }
}
impl std::error::Error for LaminarPlateFilmError {}

/// Checked coefficient and the exact physical inputs/groups used to compute it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarPlateFilm {
    geometry: LaminarPlateGeometry,
    kinematic_viscosity_m2_s: f64,
    diffusivity_m2_s: f64,
    reynolds: f64,
    schmidt: f64,
    downstream_mean_sherwood: f64,
    coefficient_m_s: f64,
}
impl LaminarPlateFilm {
    pub const fn geometry(&self) -> LaminarPlateGeometry {
        self.geometry
    }
    pub const fn kinematic_viscosity_m2_s(&self) -> f64 {
        self.kinematic_viscosity_m2_s
    }
    pub const fn diffusivity_m2_s(&self) -> f64 {
        self.diffusivity_m2_s
    }
    /// Based on the downstream distance from the original leading edge.
    pub const fn reynolds(&self) -> f64 {
        self.reynolds
    }
    pub const fn schmidt(&self) -> f64 {
        self.schmidt
    }
    /// Mean over the whole `[0, end]`, independent of the numerical strip length.
    pub const fn downstream_mean_sherwood(&self) -> f64 {
        self.downstream_mean_sherwood
    }
    /// Mean over the requested `[start, end]` strip [m/s].
    pub const fn mass_transfer_coefficient_m_s(&self) -> f64 {
        self.coefficient_m_s
    }

    /// Original thermal-owner handoff, with independently authored accommodation.
    /// The receiving owner's dilute film and gas-temperature pressure convention
    /// remain unchanged; this frozen coefficient carries no species identity.
    pub fn surface_kinetics(
        &self,
        accommodation_coefficient: f64,
    ) -> Result<SurfaceKinetics, LaminarPlateFilmError> {
        if !accommodation_coefficient.is_finite()
            || !(0.0..=1.0).contains(&accommodation_coefficient)
        {
            return Err(LaminarPlateFilmError::Invalid(
                "accommodation probability must be in [0,1]",
            ));
        }
        Ok(SurfaceKinetics {
            accommodation_coefficient,
            mass_transfer_coefficient_m_s: Some(self.coefficient_m_s),
        })
    }
}

/// Evaluate the original plate correlation with explicitly independent `nu` and `D_AB`.
///
/// The strip coefficient is the rationalized difference of cumulative transfer:
/// `k_[a,b] = k_[0,b] / (1 + sqrt(a/b))`. It remains well conditioned for narrow
/// strips and integrates back to the original whole-plate mean. Reversing flow
/// requires positions measured from the opposite physical edge. Zero speed,
/// nonpositive transport and unrepresentable positive coefficients refuse.
pub fn evaluate(
    geometry: LaminarPlateGeometry,
    kinematic_viscosity_m2_s: f64,
    diffusivity_m2_s: f64,
) -> Result<LaminarPlateFilm, LaminarPlateFilmError> {
    geometry.validate()?;
    if !positive(kinematic_viscosity_m2_s) || !positive(diffusivity_m2_s) {
        return Err(LaminarPlateFilmError::Invalid(
            "positive finite viscosity and molecular diffusivity required",
        ));
    }
    let reynolds = (geometry.speed_m_s / kinematic_viscosity_m2_s) * geometry.end_m;
    let schmidt = kinematic_viscosity_m2_s / diffusivity_m2_s;
    if !positive(reynolds) || !positive(schmidt) {
        return Err(LaminarPlateFilmError::Invalid(
            "molecular transport group is unrepresentable",
        ));
    }
    if reynolds > 5.0e5 || schmidt < 0.6 {
        return Err(LaminarPlateFilmError::OutsideCorrelation { reynolds, schmidt });
    }
    let downstream_mean_sherwood =
        0.664 * fs_math::det::sqrt(reynolds) * fs_math::det::pow(schmidt, 1.0 / 3.0);
    let mean_coefficient = downstream_mean_sherwood * (diffusivity_m2_s / geometry.end_m);
    let coefficient_m_s =
        mean_coefficient / (1.0 + fs_math::det::sqrt(geometry.start_m / geometry.end_m));
    if !positive(downstream_mean_sherwood) || !positive(coefficient_m_s) {
        return Err(LaminarPlateFilmError::Invalid(
            "positive gas-film coefficient is unrepresentable",
        ));
    }
    Ok(LaminarPlateFilm {
        geometry,
        kinematic_viscosity_m2_s,
        diffusivity_m2_s,
        reynolds,
        schmidt,
        downstream_mean_sherwood,
        coefficient_m_s,
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
#[path = "laminar_plate_tests.rs"]
mod tests;
