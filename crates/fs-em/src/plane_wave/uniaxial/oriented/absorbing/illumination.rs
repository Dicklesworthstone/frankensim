//! Partially polarized and independent spectral illumination of one crystal.
//!
//! Every borrowed field is one explicitly mutually uncorrelated peak-field
//! factor, including its square-root statistical weight. Its ordinary and
//! extraordinary waves remain coherent. Factors add as second moments and
//! material heat, never as phasors. Their original frequencies, directions,
//! spectra, loss tensors and common global optic axis remain authoritative.
//! These are finite spectral-line powers, not PSD/Hz or a bandwidth inferred
//! from unknown phase. Deterministic beat envelopes and cross-frequency
//! coherence require a different source model; no random-phase sampling is
//! performed. The source and statistical illumination remain frozen.
//!
//! Positive electromagnetic field-factor representation: Tervo, Setälä and
//! Friberg (2004), <https://doi.org/10.1364/JOSAA.21.002205>.

use super::super::OrientedCrystalFieldSide;
use super::{
    AbsorbingCrystalIlluminationPower, AbsorbingCrystalStatistics, OrientedUniaxialInterface,
    PreparedAbsorbingCrystalField, ReflectionError, kernel,
};

/// Cancellation remains distinguishable from original source/Maxwell refusal
/// so a receiving thermal owner can preserve its atomic cancellation path.
#[derive(Debug)]
pub enum CrystalIlluminationError {
    Cancelled,
    Electromagnetic(ReflectionError),
}
impl From<ReflectionError> for CrystalIlluminationError {
    fn from(error: ReflectionError) -> Self {
        Self::Electromagnetic(error)
    }
}
impl From<kernel::IlluminationError> for CrystalIlluminationError {
    fn from(error: kernel::IlluminationError) -> Self {
        match error {
            kernel::IlluminationError::Cancelled => Self::Cancelled,
            kernel::IlluminationError::Electromagnetic(message) => {
                Self::Electromagnetic(ReflectionError::Invalid(message))
            }
        }
    }
}
impl core::fmt::Display for CrystalIlluminationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => f.write_str("crystal illumination cancelled"),
            Self::Electromagnetic(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for CrystalIlluminationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cancelled => None,
            Self::Electromagnetic(error) => Some(error),
        }
    }
}
type Result<T> = core::result::Result<T, CrystalIlluminationError>;

/// A borrowed, source-bound statistical illumination. No field array is
/// copied or allocated. Preparation and each observation take O(factors)
/// work; every factor retains the fixed two-mode Maxwell/absorption query.
pub struct AbsorbingCrystalIllumination<'p, 's, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    factors: &'p [PreparedAbsorbingCrystalField<'s, 'm>],
    power: AbsorbingCrystalIlluminationPower,
}

/// Integrated ensemble dielectric heat in the original half-space volume.
/// This borrows the complete illumination and introduces no exit interface.
pub struct AbsorbingCrystalIlluminationSlabAbsorption<'i, 'p, 's, 'm> {
    illumination: &'i AbsorbingCrystalIllumination<'p, 's, 'm>,
    area_m2: f64,
    start_depth_m: f64,
    end_depth_m: f64,
    volume_m3: f64,
    absorbed_power_w: f64,
}

impl<'p, 's, 'm> AbsorbingCrystalIllumination<'p, 's, 'm> {
    /// Bind explicit uncorrelated factors to the exact original interface.
    /// Frequencies and signed incident directions may differ; each factor
    /// already passed its original spectral and light-cone admission. The
    /// source contributes no implicit field. Empty factors mean darkness.
    ///
    /// For unpolarized light at one incidence cosine c, equal transverse
    /// p/s peak amplitudes a use tangential factors [c*a,0] and [0,a]. Equal
    /// tangential Ex/Ey powers at oblique incidence would not be unpolarized.
    /// Other positive Jones factorizations admit partial polarization without
    /// inferring a phase or degree of polarization from scalar intensity.
    pub fn try_new(
        source: &'s OrientedUniaxialInterface<'m>,
        factors: &'p [PreparedAbsorbingCrystalField<'s, 'm>],
        maximum_factors: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        kernel::poll_illumination(&mut cancelled)?;
        if factors.len() > maximum_factors {
            return Err(
                ReflectionError::Invalid("crystal illumination exceeds its factor budget").into(),
            );
        }
        for factor in factors {
            kernel::poll_illumination(&mut cancelled)?;
            if !core::ptr::eq(source, factor.source()) {
                return Err(ReflectionError::Invalid(
                    "crystal illumination factors require the exact original interface",
                )
                .into());
            }
        }
        let power = kernel::illumination_power(
            factors.iter().map(|factor| factor.boundary_response()),
            &mut cancelled,
        )?;
        Ok(Self {
            source,
            factors,
            power,
        })
    }

    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }

    pub const fn factors(&self) -> &'p [PreparedAbsorbingCrystalField<'s, 'm>] {
        self.factors
    }

    pub const fn factor_count(&self) -> usize {
        self.factors.len()
    }

    /// Absolute cycle/ensemble-mean normal powers at the original interface.
    pub const fn power(&self) -> AbsorbingCrystalIlluminationPower {
        self.power
    }

    /// Physical Cartesian second moments and heat. No mean complex field is
    /// invented for a partially polarized source. The original half-space
    /// selection is checked even for darkness; cancellation publishes nothing.
    pub fn field_statistics(
        &self,
        side: OrientedCrystalFieldSide,
        position_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalStatistics> {
        kernel::illumination_field_statistics(
            self.factors.iter().map(|factor| factor.boundary_response()),
            side,
            position_m,
            cancelled,
        )
        .map_err(CrystalIlluminationError::from)
    }

    /// Signed vector energy flux without forming electric squared moments.
    /// A thermal consumer need not admit an unrelated E-coherency magnitude
    /// or recompute dielectric work merely to inspect its boundary power.
    pub fn poynting_w_m2(
        &self,
        side: OrientedCrystalFieldSide,
        position_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3]> {
        kernel::illumination_poynting(
            self.factors.iter().map(|factor| factor.boundary_response()),
            side,
            position_m,
            cancelled,
        )
        .map_err(CrystalIlluminationError::from)
    }

    /// Add the original factors' analytic volume integrals, including each
    /// factor's ordinary/extraordinary interference and wavelength-dependent
    /// loss. A finite-energy consumer binds this volume and supplies the hold
    /// interval over which these frozen mean watts are deposited exactly once.
    pub fn slab_absorption(
        &self,
        area_m2: f64,
        start_depth_m: f64,
        end_depth_m: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalIlluminationSlabAbsorption<'_, 'p, 's, 'm>> {
        let (volume_m3, absorbed_power_w) = kernel::illumination_absorption(
            self.factors.iter().map(|factor| factor.boundary_response()),
            area_m2,
            start_depth_m,
            end_depth_m,
            cancelled,
        )?;
        Ok(AbsorbingCrystalIlluminationSlabAbsorption {
            illumination: self,
            area_m2,
            start_depth_m,
            end_depth_m,
            volume_m3,
            absorbed_power_w,
        })
    }
}

impl<'i, 'p, 's, 'm> AbsorbingCrystalIlluminationSlabAbsorption<'i, 'p, 's, 'm> {
    pub const fn illumination(&self) -> &'i AbsorbingCrystalIllumination<'p, 's, 'm> {
        self.illumination
    }
    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.illumination.source()
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
    pub const fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
}
