//! One local wall-frequency law shared by every plane-wave geometry.
//!
//! Constitutive owners return specific impedance under exp(-i omega t).
//! Duct geometry supplies circumference and surface slope. A material law may
//! depend on the actual radius and frequency; it must retain its own source
//! and applicability admission. This seam does not introduce a shell solver.

use crate::{DuctError, SegmentWave};
use fs_math::c64::C64;
use fs_phs::WallPin;

/// Material-owned, rigid-backed porous absorption on a distributed wall.
mod porous;
pub use porous::PorousLining;

/// Passive locally reacting surface impedance [Pa s/m].
///
/// Radius is the actual local inner radius, including chimney and quadrature
/// radii. The law must refuse outside its material/geometric support; missing
/// data must never turn into a rigid wall or a guessed resistance.
pub trait WallLaw {
    /// Complete impedance, including inertia and every retained material arm.
    fn specific_impedance(&self, radius_m: f64, omega: f64) -> Result<C64, DuctError>;
}

/// Explicit wall binding without changing the historical constant-pin API.
#[derive(Clone, Copy)]
pub enum WallModel<'a> {
    /// No deformable wall contribution.
    Rigid,
    /// Original constant surface mass/stiffness/resistance and arithmetic.
    Pin(&'a WallPin),
    /// Source-owned frequency response evaluated at each actual local radius.
    Law(&'a dyn WallLaw),
}

impl<'a> WallModel<'a> {
    /// Preserve the original meaning of an optional constant pin.
    #[must_use]
    pub const fn from_pin(pin: Option<&'a WallPin>) -> Self {
        match pin { Some(pin) => Self::Pin(pin), None => Self::Rigid }
    }

    #[must_use]
    pub const fn is_rigid(self) -> bool { matches!(self, Self::Rigid) }

    /// An outer None denotes a general law, not a rigid wall.
    pub(crate) const fn legacy_pin(self) -> Option<Option<&'a WallPin>> {
        match self { Self::Rigid => Some(None), Self::Pin(pin) => Some(Some(pin)), Self::Law(_) => None }
    }

    pub(crate) fn apply(self, wave: SegmentWave, radius: f64, omega: f64, slant: f64)
        -> Result<SegmentWave, DuctError> {
        match self {
            Self::Rigid => Ok(wave),
            Self::Pin(pin) => crate::apply_wall_to_wave(wave, radius, omega, pin, slant),
            Self::Law(law) => {
                if !radius.is_finite() || radius <= 0.0 || !omega.is_finite() || omega <= 0.0 {
                    return Err(DuctError::BadParameter { what: "wall law needs positive finite radius and frequency" });
                }
                let z = law.specific_impedance(radius, omega)?;
                if !z.re.is_finite() || !z.im.is_finite() || z.re < 0.0 || z == C64::ZERO {
                    return Err(DuctError::BadParameter { what: "wall law impedance must be passive, finite and nonzero" });
                }
                let y = z.recip().scale(core::f64::consts::TAU * radius);
                crate::apply_wall_admittance_to_wave(wave, radius, y, slant)
            }
        }
    }
}
