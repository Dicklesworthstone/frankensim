//! Conservative one-dimensional Euler flow with original frozen NASA-9 gas laws.
//!
//! HLLC restores the contact wave; inadmissible star states fall back to an HLL
//! flux. MUSCL reconstruction and SSP RK2 are bounded by local signal speeds,
//! with transactional step reduction rather than clipping conserved variables.
//! This is an Estimated engineering baseline, not a discrete entropy theorem.
//! References: Toro, Spruce and Speares, Shock Waves 4 (1994), 25–34,
//! DOI 10.1007/BF01414629; Ketcheson, LeVeque and del Razo, *Riemann Problems
//! and Jupyter Solutions* (2020), Euler and approximate-Euler chapters.

pub(crate) mod characteristic;
mod eos;
mod evolve;
mod flux;
mod geometry;
/// Synchronized shock propagation through passages and finite rigid gas volumes.
pub mod network;
pub(crate) mod wall;

pub use eos::{CaloricDomain, FrozenGas};
pub use evolve::{Boundaries, Boundary, Euler1d, Reconstruction, StepControl, StepReport, Totals};
pub use flux::{FaceFlux, NormalConserved, NormalFaceFlux};
pub use geometry::Mesh1d;

use fs_thermochem::mixture::FrozenIdealGasMixtureErrorV1;

/// Hard state-allocation bound, checked before workspace allocation.
pub const MAX_EULER_CELLS: usize = 1_048_576;

/// Conservative cell average in coherent SI with the owner's declared energy gauge.
///
/// The physical total energy density is `energy + density * gas.energy_reference()`.
/// Evolving this exact linear change of variables avoids losing sensible energy
/// to a large formation reference. The gauge is valid only for this fixed mixture.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Conserved {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Axial momentum density [kg/(m² s)].
    pub momentum: f64,
    /// Total energy relative to the retained reference [J/m³], including kinetic energy.
    pub energy: f64,
}

impl Conserved {
    fn array(self) -> [f64; 3] {
        [self.density, self.momentum, self.energy]
    }
    fn from_array(q: [f64; 3]) -> Self {
        Self {
            density: q[0],
            momentum: q[1],
            energy: q[2],
        }
    }
    fn finite(self) -> bool {
        self.array().iter().all(|x| x.is_finite())
    }
}

/// One decoded gas state; every field refers to the same source and state.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Primitive {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Axial velocity [m/s].
    pub velocity: f64,
    /// Absolute pressure [Pa].
    pub pressure: f64,
    /// Absolute temperature [K].
    pub temperature: f64,
    /// Frozen-composition sound speed [m/s].
    pub sound_speed: f64,
}

/// Atomic source, numerical, work-budget or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// Invalid dimensions, controls, or non-finite/unrepresentable state.
    Invalid(&'static str),
    /// Original source evaluation refused.
    Thermochemistry(FrozenIdealGasMixtureErrorV1),
    /// The declared source path is unsupported or has unproven positive Cv.
    CaloricPath(&'static str),
    /// Energy is outside the declared source temperature interval.
    EnergyDomain {
        /// Requested specific energy relative to the source reference [J/kg].
        requested: f64,
        /// Largest supported relative specific internal energy [J/kg].
        maximum: f64,
    },
    /// Temperature inversion exhausted its explicit iteration budget.
    TemperatureIterations,
    /// A trial stage exceeded the explicitly admitted face-speed CFL.
    CflExceeded {
        /// Actual trial face-speed CFL.
        observed: f64,
        /// Caller-admitted maximum.
        maximum: f64,
    },
    /// All bounded time-step attempts failed; state and time are unchanged.
    StepRetries {
        /// Number of attempted steps, including the initial attempt.
        attempts: usize,
        /// Last concrete numerical/source refusal.
        last: Box<Error>,
    },
    /// A bounded workspace allocation failed.
    Allocation,
    /// Cancellation occurred before publication; state and time are unchanged.
    Cancelled,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "one-dimensional Euler gas: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<FrozenIdealGasMixtureErrorV1> for Error {
    fn from(error: FrozenIdealGasMixtureErrorV1) -> Self {
        Self::Thermochemistry(error)
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn allocated<T: Clone>(n: usize, value: T) -> Result<Vec<T>, Error> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| Error::Allocation)?;
    v.resize(n, value);
    Ok(v)
}

#[cfg(test)]
mod tests;
