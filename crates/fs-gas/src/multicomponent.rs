//! Conservative transport of an evolving neutral ideal-gas mixture.
//!
//! Every species has its own mass equation and retained NASA-9 caloric source.
//! The stored energy subtracts a separate, constant reference per species;
//! this linear change of conserved variables preserves physical total energy.
//! HLLC/MUSCL is an Estimated engineering lane. It does not establish a fully
//! discrete entropy inequality or exact pressure equilibrium across a moving
//! material interface with different caloric laws.
//!
//! Shared boundary types admit periodic, reflecting and transmissive ends here.
//! A frozen `Boundary::Reservoir` has no species-mass axis or species energy
//! references and therefore refuses. Species-bearing compact-reservoir ports
//! use [`MixtureGas::reservoir_flux`] with the original local compositions.

mod characteristic;
mod eos;
mod evolve;
mod flux;
mod normal;
mod wall;

/// Conservative shockwaves and species mixing through finite-plenum networks.
pub mod network;

/// Two-layer control volumes with species-dependent calorics and interface work.
pub mod compartment;

/// Conservative molecular diffusion, heat conduction and viscous momentum transport.
pub mod dissipation;

/// Atomic finite-rate chemistry coupled to conservative species transport.
pub mod reactive;

/// Stationary wall drag and conservative heat exchange with finite wall storage.
pub mod wall_exchange;

pub use crate::dynamics::{
    Boundaries, Boundary, CaloricDomain, Mesh1d, Primitive, Reconstruction, StepControl,
};
pub use eos::MixtureGas;
pub use evolve::{MixtureEuler1d, StepReport, Totals};
pub use flux::FaceFlux;
pub use normal::{NormalCell, NormalPrimitive, SpatialCell, SpatialPrimitive};

use fs_thermochem::{ThermochemErrorV1, kinetics::KineticsErrorV1};

/// Same bounded species axis as the original thermochemical owner.
pub const MAX_SPECIES: usize = 128;
const MAX_COMPONENTS: usize = MAX_SPECIES + 2;
const MAX_NORMAL_COMPONENTS: usize = MAX_SPECIES + 3;
const MAX_SPATIAL_COMPONENTS: usize = MAX_SPECIES + 4;

/// Owned conservative cell average, in canonical source-species order.
#[derive(Clone, Debug, PartialEq)]
pub struct Conserved {
    /// Species partial mass densities [kg/m³]; exact zero is permitted.
    pub partial_densities: Vec<f64>,
    /// Axial momentum density [kg/(m² s)].
    pub momentum: f64,
    /// Total energy less `sum(rho_i * e_i(Tmin))` [J/m³].
    pub energy: f64,
}
impl Conserved {
    /// Borrow the physical state without allocating.
    pub fn as_cell(&self) -> Cell<'_> {
        Cell {
            partial_densities: &self.partial_densities,
            momentum: self.momentum,
            energy: self.energy,
        }
    }
}

/// Borrowed conservative cell average, including all species.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cell<'a> {
    /// Species partial mass densities [kg/m³], in the owner's canonical order.
    pub partial_densities: &'a [f64],
    /// Axial momentum density [kg/(m² s)].
    pub momentum: f64,
    /// Species-referenced total energy density [J/m³].
    pub energy: f64,
}
impl Cell<'_> {
    /// Mixture density is the sum of the conserved species, never a second field.
    pub fn density(self) -> f64 {
        sum(self.partial_densities.iter().copied())
    }
}

/// Source, numerical, allocation, work-budget or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// Invalid dimensions, controls, or nonrepresentable physical state.
    Invalid(&'static str),
    /// The original NASA-9 source refused an evaluation.
    Thermochemistry(ThermochemErrorV1),
    /// The original shared positive-Cv/path admission refused the source.
    CaloricPath(KineticsErrorV1),
    /// A reflected wave leaves the admitted source path or numerical budget.
    WallCaloricPath(&'static str),
    /// A port wave leaves the admitted source path or numerical budget.
    CharacteristicPath(&'static str),
    /// Physical mesh admission refused.
    Geometry(crate::dynamics::Error),
    /// Sensible energy falls outside the source temperature interval.
    EnergyDomain,
    /// The bounded temperature inversion did not converge.
    TemperatureIterations,
    /// A stage exceeded the caller's actual face-speed CFL.
    CflExceeded { observed: f64, maximum: f64 },
    /// All atomic time-step retries failed.
    StepRetries { attempts: usize, last: Box<Error> },
    /// Bounded workspace reservation failed.
    Allocation,
    /// Cancellation occurred before publication.
    Cancelled,
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent Euler gas: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<ThermochemErrorV1> for Error {
    fn from(value: ThermochemErrorV1) -> Self {
        Self::Thermochemistry(value)
    }
}
impl From<crate::dynamics::Error> for Error {
    fn from(value: crate::dynamics::Error) -> Self {
        Self::Geometry(value)
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
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let y = value - correction;
        let next = total + y;
        correction = (next - total) - y;
        total = next;
    }
    total
}
fn cell(values: &[f64], species: usize) -> Cell<'_> {
    Cell {
        partial_densities: &values[..species],
        momentum: values[species],
        energy: values[species + 1],
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod wall_tests;

#[cfg(test)]
mod characteristic_tests;
