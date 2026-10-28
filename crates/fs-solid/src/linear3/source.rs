//! Shared material-source lowering for thermal loading of the solid owner.
//!
//! Both prescribed-temperature mechanics and spatial enthalpy coupling use
//! these exact adapters. Constitutive symmetry and thermal symmetry are
//! independent data choices; no solver or material-name dispatch lives here.
//! Closed-matter inertia is always supplied from the original reference mass.

mod expansion;
pub(super) mod quadratic_thermal;
mod stiffness;

pub use expansion::{ThermalExpansionModel, ThermalExpansionSource};
pub use quadratic_thermal::{
    NodalThermalTemperatures, QuadraticThermalSource, QuadraticThermalSources,
    ResolvedQuadraticThermalStrains,
};
pub use stiffness::{SolidStiffnessModel, ThermalStiffnessModel};

use super::TetElasticError;
use fs_material::state_point::MaterialStatePointError;

/// Original source or solid admission refusal, shared by all consumers.
#[derive(Clone, Debug, PartialEq)]
pub enum SolidMaterialError {
    /// Missing, inadmissible, ambiguous or unsupported source data.
    Material(MaterialStatePointError),
    /// Invalid constitutive tensor, frame transform or small-strain state.
    Solid(TetElasticError),
}
impl core::fmt::Display for SolidMaterialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "solid material source: {error}"),
            Self::Solid(error) => write!(f, "solid material admission: {error}"),
        }
    }
}
impl std::error::Error for SolidMaterialError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Solid(error) => Some(error),
        }
    }
}
impl From<MaterialStatePointError> for SolidMaterialError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<TetElasticError> for SolidMaterialError {
    fn from(error: TetElasticError) -> Self {
        Self::Solid(error)
    }
}
