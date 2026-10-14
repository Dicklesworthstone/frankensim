//! Shared admission of external force plus thermal buoyancy, before mutation.

/// A force field that cannot be applied to the existing thermal solver.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ThermalForcingError {
    /// The field must contain one vector per cell, including wall cells.
    Length {
        /// Required number of cell vectors.
        expected: usize,
        /// Supplied number of cell vectors.
        actual: usize,
    },
    /// An external component or its sum with buoyancy is not finite.
    NonFinite {
        /// Linear grid-cell index.
        cell: usize,
        /// Zero-based force-vector component.
        component: usize,
    },
    /// The declared wall plane cannot describe the temperature storage.
    Layout,
}
impl core::fmt::Display for ThermalForcingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Length { expected, actual } => write!(f, "thermal force field has {actual} cells, expected {expected}"),
            Self::NonFinite { cell, component } => write!(f, "nonfinite thermal force at cell {cell}, component {component}"),
            Self::Layout => write!(f, "invalid thermal force/temperature layout"),
        }
    }
}
impl std::error::Error for ThermalForcingError {}

pub(crate) fn validate_external_forces<const D: usize, const Q: usize>(
    forces: &[[f64; D]],
    temperature_populations: &[[f64; Q]],
    wall_plane: usize,
    gbeta: f64,
    reference_temperature: f64,
) -> Result<(), ThermalForcingError> {
    let n = temperature_populations.len();
    if forces.len() != n {
        return Err(ThermalForcingError::Length { expected: n, actual: forces.len() });
    }
    if D == 0 || wall_plane == 0 || wall_plane > n / 2 {
        return Err(ThermalForcingError::Layout);
    }
    for (cell, force) in forces.iter().enumerate() {
        for (component, value) in force.iter().enumerate() {
            if !value.is_finite() {
                return Err(ThermalForcingError::NonFinite { cell, component });
            }
        }
        if cell >= wall_plane && cell < n - wall_plane {
            let temperature: f64 = temperature_populations[cell].iter().sum();
            let total = gbeta * (temperature - reference_temperature) + force[D - 1];
            if !total.is_finite() {
                return Err(ThermalForcingError::NonFinite { cell, component: D - 1 });
            }
        }
    }
    Ok(())
}
