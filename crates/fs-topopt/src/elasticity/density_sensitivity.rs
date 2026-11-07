//! Analytic pullbacks to the original nodal physical-density field.
//!
//! These derivatives hold geometry, constitutive tensors and design occupancy
//! fixed. They differentiate the exact density-weighted mass, including both
//! the cell mean and normalized mass tile; no centroid-density approximation
//! or finite-difference material reconstruction is used.

use super::DensityElasticity;
use crate::control::{EvaluationStop, SolveControl};

mod moments;

/// A failed sensitivity query never returns a partial nodal gradient.
#[derive(Debug, Clone, PartialEq)]
pub enum ReferenceDensitySensitivityError {
    /// There is no explicitly authored nodal density field to differentiate.
    MissingNodalDensity,
    /// Shape, sign or finite-value requirements were not met.
    InvalidInput(&'static str),
    /// A required cell moment or accumulated nodal derivative is not finite.
    Nonfinite { cell: usize },
    /// The original shared cancellation or work control stopped evaluation.
    Stopped(EvaluationStop),
}

impl core::fmt::Display for ReferenceDensitySensitivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reference-density sensitivity: {self:?}")
    }
}

impl std::error::Error for ReferenceDensitySensitivityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Stopped(error) => Some(error),
            _ => None,
        }
    }
}

impl From<EvaluationStop> for ReferenceDensitySensitivityError {
    fn from(error: EvaluationStop) -> Self {
        Self::Stopped(error)
    }
}

impl DensityElasticity {
    /// Derivative of the physical quadratic mass form with respect to each
    /// original nodal density [kg/m³], in original mesh vertex order.
    ///
    /// `cell_mass_weights` are dimensionless and nonnegative. For a modal
    /// objective they are `mass_interp(projected_occupancy)`, with no reference
    /// density multiplier. `displacement` is a full three-component nodal field.
    /// The exact derivative is `sum_c weight_c * integral(N_a |u_h|² dV)`.
    /// Use zero displacement on fixed DOFs when differentiating reduced modes;
    /// constraint identity rows do not represent physical mass.
    ///
    /// One cell pass, fixed-size local arithmetic, and one nodal output vector;
    /// no matrix assembly or linear/eigenproblem solve. Cancellation is polled
    /// before allocation, every 256 cells, and before publication.
    pub fn reference_density_mass_gradient(
        &self,
        cell_mass_weights: &[f64],
        displacement: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ReferenceDensitySensitivityError> {
        control.checkpoint("reference-density-mass")?;
        let nodes = self.reference_density_nodes()?;
        let mut gradient = vec![0.0; nodes];
        self.accumulate_reference_density_mass_gradient(
            cell_mass_weights,
            displacement,
            1.0,
            &mut gradient,
            control,
        )?;
        control.checkpoint("reference-density-mass-publish")?;
        Ok(gradient)
    }

    /// Derivative of total physical inventory [kg] with respect to original
    /// nodal density [kg/m³]. Returned entries have units m³.
    ///
    /// Supply projected physical occupancies in [0,1]. Inventory uses linear
    /// occupancy even in the modal void regime. Original cell volume is split
    /// equally among its four P1 density basis functions, including fixed nodes.
    /// This differentiates density at fixed occupancy; the design filter's
    /// pullback remains the separate occupancy derivative.
    pub fn reference_density_inventory_gradient(
        &self,
        projected_occupancy: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ReferenceDensitySensitivityError> {
        control.checkpoint("reference-density-inventory")?;
        let nodes = self.reference_density_nodes()?;
        if projected_occupancy.len() != self.cells()
            || projected_occupancy
                .iter()
                .any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
        {
            return Err(ReferenceDensitySensitivityError::InvalidInput(
                "one projected occupancy in [0,1] per cell is required",
            ));
        }
        let mut gradient = vec![0.0; nodes];
        for (cell, tet) in self.tets.iter().enumerate() {
            if cell % 256 == 0 {
                control.checkpoint("reference-density-inventory")?;
            }
            let moment = (self.volumes[cell] * projected_occupancy[cell]) * 0.25;
            for &node in tet {
                let value = &mut gradient[node as usize];
                *value += moment;
                if !value.is_finite() {
                    return Err(ReferenceDensitySensitivityError::Nonfinite { cell });
                }
            }
        }
        control.checkpoint("reference-density-inventory-publish")?;
        Ok(gradient)
    }

    fn reference_density_nodes(&self) -> Result<usize, ReferenceDensitySensitivityError> {
        self.nodal_reference_density_kg_m3()
            .map(<[f64]>::len)
            .ok_or(ReferenceDensitySensitivityError::MissingNodalDensity)
    }

    // The modal owner accumulates all weighted mode derivatives into its single
    // unpublished candidate vector, avoiding one nodal allocation per mode.
    pub(crate) fn accumulate_reference_density_mass_gradient(
        &self,
        cell_mass_weights: &[f64],
        displacement: &[f64],
        coefficient: f64,
        gradient: &mut [f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(), ReferenceDensitySensitivityError> {
        control.checkpoint("reference-density-mass")?;
        let nodes = self.reference_density_nodes()?;
        if cell_mass_weights.len() != self.cells()
            || cell_mass_weights.iter().any(|v| !v.is_finite() || *v < 0.0)
            || displacement.len() != self.n()
            || gradient.len() != nodes
            || !coefficient.is_finite()
        {
            return Err(ReferenceDensitySensitivityError::InvalidInput(
                "finite nonnegative cell weights, full displacement and nodal output required",
            ));
        }
        for (cell, tet) in self.tets.iter().enumerate() {
            if cell % 256 == 0 {
                control.checkpoint("reference-density-mass")?;
            }
            let local =
                tet.map(|node| core::array::from_fn(|axis| displacement[3 * node as usize + axis]));
            let derivative = moments::moments(self.volumes[cell], cell_mass_weights[cell], local)
                .ok_or(ReferenceDensitySensitivityError::Nonfinite { cell })?;
            for (&node, value) in tet.iter().zip(derivative) {
                let total = &mut gradient[node as usize];
                *total = coefficient.mul_add(value, *total);
                if !total.is_finite() {
                    return Err(ReferenceDensitySensitivityError::Nonfinite { cell });
                }
            }
        }
        Ok(())
    }
}
