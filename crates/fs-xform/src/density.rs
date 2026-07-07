//! Raw SIMP density parameterization (plan §7.6): θ = per-cell densities
//! in [0, 1]. This bead owns the raw field lever + validity diagnostics;
//! Helmholtz filtering and Heaviside projection are topo-simp's,
//! downstream. Trivially linear: the design perturbation IS δθ.

use crate::XformError;

/// A per-cell density field.
#[derive(Debug, Clone)]
pub struct DensityField {
    /// Cell count.
    pub cells: usize,
}

impl DensityField {
    /// DOF count (one density per cell).
    #[must_use]
    pub fn dof(&self) -> usize {
        self.cells
    }

    /// Validate θ ∈ [0, 1]^cells; the first violation is a structured
    /// refusal naming the component (P10).
    ///
    /// # Errors
    /// [`XformError::DofMismatch`] / [`XformError::OutOfBounds`].
    pub fn validate(&self, theta: &[f64]) -> Result<(), XformError> {
        if theta.len() != self.cells {
            return Err(XformError::DofMismatch {
                expected: self.cells,
                got: theta.len(),
            });
        }
        for (index, &value) in theta.iter().enumerate() {
            if !(value.is_finite() && (0.0..=1.0).contains(&value)) {
                return Err(XformError::OutOfBounds {
                    index,
                    value,
                    bound: "density in [0, 1]",
                });
            }
        }
        Ok(())
    }

    /// The density perturbation for a design step: exactly δθ (identity
    /// Jacobian — the chain rule to charts/physics happens downstream).
    ///
    /// # Errors
    /// [`XformError::DofMismatch`] on a wrong-length δθ.
    pub fn perturbation<'a>(&self, dtheta: &'a [f64]) -> Result<&'a [f64], XformError> {
        if dtheta.len() != self.cells {
            return Err(XformError::DofMismatch {
                expected: self.cells,
                got: dtheta.len(),
            });
        }
        Ok(dtheta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validity_names_the_offender() {
        let field = DensityField { cells: 3 };
        assert!(field.validate(&[0.0, 0.5, 1.0]).is_ok());
        match field.validate(&[0.0, 1.5, 1.0]) {
            Err(XformError::OutOfBounds {
                index: 1, value, ..
            }) => {
                assert!((value - 1.5).abs() < 1e-15);
            }
            other => panic!("expected OutOfBounds at index 1, got {other:?}"),
        }
        assert!(matches!(
            field.validate(&[0.0]),
            Err(XformError::DofMismatch {
                expected: 3,
                got: 1
            })
        ));
        assert_eq!(
            field.perturbation(&[0.1, 0.2, 0.3]).unwrap(),
            &[0.1, 0.2, 0.3]
        );
    }
}
