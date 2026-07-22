//! Nodal scalar data: the one carrier every boundary value and source
//! term in this crate uses. A field is either UNIFORM (one number that
//! holds everywhere it is read) or NODAL (one number per mesh vertex).
//!
//! There is no closure-valued field. A closure could not be snapshotted,
//! replayed, or content-addressed, so spatially varying data enters as
//! evaluated nodal values whose interpolation this crate declares:
//! P₁ (linear on each element/face), matching the trial space exactly.

use crate::ConductionError;

/// Nodal or uniform scalar data over a mesh.
#[derive(Debug, Clone, PartialEq)]
pub enum ScalarField {
    /// One value, everywhere.
    Uniform(f64),
    /// One value per mesh vertex, in canonical vertex order.
    Nodal(Vec<f64>),
}

impl ScalarField {
    /// Build a nodal field, checking length and finiteness.
    ///
    /// # Errors
    /// [`ConductionError::FieldLength`] on a length mismatch and
    /// [`ConductionError::NonFinite`] for a non-finite entry.
    pub fn nodal(
        field: &'static str,
        vertex_count: usize,
        values: Vec<f64>,
    ) -> Result<ScalarField, ConductionError> {
        if values.len() != vertex_count {
            return Err(ConductionError::FieldLength {
                field,
                expected: vertex_count,
                found: values.len(),
            });
        }
        for &v in &values {
            crate::require_finite(field, v)?;
        }
        Ok(ScalarField::Nodal(values))
    }

    /// Build a uniform field, checking finiteness.
    ///
    /// # Errors
    /// [`ConductionError::NonFinite`] for a non-finite value.
    pub fn uniform(field: &'static str, value: f64) -> Result<ScalarField, ConductionError> {
        crate::require_finite(field, value)?;
        Ok(ScalarField::Uniform(value))
    }

    /// The value at a mesh vertex.
    ///
    /// # Panics
    /// If a nodal field is indexed outside its array — an internal
    /// indexing bug, never a runtime condition, because every
    /// constructor pins the length to the mesh's vertex count.
    #[must_use]
    pub fn at(&self, vertex: usize) -> f64 {
        match self {
            ScalarField::Uniform(v) => *v,
            ScalarField::Nodal(values) => values[vertex],
        }
    }

    /// Validate against a mesh's vertex count.
    ///
    /// # Errors
    /// [`ConductionError::FieldLength`] when a nodal array does not
    /// match; [`ConductionError::NonFinite`] for a non-finite entry.
    pub fn validate(
        &self,
        field: &'static str,
        vertex_count: usize,
    ) -> Result<(), ConductionError> {
        match self {
            ScalarField::Uniform(v) => {
                crate::require_finite(field, *v)?;
            }
            ScalarField::Nodal(values) => {
                if values.len() != vertex_count {
                    return Err(ConductionError::FieldLength {
                        field,
                        expected: vertex_count,
                        found: values.len(),
                    });
                }
                for &v in values {
                    crate::require_finite(field, v)?;
                }
            }
        }
        Ok(())
    }

    /// True for a spatially uniform field.
    #[must_use]
    pub const fn is_uniform(&self) -> bool {
        matches!(self, ScalarField::Uniform(_))
    }
}
