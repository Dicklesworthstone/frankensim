//! Cochain storage: one 128-byte-aligned value buffer per degree
//! (fs-soa substrate), a SINGLE dimension tag for the whole container
//! (the workspace Qty precedent: one `Dims` per field, not one Qty per
//! element), and zero-copy view descriptors for the membrane.

use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_soa::{FieldBuf, RawView};

/// A k-cochain on a tet complex: one f64 per k-cell, in canonical cell
/// order (the complex's deterministic sorted tables).
#[derive(Debug, Clone)]
pub struct Cochain {
    values: FieldBuf<f64>,
    degree: u8,
    dims: Dims,
}

/// Number of k-cells in the complex.
///
/// # Panics
/// If `degree > 3`.
#[must_use]
pub fn cell_count(complex: &TetComplex, degree: u8) -> usize {
    match degree {
        0 => complex.vertex_count,
        1 => complex.edges.len(),
        2 => complex.faces.len(),
        3 => complex.tets.len(),
        _ => panic!("cochain degree must be 0..=3"),
    }
}

impl Cochain {
    /// Zero cochain of the given degree, dimensionless.
    #[must_use]
    pub fn zeros(complex: &TetComplex, degree: u8) -> Cochain {
        let n = cell_count(complex, degree);
        let mut values = FieldBuf::with_capacity(n);
        for _ in 0..n {
            values.push(0.0);
        }
        Cochain {
            values,
            degree,
            dims: Dims([0; 5]),
        }
    }

    /// Cochain from per-cell values.
    ///
    /// # Panics
    /// If the value count does not match the complex's k-cell count.
    #[must_use]
    pub fn from_values(complex: &TetComplex, degree: u8, vals: &[f64], dims: Dims) -> Cochain {
        assert_eq!(
            vals.len(),
            cell_count(complex, degree),
            "cochain length must match the k-cell count"
        );
        let mut values = FieldBuf::with_capacity(vals.len());
        for &v in vals {
            values.push(v);
        }
        Cochain {
            values,
            degree,
            dims,
        }
    }

    /// Degree k.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.degree
    }

    /// The dimension tag (units of the represented integral quantity).
    #[must_use]
    pub const fn dims(&self) -> Dims {
        self.dims
    }

    /// Values in canonical cell order.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        self.values.as_slice()
    }

    /// Mutable values in canonical cell order.
    pub fn values_mut(&mut self) -> &mut [f64] {
        self.values.as_mut_slice()
    }

    /// Zero-copy membrane descriptor (address-free `descr()` for logs).
    #[must_use]
    pub fn view(&self) -> RawView {
        self.values.view(&format!("cochain{}", self.degree))
    }

    /// Apply the exact integer exterior derivative (the complex's
    /// incidence operator) to this cochain's f64 values: returns the
    /// (k+1)-cochain with the same dimension tag. Coefficients are
    /// ±1 so the arithmetic is pure signed addition — no rounding
    /// beyond the additions themselves, fixed row order.
    ///
    /// # Panics
    /// If `degree == 3` (no d beyond volume forms on a 3-complex).
    #[must_use]
    pub fn d(&self, complex: &TetComplex) -> Cochain {
        let inc = match self.degree {
            0 => complex.d0(),
            1 => complex.d1(),
            2 => complex.d2(),
            _ => panic!("d on a 3-cochain of a 3-complex is zero-dimensional"),
        };
        let x = self.values();
        let mut out = Vec::with_capacity(inc.rows.len());
        for row in &inc.rows {
            let mut acc = 0.0f64;
            for &(col, sign) in row {
                acc += f64::from(sign) * x[col];
            }
            out.push(acc);
        }
        let mut values = FieldBuf::with_capacity(out.len());
        for &v in &out {
            values.push(v);
        }
        Cochain {
            values,
            degree: self.degree + 1,
            dims: self.dims,
        }
    }
}
