//! FrankenNumpy interop (bead gtql item c): converters between the
//! canonical [`Csr`] and fnp's owned array. SCOUT VERDICT (the bead's
//! "find the owned array type first"): `fnp-ndarray` holds only layout
//! metadata (`NdLayout`); the owned array type in the fnp workspace is
//! `fnp_ufunc::UFuncArray` — `{ shape: Vec<usize>, values: Vec<f64>, … }`
//! with OWNED storage. Like the fnx direction, these are therefore
//! CONVERSIONS by necessity: no zero-copy view exists into another
//! crate's owned `Vec`, and the borrowed views remain [`Csr::row`] /
//! [`Csr::to_dense`] on the fs-sparse side (borrow-vs-convert decided:
//! convert, documented).
//!
//! Gated behind `fnp-interop` so the default L1 crate stays
//! dependency-lean (`cargo test -p fs-sparse --features fnp-interop`).

use crate::Csr;
use fnp_dtype::DType;
use fnp_ufunc::UFuncArray;

/// fnp interop failure — total, structured, fail-closed.
#[derive(Debug, Clone, PartialEq)]
pub enum FnpInteropError {
    /// A CSR matrix corresponds to a 2-D array; other ranks are refused.
    NotTwoDimensional {
        /// The offending rank.
        ndim: usize,
    },
    /// A non-finite dense entry cannot enter a sparse kernel (they
    /// assume finite stored values); refused with its position.
    NonFinite {
        /// Row of the offending entry.
        row: usize,
        /// Column of the offending entry.
        col: usize,
    },
    /// Densifying would overflow `usize` (nrows·ncols) — refused rather
    /// than panicking.
    DenseTooLarge {
        /// Rows requested.
        nrows: usize,
        /// Columns requested.
        ncols: usize,
    },
    /// fnp refused the constructed array (shape/length disagreement —
    /// unreachable by construction, surfaced rather than unwrapped).
    Construction(String),
}

impl std::fmt::Display for FnpInteropError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FnpInteropError::NotTwoDimensional { ndim } => {
                write!(f, "CSR interop needs a 2-D array, got rank {ndim}")
            }
            FnpInteropError::NonFinite { row, col } => write!(
                f,
                "dense entry ({row}, {col}) is non-finite; sparse kernels assume finite \
                 stored values — clean the array first (fail closed)"
            ),
            FnpInteropError::DenseTooLarge { nrows, ncols } => {
                write!(f, "densifying {nrows}x{ncols} overflows the address space")
            }
            FnpInteropError::Construction(e) => write!(f, "fnp array construction refused: {e}"),
        }
    }
}

impl std::error::Error for FnpInteropError {}

/// `Csr → UFuncArray`: a dense 2-D f64 array, shape `[nrows, ncols]`,
/// row-major, unstored entries as `+0.0`. COPIES by necessity (owned
/// target) and densifies — O(nrows·ncols) memory, the caller's explicit
/// choice. Explicit stored zeros are indistinguishable from unstored
/// zeros after this conversion (documented lossiness; the round trip
/// `Csr → dense → Csr` is the identity exactly for CSRs without
/// explicit zeros).
///
/// # Errors
/// [`FnpInteropError::DenseTooLarge`] on `nrows·ncols` overflow;
/// [`FnpInteropError::Construction`] if fnp refuses the array
/// (unreachable by construction).
pub fn csr_to_dense_array(m: &Csr) -> Result<UFuncArray, FnpInteropError> {
    if m.nrows().checked_mul(m.ncols()).is_none() {
        return Err(FnpInteropError::DenseTooLarge {
            nrows: m.nrows(),
            ncols: m.ncols(),
        });
    }
    UFuncArray::new(vec![m.nrows(), m.ncols()], m.to_dense(), DType::F64)
        .map_err(|e| FnpInteropError::Construction(format!("{e:?}")))
}

/// `UFuncArray → Csr`: the 2-D array's f64 value plane, with exact
/// zeros (either sign) dropped and everything else stored — canonical
/// CSR falls out directly because the dense walk is already row-major
/// with ascending columns. Non-f64 dtypes convert through their f64
/// value plane (fnp mirrors integer arrays there); non-finite entries
/// are REFUSED with their position.
///
/// # Errors
/// [`FnpInteropError::NotTwoDimensional`] for any rank ≠ 2;
/// [`FnpInteropError::NonFinite`] naming the first offending entry.
pub fn dense_array_to_csr(a: &UFuncArray) -> Result<Csr, FnpInteropError> {
    let shape = a.shape();
    let [nrows, ncols] = *shape else {
        return Err(FnpInteropError::NotTwoDimensional { ndim: shape.len() });
    };
    let values = a.values();
    let mut row_ptr = vec![0usize; nrows + 1];
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    for r in 0..nrows {
        for c in 0..ncols {
            let v = values[r * ncols + c];
            if !v.is_finite() {
                return Err(FnpInteropError::NonFinite { row: r, col: c });
            }
            if v != 0.0 {
                row_ptr[r + 1] += 1;
                col_idx.push(c);
                vals.push(v);
            }
        }
    }
    for r in 0..nrows {
        row_ptr[r + 1] += row_ptr[r];
    }
    Ok(Csr::from_parts(nrows, ncols, row_ptr, col_idx, vals))
}

#[cfg(test)]
mod tests {
    use super::{FnpInteropError, csr_to_dense_array, dense_array_to_csr};
    use crate::Coo;
    use fnp_dtype::DType;
    use fnp_ufunc::UFuncArray;

    fn fixture() -> crate::Csr {
        let mut coo = Coo::new(3, 4);
        coo.push(0, 1, 2.5);
        coo.push(0, 3, -1.0);
        coo.push(2, 0, 4.0);
        coo.assemble()
    }

    #[test]
    fn round_trip_is_identity_without_explicit_zeros() {
        let m = fixture();
        let dense = csr_to_dense_array(&m).expect("densify");
        assert_eq!(dense.shape(), &[3, 4]);
        assert_eq!(dense.values()[0 * 4 + 1].to_bits(), 2.5f64.to_bits());
        assert_eq!(dense.values()[2 * 4 + 0].to_bits(), 4.0f64.to_bits());
        let back = dense_array_to_csr(&dense).expect("sparsify");
        assert_eq!(back.nrows(), 3);
        assert_eq!(back.nnz(), 3);
        for r in 0..3 {
            assert_eq!(back.row(r), m.row(r), "row {r} identical (bitwise values)");
        }
        // Dense → Csr → dense is the identity for finite arrays.
        let dense2 = csr_to_dense_array(&back).expect("densify again");
        assert_eq!(dense.values(), dense2.values());
    }

    #[test]
    fn refusals_are_structured_and_fail_closed() {
        // Rank 1 and rank 3 refuse.
        let v1 = UFuncArray::new(vec![4], vec![0.0; 4], DType::F64).expect("1-d");
        assert_eq!(
            dense_array_to_csr(&v1),
            Err(FnpInteropError::NotTwoDimensional { ndim: 1 })
        );
        let v3 = UFuncArray::new(vec![2, 2, 2], vec![0.0; 8], DType::F64).expect("3-d");
        assert!(matches!(
            dense_array_to_csr(&v3),
            Err(FnpInteropError::NotTwoDimensional { ndim: 3 })
        ));
        // Non-finite entries refuse with position.
        let bad =
            UFuncArray::new(vec![2, 2], vec![1.0, 0.0, f64::NAN, 3.0], DType::F64).expect("2-d");
        assert_eq!(
            dense_array_to_csr(&bad),
            Err(FnpInteropError::NonFinite { row: 1, col: 0 })
        );
        // Negative zero is dropped like positive zero (unstored).
        let signed = UFuncArray::new(vec![1, 2], vec![-0.0, 7.0], DType::F64).expect("2-d");
        let m = dense_array_to_csr(&signed).expect("sparsify");
        assert_eq!(m.nnz(), 1);
    }
}
