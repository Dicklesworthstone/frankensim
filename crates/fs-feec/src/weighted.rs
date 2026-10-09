//! Piecewise-constant material weights for Whitney mass/Hodge operators.
//!
//! These are Galerkin pairings of primal forms, with rows representing the
//! corresponding dual test functionals. They do not change the incidence
//! operators or certify a PDE's boundary, gauge, or inf-sup choices.

use fs_ivl::Interval;
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

use crate::whitney::ElementGeometry;

/// A constant material weight on one tetrahedron, in the global Cartesian
/// frame. Tensor weights act on vector proxies of one- and two-forms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CellWeight {
    /// Strictly positive scalar, admitted for degrees zero through three.
    Scalar(f64),
    /// Exactly symmetric positive-definite tensor, admitted for degrees one
    /// and two. An inconclusive interval positivity check refuses the tensor.
    Tensor([[f64; 3]; 3]),
}

/// Explicit bounds on input cells, output rows, and staged contributions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WeightedAssemblyLimits {
    /// Maximum number of cells visited.
    pub max_cells: usize,
    /// Maximum rows, including isolated degrees of freedom.
    pub max_dofs: usize,
    /// Maximum entries before shared-cell contributions are summed.
    pub max_triplets: usize,
}

impl WeightedAssemblyLimits {
    pub(crate) fn check(
        self,
        cells: usize,
        dofs: usize,
        triplets: usize,
    ) -> Result<(), WeightedError> {
        for (resource, requested, limit) in [
            ("cells", cells, self.max_cells),
            ("degrees of freedom", dofs, self.max_dofs),
            ("triplets", triplets, self.max_triplets),
        ] {
            if requested > limit {
                return Err(WeightedError::LimitExceeded {
                    resource,
                    requested,
                    limit,
                });
            }
        }
        Ok(())
    }
}

/// A coefficient, shape, or work refusal. No partially assembled matrix is
/// returned on any error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WeightedError {
    /// The differential-form degree does not exist on this complex.
    InvalidDegree(u8),
    /// There must be exactly one material weight per cell.
    CoefficientCount {
        /// Number of cells requiring a weight.
        expected: usize,
        /// Number of supplied weights.
        actual: usize,
    },
    /// Nonfinite, nonpositive, asymmetric, unsupported, or unproved weight.
    InvalidCoefficient {
        /// Cell whose material coefficient was refused.
        cell: usize,
        /// Failed coefficient assumption.
        reason: &'static str,
    },
    /// Inconsistent geometry arrays, missing topology, or nonfinite metrics.
    InvalidGeometry {
        /// Cell with invalid metrics/topology (zero for a whole-array mismatch).
        cell: usize,
        /// Failed geometry assumption.
        reason: &'static str,
    },
    /// SI exponent arithmetic is not representable.
    DimensionOverflow,
    /// An explicitly declared work/storage extent was exceeded.
    LimitExceeded {
        /// Bounded resource whose requested extent was too large.
        resource: &'static str,
        /// Required extent before any assembly allocation.
        requested: usize,
        /// Caller-admitted maximum extent.
        limit: usize,
    },
    /// A checked array-size calculation overflowed.
    SizeOverflow,
    /// A scratch or result allocation was refused.
    Allocation,
    /// An element contribution or accumulated matrix entry is not finite.
    NonFiniteEntry,
    /// The caller's checkpoint requested cancellation.
    Cancelled,
}

impl core::fmt::Display for WeightedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDegree(degree) => write!(f, "unsupported weighted form degree {degree}"),
            Self::CoefficientCount { expected, actual } => {
                write!(f, "expected {expected} cell weights, received {actual}")
            }
            Self::InvalidCoefficient { cell, reason } => {
                write!(f, "material coefficient on cell {cell}: {reason}")
            }
            Self::InvalidGeometry { cell, reason } => {
                write!(f, "weighted geometry on cell {cell}: {reason}")
            }
            Self::DimensionOverflow => {
                f.write_str("weighted operator SI dimension exponent overflow")
            }
            Self::LimitExceeded {
                resource,
                requested,
                limit,
            } => write!(
                f,
                "weighted assembly requires {requested} {resource}, above limit {limit}"
            ),
            Self::SizeOverflow => f.write_str("weighted assembly storage extent overflow"),
            Self::Allocation => f.write_str("weighted assembly allocation refused"),
            Self::NonFiniteEntry => {
                f.write_str("weighted matrix entry is not representable as finite f64")
            }
            Self::Cancelled => f.write_str("weighted assembly cancelled before publication"),
        }
    }
}

impl std::error::Error for WeightedError {}

/// A material-weighted Galerkin mass matrix and its physical dimensions.
///
/// A degree-k basis integrates to one over its k-cell, so a 3-D entry has
/// dimensions `coefficient_dims * metre^(3-2k)`. Applying the numerical matrix
/// adds these dimensions to the input cochain dimensions. No continuum
/// coercivity or solution-error certificate follows from this matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct WeightedMass {
    matrix: Csr,
    degree: u8,
    coefficient_dims: Dims,
    entry_dims: Dims,
}

/// Galerkin Hodge-star interpretation of the same primal/dual pairing.
pub type WeightedStar = WeightedMass;

impl WeightedMass {
    pub(crate) fn from_matrix(
        matrix: Csr,
        degree: u8,
        coefficient_dims: Dims,
        entry_dims: Dims,
    ) -> Self {
        Self {
            matrix,
            degree,
            coefficient_dims,
            entry_dims,
        }
    }

    /// Canonical CSR in the complex's sorted degree-k cell order.
    #[must_use]
    pub const fn matrix(&self) -> &Csr {
        &self.matrix
    }

    /// Degree of the primal trial and test forms.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.degree
    }

    /// Dimensions shared by all per-cell material coefficients.
    #[must_use]
    pub const fn coefficient_dims(&self) -> Dims {
        self.coefficient_dims
    }

    /// Dimensions of each matrix entry after geometric integration.
    #[must_use]
    pub const fn entry_dims(&self) -> Dims {
        self.entry_dims
    }

    /// Dimensions of the dual result when this operator acts on a cochain.
    ///
    /// # Errors
    /// Returns [`WeightedError::DimensionOverflow`] on SI exponent overflow.
    pub fn output_dims(&self, input_dims: Dims) -> Result<Dims, WeightedError> {
        self.entry_dims
            .checked_plus(input_dims)
            .ok_or(WeightedError::DimensionOverflow)
    }
}

pub(crate) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WeightedError> {
    if cancelled() {
        Err(WeightedError::Cancelled)
    } else {
        Ok(())
    }
}

/// Sufficient outward-rounded Sylvester/Schur test of the supplied tensor.
/// Refusal includes positive tensors too close to singularity for this rung.
pub(crate) fn validate_tensor<const N: usize>(
    tensor: &[[f64; N]; N],
    cell: usize,
) -> Result<(), WeightedError> {
    let fail = |reason| WeightedError::InvalidCoefficient { cell, reason };
    let mut scale = 0.0f64;
    for i in 0..N {
        for j in 0..N {
            let value = tensor[i][j];
            if !value.is_finite() {
                return Err(fail("nonfinite tensor"));
            }
            if value != tensor[j][i] {
                return Err(fail("tensor must be exactly symmetric"));
            }
            scale = scale.max(value.abs());
        }
    }
    if scale == 0.0 {
        return Err(fail("zero tensor is not positive definite"));
    }
    let mut schur = [[Interval::point(0.0); N]; N];
    for i in 0..N {
        for j in 0..N {
            // Enclose division itself, rather than certify a rounded matrix.
            schur[i][j] = Interval::point(tensor[i][j]) / Interval::point(scale);
        }
    }
    for k in 0..N {
        let pivot = schur[k][k];
        if pivot.lo() <= 0.0 || !pivot.hi().is_finite() {
            return Err(fail(
                "positive definiteness not established by interval pivots",
            ));
        }
        for i in k + 1..N {
            for j in i..N {
                let entry = schur[i][j] - schur[i][k] * schur[j][k] / pivot;
                schur[i][j] = entry;
                schur[j][i] = entry;
            }
        }
    }
    Ok(())
}

/// Assemble piecewise-constant Whitney element integrals, rounded to f64,
/// in deterministic cell/insertion order. Coefficients may jump across cells.
///
/// `cancelled` is polled throughout validation, assembly, sorting, reduction,
/// and publication. A scoped caller can pass
/// `&mut || cx.checkpoint().is_err()`. All buffers remain local until the final
/// checkpoint succeeds. Arbitrary caller-supplied `ElementGeometry` is checked
/// for finite metrics and topology correspondence, not independently proved
/// to describe a particular embedded mesh.
///
/// # Errors
/// Returns [`WeightedError`] for invalid coefficients/geometry, unsupported
/// degrees, exhausted bounds, nonrepresentable entries, allocation refusal,
/// dimension overflow, or cancellation.
pub fn weighted_mass_matrix(
    complex: &TetComplex,
    geo: &ElementGeometry,
    degree: u8,
    weights: &[CellWeight],
    coefficient_dims: Dims,
    limits: WeightedAssemblyLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WeightedMass, WeightedError> {
    poll(cancelled)?;
    let local_entries = match degree {
        0 | 2 => 16,
        1 => 36,
        3 => 1,
        _ => return Err(WeightedError::InvalidDegree(degree)),
    };
    let cells = complex.tets.len();
    let dofs = crate::cochain::cell_count(complex, degree);
    let triplets = cells
        .checked_mul(local_entries)
        .ok_or(WeightedError::SizeOverflow)?;
    limits.check(cells, dofs, triplets)?;
    if weights.len() != cells {
        return Err(WeightedError::CoefficientCount {
            expected: cells,
            actual: weights.len(),
        });
    }
    if geo.vol_signed.len() != cells || geo.grads.len() != cells || geo.gram.len() != cells {
        return Err(WeightedError::InvalidGeometry {
            cell: 0,
            reason: "geometry arrays must match tetrahedra",
        });
    }
    let entry_dims = coefficient_dims
        .checked_plus(Dims([3 - 2 * degree as i8, 0, 0, 0, 0, 0]))
        .ok_or(WeightedError::DimensionOverflow)?;
    for (cell, weight) in weights.iter().enumerate() {
        poll(cancelled)?;
        match weight {
            CellWeight::Scalar(value) => {
                if !value.is_finite() || *value <= 0.0 {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "scalar must be finite and strictly positive",
                    });
                }
            }
            CellWeight::Tensor(tensor) => {
                if !matches!(degree, 1 | 2) {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "tensor weights require one- or two-forms",
                    });
                }
                validate_tensor(tensor, cell)?;
            }
        }
        validate_geometry(complex, geo, cell, degree)?;
    }
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(triplets)
        .map_err(|_| WeightedError::Allocation)?;
    for (cell, weight) in weights.iter().enumerate() {
        poll(cancelled)?;
        let normalized;
        let (scale, tensor) = match weight {
            CellWeight::Scalar(value) => (*value, None),
            CellWeight::Tensor(tensor) => {
                let scale = tensor.iter().flatten().fold(0.0f64, |s, v| s.max(v.abs()));
                normalized = tensor.map(|row| row.map(|v| v / scale));
                (scale, Some(&normalized))
            }
        };
        crate::whitney::element_mass_entries(
            complex,
            geo,
            degree,
            cell,
            tensor,
            &mut |row, col, value| {
                // Mirror one evaluated triangle so finite-precision symmetry is
                // exact, independently of reassociation in transposed formulas.
                if row > col {
                    return Ok(());
                }
                let value = scale * value;
                if !value.is_finite() {
                    return Err(WeightedError::NonFiniteEntry);
                }
                if row == col && value <= 0.0 {
                    return Err(WeightedError::InvalidGeometry {
                        cell,
                        reason: "positive mass diagonal is not representable",
                    });
                }
                entries.push((row, col, value));
                if row != col {
                    entries.push((col, row, value));
                }
                Ok(())
            },
        )?;
    }
    let matrix = assemble_triplets(dofs, entries, cancelled)?;
    Ok(WeightedMass::from_matrix(
        matrix,
        degree,
        coefficient_dims,
        entry_dims,
    ))
}

fn validate_geometry(
    complex: &TetComplex,
    geo: &ElementGeometry,
    cell: usize,
    degree: u8,
) -> Result<(), WeightedError> {
    let fail = |reason| WeightedError::InvalidGeometry { cell, reason };
    if !geo.vol_signed[cell].is_finite()
        || geo.vol_signed[cell] == 0.0
        || geo.grads[cell].iter().flatten().any(|v| !v.is_finite())
        || geo.gram[cell].iter().flatten().any(|v| !v.is_finite())
    {
        return Err(fail("nonfinite or degenerate element metrics"));
    }
    let tet = complex.tets[cell];
    for i in 0..4 {
        if tet[i] as usize >= complex.vertex_count || tet[..i].contains(&tet[i]) {
            return Err(fail("tetrahedron vertex is missing or repeated"));
        }
        if degree == 1 {
            for j in i + 1..4 {
                let edge = [tet[i].min(tet[j]), tet[i].max(tet[j])];
                if complex.edges.binary_search(&edge).is_err() {
                    return Err(fail("tetrahedron edge is missing from canonical table"));
                }
            }
        }
        if degree == 2 {
            let mut face = [0; 3];
            let mut next = 0;
            for (j, &vertex) in tet.iter().enumerate() {
                if j != i {
                    face[next] = vertex;
                    next += 1;
                }
            }
            face.sort_unstable();
            if complex.faces.binary_search(&face).is_err() {
                return Err(fail("tetrahedron face is missing from canonical table"));
            }
        }
    }
    Ok(())
}

/// Cancellable deterministic counterpart of COO assembly. Original entry
/// ordinals break sort ties and preserve element-summation order.
pub(crate) fn assemble_triplets(
    n: usize,
    entries: Vec<(usize, usize, f64)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Csr, WeightedError> {
    poll(cancelled)?;
    let mut order = Vec::new();
    order
        .try_reserve_exact(entries.len())
        .map_err(|_| WeightedError::Allocation)?;
    for (i, &(row, col, value)) in entries.iter().enumerate() {
        poll(cancelled)?;
        if row >= n || col >= n {
            return Err(WeightedError::InvalidGeometry {
                cell: i,
                reason: "mass entry lies outside its complex",
            });
        }
        if !value.is_finite() {
            return Err(WeightedError::NonFiniteEntry);
        }
        order.push(i);
    }
    // In-place heapsort: no unbounded library sort or second O(nnz) buffer.
    let key = |i: usize| (entries[i].0, entries[i].1, i);
    for start in (0..order.len() / 2).rev() {
        sift(&mut order, start, &key, cancelled)?;
    }
    for end in (1..order.len()).rev() {
        poll(cancelled)?;
        order.swap(0, end);
        sift(&mut order[..end], 0, &key, cancelled)?;
    }
    let rows = n.checked_add(1).ok_or(WeightedError::SizeOverflow)?;
    let mut row_ptr = Vec::new();
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    row_ptr
        .try_reserve_exact(rows)
        .map_err(|_| WeightedError::Allocation)?;
    col_idx
        .try_reserve_exact(entries.len())
        .map_err(|_| WeightedError::Allocation)?;
    vals.try_reserve_exact(entries.len())
        .map_err(|_| WeightedError::Allocation)?;
    for _ in 0..rows {
        poll(cancelled)?;
        row_ptr.push(0usize);
    }
    let mut i = 0;
    while i < order.len() {
        poll(cancelled)?;
        let (row, col, mut value) = entries[order[i]];
        i += 1;
        while i < order.len() && entries[order[i]].0 == row && entries[order[i]].1 == col {
            poll(cancelled)?;
            value += entries[order[i]].2;
            if !value.is_finite() {
                return Err(WeightedError::NonFiniteEntry);
            }
            i += 1;
        }
        row_ptr[row + 1] += 1;
        col_idx.push(col);
        vals.push(value);
    }
    for row in 0..n {
        poll(cancelled)?;
        row_ptr[row + 1] += row_ptr[row];
    }
    let matrix =
        Csr::try_from_parts_with_checkpoint(n, n, row_ptr, col_idx, vals, || poll(cancelled))?
            .ok_or(WeightedError::InvalidGeometry {
                cell: 0,
                reason: "noncanonical assembled mass matrix",
            })?;
    poll(cancelled)?;
    Ok(matrix)
}

fn sift(
    order: &mut [usize],
    mut root: usize,
    key: &impl Fn(usize) -> (usize, usize, usize),
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), WeightedError> {
    while root < order.len() / 2 {
        poll(cancelled)?;
        let mut child = 2 * root + 1;
        if child + 1 < order.len() && key(order[child]) < key(order[child + 1]) {
            child += 1;
        }
        if key(order[root]) >= key(order[child]) {
            break;
        }
        order.swap(root, child);
        root = child;
    }
    Ok(())
}
