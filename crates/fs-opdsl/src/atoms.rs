//! Linear atoms: the materialized building blocks expressions apply.
//! FEEC-derived atoms (exterior derivatives, Whitney masses) are
//! constructed from a complex; `External` atoms are the escape hatch
//! for operators that do not fit the derivation (they must supply a
//! transpose story and are marked `hand` in every report).

use crate::expr::Space;
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_sparse::{Csr, ops};

/// Registry index of an atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomId(pub(crate) usize);

/// How an atom's transpose is obtained.
pub enum Transpose {
    /// The matrix is symmetric: Aᵀ = A.
    Symmetric,
    /// Derive by explicit CSR transposition at registration.
    Derived,
    /// Hand-supplied transpose (must match — the consistency gate
    /// checks it like everything else).
    Explicit(Csr),
}

/// What kind of atom this is (drives dd = 0 simplification and the
/// derived/hand provenance column of the plan report).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomKind {
    /// Exterior derivative d_k (k → k+1), possibly transposed.
    D {
        /// Source degree of the UNTRANSPOSED operator.
        k: u8,
        /// Applied as dᵀ.
        transposed: bool,
    },
    /// Whitney mass / Galerkin star M_k (symmetric).
    Mass {
        /// Degree.
        k: u8,
    },
    /// Hand-supplied matrix.
    External,
}

/// A materialized linear atom.
pub struct Atom {
    /// Display name.
    pub name: String,
    /// Kind (provenance + algebraic identity hooks).
    pub kind: AtomKind,
    /// Input space.
    pub in_space: Space,
    /// Output space.
    pub out_space: Space,
    pub(crate) mat: Csr,
    pub(crate) mat_t: Option<Csr>,
    pub(crate) symmetric: bool,
}

impl Atom {
    /// (k, transposed) when this atom is an exterior derivative.
    #[must_use]
    pub fn d_signature(&self) -> Option<(u8, bool)> {
        match self.kind {
            AtomKind::D { k, transposed } => Some((k, transposed)),
            _ => None,
        }
    }

    /// The exterior derivative d_k of a complex (integer incidence
    /// materialized as ±1.0 CSR; exact).
    #[must_use]
    pub fn d(complex: &TetComplex, k: u8, dims: Dims) -> Atom {
        let inc = match k {
            0 => complex.d0(),
            1 => complex.d1(),
            2 => complex.d2(),
            _ => panic!("d_k exists for k = 0..=2 on a 3-complex"),
        };
        let mat = fs_feec::incidence_to_csr(&inc);
        let mat_t = ops::transpose(&mat);
        let in_space = Space {
            degree: k,
            n: mat.ncols(),
            dims,
        };
        let out_space = Space {
            degree: k + 1,
            n: mat.nrows(),
            dims,
        };
        Atom {
            name: format!("d{k}"),
            kind: AtomKind::D {
                k,
                transposed: false,
            },
            in_space,
            out_space,
            mat,
            mat_t: Some(mat_t),
            symmetric: false,
        }
    }

    /// The transpose dᵀ_k (maps degree k+1 → k) as a forward atom.
    #[must_use]
    pub fn d_transposed(complex: &TetComplex, k: u8, dims: Dims) -> Atom {
        let base = Atom::d(complex, k, dims);
        Atom {
            name: format!("d{k}T"),
            kind: AtomKind::D {
                k,
                transposed: true,
            },
            in_space: base.out_space,
            out_space: base.in_space,
            mat: base.mat_t.expect("built above"),
            mat_t: Some(base.mat),
            symmetric: false,
        }
    }

    /// The Whitney mass matrix M_k (Galerkin Hodge star; symmetric).
    #[must_use]
    pub fn mass(complex: &TetComplex, geo: &fs_feec::ElementGeometry, k: u8, dims: Dims) -> Atom {
        let mat = fs_feec::mass_matrix(complex, geo, k);
        let n = mat.nrows();
        let space = Space { degree: k, n, dims };
        Atom {
            name: format!("M{k}"),
            kind: AtomKind::Mass { k },
            in_space: space,
            out_space: space,
            mat,
            mat_t: None,
            symmetric: true,
        }
    }

    /// A hand-supplied matrix (the escape hatch). `degree` may be the
    /// raw-vector marker 255 for non-cochain dof spaces.
    #[must_use]
    pub fn external(
        name: &str,
        mat: Csr,
        transpose: Transpose,
        in_space: Space,
        out_space: Space,
    ) -> Atom {
        let (mat_t, symmetric) = match transpose {
            Transpose::Symmetric => (None, true),
            Transpose::Derived => (Some(ops::transpose(&mat)), false),
            Transpose::Explicit(t) => (Some(t), false),
        };
        Atom {
            name: name.to_string(),
            kind: AtomKind::External,
            in_space,
            out_space,
            mat,
            mat_t,
            symmetric,
        }
    }

    /// Forward application y = A x.
    pub(crate) fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.mat.spmv(x, y);
    }

    /// Transpose application y = Aᵀ x.
    pub(crate) fn apply_t(&self, x: &[f64], y: &mut [f64]) {
        if self.symmetric {
            self.mat.spmv(x, y);
        } else {
            self.mat_t
                .as_ref()
                .expect("non-symmetric atoms carry a transpose")
                .spmv(x, y);
        }
    }

    /// The materialized matrix (for plan folding and hand comparisons).
    #[must_use]
    pub fn matrix(&self) -> &Csr {
        &self.mat
    }

    /// True for hand-supplied atoms (the report provenance column).
    #[must_use]
    pub fn is_hand(&self) -> bool {
        matches!(self.kind, AtomKind::External)
    }
}
