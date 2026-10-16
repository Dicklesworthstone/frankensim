//! The typed expression tree and its registry. Every node carries a
//! [`Space`] (cochain degree, dof count, Qty dimensions); construction
//! is through checked builder methods that return structured
//! [`TypeError`]s — an ill-typed operator cannot be REPRESENTED, which
//! is the point of a typed IR.

use crate::atoms::{Atom, AtomId};
use crate::law::{LawId, PointwiseLaw};
use fs_qty::Dims;

/// The typed "space" of an expression value: a k-cochain vector with
/// physical dimensions. Two spaces compose only when they agree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Space {
    /// Cochain degree (0..=3) — or the raw vector-dof marker 255 for
    /// non-cochain spaces (e.g. vector elasticity dofs).
    pub degree: u8,
    /// Dof count.
    pub n: usize,
    /// Physical dimensions of the stored values.
    pub dims: Dims,
}

/// Structured type errors: what failed, where, and the two sides.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TypeError {
    /// Atom applied to an argument in the wrong space.
    ApplyMismatch {
        /// The atom's name.
        atom: String,
        /// What the atom expects.
        expected: Space,
        /// What the argument is.
        found: Space,
    },
    /// Sum of two expressions in different spaces.
    AddMismatch {
        /// Left space.
        left: Space,
        /// Right space.
        right: Space,
    },
    /// A registry id that does not exist.
    UnknownId {
        /// Which table.
        what: &'static str,
        /// The offending index.
        id: usize,
    },
    /// A law's dimension scaling left the supported i8 exponent domain
    /// (bead sj31i.11): refused typed instead of clamped.
    DimensionOverflow {
        /// The law's display name.
        law: String,
        /// The input dimensions whose scaling overflowed.
        in_dims: fs_qty::Dims,
    },
}

impl std::fmt::Display for TypeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TypeError::ApplyMismatch {
                atom,
                expected,
                found,
            } => write!(
                f,
                "atom `{atom}` expects degree {} n {} dims {:?}, got degree {} n {} dims {:?}",
                expected.degree, expected.n, expected.dims, found.degree, found.n, found.dims
            ),
            TypeError::AddMismatch { left, right } => write!(
                f,
                "cannot add degree {} n {} dims {:?} to degree {} n {} dims {:?}",
                left.degree, left.n, left.dims, right.degree, right.n, right.dims
            ),
            TypeError::UnknownId { what, id } => write!(f, "unknown {what} id {id}"),
            TypeError::DimensionOverflow { law, in_dims } => write!(
                f,
                "law `{law}` dimension scaling of {in_dims:?} leaves the supported i8 \
                 exponent domain"
            ),
        }
    }
}

/// Expression nodes. Trees are built through [`OperatorDef`] builder
/// methods only, so every constructed node is well-typed by
/// construction.
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    /// The unknown field u.
    Field(Space),
    /// Linear atom application.
    Apply(AtomId, Box<Expr>),
    /// Scalar multiple.
    Scale(f64, Box<Expr>),
    /// Sum.
    Add(Box<Expr>, Box<Expr>),
    /// Pointwise differentiable law applied value-by-value.
    Pointwise(LawId, Box<Expr>),
    /// The zero element of a space (produced by simplification —
    /// d∘d = 0 lives here).
    Zero(Space),
}

/// An operator definition: the atom/law registry plus the residual
/// expression R(u). One definition is the single source of truth for
/// everything the plan generates.
pub struct OperatorDef {
    pub(crate) atoms: Vec<Atom>,
    pub(crate) laws: Vec<Box<dyn PointwiseLaw>>,
    /// The unknown's space.
    pub field_space: Space,
}

impl OperatorDef {
    /// New definition for an unknown living in `field_space`.
    #[must_use]
    pub fn new(field_space: Space) -> OperatorDef {
        OperatorDef {
            atoms: Vec::new(),
            laws: Vec::new(),
            field_space,
        }
    }

    /// Register a linear atom, returning its id.
    pub fn add_atom(&mut self, atom: Atom) -> AtomId {
        self.atoms.push(atom);
        AtomId(self.atoms.len() - 1)
    }

    /// Register a pointwise law, returning its id.
    pub fn add_law(&mut self, law: Box<dyn PointwiseLaw>) -> LawId {
        self.laws.push(law);
        LawId(self.laws.len() - 1)
    }

    /// The unknown field expression.
    #[must_use]
    pub fn field(&self) -> Expr {
        Expr::Field(self.field_space)
    }

    /// The space an expression lives in.
    ///
    /// # Errors
    /// [`TypeError::UnknownId`] if the expression references a
    /// non-registered atom or law.
    pub fn space_of(&self, e: &Expr) -> Result<Space, TypeError> {
        match e {
            Expr::Field(s) | Expr::Zero(s) => Ok(*s),
            Expr::Apply(id, _) => {
                let atom = self.atoms.get(id.0).ok_or(TypeError::UnknownId {
                    what: "atom",
                    id: id.0,
                })?;
                Ok(atom.out_space)
            }
            Expr::Scale(_, x) => self.space_of(x),
            Expr::Add(a, _) => self.space_of(a),
            Expr::Pointwise(id, x) => {
                let law = self.laws.get(id.0).ok_or(TypeError::UnknownId {
                    what: "law",
                    id: id.0,
                })?;
                let s = self.space_of(x)?;
                let dims = law.out_dims(s.dims).ok_or(TypeError::DimensionOverflow {
                    law: law.name().to_string(),
                    in_dims: s.dims,
                })?;
                Ok(Space { dims, ..s })
            }
        }
    }

    /// Checked atom application (with d∘d = 0 simplification: applying
    /// an exterior derivative to the image of the previous one folds
    /// to [`Expr::Zero`] — the exactness identity is honored at the IR
    /// level, before any float exists).
    ///
    /// # Errors
    /// [`TypeError::ApplyMismatch`] on space disagreement.
    pub fn apply(&self, id: AtomId, x: Expr) -> Result<Expr, TypeError> {
        let atom = self.atoms.get(id.0).ok_or(TypeError::UnknownId {
            what: "atom",
            id: id.0,
        })?;
        let found = self.space_of(&x)?;
        if atom.in_space != found {
            return Err(TypeError::ApplyMismatch {
                atom: atom.name.clone(),
                expected: atom.in_space,
                found,
            });
        }
        // dd = 0: D_{k+1} ∘ D_k ≡ 0 (both untransposed, consecutive
        // degrees) and the transposed pair the other way around.
        if let Some((k_outer, t_outer)) = atom.d_signature() {
            let mut probe = &x;
            // Skip scale nodes: d(s·x) = s·d(x), and s·0 = 0.
            while let Expr::Scale(_, inner) = probe {
                probe = inner;
            }
            if let Expr::Apply(inner_id, _) = probe
                && let Some(inner_atom) = self.atoms.get(inner_id.0)
                && let Some((k_inner, t_inner)) = inner_atom.d_signature()
                && t_outer == t_inner
                && ((!t_outer && k_outer == k_inner + 1) || (t_outer && k_inner == k_outer + 1))
            {
                return Ok(Expr::Zero(atom.out_space));
            }
        }
        if matches!(x, Expr::Zero(_)) {
            return Ok(Expr::Zero(atom.out_space));
        }
        Ok(Expr::Apply(id, Box::new(x)))
    }

    /// Checked scalar multiple (folds nested scales and zero).
    #[must_use]
    pub fn scale(&self, s: f64, x: Expr) -> Expr {
        match x {
            Expr::Scale(t, inner) => Expr::Scale(s * t, inner),
            Expr::Zero(sp) => Expr::Zero(sp),
            other => Expr::Scale(s, Box::new(other)),
        }
    }

    /// Checked sum (eliminates zero summands).
    ///
    /// # Errors
    /// [`TypeError::AddMismatch`] if the spaces differ.
    pub fn add(&self, a: Expr, b: Expr) -> Result<Expr, TypeError> {
        let (sa, sb) = (self.space_of(&a)?, self.space_of(&b)?);
        if sa != sb {
            return Err(TypeError::AddMismatch {
                left: sa,
                right: sb,
            });
        }
        Ok(match (a, b) {
            (Expr::Zero(_), b) => b,
            (a, Expr::Zero(_)) => a,
            (a, b) => Expr::Add(Box::new(a), Box::new(b)),
        })
    }

    /// Checked pointwise law application.
    ///
    /// # Errors
    /// [`TypeError::UnknownId`] for unregistered laws.
    pub fn pointwise(&self, id: LawId, x: Expr) -> Result<Expr, TypeError> {
        self.laws.get(id.0).ok_or(TypeError::UnknownId {
            what: "law",
            id: id.0,
        })?;
        // Space is preserved except dims (law-declared); checked lazily
        // through space_of.
        Ok(Expr::Pointwise(id, Box::new(x)))
    }
}
