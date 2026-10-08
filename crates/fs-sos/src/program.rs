//! Sum-of-squares programs: polynomial identities over decision polynomials,
//! lowered to the block SDP of [`crate::sdp`].
//!
//! A program has decision variables of three kinds:
//!
//! - **scalars** (e.g. the bound `γ` being maximized);
//! - **free polynomials** with a declared monomial basis (equality-constraint
//!   multipliers);
//! - **SOS polynomials** `σ = m(x)ᵀ Q m(x)`, `Q ⪰ 0`, with a declared basis
//!   `m(x)` — each becomes one PSD Gram block.
//!
//! and a list of **identities** `constant + Σ_k multiplier_k · decision_k ≡ 0`
//! (coefficient-wise, as polynomials). Constants and multipliers are interval
//! polynomials ([`IPoly`]), so a derived datum (a Lie derivative, say) is
//! carried with an exact enclosure for [`crate::verify`]; the SDP sees only
//! midpoints. The objective maximizes a weighted sum of scalars.
//!
//! [`SosProgram::solve`] solves the optimization problem;
//! [`SosProgram::solve_centered`] pins chosen scalars and instead maximizes a
//! uniform margin `t` with every Gram block `Q ⪰ t·I` — the well-centred
//! solution a rigorous positivity check needs (an optimal Gram matrix sits on
//! the boundary of the PSD cone and cannot be certified positive definite).

use std::collections::BTreeMap;

use crate::mpoly::{IPoly, MPoly, Monomial};
use crate::sdp::{self, SdpProblem, SdpRow, SdpSettings, SdpSolution, SdpStatus};

/// Handle to a decision variable of an [`SosProgram`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DecisionId(pub(crate) usize);

impl DecisionId {
    /// The decision's index in the program.
    #[must_use]
    pub fn index(self) -> usize {
        self.0
    }
}

/// What a decision variable is.
#[derive(Debug, Clone, PartialEq)]
pub enum DecisionKind {
    /// A real scalar.
    Scalar,
    /// A polynomial with free coefficients over the given basis.
    Free(Vec<Monomial>),
    /// A sum-of-squares polynomial `m(x)ᵀ Q m(x)` over the given basis.
    Sos(Vec<Monomial>),
}

/// A named decision variable.
#[derive(Debug, Clone, PartialEq)]
pub struct Decision {
    /// Human-readable name (diagnostics only).
    pub name: String,
    /// Kind and basis.
    pub kind: DecisionKind,
}

/// One polynomial identity `constant + Σ multiplier·decision ≡ 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct Identity {
    pub(crate) constant: IPoly,
    pub(crate) terms: Vec<(IPoly, DecisionId)>,
}

impl Identity {
    /// Start an identity with its decision-free part.
    #[must_use]
    pub fn new(constant: impl Into<IPoly>) -> Identity {
        Identity {
            constant: constant.into(),
            terms: Vec::new(),
        }
    }

    /// Add `multiplier · decision`.
    #[must_use]
    pub fn term(mut self, multiplier: impl Into<IPoly>, decision: DecisionId) -> Identity {
        self.terms.push((multiplier.into(), decision));
        self
    }

    /// The decision-free part.
    #[must_use]
    pub fn constant(&self) -> &IPoly {
        &self.constant
    }

    /// The `(multiplier, decision)` terms.
    #[must_use]
    pub fn terms(&self) -> &[(IPoly, DecisionId)] {
        &self.terms
    }
}

impl From<MPoly> for IPoly {
    fn from(p: MPoly) -> IPoly {
        p.to_interval()
    }
}

impl From<&MPoly> for IPoly {
    fn from(p: &MPoly) -> IPoly {
        p.to_interval()
    }
}

/// Structured refusals of the SOS layer.
#[derive(Debug, Clone, PartialEq)]
pub enum SosError {
    /// A monomial of an identity's constant cannot be produced by any
    /// decision term: the identity is infeasible as stated. Repair: enlarge
    /// the bases of the decision polynomials.
    UnmatchedMonomial {
        /// Identity index.
        identity: usize,
        /// The monomial's exponents.
        monomial: Vec<u32>,
    },
    /// A basis was empty or contained a monomial of the wrong arity.
    BadBasis {
        /// Decision name.
        decision: String,
    },
    /// A decision referenced by an identity/objective is of the wrong kind
    /// or does not exist.
    BadDecision {
        /// Description.
        what: String,
    },
    /// The lowered SDP was malformed (internal invariant).
    Sdp(sdp::SdpError),
}

impl core::fmt::Display for SosError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SosError::UnmatchedMonomial { identity, monomial } => write!(
                f,
                "identity {identity}: monomial {monomial:?} of the constant is not reachable \
                 by any decision term; enlarge the decision bases"
            ),
            SosError::BadBasis { decision } => {
                write!(
                    f,
                    "decision `{decision}`: empty basis or wrong monomial arity"
                )
            }
            SosError::BadDecision { what } => write!(f, "bad decision reference: {what}"),
            SosError::Sdp(e) => write!(f, "lowered SDP rejected: {e}"),
        }
    }
}

impl std::error::Error for SosError {}

/// The numerical value of one decision in a solution.
#[derive(Debug, Clone, PartialEq)]
pub enum DecisionValue {
    /// A scalar value.
    Scalar(f64),
    /// A free polynomial.
    Free(MPoly),
    /// A Gram matrix (dense row-major, symmetric) over the decision's basis.
    Sos {
        /// The basis.
        basis: Vec<Monomial>,
        /// `basis.len()²` entries, row-major.
        gram: Vec<f64>,
    },
}

/// A numerical solution of an [`SosProgram`] (NOT a certificate — see
/// [`crate::verify::verify`]).
#[derive(Debug, Clone, PartialEq)]
pub struct SosSolution {
    /// Underlying SDP status and diagnostics.
    pub sdp: SdpSolution,
    /// One value per decision, indexed by [`DecisionId::index`].
    pub values: Vec<DecisionValue>,
    /// The uniform Gram margin `t` when solved centred.
    pub margin: Option<f64>,
    rows: Vec<(usize, Monomial)>,
}

impl SosSolution {
    /// The SDP termination status.
    #[must_use]
    pub fn status(&self) -> SdpStatus {
        self.sdp.status
    }

    /// The value of a scalar decision.
    #[must_use]
    pub fn scalar(&self, id: DecisionId) -> Option<f64> {
        match self.values.get(id.0)? {
            DecisionValue::Scalar(v) => Some(*v),
            _ => None,
        }
    }

    /// The dual multiplier of the coefficient row of `monomial` in
    /// `identity` — for the canonical `p − γ − σ₀ ≡ 0` form these are the
    /// Lasserre MOMENTS `y_α ≈ ∫ x^α dμ`.
    #[must_use]
    pub fn dual(&self, identity: usize, monomial: &Monomial) -> Option<f64> {
        self.rows
            .iter()
            .position(|(k, m)| *k == identity && m == monomial)
            .map(|r| self.sdp.y[r])
    }
}

/// One coefficient row of a lowered identity under construction.
#[derive(Default)]
struct RowAcc {
    entries: BTreeMap<(usize, usize, usize), f64>,
    free: BTreeMap<usize, f64>,
    rhs: f64,
}

/// A sum-of-squares program (see the module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct SosProgram {
    nvars: usize,
    decisions: Vec<Decision>,
    identities: Vec<Identity>,
    objective: Vec<(DecisionId, f64)>,
}

impl SosProgram {
    /// An empty program over `nvars` variables.
    #[must_use]
    pub fn new(nvars: usize) -> SosProgram {
        SosProgram {
            nvars,
            decisions: Vec::new(),
            identities: Vec::new(),
            objective: Vec::new(),
        }
    }

    /// Number of polynomial variables.
    #[must_use]
    pub fn nvars(&self) -> usize {
        self.nvars
    }

    /// The decisions.
    #[must_use]
    pub fn decisions(&self) -> &[Decision] {
        &self.decisions
    }

    /// The identities.
    #[must_use]
    pub fn identities(&self) -> &[Identity] {
        &self.identities
    }

    fn push(&mut self, name: &str, kind: DecisionKind) -> Result<DecisionId, SosError> {
        if let DecisionKind::Free(b) | DecisionKind::Sos(b) = &kind {
            let mut seen = std::collections::BTreeSet::new();
            if b.is_empty()
                || b.iter()
                    .any(|m| m.nvars() != self.nvars || !seen.insert(m.clone()))
            {
                return Err(SosError::BadBasis {
                    decision: name.to_string(),
                });
            }
        }
        self.decisions.push(Decision {
            name: name.to_string(),
            kind,
        });
        Ok(DecisionId(self.decisions.len() - 1))
    }

    /// Declare a scalar decision.
    pub fn scalar(&mut self, name: &str) -> DecisionId {
        self.push(name, DecisionKind::Scalar)
            .unwrap_or_else(|_| unreachable!("scalars have no basis"))
    }

    /// Declare a free polynomial over `basis`.
    ///
    /// # Errors
    /// [`SosError::BadBasis`] for an empty/duplicated/wrong-arity basis.
    pub fn free_poly(&mut self, name: &str, basis: Vec<Monomial>) -> Result<DecisionId, SosError> {
        self.push(name, DecisionKind::Free(basis))
    }

    /// Declare an SOS polynomial `m(x)ᵀ Q m(x)` over `basis`.
    ///
    /// # Errors
    /// [`SosError::BadBasis`] for an empty/duplicated/wrong-arity basis.
    pub fn sos(&mut self, name: &str, basis: Vec<Monomial>) -> Result<DecisionId, SosError> {
        self.push(name, DecisionKind::Sos(basis))
    }

    /// Add an identity; returns its index.
    ///
    /// # Errors
    /// [`SosError::BadDecision`] if a term references an unknown decision or
    /// a polynomial has the wrong arity.
    pub fn add_identity(&mut self, identity: Identity) -> Result<usize, SosError> {
        if identity.constant.nvars() != self.nvars {
            return Err(SosError::BadDecision {
                what: "identity constant has the wrong arity".into(),
            });
        }
        for (mult, d) in &identity.terms {
            if d.0 >= self.decisions.len() || mult.nvars() != self.nvars {
                return Err(SosError::BadDecision {
                    what: format!("term on decision {} (or its multiplier arity)", d.0),
                });
            }
        }
        self.identities.push(identity);
        Ok(self.identities.len() - 1)
    }

    /// Maximize `weight · scalar` (accumulates across calls).
    ///
    /// # Errors
    /// [`SosError::BadDecision`] if `id` is not a scalar.
    pub fn maximize(&mut self, id: DecisionId, weight: f64) -> Result<(), SosError> {
        match self.decisions.get(id.0).map(|d| &d.kind) {
            Some(DecisionKind::Scalar) => {
                self.objective.push((id, weight));
                Ok(())
            }
            _ => Err(SosError::BadDecision {
                what: format!("objective decision {} is not a scalar", id.0),
            }),
        }
    }

    /// Solve the program, maximizing the declared objective.
    ///
    /// # Errors
    /// [`SosError`] for structurally infeasible or malformed programs.
    /// Numerical non-convergence is reported in the solution status.
    pub fn solve(&self, settings: &SdpSettings) -> Result<SosSolution, SosError> {
        self.lower_and_solve(&[], false, settings)
    }

    /// Pin the listed scalars to the given values and maximize the uniform
    /// Gram margin `t` (every SOS Gram block `Q ⪰ t·I`, `t ≤ t_cap`). A
    /// positive margin is what lets [`crate::verify::verify`] prove positive
    /// definiteness with interval arithmetic.
    ///
    /// # Errors
    /// As [`SosProgram::solve`].
    pub fn solve_centered(
        &self,
        fixed: &[(DecisionId, f64)],
        settings: &SdpSettings,
    ) -> Result<SosSolution, SosError> {
        for (id, _) in fixed {
            if !matches!(
                self.decisions.get(id.0).map(|d| &d.kind),
                Some(DecisionKind::Scalar)
            ) {
                return Err(SosError::BadDecision {
                    what: format!("pinned decision {} is not a scalar", id.0),
                });
            }
        }
        self.lower_and_solve(fixed, true, settings)
    }

    #[allow(clippy::too_many_lines)]
    fn lower_and_solve(
        &self,
        fixed: &[(DecisionId, f64)],
        centered: bool,
        settings: &SdpSettings,
    ) -> Result<SosSolution, SosError> {
        // Variable layout.
        let mut block_of: Vec<Option<usize>> = vec![None; self.decisions.len()];
        let mut free_of: Vec<Option<usize>> = vec![None; self.decisions.len()];
        let mut block_sizes = Vec::new();
        let mut n_free = 0usize;
        let pinned: BTreeMap<usize, f64> = fixed.iter().map(|(d, v)| (d.0, *v)).collect();
        for (k, d) in self.decisions.iter().enumerate() {
            match &d.kind {
                DecisionKind::Scalar => {
                    if !pinned.contains_key(&k) {
                        free_of[k] = Some(n_free);
                        n_free += 1;
                    }
                }
                DecisionKind::Free(b) => {
                    free_of[k] = Some(n_free);
                    n_free += b.len();
                }
                DecisionKind::Sos(b) => {
                    block_of[k] = Some(block_sizes.len());
                    block_sizes.push(b.len());
                }
            }
        }
        let t_index = if centered {
            n_free += 1;
            Some(n_free - 1)
        } else {
            None
        };
        let cap_block = if centered {
            block_sizes.push(1);
            Some(block_sizes.len() - 1)
        } else {
            None
        };

        let mut rows: BTreeMap<(usize, Monomial), RowAcc> = BTreeMap::new();
        let mut t_cap = 1.0f64;
        for (ki, ident) in self.identities.iter().enumerate() {
            for (m, c) in ident.constant.terms() {
                let mid = c.midpoint();
                t_cap = t_cap.max(mid.abs());
                rows.entry((ki, m.clone())).or_default().rhs -= mid;
            }
            for (mult, d) in &ident.terms {
                let mult = mult.midpoint();
                match &self.decisions[d.0].kind {
                    DecisionKind::Scalar => {
                        if let Some(&v) = pinned.get(&d.0) {
                            for (m, g) in mult.terms() {
                                rows.entry((ki, m.clone())).or_default().rhs -= g * v;
                            }
                        } else if let Some(f) = free_of[d.0] {
                            for (m, g) in mult.terms() {
                                *rows
                                    .entry((ki, m.clone()))
                                    .or_default()
                                    .free
                                    .entry(f)
                                    .or_insert(0.0) += g;
                            }
                        }
                    }
                    DecisionKind::Free(basis) => {
                        let f0 =
                            free_of[d.0].unwrap_or_else(|| unreachable!("free poly has columns"));
                        for (m, g) in mult.terms() {
                            for (k, beta) in basis.iter().enumerate() {
                                *rows
                                    .entry((ki, m.mul(beta)))
                                    .or_default()
                                    .free
                                    .entry(f0 + k)
                                    .or_insert(0.0) += g;
                            }
                        }
                    }
                    DecisionKind::Sos(basis) => {
                        let b = block_of[d.0].unwrap_or_else(|| unreachable!("sos has a block"));
                        for (m, g) in mult.terms() {
                            for i in 0..basis.len() {
                                for j in i..basis.len() {
                                    let alpha = m.mul(&basis[i].mul(&basis[j]));
                                    let acc = rows.entry((ki, alpha)).or_default();
                                    *acc.entries.entry((b, i, j)).or_insert(0.0) += g;
                                    if i == j
                                        && let Some(t) = t_index
                                    {
                                        *acc.free.entry(t).or_insert(0.0) += g;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let mut problem = SdpProblem::new(block_sizes, n_free);
        let mut row_keys = Vec::new();
        for ((ki, m), acc) in rows {
            let entries: Vec<_> = acc
                .entries
                .into_iter()
                .filter(|(_, v)| *v != 0.0)
                .map(|((b, i, j), v)| (b, i, j, v))
                .collect();
            let free: Vec<_> = acc.free.into_iter().filter(|(_, v)| *v != 0.0).collect();
            if entries.is_empty() && free.is_empty() {
                if acc.rhs != 0.0 {
                    return Err(SosError::UnmatchedMonomial {
                        identity: ki,
                        monomial: m.exponents().to_vec(),
                    });
                }
                continue;
            }
            problem.add_row(SdpRow {
                entries,
                free,
                rhs: acc.rhs,
            });
            row_keys.push((ki, m));
        }
        if let (Some(t), Some(cb)) = (t_index, cap_block) {
            problem.add_row(SdpRow {
                entries: vec![(cb, 0, 0, 1.0)],
                free: vec![(t, 1.0)],
                rhs: t_cap,
            });
            problem.set_free_cost(t, -1.0);
        } else {
            for (id, w) in &self.objective {
                if let Some(f) = free_of[id.0] {
                    problem.set_free_cost(f, -w);
                }
            }
        }
        let sol = sdp::solve(&problem, settings).map_err(SosError::Sdp)?;
        let margin = t_index.map(|t| sol.x_free[t]);
        let shift = margin.unwrap_or(0.0);
        let values = self
            .decisions
            .iter()
            .enumerate()
            .map(|(k, d)| match &d.kind {
                DecisionKind::Scalar => DecisionValue::Scalar(
                    pinned
                        .get(&k)
                        .copied()
                        .or_else(|| free_of[k].map(|f| sol.x_free[f]))
                        .unwrap_or(0.0),
                ),
                DecisionKind::Free(basis) => {
                    let f0 = free_of[k].unwrap_or(0);
                    let mut p = MPoly::zero(self.nvars);
                    for (i, beta) in basis.iter().enumerate() {
                        p.add_term(beta.clone(), sol.x_free[f0 + i]);
                    }
                    DecisionValue::Free(p)
                }
                DecisionKind::Sos(basis) => {
                    let b = block_of[k].unwrap_or(0);
                    let n = basis.len();
                    let mut gram = sol.x[b].clone();
                    for i in 0..n {
                        gram[i * n + i] += shift;
                    }
                    DecisionValue::Sos {
                        basis: basis.clone(),
                        gram,
                    }
                }
            })
            .collect();
        Ok(SosSolution {
            sdp: sol,
            values,
            margin,
            rows: row_keys,
        })
    }
}
