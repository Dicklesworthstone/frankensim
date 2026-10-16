//! Lowering: one typed expression generates the primal apply, the
//! JVP, the VJP / discrete adjoint, the materialized single-matrix
//! form (a compiler optimization: constant linear chains folded
//! through deterministic spgemm), the DWR indicators, and the plan
//! report (provenance, block hints, tile-kernel metadata). The
//! transpose of every derived piece is DERIVED — never hand-kept in
//! sync — which is the whole point.

use crate::expr::{Expr, OperatorDef, Space};
use fs_sparse::{Csr, ops};
use std::fmt::Write as _;

/// A lowered operator: everything generated from one expression.
pub struct LoweredOperator<'d> {
    def: &'d OperatorDef,
    expr: Expr,
    /// Cached pointwise-law input values from the last `linearize`
    /// (the chain-rule state; one slice per Pointwise node in
    /// pre-order).
    linearization: Option<Vec<Vec<f64>>>,
}

impl OperatorDef {
    /// Lower an expression (the residual's u-dependent part).
    ///
    /// # Panics
    /// If the expression is ill-typed against this registry (cannot
    /// happen for trees built through the checked constructors).
    #[must_use]
    pub fn lower(&self, expr: Expr) -> LoweredOperator<'_> {
        self.space_of(&expr).expect("lower: ill-typed expression");
        LoweredOperator {
            def: self,
            expr,
            linearization: None,
        }
    }
}

/// Each recursion step allocates its output vector; fixture-scale v1
/// keeps this simple (allocation-light fused apply is the
/// tilelang-fusion lane's concern, recorded in the CONTRACT).
fn zeros(n: usize) -> Vec<f64> {
    vec![0.0f64; n]
}

impl LoweredOperator<'_> {
    /// The expression this operator was generated from.
    #[must_use]
    pub fn expr(&self) -> &Expr {
        &self.expr
    }

    /// Output space of the operator.
    #[must_use]
    pub fn out_space(&self) -> Space {
        self.def.space_of(&self.expr).expect("typed")
    }

    /// PRIMAL: y = R(u).
    #[must_use]
    pub fn apply(&self, u: &[f64]) -> Vec<f64> {
        self.eval(&self.expr, u, &mut None)
    }

    /// Linearize at `u0`: caches every pointwise node's input values
    /// so subsequent [`Self::jvp`]/[`Self::vjp`] differentiate THERE.
    /// Purely linear operators need no linearization.
    pub fn linearize(&mut self, u0: &[f64]) {
        let mut cache = Some(Vec::new());
        let _ = self.eval(&self.expr, u0, &mut cache);
        self.linearization = cache;
    }

    /// JVP: y = J|_(u0) · v (chain rule forward). For linear
    /// operators J = R and no linearization is needed.
    ///
    /// # Panics
    /// If the expression has pointwise nodes and `linearize` was not
    /// called.
    #[must_use]
    pub fn jvp(&self, v: &[f64]) -> Vec<f64> {
        self.jvp_rec(&self.expr, v, 0)
    }

    /// VJP / DISCRETE ADJOINT: x = Jᵀ|_(u0) · w (chain rule reverse,
    /// transposing every factor). Same linearization state as
    /// [`Self::jvp`] — the two are generated from the same tree, so
    /// transpose consistency is structural.
    ///
    /// # Panics
    /// If the expression has pointwise nodes and `linearize` was not
    /// called.
    #[must_use]
    pub fn vjp(&self, w: &[f64]) -> Vec<f64> {
        let mut grad = zeros(self.def.field_space.n);
        self.vjp_rec(&self.expr, w, &mut grad, 0);
        grad
    }

    /// Number of pointwise nodes in a subtree — the linearization
    /// cache is numbered in POST-ORDER (the order `eval` pushes), and
    /// both derivative walks compute each node's slot structurally:
    /// slot(node) = base + count(child subtree). This stays correct
    /// for nested laws where a shared mutable cursor would not (the
    /// reverse walk visits nodes in the opposite order).
    fn count_pointwise(e: &Expr) -> usize {
        match e {
            Expr::Field(_) | Expr::Zero(_) => 0,
            Expr::Apply(_, x) | Expr::Scale(_, x) => Self::count_pointwise(x),
            Expr::Add(a, b) => Self::count_pointwise(a) + Self::count_pointwise(b),
            Expr::Pointwise(_, x) => 1 + Self::count_pointwise(x),
        }
    }

    fn eval(&self, e: &Expr, u: &[f64], cache: &mut Option<Vec<Vec<f64>>>) -> Vec<f64> {
        match e {
            Expr::Field(_) => u.to_vec(),
            Expr::Zero(s) => zeros(s.n),
            Expr::Apply(id, x) => {
                let xv = self.eval(x, u, cache);
                let atom = &self.def.atoms[id.0];
                let mut y = zeros(atom.out_space.n);
                atom.apply(&xv, &mut y);
                y
            }
            Expr::Scale(s, x) => {
                let xv = self.eval(x, u, cache);
                let mut y = zeros(xv.len());
                crate::kernels::scale_k::run(&xv, *s, &mut y);
                y
            }
            Expr::Add(a, b) => {
                let av = self.eval(a, u, cache);
                let bv = self.eval(b, u, cache);
                let mut y = zeros(av.len());
                crate::kernels::add_k::run(&av, &bv, &mut y);
                y
            }
            Expr::Pointwise(id, x) => {
                let xv = self.eval(x, u, cache);
                let law = &self.def.laws[id.0];
                if let Some(c) = cache {
                    c.push(xv.clone());
                }
                xv.iter().map(|&v| law.value(v)).collect()
            }
        }
    }

    fn jvp_rec(&self, e: &Expr, v: &[f64], base: usize) -> Vec<f64> {
        match e {
            Expr::Field(_) => v.to_vec(),
            Expr::Zero(s) => zeros(s.n),
            Expr::Apply(id, x) => {
                let xv = self.jvp_rec(x, v, base);
                let atom = &self.def.atoms[id.0];
                let mut y = zeros(atom.out_space.n);
                atom.apply(&xv, &mut y);
                y
            }
            Expr::Scale(s, x) => {
                let xv = self.jvp_rec(x, v, base);
                let mut y = zeros(xv.len());
                crate::kernels::scale_k::run(&xv, *s, &mut y);
                y
            }
            Expr::Add(a, b) => {
                let av = self.jvp_rec(a, v, base);
                let bv = self.jvp_rec(b, v, base + Self::count_pointwise(a));
                let mut y = zeros(av.len());
                crate::kernels::add_k::run(&av, &bv, &mut y);
                y
            }
            Expr::Pointwise(id, x) => {
                let xv = self.jvp_rec(x, v, base);
                let slot = base + Self::count_pointwise(x);
                let u0 = &self
                    .linearization
                    .as_ref()
                    .expect("pointwise JVP requires linearize(u0)")[slot];
                let law = &self.def.laws[id.0];
                xv.iter()
                    .zip(u0)
                    .map(|(&t, &p)| law.derivative(p) * t)
                    .collect()
            }
        }
    }

    fn vjp_rec(&self, e: &Expr, w: &[f64], grad: &mut [f64], base: usize) {
        match e {
            Expr::Field(_) => {
                crate::kernels::accum_k::run(w, 1.0, grad);
            }
            Expr::Zero(_) => {}
            Expr::Apply(id, x) => {
                let atom = &self.def.atoms[id.0];
                let mut back = zeros(atom.in_space.n);
                atom.apply_t(w, &mut back);
                self.vjp_rec(x, &back, grad, base);
            }
            Expr::Scale(s, x) => {
                let mut back = zeros(w.len());
                crate::kernels::scale_k::run(w, *s, &mut back);
                self.vjp_rec(x, &back, grad, base);
            }
            Expr::Add(a, b) => {
                self.vjp_rec(a, w, grad, base);
                self.vjp_rec(b, w, grad, base + Self::count_pointwise(a));
            }
            Expr::Pointwise(id, x) => {
                let slot = base + Self::count_pointwise(x);
                let u0 = &self
                    .linearization
                    .as_ref()
                    .expect("pointwise VJP requires linearize(u0)")[slot];
                let law = &self.def.laws[id.0];
                let back: Vec<f64> = w
                    .iter()
                    .zip(u0)
                    .map(|(&t, &p)| law.derivative(p) * t)
                    .collect();
                self.vjp_rec(x, &back, grad, base);
            }
        }
    }

    /// MATERIALIZE: fold the expression into a single CSR when it is
    /// purely linear (constant chains composed through deterministic
    /// transpose/spgemm — the same association fs-feec's `stiffness`
    /// uses, so a Poisson chain reproduces it bitwise). Returns `None`
    /// when the expression contains pointwise nonlinearities.
    #[must_use]
    pub fn materialize(&self) -> Option<Csr> {
        self.mat_rec(&self.expr)
    }

    fn mat_rec(&self, e: &Expr) -> Option<Csr> {
        match e {
            Expr::Field(s) => Some(Csr::identity(s.n)),
            // A surviving Zero can only be the WHOLE expression (the
            // checked constructors fold zero summands); its matrix is
            // the empty operator.
            Expr::Zero(s) => Some(Csr::from_parts(
                s.n,
                self.def.field_space.n,
                vec![0; s.n + 1],
                vec![],
                vec![],
            )),
            Expr::Apply(id, x) => {
                let atom = &self.def.atoms[id.0];
                match x.as_ref() {
                    Expr::Field(_) => Some(atom.matrix().clone()),
                    other => {
                        let inner = self.mat_rec(other)?;
                        Some(ops::spgemm(atom.matrix(), &inner))
                    }
                }
            }
            Expr::Scale(s, x) => {
                let m = self.mat_rec(x)?;
                Some(scale_csr(&m, *s))
            }
            Expr::Add(a, b) => {
                let (ma, mb) = (self.mat_rec(a)?, self.mat_rec(b)?);
                Some(add_csr(&ma, &mb))
            }
            Expr::Pointwise(..) => None,
        }
    }

    /// Generation report: provenance (derived vs hand atoms), block
    /// hints (which spaces couple, symmetry — the single-field v1
    /// preconditioner hint), the tile-kernel metadata the executor
    /// uses, and a structural fingerprint. Deterministic: built from
    /// the registry and tree only.
    #[must_use]
    pub fn report(&self) -> PlanReport {
        let mut atoms = Vec::new();
        for a in &self.def.atoms {
            atoms.push(format!(
                "{{\"name\":\"{}\",\"kind\":\"{:?}\",\"provenance\":\"{}\",\"symmetric\":{},\
                 \"in_degree\":{},\"out_degree\":{},\"nnz\":{}}}",
                a.name,
                a.kind,
                if a.is_hand() { "hand" } else { "derived" },
                a.symmetric,
                a.in_space.degree,
                a.out_space.degree,
                a.matrix().nnz()
            ));
        }
        let laws: Vec<String> = self
            .def
            .laws
            .iter()
            .map(|l| format!("{{\"name\":\"{}\",\"provenance\":\"law\"}}", l.name()))
            .collect();
        let kernels = vec![
            crate::kernels::scale_k::META.descr(),
            crate::kernels::add_k::META.descr(),
            crate::kernels::accum_k::META.descr(),
        ];
        PlanReport {
            atoms,
            laws,
            kernels,
            structure: structure_of(&self.expr, self.def),
        }
    }
}

/// The generation report (see [`LoweredOperator::report`]).
pub struct PlanReport {
    /// One JSON line per atom (provenance column included).
    pub atoms: Vec<String>,
    /// One JSON line per pointwise law.
    pub laws: Vec<String>,
    /// Tile-kernel metadata lines (fs-tilelang META descriptors).
    pub kernels: Vec<String>,
    /// Structural fingerprint of the expression tree.
    pub structure: String,
}

impl PlanReport {
    /// The whole report as one JSON object (deterministic).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\"atoms\":[");
        s.push_str(&self.atoms.join(","));
        s.push_str("],\"laws\":[");
        s.push_str(&self.laws.join(","));
        s.push_str("],\"kernels\":[");
        s.push_str(&self.kernels.join(","));
        let _ = write!(s, "],\"structure\":\"{}\"}}", self.structure);
        s
    }
}

fn structure_of(e: &Expr, def: &OperatorDef) -> String {
    match e {
        Expr::Field(_) => "u".to_string(),
        Expr::Zero(_) => "0".to_string(),
        Expr::Apply(id, x) => {
            format!("{}({})", def.atoms[id.0].name, structure_of(x, def))
        }
        Expr::Scale(s, x) => format!("{s}*{}", structure_of(x, def)),
        Expr::Add(a, b) => {
            format!("({} + {})", structure_of(a, def), structure_of(b, def))
        }
        Expr::Pointwise(id, x) => {
            format!("{}({})", def.laws[id.0].name(), structure_of(x, def))
        }
    }
}

/// Algebraic DWR indicators: η_i = |r_i · z_i| (residual weighted by
/// the dual solution, the standard algebraic form — dwr-adaptivity's
/// food at the dof level; element-integrated forms join the
/// higher-order bead).
#[must_use]
pub fn dwr_indicators(residual: &[f64], dual: &[f64]) -> Vec<f64> {
    assert_eq!(
        residual.len(),
        dual.len(),
        "DWR needs matching residual/dual lengths"
    );
    residual
        .iter()
        .zip(dual)
        .map(|(&r, &z)| (r * z).abs())
        .collect()
}

/// α·A with the same sparsity (entry order preserved).
fn scale_csr(a: &Csr, alpha: f64) -> Csr {
    let mut row_ptr = Vec::with_capacity(a.nrows() + 1);
    let mut cols = Vec::with_capacity(a.nnz());
    let mut vals = Vec::with_capacity(a.nnz());
    row_ptr.push(0usize);
    for r in 0..a.nrows() {
        let (rc, rv) = a.row(r);
        for (&c, &v) in rc.iter().zip(rv) {
            cols.push(c);
            vals.push(alpha * v);
        }
        row_ptr.push(cols.len());
    }
    Csr::from_parts(a.nrows(), a.ncols(), row_ptr, cols, vals)
}

/// A + B through the deterministic COO accumulate path.
fn add_csr(a: &Csr, b: &Csr) -> Csr {
    assert_eq!(a.nrows(), b.nrows(), "add_csr shape mismatch");
    assert_eq!(a.ncols(), b.ncols(), "add_csr shape mismatch");
    let mut coo = fs_sparse::Coo::new(a.nrows(), a.ncols());
    for r in 0..a.nrows() {
        let (rc, rv) = a.row(r);
        for (&c, &v) in rc.iter().zip(rv) {
            coo.push(r, c, v);
        }
        let (rc, rv) = b.row(r);
        for (&c, &v) in rc.iter().zip(rv) {
            coo.push(r, c, v);
        }
    }
    coo.assemble()
}
