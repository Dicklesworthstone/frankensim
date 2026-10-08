//! Certified global polynomial optimization (plan §9.8): Lasserre/SOS
//! relaxations whose answer is a rigorous ENCLOSURE of the global minimum.
//!
//! `minimize p(x)` (optionally over `K = {x : g_k(x) ≥ 0}`) is relaxed to
//!
//! ```text
//! maximize γ  s.t.  p − γ − Σ_k g_k σ_k = σ₀,   σ_k SOS
//! ```
//!
//! (Putinar's form; unconstrained when there are no `g_k`). Any feasible `γ`
//! is a lower bound on the minimum. The returned LOWER bound is proved by
//! [`crate::verify::verify`] for a slightly backed-off `γ_c` re-solved with a
//! centred Gram margin. The returned UPPER bound is the interval-enclosed
//! value of `p` at a point extracted from the relaxation's moments (the dual
//! variables are the Lasserre moment sequence `y_α`), polished by Newton
//! steps (unconstrained) and checked feasible with interval arithmetic
//! (constrained). The pair `lower ≤ min p ≤ upper` is the theorem; the gap is
//! reported, never hidden.
//!
//! No-claim: the relaxation order is fixed by the caller (or by the degree);
//! a large gap means "increase the order or the problem is hard", not that
//! the bounds are wrong. Minimizer extraction is a heuristic (mean and
//! eigenvector candidates of the moment matrix + local polish), so the upper
//! bound may be loose when the minimizer set is not a single point.

use fs_ivl::Interval;

use crate::dense;
use crate::mpoly::{MPoly, Monomial, monomials_in_degree_range};
use crate::program::{DecisionId, DecisionValue, Identity, SosError, SosProgram};
use crate::sdp::{SdpSettings, SdpStatus};
use crate::verify::{Certificate, VerifyError, verify};

/// Options for [`minimize`].
#[derive(Debug, Clone, PartialEq)]
pub struct GlobalOptions {
    /// Relaxation half-degree `d` (the SOS multiplier σ₀ has degree `2d`).
    /// `None` uses the smallest admissible `d`.
    pub order: Option<u32>,
    /// SDP controls.
    pub sdp: SdpSettings,
    /// How many back-off attempts (×10 each) the certification may take.
    pub certify_attempts: usize,
    /// Relative back-off of the first certification attempt.
    pub initial_backoff: f64,
}

impl Default for GlobalOptions {
    fn default() -> Self {
        GlobalOptions {
            order: None,
            sdp: SdpSettings::default(),
            certify_attempts: 6,
            initial_backoff: 1e-8,
        }
    }
}

/// Why no certified enclosure was produced.
#[derive(Debug, Clone, PartialEq)]
pub enum GlobalError {
    /// The objective/constraints have inconsistent arity or are constant
    /// zero-variable problems.
    BadInput {
        /// Description.
        what: String,
    },
    /// The SOS program could not be built (structural infeasibility, e.g.
    /// odd leading degree).
    Program(SosError),
    /// The relaxation did not solve (e.g. `p − γ` is not SOS for any `γ` at
    /// this order — the Motzkin phenomenon — or `p` is unbounded below).
    RelaxationFailed {
        /// SDP status of the γ-maximization.
        status: SdpStatus,
    },
    /// The relaxation solved but no backed-off bound could be proved.
    CertificationFailed {
        /// The last verification refusal.
        last: VerifyError,
    },
}

impl core::fmt::Display for GlobalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            GlobalError::BadInput { what } => write!(f, "bad global-optimization input: {what}"),
            GlobalError::Program(e) => write!(f, "{e}"),
            GlobalError::RelaxationFailed { status } => write!(
                f,
                "SOS relaxation did not solve ({status:?}): the objective may be unbounded \
                 below or not SOS at this order; raise `order`"
            ),
            GlobalError::CertificationFailed { last } => {
                write!(f, "no backed-off bound could be proved: {last}")
            }
        }
    }
}

impl std::error::Error for GlobalError {}

/// A certified enclosure of a global minimum.
#[derive(Debug, Clone, PartialEq)]
pub struct GlobalBound {
    /// PROVED: `p(x) ≥ lower` for every `x` (in `K`).
    pub lower: f64,
    /// PROVED (when present): some feasible point attains `p ≤ upper`.
    pub upper: Option<f64>,
    /// The feasible point attaining `upper`.
    pub minimizer: Option<Vec<f64>>,
    /// The relaxation's floating optimum `γ*` (NOT certified).
    pub relaxation_value: f64,
    /// Relaxation half-degree used.
    pub order: u32,
    /// The verified certificate of `lower`.
    pub certificate: Certificate,
    /// The back-off `γ* − lower` that certification needed.
    pub backoff: f64,
}

impl GlobalBound {
    /// `upper − lower` when both exist.
    #[must_use]
    pub fn gap(&self) -> Option<f64> {
        self.upper.map(|u| u - self.lower)
    }
}

/// Basis monomials for a degree-`2d` SOS of a polynomial with support
/// `supp`: all monomials of degree ≤ `d`, pruned by sound half-Newton-
/// polytope tests along the coordinate and total-degree directions (a
/// monomial `m` with `2m` outside the Newton polytope cannot appear in any
/// SOS decomposition).
fn pruned_basis(nvars: usize, d: u32, supp: &[Monomial]) -> Vec<Monomial> {
    let all = monomials_in_degree_range(nvars, 0, d);
    if supp.is_empty() {
        return all;
    }
    let maxk: Vec<u32> = (0..nvars)
        .map(|k| supp.iter().map(|m| m.exponents()[k]).max().unwrap_or(0))
        .collect();
    let mink: Vec<u32> = (0..nvars)
        .map(|k| supp.iter().map(|m| m.exponents()[k]).min().unwrap_or(0))
        .collect();
    let maxd = supp.iter().map(Monomial::degree).max().unwrap_or(0);
    let mind = supp.iter().map(Monomial::degree).min().unwrap_or(0);
    let mut basis: Vec<Monomial> = all
        .into_iter()
        .filter(|m| {
            let e = m.exponents();
            (0..nvars).all(|k| 2 * e[k] <= maxk[k] && 2 * e[k] >= mink[k])
                && 2 * m.degree() <= maxd
                && 2 * m.degree() >= mind
        })
        .collect();
    // Iterative diagonal-consistency pruning: if x^{2m} is absent from the
    // support and is not a product of two DISTINCT basis monomials, the Gram
    // diagonal entry for m is forced to 0, hence (PSD) its whole row — drop
    // m and repeat. Sound, and it removes the forced-singular rows that would
    // make every Gram matrix unprovably definite.
    let supp_set: std::collections::BTreeSet<&Monomial> = supp.iter().collect();
    loop {
        let set: std::collections::BTreeSet<Monomial> = basis.iter().cloned().collect();
        let keep: Vec<Monomial> = basis
            .iter()
            .filter(|m| {
                let sq = m.mul(m);
                if supp_set.contains(&sq) {
                    return true;
                }
                basis.iter().any(|a| {
                    a != *m
                        && sq
                            .checked_div(a)
                            .is_some_and(|b| b != **m && b != *a && set.contains(&b))
                })
            })
            .cloned()
            .collect();
        if keep.len() == basis.len() {
            return basis;
        }
        basis = keep;
    }
}

/// Certified global minimization of `p` over `K = {x : g(x) ≥ 0 ∀ g ∈ constraints}`
/// (all of ℝⁿ when `constraints` is empty).
///
/// # Errors
/// [`GlobalError`] naming the stage that failed. An error is a refusal; it
/// never asserts anything about the true minimum.
#[allow(clippy::too_many_lines)]
pub fn minimize(
    p: &MPoly,
    constraints: &[MPoly],
    opts: &GlobalOptions,
) -> Result<GlobalBound, GlobalError> {
    let n = p.nvars();
    if n == 0 {
        return Err(GlobalError::BadInput {
            what: "objective has no variables".into(),
        });
    }
    if constraints.iter().any(|g| g.nvars() != n) {
        return Err(GlobalError::BadInput {
            what: "constraint arity differs from the objective's".into(),
        });
    }
    let half = |deg: u32| deg.div_ceil(2);
    let d_min = constraints
        .iter()
        .map(|g| half(g.degree()))
        .chain(std::iter::once(half(p.degree())))
        .max()
        .unwrap_or(1)
        .max(1);
    let d = opts.order.unwrap_or(d_min).max(d_min);
    let mut prog = SosProgram::new(n);
    let gamma = prog.scalar("gamma");
    // Unconstrained: prune σ₀'s basis by the Newton polytope of p − γ.
    let basis0 = if constraints.is_empty() {
        let mut supp: Vec<Monomial> = p.terms().map(|(m, _)| m.clone()).collect();
        supp.push(Monomial::one(n));
        pruned_basis(n, d, &supp)
    } else {
        monomials_in_degree_range(n, 0, d)
    };
    let sigma0 = prog.sos("sigma0", basis0).map_err(GlobalError::Program)?;
    let one = MPoly::constant(n, 1.0);
    let mut ident = Identity::new(p)
        .term(one.scale(-1.0), gamma)
        .term(one.scale(-1.0), sigma0);
    for (k, g) in constraints.iter().enumerate() {
        let dk = d.saturating_sub(half(g.degree()));
        let sk = prog
            .sos(
                &format!("sigma{}", k + 1),
                monomials_in_degree_range(n, 0, dk),
            )
            .map_err(GlobalError::Program)?;
        ident = ident.term(g.scale(-1.0), sk);
    }
    prog.add_identity(ident).map_err(GlobalError::Program)?;
    prog.maximize(gamma, 1.0).map_err(GlobalError::Program)?;
    let sol = prog.solve(&opts.sdp).map_err(GlobalError::Program)?;
    if !matches!(sol.status(), SdpStatus::Optimal | SdpStatus::NearOptimal) {
        return Err(GlobalError::RelaxationFailed {
            status: sol.status(),
        });
    }
    let gstar = sol.scalar(gamma).unwrap_or(f64::NAN);
    if !gstar.is_finite() {
        return Err(GlobalError::RelaxationFailed {
            status: sol.status(),
        });
    }
    let (lower, certificate, backoff) = certify_backoff(&prog, gamma, gstar, opts)
        .map_err(|last| GlobalError::CertificationFailed { last })?;

    // Upper bound: candidates from the moments, polished, feasibility-checked.
    let moment = |m: &Monomial| sol.dual(0, m);
    let y0 = moment(&Monomial::one(n)).unwrap_or(1.0);
    let mut candidates: Vec<Vec<f64>> = Vec::new();
    if y0.abs() > 1e-12 {
        let mean: Option<Vec<f64>> = (0..n)
            .map(|k| moment(&Monomial::var(n, k)).map(|v| v / y0))
            .collect();
        if let Some(c) = mean {
            candidates.push(c);
        }
    }
    // Eigenvector candidates of the moment matrix (the σ₀ dual block).
    if let (Some(DecisionValue::Sos { basis, .. }), Some(zb)) =
        (sol.values.get(sigma0.index()), sol.sdp.z.first())
    {
        let nb = basis.len();
        let pos_one = basis.iter().position(|m| m.degree() == 0);
        let pos_vars: Vec<Option<usize>> = (0..n)
            .map(|k| basis.iter().position(|m| *m == Monomial::var(n, k)))
            .collect();
        if pos_one.is_some() {
            candidates.extend(extract_atoms(basis, zb, n));
        }
        if let Some(i0) = pos_one {
            let (vals, vecs) = dense::sym_eigen(zb, nb, true);
            let top = vals.last().copied().unwrap_or(0.0);
            for col in (0..nb).rev() {
                if vals[col] < 1e-6 * top {
                    break;
                }
                let v0 = vecs[i0 * nb + col];
                if v0.abs() < 1e-9 {
                    continue;
                }
                let cand: Option<Vec<f64>> = pos_vars
                    .iter()
                    .map(|pv| pv.map(|i| vecs[i * nb + col] / v0))
                    .collect();
                if let Some(c) = cand {
                    candidates.push(c);
                }
            }
        }
    }
    let mut best: Option<(f64, Vec<f64>)> = None;
    for c in candidates {
        let x = if constraints.is_empty() {
            newton_polish(p, c)
        } else {
            c
        };
        if x.iter().any(|v| !v.is_finite()) {
            continue;
        }
        let pt: Vec<Interval> = x.iter().map(|&v| Interval::point(v)).collect();
        if constraints
            .iter()
            .any(|g| !(g.eval_interval(&pt).lo() >= 0.0))
        {
            continue;
        }
        let ub = p.eval_interval(&pt).hi();
        if best.as_ref().is_none_or(|(b, _)| ub < *b) {
            best = Some((ub, x));
        }
    }
    let (upper, minimizer) = match best {
        Some((u, x)) => (Some(u), Some(x)),
        None => (None, None),
    };
    Ok(GlobalBound {
        lower,
        upper,
        minimizer,
        relaxation_value: gstar,
        order: d,
        certificate,
        backoff,
    })
}

/// Back off from `γ*` until a centred re-solve verifies.
fn certify_backoff(
    prog: &SosProgram,
    gamma: DecisionId,
    gstar: f64,
    opts: &GlobalOptions,
) -> Result<(f64, Certificate, f64), VerifyError> {
    let mut rel = opts.initial_backoff;
    let mut last = VerifyError::ShapeMismatch {
        what: "no certification attempt was made".into(),
    };
    for _ in 0..opts.certify_attempts.max(1) {
        let backoff = rel * (1.0 + gstar.abs());
        let gc = gstar - backoff;
        match prog.solve_centered(&[(gamma, gc)], &opts.sdp) {
            Ok(cs) => match verify(prog, &cs.values) {
                Ok(cert) => return Ok((gc, cert, backoff)),
                Err(e) => last = e,
            },
            Err(e) => {
                last = VerifyError::ShapeMismatch {
                    what: format!("centred re-solve refused: {e}"),
                };
            }
        }
        rel *= 10.0;
    }
    Err(last)
}

/// Henrion–Lasserre extraction of the atoms of a (numerically) flat moment
/// matrix `M` over `basis`: factor `M = VVᵀ` at its numerical rank `r`,
/// reduce `V` to column echelon form `U` (so `b(x) = U·w(x)` on the atoms for
/// the pivot monomials `w`), read the multiplication matrices `N_i` (rows of
/// `U` at the monomials `x_i·w_j`), and diagonalize a fixed generic
/// combination `Σ c_i N_i`; the atoms' coordinates are the common
/// eigenvalues. Returns no points when the needed shifted monomials are
/// missing from the basis (raise the order) — never a guess.
#[allow(clippy::too_many_lines)] // the four extraction stages read best in sequence
fn extract_atoms(basis: &[Monomial], m: &[f64], n: usize) -> Vec<Vec<f64>> {
    let nb = basis.len();
    let (vals, vecs) = dense::sym_eigen(m, nb, true);
    let top = vals.last().copied().unwrap_or(0.0);
    if !(top > 0.0) {
        return Vec::new();
    }
    let cols: Vec<usize> = (0..nb).filter(|&c| vals[c] > 1e-6 * top).collect();
    let r = cols.len();
    if r == 0 || r > nb {
        return Vec::new();
    }
    // V (nb × r), row-major.
    let mut v = vec![0.0; nb * r];
    for (k, &c) in cols.iter().enumerate() {
        let s = vals[c].sqrt();
        for i in 0..nb {
            v[i * r + k] = vecs[i * nb + c] * s;
        }
    }
    // Column echelon form by column operations, scanning rows in graded order.
    let mut pivots: Vec<usize> = Vec::with_capacity(r); // pivot row per column
    let mut used = vec![false; r];
    let scale = v.iter().fold(0.0f64, |s, x| s.max(x.abs()));
    for row in 0..nb {
        if pivots.len() == r {
            break;
        }
        let mut best: Option<usize> = None;
        for c in 0..r {
            if !used[c] && best.is_none_or(|b| v[row * r + c].abs() > v[row * r + b].abs()) {
                best = Some(c);
            }
        }
        let Some(c) = best else { break };
        let piv = v[row * r + c];
        if piv.abs() <= 1e-8 * scale {
            continue;
        }
        for i in 0..nb {
            v[i * r + c] /= piv;
        }
        for c2 in 0..r {
            if c2 != c {
                let f = v[row * r + c2];
                if f != 0.0 {
                    for i in 0..nb {
                        v[i * r + c2] -= f * v[i * r + c];
                    }
                }
            }
        }
        used[c] = true;
        pivots.push(row);
        // Keep column order aligned with pivot order.
        let pos = pivots.len() - 1;
        if c != pos {
            for i in 0..nb {
                v.swap(i * r + c, i * r + pos);
            }
            used.swap(c, pos);
        }
    }
    if pivots.len() != r {
        return Vec::new();
    }
    let w: Vec<&Monomial> = pivots.iter().map(|&p| &basis[p]).collect();
    // Multiplication matrices N_i (r × r).
    let mut nmats: Vec<Vec<f64>> = Vec::with_capacity(n);
    for i in 0..n {
        let xi = Monomial::var(n, i);
        let mut ni = vec![0.0; r * r];
        for (j, wj) in w.iter().enumerate() {
            let target = xi.mul(wj);
            let Some(row) = basis.iter().position(|b| *b == target) else {
                return Vec::new();
            };
            for k in 0..r {
                ni[j * r + k] = v[row * r + k];
            }
        }
        nmats.push(ni);
    }
    // Generic combination (fixed, deterministic, irrational-looking weights).
    let weights: Vec<f64> = (0..n)
        .map(|i| 1.0 / (1.0 + std::f64::consts::E * (i as f64 + 0.618)))
        .collect();
    let mut comb = vec![0.0; r * r];
    for (i, ni) in nmats.iter().enumerate() {
        for k in 0..r * r {
            comb[k] += weights[i] * ni[k];
        }
    }
    let eig = real_eigenvalues(&comb, r);
    let mut atoms = Vec::new();
    for lam in eig {
        // Right eigenvector by inverse iteration; atom coordinates are the
        // Rayleigh quotients of the N_i on it.
        let mut shifted = comb.clone();
        let mag = lam.abs().max(1.0);
        for k in 0..r {
            shifted[k * r + k] -= lam + 1e-10 * mag;
        }
        let mut q = vec![1.0; r];
        for _ in 0..3 {
            let Some(sol) = dense::lu_solve(&shifted, r, &q) else {
                break;
            };
            let nrm = sol.iter().map(|x| x * x).sum::<f64>().sqrt();
            if !(nrm > 0.0) || !nrm.is_finite() {
                break;
            }
            q = sol.iter().map(|x| x / nrm).collect();
        }
        let atom: Vec<f64> = nmats
            .iter()
            .map(|ni| {
                let nq: Vec<f64> = (0..r)
                    .map(|a| (0..r).map(|b| ni[a * r + b] * q[b]).sum())
                    .collect();
                nq.iter().zip(&q).map(|(a, b)| a * b).sum::<f64>()
            })
            .collect();
        if atom.iter().all(|x| x.is_finite()) {
            atoms.push(atom);
        }
    }
    atoms
}

/// Real eigenvalues of a small nonsymmetric matrix by Wilkinson-shifted QR
/// iteration with deflation (Gram–Schmidt QR; adequate for the r ≤ ~20
/// matrices of moment extraction). Complex pairs, if any, are skipped — real
/// atoms always give real eigenvalues.
fn real_eigenvalues(a: &[f64], n: usize) -> Vec<f64> {
    let mut h = a.to_vec();
    let mut out = Vec::with_capacity(n);
    let mut size = n;
    let mut iter = 0usize;
    while size > 0 && iter < 500 * n.max(1) {
        iter += 1;
        if size == 1 {
            out.push(h[0]);
            break;
        }
        let last = size - 1;
        let (a11, a12) = (h[(last - 1) * n + (last - 1)], h[(last - 1) * n + last]);
        let (a21, a22) = (h[last * n + (last - 1)], h[last * n + last]);
        let half_tr = 0.5 * (a11 + a22);
        let disc = 0.25 * (a11 - a22) * (a11 - a22) + a12 * a21;
        if size == 2 {
            if disc >= 0.0 {
                let r = disc.sqrt();
                out.push(half_tr + r);
                out.push(half_tr - r);
            }
            break;
        }
        // Deflate when the last row's sub-diagonal part is negligible.
        let off: f64 = (0..last).map(|k| h[last * n + k].abs()).sum();
        let scale = a22.abs() + a11.abs() + 1e-300;
        if off <= 1e-13 * scale {
            out.push(a22);
            size -= 1;
            continue;
        }
        // Wilkinson shift: the trailing 2×2 eigenvalue nearer a22.
        let mut mu = if disc >= 0.0 {
            let r = disc.sqrt();
            let (e1, e2) = (half_tr + r, half_tr - r);
            if (e1 - a22).abs() <= (e2 - a22).abs() {
                e1
            } else {
                e2
            }
        } else {
            a22
        };
        if iter.is_multiple_of(11) {
            mu += 0.75 * off; // exceptional shift breaks symmetric stagnation
        }
        // QR of (H − μI) restricted to the leading `size` block.
        let mut q = vec![0.0; size * size];
        let mut rmat = vec![0.0; size * size];
        for j in 0..size {
            let mut col: Vec<f64> = (0..size)
                .map(|i| h[i * n + j] - if i == j { mu } else { 0.0 })
                .collect();
            for k in 0..j {
                let d: f64 = (0..size).map(|i| q[i * size + k] * col[i]).sum();
                rmat[k * size + j] = d;
                for i in 0..size {
                    col[i] -= d * q[i * size + k];
                }
            }
            let nrm = col.iter().map(|x| x * x).sum::<f64>().sqrt();
            rmat[j * size + j] = nrm;
            for i in 0..size {
                q[i * size + j] = if nrm > 0.0 { col[i] / nrm } else { 0.0 };
            }
        }
        // H ← R Q + μI.
        for i in 0..size {
            for j in 0..size {
                let s: f64 = (0..size)
                    .map(|k| rmat[i * size + k] * q[k * size + j])
                    .sum();
                h[i * n + j] = s + if i == j { mu } else { 0.0 };
            }
        }
    }
    out
}

/// Damped Newton on `∇p = 0`, accepting only steps that decrease `p`.
fn newton_polish(p: &MPoly, mut x: Vec<f64>) -> Vec<f64> {
    let n = p.nvars();
    let grad: Vec<MPoly> = (0..n).map(|k| p.derivative(k)).collect();
    let hess: Vec<Vec<MPoly>> = grad
        .iter()
        .map(|g| (0..n).map(|k| g.derivative(k)).collect())
        .collect();
    let mut fx = p.eval(&x);
    for _ in 0..50 {
        let g: Vec<f64> = grad.iter().map(|gk| gk.eval(&x)).collect();
        let gn = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gn < 1e-14 * (1.0 + fx.abs()) {
            break;
        }
        let h: Vec<f64> = hess
            .iter()
            .flat_map(|row| row.iter().map(|hk| hk.eval(&x)))
            .collect();
        let neg: Vec<f64> = g.iter().map(|v| -v).collect();
        // Newton first; steepest descent when the Newton step does not
        // decrease p (indefinite Hessian near a saddle).
        let mut directions = Vec::with_capacity(2);
        if let Some(s) = dense::lu_solve(&h, n, &neg) {
            directions.push(s);
        }
        directions.push(neg);
        let mut improved = false;
        'dirs: for step in &directions {
            let mut t = 1.0;
            for _ in 0..40 {
                let xn: Vec<f64> = x.iter().zip(step).map(|(a, s)| a + t * s).collect();
                let fnew = p.eval(&xn);
                if fnew < fx {
                    x = xn;
                    fx = fnew;
                    improved = true;
                    break 'dirs;
                }
                t *= 0.5;
            }
        }
        if !improved {
            break;
        }
    }
    x
}
