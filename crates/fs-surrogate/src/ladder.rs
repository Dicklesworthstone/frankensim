//! THE CERTIFIED ABSTRACTION LADDER (addendum Proposal A, bead knh1.4;
//! [F] — behind the `abstraction-ladder` feature): reduced-order VIEWS
//! at multiple abstraction levels whose fidelity to the level below is
//! itself a CERTIFIED quantity, with a LEAK ALARM — the operator
//! reasons at the concept level until a certificate says the
//! abstraction is no longer faithful HERE, and drills down only there.
//! The ladder is invisible until it leaks, which is the point.
//!
//! THE BEACHHEAD (mature technology, shipped first): certified
//! REDUCED-BASIS for the affine-parametric elliptic family
//! `−(a(x;μ) u′)′ = f` with `a = 1 + μ·χ` — offline snapshots +
//! energy-orthonormal basis + online k×k Galerkin, and the TEXTBOOK
//! a-posteriori bound `‖u−u_rb‖_a ≤ ‖r‖_{V′}/√α_LB` with the residual
//! dual norm assembled offline (Riesz representers, exact quadratic
//! form) and the coercivity floor `α_LB = min a(x;μ)` exact for the
//! affine family. The compliance QoI inherits the CLASSIC squared
//! bound `|s − s_rb| = ‖u−u_rb‖²_μ ≤ Δ²` (Galerkin symmetry).
//!
//! THE HONEST FRONTIER: non-RB concept levels carry ESTIMATED color,
//! calibrated by cross-rung discrepancy probes — the type system keeps
//! anyone from mistaking them for the RB-certified rungs.

use fs_evidence::Color;

/// The full-order "truth" model: P1 finite elements on a uniform grid
/// for `−(a u′)′ = 1`, `u(0) = u(1) = 0`, `a(x;μ) = 1 + μ·χ_{[½,1]}`.
/// (The FE model IS level 0's semantics — discretization honesty is
/// documented in the CONTRACT, not hidden in the bound.)
#[derive(Debug, Clone)]
pub struct TruthModel {
    /// Interior nodes.
    pub n: usize,
}

impl TruthModel {
    /// Assemble and solve at parameter `mu` (Thomas algorithm).
    /// Returns interior nodal values.
    #[must_use]
    pub fn solve(&self, mu: f64) -> Vec<f64> {
        let n = self.n;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        // Tridiagonal: diag[i], off[i] between i and i+1.
        let coeff = |i: usize| -> f64 {
            // Element i spans [i·h, (i+1)·h]; a is 1 or 1+mu.
            #[allow(clippy::cast_precision_loss)]
            let mid = (i as f64 + 0.5) * h;
            if mid >= 0.5 { 1.0 + mu } else { 1.0 }
        };
        let mut diag = vec![0.0f64; n];
        let mut off = vec![0.0f64; n.saturating_sub(1)];
        let rhs = vec![h; n];
        for i in 0..=n {
            let a = coeff(i) / h;
            if i < n {
                diag[i] += a;
            }
            if i > 0 {
                diag[i - 1] += a;
            }
            if i > 0 && i < n {
                off[i - 1] -= a;
            }
        }
        // Thomas solve.
        let mut c = off.clone();
        let mut d = rhs.clone();
        c[0] /= diag[0];
        d[0] /= diag[0];
        for i in 1..n {
            let m = diag[i] - off[i - 1] * c[i - 1];
            if i < n - 1 {
                c[i] = off[i] / m;
            }
            d[i] = (d[i] - off[i - 1] * d[i - 1]) / m;
        }
        let mut u = d;
        for i in (0..n - 1).rev() {
            u[i] -= c[i] * u[i + 1];
        }
        let _ = rhs;
        u
    }

    /// The energy inner product `a(u, v; mu)`.
    #[must_use]
    pub fn energy(&self, u: &[f64], v: &[f64], mu: f64) -> f64 {
        let n = self.n;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let mut acc = 0.0f64;
        for i in 0..=n {
            #[allow(clippy::cast_precision_loss)]
            let mid = (i as f64 + 0.5) * h;
            let a = if mid >= 0.5 { 1.0 + mu } else { 1.0 };
            let du = if i == 0 {
                u[0]
            } else if i == n {
                -u[n - 1]
            } else {
                u[i] - u[i - 1]
            };
            let dv = if i == 0 {
                v[0]
            } else if i == n {
                -v[n - 1]
            } else {
                v[i] - v[i - 1]
            };
            acc += a * du * dv / h;
        }
        acc
    }

    /// The compliance QoI `∫ f u = h·Σu` (f = 1).
    #[must_use]
    pub fn compliance(&self, u: &[f64]) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (self.n as f64 + 1.0);
        h * u.iter().sum::<f64>()
    }
}

/// One certified RB rung: an energy-orthonormal basis with the offline
/// residual machinery for the exact online dual-norm bound.
#[derive(Debug, Clone)]
pub struct RbLevel {
    truth: TruthModel,
    /// Basis vectors (each length n).
    basis: Vec<Vec<f64>>,
    /// μ-range the level was trained on (its declared query class).
    pub mu_range: (f64, f64),
}

impl RbLevel {
    /// OFFLINE: snapshots at `k` training parameters spread over
    /// `mu_range`, Gram–Schmidt in the reference (μ = mid) energy.
    #[must_use]
    pub fn train(truth: &TruthModel, mu_range: (f64, f64), k: usize) -> RbLevel {
        let mid = f64::midpoint(mu_range.0, mu_range.1);
        let mut basis: Vec<Vec<f64>> = Vec::with_capacity(k);
        for j in 0..k {
            #[allow(clippy::cast_precision_loss)]
            let t = if k == 1 {
                0.5
            } else {
                j as f64 / (k - 1) as f64
            };
            let mu = mu_range.0 + t * (mu_range.1 - mu_range.0);
            let mut snap = truth.solve(mu);
            // Gram–Schmidt against accepted basis (reference energy).
            for b in &basis {
                let proj = truth.energy(&snap, b, mid);
                for (s, bi) in snap.iter_mut().zip(b) {
                    *s -= proj * bi;
                }
            }
            let norm = truth.energy(&snap, &snap, mid).sqrt();
            if norm > 1e-10 {
                for s in &mut snap {
                    *s /= norm;
                }
                basis.push(snap);
            }
        }
        RbLevel {
            truth: truth.clone(),
            basis,
            mu_range,
        }
    }

    /// Basis size.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.basis.len()
    }

    /// ONLINE: Galerkin solve in the basis + the CERTIFIED energy and
    /// QoI bounds. Returns (`u_rb` in nodal form, compliance,
    /// `energy_bound`, `qoi_bound`).
    #[must_use]
    pub fn query(&self, mu: f64) -> (Vec<f64>, f64, f64, f64) {
        let k = self.basis.len();
        // RB system: A_ij = a(ξ_j, ξ_i; μ), b_i = f(ξ_i).
        let mut a = vec![vec![0.0f64; k]; k];
        let mut b = vec![0.0f64; k];
        for i in 0..k {
            for (j, aij) in a[i].iter_mut().enumerate() {
                *aij = self.truth.energy(&self.basis[j], &self.basis[i], mu);
            }
            b[i] = self.truth.compliance(&self.basis[i]);
        }
        // Dense Gaussian elimination (k is small).
        let coef = solve_dense(&mut a, &mut b);
        let n = self.truth.n;
        let mut u = vec![0.0f64; n];
        for (c, basis_vec) in coef.iter().zip(&self.basis) {
            for (ui, bi) in u.iter_mut().zip(basis_vec) {
                *ui += c * bi;
            }
        }
        let s_rb = self.truth.compliance(&u);
        // CERTIFIED bound: the residual's Riesz representer r_h solves
        // (∇r_h, ∇v) = f(v) − a(u_rb, v; μ) for all v (reference
        // Laplacian, μ = 0 metric); then ‖r‖_{V′} = ‖r_h‖_{H¹₀} and
        // ‖u − u_rb‖_a ≤ ‖r‖_{V′}/√α_LB with α_LB = min(1, 1+μ) ≥ 1
        // for μ ≥ 0 (exact for the affine family).
        let riesz = self.residual_riesz(&u, mu);
        let dual_norm = self.truth.energy(&riesz, &riesz, 0.0).sqrt();
        let alpha_lb = 1.0f64.min(1.0 + mu);
        let energy_bound = dual_norm / alpha_lb.sqrt();
        // Compliance (symmetric, f = load = QoI): s − s_rb = ‖e‖²_μ,
        // so the certified QoI bound is energy_bound².
        let qoi_bound = energy_bound * energy_bound;
        (u, s_rb, energy_bound, qoi_bound)
    }

    /// Solve the reference-Laplacian Riesz problem for the residual.
    fn residual_riesz(&self, u_rb: &[f64], mu: f64) -> Vec<f64> {
        let n = self.truth.n;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        // rhs_i = f(φ_i) − a(u_rb, φ_i; μ).
        let mut rhs = vec![h; n];
        for (i, ri) in rhs.iter_mut().enumerate() {
            // a(u_rb, φ_i; μ) over the two elements adjacent to node i.
            let mut acc = 0.0f64;
            for e in [i, i + 1] {
                #[allow(clippy::cast_precision_loss)]
                let mid = (e as f64 + 0.5) * h;
                let a = if mid >= 0.5 { 1.0 + mu } else { 1.0 };
                let du = if e == 0 {
                    u_rb[0]
                } else if e == n {
                    -u_rb[n - 1]
                } else {
                    u_rb[e] - u_rb[e - 1]
                };
                // dφ_i on element e: +1/h on [x_{i-1},x_i] side... φ_i
                // slope is +1/h on element i, −1/h on element i+1.
                let dphi = if e == i { 1.0 } else { -1.0 };
                acc += a * du * dphi / h;
            }
            *ri -= acc;
        }
        // Reference Laplacian (a ≡ 1) Thomas solve.
        let mut diag = vec![2.0 / h; n];
        let off = vec![-1.0 / h; n.saturating_sub(1)];
        let mut c = off.clone();
        let mut d = rhs;
        c[0] /= diag[0];
        d[0] /= diag[0];
        for i in 1..n {
            let m = diag[i] - off[i - 1] * c[i - 1];
            if i < n - 1 {
                c[i] = off[i] / m;
            }
            d[i] = (d[i] - off[i - 1] * d[i - 1]) / m;
        }
        for i in (0..n - 1).rev() {
            let tmp = c[i] * d[i + 1];
            d[i] -= tmp;
        }
        diag.clear();
        d
    }
}

fn solve_dense(a: &mut [Vec<f64>], b: &mut [f64]) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot (deterministic max-abs, lowest index).
        let mut piv = col;
        for r in col + 1..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let p = a[col][col];
        let pivot_row = a[col].clone();
        for r in col + 1..n {
            let f = a[r][col] / p;
            for (arc, pc) in a[r][col..].iter_mut().zip(&pivot_row[col..]) {
                *arc -= f * pc;
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0f64; n];
    for r in (0..n).rev() {
        let mut acc = b[r];
        for c in r + 1..n {
            acc -= a[r][c] * x[c];
        }
        x[r] = acc / a[r][r];
    }
    x
}

/// A NON-RB concept level: QoI lookup by linear interpolation over a
/// training grid — cheap, USEFUL, and honestly ESTIMATED. Its
/// dispersion is calibrated by CROSS-RUNG DISCREPANCY PROBES against
/// the rung below (Proposal 3's pattern pointed vertically).
#[derive(Debug, Clone)]
pub struct ConceptLevel {
    grid: Vec<(f64, f64)>,
    /// Calibrated dispersion (max |concept − rung below| over probes).
    pub dispersion: f64,
}

impl ConceptLevel {
    /// Build from a training grid and calibrate against a lower rung.
    #[must_use]
    pub fn train(rb: &RbLevel, grid_points: usize, probes: usize) -> ConceptLevel {
        let (lo, hi) = rb.mu_range;
        let grid: Vec<(f64, f64)> = (0..grid_points)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let mu = lo + (hi - lo) * i as f64 / (grid_points - 1) as f64;
                (mu, rb.query(mu).1)
            })
            .collect();
        let mut level = ConceptLevel {
            grid,
            dispersion: 0.0,
        };
        // Discrepancy probes at off-grid parameters.
        let mut disp = 0.0f64;
        for i in 0..probes {
            #[allow(clippy::cast_precision_loss)]
            let mu = lo + (hi - lo) * (i as f64 + 0.5) / probes as f64;
            disp = disp.max((level.lookup(mu) - rb.query(mu).1).abs());
        }
        level.dispersion = disp;
        level
    }

    /// Linear interpolation lookup.
    #[must_use]
    pub fn lookup(&self, mu: f64) -> f64 {
        let g = &self.grid;
        let pos = g.partition_point(|(m, _)| *m < mu).clamp(1, g.len() - 1);
        let (m0, v0) = g[pos - 1];
        let (m1, v1) = g[pos];
        let t = ((mu - m0) / (m1 - m0)).clamp(0.0, 1.0);
        v0 + t * (v1 - v0)
    }
}

/// A ladder answer: the value, its certified-or-estimated color, and
/// the drill-down forensics.
#[derive(Debug, Clone)]
pub struct LadderAnswer {
    /// The QoI value.
    pub value: f64,
    /// The evidence color: `Verified` from an RB rung (the bound is
    /// the interval half-width), `Estimated` from a concept rung.
    pub color: Color,
    /// The level that ultimately answered (0 = full order).
    pub level_used: usize,
    /// Levels that LEAKED (bound exceeded tolerance) on the way down.
    pub leaks: Vec<usize>,
}

/// The assembled ladder: level 0 = full order (truth), higher = more
/// abstract. `at_level(k)` starts at rung k and DESCENDS AUTOMATICALLY
/// on leak — invisible until it leaks.
#[derive(Debug, Clone)]
pub struct Ladder {
    truth: TruthModel,
    /// RB rungs, coarser (smaller basis) at higher indices.
    pub rb_levels: Vec<RbLevel>,
    /// Optional concept rung above the RB rungs.
    pub concept: Option<ConceptLevel>,
}

/// A view of the ladder starting at a rung.
pub struct LevelView<'a> {
    ladder: &'a Ladder,
    start: usize,
}

impl Ladder {
    /// Build: truth + RB rungs of the given basis sizes (descending
    /// fidelity) + a concept rung calibrated against the last RB rung.
    #[must_use]
    pub fn build(n: usize, mu_range: (f64, f64), rb_dims: &[usize], concept: bool) -> Ladder {
        let truth = TruthModel { n };
        let rb_levels: Vec<RbLevel> = rb_dims
            .iter()
            .map(|&k| RbLevel::train(&truth, mu_range, k))
            .collect();
        let concept = concept
            .then(|| ConceptLevel::train(rb_levels.last().expect("at least one RB rung"), 5, 7));
        Ladder {
            truth,
            rb_levels,
            concept,
        }
    }

    /// The view starting at rung `k` (0 = full order, 1.. = RB rungs,
    /// rb_levels.len()+1 = concept rung if present).
    #[must_use]
    pub fn at_level(&self, k: usize) -> LevelView<'_> {
        LevelView {
            ladder: self,
            start: k,
        }
    }

    /// Highest rung index.
    #[must_use]
    pub fn top(&self) -> usize {
        self.rb_levels.len() + usize::from(self.concept.is_some())
    }
}

impl LevelView<'_> {
    /// Query with AUTOMATIC CERTIFIED DESCENT: each rung answers only
    /// if its certificate meets `tol`; a leaking rung is recorded and
    /// the query drills down. Level 0 answers unconditionally (it IS
    /// the declared truth semantics).
    #[must_use]
    pub fn query(&self, mu: f64, tol: f64) -> LadderAnswer {
        let mut leaks = Vec::new();
        let mut level = self.start.min(self.ladder.top());
        loop {
            // Concept rung: estimated color — it can NEVER satisfy a
            // certified tolerance; it answers only if the caller's tol
            // admits estimates (tol >= dispersion is still an ESTIMATE,
            // never verified — the color says so).
            if level > self.ladder.rb_levels.len() {
                let c = self.ladder.concept.as_ref().expect("concept rung");
                if c.dispersion <= tol {
                    return LadderAnswer {
                        value: c.lookup(mu),
                        color: Color::Estimated {
                            estimator: "concept-lookup(cross-rung-probes)".to_string(),
                            dispersion: c.dispersion,
                        },
                        level_used: level,
                        leaks,
                    };
                }
                leaks.push(level);
                level -= 1;
                continue;
            }
            if level >= 1 {
                let rb = &self.ladder.rb_levels[level - 1];
                let (_, s_rb, _, qoi_bound) = rb.query(mu);
                if qoi_bound <= tol {
                    return LadderAnswer {
                        value: s_rb,
                        color: Color::Verified {
                            lo: s_rb - qoi_bound,
                            hi: s_rb + qoi_bound,
                        },
                        level_used: level,
                        leaks,
                    };
                }
                leaks.push(level);
                level -= 1;
                continue;
            }
            // Level 0: the full-order truth.
            let u = self.ladder.truth.solve(mu);
            let s = self.ladder.truth.compliance(&u);
            return LadderAnswer {
                value: s,
                color: Color::Verified { lo: s, hi: s },
                level_used: 0,
                leaks,
            };
        }
    }
}

/// THE KILL MEASUREMENT (Proposal A): the fraction of a query battery
/// answerable at an RB rung (bound ≤ tol) WITHOUT drilling to full
/// order. Below 0.2 the beachhead is too narrow — park certification.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn rb_coverage(ladder: &Ladder, mus: &[f64], tols: &[f64]) -> f64 {
    let mut covered = 0usize;
    let mut total = 0usize;
    for &mu in mus {
        for &tol in tols {
            total += 1;
            let ans = ladder.at_level(ladder.rb_levels.len()).query(mu, tol);
            if ans.level_used >= 1 {
                covered += 1;
            }
        }
    }
    covered as f64 / total.max(1) as f64
}
