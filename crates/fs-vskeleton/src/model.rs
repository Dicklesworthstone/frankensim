//! The skeleton's physics + optimization: a plate with a parametric circular
//! hole under a variable-coefficient diffusion equation.
//!
//! Problem: on the unit square with u = 0 on the boundary,
//!   -div( rho(x; r) grad u ) = 1,
//! where rho is an ersatz-material density from a SMOOTH Heaviside of the
//! hole SDF (hard Booleans have derivative discontinuities that poison shape
//! optimization — plan §7.2's R-function doctrine, here in 1-parameter form).
//!
//! Objective: J(r) = compliance + w_vol * material volume
//!   compliance = h² Σ u_i   (f ≡ 1),   volume = h² Σ rho_i.
//!
//! Gradient truth (plan §8.7 in miniature): K(r) u = f is self-adjoint, so
//! the adjoint of compliance is λ = h² u, and
//!   dJ/dr = -h² Σ_edges (∂rho_e/∂r) (u_i - u_j)² / h²  +  w_vol h² Σ ∂rho_i/∂r.
//! Both the primal apply and that sensitivity contraction are derived from
//! ONE [`EdgeLaw`] — the fs-opdsl one-source-of-truth pattern in embryo.

use crate::sexpr::{self, Sexpr};
use fs_qty::parse::parse_qty;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// Parsed study specification (the Five Explicits in miniature: seed and
/// budgets are mandatory; units flow through fs-qty; versions/capabilities
/// arrive with the real fs-ir).
#[derive(Debug, Clone)]
pub struct StudySpec {
    /// Study name.
    pub name: String,
    /// Seed (recorded provenance; this slice is deterministic without RNG).
    pub seed: u64,
    /// Nodes per side (grid x grid), boundary included.
    pub grid: usize,
    /// Initial hole radius (unit-square coordinates).
    pub initial_radius: f64,
    /// Optimizer step count.
    pub opt_steps: usize,
    /// Gradient-descent step size.
    pub step_size: f64,
    /// Volume-term weight in the objective.
    pub volume_weight: f64,
    /// Total CG-iteration budget (P4: budgets are mandatory and enforced).
    pub cg_budget: u64,
    /// Radius projection bounds.
    pub r_min: f64,
    /// Upper bound.
    pub r_max: f64,
}

impl StudySpec {
    /// Seed rendered as the ledger's hex form.
    #[must_use]
    pub fn seed_hex(&self) -> String {
        format!("{:#018x}", self.seed)
    }

    /// Parse a `(study "name" (seed ..) (grid ..) (budget (cg-iters ..)) ...)`
    /// form with teaching errors.
    ///
    /// # Errors
    /// Returns a message naming the missing/malformed field and its fix.
    pub fn parse(text: &str) -> Result<StudySpec, String> {
        let expr = sexpr::parse(text).map_err(|(at, m)| {
            format!("study parse error at byte {at}: {m}; studies look like (study \"name\" (seed 0x..) ...)")
        })?;
        let body = expr
            .as_form("study")
            .ok_or("top-level form must be (study \"name\" ...)")?;
        let name = body
            .first()
            .and_then(Sexpr::atom)
            .ok_or("study needs a name atom right after the head")?
            .to_string();

        let seed_txt = field(body, "seed")?;
        let seed = parse_seed(seed_txt)?;
        let grid: usize = num_field(body, "grid")? as usize;
        if !(5..=513).contains(&grid) {
            return Err(format!("(grid {grid}) out of the supported 5..=513 range"));
        }
        let budget_body = Sexpr::find_form(body, "budget")
            .ok_or("missing (budget (cg-iters N)) — budgets are mandatory (P4)")?;
        let cg_budget = Sexpr::find_form(budget_body, "cg-iters")
            .and_then(|b| b.first())
            .and_then(Sexpr::atom)
            .ok_or("budget needs (cg-iters N)")?
            .parse::<u64>()
            .map_err(|e| format!("cg-iters must be a positive integer: {e}"))?;

        let initial_radius = dimensionless_field(body, "hole-radius")?;
        if !(0.01..=0.45).contains(&initial_radius) {
            return Err(format!(
                "(hole-radius {initial_radius}) outside (0.01, 0.45) — the hole must fit \
                 strictly inside the unit square"
            ));
        }
        let opt_steps = num_field(body, "opt-steps")? as usize;
        let step_size = dimensionless_field(body, "step-size")?;
        let volume_weight = dimensionless_field(body, "volume-weight")?;

        Ok(StudySpec {
            name,
            seed,
            grid,
            initial_radius,
            opt_steps,
            step_size,
            volume_weight,
            cg_budget,
            r_min: 0.02,
            r_max: 0.45,
        })
    }
}

fn field<'a>(body: &'a [Sexpr], name: &str) -> Result<&'a str, String> {
    Sexpr::find_form(body, name)
        .and_then(|b| b.first())
        .and_then(Sexpr::atom)
        .ok_or_else(|| format!("missing ({name} <value>) in study"))
}

fn num_field(body: &[Sexpr], name: &str) -> Result<u64, String> {
    field(body, name)?
        .parse::<u64>()
        .map_err(|e| format!("({name} ...) must be a positive integer: {e}"))
}

/// Numeric study fields route through fs-qty so DIMENSIONED literals are
/// caught with a teaching error (the unit discipline exercised end-to-end).
fn dimensionless_field(body: &[Sexpr], name: &str) -> Result<f64, String> {
    let text = field(body, name)?;
    let q = parse_qty(text).map_err(|e| format!("({name} {text}): {e}"))?;
    if !q.dims.is_none() {
        return Err(format!(
            "({name} {text}) carries dimensions [{}] but this slice works in unit-square \
             coordinates — pass a dimensionless number",
            q.dims.unit_string()
        ));
    }
    Ok(q.value)
}

fn parse_seed(text: &str) -> Result<u64, String> {
    let t = text.trim();
    if let Some(hex) = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X")) {
        u64::from_str_radix(hex, 16).map_err(|e| format!("(seed {t}): bad hex: {e}"))
    } else {
        t.parse::<u64>()
            .map_err(|e| format!("(seed {t}): bad integer: {e}"))
    }
}

// ---------------------------------------------------------------------------
// The one-source-of-truth edge law (fs-opdsl seed).
// ---------------------------------------------------------------------------

/// Smoothed density and its radius sensitivity at a point: the SINGLE
/// definition both the primal operator and the adjoint contraction consume.
#[derive(Debug, Clone, Copy)]
pub struct EdgeLaw {
    /// Void (ersatz) density floor keeping K positive definite.
    pub eps: f64,
    /// Smoothing half-width of the Heaviside (in domain units).
    pub width: f64,
    /// Hole radius.
    pub radius: f64,
}

impl EdgeLaw {
    /// rho(x,y) = eps + (1-eps) * Hs(dist - r), Hs(d) = (1 + tanh(d/w)) / 2.
    #[must_use]
    pub fn rho(&self, x: f64, y: f64) -> f64 {
        let d = self.signed_gap(x, y);
        self.eps + (1.0 - self.eps) * 0.5 * (1.0 + (d / self.width).tanh())
    }

    /// ∂rho/∂r = (1-eps) * Hs'(d) * ∂d/∂r with ∂d/∂r = -1.
    #[must_use]
    pub fn d_rho_d_radius(&self, x: f64, y: f64) -> f64 {
        let d = self.signed_gap(x, y);
        let t = (d / self.width).tanh();
        -(1.0 - self.eps) * 0.5 * (1.0 - t * t) / self.width
    }

    /// Distance outside the hole boundary (positive = material side).
    fn signed_gap(&self, x: f64, y: f64) -> f64 {
        let (cx, cy) = (0.5, 0.5);
        ((x - cx).powi(2) + (y - cy).powi(2)).sqrt() - self.radius
    }
}

/// One primal/adjoint evaluation at a given radius.
#[derive(Debug, Clone)]
pub struct Evaluation {
    /// J = compliance + w_vol * volume.
    pub objective: f64,
    /// dJ/dr via the adjoint identity.
    pub gradient: f64,
    /// CG iterations used for this evaluation.
    pub cg_iters: u64,
}

/// Evaluate objective and ADJOINT gradient at `radius`.
///
/// # Errors
/// Returns an error if CG fails to converge within the per-solve cap or the
/// evaluation is cancelled.
pub fn evaluate(
    spec: &StudySpec,
    radius: f64,
    cancel: &Arc<AtomicBool>,
    cg_spent: &mut u64,
) -> Result<Evaluation, String> {
    let n = spec.grid;
    let h = 1.0 / (n - 1) as f64;
    let law = EdgeLaw {
        eps: 1e-3,
        width: 2.0 * h,
        radius,
    };

    // Nodal densities (fixed evaluation order: row-major — determinism by
    // construction; the parallel map below uses FIXED chunk boundaries).
    let rho = map_nodes_deterministic(n, cancel, |x, y| law.rho(x, y))?;
    let drho = map_nodes_deterministic(n, cancel, |x, y| law.d_rho_d_radius(x, y))?;

    // Solve K(rho) u = f, f ≡ 1 on interior nodes.
    let (u, iters) = cg_solve(n, h, &rho, spec, cancel)?;
    *cg_spent += iters;

    // Compliance + volume, fixed summation order (P2: deterministic
    // reductions; single accumulation pass in index order).
    let h2 = h * h;
    let mut compliance = 0.0;
    let mut volume = 0.0;
    for i in 0..n * n {
        compliance += u[i];
        volume += rho[i];
    }
    compliance *= h2;
    volume *= h2;

    // Adjoint gradient: λ = h² u (self-adjoint K, J_u = h²·1 = h²·f), so
    // dCompliance/dr = -h² Σ_edges (∂rho_e/∂r)(u_i-u_j)²/h². The edge loop
    // below is DERIVED from the same EdgeLaw as the operator apply.
    let mut d_compliance = 0.0;
    let idx = |ix: usize, iy: usize| iy * n + ix;
    for iy in 0..n {
        for ix in 0..n {
            // Right and top edges only: each edge counted once.
            if ix + 1 < n {
                let a = idx(ix, iy);
                let b = idx(ix + 1, iy);
                let d_edge = f64::midpoint(drho[a], drho[b]);
                let du = u[a] - u[b];
                d_compliance -= d_edge * du * du; // (h² and 1/h² cancel)
            }
            if iy + 1 < n {
                let a = idx(ix, iy);
                let b = idx(ix, iy + 1);
                let d_edge = f64::midpoint(drho[a], drho[b]);
                let du = u[a] - u[b];
                d_compliance -= d_edge * du * du;
            }
        }
    }
    let mut d_volume = 0.0;
    for v in &drho {
        d_volume += v;
    }
    d_volume *= h2;

    Ok(Evaluation {
        objective: compliance + spec.volume_weight * volume,
        gradient: d_compliance + spec.volume_weight * d_volume,
        cg_iters: iters,
    })
}

/// Central-difference reference for the gradient gate (2 extra solves).
///
/// # Errors
/// Propagates solver/cancellation failures.
pub fn central_difference(
    spec: &StudySpec,
    radius: f64,
    cancel: &Arc<AtomicBool>,
    cg_spent: &mut u64,
) -> Result<f64, String> {
    let dr = 1e-5;
    let plus = evaluate(spec, radius + dr, cancel, cg_spent)?;
    let minus = evaluate(spec, radius - dr, cancel, cg_spent)?;
    Ok((plus.objective - minus.objective) / (2.0 * dr))
}

// ---------------------------------------------------------------------------
// Mini deterministic executor: fixed-chunk parallel map, index-ordered
// merges, cancellation poll points. (fs-exec's two-lane executor replaces
// this pattern; chunk count is FIXED by problem size, not core count, so
// results are machine-independent.)
// ---------------------------------------------------------------------------

const FIXED_CHUNKS: usize = 4;

fn map_nodes_deterministic(
    n: usize,
    cancel: &Arc<AtomicBool>,
    f: impl Fn(f64, f64) -> f64 + Sync,
) -> Result<Vec<f64>, String> {
    let h = 1.0 / (n - 1) as f64;
    let total = n * n;
    let mut out = vec![0.0f64; total];
    let chunk = total.div_ceil(FIXED_CHUNKS);
    let cancelled = std::thread::scope(|s| {
        let mut handles = Vec::new();
        for (ci, slice) in out.chunks_mut(chunk).enumerate() {
            let f = &f;
            let cancel = Arc::clone(cancel);
            handles.push(s.spawn(move || {
                let start = ci * chunk;
                for (k, v) in slice.iter_mut().enumerate() {
                    let i = start + k;
                    // Poll point at row granularity (P7 in miniature).
                    if i.is_multiple_of(n) && cancel.load(Ordering::Relaxed) {
                        return true;
                    }
                    let (ix, iy) = (i % n, i / n);
                    *v = f(ix as f64 * h, iy as f64 * h);
                }
                false
            }));
        }
        handles.into_iter().any(|jh| jh.join().unwrap_or(true))
    });
    if cancelled {
        return Err("Cancelled: node map drained at a poll point (no torn state)".to_string());
    }
    Ok(out)
}

/// Matrix-free apply of K(rho): 5-point variable-coefficient stencil with
/// homogeneous Dirichlet boundary (boundary rows are identity·u = 0).
fn apply_k(n: usize, h: f64, rho: &[f64], u: &[f64], out: &mut [f64]) {
    let inv_h2 = 1.0 / (h * h);
    let idx = |ix: usize, iy: usize| iy * n + ix;
    for iy in 0..n {
        for ix in 0..n {
            let i = idx(ix, iy);
            if ix == 0 || iy == 0 || ix == n - 1 || iy == n - 1 {
                out[i] = u[i]; // Dirichlet row
                continue;
            }
            let mut acc = 0.0;
            for (jx, jy) in [(ix - 1, iy), (ix + 1, iy), (ix, iy - 1), (ix, iy + 1)] {
                let j = idx(jx, jy);
                let rho_e = f64::midpoint(rho[i], rho[j]);
                acc += rho_e * (u[i] - u[j]);
            }
            out[i] = acc * inv_h2;
        }
    }
}

/// Conjugate gradients with a cancellation poll each iteration and
/// deterministic (fixed-order) reductions.
fn cg_solve(
    n: usize,
    h: f64,
    rho: &[f64],
    spec: &StudySpec,
    cancel: &Arc<AtomicBool>,
) -> Result<(Vec<f64>, u64), String> {
    let total = n * n;
    let idx = |ix: usize, iy: usize| iy * n + ix;
    let mut f = vec![0.0f64; total];
    for iy in 1..n - 1 {
        for ix in 1..n - 1 {
            f[idx(ix, iy)] = 1.0;
        }
    }
    let dot = |a: &[f64], b: &[f64]| -> f64 {
        // Fixed index order: deterministic on any machine.
        let mut s = 0.0;
        for i in 0..a.len() {
            s += a[i] * b[i];
        }
        s
    };
    let mut u = vec![0.0f64; total];
    let mut r = f.clone(); // r = f - K·0
    let mut p = r.clone();
    let mut kp = vec![0.0f64; total];
    let mut rr = dot(&r, &r);
    let tol2 = 1e-22 * rr.max(1.0);
    let max_iter = 10 * total as u64;
    let mut iters = 0u64;
    while rr > tol2 && iters < max_iter {
        if cancel.load(Ordering::Relaxed) {
            return Err(format!(
                "Cancelled: CG drained after {iters} iterations at a poll point; partial \
                 state discarded, ledger untouched (request -> drain -> finalize)"
            ));
        }
        apply_k(n, h, rho, &p, &mut kp);
        let alpha = rr / dot(&p, &kp);
        for i in 0..total {
            u[i] += alpha * p[i];
            r[i] -= alpha * kp[i];
        }
        let rr_new = dot(&r, &r);
        let beta = rr_new / rr;
        rr = rr_new;
        for i in 0..total {
            p[i] = r[i] + beta * p[i];
        }
        iters += 1;
    }
    if rr > tol2 {
        return Err(format!(
            "SolverStalled: CG residual² {rr:.3e} above tolerance after {iters} iterations — \
             check conditioning (density floor eps) or raise the budget"
        ));
    }
    let _ = spec; // budget accounting is enforced by the caller across solves
    Ok((u, iters))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(grid: usize) -> StudySpec {
        StudySpec {
            name: "t".into(),
            seed: 1,
            grid,
            initial_radius: 0.25,
            opt_steps: 1,
            step_size: 0.1,
            volume_weight: 0.05,
            cg_budget: 1_000_000,
            r_min: 0.02,
            r_max: 0.45,
        }
    }

    #[test]
    fn constant_density_matches_poisson_reference() {
        // With rho ≡ 1 (radius→0 ⇒ material everywhere at ~1), the peak of u
        // for -Δu = 1 on the unit square is ≈ 0.0736713... (series solution).
        let n = 65;
        let h = 1.0 / (n - 1) as f64;
        let rho = vec![1.0; n * n];
        let cancel = Arc::new(AtomicBool::new(false));
        let (u, _) = cg_solve(n, h, &rho, &spec(n), &cancel).expect("converges");
        let peak = u.iter().copied().fold(0.0f64, f64::max);
        assert!(
            (peak - 0.07367).abs() < 5e-4,
            "peak {peak} vs series reference 0.0736713 (discretization tolerance)"
        );
    }

    #[test]
    fn adjoint_gradient_matches_central_difference() {
        let s = spec(33);
        let cancel = Arc::new(AtomicBool::new(false));
        let mut spent = 0;
        let e = evaluate(&s, 0.25, &cancel, &mut spent).expect("evaluates");
        let fd = central_difference(&s, 0.25, &cancel, &mut spent).expect("fd");
        let rel = (e.gradient - fd).abs() / e.gradient.abs().max(fd.abs()).max(1e-12);
        assert!(
            rel < 1e-4,
            "adjoint {} vs FD {fd}: rel {rel:.3e}",
            e.gradient
        );
    }

    #[test]
    fn evaluation_is_bitwise_deterministic() {
        let s = spec(33);
        let cancel = Arc::new(AtomicBool::new(false));
        let mut a_spent = 0;
        let mut b_spent = 0;
        let a = evaluate(&s, 0.3, &cancel, &mut a_spent).expect("a");
        let b = evaluate(&s, 0.3, &cancel, &mut b_spent).expect("b");
        assert_eq!(a.objective.to_bits(), b.objective.to_bits());
        assert_eq!(a.gradient.to_bits(), b.gradient.to_bits());
        assert_eq!(a_spent, b_spent);
    }

    #[test]
    fn cancellation_drains_cleanly() {
        let s = spec(129);
        let cancel = Arc::new(AtomicBool::new(true)); // pre-cancelled
        let mut spent = 0;
        let err = evaluate(&s, 0.25, &cancel, &mut spent).expect_err("must cancel");
        assert!(err.contains("Cancelled"), "{err}");
    }

    #[test]
    fn study_spec_parse_teaches_on_errors() {
        let good = r#"(study "pv" (seed 0x5EED0001) (grid 33)
            (budget (cg-iters 100000))
            (hole-radius 0.25) (opt-steps 3) (step-size 0.2) (volume-weight 0.05))"#;
        let s = StudySpec::parse(good).expect("parses");
        assert_eq!(s.grid, 33);
        assert_eq!(s.seed, 0x5EED_0001);
        // Missing budget must fail with the P4 message.
        let e = StudySpec::parse(
            r#"(study "pv" (seed 1) (grid 33) (hole-radius 0.25)
            (opt-steps 1) (step-size 0.1) (volume-weight 0.05))"#,
        )
        .expect_err("budget is mandatory");
        assert!(e.contains("budgets are mandatory"), "{e}");
        // Dimensioned radius must be rejected through the fs-qty path.
        let e = StudySpec::parse(
            r#"(study "pv" (seed 1) (grid 33) (budget (cg-iters 10))
            (hole-radius 25mm) (opt-steps 1) (step-size 0.1) (volume-weight 0.05))"#,
        )
        .expect_err("dimensioned radius");
        assert!(e.contains("dimensionless"), "{e}");
    }
}
