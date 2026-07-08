//! The nonlinear sizing pass: surviving members get cross-sections
//! from yield, compression members get EULER BUCKLING floors
//! (solid-square closed form `I = A²/12` ⇒ `A ≥ √(12|q|l²/(π²E))`),
//! joint parsimony prunes negligible members with a MANDATORY
//! least-squares equilibrium re-verification on the survivors, and
//! catalog snapping rounds areas UP (feasibility preserved by
//! construction) before the member-by-member code-check audit —
//! fs-constraint `Code` rows, every member listed, every check named.

use crate::ground::GroundStructure;
use crate::lp::LayoutLp;
use fs_constraint::ConstraintKind;
use std::fmt::Write as _;

/// One sized member.
#[derive(Debug, Clone)]
pub struct SizedMember {
    /// Member index in the ground structure.
    pub index: usize,
    /// End nodes.
    pub ends: (usize, usize),
    /// Signed axial force (positive = tension).
    pub force: f64,
    /// Length.
    pub length: f64,
    /// Continuous area from yield.
    pub area_yield: f64,
    /// Euler-governed floor (compression members; 0 for tension).
    pub area_buckling: f64,
    /// Snapped catalog area.
    pub area_catalog: f64,
}

/// The audit outcome.
#[derive(Debug, Clone, Default)]
pub struct CatalogAudit {
    /// Sized survivors.
    pub members: Vec<SizedMember>,
    /// Post-prune least-squares equilibrium residual (relative).
    pub eq_residual: f64,
    /// Members pruned by the parsimony threshold.
    pub pruned: usize,
    /// Member-by-member audit rows (JSON).
    pub rows: Vec<String>,
    /// Did every member pass every check post-snap?
    pub all_pass: bool,
}

/// Size the surviving members and snap to a catalog.
///
/// `x` is the split LP solution; members with net force below
/// `prune_frac`·max survive as pins only and are dropped. The pruned
/// topology's equilibrium is RE-VERIFIED by a least-squares refit of
/// forces on survivors (mandatory — pruning without re-verification is
/// how layouts stop being structures).
#[must_use]
#[allow(clippy::too_many_lines)] // size → prune → reverify → snap → audit
pub fn size_and_snap(
    gs: &GroundStructure,
    lp: &LayoutLp,
    x: &[f64],
    sigma_y: f64,
    youngs: f64,
    catalog: &[f64],
    prune_frac: f64,
) -> CatalogAudit {
    let m = gs.members.len();
    let forces: Vec<f64> = (0..m).map(|k| x[k] - x[m + k]).collect();
    let fmax = forces
        .iter()
        .fold(0.0f64, |a, &b| a.max(b.abs()))
        .max(1e-30);
    let survivors: Vec<usize> = (0..m)
        .filter(|&k| forces[k].abs() >= prune_frac * fmax)
        .collect();
    let pruned = m - survivors.len();
    // MANDATORY re-verification: least-squares refit of survivor
    // forces against the full load (normal equations on the reduced
    // equilibrium matrix, solved by CG on Bᵀ(B q − f) = 0 — small and
    // deterministic).
    let refit = refit_forces(gs, lp, &survivors);
    let (q_refit, eq_residual) = refit;
    let mut audit = CatalogAudit {
        eq_residual,
        pruned,
        ..CatalogAudit::default()
    };
    let code = ConstraintKind::Code {
        standard: "aisc-class-member-checks".to_string(),
    };
    let mut all_pass = true;
    for (si, &k) in survivors.iter().enumerate() {
        let q = q_refit[si];
        let l = gs.lengths[k];
        let area_yield = q.abs() / sigma_y;
        let area_buckling = if q < 0.0 {
            // Euler with pinned ends, solid square I = A²/12:
            // |q| ≤ π²EI/l² ⇒ A ≥ √(12|q|l²/(π²E)).
            (12.0 * q.abs() * l * l / (std::f64::consts::PI.powi(2) * youngs)).sqrt()
        } else {
            0.0
        };
        let need = area_yield.max(area_buckling);
        let area_catalog = catalog
            .iter()
            .copied()
            .find(|&a| a >= need)
            .unwrap_or(f64::NAN);
        // Post-snap re-verification (the checks the code row names).
        let stress_ok =
            area_catalog.is_finite() && q.abs() / area_catalog <= sigma_y * (1.0 + 1e-9);
        let buckling_ok = q >= 0.0
            || (area_catalog.is_finite()
                && q.abs()
                    <= std::f64::consts::PI.powi(2) * youngs * area_catalog * area_catalog
                        / 12.0
                        / (l * l)
                        * (1.0 + 1e-9));
        let pass = stress_ok && buckling_ok;
        all_pass &= pass;
        let mut row = String::new();
        let _ = write!(
            row,
            "{{\"member\":{k},\"force\":{q:.4e},\"len\":{l:.3},\"a_yield\":{area_yield:.4e},\
             \"a_euler\":{area_buckling:.4e},\"a_catalog\":{area_catalog:.4e},\
             \"code\":\"{}\",\"stress_ok\":{stress_ok},\"buckling_ok\":{buckling_ok}}}",
            match &code {
                ConstraintKind::Code { standard } => standard.as_str(),
                _ => unreachable!(),
            }
        );
        audit.rows.push(row);
        audit.members.push(SizedMember {
            index: k,
            ends: gs.members[k],
            force: q,
            length: l,
            area_yield,
            area_buckling,
            area_catalog,
        });
    }
    audit.all_pass = all_pass && eq_residual < 1e-6;
    audit
}

/// Least-squares force refit on a survivor subset: minimize
/// ‖B_s q − f‖ (normal equations by conjugate gradients, dense-free).
fn refit_forces(gs: &GroundStructure, lp: &LayoutLp, survivors: &[usize]) -> (Vec<f64>, f64) {
    let nrow = lp.b.len();
    let ns = survivors.len();
    // Column extraction: b_s columns of the SIGNED equilibrium matrix
    // (q⁺ columns of A, i.e. columns 0..m).
    let col = |k: usize, out: &mut Vec<f64>| {
        out.clear();
        out.resize(nrow, 0.0);
        // A's q⁺ column k = signed geometry column.
        let (a, b) = gs.members[k];
        let dx = (gs.nodes[b][0] - gs.nodes[a][0]) / gs.lengths[k];
        let dy = (gs.nodes[b][1] - gs.nodes[a][1]) / gs.lengths[k];
        for (dof, v) in [(2 * a, dx), (2 * a + 1, dy), (2 * b, -dx), (2 * b + 1, -dy)] {
            if let Some(row) = lp.dof_map[dof] {
                out[row] = v;
            }
        }
    };
    let matvec = |q: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0f64; nrow];
        let mut cbuf = Vec::new();
        for (si, &k) in survivors.iter().enumerate() {
            col(k, &mut cbuf);
            for (o, c) in out.iter_mut().zip(&cbuf) {
                *o += c * q[si];
            }
        }
        out
    };
    let rmatvec = |r: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0f64; ns];
        let mut cbuf = Vec::new();
        for (si, &k) in survivors.iter().enumerate() {
            col(k, &mut cbuf);
            out[si] = cbuf.iter().zip(r).map(|(c, r)| c * r).sum();
        }
        out
    };
    // CG on the normal equations BᵀB q = Bᵀ f.
    let bt_f = rmatvec(&lp.b);
    let mut q = vec![0.0f64; ns];
    let mut r = bt_f.clone();
    let mut p = r.clone();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    for _ in 0..4 * ns.max(32) {
        if rr.sqrt() < 1e-12 {
            break;
        }
        let bp = matvec(&p);
        let btbp = rmatvec(&bp);
        let pap: f64 = p.iter().zip(&btbp).map(|(a, b)| a * b).sum();
        if pap <= 0.0 {
            break;
        }
        let alpha = rr / pap;
        for i in 0..ns {
            q[i] += alpha * p[i];
            r[i] -= alpha * btbp[i];
        }
        let rr_new: f64 = r.iter().map(|v| v * v).sum();
        let beta = rr_new / rr;
        rr = rr_new;
        for i in 0..ns {
            p[i] = r[i] + beta * p[i];
        }
    }
    let ax = matvec(&q);
    let bnorm = lp.b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);
    let res = ax
        .iter()
        .zip(&lp.b)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt()
        / bnorm;
    (q, res)
}
