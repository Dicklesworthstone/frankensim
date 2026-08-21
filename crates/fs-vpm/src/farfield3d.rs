//! Far-field multipole aggregation + audited pruning (bead
//! wf-root-guzez.5.18.3, E4.7-iii). Mid-wake rows beyond a declared
//! age aggregate into multipole cells (monopole + dipole + quadrupole
//! about the circulation-weighted centroid); the near wake stays an
//! exact filament lattice, so the hybrid evaluation is
//! `near.induced_velocity(p) + far.eval(p)`.
//!
//! Pruning is CELL-BOUNDED and AUDITED: a cell may be dropped only
//! when its rigorous contribution bound at ALL registered probe
//! points is below the tolerance, and the audit then recomputes an
//! EXACT spot-check per pruned cell per probe — the audit distrusts
//! the bound. `pruning-certificate-failed` TERMINATES the prune (no
//! mutation), which the battery proves by executing an adversarial
//! bound-deflation driver.
//!
//! `WakeCoreEvolutionMode` is a registered enum that ENTERS IDENTITY
//! (its bytes are hashed into the far-field digest) and physically
//! feeds the pruning bounds through the spread radius.

use crate::filament3d::{FilamentWake, Refusal, segment_velocity};
use fs_blake3::hash_domain;

/// Registered core-evolution mode (enters identity).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WakeCoreEvolutionMode {
    /// Cores frozen at their shed/coarsened second moment.
    Frozen,
    /// Lamb–Oseen-class viscous spreading: a row of age `t` gains
    /// `4·nu_eff·t` of core second moment [m²].
    CoreSpreading {
        /// Effective eddy viscosity [m²/s].
        nu_eff_m2ps: f64,
    },
}

/// Spread-viscosity cap [m²/s].
pub const MAX_NU_EFF: f64 = 1.0;

impl WakeCoreEvolutionMode {
    /// Admit the mode.
    ///
    /// # Errors
    /// `core-mode-invalid` (nu_eff must be finite, > 0, <= cap —
    /// refused at cap AND cap+ulp class; Frozen always admits).
    pub fn admit(&self) -> Result<(), Refusal> {
        match self {
            WakeCoreEvolutionMode::Frozen => Ok(()),
            WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps } => {
                if nu_eff_m2ps.is_finite() && *nu_eff_m2ps > 0.0 && *nu_eff_m2ps <= MAX_NU_EFF {
                    Ok(())
                } else {
                    Err(Refusal {
                        code: "core-mode-invalid",
                        message: format!("nu_eff {nu_eff_m2ps}"),
                        ranked_repairs: vec![format!("nu_eff in (0, {MAX_NU_EFF}]")],
                    })
                }
            }
        }
    }

    /// Identity bytes (hashed into the far-field digest).
    #[must_use]
    pub fn id_bytes(&self) -> Vec<u8> {
        match self {
            WakeCoreEvolutionMode::Frozen => vec![0u8],
            WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps } => {
                let mut b = vec![1u8];
                b.extend_from_slice(&nu_eff_m2ps.to_bits().to_le_bytes());
                b
            }
        }
    }
}

/// One aggregated multipole cell.
#[derive(Clone, Debug, PartialEq)]
pub struct FarCell {
    /// Circulation-weighted centroid.
    pub centroid: [f64; 3],
    /// Bounding radius about the centroid [m] (geometry only).
    pub radius_m: f64,
    /// Monopole: Σ Γ·dl (vector).
    pub w_sum: [f64; 3],
    /// Σ |Γ|·|dl| — the rigorous bound numerator.
    pub w_abs_sum: f64,
    /// Dipole moment D[j][l] = Σ w_j·δ_l.
    pub dip: [[f64; 3]; 3],
    /// Quadrupole moment Q[j][l][m] = Σ w_j·δ_l·δ_m.
    pub quad: [[[f64; 3]; 3]; 3],
    /// Core second moment carried into the cell [m²] (mode-evolved).
    pub core2_m2: f64,
    /// Exact source segments (a, b, Γ) RETAINED for the audit's exact
    /// spot-check — the audit never trusts the expansion.
    pub segs: Vec<([f64; 3], [f64; 3], f64)>,
}

impl FarCell {
    /// Effective standoff radius: geometry + core spread.
    #[must_use]
    pub fn effective_radius_m(&self) -> f64 {
        self.radius_m + self.core2_m2.max(0.0).sqrt()
    }

    /// Rigorous contribution bound at `p`: |Γ dl × r|/|r|³ ≤ |Γ||dl|/d²
    /// with d = |p − centroid| − effective radius. Non-positive
    /// standoff ⇒ unboundable ⇒ +inf (never prunable).
    #[must_use]
    pub fn contribution_bound(&self, p: [f64; 3]) -> f64 {
        let d = dist(p, self.centroid) - self.effective_radius_m();
        if d <= 0.0 {
            return f64::INFINITY;
        }
        self.w_abs_sum / (4.0 * core::f64::consts::PI * d * d)
    }

    /// Multipole evaluation (monopole + dipole + quadrupole) at `p`.
    #[must_use]
    pub fn eval(&self, p: [f64; 3]) -> [f64; 3] {
        let r = [
            p[0] - self.centroid[0],
            p[1] - self.centroid[1],
            p[2] - self.centroid[2],
        ];
        let rho2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        let rho = rho2.sqrt();
        if rho < 1e-9 {
            return [0.0; 3];
        }
        let rho3 = rho2 * rho;
        let rho5 = rho3 * rho2;
        let rho7 = rho5 * rho2;
        // f(r) = r/ρ³ and its first/second derivatives.
        let f0 = [r[0] / rho3, r[1] / rho3, r[2] / rho3];
        let kd = |i: usize, j: usize| -> f64 { if i == j { 1.0 } else { 0.0 } };
        let jmat = |k: usize, l: usize| -> f64 { kd(k, l) / rho3 - 3.0 * r[k] * r[l] / rho5 };
        let hmat = |i: usize, j: usize, k: usize| -> f64 {
            -3.0 * (kd(i, j) * r[k] + kd(i, k) * r[j] + kd(j, k) * r[i]) / rho5
                + 15.0 * r[i] * r[j] * r[k] / rho7
        };
        // Monopole: W × f(r0).
        let mut v = cross(self.w_sum, f0);
        // Dipole: −Σ w × (J δ)  →  −ε_ijk D[j][l] J[k][l].
        for i in 0..3 {
            let (j1, k1, j2, k2) = eps_pairs(i);
            let mut acc = 0.0;
            for l in 0..3 {
                acc += self.dip[j1][l] * jmat(k1, l) - self.dip[j2][l] * jmat(k2, l);
            }
            v[i] -= acc;
        }
        // Quadrupole: ½ Σ w × (H:δδ)  →  ½ ε_ijk Q[j][l][m] H[k][l][m].
        for i in 0..3 {
            let (j1, k1, j2, k2) = eps_pairs(i);
            let mut acc = 0.0;
            for l in 0..3 {
                for m in 0..3 {
                    acc +=
                        self.quad[j1][l][m] * hmat(k1, l, m) - self.quad[j2][l][m] * hmat(k2, l, m);
                }
            }
            v[i] += 0.5 * acc;
        }
        let k = 1.0 / (4.0 * core::f64::consts::PI);
        [k * v[0], k * v[1], k * v[2]]
    }

    /// EXACT contribution at `p` from the retained source segments
    /// (audit path — bypasses the expansion entirely).
    #[must_use]
    pub fn exact(&self, p: [f64; 3]) -> [f64; 3] {
        let mut v = [0.0f64; 3];
        for (a, b, g) in &self.segs {
            let sv = segment_velocity(p, *a, *b);
            v[0] += g * sv[0];
            v[1] += g * sv[1];
            v[2] += g * sv[2];
        }
        v
    }
}

/// One pruned-cell audit row (receipt ingredient).
#[derive(Clone, Debug, PartialEq)]
pub struct PrunedCellRecord {
    /// Cell index at prune time.
    pub cell_index: usize,
    /// Worst (largest) bound across probes.
    pub bound_worst: f64,
    /// Worst EXACT spot-check magnitude across probes.
    pub exact_worst: f64,
}

/// The prune receipt (what was dropped, with its audit evidence).
#[derive(Clone, Debug, PartialEq)]
pub struct PruneReceipt {
    /// Tolerance the certificate enforced [m/s per cell per probe].
    pub tol_mps: f64,
    /// Probes registered for the certificate.
    pub n_probes: usize,
    /// Audit rows for every pruned cell.
    pub pruned: Vec<PrunedCellRecord>,
    /// Cells kept.
    pub kept: usize,
}

/// The aggregated far field.
#[derive(Clone, Debug, PartialEq)]
pub struct FarField {
    /// Registered core-evolution mode (enters identity).
    pub mode: WakeCoreEvolutionMode,
    /// Multipole cells, oldest first.
    pub cells: Vec<FarCell>,
}

impl FarField {
    /// Aggregate the OLDEST `n_rows` of `wake` into cells of
    /// `rows_per_cell` consecutive rows each, removing them from the
    /// near wake. Trailing segments are taken against the pre-removal
    /// successor rows, so no segment is lost or double-counted.
    ///
    /// # Errors
    /// `farfield-invalid` (n_rows must be a positive multiple of
    /// rows_per_cell strictly below the row count — exact multiple
    /// admits, ±1 row refuses; dt must be finite positive); mode
    /// refusals.
    pub fn aggregate(
        wake: &mut FilamentWake,
        n_rows: usize,
        rows_per_cell: usize,
        mode: WakeCoreEvolutionMode,
        dt_s: f64,
    ) -> Result<FarField, Refusal> {
        mode.admit()?;
        if n_rows == 0
            || rows_per_cell == 0
            || n_rows >= wake.rows.len()
            || n_rows % rows_per_cell != 0
            || !(dt_s.is_finite() && dt_s > 0.0)
        {
            return Err(Refusal {
                code: "farfield-invalid",
                message: format!(
                    "n_rows {n_rows} rows_per_cell {rows_per_cell} of {} rows, dt {dt_s}",
                    wake.rows.len()
                ),
                ranked_repairs: vec![
                    "aggregate a positive exact multiple of rows_per_cell, retaining >= 1 near row"
                        .into(),
                ],
            });
        }
        let n_st = wake.cert.n_stations;
        let total = wake.rows.len();
        let mut cells = Vec::with_capacity(n_rows / rows_per_cell);
        for c0 in (0..n_rows).step_by(rows_per_cell) {
            // Collect the cell's exact segments (spanwise + trailing).
            let mut segs: Vec<([f64; 3], [f64; 3], f64)> = Vec::new();
            let mut core2_num = 0.0;
            let mut core2_den = 0.0;
            for ri in c0..c0 + rows_per_cell {
                let row = &wake.rows[ri];
                for s in 0..n_st {
                    if row.gamma[s] != 0.0 {
                        segs.push((row.nodes[s], row.nodes[s + 1], row.gamma[s]));
                    }
                }
                let next_nodes = if ri + 1 < total {
                    &wake.rows[ri + 1].nodes
                } else {
                    &wake.line_nodes
                };
                for e in 0..=n_st {
                    let left = if e > 0 { row.gamma[e - 1] } else { 0.0 };
                    let right = if e < n_st { row.gamma[e] } else { 0.0 };
                    let g = left - right;
                    if g != 0.0 {
                        segs.push((row.nodes[e], next_nodes[e], g));
                    }
                }
                // Mode-evolved core second moment (age = ticks since
                // shed; oldest rows are oldest).
                let age_s = (total - ri) as f64 * dt_s;
                let spread = match mode {
                    WakeCoreEvolutionMode::Frozen => 0.0,
                    WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps } => {
                        4.0 * nu_eff_m2ps * age_s
                    }
                };
                let w_row: f64 = row.gamma.iter().map(|g| g.abs()).sum::<f64>().max(1e-12);
                core2_num += w_row * (row.core2_m2 + spread);
                core2_den += w_row;
            }
            // Circulation-weighted centroid over segment midpoints.
            let mut cw = [0.0f64; 3];
            let mut wt = 0.0;
            for (a, b, g) in &segs {
                let w = g.abs() * dist(*a, *b);
                for k in 0..3 {
                    cw[k] += w * 0.5 * (a[k] + b[k]);
                }
                wt += w;
            }
            let wt = wt.max(1e-12);
            let centroid = [cw[0] / wt, cw[1] / wt, cw[2] / wt];
            // Moments about the centroid.
            let mut w_sum = [0.0f64; 3];
            let mut w_abs_sum = 0.0;
            let mut dip = [[0.0f64; 3]; 3];
            let mut quad = [[[0.0f64; 3]; 3]; 3];
            let mut radius: f64 = 0.0;
            for (a, b, g) in &segs {
                let w = [g * (b[0] - a[0]), g * (b[1] - a[1]), g * (b[2] - a[2])];
                let mid = [
                    0.5 * (a[0] + b[0]),
                    0.5 * (a[1] + b[1]),
                    0.5 * (a[2] + b[2]),
                ];
                let d = [
                    mid[0] - centroid[0],
                    mid[1] - centroid[1],
                    mid[2] - centroid[2],
                ];
                w_abs_sum += g.abs() * dist(*a, *b);
                for j in 0..3 {
                    w_sum[j] += w[j];
                    for l in 0..3 {
                        dip[j][l] += w[j] * d[l];
                        for m in 0..3 {
                            quad[j][l][m] += w[j] * d[l] * d[m];
                        }
                    }
                }
                radius = radius.max(dist(*a, centroid)).max(dist(*b, centroid));
            }
            cells.push(FarCell {
                centroid,
                radius_m: radius,
                w_sum,
                w_abs_sum,
                dip,
                quad,
                core2_m2: core2_num / core2_den.max(1e-12),
                segs,
            });
        }
        wake.rows.drain(0..n_rows);
        Ok(FarField { mode, cells })
    }

    /// Far-field induced velocity at `p` (sum of cell expansions).
    #[must_use]
    pub fn eval(&self, p: [f64; 3]) -> [f64; 3] {
        let mut v = [0.0f64; 3];
        for c in &self.cells {
            let cv = c.eval(p);
            v[0] += cv[0];
            v[1] += cv[1];
            v[2] += cv[2];
        }
        v
    }

    /// Audited prune: drop every cell whose contribution bound at ALL
    /// probes is < `tol_mps`, with a per-cell per-probe EXACT
    /// spot-check that distrusts the bound.
    ///
    /// # Errors
    /// `prune-invalid` (no probes, or tol not finite-positive);
    /// `pruning-certificate-failed` (an exact spot-check reached the
    /// tolerance a bound promised it could not — TERMINATES with no
    /// mutation).
    pub fn prune_audited(
        &mut self,
        probes: &[[f64; 3]],
        tol_mps: f64,
    ) -> Result<PruneReceipt, Refusal> {
        self.prune_audited_scaled(probes, tol_mps, 1.0)
    }

    /// Battery falsifier driver ONLY: deflating the bound by `scale`
    /// forges over-eager prune candidacy so the audit's certificate
    /// failure path is EXECUTED, not assumed.
    #[doc(hidden)]
    pub fn prune_audited_scaled(
        &mut self,
        probes: &[[f64; 3]],
        tol_mps: f64,
        bound_scale: f64,
    ) -> Result<PruneReceipt, Refusal> {
        if probes.is_empty() || !(tol_mps.is_finite() && tol_mps > 0.0) {
            return Err(Refusal {
                code: "prune-invalid",
                message: format!("{} probes, tol {tol_mps}", probes.len()),
                ranked_repairs: vec!["register >= 1 probe and a finite positive tolerance".into()],
            });
        }
        // Phase 1: candidacy + audit (NO mutation until fully certified).
        let mut drop_flags = vec![false; self.cells.len()];
        let mut records = Vec::new();
        for (ci, cell) in self.cells.iter().enumerate() {
            let mut bound_worst = 0.0f64;
            let mut candidate = true;
            for p in probes {
                let b = cell.contribution_bound(*p) * bound_scale;
                bound_worst = bound_worst.max(b);
                if b >= tol_mps {
                    candidate = false;
                    break;
                }
            }
            if !candidate {
                continue;
            }
            // Audit: EXACT spot-check per probe — the certificate.
            let mut exact_worst = 0.0f64;
            for p in probes {
                let e = cell.exact(*p);
                let mag = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
                exact_worst = exact_worst.max(mag);
                if mag >= tol_mps {
                    return Err(Refusal {
                        code: "pruning-certificate-failed",
                        message: format!(
                            "cell {ci}: exact {mag} >= tol {tol_mps} though bound {bound_worst} \
                             promised otherwise"
                        ),
                        ranked_repairs: vec![
                            "the bound implementation is unsound for this cell; do not prune"
                                .into(),
                        ],
                    });
                }
            }
            drop_flags[ci] = true;
            records.push(PrunedCellRecord {
                cell_index: ci,
                bound_worst,
                exact_worst,
            });
        }
        // Phase 2: certified — mutate.
        let mut kept = Vec::with_capacity(self.cells.len());
        for (ci, cell) in self.cells.drain(..).enumerate() {
            if !drop_flags[ci] {
                kept.push(cell);
            }
        }
        self.cells = kept;
        Ok(PruneReceipt {
            tol_mps,
            n_probes: probes.len(),
            pruned: records,
            kept: self.cells.len(),
        })
    }

    /// Content digest (mode identity bytes + all cell moments).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut b = self.mode.id_bytes();
        for c in &self.cells {
            for v in c
                .centroid
                .iter()
                .chain(c.w_sum.iter())
                .chain([c.radius_m, c.w_abs_sum, c.core2_m2].iter())
            {
                b.extend_from_slice(&v.to_bits().to_le_bytes());
            }
            for j in 0..3 {
                for l in 0..3 {
                    b.extend_from_slice(&c.dip[j][l].to_bits().to_le_bytes());
                    for m in 0..3 {
                        b.extend_from_slice(&c.quad[j][l][m].to_bits().to_le_bytes());
                    }
                }
            }
        }
        hash_domain("org.frankensim.fs-vpm.farfield3d.v1", &b).to_hex()
    }
}

/// Hybrid induced velocity: exact near wake + far multipoles.
#[must_use]
pub fn hybrid_velocity(near: &FilamentWake, far: &FarField, p: [f64; 3]) -> [f64; 3] {
    let nv = near.induced_velocity(p);
    let fv = far.eval(p);
    [nv[0] + fv[0], nv[1] + fv[1], nv[2] + fv[2]]
}

/// V-10 receipt: hybrid wake buildup vs the fs-wakeref dense reference
/// on the overlap window. SHAPE-CLASS comparison (each series
/// normalized by its own terminal value) with honest scale notes — the
/// observables differ (wake-induced downwash buildup vs lift buildup),
/// so scale agreement is out of scope BY DECLARATION and only the
/// normalized transient shape is compared. Tier A/B KPI deltas are
/// REPORTED (V-05c), never forced to vanish.
#[derive(Clone, Debug, PartialEq)]
pub struct V10Receipt {
    /// Schema id.
    pub schema: &'static str,
    /// Core-evolution mode identity hex.
    pub mode_digest: String,
    /// Overlap window compared [ticks].
    pub overlap_ticks: usize,
    /// RMS of (hybrid_shape − reference_shape) over the window.
    pub shape_rms: f64,
    /// |terminal hybrid − terminal reference| (both ≡ 1 by
    /// construction; recorded as a non-vacuity witness).
    pub terminal_delta: f64,
    /// Tier A KPI: probe speed with the FULL exact wake [m/s].
    pub tier_a_kpi_mps: f64,
    /// Tier B KPI: probe speed with the coarsened+multipole+pruned
    /// hybrid [m/s].
    pub tier_b_kpi_mps: f64,
    /// Reported delta (B − A), never forced to vanish.
    pub kpi_delta_mps: f64,
    /// Honest scale note.
    pub scale_note: &'static str,
    /// Digest over the receipt payload.
    pub receipt_digest: String,
}

/// Declared v10 scale note.
pub const V10_SCALE_NOTE: &str = "shape-class only: fs-vpm normalized wake-downwash buildup vs \
     fs-wakeref normalized canard-lift buildup (Wagner-class); absolute \
     scales are different observables and are NOT compared";

/// Emit the V-10 receipt.
///
/// # Errors
/// `v10-invalid` (length mismatch, empty, or non-finite series).
pub fn emit_v10_receipt(
    mode: WakeCoreEvolutionMode,
    hybrid_shape: &[f64],
    reference_shape: &[f64],
    tier_a_kpi_mps: f64,
    tier_b_kpi_mps: f64,
) -> Result<V10Receipt, Refusal> {
    if hybrid_shape.is_empty()
        || hybrid_shape.len() != reference_shape.len()
        || hybrid_shape
            .iter()
            .chain(reference_shape.iter())
            .any(|v| !v.is_finite())
        || !(tier_a_kpi_mps.is_finite() && tier_b_kpi_mps.is_finite())
    {
        return Err(Refusal {
            code: "v10-invalid",
            message: format!(
                "hybrid {} vs reference {} samples",
                hybrid_shape.len(),
                reference_shape.len()
            ),
            ranked_repairs: vec!["equal-length finite normalized series over the overlap".into()],
        });
    }
    let n = hybrid_shape.len();
    let mut sum2 = 0.0;
    for i in 0..n {
        let d = hybrid_shape[i] - reference_shape[i];
        sum2 += d * d;
    }
    let shape_rms = (sum2 / n as f64).sqrt();
    let terminal_delta = (hybrid_shape[n - 1] - reference_shape[n - 1]).abs();
    let kpi_delta_mps = tier_b_kpi_mps - tier_a_kpi_mps;
    let mode_digest = hash_domain("org.frankensim.fs-vpm.core-mode.v1", &mode.id_bytes()).to_hex();
    let mut b = mode.id_bytes();
    for v in hybrid_shape.iter().chain(reference_shape.iter()).chain(
        [
            tier_a_kpi_mps,
            tier_b_kpi_mps,
            kpi_delta_mps,
            shape_rms,
            terminal_delta,
        ]
        .iter(),
    ) {
        b.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    Ok(V10Receipt {
        schema: "org.frankensim.wf.v10-receipt.v1",
        mode_digest,
        overlap_ticks: n,
        shape_rms,
        terminal_delta,
        tier_a_kpi_mps,
        tier_b_kpi_mps,
        kpi_delta_mps,
        scale_note: V10_SCALE_NOTE,
        receipt_digest: hash_domain("org.frankensim.wf.v10-receipt.v1", &b).to_hex(),
    })
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// For output index i, (w×u)_i = w[j1]·u[k1] − w[j2]·u[k2].
fn eps_pairs(i: usize) -> (usize, usize, usize, usize) {
    match i {
        0 => (1, 2, 2, 1),
        1 => (2, 0, 0, 2),
        _ => (0, 1, 1, 0),
    }
}
