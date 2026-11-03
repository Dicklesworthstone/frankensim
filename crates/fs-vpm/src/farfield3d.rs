//! Far-field multipole aggregation + audited pruning (bead
//! wf-root-guzez.5.18.3, E4.7-iii). Mid-wake rows beyond a declared
//! age aggregate into multipole cells (monopole + dipole + quadrupole
//! about the circulation-weighted centroid); the near wake stays an
//! exact filament lattice, so the hybrid evaluation is
//! `hybrid_velocity(near, far, p)`, applying the far field's declared
//! viscosity and shedding interval to retained near-row ages too.
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
//! feeds near/exact and far velocities as well as the pruning bounds.

use crate::filament3d::{FilamentWake, Refusal, segment_velocity_with_core};
use fs_blake3::hash_domain;

#[path = "farfield3d/advection.rs"]
mod advection;

/// Registered core-evolution mode (enters identity).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WakeCoreEvolutionMode {
    /// Cores frozen at their shed/coarsened width squared.
    Frozen,
    /// Reduced algebraic core spreading: a row of age `t` gains
    /// `4·nu_eff·t` of Rosenhead–Moore core width squared [m²].
    /// This is not exact Gaussian Lamb–Oseen diffusion.
    CoreSpreading {
        /// Effective kinematic viscosity [m²/s], molecular plus any declared eddy term.
        nu_eff_m2ps: f64,
    },
}

/// Spread-viscosity cap [m²/s].
pub const MAX_NU_EFF: f64 = 1.0;

impl WakeCoreEvolutionMode {
    /// Admit the mode.
    ///
    /// # Errors
    /// `core-mode-invalid` (nu_eff must be finite, >= 0, <= cap;
    /// the exact cap admits and larger values refuse; Frozen always admits).
    pub fn admit(&self) -> Result<(), Refusal> {
        match self {
            WakeCoreEvolutionMode::Frozen => Ok(()),
            WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps } => {
                if nu_eff_m2ps.is_finite() && *nu_eff_m2ps >= 0.0 && *nu_eff_m2ps <= MAX_NU_EFF {
                    Ok(())
                } else {
                    Err(Refusal {
                        code: "core-mode-invalid",
                        message: format!("nu_eff {nu_eff_m2ps}"),
                        ranked_repairs: vec![format!("nu_eff in [0, {MAX_NU_EFF}]")],
                    })
                }
            }
        }
    }

    /// Core-width-squared growth over an elapsed interval [s].
    pub(crate) fn core_growth_m2(&self, dt_s: f64) -> f64 {
        match self {
            Self::Frozen => 0.0,
            Self::CoreSpreading { nu_eff_m2ps } if *nu_eff_m2ps == 0.0 => 0.0,
            Self::CoreSpreading { nu_eff_m2ps } => (4.0 * nu_eff_m2ps) * dt_s,
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
    /// Circulation-length-weighted mean expansion core width squared [m²].
    /// Exact evaluation retains each source's own core below.
    pub core2_m2: f64,
    /// Source core width squared for each entry of `segs`, in the same order.
    pub source_core2_m2: Vec<f64>,
    /// Σ Γ dl (source_core² - core2_m2), the leading core-variation moment.
    pub core_moment: [f64; 3],
    /// Largest absolute deviation from the expansion's common core [m²].
    pub max_abs_core_delta_m2: f64,
    /// Common extra diffusion after this cell was aggregated [m²].
    pub diffusion_offset_m2: f64,
    /// Maximum base source core [m²], for O(1) advance validation.
    pub max_source_core2_m2: f64,
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
        let distance2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        // A far expansion is not a near-field solver.
        if distance2 <= self.radius_m * self.radius_m {
            return self.exact(p);
        }
        let rho2 = distance2 + self.core2_m2;
        // Expansion-domain guard only, not a certified total-velocity error.
        // Δcore² counts as spatial order two; large spreads need exact sources.
        if self.max_abs_core_delta_m2 > rho2 / 16.0 {
            return self.exact(p);
        }
        let rho = rho2.sqrt();
        if rho == 0.0 {
            return [0.0; 3];
        }
        let rho3 = rho2 * rho;
        let rho5 = rho3 * rho2;
        let rho7 = rho5 * rho2;
        // f(r) = r/(|r|² + core²)^(3/2) and its spatial derivatives.
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
        let core_correction = cross(self.core_moment, r);
        for i in 0..3 {
            v[i] -= 1.5 * core_correction[i] / rho5;
        }
        let k = 1.0 / (4.0 * core::f64::consts::PI);
        [k * v[0], k * v[1], k * v[2]]
    }

    /// EXACT contribution at `p` from the retained source segments
    /// (audit path — bypasses the expansion entirely).
    #[must_use]
    pub fn exact(&self, p: [f64; 3]) -> [f64; 3] {
        assert_eq!(
            self.segs.len(),
            self.source_core2_m2.len(),
            "one core per source segment"
        );
        let mut v = [0.0f64; 3];
        for ((a, b, g), core2) in self.segs.iter().zip(&self.source_core2_m2) {
            let sv = segment_velocity_with_core(p, *a, *b, *core2 + self.diffusion_offset_m2);
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
    /// Seconds per shedding interval used by retained near-row ages.
    pub shed_dt_s: f64,
    /// Multipole cells, oldest first.
    pub cells: Vec<FarCell>,
}

impl FarField {
    /// Aggregate the OLDEST `n_rows` of `wake` into cells of
    /// `rows_per_cell` consecutive rows each, removing them from the
    /// near wake. Physical edges come from the exact near-wake traversal.
    /// Trailing segments end at the pre-removal successors, and the last
    /// removed quad's strengths remain on the near wake's downstream
    /// interface. The shared spanwise edge is neither lost nor duplicated.
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
            || wake.downstream_gamma.as_ref().is_some_and(|g| {
                g.len() != wake.cert.n_stations || g.iter().any(|v| !v.is_finite())
            })
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
        let growth_per_tick = mode.core_growth_m2(dt_s);
        if !growth_per_tick.is_finite() {
            return Err(Refusal {
                code: "farfield-core-invalid",
                message: "viscosity × shedding interval is not representable".into(),
                ranked_repairs: vec!["use a representable viscosity and shedding interval".into()],
            });
        }
        for (ri, row) in wake.rows.iter().enumerate() {
            let grown = if growth_per_tick == 0.0 {
                row.core2_m2
            } else {
                row.core2_m2 + growth_per_tick * row.age_ticks
            };
            if !row.age_ticks.is_finite()
                || row.age_ticks < 0.0
                || !row.core2_m2.is_finite()
                || row.core2_m2 < 0.0
                || !grown.is_finite()
            {
                return Err(Refusal {
                    code: "farfield-core-invalid",
                    message: format!("row {ri}: invalid age or evolved core"),
                    ranked_repairs: vec![
                        "retain finite nonnegative row ages and core widths".into(),
                    ],
                });
            }
        }
        let mut cells = Vec::with_capacity(n_rows / rows_per_cell);
        for c0 in (0..n_rows).step_by(rows_per_cell) {
            // Collect the cell's exact segments (spanwise + trailing).
            let mut segs = Vec::new();
            let mut source_core2_m2 = Vec::new();
            for ri in c0..c0 + rows_per_cell {
                wake.visit_row_segments(ri, growth_per_tick, &mut |segment| {
                    segs.push((segment.start, segment.end, segment.circulation));
                    source_core2_m2.push(segment.core2_m2);
                });
            }
            let cell = cell_from_segs_with_cores(segs, source_core2_m2);
            let scalars = [
                cell.radius_m,
                cell.w_abs_sum,
                cell.core2_m2,
                cell.max_abs_core_delta_m2,
            ];
            let finite = cell
                .centroid
                .iter()
                .chain(cell.w_sum.iter())
                .chain(cell.core_moment.iter())
                .chain(scalars.iter())
                .chain(cell.dip.iter().flatten())
                .chain(cell.quad.iter().flatten().flatten())
                .all(|x| x.is_finite());
            if !finite {
                return Err(Refusal {
                    code: "farfield-core-invalid",
                    message: format!("cell starting at row {c0}: non-finite source/core moments"),
                    ranked_repairs: vec![
                        "reduce cell size or rescale circulation and length units".into(),
                    ],
                });
            }
            cells.push(cell);
        }
        wake.downstream_gamma = Some(wake.rows[n_rows - 1].gamma.clone());
        wake.rows.drain(0..n_rows);
        Ok(FarField {
            mode,
            shed_dt_s: dt_s,
            cells,
        })
    }

    /// Add 4ν dt to every retained source core without rebuilding moments.
    /// Spatial sources remain fixed. Callers advance geometry separately.
    /// Zero dt, Frozen, and zero viscosity are exact identity operations.
    ///
    /// # Errors
    /// Invalid elapsed time/mode, invalid cached core state, or any width
    /// overflow. Every refusal leaves all cells unchanged.
    pub fn advance_core(&mut self, dt_s: f64) -> Result<(), Refusal> {
        let growth = self.admit_core_advance(dt_s)?;
        self.apply_core_growth(growth);
        Ok(())
    }

    fn admit_core_advance(&self, dt_s: f64) -> Result<f64, Refusal> {
        if !dt_s.is_finite() || dt_s < 0.0 {
            return Err(Refusal {
                code: "core-advance-invalid",
                message: format!("elapsed time must be finite and nonnegative: {dt_s}"),
                ranked_repairs: vec!["supply elapsed seconds >= 0".into()],
            });
        }
        self.mode.admit()?;
        let growth = self.mode.core_growth_m2(dt_s);
        if !growth.is_finite() {
            return Err(Refusal {
                code: "core-advance-overflow",
                message: "viscosity × elapsed time is not representable".into(),
                ranked_repairs: vec!["reduce the elapsed interval or rescale units".into()],
            });
        }
        if growth == 0.0 {
            return Ok(0.0);
        }
        for (ci, cell) in self.cells.iter().enumerate() {
            let next_core = cell.core2_m2 + growth;
            let next_offset = cell.diffusion_offset_m2 + growth;
            let max_effective_core = cell.max_source_core2_m2 + next_offset;
            let scalars = [
                cell.core2_m2,
                cell.diffusion_offset_m2,
                cell.max_source_core2_m2,
                cell.max_abs_core_delta_m2,
            ];
            if cell.segs.len() != cell.source_core2_m2.len()
                || scalars.iter().any(|v| !v.is_finite() || *v < 0.0)
                || cell.core_moment.iter().any(|v| !v.is_finite())
                || !next_core.is_finite()
                || !next_offset.is_finite()
                || !max_effective_core.is_finite()
            {
                return Err(Refusal {
                    code: "core-advance-overflow",
                    message: format!("cell {ci}: invalid or unrepresentable evolved core"),
                    ranked_repairs: vec![
                        "use valid cell cores and a representable elapsed interval".into(),
                    ],
                });
            }
        }
        Ok(growth)
    }

    fn apply_core_growth(&mut self, growth: f64) {
        if growth != 0.0 {
            for cell in &mut self.cells {
                cell.core2_m2 += growth;
                cell.diffusion_offset_m2 += growth;
            }
        }
    }

    /// Translate cell centers and all retained exact source endpoints by one
    /// uniform displacement [m]. Relative multipole moments, radii and cores
    /// stay unchanged; no cell reconstruction or source allocation is needed.
    /// Zero displacement is an exact identity. Rebuild cached images afterward.
    ///
    /// # Errors
    /// `farfield-translation-invalid` for non-finite displacement/geometry or
    /// coordinate overflow. Every cell is checked before any mutation.
    pub fn translate(&mut self, displacement_m: [f64; 3]) -> Result<(), Refusal> {
        self.admit_translation(displacement_m)?;
        self.apply_translation(displacement_m);
        Ok(())
    }

    fn admit_translation(&self, delta: [f64; 3]) -> Result<(), Refusal> {
        let representable = delta.iter().all(|v| v.is_finite())
            && self.cells.iter().all(|cell| {
                std::iter::once(&cell.centroid)
                    .chain(cell.segs.iter().flat_map(|(a, b, _)| [a, b]))
                    .all(|p| (0..3).all(|i| (p[i] + delta[i]).is_finite()))
            });
        if !representable {
            return Err(Refusal {
                code: "farfield-translation-invalid",
                message: "uniform translation must retain finite cell/source coordinates".into(),
                ranked_repairs: vec!["use finite geometry and a representable displacement".into()],
            });
        }
        Ok(())
    }

    fn apply_translation(&mut self, delta: [f64; 3]) {
        if delta.iter().all(|v| *v == 0.0) {
            return;
        }
        for cell in &mut self.cells {
            for p in std::iter::once(&mut cell.centroid)
                .chain(cell.segs.iter_mut().flat_map(|(a, b, _)| [a, b]))
            {
                for i in 0..3 {
                    if delta[i] != 0.0 {
                        p[i] += delta[i];
                    }
                }
            }
        }
    }

    /// Advance one prescribed shedding interval as one near/far transaction.
    /// Existing sources translate by `displacement_m` [m]; existing far cores
    /// diffuse over `shed_dt_s`, while near-row ages advance one tick. The new
    /// row stays at the lifting line with zero age. This preserves the shared
    /// interface endpoints and elapsed-time schedule without cloning the wake.
    /// Call this instead of separate near shedding and far diffusion; aggregate
    /// additional rows only afterward, and rebuild any image snapshots.
    /// This is uniform prescribed convection, not free-wake self-induction.
    ///
    /// # Errors
    /// Invalid shedding interval, any core/translation refusal, or any near
    /// shedding refusal. Both input states are unchanged on every refusal.
    pub fn shed_prescribed(
        &mut self,
        near: &mut FilamentWake,
        gamma: &[f64],
        displacement_m: [f64; 3],
    ) -> Result<(), Refusal> {
        if !self.shed_dt_s.is_finite() || self.shed_dt_s <= 0.0 {
            return Err(Refusal {
                code: "farfield-shed-invalid",
                message: "prescribed shedding needs a finite positive interval".into(),
                ranked_repairs: vec!["retain the interval used for near/far aggregation".into()],
            });
        }
        let growth = self.admit_core_advance(self.shed_dt_s)?;
        if near.rows.iter().any(|row| {
            !row.core2_m2.is_finite()
                || row.core2_m2 < 0.0
                || (growth != 0.0 && !(row.core2_m2 + growth * (row.age_ticks + 1.0)).is_finite())
        }) {
            return Err(Refusal {
                code: "core-advance-overflow",
                message: "the next shedding age would give an invalid near-source core".into(),
                ranked_repairs: vec!["use finite cores and representable viscosity × age".into()],
            });
        }
        self.admit_translation(displacement_m)?;
        // Every far refusal precedes the near transaction. Once it succeeds,
        // the remaining far updates are checked arithmetic with no allocation.
        near.shed(gamma, displacement_m)?;
        self.apply_translation(displacement_m);
        self.apply_core_growth(growth);
        Ok(())
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

    /// Content digest: mode, moments, and every retained source and its core.
    #[must_use]
    pub fn digest(&self) -> String {
        let mut b = self.mode.id_bytes();
        b.extend_from_slice(&self.shed_dt_s.to_bits().to_le_bytes());
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
            for value in c.core_moment {
                b.extend_from_slice(&value.to_bits().to_le_bytes());
            }
            b.extend_from_slice(&c.max_abs_core_delta_m2.to_bits().to_le_bytes());
            b.extend_from_slice(&c.diffusion_offset_m2.to_bits().to_le_bytes());
            b.extend_from_slice(&c.max_source_core2_m2.to_bits().to_le_bytes());
            assert_eq!(
                c.segs.len(),
                c.source_core2_m2.len(),
                "one core per source segment"
            );
            b.extend_from_slice(&(c.segs.len() as u64).to_le_bytes());
            for ((a, end, gamma), core2) in c.segs.iter().zip(&c.source_core2_m2) {
                for value in a.iter().chain(end.iter()) {
                    b.extend_from_slice(&value.to_bits().to_le_bytes());
                }
                b.extend_from_slice(&gamma.to_bits().to_le_bytes());
                b.extend_from_slice(&core2.to_bits().to_le_bytes());
            }
        }
        hash_domain("org.frankensim.fs-vpm.farfield3d.v3", &b).to_hex()
    }
}

/// Hybrid induced velocity: exact near wake + far multipoles.
#[must_use]
pub fn hybrid_velocity(near: &FilamentWake, far: &FarField, p: [f64; 3]) -> [f64; 3] {
    let nv = near.induced_velocity_with_core_growth(p, far.mode.core_growth_m2(far.shed_dt_s));
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

/// Build one multipole cell from exact segments (shared by
/// aggregation AND image construction so the two paths are
/// arithmetically IDENTICAL — the conversion-symmetry battery relies
/// on this).
#[must_use]
pub fn cell_from_segs(segs: Vec<([f64; 3], [f64; 3], f64)>, core2_m2: f64) -> FarCell {
    let cores = vec![core2_m2; segs.len()];
    cell_from_segs_with_cores(segs, cores)
}

/// Build a cell while retaining the finite nonnegative core of every source.
/// `segs` and `source_core2_m2` must have the same length. The infallible
/// constructor has the existing trusted-data semantics; checked aggregation
/// additionally rejects unrepresentable moments before mutating its wake.
///
/// # Panics
/// Unequal vector lengths, or non-finite/negative source core widths.
#[must_use]
pub fn cell_from_segs_with_cores(
    segs: Vec<([f64; 3], [f64; 3], f64)>,
    source_core2_m2: Vec<f64>,
) -> FarCell {
    assert_eq!(
        segs.len(),
        source_core2_m2.len(),
        "one core per source segment"
    );
    assert!(
        source_core2_m2.iter().all(|c| c.is_finite() && *c >= 0.0),
        "finite nonnegative source cores"
    );
    let weight_scale = segs
        .iter()
        .map(|(a, b, g)| g.abs() * dist(*a, *b))
        .fold(0.0_f64, f64::max);
    let core_scale = source_core2_m2.iter().copied().fold(0.0_f64, f64::max);
    let minimum_core = source_core2_m2.iter().copied().fold(core_scale, f64::min);
    let normalized_weight = |a: [f64; 3], b: [f64; 3], g: f64| {
        if weight_scale == 0.0 {
            0.0
        } else {
            (g.abs() * dist(a, b)) / weight_scale
        }
    };
    let weight_sum: f64 = segs
        .iter()
        .map(|(a, b, g)| normalized_weight(*a, *b, *g))
        .sum();
    let mut centroid = [0.0; 3];
    let mut scaled_core_sum = 0.0;
    if weight_sum > 0.0 {
        for ((a, b, g), core2) in segs.iter().zip(&source_core2_m2) {
            let weight = normalized_weight(*a, *b, *g) / weight_sum;
            for k in 0..3 {
                centroid[k] += weight * (0.5 * a[k] + 0.5 * b[k]);
            }
            if core_scale > 0.0 {
                scaled_core_sum += weight * (*core2 / core_scale);
            }
        }
    }
    let core2_m2 = if minimum_core == core_scale {
        // A common source core remains exact, without a roundoff-only spread.
        core_scale
    } else if weight_sum == 0.0 {
        source_core2_m2.first().copied().unwrap_or(0.0)
    } else {
        // The exact normalized sum is in [0,1]; cap roundoff at its endpoint.
        core_scale * scaled_core_sum.min(1.0)
    };
    let mut w_sum = [0.0; 3];
    let mut w_abs_sum = 0.0;
    let mut dip = [[0.0; 3]; 3];
    let mut quad = [[[0.0; 3]; 3]; 3];
    let mut core_moment = [0.0; 3];
    let mut max_abs_core_delta_m2 = 0.0_f64;
    let mut radius = 0.0_f64;
    for ((a, b, g), source_core2) in segs.iter().zip(&source_core2_m2) {
        let dl = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [g * dl[0], g * dl[1], g * dl[2]];
        let d = [
            0.5 * a[0] + 0.5 * b[0] - centroid[0],
            0.5 * a[1] + 0.5 * b[1] - centroid[1],
            0.5 * a[2] + 0.5 * b[2] - centroid[2],
        ];
        let delta_core = *source_core2 - core2_m2;
        max_abs_core_delta_m2 = max_abs_core_delta_m2.max(delta_core.abs());
        w_abs_sum += g.abs() * dist(*a, *b);
        for j in 0..3 {
            w_sum[j] += w[j];
            core_moment[j] += w[j] * delta_core;
            for l in 0..3 {
                dip[j][l] += w[j] * d[l];
                for m in 0..3 {
                    // Exact line moment: mean(s²)=1/12 for s∈[-1/2,1/2].
                    quad[j][l][m] += w[j] * (d[l] * d[m] + dl[l] * dl[m] / 12.0);
                }
            }
        }
        radius = radius.max(dist(*a, centroid)).max(dist(*b, centroid));
    }
    FarCell {
        centroid,
        radius_m: radius,
        w_sum,
        w_abs_sum,
        dip,
        quad,
        core2_m2,
        source_core2_m2,
        core_moment,
        max_abs_core_delta_m2,
        diffusion_offset_m2: 0.0,
        max_source_core2_m2: core_scale,
        segs,
    }
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    d[0].hypot(d[1]).hypot(d[2])
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
