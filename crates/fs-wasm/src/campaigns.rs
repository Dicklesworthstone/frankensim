//! fs-wasm · CAMPAIGN tier (Tier IV) — ten evidence-bearing end-to-end campaign
//! pipelines in the browser.
//!
//! Each function here runs a real FrankenSim campaign targeted at
//! `wasm32-unknown-unknown`. There are no mock results; where browser-only
//! constraints require a checked transcription, that boundary is called out.
//! The
//! sum-of-squares global-optimality proofs, the numerical homogenization, the
//! Lyapunov/spectral flutter certificates, the tropical critical-path
//! scheduler, the PDHG truss layout LP, the Kalman/VoI sensor planner, the
//! interval-certified neural SDF, the MAP-Elites grammar illuminator, the
//! anytime-valid Bayesian-optimization stop, and the lattice-Boltzmann
//! credibility map are composed from the same pure-Rust numerical leaves.
//!
//! SAFETY CONTRACT (identical to the rest of the crate): `unsafe_code` is
//! forbidden, every input is clamped to a safe range, every fallible kernel
//! result is folded to `NaN` / an empty vector, and every documented panic
//! precondition of the composed crates is respected (`fs_eproc::
//! BettingEProcess::new` needs `0 < null_mean < 1`; `::observe` needs
//! `x ∈ [0, 1]`; `fs_archive::MapElites::add` needs `fitness ≥ 0` and finite;
//! `fs_lbm::plan_scaling` needs positive Reynolds and length). Nothing here
//! can trap — a wasm trap would kill the whole page.
//!
//! Where a campaign's `run_campaign` returns only final scalars but the viz
//! needs a per-iteration trajectory / geometry, the campaign's `run_campaign`
//! body is transcribed here using the same public APIs (and the campaign
//! crate's public helpers) so results stay aligned with the native conformance
//! tests. Those campaigns are: AnytimeBO, FlowCert, GrammarForge,
//! and TrussPath. CampaignSchedule now exposes its visualization fields directly
//! and is invoked through its canonical admitted API.

use fs_archive::MapElites;
use fs_bo::{Gp, Kernel, Matern, expected_improvement};
use fs_eproc::{BettingEProcess, GaussianMixtureCs};
use fs_evidence::{Color, ColorRank};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fab::min_feature_size;
use fs_lbm::{Lbm, plan_scaling, poiseuille_analytic};
use fs_rep_neural::{Layer, MlpSdf};
use fs_schedule_e2e::{ScheduleDisposition, Study};
use fs_shapeprog::{max_sdf_discrepancy, simplify};
use fs_sparse::{Coo, Csr};
use fs_truss::{LayoutCertificateLimits, LayoutCertificateProblem, PdhgSettings};
use fs_truss_e2e::{
    analyze_load_path, estimated_optimality_color, optimality_color_from_certificate,
};
use fs_viz::Grid2;
use fs_voi::{Action, ActionKind, DesignEstimate, Uncertainty};

/* ======================================================================= */
/*  Small shared helpers                                                    */
/* ======================================================================= */

/// Map an epistemic [`ColorRank`] to the wire code `2 = Verified`,
/// `1 = Validated`, `0 = Estimated`.
pub(crate) fn rank_code(r: ColorRank) -> f64 {
    match r {
        ColorRank::Verified => 2.0,
        ColorRank::Validated => 1.0,
        ColorRank::Estimated => 0.0,
    }
}

/// `(verified_flag, lo, hi)` for a [`Color`] — the `Verified` interval, or
/// `(0, NaN, NaN)` for anything weaker.
pub(crate) fn verified_bounds(c: &Color) -> (f64, f64, f64) {
    match c {
        Color::Verified { lo, hi } => (1.0, *lo, *hi),
        _ => (0.0, f64::NAN, f64::NAN),
    }
}

/// Construct the bounded deterministic context used only by cold certificate
/// work in browser campaigns. No task, thread, or context escapes this scope.
pub(crate) fn with_certificate_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x7A55_5741_534D_0001,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// A finite value passes through; `±∞` / `NaN` fold to `NaN` (log/plot-safe).
fn fon(x: f64) -> f64 {
    if x.is_finite() { x } else { f64::NAN }
}

/// `n` evenly spaced points on `[lo, hi]` inclusive (`n == 1` ⇒ just `lo`).
fn linspace(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![lo];
    }
    (0..n)
        .map(|i| lo + (hi - lo) * i as f64 / (n - 1) as f64)
        .collect()
}

/* ======================================================================= */
/*  1 · ProofRobust — fs-robustopt-e2e (fs-sos × fs-robust × fs-evidence)   */
/* ======================================================================= */

/// **ProofRobust**: three design families whose nominal global optima are each
/// PROVEN by an executable sum-of-squares certificate ([`fs_sos`]), then ranked
/// by worst-case CVaR under a `±sigma` manufacturing-tolerance grid
/// ([`fs_robust`]). The lowest-nominal family is not the robust winner — a
/// flatter family wins under CVaR — and both facts carry honest colors
/// ([`fs_evidence`]): the nominal optimum is `Verified`, the CVaR is
/// `Estimated`. Runs the real `fs_robustopt_e2e::run_campaign` on the fixed
/// `demo_families()` (curvatures 1.0 / 0.5 / 2.0, all with `x* = 2`).
///
/// `alpha` — CVaR confidence (clamped `0.5..=0.999`); `sigma` — tolerance
/// half-width (clamped `0..=5`); `n` — tolerance grid points (clamped
/// `3..=201`). Default `(0.9, 2.0, 41)`.
///
/// Output layout (length `6 + 6·F`, `F` = family count = 3):
/// - `[0]` — `F`.
/// - `[1]` — `certified_count` (families with an SOS-proven optimum).
/// - `[2]` — `reorders` (1 if the robust winner ≠ the nominal winner).
/// - `[3]` — `headline_rank` (2 Verified / 1 Validated / 0 Estimated — the
///   robust winner's weakest input color; `Estimated` for a CVaR).
/// - `[4]` — `nominal_winner_idx` (family index, lowest nominal cost).
/// - `[5]` — `robust_winner_idx` (family index, lowest robust CVaR cost).
/// - then `F` blocks of 6: `[x_star, nominal_cost, robust_cost, verified_flag,
///   cert_lo, cert_hi]` — the proven optimizer, its proven nominal cost, its
///   worst-case CVaR cost, whether the SOS proof checked, and the certified
///   `[lo, hi]` enclosure of the nominal optimum (`NaN` if unverified).
pub fn proofrobust(alpha: f64, sigma: f64, n: usize) -> Vec<f64> {
    let alpha = alpha.clamp(0.5, 0.999);
    let sigma = sigma.clamp(0.0, 5.0);
    let n = n.clamp(3, 201);

    let families = fs_robustopt_e2e::demo_families();
    let report = fs_robustopt_e2e::run_campaign(&families, alpha, sigma, n);
    let idx_of = |name: &str| -> f64 {
        report
            .families
            .iter()
            .position(|f| f.name == name)
            .map_or(-1.0, |i| i as f64)
    };

    let f = report.families.len();
    let mut out = Vec::with_capacity(6 + 6 * f);
    out.push(f as f64);
    out.push(report.certified_count as f64);
    out.push(if report.robustness_reorders { 1.0 } else { 0.0 });
    out.push(rank_code(report.headline_rank));
    out.push(idx_of(&report.nominal_winner));
    out.push(idx_of(&report.robust_winner));
    for v in &report.families {
        let (vf, lo, hi) = verified_bounds(&v.nominal_color);
        out.push(v.x_star);
        out.push(fon(v.nominal_cost));
        out.push(fon(v.robust_cost));
        out.push(vf);
        out.push(lo);
        out.push(hi);
    }
    out
}

/* ======================================================================= */
/*  2 · MetamatCert — fs-metamat-e2e (fs-lattice × fs-sos × fs-evidence)     */
/* ======================================================================= */

/// **MetamatCert**: the certified stiffness-density frontier of a holed-plate
/// metamaterial. Numerical homogenization ([`fs_lattice`]) gives each porosity
/// an effective Voigt tensor and density; every point is PROVEN PSD-stable
/// (`fs_sos::is_psd`) and Voigt-admissible (`fs_lattice::voigt_bound`). An
/// all-stable, all-admissible frontier is `Verified`. Runs the real
/// `fs_metamat_e2e::run_campaign` on `radii = linspace(0, rmax, points)`.
///
/// `n` — unit-cell resolution (clamped `6..=14`); `points` — frontier
/// resolution (clamped `2..=12`); `rmax` — max hole radius (clamped
/// `0.05..=0.45`). Default `(10, 6, 0.40)` (reproduces `default_radii()`).
///
/// Output layout (length `7 + 6·P`, `P` = frontier points):
/// - `[0]` — `P`.
/// - `[1]` — `c_solid` (solid-cell axial stiffness `C₁₁ˢᵒˡⁱᵈ`).
/// - `[2]` — `all_stable` (0/1).
/// - `[3]` — `all_admissible` (0/1).
/// - `[4]` — `stiffness_monotone` (0/1 — `C₁₁` decreasing in porosity).
/// - `[5]` — `solid_specific_optimal` (0/1 — Voigt proves the solid maximizes
///   specific stiffness).
/// - `[6]` — `color_verified` (0/1 — the frontier's stability color).
/// - then `P` blocks of 6: `[r, density, c11, specific_stiffness, stable,
///   admissible]`.
pub fn metamatcert(n: usize, points: usize, rmax: f64) -> Vec<f64> {
    let n = n.clamp(6, 14);
    let points = points.clamp(2, 12);
    let rmax = rmax.clamp(0.05, 0.45);
    let radii = linspace(0.0, rmax, points);

    let report = fs_metamat_e2e::run_campaign(n, &radii);
    let p = report.frontier.len();
    let mut out = Vec::with_capacity(7 + 6 * p);
    out.push(p as f64);
    out.push(fon(report.c_solid));
    out.push(if report.all_stable { 1.0 } else { 0.0 });
    out.push(if report.all_admissible { 1.0 } else { 0.0 });
    out.push(if report.stiffness_monotone { 1.0 } else { 0.0 });
    out.push(if report.solid_is_specific_optimal {
        1.0
    } else {
        0.0
    });
    out.push(
        if matches!(report.stability_color, Color::Verified { .. }) {
            1.0
        } else {
            0.0
        },
    );
    for c in &report.frontier {
        out.push(c.r);
        out.push(fon(c.density));
        out.push(fon(c.c11));
        out.push(fon(c.specific_stiffness));
        out.push(if c.stable { 1.0 } else { 0.0 });
        out.push(if c.admissible { 1.0 } else { 0.0 });
    }
    out
}

/* ======================================================================= */
/*  3 · FlutterCert — fs-flutter-e2e (fs-sos × fs-spectral × fs-couple)      */
/* ======================================================================= */

/// **FlutterCert**: the added-mass flutter boundary `μ* = 2` PROVEN by a
/// Lyapunov certificate (`fs_sos::lyapunov_certifies_stability`) and
/// independently cross-checked by the symmetric-part spectral abscissa
/// (`fs_spectral`); the two boundaries agreeing certifies the certifier. A
/// partitioned coupled solve ([`fs_couple`]) shows naive staggering diverges
/// early while Aitken relaxation reaches the boundary. Runs the real
/// `fs_flutter_e2e::run_campaign`.
///
/// `lo`/`hi` — μ sweep bounds (clamped `lo ∈ 0.05..=2.9`,
/// `hi ∈ lo+0.1..=3.0`); `steps` — samples (clamped `2..=200`). Default
/// `(0.55, 2.45, 20)`.
///
/// Output layout (length `9 + 5·S`, `S` = sample count):
/// - `[0]` — `S`.
/// - `[1]` — `lyapunov_boundary` (largest proven-stable μ; `NaN` if none).
/// - `[2]` — `spectral_boundary` (largest μ the abscissa calls stable).
/// - `[3]` — `boundaries_agree` (0/1).
/// - `[4]` — `naive_boundary` (largest μ the naive solve converged).
/// - `[5]` — `aitken_boundary` (largest μ Aitken converged).
/// - `[6]` — `aitken_beats_naive` (0/1).
/// - `[7]` — `witness_mu` (`NaN` if none) — a μ where the proof holds, naive
///   fails, Aitken succeeds.
/// - `[8]` — `witness_verified` (0/1).
/// - then `S` blocks of 5: `[mu, lyapunov_stable, abscissa, naive_converged,
///   aitken_converged]`.
pub fn fluttercert(lo: f64, hi: f64, steps: usize) -> Vec<f64> {
    let lo = lo.clamp(0.05, 2.9);
    let hi = hi.clamp(lo + 0.1, 3.0);
    let steps = steps.clamp(2, 200);

    let report = fs_flutter_e2e::run_campaign(lo, hi, steps);
    let s = report.samples.len();
    let mut out = Vec::with_capacity(9 + 5 * s);
    out.push(s as f64);
    out.push(fon(report.lyapunov_boundary));
    out.push(fon(report.eigen_boundary));
    out.push(if report.boundaries_agree { 1.0 } else { 0.0 });
    out.push(fon(report.naive_boundary));
    out.push(fon(report.aitken_boundary));
    out.push(if report.aitken_beats_naive { 1.0 } else { 0.0 });
    out.push(report.witness_mu.map_or(f64::NAN, fon));
    out.push(
        if matches!(report.witness_color, Some(Color::Verified { .. })) {
            1.0
        } else {
            0.0
        },
    );
    for sm in &report.samples {
        out.push(sm.mu);
        out.push(if sm.lyapunov_stable { 1.0 } else { 0.0 });
        out.push(fon(sm.spectral_abscissa));
        out.push(if sm.naive_converged { 1.0 } else { 0.0 });
        out.push(if sm.aitken_converged { 1.0 } else { 0.0 });
    }
    out
}

/* ======================================================================= */
/*  4 · CampaignSchedule — fs-schedule-e2e (fs-tropical × fs-voi)            */
/* ======================================================================= */

/// **CampaignSchedule**: the outward-bounded makespan of a design campaign as a
/// tropical (max-plus) critical path ([`fs_tropical`]) — `Verified` by an
/// enclosure —
/// plus an EVPI-driven Act/Stop recommendation over the candidate designs
/// ([`fs_voi`]). The canonical [`fs_schedule_e2e::run_campaign`] admission path
/// supplies the per-study slack vector and typed disposition for the viz.
///
/// Fixed scenario (studies `surrogate-B(2)`, `hifi-B(8,[0])`,
/// `sample-scenarios(4)`, `windtunnel-A(latency)`, `decide(1,[1,2,3])`;
/// designs A/B/C; actions hifi-B/sample-B/windtunnel-A), parametrized by the
/// long-pole latency, the contender B's mean cost, and the stop threshold.
///
/// `windtunnel_latency` — clamped `5..=15` (default 12); `design_b_mean` —
/// clamped `0.60..=1.10` (default 0.65); `stop_threshold` — clamped
/// `0..=1e3` (default 1e-6).
///
/// Output layout (length `12 + 2·N + P`, `N` = studies = 5, `P` = critical
/// path length):
/// - `[0]` — `makespan`.
/// - `[1]`,`[2]` — outward-rounded `makespan_lo`, `makespan_hi` (the
///   `Verified` interval containing `makespan`).
/// - `[3]` — `N` (study count).
/// - `[4]` — `P` (critical-path length).
/// - `[5]` — `bottleneck_idx` (study index; `-1` if none).
/// - `[6]` — `evpi`.
/// - `[7]` — `flip_risk` (top-two ranking-flip probability).
/// - `[8]` — `should_stop` (0/1).
/// - `[9]` — `leading_design_idx` (0=A,1=B,2=C — lowest cost).
/// - `[10]` — `rec_code` (0 = Act, 1 = robust Stop, 2 = expand menu).
/// - `[11]` — `value_per_cost` of the recommended action (`NaN` if Stop).
/// - then `N` study latencies; then `N` study slacks; then `P` critical-path
///   study indices (source → sink).
pub fn schedule_campaign(
    windtunnel_latency: f64,
    design_b_mean: f64,
    stop_threshold: f64,
) -> Vec<f64> {
    let wt = windtunnel_latency.clamp(5.0, 15.0);
    let b_mean = design_b_mean.clamp(0.60, 1.10);
    let thr = stop_threshold.clamp(0.0, 1.0e3);

    let studies = vec![
        Study::new("surrogate-B", 2.0, vec![]),
        Study::new("hifi-B", 8.0, vec![0]),
        Study::new("sample-scenarios", 4.0, vec![]),
        Study::new("windtunnel-A", wt, vec![]),
        Study::new("decide", 1.0, vec![1, 2, 3]),
    ];
    let latencies: Vec<f64> = studies.iter().map(|study| study.latency).collect();
    let n = studies.len();

    // Candidate designs (fs-voi MINIMIZES; lower cost is better).
    let designs = vec![
        DesignEstimate::new(
            "A",
            0.60,
            Uncertainty {
                numerical: 0.05,
                statistical: 0.05,
                model: 0.08,
            },
        ),
        DesignEstimate::new(
            "B",
            b_mean,
            Uncertainty {
                numerical: 0.08,
                statistical: 0.06,
                model: 0.10,
            },
        ),
        DesignEstimate::new(
            "C",
            0.90,
            Uncertainty {
                numerical: 0.05,
                statistical: 0.05,
                model: 0.05,
            },
        ),
    ];
    let actions = vec![
        Action {
            name: "hifi-B".into(),
            kind: ActionKind::Simulate,
            target_design: "B".into(),
            reduction: 0.9,
            cost: 8.0,
        },
        Action {
            name: "sample-B".into(),
            kind: ActionKind::Sample,
            target_design: "B".into(),
            reduction: 0.7,
            cost: 4.0,
        },
        Action {
            name: "windtunnel-A".into(),
            kind: ActionKind::Test,
            target_design: "A".into(),
            reduction: 0.8,
            cost: 12.0,
        },
    ];

    let Ok(report) = fs_schedule_e2e::run_campaign(&studies, &designs, &actions, thr) else {
        let mut out = vec![f64::NAN; 12];
        out[3] = n as f64;
        out[4] = 0.0;
        out.extend_from_slice(&latencies);
        out.extend(std::iter::repeat_n(f64::NAN, n));
        return out;
    };
    let (makespan_lo, makespan_hi) = match &report.makespan_color {
        Color::Verified { lo, hi } => (*lo, *hi),
        _ => (f64::NAN, f64::NAN),
    };
    let leading_idx = designs
        .iter()
        .position(|design| design.name == report.leading_design)
        .map_or(-1.0, |index| index as f64);
    let rec_code = match report.disposition {
        ScheduleDisposition::Act => 0.0,
        ScheduleDisposition::RobustStop => 1.0,
        ScheduleDisposition::NoEffectiveAction => 2.0,
    };

    let p = report.critical_path.len();
    let mut out = Vec::with_capacity(12 + 2 * n + p);
    out.push(report.makespan);
    out.push(makespan_lo);
    out.push(makespan_hi);
    out.push(n as f64);
    out.push(p as f64);
    out.push(report.bottleneck_index.map_or(-1.0, |index| index as f64));
    out.push(fon(report.evpi));
    out.push(fon(report.flip_risk));
    out.push(if report.should_stop { 1.0 } else { 0.0 });
    out.push(leading_idx);
    out.push(rec_code);
    out.push(fon(report
        .recommendation_value_per_cost
        .unwrap_or(f64::NAN)));
    out.extend_from_slice(&latencies);
    for &s in &report.slack {
        out.push(fon(s));
    }
    for &i in &report.critical_path {
        out.push(i as f64);
    }
    out
}

/* ======================================================================= */
/*  5 · TrussPath — fs-truss-e2e (fs-truss × fs-tropical)                    */
/* ======================================================================= */

/// A lean, **fnx-free** ground structure — exactly the fields the truss LP and
/// the critical-path logic actually read (`nodes`, `members`, `lengths`).
///
/// WHY THIS EXISTS: `fs_truss::GroundStructure::try_grid` additionally builds
/// an `fnx_classes::Graph`, whose internal compatibility-evidence path reads a
/// platform time source during construction (fine natively) that *compiles but TRAPS at
/// runtime* on `wasm32-unknown-unknown` ("time not implemented on this
/// platform") — killing the page. The truss LP never reads that graph, so the
/// grid enumeration, the numerical core of LP `try_assemble`, and the PDHG
/// `solve`/`diagnostics` are transcribed here so the wasm runtime path never
/// touches fnx. It shares the
/// same bounded tropical analysis, but cross-implementation bit identity is not
/// claimed until a retained native/WASM golden verifies it.
struct LeanGround {
    nodes: Vec<[f64; 2]>,
    members: Vec<(usize, usize)>,
    lengths: Vec<f64>,
}

/// Replicates the candidate enumeration in `GroundStructure::try_grid`
/// (ground.rs) WITHOUT the fnx
/// `Graph`: the node grid, the length-bound filter, and the collinear-through-
/// node skip. (Our rules carry an empty `angles` set, so the direction filter
/// is a no-op and is omitted.)
fn truss_grid(nx: usize, ny: usize, w: f64, h: f64, min_len: f64, max_len: f64) -> LeanGround {
    let mut nodes = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            nodes.push([
                w * i as f64 / (nx - 1) as f64,
                h * j as f64 / (ny - 1) as f64,
            ]);
        }
    }
    let n = nodes.len();
    let mut members = Vec::new();
    let mut lengths = Vec::new();
    for a in 0..n {
        for b in (a + 1)..n {
            let dx = nodes[b][0] - nodes[a][0];
            let dy = nodes[b][1] - nodes[a][1];
            let len = dx.hypot(dy);
            if len < min_len || len > max_len {
                continue;
            }
            // Skip members that pass exactly through another node.
            let mut through = false;
            for (c, node) in nodes.iter().enumerate() {
                if c == a || c == b {
                    continue;
                }
                let cx = node[0] - nodes[a][0];
                let cy = node[1] - nodes[a][1];
                let cross = cx * dy - cy * dx;
                let dot = cx * dx + cy * dy;
                if cross.abs() < 1e-9 * len && dot > 1e-12 && dot < len * len - 1e-12 {
                    through = true;
                    break;
                }
            }
            if through {
                continue;
            }
            members.push((a, b));
            lengths.push(len);
        }
    }
    LeanGround {
        nodes,
        members,
        lengths,
    }
}

/// The assembled layout LP (fs-sparse only) — a transcription of
/// `fs_truss::LayoutLp` over the lean ground structure.
struct LeanLp {
    a: Csr,
    at: Csr,
    c: Vec<f64>,
    b: Vec<f64>,
    norm_est: f64,
}

/// The PDHG solve diagnostics.
struct LeanReport {
    iters: usize,
    volume: f64,
    gap: f64,
    eq_residual: f64,
}

impl LeanLp {
    /// Numerical transcription of `fs_truss::LayoutLp::try_assemble` (lp.rs),
    /// without the native admitted-value and `Cx` boundary.
    fn assemble(
        gs: &LeanGround,
        supported: &dyn Fn(usize, usize) -> bool,
        loads: &dyn Fn(usize) -> [f64; 2],
        sigma_y: f64,
    ) -> LeanLp {
        let n = gs.nodes.len();
        let mut dof_map: Vec<Option<usize>> = Vec::with_capacity(2 * n);
        let mut nf = 0usize;
        for node in 0..n {
            for comp in 0..2 {
                if supported(node, comp) {
                    dof_map.push(None);
                } else {
                    dof_map.push(Some(nf));
                    nf += 1;
                }
            }
        }
        let m = gs.members.len();
        let mut coo = Coo::new(nf, 2 * m);
        for (k, &(a, b)) in gs.members.iter().enumerate() {
            let dx = (gs.nodes[b][0] - gs.nodes[a][0]) / gs.lengths[k];
            let dy = (gs.nodes[b][1] - gs.nodes[a][1]) / gs.lengths[k];
            let entries = [(2 * a, dx), (2 * a + 1, dy), (2 * b, -dx), (2 * b + 1, -dy)];
            for (dof, v) in entries {
                if let Some(row) = dof_map[dof] {
                    coo.push(row, k, v);
                    coo.push(row, m + k, -v);
                }
            }
        }
        let a_mat = coo.assemble();
        let at = fs_sparse::ops::transpose(&a_mat);
        let mut b_vec = vec![0.0f64; nf];
        for node in 0..n {
            let f = loads(node);
            for comp in 0..2 {
                if let Some(row) = dof_map[2 * node + comp] {
                    b_vec[row] = f[comp];
                }
            }
        }
        let mut c = Vec::with_capacity(2 * m);
        for &l in &gs.lengths {
            c.push(l / sigma_y);
        }
        for &l in &gs.lengths {
            c.push(l / sigma_y);
        }
        // Power iteration for ‖A‖ (deterministic start).
        let mut v: Vec<f64> = (0..2 * m).map(|i| 1.0 + ((i % 7) as f64) * 0.1).collect();
        let mut norm_est = 1.0;
        let mut av = vec![0.0f64; nf];
        for _ in 0..30 {
            a_mat.spmv(&v, &mut av);
            let mut atv = vec![0.0f64; 2 * m];
            at.spmv(&av, &mut atv);
            let nrm = atv.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
            norm_est = nrm.sqrt();
            for (vi, ai) in v.iter_mut().zip(&atv) {
                *vi = ai / nrm;
            }
        }
        LeanLp {
            a: a_mat,
            at,
            c,
            b: b_vec,
            norm_est,
        }
    }

    /// Faithful transcription of `fs_truss::LayoutLp::diagnostics`.
    fn diagnostics(&self, x: &[f64], y: &[f64], bnorm: f64) -> (f64, f64, f64) {
        let primal: f64 = self.c.iter().zip(x).map(|(c, x)| c * x).sum();
        let mut aty = vec![0.0f64; self.c.len()];
        self.at.spmv(y, &mut aty);
        let mut scale = 1.0f64;
        for (a, c) in aty.iter().zip(&self.c) {
            if *a < -c && *a < 0.0 {
                scale = scale.min(-c / a);
            }
        }
        let dual: f64 = -(y.iter().zip(&self.b).map(|(y, b)| y * b).sum::<f64>()) * scale.max(0.0);
        let mut ax = vec![0.0f64; self.b.len()];
        self.a.spmv(x, &mut ax);
        let eq_res = ax
            .iter()
            .zip(&self.b)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
            / bnorm;
        let gap = (primal - dual).abs() / primal.abs().max(1e-30);
        (gap, eq_res, primal)
    }

    /// Faithful transcription of `fs_truss::LayoutLp::solve` (cold start,
    /// PDHG / Chambolle–Pock; trace dropped — the viz does not draw it).
    fn solve(
        &self,
        max_iters: usize,
        gap_tol: f64,
        check_every: usize,
    ) -> (Vec<f64>, Vec<f64>, LeanReport) {
        let nvar = self.c.len();
        let nrow = self.b.len();
        let mut x = vec![0.0; nvar];
        let mut y = vec![0.0; nrow];
        let step = 0.95 / self.norm_est.max(1e-30);
        let (tau, sigma) = (step, step);
        let bnorm = self.b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);
        let mut report = LeanReport {
            iters: 0,
            volume: 0.0,
            gap: 0.0,
            eq_residual: 0.0,
        };
        let mut aty = vec![0.0f64; nvar];
        let mut ax = vec![0.0f64; nrow];
        let mut x_prev = x.clone();
        let mut xbar = vec![0.0f64; nvar];
        for it in 0..max_iters {
            self.at.spmv(&y, &mut aty);
            x_prev.copy_from_slice(&x);
            for i in 0..nvar {
                x[i] = (x[i] - tau * (self.c[i] + aty[i])).max(0.0);
            }
            for ((extrapolated, xi), previous) in xbar.iter_mut().zip(&x).zip(&x_prev) {
                *extrapolated = 2.0 * xi - previous;
            }
            self.a.spmv(&xbar, &mut ax);
            for r in 0..nrow {
                y[r] += sigma * (ax[r] - self.b[r]);
            }
            if (it + 1) % check_every == 0 || it + 1 == max_iters {
                let (gap, eq_res, primal) = self.diagnostics(&x, &y, bnorm);
                report.iters = it + 1;
                report.volume = primal;
                report.gap = gap;
                report.eq_residual = eq_res;
                if gap < gap_tol && eq_res < gap_tol {
                    break;
                }
            }
        }
        (x, y, report)
    }
}

/// **TrussPath**: a Michell ground-structure truss iterated toward minimum
/// volume and equilibrium by a first-order PDHG LP ([`fs_truss`], emitting a
/// reported objective-separation diagnostic), then an advisory tropical load
/// path through thresholded active bars. The checked endpoint/path analyzer is
/// shared with [`fs_truss_e2e`]; the solver body remains transcribed so the node
/// coordinates and per-member draw fields are available. Cantilever on
/// `[0,4]×[0,2]`, left edge supported, unit downward load at the free
/// bottom-right corner.
///
/// `nx`/`ny` — grid nodes per axis (clamped `2..=5` / `2..=4`); `gap_tol` —
/// PDHG relative-gap tolerance (clamped `1e-8..=1e-1`). Default `(4, 3, 1e-4)`
/// (43 candidate bars → 6 active, gap ≈ 7.8e-5).
///
/// Output layout (a flat array the viz slices):
/// - `[0]` — `M` (candidate members).
/// - `[1]` — `num_active`.
/// - `[2]` — `total_volume`.
/// - `[3]` — `gap` (PDHG relative primal/dual objective separation).
/// - `[4]` — `eq_residual`.
/// - `[5]` — `iters`.
/// - `[6]` — `solver_converged` (0/1; not a finite optimum certificate).
/// - `[7]` — `P` (critical-path length).
/// - `[8]` — `critical_path_volume`.
/// - `[9]` — `bottleneck_member_idx` (original member index; `-1` if none).
/// - `[10]` — `Nn` (node count).
/// - `[11]` — `load_node_idx`.
/// - then `2·Nn` node coordinates (`x, y` interleaved).
/// - then `M` blocks of 5: `[node_a, node_b, force, volume, is_active]`.
/// - then `P` critical-path member indices (original bar indices, load →
///   support).
/// - final 4 fields: `[optimality_rank, verified_flag, optimum_lo,
///   optimum_hi]`, where rank `2` and finite endpoints can come only from the
///   shared native/browser certificate-promotion gate.
pub fn trusspath(nx: usize, ny: usize, gap_tol: f64) -> Vec<f64> {
    let nx = nx.clamp(2, 5);
    let ny = ny.clamp(2, 4);
    let gap_tol = if gap_tol.is_finite() {
        gap_tol.clamp(1e-8, 1e-1)
    } else {
        1e-4
    };
    let (w, h) = (4.0f64, 2.0f64);

    // Fabrication rules (empty angle set): min/max member length.
    let min_len = 0.1;
    let max_len = (w * w + h * h).sqrt() / 1.5;
    let gs = truss_grid(nx, ny, w, h, min_len, max_len);
    let m = gs.members.len();
    let nn = gs.nodes.len();

    let support_nodes: Vec<usize> = (0..ny).map(|row| row * nx).collect();
    let supported = |node: usize, _comp: usize| node % nx == 0;
    let load_node = nx - 1;
    let loads = |node: usize| {
        if node == load_node {
            [0.0, -1.0]
        } else {
            [0.0, 0.0]
        }
    };

    // Degenerate grid: no members ⇒ emit a minimal, non-trapping header.
    if m == 0 {
        let mut out = vec![0.0, 0.0, 0.0, f64::NAN, f64::NAN, 0.0, 0.0, 0.0, 0.0, -1.0];
        out.push(nn as f64);
        out.push(load_node as f64);
        for p in &gs.nodes {
            out.push(p[0]);
            out.push(p[1]);
        }
        out.extend_from_slice(&[0.0, 0.0, f64::NAN, f64::NAN]);
        return out;
    }

    let lp = LeanLp::assemble(&gs, &supported, &loads, 1.0);
    let settings = PdhgSettings {
        max_iters: 60_000,
        gap_tol,
        check_every: 500,
    };
    let (x, y, report) = lp.solve(settings.max_iters, settings.gap_tol, settings.check_every);
    let optimality_color = match LayoutCertificateProblem::try_new(&lp.a, &lp.c, &lp.b) {
        Ok(problem) => with_certificate_cx(|cx| {
            problem
                .certify_optimum(&x, &y, settings, LayoutCertificateLimits::default(), cx)
                .map_or_else(
                    |_| estimated_optimality_color(report.gap, report.eq_residual),
                    |status| {
                        optimality_color_from_certificate(
                            &problem,
                            &x,
                            &y,
                            settings,
                            &status,
                            report.gap,
                            report.eq_residual,
                            cx,
                        )
                        .unwrap_or_else(|_| {
                            estimated_optimality_color(report.gap, report.eq_residual)
                        })
                    },
                )
        }),
        Err(_) => estimated_optimality_color(report.gap, report.eq_residual),
    };

    let force = |k: usize| x[k] - x[m + k];
    let volume = |k: usize| lp.c[k] * x[k] + lp.c[m + k] * x[m + k];
    let max_force = (0..m).map(|k| force(k).abs()).fold(0.0, f64::max);
    let active_tol = 1e-3 * max_force.max(1e-12);

    let active: Vec<usize> = (0..m).filter(|&k| force(k).abs() > active_tol).collect();
    let num_active = active.len();

    let volumes: Vec<f64> = (0..m).map(volume).collect();
    let (critical_path, critical_path_volume, bottleneck_member) = match analyze_load_path(
        &gs.nodes,
        &gs.members,
        &active,
        &volumes,
        load_node,
        &support_nodes,
    ) {
        Ok(path) => (path.members, path.weight, path.bottleneck_member),
        Err(_) => (Vec::new(), f64::NAN, None),
    };

    let solver_converged = report.gap.is_finite()
        && report.eq_residual.is_finite()
        && report.gap >= 0.0
        && report.eq_residual >= 0.0
        && report.gap < gap_tol
        && report.eq_residual < gap_tol;

    let p = critical_path.len();
    let mut out = Vec::with_capacity(12 + 2 * nn + 5 * m + p + 4);
    out.push(m as f64);
    out.push(num_active as f64);
    out.push(fon(report.volume));
    out.push(fon(report.gap));
    out.push(fon(report.eq_residual));
    out.push(report.iters as f64);
    out.push(if solver_converged { 1.0 } else { 0.0 });
    out.push(p as f64);
    out.push(fon(critical_path_volume));
    out.push(bottleneck_member.map_or(-1.0, |i| i as f64));
    out.push(nn as f64);
    out.push(load_node as f64);
    for pnode in &gs.nodes {
        out.push(pnode[0]);
        out.push(pnode[1]);
    }
    let active_set: std::collections::BTreeSet<usize> = active.iter().copied().collect();
    for k in 0..m {
        let (a, b) = gs.members[k];
        out.push(a as f64);
        out.push(b as f64);
        out.push(fon(force(k)));
        out.push(fon(volume(k)));
        out.push(if active_set.contains(&k) { 1.0 } else { 0.0 });
    }
    for &k in &critical_path {
        out.push(k as f64);
    }
    let (verified, optimum_lo, optimum_hi) = verified_bounds(&optimality_color);
    out.push(rank_code(optimality_color.rank()));
    out.push(verified);
    out.push(optimum_lo);
    out.push(optimum_hi);
    out
}

/* ======================================================================= */
/*  6 · SensorForge — fs-oed-e2e (fs-assimilate × fs-voi × fs-toleralloc)    */
/* ======================================================================= */

/// **SensorForge**: greedy value-of-information sensor placement. Each
/// candidate design is a Gaussian belief ([`fs_assimilate`]); at every step the
/// EVPI-driven `recommend` ([`fs_voi`]) places the next sensor on the candidate
/// whose measurement most sharpens the DECISION, fuses it with the exact scalar
/// Kalman update, and STOPS the instant the decision is robust. This calls the
/// checked `fs_oed_e2e::run_campaign` directly and serializes its native EVPI
/// trace, posterior summaries, and allocation. Uses the fixed `demo_candidates()`
/// (A/B/C/D), with B's prior mean and truth set to `b_prior_mean`.
///
/// `threshold` — EVPI stop threshold (clamped `1e-6..=1.0`, default 0.01);
/// `max_sensors` — placement cap (clamped `0..=64`, default 12); `b_prior_mean`
/// — B's prior mean & truth (clamped `0.60..=0.90`, default 0.65).
///
/// Output layout (a flat array the viz slices; `C` = candidates = 4,
/// `S` = sensors placed, `T = S + 1`):
/// - `[0]` — `C`.
/// - `[1]` — `S`.
/// - `[2]` — `prior_total_var`.
/// - `[3]` — `posterior_total_var`.
/// - `[4]` — `variance_reduction`.
/// - `[5]` — `initial_evpi`.
/// - `[6]` — `final_evpi`.
/// - `[7]` — `decision_robust` (0/1 — planner chose to STOP).
/// - `[8]` — `chosen_candidate_idx` (lowest-cost posterior).
/// - `[9]` — `T` (EVPI trace length).
/// - then `T` per-step EVPI values (`[0]` = initial, one per placement after).
/// - then `S` placed candidate indices.
/// - then `C` blocks of 2: `[posterior_mean, posterior_var]`.
/// - then `C` allocation tolerances (candidate order; NaN if unconstrained).
pub fn sensorforge(threshold: f64, max_sensors: usize, b_prior_mean: f64) -> Vec<f64> {
    let threshold = if threshold.is_finite() {
        threshold.clamp(1e-6, 1.0)
    } else {
        0.01
    };
    let max_sensors = max_sensors.min(64);
    let b_prior_mean = if b_prior_mean.is_finite() {
        b_prior_mean.clamp(0.60, 0.90)
    } else {
        0.65
    };

    let Ok(mut cands) = fs_oed_e2e::demo_candidates() else {
        return Vec::new();
    };
    let Some(default_b) = cands.get(1) else {
        return Vec::new();
    };
    let Ok(b) = fs_oed_e2e::Candidate::new(
        default_b.name(),
        b_prior_mean,
        b_prior_mean,
        default_b.prior_variance(),
        default_b.sensor_noise(),
        default_b.sensor_cost(),
    ) else {
        return Vec::new();
    };
    cands[1] = b;
    let c = cands.len();
    let Ok(report) = fs_oed_e2e::run_campaign(&cands, threshold, max_sensors) else {
        return Vec::new();
    };
    let Some(placements): Option<Vec<usize>> = report
        .placements
        .iter()
        .map(|name| cands.iter().position(|candidate| candidate.name() == name))
        .collect()
    else {
        return Vec::new();
    };
    let chosen_idx = cands
        .iter()
        .position(|candidate| candidate.name() == report.chosen_design)
        .map_or(-1.0, |index| index as f64);
    if report.posteriors.len() != c {
        return Vec::new();
    }
    let allocation: std::collections::BTreeMap<&str, f64> = report
        .allocation
        .iter()
        .map(|(name, tolerance)| (name.as_str(), *tolerance))
        .collect();

    let s = placements.len();
    let t = report.evpi_trace.len();
    let mut out = Vec::with_capacity(10 + t + s + 2 * c + c);
    out.push(c as f64);
    out.push(s as f64);
    out.push(fon(report.prior_total_variance));
    out.push(fon(report.posterior_total_variance));
    out.push(fon(report.variance_reduction));
    out.push(fon(report.initial_evpi));
    out.push(fon(report.final_evpi));
    out.push(if report.decision_robust { 1.0 } else { 0.0 });
    out.push(chosen_idx);
    out.push(t as f64);
    for &e in &report.evpi_trace {
        out.push(fon(e));
    }
    for &i in &placements {
        out.push(i as f64);
    }
    for posterior in &report.posteriors {
        out.push(fon(posterior.mean));
        out.push(fon(posterior.variance));
    }
    for cand in &cands {
        out.push(fon(allocation
            .get(cand.name())
            .copied()
            .unwrap_or(f64::NAN)));
    }
    out
}

/* ======================================================================= */
/*  7 · NeuroShapeCert — fs-neuroshape-e2e (fs-rep-neural × fs-viz)          */
/* ======================================================================= */

/// Build the tunable blob SDF network. `MlpSdf::new` spectral-normalizes every
/// layer to `√18`, so `L = 18`; the output bias `lift` raises the field
/// (default 6.5 reproduces `blob_sdf_net()`).
fn neuro_net(lift: f64) -> MlpSdf {
    let l1 = Layer::new(
        vec![
            vec![3.0, 0.0],
            vec![-3.0, 0.0],
            vec![0.0, 3.0],
            vec![0.0, -3.0],
        ],
        vec![-2.1, -2.1, -2.1, -2.1],
    );
    let l2 = Layer::new(vec![vec![1.0, 1.0, 1.0, 1.0]], vec![lift]);
    MlpSdf::new(vec![l1, l2], (18.0_f64).sqrt())
}

/// **NeuroShapeCert**: a PROVEN neural implicit shape. A small spectral-
/// normalized `tanh`-MLP SDF ([`fs_rep_neural`]) with certified Lipschitz
/// constant `L = 18` gives a no-tunnel sphere-trace step `|f|/L`; sound
/// Interval Bound Propagation proves a central box strictly inside (`hi < 0`)
/// and the four boundary strips of a bounding box each strictly outside
/// (`lo > 0`); tiled together (corners overlap) they wall off the interior into
/// a CLOSED frame `{f<0}` provably cannot cross — a single bounded component,
/// a proof not a mesh, strictly stronger than spot-checking discrete ring
/// boxes. A Morse cross-check ([`fs_viz`]) confirms a single interior minimum.
/// Runs the real `fs_neuroshape_e2e::run_campaign` for the certificate, then
/// renders a 64×64 SDF field for the viz.
///
/// `lift` — output bias (clamped `2.0..=12.0`, default 6.5; past ≈8.23 the
/// interior empties and the certificate flips `Verified → Estimated`); `ring_r`
/// — boundary-frame half-width (clamped `1.0..=4.0`, default 2.5); `inner` —
/// central box half-width (clamped `0.05..=1.0`, default 0.3).
///
/// Output layout (length `24 + 4096`):
/// - `[0]` — `grid_n` (64).
/// - `[1]`,`[2]` — `win_lo`, `win_hi` (the render window `[win_lo, win_hi]²`,
///   `win_lo = -(ring_r+0.5)`).
/// - `[3]` — `L` (certified Lipschitz constant, 18).
/// - `[4]` — `origin_value`.
/// - `[5]` — `safe_radius` (no-tunnel step `|f|/L`).
/// - `[6]` — `nearest_surface_radius` (`NaN` if no crossings).
/// - `[7]` — `max_crossing_radius`.
/// - `[8]`,`[9]` — `inside_lo`, `inside_hi` (IBP enclosure over the central
///   box).
/// - `[10]` — `certified_inside` (0/1 — `hi < 0`).
/// - `[11]` — `boundary_certified` (boundary strips proven strictly outside).
/// - `[12]` — `boundary_segments` (total strips forming the closed frame, 4).
/// - `[13]` — `bounded` (0/1 — closed barrier: every boundary strip certified
///   strictly outside, so the interior cannot escape).
/// - `[14]` — `single_minimum` (0/1 — Morse).
/// - `[15]` — `surface_crossings`.
/// - `[16]` — `topology_verified` (0/1).
/// - `[17]` — `component_count` (1 iff `certified_inside && bounded`, else 0).
/// - `[18]`,`[19]` — `ring_r`, `inner`.
/// - `[20..24]` — reserved (0).
/// - then `64·64` SDF field row-major (`j` outer / y, `i` inner / x) over the
///   render window.
pub fn neuroshape(lift: f64, ring_r: f64, inner: f64) -> Vec<f64> {
    let lift = lift.clamp(2.0, 12.0);
    let ring_r = ring_r.clamp(1.0, 4.0);
    let inner = inner.clamp(0.05, 1.0);
    let grid_n = 64usize;

    let net = neuro_net(lift);
    let report = fs_neuroshape_e2e::run_campaign(&net, ring_r, inner);
    let component_count = if report.certified_inside && report.bounded {
        1.0
    } else {
        0.0
    };

    let win_lo = -(ring_r + 0.5);
    let win_hi = ring_r + 0.5;
    let field = Grid2::from_fn(grid_n, grid_n, [win_lo, win_lo], [win_hi, win_hi], |p| {
        net.eval(&[p[0], p[1]])
    });

    let mut out = Vec::with_capacity(24 + grid_n * grid_n);
    out.push(grid_n as f64);
    out.push(win_lo);
    out.push(win_hi);
    out.push(fon(report.lipschitz));
    out.push(fon(report.origin_value));
    out.push(fon(report.safe_radius));
    out.push(fon(report.nearest_surface_radius));
    out.push(fon(report.max_crossing_radius));
    out.push(fon(report.inside_interval.0));
    out.push(fon(report.inside_interval.1));
    out.push(if report.certified_inside { 1.0 } else { 0.0 });
    out.push(report.boundary_certified as f64);
    out.push(report.boundary_segments as f64);
    out.push(if report.bounded { 1.0 } else { 0.0 });
    out.push(if report.single_minimum { 1.0 } else { 0.0 });
    out.push(report.surface_crossings as f64);
    out.push(if matches!(report.topology_color, Color::Verified { .. }) {
        1.0
    } else {
        0.0
    });
    out.push(component_count);
    out.push(ring_r);
    out.push(inner);
    out.extend_from_slice(&[0.0, 0.0, 0.0, 0.0]);
    for j in 0..grid_n {
        for i in 0..grid_n {
            out.push(fon(field.at(i, j)));
        }
    }
    out
}

/* ======================================================================= */
/*  8 · GrammarForge — fs-grammar-e2e (fs-shapeprog × fs-archive × fs-fab)   */
/* ======================================================================= */

/// The deterministic 7³ sample grid over `[-2, 2]³` (private in the campaign
/// crate; reconstructed here).
fn grammar_sample_points() -> Vec<[f64; 3]> {
    let mut pts = Vec::with_capacity(343);
    let n = 7;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let c = |t: usize| -2.0 + 4.0 * t as f64 / (n - 1) as f64;
                pts.push([c(i), c(j), c(k)]);
            }
        }
    }
    pts
}

/// **GrammarForge**: certified-fabricable geometric-program discovery. A
/// deterministic sweep of CSG programs ([`fs_shapeprog`]) is illuminated by
/// MAP-Elites ([`fs_archive`]) over (total material × dipole separation),
/// keeping the best-matching program in each of 6×4 niches. Each elite is
/// simplified with a fidelity certificate, then INDEPENDENTLY re-measured to
/// confirm the bound holds (certifying the certifier), and scored for
/// minimum-feature-size fabricability ([`fs_fab`]). The `run_campaign` body is
/// transcribed here so the niche grid and a representative best-program SDF can
/// be drawn.
///
/// `match_tol` — SDF-discrepancy match threshold (clamped `0.01..=1.0`,
/// default 0.2); `simplify_tol` — dropped-offset tolerance (clamped `0.0..=0.2`,
/// default 0.03; the live knob — at 0.03 the `o = 0` offsets drop 108→99).
///
/// Output layout (length `21 + 24 + 4096`):
/// - `[0]`,`[1]` — `r_bins` (6), `d_bins` (4).
/// - `[2]` — `num_elites`.
/// - `[3]` — `capacity` (24).
/// - `[4]` — `coverage`.
/// - `[5]` — `qd_score`.
/// - `[6]` — `best_discrepancy`.
/// - `[7..11]` — `best_params` `[r1, r2, d, o]`.
/// - `[11]` — `size_before` (total elite program size).
/// - `[12]` — `size_after`.
/// - `[13]` — `simplified_count`.
/// - `[14]` — `max_certified_error`.
/// - `[15]` — `simplification_sound` (0/1).
/// - `[16]` — `fab_satisfied`.
/// - `[17]` — `headline_verified` (0/1).
/// - `[18]` — `repr_grid_n` (64).
/// - `[19]`,`[20]` — `repr_lo` (-2), `repr_hi` (2).
/// - then `6·4` niche fitness grid row-major (r-bin outer, d-bin inner;
///   `-1` = empty niche).
/// - then `64·64` SDF slice (`z = 0`) of the best program over `[-2,2]²`,
///   row-major (`j` outer / y, `i` inner / x).
pub fn grammarforge(match_tol: f64, simplify_tol: f64) -> Vec<f64> {
    let match_tol = match_tol.clamp(0.01, 1.0);
    let simplify_tol = simplify_tol.clamp(0.0, 0.2);

    let target = fs_grammar_e2e::target();
    let samples = grammar_sample_points();
    let fab = min_feature_size(0.8);

    let r_vals = [0.7, 0.9, 1.0, 1.1];
    let d_vals = [0.6, 0.8, 1.0];
    let o_vals = [0.0, 0.02, 0.05];

    let mut archive = MapElites::new(vec![1.3, 0.5], vec![2.3, 1.1], vec![6, 4]);
    for &r1 in &r_vals {
        for &r2 in &r_vals {
            for &d in &d_vals {
                for &o in &o_vals {
                    let prog = fs_grammar_e2e::build_program(r1, r2, d, o);
                    let disc = max_sdf_discrepancy(&prog, &target, &samples);
                    let fitness = 1.0 / (1.0 + disc);
                    if fitness.is_finite() && fitness >= 0.0 {
                        archive.add(vec![r1, r2, d, o], vec![r1 + r2, d], fitness);
                    }
                }
            }
        }
    }

    // Niche grid (6 r-bins × 4 d-bins), row-major r-bin outer.
    let mut grid = vec![-1.0f64; 24];
    for e in archive.elites() {
        let cell = archive.cell_of(&e.descriptor);
        let (rb, db) = (cell[0], cell[1]);
        grid[rb * 4 + db] = e.fitness;
    }

    // Post-process the elites: simplification soundness + fabrication tally.
    let (mut simplified_count, mut size_before, mut size_after, mut fab_satisfied) = (0, 0, 0, 0);
    let mut max_certified_error = 0.0_f64;
    let mut simplification_sound = true;
    for e in archive.elites() {
        let prog = fs_grammar_e2e::build_program(
            e.solution[0],
            e.solution[1],
            e.solution[2],
            e.solution[3],
        );
        let before = prog.size();
        let simp = simplify(&prog, simplify_tol);
        let after = simp.program.size();
        size_before += before;
        size_after += after;
        if after < before {
            simplified_count += 1;
        }
        max_certified_error = max_certified_error.max(simp.max_error);
        let actual = max_sdf_discrepancy(&prog, &simp.program, &samples);
        if actual > simp.max_error + 1e-9 {
            simplification_sound = false;
        }
        if fab.satisfied(e.solution[0].min(e.solution[1])) {
            fab_satisfied += 1;
        }
    }

    let (best_disc, best_params, best_fab_ok) = match archive.best() {
        Some(best) => {
            let bd = 1.0 / best.fitness - 1.0;
            let bp = [
                best.solution[0],
                best.solution[1],
                best.solution[2],
                best.solution[3],
            ];
            (
                bd,
                bp,
                fab.satisfied(best.solution[0].min(best.solution[1])),
            )
        }
        None => (f64::NAN, [f64::NAN; 4], false),
    };
    let headline_verified = best_disc <= match_tol && best_fab_ok && simplification_sound;

    // Representative shape: best program z=0 slice, 64×64 over [-2,2]².
    let repr_n = 64usize;
    let repr = fs_grammar_e2e::build_program(
        best_params[0],
        best_params[1],
        best_params[2],
        best_params[3],
    );

    let mut out = Vec::with_capacity(21 + 24 + repr_n * repr_n);
    out.push(6.0);
    out.push(4.0);
    out.push(archive.num_elites() as f64);
    out.push(archive.capacity() as f64);
    out.push(fon(archive.coverage()));
    out.push(fon(archive.qd_score()));
    out.push(fon(best_disc));
    out.extend_from_slice(&best_params);
    out.push(size_before as f64);
    out.push(size_after as f64);
    out.push(simplified_count as f64);
    out.push(fon(max_certified_error));
    out.push(if simplification_sound { 1.0 } else { 0.0 });
    out.push(fab_satisfied as f64);
    out.push(if headline_verified { 1.0 } else { 0.0 });
    out.push(repr_n as f64);
    out.push(-2.0);
    out.push(2.0);
    out.extend_from_slice(&grid);
    for j in 0..repr_n {
        let y = -2.0 + 4.0 * j as f64 / (repr_n as f64 - 1.0);
        for i in 0..repr_n {
            let x = -2.0 + 4.0 * i as f64 / (repr_n as f64 - 1.0);
            out.push(fon(repr.sdf([x, y, 0.0])));
        }
    }
    out
}

/* ======================================================================= */
/*  9 · AnytimeBO — fs-adaptbo-e2e (fs-bo × fs-eproc)                        */
/* ======================================================================= */

fn bo_argmin(xs: &[Vec<f64>], ys: &[f64]) -> (f64, f64) {
    let mut bi = 0usize;
    for i in 1..ys.len() {
        if ys[i] < ys[bi] {
            bi = i;
        }
    }
    (xs[bi][0], ys[bi])
}

/// **AnytimeBO**: Bayesian optimization that provably knows when to stop. A
/// Matérn-5⁄2 GP with closed-form Expected Improvement ([`fs_bo`]) drives a
/// deterministic minimization of a tilted double well on `[0,4]`; a betting
/// e-process ([`fs_eproc`]) watches a per-iteration STALL indicator and stops
/// the search the instant its log-e-value crosses the Ville threshold
/// `−ln α` — an anytime-valid decision (no alpha-spending). The loop is
/// transcribed here to capture the per-iteration `[x, y, incumbent, log_e]`
/// trajectory.
///
/// `max_iters` — iteration cap (clamped `1..=40`, default 30); `delta` —
/// improvement threshold for "not stalled" (clamped `1e-6..=1.0`, default
/// 0.02); `alpha` — anytime level (clamped `1e-3..=0.5`, default 0.05).
///
/// Output layout (a flat array the viz slices; `I` = iterations run,
/// `G` = grid = 81):
/// - `[0]` — `iters_run` (`I`).
/// - `[1]` — `ville_threshold` (`−ln α`).
/// - `[2]` — `stopped_early` (0/1).
/// - `[3]` — `best_x`.
/// - `[4]` — `best_value`.
/// - `[5]` — `evaluations` (`3 + I`).
/// - `[6]`,`[7]` — `ci_center`, `ci_radius` (anytime confidence sequence on
///   the best-value trace).
/// - `[8]` — `n_init` (3); then 3 `[x_init, y_init]` pairs.
/// - then `I`; then `I` blocks of 4: `[x, y, incumbent, log_e]`.
/// - then `G` (81); then `G` blocks of 2: `[x, objective(x)]`.
pub fn anytimebo(max_iters: usize, delta: f64, alpha: f64) -> Vec<f64> {
    let max_iters = max_iters.clamp(1, 40);
    let delta = delta.clamp(1e-6, 1.0);
    let alpha = alpha.clamp(1e-3, 0.5);

    let grid: Vec<f64> = (0..=80).map(|i| 4.0 * f64::from(i) / 80.0).collect();
    let mut xs: Vec<Vec<f64>> = vec![vec![0.4], vec![2.6], vec![3.6]];
    let mut ys: Vec<f64> = xs.iter().map(|x| fs_adaptbo_e2e::objective(x[0])).collect();

    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.5],
    };
    let mut eproc = BettingEProcess::new(0.3);
    let mut cs = GaussianMixtureCs::new(1.0, 4.0, alpha);
    let (_, mut best_value) = bo_argmin(&xs, &ys);
    cs.observe(best_value);

    let mut iter_rows: Vec<[f64; 4]> = Vec::new();
    let mut stopped_early = false;
    for _ in 0..max_iters {
        let Some(gp) = Gp::try_fit(&xs, &ys, kernel.clone(), 1e-6) else {
            break;
        };
        let f_best = ys.iter().copied().fold(f64::INFINITY, f64::min);
        let mut best_ei = f64::NEG_INFINITY;
        let mut x_next = grid[0];
        for &g in &grid {
            let ei = expected_improvement(&gp, &[g], f_best, 0.01);
            if ei > best_ei {
                best_ei = ei;
                x_next = g;
            }
        }
        let y_next = fs_adaptbo_e2e::objective(x_next);
        xs.push(vec![x_next]);
        ys.push(y_next);

        let new_best = best_value.min(y_next);
        let improvement = best_value - new_best;
        best_value = new_best;
        cs.observe(new_best);
        let stall = if improvement < delta { 1.0 } else { 0.0 };
        eproc.observe(stall);
        iter_rows.push([x_next, y_next, new_best, eproc.log_e_value()]);
        if eproc.rejects_at(alpha) {
            stopped_early = true;
            break;
        }
    }

    let (best_x, best_val) = bo_argmin(&xs, &ys);
    let (ci_center, ci_radius) = cs.interval().unwrap_or((best_val, f64::INFINITY));
    let ville = -alpha.ln();
    let i = iter_rows.len();
    let g = grid.len();

    let mut out = Vec::with_capacity(8 + 1 + 6 + 1 + 4 * i + 1 + 2 * g);
    out.push(i as f64);
    out.push(fon(ville));
    out.push(if stopped_early { 1.0 } else { 0.0 });
    out.push(fon(best_x));
    out.push(fon(best_val));
    out.push(ys.len() as f64);
    out.push(fon(ci_center));
    out.push(fon(ci_radius));
    out.push(3.0);
    out.push(0.4);
    out.push(fon(fs_adaptbo_e2e::objective(0.4)));
    out.push(2.6);
    out.push(fon(fs_adaptbo_e2e::objective(2.6)));
    out.push(3.6);
    out.push(fon(fs_adaptbo_e2e::objective(3.6)));
    out.push(i as f64);
    for r in &iter_rows {
        out.push(fon(r[0]));
        out.push(fon(r[1]));
        out.push(fon(r[2]));
        out.push(fon(r[3]));
    }
    out.push(g as f64);
    for &x in &grid {
        out.push(x);
        out.push(fon(fs_adaptbo_e2e::objective(x)));
    }
    out
}

/* ======================================================================= */
/*  10 · FlowCert — fs-flowcert-e2e (fs-lbm × fs-archive)                    */
/* ======================================================================= */

/// One certified LBM operating point + its captured velocity profile.
struct FlowPoint {
    reynolds: f64,
    ny: usize,
    tau: f64,
    viscosity: f64,
    profile_error: f64,
    accurate: bool,
    regime_stable: bool,
    numeric: Vec<f64>,
    analytic: Vec<f64>,
}

/// Faithful transcription of `fs_flowcert_e2e::certify_point`, additionally
/// capturing the numeric & analytic velocity profiles.
fn flow_certify(reynolds: f64, ny: usize, u_lattice: f64, steps: usize, tol: f64) -> FlowPoint {
    let plan = plan_scaling(reynolds, ny as f64, u_lattice);
    let nu = plan.viscosity;
    let tau = plan.tau;
    let gx = 8.0 * nu * u_lattice / (ny as f64).powi(2);

    let mut lbm = Lbm::channel(4, ny, tau, gx);
    lbm.run(steps);
    let profile = lbm.x_velocity_profile();

    let mut peak = 0.0_f64;
    let mut max_err = 0.0_f64;
    let mut analytic = Vec::with_capacity(ny);
    for (y, &u) in profile.iter().enumerate() {
        let exact = poiseuille_analytic(gx, nu, ny, y);
        peak = peak.max(exact.abs());
        max_err = max_err.max((u - exact).abs());
        analytic.push(exact);
    }
    let profile_error = if peak > 1e-12 {
        max_err / peak
    } else {
        max_err
    };

    FlowPoint {
        reynolds,
        ny,
        tau,
        viscosity: nu,
        profile_error,
        accurate: profile_error <= tol,
        regime_stable: plan.color().rank() == ColorRank::Verified,
        numeric: profile,
        analytic,
    }
}

/// **FlowCert**: a certified credibility map for a lattice-Boltzmann channel
/// flow. Each (Reynolds × resolution) operating point is run to steady state
/// ([`fs_lbm`]) and compared to the ANALYTIC Poiseuille solution (a
/// manufactured-solution accuracy certificate); the scaling planner flags the
/// regime `Verified` only when comfortably stable. MAP-Elites ([`fs_archive`])
/// illuminates the atlas. The `certify_point` body is transcribed here so the
/// spotlight velocity profiles (numeric vs analytic) can be drawn. Fixed sweep
/// `Re ∈ {20, 60, 120}`, `ny ∈ {16, 24, 32}` (9 points); low-Re credible,
/// high-Re flagged.
///
/// `steps` — LBM steps to steady state (clamped `500..=12000`, default 12000);
/// `tol` — relative-error accuracy tolerance (clamped `1e-3..=0.5`, default
/// 0.03).
///
/// Output layout (a flat array the viz slices; `P` = 9, `S` = 2 spotlights):
/// - `[0]` — `P`.
/// - `[1]` — `coverage`.
/// - `[2]` — `qd_score`.
/// - `[3]` — `num_niches`.
/// - `[4]` — `best_error`.
/// - `[5]` — `stable_fraction`.
/// - `[6]` — `all_accurate` (0/1).
/// - `[7]` — `map_color_rank` (2/1/0 — the credibility color).
/// - then `P` blocks of 8: `[Re, ny, max_error, accurate, regime_stable,
///   fully_credible, tau, viscosity]` (`fully_credible = accurate &&
///   regime_stable`).
/// - then `S` (2); then each spotlight: `[Re, ny, then ny pairs of
///   (u_numeric, u_analytic)]` — spotlight 0 = (Re 20, ny 32) credible,
///   spotlight 1 = (Re 120, ny 32) flagged.
pub fn flowcert(steps: usize, tol: f64) -> Vec<f64> {
    let steps = steps.clamp(500, 12_000);
    let tol = tol.clamp(1e-3, 0.5);
    let u_lat = 0.05;

    let reynolds = [20.0f64, 60.0, 120.0];
    let resolutions = [16usize, 24, 32];
    let re_lo = 20.0;
    let re_hi = 120.0;
    let ny_lo = 16.0;
    let ny_hi = 32.0;
    let mut archive = MapElites::new(
        vec![re_lo - 1.0, ny_lo - 1.0],
        vec![re_hi + 1.0, ny_hi + 1.0],
        vec![3, 3],
    );

    let mut points: Vec<FlowPoint> = Vec::with_capacity(9);
    let mut spot_a: Option<usize> = None; // (20, 32)
    let mut spot_b: Option<usize> = None; // (120, 32)
    for &re in &reynolds {
        for &ny in &resolutions {
            let p = flow_certify(re, ny, u_lat, steps, tol);
            let fitness = 1.0 / (1.0 + p.profile_error);
            if fitness.is_finite() && fitness >= 0.0 {
                archive.add(vec![re, ny as f64], vec![re, ny as f64], fitness);
            }
            if (re - 20.0).abs() < 1e-9 && ny == 32 {
                spot_a = Some(points.len());
            }
            if (re - 120.0).abs() < 1e-9 && ny == 32 {
                spot_b = Some(points.len());
            }
            points.push(p);
        }
    }

    let n = points.len();
    let all_accurate = points.iter().all(|p| p.accurate);
    let stable_count = points.iter().filter(|p| p.regime_stable).count();
    let stable_fraction = stable_count as f64 / n as f64;
    let best_error = points
        .iter()
        .map(|p| p.profile_error)
        .fold(f64::INFINITY, f64::min);
    let map_rank = if all_accurate && stable_count == n {
        2.0
    } else {
        0.0
    };

    let spots: Vec<usize> = [spot_a, spot_b].into_iter().flatten().collect();
    let spot_len: usize = spots.iter().map(|&i| 2 + 2 * points[i].ny).sum();
    let mut out = Vec::with_capacity(8 + 8 * n + 1 + spot_len);
    out.push(n as f64);
    out.push(fon(archive.coverage()));
    out.push(fon(archive.qd_score()));
    out.push(archive.num_elites() as f64);
    out.push(fon(best_error));
    out.push(fon(stable_fraction));
    out.push(if all_accurate { 1.0 } else { 0.0 });
    out.push(map_rank);
    for p in &points {
        let fully = p.accurate && p.regime_stable;
        out.push(p.reynolds);
        out.push(p.ny as f64);
        out.push(fon(p.profile_error));
        out.push(if p.accurate { 1.0 } else { 0.0 });
        out.push(if p.regime_stable { 1.0 } else { 0.0 });
        out.push(if fully { 1.0 } else { 0.0 });
        out.push(fon(p.tau));
        out.push(fon(p.viscosity));
    }
    out.push(spots.len() as f64);
    for &si in &spots {
        let p = &points[si];
        out.push(p.reynolds);
        out.push(p.ny as f64);
        for y in 0..p.ny {
            out.push(fon(p.numeric.get(y).copied().unwrap_or(f64::NAN)));
            out.push(fon(p.analytic.get(y).copied().unwrap_or(f64::NAN)));
        }
    }
    out
}

/* ======================================================================= */
/*  Regression tests — the certified headline numbers must reproduce.       */
/* ======================================================================= */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn proofrobust_defaults() {
        let v = proofrobust(0.9, 2.0, 41);
        assert_eq!(v[0], 3.0, "F");
        assert_eq!(v[1], 3.0, "certified_count");
        assert_eq!(v[2], 1.0, "reorders");
        // headline rank Estimated (CVaR), winners differ.
        assert_eq!(v[3], 0.0, "headline_rank");
        assert_ne!(v[4], v[5], "nominal_winner_idx != robust_winner_idx");
        assert_eq!(v[4], 0.0, "nominal winner = champion (idx 0)");
        assert_eq!(v[5], 1.0, "robust winner = flat (idx 1)");
    }

    #[test]
    fn metamatcert_defaults() {
        let v = metamatcert(10, 6, 0.40);
        let p = v[0] as usize;
        assert_eq!(p, 6);
        assert_eq!(v[2], 1.0, "all_stable");
        assert_eq!(v[3], 1.0, "all_admissible");
        // c11 first (r=0) ~3.5, last (r=0.40) ~0.8; decreasing.
        let first_c11 = v[7 + 2];
        let last_c11 = v[7 + 6 * (p - 1) + 2];
        assert!(first_c11 > 3.0 && first_c11 < 4.0, "first c11 {first_c11}");
        assert!(last_c11 > 0.5 && last_c11 < 1.2, "last c11 {last_c11}");
        assert!(first_c11 > last_c11, "c11 decreasing");
    }

    #[test]
    fn fluttercert_defaults() {
        let v = fluttercert(0.55, 2.45, 20);
        // boundaries agree, near mu=2.
        assert_eq!(v[3], 1.0, "boundaries_agree");
        assert!((v[1] - v[2]).abs() < 1e-9, "lyapunov ~ spectral");
        assert!(v[1] > 1.7 && v[1] < 2.05, "lyapunov_boundary {}", v[1]);
    }

    #[test]
    fn schedule_defaults() {
        let v = schedule_campaign(12.0, 0.65, 1e-6);
        assert_eq!(v[0], 13.0, "makespan");
        assert!(v[1] <= v[0], "lower bound {} > {}", v[1], v[0]);
        assert!(v[0] <= v[2], "{} > upper bound {}", v[0], v[2]);
        assert_eq!(v[8], 0.0, "should_stop = false");
        assert_eq!(v[9], 0.0, "leading design A (idx 0)");
        assert_eq!(v[10], 0.0, "rec = Act");
        // critical path [3, 4]; bottleneck study idx 3 (windtunnel-A).
        assert_eq!(v[4], 2.0, "path len 2");
        assert_eq!(v[5], 3.0, "bottleneck idx 3");
    }

    #[test]
    fn trusspath_defaults() {
        let v = trusspath(4, 3, 1e-4);
        let m = v[0] as usize;
        assert_eq!(m, 43, "members");
        assert_eq!(v[1], 6.0, "active");
        assert!(v[3] < 1e-3, "gap {}", v[3]);
        assert_eq!(v[6], 1.0, "solver_converged");
        assert!(v[7] >= 2.0, "connected path length {}", v[7]);
        assert!(v[8].is_finite() && v[8] > 0.0, "path weight {}", v[8]);
        let claim = &v[v.len() - 4..];
        assert_eq!(claim[0], 2.0, "optimality rank must be Verified");
        assert_eq!(claim[1], 1.0, "verified interval flag");
        assert!(claim[2].is_finite() && claim[3].is_finite());
        assert!(claim[2] > 0.0 && claim[2] <= claim[3]);

        let native = with_certificate_cx(|cx| {
            fs_truss_e2e::run_campaign(4, 3, 4.0, 2.0, 1e-4, cx).expect("native TrussPath campaign")
        });
        let native_claim = verified_bounds(&native.optimality_color);
        assert_eq!(claim[0], rank_code(native.optimality_color.rank()));
        assert_eq!(claim[1], native_claim.0);
        assert_eq!(claim[2].to_bits(), native_claim.1.to_bits());
        assert_eq!(claim[3].to_bits(), native_claim.2.to_bits());
    }

    #[test]
    fn trusspath_non_finite_tolerance_uses_the_bounded_default() {
        let default = trusspath(4, 3, 1e-4);
        let non_finite = trusspath(4, 3, f64::NAN);
        assert_eq!(non_finite, default);
    }

    #[test]
    fn sensorforge_defaults() {
        let v = sensorforge(0.01, 12, 0.65);
        let s = v[1] as usize;
        assert_eq!(v[7], 1.0, "decision_robust");
        assert_eq!(v[8], 0.0, "chosen A (idx 0)");
        let initial = v[5];
        let final_evpi = v[6];
        assert!((initial - 0.163).abs() < 0.02, "initial_evpi {initial}");
        assert!(
            (final_evpi - 0.0097).abs() < 0.005,
            "final_evpi {final_evpi}"
        );
        // trace[0] = initial, trace[last] = final.
        let t = v[9] as usize;
        assert_eq!(t, s + 1, "trace len");
        let trace0 = v[10];
        let trace_last = v[10 + t - 1];
        assert!((trace0 - initial).abs() < 1e-12);
        assert!((trace_last - final_evpi).abs() < 1e-12);
    }

    #[test]
    fn neuroshape_defaults() {
        let v = neuroshape(6.5, 2.5, 0.3);
        assert_eq!(v[0], 64.0, "grid_n");
        assert!((v[3] - 18.0).abs() < 1e-6, "L {}", v[3]);
        assert_eq!(v[16], 1.0, "topology_verified");
        assert_eq!(v[17], 1.0, "component_count");
        assert_eq!(v.len(), 24 + 64 * 64, "total length");
    }

    #[test]
    fn grammarforge_defaults() {
        let v = grammarforge(0.2, 0.03);
        assert_eq!(v[2], 18.0, "num_elites");
        assert_eq!(v[3], 24.0, "capacity");
        assert_eq!(v[11], 108.0, "size_before");
        assert_eq!(v[12], 99.0, "size_after");
        assert_eq!(v[15], 1.0, "simplification_sound");
        assert_eq!(v[17], 1.0, "headline_verified");
        assert_eq!(v.len(), 21 + 24 + 64 * 64, "total length");
    }

    #[test]
    fn anytimebo_defaults() {
        let v = anytimebo(30, 0.02, 0.05);
        let iters = v[0] as usize;
        assert_eq!(v[2], 1.0, "stopped_early");
        assert!((v[1] - 2.9957).abs() < 0.001, "ville {}", v[1]);
        assert!((v[3] - 3.0).abs() < 0.1, "best_x {}", v[3]);
        assert!((v[4] + 0.45).abs() < 0.1, "best_value {}", v[4]);
        assert!((8..=16).contains(&iters), "iters {iters}");
    }

    #[test]
    fn flowcert_defaults() {
        let v = flowcert(12_000, 0.03);
        assert_eq!(v[0], 9.0, "P");
        // points row-major: Re=20 rows are 0,1,2 (fully_credible=1),
        // Re=120 rows are 6,7,8 (fully_credible=0).
        let block = |i: usize| &v[8 + 8 * i..8 + 8 * (i + 1)];
        for i in 0..3 {
            let b = block(i);
            assert_eq!(b[0], 20.0, "Re=20 row");
            assert_eq!(b[5], 1.0, "Re=20 fully_credible");
        }
        for i in 6..9 {
            let b = block(i);
            assert_eq!(b[0], 120.0, "Re=120 row");
            assert_eq!(b[5], 0.0, "Re=120 flagged");
        }
    }
}
