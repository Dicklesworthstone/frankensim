//! fs-flagship-e2e — the flagship E2E suite (bead mye.5): all
//! flagships as REPLAYABLE, STAGED, FORENSICALLY-LOGGED end-to-end
//! test assets. Layer: L6 (HELM).
//!
//! - STAGED FIDELITY: every flagship runs a SMOKE stage (PR-gate,
//!   minutes) here and now; MID (nightly) and FULL (weekly) stages are
//!   WIRED with envelopes and `#[ignore]` markers — their CI lanes
//!   belong to the perf-CI bead (fz2.4), and this crate does not
//!   pretend otherwise.
//! - GOLDEN-LEDGER DISCIPLINE: each smoke stage folds its metric
//!   stream into a content hash (canonical replay identity over metric bits); CI replays
//!   and compares hashes — every stage here is deterministic, so hash
//!   equality IS the gate (stochastic-labeled stages would gate on
//!   envelopes instead).
//! - SHARED-CORE AUDITS: one canonical uniform-tau D2Q9 roll hash for
//!   the collide/stream core the vessel and the ornithoid both ride
//!   (its coverage is bounded — see [`lbm_core_roll_hash`]), and one
//!   e-race audit that replays the race core and checks the ornithoid's
//!   public screening wrapper against it under the ornithoid's own
//!   declared span. `fs-vessel` exposes no public race wrapper — its
//!   convention lives in its own battery — so NO cross-flagship
//!   convention comparison is claimed here.
//! - FAILURE DRILLS: cancellation storms, budget exhaustion,
//!   ledger crash-recovery, model-form escalation — each with an
//!   EXPECTED STRUCTURED OUTCOME that is MEASURED, not a literal
//!   written into the row.
//! - FORENSIC LOGGING: every stage emits structured JSON rows
//!   (metrics, certificates, race records, timings) sufficient to
//!   diagnose failures from logs + ledger alone; the battery parses
//!   its own stream as a self-audit; the LAB NOTEBOOK artifact
//!   regenerates bitwise on replay (timings ride a separate,
//!   non-golden row — wall-clock is evidence, not identity).

use std::fmt::Write as _;

fn push_json_string(out: &mut String, value: &str) {
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch <= '\u{001f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(ch));
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
}

/// Stage fidelity tiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tier {
    /// PR-gate, minutes.
    Smoke,
    /// Nightly, hour-class.
    Mid,
    /// Weekly/on-demand, production scale.
    Full,
}

/// One stage's content-addressed artifact.
#[derive(Debug, Clone)]
pub struct StageArtifact {
    /// Flagship name.
    pub flagship: &'static str,
    /// Fidelity tier.
    pub tier: Tier,
    /// Named metrics, in emission order (the hashed content).
    pub metrics: Vec<(&'static str, f64)>,
    /// FNV-64 over the metric bit patterns (the golden identity).
    pub hash: u64,
    /// Wall-clock seconds (logged, NEVER hashed).
    pub wall_s: f64,
}

/// Metric-stream encoding version — a SEMANTIC SURFACE registered in
/// golden-couplings.json: every flagship golden depends on it, so
/// changing the encoding must deliberately re-freeze them all.
/// v1 = bare name/bits concatenation (unprefixed — "ab"+bits vs
/// "a"+... could collide); v2 = canonical replay identity (gp3.14):
/// typed, length-prefixed fs_obs::ident fields.
pub const METRIC_STREAM_ENCODING_VERSION: u32 = 2;

/// Fold a metric stream into the content hash (canonical replay
/// identity over `(name, bits)` fields, order-semantic).
#[must_use]
pub fn content_hash(metrics: &[(&'static str, f64)]) -> u64 {
    let mut b = fs_obs::ident::IdentityBuilder::new("flagship-metric-stream");
    for (name, v) in metrics {
        b = b.f64_bits(name, *v);
    }
    b.finish().root()
}

/// Build an artifact (hash computed, wall clock attached).
#[must_use]
pub fn artifact(
    flagship: &'static str,
    tier: Tier,
    metrics: Vec<(&'static str, f64)>,
    wall_s: f64,
) -> StageArtifact {
    let hash = content_hash(&metrics);
    StageArtifact {
        flagship,
        tier,
        metrics,
        hash,
        wall_s,
    }
}

/// One forensic log row: structured JSON with the suite's required
/// keys (`stage`, `kind`, `payload`). `stage` and `kind` are escaped;
/// `payload` must already be one complete JSON value.
#[must_use]
pub fn log_row(stage: &str, kind: &str, payload: &str) -> String {
    let mut row = String::from("{\"stage\":");
    push_json_string(&mut row, stage);
    row.push_str(",\"kind\":");
    push_json_string(&mut row, kind);
    row.push_str(",\"payload\":");
    row.push_str(payload);
    row.push('}');
    row
}

/// The LAB NOTEBOOK artifact: deterministic JSON over the stages'
/// golden content (hashes + metrics); timings are emitted as a
/// SEPARATE non-golden section so the notebook body replays bitwise.
#[must_use]
pub fn notebook(artifacts: &[StageArtifact]) -> String {
    let mut body = String::from("{\"suite\":\"flagship-e2e\",\"stages\":[");
    for (i, a) in artifacts.iter().enumerate() {
        if i > 0 {
            body.push(',');
        }
        body.push_str("{\"flagship\":");
        push_json_string(&mut body, a.flagship);
        let _ = write!(
            body,
            ",\"tier\":\"{:?}\",\"hash\":\"0x{:016x}\",\"metrics\":{{",
            a.tier, a.hash
        );
        for (j, (name, v)) in a.metrics.iter().enumerate() {
            if j > 0 {
                body.push(',');
            }
            push_json_string(&mut body, name);
            let _ = write!(body, ":\"0x{:016x}\"", v.to_bits());
        }
        body.push_str("}}");
    }
    body.push_str("]}");
    body
}

/// A canonical UNIFORM-TAU D2Q9 BGK roll (periodic in `x`, wall-bounded
/// in `y`) giving the vessel and ornithoid consumers ONE shared audit
/// point for the collide/stream core, with one constant to bump (with
/// justification) when that path moves.
///
/// COVERAGE, exactly: `Grid::uniform(24, 16, 0.6)` — one relaxation
/// time, no external force, no rheology — with full-width `Wall` rows at
/// `y = 0` and `y = ny-1`, 50 plain `grid.step()` calls, hashing
/// `rho`/`ux`/`uy` at 6 probe cells. It therefore pins uniform-tau BGK
/// collide + plain stream + top/bottom bounce-back on a periodic-x
/// lattice, and NOTHING ELSE.
///
/// It does NOT cover: `Rheology` (the vessel's `run_pour` passes
/// `Rheology::Newtonian`), `ContactModel`/free surface, interior
/// rasterized-obstacle bounce-back (`fs_ornith::refine`), non-periodic-x
/// with pinned inlet/outlet equilibrium columns, momentum-exchange or
/// partial-saturation variants, or D3Q19. A change confined to those
/// paths will not move this hash — the repo's own history is the
/// counterexample: the `xo2k` migration of the fs-lbm rheology `powf`
/// paths to `det::pow` changed shared machinery the vessel rides and
/// `GOLDEN_LBM_CORE` did not move for it. Extending coverage means
/// adding fixtures, not widening this sentence (bead
/// `frankensim-extreal-program-f85xj.2.36`).
#[must_use]
pub fn lbm_core_roll_hash() -> u64 {
    use fs_lbm::{Cell, Grid, Q, equilibrium};
    let (nx, ny) = (24usize, 16usize);
    let mut grid = Grid::uniform(nx, ny, 0.6);
    grid.periodic_x = true;
    grid.periodic_y = false;
    // Walls top/bottom, shear-ish init.
    for x in 0..nx {
        let b = grid.idx(x, 0);
        grid.flags[b] = Cell::Wall;
        let t = grid.idx(x, ny - 1);
        grid.flags[t] = Cell::Wall;
    }
    for y in 1..ny - 1 {
        for x in 0..nx {
            let i = grid.idx(x, y);
            let u = 0.04 * (y as f64 / ny as f64 - 0.5);
            grid.f[i] = equilibrium(1.0, u, 0.0);
        }
    }
    let mut scratch: Vec<[f64; Q]> = Vec::new();
    for _ in 0..50 {
        grid.step(&mut scratch);
    }
    let mut metrics: Vec<(&'static str, f64)> = Vec::new();
    for y in [1usize, ny / 2, ny - 2] {
        for x in [0usize, nx / 2] {
            let m = grid.moments(grid.idx(x, y));
            metrics.push(("rho", m.rho));
            metrics.push(("ux", m.u[0]));
            metrics.push(("uy", m.u[1]));
        }
    }
    content_hash(&metrics)
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
