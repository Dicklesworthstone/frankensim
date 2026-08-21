//! E4.8 — fs-bem coarse SCREENING preset + the one-shot interference
//! derivation cache (bead wf-root-guzez.5.20).
//!
//! Laws:
//!   - the screening preset carries a DECLARED panel budget and
//!     refuses on exhaustion (never a silent truncation of geometry);
//!   - the one-shot interference derivation on design commit binds
//!     the FULL cache key — geometry, operating grid, panel preset,
//!     ground mode, solver version, coefficient convention — and a
//!     result may be reused ONLY under an exactly matching key (the
//!     stale-key hostile twin refuses);
//!   - design commits bump an epoch and CANCEL stale jobs: a job
//!     minted under an older epoch refuses to deliver;
//!   - slider drags get the SCHEMATIC preview, a struct that carries
//!     NO solver numbers by construction.

use crate::{Refusal, refuse};
use fs_bem::panel3d::surface_velocity;
use fs_bem::{SpherePanels, solve_exterior};
use fs_blake3::hash_domain;

/// Screening panel budget cap (browser design-tool tier).
pub const MAX_SCREEN_PANELS: usize = 1_280;

/// Solver version pinned into every cache key.
pub const SCREEN_SOLVER_VERSION: &str = "fs-bem-exterior-gmres-v1";

/// The screening result (coarse tier, declared).
#[derive(Clone, Debug, PartialEq)]
pub struct ScreeningResult {
    /// Panels actually solved.
    pub n_panels: usize,
    /// Max surface speed ratio |u|/|U_inf| (the screening observable).
    pub max_speed_ratio: f64,
    /// Digest over the panel speeds (identity).
    pub result_digest: String,
}

/// Run the coarse screening solve on an icosphere proxy at the given
/// subdivision, under a DECLARED panel budget.
///
/// # Errors
/// `bem-screen-budget-exhausted` (panel count beyond the budget — AT
/// the budget admits, one fewer than needed refuses);
/// `bem-screen-solve` (fs-bem refusals pass through, typed).
pub fn screening_solve(
    subdivisions: u32,
    max_panels: usize,
    u_inf_mps: f64,
) -> Result<ScreeningResult, Refusal> {
    let panels = SpherePanels::icosphere(1.0, subdivisions)
        .map_err(|e| refuse("bem-screen-solve", format!("{e:?}"), "the icosphere proxy"))?;
    let n = panels.areas().len();
    if n > max_panels.min(MAX_SCREEN_PANELS) {
        return Err(refuse(
            "bem-screen-budget-exhausted",
            format!("{n} panels > budget {}", max_panels.min(MAX_SCREEN_PANELS)),
            "coarsen the preset; geometry is never silently truncated",
        ));
    }
    let u_inf = [u_inf_mps, 0.0, 0.0];
    let sol = solve_exterior(&panels, u_inf, 6, 1e-10)
        .map_err(|e| refuse("bem-screen-solve", format!("{e:?}"), "the exterior solve"))?;
    let vel = surface_velocity(&panels, &sol.sigma, u_inf, 6)
        .map_err(|e| refuse("bem-screen-solve", format!("{e:?}"), "surface velocity"))?;
    let mut max_ratio = 0.0f64;
    let mut b = Vec::with_capacity(8 * vel.len());
    for v in &vel {
        let s = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() / u_inf_mps;
        max_ratio = max_ratio.max(s);
        b.extend_from_slice(&s.to_bits().to_le_bytes());
    }
    Ok(ScreeningResult {
        n_panels: n,
        max_speed_ratio: max_ratio,
        result_digest: hash_domain("org.frankensim.wf.bem-screen.v1", &b).to_hex(),
    })
}

/// The FULL one-shot derivation cache key (every field load-bearing).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct InterferenceCacheKey {
    /// Geometry content digest.
    pub geometry_digest: String,
    /// Operating-grid digest.
    pub operating_grid_digest: String,
    /// Panel preset (subdivision tier).
    pub panel_preset: u32,
    /// Ground mode ("free-air" | "image-flat").
    pub ground_mode: &'static str,
    /// Solver version.
    pub solver_version: &'static str,
    /// Coefficient convention id.
    pub coefficient_convention: &'static str,
}

impl InterferenceCacheKey {
    /// Key digest (the cache index).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut b = self.geometry_digest.as_bytes().to_vec();
        b.extend_from_slice(self.operating_grid_digest.as_bytes());
        b.extend_from_slice(&self.panel_preset.to_le_bytes());
        b.extend_from_slice(self.ground_mode.as_bytes());
        b.extend_from_slice(self.solver_version.as_bytes());
        b.extend_from_slice(self.coefficient_convention.as_bytes());
        hash_domain("org.frankensim.wf.interference-key.v1", &b).to_hex()
    }
}

/// A derived interference result bound to the key that produced it.
#[derive(Clone, Debug, PartialEq)]
pub struct InterferenceResult {
    /// The key this result was derived under.
    pub key_digest: String,
    /// The derived interference factor (screening observable).
    pub interference_factor: f64,
}

/// The design-commit cache: one-shot derivations, epoch-guarded jobs.
#[derive(Debug, Default)]
pub struct DesignCommitCache {
    entries: Vec<(String, InterferenceResult)>,
    epoch: u64,
}

/// A job handle minted at a commit epoch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JobHandle {
    /// Epoch the job was minted under.
    pub epoch: u64,
}

/// Cache entry cap.
pub const MAX_CACHE_ENTRIES: usize = 64;

impl DesignCommitCache {
    /// Commit a design: bumps the epoch (CANCELS every outstanding
    /// job) and returns the new job handle.
    pub fn commit(&mut self) -> JobHandle {
        self.epoch += 1;
        JobHandle { epoch: self.epoch }
    }

    /// Look up a derivation by its FULL key. A miss is a miss —
    /// nothing is ever adapted from a near-key.
    #[must_use]
    pub fn get(&self, key: &InterferenceCacheKey) -> Option<&InterferenceResult> {
        let d = key.digest();
        self.entries.iter().find(|(k, _)| *k == d).map(|(_, r)| r)
    }

    /// Deliver a job's derivation into the cache.
    ///
    /// # Errors
    /// `bem-job-stale` (a job minted before the latest commit —
    /// commits cancel stale jobs); `bem-cache-key-mismatch` (the
    /// result's bound key differs from the insert key — the stale-key
    /// reuse hostile twin); `bem-cache-full` (cap AND cap+1).
    pub fn deliver(
        &mut self,
        job: JobHandle,
        key: &InterferenceCacheKey,
        result: InterferenceResult,
    ) -> Result<(), Refusal> {
        if job.epoch != self.epoch {
            return Err(refuse(
                "bem-job-stale",
                format!("job epoch {} at commit epoch {}", job.epoch, self.epoch),
                "a newer design commit cancelled this job; recompute",
            ));
        }
        let d = key.digest();
        if result.key_digest != d {
            return Err(refuse(
                "bem-cache-key-mismatch",
                "a result derived under one key may never be filed under another".into(),
                "re-derive under the full current key",
            ));
        }
        if self.entries.len() >= MAX_CACHE_ENTRIES {
            return Err(refuse(
                "bem-cache-full",
                format!("{} entries at the cap", self.entries.len()),
                "evict by design retirement, never silently",
            ));
        }
        self.entries.push((d, result));
        Ok(())
    }
}

/// The slider-drag SCHEMATIC preview: geometry proportions only — NO
/// solver numbers exist on this struct by construction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SchematicPreview {
    /// Gap-to-span ratio drawn.
    pub gap_over_span: f64,
    /// Stagger ratio drawn.
    pub stagger_over_chord: f64,
    /// Labeled schematic, never a result.
    pub label: &'static str,
}

/// Build the drag preview (cheap, labeled, numberless).
#[must_use]
pub fn schematic_preview(gap_over_span: f64, stagger_over_chord: f64) -> SchematicPreview {
    SchematicPreview {
        gap_over_span,
        stagger_over_chord,
        label: "schematic preview — commit to derive",
    }
}
