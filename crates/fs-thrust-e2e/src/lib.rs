//! fs-thrust-e2e — CertQD-Thrust: SCREENED quality-diversity discovery of
//! self-propelling point-vortex thrusters. Layer: L6 (HELM orchestration).
//!
//! "Cert" names the discipline, not the claim: every drift this crate publishes
//! is `Estimated`. Nothing here is interval-certified — see [`screen_full_drift`].
//!
//! # The campaign (the raison d'être, exercised end to end)
//!
//! Classical CFD shape optimization returns ONE design after a hand-tuned run.
//! This instead ILLUMINATES the whole diverse family of vortex configurations
//! that translate themselves through an inviscid fluid — and it does so with
//! FrankenSim's actual differentiators wired together:
//!
//! - **Physics** ([`fs_vpm`]): a "thruster" is a four-vortex quadrupole (a
//!   leading `±Γ` dipole and a trailing `±Γ·ratio` dipole). Total circulation is
//!   zero by construction, so the cluster is a generalized dipole: it TRANSLATES
//!   (the vortex mean position drifts) while its linear impulse is conserved.
//!   The two dipoles interact nonlinearly — the drift is NOT analytic and must
//!   be simulated (RK4 desingularized Biot–Savart).
//! - **Screening** ([`fs_evidence`]): a point-vortex system conserves the exact
//!   LINEAR IMPULSE `I = (Σ Γᵢ yᵢ, −Σ Γᵢ xᵢ)`. A full sim that leaked impulse is
//!   not trustworthy, so the campaign SCREENS on that residual — and reports the
//!   screen for exactly what it is. Conservation of `I` bounds nothing about the
//!   mean-`x` drift (for a zero-total-circulation quadrupole the drift is not
//!   even a component of `I`), and [`fs_vpm::simulate`] is unchecked RK4 with no
//!   step-size control, so NO drift here is `Verified`: every drift is
//!   `Estimated`, and the screened ones carry their impulse residual as the
//!   estimator's dispersion. A residual is a diagnostic, not an enclosure (bead
//!   `frankensim-extreal-program-f85xj.2.30`).
//! - **Fidelity management** ([`fs_surrogate`]): a cheap SHORT-horizon sim is a
//!   surrogate for the expensive FULL-horizon sim. A split-conformal band is
//!   calibrated on (short vs full) residuals; then `certify_or_escalate` uses the
//!   short estimate only for designs inside the calibration set's PER-AXIS
//!   support — every gene the calibration pinned or varied, not just the two
//!   descriptor axes — when the band is decision-relevant, and ESCALATES to a
//!   full sim otherwise. Whether that saves integration steps is arithmetic, not
//!   a promise: the served designs must repay the paired calibration sims (see
//!   [`CampaignReport::steps_spent`]).
//! - **Illumination** ([`fs_archive`]): a MAP-Elites archive over (circulation
//!   budget × device length) keeps the best-translating configuration in every
//!   behavioral niche — the diverse Pareto atlas, not a single optimum.
//! - **Provenance** ([`fs_report`]): the campaign emits a deterministic,
//!   content-addressed lab notebook carrying the reproducing IR.
//!
//! Everything is deterministic (a fixed design grid, no RNG) — the Five
//! Explicits: units (Γ, lengths, steps), seed, budgets, versions, capabilities.
//! Ambition: `[F]` frontier synthesis; the physics is a 2-D inviscid smoke tier.

use std::collections::BTreeMap;

use fs_archive::MapElites;
use fs_evidence::{Color, ColorRank};
use fs_report::LabNotebook;
use fs_surrogate::{Decision, certify_or_escalate, conformal_band};
use fs_vpm::{VortexParticle, simulate};

/// A four-vortex thruster design (all lengths in cell units, Γ in circulation
/// units). The leading dipole sits at `x = 0`, the trailing dipole at `x = −l`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Design {
    /// Leading-dipole circulation magnitude `Γ`.
    pub gamma: f64,
    /// Transverse spacing `d` of each `±` pair.
    pub d: f64,
    /// Longitudinal separation `l` of the two dipoles (the device length).
    pub l: f64,
    /// Trailing-dipole strength as a fraction of the leading one.
    pub ratio: f64,
}

impl Design {
    /// The four point vortices (leading `±Γ`, trailing `±Γ·ratio`).
    #[must_use]
    pub fn vortices(&self) -> Vec<VortexParticle> {
        let g2 = self.gamma * self.ratio;
        vec![
            VortexParticle::new([0.0, self.d / 2.0], self.gamma),
            VortexParticle::new([0.0, -self.d / 2.0], -self.gamma),
            VortexParticle::new([-self.l, self.d / 2.0], g2),
            VortexParticle::new([-self.l, -self.d / 2.0], -g2),
        ]
    }

    /// The circulation budget `Σ|Γᵢ|` — how much vorticity the device spends.
    #[must_use]
    pub fn budget(&self) -> f64 {
        2.0 * self.gamma * (1.0 + self.ratio)
    }

    /// The MAP-Elites behavior descriptor: `[circulation budget, device length]`.
    #[must_use]
    pub fn descriptor(&self) -> [f64; 2] {
        [self.budget(), self.l]
    }
}

/// The exact 2-D point-vortex linear impulse `I = (Σ Γᵢ yᵢ, −Σ Γᵢ xᵢ)` — a
/// conserved invariant of the inviscid dynamics.
#[must_use]
fn linear_impulse(p: &[VortexParticle]) -> [f64; 2] {
    let mut i = [0.0f64, 0.0];
    for v in p {
        i[0] += v.circulation * v.pos[1];
        i[1] -= v.circulation * v.pos[0];
    }
    i
}

fn mean_x(p: &[VortexParticle]) -> f64 {
    p.iter().map(|v| v.pos[0]).sum::<f64>() / p.len() as f64
}

/// The outcome of one thruster simulation.
#[derive(Debug, Clone, Copy)]
pub struct SimResult {
    /// Net drift of the vortex mean position in `x` over the horizon.
    pub drift: f64,
    /// `‖I(T) − I(0)‖` — the impulse-conservation error (the trust signal).
    pub impulse_error: f64,
}

/// Simulate a design for `steps` RK4 steps of size `dt` with core radius `core`.
#[must_use]
pub fn simulate_thrust(design: &Design, steps: usize, dt: f64, core: f64) -> SimResult {
    let p0 = design.vortices();
    let i0 = linear_impulse(&p0);
    let x0 = mean_x(&p0);
    let pt = simulate(&p0, dt, steps, core);
    let it = linear_impulse(&pt);
    SimResult {
        drift: mean_x(&pt) - x0,
        impulse_error: ((it[0] - i0[0]).powi(2) + (it[1] - i0[1]).powi(2)).sqrt(),
    }
}

/// One full-fidelity sim's conservation-screen verdict and the honest color it
/// earns.
#[derive(Debug, Clone, PartialEq)]
struct ScreenedDrift {
    /// The impulse-conservation screen passed (a DIAGNOSTIC on `I`, never an
    /// error bound on drift).
    screened: bool,
    /// The drift's evidence color. Always `Estimated`: see [`screen_full_drift`].
    color: Color,
    /// `drift ± impulse_error` for a screened sim — the residual band, reported
    /// so the campaign can hull it. NOT an enclosure of the true drift.
    band: Option<(f64, f64)>,
}

/// Screen a FULL-sim drift on the linear-impulse residual and color it.
///
/// The screen is `‖I(T) − I(0)‖ / Σ|Γᵢ| ≤ tol_rel`. Passing it means the
/// integration did not leak the invariant it is supposed to conserve — a
/// necessary condition for trusting the run, and nothing more. It is NOT an
/// error bound on `drift`:
///
/// - `drift` is the mean-`x` displacement, which for a zero-total-circulation
///   quadrupole is not a component of `I` at all (and `I_y ≡ 0` at `t = 0`), so
///   the residual constrains a different functional;
/// - [`fs_vpm::simulate`] takes `steps` unchecked RK4 steps with no step-size
///   control, no Richardson estimate, no interval arithmetic and no outward
///   rounding, so no executable enclosure of `drift` exists to publish.
///
/// So a screened sim earns `Estimated{estimator: "vpm-full-impulse-conserving",
/// dispersion: impulse_error}` — never `Verified`. Publishing a `Verified{lo,hi}`
/// whose half-width was `impulse_error.max(1e-9)` claimed a ±1e-9 *certificate*
/// on a value produced by unchecked RK4; the floor was chosen, not derived
/// (bead `frankensim-extreal-program-f85xj.2.30`).
///
/// Every numeric input is revalidated here because [`SimResult`] is publicly
/// constructible and the campaign tolerance is public policy state.
#[must_use]
fn screen_full_drift(res: &SimResult, impulse_scale: f64, tol_rel: f64) -> ScreenedDrift {
    let malformed = ScreenedDrift {
        screened: false,
        color: Color::Estimated {
            estimator: "vpm-full-invalid-screen-input".to_string(),
            dispersion: f64::INFINITY,
        },
        band: None,
    };
    if !res.drift.is_finite()
        || !res.impulse_error.is_finite()
        || res.impulse_error < 0.0
        || !impulse_scale.is_finite()
        || impulse_scale <= 0.0
        || !tol_rel.is_finite()
        || tol_rel < 0.0
    {
        return malformed;
    }
    let rel = res.impulse_error / impulse_scale;
    if rel.is_finite() && rel <= tol_rel {
        let (lo, hi) = (res.drift - res.impulse_error, res.drift + res.impulse_error);
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return malformed;
        }
        ScreenedDrift {
            screened: true,
            color: Color::Estimated {
                estimator: "vpm-full-impulse-conserving".to_string(),
                dispersion: res.impulse_error,
            },
            band: Some((lo, hi)),
        }
    } else {
        ScreenedDrift {
            screened: false,
            color: Color::Estimated {
                estimator: "vpm-full-nonconserving".to_string(),
                dispersion: res.impulse_error,
            },
            band: None,
        }
    }
}

/// The campaign budget + physics knobs (Five Explicits: budgets + seed).
#[derive(Debug, Clone, Copy)]
pub struct CampaignBudget {
    /// Full-fidelity horizon (RK4 steps).
    pub full_steps: usize,
    /// Surrogate (short) horizon (RK4 steps).
    pub short_steps: usize,
    /// Integration step.
    pub dt: f64,
    /// Vortex core radius (desingularization).
    pub core: f64,
    /// MAP-Elites bins per descriptor axis.
    pub bins: usize,
    /// Conformal miscoverage level.
    pub alpha: f64,
    /// Decision-relevance tolerance for certify-or-escalate (drift units).
    pub decision_tol: f64,
    /// Relative impulse-conservation tolerance for the `Verified` certificate.
    pub conserve_tol: f64,
    /// Provenance seed (the sweep is deterministic; recorded for the ledger).
    pub seed: u64,
}

impl Default for CampaignBudget {
    fn default() -> CampaignBudget {
        CampaignBudget {
            full_steps: 400,
            short_steps: 60,
            dt: 0.02,
            core: 0.05,
            bins: 8,
            // Eight calibration residuals can support at most the 8th order
            // statistic: alpha must be >= 1/(8+1). Use the binary-exact 1/8
            // margin, which retains the historical max-residual band without
            // relying on an under-covering rank clamp.
            alpha: 0.125,
            // Above the typical short-vs-full residual band, so the surrogate is
            // decision-relevant for in-domain designs; a tighter tol escalates
            // more (the certify-or-escalate knob).
            decision_tol: 1.0,
            conserve_tol: 5e-2,
            seed: 1,
        }
    }
}

/// The deterministic design sweep (a regular grid — space-filling, replayable).
#[must_use]
pub fn design_grid() -> Vec<Design> {
    let gammas = [0.6, 1.0, 1.4, 1.8];
    let ds = [0.4, 0.7, 1.0, 1.3];
    let ls = [0.6, 1.2, 1.8, 2.6];
    let ratios = [0.4, 0.7, 1.0];
    let mut out = Vec::with_capacity(gammas.len() * ds.len() * ls.len() * ratios.len());
    for &gamma in &gammas {
        for &d in &ds {
            for &l in &ls {
                for &ratio in &ratios {
                    out.push(Design { gamma, d, l, ratio });
                }
            }
        }
    }
    out
}

// A SPARSE, central CALIBRATION sub-grid (8 designs): small enough that its
// paired short+full sims are cheap overhead, yet spanning a useful slab of the
// design space that becomes the surrogate's declared validity domain. Designs
// outside that slab must escalate — no extrapolating a surrogate off its
// calibration.
fn is_calibration(d: &Design) -> bool {
    ((d.gamma - 1.0).abs() < 1e-9 || (d.gamma - 1.4).abs() < 1e-9)
        && (d.d - 0.7).abs() < 1e-9
        && ((d.l - 0.6).abs() < 1e-9 || (d.l - 1.8).abs() < 1e-9)
        && ((d.ratio - 0.4).abs() < 1e-9 || (d.ratio - 1.0).abs() < 1e-9)
}

/// The calibration designs of the default sweep — the only designs whose
/// short-vs-full residual the conformal band was ever fitted on.
#[must_use]
pub fn calibration_designs() -> Vec<Design> {
    design_grid().into_iter().filter(is_calibration).collect()
}

/// The surrogate's DECLARED validity domain: the per-axis support of the
/// calibration set.
///
/// Every axis is reported, including the ones the calibration set held FIXED.
/// A fixed axis has a degenerate support (`lo == hi`), and that is the honest
/// reading: the residuals carry no information about how the short-vs-full gap
/// behaves when that gene moves. The transverse spacing `d` is the sharp case —
/// a dipole self-advects at `~Γ/(2πd)`, so `d` is a first-order driver of the
/// very drift the surrogate extrapolates, yet it is invisible to
/// [`Design::descriptor`]. Testing only the 2-D descriptor hull served 63 of 84
/// designs at spacings the calibration never saw, each handed the `d = 0.7`
/// band as its uncertainty (bead `frankensim-extreal-program-f85xj.2.29`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalibrationSupport {
    /// `[min, max]` of `gamma` over the calibration set.
    pub gamma: (f64, f64),
    /// `[min, max]` of `d` over the calibration set (degenerate by design).
    pub d: (f64, f64),
    /// `[min, max]` of `l` over the calibration set.
    pub l: (f64, f64),
    /// `[min, max]` of `ratio` over the calibration set.
    pub ratio: (f64, f64),
    /// `[min, max]` of the circulation-budget descriptor.
    pub budget: (f64, f64),
}

impl CalibrationSupport {
    /// Is `design` inside the declared validity domain on EVERY axis?
    ///
    /// Non-finite genes and non-finite support bounds are rejected: [`Design`]
    /// is publicly constructible, so this predicate revalidates rather than
    /// trusting IEEE-754 comparisons against NaN.
    #[must_use]
    pub fn contains(&self, design: &Design) -> bool {
        fn within(v: f64, (lo, hi): (f64, f64)) -> bool {
            v.is_finite() && lo.is_finite() && hi.is_finite() && v >= lo && v <= hi
        }
        within(design.gamma, self.gamma)
            && within(design.d, self.d)
            && within(design.l, self.l)
            && within(design.ratio, self.ratio)
            && within(design.budget(), self.budget)
    }
}

/// The per-axis support of the calibration set (the declared validity domain).
///
/// # Panics
/// If the calibration sub-grid is empty — a programming error in
/// [`design_grid`]/[`is_calibration`], not a reachable input.
#[must_use]
pub fn calibration_support() -> CalibrationSupport {
    let designs = calibration_designs();
    assert!(
        !designs.is_empty(),
        "the calibration sub-grid must be non-empty"
    );
    let axis = |f: fn(&Design) -> f64| -> (f64, f64) {
        designs
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| {
                (lo.min(f(d)), hi.max(f(d)))
            })
    };
    CalibrationSupport {
        gamma: axis(|d| d.gamma),
        d: axis(|d| d.d),
        l: axis(|d| d.l),
        ratio: axis(|d| d.ratio),
        budget: axis(Design::budget),
    }
}

/// One filled niche of the illuminated archive (the QD atlas).
#[derive(Debug, Clone, Copy)]
pub struct AtlasEntry {
    /// Circulation budget descriptor of the niche's elite.
    pub budget: f64,
    /// Device-length descriptor of the niche's elite.
    pub length: f64,
    /// The elite's drift (fitness).
    pub drift: f64,
    /// Whether the elite's full sim passed the impulse-conservation SCREEN.
    /// A screened drift is still `Estimated` — the screen is a diagnostic on a
    /// different functional, not an error bound on drift.
    pub conservation_screened: bool,
    /// The elite's evidence rank. Always [`ColorRank::Estimated`] here: this
    /// campaign has no executable enclosure of drift to certify.
    pub rank: ColorRank,
}

/// The full campaign report.
#[derive(Debug, Clone)]
pub struct CampaignReport {
    /// The illuminated archive as a flat atlas (one entry per filled niche).
    pub atlas: Vec<AtlasEntry>,
    /// Fraction of behavior niches filled.
    pub coverage: f64,
    /// Quality-diversity score (Σ elite drift over the archive).
    pub qd_score: f64,
    /// Number of filled niches.
    pub num_elites: usize,
    /// The single best-translating design + its drift.
    pub best: Design,
    /// Best drift found.
    pub best_drift: f64,
    /// Elites whose full sim PASSED the impulse-conservation screen. They are
    /// still `Estimated`: passing the screen is a necessary condition for
    /// trusting the integration, not a bound on the drift it produced.
    pub conservation_screened_elites: usize,
    /// Elites that did not pass the screen — surrogate-served, non-conserving,
    /// or malformed.
    pub unscreened_elites: usize,
    /// Designs served by a full-fidelity sim.
    pub full_sims: usize,
    /// Designs served by the short surrogate.
    pub short_sims: usize,
    /// Integration steps actually spent, INCLUDING the paired calibration sims.
    /// Whether this beats [`CampaignReport::steps_all_full`] is arithmetic:
    /// `short_sims·(full_steps − short_steps) > 8·(short_steps + full_steps)`
    /// must hold for the served designs to repay the calibration overhead.
    pub steps_spent: usize,
    /// Steps an all-full-fidelity sweep would have cost.
    pub steps_all_full: usize,
    /// The calibrated conformal band half-width (short-vs-full residual).
    pub band_half_width: f64,
    /// Hull of `drift ± impulse-conservation residual` over the elites that
    /// passed the conservation screen — `None` if none did.
    ///
    /// NOT an error bound on drift, and not a certificate: it is the endpoint
    /// hull of a set of residual-width bands around `Estimated` values. Read it
    /// as "where the screened elites landed, and how much invariant each one
    /// leaked", nothing stronger.
    pub conservation_screened_drift_hull: Option<(f64, f64)>,
    /// The campaign-level claim rank (weakest elite color — no laundering).
    pub campaign_rank: ColorRank,
    /// The reproducible lab notebook (Markdown).
    pub notebook_markdown: String,
    /// Content hash of the notebook (provenance).
    pub content_hash: u64,
}

/// Run the whole CertQD-Thrust campaign under `budget`.
#[must_use]
#[allow(clippy::too_many_lines)] // one coherent campaign narrative
pub fn run_campaign(budget: &CampaignBudget) -> CampaignReport {
    let designs = design_grid();

    // --- 1. Calibrate the two-fidelity surrogate on the central designs. ---
    // Residuals = |short-horizon drift − full-horizon drift|; the PER-AXIS
    // support of the calibration set is the surrogate's declared validity
    // domain (every gene, not only the two descriptor axes).
    let support = calibration_support();
    let mut residuals = Vec::new();
    // The surrogate predicts full-horizon drift by LINEAR EXTRAPOLATION of the
    // cheap short-horizon drift (a translating cluster moves at ~constant speed);
    // the calibration residual captures the nonlinear dipole–dipole interaction
    // that makes the extrapolation imperfect.
    let extrapolate = budget.full_steps as f64 / budget.short_steps as f64;
    for d in designs.iter().filter(|d| is_calibration(d)) {
        let short = simulate_thrust(d, budget.short_steps, budget.dt, budget.core);
        let full = simulate_thrust(d, budget.full_steps, budget.dt, budget.core);
        residuals.push(short.drift * extrapolate - full.drift);
    }
    let band = conformal_band(&residuals, budget.alpha);
    let (b_lo, b_hi) = support.budget;
    let in_domain = |d: &Design| -> bool { support.contains(d) };

    // --- 2. Illuminate: certify-or-escalate each design, screen it, archive it. ---
    let mut archive = MapElites::new(
        vec![b_lo.min(1.0), 0.5],
        vec![b_hi.max(8.0), 2.8],
        vec![budget.bins, budget.bins],
    );
    // Track each elite's evidence by its niche cell (for the screen tally).
    let mut cell_evidence: BTreeMap<Vec<usize>, ScreenedDrift> = BTreeMap::new();
    // The calibration paired sims are a real (small) campaign cost.
    let (mut full_sims, mut short_sims) = (0usize, 0usize);
    let mut steps_spent = residuals.len() * (budget.short_steps + budget.full_steps);

    for d in &designs {
        let (drift, evidence) = match certify_or_escalate(&band, in_domain(d), budget.decision_tol)
        {
            Decision::UseSurrogate { band_half_width } => {
                let short = simulate_thrust(d, budget.short_steps, budget.dt, budget.core);
                short_sims += 1;
                steps_spent += budget.short_steps;
                (
                    short.drift * extrapolate,
                    ScreenedDrift {
                        screened: false,
                        color: Color::Estimated {
                            estimator: "vpm-short-surrogate".to_string(),
                            dispersion: band_half_width,
                        },
                        band: None,
                    },
                )
            }
            Decision::Escalate { .. } => {
                let full = simulate_thrust(d, budget.full_steps, budget.dt, budget.core);
                full_sims += 1;
                steps_spent += budget.full_steps;
                (
                    full.drift,
                    screen_full_drift(&full, d.budget(), budget.conserve_tol),
                )
            }
        };
        let cell = archive.cell_of(&d.descriptor());
        // A design is an elite iff it is the best drifter in its niche.
        let solution = vec![d.gamma, d.d, d.l, d.ratio];
        if archive.add(solution, d.descriptor().to_vec(), drift) {
            cell_evidence.insert(cell, evidence);
        }
    }

    // --- 3. Screen tally + the no-laundering color algebra. ---
    let conservation_screened_elites = cell_evidence.values().filter(|e| e.screened).count();
    let unscreened_elites = archive.num_elites() - conservation_screened_elites;
    // The screened drift hull: the plain endpoint hull of every screened
    // elite's residual band. This is NOT `compose`d as an interval certificate:
    // the operands are Estimated, so there is no enclosure to hull and no
    // certified rank to preserve. The hull is descriptive geometry over
    // Estimated values, and it is labelled as such.
    let conservation_screened_drift_hull = cell_evidence.values().filter_map(|e| e.band).fold(
        None,
        |acc: Option<(f64, f64)>, (lo, hi)| match acc {
            None => Some((lo, hi)),
            Some((a, b)) => Some((a.min(lo), b.max(hi))),
        },
    );
    // Campaign claim rank = the weakest elite color (min rank; no laundering).
    let campaign_rank = cell_evidence
        .values()
        .map(|e| e.color.rank())
        .min()
        .unwrap_or(ColorRank::Estimated);

    let best = archive.best().expect("archive has at least one elite");
    let best_design = Design {
        gamma: best.solution[0],
        d: best.solution[1],
        l: best.solution[2],
        ratio: best.solution[3],
    };
    let best_drift = best.fitness;
    // The naive baseline: run a full-fidelity sim on every design, no surrogate.
    let steps_all_full = designs.len() * budget.full_steps;

    // --- 4. The reproducible lab notebook (fs-report). ---
    let mut nb = LabNotebook::new(
        "CertQD-Thrust campaign",
        budget.seed,
        env!("CARGO_PKG_VERSION"),
    );
    nb.prose(
        "Quality-diversity discovery of self-propelling four-vortex thrusters: each design \
         is simulated (fs-vpm), served by a two-fidelity certify-or-escalate surrogate \
         (fs-surrogate) inside the calibration set's per-axis support, screened on the \
         linear-impulse residual and colored Estimated (fs-evidence), and illuminated into \
         a MAP-Elites archive (fs-archive). The impulse screen is a diagnostic, not a drift \
         bound: no elite claims an interval certificate.",
    )
    .metric("designs_swept", designs.len() as f64, "designs")
    .metric("coverage", archive.coverage(), "fraction")
    .metric("qd_score", archive.qd_score(), "drift")
    .metric("best_drift", best_drift, "length")
    .metric(
        "conservation_screened_elites",
        conservation_screened_elites as f64,
        "count",
    )
    .metric("unscreened_elites", unscreened_elites as f64, "count")
    .metric("full_sims", full_sims as f64, "count")
    .metric("short_sims", short_sims as f64, "count")
    .metric("steps_spent", steps_spent as f64, "steps")
    .metric("steps_all_full", steps_all_full as f64, "steps")
    .metric(
        "step_savings",
        1.0 - steps_spent as f64 / steps_all_full as f64,
        "fraction",
    )
    .metric("band_half_width", band.half_width, "length")
    .step(
        "calibrate_surrogate",
        vec![
            format!("alpha={}", budget.alpha),
            format!("short={}", budget.short_steps),
            format!("full={}", budget.full_steps),
        ],
    )
    .step(
        "illuminate",
        vec![
            format!("bins={}", budget.bins),
            format!("decision_tol={}", budget.decision_tol),
        ],
    );

    let atlas: Vec<AtlasEntry> = archive
        .elites()
        .map(|e| {
            let cell = archive.cell_of(&e.descriptor);
            let evidence = cell_evidence.get(&cell);
            AtlasEntry {
                budget: e.descriptor[0],
                length: e.descriptor[1],
                drift: e.fitness,
                conservation_screened: evidence.is_some_and(|ev| ev.screened),
                rank: evidence.map_or(ColorRank::Estimated, |ev| ev.color.rank()),
            }
        })
        .collect();

    CampaignReport {
        atlas,
        coverage: archive.coverage(),
        qd_score: archive.qd_score(),
        num_elites: archive.num_elites(),
        best: best_design,
        best_drift,
        conservation_screened_elites,
        unscreened_elites,
        full_sims,
        short_sims,
        steps_spent,
        steps_all_full,
        band_half_width: band.half_width,
        conservation_screened_drift_hull,
        campaign_rank,
        notebook_markdown: nb.render_markdown(),
        content_hash: nb.content_hash(),
    }
}

#[cfg(test)]
mod screen_and_domain_tests {
    use super::*;

    fn assert_invalid_estimated(result: SimResult, impulse_scale: f64, tol_rel: f64) {
        let out = screen_full_drift(&result, impulse_scale, tol_rel);
        assert!(!out.screened, "malformed input must not pass the screen");
        assert_eq!(out.band, None, "malformed input must publish no band");
        match out.color {
            Color::Estimated {
                estimator,
                dispersion,
            } => {
                assert_eq!(estimator, "vpm-full-invalid-screen-input");
                assert_eq!(dispersion, f64::INFINITY);
            }
            other => panic!("malformed screen input must fail closed, got {other:?}"),
        }
    }

    #[test]
    fn screen_rejects_malformed_numeric_policy_and_state() {
        // Exact old-code false-clear: IEEE-754 makes a huge finite ratio less
        // than +infinity, so the comparison-only guard passed the screen.
        assert_invalid_estimated(
            SimResult {
                drift: 1.0,
                impulse_error: 1.0e300,
            },
            1.0,
            f64::INFINITY,
        );

        for result in [
            SimResult {
                drift: f64::NAN,
                impulse_error: 0.0,
            },
            SimResult {
                drift: f64::INFINITY,
                impulse_error: 0.0,
            },
            SimResult {
                drift: 1.0,
                impulse_error: f64::NAN,
            },
            SimResult {
                drift: 1.0,
                impulse_error: f64::INFINITY,
            },
            SimResult {
                drift: 1.0,
                impulse_error: -1.0,
            },
        ] {
            assert_invalid_estimated(result, 1.0, 0.1);
        }

        let valid = SimResult {
            drift: 1.0,
            impulse_error: 0.0,
        };
        for scale in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0, -0.0, 0.0] {
            assert_invalid_estimated(valid, scale, 0.1);
        }
        for tolerance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -f64::EPSILON] {
            assert_invalid_estimated(valid, 1.0, tolerance);
        }

        // Finite inputs can still overflow the proposed interval endpoints.
        assert_invalid_estimated(
            SimResult {
                drift: f64::MAX,
                impulse_error: f64::MAX,
            },
            1.0,
            f64::MAX,
        );
    }

    #[test]
    fn screen_preserves_valid_inclusive_boundary_semantics() {
        let pass = screen_full_drift(
            &SimResult {
                drift: 2.0,
                impulse_error: 0.5,
            },
            1.0,
            0.5,
        );
        assert!(pass.screened, "the inclusive boundary must pass the screen");
        assert_eq!(
            pass.color,
            Color::Estimated {
                estimator: "vpm-full-impulse-conserving".to_string(),
                dispersion: 0.5,
            }
        );
        assert_eq!(pass.band, Some((1.5, 2.5)));

        let leaked = screen_full_drift(
            &SimResult {
                drift: 2.0,
                impulse_error: 1.0,
            },
            1.0,
            0.5,
        );
        assert!(!leaked.screened);
        assert_eq!(leaked.band, None);
        assert!(matches!(
            leaked.color,
            Color::Estimated {
                dispersion: 1.0,
                ..
            }
        ));
        // A tiny positive scale must still use the relative ratio; the prior
        // absolute-error fallback made this enormous relative error pass.
        let tiny_scale = screen_full_drift(
            &SimResult {
                drift: 2.0,
                impulse_error: 0.5,
            },
            1.0e-300,
            0.5,
        );
        assert!(!tiny_scale.screened);
        assert!(matches!(
            tiny_scale.color,
            Color::Estimated {
                dispersion: 0.5,
                ..
            }
        ));
    }

    /// REGRESSION (bead `frankensim-extreal-program-f85xj.2.30`): an
    /// impulse-conservation residual bounds nothing about drift, so passing the
    /// screen must never mint an interval certificate. The witness is the
    /// bead's own: a perfectly conserving run whose drift came out of 400
    /// unchecked RK4 steps used to be published as `Verified{2.0 ± 1e-9}` — a
    /// certified enclosure whose half-width was the `.max(1e-9)` FLOOR, i.e.
    /// chosen by the code rather than derived from the integration.
    #[test]
    fn a_conserving_full_sim_is_estimated_never_verified() {
        let out = screen_full_drift(
            &SimResult {
                drift: 2.0,
                impulse_error: 0.0,
            },
            4.0,
            5e-2,
        );
        assert!(out.screened, "a zero-residual run passes the screen");
        assert_eq!(
            out.color,
            Color::Estimated {
                estimator: "vpm-full-impulse-conserving".to_string(),
                dispersion: 0.0,
            },
            "the screen is a diagnostic; it cannot mint Verified"
        );
        assert_eq!(out.color.rank(), ColorRank::Estimated);
        // The reported band is exactly the residual — no invented floor.
        assert_eq!(out.band, Some((2.0, 2.0)));
    }

    /// REGRESSION (bead `frankensim-extreal-program-f85xj.2.29`): the declared
    /// validity domain must cover EVERY axis the calibration set pinned or
    /// varied. The descriptor hull alone admitted 63 designs at transverse
    /// spacings (`d ∈ {0.4, 1.0, 1.3}`) that no calibration residual ever saw.
    #[test]
    fn the_validity_domain_rejects_every_off_calibration_axis() {
        let support = calibration_support();
        assert_eq!(support.d, (0.7, 0.7), "the calibration set fixes d = 0.7");
        assert_eq!(support.gamma, (1.0, 1.4));
        assert_eq!(support.ratio, (0.4, 1.0));
        assert_eq!(support.l, (0.6, 1.8));

        // Every calibration design is in its own domain.
        for design in calibration_designs() {
            assert!(support.contains(&design), "calibration design {design:?}");
        }

        // The bead's reaching input: in the descriptor hull (budget 2.8 ∈
        // [2.8, 5.6], length 0.6 ∈ [0.6, 1.8]) but off-calibration in d.
        let off_axis = Design {
            gamma: 1.0,
            d: 0.4,
            l: 0.6,
            ratio: 0.4,
        };
        assert!((off_axis.budget() - 2.8).abs() < 1e-12);
        assert!(
            !support.contains(&off_axis),
            "d = 0.4 was never calibrated and must escalate"
        );

        // …and the whole sweep agrees: nothing in the domain is off-calibration.
        let admitted: Vec<Design> = design_grid()
            .into_iter()
            .filter(|d| support.contains(d))
            .collect();
        assert_eq!(admitted.len(), 18, "in-domain designs of the default sweep");
        for design in &admitted {
            assert!(
                (design.d - 0.7).abs() < 1e-12,
                "off-axis served: {design:?}"
            );
            assert!(design.gamma >= 1.0 && design.gamma <= 1.4);
            assert!(design.ratio >= 0.4 && design.ratio <= 1.0);
            assert!(design.l >= 0.6 && design.l <= 1.8);
        }
    }

    /// A publicly constructed [`Design`] cannot smuggle a non-finite gene past
    /// the domain test (the predicate revalidates instead of trusting NaN
    /// comparisons).
    #[test]
    fn the_validity_domain_refuses_non_finite_genes() {
        let support = calibration_support();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for design in [
                Design {
                    gamma: bad,
                    d: 0.7,
                    l: 0.6,
                    ratio: 0.4,
                },
                Design {
                    gamma: 1.0,
                    d: bad,
                    l: 0.6,
                    ratio: 0.4,
                },
                Design {
                    gamma: 1.0,
                    d: 0.7,
                    l: bad,
                    ratio: 0.4,
                },
                Design {
                    gamma: 1.0,
                    d: 0.7,
                    l: 0.6,
                    ratio: bad,
                },
            ] {
                assert!(
                    !support.contains(&design),
                    "{design:?} must not be admitted"
                );
            }
        }
    }

    #[test]
    fn default_conformal_policy_retains_the_supported_eighth_statistic() {
        let budget = CampaignBudget::default();
        let calibration_count = design_grid()
            .iter()
            .filter(|design| is_calibration(design))
            .count();
        assert_eq!(calibration_count, 8);
        assert!(
            budget.alpha >= 1.0 / (calibration_count as f64 + 1.0),
            "split-conformal alpha must admit an order statistic"
        );
        let rank = ((1.0 - budget.alpha) * (calibration_count as f64 + 1.0)).ceil() as usize;
        assert_eq!(rank, calibration_count);
    }
}
