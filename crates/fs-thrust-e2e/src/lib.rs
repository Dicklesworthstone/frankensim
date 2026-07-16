//! fs-thrust-e2e — CertQD-Thrust: certified quality-diversity discovery of
//! self-propelling point-vortex thrusters. Layer: L6 (HELM orchestration).
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
//! - **Certificates** ([`fs_evidence`]): a point-vortex system conserves the
//!   exact LINEAR IMPULSE `I = (Σ Γᵢ yᵢ, −Σ Γᵢ xᵢ)`. A simulation is trustworthy
//!   iff it conserved that invariant to tolerance — so a converged full sim
//!   earns a `Verified` drift band whose width is the conservation slack, while
//!   a chaotic/near-singular run that leaked impulse is honestly `Estimated`.
//!   Certificates over vibes, applied to a chaotic N-body integration.
//! - **Fidelity management** ([`fs_surrogate`]): a cheap SHORT-horizon sim is a
//!   surrogate for the expensive FULL-horizon sim. A split-conformal band is
//!   calibrated on (short vs full) residuals; then `certify_or_escalate` uses the
//!   short estimate only for designs inside the calibrated validity domain when
//!   the band is decision-relevant, and ESCALATES to a full sim otherwise. The
//!   campaign spends far fewer integration steps than an all-high-fidelity sweep
//!   at equal answer quality.
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
use fs_evidence::{Color, ColorRank, IntervalOp, compose};
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

/// Color a FULL-sim drift: `Verified` (a band from the conservation slack) when
/// the impulse invariant held to `tol_rel` of its scale, else honest `Estimated`.
#[must_use]
fn full_color(res: &SimResult, impulse_scale: f64, tol_rel: f64) -> Color {
    let rel = if impulse_scale > 1e-12 {
        res.impulse_error / impulse_scale
    } else {
        res.impulse_error
    };
    if rel <= tol_rel {
        let band = res.impulse_error.max(1e-9);
        // declared-color-ok: demo drift candidate from the local impulse-error band; admitted only at a consumer's authority boundary (6pf9)
        Color::Verified {
            lo: res.drift - band,
            hi: res.drift + band,
        }
    } else {
        Color::Estimated {
            estimator: "vpm-full-nonconserving".to_string(),
            dispersion: res.impulse_error,
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
            alpha: 0.1,
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
// paired short+full sims are cheap overhead, yet spanning a useful budget/length
// hull that becomes the surrogate's declared validity domain. Designs outside
// that hull must escalate — no extrapolating a surrogate off its calibration.
fn is_calibration(d: &Design) -> bool {
    ((d.gamma - 1.0).abs() < 1e-9 || (d.gamma - 1.4).abs() < 1e-9)
        && (d.d - 0.7).abs() < 1e-9
        && ((d.l - 0.6).abs() < 1e-9 || (d.l - 1.8).abs() < 1e-9)
        && ((d.ratio - 0.4).abs() < 1e-9 || (d.ratio - 1.0).abs() < 1e-9)
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
    /// Whether the elite's drift is impulse-conservation `Verified`.
    pub verified: bool,
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
    /// Elites whose drift is `Verified` (impulse-conservation-certified).
    pub verified_elites: usize,
    /// Elites whose drift is `Estimated` (surrogate or non-conserving).
    pub estimated_elites: usize,
    /// Designs served by a full-fidelity sim.
    pub full_sims: usize,
    /// Designs served by the short surrogate.
    pub short_sims: usize,
    /// Integration steps actually spent.
    pub steps_spent: usize,
    /// Steps an all-full-fidelity sweep would have cost.
    pub steps_all_full: usize,
    /// The calibrated conformal band half-width (short-vs-full residual).
    pub band_half_width: f64,
    /// The certified drift ENVELOPE of the archive's Verified elites (via the
    /// no-laundering `compose`/Hull) — `None` if no elite is Verified.
    pub certified_envelope: Option<(f64, f64)>,
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
    // Residuals = |short-horizon drift − full-horizon drift|; the descriptor
    // hull of the calibration set is the surrogate's declared validity domain.
    let mut residuals = Vec::new();
    let (mut b_lo, mut b_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut l_lo, mut l_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    // The surrogate predicts full-horizon drift by LINEAR EXTRAPOLATION of the
    // cheap short-horizon drift (a translating cluster moves at ~constant speed);
    // the calibration residual captures the nonlinear dipole–dipole interaction
    // that makes the extrapolation imperfect.
    let extrapolate = budget.full_steps as f64 / budget.short_steps as f64;
    for d in designs.iter().filter(|d| is_calibration(d)) {
        let short = simulate_thrust(d, budget.short_steps, budget.dt, budget.core);
        let full = simulate_thrust(d, budget.full_steps, budget.dt, budget.core);
        residuals.push(short.drift * extrapolate - full.drift);
        let [bud, len] = d.descriptor();
        b_lo = b_lo.min(bud);
        b_hi = b_hi.max(bud);
        l_lo = l_lo.min(len);
        l_hi = l_hi.max(len);
    }
    let band = conformal_band(&residuals, budget.alpha);
    let in_domain = |d: &Design| -> bool {
        let [bud, len] = d.descriptor();
        bud >= b_lo && bud <= b_hi && len >= l_lo && len <= l_hi
    };

    // --- 2. Illuminate: certify-or-escalate each design, color it, archive it. ---
    let mut archive = MapElites::new(
        vec![b_lo.min(1.0), 0.5],
        vec![b_hi.max(8.0), 2.8],
        vec![budget.bins, budget.bins],
    );
    // Track each elite's color by its niche cell (for the certificate tally).
    let mut cell_color: BTreeMap<Vec<usize>, Color> = BTreeMap::new();
    // The calibration paired sims are a real (small) campaign cost.
    let (mut full_sims, mut short_sims) = (0usize, 0usize);
    let mut steps_spent = residuals.len() * (budget.short_steps + budget.full_steps);

    for d in &designs {
        let (drift, color) = match certify_or_escalate(&band, in_domain(d), budget.decision_tol) {
            Decision::UseSurrogate { band_half_width } => {
                let short = simulate_thrust(d, budget.short_steps, budget.dt, budget.core);
                short_sims += 1;
                steps_spent += budget.short_steps;
                (
                    short.drift * extrapolate,
                    Color::Estimated {
                        estimator: "vpm-short-surrogate".to_string(),
                        dispersion: band_half_width,
                    },
                )
            }
            Decision::Escalate { .. } => {
                let full = simulate_thrust(d, budget.full_steps, budget.dt, budget.core);
                full_sims += 1;
                steps_spent += budget.full_steps;
                let color = full_color(&full, d.budget(), budget.conserve_tol);
                (full.drift, color)
            }
        };
        let cell = archive.cell_of(&d.descriptor());
        // A design is an elite iff it is the best drifter in its niche.
        let solution = vec![d.gamma, d.d, d.l, d.ratio];
        if archive.add(solution, d.descriptor().to_vec(), drift) {
            cell_color.insert(cell, color);
        }
    }

    // --- 3. Certificate tally + the no-laundering color algebra. ---
    let verified_elites = cell_color
        .values()
        .filter(|c| c.rank() == ColorRank::Verified)
        .count();
    let estimated_elites = archive.num_elites() - verified_elites;
    // The certified drift envelope: HULL-compose every Verified elite band. The
    // composed rank can never outrank Verified (the no-laundering law).
    let certified_envelope = {
        let mut acc: Option<Color> = None;
        for c in cell_color
            .values()
            .filter(|c| c.rank() == ColorRank::Verified)
        {
            acc = Some(match acc {
                None => c.clone(),
                Some(prev) => compose(&prev, c, IntervalOp::Hull),
            });
        }
        match acc {
            Some(Color::Verified { lo, hi }) => Some((lo, hi)),
            _ => None,
        }
    };
    // Campaign claim rank = the weakest elite color (min rank; no laundering).
    let campaign_rank = cell_color
        .values()
        .map(fs_evidence::Color::rank)
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
        "Certified quality-diversity discovery of self-propelling four-vortex thrusters: \
         each design is simulated (fs-vpm), served by a two-fidelity certify-or-escalate \
         surrogate (fs-surrogate), impulse-conservation-certified (fs-evidence), and \
         illuminated into a MAP-Elites archive (fs-archive).",
    )
    .metric("designs_swept", designs.len() as f64, "designs")
    .metric("coverage", archive.coverage(), "fraction")
    .metric("qd_score", archive.qd_score(), "drift")
    .metric("best_drift", best_drift, "length")
    .metric("verified_elites", verified_elites as f64, "count")
    .metric("estimated_elites", estimated_elites as f64, "count")
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
            let verified = cell_color
                .get(&cell)
                .is_some_and(|c| c.rank() == ColorRank::Verified);
            AtlasEntry {
                budget: e.descriptor[0],
                length: e.descriptor[1],
                drift: e.fitness,
                verified,
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
        verified_elites,
        estimated_elites,
        full_sims,
        short_sims,
        steps_spent,
        steps_all_full,
        band_half_width: band.half_width,
        certified_envelope,
        campaign_rank,
        notebook_markdown: nb.render_markdown(),
        content_hash: nb.content_hash(),
    }
}
