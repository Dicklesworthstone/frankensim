//! fs-wasm · FLAGSHIP tier (Tier V) — three live, end-to-end evidence-bearing
//! design pipelines in the browser, each a sized-down but genuinely real run
//! of a FrankenSim flagship campaign.
//!
//! - [`run_ornithoid`] — the ornithoid (flapping micro-flyer) pipeline:
//!   PARAMETERIZE → SCREEN (e-race) → VALIDATE (LBM) → CERTIFY (Lyapunov +
//!   conformal surrogate) → EXPLORE (NSGA-II Pareto + knee polish).
//! - [`run_vessel`] — the never-dribbling carafe: PARAMETERIZE → STABILITY
//!   (Orr–Sommerfeld min-max) → VALIDATE (free-surface pour) → ROBUSTIFY (CVaR
//!   over a fluid band) → RENDER (exact slab absorption transmittance).
//! - [`run_frame`] — the e-stopped seismic frame: LAYOUT (fnx-free truss LP) →
//!   SIZE (yield + Euler-buckling, catalog snap) → TIME HISTORY (nonlinear
//!   fiber hinge) → FRAGILITY (anytime CS + MLMC) → CVaR mass minimization.
//!
//! SAFETY CONTRACT (identical to the rest of the crate): `unsafe_code` is
//! forbidden, every input is bounded or refused, every fallible kernel result is folded to
//! `NaN` / an empty vector, and every documented panic precondition of the
//! composed crates is respected. Nothing here can trap — a wasm trap would kill
//! the whole page. In particular the LAYOUT stage of the frame flagship uses a
//! fnx-free vendor of the truss numerical core (the compatibility-evidence path
//! reached while `fs_truss::GroundStructure::try_grid` builds its `fnx` graph
//! reads a platform time source, which compiles but TRAPS at runtime on
//! `wasm32-unknown-unknown`); the trap-free racing path
//! (`fs_race::race_field` over an empty `KillRegistry`) never reads a clock.
//!
//! Determinism: no clocks, no entropy RNG. All stochastic paths are
//! counter-based Philox keyed by seed / logical identity.

use crate::campaigns::{rank_code, verified_bounds, with_certificate_cx};
use fs_bem::panel2d::{Airfoil2d, dcl_dalpha_adjoint, solve};
use fs_bem::wake2d::WakeSim;
use fs_dfo::moo::{Individual, NsgaParams, hypervolume, knee_point, nsga2};
use fs_eproc::GaussianMixtureCs;
use fs_exec::KillRegistry;
use fs_lbm::rheology::Rheology;
use fs_lbm::{Cell, ContactModel, FreeSurface, Grid, equilibrium};
use fs_ornith::certify::LdSurrogate;
use fs_ornith::param::OrnithCandidate;
use fs_ornith::screen::{PANELS, lift_to_drag};
use fs_qty::{Dims, QtyAny};
use fs_race::{LossSpan, RaceSettings, race_field};
use fs_rand::StreamKey;
use fs_render::volumes::{
    AbsorptionImageSettings, AbsorptionMixtureSettings, IndexedSpectralAbsorptionImageSettings,
    SpectralAbsorptionImageSettings, SpectralIllumination, VolumeGrid,
    prepare_density_absorption_mixtures, render_absorption_transmittance,
};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};
use fs_sos::lyapunov_certifies_stability;
use fs_sparse::{Coo, Csr};
use fs_truss::{LayoutCertificateLimits, LayoutCertificateProblem, PdhgSettings};
use fs_truss_e2e::{
    estimated_optimality_color, optimality_color_from_certificate, rescale_optimality_color,
};
use fs_vpm::{VortexParticle, advect};

/* ======================================================================= */
/*  Small shared helpers                                                    */
/* ======================================================================= */

/// A finite value passes through; `±∞` / `NaN` fold to `NaN` (plot-safe).
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

/// Ray-casting point-in-polygon test (closed polygon `poly`).
fn point_in_poly(poly: &[[f64; 2]], x: f64, y: f64) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (poly[i][0], poly[i][1]);
        let (xj, yj) = (poly[j][0], poly[j][1]);
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/* ======================================================================= */
/*  1 · ORNITHOID  (fs-ornith) — no run_campaign; composed inline           */
/* ======================================================================= */

/// The BEM lift coefficient of a candidate at trim.
fn bem_cl(c: &OrnithCandidate) -> f64 {
    solve(&c.section(PANELS), c.alpha).map_or(f64::NAN, |solution| solution.cl)
}

/// Central-difference `∂cl/∂α` (the trap-free 2-solve slope used everywhere
/// EXCEPT the single certified adjoint call).
fn slope_fd(foil: &Airfoil2d, alpha: f64) -> f64 {
    let h = 1.0e-4;
    let (Ok(plus), Ok(minus)) = (solve(foil, alpha + h), solve(foil, alpha - h)) else {
        return f64::NAN;
    };
    (plus.cl - minus.cl) / (2.0 * h)
}

/// A trap-free replica of `fs_ornith::certify::certify` that uses the 2-solve
/// slope in place of the expensive adjoint-assisted derivative (the battery certifies
/// adjoint ≈ FD to < 1e-6). Returns `(k, d, p11, p12, p22, cstar, roa,
/// maneuver, certified)`.
#[allow(clippy::type_complexity)]
fn cheap_certify(c: &OrnithCandidate) -> (f64, f64, f64, f64, f64, f64, f64, f64, bool) {
    let foil = c.section(PANELS);
    let dcl = slope_fd(&foil, c.alpha);
    let k = 0.4 * dcl;
    let d = 8.0 * c.thickness + 0.4 * c.flap_amp;
    let a = [[0.0, 1.0], [-k, -d]];
    let p12 = 1.0 / (2.0 * k);
    let p22 = (p12 + 0.5) / d;
    let p11 = d.mul_add(p12, k * p22);
    let p = [[p11, p12], [p12, p22]];
    let certified = k > 0.0 && d > 0.0 && lyapunov_certifies_stability(a, p);
    let (cstar, roa) = if certified {
        let det = p11.mul_add(p22, -(p12 * p12));
        let cstar = 0.35 * 0.35 * det / p22;
        (cstar, std::f64::consts::PI * cstar / det.sqrt())
    } else {
        (0.0, 0.0)
    };
    let maneuver = c.flap_amp * c.flap_freq / (k + 0.2);
    (k, d, p11, p12, p22, cstar, roa, maneuver, certified)
}

/// The inlet mass-flow band violation (0 inside `[0.9, 1.6]`, else the signed
/// distance outside).
fn inlet_violation(c: &OrnithCandidate) -> f64 {
    let m = c.inlet_mass_flow(PANELS);
    if m < 0.9 {
        0.9 - m
    } else if m > 1.6 {
        m - 1.6
    } else {
        0.0
    }
}

/// **ORNITHOID** — the flapping micro-flyer flagship, run end to end at reduced
/// size. Five stages, one flat `Vec<f64>`.
///
/// `seed` — the Philox seed for candidate sampling and the NSGA-II search.
///
/// OUTPUT LAYOUT (flat `Vec<f64>`; the viz slices by the header counts):
/// - HEADER, 26 values:
///   - `[0]` schema (= 1)
///   - `[1]` `Ns` — sampled candidate count (= 12)
///   - `[2]` `P`  — hero section node count
///   - `[3]` `Ne` — elimination events in the race
///   - `[4]` `nx` (= 64), `[5]` `ny` (= 32) — LBM grid
///   - `[6]` `Nrows` — NSGA-II Pareto front size (the "atlas" rows)
///   - `[7]` `atlas_pop`, `[8]` `atlas_gen` — NSGA-II population / generations
///   - `[9]` `adj_rel_err` — |adjoint − FD| / |adjoint| on the hero (STAGE 1)
///   - `[10]` `dcl_adjoint` — the one adjoint-assisted `∂cl/∂α`
///   - `[11]` `screen_winner_idx` — race winner candidate index
///   - `[12]` `evals_used`, `[13]` `fixed_n_equivalent` — race ledger
///   - `[14]` `candidates_eliminated` (= `Ne`)
///   - `[15]` `lbm_lift`, `[16]` `lbm_drag` — momentum-exchange body force
///   - `[17]` `panel_cl` — winner BEM lift coefficient
///   - `[18]` `lbm_steadiness` — relative velocity-field change near the end
///   - `[19]` `roa_volume_hero` — certified ROA proxy volume (hero)
///   - `[20]` `conformal_coverage` — surrogate band coverage on 12 fresh
///   - `[21]` `band_half_width` — surrogate conformal band half-width
///   - `[22]` `hypervolume`, `[23]` `knee_idx` — Pareto metrics
///   - `[24]` `polish_ld_before`, `[25]` `polish_ld_after` — knee polish
/// - BLOCK A (hero, STAGE 1): `2·P` section coords `[x,y]…`, then
///   `[inlet_x, dcl_adjoint, dcl_dthickness]`.
/// - BLOCK B (screen, STAGE 2): `Ns·4` `[ld, cl, thickness, alpha]`; then
///   `Ne·2` elimination events `[round, candidate_idx]`; then `W` (winner wake
///   particle count) followed by `W·3` `[x, y, circulation]`.
/// - BLOCK C (LBM, STAGE 3): `ny·nx` field, row-major `iy·nx+ix`, `|u|` for
///   fluid or `-1` for wall/obstacle; then the 4 CV-box corners
///   `[x0, y0, x1, y1]`.
/// - BLOCK D (certify, STAGE 4): `[k, d]`; `[p11, p12, p22]`;
///   `[cstar, roa_volume, maneuver, certified]`; then a 64-point ROA ellipse
///   `64·2` `[x, y]`; then `[surrogate_pred_ld, band_half_width]`.
/// - BLOCK E (explore, STAGE 5): `Nrows·11`
///   `[ld, roa, maneuver, inlet_viol, certified, surrogate_ld, g0..g4]`; then
///   the knee polish `7` `[ld_before, ld_after, polished_g0..g4]`.
pub fn run_ornithoid(seed: u32) -> Vec<f64> {
    // ---- sampling: 12 candidates from a Philox stream --------------------
    let ns = 12usize;
    let mut stream = StreamKey {
        seed: u64::from(seed),
        kernel: 0x0F1A,
        tile: 0,
    }
    .stream();
    let candidates: Vec<OrnithCandidate> = (0..ns)
        .map(|_| {
            let genes: Vec<f64> = (0..5).map(|_| stream.next_f64()).collect();
            OrnithCandidate::from_genes(&genes)
        })
        .collect();

    // hero = the highest-L/D champion (deterministic).
    let lds: Vec<f64> = candidates.iter().map(lift_to_drag).collect();
    let hero_idx = (0..ns)
        .max_by(|&a, &b| lds[a].total_cmp(&lds[b]))
        .unwrap_or(0);
    let hero = candidates[hero_idx];
    let hero_foil = hero.section(PANELS);

    // ---- STAGE 1: one adjoint-assisted derivative + the FD compare -------
    let Ok(dcl_adjoint) = dcl_dalpha_adjoint(&hero_foil, hero.alpha) else {
        return Vec::new();
    };
    let dcl_fd = slope_fd(&hero_foil, hero.alpha);
    let adj_rel_err = ((dcl_adjoint - dcl_fd) / dcl_adjoint.abs().max(1e-30)).abs();
    let dcl_dthickness = hero.cl_gradient(PANELS)[1];

    // ---- STAGE 2: e-race (mirrors fs_ornith::screen::screen_generation) --
    let base: Vec<f64> = lds.iter().map(|&ld| -ld).collect();
    let minb = base.iter().copied().fold(f64::INFINITY, f64::min);
    let spread = base.iter().copied().fold(f64::NEG_INFINITY, f64::max) - minb;
    let scale = 1.5 / spread.max(1e-9);
    let kills = KillRegistry::new();
    for candidate in 0..ns {
        let _ = kills.register(candidate as u64);
    }
    let mut loss = |i: usize, t: u64| {
        let mut h = ((i as u64) << 32) ^ t ^ u64::from(seed);
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        let jitter = ((h >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02;
        (base[i] - minb).mul_add(scale, jitter)
    };
    let race = match race_field(
        &mut loss,
        ns,
        // Normalized base lies in [0, 1.5]; jitter has total width 0.02.
        RaceSettings::new(LossSpan::new(1.52).expect("positive constant")),
        &kills,
    ) {
        Ok(race) => race,
        // The flat browser ABI uses an empty payload for a stage that cannot
        // make its claimed result; never turn an invalid e-process into a
        // winner or a wasm trap.
        Err(_) => return Vec::new(),
    };

    // winner flapping wake → vortex particles → advect.
    let winner = candidates[race.winner];
    let Ok(mut wake) = WakeSim::new(&winner.section(PANELS), winner.alpha, 0.08, 0.05) else {
        return Vec::new();
    };
    for _ in 0..24 {
        if wake.step().is_err() {
            return Vec::new();
        }
    }
    let mut particles: Vec<VortexParticle> = wake
        .wake()
        .iter()
        .map(|w| VortexParticle::new(w.pos, w.gamma))
        .collect();
    for _ in 0..12 {
        particles = advect(&particles, 0.1, 0.05);
    }

    // ---- STAGE 3: reduced LBM around the winner section ------------------
    let (nx, ny) = (64usize, 32usize);
    let (lbm_field, lbm_lift, lbm_drag, lbm_steadiness, cvbox) = ornithoid_lbm(&winner, nx, ny);
    let panel_cl = bem_cl(&winner);

    // ---- STAGE 4: certify hero + conformal surrogate ---------------------
    let (k, d, p11, p12, p22, cstar, roa, maneuver, certified) = cheap_certify(&hero);

    // surrogate: fit on 24 samples, coverage on 12 fresh (Philox-seeded).
    let mut sstream = StreamKey {
        seed: u64::from(seed),
        kernel: 0x0F1B,
        tile: 0,
    }
    .stream();
    let sample_c = |s: &mut fs_rand::Stream| {
        let genes: Vec<f64> = (0..5).map(|_| s.next_f64()).collect();
        OrnithCandidate::from_genes(&genes)
    };
    let train: Vec<(OrnithCandidate, f64)> = (0..24)
        .map(|_| {
            let c = sample_c(&mut sstream);
            let ld = lift_to_drag(&c);
            (c, ld)
        })
        .collect();
    let surrogate = LdSurrogate::fit(&train, 0.1);
    let fresh: Vec<OrnithCandidate> = (0..12).map(|_| sample_c(&mut sstream)).collect();
    let coverage = surrogate.coverage(&fresh);
    let band_hw = surrogate.band.half_width;
    let pred_ld = surrogate.predict(&hero);

    // ---- STAGE 5: NSGA-II Pareto explore + knee polish -------------------
    let (atlas_pop, atlas_gen) = (16usize, 8usize);
    let mut objective = |g: &[f64]| -> Vec<f64> {
        let c = OrnithCandidate::from_genes(g);
        let ld = lift_to_drag(&c);
        let (_, _, _, _, _, _, roa, man, _) = cheap_certify(&c);
        vec![-ld, -roa, -man, inlet_violation(&c)]
    };
    let front: Vec<Individual> = nsga2(
        &mut objective,
        5,
        (0.0, 1.0),
        &NsgaParams {
            pop: atlas_pop,
            generations: atlas_gen,
            eta_c: 15.0,
            eta_m: 20.0,
            p_mut: 0.2,
            seed: u64::from(seed) ^ 0x9E37_79B9,
        },
    );
    let nrows = front.len();
    let objs: Vec<Vec<f64>> = front.iter().map(|ind| ind.f.clone()).collect();
    // reference point (worst per objective + margin) for hypervolume.
    let mut reference = vec![f64::NEG_INFINITY; 4];
    for o in &objs {
        for j in 0..4 {
            reference[j] = reference[j].max(o[j]);
        }
    }
    for r in &mut reference {
        *r += 0.1 * r.abs() + 0.1;
    }
    let hv = if objs.len() >= 2 {
        hypervolume(&objs, &reference)
    } else {
        f64::NAN
    };
    let knee = if nrows >= 3 { knee_point(&objs) } else { 0 };

    // knee polish: backtracking FD-gradient ascent on L/D over the 5 genes.
    let (polish_before, polish_after, polished_genes) = if nrows > 0 {
        knee_polish(&front[knee].x)
    } else {
        (f64::NAN, f64::NAN, [f64::NAN; 5])
    };

    /* ---- assemble the flat output --------------------------------------- */
    let hero_p = hero_foil.nodes().len();
    let w = particles.len();
    let mut out: Vec<f64> = vec![
        1.0, // schema
        ns as f64,
        hero_p as f64,
        race.eliminated.len() as f64,
        nx as f64,
        ny as f64,
        nrows as f64,
        atlas_pop as f64,
        atlas_gen as f64,
        fon(adj_rel_err),
        fon(dcl_adjoint),
        race.winner as f64,
        race.evaluations_used as f64,
        race.fixed_n_equivalent as f64,
        race.eliminated.len() as f64,
        fon(lbm_lift),
        fon(lbm_drag),
        fon(panel_cl),
        fon(lbm_steadiness),
        fon(roa),
        fon(coverage),
        fon(band_hw),
        fon(hv),
        knee as f64,
        fon(polish_before),
        fon(polish_after),
    ];
    // BLOCK A: hero section + jacobian actions
    for nd in hero_foil.nodes() {
        out.push(nd[0]);
        out.push(nd[1]);
    }
    out.push(hero.inlet_x);
    out.push(fon(dcl_adjoint));
    out.push(fon(dcl_dthickness));
    // BLOCK B: candidate rows, elimination events, winner wake
    for (i, c) in candidates.iter().enumerate() {
        out.push(fon(lds[i]));
        out.push(fon(bem_cl(c)));
        out.push(c.thickness);
        out.push(c.alpha);
    }
    for &(round, idx) in &race.eliminated {
        out.push(f64::from(round));
        out.push(idx as f64);
    }
    out.push(w as f64);
    for p in &particles {
        out.push(p.pos[0]);
        out.push(p.pos[1]);
        out.push(p.circulation);
    }
    // BLOCK C: LBM field + CV box
    out.extend_from_slice(&lbm_field);
    out.extend_from_slice(&cvbox);
    // BLOCK D: certify
    out.push(fon(k));
    out.push(fon(d));
    out.push(fon(p11));
    out.push(fon(p12));
    out.push(fon(p22));
    out.push(fon(cstar));
    out.push(fon(roa));
    out.push(fon(maneuver));
    out.push(if certified { 1.0 } else { 0.0 });
    // ROA ellipse (64 points on { vᵀPv = cstar }).
    for i in 0..64 {
        let th = std::f64::consts::TAU * i as f64 / 64.0;
        let (cs, sn) = (th.cos(), th.sin());
        if certified {
            let quad = p11 * cs * cs + 2.0 * p12 * cs * sn + p22 * sn * sn;
            let s = (cstar / quad.max(1e-30)).sqrt();
            out.push(s * cs);
            out.push(s * sn);
        } else {
            out.push(0.0);
            out.push(0.0);
        }
    }
    out.push(fon(pred_ld));
    out.push(fon(band_hw));
    // BLOCK E: NSGA-II atlas rows + knee polish
    for ind in &front {
        let c = OrnithCandidate::from_genes(&ind.x);
        let (_, _, _, _, _, _, roa_i, man_i, cert_i) = cheap_certify(&c);
        out.push(fon(-ind.f[0])); // ld
        out.push(fon(roa_i));
        out.push(fon(man_i));
        out.push(fon(ind.f[3])); // inlet_viol
        out.push(if cert_i { 1.0 } else { 0.0 });
        out.push(fon(surrogate.predict(&c)));
        for &g in &ind.x {
            out.push(g);
        }
    }
    out.push(fon(polish_before));
    out.push(fon(polish_after));
    for &g in &polished_genes {
        out.push(g);
    }
    out
}

/// Backtracking FD-gradient ascent of L/D over the 5-gene box (8 outer
/// iterations); returns `(ld_before, ld_after, polished_genes)`.
fn knee_polish(start: &[f64]) -> (f64, f64, [f64; 5]) {
    let mut g = [0.0f64; 5];
    for (i, gi) in g.iter_mut().enumerate() {
        *gi = start.get(i).copied().unwrap_or(0.5).clamp(0.0, 1.0);
    }
    let ld_of = |g: &[f64; 5]| lift_to_drag(&OrnithCandidate::from_genes(g));
    let before = ld_of(&g);
    let mut current = before;
    let mut step = 0.05f64;
    let h = 1.0e-3;
    for _ in 0..8 {
        // central-difference gradient.
        let mut grad = [0.0f64; 5];
        for i in 0..5 {
            let mut gp = g;
            let mut gm = g;
            gp[i] = (gp[i] + h).min(1.0);
            gm[i] = (gm[i] - h).max(0.0);
            let denom = (gp[i] - gm[i]).max(1e-12);
            grad[i] = (ld_of(&gp) - ld_of(&gm)) / denom;
        }
        let gnorm = grad.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);
        // backtracking line search along the ascent direction.
        let mut accepted = false;
        for _ in 0..6 {
            let mut trial = g;
            for i in 0..5 {
                trial[i] = (g[i] + step * grad[i] / gnorm).clamp(0.0, 1.0);
            }
            let ld = ld_of(&trial);
            if ld > current {
                g = trial;
                current = ld;
                step *= 1.5;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    (before, current, g)
}

/// The reduced LBM validation: a body-force-driven periodic channel with the
/// winner airfoil rasterized as bounce-back walls, run with the real `fs_lbm`
/// D2Q9 kernel. Returns `(field, lift, drag, steadiness, cv_box)` where `field`
/// is the `ny·nx` row-major `|u|` map (`-1` at walls).
fn ornithoid_lbm(c: &OrnithCandidate, nx: usize, ny: usize) -> (Vec<f64>, f64, f64, f64, [f64; 4]) {
    let tau = 0.56;
    let u0 = 0.06;
    let nu = (tau - 0.5) / 3.0;
    let gx = 8.0 * nu * u0 / (ny as f64).powi(2);
    let steps = 600usize;

    let mut grid = Grid::uniform(nx, ny, tau);
    grid.periodic_x = true;
    grid.periodic_y = false;

    // rasterize the airfoil as an interior obstacle (bounce-back walls).
    let foil = c.section(PANELS);
    let chord = 22.0f64;
    let x0 = 13.0f64;
    let yc = ny as f64 / 2.0;
    let (ca, sa) = (c.alpha.cos(), c.alpha.sin());
    for y in 0..ny {
        for x in 0..nx {
            let uu = (x as f64 - x0) / chord;
            let vv = (y as f64 - yc) / chord;
            let qx = ca * uu + sa * vv;
            let qy = -sa * uu + ca * vv;
            if point_in_poly(foil.nodes(), qx, qy) {
                let i = grid.idx(x, y);
                grid.flags[i] = Cell::Wall;
            }
        }
    }
    // start from a uniform stream so the wake forms within the step budget.
    let f_stream = equilibrium(1.0, u0, 0.0);
    for i in 0..nx * ny {
        if matches!(grid.flags[i], Cell::Fluid) {
            grid.f[i] = f_stream;
        }
    }

    let velfield = |g: &Grid| -> Vec<f64> {
        let mut v = vec![-1.0f64; nx * ny];
        for (i, vi) in v.iter_mut().enumerate() {
            if matches!(g.flags[i], Cell::Fluid) {
                let m = g.moments(i);
                *vi = m.u[0].hypot(m.u[1]);
            }
        }
        v
    };

    let mut scratch: Vec<[f64; 9]> = Vec::new();
    grid.rebind_gravity([gx, 0.0], &mut scratch)
        .expect("bounded browser channel has a finite physical force frame");
    let airfoil_walls: Vec<bool> = grid.flags.iter().map(|&cell| cell == Cell::Wall).collect();
    let mut wall_impulse = [0.0; 2];
    let mut prev = velfield(&grid);
    for s in 0..steps {
        // Read the impulse from the same post-collision populations that
        // actually bounce back. The final stored populations have already
        // streamed and cannot reconstruct that exchange.
        if s + 1 == steps {
            wall_impulse = grid
                .step_with_wall_momentum(&mut scratch, &airfoil_walls)
                .wall_impulse;
        } else {
            grid.step(&mut scratch);
        }
        if s + 5 == steps {
            prev = velfield(&grid);
        }
    }
    let field = velfield(&grid);
    // steadiness: relative change of the speed field over the last 5 steps.
    let steadiness = {
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (a, b) in field.iter().zip(&prev) {
            if *a >= 0.0 && *b >= 0.0 {
                num += (a - b).abs();
                den += a.abs();
            }
        }
        if den > 1e-30 { num / den } else { f64::NAN }
    };

    let [fx, fy] = wall_impulse;
    let cvbox = [8.0, 8.0, 44.0, (ny - 8) as f64];
    (field, fy, fx, steadiness, cvbox)
}

/* ======================================================================= */
/*  2 · VESSEL  (fs-vessel) — drive the 5 stages; RENDER is the highlight    */
/* ======================================================================= */

/// The spectral min-max growth objective, replicated with a tunable
/// collocation size `n` (the real `growth_objective` hardcodes `n = 32`; the
/// grid/curve use a cheaper `n` for interactivity, the headline uses the real
/// one). Worst (max) growth rate over `stations` film-Reynolds stations.
fn growth_mm(
    prof: &fs_vessel::stability::VesselProfile,
    rate: f64,
    visc: f64,
    stations: usize,
    modes: usize,
    n: usize,
) -> f64 {
    let res = prof.film_reynolds(rate, visc, stations);
    res.iter()
        .map(|&re| {
            fs_cheb::orr_sommerfeld::growth_rates(re, 1.020_56, n, modes)
                .map(|v| v.iter().map(|c| c.re).fold(f64::NEG_INFINITY, f64::max))
                .unwrap_or(f64::NAN)
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

fn canonical_cvar_or_nan(losses: &[f64], beta: f64) -> f64 {
    fs_robust::cvar(losses, beta).unwrap_or(f64::NAN)
}

/// **VESSEL** — the never-dribbling carafe flagship, run end to end at reduced
/// size. Five stages, one flat `Vec<f64>`.
///
/// `lip_x1000` — lip width × 1000 (clamped so `lip ∈ [0.5, 3.0]`).
///
/// OUTPUT LAYOUT (flat `Vec<f64>`):
/// - HEADER, 20 values:
///   - `[0]` version (= 1)
///   - `[1]` `P` — profile samples (= 128)
///   - `[2]` `S` — growth-curve stations (= 8)
///   - `[3]` `M` — growth modes (= 4)
///   - `[4]` `nx` (= 32), `[5]` `ny` (= 20) — pour grid
///   - `[6]` `F` — pour frames (= 8)
///   - `[7]` `R` — render resolution (= 128)
///   - `[8]` `L` — CVaR lip samples (= 13)
///   - `[9]` `B` — fluid-band corners (= 9)
///   - `[10]` `spectral_growth_minmax` — worst SAMPLED Orr–Sommerfeld growth
///     rate at the input lip: the max over 4 film-Reynolds stations and the
///     first 4 modes at the fixed wavenumber `α = 1.02056`, `n = 32`
///     collocation (real `fs_vessel::stability::growth_objective`). A
///     floating-point eigensolve, `Estimated`: not an enclosure, no outward
///     rounding, no eigen-residual or collocation-convergence certificate, and
///     NOT a maximum over `α` — so it is not a min-max over the wavenumber axis
///     the name suggests.
///   - `[11]` `spectral_growth_offnom` — worst off-nominal sampled growth at
///     the input lip, on the same sampled (not enclosed) footing as `[10]`
///   - `[12]` `mass_ledger_residual` — worst pour mass-ledger drift
///   - `[13]` `poured_mass_neutral` — mass past the lip (Neutral contact)
///   - `[14]` `contact_poured_band` — |poured(Neutral) − poured(Wetting)|
///   - `[15]` `contact_dribble_band` — |dribble(Neutral) − dribble(Wetting)|
///   - `[16]` `fragments` — connected liquid components in this 2-D pour
///   - `[17]` `cvar_robust_offband` — off-band worst growth at the robust lip
///   - `[18]` `cvar_nominal_offband` — off-band worst growth at the nominal lip
///   - `[19]` reserved (0)
/// - BLOCK 1 PROFILE: `2·P` `[z, r(z)]`.
/// - BLOCK 2 GROWTH: `S` nominal per-station growth, then `S` off-nominal.
/// - BLOCK 3 POUR: `F·nx·ny` row-major mass frames (Neutral pour).
/// - BLOCK 4 CVaR: `3·L` `[lip, nominal_obj, cvar_obj]`; then `B` band-corner
///   losses at the robust lip.
/// - BLOCK 5 RENDER: `R·R` transmittance ∈ [0,1], row-major `py·R+px`.
/// - TAIL: `[robust_lip, nominal_lip, beta]`.
///
/// The render block is an authored lattice preview. Pass any captured mass
/// frame to [`render_vessel_transmittance`] for explicit SI absorption optics.
pub fn run_vessel(lip_x1000: u32) -> Vec<f64> {
    let lip = (f64::from(lip_x1000) / 1000.0).clamp(0.5, 3.0);
    let beta = 0.8f64;
    let prof = fs_vessel::stability::VesselProfile::carafe(lip);

    // sizes.
    let p = 128usize; // profile samples
    let s = 8usize; // growth stations
    let m = 4usize; // growth modes
    let (nx, ny) = (32usize, 20usize);
    let frames = 8usize;
    let render_res = 128usize;
    let l = 13usize; // cvar lip samples
    let band = fs_vessel::robust::fluid_band();
    let bcount = band.len();

    // ---- STAGE 1: profile r(z) ------------------------------------------
    let (z0, z1) = prof.radius.domain();
    let zs = linspace(z0, z1, p);
    let profile: Vec<(f64, f64)> = zs.iter().map(|&z| (z, prof.radius.eval(z))).collect();

    // ---- STAGE 2: spectral min-max + per-station growth curve -----------
    let growth_minmax = fs_vessel::stability::growth_objective(&prof, 1.0, 1.0, 4, 4); // real, n=32
    let re_nom = prof.film_reynolds(1.0, 1.0, s);
    let re_off = prof.film_reynolds(1.3, 0.6, s);
    let curve = |re: f64| {
        fs_cheb::orr_sommerfeld::growth_rates(re, 1.020_56, 24, m)
            .map(|v| v.iter().map(|c| c.re).fold(f64::NEG_INFINITY, f64::max))
            .unwrap_or(f64::NAN)
    };
    let growth_nom: Vec<f64> = re_nom.iter().map(|&re| curve(re)).collect();
    let growth_off: Vec<f64> = re_off.iter().map(|&re| curve(re)).collect();
    // off-nominal worst at the input lip (header [11]) — cheap n.
    let off_center = |g: &(f64, f64)| (g.0 - 1.0).abs() > 1e-12 || (g.1 - 1.0).abs() > 1e-12;
    let spectral_offnom = band
        .iter()
        .filter(|g| off_center(g))
        .map(|&(r, v)| growth_mm(&prof, r, v, 3, 3, 16))
        .fold(f64::NEG_INFINITY, f64::max);

    // ---- STAGE 3: free-surface pour (Neutral, capture frames) -----------
    let liquid = Rheology::Newtonian { nu: 0.016_667 };
    let Some(neutral) = pour_inline(nx, ny, frames, ContactModel::Neutral, liquid) else {
        return Vec::new();
    };
    let Some(wetting) = pour_inline(nx, ny, 0, ContactModel::Wetting, liquid) else {
        return Vec::new();
    };
    let contact_poured_band = (neutral.poured - wetting.poured).abs();
    let contact_dribble_band = (neutral.dribble - wetting.dribble).abs();

    // ---- STAGE 4: CVaR robustification over the fluid band --------------
    let lips = linspace(0.5, 3.0, l);
    let band_losses = |lp: f64| -> Vec<f64> {
        let p = fs_vessel::stability::VesselProfile::carafe(lp);
        band.iter()
            .map(|&(r, v)| growth_mm(&p, r, v, 3, 3, 16))
            .collect()
    };
    let mut cvar_rows: Vec<(f64, f64, f64)> = Vec::with_capacity(l);
    for &lp in &lips {
        let pr = fs_vessel::stability::VesselProfile::carafe(lp);
        let nominal = growth_mm(&pr, 1.0, 1.0, 3, 3, 16);
        let losses = band_losses(lp);
        let cvar = canonical_cvar_or_nan(&losses, beta);
        cvar_rows.push((lp, nominal, cvar));
    }
    // nominal lip minimizes nominal growth; robust lip minimizes CVaR.
    let nominal_lip = cvar_rows
        .iter()
        .filter(|r| r.1.is_finite())
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map_or(lip, |r| r.0);
    let robust_lip = cvar_rows
        .iter()
        .filter(|r| r.2.is_finite())
        .min_by(|a, b| a.2.total_cmp(&b.2))
        .map_or(lip, |r| r.0);
    let offband = |lp: f64| -> f64 {
        let pr = fs_vessel::stability::VesselProfile::carafe(lp);
        band.iter()
            .filter(|g| off_center(g))
            .map(|&(r, v)| growth_mm(&pr, r, v, 3, 3, 16))
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let cvar_nominal_offband = offband(nominal_lip);
    let cvar_robust_offband = offband(robust_lip);
    let robust_band_losses = band_losses(robust_lip);

    // ---- STAGE 5: authored exact slab view of the final pour mass field ----
    let transmittance = render_thin_slab(&neutral.mass_field, nx, ny, render_res);

    /* ---- assemble the flat output --------------------------------------- */
    let mut out: Vec<f64> = vec![
        1.0, // version
        p as f64,
        s as f64,
        m as f64,
        nx as f64,
        ny as f64,
        frames as f64,
        render_res as f64,
        l as f64,
        bcount as f64,
        fon(growth_minmax),
        fon(spectral_offnom),
        fon(neutral.drift),
        fon(neutral.poured),
        fon(contact_poured_band),
        fon(contact_dribble_band),
        neutral.fragments as f64,
        fon(cvar_robust_offband),
        fon(cvar_nominal_offband),
        0.0, // reserved
    ];
    // BLOCK 1 PROFILE
    for &(z, r) in &profile {
        out.push(z);
        out.push(fon(r));
    }
    // BLOCK 2 GROWTH
    for &g in &growth_nom {
        out.push(fon(g));
    }
    for &g in &growth_off {
        out.push(fon(g));
    }
    // BLOCK 3 POUR frames
    out.extend(neutral.frames.iter().map(|&v| fon(v)));
    // BLOCK 4 CVaR rows + robust band corners
    for &(lp, nom, cv) in &cvar_rows {
        out.push(lp);
        out.push(fon(nom));
        out.push(fon(cv));
    }
    for &v in &robust_band_losses {
        out.push(fon(v));
    }
    // BLOCK 5 RENDER
    out.extend(transmittance.iter().map(|&v| fon(v)));
    // TAIL
    out.push(robust_lip);
    out.push(nominal_lip);
    out.push(beta);
    out
}

/// One free-surface pour outcome (inline transcription of `fs_vessel::pour`).
struct PourOut {
    frames: Vec<f64>,
    mass_field: Vec<f64>,
    drift: f64,
    fragments: usize,
    poured: f64,
    dribble: f64,
}

/// Reduced free-surface geometry with the native vessel's prescribed-frame
/// forcing, capturing `capture` evenly-spaced mass frames. The browser retains
/// its fuller initial chamber and explicitly omits capillarity here. A native
/// rheology refusal discards the whole trial, including its captured frames.
fn pour_inline(
    nx: usize,
    ny: usize,
    capture: usize,
    contact: ContactModel,
    law: Rheology,
) -> Option<PourOut> {
    let g0 = 6.0e-4;
    let tilt_final = 0.7;
    let steps = 300u32;
    let lip_x = nx * 5 / 8;
    let lip_h = (ny * 2 / 5).min(ny - 2);
    let rig = fs_vessel::pour::PourRig {
        nx,
        ny,
        lip_x,
        lip_h,
        g0,
        surface_tension: 0.0,
        tilt_final,
        steps,
    };
    let pivot = [0.5 * (nx as f64 - 1.0), 0.5];

    let mut grid = Grid::uniform(nx, ny, 0.55);
    grid.periodic_x = false;
    grid.periodic_y = false;
    for i in 0..nx * ny {
        grid.flags[i] = Cell::Gas;
    }
    for x in 0..nx {
        let b = grid.idx(x, 0);
        grid.flags[b] = Cell::Wall;
        let t = grid.idx(x, ny - 1);
        grid.flags[t] = Cell::Wall;
    }
    for y in 0..ny {
        let l = grid.idx(0, y);
        grid.flags[l] = Cell::Wall;
        let r = grid.idx(nx - 1, y);
        grid.flags[r] = Cell::Wall;
    }
    for y in 1..=lip_h.min(ny - 2) {
        let w = grid.idx(lip_x, y);
        grid.flags[w] = Cell::Wall;
    }
    // Fill the left chamber ABOVE the weir crest so the (nearly) full carafe
    // pours over the lip under the tilt schedule (a live dribble test).
    let fill_h = (lip_h + 6).min(ny - 2);
    for y in 1..fill_h {
        for x in 1..lip_x {
            let i = grid.idx(x, y);
            grid.flags[i] = Cell::Fluid;
        }
    }
    let mut sim = FreeSurface::new(grid, 0.0, contact);
    let m0 = sim.ledger_mass().max(1e-30);

    let cell_mass = |s: &FreeSurface| -> Vec<f64> {
        let n = s.grid.nx * s.grid.ny;
        let mut v = vec![0.0f64; n];
        for (i, vi) in v.iter_mut().enumerate() {
            *vi = match s.grid.flags[i] {
                Cell::Fluid => s.grid.f[i].iter().sum(),
                Cell::Interface => s.mass[i],
                _ => 0.0,
            };
        }
        v
    };

    let capture_at: Vec<u32> = if capture > 0 {
        (0..capture)
            .map(|k| ((k as u32) * (steps - 1)) / (capture as u32 - 1).max(1))
            .collect()
    } else {
        Vec::new()
    };
    let mut frames: Vec<f64> = Vec::new();
    let mut drift = 0.0f64;
    for step in 0..steps {
        rig.prepare_step(&mut sim.grid, step, pivot).ok()?;
        fs_vessel::pour::advance_pour_liquid(&mut sim, law, step).ok()?;
        let dr = (sim.ledger_mass() - m0).abs() / m0;
        drift = drift.max(dr);
        if capture_at.contains(&step) {
            frames.extend_from_slice(&cell_mass(&sim));
        }
    }
    let mass_field = cell_mass(&sim);
    // poured mass / dribble past the lip column.
    let mut poured = 0.0f64;
    let mut dribble = 0.0f64;
    for y in 0..ny {
        for x in (lip_x + 1)..nx {
            let i = sim.grid.idx(x, y);
            let mm = match sim.grid.flags[i] {
                Cell::Fluid => sim.grid.f[i].iter().sum(),
                Cell::Interface => sim.mass[i],
                _ => 0.0,
            };
            poured += mm;
            if mm > 1e-6 {
                dribble += 1.0;
            }
        }
    }
    let fragments = sim.fragment_count();
    Some(PourOut {
        frames,
        mass_field,
        drift,
        fragments,
        poured,
        dribble,
    })
}

/// Exact absorption of a four-lattice-length slab of the final mass map.
/// This display mapping does not infer a calibrated liquid optical material.
fn render_thin_slab(mass: &[f64], nx: usize, ny: usize, res: usize) -> Vec<f64> {
    let Some(pixels) = res.checked_mul(res) else {
        return Vec::new();
    };
    if nx.checked_mul(ny) != Some(mass.len()) {
        return vec![f64::NAN; pixels];
    }
    let grid = VolumeGrid::new([nx, ny, 1], mass, [0.0, 0.0, 0.0], [1.0, 1.0, 4.0]);
    render_absorption_transmittance(
        &grid,
        AbsorptionImageSettings {
            resolution: [res; 2],
            max_pixels: pixels,
            max_grid_cells: mass.len(),
            max_cell_visits: u64::try_from(pixels).unwrap_or(u64::MAX),
        },
    )
    .unwrap_or_else(|_| vec![f64::NAN; pixels])
}

/// Re-render a captured pour mass frame with supplied SI absorption optics.
///
/// `mass_field` is the original row-major cell-average liquid mass, including
/// interface-cell volume fractions. Output is `resolution²` row-major direct
/// intensity transmittance at `wavelength_nm`. The caller supplies the mass
/// absorption coefficient at that wavelength; no liquid spectrum is guessed.
///
/// Square lattice cells have physical width `lattice_spacing_m`; the
/// independent `slab_depth_m` is the optical path length perpendicular to
/// the simulated plane. Density is `mass_field[i] * density_unit_kg_m3`.
/// Interface cells use the native view's homogeneous along-ray approximation.
/// Pixels average each original column's transmitted light over their exact
/// rectangular footprints, so downsampling preserves thin-stream coverage.
/// Surface reflection/refraction, wall occlusion and scattering are excluded.
///
/// Reuses captured frames without running the pour again or making an
/// intermediate extinction buffer. Refuses images or fields larger than
/// 1,048,576 samples or column-cell visits, and returns an empty vector for
/// any invalid input or exhausted budget. Nondivisible pixel footprints can
/// visit more columns than there are output pixels.
/// Physical parameters are validated, never clamped into a different material.
#[allow(clippy::too_many_arguments)]
pub fn render_vessel_transmittance(
    mass_field: &[f64],
    nx: u32,
    ny: u32,
    resolution: u32,
    wavelength_nm: f64,
    mass_absorption_m2_kg: f64,
    density_unit_kg_m3: f64,
    lattice_spacing_m: f64,
    slab_depth_m: f64,
) -> Vec<f64> {
    let (Ok(nx), Ok(ny), Ok(resolution)) = (
        usize::try_from(nx),
        usize::try_from(ny),
        usize::try_from(resolution),
    ) else {
        return Vec::new();
    };
    fs_vessel::pour::PourOpticalView {
        wavelength_nm,
        mass_absorption_m2_kg,
        density_unit_kg_m3,
        lattice_spacing_m,
        slab_depth_m,
    }
    .render_mass_field(
        mass_field,
        [nx, ny],
        AbsorptionImageSettings {
            resolution: [resolution; 2],
            max_pixels: 1_048_576,
            max_grid_cells: 1_048_576,
            max_cell_visits: 1_048_576,
        },
    )
    .unwrap_or_default()
}

/// Re-render a captured mass frame at explicitly supplied vacuum wavelengths.
///
/// Supply one to sixteen positive, strictly increasing `wavelengths_nm` and
/// equally many nonnegative `mass_absorption_m2_kg` coefficients at the observed
/// material state. They are monochromatic samples, with no interpolation,
/// wavelength-bin widths, illumination spectrum or detector response inferred.
///
/// Output contains `resolution² * wavelengths_nm.len()` direct transmittances,
/// indexed by `(py * resolution + px) * wavelengths_nm.len() + band`. Each
/// wavelength is averaged over the exact pixel footprint after attenuation.
/// The original mass cells, SI density conversion and explicit slab depth have
/// the same meaning as in [`render_vessel_transmittance`].
///
/// Reuses the captured frame and shares geometry/density traversal across all
/// wavelengths, without an intermediate extinction or full opacity/depth image.
/// Limits are 1,048,576 pixels, original cells and column-cell visits, and
/// 1,048,576 total spectral evaluations across all pixel/column overlaps.
/// Invalid inputs or exhausted budgets return an empty vector.
#[allow(clippy::too_many_arguments)]
pub fn render_vessel_spectral_transmittance(
    mass_field: &[f64],
    nx: u32,
    ny: u32,
    resolution: u32,
    wavelengths_nm: &[f64],
    mass_absorption_m2_kg: &[f64],
    density_unit_kg_m3: f64,
    lattice_spacing_m: f64,
    slab_depth_m: f64,
) -> Vec<f64> {
    let (Ok(nx), Ok(ny), Ok(resolution)) = (
        usize::try_from(nx),
        usize::try_from(ny),
        usize::try_from(resolution),
    ) else {
        return Vec::new();
    };
    fs_vessel::pour::PourSpectralOpticalView {
        wavelengths_nm,
        mass_absorption_m2_kg,
        density_unit_kg_m3,
        lattice_spacing_m,
        slab_depth_m,
    }
    .render_transmittance(
        mass_field,
        [nx, ny],
        SpectralAbsorptionImageSettings {
            image: AbsorptionImageSettings {
                resolution: [resolution; 2],
                max_pixels: 1_048_576,
                max_grid_cells: 1_048_576,
                max_cell_visits: 1_048_576,
            },
            max_spectral_evaluations: 1_048_576,
        },
    )
    .unwrap_or_default()
}

/// Prepare reusable mass-mixture spectra for captured vessel frames.
///
/// Constituent coefficients are component-major, wavelength-fastest, in m²/kg
/// of each constituent. Weights are mixture-major, component-fastest, in kg of
/// constituent per kg of the density field's declared reference mass. The
/// returned table is mixture-major, wavelength-fastest, in m²/kg of that same
/// reference mass. Its rows can be selected by the cell map passed to
/// [render_vessel_indexed_spectral_transmittance].
///
/// For total mixture density, callers supply mass fractions summing to one;
/// carrier-relative loadings need not sum to one. Weights are retained without
/// normalization. Absorption is additive within each homogeneous mixture;
/// spatial coverage is still filtered after attenuation by the renderer.
///
/// Limits are 64 constituents, 4,096 mixtures, sixteen common wavelength
/// samples and 1,048,576 weighted coefficient evaluations, including zero
/// weights. Invalid input, exceeded limits or unrepresentable coefficients
/// return an empty vector.
pub fn prepare_vessel_spectral_absorption_mixtures(
    wavelengths_nm: &[f64],
    component_absorption_m2_kg: &[f64],
    mixture_mass_weights: &[f64],
) -> Vec<f64> {
    prepare_density_absorption_mixtures(
        wavelengths_nm,
        component_absorption_m2_kg,
        mixture_mass_weights,
        AbsorptionMixtureSettings {
            max_components: 64,
            max_mixtures: 4_096,
            max_spectral_evaluations: 1_048_576,
        },
    )
    .unwrap_or_default()
}

/// Render a captured mass frame with an explicit per-cell spectral material map.
///
/// `cell_materials` has the original mass field's x-fastest layout and length.
/// Its u32 indices select material/state rows from `mass_absorption_m2_kg`;
/// the table is material-major with `wavelengths_nm.len()` coefficients per row.
/// Common positive, strictly increasing wavelengths are vacuum nanometres.
/// The caller supplies the palette and map; no material is inferred from mass.
///
/// SI scaling and homogeneous out-of-plane extrusion follow
/// [`render_vessel_transmittance`]. Each cell uses its own absorption spectrum
/// before exact pixel-footprint filtering. Output is pixel-major, wavelengths
/// fastest, with `resolution² * wavelengths_nm.len()` transmittances.
///
/// Admits at most 4,096 complete material rows and sixteen wavelengths, plus
/// 1,048,576 pixels, original cells and cell visits, and 1,048,576 total
/// cell/wavelength evaluations. The complete palette and all cell indices,
/// including vacuum cells, are checked before transport output allocation.
/// Invalid inputs or exhausted budgets return an empty vector.
#[allow(clippy::too_many_arguments)]
pub fn render_vessel_indexed_spectral_transmittance(
    mass_field: &[f64],
    cell_materials: &[u32],
    nx: u32,
    ny: u32,
    resolution: u32,
    wavelengths_nm: &[f64],
    mass_absorption_m2_kg: &[f64],
    density_unit_kg_m3: f64,
    lattice_spacing_m: f64,
    slab_depth_m: f64,
) -> Vec<f64> {
    let (Ok(nx), Ok(ny), Ok(resolution)) = (
        usize::try_from(nx),
        usize::try_from(ny),
        usize::try_from(resolution),
    ) else {
        return Vec::new();
    };
    fs_vessel::pour::PourIndexedSpectralOpticalView {
        wavelengths_nm,
        mass_absorption_m2_kg,
        cell_materials,
        density_unit_kg_m3,
        lattice_spacing_m,
        slab_depth_m,
    }
    .render_transmittance(
        mass_field,
        [nx, ny],
        IndexedSpectralAbsorptionImageSettings {
            spectral: SpectralAbsorptionImageSettings {
                image: AbsorptionImageSettings {
                    resolution: [resolution; 2],
                    max_pixels: 1_048_576,
                    max_grid_cells: 1_048_576,
                    max_cell_visits: 1_048_576,
                },
                max_spectral_evaluations: 1_048_576,
            },
            max_materials: 4_096,
        },
    )
    .unwrap_or_default()
}

/// Transmit supplied spectral illumination through a captured mass frame.
///
/// `background_w_m2_sr_nm` supplies radiance at the original wavelengths. With
/// `background_per_column == false`, it has one value per wavelength, uniform
/// across the frame. With `true`, it has `nx * ny * wavelengths_nm.len()` values,
/// indexed by `(y * nx + x) * wavelengths_nm.len() + band`. These are original
/// mass columns, independently of output resolution. Each source is constant
/// over its original column and lies behind the full extruded slab.
///
/// Each column's background is attenuated before pixel-footprint filtering;
/// downsampling retains the correlation between illumination and absorption.
/// Output is pixel-major, wavelength-fastest radiance in W/(m² sr nm).
/// Geometry and the absorption-only approximation follow
/// [`render_vessel_spectral_transmittance`]. No illuminant or exposure is inferred.
///
/// Limits are 1,048,576 pixels, original cells, cell visits and total spectral
/// evaluations. Invalid sources, field data, physical parameters or exhausted
/// limits return an empty vector, with no partial image.
#[allow(clippy::too_many_arguments)]
pub fn render_vessel_spectral_radiance(
    mass_field: &[f64],
    nx: u32,
    ny: u32,
    resolution: u32,
    wavelengths_nm: &[f64],
    mass_absorption_m2_kg: &[f64],
    density_unit_kg_m3: f64,
    lattice_spacing_m: f64,
    slab_depth_m: f64,
    background_w_m2_sr_nm: &[f64],
    background_per_column: bool,
) -> Vec<f64> {
    let (Ok(nx), Ok(ny), Ok(resolution)) = (
        usize::try_from(nx),
        usize::try_from(ny),
        usize::try_from(resolution),
    ) else {
        return Vec::new();
    };
    fs_vessel::pour::PourSpectralOpticalView {
        wavelengths_nm,
        mass_absorption_m2_kg,
        density_unit_kg_m3,
        lattice_spacing_m,
        slab_depth_m,
    }
    .render_radiance(
        mass_field,
        [nx, ny],
        vessel_illumination_settings(resolution),
        vessel_illumination(background_w_m2_sr_nm, background_per_column),
    )
    .unwrap_or_default()
}

/// Observe an illuminated captured frame as approximate visible XYZ triples.
///
/// Illumination layout and SI geometry match [`render_vessel_spectral_radiance`].
/// Each finite nonnegative `quadrature_weights_nm` value supplies the wavelength
/// integration width for its corresponding radiance-per-nanometre sample.
/// The shared analytic CIE observer admits only 380--780 nm and uses relative
/// `1 / integral(y_bar)` normalization. No spectral interpolation, illuminant,
/// exposure, tone mapping or RGB transform is inferred.
///
/// Output is `[X, Y, Z]` for each x-fastest pixel, with `resolution² * 3` values.
/// Direct projection allocates only the final XYZ image, retaining source and
/// absorption correlation at the original column resolution. It uses the same
/// work limits and empty-vector refusal as [`render_vessel_spectral_radiance`].
#[allow(clippy::too_many_arguments)]
pub fn render_vessel_spectral_xyz(
    mass_field: &[f64],
    nx: u32,
    ny: u32,
    resolution: u32,
    wavelengths_nm: &[f64],
    mass_absorption_m2_kg: &[f64],
    density_unit_kg_m3: f64,
    lattice_spacing_m: f64,
    slab_depth_m: f64,
    background_w_m2_sr_nm: &[f64],
    background_per_column: bool,
    quadrature_weights_nm: &[f64],
) -> Vec<f64> {
    let (Ok(nx), Ok(ny), Ok(resolution)) = (
        usize::try_from(nx),
        usize::try_from(ny),
        usize::try_from(resolution),
    ) else {
        return Vec::new();
    };
    fs_vessel::pour::PourSpectralOpticalView {
        wavelengths_nm,
        mass_absorption_m2_kg,
        density_unit_kg_m3,
        lattice_spacing_m,
        slab_depth_m,
    }
    .render_xyz(
        mass_field,
        [nx, ny],
        vessel_illumination_settings(resolution),
        vessel_illumination(background_w_m2_sr_nm, background_per_column),
        quadrature_weights_nm,
    )
    .map(Vec::into_flattened)
    .unwrap_or_default()
}

/// Transmit illumination through an explicit per-cell spectral material map.
///
/// Material-map and palette layouts match
/// [`render_vessel_indexed_spectral_transmittance`]; illumination layout and
/// radiance units match [`render_vessel_spectral_radiance`]. Materials are
/// evaluated before the illuminated columns are averaged. The original mass,
/// material map, palette and illumination remain borrowed.
///
/// Admits at most 4,096 material rows and the same complete-image limits as
/// [`render_vessel_spectral_radiance`]. Invalid input returns an empty vector.
#[allow(clippy::too_many_arguments)]
pub fn render_vessel_indexed_spectral_radiance(
    mass_field: &[f64],
    cell_materials: &[u32],
    nx: u32,
    ny: u32,
    resolution: u32,
    wavelengths_nm: &[f64],
    mass_absorption_m2_kg: &[f64],
    density_unit_kg_m3: f64,
    lattice_spacing_m: f64,
    slab_depth_m: f64,
    background_w_m2_sr_nm: &[f64],
    background_per_column: bool,
) -> Vec<f64> {
    let (Ok(nx), Ok(ny), Ok(resolution)) = (
        usize::try_from(nx),
        usize::try_from(ny),
        usize::try_from(resolution),
    ) else {
        return Vec::new();
    };
    fs_vessel::pour::PourIndexedSpectralOpticalView {
        wavelengths_nm,
        mass_absorption_m2_kg,
        cell_materials,
        density_unit_kg_m3,
        lattice_spacing_m,
        slab_depth_m,
    }
    .render_radiance(
        mass_field,
        [nx, ny],
        IndexedSpectralAbsorptionImageSettings {
            spectral: vessel_illumination_settings(resolution),
            max_materials: 4_096,
        },
        vessel_illumination(background_w_m2_sr_nm, background_per_column),
    )
    .unwrap_or_default()
}

/// Observe an illuminated material-indexed captured frame as visible XYZ.
///
/// The observer, explicit quadrature, flat XYZ output and absence of exposure
/// or display transforms match [`render_vessel_spectral_xyz`]. Material and
/// illumination layouts match [`render_vessel_indexed_spectral_radiance`].
/// Direct projection avoids intermediate spectral and opacity images.
/// Invalid input or exhausted work limits returns an empty vector.
#[allow(clippy::too_many_arguments)]
pub fn render_vessel_indexed_spectral_xyz(
    mass_field: &[f64],
    cell_materials: &[u32],
    nx: u32,
    ny: u32,
    resolution: u32,
    wavelengths_nm: &[f64],
    mass_absorption_m2_kg: &[f64],
    density_unit_kg_m3: f64,
    lattice_spacing_m: f64,
    slab_depth_m: f64,
    background_w_m2_sr_nm: &[f64],
    background_per_column: bool,
    quadrature_weights_nm: &[f64],
) -> Vec<f64> {
    let (Ok(nx), Ok(ny), Ok(resolution)) = (
        usize::try_from(nx),
        usize::try_from(ny),
        usize::try_from(resolution),
    ) else {
        return Vec::new();
    };
    fs_vessel::pour::PourIndexedSpectralOpticalView {
        wavelengths_nm,
        mass_absorption_m2_kg,
        cell_materials,
        density_unit_kg_m3,
        lattice_spacing_m,
        slab_depth_m,
    }
    .render_xyz(
        mass_field,
        [nx, ny],
        IndexedSpectralAbsorptionImageSettings {
            spectral: vessel_illumination_settings(resolution),
            max_materials: 4_096,
        },
        vessel_illumination(background_w_m2_sr_nm, background_per_column),
        quadrature_weights_nm,
    )
    .map(Vec::into_flattened)
    .unwrap_or_default()
}

fn vessel_illumination_settings(resolution: usize) -> SpectralAbsorptionImageSettings {
    SpectralAbsorptionImageSettings {
        image: AbsorptionImageSettings {
            resolution: [resolution; 2],
            max_pixels: 1_048_576,
            max_grid_cells: 1_048_576,
            max_cell_visits: 1_048_576,
        },
        max_spectral_evaluations: 1_048_576,
    }
}

fn vessel_illumination(background: &[f64], per_column: bool) -> SpectralIllumination<'_> {
    if per_column {
        SpectralIllumination::Columns(background)
    } else {
        SpectralIllumination::Uniform(background)
    }
}

/* ======================================================================= */
/*  3 · FRAME  (fs-frame) — fnx-free LAYOUT; real fs-frame stages 3–5        */
/* ======================================================================= */

/// A lean, fnx-free ground structure — the fields the truss LP + sizing read.
struct FrameGround {
    nodes: Vec<[f64; 2]>,
    members: Vec<(usize, usize)>,
    lengths: Vec<f64>,
}

/// Replicates the candidate enumeration in
/// `fs_truss::GroundStructure::try_grid` WITHOUT the fnx `Graph` (whose
/// construction-evidence path reads a platform time source — a wasm runtime trap).
fn frame_grid(nx: usize, ny: usize, w: f64, h: f64, min_len: f64, max_len: f64) -> FrameGround {
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
            let mut through = false;
            for (cc, node) in nodes.iter().enumerate() {
                if cc == a || cc == b {
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
    FrameGround {
        nodes,
        members,
        lengths,
    }
}

/// The assembled layout LP (fs-sparse only) — a transcription of
/// `fs_truss::LayoutLp`, additionally storing `dof_map` so the lean
/// `size_and_snap` refit can read it.
struct FrameLp {
    a: Csr,
    at: Csr,
    c: Vec<f64>,
    b: Vec<f64>,
    dof_map: Vec<Option<usize>>,
    norm_est: f64,
}

struct FrameReport {
    iters: usize,
    volume: f64,
    gap: f64,
    eq_residual: f64,
}

impl FrameLp {
    fn assemble(
        gs: &FrameGround,
        supported: &dyn Fn(usize, usize) -> bool,
        loads: &dyn Fn(usize) -> [f64; 2],
        sigma_y: f64,
    ) -> FrameLp {
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
        let _ = m;
        FrameLp {
            a: a_mat,
            at,
            c,
            b: b_vec,
            dof_map,
            norm_est,
        }
    }

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

    fn solve(
        &self,
        max_iters: usize,
        gap_tol: f64,
        check_every: usize,
    ) -> (Vec<f64>, Vec<f64>, FrameReport) {
        let nvar = self.c.len();
        let nrow = self.b.len();
        let mut x = vec![0.0; nvar];
        let mut y = vec![0.0; nrow];
        let step = 0.95 / self.norm_est.max(1e-30);
        let (tau, sigma) = (step, step);
        let bnorm = self.b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);
        let mut report = FrameReport {
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

/// A sized survivor member (lean transcription of `fs_truss::sizing`).
struct SizedRow {
    member_idx: usize,
    a: usize,
    b: usize,
    force: f64,
    area_yield: f64,
    area_euler: f64,
    area_catalog: f64,
}

struct SizeAudit {
    rows: Vec<SizedRow>,
    all_pass: bool,
    eq_residual: f64,
    pruned: usize,
}

/// Least-squares survivor force refit (transcription of `refit_forces`):
/// CG on the normal equations `BᵀB q = Bᵀb` over the survivor columns.
fn frame_refit(gs: &FrameGround, lp: &FrameLp, survivors: &[usize]) -> (Vec<f64>, f64) {
    let nrow = lp.b.len();
    let ns = survivors.len();
    let col = |k: usize, out: &mut Vec<f64>| {
        out.clear();
        out.resize(nrow, 0.0);
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

/// Lean transcription of `fs_truss::sizing::size_and_snap`.
fn frame_size_and_snap(
    gs: &FrameGround,
    lp: &FrameLp,
    x: &[f64],
    sigma_y: f64,
    youngs: f64,
    catalog: &[f64],
    prune_frac: f64,
) -> SizeAudit {
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
    let (q_refit, eq_residual) = frame_refit(gs, lp, &survivors);
    let mut rows = Vec::with_capacity(survivors.len());
    let mut all_pass = true;
    for (si, &k) in survivors.iter().enumerate() {
        let q = q_refit[si];
        let l = gs.lengths[k];
        let area_yield = q.abs() / sigma_y;
        let area_euler = if q < 0.0 {
            (12.0 * q.abs() * l * l / (std::f64::consts::PI.powi(2) * youngs)).sqrt()
        } else {
            0.0
        };
        let need = area_yield.max(area_euler);
        let area_catalog = catalog
            .iter()
            .copied()
            .find(|&a| a >= need)
            .unwrap_or(f64::NAN);
        let stress_ok =
            area_catalog.is_finite() && q.abs() / area_catalog <= sigma_y * (1.0 + 1e-9);
        let buckling_ok = q >= 0.0
            || (area_catalog.is_finite()
                && q.abs()
                    <= std::f64::consts::PI.powi(2) * youngs * area_catalog * area_catalog
                        / 12.0
                        / (l * l)
                        * (1.0 + 1e-9));
        all_pass &= stress_ok && buckling_ok;
        let (a, b) = gs.members[k];
        rows.push(SizedRow {
            member_idx: k,
            a,
            b,
            force: q,
            area_yield,
            area_euler,
            area_catalog,
        });
    }
    SizeAudit {
        rows,
        all_pass: all_pass && eq_residual < 1e-6,
        eq_residual,
        pruned,
    }
}

/// The `[m, kg, s, K, A, mol]` dimension of time (s) and rate (1/s).
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// Build the seismic ensemble (Kanai–Tajimi), sized for interactivity.
fn frame_ensemble(seed: u32, members: u32, duration: f64, dt: f64) -> StochasticEnsemble {
    StochasticEnsemble {
        name: "flagship-kt".to_string(),
        seed: u64::from(seed),
        members,
        duration: QtyAny::new(duration, TIME),
        dt: QtyAny::new(dt, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0: 0.01,
            omega_g: QtyAny::new(12.5, RATE),
            zeta_g: 0.6,
        },
    }
}

/// Native checked history, including the explicit initial row for the wire
/// format's one-to-one alignment with the original acceleration samples.
fn story_history(
    params: &fs_frame::StoryParams,
    ag: &[f64],
    dt: f64,
) -> Result<(Vec<f64>, Vec<f64>, f64), fs_frame::HistoryError> {
    let mut frame = fs_frame::StoryFrame::new(*params);
    let k0 = frame.initial_stiffness();
    let initial_x = frame.x;
    let initial_shear = frame.restoring(initial_x).0;
    let response = frame.run_from_time_zero(
        fs_frame::GroundMotion::new(ag, dt),
        fs_frame::HistoryLimits::new(ag.len(), 30, 1e-12, 1.0),
    )?;
    let (mut xs, mut shears) = (response.displacement_m, response.restoring_shear_n);
    xs.try_reserve(1)
        .map_err(|_| fs_frame::HistoryError::AllocationFailed { samples: ag.len() })?;
    shears
        .try_reserve(1)
        .map_err(|_| fs_frame::HistoryError::AllocationFailed { samples: ag.len() })?;
    xs.insert(0, initial_x);
    shears.insert(0, initial_shear);
    Ok((xs, shears, k0))
}

/// **FRAME** — the e-stopped seismic frame flagship, run end to end at reduced
/// size. Five stages, one flat `Vec<f64>` with a self-describing offset header.
///
/// `seed` — the ensemble Philox seed.
///
/// OUTPUT LAYOUT (flat `Vec<f64>`):
/// - HEADER, 12 values (offsets are f64 indices into this array):
///   `[0]` MAGIC (= 2), `[1]` version (= 2), `[2]` seed, `[3]` off_layout,
///   `[4]` off_sizing, `[5]` off_history, `[6]` off_fragility, `[7]` off_cvar,
///   `[8]` total_len, `[9..12]` reserved (0).
/// - LAYOUT block (@ off_layout): `gap, eq_residual, volume_phys, iters,
///   solver_converged, Nn, M, load_node_idx`; then `2·Nn` node coords
///   `[x,y]…`; then `M·4` `[na, nb, force_q, is_survivor]`; then
///   `optimality_rank, verified_flag, optimum_lo_phys, optimum_hi_phys`.
/// - SIZING block (@ off_sizing): `all_pass, eq_residual_postprune, pruned, Ms`;
///   then `Ms·7` `[member_idx, na, nb, force, area_yield, area_euler,
///   area_catalog]`.
/// - HISTORY block (@ off_history): `peak_drift, dt, Ns, k0, Vy`; then `Ns` ag,
///   `Ns` x, `Ns` V.
/// - FRAGILITY block (@ off_fragility): `p_hat, radius, members_used,
///   stopped_early, alpha, confidence, exceedances, drift_limit, margin,
///   mlmc_estimate, mlmc_levels, Nc`; then `Nc·3` `[member_idx,
///   cs_center_running, cs_radius_running]`.
/// - CVaR block (@ off_cvar): `scale_star, scale_snapped, cvar_star,
///   cvar_snapped, limit, mass, iters, beta, Ncv, catalog_len`; then
///   `catalog_len` catalog scales; then `Ncv·2` `[scale, cvar]`.
pub fn run_frame(seed: u32) -> Vec<f64> {
    // ---- STAGE 1: LAYOUT (fnx-free truss LP) ----------------------------
    let (nx, ny, w, h) = (4usize, 2usize, 6.0f64, 3.0f64);
    let min_len = 0.1;
    let max_len = (w * w + h * h).sqrt();
    let gs = frame_grid(nx, ny, w, h, min_len, max_len);
    let nn = gs.nodes.len();
    let m = gs.members.len();
    let supported = |node: usize, _comp: usize| gs.nodes[node][0] < 1e-9;
    let load_node = (0..nn)
        .max_by(|&a, &b| {
            (gs.nodes[a][0] + gs.nodes[a][1]).total_cmp(&(gs.nodes[b][0] + gs.nodes[b][1]))
        })
        .unwrap_or(0);
    let loads = |node: usize| {
        if node == load_node {
            [0.0, -1.0]
        } else {
            [0.0, 0.0]
        }
    };
    let sigma_y = 250.0e6;
    let lp = FrameLp::assemble(&gs, &supported, &loads, 1.0);
    let settings = PdhgSettings {
        max_iters: 60_000,
        gap_tol: 1e-4,
        check_every: 500,
    };
    let (x, y, report) = lp.solve(settings.max_iters, settings.gap_tol, settings.check_every);
    let normalized_optimality = match LayoutCertificateProblem::try_new(&lp.a, &lp.c, &lp.b) {
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
    let optimality_color = rescale_optimality_color(&normalized_optimality, sigma_y);
    let force = |k: usize| x[k] - x[m + k];
    let max_force = (0..m).map(|k| force(k).abs()).fold(0.0, f64::max);
    let active_tol = 1e-3 * max_force.max(1e-12);
    let solver_converged = report.gap < 1e-3 && report.eq_residual < 1e-3;

    // ---- STAGE 2: SIZING (yield + Euler-buckling, catalog snap) ----------
    let youngs = 200.0e9;
    let area_catalog: Vec<f64> = (1..=20).map(|k| 2.0e-4 * f64::from(k)).collect();
    let audit = frame_size_and_snap(&gs, &lp, &x, sigma_y, youngs, &area_catalog, 1e-3);

    // ---- STAGES 3–5 need the ensemble; validate it never traps ----------
    let (members, duration, dt) = (48u32, 8.0f64, 0.02f64);
    let ensemble = frame_ensemble(seed, members, duration, dt);
    let base = fs_frame::StoryParams::default();
    // guard: a realization error would panic the real stage APIs (trap on
    // wasm) — probe once and bail to a NaN body if the ensemble is malformed.
    let Ok(rep0) = ensemble.realize(0) else {
        return frame_nan_body(seed);
    };

    // ---- STAGE 3: TIME HISTORY (member 0, faithful hysteresis) ----------
    let ag = rep0.values.clone();
    let Ok((xs, vs, k0)) = story_history(&base, &ag, dt) else {
        return frame_nan_body(seed);
    };
    let ns_hist = ag.len();
    let peak = fs_frame::peak_drift(&xs, base.h);
    let vy = vs.iter().fold(0.0f64, |m, &v| m.max(v.abs()));

    // choose a drift_limit at the median member peak drift so p_hat ∈ (0,1)
    // (a thinned pre-pass — every third member — keeps this cheap).
    let history_limits = fs_frame::HistoryLimits::new(ag.len(), 30, 1e-12, 1.0);
    let mut peaks = Vec::new();
    for member in (0..members).step_by(3) {
        let Ok(record) = ensemble.realize(member) else {
            return frame_nan_body(seed);
        };
        let Ok(response) = fs_frame::StoryFrame::new(base).run_from_time_zero(
            fs_frame::GroundMotion::new(&record.values, dt),
            history_limits,
        ) else {
            return frame_nan_body(seed);
        };
        peaks.push(response.peak_abs_displacement_m / base.h);
    }
    peaks.sort_by(f64::total_cmp);
    let drift_limit = if peaks.is_empty() {
        base.h * 0.01
    } else {
        peaks[peaks.len() / 2]
    };

    // ---- STAGE 4: FRAGILITY (real e-stopped anytime CS + MLMC) ----------
    let (alpha, margin) = (0.1f64, 0.30f64);
    let Ok(frag) = with_certificate_cx(|cx| {
        fs_frame::try_e_stopped_fragility(
            &ensemble,
            base,
            fs_frame::FragilityControls {
                drift_limit,
                alpha,
                margin,
                realization_budget: fs_scenario::RealizationBudget::default(),
                history_limits,
            },
            cx,
        )
    }) else {
        return frame_nan_body(seed);
    };
    // running CS curve: replay members_used with the same GaussianMixtureCs.
    let mut cs = GaussianMixtureCs::new(0.5, 8.0, alpha);
    let mut cs_curve: Vec<(usize, f64, f64)> = Vec::new();
    for member in 0..frag.members_used {
        let Ok(record) = ensemble.realize(member) else {
            return frame_nan_body(seed);
        };
        let Ok(response) = fs_frame::StoryFrame::new(base).run_from_time_zero(
            fs_frame::GroundMotion::new(&record.values, dt),
            history_limits,
        ) else {
            return frame_nan_body(seed);
        };
        let pd = response.peak_abs_displacement_m / base.h;
        cs.observe(if pd > drift_limit { 1.0 } else { 0.0 });
        if let Some((center, radius)) = cs.interval() {
            cs_curve.push((member as usize, center, radius));
        }
    }

    // ---- STAGE 5: CVaR-vs-scale curve + mass-min design -----------------
    // The CVaR bisection re-evaluates the whole ensemble at every scale (the
    // 7-point curve plus ~12 losses() calls inside cvar_mass_min dominate the
    // whole campaign's cost), so this sub-study runs on a reduced ensemble
    // (the first 12 members at 4 s) for interactivity — enough samples for a
    // non-degenerate CVaR_0.9 while cutting the stage's work ~2×. It is fully
    // isolated from the headline: the certified fragility/history above keep
    // the full 48-member / 8 s ensemble and its p̂ / e-stop are unaffected.
    let ens_cvar = frame_ensemble(seed, 12, 4.0, dt);
    let cvar_catalog = [0.5f64, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0];
    let beta = 0.9f64;
    let mut cvar_curve: Vec<(f64, f64)> = Vec::with_capacity(cvar_catalog.len());
    for &s in &cvar_catalog {
        match fs_frame::try_ensemble_cvar(&ens_cvar, base, s, beta) {
            Ok(c) if c.is_finite() => cvar_curve.push((s, c)),
            _ => return frame_nan_body(seed),
        }
    }
    // feasible limit: the CVaR at scale 1.0 (guarantees cvar_at(4.0) ≤ limit).
    let cvar_at_1 = cvar_curve
        .iter()
        .find(|(s, _)| (*s - 1.0).abs() < 1e-9)
        .map_or(f64::NAN, |&(_, c)| c);
    let cvar_at_4 = cvar_curve
        .iter()
        .find(|(s, _)| (*s - 4.0).abs() < 1e-9)
        .map_or(f64::NAN, |&(_, c)| c);
    let limit = cvar_at_1 * 1.0001;
    let design = if limit.is_finite() && cvar_at_4.is_finite() && cvar_at_4 <= limit {
        fs_frame::try_cvar_mass_min(&ens_cvar, base, beta, limit, &cvar_catalog).ok()
    } else {
        None
    };

    /* ---- assemble the flat output (build blocks, then the offset header) */
    let mut layout: Vec<f64> = Vec::new();
    layout.push(fon(report.gap));
    layout.push(fon(report.eq_residual));
    layout.push(fon(report.volume / sigma_y));
    layout.push(report.iters as f64);
    layout.push(if solver_converged { 1.0 } else { 0.0 });
    layout.push(nn as f64);
    layout.push(m as f64);
    layout.push(load_node as f64);
    for p in &gs.nodes {
        layout.push(p[0]);
        layout.push(p[1]);
    }
    for k in 0..m {
        let (a, b) = gs.members[k];
        layout.push(a as f64);
        layout.push(b as f64);
        layout.push(fon(force(k)));
        layout.push(if force(k).abs() > active_tol {
            1.0
        } else {
            0.0
        });
    }
    let (verified, optimum_lo, optimum_hi) = verified_bounds(&optimality_color);
    layout.push(rank_code(optimality_color.rank()));
    layout.push(verified);
    layout.push(optimum_lo);
    layout.push(optimum_hi);

    let mut sizing: Vec<f64> = Vec::new();
    sizing.push(if audit.all_pass { 1.0 } else { 0.0 });
    sizing.push(fon(audit.eq_residual));
    sizing.push(audit.pruned as f64);
    sizing.push(audit.rows.len() as f64);
    for r in &audit.rows {
        sizing.push(r.member_idx as f64);
        sizing.push(r.a as f64);
        sizing.push(r.b as f64);
        sizing.push(fon(r.force));
        sizing.push(fon(r.area_yield));
        sizing.push(fon(r.area_euler));
        sizing.push(fon(r.area_catalog));
    }

    let mut history: Vec<f64> = vec![fon(peak), dt, ns_hist as f64, fon(k0), fon(vy)];
    history.extend(ag.iter().map(|&v| fon(v)));
    history.extend(xs.iter().map(|&v| fon(v)));
    history.extend(vs.iter().map(|&v| fon(v)));

    let mut fragility: Vec<f64> = Vec::new();
    fragility.push(fon(frag.p_hat));
    fragility.push(fon(frag.radius));
    fragility.push(f64::from(frag.members_used));
    fragility.push(if frag.stopped_early { 1.0 } else { 0.0 });
    fragility.push(alpha);
    fragility.push(1.0 - alpha);
    fragility.push(f64::from(frag.exceedances));
    fragility.push(fon(drift_limit));
    fragility.push(margin);
    fragility.push(fon(frag.mlmc.estimate));
    fragility.push(frag.mlmc.levels.len() as f64);
    fragility.push(cs_curve.len() as f64);
    for &(mi, center, radius) in &cs_curve {
        fragility.push(mi as f64);
        fragility.push(fon(center));
        fragility.push(fon(radius));
    }

    let mut cvar: Vec<f64> = Vec::new();
    if let Some(d) = &design {
        cvar.push(fon(d.scale_star));
        cvar.push(fon(d.scale_snapped));
        cvar.push(fon(d.cvar_star));
        cvar.push(fon(d.cvar_snapped));
        cvar.push(fon(limit));
        cvar.push(fon(d.mass));
        cvar.push(f64::from(d.iters));
        cvar.push(beta);
    } else {
        for _ in 0..8 {
            cvar.push(f64::NAN);
        }
        cvar[4] = fon(limit);
        cvar[7] = beta;
    }
    cvar.push(cvar_curve.len() as f64);
    cvar.push(cvar_catalog.len() as f64);
    for &c in &cvar_catalog {
        cvar.push(c);
    }
    for &(s, c) in &cvar_curve {
        cvar.push(s);
        cvar.push(fon(c));
    }

    // offset header.
    let header_len = 12usize;
    let off_layout = header_len;
    let off_sizing = off_layout + layout.len();
    let off_history = off_sizing + sizing.len();
    let off_fragility = off_history + history.len();
    let off_cvar = off_fragility + fragility.len();
    let total_len = off_cvar + cvar.len();

    let mut out: Vec<f64> = Vec::with_capacity(total_len);
    out.push(2.0); // MAGIC
    out.push(2.0); // version
    out.push(f64::from(seed));
    out.push(off_layout as f64);
    out.push(off_sizing as f64);
    out.push(off_history as f64);
    out.push(off_fragility as f64);
    out.push(off_cvar as f64);
    out.push(total_len as f64);
    out.push(0.0);
    out.push(0.0);
    out.push(0.0);
    out.extend_from_slice(&layout);
    out.extend_from_slice(&sizing);
    out.extend_from_slice(&history);
    out.extend_from_slice(&fragility);
    out.extend_from_slice(&cvar);
    out
}

/// A minimal, non-trapping frame body used only if the ensemble is malformed.
fn frame_nan_body(seed: u32) -> Vec<f64> {
    let mut out = vec![2.0, 2.0, f64::from(seed)];
    out.extend(std::iter::repeat_n(f64::NAN, 9));
    out
}

/* ======================================================================= */
/*  The JavaScript boundary (wasm32 only)                                   */
/* ======================================================================= */

#[cfg(target_arch = "wasm32")]
mod wasm {
    use wasm_bindgen::prelude::*;

    #[wasm_bindgen]
    pub fn run_ornithoid(seed: u32) -> Vec<f64> {
        super::run_ornithoid(seed)
    }

    #[wasm_bindgen]
    pub fn run_vessel(lip_x1000: u32) -> Vec<f64> {
        super::run_vessel(lip_x1000)
    }


    /// Monochromatic SI absorption of an existing captured mass frame.
    #[wasm_bindgen]
    #[allow(clippy::too_many_arguments)]
    pub fn render_vessel_transmittance(
        mass_field: &[f64],
        nx: u32,
        ny: u32,
        resolution: u32,
        wavelength_nm: f64,
        mass_absorption_m2_kg: f64,
        density_unit_kg_m3: f64,
        lattice_spacing_m: f64,
        slab_depth_m: f64,
    ) -> Vec<f64> {
        super::render_vessel_transmittance(
            mass_field,
            nx,
            ny,
            resolution,
            wavelength_nm,
            mass_absorption_m2_kg,
            density_unit_kg_m3,
            lattice_spacing_m,
            slab_depth_m,
        )
    }

    /// Spectral SI absorption of an existing captured mass frame.
    #[wasm_bindgen]
    #[allow(clippy::too_many_arguments)]
    pub fn render_vessel_spectral_transmittance(
        mass_field: &[f64],
        nx: u32,
        ny: u32,
        resolution: u32,
        wavelengths_nm: &[f64],
        mass_absorption_m2_kg: &[f64],
        density_unit_kg_m3: f64,
        lattice_spacing_m: f64,
        slab_depth_m: f64,
    ) -> Vec<f64> {
        super::render_vessel_spectral_transmittance(
            mass_field,
            nx,
            ny,
            resolution,
            wavelengths_nm,
            mass_absorption_m2_kg,
            density_unit_kg_m3,
            lattice_spacing_m,
            slab_depth_m,
        )
    }

    /// Prepare mass-mixture coefficient rows for the indexed vessel renderer.
    #[wasm_bindgen]
    pub fn prepare_vessel_spectral_absorption_mixtures(
        wavelengths_nm: &[f64],
        component_absorption_m2_kg: &[f64],
        mixture_mass_weights: &[f64],
    ) -> Vec<f64> {
        super::prepare_vessel_spectral_absorption_mixtures(
            wavelengths_nm,
            component_absorption_m2_kg,
            mixture_mass_weights,
        )
    }

    /// Spectral SI absorption using explicit original cell material indices.
    #[wasm_bindgen]
    #[allow(clippy::too_many_arguments)]
    pub fn render_vessel_indexed_spectral_transmittance(
        mass_field: &[f64],
        cell_materials: &[u32],
        nx: u32,
        ny: u32,
        resolution: u32,
        wavelengths_nm: &[f64],
        mass_absorption_m2_kg: &[f64],
        density_unit_kg_m3: f64,
        lattice_spacing_m: f64,
        slab_depth_m: f64,
    ) -> Vec<f64> {
        super::render_vessel_indexed_spectral_transmittance(
            mass_field,
            cell_materials,
            nx,
            ny,
            resolution,
            wavelengths_nm,
            mass_absorption_m2_kg,
            density_unit_kg_m3,
            lattice_spacing_m,
            slab_depth_m,
        )
    }

    /// Physical radiance through captured mass, with explicit illumination layout.
    #[wasm_bindgen]
    #[allow(clippy::too_many_arguments)]
    pub fn render_vessel_spectral_radiance(
        mass_field: &[f64],
        nx: u32,
        ny: u32,
        resolution: u32,
        wavelengths_nm: &[f64],
        mass_absorption_m2_kg: &[f64],
        density_unit_kg_m3: f64,
        lattice_spacing_m: f64,
        slab_depth_m: f64,
        background_w_m2_sr_nm: &[f64],
        background_per_column: bool,
    ) -> Vec<f64> {
        super::render_vessel_spectral_radiance(
            mass_field,
            nx,
            ny,
            resolution,
            wavelengths_nm,
            mass_absorption_m2_kg,
            density_unit_kg_m3,
            lattice_spacing_m,
            slab_depth_m,
            background_w_m2_sr_nm,
            background_per_column,
        )
    }

    /// Approximate visible XYZ with explicit source radiance and quadrature.
    #[wasm_bindgen]
    #[allow(clippy::too_many_arguments)]
    pub fn render_vessel_spectral_xyz(
        mass_field: &[f64],
        nx: u32,
        ny: u32,
        resolution: u32,
        wavelengths_nm: &[f64],
        mass_absorption_m2_kg: &[f64],
        density_unit_kg_m3: f64,
        lattice_spacing_m: f64,
        slab_depth_m: f64,
        background_w_m2_sr_nm: &[f64],
        background_per_column: bool,
        quadrature_weights_nm: &[f64],
    ) -> Vec<f64> {
        super::render_vessel_spectral_xyz(
            mass_field,
            nx,
            ny,
            resolution,
            wavelengths_nm,
            mass_absorption_m2_kg,
            density_unit_kg_m3,
            lattice_spacing_m,
            slab_depth_m,
            background_w_m2_sr_nm,
            background_per_column,
            quadrature_weights_nm,
        )
    }

    /// Physical radiance through captured mass and an explicit material map.
    #[wasm_bindgen]
    #[allow(clippy::too_many_arguments)]
    pub fn render_vessel_indexed_spectral_radiance(
        mass_field: &[f64],
        cell_materials: &[u32],
        nx: u32,
        ny: u32,
        resolution: u32,
        wavelengths_nm: &[f64],
        mass_absorption_m2_kg: &[f64],
        density_unit_kg_m3: f64,
        lattice_spacing_m: f64,
        slab_depth_m: f64,
        background_w_m2_sr_nm: &[f64],
        background_per_column: bool,
    ) -> Vec<f64> {
        super::render_vessel_indexed_spectral_radiance(
            mass_field,
            cell_materials,
            nx,
            ny,
            resolution,
            wavelengths_nm,
            mass_absorption_m2_kg,
            density_unit_kg_m3,
            lattice_spacing_m,
            slab_depth_m,
            background_w_m2_sr_nm,
            background_per_column,
        )
    }

    /// Approximate visible XYZ through captured mass and an explicit material map.
    #[wasm_bindgen]
    #[allow(clippy::too_many_arguments)]
    pub fn render_vessel_indexed_spectral_xyz(
        mass_field: &[f64],
        cell_materials: &[u32],
        nx: u32,
        ny: u32,
        resolution: u32,
        wavelengths_nm: &[f64],
        mass_absorption_m2_kg: &[f64],
        density_unit_kg_m3: f64,
        lattice_spacing_m: f64,
        slab_depth_m: f64,
        background_w_m2_sr_nm: &[f64],
        background_per_column: bool,
        quadrature_weights_nm: &[f64],
    ) -> Vec<f64> {
        super::render_vessel_indexed_spectral_xyz(
            mass_field,
            cell_materials,
            nx,
            ny,
            resolution,
            wavelengths_nm,
            mass_absorption_m2_kg,
            density_unit_kg_m3,
            lattice_spacing_m,
            slab_depth_m,
            background_w_m2_sr_nm,
            background_per_column,
            quadrature_weights_nm,
        )
    }

    #[wasm_bindgen]
    pub fn run_frame(seed: u32) -> Vec<f64> {
        super::run_frame(seed)
    }
}

/* ======================================================================= */
/*  Regression tests — the certified headline numbers must reproduce.       */
/* ======================================================================= */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pour_slab_uses_exact_absorption_and_refuses_nonphysical_coefficients() {
        let mass = [0.0, 1e-9, 0.5, 1.0];
        let image = render_thin_slab(&mass, 2, 2, 2);
        for (actual, coefficient) in image.iter().zip(mass) {
            assert!((actual - (-coefficient * 4.0).exp()).abs() < 2e-15);
        }
        assert!(image[1] < 1.0);
        assert!(render_thin_slab(&[-0.1], 1, 1, 1)[0].is_nan());
        assert!(render_thin_slab(&[f64::NAN], 1, 1, 1)[0].is_nan());
    }

    #[test]
    fn browser_physical_pour_view_matches_native_optics_and_source_wavelength() {
        let mass = [0.0, 0.25, 0.5, 1.0];
        let optics = fs_vessel::pour::PourOpticalView {
            wavelength_nm: 550.0,
            mass_absorption_m2_kg: 0.4,
            density_unit_kg_m3: 1000.0,
            lattice_spacing_m: 1e-3,
            slab_depth_m: 2e-3,
        };
        let native = optics
            .render_mass_field(
                &mass,
                [2, 2],
                AbsorptionImageSettings {
                    resolution: [2; 2],
                    max_pixels: 4,
                    max_grid_cells: 4,
                    max_cell_visits: 4,
                },
            )
            .unwrap();
        let browser = render_vessel_transmittance(
            &mass, 2, 2, 2, 550.0, 0.4, 1000.0, 1e-3, 2e-3,
        );
        assert_eq!(browser, native);
        for (&actual, &density) in browser.iter().zip(&mass) {
            assert!((actual - (-density * 0.8).exp()).abs() < 8e-15);
        }
        let red = render_vessel_transmittance(
            &mass, 2, 2, 2, 700.0, 0.1, 1000.0, 1e-3, 2e-3,
        );
        assert!((red[3] - (-0.2_f64).exp()).abs() < 8e-15);
    }

    #[test]
    fn browser_physical_pour_keeps_subpixel_liquid_coverage() {
        let mass = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        assert_eq!(
            render_vessel_transmittance(&mass, 4, 2, 1, 550.0, 1000.0, 1.0, 1e-3, 1.0),
            vec![0.75]
        );
        let transmission = render_vessel_transmittance(
            &[0.0, 0.25, 0.5, 1.0], 2, 2, 1, 550.0, 0.4, 1000.0, 1e-3, 2e-3,
        );
        let expected = [0.0, 0.25, 0.5, 1.0].iter()
            .map(|rho| (-0.8_f64 * rho).exp()).sum::<f64>() / 4.0;
        assert!((transmission[0] - expected).abs() < 8e-15);
    }

    #[test]
    fn browser_physical_pour_view_refuses_invalid_or_oversized_input() {
        let mass = [0.0, 0.25, 0.5, 1.0];
        for (nx, ny, resolution, wavelength, coefficient, depth) in [
            (3, 2, 2, 550.0, 0.4, 2e-3),
            (2, 2, 1_000_000, 550.0, 0.4, 2e-3),
            (2, 2, 2, f64::NAN, 0.4, 2e-3),
            (2, 2, 2, 550.0, -0.4, 2e-3),
            (2, 2, 2, 550.0, 0.4, 0.0),
        ] {
            assert!(
                render_vessel_transmittance(
                    &mass, nx, ny, resolution, wavelength, coefficient, 1000.0, 1e-3, depth,
                )
                .is_empty()
            );
        }
    }

    /// G1/G3: the browser consumes the native spectrum in pixel-major order.
    #[test]
    fn browser_spectral_pour_matches_native_and_preserves_pixel_band_layout() {
        let mass = [0.0, 0.25, 0.5, 1.0];
        let wavelengths = [450.0, 550.0, 650.0];
        let coefficients = [0.0, 0.4, 1000.0];
        let optics = fs_vessel::pour::PourSpectralOpticalView {
            wavelengths_nm: &wavelengths,
            mass_absorption_m2_kg: &coefficients,
            density_unit_kg_m3: 1000.0,
            lattice_spacing_m: 1e-3,
            slab_depth_m: 2e-3,
        };
        let native = optics
            .render_transmittance(
                &mass,
                [2, 2],
                SpectralAbsorptionImageSettings {
                    image: AbsorptionImageSettings {
                        resolution: [2; 2],
                        max_pixels: 4,
                        max_grid_cells: 4,
                        max_cell_visits: 4,
                    },
                    max_spectral_evaluations: 12,
                },
            )
            .unwrap();
        let browser = render_vessel_spectral_transmittance(
            &mass,
            2,
            2,
            2,
            &wavelengths,
            &coefficients,
            1000.0,
            1e-3,
            2e-3,
        );
        assert_eq!(browser, native);
        assert_eq!(browser.len(), 12);
        for (pixel, &density) in browser.chunks_exact(3).zip(&mass) {
            for (&actual, &coefficient) in pixel.iter().zip(&coefficients) {
                let expected = (-density * 2.0 * coefficient).exp();
                assert!((actual - expected).abs() < 8e-15);
            }
        }
        let coarse = render_vessel_spectral_transmittance(
            &mass,
            2,
            2,
            1,
            &wavelengths,
            &coefficients,
            1000.0,
            1e-3,
            2e-3,
        );
        for (band, &actual) in coarse.iter().enumerate() {
            let expected = mass
                .iter()
                .map(|density| (-density * 2.0 * coefficients[band]).exp())
                .sum::<f64>()
                / 4.0;
            assert!((actual - expected).abs() < 8e-15);
        }
        assert_eq!(coarse[0], 1.0);
        assert_eq!(coarse[2], 0.25);
        let single = render_vessel_spectral_transmittance(
            &mass,
            2,
            2,
            2,
            &[550.0],
            &[0.4],
            1000.0,
            1e-3,
            2e-3,
        );
        let monochromatic =
            render_vessel_transmittance(&mass, 2, 2, 2, 550.0, 0.4, 1000.0, 1e-3, 2e-3);
        assert_eq!(single.len(), 4);
        for (&spectral, &monochromatic) in single.iter().zip(&monochromatic) {
            assert!((spectral - monochromatic).abs() < 8e-15);
        }
    }

    /// G0: the browser limit applies to all wavelengths together.
    #[test]
    fn browser_spectral_pour_refuses_invalid_spectra_and_total_work() {
        for (wavelengths, coefficients) in [
            (&[][..], &[][..]),
            (&[450.0][..], &[0.4, 0.1][..]),
            (&[650.0, 450.0][..], &[0.4, 0.1][..]),
            (&[450.0, 650.0][..], &[0.4, -0.1][..]),
        ] {
            assert!(
                render_vessel_spectral_transmittance(
                    &[1.0],
                    1,
                    1,
                    1,
                    wavelengths,
                    coefficients,
                    1000.0,
                    1e-3,
                    2e-3,
                )
                .is_empty()
            );
        }
        // Pixel and cell-visit counts alone admit this image, but two bands
        // require 2,097,152 spectral evaluations. Refuse before allocation.
        assert!(
            render_vessel_spectral_transmittance(
                &[1.0],
                1,
                1,
                1024,
                &[450.0, 650.0],
                &[0.4, 0.1],
                1000.0,
                1e-3,
                2e-3,
            )
            .is_empty()
        );
    }

    /// G1: mass mixtures and spatial coverage retain their different physics.
    #[test]
    fn browser_prepared_mass_mixtures_feed_the_original_indexed_renderer() {
        let wavelengths = [450.0, 650.0];
        let a = core::f64::consts::LN_2;
        let components = [0.0, 0.0, 2.0 * a, 4.0 * a];
        let palette = prepare_vessel_spectral_absorption_mixtures(
            &wavelengths,
            &components,
            &[1.0, 0.0, 0.0, 1.0, 0.5, 0.5],
        );
        assert_eq!(palette, vec![0.0, 0.0, 2.0 * a, 4.0 * a, a, 2.0 * a]);
        assert_eq!(
            prepare_vessel_spectral_absorption_mixtures(&wavelengths, &components, &[0.0, 2.0]),
            vec![4.0 * a, 8.0 * a],
        );
        let mixed = render_vessel_indexed_spectral_transmittance(
            &[1.0], &[2], 1, 1, 1, &wavelengths, &palette, 1000.0, 1e-3, 1e-3,
        );
        let covered = render_vessel_indexed_spectral_transmittance(
            &[1.0, 1.0], &[0, 1], 2, 1, 1, &wavelengths, &palette, 1000.0, 1e-3, 1e-3,
        );
        assert_eq!(mixed.len(), 2);
        assert_eq!(covered.len(), 2);
        for ((mixed, covered), (expected_mixture, expected_coverage)) in mixed
            .into_iter()
            .zip(covered)
            .zip([(0.5, 0.625), (0.25, 0.53125)])
        {
            assert!((mixed - expected_mixture).abs() < 8e-15);
            assert!((covered - expected_coverage).abs() < 8e-15);
        }
    }

    /// G0: preparation bounds complete inputs, including unused constituents.
    #[test]
    fn browser_mass_mixture_preparation_refuses_invalid_inputs_and_excess_work() {
        let wavelengths = [450.0, 650.0];
        for (coefficients, weights) in [
            (&[0.1, 0.2, 0.3][..], &[1.0, 0.0][..]),
            (&[0.1, 0.2, 0.3, 0.4][..], &[1.0][..]),
            (&[0.1, 0.2, 0.3, 0.4][..], &[1.0, -1.0][..]),
            (&[0.1, 0.2, f64::NAN, 0.4][..], &[1.0, 0.0][..]),
        ] {
            assert!(
                prepare_vessel_spectral_absorption_mixtures(&wavelengths, coefficients, weights)
                    .is_empty()
            );
        }
        assert!(
            prepare_vessel_spectral_absorption_mixtures(
                &[550.0], &[f64::from_bits(1)], &[0.25],
            )
            .is_empty()
        );
        assert!(
            prepare_vessel_spectral_absorption_mixtures(
                &wavelengths, &vec![0.0; 65 * 2], &vec![0.0; 65],
            )
            .is_empty()
        );
        assert!(
            prepare_vessel_spectral_absorption_mixtures(
                &wavelengths, &[0.1, 0.2], &vec![0.0; 4_097],
            )
            .is_empty()
        );
        let nodes: [f64; 16] = std::array::from_fn(|lane| 400.0 + 20.0 * lane as f64);
        let components = vec![0.0; 64 * 16];
        let recipes = vec![0.0; 1_025 * 64];
        assert_eq!(
            prepare_vessel_spectral_absorption_mixtures(&nodes, &components, &recipes[..1_024 * 64])
                .len(),
            1_024 * 16,
        );
        assert!(
            prepare_vessel_spectral_absorption_mixtures(&nodes, &components, &recipes).is_empty()
        );
    }

    /// G1/G3: explicit neighboring materials retain their spectra and coverage.
    #[test]
    fn browser_indexed_spectral_pour_matches_native_materials_and_pixel_layout() {
        let mass = [1.0, 1.0];
        let materials = [0, 1];
        let wavelengths = [450.0, 650.0];
        let coefficients = [
            core::f64::consts::LN_2,
            0.0,
            0.0,
            2.0 * core::f64::consts::LN_2,
        ];
        let native = fs_vessel::pour::PourIndexedSpectralOpticalView {
            wavelengths_nm: &wavelengths,
            mass_absorption_m2_kg: &coefficients,
            cell_materials: &materials,
            density_unit_kg_m3: 1000.0,
            lattice_spacing_m: 1e-3,
            slab_depth_m: 1e-3,
        }
        .render_transmittance(
            &mass,
            [2, 1],
            IndexedSpectralAbsorptionImageSettings {
                spectral: SpectralAbsorptionImageSettings {
                    image: AbsorptionImageSettings {
                        resolution: [2, 2],
                        max_pixels: 4,
                        max_grid_cells: 2,
                        max_cell_visits: 4,
                    },
                    max_spectral_evaluations: 8,
                },
                max_materials: 2,
            },
        )
        .unwrap();
        let browser = render_vessel_indexed_spectral_transmittance(
            &mass,
            &materials,
            2,
            1,
            2,
            &wavelengths,
            &coefficients,
            1000.0,
            1e-3,
            1e-3,
        );
        assert_eq!(browser, native);
        assert_eq!(browser.len(), 8);
        for (&actual, expected) in browser
            .iter()
            .zip([0.5, 1.0, 1.0, 0.25, 0.5, 1.0, 1.0, 0.25])
        {
            assert!((actual - expected).abs() < 8e-15);
        }
        let coarse = render_vessel_indexed_spectral_transmittance(
            &mass,
            &materials,
            2,
            1,
            1,
            &wavelengths,
            &coefficients,
            1000.0,
            1e-3,
            1e-3,
        );
        assert_eq!(coarse.len(), 2);
        for (&actual, expected) in coarse.iter().zip([0.75, 0.625]) {
            assert!((actual - expected).abs() < 8e-15);
        }
    }

    /// G0: browser material maps and complete palette/work budgets are explicit.
    #[test]
    fn browser_indexed_spectral_pour_refuses_invalid_maps_and_excess_work() {
        for materials in [&[0][..], &[u32::MAX, 0][..]] {
            assert!(
                render_vessel_indexed_spectral_transmittance(
                    &[0.0, 1.0],
                    materials,
                    2,
                    1,
                    1,
                    &[450.0, 650.0],
                    &[0.4, 0.1, 0.2, 0.7],
                    1000.0,
                    1e-3,
                    1e-3,
                )
                .is_empty()
            );
        }
        let oversized_palette = vec![0.0; 4_097 * 2];
        assert!(
            render_vessel_indexed_spectral_transmittance(
                &[1.0],
                &[0],
                1,
                1,
                1,
                &[450.0, 650.0],
                &oversized_palette,
                1000.0,
                1e-3,
                1e-3,
            )
            .is_empty()
        );
        // Scalar pixel/cell limits admit this image; two wavelengths exceed
        // the total 1,048,576 cell/wavelength evaluation cap.
        assert!(
            render_vessel_indexed_spectral_transmittance(
                &[1.0],
                &[0],
                1,
                1,
                1024,
                &[450.0, 650.0],
                &[0.4, 0.1],
                1000.0,
                1e-3,
                1e-3,
            )
            .is_empty()
        );
    }

    #[test]
    fn browser_illuminated_pour_retains_original_column_correlation_and_layout() {
        let wavelengths = [450.0, 650.0];
        let coefficients = [core::f64::consts::LN_2, 2.0 * core::f64::consts::LN_2];
        let background = [2.0, 4.0, 10.0, 20.0];
        let render = |resolution, source: &[f64], per_column| {
            render_vessel_spectral_radiance(
                &[0.0, 1.0],
                2,
                1,
                resolution,
                &wavelengths,
                &coefficients,
                1000.0,
                1e-3,
                1e-3,
                source,
                per_column,
            )
        };
        let coarse = render(1, &background, true);
        assert_eq!(coarse.len(), 2);
        for (&actual, expected) in coarse.iter().zip([3.5, 4.5]) {
            assert!((actual - expected).abs() < 8e-15);
        }
        let uniform = render(1, &[6.0, 12.0], false);
        for (&actual, expected) in uniform.iter().zip([4.5, 7.5]) {
            assert!((actual - expected).abs() < 8e-15);
        }
        // The source remains two original columns while output has four pixels.
        let fine = render(2, &background, true);
        assert_eq!(fine.len(), 8);
        for (&actual, expected) in fine.iter().zip([2.0, 4.0, 5.0, 5.0, 2.0, 4.0, 5.0, 5.0]) {
            assert!((actual - expected).abs() < 8e-15);
        }
        let palette = [0.0, 0.0, coefficients[0], coefficients[1]];
        let indexed = render_vessel_indexed_spectral_radiance(
            &[1.0, 1.0],
            &[0, 1],
            2,
            1,
            2,
            &wavelengths,
            &palette,
            1000.0,
            1e-3,
            1e-3,
            &background,
            true,
        );
        assert_eq!(indexed.len(), fine.len());
        for (&actual, expected) in indexed.iter().zip(fine) {
            assert!((actual - expected).abs() < 8e-15);
        }
    }

    #[test]
    fn browser_illuminated_pour_xyz_matches_native_observation_and_flat_layout() {
        use fs_render::spectral::{cie_x, cie_y, cie_z, y_integral};

        let wavelengths = [450.0, 650.0];
        let coefficients = [core::f64::consts::LN_2, 2.0 * core::f64::consts::LN_2];
        let background = [2.0, 4.0, 10.0, 20.0];
        let weights = [20.0, 40.0];
        let optics = fs_vessel::pour::PourSpectralOpticalView {
            wavelengths_nm: &wavelengths,
            mass_absorption_m2_kg: &coefficients,
            density_unit_kg_m3: 1000.0,
            lattice_spacing_m: 1e-3,
            slab_depth_m: 1e-3,
        };
        let native = optics
            .render_xyz(
                &[0.0, 1.0],
                [2, 1],
                SpectralAbsorptionImageSettings {
                    image: AbsorptionImageSettings {
                        resolution: [2, 2],
                        max_pixels: 4,
                        max_grid_cells: 2,
                        max_cell_visits: 4,
                    },
                    max_spectral_evaluations: 8,
                },
                SpectralIllumination::Columns(&background),
                &weights,
            )
            .unwrap();
        let browser = render_vessel_spectral_xyz(
            &[0.0, 1.0],
            2,
            1,
            2,
            &wavelengths,
            &coefficients,
            1000.0,
            1e-3,
            1e-3,
            &background,
            true,
            &weights,
        );
        assert_eq!(browser, native.into_flattened());
        assert_eq!(browser.len(), 12);
        for (pixel, radiance) in
            browser
                .chunks_exact(3)
                .zip([[2.0, 4.0], [5.0, 5.0], [2.0, 4.0], [5.0, 5.0]])
        {
            for (channel, cmf) in [cie_x, cie_y, cie_z].into_iter().enumerate() {
                let expected = (cmf(450.0) * weights[0] * radiance[0]
                    + cmf(650.0) * weights[1] * radiance[1])
                    / y_integral();
                assert!((pixel[channel] - expected).abs() < 8e-15);
            }
        }
        let palette = [0.0, 0.0, coefficients[0], coefficients[1]];
        let indexed = render_vessel_indexed_spectral_xyz(
            &[1.0, 1.0],
            &[0, 1],
            2,
            1,
            2,
            &wavelengths,
            &palette,
            1000.0,
            1e-3,
            1e-3,
            &background,
            true,
            &weights,
        );
        assert_eq!(indexed.len(), browser.len());
        for (&actual, expected) in indexed.iter().zip(browser) {
            assert!((actual - expected).abs() < 8e-15);
        }
    }

    #[test]
    fn browser_illuminated_pour_refuses_invalid_sources_observers_and_work() {
        let mass = [0.0, 1.0];
        let wavelengths = [450.0, 650.0];
        let coefficients = [0.4, 0.1];
        for (source, per_column) in [
            (&[1.0][..], false),
            (&[1.0, 2.0][..], true),
            (&[1.0, 2.0, -1.0, 4.0][..], true),
            (&[1.0, 2.0, f64::NAN, 4.0][..], true),
        ] {
            assert!(
                render_vessel_spectral_radiance(
                    &mass,
                    2,
                    1,
                    1,
                    &wavelengths,
                    &coefficients,
                    1000.0,
                    1e-3,
                    1e-3,
                    source,
                    per_column,
                )
                .is_empty()
            );
            assert!(
                render_vessel_spectral_xyz(
                    &mass,
                    2,
                    1,
                    1,
                    &wavelengths,
                    &coefficients,
                    1000.0,
                    1e-3,
                    1e-3,
                    source,
                    per_column,
                    &[20.0, 40.0],
                )
                .is_empty()
            );
        }
        for weights in [&[1.0][..], &[1.0, -1.0][..], &[1.0, f64::INFINITY][..]] {
            assert!(
                render_vessel_spectral_xyz(
                    &mass,
                    2,
                    1,
                    1,
                    &wavelengths,
                    &coefficients,
                    1000.0,
                    1e-3,
                    1e-3,
                    &[1.0, 2.0],
                    false,
                    weights,
                )
                .is_empty()
            );
        }
        for materials in [&[0][..], &[u32::MAX, 0][..]] {
            assert!(
                render_vessel_indexed_spectral_radiance(
                    &mass,
                    materials,
                    2,
                    1,
                    1,
                    &wavelengths,
                    &coefficients,
                    1000.0,
                    1e-3,
                    1e-3,
                    &[1.0, 2.0],
                    false,
                )
                .is_empty()
            );
            assert!(
                render_vessel_indexed_spectral_xyz(
                    &mass,
                    materials,
                    2,
                    1,
                    1,
                    &wavelengths,
                    &coefficients,
                    1000.0,
                    1e-3,
                    1e-3,
                    &[1.0, 2.0],
                    false,
                    &[20.0, 40.0],
                )
                .is_empty()
            );
        }
        // Each image would need 2,097,152 spectral evaluations: refuse before allocation.
        assert!(
            render_vessel_spectral_radiance(
                &[1.0],
                1,
                1,
                1024,
                &wavelengths,
                &coefficients,
                1000.0,
                1e-3,
                1e-3,
                &[1.0, 2.0],
                false,
            )
            .is_empty()
        );
        assert!(
            render_vessel_indexed_spectral_xyz(
                &[1.0],
                &[0],
                1,
                1,
                1024,
                &wavelengths,
                &coefficients,
                1000.0,
                1e-3,
                1e-3,
                &[1.0, 2.0],
                false,
                &[20.0, 40.0],
            )
            .is_empty()
        );
        assert!(
            render_vessel_spectral_xyz(
                &[1.0],
                1,
                1,
                1,
                &[900.0, 1000.0],
                &coefficients,
                1000.0,
                1e-3,
                1e-3,
                &[1.0, 2.0],
                false,
                &[20.0, 40.0],
            )
            .is_empty()
        );
    }

    #[test]
    fn standalone_qty_dependency_round_trips_nonzero_mole_axis() {
        let amount = QtyAny::new(2.5, Dims([0, 0, 0, 0, 0, 1]));
        let canonical = fs_qty::json::to_json(amount).expect("finite amount encodes");
        assert_eq!(
            canonical,
            r#"{"schema_version":2,"value":2.5,"dims":[0,0,0,0,0,1]}"#
        );
        let decoded = fs_qty::json::decode_json(&canonical).expect("canonical v2 decodes");
        assert_eq!(decoded.qty(), amount);
        assert!(
            decoded.migration().is_none(),
            "canonical six-base input must not mint a legacy crosswalk receipt"
        );
    }

    #[test]
    fn ornithoid_headline() {
        let v = run_ornithoid(1);
        let adj_rel_err = v[9];
        let nrows = v[6] as usize;
        let polish_before = v[24];
        let polish_after = v[25];
        eprintln!(
            "ornithoid: adj_rel_err={adj_rel_err:e} nrows={nrows} polish {polish_before:.4}→{polish_after:.4} \
             dcl_adjoint={} race_winner={} evals={} fixed_n={} elim={} \
             lbm_lift={} lbm_drag={} steadiness={:e} roa={} coverage={} hv={} knee={}",
            v[10], v[11], v[12], v[13], v[14], v[15], v[16], v[18], v[19], v[20], v[22], v[23]
        );
        assert!(adj_rel_err < 1e-5, "adj_rel_err {adj_rel_err:e}");
        assert!(nrows >= 12, "atlas rows {nrows}");
        assert!(
            polish_after > polish_before,
            "polish {polish_before} !> {polish_after}"
        );
    }

    #[test]
    fn vessel_headline() {
        let v = run_vessel(1000);
        let robust_offband = v[17];
        let nominal_offband = v[18];
        let profile_samples = v[1] as usize;
        let growth_stations = v[2] as usize;
        let nx = v[4] as usize;
        let ny = v[5] as usize;
        let frames = v[6] as usize;
        let render_res = v[7] as usize;
        let cvar_rows = v[8] as usize;
        let cvar_start = 20 + 2 * profile_samples + 2 * growth_stations + frames * nx * ny;
        let tail = v.len();
        let robust_lip = v[tail - 3];
        let nominal_lip = v[tail - 2];
        // render block precedes the 3-value tail.
        let render_start = tail - 3 - render_res * render_res;
        let mut tmin = f64::INFINITY;
        let mut tmax = f64::NEG_INFINITY;
        for &t in &v[render_start..render_start + render_res * render_res] {
            if t.is_finite() {
                tmin = tmin.min(t);
                tmax = tmax.max(t);
            }
        }
        eprintln!(
            "vessel: growth_minmax={} offnom={} robust_offband={robust_offband} nominal_offband={nominal_offband} \
             robust_lip={robust_lip} nominal_lip={nominal_lip} drift={:e} poured_neutral={} \
             contact_band={} fragments={} transmittance∈[{tmin},{tmax}]",
            v[10], v[11], v[12], v[13], v[14], v[16]
        );
        assert!(
            robust_offband < nominal_offband,
            "robust_offband {robust_offband} !< nominal_offband {nominal_offband}"
        );
        assert!(
            v[cvar_start..cvar_start + 3 * cvar_rows]
                .as_chunks::<3>()
                .0
                .iter()
                .all(|row| row[2].is_finite()),
            "canonical CVaR returned a non-finite value for a finite vessel fluid band"
        );
        assert!(
            tmin >= 0.0 && tmax <= 1.0,
            "transmittance out of [0,1]: [{tmin},{tmax}]"
        );
        assert!(
            tmax - tmin > 0.05,
            "transmittance range too small: [{tmin},{tmax}]"
        );
    }

    #[test]
    fn vessel_cvar_validation_refuses_without_trapping() {
        assert!(canonical_cvar_or_nan(&[], 0.8).is_nan());
        assert!(canonical_cvar_or_nan(&[1.0], f64::NAN).is_nan());
        assert!(canonical_cvar_or_nan(&[1.0, f64::INFINITY], 0.8).is_nan());
    }

    #[test]
    fn browser_pour_frames_come_from_the_physical_native_trajectory() {
        let (nx, ny, captures) = (32, 20, 5);
        let liquid = Rheology::Newtonian { nu: 0.016_667 };
        let recorded = pour_inline(nx, ny, captures, ContactModel::Neutral, liquid)
            .expect("bounded physical pour must remain representable");
        let unrecorded = pour_inline(nx, ny, 0, ContactModel::Neutral, liquid)
            .expect("capture policy cannot change the physical trajectory");
        assert_eq!(recorded.frames.len(), captures * nx * ny);
        assert_eq!(
            &recorded.frames[(captures - 1) * nx * ny..],
            &recorded.mass_field
        );
        assert_eq!(recorded.mass_field, unrecorded.mass_field);
        assert!(recorded.mass_field.iter().all(|mass| mass.is_finite()));
        assert!(recorded.drift < 1e-10, "mass drift {}", recorded.drift);
        assert!(recorded.poured > 0.0);
    }

    #[test]
    fn browser_pour_refuses_clipped_fluids_without_publishing_frames() {
        for nu in [1e-5, 2.0] {
            assert!(
                pour_inline(12, 8, 2, ContactModel::Neutral, Rheology::Newtonian { nu }).is_none()
            );
        }
    }

    #[test]
    fn flagship_determinism() {
        // No clocks, no entropy RNG: identical inputs → bit-identical output.
        let a = run_ornithoid(7);
        let b = run_ornithoid(7);
        assert_eq!(a.len(), b.len());
        assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
        let a = run_vessel(1500);
        let b = run_vessel(1500);
        assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
        let a = run_frame(1234);
        let b = run_frame(1234);
        assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
    }

    #[test]
    fn flagship_timings() {
        use std::time::Instant;
        for &seed in &[1u32, 42, 90210] {
            let t = Instant::now();
            let o = run_ornithoid(seed);
            let d_o = t.elapsed();
            let t = Instant::now();
            let ve = run_vessel(1000);
            let d_v = t.elapsed();
            let t = Instant::now();
            let f = run_frame(seed);
            let d_f = t.elapsed();
            eprintln!(
                "TIMING seed={seed}: ornithoid {:.1}ms (len {}), vessel {:.1}ms (len {}), frame {:.1}ms (len {})",
                d_o.as_secs_f64() * 1e3,
                o.len(),
                d_v.as_secs_f64() * 1e3,
                ve.len(),
                d_f.as_secs_f64() * 1e3,
                f.len()
            );
        }
    }

    #[test]
    fn frame_history_preserves_initial_rows_and_native_checked_response() {
        let params = fs_frame::StoryParams::default();
        let acceleration = [0.8, 0.2, -0.4, 0.1, 0.3, -0.2];
        let dt = 0.01;
        let native = fs_frame::StoryFrame::new(params)
            .run_from_time_zero(
                fs_frame::GroundMotion::new(&acceleration, dt),
                fs_frame::HistoryLimits::new(acceleration.len(), 30, 1e-12, 1.0),
            )
            .unwrap();
        let (displacement, shear, stiffness) = story_history(&params, &acceleration, dt).unwrap();
        assert_eq!(displacement.len(), acceleration.len());
        assert_eq!(shear.len(), acceleration.len());
        assert_eq!(displacement[0], 0.0);
        assert_eq!(shear[0], 0.0);
        assert_eq!(&displacement[1..], native.displacement_m.as_slice());
        assert_eq!(&shear[1..], native.restoring_shear_n.as_slice());
        assert_eq!(
            stiffness,
            fs_frame::StoryFrame::new(params).initial_stiffness()
        );
        assert!(story_history(&params, &[0.0, f64::NAN], dt).is_err());
    }

    #[test]
    fn frame_headline() {
        let v = run_frame(90210);
        assert_eq!(v[0], 2.0, "MAGIC");
        assert_eq!(v[1], 2.0, "wire version");
        let off_layout = v[3] as usize;
        let off_sizing = v[4] as usize;
        let off_fragility = v[6] as usize;
        let gap = v[off_layout];
        let eq_residual = v[off_layout + 1];
        let solver_converged = v[off_layout + 4];
        let optimality = &v[off_sizing - 4..off_sizing];
        let all_pass = v[off_sizing];
        let sizing_eq = v[off_sizing + 1];
        let p_hat = v[off_fragility];
        let radius = v[off_fragility + 1];
        let members_used = v[off_fragility + 2];
        let stopped_early = v[off_fragility + 3];
        let exceedances = v[off_fragility + 6];
        let mlmc_estimate = v[off_fragility + 9];
        eprintln!(
            "frame: gap={gap:e} eq_residual={eq_residual:e} solver_converged={solver_converged} \
             all_pass={all_pass} sizing_eq={sizing_eq:e} p_hat={p_hat} radius={radius} \
             members_used={members_used} stopped_early={stopped_early} exceedances={exceedances} \
             mlmc_estimate={mlmc_estimate}"
        );
        assert!(gap < 1e-3, "layout gap {gap:e}");
        assert_eq!(solver_converged, 1.0, "layout solver convergence");
        assert_eq!(optimality[0], 2.0, "layout optimality rank");
        assert_eq!(optimality[1], 1.0, "layout verified flag");
        assert!(optimality[2].is_finite() && optimality[3].is_finite());
        assert!(optimality[2] > 0.0 && optimality[2] <= optimality[3]);
        assert_eq!(all_pass, 1.0, "sizing all_pass");
        assert!(p_hat > 0.0 && p_hat < 1.0, "fragility p_hat {p_hat}");
        assert!(members_used < 48.0, "members_used {members_used}");
        assert_eq!(v[off_fragility + 10], 1.0, "original-spectrum level count");
        let off_history = v[5] as usize;
        let history_samples = v[off_history + 2] as usize;
        assert_eq!(
            v[off_history + 5 + history_samples],
            0.0,
            "initial displacement"
        );
        assert_eq!(
            v[off_history + 5 + 2 * history_samples],
            0.0,
            "initial restoring shear"
        );
        let curve_samples = v[off_fragility + 11] as usize;
        let last_curve = off_fragility + 12 + 3 * (curve_samples - 1);
        assert_eq!(
            v[last_curve + 1].to_bits(),
            p_hat.to_bits(),
            "same CS replay"
        );
        assert_eq!(
            v[last_curve + 2].to_bits(),
            radius.to_bits(),
            "same CS radius"
        );

        let off_cvar = v[7] as usize;
        let scale_star = v[off_cvar];
        let scale_snapped = v[off_cvar + 1];
        let cvar_star = v[off_cvar + 2];
        let cvar_snapped = v[off_cvar + 3];
        let limit = v[off_cvar + 4];
        let mass = v[off_cvar + 5];
        let beta = v[off_cvar + 7];
        let cvar_curve_len = v[off_cvar + 8] as usize;
        let cvar_cat_len = v[off_cvar + 9] as usize;
        assert!(
            scale_star.is_finite() && scale_star > 0.0,
            "scale_star {scale_star}"
        );
        assert!(
            scale_snapped.is_finite() && scale_snapped >= scale_star,
            "scale_snapped {scale_snapped}"
        );
        assert!(
            cvar_star.is_finite() && cvar_star > 0.0,
            "cvar_star {cvar_star}"
        );
        assert!(
            cvar_snapped.is_finite() && cvar_snapped <= limit,
            "cvar_snapped {cvar_snapped} <= limit {limit}"
        );
        assert!(mass.is_finite() && mass > 0.0, "mass {mass}");
        assert_eq!(beta, 0.9, "beta");
        assert_eq!(cvar_curve_len, 7, "cvar_curve_len");
        assert_eq!(cvar_cat_len, 7, "cvar_cat_len");
    }

    #[test]
    fn frame_nan_body_on_malformed_ensemble() {
        let nan_body = frame_nan_body(42);
        assert_eq!(nan_body[0], 2.0, "MAGIC");
        assert_eq!(nan_body[1], 2.0, "wire version");
        assert_eq!(nan_body[2], 42.0, "seed");
        // Remaining payload should be NaNs
        assert!(nan_body[12..].iter().all(|x| x.is_nan()));
    }
}
