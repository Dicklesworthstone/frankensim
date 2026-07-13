//! fs-wasm · FLAGSHIP tier (Tier V) — three live, end-to-end evidence-bearing
//! design pipelines in the browser, each a sized-down but genuinely real run
//! of a FrankenSim flagship campaign.
//!
//! - [`run_ornithoid`] — the ornithoid (flapping micro-flyer) pipeline:
//!   PARAMETERIZE → SCREEN (e-race) → VALIDATE (LBM) → CERTIFY (Lyapunov +
//!   conformal surrogate) → EXPLORE (NSGA-II Pareto + knee polish).
//! - [`run_vessel`] — the never-dribbling carafe: PARAMETERIZE → STABILITY
//!   (Orr–Sommerfeld min-max) → VALIDATE (free-surface pour) → ROBUSTIFY (CVaR
//!   over a fluid band) → RENDER (Woodcock volume transmittance).
//! - [`run_frame`] — the e-stopped seismic frame: LAYOUT (fnx-free truss LP) →
//!   SIZE (yield + Euler-buckling, catalog snap) → TIME HISTORY (nonlinear
//!   fiber hinge) → FRAGILITY (anytime CS + MLMC) → CVaR mass minimization.
//!
//! SAFETY CONTRACT (identical to the rest of the crate): `unsafe_code` is
//! forbidden, every input is clamped, every fallible kernel result is folded to
//! `NaN` / an empty vector, and every documented panic precondition of the
//! composed crates is respected. Nothing here can trap — a wasm trap would kill
//! the whole page. In particular the LAYOUT stage of the frame flagship uses a
//! fnx-free vendor of the truss numerical core (the compatibility-evidence path
//! reached while `fs_truss::GroundStructure::try_grid` builds its `fnx` graph
//! reads a wall clock, which compiles but TRAPS at runtime on
//! `wasm32-unknown-unknown`); the trap-free racing path
//! (`fs_race::race_field` over an empty `KillRegistry`) never reads a clock.
//!
//! Determinism: no clocks, no entropy RNG. All stochastic paths are
//! counter-based Philox keyed by seed / logical identity.

use fs_bem::panel2d::{Airfoil2d, dcl_dalpha_adjoint, solve};
use fs_bem::wake2d::WakeSim;
use fs_dfo::moo::{Individual, NsgaParams, hypervolume, knee_point, nsga2};
use fs_eproc::GaussianMixtureCs;
use fs_exec::KillRegistry;
use fs_lbm::rheology::{Rheology, update_tau};
use fs_lbm::{Cell, ContactModel, FreeSurface, Grid, equilibrium};
use fs_ornith::certify::LdSurrogate;
use fs_ornith::param::OrnithCandidate;
use fs_ornith::screen::{PANELS, lift_to_drag};
use fs_qty::{Dims, QtyAny};
use fs_race::{LossSpan, RaceSettings, race_field};
use fs_rand::StreamKey;
use fs_render::volumes::{MajorantGrid, VolumeGrid, render_transmittance};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};
use fs_solid::fiber::{Section, rc_section};
use fs_sos::lyapunov_certifies_stability;
use fs_sparse::{Coo, Csr};
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

/// The D2Q9 lattice velocities, in `fs_lbm`'s canonical order (its own `E`
/// table is `pub(crate)`, so it is mirrored here for the momentum-exchange
/// boundary force — the value read from `Grid::moments` needs no table).
const E9: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

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
    grid.g = [gx, 0.0];

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
    let mut prev = velfield(&grid);
    for s in 0..steps {
        grid.step(&mut scratch);
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

    // momentum-exchange body force (Ladd): sum over fluid→wall links.
    let (mut fx, mut fy) = (0.0f64, 0.0f64);
    for y in 0..ny {
        for x in 0..nx {
            let i = grid.idx(x, y);
            if !matches!(grid.flags[i], Cell::Fluid) {
                continue;
            }
            for (q, &(ex, ey)) in E9.iter().enumerate() {
                let nxx = x as i64 + ex as i64;
                let nyy = y as i64 + ey as i64;
                if nxx < 0 || nyy < 0 || nxx >= nx as i64 || nyy >= ny as i64 {
                    continue;
                }
                let j = grid.idx(nxx as usize, nyy as usize);
                if matches!(grid.flags[j], Cell::Wall) {
                    let f = 2.0 * grid.f[i][q];
                    fx += f64::from(ex) * f;
                    fy += f64::from(ey) * f;
                }
            }
        }
    }
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
///   - `[10]` `spectral_growth_minmax` — certified min-max at the input lip
///     (real `growth_objective`, `n = 32`)
///   - `[11]` `spectral_growth_offnom` — worst off-nominal growth at input lip
///   - `[12]` `mass_ledger_residual` — worst pour mass-ledger drift
///   - `[13]` `poured_mass_neutral` — mass past the lip (Neutral contact)
///   - `[14]` `contact_poured_band` — |poured(Neutral) − poured(Wetting)|
///   - `[15]` `contact_dribble_band` — |dribble(Neutral) − dribble(Wetting)|
///   - `[16]` `fragments` — Plateau–Rayleigh fragment count (Neutral)
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
    let neutral = pour_inline(nx, ny, frames, ContactModel::Neutral, true);
    let wetting = pour_inline(nx, ny, 0, ContactModel::Wetting, false);
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
        let cvar = if losses.iter().all(|v| v.is_finite()) {
            fs_vessel::robust::empirical_cvar(&losses, beta)
        } else {
            f64::NAN
        };
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

    // ---- STAGE 5: Woodcock volume render of the final pour mass field ----
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

/// Inline free-surface pour on an `nx×ny` lattice (transcribes the setup of
/// `fs_vessel::pour::run_pour`), capturing `capture` evenly-spaced mass frames.
fn pour_inline(
    nx: usize,
    ny: usize,
    capture: usize,
    contact: ContactModel,
    _keep_frames: bool,
) -> PourOut {
    let g0 = 6.0e-4;
    let tilt_final = 0.7;
    let steps = 300u32;
    let lip_x = nx * 5 / 8;
    let lip_h = (ny * 2 / 5).min(ny - 2);

    let mut grid = Grid::uniform(nx, ny, 0.55);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = [0.0, -g0];
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
    let law = Rheology::Newtonian { nu: 0.016_667 };
    for step in 0..steps {
        let theta = tilt_final * f64::from(step) / f64::from(steps);
        sim.grid.g = [g0 * theta.sin(), -g0 * theta.cos()];
        let _ = update_tau(&mut sim.grid, law);
        sim.step();
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
    PourOut {
        frames,
        mass_field,
        drift,
        fragments,
        poured,
        dribble,
    }
}

/// Woodcock volume transmittance of a thin z-slab whose extinction field is the
/// final pour mass map (`fs_render::volumes`).
fn render_thin_slab(mass: &[f64], nx: usize, ny: usize, res: usize) -> Vec<f64> {
    if mass.len() != nx * ny {
        return vec![f64::NAN; res * res];
    }
    let grid = VolumeGrid::new([nx, ny, 1], mass, [0.0, 0.0, 0.0], [1.0, 1.0, 4.0]);
    let majorant = MajorantGrid::build(&grid, 8);
    render_transmittance(&grid, &majorant, res, 16, 0x7E55_E1D0)
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
/// construction-evidence path reads a wall clock — a wasm runtime trap).
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

    fn solve(&self, max_iters: usize, gap_tol: f64, check_every: usize) -> (Vec<f64>, FrameReport) {
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
        (x, report)
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

/// The `[m, kg, s, K, A]` dimension of time (s) and rate (1/s).
const TIME: Dims = Dims([0, 0, 1, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0]);

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

/// A faithful capture of the nonlinear story time history (transcribes
/// `fs_frame::history::StoryFrame::run` with a locally-owned fiber section so
/// the true path-dependent shear `V(t)` — the hysteresis loop — is recorded).
/// Returns `(x_history, v_history, k0)`.
fn story_history(params: &fs_frame::StoryParams, ag: &[f64], dt: f64) -> (Vec<f64>, Vec<f64>, f64) {
    let mut hinge: Section = rc_section(0.5, 0.35, 12, 0.002);
    for f in &mut hinge.fibers {
        f.area *= params.scale;
    }
    let restoring = |hinge: &Section, x: f64| -> (f64, f64) {
        let kappa = x / (params.h * params.lp);
        let st = hinge.respond(0.0, kappa);
        let v = 2.0 * st.m / params.h;
        let dv_dx = 2.0 * st.tangent[1][1] / (params.h * params.h * params.lp);
        (v, dv_dx)
    };
    let k0 = restoring(&hinge, 1e-9).1;
    let mm = params.mass;
    let c = 2.0 * params.zeta * (k0 * mm).sqrt();
    let (beta, gamma) = (0.25f64, 0.5f64);
    let (mut cx, mut cv, mut ca) = (0.0f64, 0.0f64, 0.0f64);
    let mut xs = Vec::with_capacity(ag.len());
    let mut vs = Vec::with_capacity(ag.len());
    for &agi in ag {
        let p_ext = -mm * agi;
        let (x0, v0, a0) = (cx, cv, ca);
        let mut x = x0;
        for _ in 0..30 {
            let a_new = (x - x0 - dt * v0) / (beta * dt * dt) - (0.5 - beta) / beta * a0;
            let v_new = v0 + dt * ((1.0 - gamma) * a0 + gamma * a_new);
            let (fs, kt) = restoring(&hinge, x);
            let r = mm * a_new + c * v_new + fs - p_ext;
            let kdyn = mm / (beta * dt * dt) + c * gamma / (beta * dt) + kt;
            let dx = -r / kdyn;
            x += dx;
            if dx.abs() < 1e-12 {
                break;
            }
        }
        let a_new = (x - x0 - dt * v0) / (beta * dt * dt) - (0.5 - beta) / beta * a0;
        let v_new = v0 + dt * ((1.0 - gamma) * a0 + gamma * a_new);
        let (shear, _) = restoring(&hinge, x);
        let kappa = x / (params.h * params.lp);
        hinge.commit(0.0, kappa);
        cx = x;
        cv = v_new;
        ca = a_new;
        xs.push(x);
        vs.push(shear);
    }
    (xs, vs, k0)
}

/// **FRAME** — the e-stopped seismic frame flagship, run end to end at reduced
/// size. Five stages, one flat `Vec<f64>` with a self-describing offset header.
///
/// `seed` — the ensemble Philox seed.
///
/// OUTPUT LAYOUT (flat `Vec<f64>`):
/// - HEADER, 12 values (offsets are f64 indices into this array):
///   `[0]` MAGIC (= 2), `[1]` version, `[2]` seed, `[3]` off_layout,
///   `[4]` off_sizing, `[5]` off_history, `[6]` off_fragility, `[7]` off_cvar,
///   `[8]` total_len, `[9..12]` reserved (0).
/// - LAYOUT block (@ off_layout): `gap, eq_residual, volume_phys, iters,
///   solver_converged, Nn, M, load_node_idx`; then `2·Nn` node coords
///   `[x,y]…`; then `M·4` `[na, nb, force_q, is_survivor]`.
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
    let lp = FrameLp::assemble(&gs, &supported, &loads, 1.0);
    let (x, report) = lp.solve(60_000, 1e-4, 500);
    let force = |k: usize| x[k] - x[m + k];
    let max_force = (0..m).map(|k| force(k).abs()).fold(0.0, f64::max);
    let active_tol = 1e-3 * max_force.max(1e-12);
    let solver_converged = report.gap < 1e-3 && report.eq_residual < 1e-3;

    // ---- STAGE 2: SIZING (yield + Euler-buckling, catalog snap) ----------
    let sigma_y = 250.0e6;
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
    let (xs, vs, k0) = story_history(&base, &ag, dt);
    let ns_hist = ag.len();
    let peak = fs_frame::peak_drift(&xs, base.h);
    let vy = vs.iter().fold(0.0f64, |m, &v| m.max(v.abs()));

    // choose a drift_limit at the median member peak drift so p_hat ∈ (0,1)
    // (a thinned pre-pass — every third member — keeps this cheap).
    let mut peaks: Vec<f64> = (0..members)
        .step_by(3)
        .filter_map(|mm| ensemble.realize(mm).ok())
        .map(|r| {
            let mut fr = fs_frame::StoryFrame::new(base);
            fs_frame::peak_drift(&fr.run(&r.values, dt), base.h)
        })
        .collect();
    peaks.sort_by(f64::total_cmp);
    let drift_limit = if peaks.is_empty() {
        base.h * 0.01
    } else {
        peaks[peaks.len() / 2]
    };

    // ---- STAGE 4: FRAGILITY (real e-stopped anytime CS + MLMC) ----------
    let (alpha, margin) = (0.1f64, 0.30f64);
    let frag = fs_frame::e_stopped_fragility(&ensemble, base, drift_limit, alpha, margin);
    // running CS curve: replay members_used with the same GaussianMixtureCs.
    let mut cs = GaussianMixtureCs::new(0.5, 8.0, alpha);
    let mut cs_curve: Vec<(usize, f64, f64)> = Vec::new();
    for member in 0..frag.members_used {
        if let Ok(r) = ensemble.realize(member) {
            let mut fr = fs_frame::StoryFrame::new(base);
            let pd = fs_frame::peak_drift(&fr.run(&r.values, dt), base.h);
            cs.observe(if pd > drift_limit { 1.0 } else { 0.0 });
            if let Some((center, radius)) = cs.interval() {
                cs_curve.push((member as usize, center, radius));
            }
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
    let cvar_curve: Vec<(f64, f64)> = cvar_catalog
        .iter()
        .map(|&s| (s, fs_frame::ensemble_cvar(&ens_cvar, base, s, beta)))
        .collect();
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
        Some(fs_frame::cvar_mass_min(
            &ens_cvar,
            base,
            beta,
            limit,
            &cvar_catalog,
        ))
    } else {
        None
    };

    /* ---- assemble the flat output (build blocks, then the offset header) */
    let mut layout: Vec<f64> = Vec::new();
    layout.push(fon(report.gap));
    layout.push(fon(report.eq_residual));
    layout.push(fon(report.volume)); // volume_phys (sigma_y = 1 in the LP)
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
    out.push(1.0); // version
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
    let mut out = vec![2.0, 1.0, f64::from(seed)];
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
        let render_res = v[7] as usize;
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
            tmin >= 0.0 && tmax <= 1.0,
            "transmittance out of [0,1]: [{tmin},{tmax}]"
        );
        assert!(
            tmax - tmin > 0.05,
            "transmittance range too small: [{tmin},{tmax}]"
        );
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
    fn frame_headline() {
        let v = run_frame(90210);
        assert_eq!(v[0], 2.0, "MAGIC");
        let off_layout = v[3] as usize;
        let off_sizing = v[4] as usize;
        let off_fragility = v[6] as usize;
        let gap = v[off_layout];
        let eq_residual = v[off_layout + 1];
        let solver_converged = v[off_layout + 4];
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
        assert_eq!(all_pass, 1.0, "sizing all_pass");
        assert!(p_hat > 0.0 && p_hat < 1.0, "fragility p_hat {p_hat}");
        assert!(members_used < 48.0, "members_used {members_used}");
    }
}
