//! The compliance descent loop — THE marquee coupling: physics by
//! fs-solid's CutFEM elasticity DIRECTLY on the evolving [`GridSdf`]
//! (zero meshing, ever), shape velocity `v_n = w − ℓ` from the energy
//! density (grow where strain energy exceeds the volume multiplier,
//! shrink where it doesn't), Sobolev-smoothed through fs-adjoint's
//! Riesz step, extended off the interface, advected by WENO5 normal
//! flow, redistanced with drift audits, and — on schedule — hole
//! nucleation by topological derivative. Volume rides an augmented
//! Lagrangian multiplier; every iteration ledgers compliance, volume,
//! multiplier, audits, events, and an FNV snapshot hash of φ.

use crate::fim::{RedistanceAudit, redistance};
use crate::gridsdf::GridSdf;
use crate::topder::{NucleationEvent, nucleate, topological_derivative};
use crate::veloext::extend_velocity;
use crate::weno::{Velocity, advect, build_band};
use fs_cutfem::quad::cut_cell_rules;
use fs_cutfem::{CutSdf, Quadtree};
use fs_solid::linear::lame;
use fs_solid::{CutElasticity, PlaneKind, SolidError};
use std::fmt::Write as _;

/// Optimizer controls.
#[derive(Debug, Clone, Copy)]
pub struct OptimizeSettings {
    /// Grid level (cells per side = 2^level; GridSdf n must match).
    pub level: u32,
    /// Target volume fraction of the design box.
    pub volfrac: f64,
    /// Descent iterations.
    pub iterations: usize,
    /// Narrow band half-width in cells.
    pub band_cells: f64,
    /// Interface travel per iteration, in cells.
    pub move_cells: f64,
    /// Initial volume multiplier ℓ.
    pub ell0: f64,
    /// Augmented-Lagrangian multiplier gain.
    pub mu_al: f64,
    /// Sobolev smoothing α (≈ h²·scale).
    pub sobolev_alpha: f64,
    /// Nucleation period (0 disables).
    pub nucleation_period: usize,
    /// Nucleation hole radius (in cells).
    pub hole_radius_cells: f64,
    /// Young's modulus.
    pub youngs: f64,
    /// Poisson ratio.
    pub poisson: f64,
}

impl Default for OptimizeSettings {
    fn default() -> Self {
        OptimizeSettings {
            level: 5,
            volfrac: 0.5,
            iterations: 25,
            band_cells: 6.0,
            move_cells: 0.5,
            ell0: 0.0,
            mu_al: 4.0,
            sobolev_alpha: 2.0,
            nucleation_period: 8,
            hole_radius_cells: 2.5,
            youngs: 1.0,
            poisson: 0.3,
        }
    }
}

/// The ledgered trajectory.
#[derive(Debug, Clone, Default)]
pub struct OptimizeReport {
    /// Compliance per iteration.
    pub compliance: Vec<f64>,
    /// Material volume per iteration.
    pub volume: Vec<f64>,
    /// Multiplier per iteration.
    pub ell: Vec<f64>,
    /// Redistancing audits.
    pub audits: Vec<RedistanceAudit>,
    /// Nucleation events.
    pub events: Vec<NucleationEvent>,
    /// FNV-64 hashes of the φ bits per iteration (evolution snapshots).
    pub snapshots: Vec<u64>,
    /// Ledger rows.
    pub rows: Vec<String>,
}

fn fnv(phi: &GridSdf) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for v in phi.nodes() {
        for b in v.to_bits().to_le_bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

/// Material area of `{φ < 0}` by certified cut quadrature.
#[must_use]
pub fn material_volume(grid: &Quadtree, phi: &GridSdf) -> f64 {
    let mut vol = 0.0;
    for c in grid.leaves() {
        let (lo, hi) = grid.rect(c);
        let iv = phi.enclose(lo, hi);
        if iv.hi() < 0.0 {
            vol += (hi[0] - lo[0]) * (hi[1] - lo[1]);
        } else if iv.lo() <= 0.0 {
            vol += cut_cell_rules(phi, lo, hi, 2)
                .bulk
                .iter()
                .map(|&(_, w)| w)
                .sum::<f64>();
        }
    }
    vol
}

/// The cantilever fixture: clamp on the left box edge, downward
/// traction band on the right edge around mid-height.
#[derive(Debug, Clone, Copy)]
pub struct Cantilever {
    /// Traction magnitude.
    pub load: f64,
    /// Load band half-width.
    pub band: f64,
}

/// Uniform Q1 mass/stiffness on the full node lattice (Sobolev step).
fn mass_stiffness(n: usize) -> (fs_sparse::Csr, fs_sparse::Csr) {
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / n as f64;
    let stride = n + 1;
    let nn = stride * stride;
    let mut mc = fs_sparse::Coo::new(nn, nn);
    let mut kc = fs_sparse::Coo::new(nn, nn);
    // Q1 element matrices on a square of side h (standard closed
    // forms): lumped mass h²/4 per corner; stiffness pattern of the
    // Laplacian.
    let ke = [
        [2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0],
        [-1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0],
    ];
    for cj in 0..n {
        for ci in 0..n {
            let ids = [
                ci + cj * stride,
                ci + 1 + cj * stride,
                ci + 1 + (cj + 1) * stride,
                ci + (cj + 1) * stride,
            ];
            for (a, &ia) in ids.iter().enumerate() {
                mc.push(ia, ia, 0.25 * h * h);
                for (b, &ib) in ids.iter().enumerate() {
                    kc.push(ia, ib, ke[a][b]);
                }
            }
        }
    }
    (mc.assemble(), kc.assemble())
}

/// Run the level-set compliance descent. Returns the report; the
/// level set evolves in place.
///
/// # Errors
/// fs-solid solve refusals propagate.
///
/// # Panics
/// If `phi.n() != 2^level` (the SDF lattice must match the CutFEM
/// grid so cells align).
#[allow(clippy::too_many_lines)] // the descent loop is one narrative
pub fn optimize_compliance(
    phi: &mut GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
) -> Result<OptimizeReport, SolidError> {
    let n = 1usize << settings.level;
    assert_eq!(phi.n(), n, "SDF lattice must match the CutFEM grid");
    let grid = Quadtree::uniform(settings.level);
    let h = phi.h();
    let stride = n + 1;
    let (lambda, mu) = lame(settings.youngs, settings.poisson, PlaneKind::Strain);
    let (mass, stiffness) = mass_stiffness(n);
    let mut ell = settings.ell0;
    let mut report = OptimizeReport::default();
    let clamp = |x: f64, _y: f64| x < 1e-9;
    let load_band = fixture.band;
    let load = fixture.load;
    let traction = move |x: f64, y: f64| -> [f64; 2] {
        if x > 1.0 - 1e-9 && (y - 0.5).abs() <= load_band {
            [0.0, -load]
        } else {
            [0.0, 0.0]
        }
    };
    for iter in 0..settings.iterations {
        // 1. Physics on the level set (zero meshing).
        let solver = CutElasticity {
            grid: &grid,
            sdf: phi,
            youngs: settings.youngs,
            poisson: settings.poisson,
            nitsche_beta: 20.0,
            ghost_gamma: 0.5,
            quad_depth: 2,
            clamp: Some(&clamp),
            boundary_traction: Some(&traction),
            traction_free_interface: true,
        };
        let sol = solver.solve(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])?;
        // 2. Compliance J = ∫ t·u over the loaded edge.
        let mut compliance = 0.0;
        for (&(gi, gj), u) in sol.nodal() {
            if gi as usize == n {
                let y = f64::from(gj) * h;
                let t = traction(1.0, y);
                // Nodal quadrature of the edge integral (weight h,
                // halved at the ends).
                let w = if gj == 0 || gj as usize == n {
                    0.5 * h
                } else {
                    h
                };
                compliance += w * (t[0] * u[0] + t[1] * u[1]);
            }
        }
        // t·u is already positive: t = (0,−P) against uy < 0.
        // 3. Nodal strain-energy density w(x) = ½ σ:ε from adjacent
        // material; seeds for the extension are band nodes.
        let mut energy = vec![0.0f64; stride * stride];
        let mut seeded = vec![false; stride * stride];
        for j in 0..=n {
            for i in 0..=n {
                let k = i + j * stride;
                let p = phi.pos(i, j);
                // Sample slightly inside material along −∇φ.
                let g = phi.gradient_at(p);
                let gn = g[0].hypot(g[1]).max(1e-12);
                let q = [
                    (p[0] - 0.75 * h * g[0] / gn).clamp(0.0, 1.0),
                    (p[1] - 0.75 * h * g[1] / gn).clamp(0.0, 1.0),
                ];
                if phi.value_at(q) > 0.0 {
                    continue;
                }
                let (eps, ok) = strain_at(&grid, phi, &sol, q);
                if !ok {
                    continue;
                }
                let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
                let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
                let sxy = 2.0 * mu * eps[2];
                energy[k] = 0.5 * (sxx * eps[0] + syy * eps[1] + 2.0 * sxy * eps[2]);
                seeded[k] = phi.node(i, j).abs() <= 2.0 * h;
            }
        }
        // 4. Extend off the interface, smooth, form v_n = w − ℓ.
        extend_velocity(phi, &mut energy, &seeded);
        let (smooth, _iters) = fs_adjoint::sobolev::sobolev_smooth(
            &mass,
            &stiffness,
            settings.sobolev_alpha * h * h,
            &energy,
            1e-10,
        );
        // The multiplier lives on the ENERGY-DENSITY scale: its
        // update is normalized by the mean band energy so the volume
        // feedback competes with the shape term instead of drowning it
        // (an O(1) multiplier against O(J) energies shrinks the
        // structure to nothing at full speed — measured failure mode).
        let vn: Vec<f64> = smooth.iter().map(|w| w - ell).collect();
        // 5. Advect one interface move, on the band.
        let band = build_band(phi, settings.band_cells);
        let vmax = vn.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-12);
        advect(
            phi,
            &band,
            &Velocity::Normal(&vn),
            settings.move_cells * h / vmax,
            0.45,
        );
        // 6. Redistance + audit.
        let audit = redistance(phi, settings.band_cells);
        // 7. Volume + multiplier update.
        let volume = material_volume(&grid, phi);
        #[allow(clippy::cast_precision_loss)]
        let w_mean = smooth.iter().sum::<f64>() / smooth.len() as f64;
        ell = (ell
            + settings.mu_al * w_mean.abs().max(1e-30) * (volume - settings.volfrac)
                / settings.volfrac)
            .max(0.0);
        // 8. Scheduled nucleation by topological derivative.
        if settings.nucleation_period > 0 && iter > 0 && iter % settings.nucleation_period == 0 {
            let mut dt_field = vec![f64::INFINITY; stride * stride];
            for j in 0..=n {
                for i in 0..=n {
                    let k = i + j * stride;
                    let p = phi.pos(i, j);
                    if phi.value_at(p) > -2.0 * h {
                        continue;
                    }
                    let (eps, ok) = strain_at(&grid, phi, &sol, p);
                    if !ok {
                        continue;
                    }
                    let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
                    let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
                    let sxy = 2.0 * mu * eps[2];
                    dt_field[k] = topological_derivative(lambda, mu, [sxx, syy, sxy], eps);
                }
            }
            let events = nucleate(
                phi,
                &dt_field,
                ell,
                settings.hole_radius_cells * h,
                6.0 * settings.hole_radius_cells * h,
                2,
            );
            if !events.is_empty() {
                let _ = redistance(phi, settings.band_cells);
            }
            report.events.extend(events);
        }
        // 9. Ledger.
        let snap = fnv(phi);
        let mut row = String::new();
        let _ = write!(
            row,
            "{{\"iter\":{iter},\"compliance\":{compliance:.6e},\"volume\":{volume:.4},\
             \"ell\":{ell:.4e},\"drift_h\":{:.2e},\"snapshot\":\"{snap:#018x}\"}}",
            audit.interface_drift_h
        );
        report.rows.push(row);
        report.compliance.push(compliance);
        report.volume.push(volume);
        report.ell.push(ell);
        report.audits.push(audit);
        report.snapshots.push(snap);
    }
    Ok(report)
}

/// Strain at a point from the CutFEM solution (bilinear gradient on
/// the containing cell); `ok = false` outside the active mesh.
fn strain_at(
    grid: &Quadtree,
    phi: &GridSdf,
    sol: &fs_solid::CutSolution,
    p: [f64; 2],
) -> ([f64; 3], bool) {
    let level = grid.max_level();
    let nf = f64::from(1u32 << level);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let ci = ((p[0] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let cj = ((p[1] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    let cell = (level, ci, cj);
    let (lo, hi) = grid.rect(cell);
    let corners = grid.corner_nodes(cell);
    let nodal = sol.nodal();
    let mut vals = [[0.0f64; 2]; 4];
    for (a, c) in corners.iter().enumerate() {
        match nodal.get(c) {
            Some(u) => vals[a] = *u,
            None => return ([0.0; 3], false),
        }
    }
    let _ = phi;
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = ((p[0] - lo[0]) / hx).clamp(0.0, 1.0);
    let et = ((p[1] - lo[1]) / hy).clamp(0.0, 1.0);
    let g = [
        [-(1.0 - et) / hx, -(1.0 - xi) / hy],
        [(1.0 - et) / hx, -xi / hy],
        [et / hx, xi / hy],
        [-et / hx, (1.0 - xi) / hy],
    ];
    let mut gu = [[0.0f64; 2]; 2];
    for a in 0..4 {
        for c in 0..2 {
            gu[c][0] += g[a][0] * vals[a][c];
            gu[c][1] += g[a][1] * vals[a][c];
        }
    }
    (
        [gu[0][0], gu[1][1], f64::midpoint(gu[0][1], gu[1][0])],
        true,
    )
}
