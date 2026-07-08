//! Free-surface LBM (bead tfz.19): Körner-lineage mass-tracking VOF.
//! Interface cells carry a partial mass m (fill fraction κ = m/ρ);
//! mass moves along streaming links with PAIRWISE-ANTISYMMETRIC
//! exchange terms (so the ledger balances exactly by construction),
//! missing populations from gas neighbors are reconstructed with the
//! atmospheric-pressure anti-population closure (optionally shifted
//! by a Laplace term from fill-field curvature), and cell conversions
//! (fluid ↔ interface ↔ gas) redistribute their excess/deficit mass
//! conservatively — with a carry accumulator so NOTHING is ever
//! silently dropped. The battery gates the ledger at 1e-10 relative
//! EVERY step; contact-line physics is deliberately MODEL-BRACKETED
//! (neutral vs wetting wall ghosts), not pretended-certain.

use crate::core2::{Cell, Grid};
use crate::{CS2, E, OPP, Q, equilibrium};

/// Wall wetting model for the curvature fill-field ghost — the
/// contact-line bracket per the plan's honesty clause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactModel {
    /// Neutral (≈90°): wall ghost copies the adjacent cell's fill.
    Neutral,
    /// Wetting: wall ghost is full (φ = 1) — the fluid is drawn along
    /// the wall.
    Wetting,
}

/// Free-surface simulation state.
pub struct FreeSurface {
    /// The lattice (flags: Fluid/Interface/Gas/Wall).
    pub grid: Grid,
    /// Per-cell tracked mass (meaningful for Interface; Fluid uses
    /// Σf, Gas is 0).
    pub mass: Vec<f64>,
    /// Surface tension coefficient (0 = off).
    pub sigma: f64,
    /// Contact-line model (curvature ghost at walls).
    pub contact: ContactModel,
    /// Conversion-conservation carry (mass awaiting redistribution).
    pub carry: f64,
    /// Cell-conversion statistics (cumulative).
    pub conversions: ConversionStats,
    scratch: Vec<[f64; Q]>,
    fill_smooth: Vec<f64>,
}

/// Cumulative conversion ledger.
#[derive(Debug, Clone, Copy, Default)]
pub struct ConversionStats {
    /// Interface → fluid events.
    pub to_fluid: u64,
    /// Interface → gas events.
    pub to_gas: u64,
    /// Gas → interface events (closure repair).
    pub gas_to_interface: u64,
    /// Fluid → interface events (closure repair).
    pub fluid_to_interface: u64,
}

const EPS: f64 = 1e-3;

impl FreeSurface {
    /// Build from a grid whose flags are already set (Fluid regions,
    /// Gas elsewhere, Wall boundaries). Interface cells are inserted
    /// automatically between fluid and gas; masses initialized
    /// (fluid = ρ, interface = ρ/2).
    ///
    /// # Panics
    /// If a fluid cell touches a gas cell after interface insertion
    /// (impossible by construction).
    #[must_use]
    pub fn new(mut grid: Grid, sigma: f64, contact: ContactModel) -> FreeSurface {
        let (nx, ny) = (grid.nx, grid.ny);
        // Insert interface cells: any fluid cell with a gas neighbor.
        let mut promote = Vec::new();
        for y in 0..ny {
            for x in 0..nx {
                let i = grid.idx(x, y);
                if grid.flags[i] != Cell::Fluid {
                    continue;
                }
                for q in 1..Q {
                    if let Some(nb) = neighbor(&grid, x, y, q)
                        && grid.flags[nb] == Cell::Gas
                    {
                        promote.push(i);
                        break;
                    }
                }
            }
        }
        for i in promote {
            grid.flags[i] = Cell::Interface;
        }
        let mut mass = vec![0.0f64; nx * ny];
        for (i, m) in mass.iter_mut().enumerate().take(nx * ny) {
            match grid.flags[i] {
                Cell::Fluid => *m = grid.f[i].iter().sum(),
                Cell::Interface => *m = 0.5 * grid.f[i].iter().sum::<f64>(),
                _ => {}
            }
        }
        let fs = FreeSurface {
            grid,
            mass,
            sigma,
            contact,
            carry: 0.0,
            conversions: ConversionStats::default(),
            scratch: Vec::new(),
            fill_smooth: vec![0.0; nx * ny],
        };
        fs.assert_closure();
        fs
    }

    fn assert_closure(&self) {
        for y in 0..self.grid.ny {
            for x in 0..self.grid.nx {
                let i = self.grid.idx(x, y);
                if self.grid.flags[i] != Cell::Fluid {
                    continue;
                }
                for q in 1..Q {
                    if let Some(nb) = neighbor(&self.grid, x, y, q) {
                        assert!(
                            self.grid.flags[nb] != Cell::Gas,
                            "closure violated: fluid touches gas at ({x},{y})"
                        );
                    }
                }
            }
        }
    }

    /// Fill fraction of a cell (fluid 1, gas 0, interface m/ρ).
    #[must_use]
    pub fn fill(&self, i: usize) -> f64 {
        match self.grid.flags[i] {
            Cell::Fluid => 1.0,
            Cell::Interface => {
                let rho: f64 = self.grid.f[i].iter().sum();
                (self.mass[i] / rho.max(1e-12)).clamp(0.0, 1.0)
            }
            _ => 0.0,
        }
    }

    /// The strict ledger: Σ_fluid Σf + Σ_interface m + carry.
    #[must_use]
    pub fn ledger_mass(&self) -> f64 {
        let mut total = self.carry;
        for i in 0..self.grid.nx * self.grid.ny {
            match self.grid.flags[i] {
                Cell::Fluid => total += self.grid.f[i].iter().sum::<f64>(),
                Cell::Interface => total += self.mass[i],
                _ => {}
            }
        }
        total
    }

    /// Count of connected fluid+interface components (4-connectivity)
    /// — the breaking-jet fragment counter.
    #[must_use]
    pub fn fragment_count(&self) -> usize {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let wet = |i: usize| matches!(self.grid.flags[i], Cell::Fluid | Cell::Interface);
        let mut seen = vec![false; nx * ny];
        let mut count = 0;
        for start in 0..nx * ny {
            if !wet(start) || seen[start] {
                continue;
            }
            count += 1;
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(i) = stack.pop() {
                let (x, y) = (i % nx, i / nx);
                for q in [1usize, 2, 3, 4] {
                    if let Some(nb) = neighbor(&self.grid, x, y, q)
                        && wet(nb)
                        && !seen[nb]
                    {
                        seen[nb] = true;
                        stack.push(nb);
                    }
                }
            }
        }
        count
    }

    /// Smoothed fill field + curvature-adjusted reference density for
    /// interface reconstruction at cell i.
    fn reference_density(&self, x: usize, y: usize) -> f64 {
        if self.sigma == 0.0 {
            return 1.0;
        }
        let kappa = self.curvature(x, y);
        self.sigma.mul_add(kappa / CS2, 1.0)
    }

    /// Fill-field value with the wall ghost per contact model.
    fn phi_at(&self, i: usize, from: usize) -> f64 {
        match self.grid.flags[i] {
            Cell::Wall => match self.contact {
                ContactModel::Neutral => self.fill(from),
                ContactModel::Wetting => 1.0,
            },
            _ => self.fill_smooth[i],
        }
    }

    /// Curvature of the smoothed fill field at (x, y): div(n̂) with
    /// n̂ = −∇φ/|∇φ| (outward from fluid), central differences.
    fn curvature(&self, x: usize, y: usize) -> f64 {
        let g = &self.grid;
        let i = g.idx(x, y);
        // Normal components at the four face neighbors via one-sided
        // gradients of the smoothed field; divergence by differencing.
        let phi = |dx: i32, dy: i32| -> f64 {
            let xx = offset_coord(x, dx, g.nx, g.periodic_x);
            let yy = offset_coord(y, dy, g.ny, g.periodic_y);
            match (xx, yy) {
                (Some(a), Some(b)) => self.phi_at(g.idx(a, b), i),
                _ => self.fill_smooth[i],
            }
        };
        let nhat = |dx: i32, dy: i32| -> [f64; 2] {
            let gx = (phi(dx + 1, dy) - phi(dx - 1, dy)) / 2.0;
            let gy = (phi(dx, dy + 1) - phi(dx, dy - 1)) / 2.0;
            let m = gx.hypot(gy).max(1e-9);
            [-gx / m, -gy / m]
        };
        let div = (nhat(1, 0)[0] - nhat(-1, 0)[0]) / 2.0 + (nhat(0, 1)[1] - nhat(0, -1)[1]) / 2.0;
        div.clamp(-1.0, 1.0)
    }

    fn refresh_fill(&mut self) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let raw: Vec<f64> = (0..nx * ny).map(|i| self.fill(i)).collect();
        for y in 0..ny {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                if self.grid.flags[i] == Cell::Wall {
                    self.fill_smooth[i] = raw[i];
                    continue;
                }
                let mut acc = raw[i];
                let mut count = 1.0;
                for q in 1..Q {
                    if let Some(nb) = neighbor(&self.grid, x, y, q)
                        && self.grid.flags[nb] != Cell::Wall
                    {
                        acc += raw[nb];
                        count += 1.0;
                    }
                }
                self.fill_smooth[i] = acc / count;
            }
        }
    }

    /// One free-surface step: collide, mass exchange + reconstruction
    /// + stream, conversion cascade with conservative redistribution.
    pub fn step(&mut self) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        self.refresh_fill();
        self.grid.collide_into(&mut self.scratch);
        let post = std::mem::take(&mut self.scratch);
        // Mass exchange (pairwise antisymmetric, from POST populations)
        // and streaming with gas reconstruction.
        let fills: Vec<f64> = (0..nx * ny).map(|i| self.fill(i)).collect();
        let mut new_f = self.grid.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                let flag = self.grid.flags[i];
                if !matches!(flag, Cell::Fluid | Cell::Interface) {
                    continue;
                }
                // Mass exchange along all links (interface only; fluid
                // cells' Σf tracks their mass through plain streaming).
                if flag == Cell::Interface {
                    let mut dm = 0.0f64;
                    for q in 1..Q {
                        if let Some(nb) = neighbor(&self.grid, x, y, q) {
                            let w = match self.grid.flags[nb] {
                                Cell::Fluid => 1.0,
                                Cell::Interface => f64::midpoint(fills[i], fills[nb]),
                                _ => 0.0,
                            };
                            if w > 0.0 {
                                dm += w * (post[nb][OPP[q]] - post[i][q]);
                            }
                        }
                    }
                    self.mass[i] += dm;
                }
                // Stream (pull) with reconstruction from gas sources.
                let mm = self.grid.moments(i);
                for q in 0..Q {
                    let src = self.grid.source(x, y, q);
                    new_f[i][q] = match src {
                        Some(s) if self.grid.flags[s] == Cell::Gas => {
                            let rho_ref = self.reference_density(x, y);
                            let eq = equilibrium(rho_ref, mm.u[0], mm.u[1]);
                            eq[q] + eq[OPP[q]] - post[i][OPP[q]]
                        }
                        Some(s) if self.grid.flags[s] == Cell::Wall => post[i][OPP[q]],
                        Some(s) => post[s][q],
                        None => post[i][OPP[q]],
                    };
                }
            }
        }
        self.grid.f = new_f;
        self.scratch = post;
        self.apply_conversions();
    }

    fn apply_conversions(&mut self) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        // Fluid mass is Σf; interface mass tracked. Conversions:
        let mut excess_pool = std::mem::take(&mut self.carry);
        let mut to_fluid = Vec::new();
        let mut to_gas = Vec::new();
        for i in 0..nx * ny {
            if self.grid.flags[i] != Cell::Interface {
                continue;
            }
            let rho: f64 = self.grid.f[i].iter().sum();
            if self.mass[i] > (1.0 + EPS) * rho {
                to_fluid.push(i);
            } else if self.mass[i] < -EPS * rho {
                to_gas.push(i);
            }
        }
        // Interface → fluid: excess to the pool; gas neighbors become
        // interface (closure).
        for &i in &to_fluid {
            let rho: f64 = self.grid.f[i].iter().sum();
            excess_pool += self.mass[i] - rho;
            self.grid.flags[i] = Cell::Fluid;
            self.mass[i] = rho;
            self.conversions.to_fluid += 1;
            let (x, y) = (i % nx, i / nx);
            for q in 1..Q {
                if let Some(nb) = neighbor(&self.grid, x, y, q)
                    && self.grid.flags[nb] == Cell::Gas
                {
                    // Initialize from the average of wet neighbors.
                    let (nbx, nby) = (nb % nx, nb / nx);
                    let mut rho_avg = 0.0;
                    let mut u_avg = [0.0f64; 2];
                    let mut cnt = 0.0;
                    for q2 in 1..Q {
                        if let Some(nn) = neighbor(&self.grid, nbx, nby, q2)
                            && matches!(self.grid.flags[nn], Cell::Fluid | Cell::Interface)
                        {
                            let m2 = self.grid.moments(nn);
                            rho_avg += m2.rho;
                            u_avg[0] += m2.u[0];
                            u_avg[1] += m2.u[1];
                            cnt += 1.0;
                        }
                    }
                    if cnt > 0.0 {
                        rho_avg /= cnt;
                        u_avg[0] /= cnt;
                        u_avg[1] /= cnt;
                    } else {
                        rho_avg = 1.0;
                    }
                    self.grid.f[nb] = equilibrium(rho_avg, u_avg[0], u_avg[1]);
                    self.grid.flags[nb] = Cell::Interface;
                    self.mass[nb] = 0.0;
                    self.conversions.gas_to_interface += 1;
                }
            }
        }
        // Interface → gas: deficit to the pool; fluid neighbors become
        // interface (their Σf IS their mass — ledger unchanged).
        for &i in &to_gas {
            if self.grid.flags[i] != Cell::Interface {
                continue; // may have been re-flagged by the cascade
            }
            excess_pool += self.mass[i];
            self.grid.flags[i] = Cell::Gas;
            self.mass[i] = 0.0;
            self.conversions.to_gas += 1;
            let (x, y) = (i % nx, i / nx);
            for q in 1..Q {
                if let Some(nb) = neighbor(&self.grid, x, y, q)
                    && self.grid.flags[nb] == Cell::Fluid
                {
                    self.grid.flags[nb] = Cell::Interface;
                    self.mass[nb] = self.grid.f[nb].iter().sum();
                    self.conversions.fluid_to_interface += 1;
                }
            }
        }
        // Conservative redistribution of the pool over interface cells.
        let interfaces: Vec<usize> = (0..nx * ny)
            .filter(|&i| self.grid.flags[i] == Cell::Interface)
            .collect();
        if interfaces.is_empty() {
            self.carry = excess_pool;
        } else {
            let share = excess_pool / interfaces.len() as f64;
            for &i in &interfaces {
                self.mass[i] += share;
            }
        }
    }
}

fn offset_coord(coord: usize, delta: i32, n: usize, periodic: bool) -> Option<usize> {
    assert!(n > 0, "grid dimension must be positive");
    let mut out = coord;
    if delta >= 0 {
        for _ in 0..delta {
            if out + 1 == n {
                if periodic {
                    out = 0;
                } else {
                    return None;
                }
            } else {
                out += 1;
            }
        }
    } else {
        for _ in 0..(-delta) {
            if out == 0 {
                if periodic {
                    out = n - 1;
                } else {
                    return None;
                }
            } else {
                out -= 1;
            }
        }
    }
    Some(out)
}

/// Neighbor cell index in direction q (None across non-periodic
/// boundaries).
fn neighbor(grid: &Grid, x: usize, y: usize, q: usize) -> Option<usize> {
    let (ex, ey) = E[q];
    let xx = offset_coord(x, ex, grid.nx, grid.periodic_x)?;
    let yy = offset_coord(y, ey, grid.ny, grid.periodic_y)?;
    Some(grid.idx(xx, yy))
}

/// A closed-box dam-break fixture: walls on all four sides, a fluid
/// column of `a × 2a` cells in the lower-left corner, gas elsewhere,
/// gravity `g` pointing down.
#[must_use]
pub fn dam_break(
    nx: usize,
    ny: usize,
    a: usize,
    g: f64,
    sigma: f64,
    contact: ContactModel,
) -> FreeSurface {
    let mut grid = Grid::uniform(nx, ny, 0.55);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = [0.0, -g];
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
    for y in 1..=(2 * a).min(ny - 2) {
        for x in 1..=a.min(nx - 2) {
            let i = grid.idx(x, y);
            grid.flags[i] = Cell::Fluid;
        }
    }
    FreeSurface::new(grid, sigma, contact)
}

/// Surge-front x position (rightmost wet cell in the bottom fluid
/// row), in cells from the left wall.
#[must_use]
pub fn surge_front(fs: &FreeSurface) -> usize {
    let mut front = 0;
    for x in 1..fs.grid.nx - 1 {
        let i = fs.grid.idx(x, 1);
        if matches!(fs.grid.flags[i], Cell::Fluid | Cell::Interface) {
            front = x;
        }
    }
    front
}
