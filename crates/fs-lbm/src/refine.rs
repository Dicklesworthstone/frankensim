//! Grid refinement with rescaling (bead tfz.19): a two-level 1:2
//! coupled channel — coarse below the split plane, fine (2× in space,
//! 2× in time) above it — with the Dupuis–Chopard level-transition
//! physics: equilibria transfer as-is under convective scaling (the
//! lattice velocity is scale-invariant), NON-equilibria rescale by
//! the τ·dt ratio (fine←coarse: ×τ_f/(2τ_c); coarse←fine: the
//! inverse), and the viscosity match fixes τ_f = 2τ_c − ½. Ghost
//! layers carry PRE-collision states and participate in the receiving
//! grid's collision, so transferred populations behave exactly like
//! native fluid for the streaming step; fine sub-step ghosts use
//! linear temporal interpolation between coarse times.

use crate::core2::{Cell, Grid};
use crate::{Q, equilibrium};

/// A channel refined above the split: coarse grid owns the lower half
/// (wall at the bottom), fine grid owns the upper half (wall at the
/// top), periodic in x.
pub struct RefinedChannel {
    /// Coarse grid: row 0 wall, rows 1..=own fluid, row own+1 ghost.
    pub coarse: Grid,
    /// Fine grid: rows 0..2 ghosts, rows 2..2+own fluid, top wall.
    pub fine: Grid,
    /// Coarse relaxation time.
    pub tau_c: f64,
    /// Fine relaxation time (2τ_c − ½).
    pub tau_f: f64,
    /// Coarse fluid rows owned.
    pub own_c: usize,
    /// Fine fluid rows owned.
    pub own_f: usize,
    scratch: Vec<[f64; Q]>,
}

/// Rescale one population set between levels: keep the equilibrium,
/// scale the non-equilibrium by `factor`.
fn rescale(f: &[f64; Q], factor: f64) -> [f64; Q] {
    let rho: f64 = f.iter().sum();
    let mut m = [0.0f64; 2];
    for (q, fi) in f.iter().enumerate() {
        m[0] += f64::from(crate::E[q].0) * fi;
        m[1] += f64::from(crate::E[q].1) * fi;
    }
    let eq = equilibrium(rho, m[0] / rho, m[1] / rho);
    let mut out = [0.0f64; Q];
    for q in 0..Q {
        out[q] = factor.mul_add(f[q] - eq[q], eq[q]);
    }
    out
}

impl RefinedChannel {
    /// Build with `nx_c` coarse columns, `own_c` coarse fluid rows
    /// below the split and `own_f` fine fluid rows above it, gravity
    /// `gx` along the channel.
    #[must_use]
    pub fn new(nx_c: usize, own_c: usize, own_f: usize, tau_c: f64, gx: f64) -> RefinedChannel {
        let tau_f = 2.0f64.mul_add(tau_c, -0.5);
        // Coarse: wall + own + ghost.
        let mut coarse = Grid::uniform(nx_c, own_c + 2, tau_c);
        coarse.periodic_y = false;
        coarse.g = [gx, 0.0];
        for x in 0..nx_c {
            let w = coarse.idx(x, 0);
            coarse.flags[w] = Cell::Wall;
        }
        // Fine: 2 ghosts + own + wall. Convective scaling halves the
        // body force per fine step? No: acceleration in lattice units
        // g_lat = g_phys·dt²/dx; dt and dx both halve, so g_lat
        // halves on the fine grid.
        let mut fine = Grid::uniform(2 * nx_c, own_f + 3, tau_f);
        fine.periodic_y = false;
        fine.g = [gx / 2.0, 0.0];
        for x in 0..2 * nx_c {
            let w = fine.idx(x, own_f + 2);
            fine.flags[w] = Cell::Wall;
        }
        RefinedChannel {
            coarse,
            fine,
            tau_c,
            tau_f,
            own_c,
            own_f,
            scratch: Vec::new(),
        }
    }

    /// Prolong coarse → fine ghost populations (rows 0 and 1),
    /// linear in x and y, non-equilibrium rescaled.
    fn prolong(&self) -> Vec<[f64; Q]> {
        let nxf = self.fine.nx;
        let nxc = self.coarse.nx;
        let factor = self.tau_f / (2.0 * self.tau_c);
        // Coarse rows: own_c (top owned, center Ys−dx_f in fine
        // units), own_c−1 (center Ys−3dx_f), ghost row own_c+1
        // (center Ys+dx_f).
        let rows = [self.own_c - 1, self.own_c, self.own_c + 1];
        let mut out = vec![[0.0f64; Q]; 2 * nxf];
        for xf in 0..nxf {
            let xc = xf / 2;
            // x-linear neighbor with weight 0.25.
            let xn = if xf % 2 == 0 {
                (xc + nxc - 1) % nxc
            } else {
                (xc + 1) % nxc
            };
            let blend_x = |row: usize, buf: &mut [f64; Q]| {
                let a = &self.coarse.f[self.coarse.idx(xc, row)];
                let b = &self.coarse.f[self.coarse.idx(xn, row)];
                for q in 0..Q {
                    buf[q] = 0.75f64.mul_add(a[q], 0.25 * b[q]);
                }
            };
            let mut low = [0.0f64; Q];
            let mut mid = [0.0f64; Q];
            let mut high = [0.0f64; Q];
            blend_x(rows[0], &mut low);
            blend_x(rows[1], &mut mid);
            blend_x(rows[2], &mut high);
            // Fine ghost k=0 (center Ys−1.5dx_f): 0.25 low + 0.75 mid.
            // Fine ghost k=1 (center Ys−0.5dx_f): 0.75 mid + 0.25 high.
            let mut g0 = [0.0f64; Q];
            let mut g1 = [0.0f64; Q];
            for q in 0..Q {
                g0[q] = 0.25f64.mul_add(low[q], 0.75 * mid[q]);
                g1[q] = 0.75f64.mul_add(mid[q], 0.25 * high[q]);
            }
            out[xf] = rescale(&g0, factor);
            out[nxf + xf] = rescale(&g1, factor);
        }
        out
    }

    fn set_fine_ghosts(&mut self, ghosts: &[[f64; Q]]) {
        let nxf = self.fine.nx;
        for xf in 0..nxf {
            let i0 = self.fine.idx(xf, 0);
            self.fine.f[i0] = ghosts[xf];
            let i1 = self.fine.idx(xf, 1);
            self.fine.f[i1] = ghosts[nxf + xf];
        }
    }

    /// Restrict fine rows 2, 3 → the coarse ghost row (average of the
    /// 2×2 block, non-equilibrium rescaled).
    fn restrict(&mut self) {
        let nxc = self.coarse.nx;
        let factor = 2.0 * self.tau_c / self.tau_f;
        for xc in 0..nxc {
            let mut acc = [0.0f64; Q];
            for (dx, row) in [(0usize, 2usize), (1, 2), (0, 3), (1, 3)] {
                let f = &self.fine.f[self.fine.idx(2 * xc + dx, row)];
                let r = rescale(f, factor);
                for q in 0..Q {
                    acc[q] += 0.25 * r[q];
                }
            }
            let gi = self.coarse.idx(xc, self.own_c + 1);
            self.coarse.f[gi] = acc;
        }
    }

    /// One coarse step (two fine sub-steps) with the full transfer
    /// cycle.
    pub fn step(&mut self) {
        let g_t = self.prolong();
        // The coarse ghost row streams against a phantom ceiling
        // during the coarse step (bounce-back garbage); snapshot and
        // restore it so the t+1 prolongation sees the last RESTRICTED
        // value (first-order lag) instead — without this the 0.25
        // ghost weight injects wall-like drag at the interface every
        // cycle (measured as a 12% global Poiseuille deficit).
        let ghost_backup: Vec<[f64; Q]> = (0..self.coarse.nx)
            .map(|x| self.coarse.f[self.coarse.idx(x, self.own_c + 1)])
            .collect();
        self.coarse.step(&mut self.scratch);
        for (x, g) in ghost_backup.into_iter().enumerate() {
            let i = self.coarse.idx(x, self.own_c + 1);
            self.coarse.f[i] = g;
        }
        let g_t1 = self.prolong();
        let g_half: Vec<[f64; Q]> = g_t
            .iter()
            .zip(&g_t1)
            .map(|(a, b)| {
                let mut m = [0.0f64; Q];
                for q in 0..Q {
                    m[q] = f64::midpoint(a[q], b[q]);
                }
                m
            })
            .collect();
        self.set_fine_ghosts(&g_t);
        self.fine.step(&mut self.scratch);
        self.set_fine_ghosts(&g_half);
        self.fine.step(&mut self.scratch);
        self.restrict();
    }

    /// The x-velocity profile across the WHOLE channel at physical
    /// resolution dx_f: coarse rows report their center value twice.
    /// Returns (y_phys in fine units, u_x) pairs, walls excluded.
    #[must_use]
    pub fn profile(&self) -> Vec<(f64, f64)> {
        let mut out = Vec::new();
        for j in 1..=self.own_c {
            let mm = self.coarse.moments(self.coarse.idx(0, j));
            // Coarse row j center in fine units: 2(j − 0.5) = 2j − 1.
            out.push(((2 * j) as f64 - 1.0, mm.u[0]));
        }
        for k in 2..2 + self.own_f {
            let mm = self.fine.moments(self.fine.idx(0, k));
            // Fine row k center: split + (k − 2) + 0.5.
            let split = 2.0 * self.own_c as f64;
            out.push((split + (k as f64 - 2.0) + 0.5, mm.u[0]));
        }
        out
    }

    /// Total kinetic energy, area-weighted (coarse cells count 4×).
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0f64;
        for j in 1..=self.own_c {
            for x in 0..self.coarse.nx {
                let mm = self.coarse.moments(self.coarse.idx(x, j));
                ke += 4.0 * 0.5 * mm.rho * mm.u[1].mul_add(mm.u[1], mm.u[0] * mm.u[0]);
            }
        }
        for k in 2..2 + self.own_f {
            for x in 0..self.fine.nx {
                let mm = self.fine.moments(self.fine.idx(x, k));
                ke += 0.5 * mm.rho * mm.u[1].mul_add(mm.u[1], mm.u[0] * mm.u[0]);
            }
        }
        ke
    }

    /// Seed a transverse shear wave u_y = a·sin(2πx/L) on both grids
    /// (same PHYSICAL field; lattice velocity is scale-invariant under
    /// convective scaling).
    pub fn seed_shear_wave(&mut self, a: f64) {
        let nxc = self.coarse.nx;
        for j in 0..self.coarse.ny {
            for x in 0..nxc {
                let i = self.coarse.idx(x, j);
                if self.coarse.flags[i] != Cell::Wall {
                    let uy = a * fs_math::det::sin(
                        std::f64::consts::TAU * (x as f64 + 0.5) / nxc as f64,
                    );
                    self.coarse.f[i] = equilibrium(1.0, 0.0, uy);
                }
            }
        }
        let nxf = self.fine.nx;
        for k in 0..self.fine.ny {
            for x in 0..nxf {
                let i = self.fine.idx(x, k);
                if self.fine.flags[i] != Cell::Wall {
                    // Same physical phase: fine cell x covers phase
                    // (x + 0.5)/nxf of the same period.
                    let uy = a * fs_math::det::sin(
                        std::f64::consts::TAU * (x as f64 + 0.5) / nxf as f64,
                    );
                    self.fine.f[i] = equilibrium(1.0, 0.0, uy);
                }
            }
        }
    }
}
