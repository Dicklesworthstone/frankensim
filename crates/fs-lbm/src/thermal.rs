//! Thermal LBM (bead tfz.19): double population — the D2Q9 flow field
//! plus a D2Q5 temperature distribution advected by it, coupled back
//! through a Boussinesq buoyancy force. Fixed-temperature walls use
//! anti-bounce-back. The battery gates the Rayleigh–Bénard onset
//! bracket (decay below Ra_c ≈ 1708, growth above) — physics the
//! scheme cannot fake.

use crate::core2::{Cell, Grid};
use crate::{CS2, Q};

/// D2Q5 velocities (rest + 4 axis directions).
const E5: [(i32, i32); 5] = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)];
/// D2Q5 weights.
const W5: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];
/// D2Q5 opposites.
const OPP5: [usize; 5] = [0, 3, 4, 1, 2];

/// D2Q5 advection-diffusion equilibrium.
fn geq(t: f64, ux: f64, uy: f64) -> [f64; 5] {
    let mut g = [0.0f64; 5];
    for q in 0..5 {
        let eu = f64::from(E5[q].0) * ux + f64::from(E5[q].1) * uy;
        g[q] = W5[q] * t * (1.0 + eu / CS2);
    }
    g
}

/// Double-population thermal lattice: rigid walls top/bottom at fixed
/// temperatures, periodic in x, Boussinesq buoyancy on the flow.
pub struct ThermalLbm {
    /// The flow grid (walls in rows 0 and ny+1).
    pub grid: Grid,
    /// Temperature populations.
    pub gpop: Vec<[f64; 5]>,
    /// Thermal relaxation time (α = c_s²(τ_g − ½)).
    pub tau_g: f64,
    /// Bottom-wall temperature.
    pub t_bottom: f64,
    /// Top-wall temperature.
    pub t_top: f64,
    /// Buoyancy coefficient g·β (force = gβ(T − T_ref) ŷ).
    pub gbeta: f64,
    /// Reference temperature.
    pub t_ref: f64,
    scratch: Vec<[f64; Q]>,
}

impl ThermalLbm {
    /// A quiescent conducting slab: `nx × ny` fluid rows between two
    /// wall rows, linear conduction profile as the initial state.
    #[must_use]
    pub fn slab(nx: usize, ny: usize, tau_f: f64, tau_g: f64, gbeta: f64) -> ThermalLbm {
        assert!(
            nx > 0 && ny >= 2,
            "thermal slab needs positive width and at least two fluid rows"
        );
        assert!(
            tau_g.is_finite() && tau_g > 0.5,
            "thermal relaxation time tau_g must be finite and greater than 0.5"
        );
        assert!(gbeta.is_finite(), "buoyancy coefficient must be finite");
        let grid_ny = ny + 2;
        let mut grid = Grid::uniform(nx, grid_ny, tau_f);
        grid.periodic_y = false;
        for x in 0..nx {
            let b = grid.idx(x, 0);
            grid.flags[b] = Cell::Wall;
            let t = grid.idx(x, grid_ny - 1);
            grid.flags[t] = Cell::Wall;
        }
        let (t_bottom, t_top) = (1.0f64, 0.0f64);
        let mut gpop = vec![[0.0f64; 5]; nx * grid_ny];
        for x in 0..nx {
            gpop[grid.idx(x, 0)] = geq(t_bottom, 0.0, 0.0);
            gpop[grid.idx(x, grid_ny - 1)] = geq(t_top, 0.0, 0.0);
        }
        for y in 1..=ny {
            // Linear conduction profile between the halfway wall planes.
            let t = t_bottom + (t_top - t_bottom) * ((y as f64 - 0.5) / ny as f64);
            for x in 0..nx {
                gpop[y * nx + x] = geq(t, 0.0, 0.0);
            }
        }
        ThermalLbm {
            grid,
            gpop,
            tau_g,
            t_bottom,
            t_top,
            gbeta,
            t_ref: 0.5,
            scratch: Vec::new(),
        }
    }

    /// Temperature of cell (x, y).
    #[must_use]
    pub fn temperature(&self, x: usize, y: usize) -> f64 {
        self.gpop[self.grid.idx(x, y)].iter().sum()
    }

    /// Thermal diffusivity α = c_s²(τ_g − ½).
    #[must_use]
    pub fn diffusivity(&self) -> f64 {
        CS2 * (self.tau_g - 0.5)
    }

    /// Seed a sinusoidal vertical-velocity perturbation (onset mode).
    pub fn perturb(&mut self, amplitude: f64) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        for y in 1..ny - 1 {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                let s = (std::f64::consts::TAU * x as f64 / nx as f64).sin()
                    * (std::f64::consts::PI * (y as f64 - 0.5) / (ny as f64 - 2.0)).sin();
                let mm = self.grid.moments(i);
                self.grid.f[i] = crate::equilibrium(mm.rho, mm.u[0], amplitude.mul_add(s, mm.u[1]));
            }
        }
    }

    /// One coupled step: buoyancy from T, flow step, then temperature
    /// collide + stream (anti-bounce-back at the fixed-T walls).
    pub fn step(&mut self) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        // Boussinesq: per-cell force via the external-force field.
        for y in 1..ny - 1 {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                let t: f64 = self.gpop[i].iter().sum();
                self.grid.fext[i] = [0.0, self.gbeta * (t - self.t_ref)];
            }
        }
        self.grid.step(&mut self.scratch);
        // Temperature collide.
        let mut post = vec![[0.0f64; 5]; nx * ny];
        for y in 1..ny - 1 {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                let mm = self.grid.moments(i);
                let t: f64 = self.gpop[i].iter().sum();
                let eq = geq(t, mm.u[0], mm.u[1]);
                for q in 0..5 {
                    post[i][q] = self.gpop[i][q] + (eq[q] - self.gpop[i][q]) / self.tau_g;
                }
            }
        }
        // Temperature stream with anti-bounce-back walls.
        for y in 1..ny - 1 {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                for q in 0..5 {
                    let (ex, ey) = E5[q];
                    let sx = match ex {
                        1 => (x + nx - 1) % nx,
                        -1 => (x + 1) % nx,
                        _ => x,
                    };
                    let sy_opt = match ey {
                        1 => Some(y - 1),
                        -1 => Some(y + 1),
                        _ => Some(y),
                    };
                    let s = sy_opt.map(|sy| self.grid.idx(sx, sy));
                    self.gpop[i][q] = match s {
                        Some(si) if self.grid.flags[si] == Cell::Wall => {
                            let tw = if y == 1 && ey == 1 {
                                self.t_bottom
                            } else {
                                self.t_top
                            };
                            // Anti-bounce-back: fixed halfway temperature.
                            2.0f64.mul_add(W5[q] * tw, -post[i][OPP5[q]])
                        }
                        Some(si) => post[si][q],
                        None => post[i][OPP5[q]],
                    };
                }
            }
        }
    }

    /// Total kinetic energy of the fluid.
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0f64;
        for i in 0..self.grid.nx * self.grid.ny {
            if self.grid.flags[i] == Cell::Fluid {
                let mm = self.grid.moments(i);
                ke += 0.5 * mm.rho * mm.u[1].mul_add(mm.u[1], mm.u[0] * mm.u[0]);
            }
        }
        ke
    }

    /// Nusselt number: 1 + ⟨u_y·T⟩·H / (α·ΔT).
    #[must_use]
    pub fn nusselt(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let h = (ny - 2) as f64;
        let mut adv = 0.0f64;
        let mut count = 0usize;
        for y in 1..ny - 1 {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                let t: f64 = self.gpop[i].iter().sum();
                adv += self.grid.moments(i).u[1] * t;
                count += 1;
            }
        }
        adv /= count as f64;
        1.0 + adv * h / (self.diffusivity() * (self.t_bottom - self.t_top))
    }

    /// The Rayleigh number of the current configuration.
    #[must_use]
    pub fn rayleigh(&self) -> f64 {
        let h = (self.grid.ny - 2) as f64;
        let nu = CS2 * (self.grid.tau[self.grid.idx(0, 1)] - 0.5);
        self.gbeta * (self.t_bottom - self.t_top) * h * h * h / (nu * self.diffusivity())
    }
}

/// gβ needed for a target Rayleigh number at the given lattice setup.
#[must_use]
pub fn gbeta_for_rayleigh(ra: f64, ny: usize, tau_f: f64, tau_g: f64) -> f64 {
    assert!(ra.is_finite(), "Rayleigh number must be finite");
    assert!(ny > 0, "Rayleigh height must be positive");
    assert!(
        tau_f.is_finite() && tau_f > 0.5,
        "flow relaxation time tau_f must be finite and greater than 0.5"
    );
    assert!(
        tau_g.is_finite() && tau_g > 0.5,
        "thermal relaxation time tau_g must be finite and greater than 0.5"
    );
    let h = ny as f64;
    let nu = CS2 * (tau_f - 0.5);
    let alpha = CS2 * (tau_g - 0.5);
    ra * nu * alpha / (h * h * h)
}
