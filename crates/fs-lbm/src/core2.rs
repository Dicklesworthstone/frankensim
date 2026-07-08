//! The general D2Q9 core (bead tfz.19): cell flags, VECTOR gravity
//! (tilt schedules rotate it), per-cell relaxation time (the
//! non-Newtonian hook), Guo forcing, pull streaming with halfway
//! bounce-back at walls — the substrate the thermal, rheology,
//! refinement, and free-surface extensions all share. Deterministic:
//! fixed row-major cell order, no RNG.

use crate::{CS2, E, OPP, Q, W};

/// Cell classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cell {
    /// Bulk fluid.
    Fluid,
    /// Solid wall (halfway bounce-back).
    Wall,
    /// Free-surface interface cell (carries partial mass).
    Interface,
    /// Gas cell (no populations).
    Gas,
}

/// The general D2Q9 lattice.
#[derive(Debug, Clone)]
pub struct Grid {
    /// Cells in x.
    pub nx: usize,
    /// Cells in y.
    pub ny: usize,
    /// Cell flags.
    pub flags: Vec<Cell>,
    /// Populations.
    pub f: Vec<[f64; Q]>,
    /// Per-cell relaxation time.
    pub tau: Vec<f64>,
    /// Gravity vector (lattice units).
    pub g: [f64; 2],
    /// Per-cell external force (Boussinesq buoyancy etc.), added to
    /// ρ·g.
    pub fext: Vec<[f64; 2]>,
    /// Periodic in x?
    pub periodic_x: bool,
    /// Periodic in y?
    pub periodic_y: bool,
}

/// Macroscopic moments of one cell.
#[derive(Debug, Clone, Copy)]
pub struct Moments {
    /// Density.
    pub rho: f64,
    /// Velocity (force-corrected).
    pub u: [f64; 2],
}

impl Grid {
    /// A grid of fluid at rest (unit density), uniform `tau`.
    #[must_use]
    pub fn uniform(nx: usize, ny: usize, tau: f64) -> Grid {
        assert!(nx > 0 && ny > 0, "grid dimensions must be positive");
        assert!(
            tau.is_finite() && tau > 0.5,
            "relaxation time tau must be finite and greater than 0.5"
        );
        let f0 = crate::equilibrium(1.0, 0.0, 0.0);
        Grid {
            nx,
            ny,
            flags: vec![Cell::Fluid; nx * ny],
            f: vec![f0; nx * ny],
            tau: vec![tau; nx * ny],
            g: [0.0, 0.0],
            fext: vec![[0.0; 2]; nx * ny],
            periodic_x: true,
            periodic_y: true,
        }
    }

    /// Row-major index.
    #[must_use]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Moments of cell `i` (Guo half-force correction).
    #[must_use]
    pub fn moments(&self, i: usize) -> Moments {
        let f = &self.f[i];
        let rho: f64 = f.iter().sum();
        assert!(
            rho.is_finite() && rho > 0.0,
            "moments require positive finite density"
        );
        let mut m = [0.0f64; 2];
        for (q, fi) in f.iter().enumerate() {
            m[0] += f64::from(E[q].0) * fi;
            m[1] += f64::from(E[q].1) * fi;
        }
        let fx = self.g[0].mul_add(rho, self.fext[i][0]);
        let fy = self.g[1].mul_add(rho, self.fext[i][1]);
        Moments {
            rho,
            u: [(m[0] + 0.5 * fx) / rho, (m[1] + 0.5 * fy) / rho],
        }
    }

    /// Total mass over non-gas cells.
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.f
            .iter()
            .zip(&self.flags)
            .filter(|&(_, &fl)| fl != Cell::Gas && fl != Cell::Wall)
            .map(|(c, _)| c.iter().sum::<f64>())
            .sum()
    }

    /// Collide (per-cell tau, vector Guo forcing) into `post`.
    pub fn collide_into(&self, post: &mut Vec<[f64; Q]>) {
        post.clear();
        post.resize(self.nx * self.ny, [0.0; Q]);
        for (i, out) in post.iter_mut().enumerate().take(self.nx * self.ny) {
            if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                *out = self.f[i];
                continue;
            }
            let mm = self.moments(i);
            let (rho, ux, uy) = (mm.rho, mm.u[0], mm.u[1]);
            let feq = crate::equilibrium(rho, ux, uy);
            let tau = self.tau[i];
            assert!(
                tau.is_finite() && tau > 0.5,
                "cell relaxation time tau must be finite and greater than 0.5"
            );
            let coef = 1.0 - 0.5 / tau;
            let (gx, gy) = (
                self.g[0].mul_add(rho, self.fext[i][0]),
                self.g[1].mul_add(rho, self.fext[i][1]),
            );
            for q in 0..Q {
                let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
                let eu = ex * ux + ey * uy;
                // Guo forcing, vector form.
                let fx = (ex - ux) / CS2 + eu * ex / (CS2 * CS2);
                let fy = (ey - uy) / CS2 + eu * ey / (CS2 * CS2);
                let force = coef * W[q] * (fx * gx + fy * gy);
                out[q] = self.f[i][q] + (feq[q] - self.f[i][q]) / tau + force;
            }
        }
    }

    /// Source cell for pull-streaming direction `q` into (x, y);
    /// `None` when the pull crosses a non-periodic boundary (treated
    /// as wall bounce-back).
    #[must_use]
    pub fn source(&self, x: usize, y: usize, q: usize) -> Option<usize> {
        let (ex, ey) = E[q];
        let sx = match ex {
            1 => {
                if x == 0 {
                    if self.periodic_x {
                        self.nx - 1
                    } else {
                        return None;
                    }
                } else {
                    x - 1
                }
            }
            -1 => {
                if x + 1 == self.nx {
                    if self.periodic_x {
                        0
                    } else {
                        return None;
                    }
                } else {
                    x + 1
                }
            }
            _ => x,
        };
        let sy = match ey {
            1 => {
                if y == 0 {
                    if self.periodic_y {
                        self.ny - 1
                    } else {
                        return None;
                    }
                } else {
                    y - 1
                }
            }
            -1 => {
                if y + 1 == self.ny {
                    if self.periodic_y {
                        0
                    } else {
                        return None;
                    }
                } else {
                    y + 1
                }
            }
            _ => y,
        };
        Some(self.idx(sx, sy))
    }

    /// Stream `post` into `self.f` (fluid pull; wall and out-of-domain
    /// pulls bounce back).
    pub fn stream_from(&mut self, post: &[[f64; Q]]) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                for q in 0..Q {
                    let pulled = match self.source(x, y, q) {
                        Some(s) if matches!(self.flags[s], Cell::Wall | Cell::Gas) => {
                            post[i][OPP[q]]
                        }
                        Some(s) => post[s][q],
                        None => post[i][OPP[q]],
                    };
                    self.f[i][q] = pulled;
                }
            }
        }
    }

    /// One plain step (no free-surface bookkeeping).
    pub fn step(&mut self, scratch: &mut Vec<[f64; Q]>) {
        self.collide_into(scratch);
        let post = std::mem::take(scratch);
        self.stream_from(&post);
        *scratch = post;
    }
}

/// Strain-rate magnitude (sqrt(2 S:S)) of one cell from its
/// non-equilibrium populations — the LOCAL quantity non-Newtonian
/// relaxation adapts to. `feq` must match the cell's moments.
#[must_use]
pub fn shear_rate(f: &[f64; Q], feq: &[f64; Q], rho: f64, tau: f64) -> f64 {
    // S_ab = −(3 / (2 ρ τ)) Σ_q e_qa e_qb (f_q − feq_q)   (c_s² = 1/3)
    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    let mut syy = 0.0f64;
    for q in 0..Q {
        let neq = f[q] - feq[q];
        let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
        sxx += ex * ex * neq;
        sxy += ex * ey * neq;
        syy += ey * ey * neq;
    }
    let c = -3.0 / (2.0 * rho * tau);
    let (sxx, sxy, syy) = (c * sxx, c * sxy, c * syy);
    let ss = 2.0f64.mul_add(sxy * sxy, sxx.mul_add(sxx, syy * syy));
    (2.0 * ss).sqrt()
}
