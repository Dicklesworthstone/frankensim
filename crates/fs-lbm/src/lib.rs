//! fs-lbm — lattice Boltzmann core (D2Q9). Layer: L3.
//!
//! Lattice Boltzmann is the many-core queen. This v0 is the D2Q9 BGK
//! stream-and-collide kernel with Guo body forcing and halfway bounce-back
//! walls — enough to reproduce analytic Poiseuille channel flow (a parabolic
//! velocity profile) from first principles — plus the LATTICE-SCALING
//! ASSISTANT.
//!
//! The scaling assistant ([`ScalingPlan`]) automates the `dx`/`dt`/`τ`/`Mach`
//! bookkeeping that is a chronic source of human LBM error: given a Reynolds
//! number, a resolution, and a lattice velocity, it derives the relaxation time
//! `τ = 3ν + ½` and checks the stability constraints (`τ > ½`, low Mach) with
//! explicit margins — an Evidence-typed plan agents consume instead of touching
//! the raw knobs. Deterministic (fixed cell order).

pub use fs_evidence::Color;

pub mod core2;
pub mod freesurface;
pub mod refine;
pub mod rheology;
pub mod thermal;

pub use core2::{Cell, Grid};
pub use freesurface::{ContactModel, FreeSurface};
pub use refine::RefinedChannel;
pub use rheology::Rheology;
pub use thermal::ThermalLbm;

/// The D2Q9 population count.
pub const Q: usize = 9;

/// The D2Q9 lattice velocities.
pub(crate) const E: [(i32, i32); Q] = [
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

/// The D2Q9 lattice weights.
pub(crate) const W: [f64; Q] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Opposite-direction indices (for bounce-back).
pub(crate) const OPP: [usize; Q] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// The lattice sound speed squared `c_s² = 1/3`.
pub const CS2: f64 = 1.0 / 3.0;

/// The D2Q9 equilibrium distribution at density `rho` and velocity `(ux, uy)`.
#[must_use]
pub fn equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; Q] {
    let usq = ux * ux + uy * uy;
    let mut f = [0.0; Q];
    for i in 0..Q {
        let (ex, ey) = (f64::from(E[i].0), f64::from(E[i].1));
        let eu = ex * ux + ey * uy;
        f[i] = W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq);
    }
    f
}

/// A D2Q9 lattice Boltzmann channel: periodic in x, walls (bounce-back) in y,
/// driven by a body force in x.
#[derive(Debug, Clone)]
pub struct Lbm {
    nx: usize,
    ny: usize,
    tau: f64,
    gx: f64,
    f: Vec<[f64; Q]>,
}

impl Lbm {
    /// A channel at rest (unit density) with relaxation time `tau` and body
    /// force `gx`.
    #[must_use]
    pub fn channel(nx: usize, ny: usize, tau: f64, gx: f64) -> Lbm {
        let f0 = equilibrium(1.0, 0.0, 0.0);
        Lbm {
            nx,
            ny,
            tau,
            gx,
            f: vec![f0; nx * ny],
        }
    }

    fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// The kinematic viscosity `ν = (τ − ½)/3`.
    #[must_use]
    pub fn viscosity(&self) -> f64 {
        (self.tau - 0.5) / 3.0
    }

    /// The macroscopic density at `(x, y)`.
    #[must_use]
    pub fn density(&self, x: usize, y: usize) -> f64 {
        self.f[self.idx(x, y)].iter().sum()
    }

    /// The macroscopic velocity at `(x, y)` (with the body-force correction).
    #[must_use]
    pub fn velocity(&self, x: usize, y: usize) -> (f64, f64) {
        let f = &self.f[self.idx(x, y)];
        let rho: f64 = f.iter().sum();
        let mut ux = 0.0;
        let mut uy = 0.0;
        for i in 0..Q {
            ux += f64::from(E[i].0) * f[i];
            uy += f64::from(E[i].1) * f[i];
        }
        // Guo: half the force is added to the momentum.
        ((ux + 0.5 * self.gx) / rho, uy / rho)
    }

    /// Total mass (conserved).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.f.iter().map(|c| c.iter().sum::<f64>()).sum()
    }

    /// One collide-force-stream step.
    pub fn step(&mut self) {
        let mut post = vec![[0.0; Q]; self.nx * self.ny];
        // collide + Guo forcing.
        for y in 0..self.ny {
            for x in 0..self.nx {
                let idx = self.idx(x, y);
                let f = &self.f[idx];
                let rho: f64 = f.iter().sum();
                let mut mx = 0.0;
                let mut my = 0.0;
                for i in 0..Q {
                    mx += f64::from(E[i].0) * f[i];
                    my += f64::from(E[i].1) * f[i];
                }
                let ux = (mx + 0.5 * self.gx) / rho;
                let uy = my / rho;
                let feq = equilibrium(rho, ux, uy);
                let coef = 1.0 - 0.5 / self.tau;
                for i in 0..Q {
                    let (ex, ey) = (f64::from(E[i].0), f64::from(E[i].1));
                    let eu = ex * ux + ey * uy;
                    // Guo forcing term (force = (gx, 0)).
                    let force = coef * W[i] * (3.0 * (ex - ux) + 9.0 * eu * ex) * self.gx;
                    post[idx][i] = f[i] + (feq[i] - f[i]) / self.tau + force;
                }
            }
        }
        // pull-streaming (source = x − eᵢ) with halfway bounce-back at the
        // y-walls and periodic x. Offsets are ±1, so plain usize arithmetic.
        for y in 0..self.ny {
            for x in 0..self.nx {
                let idx = self.idx(x, y);
                for i in 0..Q {
                    // periodic source in x.
                    let sx = match E[i].0 {
                        1 => (x + self.nx - 1) % self.nx,
                        -1 => (x + 1) % self.nx,
                        _ => x,
                    };
                    // source in y, flagging a wall crossing.
                    let (in_domain, sy) = match E[i].1 {
                        1 if y == 0 => (false, 0),
                        1 => (true, y - 1),
                        -1 if y + 1 == self.ny => (false, 0),
                        -1 => (true, y + 1),
                        _ => (true, y),
                    };
                    self.f[idx][i] = if in_domain {
                        post[self.idx(sx, sy)][i]
                    } else {
                        post[idx][OPP[i]] // wall: bounce back
                    };
                }
            }
        }
    }

    /// Run `steps` steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// The mid-channel `x`-velocity profile along `y` (column `x = 0`).
    #[must_use]
    pub fn x_velocity_profile(&self) -> Vec<f64> {
        (0..self.ny).map(|y| self.velocity(0, y).0).collect()
    }
}

// -- The lattice-scaling assistant ------------------------------------------

/// A lattice-scaling plan: the derived lattice parameters + stability checks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScalingPlan {
    /// The relaxation time `τ = 3ν + ½`.
    pub tau: f64,
    /// The lattice kinematic viscosity `ν = u·L / Re`.
    pub viscosity: f64,
    /// The lattice velocity (Mach-limited).
    pub u_lattice: f64,
    /// The Mach number `u / c_s`.
    pub mach: f64,
    /// The stability margin on `τ` (`τ − ½`, must be `> 0`).
    pub tau_margin: f64,
    /// Is the plan stable (positive viscosity AND low Mach)?
    pub stable: bool,
}

/// The maximum Mach number for the low-Mach (incompressible) regime.
pub const MACH_LIMIT: f64 = 0.3;

impl ScalingPlan {
    /// The Evidence color of the plan: verified when comfortably stable,
    /// estimated when near a stability boundary.
    #[must_use]
    pub fn color(&self) -> Color {
        if self.stable && self.tau_margin > 0.05 && self.mach < 0.5 * MACH_LIMIT {
            Color::Verified {
                lo: 0.0,
                hi: self.tau_margin,
            }
        } else {
            Color::Estimated {
                estimator: "lbm-scaling".to_string(),
                dispersion: (MACH_LIMIT - self.mach).abs(),
            }
        }
    }
}

/// Plan the lattice scaling for a target Reynolds number, a characteristic
/// length in lattice units, and a chosen lattice velocity. Derives `ν`, `τ`,
/// and the Mach number, and flags stability (`τ > ½`, `Mach < MACH_LIMIT`) —
/// the bookkeeping agents should never do by hand.
///
/// # Panics
/// If `reynolds <= 0` or `char_length_lu <= 0` (a nonsensical request).
#[must_use]
pub fn plan_scaling(reynolds: f64, char_length_lu: f64, u_lattice: f64) -> ScalingPlan {
    assert!(
        reynolds > 0.0 && char_length_lu > 0.0,
        "Reynolds number and characteristic length must be positive"
    );
    let viscosity = u_lattice * char_length_lu / reynolds;
    let tau = 3.0 * viscosity + 0.5;
    let mach = u_lattice / CS2.sqrt();
    let tau_margin = tau - 0.5;
    let stable = tau_margin > 0.0 && mach < MACH_LIMIT;
    ScalingPlan {
        tau,
        viscosity,
        u_lattice,
        mach,
        tau_margin,
        stable,
    }
}

/// The analytic steady Poiseuille `x`-velocity at lattice row `y` for a channel
/// of `ny` rows under body force `gx`, with halfway bounce-back walls at
/// `y = −½` and `y = ny − ½`: `u(y) = (gx / 2ν)·(y + ½)·(ny − ½ − y)`.
#[must_use]
pub fn poiseuille_analytic(gx: f64, viscosity: f64, ny: usize, y: usize) -> f64 {
    let yf = y as f64;
    let h = ny as f64;
    // walls at y = -1/2 and y = ny - 1/2.
    (gx / (2.0 * viscosity)) * (yf + 0.5) * ((h - 0.5) - yf)
}
