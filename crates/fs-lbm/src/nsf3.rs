//! Hybrid compressible Navier-Stokes-Fourier dynamics on the D3Q27 lattice (bead pexfp.133, M3.2-W3a).
//!
//! FORMULATION (disclosed conservation class). Mass and momentum are carried by the D3Q27 populations of an
//! ISOTHERMAL lattice (equilibrium at the reference lattice temperature `theta0 = 1/3`, the lattice `c_s^2`,
//! regularised collision with the shear rate `1 / (mu(T) / (rho theta0) + 1/2)`, a DECLARED trace (bulk) rate and
//! the third-order non-equilibrium part rebuilt recursively from the deviatoric stress). The thermodynamic
//! pressure difference `p' = rho (theta - theta0)` of the ideal gas `p = rho theta` enters the momentum as the
//! Guo body force `F = -grad p'`, in conservative face form (fourth-order face values away from walls, a
//! zero-gradient wall face), with the half-force velocity `u = (sum c f + F/2) / rho`. The equilibrium therefore
//! never sees the temperature. Total energy `E = rho (c_v theta + |u|^2 / 2)` is a conservative finite-volume
//! field on the same cells with the face flux `(E + p) u - sigma . u + q`, `q = -kappa grad T`; `sigma = -(1 -
//! w/2) Pi_neq` comes from the same non-equilibrium moments that drive the momentum. Interior face values are
//! single-valued, so the sum of E over the fluid cells changes only by the wall heat flux: ENERGY IS CONSERVED
//! TO ROUND-OFF in a closed adiabatic box by flux telescoping (class: total-energy conservative). The internal
//! energy is `E` minus the kinetic energy of the force-corrected velocity, so kinetic energy removed by
//! viscosity or by a wall becomes internal energy; `theta = (gamma - 1) e` feeds back only through the force.
//! Time coupling: the trapezoid rule on the flux divergences of the two time levels with a fixed number of
//! predictor-corrector passes (`COUPLING_PASSES`).
//!
//! WHY NOT THE TEMPERATURE-CORRECTED EQUILIBRIUM. The first formulation put `theta - theta0` into the
//! equilibrium (Hermite expansion through third order). Its Bloch-linearised one-step matrix has a real
//! eigenvalue of 1.4 to 2.0 at oblique wavevectors near the Nyquist limit at air viscosity (eigenvector: the rest
//! population and the energy): at fixed total energy `p = (gamma - 1)(E - KE)` does not depend on `rho`, so a
//! density perturbation enters the equilibrium with +2.5 on the rest population and -2 on the corners, and the
//! trace relaxation amplifies it. The energy-free equilibrium used here has |lambda| < 1 on the same scan.
//!
//! LEDGER SPLIT. `pressure_work_j = -integral(p div u) dV dt` is the reversible kinetic-to-internal exchange
//! (trapezoid of the two time levels, fourth-order central divergence). `viscous_heating_j` is DEFINED as
//! `internal change - wall heat - pressure work`; total-energy conservation then makes it equal the kinetic
//! energy loss net of pressure work TERM BY TERM to round-off. That identity is bookkeeping, not evidence; the
//! physical content (the heating equals the continuum dissipation) is tested on a Taylor-Green vortex.
//!
//! WHAT THIS IS NOT: not Boussinesq (`thermal.rs` and `d3q19/coupled.rs` are); not a thermal-LBM with an energy
//! population; the thermal pressure enters as a force, so the lattice resolves the sound speed `sqrt(gamma
//! theta)` only while it stays below the lattice speed (declared Mach and temperature-ratio bounds, refused
//! beyond); no turbulence model, open boundary or moving wall. Sound-speed dispersion is second order; see the
//! battery.

use crate::lattice3::Lattice3;
use core::fmt;
use fs_material::gas::{GasSpec, GasState};

/// Reference lattice temperature `theta0 = r T0` in lattice units; equals the lattice `c_s^2`.
pub const THETA0: f64 = 1.0 / 3.0;
const Q: usize = 27;
/// Predictor-corrector passes of the energy/force coupling per step.
pub const COUPLING_PASSES: usize = 3;
/// Version of the arithmetic order, coefficients and coupling (bump on any change that moves bits).
pub const NSF3_BIT_SEMANTICS_VERSION: u32 = 2;

/// Refusals of the NSF grid. A failed step leaves the published state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum NsfError {
    /// A malformed input, with the reason.
    Input(&'static str),
    /// The gas receipt refused the state (outside its validity window or non-physical spec).
    Gas(String),
    /// A temperature outside the declared run window.
    Window {
        /// Offending cell index.
        cell: usize,
        /// Its temperature [K].
        temperature_k: f64,
        /// Window ends [K].
        low_k: f64,
        /// Window ends [K].
        high_k: f64,
    },
    /// More cells than the declared cap.
    Cap {
        /// Requested cell count.
        requested: usize,
        /// Cap.
        maximum: usize,
    },
    /// Local Mach number above the declared bound.
    Mach {
        /// Offending cell index.
        cell: usize,
        /// Observed Mach number.
        observed: f64,
        /// Declared bound.
        maximum: f64,
    },
    /// `|theta / theta0 - 1|` above the declared validity bound.
    TemperatureRatio {
        /// Offending cell index.
        cell: usize,
        /// Observed ratio.
        observed: f64,
        /// Declared bound.
        maximum: f64,
    },
    /// A non-finite or non-positive density or energy.
    NonFinite {
        /// Offending cell index.
        cell: usize,
    },
}

impl fmt::Display for NsfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for NsfError {}

/// Thermal condition of every wall face (solid cells and non-periodic domain faces).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WallThermal {
    /// Zero heat flux.
    Adiabatic,
    /// Fixed wall temperature [K], flux through the half cell.
    Isothermal(f64),
}

/// Everything the grid needs, in physical units. Nothing is defaulted.
#[derive(Clone, Debug)]
pub struct NsfSpec {
    /// Cell counts per axis.
    pub dims: [usize; 3],
    /// Periodic axes; a non-periodic face is a halfway bounce-back wall.
    pub periodic: [bool; 3],
    /// Gas law (receipt source of `gamma`, `R`, `mu(T)`, `kappa(T)`).
    pub gas: GasSpec,
    /// Reference temperature [K] and pressure [Pa] (lattice density 1, `theta0 = 1/3`).
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// Cell size [m]; the time step follows from `theta0` and the reference gas constant.
    pub dx_m: f64,
    /// Declared temperature window of the run [K]; transport is tabulated across it from the gas receipt.
    pub temperature_window_k: [f64; 2],
    /// Local Mach bound `|u| / sqrt(gamma theta)`.
    pub maximum_mach: f64,
    /// Bound of `|theta / theta0 - 1|`.
    pub maximum_temperature_ratio: f64,
    /// Trace (bulk) relaxation rate in (0, 2); the bulk viscosity is `(2/3) theta0 (1/w - 1/2)` per density.
    pub bulk_rate: f64,
    /// Cap on the number of cells (cap admitted, cap + 1 refused).
    pub maximum_cells: usize,
    /// Thermal condition of every wall face.
    pub wall: WallThermal,
    /// Threads for the cell passes; the result is bitwise independent of this number.
    pub workers: usize,
    /// Table points of the transport laws (>= 2).
    pub transport_points: usize,
}

/// Physical <-> lattice scaling and the receipt-derived transport laws.
#[derive(Clone, Debug)]
pub struct NsfLaw {
    /// Heat-capacity ratio of the reference state.
    pub gamma: f64,
    /// Specific gas constant [J/(kg K)].
    pub gas_constant: f64,
    /// Reference density [kg/m^3] (lattice density 1).
    pub rho_ref_kg_m3: f64,
    /// Cell size [m].
    pub dx_m: f64,
    /// Time step [s]: `dx / sqrt(3 R T0)`.
    pub dt_s: f64,
    /// Reference temperature [K] (lattice `theta0 = 1/3`).
    pub temperature_ref_k: f64,
    t_low: f64,
    t_high: f64,
    mu: Vec<f64>,
    kappa: Vec<f64>,
}

impl NsfLaw {
    /// Lattice temperature of a physical temperature: `theta = T / (3 T0)`.
    #[must_use]
    pub fn theta_of_kelvin(&self, t: f64) -> f64 {
        t / (3.0 * self.temperature_ref_k)
    }
    /// Physical temperature of a lattice temperature.
    #[must_use]
    pub fn kelvin_of_theta(&self, theta: f64) -> f64 {
        theta * 3.0 * self.temperature_ref_k
    }
    /// Energy unit [J] of one lattice energy unit: `rho_ref (dx/dt)^2 dx^3`.
    #[must_use]
    pub fn energy_unit_j(&self) -> f64 {
        self.rho_ref_kg_m3 * (self.dx_m / self.dt_s).powi(2) * self.dx_m.powi(3)
    }
    fn table(&self, tab: &[f64], theta: f64) -> f64 {
        let t = self.kelvin_of_theta(theta);
        let x = ((t - self.t_low) / (self.t_high - self.t_low)).clamp(0.0, 1.0) * (tab.len() - 1) as f64;
        let i = (x.floor() as usize).min(tab.len() - 2);
        let w = x - i as f64;
        (1.0 - w) * tab[i] + w * tab[i + 1]
    }
    /// Lattice dynamic viscosity at lattice temperature `theta`, from the gas receipt.
    #[must_use]
    pub fn mu_lattice(&self, theta: f64) -> f64 {
        self.table(&self.mu, theta)
    }
    /// Lattice thermal conductivity at lattice temperature `theta`, from the gas receipt.
    #[must_use]
    pub fn kappa_lattice(&self, theta: f64) -> f64 {
        self.table(&self.kappa, theta)
    }
}

/// Energy ledger in joules for one step and cumulative. `residual_rel` closes the total-energy budget.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnergyLedger {
    /// Internal energy after the step [J].
    pub internal_j: f64,
    /// Kinetic energy after the step [J].
    pub kinetic_j: f64,
    /// Viscous heating of the step [J] (defined from the budget, see the module note).
    pub viscous_heating_j: f64,
    /// Reversible kinetic-to-internal exchange `-integral(p div u)` of the step [J].
    pub pressure_work_j: f64,
    /// Heat entering through the walls during the step [J].
    pub wall_heat_j: f64,
    /// `|change of total energy - wall heat| / total energy` (round-off for a conservative step).
    pub residual_rel: f64,
}

/// Equilibrium populations (third-order Hermite about `THETA0`; `theta` corrects the second and third moments).
#[must_use]
pub fn equilibrium<L: Lattice3<27>>(rho: f64, u: [f64; 3], theta: f64) -> [f64; Q] {
    let dth = theta - THETA0;
    let u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let mut out = [0.0; Q];
    for q in 0..Q {
        let c = L::E[q];
        let cf = [f64::from(c[0]), f64::from(c[1]), f64::from(c[2])];
        let cu = cf[0] * u[0] + cf[1] * u[1] + cf[2] * u[2];
        let c2 = cf[0] * cf[0] + cf[1] * cf[1] + cf[2] * cf[2];
        let t2 = (cu * cu - THETA0 * u2 + dth * (c2 - 3.0 * THETA0)) / (2.0 * THETA0 * THETA0);
        let t3 = (cu * cu * cu - 3.0 * THETA0 * cu * u2 + 3.0 * dth * cu * (c2 - 5.0 * THETA0))
            / (6.0 * THETA0 * THETA0 * THETA0);
        out[q] = L::W[q] * rho * (1.0 + cu / THETA0 + t2 + t3);
    }
    out
}

/// Density and first moment `sum c f`.
fn raw_moments<L: Lattice3<27>>(f: &[f64; Q]) -> (f64, [f64; 3]) {
    let mut rho = 0.0;
    let mut m = [0.0; 3];
    for q in 0..Q {
        rho += f[q];
        for a in 0..3 {
            m[a] += f64::from(L::E[q][a]) * f[q];
        }
    }
    (rho, m)
}

/// Second moment of the non-equilibrium part: `[xx, yy, zz, xy, xz, yz]`.
fn neq_second<L: Lattice3<27>>(f: &[f64; Q], feq: &[f64; Q]) -> [f64; 6] {
    let mut p = [0.0; 6];
    for q in 0..Q {
        let d = f[q] - feq[q];
        let c = L::E[q];
        let (x, y, z) = (f64::from(c[0]), f64::from(c[1]), f64::from(c[2]));
        p[0] += x * x * d;
        p[1] += y * y * d;
        p[2] += z * z * d;
        p[3] += x * y * d;
        p[4] += x * z * d;
        p[5] += y * z * d;
    }
    p
}

fn half_force_velocity(rho: f64, mom: [f64; 3], force: [f64; 3]) -> [f64; 3] {
    [(mom[0] + 0.5 * force[0]) / rho, (mom[1] + 0.5 * force[1]) / rho, (mom[2] + 0.5 * force[2]) / rho]
}

/// Isothermal regularised collision with Guo forcing: shear relaxation `ws`, trace relaxation `wb`. The first-moment
/// non-equilibrium `-F/2` of Guo's scheme, relaxed at `ws`, is added back because the regularised projection drops it
/// (without it only `(3 - ws) / 2` of the force would reach the momentum).
fn collide<L: Lattice3<27>>(f: &[f64; Q], force: [f64; 3], ws: f64, wb: f64) -> [f64; Q] {
    let (rho, mom) = raw_moments::<L>(f);
    let u = half_force_velocity(rho, mom, force);
    let feq = equilibrium::<L>(rho, u, THETA0);
    let p = neq_second::<L>(f, &feq);
    let tr = p[0] + p[1] + p[2];
    let d = [p[0] - tr / 3.0, p[1] - tr / 3.0, p[2] - tr / 3.0, p[3], p[4], p[5]];
    let du = [
        d[0] * u[0] + d[3] * u[1] + d[4] * u[2],
        d[3] * u[0] + d[1] * u[1] + d[5] * u[2],
        d[4] * u[0] + d[5] * u[1] + d[2] * u[2],
    ];
    let uf = u[0] * force[0] + u[1] * force[1] + u[2] * force[2];
    let mut out = [0.0; Q];
    for q in 0..Q {
        let c = L::E[q];
        let (x, y, z) = (f64::from(c[0]), f64::from(c[1]), f64::from(c[2]));
        let c2 = x * x + y * y + z * z;
        let cdc = x * x * d[0] + y * y * d[1] + z * z * d[2] + 2.0 * (x * y * d[3] + x * z * d[4] + y * z * d[5]);
        let cu = x * u[0] + y * u[1] + z * u[2];
        let cdu = x * du[0] + y * du[1] + z * du[2];
        let cf = x * force[0] + y * force[1] + z * force[2];
        let h2 = (1.0 - ws) * cdc + (1.0 - wb) * (tr / 3.0) * c2 - THETA0 * (1.0 - wb) * tr;
        let h3 = (1.0 - ws) * (3.0 * cu * cdc - 6.0 * THETA0 * cdu);
        let source = (cf - uf) / THETA0 + cu * cf / (THETA0 * THETA0);
        let first = (1.0 - ws) * (-0.5) * cf / THETA0;
        out[q] = feq[q]
            + L::W[q]
                * (h2 / (2.0 * THETA0 * THETA0) + h3 / (6.0 * THETA0 * THETA0 * THETA0)
                    + (1.0 - 0.5 * ws) * source
                    + first);
    }
    out
}

/// Viscous stress `sigma_ab = -(1 - ws/2) dev(Pi_neq) - (1 - wb/2) (tr/3) delta`: `[xx, yy, zz, xy, xz, yz]`.
fn stress<L: Lattice3<27>>(f: &[f64; Q], force: [f64; 3], ws: f64, wb: f64) -> [f64; 6] {
    let (rho, mom) = raw_moments::<L>(f);
    let u = half_force_velocity(rho, mom, force);
    let feq = equilibrium::<L>(rho, u, THETA0);
    let p = neq_second::<L>(f, &feq);
    let tr = p[0] + p[1] + p[2];
    let a = 1.0 - ws / 2.0;
    let b = (1.0 - wb / 2.0) * tr / 3.0;
    [-a * (p[0] - tr / 3.0) - b, -a * (p[1] - tr / 3.0) - b, -a * (p[2] - tr / 3.0) - b, -a * p[3], -a * p[4], -a * p[5]]
}

/// Fill `out[i] = f(i)` over `workers` scoped threads with fixed contiguous chunks; the value of every cell is a
/// pure function of read-only inputs, so the bits do not depend on the worker count.
fn par_fill<T: Send>(out: &mut [T], workers: usize, f: &(impl Fn(usize) -> T + Sync)) {
    let n = out.len();
    if workers <= 1 || n < 2 {
        for (i, o) in out.iter_mut().enumerate() {
            *o = f(i);
        }
        return;
    }
    let chunk = n.div_ceil(workers);
    std::thread::scope(|s| {
        for (k, part) in out.chunks_mut(chunk).enumerate() {
            s.spawn(move || {
                for (j, o) in part.iter_mut().enumerate() {
                    *o = f(k * chunk + j);
                }
            });
        }
    });
}

/// The hybrid grid. See the module documentation for the formulation and its conservation class.
pub struct NsfGrid3<L: Lattice3<27>> {
    spec: NsfSpec,
    law: NsfLaw,
    n: usize,
    solid: Vec<bool>,
    /// Pre-collision populations (first moment `sum c f = rho u - F/2`).
    f: Vec<[f64; Q]>,
    e: Vec<f64>,
    theta: Vec<f64>,
    force: Vec<[f64; 3]>,
    div: Vec<f64>,
    wall_out: f64,
    steps: usize,
    cumulative: EnergyLedger,
    marker: core::marker::PhantomData<L>,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Build the grid. `initial(cell) -> (rho_lattice, velocity_lattice, temperature_K)`; solid cells are given by
    /// `solid`. Every refusal (cap, window, Mach, temperature ratio, non-finite) happens here and in `step`.
    ///
    /// # Errors
    /// [`NsfError`] for any refused input or state.
    pub fn new(
        spec: NsfSpec,
        solid: impl Fn([usize; 3]) -> bool,
        initial: impl Fn([usize; 3]) -> (f64, [f64; 3], f64),
    ) -> Result<Self, NsfError> {
        let n = spec.dims.iter().try_fold(1usize, |a, &d| a.checked_mul(d)).ok_or(NsfError::Input("dimension overflow"))?;
        if spec.dims.contains(&0) {
            return Err(NsfError::Input("dimensions must be positive"));
        }
        if n > spec.maximum_cells {
            return Err(NsfError::Cap { requested: n, maximum: spec.maximum_cells });
        }
        let ok = spec.dx_m.is_finite() && spec.dx_m > 0.0 && spec.maximum_mach.is_finite() && spec.maximum_mach > 0.0
            && spec.maximum_temperature_ratio.is_finite() && spec.maximum_temperature_ratio > 0.0
            && spec.bulk_rate > 0.0 && spec.bulk_rate < 2.0 && spec.workers >= 1 && spec.transport_points >= 2
            && spec.temperature_window_k[0] < spec.temperature_window_k[1];
        if !ok {
            return Err(NsfError::Input("dx, Mach and temperature bounds, bulk rate in (0,2), workers, table points, window"));
        }
        if let WallThermal::Isothermal(t) = spec.wall {
            if !(t.is_finite() && t >= spec.temperature_window_k[0] && t <= spec.temperature_window_k[1]) {
                return Err(NsfError::Input("isothermal wall temperature must lie in the declared window"));
            }
        }
        let state = GasState::try_new(&spec.gas, spec.temperature_k, spec.pressure_pa).map_err(|e| NsfError::Gas(e.to_string()))?;
        let r = state.specific_gas_constant;
        let dt_s = spec.dx_m / (3.0 * r * state.temperature).sqrt();
        let sc = dt_s / (state.density * spec.dx_m * spec.dx_m);
        let mut mu = Vec::with_capacity(spec.transport_points);
        let mut kappa = Vec::with_capacity(spec.transport_points);
        for i in 0..spec.transport_points {
            let t = spec.temperature_window_k[0]
                + (spec.temperature_window_k[1] - spec.temperature_window_k[0]) * i as f64 / (spec.transport_points - 1) as f64;
            let g = GasState::try_new(&spec.gas, t, spec.pressure_pa).map_err(|e| NsfError::Gas(e.to_string()))?;
            mu.push(g.dynamic_viscosity * sc);
            kappa.push(g.thermal_conductivity * sc / r);
        }
        let law = NsfLaw {
            gamma: state.gamma,
            gas_constant: r,
            rho_ref_kg_m3: state.density,
            dx_m: spec.dx_m,
            dt_s,
            temperature_ref_k: state.temperature,
            t_low: spec.temperature_window_k[0],
            t_high: spec.temperature_window_k[1],
            mu,
            kappa,
        };
        let cv = 1.0 / (law.gamma - 1.0);
        let at = |i: usize| [i % spec.dims[0], (i / spec.dims[0]) % spec.dims[1], i / (spec.dims[0] * spec.dims[1])];
        let mut rho = vec![1.0; n];
        let mut vel = vec![[0.0; 3]; n];
        let mut theta = vec![THETA0; n];
        let mut solid_v = vec![false; n];
        let mut e = vec![cv * THETA0; n];
        for i in 0..n {
            let x = at(i);
            if solid(x) {
                solid_v[i] = true;
                continue;
            }
            let (r0, u, t_k) = initial(x);
            if !(r0.is_finite() && r0 > 0.0 && u.iter().all(|v| v.is_finite()) && t_k.is_finite()) {
                return Err(NsfError::NonFinite { cell: i });
            }
            if t_k < spec.temperature_window_k[0] || t_k > spec.temperature_window_k[1] {
                return Err(NsfError::Window { cell: i, temperature_k: t_k, low_k: spec.temperature_window_k[0], high_k: spec.temperature_window_k[1] });
            }
            let th = law.theta_of_kelvin(t_k);
            check_state(&spec, &law, i, r0, u, th)?;
            rho[i] = r0;
            vel[i] = u;
            theta[i] = th;
            e[i] = r0 * (cv * th + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
        }
        let mut grid = NsfGrid3 {
            spec,
            law,
            n,
            solid: solid_v,
            f: vec![[0.0; Q]; n],
            e,
            theta,
            force: vec![[0.0; 3]; n],
            div: vec![0.0; n],
            wall_out: 0.0,
            steps: 0,
            cumulative: EnergyLedger::default(),
            marker: core::marker::PhantomData,
        };
        grid.force = grid.force_field(&rho, &grid.theta);
        // first moment sum c f = rho u - F/2 so that the half-force velocity equals the requested one
        for i in 0..n {
            grid.f[i] = if grid.solid[i] {
                equilibrium::<L>(1.0, [0.0; 3], THETA0)
            } else {
                let f = grid.force[i];
                let u = [vel[i][0] - 0.5 * f[0] / rho[i], vel[i][1] - 0.5 * f[1] / rho[i], vel[i][2] - 0.5 * f[2] / rho[i]];
                equilibrium::<L>(rho[i], u, THETA0)
            };
        }
        let (rho_g, mom) = grid.macro_fields(&grid.f);
        let u = grid.velocity_field(&rho_g, &mom, &grid.force);
        let sig = grid.stress_field(&grid.f, &grid.force, &rho_g, &grid.theta);
        let (div, wall) = grid.flux_divergence(&grid.e, &rho_g, &u, &grid.theta, &sig);
        grid.div = div;
        grid.wall_out = wall;
        Ok(grid)
    }

    fn at(&self, i: usize) -> [usize; 3] {
        let d = self.spec.dims;
        [i % d[0], (i / d[0]) % d[1], i / (d[0] * d[1])]
    }
    fn index(&self, x: [usize; 3]) -> usize {
        let d = self.spec.dims;
        x[0] + d[0] * (x[1] + d[1] * x[2])
    }
    /// Neighbour of cell `x` shifted by `s` along `axis`; `None` beyond a non-periodic face.
    fn shifted(&self, x: [usize; 3], axis: usize, s: i64) -> Option<usize> {
        let d = self.spec.dims[axis] as i64;
        let v = x[axis] as i64 + s;
        let w = if (0..d).contains(&v) {
            v
        } else if self.spec.periodic[axis] {
            v.rem_euclid(d)
        } else {
            return None;
        };
        let mut y = x;
        y[axis] = w as usize;
        Some(self.index(y))
    }
    fn fluid_at(&self, x: [usize; 3], axis: usize, s: i64) -> Option<usize> {
        self.shifted(x, axis, s).filter(|&j| !self.solid[j])
    }

    fn macro_fields(&self, f: &[[f64; Q]]) -> (Vec<f64>, Vec<[f64; 3]>) {
        let mut rm = vec![(1.0, [0.0; 3]); self.n];
        let solid = &self.solid;
        par_fill(&mut rm, self.spec.workers, &|i| if solid[i] { (1.0, [0.0; 3]) } else { raw_moments::<L>(&f[i]) });
        (rm.iter().map(|v| v.0).collect(), rm.iter().map(|v| v.1).collect())
    }

    fn velocity_field(&self, rho: &[f64], mom: &[[f64; 3]], force: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut out = vec![[0.0; 3]; self.n];
        let solid = &self.solid;
        par_fill(&mut out, self.spec.workers, &|i| if solid[i] { [0.0; 3] } else { half_force_velocity(rho[i], mom[i], force[i]) });
        out
    }

    fn rates(&self, rho: f64, theta: f64) -> (f64, f64) {
        let mu = self.law.mu_lattice(theta);
        (1.0 / (mu / (rho * THETA0) + 0.5), self.spec.bulk_rate)
    }

    /// `F = -grad p'`, `p' = rho (theta - theta0)`, in conservative face form: fourth-order face values where the four
    /// stencil cells are fluid, the two-point value beside a solid cell, the cell's own value at a wall face.
    fn force_field(&self, rho: &[f64], theta: &[f64]) -> Vec<[f64; 3]> {
        let pp: Vec<f64> = (0..self.n).map(|i| if self.solid[i] { 0.0 } else { rho[i] * (theta[i] - THETA0) }).collect();
        let solid = &self.solid;
        let face = |i: usize, a: usize, s: i64| -> f64 {
            // face value on the +a side (s = 1) or -a side (s = -1) of cell i
            let x = self.at(i);
            let Some(j) = self.fluid_at(x, a, s) else { return pp[i] };
            let xj = self.at(j);
            let far_i = self.fluid_at(x, a, -s);
            let far_j = self.fluid_at(xj, a, s);
            match (far_i, far_j) {
                (Some(m), Some(p)) => (7.0 / 12.0) * (pp[i] + pp[j]) - (1.0 / 12.0) * (pp[m] + pp[p]),
                _ => 0.5 * (pp[i] + pp[j]),
            }
        };
        let mut out = vec![[0.0; 3]; self.n];
        par_fill(&mut out, self.spec.workers, &|i| {
            if solid[i] {
                return [0.0; 3];
            }
            [-(face(i, 0, 1) - face(i, 0, -1)), -(face(i, 1, 1) - face(i, 1, -1)), -(face(i, 2, 1) - face(i, 2, -1))]
        });
        out
    }

    fn stress_field(&self, f: &[[f64; Q]], force: &[[f64; 3]], rho: &[f64], theta: &[f64]) -> Vec<[f64; 6]> {
        let mut out = vec![[0.0; 6]; self.n];
        let solid = &self.solid;
        par_fill(&mut out, self.spec.workers, &|i| {
            if solid[i] {
                [0.0; 6]
            } else {
                let (ws, wb) = self.rates(rho[i], theta[i]);
                stress::<L>(&f[i], force[i], ws, wb)
            }
        });
        out
    }

    /// Divergence of the total-energy flux at ONE time level and the total heat leaving through the walls.
    fn flux_divergence(&self, et: &[f64], rho: &[f64], u: &[[f64; 3]], theta: &[f64], sig: &[[f64; 6]]) -> (Vec<f64>, f64) {
        let n = self.n;
        let mut g = vec![[0.0; 3]; n];
        let solid = &self.solid;
        par_fill(&mut g, self.spec.workers, &|i| {
            if solid[i] {
                return [0.0; 3];
            }
            let s = sig[i];
            let p = rho[i] * theta[i];
            let su = [
                s[0] * u[i][0] + s[3] * u[i][1] + s[4] * u[i][2],
                s[3] * u[i][0] + s[1] * u[i][1] + s[5] * u[i][2],
                s[4] * u[i][0] + s[5] * u[i][1] + s[2] * u[i][2],
            ];
            [(et[i] + p) * u[i][0] - su[0], (et[i] + p) * u[i][1] - su[1], (et[i] + p) * u[i][2] - su[2]]
        });
        let wall_theta = match self.spec.wall {
            WallThermal::Adiabatic => None,
            WallThermal::Isothermal(t) => Some(self.law.theta_of_kelvin(t)),
        };
        let face = |i: usize, a: usize| -> Option<f64> {
            let x = self.at(i);
            let j = self.fluid_at(x, a, 1)?;
            let kf = self.law.kappa_lattice(0.5 * (theta[i] + theta[j]));
            let wide = match (self.fluid_at(x, a, -1), self.fluid_at(self.at(j), a, 1)) {
                (Some(m), Some(p)) => Some((m, p)),
                _ => None,
            };
            Some(match wide {
                Some((m, p)) => {
                    let gf = (7.0 / 12.0) * (g[i][a] + g[j][a]) - (1.0 / 12.0) * (g[m][a] + g[p][a]);
                    let dt4 = (27.0 * (theta[j] - theta[i]) - (theta[p] - theta[m])) / 24.0;
                    gf - kf * dt4
                }
                None => 0.5 * (g[i][a] + g[j][a]) - kf * (theta[j] - theta[i]),
            })
        };
        let mut plus = vec![[0.0; 3]; n];
        par_fill(&mut plus, self.spec.workers, &|i| {
            if solid[i] {
                return [0.0; 3];
            }
            [face(i, 0).unwrap_or(0.0), face(i, 1).unwrap_or(0.0), face(i, 2).unwrap_or(0.0)]
        });
        // wall face: flux in the +axis direction (heat only; u = 0 at a bounce-back wall)
        let wall_flux = |i: usize, s: i64| -> f64 {
            match wall_theta {
                None => 0.0,
                Some(tw) => {
                    let kf = self.law.kappa_lattice(0.5 * (theta[i] + tw));
                    if s > 0 { -kf * (tw - theta[i]) / 0.5 } else { -kf * (theta[i] - tw) / 0.5 }
                }
            }
        };
        let mut div = vec![0.0; n];
        par_fill(&mut div, self.spec.workers, &|i| {
            if solid[i] {
                return 0.0;
            }
            let x = self.at(i);
            let mut d = 0.0;
            for a in 0..3 {
                d += match self.fluid_at(x, a, 1) {
                    Some(_) => plus[i][a],
                    None => wall_flux(i, 1),
                };
                d -= match self.fluid_at(x, a, -1) {
                    Some(m) => plus[m][a],
                    None => wall_flux(i, -1),
                };
            }
            d
        });
        let mut out = 0.0;
        for i in 0..n {
            if solid[i] {
                continue;
            }
            let x = self.at(i);
            for a in 0..3 {
                if self.fluid_at(x, a, 1).is_none() {
                    out += wall_flux(i, 1);
                }
                if self.fluid_at(x, a, -1).is_none() {
                    out -= wall_flux(i, -1);
                }
            }
        }
        (div, out)
    }

    fn stream(&self, post: &[[f64; Q]]) -> Vec<[f64; Q]> {
        let mut next = vec![[0.0; Q]; self.n];
        let solid = &self.solid;
        par_fill(&mut next, self.spec.workers, &|i| {
            if solid[i] {
                return post[i];
            }
            let x = self.at(i);
            let mut out = [0.0; Q];
            for q in 0..Q {
                let c = L::E[q];
                let mut y = Some(x);
                for a in 0..3 {
                    if c[a] != 0 {
                        y = y.and_then(|yy| self.shifted(yy, a, -i64::from(c[a])).map(|j| self.at(j)));
                    }
                }
                out[q] = match y.map(|yy| self.index(yy)) {
                    Some(j) if !self.solid[j] => post[j][q],
                    _ => post[i][L::OPP[q]],
                };
            }
            out
        });
        next
    }

    /// Advance one lattice step and return the ledger of that step in joules.
    ///
    /// # Errors
    /// [`NsfError`] when the new state leaves the declared Mach, temperature-ratio, window or finite domain; the
    /// published state is unchanged.
    pub fn step(&mut self) -> Result<EnergyLedger, NsfError> {
        let n = self.n;
        let cv = 1.0 / (self.law.gamma - 1.0);
        let (rho0, mom0) = self.macro_fields(&self.f);
        let u0 = self.velocity_field(&rho0, &mom0, &self.force);
        let mut post = vec![[0.0; Q]; n];
        {
            let solid = &self.solid;
            let f = &self.f;
            let force = &self.force;
            let theta = &self.theta;
            let rho = &rho0;
            par_fill(&mut post, self.spec.workers, &|i| {
                if solid[i] {
                    f[i]
                } else {
                    let (ws, wb) = self.rates(rho[i], theta[i]);
                    collide::<L>(&f[i], force[i], ws, wb)
                }
            });
        }
        let next = self.stream(&post);
        let (rho1, mom1) = self.macro_fields(&next);
        let mut th1 = self.theta.clone();
        let mut result = None;
        for _ in 0..COUPLING_PASSES {
            let f1 = self.force_field(&rho1, &th1);
            let u1 = self.velocity_field(&rho1, &mom1, &f1);
            let et_trial: Vec<f64> = (0..n)
                .map(|i| rho1[i] * (cv * th1[i] + 0.5 * (u1[i][0] * u1[i][0] + u1[i][1] * u1[i][1] + u1[i][2] * u1[i][2])))
                .collect();
            let sig = self.stress_field(&next, &f1, &rho1, &th1);
            let (div1, wall1) = self.flux_divergence(&et_trial, &rho1, &u1, &th1, &sig);
            let mut e1 = self.e.clone();
            let mut th_new = th1.clone();
            for i in 0..n {
                if self.solid[i] {
                    continue;
                }
                e1[i] = self.e[i] - 0.5 * (div1[i] + self.div[i]);
                let ke = 0.5 * rho1[i] * (u1[i][0] * u1[i][0] + u1[i][1] * u1[i][1] + u1[i][2] * u1[i][2]);
                th_new[i] = (e1[i] - ke) / (rho1[i] * cv);
            }
            th1 = th_new;
            result = Some((f1, u1, e1, div1, wall1));
        }
        let (force1, u1, e1, div1, wall1) = result.expect("at least one coupling pass");
        for i in 0..n {
            if self.solid[i] {
                continue;
            }
            if !(rho1[i].is_finite() && rho1[i] > 0.0 && e1[i].is_finite() && th1[i].is_finite() && th1[i] > 0.0) {
                return Err(NsfError::NonFinite { cell: i });
            }
            check_state(&self.spec, &self.law, i, rho1[i], u1[i], th1[i])?;
            let t_k = self.law.kelvin_of_theta(th1[i]);
            if t_k < self.spec.temperature_window_k[0] || t_k > self.spec.temperature_window_k[1] {
                return Err(NsfError::Window { cell: i, temperature_k: t_k, low_k: self.spec.temperature_window_k[0], high_k: self.spec.temperature_window_k[1] });
            }
        }
        let unit = self.law.energy_unit_j();
        let (mut ke0, mut ke1, mut e_old, mut e_new, mut pw) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            if self.solid[i] {
                continue;
            }
            ke0 += 0.5 * rho0[i] * (u0[i][0] * u0[i][0] + u0[i][1] * u0[i][1] + u0[i][2] * u0[i][2]);
            ke1 += 0.5 * rho1[i] * (u1[i][0] * u1[i][0] + u1[i][1] * u1[i][1] + u1[i][2] * u1[i][2]);
            e_old += self.e[i];
            e_new += e1[i];
        }
        let p0 = self.pressure_work_density(&rho0, &u0, &self.theta);
        let p1 = self.pressure_work_density(&rho1, &u1, &th1);
        for i in 0..n {
            if !self.solid[i] {
                pw += 0.5 * (p0[i] + p1[i]);
            }
        }
        let wall_in = -0.5 * (wall1 + self.wall_out);
        let d_total = e_new - e_old;
        let d_kin = ke1 - ke0;
        let d_int = d_total - d_kin;
        let pressure_work = -pw;
        let viscous = d_int - wall_in - pressure_work;
        let ledger = EnergyLedger {
            internal_j: unit * (e_new - ke1),
            kinetic_j: unit * ke1,
            viscous_heating_j: unit * viscous,
            pressure_work_j: unit * pressure_work,
            wall_heat_j: unit * wall_in,
            residual_rel: (d_total - wall_in).abs() / e_new.abs().max(f64::MIN_POSITIVE),
        };
        self.f = next;
        self.e = e1;
        self.theta = th1;
        self.force = force1;
        self.div = div1;
        self.wall_out = wall1;
        self.steps += 1;
        self.cumulative.viscous_heating_j += ledger.viscous_heating_j;
        self.cumulative.pressure_work_j += ledger.pressure_work_j;
        self.cumulative.wall_heat_j += ledger.wall_heat_j;
        self.cumulative.internal_j = ledger.internal_j;
        self.cumulative.kinetic_j = ledger.kinetic_j;
        Ok(ledger)
    }

    /// `p div u` per fluid cell with the fourth-order central divergence (second order beside walls).
    fn pressure_work_density(&self, rho: &[f64], u: &[[f64; 3]], theta: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.n];
        let solid = &self.solid;
        par_fill(&mut out, self.spec.workers, &|i| {
            if solid[i] {
                return 0.0;
            }
            let x = self.at(i);
            let mut div = 0.0;
            for a in 0..3 {
                let (p1, m1) = (self.fluid_at(x, a, 1), self.fluid_at(x, a, -1));
                div += match (m1, p1) {
                    (Some(m), Some(p)) => match (self.fluid_at(self.at(m), a, -1), self.fluid_at(self.at(p), a, 1)) {
                        (Some(mm), Some(pp)) => (8.0 * (u[p][a] - u[m][a]) - (u[pp][a] - u[mm][a])) / 12.0,
                        _ => 0.5 * (u[p][a] - u[m][a]),
                    },
                    (None, Some(p)) => 0.5 * (u[p][a] + u[i][a]),
                    (Some(m), None) => -0.5 * (u[m][a] + u[i][a]),
                    (None, None) => 0.0,
                };
            }
            rho[i] * theta[i] * div
        });
        out
    }

    /// Cumulative ledger since construction (internal and kinetic are the latest values).
    #[must_use]
    pub fn cumulative_ledger(&self) -> EnergyLedger {
        self.cumulative
    }
    /// Kinetic energy `sum rho |u|^2 / 2` in joules over the fluid cells in fixed cell order.
    #[must_use]
    pub fn kinetic_energy_j(&self) -> f64 {
        let (rho, mom) = self.macro_fields(&self.f);
        let u = self.velocity_field(&rho, &mom, &self.force);
        let mut s = 0.0;
        for i in 0..self.n {
            if !self.solid[i] {
                s += 0.5 * rho[i] * (u[i][0] * u[i][0] + u[i][1] * u[i][1] + u[i][2] * u[i][2]);
            }
        }
        s * self.law.energy_unit_j()
    }
    /// Total energy in joules summed over the fluid cells in fixed cell order.
    #[must_use]
    pub fn total_energy_j(&self) -> f64 {
        let mut s = 0.0;
        for i in 0..self.n {
            if !self.solid[i] {
                s += self.e[i];
            }
        }
        s * self.law.energy_unit_j()
    }
    /// Density, lattice velocity (half-force corrected) and physical temperature of a cell.
    #[must_use]
    pub fn cell(&self, x: [usize; 3]) -> (f64, [f64; 3], f64) {
        let i = self.index(x);
        let (rho, mom) = raw_moments::<L>(&self.f[i]);
        (rho, half_force_velocity(rho, mom, self.force[i]), self.law.kelvin_of_theta(self.theta[i]))
    }
    /// The scaling and transport laws in force.
    #[must_use]
    pub fn law(&self) -> &NsfLaw {
        &self.law
    }
    /// Steps completed.
    #[must_use]
    pub fn steps(&self) -> usize {
        self.steps
    }
    /// FNV-1a hash of every population and energy bit (bitwise invariance checks).
    #[must_use]
    pub fn state_hash(&self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |v: f64| {
            for b in v.to_bits().to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0100_0000_01b3);
            }
        };
        for i in 0..self.n {
            for q in 0..Q {
                feed(self.f[i][q]);
            }
            feed(self.e[i]);
        }
        h
    }
}

fn check_state(spec: &NsfSpec, law: &NsfLaw, cell: usize, rho: f64, u: [f64; 3], theta: f64) -> Result<(), NsfError> {
    if !(theta.is_finite() && theta > 0.0 && rho.is_finite() && rho > 0.0) {
        return Err(NsfError::NonFinite { cell });
    }
    let ratio = theta / THETA0 - 1.0;
    if ratio.abs() > spec.maximum_temperature_ratio {
        return Err(NsfError::TemperatureRatio { cell, observed: ratio.abs(), maximum: spec.maximum_temperature_ratio });
    }
    let mach = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt() / (law.gamma * theta).sqrt();
    if !mach.is_finite() || mach > spec.maximum_mach {
        return Err(NsfError::Mach { cell, observed: mach, maximum: spec.maximum_mach });
    }
    Ok(())
}
