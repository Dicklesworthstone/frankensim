//! Hybrid compressible Navier-Stokes-Fourier dynamics on the D3Q27 lattice (bead pexfp.133, M3.2-W3a).
//!
//! FORMULATION (disclosed conservation class). Mass and momentum are carried by the D3Q27 populations of an
//! ISOTHERMAL lattice (equilibrium at the reference lattice temperature `theta0 = 1/3`, the lattice `c_s^2`,
//! regularised collision with the shear rate `1 / (mu(T) / (rho theta0) + 1/2)`, the trace (bulk) rate derived from a
//! SOURCED bulk-to-shear viscosity ratio and the third-order non-equilibrium part rebuilt recursively from the deviatoric stress). The thermodynamic
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
//! MATERIAL-OWNED GAS. NsfSpec also accepts ThermochemicalNsfGas. The same
//! collision, pressure force, conservative face flux and fixed coupling passes
//! then use source-owned NASA9 u(T), cv(T), gamma(T), viscosity and conductivity
//! at each actual density-derived pressure. Temperature is recovered by the
//! original bounded material energy solve, not by a constant-cv division.
//! The conserved internal energy uses one fixed sensible offset, disclosed by
//! ThermochemicalNsfGas::energy_reference_j_kg, so an arbitrary formation-energy
//! datum cannot drive this hybrid mass/energy discretization. Composition is
//! frozen; no reaction, species diffusion or phase conversion is implied.
//! The original GasSpec arithmetic and qualification remain separate from
//! physical validation of this source-driven extension.
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
//! `internal change - wall heat - pressure work - source energy` for heat-only volume sources.
//! With no volume source, total-energy conservation makes it equal the kinetic energy loss net of pressure
//! work TERM BY TERM to round-off. Manufactured mass/momentum sources also do mechanical work; without a
//! separately qualified source-work split, that signed remainder is reported as
//! `unresolved_internal_exchange_j`, never as viscous heating. The exact total source supply is still
//! independently audited. These identities are bookkeeping, not evidence; the physical content of
//! viscous heating is tested on a Taylor-Green vortex.
//!
//! WHAT THIS IS NOT: not Boussinesq (`thermal.rs` and `d3q19/coupled.rs` are); not a thermal-LBM with an energy
//! population; the thermal pressure enters as a force, so the lattice resolves the sound speed `sqrt(gamma
//! theta)` only while it stays below the lattice speed (declared Mach and temperature-ratio bounds, refused
//! beyond). Explicit exploratory Smagorinsky/WALE transport is available through
//! [`NsfGrid3::with_sgs`], including conservative turbulent heat transport; it is
//! never enabled by default. The optional open-axis boundary is a separately
//! bounded local linear characteristic rung, not a qualified nonreflecting jet interface. Sound-speed dispersion is second order; see the
//! battery.
//! Prescribed reservoir ghost states are also available as a separate, mutually
//! exclusive boundary mode; they do not preserve outgoing acoustic characteristics.

use crate::lattice3::Lattice3;
use core::fmt;
use fs_material::gas::GasSpec;
use material::NsfCellMaterial;

mod material;
pub use material::{NsfGasMedium, ThermochemicalNsfGas};
mod sgs;
pub use sgs::{NsfSgsConfig, NsfSgsTransport};
mod open;
pub use open::NsfOpenState;
mod open_boundary;
use open_boundary::{BoundaryFlux, IncomingPopulations, OpenFaces};
pub use open_boundary::{
    NsfBoundaryExchange, NsfBoundaryTrace, NsfCharacteristicAxis, NsfCharacteristicBoundary,
};
mod characteristic_trial;
/// Unpublished original NSF interval for reciprocal characteristic coupling.
pub use characteristic_trial::NsfCharacteristicTrial;

/// Exact linearized Navier-Stokes-Fourier modes (test support of the qualification, shipped for reuse).
pub mod dispersion;

/// Reference lattice temperature `theta0 = r T0` in lattice units; equals the lattice `c_s^2`.
pub const THETA0: f64 = 1.0 / 3.0;
const Q: usize = 27;
/// Predictor-corrector passes of the energy/force coupling per step.
pub const COUPLING_PASSES: usize = 3;
/// Largest lattice viscosity `mu / rho0` of the transport tables for which the Bloch scan of this scheme was run and found
/// stable at the sourced bulk ratio (|lambda| <= 1.0000 on eight oblique modes for 1e-5 <= nu <= 0.03; 1.146 at nu = 0.1). A
/// larger one is refused, not adjusted.
pub const MAX_VERIFIED_LATTICE_VISCOSITY: f64 = 0.03;
/// Largest Mach number `|u| / sqrt(gamma theta)` the formulation is qualified for: the Galilean-invariance map of
/// `tests/nsf3_qualification.rs` (advected entropy and shear waves at Mach 0.05, 0.1 and 0.2) keeps every Mach-dependent defect
/// of a decay rate below 5e-3 of the quiescent rate up to this value, and nothing was measured above it. A larger declared
/// `maximum_mach` is refused ([`NsfError::MachBound`]), not extrapolated.
pub const MAX_QUALIFIED_MACH: f64 = 0.2;
/// Version of the arithmetic order, coefficients and coupling (bump on any change that moves bits). 5: initialization
/// and reservoir updates preserve the full Guo half-source chart, including its physical stress.
pub const NSF3_BIT_SEMANTICS_VERSION: u32 = 5;

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
    /// A declared Mach bound above [`MAX_QUALIFIED_MACH`].
    MachBound {
        /// Requested bound.
        requested: f64,
        /// `MAX_QUALIFIED_MACH`.
        qualified: f64,
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
    /// Lattice viscosity above the verified stable envelope for this dx (refine fewer cells per metre, i.e. enlarge dx).
    Viscosity {
        /// Largest reference-table or actual local total kinematic viscosity.
        observed: f64,
        /// `MAX_VERIFIED_LATTICE_VISCOSITY`.
        maximum: f64,
    },
    /// Explicit SGS controls or the local modeled/molecular viscosity ratio refused.
    Sgs(crate::sgs3::SgsError),
    /// A non-finite state or non-positive density/absolute temperature.
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
pub struct NsfSpec<G = GasSpec, S = &'static str> {
    /// Cell counts per axis.
    pub dims: [usize; 3],
    /// Periodic axes; a non-periodic face initially has a halfway bounce-back
    /// wall. `with_characteristic_axis` can open one complete axis cold.
    pub periodic: [bool; 3],
    /// Gas law: legacy GasSpec or source-owned ThermochemicalNsfGas.
    pub gas: G,
    /// Reference temperature [K] and pressure [Pa] (lattice density 1, `theta0 = 1/3`).
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// Cell size [m]; the time step follows from `theta0` and the reference gas constant.
    pub dx_m: f64,
    /// Declared temperature window [K]; source caloric inversion uses this exact bracket.
    pub temperature_window_k: [f64; 2],
    /// Local Mach bound `|u| / sqrt(gamma theta)`.
    pub maximum_mach: f64,
    /// Bound of `|theta / theta0 - 1|`.
    pub maximum_temperature_ratio: f64,
    /// Bulk-to-shear viscosity ratio `mu_b / mu` of the gas: a SOURCED input (dry air near 293 K: 0.73, Cramer 2012 for N2,
    /// range 0.67 to 0.84; see the CONTRACT). The trace rate follows from `(2/3) theta0 (1/w_b - 1/2) = ratio * mu / rho`.
    pub bulk_to_shear_ratio: f64,
    /// Citation of the ratio; must be non-empty so that no ratio is anonymous.
    /// Construction retains an owned string, accepting both string literals
    /// and source text loaded from a physical description without leaking it.
    pub bulk_ratio_source: S,
    /// Cap on the number of cells (cap admitted, cap + 1 refused).
    pub maximum_cells: usize,
    /// Thermal condition of every wall face.
    pub wall: WallThermal,
    /// Threads for the cell passes; the result is bitwise independent of this number.
    pub workers: usize,
    /// Table points of the transport laws (>= 2).
    pub transport_points: usize,
}

impl<G, S: Into<String>> NsfSpec<G, S> {
    fn map_gas<H>(self, map: impl FnOnce(G) -> H) -> NsfSpec<H, String> {
        NsfSpec {
            dims: self.dims,
            periodic: self.periodic,
            gas: map(self.gas),
            temperature_k: self.temperature_k,
            pressure_pa: self.pressure_pa,
            dx_m: self.dx_m,
            temperature_window_k: self.temperature_window_k,
            maximum_mach: self.maximum_mach,
            maximum_temperature_ratio: self.maximum_temperature_ratio,
            bulk_to_shear_ratio: self.bulk_to_shear_ratio,
            bulk_ratio_source: self.bulk_ratio_source.into(),
            maximum_cells: self.maximum_cells,
            wall: self.wall,
            workers: self.workers,
            transport_points: self.transport_points,
        }
    }
}

/// Physical <-> lattice scaling and the receipt-derived transport laws.
#[derive(Clone, Debug)]
pub struct NsfLaw {
    /// Heat-capacity ratio of the reference state; source-driven evolution uses local gamma(T).
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
    /// The declared bulk-to-shear viscosity ratio in force.
    pub bulk_to_shear_ratio: f64,
    t_low: f64,
    t_high: f64,
    mu: Vec<f64>,
    kappa: Vec<f64>,
    /// Exact source reference coefficients for uniform-state linearization.
    source_reference_transport: Option<[f64; 2]>,
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
        let x = ((t - self.t_low) / (self.t_high - self.t_low)).clamp(0.0, 1.0)
            * (tab.len() - 1) as f64;
        let i = (x.floor() as usize).min(tab.len() - 2);
        let w = x - i as f64;
        (1.0 - w) * tab[i] + w * tab[i + 1]
    }
    /// Reference-pressure lattice dynamic viscosity table at lattice temperature.
    /// Source-driven evolution uses direct local material queries instead of this diagnostic table.
    #[must_use]
    pub fn mu_lattice(&self, theta: f64) -> f64 {
        self.table(&self.mu, theta)
    }
    /// Reference-pressure lattice conductivity table at lattice temperature.
    /// Source-driven evolution uses direct local material queries instead of this diagnostic table.
    #[must_use]
    pub fn kappa_lattice(&self, theta: f64) -> f64 {
        self.table(&self.kappa, theta)
    }
}

/// Energy ledger in joules for one step and cumulative. `residual_rel` closes the total-energy budget.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnergyLedger {
    /// Internal energy after the step [J], in the medium's fixed sensible-energy chart.
    pub internal_j: f64,
    /// Kinetic energy after the step [J].
    pub kinetic_j: f64,
    /// Viscous heating of the step [J] (defined from the budget, see the module note).
    /// With prescribed reservoir faces this remainder also includes kinetic-energy
    /// transport and is not an independent viscous-heating measurement. Zero for
    /// a step with manufactured mass/momentum sources; see the unresolved account.
    /// With SGS this remainder also includes modeled and numerical dissipation;
    /// it is not an independent measurement of molecular viscous heating.
    pub viscous_heating_j: f64,
    /// Signed internal-energy remainder for manufactured mass/momentum sources [J].
    /// Includes physical/numerical dissipation and the unseparated mechanical
    /// action of those sources after their total energy supply is removed.
    /// It can be negative. It is neither heat nor a total-energy balancing credit.
    /// Zero for source-free and heat-only steps.
    pub unresolved_internal_exchange_j: f64,
    /// Reversible kinetic-to-internal exchange `-integral(p div u)` of the step [J].
    pub pressure_work_j: f64,
    /// Heat entering through the walls during the step [J].
    pub wall_heat_j: f64,
    /// Net total energy entering the open faces [J], including enthalpy,
    /// conduction and viscous work. Zero without characteristic boundaries.
    pub open_boundary_energy_j: f64,
    /// Total energy entering prescribed reservoir ghost faces [J].
    /// Separate from characteristic-face energy; the two modes cannot coexist.
    pub open_energy_j: f64,
    /// Total energy supplied by the volume source during this step [J].
    /// The exact trapezoid used by the energy update, in the medium's fixed
    /// sensible chart. Heat-only sources supply heat; arbitrary manufactured
    /// mass/momentum sources need not. Negative values remove energy.
    pub source_energy_j: f64,
    /// `|change of total energy - wall heat - open energy - source energy| / total energy`.
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
    [
        (mom[0] + 0.5 * force[0]) / rho,
        (mom[1] + 0.5 * force[1]) / rho,
        (mom[2] + 0.5 * force[2]) / rho,
    ]
}

/// Add half the unrelaxed Guo source for `difference = F_old - F_new`.
/// The physical populations are `f + S(u,F)/2`: preserving only their first
/// moment would add `(u delta_F + delta_F u)/2` to the physical second moment
/// and turn a force-chart change into viscous stress and energy transport.
/// Guo, Zheng and Shi (2002), Phys. Rev. E 65, 046308,
/// https://doi.org/10.1103/PhysRevE.65.046308.
fn rebase_force<L: Lattice3<27>>(f: &mut [f64; Q], u: [f64; 3], difference: [f64; 3]) {
    if difference == [0.0; 3] {
        return;
    }
    let uf = (0..3).map(|a| u[a] * difference[a]).sum::<f64>();
    for q in 0..Q {
        let c = L::E[q];
        let cf = (0..3).map(|a| f64::from(c[a]) * difference[a]).sum::<f64>();
        let cu = (0..3).map(|a| f64::from(c[a]) * u[a]).sum::<f64>();
        f[q] += 0.5 * L::W[q] * ((cf - uf) / THETA0 + cu * cf / (THETA0 * THETA0));
    }
}

/// Isothermal regularised collision with Guo forcing: shear relaxation `ws`, trace relaxation `wb`. The first-moment
/// non-equilibrium `-F/2` of Guo's scheme, relaxed at `ws`, is added back because the regularised projection drops it
/// (without it only `(3 - ws) / 2` of the force would reach the momentum).
/// Each second-moment source mode carries its OWN `(1 - rate/2)` factor:
/// Premnath, Pattison and Banerjee (2009), Eq. (1), https://arxiv.org/abs/0901.2786.
fn collide<L: Lattice3<27>>(f: &[f64; Q], force: [f64; 3], ws: f64, wb: f64) -> [f64; Q] {
    let (rho, mom) = raw_moments::<L>(f);
    let u = half_force_velocity(rho, mom, force);
    let feq = equilibrium::<L>(rho, u, THETA0);
    let p = neq_second::<L>(f, &feq);
    let tr = p[0] + p[1] + p[2];
    let d = [
        p[0] - tr / 3.0,
        p[1] - tr / 3.0,
        p[2] - tr / 3.0,
        p[3],
        p[4],
        p[5],
    ];
    let du = [
        d[0] * u[0] + d[3] * u[1] + d[4] * u[2],
        d[3] * u[0] + d[1] * u[1] + d[5] * u[2],
        d[4] * u[0] + d[5] * u[1] + d[2] * u[2],
    ];
    let uf = u[0] * force[0] + u[1] * force[1] + u[2] * force[2];
    // The raw Guo source has second moment uF + Fu and trace 2 u.F.
    // Its isotropic Hermite projection has zero mass/first moment and
    // second moment (2 u.F / 3) I. Replace only its shear-rate prefactor
    // by the bulk-rate prefactor; the deviatoric source stays unchanged.
    let trace_force = 0.5 * (ws - wb) * uf / (3.0 * THETA0 * THETA0);
    let mut out = [0.0; Q];
    for q in 0..Q {
        let c = L::E[q];
        let (x, y, z) = (f64::from(c[0]), f64::from(c[1]), f64::from(c[2]));
        let c2 = x * x + y * y + z * z;
        let cdc = x * x * d[0]
            + y * y * d[1]
            + z * z * d[2]
            + 2.0 * (x * y * d[3] + x * z * d[4] + y * z * d[5]);
        let cu = x * u[0] + y * u[1] + z * u[2];
        let cdu = x * du[0] + y * du[1] + z * du[2];
        let cf = x * force[0] + y * force[1] + z * force[2];
        let h2 = (1.0 - ws) * cdc + (1.0 - wb) * (tr / 3.0) * c2 - THETA0 * (1.0 - wb) * tr;
        let h3 = (1.0 - ws) * (3.0 * cu * cdc - 6.0 * THETA0 * cdu);
        let source = (cf - uf) / THETA0 + cu * cf / (THETA0 * THETA0);
        let first = (1.0 - ws) * (-0.5) * cf / THETA0;
        out[q] = feq[q]
            + L::W[q]
                * (h2 / (2.0 * THETA0 * THETA0)
                    + h3 / (6.0 * THETA0 * THETA0 * THETA0)
                    + (1.0 - 0.5 * ws) * source
                    + trace_force * (c2 - 3.0 * THETA0)
                    + first);
    }
    out
}

/// Viscous stress `sigma_ab = -(1 - ws/2) dev(Pi_neq) - (1 - wb/2) (tr/3) delta`: `[xx, yy, zz, xy, xz, yz]`, with `Pi_neq` the
/// second moment of the pre-collision populations about `feq` PLUS `(F u + u F) / 2`: Guo's scheme stores the populations
/// half a force source behind the physical ones (`f = f_phys - S/2`, and `S` has second moment `F u + u F`), so without the
/// shift a thermal-pressure force `F` acting on a flow `u` would appear as a viscous stress. Measured on an entropy wave advected
/// at Mach 0.2: the shift takes the decay-rate defect from -2.2e-2 to +2.3e-3 relative to the quiescent wave.
fn stress<L: Lattice3<27>>(f: &[f64; Q], force: [f64; 3], ws: f64, wb: f64) -> [f64; 6] {
    let (rho, mom) = raw_moments::<L>(f);
    let u = half_force_velocity(rho, mom, force);
    let feq = equilibrium::<L>(rho, u, THETA0);
    let mut p = neq_second::<L>(f, &feq);
    p[0] += force[0] * u[0];
    p[1] += force[1] * u[1];
    p[2] += force[2] * u[2];
    p[3] += 0.5 * (force[0] * u[1] + force[1] * u[0]);
    p[4] += 0.5 * (force[0] * u[2] + force[2] * u[0]);
    p[5] += 0.5 * (force[1] * u[2] + force[2] * u[1]);
    let tr = p[0] + p[1] + p[2];
    let a = 1.0 - ws / 2.0;
    let b = (1.0 - wb / 2.0) * tr / 3.0;
    [
        -a * (p[0] - tr / 3.0) - b,
        -a * (p[1] - tr / 3.0) - b,
        -a * (p[2] - tr / 3.0) - b,
        -a * p[3],
        -a * p[4],
        -a * p[5],
    ]
}

/// Fill `out[i] = f(i)` over `workers` scoped threads with fixed contiguous chunks; the value of every cell is a
/// pure function of read-only inputs, so the bits do not depend on the worker count.
fn par_fill<T: Send>(out: &mut [T], workers: usize, f: &(impl Fn(usize) -> T + Sync)) {
    let n = out.len();
    // wasm32 cannot spawn threads: the serial path gives the same bits (the pass is order-independent)
    if workers <= 1 || n < 2 || cfg!(target_arch = "wasm32") {
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

/// Sources of the method of manufactured solutions, in lattice units, for the cell with integer index `x` at lattice time `t`
/// (steps): `(S_rho, S_m, S_E)` of `d rho / dt + div(rho u) = S_rho`, `d(rho u) / dt + div(rho u u + p I - sigma) = S_m` and
/// `dE / dt + div((E + p) u - sigma . u - kappa grad theta) = S_E`. The mass source is applied at the half step (added to the
/// populations with the lattice weights, half before and half after the stream, so it carries no momentum), the momentum source joins the force at both time levels and
/// the energy source enters the trapezoid with the flux divergence. A grid without sources is bitwise unchanged.
/// For ThermochemicalNsfGas, E and S_E use its fixed sensible chart. Converting an absolute NASA9 source requires
/// subtracting energy_reference_j_kg * S_rho / specific_energy_unit_j_kg from S_E; composition remains frozen.
pub type SourceFn = Box<dyn Fn([usize; 3], f64) -> (f64, [f64; 3], f64) + Send + Sync>;

/// The hybrid grid. See the module documentation for the formulation and its conservation class.
pub struct NsfGrid3<L: Lattice3<27>> {
    spec: NsfSpec<NsfGasMedium, String>,
    law: NsfLaw,
    /// Local source response only; empty for the bit-preserving GasSpec route.
    material: Vec<Option<NsfCellMaterial>>,
    sgs: Option<NsfSgsConfig>,
    /// Accepted-state closure; empty in the unchanged molecular route.
    sgs_transport: Vec<NsfSgsTransport>,
    n: usize,
    solid: Vec<bool>,
    /// Pre-collision populations (first moment `sum c f = rho u - F/2`).
    f: Vec<[f64; Q]>,
    e: Vec<f64>,
    theta: Vec<f64>,
    force: Vec<[f64; 3]>,
    div: Vec<f64>,
    wall_out: f64,
    open: Option<NsfCharacteristicAxis>,
    open_faces: Option<OpenFaces>,
    reservoir_faces: [[Option<open::OpenGhost>; 2]; 3],
    open_out: BoundaryFlux,
    last_boundary_exchange: Option<[NsfBoundaryExchange; 2]>,
    steps: usize,
    cumulative: EnergyLedger,
    /// Velocity of the bounce-back wall at the low (index 0) and high (index 1) edge of each non-periodic axis.
    edge_velocity: [[[f64; 3]; 2]; 3],
    moving_walls: bool,
    sources: Option<SourceFn>,
    /// Sources at the current time level (empty without sources).
    src_now: Vec<(f64, [f64; 3], f64)>,
    marker: core::marker::PhantomData<L>,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Build the grid. `initial(cell) -> (rho_lattice, velocity_lattice, temperature_K)`; solid cells are given by
    /// `solid`. Every refusal (cap, window, Mach, temperature ratio, non-finite) happens here and in `step`.
    ///
    /// # Errors
    /// [`NsfError`] for any refused input or state.
    pub fn new<G: Into<NsfGasMedium>, S: Into<String>>(
        spec: NsfSpec<G, S>,
        solid: impl Fn([usize; 3]) -> bool,
        initial: impl Fn([usize; 3]) -> (f64, [f64; 3], f64),
    ) -> Result<Self, NsfError> {
        Self::new_with_sources(spec, solid, initial, None)
    }

    /// As [`NsfGrid3::new`] with manufactured-solution sources ([`SourceFn`]); the momentum source of time 0 is part of the
    /// initial force, so the half-force velocity of the initial state is the requested one.
    ///
    /// # Errors
    /// [`NsfError`] for any refused input or state.
    pub fn new_with_sources<G: Into<NsfGasMedium>, S: Into<String>>(
        spec: NsfSpec<G, S>,
        solid: impl Fn([usize; 3]) -> bool,
        initial: impl Fn([usize; 3]) -> (f64, [f64; 3], f64),
        sources: Option<SourceFn>,
    ) -> Result<Self, NsfError> {
        let spec = spec.map_gas(Into::into);
        let n = spec
            .dims
            .iter()
            .try_fold(1usize, |a, &d| a.checked_mul(d))
            .ok_or(NsfError::Input("dimension overflow"))?;
        if spec.dims.contains(&0) {
            return Err(NsfError::Input("dimensions must be positive"));
        }
        if n > spec.maximum_cells {
            return Err(NsfError::Cap {
                requested: n,
                maximum: spec.maximum_cells,
            });
        }
        let ok = spec.dx_m.is_finite()
            && spec.dx_m > 0.0
            && spec.maximum_mach.is_finite()
            && spec.maximum_mach > 0.0
            && spec.maximum_temperature_ratio.is_finite()
            && spec.maximum_temperature_ratio > 0.0
            && spec.bulk_to_shear_ratio.is_finite()
            && spec.bulk_to_shear_ratio > 0.0
            && !spec.bulk_ratio_source.is_empty()
            && spec.workers >= 1
            && spec.transport_points >= 2
            && spec.temperature_window_k[0] < spec.temperature_window_k[1];
        if spec.maximum_mach.is_finite() && spec.maximum_mach > MAX_QUALIFIED_MACH {
            return Err(NsfError::MachBound {
                requested: spec.maximum_mach,
                qualified: MAX_QUALIFIED_MACH,
            });
        }
        if !ok {
            return Err(NsfError::Input(
                "dx, Mach and temperature bounds, a cited positive bulk-to-shear ratio, workers, table points, window",
            ));
        }
        if let WallThermal::Isothermal(t) = spec.wall {
            if !(t.is_finite()
                && t >= spec.temperature_window_k[0]
                && t <= spec.temperature_window_k[1])
            {
                return Err(NsfError::Input(
                    "isothermal wall temperature must lie in the declared window",
                ));
            }
        }
        let state = spec.gas.reference(&spec)?;
        let r = state.specific_gas_constant;
        let dt_s = spec.dx_m / (3.0 * r * state.temperature).sqrt();
        let sc = dt_s / (state.density * spec.dx_m * spec.dx_m);
        let mut mu = Vec::with_capacity(spec.transport_points);
        let mut kappa = Vec::with_capacity(spec.transport_points);
        for i in 0..spec.transport_points {
            let t = spec.temperature_window_k[0]
                + (spec.temperature_window_k[1] - spec.temperature_window_k[0]) * i as f64
                    / (spec.transport_points - 1) as f64;
            let g = spec.gas.at(t, spec.pressure_pa)?;
            mu.push(g.dynamic_viscosity * sc);
            kappa.push(g.thermal_conductivity * sc / r);
        }
        let max_nu = mu.iter().fold(0.0f64, |a, &v| a.max(v)); // rho0 = 1 in lattice units
        if max_nu > MAX_VERIFIED_LATTICE_VISCOSITY {
            return Err(NsfError::Viscosity {
                observed: max_nu,
                maximum: MAX_VERIFIED_LATTICE_VISCOSITY,
            });
        }
        let law = NsfLaw {
            gamma: state.gamma,
            gas_constant: r,
            rho_ref_kg_m3: state.density,
            dx_m: spec.dx_m,
            dt_s,
            temperature_ref_k: state.temperature,
            bulk_to_shear_ratio: spec.bulk_to_shear_ratio,
            t_low: spec.temperature_window_k[0],
            t_high: spec.temperature_window_k[1],
            mu,
            kappa,
            source_reference_transport: match &spec.gas {
                NsfGasMedium::CaloricallyPerfect(_) => None,
                NsfGasMedium::Thermochemical(_) => Some([
                    state.dynamic_viscosity * sc,
                    state.thermal_conductivity * sc / r,
                ]),
            },
        };
        let cv = 1.0 / (law.gamma - 1.0);
        let at = |i: usize| {
            [
                i % spec.dims[0],
                (i / spec.dims[0]) % spec.dims[1],
                i / (spec.dims[0] * spec.dims[1]),
            ]
        };
        let mut rho = vec![1.0; n];
        let mut vel = vec![[0.0; 3]; n];
        let mut theta = vec![THETA0; n];
        let mut solid_v = vec![false; n];
        let mut e = vec![cv * THETA0; n];
        let mut material = match &spec.gas {
            NsfGasMedium::CaloricallyPerfect(_) => Vec::new(),
            NsfGasMedium::Thermochemical(_) => vec![None; n],
        };
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
                return Err(NsfError::Window {
                    cell: i,
                    temperature_k: t_k,
                    low_k: spec.temperature_window_k[0],
                    high_k: spec.temperature_window_k[1],
                });
            }
            let th = law.theta_of_kelvin(t_k);
            let local = spec.gas.resolve(&law, r0, th)?;
            let gamma = local.as_ref().map_or(law.gamma, |m| m.state.gamma);
            check_state(&spec, gamma, i, r0, u, th)?;
            rho[i] = r0;
            vel[i] = u;
            theta[i] = th;
            let internal = local
                .as_ref()
                .map_or(cv * th, |m| m.sensible_energy_lattice);
            e[i] = r0 * (internal + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
            if let Some(slot) = material.get_mut(i) {
                *slot = local;
            }
        }
        let mut grid = NsfGrid3 {
            spec,
            law,
            material,
            n,
            solid: solid_v,
            sgs: None,
            sgs_transport: Vec::new(),
            f: vec![[0.0; Q]; n],
            e,
            theta,
            force: vec![[0.0; 3]; n],
            div: vec![0.0; n],
            wall_out: 0.0,
            open: None,
            open_faces: None,
            reservoir_faces: [[None; 2]; 3],
            open_out: BoundaryFlux::default(),
            last_boundary_exchange: None,
            steps: 0,
            cumulative: EnergyLedger::default(),
            edge_velocity: [[[0.0; 3]; 2]; 3],
            moving_walls: false,
            sources,
            src_now: Vec::new(),
            marker: core::marker::PhantomData,
        };
        grid.force = grid.force_field(&rho, &grid.theta, None);
        if grid.sources.is_some() {
            grid.src_now = grid.source_field(0.0);
            for i in 0..n {
                if !grid.solid[i] {
                    for a in 0..3 {
                        grid.force[i][a] += grid.src_now[i].1[a];
                    }
                }
            }
        }
        // Initialize local physical equilibrium in the full Guo chart.
        // feq(rho, u - F/(2 rho)) has the correct first moment but leaves a
        // spurious physical second moment F F / (4 rho). Subtract S(u,F)/2
        // instead, keeping both the requested velocity and zero initial stress.
        for i in 0..n {
            grid.f[i] = if grid.solid[i] {
                equilibrium::<L>(1.0, [0.0; 3], THETA0)
            } else {
                let mut f = equilibrium::<L>(rho[i], vel[i], THETA0);
                rebase_force::<L>(&mut f, vel[i], grid.force[i].map(|force| -force));
                f
            };
        }
        let (rho_g, mom) = grid.macro_fields(&grid.f);
        let u = grid.velocity_field(&rho_g, &mom, &grid.force);
        grid.material = grid.material_field(&rho_g, &grid.theta)?;
        let sig = grid.stress_field(
            &grid.f,
            &grid.force,
            &rho_g,
            &grid.theta,
            &grid.material,
            &[],
        );
        let (div, wall, _) =
            grid.flux_divergence(&grid.e, &rho_g, &u, &grid.theta, &sig, &[], None)?;
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

    /// The sources at lattice time `t` for every fluid cell (zero without sources).
    fn source_field(&self, t: f64) -> Vec<(f64, [f64; 3], f64)> {
        let mut out = vec![(0.0, [0.0; 3], 0.0); self.n];
        if let Some(src) = &self.sources {
            let solid = &self.solid;
            par_fill(&mut out, self.spec.workers, &|i| {
                if solid[i] {
                    (0.0, [0.0; 3], 0.0)
                } else {
                    src(self.at(i), t)
                }
            });
        }
        out
    }

    /// Give the bounce-back wall at the low (`high = false`) or high edge of the non-periodic `axis` a velocity [lattice units,
    /// tangential to the wall]. A link that crosses the wall gains the momentum `2 w_q rho (c_q . u_w) / theta0` (Ladd), and
    /// the energy flux through the wall face gains the viscous work `-sigma_{a b} u_{w b}` with `sigma` of the cell beside the
    /// wall (first order in dx for a developing flow, exact for a uniform stress). Solid cells inside the domain stay at
    /// rest. The flux history is restarted at the current state.
    ///
    /// # Errors
    /// [`NsfError::Input`] for a periodic axis, a non-finite or wall-normal velocity, or a speed above the declared Mach bound.
    pub fn set_edge_velocity(
        &mut self,
        axis: usize,
        high: bool,
        velocity: [f64; 3],
    ) -> Result<(), NsfError> {
        if axis >= 3
            || self.spec.periodic[axis]
            || self.open.as_ref().is_some_and(|open| open.axis == axis)
            || self.reservoir_faces[axis][usize::from(high)].is_some()
        {
            return Err(NsfError::Input(
                "an edge wall exists only on a non-periodic axis",
            ));
        }
        if !velocity.iter().all(|v| v.is_finite()) || velocity[axis] != 0.0 {
            return Err(NsfError::Input(
                "a wall velocity must be finite and tangential to the wall",
            ));
        }
        let speed =
            (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2])
                .sqrt();
        if speed / (self.law.gamma * THETA0).sqrt() > self.spec.maximum_mach {
            return Err(NsfError::Input("wall speed above the declared Mach bound"));
        }
        let old_edge_velocity = self.edge_velocity;
        let old_moving_walls = self.moving_walls;
        self.edge_velocity[axis][usize::from(high)] = velocity;
        self.moving_walls = self
            .edge_velocity
            .iter()
            .flatten()
            .any(|u| u.iter().any(|&v| v != 0.0));
        let (rho, mom) = self.macro_fields(&self.f);
        let u = self.velocity_field(&rho, &mom, &self.force);
        let transport = match self.sgs_field(&rho, &u, &self.theta, &self.material) {
            Ok(transport) => transport,
            Err(error) => {
                self.edge_velocity = old_edge_velocity;
                self.moving_walls = old_moving_walls;
                return Err(error);
            }
        };
        let sig = self.stress_field(
            &self.f,
            &self.force,
            &rho,
            &self.theta,
            &self.material,
            &transport,
        );
        let (div, wall, open) = match self.flux_divergence(
            &self.e,
            &rho,
            &u,
            &self.theta,
            &sig,
            &transport,
            self.open_faces.as_ref(),
        ) {
            Ok(history) => history,
            Err(error) => {
                self.edge_velocity = old_edge_velocity;
                self.moving_walls = old_moving_walls;
                return Err(error);
            }
        };
        self.sgs_transport = transport;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(())
    }

    fn macro_fields(&self, f: &[[f64; Q]]) -> (Vec<f64>, Vec<[f64; 3]>) {
        let mut rm = vec![(1.0, [0.0; 3]); self.n];
        let solid = &self.solid;
        par_fill(&mut rm, self.spec.workers, &|i| {
            if solid[i] {
                (1.0, [0.0; 3])
            } else {
                raw_moments::<L>(&f[i])
            }
        });
        (
            rm.iter().map(|v| v.0).collect(),
            rm.iter().map(|v| v.1).collect(),
        )
    }

    fn velocity_field(&self, rho: &[f64], mom: &[[f64; 3]], force: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut out = vec![[0.0; 3]; self.n];
        let solid = &self.solid;
        par_fill(&mut out, self.spec.workers, &|i| {
            if solid[i] {
                [0.0; 3]
            } else {
                half_force_velocity(rho[i], mom[i], force[i])
            }
        });
        out
    }

    fn rates(&self, rho: f64, theta: f64, material: Option<&NsfCellMaterial>) -> (f64, f64) {
        let mu = material.map_or_else(|| self.law.mu_lattice(theta), |m| m.mu_lattice);
        let nu = mu / (rho * THETA0);
        (
            1.0 / (nu + 0.5),
            1.0 / (0.5 + 1.5 * self.spec.bulk_to_shear_ratio * nu),
        )
    }

    fn material_field(
        &self,
        rho: &[f64],
        theta: &[f64],
    ) -> Result<Vec<Option<NsfCellMaterial>>, NsfError> {
        let resolve = |i: usize| {
            if self.solid[i] {
                return Ok(None);
            }
            if !(rho[i].is_finite() && rho[i] > 0.0 && theta[i].is_finite() && theta[i] > 0.0) {
                return Err(NsfError::NonFinite { cell: i });
            }
            let temperature_k = self.law.kelvin_of_theta(theta[i]);
            let [low_k, high_k] = self.spec.temperature_window_k;
            if temperature_k < low_k || temperature_k > high_k {
                return Err(NsfError::Window {
                    cell: i,
                    temperature_k,
                    low_k,
                    high_k,
                });
            }
            self.spec.gas.resolve(&self.law, rho[i], theta[i])
        };
        if matches!(&self.spec.gas, NsfGasMedium::CaloricallyPerfect(_)) {
            // Keep the legacy path allocation-free while admitting its
            // actual local kinematic viscosity at each corrector state.
            for i in 0..self.n {
                resolve(i)?;
            }
            return Ok(Vec::new());
        }
        let mut values = vec![Ok(None); self.n];
        par_fill(&mut values, self.spec.workers, &resolve);
        values.into_iter().collect()
    }

    /// `F = -grad p'`, `p' = rho (theta - theta0)`, in conservative face form: fourth-order face values where the four
    /// stencil cells are fluid, the two-point value beside a solid cell, the cell's own value at a wall face.
    fn force_field(&self, rho: &[f64], theta: &[f64], faces: Option<&OpenFaces>) -> Vec<[f64; 3]> {
        let pp: Vec<f64> = (0..self.n)
            .map(|i| {
                if self.solid[i] {
                    0.0
                } else {
                    rho[i] * (theta[i] - THETA0)
                }
            })
            .collect();
        let solid = &self.solid;
        let face = |i: usize, a: usize, s: i64| -> f64 {
            // face value on the +a side (s = 1) or -a side (s = -1) of cell i
            let x = self.at(i);
            let Some(j) = self.fluid_at(x, a, s) else {
                if let Some(face) = self.characteristic_face(faces, i, a, s) {
                    return face.rho * (face.theta - THETA0);
                }
                if let Some(b) = self.open_at(i, a, s) {
                    return 0.5 * (pp[i] + b.rho * (b.theta - THETA0));
                }
                return pp[i];
            };
            let xj = self.at(j);
            let far_i = self.fluid_at(x, a, -s);
            let far_j = self.fluid_at(xj, a, s);
            match (far_i, far_j) {
                (Some(m), Some(p)) => {
                    (7.0 / 12.0) * (pp[i] + pp[j]) - (1.0 / 12.0) * (pp[m] + pp[p])
                }
                _ => 0.5 * (pp[i] + pp[j]),
            }
        };
        let mut out = vec![[0.0; 3]; self.n];
        par_fill(&mut out, self.spec.workers, &|i| {
            if solid[i] {
                return [0.0; 3];
            }
            [
                -(face(i, 0, 1) - face(i, 0, -1)),
                -(face(i, 1, 1) - face(i, 1, -1)),
                -(face(i, 2, 1) - face(i, 2, -1)),
            ]
        });
        out
    }

    fn stress_field(
        &self,
        f: &[[f64; Q]],
        force: &[[f64; 3]],
        rho: &[f64],
        theta: &[f64],
        material: &[Option<NsfCellMaterial>],
        transport: &[NsfSgsTransport],
    ) -> Vec<[f64; 6]> {
        let mut out = vec![[0.0; 6]; self.n];
        let solid = &self.solid;
        par_fill(&mut out, self.spec.workers, &|i| {
            if solid[i] {
                [0.0; 6]
            } else {
                let (ws, wb) = self.sgs_rates(
                    rho[i],
                    theta[i],
                    material.get(i).and_then(Option::as_ref),
                    transport.get(i),
                );
                stress::<L>(&f[i], force[i], ws, wb)
            }
        });
        out
    }

    /// Divergence of the total-energy flux at ONE time level and the total heat leaving through the walls.
    fn flux_divergence(
        &self,
        et: &[f64],
        rho: &[f64],
        u: &[[f64; 3]],
        theta: &[f64],
        sig: &[[f64; 6]],
        transport: &[NsfSgsTransport],
        faces: Option<&OpenFaces>,
    ) -> Result<(Vec<f64>, f64, BoundaryFlux), NsfError> {
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
            [
                (et[i] + p) * u[i][0] - su[0],
                (et[i] + p) * u[i][1] - su[1],
                (et[i] + p) * u[i][2] - su[2],
            ]
        });
        let wall_theta = match self.spec.wall {
            WallThermal::Adiabatic => None,
            WallThermal::Isothermal(t) => Some(self.law.theta_of_kelvin(t)),
        };
        let face = |i: usize, a: usize| -> Result<f64, NsfError> {
            let x = self.at(i);
            let Some(j) = self.fluid_at(x, a, 1) else {
                return Ok(0.0);
            };
            let mut kf = self.spec.gas.face_conductivity(
                &self.law,
                0.5 * (theta[i] + theta[j]),
                0.5 * (self.law.pressure_pa(rho[i], theta[i])
                    + self.law.pressure_pa(rho[j], theta[j])),
            )?;
            // One face value shared by both cells preserves the total-energy
            // telescoping ledger. Never add modeled heat flux through an
            // impermeable wall; its molecular half-cell law below is retained.
            if let (Some(left), Some(right)) = (transport.get(i), transport.get(j)) {
                kf += 0.5 * (left.turbulent_conductivity_lu + right.turbulent_conductivity_lu);
            }
            let wide = match (self.fluid_at(x, a, -1), self.fluid_at(self.at(j), a, 1)) {
                (Some(m), Some(p)) => Some((m, p)),
                _ => None,
            };
            Ok(match wide {
                Some((m, p)) => {
                    let gf =
                        (7.0 / 12.0) * (g[i][a] + g[j][a]) - (1.0 / 12.0) * (g[m][a] + g[p][a]);
                    let dt4 = (27.0 * (theta[j] - theta[i]) - (theta[p] - theta[m])) / 24.0;
                    gf - kf * dt4
                }
                None => 0.5 * (g[i][a] + g[j][a]) - kf * (theta[j] - theta[i]),
            })
        };
        let mut plus: Vec<Result<[f64; 3], NsfError>> = vec![Ok([0.0; 3]); n];
        par_fill(&mut plus, self.spec.workers, &|i| {
            if solid[i] {
                return Ok([0.0; 3]);
            }
            Ok([face(i, 0)?, face(i, 1)?, face(i, 2)?])
        });
        let plus: Vec<[f64; 3]> = plus.into_iter().collect::<Result<_, _>>()?;
        // wall face: flux in the +axis direction. Heat through the half cell; the viscous work `-sigma_{a b} u_{w b}` of a moving
        // edge wall (zero normal velocity, so no convective part; solid cells and still walls do no work).
        let wall_flux = |i: usize, a: usize, s: i64| -> Result<f64, NsfError> {
            if let Some(face) = self.characteristic_face(faces, i, a, s) {
                return self
                    .open_energy_flux(face, i, a, s, rho, theta, sig)
                    .map(|f| f.0);
            }
            if let Some(b) = self.open_at(i, a, s) {
                let sg = sig[i];
                let row = match a {
                    0 => [sg[0], sg[3], sg[4]],
                    1 => [sg[3], sg[1], sg[5]],
                    _ => [sg[4], sg[5], sg[2]],
                };
                let stress_work = row[0] * b.u[0] + row[1] * b.u[1] + row[2] * b.u[2];
                let gb = (b.energy + b.rho * b.theta) * b.u[a] - stress_work;
                let k = self.spec.gas.face_conductivity(
                    &self.law,
                    0.5 * (theta[i] + b.theta),
                    0.5 * (self.law.pressure_pa(rho[i], theta[i])
                        + self.law.pressure_pa(b.rho, b.theta)),
                )?;
                return Ok(0.5 * (g[i][a] + gb) - (s as f64) * k * (b.theta - theta[i]));
            }
            let heat = match wall_theta {
                None => 0.0,
                Some(tw) => {
                    // Same local pressure as the adjacent cell; only the wall
                    // temperature changes across its half-cell resistance.
                    let kf = self.spec.gas.face_conductivity(
                        &self.law,
                        0.5 * (theta[i] + tw),
                        self.law.pressure_pa(rho[i], theta[i]),
                    )?;
                    if s > 0 {
                        -kf * (tw - theta[i]) / 0.5
                    } else {
                        -kf * (theta[i] - tw) / 0.5
                    }
                }
            };
            if !self.moving_walls || self.shifted(self.at(i), a, s).is_some() {
                return Ok(heat);
            }
            let uw = self.edge_velocity[a][usize::from(s > 0)];
            let sg = sig[i];
            let row = match a {
                0 => [sg[0], sg[3], sg[4]],
                1 => [sg[3], sg[1], sg[5]],
                _ => [sg[4], sg[5], sg[2]],
            };
            Ok(heat - (row[0] * uw[0] + row[1] * uw[1] + row[2] * uw[2]))
        };
        let mut div: Vec<Result<f64, NsfError>> = vec![Ok(0.0); n];
        par_fill(&mut div, self.spec.workers, &|i| {
            if solid[i] {
                return Ok(0.0);
            }
            let x = self.at(i);
            let mut d = 0.0;
            for a in 0..3 {
                d += match self.fluid_at(x, a, 1) {
                    Some(_) => plus[i][a],
                    None => wall_flux(i, a, 1)?,
                };
                d -= match self.fluid_at(x, a, -1) {
                    Some(m) => plus[m][a],
                    None => wall_flux(i, a, -1)?,
                };
            }
            Ok(d)
        });
        let div = div.into_iter().collect::<Result<Vec<_>, _>>()?;
        let mut out = 0.0;
        let mut open = BoundaryFlux::default();
        for i in 0..n {
            if solid[i] {
                continue;
            }
            let x = self.at(i);
            for a in 0..3 {
                if self.fluid_at(x, a, 1).is_none() {
                    if let Some(face) = self.characteristic_face(faces, i, a, 1) {
                        let (total, internal) =
                            self.open_energy_flux(face, i, a, 1, rho, theta, sig)?;
                        open.total_out[1] += total;
                        open.internal_out[1] += internal;
                    } else if self.open_at(i, a, 1).is_some() {
                        open.prescribed_out += wall_flux(i, a, 1)?;
                    } else {
                        out += wall_flux(i, a, 1)?;
                    }
                }
                if self.fluid_at(x, a, -1).is_none() {
                    if let Some(face) = self.characteristic_face(faces, i, a, -1) {
                        let (total, internal) =
                            self.open_energy_flux(face, i, a, -1, rho, theta, sig)?;
                        open.total_out[0] -= total;
                        open.internal_out[0] -= internal;
                    } else if self.open_at(i, a, -1).is_some() {
                        open.prescribed_out -= wall_flux(i, a, -1)?;
                    } else {
                        out -= wall_flux(i, a, -1)?;
                    }
                }
            }
        }
        Ok((div, out, open))
    }

    fn stream(&self, post: &[[f64; Q]], incoming: Option<&IncomingPopulations>) -> Vec<[f64; Q]> {
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
                let mut y = x;
                // the first axis whose step leaves a non-periodic domain: the edge wall that the link bounces from
                let mut edge: Option<(usize, usize)> = None;
                for a in 0..3 {
                    if c[a] != 0 && edge.is_none() {
                        match self.shifted(y, a, -i64::from(c[a])) {
                            Some(j) => y = self.at(j),
                            None => edge = Some((a, usize::from(c[a] < 0))),
                        }
                    }
                }
                out[q] = match edge {
                    Some((a, side))
                        if self.open.as_ref().is_some_and(|open| open.axis == a)
                            && self.characteristic_link_is_open(x, c, -1) =>
                    {
                        incoming
                            .expect("open stream has prepared incoming populations")
                            .at(x, self.spec.dims, side)[q]
                    }
                    None if !self.solid[self.index(y)] => post[self.index(y)][q],
                    Some((a, side)) if self.reservoir_faces[a][side].is_some() => {
                        let b = self.reservoir_faces[a][side].expect("guarded reservoir face");
                        let (rho, mom) = raw_moments::<L>(&post[i]);
                        let local_u = mom.map(|m| m / rho);
                        // Post-collision populations carry rho*u + F/2. Extrapolate
                        // that half-force shift as well as the non-equilibrium part.
                        let ghost_u =
                            core::array::from_fn(|k| b.u[k] + 0.5 * self.force[i][k] / b.rho);
                        equilibrium::<L>(b.rho, ghost_u, THETA0)[q] + post[i][q]
                            - equilibrium::<L>(rho, local_u, THETA0)[q]
                    }
                    _ => {
                        let mut v = post[i][L::OPP[q]];
                        if let (true, Some((a, side))) = (self.moving_walls, edge) {
                            let uw = self.edge_velocity[a][side];
                            let cu = f64::from(c[0]) * uw[0]
                                + f64::from(c[1]) * uw[1]
                                + f64::from(c[2]) * uw[2];
                            let rho: f64 = post[i].iter().sum();
                            v += 2.0 * L::W[q] * rho * cu / THETA0;
                        }
                        v
                    }
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
                    let (ws, wb) = self.sgs_rates(
                        rho[i],
                        theta[i],
                        self.material.get(i).and_then(Option::as_ref),
                        self.sgs_transport.get(i),
                    );
                    collide::<L>(&f[i], force[i], ws, wb)
                }
            });
        }
        let have_sources = self.sources.is_some();
        // The mass increment of the step is injected half before and half after the stream: the stream's pressure gradient then
        // sees the mid-step density, as the momentum source (the average of the two levels) does. All of it before the stream
        // would leave a first-order mismatch of half theta0 grad(S) in the momentum.
        let mass_half = if have_sources {
            self.source_field(self.steps as f64 + 0.5)
        } else {
            Vec::new()
        };
        if have_sources {
            for i in 0..n {
                if !self.solid[i] {
                    for q in 0..Q {
                        post[i][q] += 0.5 * L::W[q] * mass_half[i].0;
                    }
                }
            }
        }
        let src_next = self.source_field(self.steps as f64 + 1.0);
        let incoming = self.open_incoming_populations(&post, &rho0, &u0);
        let boundary_links = self.open_link_exchange(&post, incoming.as_ref());
        let mut next = self.stream(&post, incoming.as_ref());
        if have_sources {
            for i in 0..n {
                if !self.solid[i] {
                    for q in 0..Q {
                        next[i][q] += 0.5 * L::W[q] * mass_half[i].0;
                    }
                }
            }
        }
        let (rho1, mom1) = self.macro_fields(&next);
        let mut th1 = self.theta.clone();
        let mut u_guess = if self.open.is_some() {
            self.velocity_field(&rho1, &mom1, &self.force)
        } else {
            Vec::new()
        };
        let mut result = None;
        for _ in 0..COUPLING_PASSES {
            let faces1 = self.resolve_open_faces(self.open.as_ref(), &rho1, &u_guess, &th1)?;
            let mut f1 = self.force_field(&rho1, &th1, faces1.as_ref());
            if have_sources {
                for i in 0..n {
                    if !self.solid[i] {
                        for a in 0..3 {
                            f1[i][a] += src_next[i].1[a];
                        }
                    }
                }
            }
            let u1 = self.velocity_field(&rho1, &mom1, &f1);
            let trial_material = self.material_field(&rho1, &th1)?;
            let trial_transport = self.sgs_field(&rho1, &u1, &th1, &trial_material)?;
            let et_trial: Vec<f64> = (0..n)
                .map(|i| {
                    let internal = trial_material
                        .get(i)
                        .and_then(Option::as_ref)
                        .map_or(cv * th1[i], |m| m.sensible_energy_lattice);
                    rho1[i]
                        * (internal
                            + 0.5
                                * (u1[i][0] * u1[i][0] + u1[i][1] * u1[i][1] + u1[i][2] * u1[i][2]))
                })
                .collect();
            let sig = self.stress_field(&next, &f1, &rho1, &th1, &trial_material, &trial_transport);
            let (div1, wall1, open1) = self.flux_divergence(
                &et_trial,
                &rho1,
                &u1,
                &th1,
                &sig,
                &trial_transport,
                faces1.as_ref(),
            )?;
            let mut material_new = trial_material;
            let mut e1 = self.e.clone();
            let mut th_new = th1.clone();
            for i in 0..n {
                if self.solid[i] {
                    continue;
                }
                e1[i] = self.e[i] - 0.5 * (div1[i] + self.div[i]);
                if have_sources {
                    e1[i] += 0.5 * (src_next[i].2 + self.src_now[i].2);
                }
                let ke = 0.5
                    * rho1[i]
                    * (u1[i][0] * u1[i][0] + u1[i][1] * u1[i][1] + u1[i][2] * u1[i][2]);
                match &self.spec.gas {
                    NsfGasMedium::CaloricallyPerfect(_) => {
                        // Preserve the original division/order exactly.
                        th_new[i] = (e1[i] - ke) / (rho1[i] * cv);
                    }
                    NsfGasMedium::Thermochemical(source) => {
                        let (theta, local) =
                            source.reconstruct(&self.law, rho1[i], (e1[i] - ke) / rho1[i])?;
                        th_new[i] = theta;
                        material_new[i] = Some(local);
                    }
                }
            }
            th1 = th_new;
            if self.open.is_some() {
                u_guess.clone_from(&u1);
            }
            result = Some((f1, u1, e1, div1, wall1, open1, faces1, material_new));
        }
        let (force1, u1, e1, div1, wall1, open1, faces1, material_new) =
            result.expect("at least one coupling pass");
        for i in 0..n {
            if self.solid[i] {
                continue;
            }
            if !(rho1[i].is_finite()
                && rho1[i] > 0.0
                && e1[i].is_finite()
                && th1[i].is_finite()
                && th1[i] > 0.0)
            {
                return Err(NsfError::NonFinite { cell: i });
            }
            let gamma = material_new
                .get(i)
                .and_then(Option::as_ref)
                .map_or(self.law.gamma, |m| m.state.gamma);
            check_state(&self.spec, gamma, i, rho1[i], u1[i], th1[i])?;
            let t_k = self.law.kelvin_of_theta(th1[i]);
            if t_k < self.spec.temperature_window_k[0] || t_k > self.spec.temperature_window_k[1] {
                return Err(NsfError::Window {
                    cell: i,
                    temperature_k: t_k,
                    low_k: self.spec.temperature_window_k[0],
                    high_k: self.spec.temperature_window_k[1],
                });
            }
            if material_new.is_empty() {
                // The final caloric corrector changes temperature after the
                // trial material check. Admit the actual published legacy
                // endpoint too; source reconstruction already does this.
                self.spec.gas.resolve(&self.law, rho1[i], th1[i])?;
            }
        }
        // The last corrector changes temperature after faces1 was evaluated.
        // Admit the actual endpoint too, without replacing the common
        // predictor faces that supplied its force and conservative FV flux.
        let _ = self.resolve_open_faces(self.open.as_ref(), &rho1, &u1, &th1)?;
        // The next collision must use the actual accepted force-corrected
        // velocity and source state, after the final temperature corrector.
        let transport_new = self.sgs_field(&rho1, &u1, &th1, &material_new)?;
        let unit = self.law.energy_unit_j();
        let (mut ke0, mut ke1, mut e_old, mut e_new, mut pw) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            if self.solid[i] {
                continue;
            }
            ke0 +=
                0.5 * rho0[i] * (u0[i][0] * u0[i][0] + u0[i][1] * u0[i][1] + u0[i][2] * u0[i][2]);
            ke1 +=
                0.5 * rho1[i] * (u1[i][0] * u1[i][0] + u1[i][1] * u1[i][1] + u1[i][2] * u1[i][2]);
            e_old += self.e[i];
            e_new += e1[i];
        }
        let p0 = self.pressure_work_density(&rho0, &u0, &self.theta, self.open_faces.as_ref());
        let p1 = self.pressure_work_density(&rho1, &u1, &th1, faces1.as_ref());
        for i in 0..n {
            if !self.solid[i] {
                pw += 0.5 * (p0[i] + p1[i]);
            }
        }
        let wall_in = -0.5 * (wall1 + self.wall_out);
        let open_in = -0.5
            * (open1.total_out.iter().sum::<f64>() + self.open_out.total_out.iter().sum::<f64>());
        let open_internal_in = -0.5
            * (open1.internal_out.iter().sum::<f64>()
                + self.open_out.internal_out.iter().sum::<f64>());
        let prescribed_in = -0.5 * (open1.prescribed_out + self.open_out.prescribed_out);
        // Use the exact source quadrature from the native energy update and
        // the already-sampled fields; never query a callback again for the audit.
        let mut source_in = 0.0;
        let mut source_carries_mechanics = false;
        if have_sources {
            for i in 0..n {
                if !self.solid[i] {
                    source_in += 0.5 * (src_next[i].2 + self.src_now[i].2);
                    source_carries_mechanics |= mass_half[i].0 != 0.0
                        || self.src_now[i].1.iter().any(|&value| value != 0.0)
                        || src_next[i].1.iter().any(|&value| value != 0.0);
                }
            }
        }
        let d_total = e_new - e_old;
        let d_kin = ke1 - ke0;
        let d_int = d_total - d_kin;
        let pressure_work = -pw;
        let viscous = if self.open.is_some() {
            d_int - wall_in - open_internal_in - pressure_work
        } else {
            let remainder = d_int - wall_in - pressure_work;
            if prescribed_in == 0.0 {
                remainder
            } else {
                remainder - prescribed_in
            }
        };
        let remainder = if source_in == 0.0 {
            viscous
        } else {
            viscous - source_in
        };
        let (viscous, unresolved) = if source_carries_mechanics {
            (0.0, remainder)
        } else {
            (remainder, 0.0)
        };
        let balance = if self.open.is_some() {
            d_total - wall_in - open_in
        } else {
            d_total - wall_in
        };
        // Retain closed/characteristic arithmetic, including signed zero.
        let balance = if prescribed_in == 0.0 {
            balance
        } else {
            balance - prescribed_in
        };
        let balance = if source_in == 0.0 {
            balance
        } else {
            balance - source_in
        };
        let ledger = EnergyLedger {
            internal_j: unit * (e_new - ke1),
            kinetic_j: unit * ke1,
            viscous_heating_j: unit * viscous,
            unresolved_internal_exchange_j: unit * unresolved,
            pressure_work_j: unit * pressure_work,
            wall_heat_j: unit * wall_in,
            open_boundary_energy_j: unit * open_in,
            open_energy_j: unit * prescribed_in,
            source_energy_j: unit * source_in,
            residual_rel: balance.abs() / e_new.abs().max(f64::MIN_POSITIVE),
        };
        let boundary_exchange = self.finish_boundary_exchange(
            boundary_links,
            self.open_faces.as_ref(),
            faces1.as_ref(),
            open1,
        );
        self.f = next;
        self.e = e1;
        self.material = material_new;
        self.sgs_transport = transport_new;
        if have_sources {
            self.src_now = src_next;
        }
        self.theta = th1;
        self.force = force1;
        self.div = div1;
        self.wall_out = wall1;
        self.open_out = open1;
        self.open_faces = faces1;
        self.last_boundary_exchange = boundary_exchange;
        self.steps += 1;
        self.cumulative.viscous_heating_j += ledger.viscous_heating_j;
        self.cumulative.unresolved_internal_exchange_j += ledger.unresolved_internal_exchange_j;
        self.cumulative.source_energy_j += ledger.source_energy_j;
        self.cumulative.pressure_work_j += ledger.pressure_work_j;
        self.cumulative.wall_heat_j += ledger.wall_heat_j;
        self.cumulative.open_boundary_energy_j += ledger.open_boundary_energy_j;
        self.cumulative.open_energy_j += ledger.open_energy_j;
        self.cumulative.internal_j = ledger.internal_j;
        self.cumulative.kinetic_j = ledger.kinetic_j;
        Ok(ledger)
    }

    /// `p div u` per fluid cell with the fourth-order central divergence (second order beside walls).
    fn pressure_work_density(
        &self,
        rho: &[f64],
        u: &[[f64; 3]],
        theta: &[f64],
        faces: Option<&OpenFaces>,
    ) -> Vec<f64> {
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
                    (Some(m), Some(p)) => match (
                        self.fluid_at(self.at(m), a, -1),
                        self.fluid_at(self.at(p), a, 1),
                    ) {
                        (Some(mm), Some(pp)) => {
                            (8.0 * (u[p][a] - u[m][a]) - (u[pp][a] - u[mm][a])) / 12.0
                        }
                        _ => 0.5 * (u[p][a] - u[m][a]),
                    },
                    (None, Some(p)) => match self.characteristic_face(faces, i, a, -1) {
                        Some(f) => 0.5 * (u[p][a] + u[i][a]) - f.u[a],
                        None => self
                            .open_at(i, a, -1)
                            .map_or(0.5 * (u[p][a] + u[i][a]), |b| 0.5 * (u[p][a] - b.u[a])),
                    },
                    (Some(m), None) => match self.characteristic_face(faces, i, a, 1) {
                        Some(f) => f.u[a] - 0.5 * (u[m][a] + u[i][a]),
                        None => self
                            .open_at(i, a, 1)
                            .map_or(-0.5 * (u[m][a] + u[i][a]), |b| 0.5 * (b.u[a] - u[m][a])),
                    },
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
    /// Total energy in joules in the medium's fixed sensible chart, summed over
    /// fluid cells in fixed order. NASA9 absolute formation energy is excluded.
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
        (
            rho,
            half_force_velocity(rho, mom, self.force[i]),
            self.law.kelvin_of_theta(self.theta[i]),
        )
    }
    /// Retained source owner or original GasSpec.
    #[must_use]
    pub fn medium(&self) -> &NsfGasMedium {
        &self.spec.gas
    }

    /// Re-evaluate the original source at this accepted cell's actual density
    /// and temperature. The returned receipt includes the same local EOS,
    /// heat capacities, sound speed and transport used by the flow.
    ///
    /// Returns None for the legacy GasSpec path. Solid/out-of-grid cells refuse.
    /// No material selection, composition or fixed query coordinate is changed.
    pub fn thermochemical_state(
        &self,
        x: [usize; 3],
    ) -> Result<Option<fs_material::thermochemical_gas::ResolvedThermochemicalGasState>, NsfError>
    {
        if x.iter().zip(self.spec.dims).any(|(&i, n)| i >= n) {
            return Err(NsfError::Input("material state cell outside the grid"));
        }
        let i = self.index(x);
        if self.solid[i] {
            return Err(NsfError::Input("solid cell has no gas material state"));
        }
        let (rho, _) = raw_moments::<L>(&self.f[i]);
        self.spec.gas.source_state(&self.law, rho, self.theta[i])
    }

    /// The scaling and reference-pressure diagnostic transport laws in force.
    #[must_use]
    pub fn law(&self) -> &NsfLaw {
        &self.law
    }
    /// Retained citation of the bulk/shear viscosity ratio used by this grid.
    /// Ownership is independent of the caller's input buffer lifetime.
    #[must_use]
    pub fn bulk_ratio_source(&self) -> &str {
        &self.spec.bulk_ratio_source
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

fn check_state(
    spec: &NsfSpec<NsfGasMedium, String>,
    gamma: f64,
    cell: usize,
    rho: f64,
    u: [f64; 3],
    theta: f64,
) -> Result<(), NsfError> {
    if !(theta.is_finite() && theta > 0.0 && rho.is_finite() && rho > 0.0) {
        return Err(NsfError::NonFinite { cell });
    }
    let ratio = theta / THETA0 - 1.0;
    if ratio.abs() > spec.maximum_temperature_ratio {
        return Err(NsfError::TemperatureRatio {
            cell,
            observed: ratio.abs(),
            maximum: spec.maximum_temperature_ratio,
        });
    }
    let mach = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt() / (gamma * theta).sqrt();
    if !mach.is_finite() || mach > spec.maximum_mach {
        return Err(NsfError::Mach {
            cell,
            observed: mach,
            maximum: spec.maximum_mach,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q27::D3q27;

    fn raw_source(u: [f64; 3], force: [f64; 3]) -> [f64; Q] {
        let uf: f64 = u.iter().zip(force).map(|(u, force)| u * force).sum();
        std::array::from_fn(|q| {
            let c = D3q27::E[q].map(f64::from);
            let cu: f64 = c.iter().zip(u).map(|(c, u)| c * u).sum();
            let cf: f64 = c.iter().zip(force).map(|(c, force)| c * force).sum();
            D3q27::W[q] * (cf / THETA0 + (cu * cf - THETA0 * uf) / THETA0.powi(2))
        })
    }

    fn half_source_second(u: [f64; 3], force: [f64; 3]) -> [f64; 6] {
        [
            u[0] * force[0],
            u[1] * force[1],
            u[2] * force[2],
            0.5 * (u[0] * force[1] + u[1] * force[0]),
            0.5 * (u[0] * force[2] + u[2] * force[0]),
            0.5 * (u[1] * force[2] + u[2] * force[1]),
        ]
    }

    #[test]
    fn g0_initial_forced_gas_has_requested_velocity_and_no_spurious_stress() {
        let velocity = [0.07, -0.03, 0.02];
        let source_force = [0.001, -0.0017, 0.0004];
        for workers in [1, 3] {
            let spec = NsfSpec {
                dims: [8, 2, 2],
                periodic: [true; 3],
                gas: GasSpec::dry_air_ussa1976(),
                temperature_k: 293.15,
                pressure_pa: 101325.0,
                dx_m: 5e-5,
                temperature_window_k: [273.15, 333.15],
                maximum_mach: MAX_QUALIFIED_MACH,
                maximum_temperature_ratio: 0.2,
                bulk_to_shear_ratio: 0.73,
                bulk_ratio_source: "Cramer 2012, N2 at 293 K",
                maximum_cells: 32,
                wall: WallThermal::Adiabatic,
                workers,
                transport_points: 17,
            };
            let initial = |x: [usize; 3]| {
                let temperature = 293.15 + 3.0 * x[0] as f64;
                (1.2, velocity, temperature)
            };
            let grid = NsfGrid3::<D3q27>::new_with_sources(
                spec,
                |_| false,
                initial,
                Some(Box::new(move |_, _| (0.0, source_force, 0.0))),
            )
            .unwrap();
            // Both a nonuniform thermal-pressure force and a nonparallel
            // manufactured force enter the actual constructor's chart.
            for i in 0..grid.n {
                let (rho, u, temperature) = grid.cell(grid.at(i));
                assert!((rho - 1.2).abs() < 2e-14);
                for a in 0..3 {
                    assert!((u[a] - velocity[a]).abs() < 2e-14);
                }
                assert!((temperature - initial(grid.at(i)).2).abs() < 1e-12);
                let (ws, wb) = grid.rates(rho, grid.theta[i], None);
                let sigma = stress::<D3q27>(&grid.f[i], grid.force[i], ws, wb);
                assert!(
                    sigma.iter().all(|component| component.abs() < 2e-14),
                    "cell={i} workers={workers} initial stress={sigma:?}"
                );
                // Independently check the physical second moment itself,
                // without allowing a near-two relaxation rate to hide error.
                let eq = equilibrium::<D3q27>(rho, velocity, THETA0);
                let second = neq_second::<D3q27>(&grid.f[i], &eq);
                let half_source = half_source_second(velocity, grid.force[i]);
                for component in 0..6 {
                    assert!((second[component] + half_source[component]).abs() < 2e-14);
                }
            }
        }
    }

    #[test]
    fn guo_source_respects_bulk_relaxation_for_uniform_acceleration() {
        // G0: the transformed populations of a locally unstrained flow are
        // feq - S/2. Collision must produce second moment Pi_eq + (uF+Fu)/2
        // independently of both viscosities, while supplying exactly F to
        // momentum. A single source prefactor leaves a spurious isotropic
        // stress of (wb-ws) (u.F)/3 when the two rates differ.
        let rho = 1.3;
        let u = [0.07, -0.03, 0.02];
        let force = [0.001, -0.0017, 0.0004];
        let eq = equilibrium::<D3q27>(rho, u, THETA0);
        let source = raw_source(u, force);
        let initial = std::array::from_fn(|q| eq[q] - 0.5 * source[q]);
        let wanted_second = half_source_second(u, force);
        for (ws, wb) in [(1.82, 1.15), (0.75, 1.95), (1.2, 1.2)] {
            let initial_stress = stress::<D3q27>(&initial, force, ws, wb);
            assert!(initial_stress.iter().all(|v| v.abs() < 2e-14));
            let next = collide::<D3q27>(&initial, force, ws, wb);
            let (density, momentum) = raw_moments::<D3q27>(&next);
            assert!((density - rho).abs() < 2e-14);
            for axis in 0..3 {
                assert!((momentum[axis] - rho * u[axis] - 0.5 * force[axis]).abs() < 2e-14);
            }
            let second = neq_second::<D3q27>(&next, &eq);
            for component in 0..6 {
                assert!(
                    (second[component] - wanted_second[component]).abs() < 2e-14,
                    "mode {component}, ws={ws}, wb={wb}: {} vs {}",
                    second[component],
                    wanted_second[component]
                );
            }
        }
    }

    #[test]
    fn bulk_and_shear_stress_relax_independently_in_forced_flow() {
        // G0: a supplied physical second moment must relax in its own trace
        // and deviatoric eigenspaces. Nonparallel u and F exercise all source
        // components; an isotropic-only patch must not change the shear modes.
        let rho = 0.9;
        let u = [-0.04, 0.03, 0.05];
        let force = [0.0006, -0.0002, 0.0007];
        let physical = [2e-4, -1e-4, 3e-4, 1.5e-4, -2.5e-4, 5e-5];
        let trace = physical[0] + physical[1] + physical[2];
        let eq = equilibrium::<D3q27>(rho, u, THETA0);
        let source = raw_source(u, force);
        let initial = std::array::from_fn(|q| {
            let [x, y, z] = D3q27::E[q].map(f64::from);
            let cpc = x * x * physical[0]
                + y * y * physical[1]
                + z * z * physical[2]
                + 2.0 * (x * y * physical[3] + x * z * physical[4] + y * z * physical[5]);
            eq[q] + D3q27::W[q] * (cpc - THETA0 * trace) / (2.0 * THETA0.powi(2)) - 0.5 * source[q]
        });
        let half_source = half_source_second(u, force);
        for (ws, wb) in [(1.87, 0.64), (0.72, 1.91), (1.3, 1.3)] {
            let next = collide::<D3q27>(&initial, force, ws, wb);
            let second = neq_second::<D3q27>(&next, &eq);
            for component in 0..6 {
                let isotropic = if component < 3 { trace / 3.0 } else { 0.0 };
                let expected = (1.0 - ws) * (physical[component] - isotropic)
                    + (1.0 - wb) * isotropic
                    + half_source[component];
                assert!(
                    (second[component] - expected).abs() < 2e-14,
                    "mode {component}, ws={ws}, wb={wb}: {} vs {expected}",
                    second[component]
                );
            }
        }
    }
}
