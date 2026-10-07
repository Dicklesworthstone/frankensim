//! Steady incompressible Navier–Stokes by finite volumes on the staggered
//! (MAC) faces of a [`VoxelDomain`], coupled by SIMPLEC (Van Doormaal &
//! Raithby, *Numer. Heat Transfer* 7, 1984) — the Cartesian-voxel method of
//! the electronics-cooling class of tools, sharing the grid, the face-flux
//! layout and the convection schemes of the conjugate energy equation.
//!
//! # Discretization
//!
//! Velocity component `a` lives on the faces normal to axis `a`; pressure
//! lives on cells. Each face velocity has a cell-sized staggered control
//! volume: convective fluxes through its sides interpolate the neighbouring
//! face velocities, diffusion is `mu A / dx` (`2 mu A / dx` to a no-slip wall
//! half a cell away), and convection uses Patankar's power-law (or upwind)
//! coefficients. Faces touching a solid voxel are blocked (zero velocity);
//! a blocked transverse neighbour is a no-slip wall half a cell away
//! (staircase geometry).
//!
//! Domain faces are [`FvBoundary`]: `Wall` (no-slip, optionally moving
//! tangentially), `Symmetry` (free slip), `Inlet` (prescribed velocity) and
//! `Outlet` (pressure zero, zero normal gradient of velocity). An outlet
//! face is a momentum unknown like an interior face: the cell beyond the
//! boundary mirrors the inside cell's velocities and carries the ghost
//! pressure `-p_inside`, so the face holds pressure zero, and the pressure
//! correction uses the same linearization (`u' = 2 d p'_inside`). The steady
//! state is therefore one fixed point of the discrete equations, independent
//! of the relaxation, also where the outflow is not fully developed.
//!
//! # Fans and resistances
//!
//! [`InternalFan`]s raise the static pressure across a planar patch of
//! interior faces by their curve's value at the flow through the patch; the
//! source is linearized about each sweep's flow (`S(u) = s A dp(Q0) + A
//! dp'(Q0) A_fan (u - u0)`, exact at the fixed point) so the fan-system
//! loop is damped by the curve's own slope. [`FlowResistance::Planar`]
//! drops it by `1/2 rho K |u| u` across a patch (grilles, perforated
//! plates); [`FlowResistance::Volume`] is a Darcy–Forchheimer porous block.
//! Both are implicit in `a_P` (Picard in `|u|`), so they also enter
//! SIMPLEC's `d`. They act on momentum only (no heat).
//!
//! # SIMPLEC iteration
//!
//! Momentum is under-relaxed (`a_P / alpha`) and solved per component with
//! ILU(0)-BiCGStab; the pressure correction uses `d = A / (a_P / alpha -
//! sum a_nb)` and is solved with AMG-preconditioned CG (or ILU(0)-BiCGStab)
//! on the Jacobi-scaled system (components without an outlet are pinned);
//! velocities and pressure are
//! corrected with no pressure under-relaxation. Convergence requires two steady residuals of the same
//! iterate below the tolerance: the largest cell mass imbalance (before
//! correction) over the largest face mass flux, and, per velocity component,
//! the residual `||b - A u||` of the Jacobi-scaled momentum system at the
//! velocities entering the iteration (under-relaxation cancels there, so it
//! is the residual of the unrelaxed equations), relative to the flow's
//! momentum scale `sqrt(rows) max_c rms(b_c)` (a component that vanishes by
//! symmetry has a round-off `b`; its own norm is no scale). Inner solves
//! only reduce their entry residuals (by `momentum_tolerance` and
//! `pressure_tolerance`), so the outer residuals carry the convergence claim
//! and a converged report never rests on a skipped inner solve.
//!
//! # No-claim boundaries
//!
//! Constant-property incompressible flow, steady here and unsteady through
//! [`super::simple_unsteady`]; turbulence only through the algebraic LVEL
//! closure (a laminar solution above transition is a laminar idealization,
//! not a prediction); buoyancy only as the body force the natural-convection
//! drivers set; staircase voxel walls; power-law convection (first order at
//! high cell Péclet numbers: numerical diffusion is not bounded here). One
//! run makes no mesh-convergence claim.

use fs_exec::CancelGate;
use fs_sparse::Coo;
use fs_sparse::precond::{SaAmg, pcg};

use super::domain::{FluidProperties, VoxelDomain};
use super::energy::ConvectionScheme;
use super::flow::{FlowField, scale_rows};
use super::krylov::bicgstab_ilu0;
use super::turbulence::{TURBULENT_PRANDTL, law_of_the_wall_ratios, wall_distance};
use super::{ChtError, finite, finite_positive, poll};
use crate::d3q19::Face3;

/// Flow rule on one domain face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FvBoundary {
    /// No-slip wall moving with the given (tangential) velocity, m/s.
    Wall {
        /// Wall velocity; the normal component must be zero.
        velocity: [f64; 3],
    },
    /// Free slip: zero normal velocity and zero tangential shear.
    Symmetry,
    /// Prescribed inflow velocity, m/s (the normal component must point into
    /// the domain).
    Inlet {
        /// Inflow velocity vector.
        velocity: [f64; 3],
    },
    /// Prescribed (zero) pressure with zero-gradient velocity.
    Outlet,
}

impl FvBoundary {
    /// A stationary no-slip wall.
    #[must_use]
    pub const fn wall() -> Self {
        Self::Wall { velocity: [0.0; 3] }
    }
}

/// Piecewise-linear fan characteristic: static pressure rise (Pa) against
/// volume flow (m^3/s), 2 to 8 points with strictly increasing flow and
/// non-increasing pressure; linear extrapolation beyond the ends.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FanCurve {
    points: [(f64, f64); 8],
    len: usize,
}

impl FanCurve {
    /// Admit a characteristic.
    ///
    /// # Errors
    /// [`ChtError::InvalidInput`] for fewer than 2 or more than 8 points,
    /// non-finite values, negative or non-increasing flows, or a rising
    /// pressure.
    pub fn new(points: &[(f64, f64)]) -> Result<Self, ChtError> {
        let refuse = |reason: &str| ChtError::InvalidInput {
            field: "fan.curve",
            reason: reason.to_string(),
        };
        if !(2..=8).contains(&points.len()) {
            return Err(refuse("needs 2 to 8 (flow, pressure) points"));
        }
        for (i, &(q, p)) in points.iter().enumerate() {
            if !(q.is_finite() && p.is_finite()) || q < 0.0 {
                return Err(refuse(
                    "flows must be finite and non-negative, pressures finite",
                ));
            }
            if i > 0 && (q <= points[i - 1].0 || p > points[i - 1].1) {
                return Err(refuse("flow must increase and pressure must not rise"));
            }
        }
        let mut stored = [(0.0, 0.0); 8];
        stored[..points.len()].copy_from_slice(points);
        Ok(Self {
            points: stored,
            len: points.len(),
        })
    }

    fn segment(&self, q: f64) -> ((f64, f64), (f64, f64)) {
        let pts = &self.points[..self.len];
        let i = pts[1..self.len - 1]
            .iter()
            .take_while(|(qi, _)| *qi <= q)
            .count();
        (pts[i], pts[i + 1])
    }

    /// Static pressure rise at flow `q`, Pa.
    #[must_use]
    pub fn pressure(&self, q: f64) -> f64 {
        let ((q0, p0), (q1, p1)) = self.segment(q);
        p0 + (p1 - p0) * (q - q0) / (q1 - q0)
    }

    /// `d(pressure)/dq` on the segment containing `q` (non-positive).
    #[must_use]
    pub fn slope(&self, q: f64) -> f64 {
        let ((q0, p0), (q1, p1)) = self.segment(q);
        (p1 - p0) / (q1 - q0)
    }
}

/// A fan on one inlet face: the face's uniform normal inflow velocity is
/// solved for so that the mean pressure of the fluid layer behind the face
/// equals the fan's static pressure rise from ambient (pressure zero) at
/// the delivered flow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FanInlet {
    /// The face; its rule must be [`FvBoundary::Inlet`] (the declared
    /// velocity is the initial guess and fixes the direction).
    pub face: Face3,
    /// The characteristic.
    pub curve: FanCurve,
}

/// A planar patch of interior faces: the faces normal to `axis` on the face
/// plane `index` (between cells `index - 1` and `index` along `axis`, so
/// `1..n[axis]`), over the cells `lo[t]..hi[t]` of the two transverse axes
/// `t` (the `axis` entries of `lo`/`hi` are ignored). Faces touching a solid
/// voxel are skipped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FacePatch {
    /// Normal axis (0, 1, 2).
    pub axis: usize,
    /// Face-plane index along `axis`.
    pub index: usize,
    /// First transverse cells (inclusive).
    pub lo: [usize; 3],
    /// Last transverse cells (exclusive).
    pub hi: [usize; 3],
}

impl FacePatch {
    /// Face coordinates of the patch whose two adjacent cells are fluid.
    fn open_faces(&self, domain: &VoxelDomain) -> Vec<[usize; 3]> {
        let n = domain.dims();
        let (u, w) = ((self.axis + 1) % 3, (self.axis + 2) % 3);
        let mut faces = Vec::new();
        for j in self.lo[w]..self.hi[w].min(n[w]) {
            for i in self.lo[u]..self.hi[u].min(n[u]) {
                let mut f = [0usize; 3];
                f[self.axis] = self.index;
                f[u] = i;
                f[w] = j;
                let mut minus = f;
                minus[self.axis] -= 1;
                if domain.is_fluid(cell_of(n, minus)) && domain.is_fluid(cell_of(n, f)) {
                    faces.push(f);
                }
            }
        }
        faces
    }

    fn admit(&self, domain: &VoxelDomain, field: &'static str) -> Result<(), ChtError> {
        let n = domain.dims();
        let refuse = |reason: String| Err(ChtError::InvalidInput { field, reason });
        if self.axis > 2 {
            return refuse(format!("axis {} is not 0, 1 or 2", self.axis));
        }
        if self.index == 0 || self.index >= n[self.axis] {
            return refuse(format!(
                "face plane {} must be interior (1..{} along axis {})",
                self.index, n[self.axis], self.axis
            ));
        }
        for t in [(self.axis + 1) % 3, (self.axis + 2) % 3] {
            if self.lo[t] >= self.hi[t] || self.hi[t] > n[t] {
                return refuse(format!(
                    "transverse cells {}..{} on axis {t} must be non-empty within 0..{}",
                    self.lo[t], self.hi[t], n[t]
                ));
            }
        }
        Ok(())
    }
}

/// A flow resistance inside the fluid (a sink of momentum; no heat).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FlowResistance {
    /// Thin planar resistance (grille, perforated plate, filter sheet): the
    /// static pressure drops by `1/2 rho K |u| u` across each face of the
    /// patch, `u` the face (approach) velocity.
    Planar {
        /// The faces.
        patch: FacePatch,
        /// Loss coefficient `K` (non-negative).
        loss_coefficient: f64,
    },
    /// Porous block over the cells `lo..hi` (a heat-exchanger core, filter,
    /// dense component array): per axis `a`,
    /// `-dp/dx_a = mu / kappa_a u_a + 1/2 rho C_a |u_a| u_a`
    /// (Darcy–Forchheimer on the superficial velocity, orthotropic).
    Volume {
        /// First cells (inclusive).
        lo: [usize; 3],
        /// Last cells (exclusive).
        hi: [usize; 3],
        /// Permeability per axis, m^2 (`f64::INFINITY`: no viscous term).
        permeability_m2: [f64; 3],
        /// Inertial (Forchheimer) coefficient per axis, 1/m.
        inertial_per_m: [f64; 3],
    },
}

impl FlowResistance {
    /// Idelchik's loss coefficient of a thin sharp-edged perforated plate
    /// with free-area ratio `f` (Handbook of Hydraulic Resistance, diagram
    /// 8-1, turbulent): `K = (1 + 0.707 sqrt(1 - f) - f)^2 / f^2` on the
    /// approach velocity; `None` outside `(0, 1]`.
    #[must_use]
    pub fn perforated_plate_loss(free_area_ratio: f64) -> Option<f64> {
        let f = free_area_ratio;
        if !(f > 0.0 && f <= 1.0) {
            return None;
        }
        let k = 0.707f64.mul_add(fs_math::det::sqrt(1.0 - f), 1.0 - f) / f;
        Some(k * k)
    }
}

/// A fan inside the domain (an axial fan in an enclosure): across each face
/// of the patch the static pressure rises by the curve's value at the flow
/// through the whole patch, in the blowing direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InternalFan {
    /// The faces (at least one must be open).
    pub patch: FacePatch,
    /// Blowing toward `+axis` (true) or `-axis` (false).
    pub blows_positive: bool,
    /// The characteristic (flow in the blowing direction).
    pub curve: FanCurve,
}

/// Linear solver of the (symmetric positive definite) pressure correction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PressureSolver {
    /// Conjugate gradients preconditioned by fs-sparse's smoothed-aggregation
    /// AMG V-cycle; the hierarchy is rebuilt every
    /// [`AMG_REBUILD_SWEEPS`] sweeps (the matrix changes slowly), and a
    /// failed solve falls back to ILU(0)-BiCGStab.
    #[default]
    AmgCg,
    /// ILU(0)-preconditioned BiCGStab on the Jacobi-scaled system.
    IluBicgstab,
}

/// Turbulence closure of the momentum equations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Turbulence {
    /// Laminar (molecular viscosity only).
    #[default]
    Laminar,
    /// LVEL algebraic eddy viscosity from the wall distance and local speed
    /// (see [`super::turbulence`]).
    Lvel,
}

/// Implicit time discretization of the unsteady momentum equations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TimeScheme {
    /// First-order backward Euler.
    BackwardEuler,
    /// Second-order backward differentiation (the first step is backward
    /// Euler).
    #[default]
    Bdf2,
}

/// Sweeps between AMG hierarchy rebuilds for [`PressureSolver::AmgCg`].
pub const AMG_REBUILD_SWEEPS: usize = 10;

/// SIMPLEC controls.
#[derive(Debug, Clone, PartialEq)]
pub struct SimpleConfig {
    /// One rule per domain face in [`Face3::ALL`] order.
    pub faces: [FvBoundary; 6],
    /// Convection scheme for momentum.
    pub scheme: ConvectionScheme,
    /// Momentum under-relaxation in `(0, 1)`.
    pub velocity_relaxation: f64,
    /// Outer iteration budget.
    pub max_iterations: usize,
    /// Convergence tolerance (mass imbalance and velocity change, relative).
    pub tolerance: f64,
    /// Factor by which each inner momentum solve reduces its entry residual.
    pub momentum_tolerance: f64,
    /// Factor by which each pressure-correction solve reduces its entry
    /// residual (the correction starts from zero, so this is its residual
    /// relative to the mass imbalance it corrects).
    pub pressure_tolerance: f64,
    /// Optional fan on one inlet face (operating point solved).
    pub fan: Option<FanInlet>,
    /// Pressure-correction linear solver.
    pub pressure_solver: PressureSolver,
    /// Turbulence closure.
    pub turbulence: Turbulence,
    /// Flow resistances (grilles, porous blocks).
    pub resistances: Vec<FlowResistance>,
    /// Fans inside the domain.
    pub internal_fans: Vec<InternalFan>,
}

impl SimpleConfig {
    /// Defaults with the given face rules.
    #[must_use]
    pub const fn new(faces: [FvBoundary; 6]) -> Self {
        Self {
            faces,
            scheme: ConvectionScheme::PowerLaw,
            velocity_relaxation: 0.7,
            max_iterations: 3000,
            tolerance: 1e-6,
            momentum_tolerance: 1e-1,
            pressure_tolerance: 1e-2,
            fan: None,
            pressure_solver: PressureSolver::AmgCg,
            turbulence: Turbulence::Laminar,
            resistances: Vec::new(),
            internal_fans: Vec::new(),
        }
    }
}

/// Run evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct SimpleReport {
    /// Outer iterations.
    pub iterations: usize,
    /// Largest cell mass imbalance before the last correction over the
    /// largest face mass flux.
    pub mass_residual: f64,
    /// Largest steady momentum residual over the three components at the
    /// last iteration (see [`simple_flow`]).
    pub momentum_residual: f64,
    /// Largest per-cell |net outflow| of the returned fluxes, m^3/s.
    pub max_divergence_m3_s: f64,
    /// Total inflow through inlet faces and backflow through outlets
    /// (openings), m^3/s.
    pub inflow_m3_s: f64,
    /// Total outflow through outlet faces (openings), m^3/s.
    pub outflow_m3_s: f64,
    /// Largest cell Reynolds number `|u| dx / nu`.
    pub max_cell_reynolds: f64,
    /// Total Krylov iterations of the momentum solves.
    pub momentum_krylov_iterations: usize,
    /// Total Krylov iterations of the pressure-correction solves.
    pub pressure_krylov_iterations: usize,
    /// Fan operating point `(flow m^3/s, static pressure rise Pa, relative
    /// mismatch between the curve and the mean inlet-layer pressure)`, when a
    /// fan is declared.
    pub fan: Option<(f64, f64, f64)>,
    /// Operating point `(flow m^3/s in the blowing direction, static
    /// pressure rise Pa)` of each internal fan, in declaration order.
    pub internal_fans: Vec<(f64, f64)>,
}

/// Steady finite-volume flow.
#[derive(Debug, Clone, PartialEq)]
pub struct FvFlow {
    /// Face fluxes for the conjugate energy equation.
    pub field: FlowField,
    /// Cell-centred velocity (mean of each cell's two face velocities per
    /// axis), m/s; zero in solids.
    pub velocity_m_s: Vec<[f64; 3]>,
    /// Cell pressure, Pa (outlet reference zero); zero in solids.
    pub pressure_pa: Vec<f64>,
    /// Eddy viscosity per cell, m^2/s (zeros when laminar).
    pub eddy_viscosity_m2_s: Vec<f64>,
    /// The eddy viscosity the energy equation should use per cell, m^2/s:
    /// the law-of-the-wall secant value in wall-adjacent cells, the LVEL
    /// value elsewhere (zeros when laminar).
    pub heat_eddy_viscosity_m2_s: Vec<f64>,
    /// Run evidence.
    pub report: SimpleReport,
}

impl FvFlow {
    /// Turbulent conductivity `rho c_p nu_t / Pr_t` per cell, W/(m K), for
    /// `ThermalSetup::eddy_conductivity_w_m_k` (zeros when laminar).
    #[must_use]
    pub fn eddy_conductivity(&self, fluid: &FluidProperties) -> Vec<f64> {
        self.heat_eddy_viscosity_m2_s
            .iter()
            .map(|nu_t| fluid.volumetric_heat_capacity() * nu_t / TURBULENT_PRANDTL)
            .collect()
    }

    /// Mean pressure over the fluid cells of layer `index` along `axis`, Pa.
    #[must_use]
    pub fn mean_pressure(&self, domain: &VoxelDomain, axis: usize, index: usize) -> Option<f64> {
        let (mut sum, mut count) = (0.0f64, 0usize);
        for c in 0..domain.cell_count() {
            if domain.is_fluid(c) && domain.coords(c)[axis] == index {
                sum += self.pressure_pa[c];
                count += 1;
            }
        }
        (count > 0).then(|| sum / count as f64)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Kind {
    /// Solved by the momentum equation (row index).
    Unknown(usize),
    /// Prescribed value (blocked faces carry zero).
    Fixed(f64),
    /// Outlet face (row index): solved by its own momentum equation with
    /// zero normal gradient and a ghost pressure `-p_inside`, so the face
    /// carries the boundary pressure zero.
    Outlet(usize),
}

/// Face lattice of one velocity component.
#[derive(Debug, Clone)]
struct Component {
    dims: [usize; 3],
    kind: Vec<Kind>,
    unknowns: Vec<usize>,
}

impl Component {
    fn index(&self, f: [usize; 3]) -> usize {
        (f[2] * self.dims[1] + f[1]) * self.dims[0] + f[0]
    }

    fn coords(&self, index: usize) -> [usize; 3] {
        let x = index % self.dims[0];
        let yz = index / self.dims[0];
        [x, yz % self.dims[1], yz / self.dims[1]]
    }
}

pub(super) struct Solver<'a> {
    domain: &'a VoxelDomain,
    config: &'a SimpleConfig,
    n: [usize; 3],
    rho: f64,
    mu: f64,
    area: f64,
    comps: [Component; 3],
    /// Face velocities per component, m/s.
    vel: [Vec<f64>; 3],
    pressure: Vec<f64>,
    /// SIMPLEC `d` per face (0 for non-unknown faces until set).
    d: [Vec<f64>; 3],
    /// Krylov iterations spent on momentum and pressure correction.
    momentum_krylov: usize,
    pressure_krylov: usize,
    /// Body force per cell, N/m^3 (empty: none).
    force: Vec<[f64; 3]>,
    /// Face rules in force (a fan updates its inlet velocity).
    faces: [FvBoundary; 6],
    /// Fan state, when declared.
    fan: Option<FanState>,
    /// SIMPLEC iterations taken.
    sweeps: usize,
    /// AMG hierarchy of the pressure correction and the sweep it was built at.
    amg: Option<(SaAmg, usize)>,
    /// LVEL: wall distance per cell (m), eddy viscosity per cell (m^2/s),
    /// the law-of-the-wall secant viscosity ratio `y+ / u+` per cell, and
    /// whether a cell borders a wall; all empty when laminar.
    wall_distance: Vec<f64>,
    eddy_viscosity: Vec<f64>,
    wall_ratio: Vec<f64>,
    wall_adjacent: Vec<bool>,
    /// Planar loss coefficient per face and component (empty: none).
    planar_loss: [Vec<f64>; 3],
    /// Porous `(mu / kappa, rho C / 2)` per cell and axis (empty: none).
    porous: Vec<[(f64, f64); 3]>,
    /// Internal fan index per face and component (`usize::MAX`: none;
    /// empty: no internal fans).
    fan_face: [Vec<usize>; 3],
    internal_fans: Vec<InternalFanState>,
    /// Last RMS momentum right-hand side per component (residual scale).
    momentum_rms: [f64; 3],
    /// Unsteady mode: step, scheme, and the face velocities at the previous
    /// (and, for BDF2, the one before) time level.
    time: Option<TimeLevels>,
}

#[derive(Debug, Clone)]
struct TimeLevels {
    dt: f64,
    scheme: TimeScheme,
    previous: [Vec<f64>; 3],
    older: Option<[Vec<f64>; 3]>,
}

/// Internal fan state: the flow through the patch and the curve's rise and
/// slope there, refreshed at the start of each sweep.
#[derive(Debug, Clone)]
struct InternalFanState {
    axis: usize,
    /// +1 blowing toward +axis, -1 toward -axis.
    sign: f64,
    faces: Vec<usize>,
    /// Open area of the patch, m^2.
    area: f64,
    curve: FanCurve,
    flow: f64,
    rise: f64,
    slope: f64,
}

/// Fan operating-point state.
#[derive(Debug, Clone)]
struct FanState {
    side: usize,
    curve: FanCurve,
    /// Open (fluid) area of the fan face, m^2.
    area: f64,
    /// Fluid cells behind the face.
    cells: Vec<usize>,
    residual: f64,
    flow: f64,
    /// The previous settled `(flow, mean inlet pressure)`, for the secant.
    previous: Option<(f64, f64)>,
}

fn norm(v: &[f64]) -> f64 {
    fs_math::det::sqrt(v.iter().map(|x| x * x).sum::<f64>())
}

fn cell_of(n: [usize; 3], c: [usize; 3]) -> usize {
    (c[2] * n[1] + c[1]) * n[0] + c[0]
}

impl<'a> Solver<'a> {
    #[allow(clippy::too_many_lines)] // one constructor of the solver state
    pub(super) fn new(
        domain: &'a VoxelDomain,
        fluid: &FluidProperties,
        config: &'a SimpleConfig,
    ) -> Self {
        let n = domain.dims();
        let comps = [0, 1, 2].map(|axis| {
            let mut dims = n;
            dims[axis] += 1;
            let count = dims[0] * dims[1] * dims[2];
            let mut kind = vec![Kind::Fixed(0.0); count];
            let mut unknowns = Vec::new();
            let mut comp = Component {
                dims,
                kind: Vec::new(),
                unknowns: Vec::new(),
            };
            for (index, slot) in kind.iter_mut().enumerate() {
                let f = comp.coords(index);
                let along = f[axis];
                let mut minus = f;
                let fluid_at = |c: [usize; 3]| domain.is_fluid(cell_of(n, c));
                if along == 0 || along == n[axis] {
                    let cell = if along == 0 {
                        f
                    } else {
                        minus[axis] -= 1;
                        minus
                    };
                    if !fluid_at(cell) {
                        continue;
                    }
                    let rule = config.faces[2 * axis + usize::from(along != 0)];
                    *slot = match rule {
                        FvBoundary::Wall { .. } | FvBoundary::Symmetry => Kind::Fixed(0.0),
                        FvBoundary::Inlet { velocity } => Kind::Fixed(velocity[axis]),
                        FvBoundary::Outlet => {
                            unknowns.push(index);
                            Kind::Outlet(unknowns.len() - 1)
                        }
                    };
                } else {
                    minus[axis] -= 1;
                    if fluid_at(minus) && fluid_at(f) {
                        *slot = Kind::Unknown(unknowns.len());
                        unknowns.push(index);
                    }
                }
            }
            comp.kind = kind;
            comp.unknowns = unknowns;
            comp
        });
        let vel = [0, 1, 2].map(|a| {
            comps[a]
                .kind
                .iter()
                .map(|k| match *k {
                    Kind::Fixed(v) => v,
                    _ => 0.0,
                })
                .collect::<Vec<f64>>()
        });
        let d = [0, 1, 2].map(|a| vec![0.0; comps[a].kind.len()]);
        let planar_loss = {
            let mut loss = [0, 1, 2].map(|_| Vec::new());
            for resistance in &config.resistances {
                if let FlowResistance::Planar {
                    patch,
                    loss_coefficient,
                } = *resistance
                {
                    let comp = &comps[patch.axis];
                    let slot = &mut loss[patch.axis];
                    if slot.is_empty() {
                        *slot = vec![0.0; comp.kind.len()];
                    }
                    for f in patch.open_faces(domain) {
                        slot[comp.index(f)] += loss_coefficient;
                    }
                }
            }
            loss
        };
        let porous = {
            let mut porous = Vec::new();
            let mu = fluid.density_kg_m3 * fluid.kinematic_viscosity_m2_s;
            for resistance in &config.resistances {
                if let FlowResistance::Volume {
                    lo,
                    hi,
                    permeability_m2,
                    inertial_per_m,
                } = *resistance
                {
                    if porous.is_empty() {
                        porous = vec![[(0.0, 0.0); 3]; domain.cell_count()];
                    }
                    for z in lo[2]..hi[2] {
                        for y in lo[1]..hi[1] {
                            for x in lo[0]..hi[0] {
                                let cell: &mut [(f64, f64); 3] = &mut porous[domain.index(x, y, z)];
                                for a in 0..3 {
                                    cell[a].0 += mu / permeability_m2[a];
                                    cell[a].1 += 0.5 * fluid.density_kg_m3 * inertial_per_m[a];
                                }
                            }
                        }
                    }
                }
            }
            porous
        };
        let fan_face = if config.internal_fans.is_empty() {
            [Vec::new(), Vec::new(), Vec::new()]
        } else {
            let mut map = [0, 1, 2].map(|a| vec![usize::MAX; comps[a].kind.len()]);
            for (i, fan) in config.internal_fans.iter().enumerate() {
                let comp = &comps[fan.patch.axis];
                for f in fan.patch.open_faces(domain) {
                    map[fan.patch.axis][comp.index(f)] = i;
                }
            }
            map
        };
        let internal_fans: Vec<InternalFanState> = config
            .internal_fans
            .iter()
            .map(|fan| {
                let comp = &comps[fan.patch.axis];
                let faces: Vec<usize> = fan
                    .patch
                    .open_faces(domain)
                    .into_iter()
                    .map(|f| comp.index(f))
                    .collect();
                InternalFanState {
                    axis: fan.patch.axis,
                    sign: if fan.blows_positive { 1.0 } else { -1.0 },
                    area: domain.dx() * domain.dx() * faces.len() as f64,
                    faces,
                    curve: fan.curve,
                    flow: 0.0,
                    rise: fan.curve.pressure(0.0),
                    slope: fan.curve.slope(0.0),
                }
            })
            .collect();
        Self {
            domain,
            config,
            n,
            rho: fluid.density_kg_m3,
            mu: fluid.density_kg_m3 * fluid.kinematic_viscosity_m2_s,
            area: domain.dx() * domain.dx(),
            comps,
            vel,
            pressure: vec![0.0; domain.cell_count()],
            d,
            momentum_krylov: 0,
            pressure_krylov: 0,
            force: Vec::new(),
            faces: config.faces,
            sweeps: 0,
            amg: None,
            wall_distance: match config.turbulence {
                Turbulence::Laminar => Vec::new(),
                Turbulence::Lvel => wall_distance(
                    domain,
                    config
                        .faces
                        .map(|rule| matches!(rule, FvBoundary::Wall { .. })),
                ),
            },
            eddy_viscosity: match config.turbulence {
                Turbulence::Laminar => Vec::new(),
                Turbulence::Lvel => vec![0.0; domain.cell_count()],
            },
            wall_ratio: match config.turbulence {
                Turbulence::Laminar => Vec::new(),
                Turbulence::Lvel => vec![1.0; domain.cell_count()],
            },
            wall_adjacent: match config.turbulence {
                Turbulence::Laminar => Vec::new(),
                Turbulence::Lvel => (0..domain.cell_count())
                    .map(|c| {
                        domain.is_fluid(c)
                            && (0..6).any(|f| match domain.neighbor(c, f) {
                                Some(n) => !domain.is_fluid(n),
                                None => matches!(config.faces[f], FvBoundary::Wall { .. }),
                            })
                    })
                    .collect(),
            },
            planar_loss,
            porous,
            fan_face,
            internal_fans,
            momentum_rms: [0.0; 3],
            time: None,
            fan: config.fan.map(|fan| {
                let side = fan.face as usize;
                let axis = side / 2;
                let cells: Vec<usize> = (0..domain.cell_count())
                    .filter(|&c| domain.is_fluid(c) && domain.neighbor(c, side).is_none())
                    .collect();
                let area = domain.dx() * domain.dx() * cells.len() as f64;
                let speed = match config.faces[side] {
                    FvBoundary::Inlet { velocity } => velocity[axis].abs(),
                    _ => 0.0,
                };
                FanState {
                    side,
                    curve: fan.curve,
                    area,
                    cells,
                    // Unmeasured: a full mismatch until the first settled
                    // state (finite, so the divergence guard does not trip).
                    residual: 1.0,
                    flow: speed * area,
                    previous: None,
                }
            }),
        }
    }

    /// Effective dynamic viscosity of `cell` (molecular plus eddy).
    fn mu_eff(&self, cell: [usize; 3]) -> f64 {
        if self.eddy_viscosity.is_empty() {
            self.mu
        } else {
            self.rho
                .mul_add(self.eddy_viscosity[cell_of(self.n, cell)], self.mu)
        }
    }

    /// Wall viscosity of the half cell between `minus`/`plus` and a wall:
    /// the law-of-the-wall secant (molecular when laminar).
    fn wall_mu(&self, minus: [usize; 3], plus: [usize; 3]) -> f64 {
        if self.wall_ratio.is_empty() {
            return self.mu;
        }
        let ratio = 0.5
            * (self.wall_ratio[cell_of(self.n, minus)] + self.wall_ratio[cell_of(self.n, plus)]);
        self.mu * ratio
    }

    /// Mean effective viscosity at the CV edge shared by `minus` and
    /// `plus` across direction `d` (`up`: toward +d): the fluid cells among
    /// the two and their neighbours across the edge.
    fn edge_mu(&self, minus: [usize; 3], plus: [usize; 3], d: usize, up: bool) -> f64 {
        if self.eddy_viscosity.is_empty() {
            return self.mu;
        }
        let (mut sum, mut count) = (self.mu_eff(minus) + self.mu_eff(plus), 2.0);
        for c in [minus, plus] {
            let across = if up {
                (c[d] + 1 < self.n[d]).then(|| {
                    let mut g = c;
                    g[d] += 1;
                    g
                })
            } else {
                (c[d] > 0).then(|| {
                    let mut g = c;
                    g[d] -= 1;
                    g
                })
            };
            if let Some(g) = across
                && self.domain.is_fluid(cell_of(self.n, g))
            {
                sum += self.mu_eff(g);
                count += 1.0;
            }
        }
        sum / count
    }

    /// LVEL: refresh the eddy viscosity from the current cell-centred
    /// speeds (under-relaxed by one half after the first sweep).
    fn update_turbulence(&mut self) {
        if self.eddy_viscosity.is_empty() {
            return;
        }
        let nu = self.mu / self.rho;
        let first = self.sweeps <= 1;
        for c in 0..self.domain.cell_count() {
            if !self.domain.is_fluid(c) {
                continue;
            }
            let at = self.domain.coords(c);
            let mut speed2 = 0.0;
            for a in 0..3 {
                let v = 0.5
                    * (self.vel[a][self.cell_face(a, at, false)]
                        + self.vel[a][self.cell_face(a, at, true)]);
                speed2 += v * v;
            }
            let reynolds = fs_math::det::sqrt(speed2) * self.wall_distance[c].min(1e6) / nu;
            let (tangent, secant) = law_of_the_wall_ratios(reynolds);
            let target = tangent * nu;
            if first {
                self.eddy_viscosity[c] = target;
                self.wall_ratio[c] = secant;
            } else {
                self.eddy_viscosity[c] = 0.5 * (self.eddy_viscosity[c] + target);
                self.wall_ratio[c] = 0.5 * (self.wall_ratio[c] + secant);
            }
        }
    }

    /// Eddy viscosity per cell, m^2/s (zeros when laminar).
    pub(super) fn eddy_viscosity(&self) -> Vec<f64> {
        if self.eddy_viscosity.is_empty() {
            vec![0.0; self.domain.cell_count()]
        } else {
            self.eddy_viscosity.clone()
        }
    }

    /// The eddy viscosity the energy equation uses per cell, m^2/s: the
    /// law-of-the-wall secant `nu (y+/u+ - 1)` in wall-adjacent cells (their
    /// conduction to the wall spans the half cell), LVEL's value elsewhere.
    pub(super) fn heat_eddy_viscosity(&self) -> Vec<f64> {
        if self.eddy_viscosity.is_empty() {
            return vec![0.0; self.domain.cell_count()];
        }
        let nu = self.mu / self.rho;
        (0..self.domain.cell_count())
            .map(|c| {
                if self.wall_adjacent[c] {
                    nu * (self.wall_ratio[c] - 1.0)
                } else {
                    self.eddy_viscosity[c]
                }
            })
            .collect()
    }

    /// Set the uniform normal inflow speed of inlet face `side`.
    fn set_inlet_speed(&mut self, side: usize, speed: f64) {
        let axis = side / 2;
        let inward = if side % 2 == 0 { 1.0 } else { -1.0 };
        if let FvBoundary::Inlet { velocity } = &mut self.faces[side] {
            velocity[axis] = inward * speed;
        }
        let along = if side % 2 == 0 { 0 } else { self.n[axis] };
        for index in 0..self.comps[axis].kind.len() {
            let f = self.comps[axis].coords(index);
            if f[axis] != along || !matches!(self.comps[axis].kind[index], Kind::Fixed(_)) {
                continue;
            }
            let mut cell = f;
            if along != 0 {
                cell[axis] -= 1;
            }
            if self.domain.is_fluid(cell_of(self.n, cell)) {
                self.comps[axis].kind[index] = Kind::Fixed(inward * speed);
                self.vel[axis][index] = inward * speed;
            }
        }
    }

    /// Scale every declared inlet's velocity by `scale` (an inflow schedule
    /// of the unsteady march).
    pub(super) fn scale_inlets(&mut self, scale: f64) {
        for side in 0..6 {
            if let FvBoundary::Inlet { velocity } = self.config.faces[side] {
                let axis = side / 2;
                let inward = if side % 2 == 0 { 1.0 } else { -1.0 };
                self.set_inlet_speed(side, inward * velocity[axis] * scale);
                self.faces[side] = FvBoundary::Inlet {
                    velocity: velocity.map(|v| v * scale),
                };
            }
        }
    }

    /// Secant step of the fan operating point, taken only on a settled
    /// flow: the mean inlet-layer pressure just after an inflow change is
    /// dominated by the pressure-correction transient, so the flow at the
    /// current delivery must first converge to a gate that tightens with the
    /// fan mismatch (finally the solver tolerance). The system-curve slope is
    /// the secant through the last two settled states (`p / Q` at the
    /// first). A passing mismatch leaves the delivery unchanged.
    fn update_fan(&mut self, mass: f64, momentum: f64) {
        let tolerance = self.config.tolerance;
        let Some(fan) = self.fan.as_mut() else {
            return;
        };
        let gate = tolerance.max(1e-3 * fan.residual.min(1.0));
        if mass > gate || momentum > gate {
            return;
        }
        let mean = fan.cells.iter().map(|&c| self.pressure[c]).sum::<f64>()
            / fan.cells.len().max(1) as f64;
        let mismatch = fan.curve.pressure(fan.flow) - mean;
        let scale = fan
            .curve
            .pressure(0.0)
            .abs()
            .max(mean.abs())
            .max(f64::MIN_POSITIVE);
        fan.residual = mismatch.abs() / scale;
        if fan.residual <= tolerance {
            return;
        }
        let system = match fan.previous {
            Some((q, p)) if (fan.flow - q).abs() > f64::EPSILON * fan.flow => {
                (mean - p) / (fan.flow - q)
            }
            _ => mean / fan.flow.max(f64::MIN_POSITIVE),
        }
        .max(0.0);
        fan.previous = Some((fan.flow, mean));
        let step = mismatch / (system - fan.curve.slope(fan.flow)).max(f64::MIN_POSITIVE);
        // Never stop the fan outright: a zero delivery has no system slope.
        fan.flow = (fan.flow + step).max(1e-3 * fan.flow);
        let (side, speed) = (fan.side, fan.flow / fan.area.max(f64::MIN_POSITIVE));
        self.set_inlet_speed(side, speed);
    }

    /// Replace the per-cell body force (N/m^3; zero in solids).
    pub(super) fn set_force(&mut self, force: Vec<[f64; 3]>) {
        debug_assert_eq!(force.len(), self.domain.cell_count());
        self.force = force;
    }

    /// One SIMPLEC iteration: three momentum solves and a pressure
    /// correction. Returns `(mass_residual, momentum_residual)` of the
    /// velocities entering it (see the module docs).
    pub(super) fn sweep(&mut self, gate: &CancelGate) -> Result<(f64, f64), ChtError> {
        self.sweeps += 1;
        self.update_turbulence();
        for fan in &mut self.internal_fans {
            let sum: f64 = fan.faces.iter().map(|&f| self.vel[fan.axis][f]).sum();
            fan.flow = fan.sign * self.area * sum;
            fan.rise = fan.curve.pressure(fan.flow);
            fan.slope = fan.curve.slope(fan.flow);
        }
        let mut steady = 0.0f64;
        for a in 0..3 {
            steady = steady.max(self.momentum(a, gate)?);
            self.guard_finite()?;
        }
        let imbalance = self.correct(self.config.pressure_tolerance, gate)?;
        self.guard_finite()?;
        // The fan residual joins the momentum residual: the operating point
        // is part of the steady state.
        let largest_velocity = self
            .vel
            .iter()
            .flat_map(|v| v.iter())
            .fold(0.0f64, |m, v| m.max(v.abs()))
            .max(f64::MIN_POSITIVE);
        let mass = imbalance / (self.rho * self.area * largest_velocity);
        self.update_fan(mass, steady);
        let fan = self.fan.as_ref().map_or(0.0, |fan| fan.residual);
        Ok((mass, steady.max(fan)))
    }

    /// Enter unsteady mode at the current velocities (time level 0).
    pub(super) fn start_time(&mut self, dt: f64, scheme: TimeScheme) {
        self.time = Some(TimeLevels {
            dt,
            scheme,
            previous: self.vel.clone(),
            older: None,
        });
    }

    /// Close the current time step: its velocities become the previous
    /// level.
    pub(super) fn advance_time(&mut self) {
        if let Some(time) = &mut self.time {
            let previous = std::mem::replace(&mut time.previous, self.vel.clone());
            time.older = Some(previous);
        }
    }

    /// Remove the residual divergence of the current iterate by one tight
    /// pressure correction (1e-8 of its imbalance).
    pub(super) fn project(&mut self, gate: &CancelGate) -> Result<(), ChtError> {
        self.correct(1e-8, gate)?;
        self.guard_finite()
    }

    /// Cell-centred velocity of `cell` (zero in solids).
    pub(super) fn cell_velocity(&self, cell: usize) -> [f64; 3] {
        if !self.domain.is_fluid(cell) {
            return [0.0; 3];
        }
        let at = self.domain.coords(cell);
        [0, 1, 2].map(|a| {
            0.5 * (self.vel[a][self.cell_face(a, at, false)]
                + self.vel[a][self.cell_face(a, at, true)])
        })
    }

    /// A diverging iteration refuses here, before non-finite coefficients
    /// reach an incomplete factorization.
    fn guard_finite(&self) -> Result<(), ChtError> {
        let finite = self.vel.iter().all(|v| v.iter().all(|u| u.is_finite()))
            && self.pressure.iter().all(|p| p.is_finite());
        if finite {
            Ok(())
        } else {
            Err(ChtError::FlowDiverged { step: self.sweeps })
        }
    }

    /// The current face fluxes.
    pub(super) fn field(&self) -> FlowField {
        let area = self.area;
        let [fx, fy, fz] = self
            .vel
            .clone()
            .map(|v| v.into_iter().map(|u| u * area).collect::<Vec<f64>>());
        FlowField::from_face_arrays(self.domain, fx, fy, fz)
    }

    fn weight(&self, peclet: f64) -> f64 {
        match self.config.scheme {
            ConvectionScheme::Upwind => 1.0,
            ConvectionScheme::PowerLaw => {
                let t = 0.1f64.mul_add(-peclet.abs(), 1.0).max(0.0);
                let t2 = t * t;
                t2 * t2 * t
            }
        }
    }

    /// Face of component `d` on side `s` of cell `c` (`s = 1`: plus side).
    fn cell_face(&self, d: usize, c: [usize; 3], plus: bool) -> usize {
        let mut f = c;
        f[d] += usize::from(plus);
        self.comps[d].index(f)
    }

    /// Momentum system of component `a`: assemble, under-relax, solve, and
    /// record SIMPLEC `d`. Returns the steady momentum residual of the
    /// velocities on entry: `||b - A u||` of the row-scaled system over the
    /// flow's momentum scale (under-relaxation cancels at `u = u_old`, so
    /// this is the residual of the unrelaxed steady equations).
    #[allow(clippy::too_many_lines)] // one staggered control-volume assembly
    fn momentum(&mut self, a: usize, gate: &CancelGate) -> Result<f64, ChtError> {
        let comp = &self.comps[a];
        let rows = comp.unknowns.len();
        if rows == 0 {
            return Ok(0.0);
        }
        let alpha = self.config.velocity_relaxation;
        let dx = self.domain.dx();
        let conductance = self.area / dx;
        let mut coo = Coo::new(rows, rows);
        let mut b = vec![0.0f64; rows];
        let mut a_p_relaxed = vec![0.0f64; rows];
        let mut neighbour_sum = vec![0.0f64; rows];
        for (row, &index) in comp.unknowns.iter().enumerate() {
            if row % 4096 == 0 {
                poll(gate)?;
            }
            let f = comp.coords(index);
            // The two cells the face separates; an outlet face has one
            // inside, and its outside twin mirrors it (zero normal gradient).
            let minus_cell = (f[a] > 0).then(|| {
                let mut c = f;
                c[a] -= 1;
                c
            });
            let plus_cell = (f[a] < self.n[a]).then_some(f);
            let (minus_cell, plus_cell) = match (minus_cell, plus_cell) {
                (Some(m), Some(p)) => (m, p),
                (Some(m), None) => (m, m),
                (None, Some(p)) => (p, p),
                (None, None) => unreachable!("every face touches a cell"),
            };
            let mut sum_nb = 0.0f64;
            let mut a_p = 0.0f64;
            let mut rhs = 0.0f64;
            let mut net_out = 0.0f64;
            for d in 0..3 {
                for plus in [false, true] {
                    let s = if plus { 1.0 } else { -1.0 };
                    if d == a {
                        // CV side at the centre of the minus/plus cell.
                        let beyond = if plus { f[a] == self.n[a] } else { f[a] == 0 };
                        if beyond {
                            // Outlet: the mirrored outside centre moves with
                            // the face itself, so only the flux remains.
                            net_out += s * self.rho * self.area * self.vel[a][index];
                            continue;
                        }
                        let nf = {
                            let mut g = f;
                            if plus {
                                g[a] += 1;
                            } else {
                                g[a] -= 1;
                            }
                            comp.index(g)
                        };
                        let flux =
                            s * self.rho * self.area * 0.5 * (self.vel[a][index] + self.vel[a][nf]);
                        net_out += flux;
                        let diff =
                            conductance * self.mu_eff(if plus { plus_cell } else { minus_cell });
                        let coef = diff.mul_add(self.weight(flux / diff), (-flux).max(0.0));
                        match comp.kind[nf] {
                            Kind::Unknown(col) | Kind::Outlet(col) => {
                                coo.push(row, col, -coef);
                                sum_nb += coef;
                            }
                            Kind::Fixed(v) => rhs = coef.mul_add(v, rhs),
                        }
                        a_p += coef;
                        continue;
                    }
                    // Transverse side at the edge shared with the next row.
                    let edge_flux = {
                        let vm = self.vel[d][self.cell_face(d, minus_cell, plus)];
                        let vp = self.vel[d][self.cell_face(d, plus_cell, plus)];
                        s * self.rho * self.area * 0.5 * (vm + vp)
                    };
                    let outside = if plus {
                        f[d] + 1 >= self.n[d]
                    } else {
                        f[d] == 0
                    };
                    let diff = conductance * self.edge_mu(minus_cell, plus_cell, d, plus);
                    if outside {
                        match self.faces[2 * d + usize::from(plus)] {
                            FvBoundary::Wall { velocity } => {
                                let wall = conductance * self.wall_mu(minus_cell, plus_cell);
                                a_p += 2.0 * wall;
                                rhs = (2.0 * wall).mul_add(velocity[a], rhs);
                            }
                            FvBoundary::Symmetry => {}
                            FvBoundary::Inlet { velocity } => {
                                let flux = s * self.rho * self.area * velocity[d];
                                net_out += flux;
                                let coef = (2.0 * diff)
                                    .mul_add(self.weight(flux / (2.0 * diff)), (-flux).max(0.0));
                                a_p += coef;
                                rhs = coef.mul_add(velocity[a], rhs);
                            }
                            FvBoundary::Outlet => {
                                net_out += edge_flux;
                                if edge_flux < 0.0 {
                                    // Backflow carries zero tangential momentum.
                                    a_p += -edge_flux;
                                }
                            }
                        }
                        continue;
                    }
                    let mut g = f;
                    if plus {
                        g[d] += 1;
                    } else {
                        g[d] -= 1;
                    }
                    let nf = comp.index(g);
                    match comp.kind[nf] {
                        Kind::Unknown(col) | Kind::Outlet(col) => {
                            net_out += edge_flux;
                            let coef =
                                diff.mul_add(self.weight(edge_flux / diff), (-edge_flux).max(0.0));
                            coo.push(row, col, -coef);
                            sum_nb += coef;
                            a_p += coef;
                        }
                        // A blocked transverse neighbour touches a solid:
                        // no-slip wall half a cell away; any half-face
                        // inflow through the edge brings zero momentum.
                        _ => {
                            net_out += edge_flux;
                            let wall = conductance * self.wall_mu(minus_cell, plus_cell);
                            a_p += 2.0f64.mul_add(wall, (-edge_flux).max(0.0));
                        }
                    }
                }
            }
            // Continuity is satisfied only at convergence: keep a_P >= sum.
            a_p += net_out.max(0.0);
            let u_old = self.vel[a][index];
            if let Some(&loss) = self.planar_loss[a].get(index) {
                // 1/2 rho K |u| u across the face (Picard in |u|).
                a_p += 0.5 * self.rho * loss * u_old.abs() * self.area;
            }
            if !self.porous.is_empty() {
                // Half of the staggered volume lies in each adjacent cell.
                let (m, p) = (
                    self.porous[cell_of(self.n, minus_cell)][a],
                    self.porous[cell_of(self.n, plus_cell)][a],
                );
                let per_volume = (0.5 * (m.1 + p.1)).mul_add(u_old.abs(), 0.5 * (m.0 + p.0));
                a_p += per_volume * self.area * dx;
            }
            if let Some(time) = &self.time {
                // rho V du/dt over the cell-sized staggered volume.
                let c = self.rho * self.area * dx / time.dt;
                match (time.scheme, &time.older) {
                    (TimeScheme::Bdf2, Some(older)) => {
                        a_p += 1.5 * c;
                        rhs = c.mul_add(
                            (-0.5f64).mul_add(older[a][index], 2.0 * time.previous[a][index]),
                            rhs,
                        );
                    }
                    _ => {
                        a_p += c;
                        rhs = c.mul_add(time.previous[a][index], rhs);
                    }
                }
            }
            if let Some(&i) = self.fan_face[a].get(index)
                && i != usize::MAX
            {
                // Rise dp(Q) on the face, linearized about the sweep's flow
                // as if the whole patch moved with this face:
                // S(u) = s A dp(Q0) + A slope A_fan (u - u0); exact at the
                // fixed point, and the stiffness damps the fan-system loop.
                let fan = &self.internal_fans[i];
                let stiffness = -fan.slope * fan.area * self.area;
                a_p += stiffness;
                rhs += (fan.sign * self.area).mul_add(fan.rise, stiffness * u_old);
            }
            // Interior faces separate two fluid cells; an outlet face sees
            // the ghost pressure -p_inside (boundary pressure zero).
            let drop = match comp.kind[index] {
                Kind::Outlet(_) if f[a] == 0 => -2.0 * self.pressure[cell_of(self.n, plus_cell)],
                Kind::Outlet(_) => 2.0 * self.pressure[cell_of(self.n, minus_cell)],
                _ => {
                    self.pressure[cell_of(self.n, minus_cell)]
                        - self.pressure[cell_of(self.n, plus_cell)]
                }
            };
            rhs = self.area.mul_add(drop, rhs);
            if !self.force.is_empty() {
                // Cell-sized staggered volume: half in each adjacent cell
                // (an outlet face's mirrored twin repeats the inside cell).
                let mean = 0.5
                    * (self.force[cell_of(self.n, minus_cell)][a]
                        + self.force[cell_of(self.n, plus_cell)][a]);
                rhs = (self.area * dx).mul_add(mean, rhs);
            }
            let relaxed = a_p / alpha;
            rhs = ((1.0 - alpha) * relaxed).mul_add(self.vel[a][index], rhs);
            coo.push(row, row, relaxed);
            b[row] = rhs;
            a_p_relaxed[row] = relaxed;
            neighbour_sum[row] = sum_nb;
        }
        let matrix = scale_rows(&coo, &mut b);
        let mut x: Vec<f64> = comp.unknowns.iter().map(|&i| self.vel[a][i]).collect();
        let mut r = vec![0.0f64; rows];
        matrix.spmv(&x, &mut r);
        // Relative to the flow's momentum scale: the largest RMS right-hand
        // side over the components (a component that vanishes by symmetry
        // has a round-off `b`, and its own ratio would be noise over noise).
        let rows_root = fs_math::det::sqrt(rows as f64);
        self.momentum_rms[a] = norm(&b) / rows_root;
        let scale = rows_root * self.momentum_rms.iter().fold(0.0f64, |m, v| m.max(*v));
        let residual = fs_math::det::sqrt(
            r.iter()
                .zip(&b)
                .map(|(ax, bi)| (bi - ax) * (bi - ax))
                .sum::<f64>(),
        );
        let steady = if scale > 0.0 {
            residual / scale
        } else {
            residual
        };
        // Inner solves reduce the entry residual by `momentum_tolerance`.
        let inner = (self.config.momentum_tolerance * steady).max(1e-15);
        if steady > 1e-15 {
            let outcome = bicgstab_ilu0("momentum", &matrix, &b, &mut x, inner, 20_000, gate)?;
            self.momentum_krylov += outcome.iterations;
        }
        for (row, &index) in self.comps[a].unknowns.iter().enumerate() {
            self.vel[a][index] = x[row];
            self.d[a][index] =
                self.area / (a_p_relaxed[row] - neighbour_sum[row]).max(f64::MIN_POSITIVE);
        }
        Ok(steady)
    }

    /// Outward mass flux of cell `c` through local face `side` (Face3 order).
    fn outward_mass(&self, c: [usize; 3], side: usize) -> f64 {
        let axis = side / 2;
        let plus = side % 2 == 1;
        let v = self.vel[axis][self.cell_face(axis, c, plus)];
        (if plus { v } else { -v }) * self.rho * self.area
    }

    /// Pressure correction solved to `tolerance` (relative to the imbalance
    /// it corrects); returns the largest cell mass imbalance before
    /// correction.
    fn correct(&mut self, tolerance: f64, gate: &CancelGate) -> Result<f64, ChtError> {
        let cells: Vec<usize> = (0..self.domain.cell_count())
            .filter(|&c| self.domain.is_fluid(c))
            .collect();
        let mut row_of = vec![usize::MAX; self.domain.cell_count()];
        for (row, &c) in cells.iter().enumerate() {
            row_of[c] = row;
        }
        let rows = cells.len();
        let mut coo = Coo::new(rows, rows);
        let mut b = vec![0.0f64; rows];
        let mut imbalance = 0.0f64;
        let mut drains = vec![false; rows];
        for (row, &c) in cells.iter().enumerate() {
            if row % 4096 == 0 {
                poll(gate)?;
            }
            let at = self.domain.coords(c);
            let mut diag = 0.0;
            let mut net = 0.0;
            for side in 0..6 {
                net += self.outward_mass(at, side);
                let axis = side / 2;
                let plus = side % 2 == 1;
                let face = self.cell_face(axis, at, plus);
                match self.comps[axis].kind[face] {
                    Kind::Unknown(_) => {
                        let coef = self.rho * self.area * self.d[axis][face];
                        diag += coef;
                        let mut nc = at;
                        if plus {
                            nc[axis] += 1;
                        } else {
                            nc[axis] -= 1;
                        }
                        coo.push(row, row_of[cell_of(self.n, nc)], -coef);
                    }
                    Kind::Outlet(_) => {
                        // u' = 2 d p'_inside against the ghost pressure.
                        diag += 2.0 * self.rho * self.area * self.d[axis][face];
                        drains[row] = true;
                    }
                    Kind::Fixed(_) => {}
                }
            }
            imbalance = imbalance.max(net.abs());
            coo.push(row, row, diag);
            b[row] = -net;
        }
        // Components that reach no outlet are singular: pin one cell each.
        let mut component = vec![usize::MAX; rows];
        let mut pins = Vec::new();
        for start in 0..rows {
            if component[start] != usize::MAX {
                continue;
            }
            let mut stack = vec![start];
            component[start] = start;
            let mut drained = false;
            while let Some(row) = stack.pop() {
                drained |= drains[row];
                let at = self.domain.coords(cells[row]);
                for side in 0..6 {
                    let axis = side / 2;
                    let plus = side % 2 == 1;
                    if matches!(
                        self.comps[axis].kind[self.cell_face(axis, at, plus)],
                        Kind::Unknown(_)
                    ) {
                        let mut nc = at;
                        if plus {
                            nc[axis] += 1;
                        } else {
                            nc[axis] -= 1;
                        }
                        let r = row_of[cell_of(self.n, nc)];
                        if component[r] == usize::MAX {
                            component[r] = start;
                            stack.push(r);
                        }
                    }
                }
            }
            if !drained {
                pins.push(start);
            }
        }
        for &row in &pins {
            coo.push(row, row, self.rho * self.area * self.domain.dx() / self.mu);
        }
        let mut correction = vec![0.0f64; rows];
        if b.iter().any(|v| *v != 0.0) {
            let mut solved = false;
            if self.config.pressure_solver == PressureSolver::AmgCg {
                poll(gate)?;
                let matrix = coo.assemble();
                let stale = self
                    .amg
                    .as_ref()
                    .is_none_or(|(_, built)| self.sweeps >= built + AMG_REBUILD_SWEEPS);
                if stale {
                    self.amg = Some((SaAmg::new(&matrix, 0.08, 3), self.sweeps));
                }
                let (amg, _) = self.amg.as_ref().expect("built above");
                let report = pcg(&matrix, &b, &mut correction, amg, tolerance, 500);
                self.pressure_krylov += report.iters;
                solved = report.converged && correction.iter().all(|v| v.is_finite());
                if !solved {
                    correction.fill(0.0);
                }
            }
            if !solved {
                let matrix = scale_rows(&coo, &mut b);
                let outcome = bicgstab_ilu0(
                    "pressure",
                    &matrix,
                    &b,
                    &mut correction,
                    tolerance,
                    20_000,
                    gate,
                )?;
                self.pressure_krylov += outcome.iterations;
            }
        }
        let p_of = |cell: [usize; 3]| -> f64 {
            let c = cell_of(self.n, cell);
            if self.domain.is_fluid(c) {
                correction[row_of[c]]
            } else {
                0.0
            }
        };
        for a in 0..3 {
            for index in 0..self.comps[a].kind.len() {
                let f = self.comps[a].coords(index);
                match self.comps[a].kind[index] {
                    Kind::Unknown(_) => {
                        let mut minus = f;
                        minus[a] -= 1;
                        self.vel[a][index] += self.d[a][index] * (p_of(minus) - p_of(f));
                    }
                    Kind::Outlet(_) => {
                        if f[a] == 0 {
                            self.vel[a][index] -= 2.0 * self.d[a][index] * p_of(f);
                        } else {
                            let mut minus = f;
                            minus[a] -= 1;
                            self.vel[a][index] += 2.0 * self.d[a][index] * p_of(minus);
                        }
                    }
                    Kind::Fixed(_) => {}
                }
            }
        }
        for (row, &c) in cells.iter().enumerate() {
            self.pressure[c] += correction[row];
        }
        Ok(imbalance)
    }
}

/// Admission shared by every SIMPLEC driver.
#[allow(clippy::too_many_lines)] // one check per declared input
pub(super) fn admit(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
) -> Result<(), ChtError> {
    fluid.validate()?;
    // SIMPLEC's d = A / (a_P / alpha - sum a_nb) degenerates at alpha = 1.
    if !(config.velocity_relaxation > 0.0 && config.velocity_relaxation < 1.0) {
        return Err(ChtError::InvalidInput {
            field: "simple.velocity_relaxation",
            reason: format!("must lie in (0, 1), got {}", config.velocity_relaxation),
        });
    }
    finite_positive("simple.tolerance", config.tolerance)?;
    finite_positive("simple.momentum_tolerance", config.momentum_tolerance)?;
    finite_positive("simple.pressure_tolerance", config.pressure_tolerance)?;
    for (side, rule) in config.faces.iter().enumerate() {
        let axis = side / 2;
        let inward = if side % 2 == 0 { 1.0 } else { -1.0 };
        match *rule {
            FvBoundary::Wall { velocity } => {
                for v in velocity {
                    finite("simple.wall_velocity", v)?;
                }
                if velocity[axis] != 0.0 {
                    return Err(ChtError::InvalidInput {
                        field: "simple.wall_velocity",
                        reason: format!("wall on face {side} has a normal velocity"),
                    });
                }
            }
            FvBoundary::Inlet { velocity } => {
                for v in velocity {
                    finite("simple.inlet_velocity", v)?;
                }
                if velocity[axis] * inward <= 0.0 {
                    return Err(ChtError::InvalidInput {
                        field: "simple.inlet_velocity",
                        reason: format!("inlet on face {side} must point into the domain"),
                    });
                }
            }
            FvBoundary::Symmetry | FvBoundary::Outlet => {}
        }
    }
    if let Some(fan) = config.fan
        && !matches!(config.faces[fan.face as usize], FvBoundary::Inlet { .. })
    {
        return Err(ChtError::InvalidInput {
            field: "simple.fan",
            reason: format!("the fan face {:?} must be declared an Inlet", fan.face),
        });
    }
    for resistance in &config.resistances {
        match *resistance {
            FlowResistance::Planar {
                patch,
                loss_coefficient,
            } => {
                patch.admit(domain, "simple.resistance.patch")?;
                finite("simple.resistance.loss_coefficient", loss_coefficient)?;
                if loss_coefficient < 0.0 {
                    return Err(ChtError::InvalidInput {
                        field: "simple.resistance.loss_coefficient",
                        reason: format!("must be non-negative, got {loss_coefficient}"),
                    });
                }
            }
            FlowResistance::Volume {
                lo,
                hi,
                permeability_m2,
                inertial_per_m,
            } => {
                let n = domain.dims();
                if (0..3).any(|a| lo[a] >= hi[a] || hi[a] > n[a]) {
                    return Err(ChtError::InvalidInput {
                        field: "simple.resistance.cells",
                        reason: format!("cells {lo:?}..{hi:?} must be non-empty within {n:?}"),
                    });
                }
                for a in 0..3 {
                    if !(permeability_m2[a] > 0.0) {
                        return Err(ChtError::InvalidInput {
                            field: "simple.resistance.permeability_m2",
                            reason: format!(
                                "must be positive (infinite: no viscous term), got {}",
                                permeability_m2[a]
                            ),
                        });
                    }
                    finite("simple.resistance.inertial_per_m", inertial_per_m[a])?;
                    if inertial_per_m[a] < 0.0 {
                        return Err(ChtError::InvalidInput {
                            field: "simple.resistance.inertial_per_m",
                            reason: format!("must be non-negative, got {}", inertial_per_m[a]),
                        });
                    }
                }
            }
        }
    }
    for fan in &config.internal_fans {
        fan.patch.admit(domain, "simple.internal_fan.patch")?;
        if fan.patch.open_faces(domain).is_empty() {
            return Err(ChtError::InvalidInput {
                field: "simple.internal_fan.patch",
                reason: "the fan covers no face between two fluid cells".into(),
            });
        }
    }
    if domain.fluid_count() == 0 {
        return Err(ChtError::InvalidDomain {
            reason: "no fluid cell".into(),
        });
    }
    Ok(())
}

impl Solver<'_> {
    /// Project tightly and assemble the public result.
    pub(super) fn finish(
        mut self,
        fluid: &FluidProperties,
        iterations: usize,
        (mass_residual, momentum_residual): (f64, f64),
        gate: &CancelGate,
    ) -> Result<FvFlow, ChtError> {
        // Hand over an exactly projected field: one tight correction removes
        // the residual divergence the loose inner solves leave (its size is
        // bounded by the converged imbalance, so the flow is unchanged to
        // tolerance).
        self.correct(1e-8, gate)?;
        let domain = self.domain;
        let field = self.field();
        let mut velocity_m_s = vec![[0.0f64; 3]; domain.cell_count()];
        let mut max_cell_reynolds = 0.0f64;
        for (c, slot) in velocity_m_s.iter_mut().enumerate() {
            if !domain.is_fluid(c) {
                self.pressure[c] = 0.0;
                continue;
            }
            let at = domain.coords(c);
            for a in 0..3 {
                slot[a] = 0.5
                    * (self.vel[a][self.cell_face(a, at, false)]
                        + self.vel[a][self.cell_face(a, at, true)]);
            }
            let speed = fs_math::det::sqrt(slot.iter().map(|v| v * v).sum::<f64>());
            max_cell_reynolds =
                max_cell_reynolds.max(speed * domain.dx() / fluid.kinematic_viscosity_m2_s);
        }
        let (mut inflow, mut outflow) = (0.0f64, 0.0f64);
        for face in Face3::ALL {
            let net = field.boundary_outflow(domain, face);
            match self.faces[face as usize] {
                FvBoundary::Inlet { .. } => inflow -= net,
                // Openings pass flow both ways: count each direction.
                FvBoundary::Outlet => {
                    for c in 0..domain.cell_count() {
                        if domain.neighbor(c, face as usize).is_none() && domain.is_fluid(c) {
                            let out = field.outward(domain, c, face as usize);
                            if out > 0.0 {
                                outflow += out;
                            } else {
                                inflow -= out;
                            }
                        }
                    }
                }
                _ => {}
            }
        }
        let report = SimpleReport {
            iterations,
            mass_residual,
            momentum_residual,
            max_divergence_m3_s: field.max_divergence(domain),
            inflow_m3_s: inflow,
            outflow_m3_s: outflow,
            max_cell_reynolds,
            momentum_krylov_iterations: self.momentum_krylov,
            pressure_krylov_iterations: self.pressure_krylov,
            fan: self
                .fan
                .as_ref()
                .map(|fan| (fan.flow, fan.curve.pressure(fan.flow), fan.residual)),
            internal_fans: self
                .internal_fans
                .iter()
                .map(|fan| {
                    let sum: f64 = fan.faces.iter().map(|&f| self.vel[fan.axis][f]).sum();
                    let flow = fan.sign * self.area * sum;
                    (flow, fan.curve.pressure(flow))
                })
                .collect(),
        };
        let eddy_viscosity_m2_s = self.eddy_viscosity();
        let heat_eddy_viscosity_m2_s = self.heat_eddy_viscosity();
        Ok(FvFlow {
            field,
            velocity_m_s,
            pressure_pa: self.pressure,
            eddy_viscosity_m2_s,
            heat_eddy_viscosity_m2_s,
            report,
        })
    }
}

/// Solve steady incompressible flow on `domain` by SIMPLEC.
///
/// # Errors
/// Input refusals (non-finite rules, an inlet pointing outward, a wall with
/// normal velocity), [`ChtError::FlowNotSteady`] when the iteration budget
/// ends first, solver refusals, or [`ChtError::Cancelled`].
pub fn simple_flow(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    gate: &CancelGate,
) -> Result<FvFlow, ChtError> {
    admit(domain, fluid, config)?;
    let mut solver = Solver::new(domain, fluid, config);
    let mut iterations = 0usize;
    let mut residuals = (f64::INFINITY, f64::INFINITY);
    while iterations < config.max_iterations {
        poll(gate)?;
        iterations += 1;
        residuals = solver.sweep(gate)?;
        if !(residuals.0.is_finite() && residuals.1.is_finite()) {
            return Err(ChtError::FlowDiverged { step: iterations });
        }
        if residuals.0 <= config.tolerance && residuals.1 <= config.tolerance {
            break;
        }
    }
    if residuals.0 > config.tolerance || residuals.1 > config.tolerance {
        return Err(ChtError::FlowNotSteady {
            steps: iterations,
            last_change: residuals.0.max(residuals.1),
            tolerance: config.tolerance,
        });
    }
    solver.finish(fluid, iterations, residuals, gate)
}
