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
//!
//! With spatially varying molecular or eddy viscosity, momentum uses the
//! complete symmetric stress `div(mu_eff (grad(u) + grad(u)^T))`. The isotropic
//! Reynolds stress is absorbed into the incompressible pressure; LVEL does not
//! supply a separately resolved turbulent kinetic energy.

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

pub mod sgs;
pub mod viscous_heat;
#[cfg(feature = "enthalpy-coupling")]
pub(in crate::conjugate) mod phase_capture;

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
    /// LVEL uses its wall-law secant value in wall-adjacent cells. Explicit
    /// transient SGS uses its resolved-gradient value (zeros when laminar).
    pub heat_eddy_viscosity_m2_s: Vec<f64>,
    /// Run evidence.
    pub report: SimpleReport,
}

impl FvFlow {
    /// Exact empty/rest field, without inventing a fluid transport tuple.
    #[cfg(feature = "enthalpy-coupling")]
    pub(in crate::conjugate) fn stationary(domain: &VoxelDomain) -> Self {
        Self {
            field: FlowField::quiescent(domain),
            velocity_m_s: vec![[0.0; 3]; domain.cell_count()],
            pressure_pa: vec![0.0; domain.cell_count()],
            eddy_viscosity_m2_s: vec![0.0; domain.cell_count()],
            heat_eddy_viscosity_m2_s: vec![0.0; domain.cell_count()],
            report: SimpleReport {
                iterations: 0,
                mass_residual: 0.0,
                momentum_residual: 0.0,
                max_divergence_m3_s: 0.0,
                inflow_m3_s: 0.0,
                outflow_m3_s: 0.0,
                max_cell_reynolds: 0.0,
                momentum_krylov_iterations: 0,
                pressure_krylov_iterations: 0,
                fan: None,
                internal_fans: Vec::new(),
            },
        }
    }

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
        // Area-weighted over the layer's fluid cells (graded grids).
        let (mut sum, mut area) = (0.0f64, 0.0f64);
        for c in 0..domain.cell_count() {
            if domain.is_fluid(c) && domain.coords(c)[axis] == index {
                let a = domain.face_area(c, axis);
                sum = a.mul_add(self.pressure_pa[c], sum);
                area += a;
            }
        }
        (area > 0.0).then(|| sum / area)
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
    /// Admitted laminar molecular viscosity per cell. Either a molecular field
    /// or an eddy field selects the complete symmetric stress.
    molecular_viscosity: Vec<f64>,
    /// Explicit Stokes reduction; all ordinary SIMPLE entries keep convection.
    creeping: bool,
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
    /// Explicit transient SGS closure and reusable gradient workspace.
    sgs: Option<sgs::State>,
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

/// Exact physical state of the native staggered momentum owner. Iteration
/// caches are rebuilt when source transport changes; face momentum is retained.
#[cfg(feature = "enthalpy-coupling")]
#[derive(Debug, Clone)]
pub(super) struct MomentumState {
    velocity: [Vec<f64>; 3],
    pressure: Vec<f64>,
}

#[cfg(feature = "enthalpy-coupling")]
impl MomentumState {
    /// Closed phase domains can change occupancy without a kinetic handoff
    /// only when every newly constrained native velocity is already zero.
    pub(super) fn admit_phase_domain(
        &self,
        domain: &VoxelDomain,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        let n = domain.dims();
        for (axis, velocity) in self.velocity.iter().enumerate() {
            let mut dims = n;
            dims[axis] += 1;
            if velocity.len() != dims.iter().product::<usize>()
                || self.pressure.len() != domain.cell_count()
            {
                return Err(ChtError::InvalidInput {
                    field: "simple.momentum_state",
                    reason: "phase momentum must retain its original staggered grid".into(),
                });
            }
            for (index, &speed) in velocity.iter().enumerate() {
                if index % 256 == 0 {
                    poll(gate)?;
                }
                let plus = [
                    index % dims[0],
                    (index / dims[0]) % dims[1],
                    index / (dims[0] * dims[1]),
                ];
                let along = plus[axis];
                let constrained = if along == 0 || along == n[axis] {
                    true
                } else {
                    let mut minus = plus;
                    minus[axis] -= 1;
                    !domain.is_fluid(cell_of(n, minus)) || !domain.is_fluid(cell_of(n, plus))
                };
                if constrained && speed != 0.0 {
                    return Err(ChtError::InvalidInput {
                        field: "simple.momentum_state",
                        reason: "freezing a moving native face requires an explicit momentum/kinetic-energy handoff".into(),
                    });
                }
            }
        }
        Ok(())
    }

    pub(super) fn clear_solid_pressure(&mut self, domain: &VoxelDomain) {
        for (cell, pressure) in self.pressure.iter_mut().enumerate() {
            if !domain.is_fluid(cell) {
                *pressure = 0.0;
            }
        }
    }
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
    /// Area of each face, m^2.
    face_areas: Vec<f64>,
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
                let open = fan.patch.open_faces(domain);
                let face_areas: Vec<f64> = open
                    .iter()
                    .map(|f| {
                        let (u, w) = ((fan.patch.axis + 1) % 3, (fan.patch.axis + 2) % 3);
                        domain.width(u, f[u]) * domain.width(w, f[w])
                    })
                    .collect();
                let faces: Vec<usize> = open.into_iter().map(|f| comp.index(f)).collect();
                InternalFanState {
                    axis: fan.patch.axis,
                    sign: if fan.blows_positive { 1.0 } else { -1.0 },
                    area: face_areas.iter().sum(),
                    face_areas,
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
            molecular_viscosity: Vec::new(),
            creeping: false,
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
            sgs: None,
            fan: config.fan.map(|fan| {
                let side = fan.face as usize;
                let axis = side / 2;
                let cells: Vec<usize> = (0..domain.cell_count())
                    .filter(|&c| domain.is_fluid(c) && domain.neighbor(c, side).is_none())
                    .collect();
                let area = cells.iter().map(|&c| domain.face_area(c, axis)).sum();
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
        if !self.molecular_viscosity.is_empty() {
            return self.molecular_viscosity[cell_of(self.n, cell)];
        }
        if self.eddy_viscosity.is_empty() {
            self.mu
        } else {
            self.rho
                .mul_add(self.eddy_viscosity[cell_of(self.n, cell)], self.mu)
        }
    }

    /// `div(mu grad(u)^T)` vanishes only for constant viscosity and a
    /// divergence-free velocity. Eddy viscosity is a spatial field too.
    fn symmetric_stress(&self) -> bool {
        !self.molecular_viscosity.is_empty() || !self.eddy_viscosity.is_empty()
    }

    /// Wall viscosity of the half cell between `minus`/`plus` and a wall:
    /// the law-of-the-wall secant (molecular when laminar).
    fn wall_mu(&self, minus: [usize; 3], plus: [usize; 3]) -> f64 {
        if !self.molecular_viscosity.is_empty() {
            let Some(a) = (0..3).find(|&a| minus[a] != plus[a]) else {
                return self.mu_eff(minus);
            };
            let wm = self.domain.width(a, minus[a]);
            let wp = self.domain.width(a, plus[a]);
            // Parallel half-cell wall patches, weighted by their actual areas.
            return (wm / (wm + wp)) * self.mu_eff(minus) + (wp / (wm + wp)) * self.mu_eff(plus);
        }
        if self.wall_ratio.is_empty() {
            return self.mu;
        }
        let rm = self.wall_ratio[cell_of(self.n, minus)];
        let rp = self.wall_ratio[cell_of(self.n, plus)];
        let ratio = if let Some(a) = (0..3).find(|&a| minus[a] != plus[a]) {
            let wm = self.domain.width(a, minus[a]);
            let wp = self.domain.width(a, plus[a]);
            (wm / (wm + wp)) * rm + (wp / (wm + wp)) * rp
        } else {
            rm
        };
        self.mu * ratio
    }

    /// Original wall conductance over one full row width. The no-slip
    /// half-cell boundary coefficient is TWICE this value. Momentum and heat
    /// share both its source viscosity and geometry; the original arithmetic
    /// order and momentum FMA remain unchanged. `minus` and `plus` border an
    /// active component-`a` face; `d != a` is the wall-normal direction.
    fn wall_shear_base_conductance(
        &self,
        minus: [usize; 3],
        plus: [usize; 3],
        a: usize,
        d: usize,
    ) -> f64 {
        let length = 0.5 * (self.domain.width(a, minus[a]) + self.domain.width(a, plus[a]));
        let side_area = length * self.domain.width(3 - a - d, plus[3 - a - d]);
        (side_area / self.domain.width(d, plus[d])) * self.wall_mu(minus, plus)
    }

    /// Mean effective viscosity at the CV edge shared by `minus` and
    /// `plus` across direction `d` (`up`: toward +d): the fluid cells among
    /// the two and their neighbours across the edge.
    fn edge_mu(&self, minus: [usize; 3], plus: [usize; 3], d: usize, up: bool) -> f64 {
        if self.symmetric_stress() {
            // One positive coefficient at the MAC shear edge, shared by
            // tau_ad and tau_da. Material jumps use harmonic transmission;
            // smooth eddy fields retain the native arithmetic interpolation.
            // The same area weights also handle graded cells.
            let a = (0..3).find(|&a| minus[a] != plus[a]).unwrap_or((d + 1) % 3);
            let mut weights = [(0usize, 0.0, 1.0); 4];
            let mut count = 0;
            for c in [minus, plus] {
                weights[count] = (
                    cell_of(self.n, c),
                    self.domain.width(a, c[a]) * self.domain.width(d, c[d]),
                    self.mu_eff(c),
                );
                count += 1;
                let mut g = c;
                let exists = if up { c[d] + 1 < self.n[d] } else { c[d] > 0 };
                if exists {
                    g[d] = if up { c[d] + 1 } else { c[d] - 1 };
                    if self.domain.is_fluid(cell_of(self.n, g)) {
                        weights[count] = (
                            cell_of(self.n, g),
                            self.domain.width(a, g[a]) * self.domain.width(d, g[d]),
                            self.mu_eff(g),
                        );
                        count += 1;
                    }
                }
            }
            // Opposite CV sides and reciprocal components must visit the
            // same four source cells in the same order, including roundoff.
            weights[..count].sort_unstable_by_key(|(cell, _, _)| *cell);
            let area: f64 = weights[..count].iter().map(|(_, w, _)| w).sum();
            if self.molecular_viscosity.is_empty() {
                return weights[..count]
                    .iter()
                    .map(|(_, w, mu)| (w / area) * mu)
                    .sum();
            }
            let minimum = weights[..count]
                .iter()
                .map(|(_, _, mu)| *mu)
                .fold(f64::INFINITY, f64::min);
            let inverse: f64 = weights[..count]
                .iter()
                .map(|(_, w, mu)| (w / area) * (minimum / mu))
                .sum();
            return minimum / inverse;
        }
        self.mu
    }

    /// LVEL: refresh the eddy viscosity from the current cell-centred
    /// speeds (under-relaxed by one half after the first sweep).
    fn update_turbulence(&mut self) {
        if self.config.turbulence != Turbulence::Lvel {
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
            let reynolds = fs_math::det::sqrt(speed2) * self.wall_distance[c] / nu;
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
    /// law-of-the-wall secant `nu (y+/u+ - 1)` in LVEL wall-adjacent cells
    /// (their conduction spans the half cell), LVEL's value elsewhere. SGS
    /// uses its own resolved-gradient viscosity without a wall-law overlay.
    pub(super) fn heat_eddy_viscosity(&self) -> Vec<f64> {
        if self.eddy_viscosity.is_empty() {
            return vec![0.0; self.domain.cell_count()];
        }
        if self.sgs.is_some() {
            return self.eddy_viscosity.clone();
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
        // Area-weighted mean pressure of the layer behind the fan face.
        let axis = fan.side / 2;
        let mean = fan
            .cells
            .iter()
            .map(|&c| self.pressure[c] * self.domain.face_area(c, axis))
            .sum::<f64>()
            / fan.area.max(f64::MIN_POSITIVE);
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
        self.update_sgs(gate)?;
        for fan in &mut self.internal_fans {
            let sum: f64 = fan
                .faces
                .iter()
                .zip(&fan.face_areas)
                .map(|(&f, area)| self.vel[fan.axis][f] * area)
                .sum();
            fan.flow = fan.sign * sum;
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
        let largest_flux = (0..3)
            .flat_map(|a| (0..self.vel[a].len()).map(move |i| (a, i)))
            .fold(0.0f64, |m, (a, i)| {
                m.max((self.vel[a][i] * self.comp_face_area(a, i)).abs())
            })
            .max(f64::MIN_POSITIVE);
        let mass = imbalance / (self.rho * largest_flux);
        self.update_fan(mass, steady);
        let fan = self.fan.as_ref().map_or(0.0, |fan| fan.residual);
        Ok((mass, steady.max(fan)))
    }

    /// Checkpoint native faces directly, without a cell-centred reconstruction.
    #[cfg(feature = "enthalpy-coupling")]
    pub(super) fn momentum_state(&self) -> MomentumState {
        MomentumState {
            velocity: self.vel.clone(),
            pressure: self.pressure.clone(),
        }
    }

    /// Restore a physical state into this same domain's material operator.
    #[cfg(feature = "enthalpy-coupling")]
    pub(super) fn restore_momentum(&mut self, state: &MomentumState) -> Result<(), ChtError> {
        if self.pressure.len() != state.pressure.len()
            || self
                .vel
                .iter()
                .zip(&state.velocity)
                .any(|(a, b)| a.len() != b.len())
        {
            return Err(ChtError::InvalidInput {
                field: "simple.momentum_state",
                reason: "retained momentum must match the native staggered domain".into(),
            });
        }
        // Opening a formerly constrained face retains its exact zero speed.
        // Closing one is admissible only if it already satisfies the new
        // constraint. Never clip a live velocity and silently delete its
        // momentum or native face-mass kinetic energy.
        for (component, velocity) in self.comps.iter().zip(&state.velocity) {
            for (kind, &speed) in component.kind.iter().zip(velocity) {
                if let Kind::Fixed(required) = kind {
                    if speed != *required {
                        return Err(ChtError::InvalidInput {
                            field: "simple.momentum_state",
                            reason: "a new fixed face requires an explicit momentum/kinetic-energy handoff".into(),
                        });
                    }
                }
            }
        }
        self.vel.clone_from(&state.velocity);
        self.pressure.clone_from(&state.pressure);
        self.guard_finite()
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
        let [fx, fy, fz] = [0, 1, 2].map(|a| {
            self.vel[a]
                .iter()
                .enumerate()
                .map(|(i, u)| u * self.comp_face_area(a, i))
                .collect::<Vec<f64>>()
        });
        FlowField::from_face_arrays(self.domain, fx, fy, fz)
    }

    /// Area of face `index` of component `a`, m^2.
    fn comp_face_area(&self, a: usize, index: usize) -> f64 {
        let f = self.comps[a].coords(index);
        let (u, w) = ((a + 1) % 3, (a + 2) % 3);
        self.domain.width(u, f[u]) * self.domain.width(w, f[w])
    }

    /// Native inertia of an interior staggered face, shared by kinetic storage
    /// and the explicit stationary-support capture law. Both adjacent cells
    /// must belong to the original mobile domain.
    #[cfg(feature = "enthalpy-coupling")]
    fn interior_face_mass_kg(&self, axis: usize, index: usize) -> f64 {
        let plus = self.comps[axis].coords(index);
        self.rho
            * self.comp_face_area(axis, index)
            * (0.5
                * (self.domain.width(axis, plus[axis] - 1)
                    + self.domain.width(axis, plus[axis])))
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
        let advective_density = if self.creeping { 0.0 } else { self.rho };
        let width = |axis: usize, i: usize| self.domain.width(axis, i);
        let mut coo = Coo::new(rows, rows);
        let mut b = vec![0.0f64; rows];
        let mut a_p_relaxed = vec![0.0f64; rows];
        let mut neighbour_sum = vec![0.0f64; rows];
        let mut areas = vec![0.0f64; rows];
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
            // Staggered control volume: the face's area across, from the
            // minus to the plus cell centre along (an outlet's mirrored twin
            // repeats the inside cell's width).
            let face_area = width((a + 1) % 3, f[(a + 1) % 3]) * width((a + 2) % 3, f[(a + 2) % 3]);
            let (w_minus, w_plus) = (width(a, minus_cell[a]), width(a, plus_cell[a]));
            let length = 0.5 * (w_minus + w_plus);
            areas[row] = face_area;
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
                            net_out += s * advective_density * face_area * self.vel[a][index];
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
                        let flux = s
                            * advective_density
                            * face_area
                            * 0.5
                            * (self.vel[a][index] + self.vel[a][nf]);
                        net_out += flux;
                        // Diffusion across the cell between the two faces.
                        let side_cell = if plus { plus_cell } else { minus_cell };
                        // A variable molecular OR eddy viscosity requires the
                        // full stress: normal stress is 2 mu_eff du_a/dx_a.
                        // The legacy constant-viscosity Laplacian keeps its
                        // arithmetic; div(grad(u)^T)=0 for incompressible flow.
                        let normal = if self.symmetric_stress() { 2.0 } else { 1.0 };
                        let diff =
                            normal * self.mu_eff(side_cell) * face_area / width(a, side_cell[a]);
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
                    // Transverse side at the edge shared with the next row:
                    // half of it beside each cell, area-weighted.
                    let third = 3 - a - d;
                    let w3 = width(third, f[third]);
                    let side_area = length * w3;
                    let edge_flux = {
                        let vm = self.vel[d][self.cell_face(d, minus_cell, plus)];
                        let vp = self.vel[d][self.cell_face(d, plus_cell, plus)];
                        s * advective_density * w3 * (0.5 * w_minus).mul_add(vm, 0.5 * w_plus * vp)
                    };
                    let outside = if plus {
                        f[d] + 1 >= self.n[d]
                    } else {
                        f[d] == 0
                    };
                    // Node-to-node distance across the edge (to a wall: the
                    // half row width).
                    let row_width = width(d, f[d]);
                    let distance = if outside {
                        row_width
                    } else {
                        let next = if plus { f[d] + 1 } else { f[d] - 1 };
                        0.5 * (row_width + width(d, next))
                    };
                    let edge_mu = self.edge_mu(minus_cell, plus_cell, d, plus);
                    let diff = edge_mu * side_area / distance;
                    // Reciprocal shear mu du_d/dx_a on this SAME staggered
                    // control-volume side. Keeping only mu du_a/dx_d would
                    // give a scalar diffusion model for a variable-mu fluid.
                    // This term is iterated with the original coupled velocity
                    // solve; it is not an externally supplied body force.
                    let cross_shear = if self.symmetric_stress() {
                        let vm = self.vel[d][self.cell_face(d, minus_cell, plus)];
                        let vp = self.vel[d][self.cell_face(d, plus_cell, plus)];
                        s * edge_mu * w3 * (vp - vm)
                    } else {
                        0.0
                    };
                    if outside {
                        match self.faces[2 * d + usize::from(plus)] {
                            FvBoundary::Wall { velocity } => {
                                let wall =
                                    self.wall_shear_base_conductance(minus_cell, plus_cell, a, d);
                                a_p += 2.0 * wall;
                                rhs = (2.0 * wall).mul_add(velocity[a], rhs);
                            }
                            FvBoundary::Symmetry => {}
                            FvBoundary::Inlet { velocity } => {
                                let flux = s * advective_density * side_area * velocity[d];
                                net_out += flux;
                                let coef = (2.0 * diff)
                                    .mul_add(self.weight(flux / (2.0 * diff)), (-flux).max(0.0));
                                a_p += coef;
                                rhs = coef.mul_add(velocity[a], rhs);
                            }
                            FvBoundary::Outlet => {
                                rhs += cross_shear;
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
                            rhs += cross_shear;
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
                            let wall =
                                self.wall_shear_base_conductance(minus_cell, plus_cell, a, d);
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
                a_p += 0.5 * self.rho * loss * u_old.abs() * face_area;
            }
            if !self.porous.is_empty() {
                // Half of each adjacent cell lies in the staggered volume.
                let (m, p) = (
                    self.porous[cell_of(self.n, minus_cell)][a],
                    self.porous[cell_of(self.n, plus_cell)][a],
                );
                let coef = |c: (f64, f64)| c.1.mul_add(u_old.abs(), c.0);
                a_p += face_area * (0.5 * w_minus).mul_add(coef(m), 0.5 * w_plus * coef(p));
            }
            if let Some(time) = &self.time {
                // rho V du/dt over the staggered volume.
                let c = self.rho * face_area * length / time.dt;
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
                let stiffness = -fan.slope * fan.area * face_area;
                a_p += stiffness;
                rhs += (fan.sign * face_area).mul_add(fan.rise, stiffness * u_old);
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
            rhs = face_area.mul_add(drop, rhs);
            if !self.force.is_empty() {
                // Half of each adjacent cell lies in the staggered volume
                // (an outlet face's mirrored twin repeats the inside cell).
                let (fm, fp) = (
                    self.force[cell_of(self.n, minus_cell)][a],
                    self.force[cell_of(self.n, plus_cell)][a],
                );
                rhs = face_area.mul_add((0.5 * w_minus).mul_add(fm, 0.5 * w_plus * fp), rhs);
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
        // The relaxed momentum operator is a strictly diagonally dominant
        // M-matrix for finite coefficients, so a non-finite residual or an
        // incomplete-factorization breakdown means the iteration diverged.
        if !steady.is_finite() {
            return Err(ChtError::FlowDiverged { step: self.sweeps });
        }
        if steady > 1e-15 {
            let outcome = bicgstab_ilu0("momentum", &matrix, &b, &mut x, inner, 20_000, gate)
                .map_err(|error| match error {
                    ChtError::PreconditionerBreakdown { .. } => {
                        ChtError::FlowDiverged { step: self.sweeps }
                    }
                    other => other,
                })?;
            self.momentum_krylov += outcome.iterations;
        }
        for (row, &index) in self.comps[a].unknowns.iter().enumerate() {
            self.vel[a][index] = x[row];
            self.d[a][index] =
                areas[row] / (a_p_relaxed[row] - neighbour_sum[row]).max(f64::MIN_POSITIVE);
        }
        Ok(steady)
    }

    /// Outward mass flux of cell `c` through local face `side` (Face3 order).
    fn outward_mass(&self, c: [usize; 3], side: usize) -> f64 {
        let axis = side / 2;
        let plus = side % 2 == 1;
        let v = self.vel[axis][self.cell_face(axis, c, plus)];
        (if plus { v } else { -v }) * self.rho * self.domain.face_area(cell_of(self.n, c), axis)
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
                        let coef = self.rho * self.domain.face_area(c, axis) * self.d[axis][face];
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
                        diag +=
                            2.0 * self.rho * self.domain.face_area(c, axis) * self.d[axis][face];
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
            // Any positive scale pins the component (its correction is a
            // free constant); the cell's own A dx / nu keeps it conditioned.
            let c = cells[row];
            coo.push(
                row,
                row,
                self.rho * self.domain.face_area(c, 0) * self.domain.widths(c)[0] / self.mu,
            );
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
    if config.turbulence == Turbulence::Lvel
        && domain.fluid_count() == domain.cell_count()
        && !config
            .faces
            .iter()
            .any(|face| matches!(face, FvBoundary::Wall { .. }))
    {
        return Err(ChtError::InvalidInput {
            field: "simple.turbulence",
            reason: "LVEL requires a no-slip wall or solid voxel; open and symmetry faces do not define a wall distance".into(),
        });
    }
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
        self.update_sgs(gate)?;
        Ok(self.finish_projected(fluid, iterations, (mass_residual, momentum_residual)))
    }

    /// Report an already projected field, keeping a checkpoint and its receipt
    /// on exactly the same native face/pressure state.
    pub(super) fn finish_projected(
        mut self,
        fluid: &FluidProperties,
        iterations: usize,
        (mass_residual, momentum_residual): (f64, f64),
    ) -> FvFlow {
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
            max_cell_reynolds = max_cell_reynolds.max(
                speed * domain.widths(c).iter().fold(0.0f64, |m, w| m.max(*w))
                    / if self.molecular_viscosity.is_empty() {
                        fluid.kinematic_viscosity_m2_s
                    } else {
                        self.molecular_viscosity[c] / self.rho
                    },
            );
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
                    let sum: f64 = fan
                        .faces
                        .iter()
                        .zip(&fan.face_areas)
                        .map(|(&f, area)| self.vel[fan.axis][f] * area)
                        .sum();
                    let flow = fan.sign * sum;
                    (flow, fan.curve.pressure(flow))
                })
                .collect(),
        };
        let eddy_viscosity_m2_s = self.eddy_viscosity();
        let heat_eddy_viscosity_m2_s = self.heat_eddy_viscosity();
        FvFlow {
            field,
            velocity_m_s,
            pressure_pa: self.pressure,
            eddy_viscosity_m2_s,
            heat_eddy_viscosity_m2_s,
            report,
        }
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
    solve_steady(Solver::new(domain, fluid, config), fluid, config, gate)
}

/// Solve constant-density laminar flow with local dynamic viscosity [Pa s].
///
/// The existing MAC momentum owner assembles the full symmetric Newtonian
/// stress `div(mu * (grad(u) + grad(u)^T))`. Normal and reciprocal shear terms
/// share the same local coefficients. Merely substituting a field into the
/// constant-mu Laplacian would omit a physical force when viscosity varies.
/// The original pressure correction, convection, fans, porous resistance,
/// boundaries and convergence/refusal owners remain in charge.
///
/// One positive value is required for every voxel, including unused solids.
/// Density and geometry remain fixed. No variable-density continuity, phase
/// mobility, non-Newtonian rheology, viscous heating or turbulent closure is
/// inferred. Source identity and the temperature epoch belong to the calling
/// material owner. This is a steady solve; it retains no transient momentum.
///
/// # Errors
/// Refuses a malformed/nonpositive field or LVEL, and propagates the original
/// numerical and cancellation refusals without publishing a partial field.
pub fn simple_flow_with_viscosity(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    dynamic_viscosity_pa_s: &[f64],
    config: &SimpleConfig,
    gate: &CancelGate,
) -> Result<FvFlow, ChtError> {
    material_flow(
        domain,
        fluid,
        dynamic_viscosity_pa_s,
        None,
        None,
        config,
        gate,
    )
}

/// Solve the same variable-viscosity momentum system with a cell body force.
///
/// Forces are physical force densities [N/m³], with one finite vector per
/// voxel and zero on solids. They enter the existing staggered momentum
/// volumes together with pressure, fans, porous resistance and full molecular
/// stress. No separate projection or cell-velocity interpolation is added.
/// For Boussinesq flow supply `(rho(T)-rho0)*g`, with `fluid` retaining rho0;
/// pressure is then relative to that removed hydrostatic reference. The caller
/// owns source selection, the force epoch and any mechanical-energy account.
pub fn simple_flow_with_viscosity_and_force(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    dynamic_viscosity_pa_s: &[f64],
    force_density_n_m3: &[[f64; 3]],
    config: &SimpleConfig,
    gate: &CancelGate,
) -> Result<FvFlow, ChtError> {
    material_flow(
        domain,
        fluid,
        dynamic_viscosity_pa_s,
        Some(force_density_n_m3),
        None,
        config,
        gate,
    )
}

/// Add declared local isotropic Darcy resistance to the same momentum owner.
///
/// `drag_pa_s_m2` is the coefficient of `-drag*u` [N/m3], one nonnegative
/// finite value per voxel, zero on solids. It is implicit in momentum and
/// SIMPLEC's pressure correction, and adds to authored porous resistances.
/// Zero restores the existing fluid equations. The calling material owner
/// supplies its phase/permeability closure and actual h/T epoch; this solver
/// does not infer a solid-locking threshold or turn drag work into heat.
pub fn simple_flow_with_material_resistance(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    dynamic_viscosity_pa_s: &[f64],
    drag_pa_s_m2: &[f64],
    force_density_n_m3: Option<&[[f64; 3]]>,
    config: &SimpleConfig,
    gate: &CancelGate,
) -> Result<FvFlow, ChtError> {
    material_flow(
        domain,
        fluid,
        dynamic_viscosity_pa_s,
        force_density_n_m3,
        Some(drag_pa_s_m2),
        config,
        gate,
    )
}

fn material_flow(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    dynamic_viscosity_pa_s: &[f64],
    force_density_n_m3: Option<&[[f64; 3]]>,
    drag_pa_s_m2: Option<&[f64]>,
    config: &SimpleConfig,
    gate: &CancelGate,
) -> Result<FvFlow, ChtError> {
    let solver = material_solver(
        domain,
        fluid,
        dynamic_viscosity_pa_s,
        force_density_n_m3,
        drag_pa_s_m2,
        config,
        gate,
    )?;
    solve_steady(solver, fluid, config, gate)
}

pub(super) fn material_solver<'a>(
    domain: &'a VoxelDomain,
    fluid: &FluidProperties,
    dynamic_viscosity_pa_s: &[f64],
    force_density_n_m3: Option<&[[f64; 3]]>,
    drag_pa_s_m2: Option<&[f64]>,
    config: &'a SimpleConfig,
    gate: &CancelGate,
) -> Result<Solver<'a>, ChtError> {
    poll(gate)?;
    admit(domain, fluid, config)?;
    if dynamic_viscosity_pa_s.len() != domain.cell_count()
        || config.turbulence != Turbulence::Laminar
    {
        return Err(ChtError::InvalidInput {
            field: "simple.dynamic_viscosity_pa_s",
            reason: "one molecular viscosity per voxel and the laminar closure are required".into(),
        });
    }
    for (c, &mu) in dynamic_viscosity_pa_s.iter().enumerate() {
        if c % 256 == 0 {
            poll(gate)?;
        }
        finite_positive("simple.dynamic_viscosity_pa_s", mu)?;
        finite_positive("simple.kinematic_viscosity_m2_s", mu / fluid.density_kg_m3)?;
    }
    if let Some(force) = force_density_n_m3 {
        if force.len() != domain.cell_count() {
            return Err(ChtError::InvalidInput {
                field: "simple.force_density_n_m3",
                reason: "one body-force vector per voxel is required".into(),
            });
        }
        for (c, components) in force.iter().enumerate() {
            if c % 256 == 0 {
                poll(gate)?;
            }
            for &value in components {
                finite("simple.force_density_n_m3", value)?;
                if !domain.is_fluid(c) && value != 0.0 {
                    return Err(ChtError::InvalidInput {
                        field: "simple.force_density_n_m3",
                        reason: "body force must be zero in solid voxels".into(),
                    });
                }
            }
        }
    }
    let mut solver = Solver::new(domain, fluid, config);
    // Existing Darcy coefficients are linear in molecular viscosity. Keep
    // their geometric/permeability owner and update only that coefficient.
    for (c, axes) in solver.porous.iter_mut().enumerate() {
        for (viscous, _) in axes {
            *viscous = (*viscous / solver.mu) * dynamic_viscosity_pa_s[c];
            finite("simple.material_darcy_resistance", *viscous)?;
        }
    }
    if let Some(drag) = drag_pa_s_m2 {
        if drag.len() != domain.cell_count() {
            return Err(ChtError::InvalidInput {
                field: "simple.drag_pa_s_m2",
                reason: "one isotropic drag coefficient per voxel is required".into(),
            });
        }
        if solver.porous.is_empty() {
            solver.porous = vec![[(0.0, 0.0); 3]; domain.cell_count()];
        }
        for (c, &coefficient) in drag.iter().enumerate() {
            if c % 256 == 0 {
                poll(gate)?;
            }
            finite("simple.drag_pa_s_m2", coefficient)?;
            if coefficient < 0.0 || (!domain.is_fluid(c) && coefficient != 0.0) {
                return Err(ChtError::InvalidInput {
                    field: "simple.drag_pa_s_m2",
                    reason: "drag must be nonnegative and zero in solid voxels".into(),
                });
            }
            for (viscous, _) in &mut solver.porous[c] {
                *viscous += coefficient;
                finite("simple.material_darcy_resistance", *viscous)?;
            }
        }
    }
    solver.molecular_viscosity = dynamic_viscosity_pa_s.to_vec();
    if let Some(force) = force_density_n_m3 {
        solver.set_force(force.to_vec());
    }
    Ok(solver)
}

fn solve_steady(
    mut solver: Solver<'_>,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    gate: &CancelGate,
) -> Result<FvFlow, ChtError> {
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

#[cfg(test)]
#[path = "simple/material_viscosity_tests.rs"]
mod material_viscosity_tests;

#[cfg(test)]
#[path = "simple/turbulence_tests.rs"]
mod turbulence_tests;
