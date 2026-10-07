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
//! `Outlet` (prescribed pressure zero, zero-gradient velocity; its normal
//! velocity is corrected by the pressure correction like an interior face,
//! so the outlet absorbs the mass imbalance).
//!
//! # SIMPLEC iteration
//!
//! Momentum is under-relaxed (`a_P / alpha`) and solved per component with
//! ILU(0)-BiCGStab; the pressure correction uses `d = A / (a_P / alpha -
//! sum a_nb)` and is solved with ILU(0)-PCG (components without an outlet are
//! pinned); velocities and pressure are corrected with no pressure
//! under-relaxation. Convergence requires two steady residuals of the same
//! iterate below the tolerance: the largest cell mass imbalance (before
//! correction) over the largest face mass flux, and, per velocity component,
//! the relative residual `||b - A u|| / ||b||` of the Jacobi-scaled momentum
//! system at the velocities entering the iteration (under-relaxation cancels
//! there, so it is the residual of the unrelaxed equations). Inner solves
//! only reduce that residual by `momentum_tolerance`, so a converged report
//! never rests on a skipped inner solve.
//!
//! # No-claim boundaries
//!
//! Steady, laminar, constant-property, incompressible flow; no turbulence
//! model (a laminar solution above transition is a laminar idealization, not
//! a prediction), no buoyancy, staircase voxel walls, power-law convection
//! (first order at high cell Péclet numbers: numerical diffusion is not
//! bounded here). One run makes no mesh-convergence claim.

use fs_exec::CancelGate;
use fs_sparse::Coo;

use super::domain::{FluidProperties, VoxelDomain};
use super::energy::ConvectionScheme;
use super::flow::{FlowField, scale_rows};
use super::krylov::bicgstab_ilu0;
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

/// SIMPLEC controls.
#[derive(Debug, Clone, Copy, PartialEq)]
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
    /// Relative residual of each pressure-correction solve.
    pub pressure_tolerance: f64,
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
            momentum_tolerance: 1e-2,
            pressure_tolerance: 1e-8,
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
    /// Total inflow through inlet faces, m^3/s.
    pub inflow_m3_s: f64,
    /// Total outflow through outlet faces, m^3/s.
    pub outflow_m3_s: f64,
    /// Largest cell Reynolds number `|u| dx / nu`.
    pub max_cell_reynolds: f64,
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
    /// Run evidence.
    pub report: SimpleReport,
}

impl FvFlow {
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
    /// Outlet face: extrapolated, then pressure-corrected.
    Outlet,
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

struct Solver<'a> {
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
}

fn norm(v: &[f64]) -> f64 {
    fs_math::det::sqrt(v.iter().map(|x| x * x).sum::<f64>())
}

fn cell_of(n: [usize; 3], c: [usize; 3]) -> usize {
    (c[2] * n[1] + c[1]) * n[0] + c[0]
}

impl<'a> Solver<'a> {
    fn new(domain: &'a VoxelDomain, fluid: &FluidProperties, config: &'a SimpleConfig) -> Self {
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
                        FvBoundary::Outlet => Kind::Outlet,
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
        }
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
    /// velocities on entry: `||b - A u|| / ||b||` of the row-scaled system
    /// (under-relaxation cancels at `u = u_old`, so this is the residual of
    /// the unrelaxed steady equations).
    #[allow(clippy::too_many_lines)] // one staggered control-volume assembly
    fn momentum(&mut self, a: usize, gate: &CancelGate) -> Result<f64, ChtError> {
        let comp = &self.comps[a];
        let rows = comp.unknowns.len();
        if rows == 0 {
            return Ok(0.0);
        }
        let alpha = self.config.velocity_relaxation;
        let dx = self.domain.dx();
        let diff = self.mu * self.area / dx;
        let mut coo = Coo::new(rows, rows);
        let mut b = vec![0.0f64; rows];
        let mut a_p_relaxed = vec![0.0f64; rows];
        let mut neighbour_sum = vec![0.0f64; rows];
        for (row, &index) in comp.unknowns.iter().enumerate() {
            if row % 4096 == 0 {
                poll(gate)?;
            }
            let f = comp.coords(index);
            let mut minus_cell = f;
            minus_cell[a] -= 1;
            let plus_cell = f;
            let mut sum_nb = 0.0f64;
            let mut a_p = 0.0f64;
            let mut rhs = 0.0f64;
            let mut net_out = 0.0f64;
            for d in 0..3 {
                for plus in [false, true] {
                    let s = if plus { 1.0 } else { -1.0 };
                    if d == a {
                        // CV side at the centre of the minus/plus cell.
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
                        let coef = diff.mul_add(self.weight(flux / diff), (-flux).max(0.0));
                        match comp.kind[nf] {
                            Kind::Unknown(col) => {
                                coo.push(row, col, -coef);
                                sum_nb += coef;
                            }
                            Kind::Fixed(v) => rhs = coef.mul_add(v, rhs),
                            Kind::Outlet => rhs = coef.mul_add(self.vel[a][nf], rhs),
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
                    if outside {
                        match self.config.faces[2 * d + usize::from(plus)] {
                            FvBoundary::Wall { velocity } => {
                                a_p += 2.0 * diff;
                                rhs = (2.0 * diff).mul_add(velocity[a], rhs);
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
                        Kind::Unknown(col) => {
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
                            a_p += 2.0f64.mul_add(diff, (-edge_flux).max(0.0));
                        }
                    }
                }
            }
            // Continuity is satisfied only at convergence: keep a_P >= sum.
            a_p += net_out.max(0.0);
            // Unknown faces separate two fluid cells.
            let drop = self.pressure[cell_of(self.n, minus_cell)]
                - self.pressure[cell_of(self.n, plus_cell)];
            rhs = self.area.mul_add(drop, rhs);
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
        let b_norm = norm(&b);
        let steady = if b_norm > 0.0 {
            fs_math::det::sqrt(
                r.iter()
                    .zip(&b)
                    .map(|(ax, bi)| (bi - ax) * (bi - ax))
                    .sum::<f64>(),
            ) / b_norm
        } else {
            norm(&r)
        };
        // Inner solves reduce the entry residual by `momentum_tolerance`.
        let inner = (self.config.momentum_tolerance * steady).max(1e-15);
        if steady > 1e-15 {
            bicgstab_ilu0("momentum", &matrix, &b, &mut x, inner, 20_000, gate)?;
        }
        for (row, &index) in self.comps[a].unknowns.iter().enumerate() {
            self.vel[a][index] = x[row];
            self.d[a][index] =
                self.area / (a_p_relaxed[row] - neighbour_sum[row]).max(f64::MIN_POSITIVE);
        }
        // Outlet faces: zero-gradient extrapolation from the interior face
        // upstream; their d follows the same face.
        for index in 0..self.comps[a].kind.len() {
            if self.comps[a].kind[index] != Kind::Outlet {
                continue;
            }
            let f = self.comps[a].coords(index);
            let mut g = f;
            if f[a] == 0 {
                g[a] += 1;
            } else {
                g[a] -= 1;
            }
            let upstream = self.comps[a].index(g);
            let (value, d) = match self.comps[a].kind[upstream] {
                Kind::Unknown(_) => (self.vel[a][upstream], self.d[a][upstream]),
                _ => (self.vel[a][index], 0.5 * dx / self.mu),
            };
            // An outlet admits no inflow from outside in the extrapolation.
            let outward = if f[a] == 0 { -value } else { value };
            self.vel[a][index] = if outward < 0.0 { 0.0 } else { value };
            self.d[a][index] = d;
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

    /// Pressure correction; returns the largest cell mass imbalance before
    /// correction.
    fn correct(&mut self, gate: &CancelGate) -> Result<f64, ChtError> {
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
                    Kind::Outlet => {
                        diag += self.rho * self.area * self.d[axis][face];
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
        let matrix = scale_rows(&coo, &mut b);
        let mut correction = vec![0.0f64; rows];
        if b.iter().any(|v| *v != 0.0) {
            bicgstab_ilu0(
                "pressure",
                &matrix,
                &b,
                &mut correction,
                self.config.pressure_tolerance,
                20_000,
                gate,
            )?;
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
                    Kind::Outlet => {
                        if f[a] == 0 {
                            self.vel[a][index] -= self.d[a][index] * p_of(f);
                        } else {
                            let mut minus = f;
                            minus[a] -= 1;
                            self.vel[a][index] += self.d[a][index] * p_of(minus);
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

/// Solve steady incompressible flow on `domain` by SIMPLEC.
///
/// # Errors
/// Input refusals (non-finite rules, an inlet pointing outward, a wall with
/// normal velocity), [`ChtError::FlowNotSteady`] when the iteration budget
/// ends first, solver refusals, or [`ChtError::Cancelled`].
#[allow(clippy::too_many_lines)] // admission, iteration, flux handover, report
pub fn simple_flow(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    gate: &CancelGate,
) -> Result<FvFlow, ChtError> {
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
    if domain.fluid_count() == 0 {
        return Err(ChtError::InvalidDomain {
            reason: "no fluid cell".into(),
        });
    }
    let mut solver = Solver::new(domain, fluid, config);
    let mut iterations = 0usize;
    let mut mass_residual = f64::INFINITY;
    let mut momentum_residual = f64::INFINITY;
    while iterations < config.max_iterations {
        poll(gate)?;
        iterations += 1;
        let mut steady = 0.0f64;
        for a in 0..3 {
            steady = steady.max(solver.momentum(a, gate)?);
        }
        let imbalance = solver.correct(gate)?;
        let largest_velocity = solver
            .vel
            .iter()
            .flat_map(|v| v.iter())
            .fold(0.0f64, |m, v| m.max(v.abs()))
            .max(f64::MIN_POSITIVE);
        let largest_flux = solver.rho * solver.area * largest_velocity;
        mass_residual = imbalance / largest_flux;
        momentum_residual = steady;
        if !(mass_residual.is_finite() && momentum_residual.is_finite()) {
            return Err(ChtError::FlowDiverged { step: iterations });
        }
        if mass_residual <= config.tolerance && momentum_residual <= config.tolerance {
            break;
        }
    }
    if mass_residual > config.tolerance || momentum_residual > config.tolerance {
        return Err(ChtError::FlowNotSteady {
            steps: iterations,
            last_change: mass_residual.max(momentum_residual),
            tolerance: config.tolerance,
        });
    }
    let area = solver.area;
    let fluxes = solver
        .vel
        .clone()
        .map(|v| v.into_iter().map(|u| u * area).collect::<Vec<f64>>());
    let [fx, fy, fz] = fluxes;
    let field = FlowField::from_face_arrays(domain, fx, fy, fz);
    let mut velocity_m_s = vec![[0.0f64; 3]; domain.cell_count()];
    let mut max_cell_reynolds = 0.0f64;
    for (c, slot) in velocity_m_s.iter_mut().enumerate() {
        if !domain.is_fluid(c) {
            solver.pressure[c] = 0.0;
            continue;
        }
        let at = domain.coords(c);
        for a in 0..3 {
            slot[a] = 0.5
                * (solver.vel[a][solver.cell_face(a, at, false)]
                    + solver.vel[a][solver.cell_face(a, at, true)]);
        }
        let speed = fs_math::det::sqrt(slot.iter().map(|v| v * v).sum::<f64>());
        max_cell_reynolds =
            max_cell_reynolds.max(speed * domain.dx() / fluid.kinematic_viscosity_m2_s);
    }
    let (mut inflow, mut outflow) = (0.0f64, 0.0f64);
    for face in Face3::ALL {
        let net = field.boundary_outflow(domain, face);
        match config.faces[face as usize] {
            FvBoundary::Inlet { .. } => inflow -= net,
            FvBoundary::Outlet => outflow += net,
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
    };
    Ok(FvFlow {
        field,
        velocity_m_s,
        pressure_pa: solver.pressure,
        report,
    })
}
