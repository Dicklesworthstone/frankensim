//! Staggered volumetric face fluxes, their discrete projection onto the
//! divergence-free subspace, and the steady D3Q19 duct-flow producer.

// Local face numbers index both the neighbour stencil and the face-rule
// table; iterating them as integers is the clearest form.
#![allow(clippy::needless_range_loop)]

use fs_exec::{CancelGate, TilePool};
use fs_sparse::Coo;

use super::domain::{FluidProperties, VoxelDomain};
use super::krylov::bicgstab_ilu0;
use super::{ChtError, finite, finite_positive, poll};
use crate::d3q19::equilibrium3;
use crate::d3q19::{
    BoundaryGrid3, BoundarySpec3, BoundaryStepError3, CollisionModel3, E3, Face3, FaceBoundary3,
    TILE,
};

/// Mass-flow role of one domain face for flux construction and projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowFace {
    /// No flow crosses the face.
    Wall,
    /// Open face whose supplied flux is held fixed by the projection (an
    /// analytic or measured inflow profile).
    Fixed,
    /// Open face where the projection's correction potential vanishes, so
    /// its flux adjusts to the interior field (an outflow, or an LBM on-site
    /// boundary whose cell velocity is not the transported mass flux).
    Free,
}

/// Volumetric fluxes (m^3/s, positive along `+axis`) through every cell face
/// of a [`VoxelDomain`].
///
/// Layout: `x` faces `(z ny + y)(nx + 1) + i`, `y` faces `(z (ny + 1) + j) nx
/// + x`, `z` faces `(k ny + y) nx + x`, where `i`, `j`, `k` index the face
/// planes `0..=n`. Faces touching a solid cell carry exactly zero flux.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowField {
    dims: [usize; 3],
    dx: f64,
    fx: Vec<f64>,
    fy: Vec<f64>,
    fz: Vec<f64>,
}

/// Evidence retained by the divergence-free projection.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectionReport {
    /// Fluid cells (unknowns of the correction-potential system).
    pub fluid_cells: usize,
    /// Fluid connected components with no outflow face, each pinned at one
    /// cell (their net imbalance cannot be removed and stays in the report).
    pub pinned_components: usize,
    /// BiCGStab iterations.
    pub iterations: usize,
    /// Recomputed relative residual of the row-scaled potential system.
    pub relative_residual: f64,
    /// Largest |face flux| after projection, m^3/s.
    pub max_flux_m3_s: f64,
    /// Largest per-cell |net outflow| before projection, m^3/s.
    pub max_divergence_before_m3_s: f64,
    /// Largest per-cell |net outflow| after projection, m^3/s.
    pub max_divergence_after_m3_s: f64,
    /// Largest |flux correction| applied to any face, m^3/s.
    pub max_correction_m3_s: f64,
}

impl FlowField {
    /// No flow anywhere (pure conduction).
    #[must_use]
    pub fn quiescent(domain: &VoxelDomain) -> Self {
        let [nx, ny, nz] = domain.dims();
        Self {
            dims: [nx, ny, nz],
            dx: domain.min_width(),
            fx: vec![0.0; (nx + 1) * ny * nz],
            fy: vec![0.0; nx * (ny + 1) * nz],
            fz: vec![0.0; nx * ny * (nz + 1)],
        }
    }

    /// Adopt face fluxes (m^3/s) already laid out as this type stores them.
    pub(crate) fn from_face_arrays(
        domain: &VoxelDomain,
        fx: Vec<f64>,
        fy: Vec<f64>,
        fz: Vec<f64>,
    ) -> Self {
        let [nx, ny, nz] = domain.dims();
        debug_assert_eq!(fx.len(), (nx + 1) * ny * nz);
        debug_assert_eq!(fy.len(), nx * (ny + 1) * nz);
        debug_assert_eq!(fz.len(), nx * ny * (nz + 1));
        Self {
            dims: [nx, ny, nz],
            dx: domain.min_width(),
            fx,
            fy,
            fz,
        }
    }

    /// The x, y and z face-flux arrays, m^3/s (the layout of
    /// [`Self::from_face_arrays`]).
    pub(crate) fn face_arrays(&self) -> [&[f64]; 3] {
        [&self.fx, &self.fy, &self.fz]
    }

    /// Sample a face-normal velocity (m/s) at every face centre that carries
    /// flow: interior faces between two fluid cells, and boundary faces of
    /// fluid cells on [`FlowFace::Fixed`] / [`FlowFace::Free`] domain
    /// faces. `velocity(axis, centre)` returns the `axis` component. No
    /// projection is applied; the caller owns the divergence of the result,
    /// which [`FlowField::max_divergence`] measures.
    ///
    /// # Errors
    /// [`ChtError::InvalidInput`] for a non-finite sample.
    pub fn from_face_velocity(
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        mut velocity: impl FnMut(usize, [f64; 3]) -> f64,
    ) -> Result<Self, ChtError> {
        let mut field = Self::quiescent(domain);
        for c in 0..domain.cell_count() {
            if !domain.is_fluid(c) {
                continue;
            }
            let [x, y, z] = domain.coords(c);
            for f in [1usize, 3, 5, 0, 2, 4] {
                let carries = match domain.neighbor(c, f) {
                    Some(n) => f % 2 == 1 && domain.is_fluid(n),
                    None => faces[f] != FlowFace::Wall,
                };
                if !carries {
                    continue;
                }
                let axis = f / 2;
                let mut centre = domain.center(x, y, z);
                let half = if f % 2 == 1 { 0.5 } else { -0.5 };
                centre[axis] += half * domain.widths(c)[axis];
                let u = velocity(axis, centre);
                finite("flow.face_velocity", u)?;
                let (slot, _) = field.slot(domain, c, f);
                *field.face_mut(axis, slot) = u * domain.face_area(c, axis);
            }
        }
        Ok(field)
    }

    /// Interpolate cell-centred velocities (m/s, one per cell; solid entries
    /// are ignored) to faces and project the result onto discretely
    /// divergence-free fluxes. Interior fluid faces take the arithmetic mean
    /// of their two cells; open boundary faces take the boundary cell's
    /// velocity. [`FlowFace::Fixed`] fluxes are held fixed and the correction
    /// potential vanishes on [`FlowFace::Free`] faces.
    ///
    /// # Errors
    /// [`ChtError::InvalidInput`] for a wrong-length or non-finite velocity
    /// array; solver refusals from the potential solve; [`ChtError::Cancelled`].
    pub fn from_cell_velocities(
        domain: &VoxelDomain,
        velocities: &[[f64; 3]],
        faces: [FlowFace; 6],
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<(Self, ProjectionReport), ChtError> {
        if velocities.len() != domain.cell_count() {
            return Err(ChtError::InvalidInput {
                field: "flow.velocities",
                reason: format!(
                    "expected {} cell velocities, got {}",
                    domain.cell_count(),
                    velocities.len()
                ),
            });
        }
        finite_positive("flow.projection_tolerance", tolerance)?;
        let dx = domain.require_uniform("the cell-velocity projection")?;
        let mut field = Self::quiescent(domain);
        let area = dx * dx;
        for c in 0..domain.cell_count() {
            if !domain.is_fluid(c) {
                continue;
            }
            for component in velocities[c] {
                finite("flow.velocities", component)?;
            }
            for f in [1usize, 3, 5, 0, 2, 4] {
                let axis = f / 2;
                let value = match domain.neighbor(c, f) {
                    Some(n) if f % 2 == 1 && domain.is_fluid(n) => {
                        0.5 * (velocities[c][axis] + velocities[n][axis])
                    }
                    None if faces[f] != FlowFace::Wall => velocities[c][axis],
                    // Plus-side copies of interior faces, solid neighbours,
                    // and walls carry no flux from this cell.
                    _ => continue,
                };
                let (slot, _) = field.slot(domain, c, f);
                *field.face_mut(axis, slot) = value * area;
            }
        }
        let report = field.project(domain, faces, tolerance, gate)?;
        Ok((field, report))
    }

    /// Domain dimensions this field was built for.
    #[must_use]
    pub const fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Flux through the `x` face plane `i` (`0..=nx`) at `(y, z)`, m^3/s.
    #[must_use]
    pub fn flux_x(&self, i: usize, y: usize, z: usize) -> f64 {
        let [nx, ny, _] = self.dims;
        self.fx[(z * ny + y) * (nx + 1) + i]
    }

    /// Flux through the `y` face plane `j` (`0..=ny`) at `(x, z)`, m^3/s.
    #[must_use]
    pub fn flux_y(&self, x: usize, j: usize, z: usize) -> f64 {
        let [nx, ny, _] = self.dims;
        self.fy[(z * (ny + 1) + j) * nx + x]
    }

    /// Flux through the `z` face plane `k` (`0..=nz`) at `(x, y)`, m^3/s.
    #[must_use]
    pub fn flux_z(&self, x: usize, y: usize, k: usize) -> f64 {
        let [nx, ny, _] = self.dims;
        self.fz[(k * ny + y) * nx + x]
    }

    /// Outward flux of cell `c` through local face `f` ([`Face3::ALL`]
    /// order), m^3/s.
    #[must_use]
    pub fn outward(&self, domain: &VoxelDomain, c: usize, f: usize) -> f64 {
        let (slot, sign) = self.slot(domain, c, f);
        sign * self.face(f / 2, slot)
    }

    /// Net outward flux of cell `c`, m^3/s.
    #[must_use]
    pub fn net_outflow(&self, domain: &VoxelDomain, c: usize) -> f64 {
        (0..6).map(|f| self.outward(domain, c, f)).sum()
    }

    /// Largest per-cell |net outflow|, m^3/s.
    #[must_use]
    pub fn max_divergence(&self, domain: &VoxelDomain) -> f64 {
        (0..domain.cell_count())
            .map(|c| self.net_outflow(domain, c).abs())
            .fold(0.0, f64::max)
    }

    /// Largest |face flux|, m^3/s.
    #[must_use]
    pub fn max_flux(&self) -> f64 {
        self.fx
            .iter()
            .chain(&self.fy)
            .chain(&self.fz)
            .fold(0.0f64, |m, v| m.max(v.abs()))
    }

    /// Net outward flux through one whole domain face, m^3/s.
    #[must_use]
    pub fn boundary_outflow(&self, domain: &VoxelDomain, face: Face3) -> f64 {
        let f = face as usize;
        let mut total = 0.0;
        for c in 0..domain.cell_count() {
            if domain.neighbor(c, f).is_none() {
                total += self.outward(domain, c, f);
            }
        }
        total
    }

    /// Flux slot of cell `c`'s local face `f` and the outward sign.
    pub(crate) fn slot(&self, domain: &VoxelDomain, c: usize, f: usize) -> (usize, f64) {
        let [nx, ny, _] = self.dims;
        let [x, y, z] = domain.coords(c);
        let plus = f % 2 == 1;
        let sign = if plus { 1.0 } else { -1.0 };
        let o = usize::from(plus);
        let slot = match f / 2 {
            0 => (z * ny + y) * (nx + 1) + x + o,
            1 => (z * (ny + 1) + y + o) * nx + x,
            _ => ((z + o) * ny + y) * nx + x,
        };
        (slot, sign)
    }

    fn face(&self, axis: usize, slot: usize) -> f64 {
        match axis {
            0 => self.fx[slot],
            1 => self.fy[slot],
            _ => self.fz[slot],
        }
    }

    fn face_mut(&mut self, axis: usize, slot: usize) -> &mut f64 {
        match axis {
            0 => &mut self.fx[slot],
            1 => &mut self.fy[slot],
            _ => &mut self.fz[slot],
        }
    }

    /// Remove the discrete divergence with one correction potential `phi` on
    /// the fluid cells: corrected flux `F - (A/dx)(phi_N - phi_P)` on interior
    /// faces and `F + (2A/dx) phi_P` on outflow faces (`phi = 0` there).
    #[allow(clippy::too_many_lines)] // components, assembly, solve, correction
    fn project(
        &mut self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<ProjectionReport, ChtError> {
        let max_divergence_before_m3_s = self.max_divergence(domain);
        let fluid: Vec<usize> = (0..domain.cell_count())
            .filter(|&c| domain.is_fluid(c))
            .collect();
        let mut row_of = vec![usize::MAX; domain.cell_count()];
        for (row, &c) in fluid.iter().enumerate() {
            row_of[c] = row;
        }
        let n = fluid.len();
        // Components without an outflow face are pinned at their first cell.
        let mut component = vec![usize::MAX; n];
        let mut pins = Vec::new();
        for start in 0..n {
            if component[start] != usize::MAX {
                continue;
            }
            let id = pins.len();
            let mut stack = vec![start];
            component[start] = id;
            let mut drains = false;
            while let Some(row) = stack.pop() {
                let c = fluid[row];
                for f in 0..6 {
                    match domain.neighbor(c, f) {
                        Some(nb) if domain.is_fluid(nb) => {
                            let r = row_of[nb];
                            if component[r] == usize::MAX {
                                component[r] = id;
                                stack.push(r);
                            }
                        }
                        Some(_) => {}
                        None => drains |= faces[f] == FlowFace::Free,
                    }
                }
            }
            pins.push((start, drains));
        }
        let mut is_pinned = vec![false; n];
        let mut pinned_components = 0usize;
        for &(row, drains) in &pins {
            if !drains {
                is_pinned[row] = true;
                pinned_components += 1;
            }
        }
        let g = domain.dx(); // A / dx
        let mut coo = Coo::new(n, n);
        let mut b = vec![0.0f64; n];
        for (row, &c) in fluid.iter().enumerate() {
            poll_every(gate, row)?;
            let mut diag = 0.0;
            let rhs = -self.net_outflow(domain, c);
            for f in 0..6 {
                match domain.neighbor(c, f) {
                    Some(nb) if domain.is_fluid(nb) => {
                        diag += g;
                        coo.push(row, row_of[nb], -g);
                    }
                    None if faces[f] == FlowFace::Free => diag += 2.0 * g,
                    _ => {}
                }
            }
            // A pinned cell regularizes its otherwise singular Neumann
            // component; a component with zero net boundary flux then
            // solves to phi_pin = 0 exactly and loses no equation.
            if is_pinned[row] {
                diag += g;
            }
            coo.push(row, row, diag);
            b[row] = rhs;
        }
        let a = scale_rows(&coo, &mut b);
        let mut phi = vec![0.0f64; n];
        let outcome = bicgstab_ilu0(
            "projection",
            &a,
            &b,
            &mut phi,
            tolerance,
            50 * n.max(20),
            gate,
        )?;
        let mut max_correction_m3_s = 0.0f64;
        for (row, &c) in fluid.iter().enumerate() {
            for f in 0..6 {
                let correction = match domain.neighbor(c, f) {
                    // Each interior face once, from its minus-side cell.
                    Some(nb) if f % 2 == 1 && domain.is_fluid(nb) => {
                        -g * (phi[row_of[nb]] - phi[row])
                    }
                    None if faces[f] == FlowFace::Free => 2.0 * g * phi[row],
                    _ => continue,
                };
                let (slot, sign) = self.slot(domain, c, f);
                *self.face_mut(f / 2, slot) += sign * correction;
                max_correction_m3_s = max_correction_m3_s.max(correction.abs());
            }
        }
        Ok(ProjectionReport {
            fluid_cells: n,
            pinned_components,
            iterations: outcome.iterations,
            relative_residual: outcome.relative_residual,
            max_flux_m3_s: self.max_flux(),
            max_divergence_before_m3_s,
            max_divergence_after_m3_s: self.max_divergence(domain),
            max_correction_m3_s,
        })
    }
}

fn poll_every(gate: &CancelGate, i: usize) -> Result<(), ChtError> {
    if i.is_multiple_of(4096) {
        poll(gate)
    } else {
        Ok(())
    }
}

/// Jacobi row scaling `D^-1 A x = D^-1 b`, so residuals compare rows in the
/// unknown's own units regardless of the row's conductance magnitude.
pub(crate) fn scale_rows(coo: &Coo, b: &mut [f64]) -> fs_sparse::Csr {
    let a = coo.assemble();
    let n = a.nrows();
    let mut row_ptr = Vec::with_capacity(n + 1);
    let mut cols = Vec::with_capacity(a.nnz());
    let mut vals = Vec::with_capacity(a.nnz());
    row_ptr.push(0);
    for (r, br) in b.iter_mut().enumerate().take(n) {
        let (c, v) = a.row(r);
        let d = a.get(r, r);
        let s = if d == 0.0 { 1.0 } else { 1.0 / d };
        cols.extend_from_slice(c);
        vals.extend(v.iter().map(|x| x * s));
        *br *= s;
        row_ptr.push(cols.len());
    }
    fs_sparse::Csr::from_parts(n, n, row_ptr, cols, vals)
}

/// Collision operator selection for [`lbm_duct_flow`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LbmCollisionChoice {
    /// BGK. Measured on the plate-fin heatsink example at `tau = 0.519`:
    /// the central-moment operator diverged at step 85 while BGK stayed
    /// bounded through 2500 steps, so `Auto` does not switch operators
    /// near `tau = 1/2`; `min_tau` is the admission floor instead.
    #[default]
    Auto,
    /// Single-relaxation-time BGK.
    Bgk,
    /// Central-moment relaxation: second-order rate `1/tau`, higher rate 1.
    CentralMoment,
}

/// Steady duct-flow run configuration. Inflow enters through `x-min` with a
/// uniform velocity; `x-max` is a constant-pressure outlet; the four other
/// domain faces and every solid voxel are no-slip walls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LbmFlowConfig {
    /// Physical inlet velocity along `+x`, m/s.
    pub inlet_velocity_m_s: f64,
    /// Lattice inlet velocity (Mach `= u sqrt(3)`); fixes the time step.
    pub lattice_inlet_velocity: f64,
    /// Collision operator.
    pub collision: LbmCollisionChoice,
    /// Smallest admitted relaxation time.
    pub min_tau: f64,
    /// Step budget.
    pub max_steps: usize,
    /// Steps between steady checks.
    pub check_interval: usize,
    /// Steady criterion: relative L2 velocity change between checks.
    pub steady_tolerance: f64,
    /// Relative residual for the divergence-free projection.
    pub projection_tolerance: f64,
    /// Pool workers for the collide/stream passes; `0` means every core the
    /// host reports. Results are bit-identical for any value.
    pub workers: usize,
}

impl Default for LbmFlowConfig {
    fn default() -> Self {
        Self {
            inlet_velocity_m_s: 1.0,
            lattice_inlet_velocity: 0.05,
            collision: LbmCollisionChoice::Auto,
            min_tau: 0.505,
            max_steps: 400_000,
            check_interval: 100,
            steady_tolerance: 1e-7,
            projection_tolerance: 1e-12,
            workers: 0,
        }
    }
}

/// Evidence retained by one steady LBM duct-flow run.
#[derive(Debug, Clone, PartialEq)]
pub struct LbmFlowReport {
    /// Steps executed.
    pub steps: usize,
    /// Relative L2 velocity change at the last check.
    pub last_change: f64,
    /// Relaxation time `3 nu_lat + 1/2`.
    pub tau: f64,
    /// Lattice kinematic viscosity.
    pub lattice_viscosity: f64,
    /// Lattice inlet velocity.
    pub lattice_inlet_velocity: f64,
    /// Inlet Mach number `u_lat sqrt(3)`.
    pub inlet_mach: f64,
    /// Largest lattice speed in the converged field.
    pub max_lattice_speed: f64,
    /// Cell Reynolds number `U dx / nu`.
    pub cell_reynolds: f64,
    /// Nominal inflow `U` times the open inlet area, m^3/s.
    pub nominal_inflow_m3_s: f64,
    /// Volumetric inflow the converged lattice actually admits (mass flux
    /// over reference density through the inlet face), m^3/s. Inlet rim
    /// cells share bounce-back links with the walls, so this sits below the
    /// nominal value on coarse grids; downstream physics uses this value.
    pub realized_inflow_m3_s: f64,
    /// Metres per second represented by one lattice velocity unit.
    pub velocity_scale_m_s: f64,
    /// Mean gauge pressure of the first interior layer minus that of the
    /// last interior layer, Pa (the on-site boundary layers are excluded).
    pub pressure_drop_pa: f64,
    /// Collision operator actually used.
    pub collision: CollisionModel3,
    /// The projection that produced the returned fluxes.
    pub projection: ProjectionReport,
}

/// Steady duct flow in physical units.
#[derive(Debug, Clone, PartialEq)]
pub struct LbmFlow {
    /// Projected, discretely divergence-free face fluxes.
    pub field: FlowField,
    /// Cell mass flux over reference density, m/s (zero in solids). This is
    /// the incompressible velocity the lattice transports.
    pub velocity_m_s: Vec<[f64; 3]>,
    /// Gauge pressure `c_s^2 (rho_lat - 1) rho S^2`, Pa, relative to the
    /// outlet reference density (zero in solids), where `S` is the velocity
    /// scale.
    pub pressure_pa: Vec<f64>,
    /// Run evidence.
    pub report: LbmFlowReport,
}

impl LbmFlow {
    /// Mean gauge pressure over the fluid cells of layer `x`, Pa.
    #[must_use]
    pub fn mean_pressure_x(&self, domain: &VoxelDomain, x: usize) -> Option<f64> {
        let [_, ny, nz] = domain.dims();
        let (mut sum, mut count) = (0.0f64, 0usize);
        for z in 0..nz {
            for y in 0..ny {
                let c = domain.index(x, y, z);
                if domain.is_fluid(c) {
                    sum += self.pressure_pa[c];
                    count += 1;
                }
            }
        }
        (count > 0).then(|| sum / count as f64)
    }
}

/// Run the D3Q19 duct flow to a steady state and return projected physical
/// face fluxes, physical cell velocities and gauge pressures (zero in
/// solids), and the run report.
///
/// Requirements: every dimension a positive multiple of four (the D3Q19 tile
/// edge); every fluid cell on the `x-min`/`x-max` faces has a fluid inward
/// neighbour; at least one fluid cell.
///
/// # Errors
/// [`ChtError::InvalidDomain`] for an inadmissible lattice geometry,
/// [`ChtError::LatticeResolution`] when `tau < min_tau`,
/// [`ChtError::FlowDiverged`], [`ChtError::FlowNotSteady`], projection
/// refusals, or [`ChtError::Cancelled`].
#[allow(clippy::too_many_lines)] // one linear admission -> run -> projection story
pub fn lbm_duct_flow(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &LbmFlowConfig,
    gate: &CancelGate,
) -> Result<LbmFlow, ChtError> {
    domain.require_uniform("the lattice-Boltzmann duct flow")?;
    fluid.validate()?;
    finite_positive("lbm.inlet_velocity_m_s", config.inlet_velocity_m_s)?;
    finite_positive("lbm.lattice_inlet_velocity", config.lattice_inlet_velocity)?;
    finite_positive("lbm.steady_tolerance", config.steady_tolerance)?;
    if config.lattice_inlet_velocity > 0.15 {
        return Err(ChtError::InvalidInput {
            field: "lbm.lattice_inlet_velocity",
            reason: "must not exceed 0.15 (inlet Mach about 0.26)".into(),
        });
    }
    if !(config.min_tau.is_finite() && config.min_tau > 0.5) {
        return Err(ChtError::InvalidInput {
            field: "lbm.min_tau",
            reason: format!(
                "must be finite and greater than 0.5, got {}",
                config.min_tau
            ),
        });
    }
    if config.check_interval == 0 {
        return Err(ChtError::InvalidInput {
            field: "lbm.check_interval",
            reason: "must be positive".into(),
        });
    }
    let [nx, ny, nz] = domain.dims();
    if [nx, ny, nz].iter().any(|n| !n.is_multiple_of(TILE)) {
        return Err(ChtError::InvalidDomain {
            reason: format!("LBM dimensions must be multiples of {TILE}, got {nx}x{ny}x{nz}"),
        });
    }
    if domain.fluid_count() == 0 {
        return Err(ChtError::InvalidDomain {
            reason: "no fluid cell".into(),
        });
    }
    if nx < 2 {
        return Err(ChtError::InvalidDomain {
            reason: "need at least two x layers".into(),
        });
    }
    for z in 0..nz {
        for y in 0..ny {
            for (x, inner) in [(0, 1), (nx - 1, nx - 2)] {
                if domain.is_fluid(domain.index(x, y, z))
                    && !domain.is_fluid(domain.index(inner, y, z))
                {
                    return Err(ChtError::InvalidDomain {
                        reason: format!(
                            "open-face fluid cell ({x},{y},{z}) has a solid inward neighbour; pad the inlet/outlet with fluid"
                        ),
                    });
                }
            }
        }
    }

    let u_lat = config.lattice_inlet_velocity;
    let velocity_scale_m_s = config.inlet_velocity_m_s / u_lat;
    let lattice_viscosity = fluid.kinematic_viscosity_m2_s / (velocity_scale_m_s * domain.dx());
    let tau = 3.0f64.mul_add(lattice_viscosity, 0.5);
    let cell_reynolds = config.inlet_velocity_m_s * domain.dx() / fluid.kinematic_viscosity_m2_s;
    if tau < config.min_tau {
        let dx_needed = fluid.kinematic_viscosity_m2_s * u_lat
            / (config.inlet_velocity_m_s * (config.min_tau - 0.5) / 3.0);
        return Err(ChtError::LatticeResolution {
            tau,
            cell_reynolds,
            remedy: format!(
                "refine the voxel size to at most {dx_needed:.3e} m at this lattice velocity, or raise lattice_inlet_velocity (Mach permitting)"
            ),
        });
    }
    let bgk = match config.collision {
        LbmCollisionChoice::Bgk | LbmCollisionChoice::Auto => true,
        LbmCollisionChoice::CentralMoment => false,
    };
    let collision = if bgk {
        CollisionModel3::Bgk { tau }
    } else {
        CollisionModel3::CentralMoment {
            second_order_rate: 1.0 / tau,
            higher_order_rate: 1.0,
        }
    };
    let wall = FaceBoundary3::stationary_wall();
    let spec = BoundarySpec3::new([
        FaceBoundary3::Velocity {
            velocity: [u_lat, 0.0, 0.0],
        },
        FaceBoundary3::Pressure { density: 1.0 },
        wall,
        wall,
        wall,
        wall,
    ]);
    let mut grid = BoundaryGrid3::with_collision_model(nx, ny, nz, collision, [0.0; 3], spec);
    grid.voxelize_sdf(|p| {
        let (x, y, z) = (p[0] as usize, p[1] as usize, p[2] as usize);
        if domain.is_fluid(domain.index(x, y, z)) {
            1.0
        } else {
            -1.0
        }
    });
    // Plug-flow start shortens the transient; it is not part of the answer.
    let start = equilibrium3(1.0, [u_lat, 0.0, 0.0]);
    let fluid_cells: Vec<usize> = (0..domain.cell_count())
        .filter(|&c| domain.is_fluid(c))
        .collect();
    for &c in &fluid_cells {
        let [x, y, z] = domain.coords(c);
        grid.set_populations(x, y, z, &start);
    }
    // The weakly compressible lattice conserves MASS: density falls along
    // the pressure drop and the velocity rises with it. The incompressible
    // volumetric flux is therefore the momentum `sum e_q f_q` (force-free
    // grid) over the reference density 1, which is (to steady tolerance)
    // constant along a duct. Moments are taken from the populations
    // directly, so a diverging state is a refusal, never a panic.
    let sample = |grid: &BoundaryGrid3, step: usize| -> Result<Vec<[f64; 4]>, ChtError> {
        fluid_cells
            .iter()
            .map(|&c| {
                let [x, y, z] = domain.coords(c);
                let mut moments = [0.0f64; 4];
                for (q, fq) in grid.populations(x, y, z).into_iter().enumerate() {
                    let e = E3[q];
                    moments[0] += f64::from(e.0) * fq;
                    moments[1] += f64::from(e.1) * fq;
                    moments[2] += f64::from(e.2) * fq;
                    moments[3] += fq;
                }
                if moments.iter().all(|m| m.is_finite()) && moments[3] > 0.0 {
                    Ok(moments)
                } else {
                    Err(ChtError::FlowDiverged { step })
                }
            })
            .collect()
    };
    let workers = if config.workers == 0 {
        std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
    } else {
        config.workers
    };
    let pool = TilePool::for_host(workers, 0);
    let (steps, last_change, previous) = pool.with_parked_crew_local(|parked| {
        let mut previous = sample(&grid, 0)?;
        let mut steps = 0usize;
        let mut last_change = f64::INFINITY;
        while steps < config.max_steps {
            poll(gate)?;
            let batch = config.check_interval.min(config.max_steps - steps);
            for _ in 0..batch {
                grid.step_pooled(parked, gate)
                    .map_err(|error| match error {
                        BoundaryStepError3::Cancelled => ChtError::Cancelled,
                        BoundaryStepError3::Collision { .. }
                        | BoundaryStepError3::Unphysical { .. } => {
                            ChtError::FlowDiverged { step: steps }
                        }
                        BoundaryStepError3::Pool(detail) => ChtError::Executor { detail },
                    })?;
                steps += 1;
            }
            let current = sample(&grid, steps)?;
            let (mut diff, mut norm) = (0.0f64, 0.0f64);
            for (a, b) in current.iter().zip(&previous) {
                for k in 0..3 {
                    diff = (a[k] - b[k]).mul_add(a[k] - b[k], diff);
                    norm = a[k].mul_add(a[k], norm);
                }
            }
            last_change = fs_math::det::sqrt(diff / norm.max(f64::MIN_POSITIVE));
            previous = current;
            if last_change <= config.steady_tolerance {
                break;
            }
        }
        Ok::<_, ChtError>((steps, last_change, previous))
    })?;
    if last_change > config.steady_tolerance {
        return Err(ChtError::FlowNotSteady {
            steps,
            last_change,
            tolerance: config.steady_tolerance,
        });
    }
    let mut max_lattice_speed = 0.0f64;
    let mut velocities = vec![[0.0f64; 3]; domain.cell_count()];
    let mut pressure_pa = vec![0.0f64; domain.cell_count()];
    let pressure_scale = crate::CS2 * fluid.density_kg_m3 * velocity_scale_m_s * velocity_scale_m_s;
    for (&c, m) in fluid_cells.iter().zip(&previous) {
        let speed = fs_math::det::sqrt(m[0].mul_add(m[0], m[1].mul_add(m[1], m[2] * m[2])));
        max_lattice_speed = max_lattice_speed.max(speed);
        velocities[c] = [m[0], m[1], m[2]].map(|v| v * velocity_scale_m_s);
        pressure_pa[c] = (m[3] - 1.0) * pressure_scale;
    }
    // Both on-site faces are Free: the interior lattice layers conserve mass
    // exactly, while the boundary cells' momentum is a reconstruction, so the
    // interior field decides the transported flux.
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let (field, projection) = FlowField::from_cell_velocities(
        domain,
        &velocities,
        faces,
        config.projection_tolerance,
        gate,
    )?;
    let open_cells = (0..ny * nz).filter(|&i| domain.is_fluid(i * nx)).count();
    let nominal_inflow_m3_s =
        config.inlet_velocity_m_s * open_cells as f64 * domain.dx() * domain.dx();
    let realized_inflow_m3_s = -field.boundary_outflow(domain, Face3::XMin);
    let report = LbmFlowReport {
        steps,
        nominal_inflow_m3_s,
        realized_inflow_m3_s,
        last_change,
        tau,
        lattice_viscosity,
        lattice_inlet_velocity: u_lat,
        inlet_mach: u_lat * fs_math::det::sqrt(3.0),
        max_lattice_speed,
        cell_reynolds,
        velocity_scale_m_s,
        pressure_drop_pa: f64::NAN,
        collision,
        projection,
    };
    let mut flow = LbmFlow {
        field,
        velocity_m_s: velocities,
        pressure_pa,
        report,
    };
    flow.report.pressure_drop_pa = match (
        flow.mean_pressure_x(domain, 1),
        flow.mean_pressure_x(domain, nx - 2),
    ) {
        (Some(inlet), Some(outlet)) => inlet - outlet,
        _ => f64::NAN,
    };
    Ok(flow)
}
