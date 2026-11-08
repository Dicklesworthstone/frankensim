//! Mass-consistent variable-density transport on the native BDM1/P0 solver.
//!
//! Density uses the original positive donor-cell step. Its face mass flux is
//! reconstructed in H(div) with constant cell divergence and quadratic face
//! moments. The SAME flux then transports the complete affine BDM velocity:
//! `M(rho_new) u_new/dt + C(j) u_new + A(mu) u_new + G p
//!  = M(rho_old) u_old/dt + f`, with the original continuity equations.
//! The advecting volume velocity is explicitly lagged. Momentum, viscosity and
//! pressure are implicit; this is a first-order fixed-mesh IMEX discretization,
//! not a fully implicit density/velocity iteration or an EOS closure.
//!
//! Constant velocity is compatible with the density update because div(j)
//! equals the cell mass loss, up to reported density-conversion roundoff.
//! Quadratic face moments also match |u|², so the
//! conservative kinetic identity retains old-density time dissipation and
//! nonnegative upwind loss, including faces with simultaneous counterflow.
//! Spatial viscous/SIP loss and external work retain their native meanings.

use super::reaction::{ReactionAssembly, velocity};
use super::transport::Sums;
use super::{
    AdmittedBoundary, OpenBoundaryPolicy, StokesAssembly, StokesBoundaryError, StokesSolution,
};
use crate::bdm::{eval_basis, tri_quad};
use crate::material::{
    MaterialFluxNumerics, MaterialNavierStokesBoundary, with_kinematic_boundaries,
};
use crate::ns::transport::density::{
    DensityTransportCell, DensityTransportError, DensityTransportReport,
};
use crate::ns::transport::mass_flux::PreparedDensityMassFlux;
use crate::ns::{FluxParams, FluxSystem, convection_boundary_trace, momentum_entry};
use fs_sparse::Coo;

/// Bounded density subcycling with one original native macro momentum solve.
pub mod subcycle;

/// Conservative SSPRK2 density stages with one original BDF1 momentum solve.
pub mod ssprk2;

/// Implicit SGS and molecular rheology on the same transported-density step.
#[cfg(feature = "experimental-sgs")]
pub mod sgs;

/// Physical normalization and explicit density-transport allowances.
#[derive(Clone, Copy, Debug)]
pub struct DensityMomentumControls {
    /// One positive reference density [kg/m³] for conditioning the entire
    /// momentum/stress equation. This is never a replacement cell density.
    pub reference_density_kg_m3: f64,
    /// Positive out-of-plane depth [m] for physical inventory observations.
    pub depth_m: f64,
    /// Outgoing-volume Courant cap, strictly between zero and one.
    /// Bounded linear density applies one third of this sufficient allowance.
    pub maximum_courant: f64,
    /// Admitted relative cell volume divergence, from zero through 1e-12.
    pub relative_volume_divergence_tolerance: f64,
    /// Per-owner cap for the density plan, its workspace, the mass-flux
    /// reconstruction, and this step's transport buffers. It does not cap the
    /// existing native sparse assembly or dense/sparse linear solver storage.
    pub maximum_transport_storage_bytes: usize,
}

/// Refusal before a coupled density/velocity candidate is accepted.
#[derive(Debug, Clone, PartialEq)]
pub enum DensityMomentumError {
    /// Original conservative density or reconstructed mass-flux refusal.
    Transport(DensityTransportError),
    /// Original native boundary, viscosity, pressure or linear-solve refusal.
    Flow(StokesBoundaryError),
    /// Invalid units, dimensions, controls or unrepresentable physical output.
    Invalid(&'static str),
    /// A bounded transport allocation could not be reserved.
    Allocation,
    /// Cancellation before or after the original synchronous solver.
    Cancelled,
}

impl core::fmt::Display for DensityMomentumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "conservative density/momentum step: {self:?}")
    }
}
impl std::error::Error for DensityMomentumError {}
impl From<DensityTransportError> for DensityMomentumError {
    fn from(value: DensityTransportError) -> Self {
        Self::Transport(value)
    }
}
impl From<StokesBoundaryError> for DensityMomentumError {
    fn from(value: StokesBoundaryError) -> Self {
        Self::Flow(value)
    }
}
type Result<T> = core::result::Result<T, DensityMomentumError>;

/// Actual physical momentum and kinetic storage of the native BDM fields.
/// Numerical losses have no automatic thermal destination.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityMomentumObservations {
    /// Domain momentum at the old and new density/velocity states [kg m/s].
    pub momentum_before_kg_m_s: [f64; 2],
    pub momentum_after_kg_m_s: [f64; 2],
    /// Consistent kinetic storage, including subcell affine velocity [J].
    pub kinetic_energy_before_j: f64,
    pub kinetic_energy_after_j: f64,
    /// `0.5 * integral rho_old |u_new-u_old|² dV` [J].
    pub temporal_dissipation_j: f64,
    /// `0.5 * dt * depth * sum_faces integral |j.n| |jump u_new|²` [J].
    /// Boundary upwind loss is retained in the native boundary transport view.
    pub interior_upwind_dissipation_j: f64,
}

/// Accepted physical inventory and complete native pressure/velocity field.
pub struct DensityMomentumStep<const N: usize> {
    /// Native BDM velocity and P0 pressure, retaining boundary reactions.
    /// Pressure/stress/momentum diagnostics in this native view are divided by
    /// `reference_density_kg_m3`; volume rates and velocities are physical.
    pub flow: StokesSolution,
    /// Mechanical pressure [Pa] in the original component-local pressure datum.
    /// This is not an absolute thermodynamic/EOS pressure.
    pub pressure_pa: Vec<f64>,
    /// Accepted positive P0 density and transported specific inventories.
    pub cells: Vec<DensityTransportCell<N>>,
    /// Same accepted donor-cell mass and carried-inventory accounting.
    pub transport: DensityTransportReport<N>,
    /// Consistent physical inertia and separated numerical dissipation.
    pub observations: DensityMomentumObservations,
    /// The single conditioning reference used for native normalized outputs.
    pub reference_density_kg_m3: f64,
}

/// Coupled candidate bound to immutable old cell and velocity states.
/// The original matrix/RHS remain available for native linear transpose solves;
/// this is not a Jacobian of the complete density/velocity transport step.
pub struct DensityMomentumAssembly<'a, 'm, const N: usize> {
    momentum: StokesAssembly<'a, 'm>,
    cells: Vec<DensityTransportCell<N>>,
    transport: DensityTransportReport<N>,
    old_cells: &'a [DensityTransportCell<N>],
    old_velocity: &'a [f64],
    faces: Vec<RelativeMassFluxFace>,
    controls: DensityMomentumControls,
}

// Transport, physical force and boundary data are prepared once. Only the
// viscous stress depends on the candidate in a nonlinear constitutive solve.
struct PreparedDensityMomentum<'a, 'm, const N: usize> {
    system: &'a FluxSystem<'m>,
    params: FluxParams,
    admitted: AdmittedBoundary,
    inv_dt: f64,
    current_weights: Vec<[f64; 6]>,
    conservative: ConservativeMomentum,
    cells: Vec<DensityTransportCell<N>>,
    transport: DensityTransportReport<N>,
    old_cells: &'a [DensityTransportCell<N>],
    old_velocity: &'a [f64],
    controls: DensityMomentumControls,
}

impl<'a, 'm, const N: usize> PreparedDensityMomentum<'a, 'm, N> {
    fn assemble(
        &self,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        current_weights: Vec<[f64; 6]>,
    ) -> core::result::Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        self.system.assemble_stokes_admitted_conservative(
            self.params,
            &|cell, point| dynamic_viscosity(cell, point) / self.controls.reference_density_kg_m3,
            &|_| [0.0; 2],
            &self.admitted,
            None,
            self.inv_dt,
            None,
            None,
            None,
            Some(current_weights),
            Some(&self.conservative),
        )
    }

    fn into_assembly(self, momentum: StokesAssembly<'a, 'm>) -> DensityMomentumAssembly<'a, 'm, N> {
        DensityMomentumAssembly {
            momentum,
            cells: self.cells,
            transport: self.transport,
            old_cells: self.old_cells,
            old_velocity: self.old_velocity,
            faces: self.conservative.faces,
            controls: self.controls,
        }
    }
}

impl<'a, 'm, const N: usize> DensityMomentumAssembly<'a, 'm, N> {
    /// Native mass, conservative momentum, viscosity and pressure operator.
    pub fn momentum(&self) -> &StokesAssembly<'a, 'm> {
        &self.momentum
    }

    /// Transported candidate cells; the coupled velocity solve is still pending.
    pub fn transported_cells(&self) -> &[DensityTransportCell<N>] {
        &self.cells
    }

    /// Solve and observe the complete native field before publishing any state.
    /// Cancellation is checked before and after the existing synchronous solve,
    /// and per cell/edge in the observations; no mid-solve cancellation is claimed.
    pub fn solve(&self, mut cancelled: impl FnMut() -> bool) -> Result<DensityMomentumStep<N>> {
        poll(&mut cancelled)?;
        let flow = self.momentum.solve()?;
        observe_step(
            self.momentum.system,
            flow,
            &self.cells,
            self.transport,
            self.old_cells,
            self.old_velocity,
            &self.faces,
            self.controls,
            cancelled,
        )
    }
}

// Linear and nonlinear endpoints share the same physical inventory and
// numerical-dissipation observations. This never queries a constitutive law.
#[allow(clippy::too_many_arguments)]
fn observe_step<const N: usize>(
    system: &FluxSystem<'_>,
    flow: StokesSolution,
    cells: &[DensityTransportCell<N>],
    transport: DensityTransportReport<N>,
    old_cells: &[DensityTransportCell<N>],
    old_velocity: &[f64],
    faces: &[RelativeMassFluxFace],
    controls: DensityMomentumControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DensityMomentumStep<N>> {
    poll(&mut cancelled)?;
    let before = observe(
        system,
        old_cells,
        old_velocity,
        controls.depth_m,
        &mut cancelled,
    )?;
    let after = observe(
        system,
        cells,
        &flow.field.x,
        controls.depth_m,
        &mut cancelled,
    )?;
    let mut temporal = Sums::<1>::new();
    for (cell, nodes) in system.mesh.tris.iter().enumerate() {
        poll(&mut cancelled)?;
        let points = nodes.map(|node| system.mesh.verts[node]);
        for (point, weight) in tri_quad(points, system.mesh.areas[cell]) {
            let old = velocity(system, cell, point, old_velocity);
            let new = velocity(system, cell, point, &flow.field.x);
            let mass = weight * old_cells[cell].density_kg_m3 * controls.depth_m;
            temporal.add([half_square(mass, new[0] - old[0]) + half_square(mass, new[1] - old[1])]);
        }
    }
    let mut upwind = Sums::<1>::new();
    for (index, edge) in system.mesh.edges.iter().enumerate() {
        poll(&mut cancelled)?;
        let (a, b) = edge.tris;
        if b == usize::MAX {
            continue;
        }
        for sample in faces[index].points() {
            let left = velocity(system, a, sample.point_m, &flow.field.x);
            let right = velocity(system, b, sample.point_m, &flow.field.x);
            let mass = sample.mass_rate_kg_m_s.abs() * transport.timestep_s * controls.depth_m;
            upwind.add([
                half_square(mass, left[0] - right[0]) + half_square(mass, left[1] - right[1])
            ]);
        }
    }
    let temporal = temporal.values()[0];
    let upwind = upwind.values()[0];
    if !temporal.is_finite() || temporal < 0.0 || !upwind.is_finite() || upwind < 0.0 {
        return Err(invalid("unrepresentable numerical kinetic dissipation"));
    }
    let mut pressure_pa = allocate(flow.pressure.len(), 0.0)?;
    for (target, &pressure) in pressure_pa.iter_mut().zip(&flow.pressure) {
        *target = pressure * controls.reference_density_kg_m3;
        if !target.is_finite() {
            return Err(invalid("unrepresentable physical pressure"));
        }
    }
    let mut accepted_cells = allocate(cells.len(), DensityTransportCell::default())?;
    accepted_cells.copy_from_slice(cells);
    poll(&mut cancelled)?;
    Ok(DensityMomentumStep {
        flow,
        pressure_pa,
        cells: accepted_cells,
        transport,
        observations: DensityMomentumObservations {
            momentum_before_kg_m_s: [before[0], before[1]],
            momentum_after_kg_m_s: [after[0], after[1]],
            kinetic_energy_before_j: before[2],
            kinetic_energy_after_j: after[2],
            temporal_dissipation_j: temporal,
            interior_upwind_dissipation_j: upwind,
        },
        reference_density_kg_m3: controls.reference_density_kg_m3,
    })
}

impl<'m> FluxSystem<'m> {
    /// Prepare first-order variable-density incompressible momentum on the
    /// original BDM1/P0 mesh. Density and N mass-specific quantities use the
    /// same selected directional mass transfers. The N slots are passive specific
    /// quantities (e.g. enthalpy or species); forces act on the native velocity,
    /// not on passive slots that happen to contain a velocity-like quantity.
    /// Native momentum retains ALL affine
    /// velocity modes and the consistent mass, rather than projecting cell
    /// average momentum into a different velocity representation.
    ///
    /// `advector` is an explicitly frozen divergence-free volume velocity [m/s]
    /// in native BDM moments; normally supply the previous accepted velocity.
    /// `previous_velocity` is the immutable momentum history. Boundary inflow
    /// densities/specifics follow original edge order. Momentum inflow follows
    /// the existing physical boundary owner (velocity, slip, literal traction,
    /// reservoir velocity or explicitly dissipative backflow).
    ///
    /// `dynamic_viscosity(cell, point, candidate_cell)` returns mu [Pa s], and
    /// `force_density(...)` returns physical force density [N/m³]. Both are
    /// frozen endpoint laws; no EOS, temperature inversion, diffusion, phase
    /// source or nonlinear rheology is inferred from the transported specifics.
    /// All stresses, pressure and inertia share ONE conditioning reference.
    ///
    /// For homogeneous closed velocity data, the conservative temporal/advective
    /// kinetic identity has old-density time loss and nonnegative upwind loss.
    /// A complete energy balance also includes the original viscous/SIP form,
    /// source work, boundary kinetic flux and any explicit backflow resistance.
    /// These numerical losses are not silently deposited as physical heat.
    ///
    /// # Errors
    /// Original density/CFL/boundary/viscosity refusals, invalid physical data,
    /// exhausted transport storage or cancellation. No caller state is mutated.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn assemble_density_momentum_bdf1<'a, const N: usize>(
        &'a self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &'a [f64],
        cells: &'a [DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumAssembly<'a, 'm, N>> {
        let mut prepared = self.prepare_density_momentum_bdf1(
            numerics,
            force_density,
            boundaries,
            advector,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            &mut cancelled,
        )?;
        let weights = core::mem::take(&mut prepared.current_weights);
        let momentum = prepared.assemble(
            &|cell, point| dynamic_viscosity(cell, point, &prepared.cells[cell]),
            weights,
        )?;
        poll(&mut cancelled)?;
        Ok(prepared.into_assembly(momentum))
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_density_momentum_bdf1<'a, const N: usize>(
        &'a self,
        numerics: MaterialFluxNumerics,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &'a [f64],
        cells: &'a [DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDensityMomentum<'a, 'm, N>> {
        poll(&mut cancelled)?;
        if !positive(controls.reference_density_kg_m3) || !positive(controls.depth_m) {
            return Err(invalid(
                "positive physical reference density and depth required",
            ));
        }
        let inv_dt = self.admit_mixed_bdf1_history(previous_velocity, timestep_s)?;
        let params = FluxParams {
            // Unused: dynamic viscosity is sampled and normalized explicitly.
            nu: 1.0,
            sigma: numerics.penalty,
            tol: numerics.linear_tolerance,
            max_iters: numerics.max_linear_iterations,
        };
        let nc = self.mesh.tris.len();
        let ne = self.mesh.edges.len();
        admit_transport_buffers::<N>(nc, ne, controls.maximum_transport_storage_bytes)?;
        let plan = self.prepare_density_transport(
            advector,
            controls.depth_m,
            controls.relative_volume_divergence_tolerance,
            controls.maximum_transport_storage_bytes,
            &mut cancelled,
        )?;
        let mut scratch = plan.workspace::<N>(controls.maximum_transport_storage_bytes)?;
        let mut next = allocate(nc, DensityTransportCell::default())?;
        let transport = plan.advance_into(
            cells,
            inflow,
            timestep_s,
            controls.maximum_courant,
            &mut next,
            &mut scratch,
            &mut cancelled,
        )?;
        drop(scratch);
        let mass_flux = plan.prepare_mass_flux(
            cells,
            inflow,
            controls.maximum_transport_storage_bytes,
            &mut cancelled,
        )?;
        drop(plan);
        self.prepare_transported_density_momentum_bdf1(
            params,
            force_density,
            boundaries,
            previous_velocity,
            cells,
            next,
            transport,
            mass_flux,
            inv_dt,
            controls,
            &mut cancelled,
        )
    }

    // The single-step and density-subcycled paths share this exact original
    // momentum owner. No density stage may rebuild or solve it independently.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn assemble_transported_density_momentum_bdf1<'a, const N: usize>(
        &'a self,
        params: FluxParams,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &'a [f64],
        cells: &'a [DensityTransportCell<N>],
        next: Vec<DensityTransportCell<N>>,
        transport: DensityTransportReport<N>,
        mass_flux: PreparedDensityMassFlux<'m>,
        inv_dt: f64,
        controls: DensityMomentumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumAssembly<'a, 'm, N>> {
        let mut prepared = self.prepare_transported_density_momentum_bdf1(
            params,
            force_density,
            boundaries,
            previous_velocity,
            cells,
            next,
            transport,
            mass_flux,
            inv_dt,
            controls,
            &mut cancelled,
        )?;
        let weights = core::mem::take(&mut prepared.current_weights);
        let momentum = prepared.assemble(
            &|cell, point| dynamic_viscosity(cell, point, &prepared.cells[cell]),
            weights,
        )?;
        poll(&mut cancelled)?;
        Ok(prepared.into_assembly(momentum))
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn prepare_transported_density_momentum_bdf1<'a, const N: usize>(
        &'a self,
        params: FluxParams,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &'a [f64],
        cells: &'a [DensityTransportCell<N>],
        next: Vec<DensityTransportCell<N>>,
        transport: DensityTransportReport<N>,
        mass_flux: PreparedDensityMassFlux<'m>,
        inv_dt: f64,
        controls: DensityMomentumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDensityMomentum<'a, 'm, N>> {
        let nc = self.mesh.tris.len();
        let ne = self.mesh.edges.len();
        let admitted =
            with_kinematic_boundaries(controls.reference_density_kg_m3, boundaries, |native| {
                self.admit_navier_stokes_boundaries(params, native)
            })?;
        let reference = controls.reference_density_kg_m3;
        let mut current_weights = allocate(nc, [0.0; 6])?;
        let mut volume = allocate(nc, [[0.0; 2]; 6])?;
        let mut loads = allocate(nc, [0.0; 6])?;
        for (cell, nodes) in self.mesh.tris.iter().enumerate() {
            poll(&mut cancelled)?;
            let old_ratio = positive_ratio(cells[cell].density_kg_m3, reference)?;
            let new_ratio = positive_ratio(next[cell].density_kg_m3, reference)?;
            let points = nodes.map(|node| self.mesh.verts[node]);
            for (q, (point, weight)) in tri_quad(points, self.mesh.areas[cell])
                .into_iter()
                .enumerate()
            {
                current_weights[cell][q] = weight * new_ratio;
                if !positive(current_weights[cell][q]) {
                    return Err(invalid("unrepresentable current mass quadrature"));
                }
                let j = mass_flux.cell_value(cell, point)?;
                volume[cell][q] = [ratio(j[0], reference)?, ratio(j[1], reference)?];
                let old = velocity(self, cell, point, previous_velocity);
                if old.iter().any(|value| !value.is_finite()) {
                    return Err(invalid("unrepresentable old velocity trace"));
                }
                let old_weight = weight * old_ratio;
                if !positive(old_weight) {
                    return Err(invalid("unrepresentable old mass quadrature"));
                }
                let force = force_density(cell, point, &next[cell]);
                if force.iter().any(|value| !value.is_finite()) {
                    return Err(StokesBoundaryError::NonFiniteForce.into());
                }
                let force = [ratio(force[0], reference)?, ratio(force[1], reference)?];
                for i in 0..6 {
                    let value = eval_basis(&self.bases[cell], i, point);
                    // Form the OLD mass history directly. Subtracting a new-mass
                    // history correction would erase a light donor at high contrast.
                    loads[cell][i] += inv_dt * old_weight * (value[0] * old[0] + value[1] * old[1]);
                    loads[cell][i] += weight * (value[0] * force[0] + value[1] * force[1]);
                }
            }
            if loads[cell].iter().any(|value| !value.is_finite()) {
                return Err(invalid("unrepresentable physical momentum history or load"));
            }
        }
        let mut faces = allocate(ne, RelativeMassFluxFace::default())?;
        for (edge, face) in faces.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            let points = mass_flux
                .edge_quadrature(edge)
                .ok_or_else(|| invalid("missing mass-flux edge"))?;
            if points.len() > face.samples.len() {
                return Err(invalid("mass-flux quadrature exceeds native face capacity"));
            }
            face.count = points.len();
            for (target, source) in face.samples.iter_mut().zip(points) {
                *target = RelativeMassFluxPoint {
                    point_m: source.point_m,
                    volume_rate_m2_s: source.volume_rate_m2_s,
                    mass_rate_kg_m_s: source.mass_rate_kg_m_s,
                    mass_rate_per_reference: ratio(source.mass_rate_kg_m_s, reference)?,
                    boundary_fraction: source.boundary_fraction,
                };
            }
        }
        drop(mass_flux);
        let conservative = ConservativeMomentum {
            volume,
            loads,
            faces,
        };
        poll(&mut cancelled)?;
        Ok(PreparedDensityMomentum {
            system: self,
            params,
            admitted,
            inv_dt,
            current_weights,
            conservative,
            cells: next,
            transport,
            old_cells: cells,
            old_velocity: previous_velocity,
            controls,
        })
    }

    /// Prepare and solve [`Self::assemble_density_momentum_bdf1`] atomically.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_bdf1_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumStep<N>> {
        self.assemble_density_momentum_bdf1(
            numerics,
            dynamic_viscosity,
            force_density,
            boundaries,
            advector,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            &mut cancelled,
        )?
        .solve(&mut cancelled)
    }
}

fn admit_transport_buffers<const N: usize>(nc: usize, ne: usize, maximum: usize) -> Result<()> {
    let bytes = nc
        .checked_mul(
            size_of::<DensityTransportCell<N>>()
                + size_of::<[[f64; 2]; 6]>()
                + 2 * size_of::<[f64; 6]>(),
        )
        .and_then(|value| {
            ne.checked_mul(size_of::<RelativeMassFluxFace>())
                .and_then(|n| value.checked_add(n))
        })
        .ok_or_else(|| invalid("transport buffer size overflow"))?;
    if bytes > maximum {
        Err(DensityTransportError::Storage {
            required: bytes,
            maximum,
        }
        .into())
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
pub(super) struct RelativeMassFluxPoint {
    pub(super) point_m: [f64; 2],
    pub(super) volume_rate_m2_s: f64,
    pub(super) mass_rate_kg_m_s: f64,
    pub(super) mass_rate_per_reference: f64,
    pub(super) boundary_fraction: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct RelativeMassFluxFace {
    samples: [RelativeMassFluxPoint; crate::ns::CONVECTION_EDGE_QUADRATURE_CAPACITY],
    count: usize,
}
impl RelativeMassFluxFace {
    pub(super) fn points(&self) -> &[RelativeMassFluxPoint] {
        &self.samples[..self.count]
    }
}

/// Prevalidated normalized transport; used before the existing native owner
/// captures reaction loads and replaces strong normal-velocity rows.
pub(in crate::ns) struct ConservativeMomentum {
    volume: Vec<[[f64; 2]; 6]>,
    loads: Vec<[f64; 6]>,
    pub(super) faces: Vec<RelativeMassFluxFace>,
}

impl ConservativeMomentum {
    pub(in crate::ns) fn add_operator(
        &self,
        system: &FluxSystem<'_>,
        matrix: &mut Coo,
        rhs: &mut [f64],
        constraints: &[Option<f64>],
        boundary: &AdmittedBoundary,
        mut reactions: Option<&mut ReactionAssembly>,
    ) {
        for (cell, nodes) in system.mesh.tris.iter().enumerate() {
            let points = nodes.map(|node| system.mesh.verts[node]);
            for i in 0..6 {
                rhs[system.dof(cell, i)] += self.loads[cell][i];
            }
            for (q, (point, weight)) in tri_quad(points, system.mesh.areas[cell])
                .into_iter()
                .enumerate()
            {
                let flux = self.volume[cell][q];
                let values: [[f64; 2]; 6] =
                    core::array::from_fn(|i| eval_basis(&system.bases[cell], i, point));
                for i in 0..6 {
                    let gradient = system.bases[cell].grad[i];
                    let derivative = [
                        flux[0] * gradient[0][0] + flux[1] * gradient[0][1],
                        flux[0] * gradient[1][0] + flux[1] * gradient[1][1],
                    ];
                    for (j, value) in values.iter().enumerate() {
                        momentum_entry(
                            matrix,
                            &mut reactions,
                            constraints,
                            system.dof(cell, i),
                            system.dof(cell, j),
                            -weight * (value[0] * derivative[0] + value[1] * derivative[1]),
                        );
                    }
                }
            }
        }
        for (index, edge) in system.mesh.edges.iter().enumerate() {
            let (left, right) = edge.tris;
            for sample in self.faces[index].points() {
                let rate = sample.mass_rate_per_reference;
                let point = sample.point_m;
                if right == usize::MAX {
                    let reservoir = boundary.reservoir_velocity(index);
                    let slip = boundary.tangential_resistance(index).is_some();
                    let policy = boundary.open_policy(index);
                    let prescribed = if rate < 0.0 && !slip {
                        reservoir.or_else(|| {
                            policy
                                .is_none()
                                .then(|| [boundary.velocity(index, 0), boundary.velocity(index, 1)])
                        })
                    } else {
                        None
                    };
                    for i in 0..6 {
                        let vi = eval_basis(&system.bases[left], i, point);
                        let row = system.dof(left, i);
                        if let Some(samples) = prescribed {
                            let value =
                                convection_boundary_trace(samples, sample.boundary_fraction);
                            rhs[row] -= rate * (value[0] * vi[0] + value[1] * vi[1]);
                        } else {
                            let rate = if policy == Some(OpenBoundaryPolicy::DissipativeBackflow) {
                                rate.max(0.0)
                            } else {
                                rate
                            };
                            for j in 0..6 {
                                let uj = eval_basis(&system.bases[left], j, point);
                                momentum_entry(
                                    matrix,
                                    &mut reactions,
                                    constraints,
                                    row,
                                    system.dof(left, j),
                                    rate * (vi[0] * uj[0] + vi[1] * uj[1]),
                                );
                            }
                        }
                    }
                    continue;
                }
                let donor = if rate >= 0.0 { left } else { right };
                for (cell, sign) in [(left, 1.0), (right, -1.0)] {
                    for i in 0..6 {
                        let vi = eval_basis(&system.bases[cell], i, point);
                        for j in 0..6 {
                            let uj = eval_basis(&system.bases[donor], j, point);
                            momentum_entry(
                                matrix,
                                &mut reactions,
                                constraints,
                                system.dof(cell, i),
                                system.dof(donor, j),
                                sign * rate * (vi[0] * uj[0] + vi[1] * uj[1]),
                            );
                        }
                    }
                }
            }
        }
    }
}

fn observe<const N: usize>(
    system: &FluxSystem<'_>,
    cells: &[DensityTransportCell<N>],
    field: &[f64],
    depth: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3]> {
    let mut sums = Sums::<3>::new();
    for (cell, nodes) in system.mesh.tris.iter().enumerate() {
        poll(cancelled)?;
        let points = nodes.map(|node| system.mesh.verts[node]);
        for (point, weight) in tri_quad(points, system.mesh.areas[cell]) {
            let u = velocity(system, cell, point, field);
            let mass = weight * cells[cell].density_kg_m3 * depth;
            sums.add([
                mass * u[0],
                mass * u[1],
                half_square(mass, u[0]) + half_square(mass, u[1]),
            ]);
        }
    }
    let result = sums.values();
    if result.iter().any(|value| !value.is_finite()) || result[2] < 0.0 {
        return Err(invalid(
            "unrepresentable physical momentum or kinetic storage",
        ));
    }
    Ok(result)
}

fn half_square(weight: f64, value: f64) -> f64 {
    (0.5 * value * weight) * value
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn ratio(value: f64, reference: f64) -> Result<f64> {
    let result = value / reference;
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        return Err(invalid("unrepresentable common-reference normalization"));
    }
    Ok(result)
}
fn positive_ratio(value: f64, reference: f64) -> Result<f64> {
    let result = ratio(value, reference)?;
    if !positive(result) {
        return Err(invalid("unrepresentable positive relative density"));
    }
    Ok(result)
}
fn allocate<T: Clone>(length: usize, value: T) -> Result<Vec<T>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(length)
        .map_err(|_| DensityMomentumError::Allocation)?;
    result.resize(length, value);
    Ok(result)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(DensityMomentumError::Cancelled)
    } else {
        Ok(())
    }
}
fn invalid(reason: &'static str) -> DensityMomentumError {
    DensityMomentumError::Invalid(reason)
}
