//! One accepted density transport and conservative native IMEX momentum step.
//!
//! The original old velocity supplies BOTH mass and momentum transport. Old
//! and new density weight their own consistent inertia; pressure, viscous
//! stress and velocity constraints remain with the native mixed operator.

use super::{MaterialStokesBoundary, with_physical_boundaries};
use crate::bdm::{eval_basis, tri_quad};
use crate::ns::boundary::{StokesAssembly, StokesBoundaryError, StokesPressureReference};
use crate::ns::material::MaterialFluxNumerics;
use crate::ns::transport::density::{
    DensityTransportCell, DensityTransportError, DensityTransportReport,
};
use crate::ns::transport::mass_flux::PreparedDensityMassFlux;
use crate::ns::{FluxParams, FluxSolution, FluxSystem};
use fs_math::dd::Dd;
use fs_sparse::Csr;

/// Immutable history and exterior material on the original mesh.
pub struct ConservativeMomentumHistory<'a> {
    /// One physical cell-average density [kg/m³] per original triangle.
    pub density: &'a [DensityTransportCell<0>],
    /// Original divergence-free BDM1 velocity moments [m/s]. This same field
    /// advects density and supplies the explicitly transported momentum.
    pub velocity: &'a [f64],
    /// Original edge order: density on each incoming boundary stream; None
    /// on interior edges. The original density owner performs admission.
    pub inflow_density: &'a [Option<DensityTransportCell<0>>],
    /// Affine exterior velocity at original boundary Gauss sites [m/s].
    /// Needed on incoming traction/Navier-slip boundaries. Full Velocity
    /// boundaries use their ONE admitted velocity trace, so these entries
    /// must be None there. Interior entries must also be None.
    pub inflow_velocity: &'a [Option<[[f64; 2]; 2]>],
}

/// Physical time/depth and explicit transport admission, separate from viscosity.
#[derive(Clone, Copy, Debug)]
pub struct ConservativeMomentumControls {
    pub timestep_s: f64,
    /// Fixed out-of-plane depth for the density owner's total mass [m].
    pub depth_m: f64,
    /// One positive common normalization for inertia, viscosity, pressure and
    /// force. This is NOT a replacement for the evolving local density.
    pub reference_density_kg_m3: f64,
    /// Strictly between zero and one; applies to the original outgoing-volume
    /// density CFL. It does not certify BDM momentum time-step stability.
    pub maximum_courant: f64,
    /// Cap for EACH prepared transport owner or scratch vector bundle, as in
    /// their original APIs. Not a total step cap: native sparse assembly,
    /// solver storage and this wrapper's cell/velocity arrays are additional.
    pub maximum_transport_storage_bytes: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ConservativeMomentumError {
    Transport(DensityTransportError),
    Flow(StokesBoundaryError),
}

impl core::fmt::Display for ConservativeMomentumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "conservative density/momentum: {self:?}")
    }
}
impl std::error::Error for ConservativeMomentumError {}
impl From<DensityTransportError> for ConservativeMomentumError {
    fn from(value: DensityTransportError) -> Self {
        Self::Transport(value)
    }
}
impl From<StokesBoundaryError> for ConservativeMomentumError {
    fn from(value: StokesBoundaryError) -> Self {
        Self::Flow(value)
    }
}
type Result<T> = core::result::Result<T, ConservativeMomentumError>;

/// Prepared native saddle operator with one immutable density/transport pair.
/// Neither the supplied density nor velocity history is changed by preparation.
pub struct ConservativeMomentumAssembly<'s, 'm> {
    system: &'s FluxSystem<'m>,
    native: StokesAssembly<'s, 'm>,
    rhs: Vec<f64>,
    params: FluxParams,
    density: Vec<DensityTransportCell<0>>,
    mass_transport: DensityTransportReport<0>,
    mass_flux: PreparedDensityMassFlux<'m>,
    reference_density: f64,
}

/// Joint accepted density, velocity and pressure from the same native step.
/// The retained mass flux is the old-state transport used by BOTH equations.
pub struct ConservativeMomentumStep<'m> {
    field: FluxSolution,
    density: Vec<DensityTransportCell<0>>,
    pressure_pa: Vec<f64>,
    pressure_reference: Vec<StokesPressureReference>,
    kinetic_energy_j_per_m: f64,
    mass_transport: DensityTransportReport<0>,
    mass_flux: PreparedDensityMassFlux<'m>,
}

impl<'m> FluxSystem<'m> {
    /// Prepare one conservative incompressible IMEX BDF1 step in coherent SI:
    /// `(rho_new u_new-rho_old u_old)/dt + C(m_old)u_old
    ///  - div(2 mu D(u_new)) + grad(p_new) = force_density`.
    ///
    /// The method itself prepares the original divergence-free BDM volume
    /// flow, advances density and reconstructs mass flux from the SAME old
    /// cell/inlet densities. Callers cannot supply an independently edited
    /// new density or a different momentum mass flux. The original degree-four
    /// quadrature weights old and new inertia separately; no cell-average
    /// velocity projection or nonconservative `rho_new (u_new-u_old)` is used.
    /// Old velocity must pass the original 1e-12 relative volume-divergence
    /// allowance. An insufficiently converged history is refused, not projected.
    ///
    /// `dynamic_viscosity(cell, point)` is positive mu [Pa s], force is per
    /// volume [N/m³], and the existing physical boundary type supplies velocity
    /// [m/s], total traction [Pa] and slip resistance [Pa s/m]. One fixed
    /// reference density normalizes ALL equations. Callbacks must be pure and
    /// frozen during preparation. Full velocity boundaries supply the exact
    /// same admitted samples to momentum inflow and the native constraints.
    ///
    /// Explicit convection is first order in time; viscosity/pressure are
    /// implicit. This density CFL preserves positive density but is not a
    /// kinetic-energy stability certificate for continuous BDM momentum. The
    /// original synchronous assembly/solve cancellation limitations remain;
    /// wrapper loops poll at cells/edges and before publication. No EOS,
    /// composition/enthalpy, moving mesh, support-reaction or physical-heating
    /// output is added by this entry. Total traction fixes the dynamic pressure
    /// datum; it does not determine thermodynamic absolute pressure.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_conservative_momentum_bdf1<'s>(
        &'s self,
        history: ConservativeMomentumHistory<'_>,
        controls: ConservativeMomentumControls,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConservativeMomentumAssembly<'s, 'm>> {
        poll(&mut cancelled)?;
        let reference = controls.reference_density_kg_m3;
        if !reference.is_finite()
            || reference <= 0.0
            || history.density.len() != self.mesh.tris.len()
            || history.inflow_density.len() != self.mesh.edges.len()
            || history.inflow_velocity.len() != self.mesh.edges.len()
        {
            return Err(invalid(
                "positive reference density and original density/inflow dimensions required",
            ));
        }
        let inv_dt = self.admit_mixed_bdf1_history(history.velocity, controls.timestep_s)?;
        let params = FluxParams {
            nu: 1.0, // unused: original variable-viscosity operator owns mu/rho_ref
            sigma: numerics.penalty,
            tol: numerics.linear_tolerance,
            max_iters: numerics.max_linear_iterations,
        };
        let anchors = vec![true; self.pressure_component_count()];
        let admitted = with_physical_boundaries(boundaries, reference, |native| {
            self.admit_stokes_problem_with_velocity_anchors(params, native, &anchors)
        })?;
        poll(&mut cancelled)?;
        let cap = controls.maximum_transport_storage_bytes;
        let plan = self.prepare_density_transport(
            history.velocity,
            controls.depth_m,
            1e-12,
            cap,
            &mut cancelled,
        )?;
        let mass_flux =
            plan.prepare_mass_flux(history.density, history.inflow_density, cap, &mut cancelled)?;
        let mut density = vec![DensityTransportCell::default(); history.density.len()];
        let mut density_work = plan.workspace::<0>(cap)?;
        let mass_transport = plan.advance_into(
            history.density,
            history.inflow_density,
            controls.timestep_s,
            controls.maximum_courant,
            &mut density,
            &mut density_work,
            &mut cancelled,
        )?;
        let mut exterior = Vec::with_capacity(self.mesh.edges.len());
        for (edge, geometry) in self.mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let prescribed = geometry.tris.1 == usize::MAX
                && !admitted.is_traction(edge)
                && admitted.tangential_resistance(edge).is_none();
            if prescribed {
                if history.inflow_velocity[edge].is_some() {
                    return Err(invalid(
                        "full velocity boundaries have one admitted inflow velocity owner",
                    ));
                }
                exterior.push(Some([
                    admitted.velocity(edge, 0),
                    admitted.velocity(edge, 1),
                ]));
            } else {
                exterior.push(history.inflow_velocity[edge]);
            }
        }
        let mut convection = vec![0.0; self.n_u];
        let mut momentum_work = mass_flux.momentum_workspace(cap)?;
        mass_flux.apply_momentum_into(
            self,
            history.velocity,
            &exterior,
            &mut convection,
            &mut momentum_work,
            &mut cancelled,
        )?;
        // Explicit transport and old momentum are integrated independently of
        // the new metric, avoiding cancellation of large new/old-density loads.
        let mut load = Vec::with_capacity(self.n_u);
        for chunk in convection.chunks(1024) {
            poll(&mut cancelled)?;
            for &value in chunk {
                load.push(Dd::from_f64(normalize(-value, reference)?));
            }
        }
        let mut weights = Vec::with_capacity(density.len());
        for (cell, nodes) in self.mesh.tris.iter().enumerate() {
            poll(&mut cancelled)?;
            let old = normalize(history.density[cell].density_kg_m3, reference)?;
            let new = normalize(density[cell].density_kg_m3, reference)?;
            let points = nodes.map(|node| self.mesh.verts[node]);
            let mut current_weights = [0.0; 6];
            for (q, (point, weight)) in tri_quad(points, self.mesh.areas[cell])
                .into_iter()
                .enumerate()
            {
                current_weights[q] = weight * new;
                let old_weight = weight * old * inv_dt;
                if !current_weights[q].is_finite()
                    || current_weights[q] <= 0.0
                    || !old_weight.is_finite()
                    || old_weight <= 0.0
                {
                    return Err(invalid("old/new density inertia is not representable"));
                }
                let values: [[f64; 2]; 6] =
                    core::array::from_fn(|i| eval_basis(&self.bases[cell], i, point));
                let mut old_velocity = [0.0; 2];
                for (i, value) in values.iter().enumerate() {
                    let coefficient = history.velocity[self.dof(cell, i)];
                    old_velocity[0] += coefficient * value[0];
                    old_velocity[1] += coefficient * value[1];
                }
                for (i, test) in values.iter().enumerate() {
                    let value =
                        old_weight * (old_velocity[0] * test[0] + old_velocity[1] * test[1]);
                    if !value.is_finite() {
                        return Err(invalid("old momentum load is not representable"));
                    }
                    let row = self.dof(cell, i);
                    load[row] = load[row] + Dd::from_f64(value);
                }
            }
            weights.push(current_weights);
        }
        let native = self.assemble_stokes_admitted_cell_resistance(
            params,
            &|cell, point| normalize(dynamic_viscosity(cell, point), reference).unwrap_or(f64::NAN),
            &|point| {
                force_density(point).map(|value| normalize(value, reference).unwrap_or(f64::NAN))
            },
            &admitted,
            None,
            inv_dt,
            None,
            None,
            None,
            Some(weights),
        )?;
        let mut rhs = native.rhs().to_vec();
        for (edge, geometry) in self.mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let constrained = geometry.tris.1 == usize::MAX && !admitted.is_traction(edge);
            for row in [2 * edge, 2 * edge + 1] {
                let value = load[row].to_f64();
                if !value.is_finite() {
                    return Err(invalid("conservative momentum load is not representable"));
                }
                if !constrained {
                    rhs[row] = (Dd::from_f64(rhs[row]) + load[row]).to_f64();
                    if !rhs[row].is_finite() {
                        return Err(invalid("total momentum load is not representable"));
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(ConservativeMomentumAssembly {
            system: self,
            native,
            rhs,
            params,
            density,
            mass_transport,
            mass_flux,
            reference_density: reference,
        })
    }
}

impl<'m> ConservativeMomentumAssembly<'_, 'm> {
    pub fn matrix(&self) -> &Csr {
        self.native.matrix()
    }
    /// Complete normalized RHS, including old momentum and explicit transport.
    pub fn rhs(&self) -> &[f64] {
        &self.rhs
    }
    /// New density obtained by the original transport operation, not a caller
    /// supplied prediction. It is published jointly with velocity by solve.
    pub fn density(&self) -> &[DensityTransportCell<0>] {
        &self.density
    }
    pub fn mass_transport(&self) -> &DensityTransportReport<0> {
        &self.mass_transport
    }
    pub fn mass_flux(&self) -> &PreparedDensityMassFlux<'m> {
        &self.mass_flux
    }

    /// Solve with the original refined dense-LU/sparse-block-FGMRES selection.
    /// Both original and row-scaled residuals must converge. All physical
    /// history is immutable; refusal or final cancellation publishes no step.
    /// Native factorization/Krylov work itself is synchronous and unpolled.
    pub fn solve(
        self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConservativeMomentumStep<'m>> {
        poll(&mut cancelled)?;
        let mut field =
            self.system
                .solve_linear_result(self.native.matrix(), &self.rhs, self.params)?;
        let mut residual = vec![0.0; self.rhs.len()];
        let (relative_residual, scaled_relative_residual) =
            crate::ns::sparse::residuals(self.native.matrix(), &self.rhs, &field.x, &mut residual);
        field.rel_residual = relative_residual;
        field.scaled_rel_residual = scaled_relative_residual;
        let mut finite = true;
        for chunk in field.x.chunks(1024) {
            poll(&mut cancelled)?;
            finite &= chunk.iter().all(|value| value.is_finite());
        }
        if !finite || !field.converged(self.params.tol) {
            return Err(StokesBoundaryError::NotConverged {
                iterations: field.iters,
                relative_residual,
                scaled_relative_residual,
            }
            .into());
        }
        let mut has_traction = Vec::with_capacity(self.system.pressure_component_count());
        for &cell in &self.system.pressure_roots {
            poll(&mut cancelled)?;
            has_traction.push(
                self.native.pressure_reference_at(cell)
                    == Some(StokesPressureReference::TractionDatum),
            );
        }
        let mut pressure_pa = self.system.pressure_in_gauges(&field.x, &has_traction);
        for chunk in pressure_pa.chunks_mut(1024) {
            poll(&mut cancelled)?;
            for value in chunk {
                *value *= self.reference_density;
                if !value.is_finite() {
                    return Err(invalid("physical pressure is not representable"));
                }
            }
        }
        let kinetic_energy_j_per_m =
            self.reference_density * self.native.kinetic_energy(&field.x)?;
        if !kinetic_energy_j_per_m.is_finite() {
            return Err(invalid("physical kinetic energy is not representable"));
        }
        let mut pressure_reference = Vec::with_capacity(self.density.len());
        for cell in 0..self.density.len() {
            poll(&mut cancelled)?;
            pressure_reference.push(
                self.native
                    .pressure_reference_at(cell)
                    .expect("original cell"),
            );
        }
        poll(&mut cancelled)?;
        Ok(ConservativeMomentumStep {
            field,
            density: self.density,
            pressure_pa,
            pressure_reference,
            kinetic_energy_j_per_m,
            mass_transport: self.mass_transport,
            mass_flux: self.mass_flux,
        })
    }
}

impl<'m> ConservativeMomentumStep<'m> {
    /// Original BDM1 velocity moments [m/s], ready to use as the next history.
    pub fn velocity_m_s(&self) -> &[f64] {
        &self.field.x[..self.field.n_u]
    }

    /// Native velocity moments [m/s] followed by pressure/rho_reference, with
    /// original solver diagnostics. Use pressure_pa for physical pressure.
    pub fn field(&self) -> &FluxSolution {
        &self.field
    }
    pub fn density(&self) -> &[DensityTransportCell<0>] {
        &self.density
    }
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    pub fn pressure_reference_at(&self, cell: usize) -> Option<StokesPressureReference> {
        self.pressure_reference.get(cell).copied()
    }
    pub fn kinetic_energy_j_per_m(&self) -> f64 {
        self.kinetic_energy_j_per_m
    }
    pub fn mass_transport(&self) -> &DensityTransportReport<0> {
        &self.mass_transport
    }
    pub fn mass_flux(&self) -> &PreparedDensityMassFlux<'m> {
        &self.mass_flux
    }
}

fn normalize(value: f64, reference: f64) -> Result<f64> {
    let scaled = value / reference;
    if !value.is_finite() || !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
        Err(invalid(
            "physical momentum normalization is not representable",
        ))
    } else {
        Ok(scaled)
    }
}
fn invalid(message: &'static str) -> ConservativeMomentumError {
    StokesBoundaryError::InvalidInput(message).into()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(DensityTransportError::Cancelled.into())
    } else {
        Ok(())
    }
}
