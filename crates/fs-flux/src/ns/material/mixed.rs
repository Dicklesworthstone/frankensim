//! Source-owned, shear-dependent steady flow with physical traction boundaries.

use super::{
    MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxRheologySample, MaterialFluxSystem,
    rheology::RheologyField,
};
use crate::ns::FluxSolution;
use crate::ns::boundary::{
    StokesBoundary, StokesPressureReference, StokesRheologyControls, StokesRheologySolution,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// Exactly one physical velocity or total-traction owner per boundary edge.
/// Geometry is in metres. Callbacks must be pure and fixed during the solve.
#[derive(Clone, Copy)]
pub enum MaterialStokesBoundary<'a> {
    /// Prescribed velocity, with the existing strong-normal/SIP-tangent rule.
    Velocity {
        /// Mesh boundary edge index.
        edge: usize,
        /// Physical velocity [m/s].
        value: &'a dyn Fn([f64; 2]) -> [f64; 2],
    },
    /// Physical total stress vector (2 mu D(u) - p I)n.
    TotalTraction {
        /// Mesh boundary edge index.
        edge: usize,
        /// Traction [Pa], evaluated at position [m] with outward unit normal.
        value: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
    },
}

/// Accepted source-backed steady flow and its exact final material responses.
pub struct MaterialStokesRheologySolution<G> {
    flow: StokesRheologySolution,
    pressure_pa: Vec<f64>,
    kinetic_energy_j_per_m: f64,
    samples: Vec<MaterialFluxRheologySample<G>>,
}

impl<G> MaterialStokesRheologySolution<G> {
    /// Primitive velocity/kinematic-pressure coefficients and independently
    /// reassembled original and row-scaled nonlinear endpoint residuals.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        &self.flow.solution.field
    }
    /// Physical velocity coefficients [m/s].
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        let x = &self.flow.solution.field.x;
        &x[..x.len() - self.pressure_pa.len()]
    }
    /// Physical pressure [Pa] in the retained traction datum, or mean-zero
    /// gauge for all-velocity boundaries; not an EOS pressure declaration.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    /// Declared pressure interpretation, with no hidden shift at open ends.
    #[must_use]
    pub const fn pressure_reference(&self) -> StokesPressureReference {
        self.flow.solution.pressure_reference
    }
    /// Physical kinetic energy per metre of out-of-plane depth.
    #[must_use]
    pub const fn kinetic_energy_j_per_m(&self) -> f64 {
        self.kinetic_energy_j_per_m
    }
    /// Exact final source owners, in triangle order, without a reporting requery.
    #[must_use]
    pub fn samples(&self) -> &[MaterialFluxRheologySample<G>] {
        &self.samples
    }
    /// Completed linearized constitutive/flow solves.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.flow.iterations
    }
    /// Final undamped relative velocity correction.
    #[must_use]
    pub const fn relative_update(&self) -> f64 {
        self.flow.relative_update
    }
    /// Larger of the two independently reassembled nonlinear residuals.
    #[must_use]
    pub const fn nonlinear_residual(&self) -> f64 {
        self.flow.nonlinear_residual
    }
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Solve steady generalized-Newtonian mixed Stokes with physical body
    /// force [N/m³], velocity [m/s] and TOTAL traction [Pa].
    ///
    /// Each constitutive response must retain the reference density. Dynamic
    /// viscosity, force and traction are divided by that SAME density; pressure
    /// and kinetic energy are converted back together. A traction pressure
    /// datum is retained rather than shifted to zero mean.
    ///
    /// resolve receives the actual cell centroid and shear invariant [1/s].
    /// Captured temperature/composition remain fixed during this solve. One
    /// cellwise response supplies volume and both one-sided interface traces.
    /// The shared BDF1 material-query owner caches unchanged cell/rate queries
    /// and retains the exact source owners used by the accepted final residual.
    /// It adds no fit, viscosity floor, thermal evolution or viscoelastic state.
    ///
    /// # Errors
    /// Preserves original material refusals and rejects changed density,
    /// invalid boundaries/conversions, nonfinite states and unconverged solves.
    /// Refusal publishes neither fallback properties nor a partial flow, and
    /// never mutates initial_velocity.
    pub fn solve_stokes_mixed_rheology<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialStokesRheologySolution<G>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        // Homogeneous closure storage avoids a boxed callback per edge.
        // Admission still happens before any force, traction or source query.
        let scaled_tractions: Vec<_> = boundaries
            .iter()
            .map(|boundary| {
                move |point: [f64; 2], normal: [f64; 2]| match boundary {
                    MaterialStokesBoundary::TotalTraction { value, .. } => {
                        let physical = value(point, normal);
                        [physical[0] / density, physical[1] / density]
                    }
                    MaterialStokesBoundary::Velocity { .. } => {
                        unreachable!("only traction owners request stress samples")
                    }
                }
            })
            .collect();
        let primitive_boundaries: Vec<_> = boundaries
            .iter()
            .zip(&scaled_tractions)
            .map(|(boundary, scaled)| match boundary {
                MaterialStokesBoundary::Velocity { edge, value } => StokesBoundary::Velocity {
                    edge: *edge,
                    value: *value,
                },
                MaterialStokesBoundary::TotalTraction { edge, .. } => {
                    StokesBoundary::TotalTraction {
                        edge: *edge,
                        value: scaled,
                    }
                }
            })
            .collect();
        let field = RheologyField::new(&self.system.mesh.centroids, &self.transport, resolve);
        let result = self
            .system
            .solve_stokes_mixed_rheology_checked(
                self.params,
                &|cell, _, rate| field.viscosity(cell, rate),
                &|point| {
                    let physical = force_density(point);
                    [physical[0] / density, physical[1] / density]
                },
                &primitive_boundaries,
                initial_velocity,
                controls,
                &|| field.is_valid(),
            )
            .map_err(MaterialFluxError::Boundary);
        let (flow, samples) = field.finish(result)?;
        let pressure_pa: Vec<_> = flow.solution.pressure.iter().map(|p| density * p).collect();
        let kinetic_energy_j_per_m = density * self.system.kinetic_energy(&flow.solution.field.x);
        if pressure_pa.iter().any(|value| !value.is_finite()) || !kinetic_energy_j_per_m.is_finite()
        {
            return Err(MaterialFluxError::Invalid(
                "physical pressure or energy is not representable",
            ));
        }
        Ok(MaterialStokesRheologySolution {
            flow,
            pressure_pa,
            kinetic_energy_j_per_m,
            samples,
        })
    }
}
