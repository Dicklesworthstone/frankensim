//! Shared-material entry to the existing BDM1 incompressible solver.
//!
//! The homogeneous entry advances `du/dt + u.grad(u) - nu*lap(u) + grad(pi)
//! = a`. The rheological entry uses `-div(2 nu D(u))`, including one-sided
//! viscous tractions at material interfaces. Both use `nu = mu/rho`,
//! `pi = p/rho`, and `a = force_density/rho` with the SAME constant density.
//! This binding owns those conversions together. It retains the original
//! gas/liquid evaluation (and hence its source receipts) and a single immutable
//! transport snapshot. A material swap never requires a solver-specific fit.
//!
//! Geometry is in metres, velocity in m/s, force density in N/m³, and time
//! in seconds. The base binding is a homogeneous Newtonian incompressible
//! snapshot. Its rheological BDF1 and mixed Stokes entries resolve cellwise
//! material responses at the evolving strain, retaining their exact owners.
//! Neither entry infers gas compressibility, phase validity or thermal evolution.

use std::cell::Cell;

use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_sparse::Csr;

use crate::ns::{FluxParams, FluxRheologyError, FluxSolution, FluxSystem};

#[path = "material/mixed.rs"]
pub mod mixed;
#[path = "material/rheology.rs"]
mod rheology;
pub use mixed::{MaterialStokesBoundary, MaterialStokesRheologySolution};

/// Numerical choices, separate from the material's physical viscosity.
#[derive(Clone, Copy, Debug)]
pub struct MaterialFluxNumerics {
    /// Dimensionless symmetric interior-penalty coefficient.
    pub penalty: f64,
    /// Required relative residual of the final linear solve.
    pub linear_tolerance: f64,
    /// Krylov iteration budget inherited by the primitive solver.
    pub max_linear_iterations: usize,
}

impl Default for MaterialFluxNumerics {
    fn default() -> Self {
        let params = FluxParams::default();
        Self {
            penalty: params.sigma,
            linear_tolerance: params.tol,
            max_linear_iterations: params.max_iters,
        }
    }
}

/// Refusal before a physical result is published.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialFluxError {
    /// The original material owner refused its transport state.
    Material(MaterialError),
    /// Invalid physical/numerical input or an unrepresentable conversion.
    Invalid(&'static str),
    /// The primitive solver did not satisfy the requested linear residual.
    LinearConvergence {
        /// Measured relative residual.
        residual: f64,
        /// Requested relative tolerance.
        tolerance: f64,
    },
    /// The steady nonlinear iteration exhausted its update criterion.
    NonlinearConvergence {
        /// Completed Picard iterations.
        iterations: usize,
        /// Last relative velocity update.
        update: f64,
    },
    /// The fully implicit constitutive/momentum step refused its candidate.
    Rheology(FluxRheologyError),
    /// Mixed-boundary admission or a checked constitutive flow solve refused.
    Boundary(crate::ns::boundary::StokesBoundaryError),
}

impl core::fmt::Display for MaterialFluxError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "material flow: {error}"),
            Self::Invalid(message) => write!(f, "material flow: {message}"),
            Self::LinearConvergence {
                residual,
                tolerance,
            } => write!(
                f,
                "material flow: linear residual {residual} exceeds {tolerance}"
            ),
            Self::NonlinearConvergence { iterations, update } => write!(
                f,
                "material flow: Picard did not converge in {iterations} iterations (update {update})"
            ),
            Self::Rheology(error) => write!(f, "material flow: {error}"),
            Self::Boundary(error) => write!(f, "material flow: {error}"),
        }
    }
}
impl std::error::Error for MaterialFluxError {}
impl From<MaterialError> for MaterialFluxError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}

/// One mesh/operator owner bound to one homogeneous material evaluation.
/// The reference preserves the source owner; a raw `FluidTransportState`
/// remains an explicitly unsourced numerical input.
pub struct MaterialFluxSystem<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    system: &'s FluxSystem<'m>,
    material: &'a F,
    transport: FluidTransportState,
    params: FluxParams,
}

/// A solved velocity and pressure with physical outputs from the SAME material
/// snapshot. The primitive pressure coordinates remain `p/rho`; use
/// `pressure_pa()` for physical pressure differences.
pub struct MaterialFluxSolution {
    discrete: FluxSolution,
    transport: FluidTransportState,
    pressure_pa: Vec<f64>,
    kinetic_energy_j_per_m: f64,
}

/// Physical coordinates for a cellwise generalized-Newtonian material query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFluxRheologyPoint {
    /// Triangle index in the bound mesh.
    pub cell: usize,
    /// Cell centroid [m]. The endpoint material field is piecewise constant.
    pub position_m: [f64; 2],
    /// Objective strain invariant `sqrt(2 D:D)` [1/s] at this iterate.
    pub shear_rate_per_s: f64,
}

/// Exact retained material response used at the accepted endpoint strain.
pub struct MaterialFluxRheologySample<F> {
    point: MaterialFluxRheologyPoint,
    material: F,
    transport: FluidTransportState,
}

impl<F> MaterialFluxRheologySample<F> {
    /// Actual material query, including the endpoint shear rate.
    #[must_use]
    pub const fn point(&self) -> MaterialFluxRheologyPoint {
        self.point
    }
    /// Original material owner and its source/property receipts.
    #[must_use]
    pub const fn material(&self) -> &F {
        &self.material
    }
    /// Validated transport response read once from this owner.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
}

/// Accepted fully implicit flow and the exact endpoint material field.
/// Pressure and energy use the immutable reference density; viscosity comes
/// from `samples()`, with one retained owner per triangle in mesh order.
pub struct MaterialFluxRheologyStep<F> {
    flow: MaterialFluxSolution,
    samples: Vec<MaterialFluxRheologySample<F>>,
    /// Coupled constitutive/convection iterations.
    pub iterations: usize,
    /// Last relative velocity update.
    pub relative_update: f64,
    /// Residual reassembled at the accepted velocity and material field.
    pub nonlinear_residual: f64,
}

impl<F> MaterialFluxRheologyStep<F> {
    /// Primitive discrete result, including its linear residual.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        self.flow.discrete()
    }
    /// BDM1 velocity moments [m/s] for the next physical timestep.
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        self.flow.velocity_dofs()
    }
    /// Mean-zero cell pressure differences [Pa].
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        self.flow.pressure_pa()
    }
    /// Kinetic energy per out-of-plane depth [J/m].
    #[must_use]
    pub const fn kinetic_energy_j_per_m(&self) -> f64 {
        self.flow.kinetic_energy_j_per_m()
    }
    /// Endpoint material responses, indexed by triangle.
    #[must_use]
    pub fn samples(&self) -> &[MaterialFluxRheologySample<F>] {
        &self.samples
    }
}

impl MaterialFluxSolution {
    /// Original discrete result, including residual and iteration count.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        &self.discrete
    }
    /// BDM1 normal-velocity moments [m/s], suitable for the next BDF1 step.
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        &self.discrete.x[..self.discrete.n_u]
    }
    /// Cell pressure differences [Pa], shifted to zero area-weighted mean.
    /// The incompressible pressure gauge is not thermodynamic absolute pressure.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    /// Kinetic energy per unit out-of-plane depth [J/m] for this 2-D mesh.
    #[must_use]
    pub const fn kinetic_energy_j_per_m(&self) -> f64 {
        self.kinetic_energy_j_per_m
    }
    /// The exact material snapshot used for assembly and physical outputs.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> MaterialFluxSystem<'s, 'm, 'a, F> {
    /// Bind an existing FEEC operator to the shared gas/liquid transport owner.
    ///
    /// # Errors
    /// Refuses missing material data, an empty mesh, or invalid solver settings.
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
    ) -> Result<Self, MaterialFluxError> {
        let transport = material.transport_state()?;
        if system.mesh.tris.is_empty()
            || !positive(numerics.penalty)
            || !positive(numerics.linear_tolerance)
            || numerics.max_linear_iterations == 0
        {
            return Err(MaterialFluxError::Invalid(
                "nonempty mesh and positive finite numerical controls required",
            ));
        }
        let params = FluxParams {
            nu: transport.kinematic_viscosity_m2_s(),
            sigma: numerics.penalty,
            tol: numerics.linear_tolerance,
            max_iters: numerics.max_linear_iterations,
        };
        Ok(Self {
            system,
            material,
            transport,
            params,
        })
    }

    /// Original material owner, retaining its claims and query conditions.
    #[must_use]
    pub const fn material(&self) -> &F {
        self.material
    }
    /// One snapshot supplies both momentum transport and physical output units.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }

    /// Assemble steady Stokes with volumetric force [N/m³] and velocity [m/s].
    /// The returned saddle operator uses the primitive `p/rho` coordinate and
    /// can be consumed by its existing adjoint/linear-solver interfaces.
    ///
    /// # Errors
    /// Refuses nonfinite boundary/load samples or overflow in force/rho.
    pub fn assemble_stokes(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
    ) -> Result<(Csr, Vec<f64>), MaterialFluxError> {
        self.with_inputs(
            force_density,
            boundary_velocity,
            |acceleration, velocity| {
                self.system
                    .assemble(self.params, acceleration, velocity, None)
            },
        )
    }

    /// Solve steady Stokes through the existing BDM1 operator.
    ///
    /// # Errors
    /// Refuses invalid inputs, unconverged solves, and nonfinite physical outputs.
    /// # Panics
    /// Inherits the primitive solver's singular-mesh/factorization contract.
    pub fn solve_stokes(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
    ) -> Result<MaterialFluxSolution, MaterialFluxError> {
        let (operator, rhs) = self.assemble_stokes(force_density, boundary_velocity)?;
        self.physical_solution(self.system.solve_linear(&operator, &rhs, self.params))
    }

    /// Steady Navier–Stokes using the existing Picard iteration. Loads are
    /// N/m³ and boundary velocities m/s. Returns the physical solution, count,
    /// and last relative velocity update; no unconverged result is published.
    ///
    /// # Errors
    /// Refuses invalid inputs, invalid budgets, failed convergence or overflow.
    /// # Panics
    /// Inherits the primitive solver's singular-mesh/factorization contract.
    pub fn solve_picard(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        max_iterations: usize,
        relative_tolerance: f64,
    ) -> Result<(MaterialFluxSolution, usize, f64), MaterialFluxError> {
        if max_iterations == 0 || !positive(relative_tolerance) {
            return Err(MaterialFluxError::Invalid(
                "Picard needs a positive iteration budget and finite positive tolerance",
            ));
        }
        let (solution, iterations, update) = self.with_inputs(
            force_density,
            boundary_velocity,
            |acceleration, velocity| {
                self.system.picard(
                    self.params,
                    acceleration,
                    velocity,
                    max_iterations,
                    relative_tolerance,
                )
            },
        )?;
        if !update.is_finite() || update >= relative_tolerance {
            return Err(MaterialFluxError::NonlinearConvergence { iterations, update });
        }
        Ok((self.physical_solution(solution)?, iterations, update))
    }

    /// One IMEX BDF1 step in physical seconds. `previous_velocity` contains
    /// exactly the mesh's BDM1 velocity moments [m/s]. This same homogeneous
    /// density applies to both time levels; changing material density during a
    /// trajectory requires a separate conservative state-transfer model.
    ///
    /// # Errors
    /// Refuses invalid velocity dimensions/values, timestep, loads or output.
    /// # Panics
    /// Inherits the primitive solver's singular-mesh/factorization contract.
    pub fn bdf1_step(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<MaterialFluxSolution, MaterialFluxError> {
        if previous_velocity.len() != self.system.n_u
            || !previous_velocity.iter().all(|value| value.is_finite())
            || !positive(timestep_s)
            || !positive(timestep_s.recip())
        {
            return Err(MaterialFluxError::Invalid(
                "BDF1 needs finite velocity moments and a representable positive timestep in seconds",
            ));
        }
        let solution = self.with_inputs(
            force_density,
            boundary_velocity,
            |acceleration, velocity| {
                self.system.bdf1_step(
                    self.params,
                    acceleration,
                    velocity,
                    previous_velocity,
                    timestep_s,
                )
            },
        )?;
        self.physical_solution(solution)
    }

    /// Fully implicit BDF1 with source-backed generalized-Newtonian viscosity.
    /// Loads are N/m³, velocity moments m/s, and the timestep seconds.
    ///
    /// `resolve` owns the endpoint temperature, composition and other material
    /// context. Those fields must remain fixed during this solve. It receives
    /// each triangle's centroid and current strain rate, and the immutable
    /// reference transport snapshot. The returned viscosity is cellwise
    /// constant; volume and both one-sided face traces use that SAME response.
    /// Responses at unchanged cell/rate coordinates are reused, avoiding source
    /// queries inside quadrature loops. No constitutive fit or floor is added.
    ///
    /// Every response must retain the reference density, since momentum,
    /// pressure, energy and both time levels share one incompressible inertia.
    /// A changed density requires a separate conservative state-transfer model.
    /// Material owners at the accepted endpoint are retained without requery.
    /// This method does not evolve temperature or add viscoelastic memory.
    ///
    /// # Errors
    /// Refuses material/source errors, changed density, invalid inputs and
    /// unconverged implicit steps. A source refusal stops before the next linear
    /// solve; no fallback transport or partial solution is published. The
    /// caller's previous velocity is never mutated.
    /// # Panics
    /// Inherits the primitive solver's singular-mesh/factorization contract.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_rheology_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
        maximum_iterations: usize,
        relative_tolerance: f64,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialFluxRheologyStep<G>, MaterialFluxError> {
        let field =
            rheology::RheologyField::new(&self.system.mesh.centroids, &self.transport, resolve);
        let result = self.with_inputs(
            force_density,
            boundary_velocity,
            |acceleration, velocity| {
                self.system.bdf1_rheology_step_checked(
                    self.params,
                    &|cell, _, rate| field.viscosity(cell, rate),
                    acceleration,
                    velocity,
                    previous_velocity,
                    timestep_s,
                    maximum_iterations,
                    relative_tolerance,
                    &|| field.is_valid(),
                )
            },
        );
        let (step, samples) =
            field.finish(result.and_then(|step| step.map_err(MaterialFluxError::Rheology)))?;
        Ok(MaterialFluxRheologyStep {
            flow: self.physical_solution(step.solution)?,
            samples,
            iterations: step.iterations,
            relative_update: step.relative_update,
            nonlinear_residual: step.nonlinear_residual,
        })
    }

    fn with_inputs<T>(
        &self,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        run: impl FnOnce(&dyn Fn([f64; 2]) -> [f64; 2], &dyn Fn([f64; 2]) -> [f64; 2]) -> T,
    ) -> Result<T, MaterialFluxError> {
        let invalid = Cell::new(false);
        let density = self.transport.density_kg_m3();
        let acceleration = |point| {
            let load = force(point);
            let value = [load[0] / density, load[1] / density];
            if value.iter().all(|component| component.is_finite()) {
                value
            } else {
                invalid.set(true);
                [0.0; 2]
            }
        };
        let velocity = |point| {
            let value = boundary(point);
            if value.iter().all(|component| component.is_finite()) {
                value
            } else {
                invalid.set(true);
                [0.0; 2]
            }
        };
        let result = run(&acceleration, &velocity);
        if invalid.get() {
            Err(MaterialFluxError::Invalid(
                "nonfinite force-density/boundary sample or force/rho conversion",
            ))
        } else {
            Ok(result)
        }
    }

    fn physical_solution(
        &self,
        discrete: FluxSolution,
    ) -> Result<MaterialFluxSolution, MaterialFluxError> {
        if !discrete.converged(self.params.tol) {
            return Err(MaterialFluxError::LinearConvergence {
                residual: discrete.rel_residual.max(discrete.scaled_rel_residual),
                tolerance: self.params.tol,
            });
        }
        if discrete.x.len() != self.system.n || !discrete.x.iter().all(|value| value.is_finite()) {
            return Err(MaterialFluxError::Invalid(
                "solver returned an invalid discrete state",
            ));
        }
        let density = self.transport.density_kg_m3();
        let pressure_pa: Vec<_> = self
            .system
            .pressure(&discrete.x)
            .into_iter()
            .map(|value| density * value)
            .collect();
        let kinetic_energy_j_per_m = density * self.system.kinetic_energy(&discrete.x);
        if !pressure_pa.iter().all(|value| value.is_finite()) || !kinetic_energy_j_per_m.is_finite()
        {
            return Err(MaterialFluxError::Invalid(
                "physical pressure or energy is not representable",
            ));
        }
        Ok(MaterialFluxSolution {
            discrete,
            transport: self.transport,
            pressure_pa,
            kinetic_energy_j_per_m,
        })
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
