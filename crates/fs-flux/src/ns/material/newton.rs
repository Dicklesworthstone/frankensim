//! Source-owned constitutive tangents for the existing damped Newton solver.

use super::{
    MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxRheologyStep, MaterialFluxSystem,
    MaterialStokesBoundary, MaterialStokesRheologySolution, mixed::with_physical_boundaries,
    rheology::RheologyField,
};
use crate::ns::boundary::StokesRheologySolution;
use crate::ns::linearization::FluxViscosityTangent;
use crate::ns::newton::{FluxNewtonControls, FluxNewtonError};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::RefCell;

// This file is itself included by path, so a bare `mod mixed;` would resolve
// beside it and compile `material/mixed.rs` a second time.
#[path = "newton/mixed.rs"]
mod mixed;
pub use mixed::adaptive;
pub use mixed::{MaterialNavierStokesNewtonStep, MaterialNavierStokesSdirk2Step};
#[path = "newton/sdirk2.rs"]
mod sdirk2;
pub use sdirk2::MaterialFluxSdirk2Step;

/// Original material response and its analytic derivative at the SAME query.
/// Temperature, composition, position and reference density are held fixed.
pub struct MaterialFluxViscosityTangent<G> {
    /// Original owner supplying the dynamic viscosity and transport snapshot.
    pub material: G,
    /// `d mu / d gamma` [Pa s²], where `gamma = sqrt(2 D:D)` [1/s].
    /// Must be finite and zero at zero strain for this smooth radial law.
    pub derivative_pa_s2: f64,
}

/// Material-bound Newton endpoint, with the accepted original cell sources.
pub struct MaterialFluxNewtonStep<G> {
    /// Physical velocity, pressure, kinetic energy, molecular heat and sources.
    /// Iterations count Newton corrections; the residual belongs to the actual
    /// nonlinear endpoint, not a frozen-viscosity linear solve.
    pub step: MaterialFluxRheologyStep<G>,
    /// Total line-search candidates rejected by the original Newton owner.
    pub backtracks: usize,
    derivatives_pa_s2: Vec<f64>,
}

impl<G> MaterialFluxNewtonStep<G> {
    /// Exact analytic dynamic-viscosity derivatives at the accepted cell rates,
    /// in the same order as `step.samples()`. No source is evaluated again.
    #[must_use]
    pub fn viscosity_derivatives_pa_s2(&self) -> &[f64] {
        &self.derivatives_pa_s2
    }
}

/// Material-bound mixed Stokes Newton endpoint in physical units.
pub struct MaterialStokesRheologyNewtonSolution<G> {
    /// Accepted flow, original cell sources, pressure datum, bulk heat and
    /// physical boundary reactions, using the existing mixed-flow owner.
    pub flow: MaterialStokesRheologySolution<G>,
    /// Total line-search candidates rejected by the original Newton owner.
    pub backtracks: usize,
    derivatives_pa_s2: Vec<f64>,
}

impl<G> MaterialStokesRheologyNewtonSolution<G> {
    /// Exact dynamic-viscosity derivatives [Pa s²] at the accepted cell rates,
    /// in the same order as `flow.samples()`. No source is evaluated again.
    #[must_use]
    pub fn viscosity_derivatives_pa_s2(&self) -> &[f64] {
        &self.derivatives_pa_s2
    }
}

// Both velocity-only and mixed-boundary Newton entries retain the SAME
// original-source derivative and reject unrepresentable density scaling.
fn retain_tangent<G>(
    response: MaterialFluxViscosityTangent<G>,
    point: MaterialFluxRheologyPoint,
    density: f64,
    derivatives: &RefCell<Vec<f64>>,
) -> Result<G, MaterialError> {
    let derivative = response.derivative_pa_s2;
    let normalized = derivative / density;
    if !derivative.is_finite()
        || !normalized.is_finite()
        || (derivative != 0.0 && normalized == 0.0)
        || (point.shear_rate_per_s == 0.0 && derivative != 0.0)
    {
        return Err(MaterialError::State {
            what: "material Newton viscosity derivative must be finite, representable after density scaling, and zero at zero strain".into(),
        });
    }
    derivatives.borrow_mut()[point.cell] = derivative;
    Ok(response.material)
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Solve steady mixed velocity, total-traction and Navier-slip Stokes flow
    /// using original material viscosities and their analytic rate derivatives.
    ///
    /// Force density is N/m³, total traction Pa, wall velocity m/s and Navier
    /// resistance Pa s/m. Viscosity, its derivative, force, traction and slip
    /// resistance all use the SAME fixed reference density for normalization.
    /// Each cell/rate retains one original material response and its derivative;
    /// changing that density or refusing any source query terminates the solve,
    /// including during line search and the final reaction assembly.
    ///
    /// The existing mixed Newton owner differentiates the original symmetric
    /// stress and viscosity-weighted interfaces. Traction and wall resistance
    /// remain fixed physical data. Pressure retains the original traction datum
    /// or componentwise closed-domain gauge. Accepted heat is `mu * gamma²`;
    /// wall-slip loss stays in the physical boundary loads. Reporting reuses
    /// the accepted sources without reevaluating their transport properties.
    ///
    /// Temperature and composition are fixed constitutive inputs. This is
    /// steady creeping flow with a smooth finite radial viscosity law; it adds
    /// no convection, time history, thermal evolution or viscoelastic memory.
    ///
    /// # Errors
    /// Preserves original material errors and refuses invalid tangents,
    /// density changes, inadmissible boundaries and unconverged nonlinear flow.
    pub fn solve_stokes_mixed_rheology_newton<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: FluxNewtonControls,
        mut resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialStokesRheologyNewtonSolution<G>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let derivatives = RefCell::new(vec![0.0; self.system.mesh.tris.len()]);
        let field = RheologyField::new(
            &self.system.mesh.centroids,
            &self.transport,
            |point, reference| {
                retain_tangent(resolve(point, reference)?, point, density, &derivatives)
            },
        );
        let result = with_physical_boundaries(boundaries, density, |primitive| {
            let tangent = |cell, _, rate| {
                let nu = field.viscosity(cell, rate);
                FluxViscosityTangent {
                    nu,
                    derivative: if field.is_valid() {
                        derivatives.borrow()[cell] / density
                    } else {
                        0.0
                    },
                }
            };
            self.system.solve_stokes_mixed_rheology_newton_checked(
                self.params,
                &tangent,
                &|point| {
                    let physical = force_density(point);
                    [physical[0] / density, physical[1] / density]
                },
                primitive,
                initial_velocity,
                controls,
                &|| field.is_valid(),
            )
        })
        .map_err(MaterialFluxError::Newton);
        let (solved, samples) = field.finish(result)?;
        let flow = self.material_stokes_rheology_solution(
            StokesRheologySolution {
                solution: solved.solution,
                iterations: solved.iterations,
                relative_update: solved.relative_update,
                nonlinear_residual: solved.nonlinear_residual,
            },
            samples,
        )?;
        Ok(MaterialStokesRheologyNewtonSolution {
            flow,
            backtracks: solved.backtracks,
            derivatives_pa_s2: derivatives.into_inner(),
        })
    }

    /// Solve steady shear-dependent flow using the original material
    /// viscosities and analytic tangents. `initial_velocity` contains an
    /// immutable guess in m/s; no time history or mass term is introduced.
    /// Set `convection` for steady Navier–Stokes, or leave it false for creeping
    /// flow. Forces are N/m³, and viscosity derivatives are Pa s².
    ///
    /// Uses the same source retention, density scaling, bounded line search
    /// and independently checked nonlinear endpoint as
    /// [`Self::bdf1_rheology_newton_step`]. The returned heat is the physical
    /// `mu * gamma²` at the accepted material state. Temperature/composition
    /// remain fixed constitutive inputs, and source refusals remain terminal.
    ///
    /// # Errors
    /// Preserves original material errors and refuses invalid data,
    /// incompatible boundary flux or unconverged nonlinear equations.
    #[allow(clippy::too_many_arguments)]
    pub fn steady_rheology_newton<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        initial_velocity: &[f64],
        convection: bool,
        controls: FluxNewtonControls,
        resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialFluxNewtonStep<G>, MaterialFluxError> {
        self.rheology_newton_solve(
            force_density,
            boundary_velocity,
            initial_velocity,
            None,
            None,
            None,
            convection,
            controls,
            resolve,
        )
    }

    /// Fully implicit BDF1 using original material viscosities and consistent
    /// analytic tangents in the existing sparse damped Newton solver.
    ///
    /// Forces are N/m³, velocities m/s, time seconds and `dmu/dgamma` Pa s².
    /// Both viscosity and its derivative are divided by the SAME immutable
    /// reference density. The original Jacobian differentiates volume stress,
    /// one-sided interface viscosity, SIP weights, wall loads and convection.
    /// No finite differencing, second constitutive model or viscosity floor is
    /// introduced. Trial residuals reuse the original primal assembler.
    ///
    /// `resolve` supplies one material owner and derivative per cell/rate;
    /// repeated volume/face quadrature queries reuse that response. Its endpoint
    /// temperature/composition must stay fixed throughout the solve. All source
    /// responses must retain the reference density. The final physical pressure,
    /// energy and molecular heat use those exact accepted owners.
    ///
    /// A material refusal is terminal, including in line search: this adapter
    /// never interprets a failed source query as permission to extrapolate or
    /// replace its domain. Pure numerical trial refusals retain the native
    /// Newton backtracking behavior. No thermal evolution, viscoelastic memory,
    /// mixed pressure boundary or variable-density transfer is inferred.
    ///
    /// # Errors
    /// Preserves original source errors; refuses invalid/unrepresentable
    /// derivatives, density changes, incompatible boundaries and nonlinear
    /// nonconvergence. Neither previous velocity nor accepted sources mutate
    /// on refusal. No partial or fallback-coefficient endpoint is published.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_rheology_newton_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
        resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialFluxNewtonStep<G>, MaterialFluxError> {
        self.rheology_newton_solve(
            force_density,
            boundary_velocity,
            previous_velocity,
            Some(timestep_s),
            None,
            None,
            true,
            controls,
            resolve,
        )
    }

    /// Equal-step BDF2 with original material viscosities and analytic
    /// tangents, evaluated at the fully implicit endpoint. `older_velocity`
    /// is one interval before `previous_velocity`; `None` is exact BDF1
    /// startup. Both histories remain unchanged on success or refusal.
    ///
    /// # Errors
    /// Preserves the source, history and Newton refusals of
    /// [`Self::bdf2_rheology_newton_variable_step`].
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_rheology_newton_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        older_velocity: Option<&[f64]>,
        timestep_s: f64,
        controls: FluxNewtonControls,
        resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialFluxNewtonStep<G>, MaterialFluxError> {
        self.bdf2_rheology_newton_variable_step(
            force_density,
            boundary_velocity,
            previous_velocity,
            older_velocity.map(|older| (older, timestep_s)),
            timestep_s,
            controls,
            resolve,
        )
    }

    /// Unequal-step, fully implicit material BDF2 using the consistent sparse
    /// Newton solver. `older_history` pairs the older velocity with the elapsed
    /// time from it to `previous_velocity`. The current/previous step ratio
    /// must lie in `(0, 2]`; absent older history selects exact BDF1 startup.
    ///
    /// Forces are N/m³, velocities m/s, intervals seconds and tangents Pa s².
    /// The actual previous velocity initializes Newton; only the inertial load
    /// uses the BDF2 history combination. Viscosity, its tangent, pressure and
    /// `mu * gamma²` heat use the accepted source owners at the new endpoint.
    /// This heat is an instantaneous rate, not a BDF2 integrated energy loss.
    /// Temperature/composition are fixed endpoint inputs and density remains
    /// the same reference value. A refused source is terminal even in a trial.
    ///
    /// # Errors
    /// Preserves original material errors and refuses invalid histories,
    /// intervals/ratios, density changes, tangents, boundaries or nonlinear
    /// nonconvergence. No partial endpoint or mutated history is published.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_rheology_newton_variable_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: FluxNewtonControls,
        resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialFluxNewtonStep<G>, MaterialFluxError> {
        self.rheology_newton_solve(
            force_density,
            boundary_velocity,
            previous_velocity,
            Some(timestep_s),
            older_history,
            None,
            true,
            controls,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn rheology_newton_solve<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn([f64; 2]) -> [f64; 2],
        initial_velocity: &[f64],
        timestep_s: Option<f64>,
        older_history: Option<(&[f64], f64)>,
        stage_mass: Option<(f64, &[f64])>,
        convection: bool,
        controls: FluxNewtonControls,
        mut resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialFluxNewtonStep<G>, MaterialFluxError> {
        // Admit and prepare the two-history mass load before with_inputs
        // queries endpoint boundary data. Keep u_n as the constitutive/Newton
        // initial state and reuse this single preparation for every iteration.
        let bdf2 = if let (Some(timestep), Some(_)) = (timestep_s, older_history) {
            Some(
                self.system
                    .mixed_bdf2_history(initial_velocity, older_history, timestep)
                    .map_err(|error| MaterialFluxError::Newton(FluxNewtonError::History(error)))?,
            )
        } else {
            None
        };
        let density = self.transport.density_kg_m3();
        let derivatives = RefCell::new(vec![0.0; self.system.mesh.tris.len()]);
        let field = RheologyField::new(
            &self.system.mesh.centroids,
            &self.transport,
            |point, reference| {
                retain_tangent(resolve(point, reference)?, point, density, &derivatives)
            },
        );
        let result = self.with_inputs(
            force_density,
            boundary_velocity,
            |acceleration, velocity| {
                let tangent = |cell, _, rate| {
                    let nu = field.viscosity(cell, rate);
                    FluxViscosityTangent {
                        nu,
                        // A source refusal only finishes the already-started
                        // assembly. Its fallback cannot be solved or accepted.
                        derivative: if field.is_valid() {
                            derivatives.borrow()[cell] / density
                        } else {
                            0.0
                        },
                    }
                };
                if let Some((stage_timestep, history)) = stage_mass {
                    self.system.rheology_newton_solve_checked(
                        self.params,
                        &tangent,
                        acceleration,
                        velocity,
                        initial_velocity,
                        convection,
                        stage_timestep.recip(),
                        Some(history),
                        controls,
                        &|| field.is_valid(),
                    )
                } else if let Some((history, effective_timestep)) = &bdf2 {
                    self.system.rheology_newton_solve_checked(
                        self.params,
                        &tangent,
                        acceleration,
                        velocity,
                        initial_velocity,
                        true,
                        effective_timestep.recip(),
                        Some(history.as_deref().unwrap_or(initial_velocity)),
                        controls,
                        &|| field.is_valid(),
                    )
                } else if let Some(timestep) = timestep_s {
                    self.system.bdf1_rheology_newton_step_checked(
                        self.params,
                        &tangent,
                        acceleration,
                        velocity,
                        initial_velocity,
                        timestep,
                        controls,
                        &|| field.is_valid(),
                    )
                } else {
                    self.system.rheology_newton_solve_checked(
                        self.params,
                        &tangent,
                        acceleration,
                        velocity,
                        initial_velocity,
                        convection,
                        0.0,
                        None,
                        controls,
                        &|| field.is_valid(),
                    )
                }
            },
        );
        let (step, samples) =
            field.finish(result.and_then(|step| step.map_err(MaterialFluxError::Newton)))?;
        let flow = self.physical_solution_with_viscosity_tolerance(
            step.solution,
            &|cell| samples[cell].transport().dynamic_viscosity_pa_s(),
            controls.residual_tolerance,
        )?;
        Ok(MaterialFluxNewtonStep {
            step: MaterialFluxRheologyStep {
                flow,
                samples,
                iterations: step.iterations,
                relative_update: step.relative_update,
                nonlinear_residual: step.nonlinear_residual,
            },
            backtracks: step.backtracks,
            derivatives_pa_s2: derivatives.into_inner(),
        })
    }
}
