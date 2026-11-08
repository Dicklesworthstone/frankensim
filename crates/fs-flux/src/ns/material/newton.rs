//! Source-owned constitutive tangents for the existing damped Newton solver.

use super::{
    MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxRheologyStep, MaterialFluxSystem,
    rheology::RheologyField,
};
use crate::ns::linearization::FluxViscosityTangent;
use crate::ns::newton::FluxNewtonControls;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::RefCell;

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

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
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
        convection: bool,
        controls: FluxNewtonControls,
        mut resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialFluxNewtonStep<G>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let derivatives = RefCell::new(vec![0.0; self.system.mesh.tris.len()]);
        let field = RheologyField::new(
            &self.system.mesh.centroids,
            &self.transport,
            |point, reference| {
                let response = resolve(point, reference)?;
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
                if let Some(timestep) = timestep_s {
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
