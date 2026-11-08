//! Stage-specific original material owners on the shared SDIRK2 time method.

use super::{MaterialFluxNewtonStep, MaterialFluxViscosityTangent};
use crate::ns::material::{MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxSystem};
use crate::ns::newton::{FluxNewtonControls, FluxSdirk2Error, FluxSdirk2Stage};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// Accepted SDIRK2 stages with each stage's original constitutive responses.
pub struct MaterialFluxSdirk2Step<G> {
    /// Internal-stage velocity, pressure, instantaneous heat and material owners.
    pub first_stage: MaterialFluxNewtonStep<G>,
    /// Final velocity, pressure, instantaneous heat and material owners.
    pub endpoint: MaterialFluxNewtonStep<G>,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Advance second-order, stiffly accurate SDIRK2 with a separate material
    /// law, physical body force and boundary velocity at each stage. Forces are
    /// N/m³, velocities m/s, time seconds, and viscosity derivatives Pa s².
    ///
    /// Stage time is `t_n + stage.time_fraction() * timestep_s`. `resolve`
    /// retains one original source per cell/rate in each stage, including its
    /// supplied temperature/composition. Those external fields remain fixed
    /// within that stage's Newton solve; this method adds no thermal evolution.
    /// Both stages retain the SAME reference density for force, viscosity,
    /// tangent, pressure and energy conversions. `convection` selects whether
    /// advective momentum is included, while both modes retain physical inertia.
    ///
    /// The source cache and exact constitutive Jacobian are the existing
    /// material Newton owners. Stage two starts from the first-stage velocity,
    /// never the extrapolated inertial load. Each stage's heat is the original
    /// `mu * gamma²` instantaneous rate; neither it nor a kinetic-energy
    /// difference alone is a discrete heat balance. Work limits apply per stage.
    ///
    /// # Errors
    /// Preserves the failing stage and original material/numerical error.
    /// Invalid old history/timestep refuses before callbacks. Material source
    /// refusals are terminal even during a line search. No partial stage pair
    /// is returned and the caller's old velocity is never modified.
    #[allow(clippy::too_many_arguments)]
    pub fn sdirk2_rheology_newton_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn(FluxSdirk2Stage, [f64; 2]) -> [f64; 2],
        boundary_velocity: &dyn Fn(FluxSdirk2Stage, [f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
        convection: bool,
        controls: FluxNewtonControls,
        mut resolve: impl FnMut(
            FluxSdirk2Stage,
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialFluxSdirk2Step<G>, FluxSdirk2Error<MaterialFluxError>> {
        let stage_dt = self
            .system
            .sdirk2_stage_timestep(previous_velocity, timestep_s)
            .map_err(|error| FluxSdirk2Error {
                stage: None,
                error: MaterialFluxError::Newton(error),
            })?;
        let mut solve = |stage, initial: &[f64], history: &[f64]| {
            self.rheology_newton_solve(
                &|point| force_density(stage, point),
                &|point| boundary_velocity(stage, point),
                initial,
                None,
                None,
                Some((stage_dt, history)),
                convection,
                controls,
                |point, reference| resolve(stage, point, reference),
            )
            .map_err(|error| FluxSdirk2Error {
                stage: Some(stage),
                error,
            })
        };
        let first_stage = solve(FluxSdirk2Stage::First, previous_velocity, previous_velocity)?;
        let first_velocity = &first_stage.step.discrete().x[..self.system.n_u];
        let history = self
            .system
            .sdirk2_mass_history(previous_velocity, first_velocity)
            .map_err(|error| FluxSdirk2Error {
                stage: Some(FluxSdirk2Stage::Endpoint),
                error: MaterialFluxError::Newton(error),
            })?;
        let endpoint = solve(FluxSdirk2Stage::Endpoint, first_velocity, &history)?;
        Ok(MaterialFluxSdirk2Step {
            first_stage,
            endpoint,
        })
    }
}
