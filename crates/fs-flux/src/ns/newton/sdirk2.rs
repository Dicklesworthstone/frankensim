//! Two-stage stiffly accurate SDIRK2 on the original Newton momentum owner.

use super::{FluxNewtonControls, FluxNewtonError, FluxNewtonStep};
use crate::ns::linearization::FluxViscosityTangent;
use crate::ns::{FluxParams, FluxSystem};

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

/// The physical time at which one SDIRK2 stage's external data are supplied.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FluxSdirk2Stage {
    /// The internal stage at `t_n + (1 - 1/sqrt(2)) * dt`.
    First,
    /// The stiffly accurate endpoint at `t_n + dt`.
    Endpoint,
}

impl FluxSdirk2Stage {
    /// Fraction of the requested physical timestep for this stage's data.
    #[must_use]
    pub const fn time_fraction(self) -> f64 {
        match self {
            Self::First => GAMMA,
            Self::Endpoint => 1.0,
        }
    }
}

/// Both accepted stage fields, published only after the whole step succeeds.
pub struct FluxSdirk2Step {
    /// Accepted internal stage, including its actual nonlinear work/residual.
    pub first_stage: FluxNewtonStep,
    /// Accepted final velocity and pressure at the requested new time.
    pub endpoint: FluxNewtonStep,
}

/// The original refusal together with the stage that failed.
#[derive(Clone, Debug, PartialEq)]
pub struct FluxSdirk2Error<E> {
    /// `None` means the input history or timestep failed before either stage.
    pub stage: Option<FluxSdirk2Stage>,
    /// Original numerical or material refusal; no partial flow is returned.
    pub error: E,
}

impl<E: core::fmt::Display> core::fmt::Display for FluxSdirk2Error<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self.stage {
            Some(stage) => write!(f, "SDIRK2 flow {stage:?}: {}", self.error),
            None => write!(f, "SDIRK2 flow input: {}", self.error),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for FluxSdirk2Error<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
    }
}

impl FluxSystem<'_> {
    /// Advance fully implicit second-order SDIRK2 with the exact constitutive
    /// Jacobian at BOTH stages. This is the stiffly accurate, L-stable method
    /// with `gamma = 1 - 1/sqrt(2)` and tableau rows `[gamma, 0]` and
    /// `[1-gamma, gamma]`. No previous interval or two-step startup is needed.
    ///
    /// Each callback receives the stage whose physical time is
    /// `t_n + stage.time_fraction() * timestep_s`. Force, boundary velocity,
    /// temperature/composition captured by the viscosity law and its analytic
    /// tangent must belong to that time and remain fixed during its Newton
    /// solve. Set `convection` for Navier–Stokes; false retains inertia and
    /// nonlinear viscous stress but omits convective momentum.
    ///
    /// The first stage solves from `u_n`. The second mass history is
    /// `u_n + (1-gamma)/gamma * (u_first-u_n)`, while its Newton initial guess
    /// is the accepted first-stage velocity. This mass combination is never
    /// treated as a constitutive state. Both stages use the existing primal
    /// residual, consistent Jacobian, line search and endpoint acceptance.
    /// The supplied work limits apply separately to each stage.
    ///
    /// Second-order accuracy requires resolved spatial/nonlinear errors and
    /// sufficiently smooth stage data; L-stability is the linear decay
    /// property, not a general nonlinear positivity or energy theorem. All
    /// boundaries prescribe velocity. Input history is immutable, and no
    /// partial step is returned when either stage refuses.
    ///
    /// # Errors
    /// Reports the failing stage and original Newton refusal, including
    /// invalid/unrepresentable history, timestep or second-stage mass load.
    #[allow(clippy::too_many_arguments)]
    pub fn sdirk2_rheology_newton_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(FluxSdirk2Stage, usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn(FluxSdirk2Stage, [f64; 2]) -> [f64; 2],
        boundary: &dyn Fn(FluxSdirk2Stage, [f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
        convection: bool,
        controls: FluxNewtonControls,
    ) -> Result<FluxSdirk2Step, FluxSdirk2Error<FluxNewtonError>> {
        let stage_dt = self
            .sdirk2_stage_timestep(previous_velocity, timestep_s)
            .map_err(|error| FluxSdirk2Error { stage: None, error })?;
        let solve = |stage, initial: &[f64], history: &[f64]| {
            self.rheology_newton_solve_checked(
                params,
                &|cell, point, rate| law(stage, cell, point, rate),
                &|point| force(stage, point),
                &|point| boundary(stage, point),
                initial,
                convection,
                stage_dt.recip(),
                Some(history),
                controls,
                &|| true,
            )
            .map_err(|error| FluxSdirk2Error {
                stage: Some(stage),
                error,
            })
        };
        let first_stage = solve(FluxSdirk2Stage::First, previous_velocity, previous_velocity)?;
        let first_velocity = &first_stage.solution.x[..self.n_u];
        let history = self
            .sdirk2_mass_history(previous_velocity, first_velocity)
            .map_err(|error| FluxSdirk2Error {
                stage: Some(FluxSdirk2Stage::Endpoint),
                error,
            })?;
        let endpoint = solve(FluxSdirk2Stage::Endpoint, first_velocity, &history)?;
        Ok(FluxSdirk2Step {
            first_stage,
            endpoint,
        })
    }

    // Shared by native and source-owned stages, before querying stage data.
    pub(crate) fn sdirk2_stage_timestep(
        &self,
        previous: &[f64],
        timestep_s: f64,
    ) -> Result<f64, FluxNewtonError> {
        let stage_dt = GAMMA * timestep_s;
        if self.mesh.tris.is_empty()
            || previous.len() != self.n_u
            || previous.iter().any(|value| !value.is_finite())
            || !timestep_s.is_finite()
            || timestep_s <= 0.0
            || !stage_dt.is_finite()
            || stage_dt <= 0.0
            || !stage_dt.recip().is_finite()
            || stage_dt.recip() <= 0.0
        {
            return Err(FluxNewtonError::Invalid(
                "SDIRK2 needs finite velocity history and a representable positive stage timestep",
            ));
        }
        Ok(stage_dt)
    }

    pub(crate) fn sdirk2_mass_history(
        &self,
        previous: &[f64],
        first: &[f64],
    ) -> Result<Vec<f64>, FluxNewtonError> {
        if previous.len() != self.n_u || first.len() != self.n_u {
            return Err(FluxNewtonError::Invalid(
                "SDIRK2 stage velocity count mismatch",
            ));
        }
        let ratio = (1.0 - GAMMA) / GAMMA;
        let history: Vec<_> = previous
            .iter()
            .zip(first)
            // Fusing the update avoids an overflowing intermediate product
            // when its sum with the old velocity is still representable.
            .map(|(&old, &stage)| ratio.mul_add(stage - old, old))
            .collect();
        if history.iter().any(|value| !value.is_finite()) {
            return Err(FluxNewtonError::Invalid(
                "SDIRK2 second-stage mass history is not representable",
            ));
        }
        Ok(history)
    }
}
