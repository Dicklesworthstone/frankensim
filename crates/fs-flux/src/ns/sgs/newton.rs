//! Consistent Smagorinsky stress on the existing damped Newton owner.
//!
//! For nu_eff = nu + (C_s Delta)^2 gamma, the positive-strain tangent is
//! d(nu_eff)/d(gamma) = (C_s Delta)^2. The original residual/Jacobian owner
//! differentiates the complete cell stress, harmonic SIP weights, penalty,
//! boundary loads and convection. The constitutive law is neither frozen nor
//! regularized to obtain convergence.
//!
//! At exactly zero strain, use the zero generalized derivative of the strain
//! norm. The modeled volume stress |D|D has zero derivative there; viscosity-
//! dependent DG penalties can remain nonsmooth when a tangential jump is
//! present. This is a semismooth selection, not a classical derivative claim
//! at those states. WALE and Vreman depend on the full gradient and cannot use
//! a scalar strain tangent; they retain their existing Picard entry points.

use super::{PlanarSgs, PlanarSgsError, PlanarSgsField, SgsModel};
use crate::ns::FluxParams;
use crate::ns::linearization::FluxViscosityTangent;
use crate::ns::newton::{FluxNewtonControls, FluxNewtonError, FluxNewtonStep};

/// Accepted native Newton flow and SGS diagnostics at the same endpoint.
pub struct PlanarSgsNewtonStep {
    /// Original velocity/pressure, consistent residual and nonlinear work.
    pub step: FluxNewtonStep,
    /// Instantaneous modeled transfer; no extra heat is added to the flow.
    pub sgs: PlanarSgsField,
}

impl From<FluxNewtonError> for PlanarSgsError {
    fn from(error: FluxNewtonError) -> Self {
        Self::Newton(error)
    }
}

impl PlanarSgs<'_, '_> {
    /// Solve true steady Smagorinsky Navier--Stokes flow with damped Newton.
    ///
    /// Every boundary prescribes velocity. Initial velocity is only a nonlinear
    /// guess: no artificial timestep or mass is introduced. The existing sparse
    /// constitutive/convection Jacobian and Armijo line search reduce the actual
    /// residual; both original and row-scaled residuals must converge. Params.nu
    /// is unused: the prepared molecular field remains authoritative.
    ///
    /// The analytic Smagorinsky tangent includes the model's influence on SIP
    /// fluxes and boundary loads. Zero strain uses the documented semismooth
    /// selection. Force and boundary data must be pure and fixed during the
    /// solve. This retains the Estimated planar screen's physical limitations;
    /// a better nonlinear solve does not establish a three-dimensional cascade.
    ///
    /// # Errors
    /// Refuses WALE/Vreman (no scalar tangent), unrepresentable tangents,
    /// invalid boundary/control data and the native Newton/endpoint failures.
    pub fn solve_steady_newton(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        velocity_boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        initial_velocity: &[f64],
        controls: FluxNewtonControls,
    ) -> Result<PlanarSgsNewtonStep, PlanarSgsError> {
        let tangents = self.smagorinsky_tangents()?;
        let step = self.system.steady_rheology_newton(
            params,
            &|cell, _, rate| self.newton_response(&tangents, cell, rate),
            force,
            velocity_boundary,
            initial_velocity,
            true,
            controls,
        )?;
        self.newton_endpoint(step)
    }

    /// Advance fully implicit Smagorinsky momentum with damped Newton BDF1.
    ///
    /// Stress and convection belong to the accepted endpoint. Every boundary
    /// prescribes velocity; data describe the requested new time. The caller's
    /// previous velocity is immutable on success and refusal. The constitutive
    /// and boundary qualifications of [`Self::solve_steady_newton`] apply.
    pub fn bdf1_newton_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        velocity_boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<PlanarSgsNewtonStep, PlanarSgsError> {
        self.bdf2_newton_variable_step(
            params,
            force,
            velocity_boundary,
            previous_velocity,
            None,
            timestep_s,
            controls,
        )
    }

    /// Advance fully implicit Smagorinsky momentum with unequal-step BDF2.
    ///
    /// The older velocity is paired with its ACTUAL interval to the previous
    /// velocity [s]. The new/previous ratio must lie in (0,2]; an equal interval
    /// gives ordinary BDF2, and None is exactly BDF1 startup. The native history
    /// owner combines only the mass load: the Newton initial iterate is the
    /// actual previous velocity. Stress and convection converge at the endpoint.
    ///
    /// Smooth resolved motion retains second-order time accuracy when spatial
    /// and nonlinear errors are smaller. Endpoint modeled transfer remains an
    /// instantaneous rate, not BDF2's two-history energy or accumulated heat.
    /// Prescribed-velocity boundaries and the Smagorinsky-only tangent scope
    /// follow [`Self::solve_steady_newton`]; neither history is modified.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_newton_variable_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        velocity_boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<PlanarSgsNewtonStep, PlanarSgsError> {
        let tangents = self.smagorinsky_tangents()?;
        let step = self.system.bdf2_rheology_newton_variable_step(
            params,
            &|cell, _, rate| self.newton_response(&tangents, cell, rate),
            force,
            velocity_boundary,
            previous_velocity,
            older_history,
            timestep_s,
            controls,
        )?;
        self.newton_endpoint(step)
    }

    // Prepared once per solve, not recomputed at every DG quadrature point.
    fn smagorinsky_tangents(&self) -> Result<Vec<f64>, PlanarSgsError> {
        if !matches!(self.model, SgsModel::Smagorinsky { .. }) {
            return Err(PlanarSgsError::InvalidInput(
                "Newton SGS requires Smagorinsky's scalar strain tangent",
            ));
        }
        self.closure_lengths_m
            .iter()
            .enumerate()
            .map(|(cell, &length)| {
                let tangent = length * length;
                if !tangent.is_finite() || (length > 0.0 && tangent == 0.0) {
                    return Err(PlanarSgsError::UnrepresentableTransport { cell });
                }
                Ok(tangent)
            })
            .collect()
    }

    fn newton_response(&self, tangents: &[f64], cell: usize, rate: f64) -> FluxViscosityTangent {
        // Invalid trial arithmetic is a native constitutive refusal, allowing
        // the ORIGINAL line search to backtrack. Never poison the whole solve
        // merely because a rejected Newton trial overflowed its SGS response.
        let nu = self
            .viscosities(cell, rate)
            .map_or(f64::NAN, |(_, effective)| effective);
        FluxViscosityTangent {
            nu,
            derivative: if rate == 0.0 { 0.0 } else { tangents[cell] },
        }
    }

    fn newton_endpoint(&self, step: FluxNewtonStep) -> Result<PlanarSgsNewtonStep, PlanarSgsError> {
        let sgs = self.field(&step.solution.x[..self.system.n_u])?;
        Ok(PlanarSgsNewtonStep { step, sgs })
    }
}
