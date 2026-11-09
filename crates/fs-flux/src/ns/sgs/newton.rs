//! Consistent Smagorinsky, WALE and Vreman stress on the original Newton owner.
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
//! at those states. WALE and Vreman use the shared analytic full-gradient
//! derivative, including rotation, prepared once per triangle per assembly.
//! Their zero-gradient and Vreman rank-one selections follow the shared kernel;
//! they do not establish a classical derivative at nonsmooth states.

use super::{PlanarSgs, PlanarSgsError, PlanarSgsField, SgsModel};
use crate::ns::FluxParams;
use crate::ns::boundary::NavierStokesBoundary;
use crate::ns::linearization::{FluxViscosityGradientTangent, FluxViscosityTangent, ViscosityLaw};
use crate::ns::newton::{FluxMixedNewtonStep, FluxNewtonControls, FluxNewtonError, FluxNewtonStep};
use fs_lbm::sgs3::eddy_viscosity_tangent;

/// Accepted native Newton flow and SGS diagnostics at the same endpoint.
pub struct PlanarSgsNewtonStep {
    /// Original velocity/pressure, consistent residual and nonlinear work.
    pub step: FluxNewtonStep,
    /// Instantaneous modeled transfer; no extra heat is added to the flow.
    pub sgs: PlanarSgsField,
}

/// Accepted mixed-boundary Newton flow and its instantaneous SGS diagnostics.
pub struct PlanarSgsMixedNewtonStep {
    /// Original pressure datum, physical reactions, boundary policies and work.
    pub step: FluxMixedNewtonStep,
    /// Modeled transfer evaluated at this same accepted velocity.
    pub sgs: PlanarSgsField,
}

impl From<FluxNewtonError> for PlanarSgsError {
    fn from(error: FluxNewtonError) -> Self {
        Self::Newton(error)
    }
}

impl PlanarSgs<'_, '_> {
    /// Solve true steady algebraic SGS Navier--Stokes flow with damped Newton.
    ///
    /// Every boundary prescribes velocity. Initial velocity is only a nonlinear
    /// guess: no artificial timestep or mass is introduced. The existing sparse
    /// constitutive/convection Jacobian and Armijo line search reduce the actual
    /// residual; both original and row-scaled residuals must converge. Params.nu
    /// is unused: the prepared molecular field remains authoritative.
    ///
    /// The analytic tangent includes the model's influence on SIP fluxes and
    /// boundary loads. WALE/Vreman use the full velocity gradient, not a fitted
    /// strain-only law. Degenerate states use the shared generalized selections.
    /// Force and boundary data must be pure and fixed during the
    /// solve. This retains the Estimated planar screen's physical limitations;
    /// a better nonlinear solve does not establish a three-dimensional cascade.
    ///
    /// # Errors
    /// Refuses unrepresentable tangents, invalid boundary/control data and
    /// the native Newton/endpoint failures. All three local models are supported.
    pub fn solve_steady_newton(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        velocity_boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        initial_velocity: &[f64],
        controls: FluxNewtonControls,
    ) -> Result<PlanarSgsNewtonStep, PlanarSgsError> {
        let step = self.with_newton_law(|law| {
            self.system.constitutive_newton_solve_admitted_checked(
                params,
                law,
                force,
                velocity_boundary,
                initial_velocity,
                true,
                0.0,
                None,
                controls,
                &|| true,
                None,
            )
        })?;
        self.newton_endpoint(step)
    }

    /// Advance fully implicit local SGS momentum with damped Newton BDF1.
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

    /// Advance fully implicit local SGS momentum with unequal-step BDF2.
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
    /// Prescribed-velocity boundaries and the constitutive tangent scope
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
        let step = self.with_newton_law(|law| {
            self.system.bdf2_constitutive_newton_variable_step_checked(
                params,
                law,
                force,
                velocity_boundary,
                previous_velocity,
                older_history,
                timestep_s,
                controls,
                &|| true,
            )
        })?;
        self.newton_endpoint(step)
    }

    /// Resolve steady local SGS flow with velocity, slip and open boundaries.
    ///
    /// Reuses the SAME exact stress/convection Newton owner. Literal traction,
    /// dissipative backflow and reservoir inflow retain their original distinct
    /// meanings. Boundary reactions use accepted effective viscosity; no heat
    /// source is introduced. Full-gradient WALE/Vreman tangents retain the same
    /// closure values and Estimated planar scope; generalized selections at
    /// degenerate gradients do not imply a unique classical derivative.
    pub fn solve_steady_mixed_newton(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: FluxNewtonControls,
    ) -> Result<PlanarSgsMixedNewtonStep, PlanarSgsError> {
        let step = self.with_newton_law(|law| {
            self.system.mixed_constitutive_newton_endpoint(
                params,
                law,
                force,
                boundaries,
                initial_velocity,
                0.0,
                None,
                controls,
                &|| true,
            )
        })?;
        self.mixed_newton_endpoint(step)
    }

    /// Fully implicit BDF1 local SGS momentum on admitted mixed boundaries.
    /// Previous velocity is immutable; all external data belong to the endpoint.
    pub fn bdf1_mixed_newton_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<PlanarSgsMixedNewtonStep, PlanarSgsError> {
        self.bdf2_mixed_newton_variable_step(
            params,
            force,
            boundaries,
            previous_velocity,
            None,
            timestep_s,
            controls,
        )
    }

    /// Unequal-step BDF2 local SGS momentum with original mixed-boundary work.
    ///
    /// The older velocity is paired with its actual previous interval. The
    /// native owner admits ratios in (0,2], and None is BDF1 startup. The mass
    /// history is combined once; velocity, viscosity and convection converge
    /// together at the endpoint. Neither history nor boundary law is modified.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_mixed_newton_variable_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<PlanarSgsMixedNewtonStep, PlanarSgsError> {
        let step = self.with_newton_law(|law| {
            self.system.bdf2_mixed_constitutive_newton_variable_step(
                params,
                law,
                force,
                boundaries,
                previous_velocity,
                older_history,
                timestep_s,
                controls,
                &|| true,
            )
        })?;
        self.mixed_newton_endpoint(step)
    }

    fn mixed_newton_endpoint(
        &self,
        step: FluxMixedNewtonStep,
    ) -> Result<PlanarSgsMixedNewtonStep, PlanarSgsError> {
        let sgs = self.field(&step.solution.field.x[..self.system.n_u])?;
        Ok(PlanarSgsMixedNewtonStep { step, sgs })
    }

    fn with_newton_law<T>(
        &self,
        solve: impl FnOnce(ViscosityLaw<'_>) -> Result<T, FluxNewtonError>,
    ) -> Result<T, PlanarSgsError> {
        if matches!(self.model, SgsModel::Smagorinsky { .. }) {
            let tangents = self.smagorinsky_tangents()?;
            Ok(solve(ViscosityLaw::Rate(&|cell, _, rate| {
                self.newton_response(&tangents, cell, rate)
            }))?)
        } else {
            Ok(solve(ViscosityLaw::Gradient(&|cell, gradient| {
                self.gradient_newton_response(cell, gradient)
            }))?)
        }
    }

    fn gradient_newton_response(
        &self,
        cell: usize,
        gradient: [[f64; 2]; 2],
    ) -> FluxViscosityGradientTangent {
        let embedded = [
            [gradient[0][0], gradient[0][1], 0.0],
            [gradient[1][0], gradient[1][1], 0.0],
            [0.0; 3],
        ];
        match eddy_viscosity_tangent(self.model, self.closure_lengths_m[cell], embedded) {
            Ok(value) => FluxViscosityGradientTangent {
                nu: self.molecular_viscosity_m2_s[cell] + value.viscosity,
                derivative: core::array::from_fn(|i| {
                    core::array::from_fn(|j| value.derivative[i][j])
                }),
            },
            // Rejected trial arithmetic stays inside the original line search.
            // No shared error latch may poison a later admissible damped trial.
            Err(_) => FluxViscosityGradientTangent {
                nu: f64::NAN,
                derivative: [[f64::NAN; 2]; 2],
            },
        }
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
