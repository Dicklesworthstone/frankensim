//! Original material tangents in mixed-boundary convective Newton flow.

use super::{
    FluxNewtonControls, MaterialError, MaterialFluxError, MaterialFluxRheologyPoint,
    MaterialFluxSystem, MaterialFluxViscosityTangent, RheologyField, retain_tangent,
};
use crate::material::{
    MaterialNavierStokesBoundary, MaterialNavierStokesRheologyStep, with_kinematic_boundaries,
};
use crate::ns::boundary::NavierStokesStep;
use crate::ns::linearization::FluxViscosityTangent;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::RefCell;

/// Accepted mixed Navier--Stokes Newton endpoint and original material owners.
pub struct MaterialNavierStokesNewtonStep<G> {
    /// Physical flow, endpoint sources, pressure datums, reactions and bulk heat.
    /// `flow.step.iterations` counts Newton corrections, not Picard updates.
    /// Viscous/slip losses are endpoint powers, not integrated heat deposition.
    pub flow: MaterialNavierStokesRheologyStep<G>,
    /// Line-search candidates rejected by the native Newton owner.
    pub backtracks: usize,
    derivatives_pa_s2: Vec<f64>,
}

impl<G> MaterialNavierStokesNewtonStep<G> {
    /// Original analytic `dmu/dgamma` [Pa s²] at the accepted cell rates, in
    /// the same order as `flow.samples()`. This does not requery any source.
    #[must_use]
    pub fn viscosity_derivatives_pa_s2(&self) -> &[f64] {
        &self.derivatives_pa_s2
    }
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Solve steady mixed-boundary Navier--Stokes with original viscosity and
    /// analytic rate derivatives. The initial velocity is only a guess; no
    /// transient mass or fabricated history enters the steady equations.
    ///
    /// Force density is N/m³, total traction Pa, velocity m/s and Navier
    /// resistance Pa s/m. The SAME fixed source density scales viscosity,
    /// its derivative, convection, force, stress, pressure and wall resistance.
    /// Natural traction, dissipative backflow and reservoir inflow preserve
    /// their original separate policies. No reservoir imposes its flow rate.
    ///
    /// Captured temperature/composition and boundary data are fixed and pure.
    /// Original sources supply `mu` and `dmu/dgamma` at the same cell/rate;
    /// the cache retains the accepted owners for physical heat and reporting.
    /// Source failures remain terminal during line search and final reaction
    /// assembly. No fitted or fallback-viscosity endpoint can be published.
    ///
    /// # Errors
    /// Preserves the original material, tangent, density, steady anchoring,
    /// boundary and Newton work/convergence refusals. Initial state is immutable.
    pub fn solve_navier_stokes_mixed_rheology_newton<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: FluxNewtonControls,
        resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialNavierStokesNewtonStep<G>, MaterialFluxError> {
        self.mixed_rheology_newton(
            force_density,
            boundaries,
            initial_velocity,
            None,
            None,
            controls,
            resolve,
        )
    }

    /// Fully implicit physical BDF1 with mixed boundaries and exact material
    /// tangents. Endpoint stress and convection converge in the SAME original
    /// Newton solve. Finite mass anchors velocity even on all-open components.
    ///
    /// The supplied time is in seconds and previous velocity in m/s. All
    /// external fields describe the new time, with the same constant reference
    /// density at both levels. Accepted original sources determine endpoint
    /// pressure, reactions and `mu*gamma²` bulk heating; no caloric update or
    /// interpretation of mechanical pressure as an EOS state is inferred.
    ///
    /// # Errors
    /// Preserves the source and Newton refusals of the steady entry and the
    /// native history/timestep checks. Failure changes neither history nor any
    /// accepted source. No timestep clipping or hidden substepping is applied.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_rheology_newton_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
        resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialNavierStokesNewtonStep<G>, MaterialFluxError> {
        self.bdf2_navier_stokes_mixed_rheology_newton_variable_step(
            force_density,
            boundaries,
            previous_velocity,
            None,
            timestep_s,
            controls,
            resolve,
        )
    }

    /// Unequal-step BDF2 for the same source-owned mixed Newton flow.
    ///
    /// `older_history` pairs the older velocity with its interval to the
    /// previous velocity. The native BDF owner admits new/previous step ratios
    /// in (0,2]; `None` is exact BDF1 startup. It combines only the inertial
    /// mass load. The actual previous velocity remains the Newton/material
    /// initial state, rather than an artificial combined-history shear field.
    ///
    /// Both histories and reference density remain fixed. Returned heat and
    /// reactions are instantaneous endpoint rates; neither a two-history BDF
    /// energy identity nor thermal deposition is claimed. Temperature evolution,
    /// variable-density flow and viscoelastic memory remain outside this API.
    ///
    /// # Errors
    /// Invalid histories/intervals refuse before material or load callbacks.
    /// Original source failures take precedence over numerical refusals; all
    /// steady/BDF1 constitutive, boundary and Newton work bounds remain active.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_newton_variable_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: FluxNewtonControls,
        resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialNavierStokesNewtonStep<G>, MaterialFluxError> {
        self.mixed_rheology_newton(
            force_density,
            boundaries,
            previous_velocity,
            Some(timestep_s),
            older_history,
            controls,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn mixed_rheology_newton<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        timestep_s: Option<f64>,
        older_history: Option<(&[f64], f64)>,
        controls: FluxNewtonControls,
        mut resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialNavierStokesNewtonStep<G>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let derivatives = RefCell::new(vec![0.0; self.system.mesh.tris.len()]);
        let field = RheologyField::new(
            &self.system.mesh.centroids,
            &self.transport,
            |point, reference| {
                retain_tangent(resolve(point, reference)?, point, density, &derivatives)
            },
        );
        let result = with_kinematic_boundaries(density, boundaries, |primitive| {
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
            let acceleration = |point| {
                let physical = force_density(point);
                [physical[0] / density, physical[1] / density]
            };
            if let Some(timestep) = timestep_s {
                self.system
                    .bdf2_navier_stokes_mixed_rheology_newton_variable_step_checked(
                        self.params,
                        &tangent,
                        &acceleration,
                        primitive,
                        initial_velocity,
                        older_history,
                        timestep,
                        controls,
                        &|| field.is_valid(),
                    )
            } else {
                self.system
                    .solve_navier_stokes_mixed_rheology_newton_checked(
                        self.params,
                        &tangent,
                        &acceleration,
                        primitive,
                        initial_velocity,
                        controls,
                        &|| field.is_valid(),
                    )
            }
        })
        .map_err(MaterialFluxError::Newton);
        let (solved, samples) = field.finish(result)?;
        let flow = self.material_navier_stokes_rheology_step(
            NavierStokesStep {
                solution: solved.solution,
                iterations: solved.iterations,
                relative_update: solved.relative_update,
                nonlinear_residual: solved.nonlinear_residual,
                open_boundary_policies: solved.open_boundary_policies,
            },
            samples,
        )?;
        Ok(MaterialNavierStokesNewtonStep {
            flow,
            backtracks: solved.backtracks,
            derivatives_pa_s2: derivatives.into_inner(),
        })
    }
}
