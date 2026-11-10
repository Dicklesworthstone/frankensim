//! Original material tangents in mixed-boundary convective Newton flow.

use super::{
    FluxNewtonControls, MaterialError, MaterialFluxError, MaterialFluxRheologyPoint,
    MaterialFluxSystem, MaterialFluxViscosityTangent, RheologyField, retain_tangent,
};
use crate::material::{
    MaterialDissipatedHeat, MaterialNavierStokesBoundary, MaterialNavierStokesRheologyStep,
    with_kinematic_boundaries,
};
use crate::ns::boundary::NavierStokesStep;
use crate::ns::linearization::{FluxViscosityTangent, ViscosityLaw};
use crate::ns::material::heat::DissipatedHeatAccumulator;
use crate::ns::newton::{FluxSdirk2Error, FluxSdirk2Stage};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::dd::Dd;
use std::cell::RefCell;

/// Adaptive physical-time advancement on the original mixed SDIRK2 owner.
#[path = "mixed/adaptive.rs"]
pub mod adaptive;

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

/// Accepted mixed-boundary SDIRK2 stages with their original material sources.
pub struct MaterialNavierStokesSdirk2Step<G> {
    /// Internal-stage flow, physical reactions and instantaneous molecular heat.
    pub first_stage: MaterialNavierStokesNewtonStep<G>,
    /// Stiffly accurate endpoint, including its own original material sources.
    pub endpoint: MaterialNavierStokesNewtonStep<G>,
    /// Original-stage bulk and Navier-slip production integrated over this
    /// physical step [J/m], with positive SDIRK2 weights (1−gamma, gamma).
    /// No thermal deposition or wall/fluid slip-heat partition is inferred.
    pub dissipated_heat: MaterialDissipatedHeat,
    /// Signed kinetic quadrature defect of this SDIRK2 time step [J/m].
    /// With d1=Y1−Un and d2=Y2−Un−((1−gamma)/gamma)*d1, this is
    /// rho*(K(d2)−K(d1)), using the original finite-element mass quadrature.
    /// It can have EITHER sign. It is not physical heat, nonnegative numerical
    /// damping, spatial discretization loss or a full energy/passivity balance.
    pub temporal_kinetic_defect_j_per_m: f64,
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
            None,
            controls,
            resolve,
        )
    }

    /// Advance mixed-boundary Navier--Stokes with the original two-stage,
    /// second-order, stiffly accurate SDIRK2 method and constitutive Newton
    /// owner. No older interval or first-order startup step is required.
    ///
    /// The first boundary slice belongs to `t_n + gamma * timestep_s`, the
    /// second to `t_n + timestep_s`, with `gamma = 1 - 1/sqrt(2)`. Each may
    /// supply its own physical traction, reservoir velocity, wall velocity and
    /// Navier resistance. Force and material callbacks receive the same stage.
    /// All external fields must remain pure and fixed within each Newton solve.
    /// Forces are N/m³, traction Pa, velocity m/s, resistance Pa s/m, time s,
    /// and the viscosity derivative Pa s².
    ///
    /// Both stages use the same reference density, original material cache,
    /// analytic constitutive/advection Jacobian, bounded line search and mixed
    /// reaction owner. Finite stage mass anchors all-open components. Stage two
    /// starts from the actual first-stage velocity; the SDIRK mass combination
    /// enters only its inertial load and is never a constitutive initial state.
    ///
    /// Each accepted stage retains its own physical pressure, source owners,
    /// tangents, reactions and `mu*gamma²` heating. Positive SDIRK2 weights
    /// (1−gamma, gamma) integrate those original stage powers into separate
    /// cell bulk and boundary-slip heat [J/m]. No source is queried again.
    /// Kinetic transport, wall work and numerical damping remain separate.
    /// This adds no caloric evolution, slip-heat partition or full energy balance.
    /// Temperature/composition may differ between stages but density remains
    /// fixed. Work limits apply per stage; temporal order requires smooth stage
    /// data and resolved spatial/nonlinear errors. Linear L-stability is not a
    /// nonlinear positivity or passivity theorem.
    ///
    /// # Errors
    /// Invalid history/timestep refuses before any stage callback. Errors retain
    /// the failing stage and original material/numerical cause, including source
    /// loss during a trial or final reaction assembly. A second-stage refusal
    /// publishes no partial pair; caller history and accepted sources are never
    /// modified. There is no hidden timestep clipping or fallback boundary law.
    #[allow(clippy::too_many_arguments)]
    pub fn sdirk2_navier_stokes_mixed_rheology_newton_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn(FluxSdirk2Stage, [f64; 2]) -> [f64; 2],
        first_stage_boundaries: &[MaterialNavierStokesBoundary<'_>],
        endpoint_boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
        mut resolve: impl FnMut(
            FluxSdirk2Stage,
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
    ) -> Result<MaterialNavierStokesSdirk2Step<G>, FluxSdirk2Error<MaterialFluxError>> {
        let stage_dt = self
            .system
            .sdirk2_stage_timestep(previous_velocity, timestep_s)
            .map_err(|error| FluxSdirk2Error {
                stage: None,
                error: MaterialFluxError::Newton(error),
            })?;
        let mut solve = |stage, initial: &[f64], history: &[f64]| {
            let boundaries = match stage {
                FluxSdirk2Stage::First => first_stage_boundaries,
                FluxSdirk2Stage::Endpoint => endpoint_boundaries,
            };
            self.mixed_rheology_newton(
                &|point| force_density(stage, point),
                boundaries,
                initial,
                None,
                None,
                Some((stage_dt, history)),
                controls,
                |point, reference| resolve(stage, point, reference),
            )
            .map_err(|error| FluxSdirk2Error {
                stage: Some(stage),
                error,
            })
        };
        let first_stage = solve(FluxSdirk2Stage::First, previous_velocity, previous_velocity)?;
        let mut heat = DissipatedHeatAccumulator::new(self.system.mesh.areas.len());
        heat.add_stage(
            &self.system.mesh.areas,
            &first_stage.flow.step.solution,
            timestep_s - stage_dt,
        )
        .map_err(|error| FluxSdirk2Error {
            stage: Some(FluxSdirk2Stage::First),
            error,
        })?;
        let first_velocity = &first_stage.flow.step.solution.discrete().x[..self.system.n_u];
        let mut history = self
            .system
            .sdirk2_mass_history(previous_velocity, first_velocity)
            .map_err(|error| FluxSdirk2Error {
                stage: Some(FluxSdirk2Stage::Endpoint),
                error: MaterialFluxError::Newton(error),
            })?;
        let endpoint = solve(FluxSdirk2Stage::Endpoint, first_velocity, &history)?;
        heat.add_stage(
            &self.system.mesh.areas,
            &endpoint.flow.step.solution,
            stage_dt,
        )
        .map_err(|error| FluxSdirk2Error {
            stage: Some(FluxSdirk2Stage::Endpoint),
            error,
        })?;
        let dissipated_heat = heat.finish(timestep_s).map_err(|error| FluxSdirk2Error {
            stage: Some(FluxSdirk2Stage::Endpoint),
            error,
        })?;
        // Reuse the consumed mass-history allocation for the difference of
        // increments. Form the signed quadratic difference as a mass pairing
        // instead of subtracting two nearly equal kinetic energies.
        let endpoint_velocity = endpoint.flow.step.solution.velocity_dofs();
        let mut sum = Vec::with_capacity(self.system.n_u);
        let ratio = Dd::from_f64(
            (1.0 - FluxSdirk2Stage::First.time_fraction()) / FluxSdirk2Stage::First.time_fraction(),
        );
        for i in 0..self.system.n_u {
            let first = Dd::from_f64(first_velocity[i]) - Dd::from_f64(previous_velocity[i]);
            let second = Dd::from_f64(endpoint_velocity[i])
                - Dd::from_f64(previous_velocity[i])
                - ratio * first;
            history[i] = (second - first).to_f64();
            sum.push((second + first).to_f64());
        }
        let pairing = self.system.mass_inner_product(&history, &sum);
        let temporal_kinetic_defect_j_per_m = 0.5 * self.transport.density_kg_m3() * pairing;
        if !temporal_kinetic_defect_j_per_m.is_finite()
            || (pairing != 0.0 && temporal_kinetic_defect_j_per_m == 0.0)
        {
            return Err(FluxSdirk2Error {
                stage: Some(FluxSdirk2Stage::Endpoint),
                error: MaterialFluxError::Invalid("unrepresentable SDIRK2 kinetic time defect"),
            });
        }
        Ok(MaterialNavierStokesSdirk2Step {
            first_stage,
            endpoint,
            dissipated_heat,
            temporal_kinetic_defect_j_per_m,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn mixed_rheology_newton<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        timestep_s: Option<f64>,
        older_history: Option<(&[f64], f64)>,
        stage_mass: Option<(f64, &[f64])>,
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
            if let Some((stage_timestep, history)) = stage_mass {
                self.system.mixed_constitutive_newton_endpoint(
                    self.params,
                    ViscosityLaw::Rate(&tangent),
                    &acceleration,
                    primitive,
                    initial_velocity,
                    stage_timestep.recip(),
                    Some(history),
                    controls,
                    &|| field.is_valid(),
                )
            } else if let Some(timestep) = timestep_s {
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
