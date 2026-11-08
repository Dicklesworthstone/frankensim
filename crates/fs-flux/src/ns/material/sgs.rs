//! Material rheology and algebraic SGS on one native nonlinear endpoint.
//!
//! Molecular viscosity remains owned by the original source at the actual
//! symmetric strain; Smagorinsky/WALE/Vreman add only modeled stress. The
//! source's physical molecular heating is never relabeled as turbulent heat.

use super::mixed::physical_control;
use super::rheology::RheologyField;
use super::{
    MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxSystem, MaterialNavierStokesBoundary,
    MaterialNavierStokesRheologyStep, MaterialStokesFlowControlSolution, with_kinematic_boundaries,
};
use crate::ns::boundary::{NavierStokesControls, StokesBoundaryError, StokesFlowControlSolution};
use crate::ns::sgs::{PlanarSgs, PlanarSgsError, PlanarSgsField, SgsModel};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::RefCell;

#[path = "sgs/flow_control.rs"]
mod flow_control;
use flow_control::SgsPortControl;

/// Explicit algebraic model and isotropic filter length for each native cell.
/// No width, model coefficient, material parameter or wall correction is inferred.
#[derive(Clone, Copy, Debug)]
pub struct MaterialSgsConfig<'a> {
    pub model: SgsModel,
    /// Positive physical filter widths [m], in the original triangle order.
    pub filter_widths_m: &'a [f64],
}

/// Original source/native failure, or refusal of the shared modeled transport.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialSgsError {
    Material(MaterialFluxError),
    Subgrid(PlanarSgsError),
}

impl core::fmt::Display for MaterialSgsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material SGS flow: {self:?}")
    }
}
impl std::error::Error for MaterialSgsError {}
impl From<MaterialFluxError> for MaterialSgsError {
    fn from(error: MaterialFluxError) -> Self {
        Self::Material(error)
    }
}
impl From<PlanarSgsError> for MaterialSgsError {
    fn from(error: PlanarSgsError) -> Self {
        Self::Subgrid(error)
    }
}

/// One accepted native flow with its exact molecular owners and separate SGS.
pub struct MaterialSgsStep<G> {
    /// Physical pressure/velocity, convergence evidence and endpoint source
    /// samples. Its viscous heating contains ONLY original molecular viscosity.
    /// Boundary reactions, however, belong to the full solved molecular+SGS stress.
    pub flow: MaterialNavierStokesRheologyStep<G>,
    /// Current endpoint strain, eddy/effective viscosities and modeled transfer.
    /// Effective viscosity includes the actual endpoint molecular response.
    pub sgs: PlanarSgsField,
    /// Area-integrated rho*nu_t*gamma^2 [W/m], per out-of-plane span.
    /// This is resolved-to-modeled transfer, not an added heat source or an
    /// accumulated energy. SIP/upwind and time-discretization losses are separate.
    pub modeled_transfer_w_per_m: f64,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Solve steady flow with original molecular rheology plus explicit SGS.
    ///
    /// Both contributions and convection are evaluated at the SAME nonlinear
    /// candidate. The resolver receives the objective cell shear rate, not an
    /// eddy-viscosity-adjusted rate; its density must retain the original fixed
    /// reference. Original Cp, conductivity, provenance and molecular viscosity
    /// are retained on the accepted samples. Temperature/composition captured
    /// by the resolver remain prescribed; mechanical pressure is not an EOS.
    ///
    /// Mixed velocity/physical total traction/Navier-slip boundaries, pressure
    /// datums, backflow policy and endpoint residual gates remain native. The
    /// initial velocity is only a guess; there is no artificial mass term.
    /// No caller history or partial source set is published on failure.
    ///
    /// This is the existing exploratory planar screen. A 2-D velocity field
    /// does not resolve a 3-D turbulent cascade; no transition, wall model,
    /// thermal evolution or whole-scheme energy balance is supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_navier_stokes_mixed_rheology_sgs<G: FluidTransport>(
        &self,
        sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: NavierStokesControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialSgsStep<G>, MaterialSgsError> {
        self.sgs_endpoint(
            sgs,
            force_density,
            boundaries,
            initial_velocity,
            None,
            None,
            controls,
            None,
            resolve,
        )
        .map(|controlled| controlled.flow)
    }

    /// BDF1 with source molecular rheology and SGS converged at the new endpoint.
    /// Uses the same SI conversions, material cache and refusal as the steady
    /// entry. Input velocity is the original physical history and remains intact.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_rheology_sgs_step<G: FluidTransport>(
        &self,
        sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialSgsStep<G>, MaterialSgsError> {
        self.bdf2_navier_stokes_mixed_rheology_sgs_step(
            sgs,
            force_density,
            boundaries,
            previous_velocity,
            None,
            timestep_s,
            controls,
            resolve,
        )
    }

    /// Equal-step BDF2, or original BDF1 startup with no older history.
    /// Only inertia uses the two-history combination; both constitutive models
    /// use the actual candidate. Endpoint powers are rates, not integrated heat.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_sgs_step<G: FluidTransport>(
        &self,
        sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_velocity: Option<&[f64]>,
        timestep_s: f64,
        controls: NavierStokesControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialSgsStep<G>, MaterialSgsError> {
        self.bdf2_navier_stokes_mixed_rheology_sgs_variable_step(
            sgs,
            force_density,
            boundaries,
            previous_velocity,
            older_velocity.map(|older| (older, timestep_s)),
            timestep_s,
            controls,
            resolve,
        )
    }

    /// Unequal-interval BDF2 with the actual older physical interval.
    /// The original time owner admits both histories and its step-growth bound;
    /// this closure does not alter time coefficients or extrapolate materials.
    /// Neither velocity history nor the reference material is mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_sgs_variable_step<G: FluidTransport>(
        &self,
        sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_velocity: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: NavierStokesControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialSgsStep<G>, MaterialSgsError> {
        self.sgs_endpoint(
            sgs,
            force_density,
            boundaries,
            previous_velocity,
            older_velocity,
            Some(timestep_s),
            controls,
            None,
            resolve,
        )
        .map(|controlled| controlled.flow)
    }

    #[allow(clippy::too_many_arguments)]
    fn sgs_endpoint<G: FluidTransport>(
        &self,
        controls_sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        older_velocity: Option<(&[f64], f64)>,
        timestep_s: Option<f64>,
        controls: NavierStokesControls,
        port_control: Option<SgsPortControl<'_>>,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialSgsStep<G>>, MaterialSgsError> {
        let count = self.system.mesh.tris.len();
        let density = self.transport.density_kg_m3();
        // Setup-only storage for the existing authored model. Its frozen
        // reference viscosity never enters the nonlinear composition below.
        let closure = PlanarSgs::with_model(
            self.system,
            controls_sgs.model,
            controls_sgs.filter_widths_m,
            &vec![self.transport.kinematic_viscosity_m2_s(); count],
        )?;
        let field = RheologyField::new(&self.system.mesh.centroids, &self.transport, resolve);
        let subgrid_error = RefCell::new(None);
        let result = {
            let prepare = |velocity: &[f64]| {
                let mut effective = Vec::with_capacity(count);
                for cell in 0..count {
                    let (rate, eddy) =
                        closure.modeled_response(cell, velocity).map_err(|error| {
                            *subgrid_error.borrow_mut() = Some(error);
                            StokesBoundaryError::InvalidInput("invalid material SGS response")
                        })?;
                    let molecular = field.viscosity(cell, rate);
                    if !field.is_valid() {
                        return Err(StokesBoundaryError::InvalidInput(
                            "invalid molecular material response",
                        ));
                    }
                    let total = molecular + eddy;
                    if !total.is_finite() {
                        *subgrid_error.borrow_mut() =
                            Some(PlanarSgsError::UnrepresentableTransport { cell });
                        return Err(StokesBoundaryError::InvalidInput(
                            "unrepresentable combined material SGS viscosity",
                        ));
                    }
                    effective.push(total);
                }
                Ok(effective)
            };
            let force = |point| {
                let value = force_density(point);
                [value[0] / density, value[1] / density]
            };
            with_kinematic_boundaries(density, boundaries, |native| {
                let source_is_valid = || field.is_valid() && subgrid_error.borrow().is_none();
                let uncontrolled = |flow| StokesFlowControlSolution {
                    flow,
                    pressure_offsets_m2_s2: Vec::new(),
                    outward_volume_rates_m2_s: Vec::new(),
                    maximum_flow_residual_m2_s: 0.0,
                };
                if let Some(control) = port_control {
                    self.system
                        .solve_navier_stokes_mixed_prepared_flow_control_checked(
                            self.params,
                            &prepare,
                            &|_, _, value| value,
                            &force,
                            native,
                            initial_velocity,
                            older_velocity,
                            timestep_s,
                            controls,
                            control.ports,
                            control.targets,
                            control.absolute_tolerance_m2_s,
                            &source_is_valid,
                        )
                } else if let Some(dt) = timestep_s {
                    self.system
                        .bdf2_navier_stokes_mixed_prepared_rheology_variable_step_checked(
                            self.params,
                            &prepare,
                            &|_, _, value| value,
                            &force,
                            native,
                            initial_velocity,
                            older_velocity,
                            dt,
                            controls,
                            &source_is_valid,
                        )
                        .map(uncontrolled)
                } else {
                    self.system
                        .solve_navier_stokes_mixed_prepared_rheology_checked(
                            self.params,
                            &prepare,
                            &|_, _, value| value,
                            &force,
                            native,
                            initial_velocity,
                            controls,
                            &source_is_valid,
                        )
                        .map(uncontrolled)
                }
            })
            .map_err(MaterialFluxError::Boundary)
        };
        // Native acceptance reassembles the endpoint. Keep those original
        // source samples and do not query a provider again for diagnostics.
        let resolved = field.finish(result);
        if let Some(error) = subgrid_error.into_inner() {
            return Err(error.into());
        }
        let (endpoint, samples) = resolved?;
        let sgs = closure
            .field_with_molecular(&endpoint.flow.solution.field.x[..self.system.n_u], &|cell| {
                samples[cell].transport().kinematic_viscosity_m2_s()
            })?;
        let modeled_transfer_w_per_m = density * sgs.integrated_dissipation_m4_s3;
        if !modeled_transfer_w_per_m.is_finite()
            || (sgs.integrated_dissipation_m4_s3 > 0.0 && modeled_transfer_w_per_m == 0.0)
            || sgs.cells.iter().any(|cell| {
                let power = density * cell.dissipation_m2_s3;
                !power.is_finite() || (cell.dissipation_m2_s3 > 0.0 && power == 0.0)
            })
        {
            return Err(MaterialFluxError::Invalid("unrepresentable physical SGS transfer").into());
        }
        let controlled = physical_control(endpoint, density, |endpoint| {
            self.material_navier_stokes_rheology_step(endpoint, samples)
        })?;
        Ok(MaterialStokesFlowControlSolution {
            flow: MaterialSgsStep {
                flow: controlled.flow,
                sgs,
                modeled_transfer_w_per_m,
            },
            pressure_offsets_pa: controlled.pressure_offsets_pa,
            outward_volume_rates_m2_s: controlled.outward_volume_rates_m2_s,
            maximum_flow_residual_m2_s: controlled.maximum_flow_residual_m2_s,
        })
    }
}
