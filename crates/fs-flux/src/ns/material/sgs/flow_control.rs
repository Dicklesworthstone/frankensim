//! Required physical port pressures with native molecular and SGS stress.
//!
//! The integral constraints, pressure responses, field relaxation and endpoint
//! residuals belong to the existing native nonlinear flow-control owner. This
//! entry only composes its prepared constitutive response with the original
//! material source and existing full-gradient SGS kernels.

use super::{MaterialSgsConfig, MaterialSgsError, MaterialSgsStep};
use crate::material::{
    MaterialFluxRheologyPoint, MaterialFluxSystem, MaterialNavierStokesBoundary,
    MaterialStokesFlowControlSolution,
};
use crate::ns::boundary::NavierStokesControls;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

#[derive(Clone, Copy)]
pub(super) struct SgsPortControl<'a> {
    pub(super) ports: &'a [&'a [usize]],
    pub(super) targets: &'a [f64],
    pub(super) absolute_tolerance_m2_s: f64,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Solve the steady pressures required by integral port flows with SGS.
    ///
    /// Each target is a signed OUTWARD volume rate per out-of-plane depth
    /// [m²/s]; multiply by the physical depth for m³/s. Ports contain disjoint
    /// original traction edges, and each controlled component must retain an
    /// uncontrolled traction datum. One uniform pressure offset is added to
    /// each port's original stress. The velocity profile remains a resolved
    /// unknown, including the original tangential stress and backflow policy.
    ///
    /// Native Picard updates converge convection, molecular rheology, SGS stress
    /// and pressure offsets together. Every pressure response in an update uses
    /// that same prepared viscosity. Pressure offsets relax with velocity and
    /// pressure; final acceptance independently checks both pressure-inclusive
    /// equation residuals and the actual integral flows.
    ///
    /// Force is N/m³, traction and returned pressure offsets are Pa, coordinates
    /// metres and velocity m/s. The source's fixed density normalizes all native
    /// quantities consistently. Retained sources and separate molecular heating
    /// and modeled transfer belong to the accepted endpoint, without a reporting
    /// source query. The initial velocity is only a guess: no artificial mass
    /// or timestep is introduced into this steady solve.
    ///
    /// # Errors
    /// Preserves source, boundary, port/datum, physical-unit and native nonlinear
    /// refusals. Neither a partial pressure/flow result nor a partial source set
    /// is published. This retains the existing Estimated planar SGS limits.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_navier_stokes_mixed_rheology_sgs_flow_control<G: FluidTransport>(
        &self,
        sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        initial_velocity: &[f64],
        controls: NavierStokesControls,
        absolute_tolerance_m2_s: f64,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialSgsStep<G>>, MaterialSgsError> {
        self.sgs_endpoint(
            sgs,
            force_density,
            boundaries,
            initial_velocity,
            None,
            None,
            controls,
            Some(SgsPortControl {
                ports,
                targets: outward_volume_rates_m2_s,
                absolute_tolerance_m2_s,
            }),
            resolve,
        )
    }

    /// Advance backward Euler with SGS and prescribed endpoint port flows.
    ///
    /// The pressure includes the physical acceleration of the original previous
    /// velocity, as well as molecular/SGS resistance and native convection.
    /// Sources, boundary data and targets describe the new endpoint. Units,
    /// pressure datum and acceptance follow the steady flow-control entry.
    /// History remains immutable on success and refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_rheology_sgs_flow_control<G: FluidTransport>(
        &self,
        sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_tolerance_m2_s: f64,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialSgsStep<G>>, MaterialSgsError> {
        self.bdf2_navier_stokes_mixed_rheology_sgs_variable_step_flow_control(
            sgs,
            force_density,
            boundaries,
            ports,
            outward_volume_rates_m2_s,
            previous_velocity,
            None,
            timestep_s,
            controls,
            absolute_tolerance_m2_s,
            resolve,
        )
    }

    /// Advance variable-interval BDF2 with SGS and unknown port pressures.
    ///
    /// The older velocity is paired with its ACTUAL interval to the previous
    /// velocity [s]. An equal interval gives ordinary BDF2; None gives native
    /// BDF1 startup. The original time owner requires the new/previous interval
    /// ratio in (0, 2]. Both histories enter only its consistent inertial load;
    /// the physical previous velocity initializes the nonlinear solve.
    ///
    /// Required pressure, convection and both constitutive contributions close
    /// at the same endpoint. Pressure responses retain the exact same original
    /// BDF mass/history, without extrapolating viscosity or SGS from either old
    /// velocity. Molecular heating and modeled transfer are instantaneous rates,
    /// not time-integrated heat or a nonlinear passivity guarantee.
    ///
    /// # Errors
    /// Preserves temporal admission before source queries, original material
    /// errors, port-flow/datum checks and pressure-inclusive nonlinear residuals.
    /// Neither caller history nor the reference material is changed.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_sgs_variable_step_flow_control<G: FluidTransport>(
        &self,
        sgs: MaterialSgsConfig<'_>,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_tolerance_m2_s: f64,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialSgsStep<G>>, MaterialSgsError> {
        self.sgs_endpoint(
            sgs,
            force_density,
            boundaries,
            previous_velocity,
            older_history,
            Some(timestep_s),
            controls,
            Some(SgsPortControl {
                ports,
                targets: outward_volume_rates_m2_s,
                absolute_tolerance_m2_s,
            }),
            resolve,
        )
    }
}
