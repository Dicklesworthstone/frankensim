//! Molecular endpoint rheology with the original integral-flow and RCR owners.

use super::super::windkessel::DensityMomentumWindkesselError;
use super::super::{
    CoupledViscosity, DensityMomentumControls, DensityTransportCell, FluxSystem,
    MaterialFluxNumerics, MaterialNavierStokesBoundary, Result, StokesRheologyControls,
    physical_flow_control,
};
use super::{CellRheology, CoupledDensityMomentumRheologyStep};
use crate::material::{
    MaterialStokesWindkesselHistory, MaterialStokesWindkesselPort, MaterialStokesWindkesselSolution,
};
use crate::ns::boundary::flow_control::FlowControl;

/// One accepted density/rheology endpoint with physical integral-flow pressures.
pub struct CoupledDensityMomentumRheologyFlowControlStep<const N: usize> {
    /// Retained transport, source samples, flow, reactions and molecular heating.
    pub step: CoupledDensityMomentumRheologyStep<N>,
    /// Uniform pressures added to the ports' original traction as `-p*n` [Pa].
    pub pressure_offsets_pa: Vec<f64>,
    /// Accepted signed outward volume flow per unit depth [m²/s].
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Largest accepted absolute error in the requested volume rates [m²/s].
    pub maximum_flow_residual_m2_s: f64,
}

impl FluxSystem<'_> {
    /// Resolve integral port pressure, transported density and molecular rheology
    /// at one native BDF1 endpoint.
    ///
    /// Every trial rebuilds native pressure responses using that same trial's
    /// transported density, objective shear and retained cellwise source stress.
    /// Accepted port increments also enter pressure, reactions and boundary work.
    /// The original uncontrolled traction datum, disjoint-port admission, volume
    /// units and donor/reversal rules of
    /// [`Self::density_momentum_coupled_bdf1_flow_control_step`] apply unchanged.
    /// Source purity, original errors, convergence, molecular-heating separation
    /// and atomic refusal follow [`Self::density_momentum_coupled_bdf1_rheology_step`].
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_coupled_bdf1_rheology_flow_control_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &CellRheology<'_, N>,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        absolute_tolerance_m2_s: f64,
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        nonlinear: StokesRheologyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledDensityMomentumRheologyFlowControlStep<N>> {
        let (solved, response) = self.coupled_density_momentum_bdf1_viscosity(
            numerics,
            CoupledViscosity::Rheology(dynamic_viscosity),
            force_density,
            boundaries,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            nonlinear,
            Some(FlowControl::targets(
                ports,
                outward_volume_rates_m2_s,
                absolute_tolerance_m2_s,
            )),
            &mut cancelled,
        )?;
        let physical = physical_flow_control(solved, controls, &mut cancelled)?;
        Ok(CoupledDensityMomentumRheologyFlowControlStep {
            step: self.finish_coupled_rheology(
                physical.step,
                response,
                controls,
                &mut cancelled,
            )?,
            pressure_offsets_pa: physical.pressure_offsets_pa,
            outward_volume_rates_m2_s: physical.outward_volume_rates_m2_s,
            maximum_flow_residual_m2_s: physical.maximum_flow_residual_m2_s,
        })
    }

    /// Resolve source molecular stress, transported density and finite RCR load
    /// histories at one accepted native BDF1 endpoint.
    ///
    /// The original RCR reduction is formed once from immutable old velocity and
    /// compliance. Its pressure responses are rebuilt at every density/shear
    /// trial; accepted original pressure and storage laws are checked with the
    /// measured endpoint flow. Reporting retains the same source samples without
    /// querying again. Source refusals stay intact in the momentum error variant.
    ///
    /// All physical units, history, donor, cancellation and reservoir-reference
    /// semantics of [`Self::density_momentum_coupled_bdf1_windkessel_step`] apply.
    /// Molecular heat in `flow` is depth-integrated W, while external RCR powers
    /// remain W/m. External resistor loss, BDF defect and molecular heating are
    /// separate outputs; none is deposited in a transported caloric inventory.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_coupled_bdf1_rheology_windkessel_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &CellRheology<'_, N>,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesWindkesselPort],
        previous: MaterialStokesWindkesselHistory<'_>,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        nonlinear: StokesRheologyControls,
        absolute_pressure_tolerance_pa: f64,
        absolute_flow_tolerance_m2_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        MaterialStokesWindkesselSolution<CoupledDensityMomentumRheologyStep<N>>,
        DensityMomentumWindkesselError,
    > {
        self.coupled_density_momentum_bdf1_windkessel_viscosity(
            numerics,
            CoupledViscosity::Rheology(dynamic_viscosity),
            force_density,
            boundaries,
            ports,
            loads,
            previous,
            cells,
            inflow,
            timestep_s,
            controls,
            nonlinear,
            absolute_pressure_tolerance_pa,
            absolute_flow_tolerance_m2_s,
            |step, response, cancelled| {
                self.finish_coupled_rheology(step, response, controls, cancelled)
            },
            cancelled,
        )
    }
}
