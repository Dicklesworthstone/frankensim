//! Source-owned physical pressure required by resolved integral port flows.

use fs_material::fluid_transport::FluidTransport;

use super::super::{MaterialFluxError, MaterialFluxSystem};
use super::{MaterialStokesBoundary, MaterialStokesSolution, with_physical_boundaries};
use crate::ns::boundary::{StokesFlowControlSolution, StokesSolution};

/// Physical source-backed result of an integral port-flow constraint.
pub struct MaterialStokesFlowControlSolution<S> {
    /// Homogeneous or cell-material resolved flow, including physical reactions.
    pub flow: S,
    /// Uniform pressure added to each controlled port's original traction [Pa].
    pub pressure_offsets_pa: Vec<f64>,
    /// Actual signed outward flow per unit out-of-plane depth [m²/s].
    /// Multiply by the physical depth to obtain volume per time [m³/s].
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Largest absolute error in the requested port flows [m²/s].
    pub maximum_flow_residual_m2_s: f64,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Determine pressure increments required for total port flows in steady Stokes.
    ///
    /// Ports are disjoint original traction-edge groups. At least one
    /// uncontrolled traction edge must retain the pressure datum in each
    /// controlled component. Input body force is N/m³ and total traction Pa;
    /// flow targets and their absolute allowance are m²/s per unit depth.
    /// Velocity profiles remain unknown. The original source supplies viscosity
    /// and the shared reference density for pressure, loads and energy.
    ///
    /// # Errors
    /// Preserves native pressure-datum, flow-constraint and equation checks,
    /// original material errors and physical-output representability.
    pub fn solve_stokes_mixed_flow_control(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        absolute_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesSolution>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        with_physical_boundaries(boundaries, density, |primitive| {
            let prepared = self
                .system
                .assemble_stokes_mixed(
                    self.params,
                    &|point| {
                        let force = force_density(point);
                        [force[0] / density, force[1] / density]
                    },
                    primitive,
                )
                .and_then(|assembly| assembly.prepare_flow_control(ports))
                .map_err(MaterialFluxError::Boundary)?;
            let result = prepared
                .solve(outward_volume_rates_m2_s, absolute_tolerance_m2_s)
                .map_err(MaterialFluxError::Boundary)?;
            physical_control(result, density, |flow| self.material_stokes_solution(flow))
        })
    }
}

pub(super) fn physical_control<S>(
    native: StokesFlowControlSolution,
    density: f64,
    recover: impl FnOnce(StokesSolution) -> Result<S, MaterialFluxError>,
) -> Result<MaterialStokesFlowControlSolution<S>, MaterialFluxError> {
    let pressure_offsets_pa: Vec<_> = native
        .pressure_offsets_m2_s2
        .iter()
        .map(|pressure| density * pressure)
        .collect();
    if pressure_offsets_pa
        .iter()
        .any(|pressure| !pressure.is_finite())
    {
        return Err(MaterialFluxError::Invalid(
            "unrepresentable controlled physical pressure",
        ));
    }
    Ok(MaterialStokesFlowControlSolution {
        flow: recover(native.flow)?,
        pressure_offsets_pa,
        outward_volume_rates_m2_s: native.outward_volume_rates_m2_s,
        maximum_flow_residual_m2_s: native.maximum_flow_residual_m2_s,
    })
}
