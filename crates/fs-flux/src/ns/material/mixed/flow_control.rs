//! Source-owned physical pressure required by resolved integral port flows.

use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

use super::super::{
    MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxSystem, rheology::RheologyField,
};
use super::{
    MaterialStokesBoundary, MaterialStokesRheologySolution, MaterialStokesSolution,
    with_physical_boundaries,
};
use crate::ns::boundary::{StokesAssembly, StokesFlowControlSolution, StokesRheologyControls};

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
        let assembly = with_physical_boundaries(boundaries, density, |primitive| {
            self.system.assemble_stokes_mixed(
                self.params,
                &|point| {
                    let force = force_density(point);
                    [force[0] / density, force[1] / density]
                },
                primitive,
            )
        })
        .map_err(MaterialFluxError::Boundary)?;
        self.solve_controlled_assembly(
            assembly,
            ports,
            outward_volume_rates_m2_s,
            absolute_tolerance_m2_s,
        )
    }

    /// Advance backward-Euler Stokes to prescribed endpoint port flows.
    ///
    /// This solves the pressure needed to accelerate the previous velocity,
    /// retaining its consistent inertia and the material's viscosity. Force
    /// density is N/m³, total traction and required pressure are Pa, time is
    /// seconds, and outward flow targets and absolute tolerance are m²/s per
    /// unit depth. Boundary callbacks and targets describe the new time level.
    /// Port profiles remain resolved unknowns. Each controlled component needs
    /// an uncontrolled traction datum, as in the steady flow-control entry.
    ///
    /// # Errors
    /// Preserves history, boundary, equation and port-flow checks and SI
    /// representability. No history is modified; the returned velocity may be
    /// accepted only after this complete endpoint succeeds.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_stokes_mixed_flow_control(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        previous_velocity: &[f64],
        timestep_s: f64,
        absolute_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesSolution>, MaterialFluxError> {
        self.bdf2_stokes_mixed_variable_step_flow_control(
            force_density,
            boundaries,
            ports,
            outward_volume_rates_m2_s,
            previous_velocity,
            None,
            timestep_s,
            absolute_tolerance_m2_s,
        )
    }

    /// Advance variable-step BDF2 Stokes to prescribed endpoint port flows.
    ///
    /// `older_history` pairs u_(n-1) with its actual interval to u_n in seconds;
    /// `None` performs native backward-Euler startup. The new/previous interval
    /// ratio lies in (0, 2]. Both histories enter the same native momentum
    /// assembly as ordinary material Stokes, so the solved pressure retains
    /// their inertia. Density and viscosity stay frozen on the fixed mesh.
    /// Force, traction, flow units and port-datum requirements are those of
    /// [`Self::bdf1_stokes_mixed_flow_control`]. Required pressures and reactions
    /// describe this discrete endpoint; targets do not prescribe port profiles.
    ///
    /// # Errors
    /// Refuses invalid histories/intervals, incompatible or unresolved port
    /// flows, failed equations and unrepresentable SI outputs. Both caller
    /// histories remain unchanged on success or refusal; no partial endpoint
    /// is published. This adds no convection or material transport.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_stokes_mixed_variable_step_flow_control(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        absolute_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesSolution>, MaterialFluxError> {
        let assembly = self.assemble_material_stokes_bdf2(
            force_density,
            boundaries,
            previous_velocity,
            older_history,
            timestep_s,
        )?;
        self.solve_controlled_assembly(
            assembly,
            ports,
            outward_volume_rates_m2_s,
            absolute_tolerance_m2_s,
        )
    }

    /// Solve the material pressure needed for shear-dependent steady port flows.
    ///
    /// Each cell's source is queried at its actual resolved shear, and the
    /// pressure needed for the requested integral flow changes with that
    /// viscosity. The existing coupled Picard loop accepts only the final
    /// viscosity's nonlinear equations and the original port-flow tolerances.
    /// Pressure offsets are relaxed with the complete field. Returned source
    /// owners, pressure, reactions, kinetic energy and viscous heat all belong
    /// to that accepted endpoint.
    ///
    /// Force density is N/m³, total traction and required pressures are Pa,
    /// velocity is m/s, and signed outward flow and its absolute tolerance are
    /// m²/s per unit depth. Ports retain the original traction-datum rules;
    /// their velocity profiles remain unknown. Sources must retain the bound
    /// reference density; captured temperature/composition remain fixed. This
    /// is generalized-Newtonian creeping flow, without convection, thermal
    /// evolution, material transport or viscoelastic memory.
    ///
    /// # Errors
    /// Preserves typed material refusals, native boundary/flow admission,
    /// coupled nonlinear convergence and SI representability. The initial
    /// velocity is immutable; refusal publishes no partial flow or source set.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_stokes_mixed_rheology_flow_control<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        absolute_tolerance_m2_s: f64,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<
        MaterialStokesFlowControlSolution<MaterialStokesRheologySolution<G>>,
        MaterialFluxError,
    > {
        let density = self.transport.density_kg_m3();
        let field = RheologyField::new(&self.system.mesh.centroids, &self.transport, resolve);
        let result = with_physical_boundaries(boundaries, density, |primitive| {
            self.system
                .solve_stokes_mixed_rheology_flow_control_checked(
                    self.params,
                    &|cell, _, rate| field.viscosity(cell, rate),
                    &|point| {
                        let physical = force_density(point);
                        [physical[0] / density, physical[1] / density]
                    },
                    primitive,
                    ports,
                    outward_volume_rates_m2_s,
                    initial_velocity,
                    controls,
                    absolute_tolerance_m2_s,
                    &|| field.is_valid(),
                )
        })
        .map_err(MaterialFluxError::Boundary);
        let (native, samples) = field.finish(result)?;
        physical_control(native, density, |flow| {
            self.material_stokes_rheology_solution(flow, samples)
        })
    }

    fn solve_controlled_assembly(
        &self,
        assembly: StokesAssembly<'_, '_>,
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        absolute_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesSolution>, MaterialFluxError> {
        let prepared = assembly
            .prepare_flow_control(ports)
            .map_err(MaterialFluxError::Boundary)?;
        let result = prepared
            .solve(outward_volume_rates_m2_s, absolute_tolerance_m2_s)
            .map_err(MaterialFluxError::Boundary)?;
        physical_control(result, self.transport.density_kg_m3(), |flow| {
            self.material_stokes_solution(flow)
        })
    }
}

pub(in crate::ns::material) fn physical_control<N, S>(
    native: StokesFlowControlSolution<N>,
    density: f64,
    recover: impl FnOnce(N) -> Result<S, MaterialFluxError>,
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
