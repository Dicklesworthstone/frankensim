//! Frozen source-material fields for density-weighted mixed Stokes inertia.

use std::sync::Arc;

use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

use super::super::{
    MaterialFluxError, MaterialFluxSystem, MaterialViscousDissipation, physical_viscous_dissipation,
};
use super::windkessel::{
    MaterialStokesWindkesselHistory, MaterialStokesWindkesselPort,
    MaterialStokesWindkesselSolution, solve_material_windkessel_assembly,
};
use super::{
    MaterialBoundaryLoad, MaterialStokesBoundary, MaterialStokesFlowControlSolution,
    MaterialStokesResistancePort, MaterialStokesResistanceSolution, flow_control::physical_control,
    physical_boundary_load, resistance::solve_resistance_assembly, with_physical_boundaries,
};
use crate::ns::boundary::{StokesAssembly, StokesPressureReference, StokesSolution};
use crate::ns::{FluxParams, FluxSolution, FluxSystem};

/// One original material owner and its admitted, immutable cell snapshot.
/// Cell order is the original triangle order. A raw `FluidTransportState`
/// input remains an explicitly unsourced numerical fixture.
pub struct MaterialStokesCell<G> {
    material: G,
    transport: FluidTransportState,
    relative_density: f64,
    normalized_viscosity: f64,
}

impl<G> MaterialStokesCell<G> {
    /// Original source, including its query conditions and source receipts.
    #[must_use]
    pub const fn material(&self) -> &G {
        &self.material
    }

    /// The single source evaluation used for inertia, viscosity and reporting.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
}

/// Prepared, fixed-mesh material field for linear incompressible Stokes.
///
/// Each original triangle owns one density and dynamic viscosity. All source
/// evaluations and reference-density conversions occur once during binding;
/// repeated timesteps and physical reporting never re-query a source. Density
/// is frozen at both time levels. This represents stationary stratification
/// or shear parallel to material layers, not material transport or mixing.
pub struct MaterialStokesField<'s, 'm, G> {
    system: &'s FluxSystem<'m>,
    params: FluxParams,
    reference_density_kg_m3: f64,
    samples: Arc<[MaterialStokesCell<G>]>,
}

/// Physical result retaining the exact spatial material owners used to solve.
/// Source snapshots are shared with the prepared field without copying owners
/// or querying their laws again.
pub struct MaterialStokesFieldSolution<G> {
    flow: StokesSolution,
    reference_density_kg_m3: f64,
    samples: Arc<[MaterialStokesCell<G>]>,
    pressure_pa: Vec<f64>,
    kinetic_energy_j_per_m: f64,
    viscous_dissipation: MaterialViscousDissipation,
}

impl<'s, 'm, F: FluidTransport + ?Sized> MaterialFluxSystem<'s, 'm, '_, F> {
    /// Bind one original transport owner per triangle to weighted Stokes.
    ///
    /// The existing binding supplies numerical controls and ONE fixed
    /// reference density, used only as a common normalization. Each cell's
    /// original source supplies its actual `rho` and `mu`; no homogeneous
    /// reference viscosity enters this field's equations. Heterogeneous
    /// sources may be passed as `Vec<&dyn FluidTransport>`.
    ///
    /// The resulting equation is `rho(x) du/dt - div(2 mu(x) D(u)) + grad(p)
    /// = force_density`, with div(u)=0. Internally `rho/rho_reference`,
    /// `mu/rho_reference`, pressure, force and traction all share that same
    /// reference; local `mu/rho(x)` or `p/rho(x)` would be a different PDE.
    ///
    /// # Errors
    /// Refuses the wrong number of sources before any source queries, keeps
    /// original material errors intact, and refuses unrepresentable positive
    /// normalization. Failed binding publishes no partial material field.
    pub fn bind_stokes_materials<G: FluidTransport>(
        &self,
        materials: Vec<G>,
    ) -> Result<MaterialStokesField<'s, 'm, G>, MaterialFluxError> {
        if materials.len() != self.system.mesh.tris.len() {
            return Err(MaterialFluxError::Invalid(
                "one frozen Stokes material per original triangle is required",
            ));
        }
        let reference_density_kg_m3 = self.transport.density_kg_m3();
        let mut samples = Vec::with_capacity(materials.len());
        for material in materials {
            let transport = material.transport_state()?;
            let relative_density = transport.density_kg_m3() / reference_density_kg_m3;
            let normalized_viscosity = transport.dynamic_viscosity_pa_s() / reference_density_kg_m3;
            if ![relative_density, normalized_viscosity]
                .into_iter()
                .all(|value| value.is_finite() && value > 0.0)
            {
                return Err(MaterialFluxError::Invalid(
                    "spatial Stokes density or viscosity normalization is not representable",
                ));
            }
            samples.push(MaterialStokesCell {
                material,
                transport,
                relative_density,
                normalized_viscosity,
            });
        }
        Ok(MaterialStokesField {
            system: self.system,
            params: self.params,
            reference_density_kg_m3,
            samples: samples.into(),
        })
    }
}

impl<'s, 'm, G> MaterialStokesField<'s, 'm, G> {
    /// Exact original owners and frozen snapshots, in triangle order.
    #[must_use]
    pub fn samples(&self) -> &[MaterialStokesCell<G>] {
        &self.samples
    }

    /// Common numerical normalization [kg/m³], not a local material density.
    #[must_use]
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density_kg_m3
    }

    /// Advance frozen source-material layers and external RCR/inertance loads.
    ///
    /// Every cell retains its original rho/rho_reference mass and
    /// mu/rho_reference stress; the SAME reference normalizes all port pressures,
    /// R,L and compliance. Fluid and storage share BDF1 startup or the supplied
    /// paired BDF2 histories and interval. Inertance Q histories come directly
    /// from those native velocities. No source is re-queried during the solve.
    ///
    /// External storage, resistor heat and signed numerical BDF defect remain
    /// separate from bulk fluid heating. External coefficients and compliance
    /// reference stay fixed through the histories; reservoirs are held per
    /// interval. This fixed-layer model does not advect density or material
    /// interfaces, update temperature, or infer compliance from transport data.
    /// Every traction port may be controlled; reservoir laws select pressure.
    ///
    /// # Errors
    /// Preserves source snapshots and immutable histories on every refusal;
    /// checks original continuum and both dynamic port laws before publication.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_stokes_mixed_variable_step_windkessel_ports(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesWindkesselPort],
        previous: MaterialStokesWindkesselHistory<'_>,
        older_history: Option<(MaterialStokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        absolute_pressure_tolerance_pa: f64,
        absolute_flow_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesWindkesselSolution<MaterialStokesFieldSolution<G>>, MaterialFluxError>
    {
        let assembly = self.assemble_material_stokes_bdf2(
            force_density,
            boundaries,
            previous.velocity,
            older_history.map(|(history, dt)| (history.velocity, dt)),
            timestep_s,
        )?;
        solve_material_windkessel_assembly(
            assembly,
            self.reference_density_kg_m3,
            ports,
            loads,
            previous,
            older_history,
            timestep_s,
            absolute_pressure_tolerance_pa,
            absolute_flow_tolerance_m2_s,
            |flow, energy| self.recover_solution(flow, energy),
        )
    }

    /// Resolve steady laminar port resistances through frozen material layers.
    ///
    /// Reservoir increments and tolerance are Pa, force density is N/m³ and
    /// resistance is Pa s/m² for flow per depth. The same reference density
    /// normalizes every port, while each cell retains its OWN dynamic viscosity.
    /// Source snapshots, bulk heat and reactions are retained without re-query.
    /// External R*Q² heat remains separate. Each controlled component retains
    /// reservoir-selected pressure; this does not advect material layers.
    ///
    /// # Errors
    /// Preserves native port/equation refusals and physical output checks.
    pub fn solve_stokes_mixed_resistance_ports(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesResistancePort],
        absolute_pressure_tolerance_pa: f64,
    ) -> Result<MaterialStokesResistanceSolution<MaterialStokesFieldSolution<G>>, MaterialFluxError>
    {
        let reference = self.reference_density_kg_m3;
        let assembly = with_physical_boundaries(boundaries, reference, |primitive| {
            self.system.assemble_stokes_mixed_variable_density(
                self.params,
                &|cell, _| self.samples[cell].relative_density,
                &|cell, _| self.samples[cell].normalized_viscosity,
                &|point| force_density(point).map(|f| f / reference),
                primitive,
            )
        })
        .map_err(MaterialFluxError::Boundary)?;
        solve_resistance_assembly(
            assembly,
            reference,
            ports,
            loads,
            absolute_pressure_tolerance_pa,
            |flow, energy| self.recover_solution(flow, energy),
        )
    }

    /// Couple endpoint reservoir resistances to density-weighted BDF2 Stokes.
    ///
    /// `older_history=None` performs BDF1 startup. Source snapshots and the
    /// density metric stay frozen; both histories, endpoint pressure and all
    /// port interactions use the original momentum operator. Physical units
    /// and external-model scope match the steady resistance entry.
    ///
    /// # Errors
    /// Preserves native history/ratio, port/equation and physical output checks.
    /// Both caller histories and all material owners remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_stokes_mixed_variable_step_resistance_ports(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesResistancePort],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        absolute_pressure_tolerance_pa: f64,
    ) -> Result<MaterialStokesResistanceSolution<MaterialStokesFieldSolution<G>>, MaterialFluxError>
    {
        let assembly = self.assemble_material_stokes_bdf2(
            force_density,
            boundaries,
            previous_velocity,
            older_history,
            timestep_s,
        )?;
        solve_resistance_assembly(
            assembly,
            self.reference_density_kg_m3,
            ports,
            loads,
            absolute_pressure_tolerance_pa,
            |flow, energy| self.recover_solution(flow, energy),
        )
    }

    /// Resolve the pressures needed for total port flows through frozen materials.
    ///
    /// Input force density is N/m³, total traction is Pa, and signed outward
    /// target flows are m²/s per unit depth. Each controlled port contains only
    /// original traction edges; every controlled component retains a separate
    /// uncontrolled traction datum. No normal velocity profile is invented.
    /// Original one-sided viscosities, density-weighted energy, source owners
    /// and pressure/load normalization are retained without source re-query.
    ///
    /// # Errors
    /// Preserves native flow/datum/residual admission and physical SI recovery.
    pub fn solve_stokes_mixed_flow_control(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        absolute_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesFieldSolution<G>>, MaterialFluxError>
    {
        let reference = self.reference_density_kg_m3;
        let assembly = with_physical_boundaries(boundaries, reference, |primitive| {
            self.system.assemble_stokes_mixed_variable_density(
                self.params,
                &|cell, _| self.samples[cell].relative_density,
                &|cell, _| self.samples[cell].normalized_viscosity,
                &|point| {
                    let force = force_density(point);
                    [force[0] / reference, force[1] / reference]
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

    /// Resolve the pressure for one backward-Euler endpoint flow constraint.
    ///
    /// Each frozen cell's original density weights inertia and history, and
    /// its dynamic viscosity supplies the one-sided stress. Input force is
    /// N/m³, total traction and required pressure are Pa, velocity is m/s and
    /// time is seconds. Outward flow targets and their absolute tolerance are
    /// m²/s per depth. Callbacks and targets specify the new time level.
    /// Port profiles remain unknown; every controlled component retains an
    /// uncontrolled traction datum. No material source is queried again.
    ///
    /// # Errors
    /// Retains native history, material, port-flow and equation checks and SI
    /// recovery. Previous velocity and frozen sources remain unchanged; only
    /// a completely accepted endpoint supplies the next velocity history.
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
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesFieldSolution<G>>, MaterialFluxError>
    {
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

    /// Resolve variable-step BDF2 pressure for prescribed endpoint port flows.
    ///
    /// `older_history` pairs u_(n-1) with the actual preceding interval; `None`
    /// performs native backward-Euler startup. The new/previous step ratio is
    /// in (0, 2]. The same frozen cell density metric weights both histories,
    /// endpoint inertia, reactions and physical kinetic energy. Dynamic
    /// viscosities and source owners remain those bound to the original cells.
    /// Force, traction, flow units and pressure datums follow
    /// [`Self::bdf1_stokes_mixed_flow_control`]. No convection, density transport
    /// or material-interface evolution is added.
    ///
    /// # Errors
    /// Retains native history/ratio, flow/datum and equation checks and SI
    /// representability. Histories and source snapshots are immutable; a
    /// refused candidate returns no partial velocity or pressure endpoint.
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
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesFieldSolution<G>>, MaterialFluxError>
    {
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

    fn solve_controlled_assembly(
        &self,
        assembly: StokesAssembly<'_, '_>,
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        absolute_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesFlowControlSolution<MaterialStokesFieldSolution<G>>, MaterialFluxError>
    {
        let prepared = assembly
            .prepare_flow_control(ports)
            .map_err(MaterialFluxError::Boundary)?;
        let result = prepared
            .solve(outward_volume_rates_m2_s, absolute_tolerance_m2_s)
            .map_err(MaterialFluxError::Boundary)?;
        let energy = self.reference_density_kg_m3
            * prepared
                .kinetic_energy(&result.flow.field.x)
                .map_err(MaterialFluxError::Boundary)?;
        physical_control(result, self.reference_density_kg_m3, |flow| {
            self.recover_solution(flow, energy)
        })
    }

    /// Advance weighted Stokes by one backward-Euler step in coherent SI.
    ///
    /// Force is per volume [N/m³], traction is total stress [Pa], velocity is
    /// in m/s and time in seconds. Callbacks specify the endpoint and must be
    /// pure during admission/solution. Previous velocity is in the same fixed
    /// mesh and frozen material field. Changing density requires a separate
    /// conservative state transfer; this entry does not transport density,
    /// advance interfaces, add convection or infer thermodynamic pressure.
    ///
    /// Density-weighted consistent inertia, history, support reactions and
    /// kinetic energy use the same native quadrature. Dynamic viscosity is
    /// one-sided at material interfaces; bulk heat uses that actual viscosity.
    ///
    /// # Errors
    /// Preserves native admission and convergence errors and refuses
    /// unrepresentable physical outputs. History and sources remain unchanged.
    pub fn bdf1_stokes_mixed_step(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<MaterialStokesFieldSolution<G>, MaterialFluxError> {
        self.bdf2_stokes_mixed_step(
            force_density,
            boundaries,
            previous_velocity,
            None,
            timestep_s,
        )
    }

    /// Advance the frozen material field with equal-step second-order BDF2.
    ///
    /// `previous_velocity` is u_n. Some(u_(n-1)) supplies the second immutable
    /// history; BOTH time intervals must equal timestep_s. None performs the
    /// exact original BDF1 startup. Source owners are not queried again: one
    /// density-weighted metric acts on all three velocity levels and supplies
    /// reactions and physical kinetic energy. Dynamic viscosity, pressure,
    /// force and traction retain the same fixed reference normalization.
    ///
    /// Equal-step G-energy stability has the scope of the native
    /// [`FluxSystem::assemble_stokes_mixed_bdf2`] entry. Single-state kinetic
    /// energy need not decrease at every unforced step. Use
    /// [`Self::bdf2_stokes_mixed_variable_step`] when intervals change. Neither
    /// density transport nor nonlinear dynamics is included.
    ///
    /// # Errors
    /// Preserves native two-history admission, convergence and physical-output
    /// refusals. Neither history nor the frozen source snapshots are changed.
    pub fn bdf2_stokes_mixed_step(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_velocity: Option<&[f64]>,
        timestep_s: f64,
    ) -> Result<MaterialStokesFieldSolution<G>, MaterialFluxError> {
        self.bdf2_stokes_mixed_variable_step(
            force_density,
            boundaries,
            previous_velocity,
            older_velocity.map(|older| (older, timestep_s)),
            timestep_s,
        )
    }

    /// Advance frozen density/viscosity layers with unequal BDF2 timesteps.
    ///
    /// `older_history` pairs u_(n-1) with its actual interval to u_n in seconds.
    /// `None` starts with BDF1; otherwise the new/previous step ratio lies in
    /// (0, 2]. One unchanged spatial density metric weights both histories and
    /// recovers physical energy and reactions. Source snapshots are not queried
    /// again. All endpoint loads and recovered outputs retain physical SI units.
    /// See [`FluxSystem::assemble_stokes_mixed_bdf2_variable_step`] for stability
    /// scope. Neither density transport nor nonlinear advection is introduced.
    ///
    /// # Errors
    /// Preserves native history/ratio admission, convergence and SI recovery;
    /// caller-owned histories and the frozen material field stay unchanged.
    pub fn bdf2_stokes_mixed_variable_step(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
    ) -> Result<MaterialStokesFieldSolution<G>, MaterialFluxError> {
        let assembly = self.assemble_material_stokes_bdf2(
            force_density,
            boundaries,
            previous_velocity,
            older_history,
            timestep_s,
        )?;
        let flow = assembly.solve().map_err(MaterialFluxError::Boundary)?;
        let kinetic_energy_j_per_m = self.reference_density_kg_m3
            * assembly
                .kinetic_energy(&flow.field.x)
                .map_err(MaterialFluxError::Boundary)?;
        self.recover_solution(flow, kinetic_energy_j_per_m)
    }

    // Share the exact sampled material metric and native BDF history between
    // prescribed-traction and prescribed-flow endpoints.
    fn assemble_material_stokes_bdf2(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
    ) -> Result<StokesAssembly<'s, 'm>, MaterialFluxError> {
        let reference = self.reference_density_kg_m3;
        with_physical_boundaries(boundaries, reference, |primitive_boundaries| {
            self.system
                .assemble_stokes_mixed_bdf2_variable_density_variable_step(
                    self.params,
                    &|cell, _| self.samples[cell].relative_density,
                    &|cell, _| self.samples[cell].normalized_viscosity,
                    &|point| {
                        let physical = force_density(point);
                        [physical[0] / reference, physical[1] / reference]
                    },
                    primitive_boundaries,
                    previous_velocity,
                    older_history,
                    timestep_s,
                )
                .map_err(MaterialFluxError::Boundary)
        })
    }

    fn recover_solution(
        &self,
        flow: StokesSolution,
        kinetic_energy_j_per_m: f64,
    ) -> Result<MaterialStokesFieldSolution<G>, MaterialFluxError> {
        let reference = self.reference_density_kg_m3;
        let pressure_pa: Vec<_> = flow.pressure.iter().map(|p| reference * p).collect();
        let viscous_dissipation =
            physical_viscous_dissipation(self.system, &flow.field.x, &|cell| {
                self.samples[cell].transport.dynamic_viscosity_pa_s()
            })?;
        if pressure_pa.iter().any(|value| !value.is_finite()) || !kinetic_energy_j_per_m.is_finite()
        {
            return Err(MaterialFluxError::Invalid(
                "physical pressure or energy is not representable",
            ));
        }
        Ok(MaterialStokesFieldSolution {
            flow,
            reference_density_kg_m3: reference,
            samples: Arc::clone(&self.samples),
            pressure_pa,
            kinetic_energy_j_per_m,
            viscous_dissipation,
        })
    }
}

impl<G> MaterialStokesFieldSolution<G> {
    /// Original discrete result, including checked linear residuals.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        &self.flow.field
    }

    /// Physical velocity moments [m/s] for the next step of this material field.
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        &self.flow.field.x[..self.flow.field.n_u]
    }

    /// Physical pressure [Pa] in the original component pressure datums.
    /// This is incompressible mechanical pressure, not an EOS state update.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }

    /// Declared pressure interpretation, preserving natural traction datums.
    #[must_use]
    pub const fn pressure_reference(&self) -> StokesPressureReference {
        self.flow.pressure_reference
    }

    /// Local datum in original cell order; None outside this mesh.
    #[must_use]
    pub fn pressure_reference_at(&self, cell: usize) -> Option<StokesPressureReference> {
        self.flow.pressure_reference_per_cell.get(cell).copied()
    }

    /// Physical kinetic energy per out-of-plane depth [J/m].
    #[must_use]
    pub const fn kinetic_energy_j_per_m(&self) -> f64 {
        self.kinetic_energy_j_per_m
    }

    /// Physical bulk heat from the same cell viscosities and accepted velocity.
    #[must_use]
    pub const fn viscous_dissipation(&self) -> &MaterialViscousDissipation {
        &self.viscous_dissipation
    }

    /// Exact original source owners and snapshots, retained after the solve.
    #[must_use]
    pub fn samples(&self) -> &[MaterialStokesCell<G>] {
        &self.samples
    }

    /// Stress forces [N/m], moments [N] and slip power [W/m] about origin_m.
    /// Offsets for selected mean-zero components are explicit physical Pa.
    /// Advective transport is zero in this linear Stokes entry.
    ///
    /// # Errors
    /// Preserves original edge/datum admission and refuses nonfinite outputs.
    pub fn boundary_loads(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryLoad, MaterialFluxError> {
        physical_boundary_load(
            &self.flow,
            self.reference_density_kg_m3,
            edges,
            origin_m,
            pressure_offsets_pa,
        )
    }
}
