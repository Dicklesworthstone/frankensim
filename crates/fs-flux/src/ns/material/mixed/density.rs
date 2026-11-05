//! Frozen source-material fields for density-weighted mixed Stokes inertia.

use std::sync::Arc;

use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

use super::super::{
    MaterialFluxError, MaterialFluxSystem, MaterialViscousDissipation, physical_viscous_dissipation,
};
use super::{
    MaterialBoundaryLoad, MaterialStokesBoundary, physical_boundary_load, with_physical_boundaries,
};
use crate::ns::boundary::{StokesPressureReference, StokesSolution};
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

impl<G> MaterialStokesField<'_, '_, G> {
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
        let reference = self.reference_density_kg_m3;
        with_physical_boundaries(boundaries, reference, |primitive_boundaries| {
            let assembly = self
                .system
                .assemble_stokes_mixed_bdf1_variable_density(
                    self.params,
                    &|cell, _| self.samples[cell].relative_density,
                    &|cell, _| self.samples[cell].normalized_viscosity,
                    &|point| {
                        let physical = force_density(point);
                        [physical[0] / reference, physical[1] / reference]
                    },
                    primitive_boundaries,
                    previous_velocity,
                    timestep_s,
                )
                .map_err(MaterialFluxError::Boundary)?;
            let flow = assembly.solve().map_err(MaterialFluxError::Boundary)?;
            let pressure_pa: Vec<_> = flow.pressure.iter().map(|p| reference * p).collect();
            let kinetic_energy_j_per_m = reference
                * assembly
                    .kinetic_energy(&flow.field.x)
                    .map_err(MaterialFluxError::Boundary)?;
            let viscous_dissipation =
                physical_viscous_dissipation(self.system, &flow.field.x, &|cell| {
                    self.samples[cell].transport.dynamic_viscosity_pa_s()
                })?;
            if pressure_pa.iter().any(|value| !value.is_finite())
                || !kinetic_energy_j_per_m.is_finite()
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
