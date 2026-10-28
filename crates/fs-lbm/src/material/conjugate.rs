//! One material snapshot drives the existing conjugate flow/heat pipeline.
//!
//! Gas composition, pressure and temperature, or a liquid material selection,
//! are evaluated by the original `FluidTransport` owner. This adapter neither
//! chooses a fluid by name nor invents another equation of state. Momentum,
//! heat-capacity flux and molecular/eddy conductivity use that SAME snapshot.
//! The forced-flow entry point freezes properties. [`buoyant`] explicitly closes
//! source-resolved steady temperature feedback through these same numerical
//! owners; neither path infers compressibility, chemistry or phase change.

/// Closed, source-resolved steady buoyancy, viscosity and conductivity feedback.
pub mod buoyant;
/// Original solid material sources on the coupled fluid/solid temperature field.
pub mod solid;

use fs_exec::CancelGate;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use crate::conjugate::{ChtError, EnergyConfig, EnergySolution, FluidProperties,
    FvFlow, SimpleConfig, SolidMaterial, ThermalSetup, VoxelDomain, simple_flow, solve_energy};
use super::{MaterialLbmError, MaterialScaling, positive};

/// A fluid's original source owner, immutable transport snapshot and CHT view.
pub struct MaterialChtFluid<'a, F: FluidTransport + ?Sized> {
    material: &'a F,
    transport: FluidTransportState,
    properties: FluidProperties,
}

impl<'a, F: FluidTransport + ?Sized> MaterialChtFluid<'a, F> {
    /// Evaluate the original source exactly once, with no dry-air fallback.
    /// Source receipts remain on `material()`, not on a fabricated CHT card.
    ///
    /// # Errors
    /// Propagates source refusals and rejects derived products required by CHT
    /// which cannot be represented even if the factored transport ratios can.
    pub fn try_new(material: &'a F) -> Result<Self, MaterialLbmError> {
        Self::from_snapshot(material, material.transport_state()?)
    }

    fn from_snapshot(material: &'a F, transport: FluidTransportState) -> Result<Self, MaterialLbmError> {
        let properties = FluidProperties {
            density_kg_m3: transport.density_kg_m3(),
            specific_heat_j_kg_k: transport.specific_heat_j_kg_k(),
            conductivity_w_m_k: transport.thermal_conductivity_w_m_k(),
            kinematic_viscosity_m2_s: transport.kinematic_viscosity_m2_s(),
        };
        positive(properties.volumetric_heat_capacity(), "CHT volumetric heat capacity is unrepresentable")?;
        positive(properties.diffusivity(), "CHT thermal diffusivity is unrepresentable")?;
        positive(properties.prandtl(), "CHT Prandtl number is unrepresentable")?;
        positive(properties.density_kg_m3 * properties.kinematic_viscosity_m2_s,
            "CHT dynamic viscosity is unrepresentable")?;
        Ok(Self { material, transport, properties })
    }

    /// Original material evaluation, including its source receipts when present.
    #[must_use]
    pub const fn material(&self) -> &'a F { self.material }
    /// Exact frozen transport state; no solver samples the owner again.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// The ordinary CHT parameter type for lower-level existing solver calls.
    #[must_use]
    pub const fn properties(&self) -> &FluidProperties { &self.properties }

    /// Solve momentum with the retained material snapshot and keep its actual
    /// conservative face fluxes together with their domain and source owner.
    /// The result can drive the shared enthalpy timeline through
    /// `MaterialChtFlow::advect_enthalpy` when `enthalpy-coupling` is enabled.
    /// No constant-Cp energy solve is required to obtain the flow.
    pub fn solve_flow(
        &self,
        domain: &'a VoxelDomain,
        config: &SimpleConfig,
        gate: &CancelGate,
    ) -> Result<MaterialChtFlow<'a, F>, ChtError> {
        let flow = simple_flow(domain, &self.properties, config, gate)?;
        Ok(MaterialChtFlow {
            material: self.material, transport: self.transport, domain, flow,
        })
    }

    /// Run existing SIMPLEC flow, then the conservative conjugate energy solve.
    ///
    /// The solved volume flux is converted to heat-capacity flux using this
    /// same fluid. LVEL's wall-corrected eddy conductivity is forwarded exactly
    /// once through the existing `FvFlow::eddy_conductivity` owner. Molecular
    /// conductivity, solid tables, interfaces, sources and solver controls are
    /// unchanged. This is a one-way, frozen-fluid forced-convection solve, NOT
    /// feedback of the computed temperature into fluid density or viscosity.
    ///
    /// # Errors
    /// Rejects an independently authored eddy-conductivity field instead of
    /// silently replacing or double-counting it. Propagates normal solver,
    /// boundary, material-index and cancellation refusals; no partial result is
    /// returned and no caller-owned input is changed.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_forced_convection(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        gate: &CancelGate,
    ) -> Result<MaterialChtSolution<'a, F>, ChtError> {
        if !thermal.eddy_conductivity_w_m_k.is_empty() {
            return Err(ChtError::InvalidInput {
                field: "thermal.eddy_conductivity_w_m_k",
                reason: "material-coupled flow owns eddy conductivity; an independent field would overwrite or double-count it".into(),
            });
        }
        let flow = simple_flow(domain, &self.properties, flow_config, gate)?;
        let mut setup = thermal.clone();
        setup.eddy_conductivity_w_m_k = flow.eddy_conductivity(&self.properties);
        // Fail before conductivity addition if a representable eddy term and
        // molecular term together overflow the energy kernel's coefficient.
        for (cell, &eddy) in setup.eddy_conductivity_w_m_k.iter().enumerate() {
            if !eddy.is_finite() || eddy < 0.0
                || (domain.is_fluid(cell) && !(self.properties.conductivity_w_m_k + eddy).is_finite())
            {
                return Err(ChtError::InvalidInput {
                    field: "thermal.eddy_conductivity_w_m_k",
                    reason: "flow-derived thermal conductivity is unrepresentable".into(),
                });
            }
        }
        let energy = solve_energy(domain, &self.properties, solids, &flow.field, &setup, energy_config, gate)?;
        Ok(MaterialChtSolution { material: self.material, transport: self.transport, flow, energy })
    }
}

/// A material-driven momentum solve, with the original domain and volume flux.
/// Molecular viscosity and density are the same snapshot later used to turn
/// volume flux into advected mass. Caloric storage remains a separate law.
pub struct MaterialChtFlow<'a, F: FluidTransport + ?Sized> {
    material: &'a F,
    transport: FluidTransportState,
    domain: &'a VoxelDomain,
    flow: FvFlow,
}

impl<'a, F: FluidTransport + ?Sized> MaterialChtFlow<'a, F> {
    /// Original material evaluation, including its source receipts.
    #[must_use]
    pub const fn material(&self) -> &'a F { self.material }
    /// Exact snapshot used by the momentum solve.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// Original geometry, immutably borrowed for this solution's lifetime.
    #[must_use]
    pub const fn domain(&self) -> &'a VoxelDomain { self.domain }
    /// Existing SIMPLEC result, including its conservative volume face fluxes.
    #[must_use]
    pub const fn flow(&self) -> &FvFlow { &self.flow }
}

impl<'a, F: FluidTransport + ?Sized> MaterialScaling<'a, F> {
    /// Use the lattice plan's EXACT material snapshot for the CHT solver too.
    /// This never re-evaluates an interior-mutable material owner.
    ///
    /// # Errors
    /// Refuses unrepresentable CHT products, as `MaterialChtFluid::try_new`.
    pub fn conjugate_fluid(&self) -> Result<MaterialChtFluid<'a, F>, MaterialLbmError> {
        MaterialChtFluid::from_snapshot(self.material, self.transport)
    }
}

/// Existing solver results with the original material and exact snapshot kept.
pub struct MaterialChtSolution<'a, F: FluidTransport + ?Sized> {
    material: &'a F,
    transport: FluidTransportState,
    flow: FvFlow,
    energy: EnergySolution,
}
impl<'a, F: FluidTransport + ?Sized> MaterialChtSolution<'a, F> {
    /// Original material owner, not an inferred fluid name.
    #[must_use]
    pub const fn material(&self) -> &'a F { self.material }
    /// One source snapshot used by every coefficient in both solves.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// Existing flow, including fluxes, pressures and eddy-viscosity fields.
    #[must_use]
    pub const fn flow(&self) -> &FvFlow { &self.flow }
    /// Existing temperatures, interface heat-flow observers and energy balance.
    #[must_use]
    pub const fn energy(&self) -> &EnergySolution { &self.energy }
}

