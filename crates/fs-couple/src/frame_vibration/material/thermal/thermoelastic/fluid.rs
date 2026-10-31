//! Source thermoelasticity and finite heat with the ORIGINAL exterior fluid.
//!
//! This changes coordinates once, not physics or material laws. The existing
//! complete-frame thermoelastic owner keeps its per-layer E(T), integrated
//! expansion, Psi+TS, entropy, capacity and strain checks. Its SAME simultaneous
//! mechanical assembly now includes every original fluid momentum and dashpot.
//! Exterior-fluid loss is exported, never reclassified as finite-solid heat.

use super::*;
use crate::air_path::cylinder::frame::{FrameFluidError, MaterialFluidFrame};
use crate::air_path::cylinder::frame::modal::{FrameFluidCoordinates, FrameFluidModalParts};
use crate::modal_acoustic_time::ModalAcousticTimeBudget;
use crate::material_circuit::MaterialCircuit;
use fs_material::conductor::ElectricalPropertyBasis;
use crate::circuit_enthalpy::contacts::heat::{HeatAssembly, HeatNetworkError};

mod assembly;
pub use assembly::{FrameFluidThermoelasticBinding, FrameFluidThermoelasticParts};

/// Cold preparation controls. The native budgets cover ALL structural modes,
/// fluid-memory coordinates (including free position gauges), and source ports.
pub struct FrameFluidThermoelasticConfig {
    pub sample_rate_hz: u32,
    pub modal: ModalAcousticTimeBudget,
    pub coupling: ModalCouplingConfig,
    pub prepared: Option<PreparedModalBudget>,
    pub maximum_projection_terms: usize,
    pub maximum_absolute_strain: f64,
    /// Independent original-vs-composed initial mechanical free-energy check.
    pub initial_energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum FrameFluidThermoelasticError<E = HeatNetworkError> {
    Thermal(FrameThermalError),
    Fluid(FrameFluidError),
    Heat(FrameHeatError<E>),
}
impl<E: core::fmt::Debug> core::fmt::Display for FrameFluidThermoelasticError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fluid-loaded frame thermodynamics: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for FrameFluidThermoelasticError<E> {}
impl<E> From<FrameThermalError> for FrameFluidThermoelasticError<E> {
    fn from(e: FrameThermalError) -> Self { Self::Thermal(e) }
}
impl<E> From<FrameFluidError> for FrameFluidThermoelasticError<E> {
    fn from(e: FrameFluidError) -> Self { Self::Fluid(e) }
}
impl<E> From<FrameHeatError<E>> for FrameFluidThermoelasticError<E> {
    fn from(e: FrameHeatError<E>) -> Self {
        match e { FrameHeatError::Frame(error) => Self::Thermal(error), other => Self::Heat(other) }
    }
}
type FluidResult<T, E = HeatNetworkError> = core::result::Result<T, FrameFluidThermoelasticError<E>>;
fn frame_error(e: FrameVibrationError) -> FrameFluidThermoelasticError {
    match e {
        FrameVibrationError::Cancelled => FrameThermalError::Cancelled.into(),
        other => FrameThermalError::Source(FrameThermalSourceError::Frame(other)).into(),
    }
}

/// One original source model, one original thermoelastic/enthalpy transaction,
/// and one constant full-inertia chart. The underlying owner's dry basis is
/// intentionally not exposed as a map for native normalized states: use
/// `modal_states` or `reconstruct_into` for physical observations.
///
/// Complete structural modes, fixed reference mass/geometry and held exterior
/// transport are required. Temperature can change solid E, free strain and eta;
/// it cannot change the fluid density, reference radius or kinematic chart here.
/// This is linear reference-geometry thermoelasticity, not geometric buckling,
/// finite-deformation moving-boundary flow, melting or liquid self-heating.
pub struct FrameFluidThermoelasticSystem<'m, 'f, 'n, 'a, F> {
    coordinates: FrameFluidCoordinates<'m, 'f, F>,
    inner: FrameThermoelasticSystem<'n, 'a>,
    dry_force: Vec<f64>,
    native_force: Vec<f64>,
}

impl<'m, 'f, 'n, 'a, F> FrameFluidThermoelasticSystem<'m, 'f, 'n, 'a, F> {
    /// Cold-bind the original full frame/fluid state and every solid caloric
    /// layer through the same parts/validation path used for shared assemblies.
    /// Source pins, reference state and ALL fluid momenta survive. No eigensolve
    /// or Bessel quadrature is repeated. This is not a warm material exchange.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        model: &'m MaterialFluidFrame<'f, F>,
        initial_state: &[f64],
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        destinations: &[FrameLayerCell],
        config: FrameFluidThermoelasticConfig,
        gate: &CancelGate,
    ) -> FluidResult<Self> {
        Self::new_bound(model, initial_state, thermal, source, None, destinations, config, gate)
    }

    /// Original complete wet-frame preparation alongside explicitly owned
    /// finite coolant or other external thermal links. The original material
    /// graph validates their disjoint coverage; subsequent complete assembly
    /// steps must supply fresh coefficients from the original physical sources.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_link_conductances(
        model: &'m MaterialFluidFrame<'f, F>,
        initial_state: &[f64],
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        additional_links: &BTreeMap<u64, f64>,
        destinations: &[FrameLayerCell],
        config: FrameFluidThermoelasticConfig,
        gate: &CancelGate,
    ) -> FluidResult<Self> {
        Self::new_bound(model, initial_state, thermal, source, Some(additional_links),
            destinations, config, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn new_bound(
        model: &'m MaterialFluidFrame<'f, F>,
        initial_state: &[f64],
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        additional_links: Option<&BTreeMap<u64, f64>>,
        destinations: &[FrameLayerCell],
        config: FrameFluidThermoelasticConfig,
        gate: &CancelGate,
    ) -> FluidResult<Self> {
        check(gate)?;
        let count = model.basis().modes().len();
        let coupling = config.coupling;
        let prepared = config.prepared;
        let tolerance = config.initial_energy_tolerance_j;
        let FrameFluidThermoelasticParts { placement, binding, components, connections } =
            model.thermoelastic_parts_bound(initial_state, thermal, source, additional_links, destinations,
                config, 0, 0, gate)?;
        let mechanical = match prepared {
            Some(budget) => CoupledModalSystem::new_prepared(components, connections,
                coupling, budget, gate),
            None => CoupledModalSystem::new(components, connections, coupling, gate),
        }.map_err(FrameThermalError::from)?;
        let original = model.system().hamiltonian(initial_state);
        let composed = mechanical.total_energy_j().map_err(FrameThermalError::from)?;
        if finite(composed-original)?.abs() > tolerance {
            return Err(invalid("complete fluid/thermoelastic binding changed original initial storage").into());
        }
        let total = mechanical.mode_count();
        let inner = binding.into_system(mechanical, thermal, gate)?;
        check(gate)?;
        Ok(Self { coordinates: placement.into_coordinates(), inner,
            dry_force: vec![0.0; count], native_force: vec![0.0; total] })
    }

    pub const fn source(&self) -> &'m MaterialFluidFrame<'f, F> { self.coordinates.source() }
    /// Actual transformed native components; physical maps use `modal_states`.
    pub fn mechanical(&self) -> &CoupledModalSystem { self.inner.mechanical() }
    pub fn members(&self) -> &[FrameThermoelasticMemberState] { self.inner.members() }
    pub fn mechanical_internal_energy_j(&self) -> FluidResult<f64> {
        Ok(self.inner.mechanical_internal_energy_j()?)
    }
    /// Original dry-basis displacements/velocities, including thermal expansion.
    pub fn modal_states(&self) -> FluidResult<Vec<ModalAcousticState>> {
        Ok(self.coordinates.original_states(self.inner.mechanical().components()[0].states())?)
    }
    /// Reconstruct through the ORIGINAL source frame. Buffers are observation
    /// scratch; refusal never changes any accepted mechanics or thermal history.
    pub fn reconstruct_into(&self, displacement: &mut [f64], velocity: &mut [f64],
        gate: &CancelGate) -> FluidResult<()> {
        let state = self.modal_states()?;
        self.source().basis().reconstruct_into(&state, displacement, velocity,
            || gate.is_requested()).map_err(frame_error)
    }
    /// Complete original [dry q, loaded momentum, spectral fluid momenta] chart.
    /// Native free fluid positions are unobservable gauges, not discarded heat
    /// or mechanical memory. Do not re-use the source's FROZEN stiffness to
    /// evolve this state; continuing motion belongs to this live source owner.
    pub fn physical_state(&self) -> FluidResult<Vec<f64>> {
        let motion = self.modal_states()?;
        let masses = self.source().memory_masses_kg();
        let mut velocity = Vec::with_capacity(masses.len());
        if !masses.is_empty() {
            for (s, mass) in self.inner.mechanical().components()[1].states().iter().zip(masses) {
                velocity.push(finite(s.velocity_m_sqrt_kg_per_s/mass.sqrt())?);
            }
        }
        Ok(self.source().state_from_motion(&motion, &velocity)?)
    }

    fn project_forces(&mut self, nodal: &[f64], gate: &CancelGate) -> FluidResult<()> {
        check(gate)?;
        self.coordinates.source().basis().project_force_into(nodal, [0.0; 2],
            &mut self.dry_force, || gate.is_requested()).map_err(frame_error)?;
        let count = self.dry_force.len();
        self.coordinates.project_column_into(&self.dry_force, &mut self.native_force[..count])?;
        // No authored force acts on the fluid's private memory/gauge coordinates.
        self.native_force[count..].fill(0.0);
        check(gate)?;
        Ok(())
    }

    /// ONE exact frozen-start-temperature mechanical sample followed by the
    /// original joint caloric/source admission. E(T), expansion and eta are
    /// installed at unchanged accepted coordinates for the next sample.
    /// This retains the source owner's first-order thermoelastic split; no new
    /// integrator, residual-energy heat correction or constitutive law is added.
    /// Fluid loss remains an explicit exterior sink. On Err, BOTH histories and
    /// source coefficients remain at their accepted state; only scratch changes.
    pub fn step(&mut self, thermal: &mut DepositionSimulation<'_>, nodal_force: &[f64],
        controls: FrameThermalSolve, energy_guess_j: &[f64], gate: &CancelGate)
        -> FluidResult<FrameThermoelasticReport> {
        self.project_forces(nodal_force, gate)?;
        Ok(self.inner.step(thermal, &self.native_force, controls, energy_guess_j, gate)?)
    }

    /// Reuse the original same-cell electrical/thermal binding. Fixed electrical
    /// geometry and C/L, no motional EMF or field-force claim. A source swap
    /// changes R(T), heat, expansion/stiffness and the SAME wet-frame response.
    pub fn with_material_circuit(mut self, thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>, circuit: MaterialCircuit,
        bases: &[(usize, ElectricalPropertyBasis)], gate: &CancelGate) -> FluidResult<Self> {
        self.inner = self.inner.with_material_circuit(thermal, source, circuit, bases, gate)?;
        Ok(self)
    }
    /// Preserve the full exterior memory while admitting source-bound electrical
    /// interfaces and their original two-sided thermal receivers.
    pub fn with_material_contact_circuit(mut self, thermal: &DepositionSimulation<'_>,
        source: &'n crate::circuit_enthalpy::contacts::material::ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit, gate: &CancelGate) -> FluidResult<Self> {
        self.inner = self.inner.with_material_contact_circuit(thermal, source, circuit, gate)?;
        Ok(self)
    }

    /// Material C(T)/L(T) uses the same full-inertia frame/circuit/heat owner.
    /// This changes material response, not geometry or exterior transport.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_field_circuit(mut self, thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>, circuit: MaterialCircuit,
        bases: &[(usize, ElectricalPropertyBasis)],
        fields: Vec<crate::circuit_enthalpy::field_storage::FieldThermalBinding>,
        gate: &CancelGate) -> FluidResult<Self> {
        self.inner = self.inner.with_material_field_circuit(thermal, source, circuit, bases, fields, gate)?;
        Ok(self)
    }

    /// Same source-dependent field storage with original electrical interfaces.
    pub fn with_material_contact_field_circuit(mut self, thermal: &DepositionSimulation<'_>,
        source: &'n crate::circuit_enthalpy::contacts::material::ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        fields: Vec<crate::circuit_enthalpy::field_storage::FieldThermalBinding>,
        gate: &CancelGate) -> FluidResult<Self> {
        self.inner = self.inner.with_material_contact_field_circuit(thermal, source, circuit, fields, gate)?;
        Ok(self)
    }

    /// Original bulk circuit with explicit external thermal-link ownership.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_circuit_and_link_conductances(mut self, thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>, circuit: MaterialCircuit,
        bases: &[(usize, ElectricalPropertyBasis)], additional_links: &BTreeMap<u64, f64>,
        gate: &CancelGate) -> FluidResult<Self> {
        self.inner = self.inner.with_material_circuit_and_link_conductances(
            thermal, source, circuit, bases, additional_links, gate)?;
        Ok(self)
    }

    /// Original contact circuit with explicit external thermal-link ownership.
    pub fn with_material_contact_circuit_and_link_conductances(mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n crate::circuit_enthalpy::contacts::material::ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit, additional_links: &BTreeMap<u64, f64>, gate: &CancelGate,
    ) -> FluidResult<Self> {
        self.inner = self.inner.with_material_contact_circuit_and_link_conductances(
            thermal, source, circuit, additional_links, gate)?;
        Ok(self)
    }

    /// Original material C(T)/L(T) circuit with explicit external thermal links.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_field_circuit_and_link_conductances(mut self,
        thermal: &DepositionSimulation<'_>, source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit, bases: &[(usize, ElectricalPropertyBasis)],
        fields: Vec<crate::circuit_enthalpy::field_storage::FieldThermalBinding>,
        additional_links: &BTreeMap<u64, f64>, gate: &CancelGate,
    ) -> FluidResult<Self> {
        self.inner = self.inner.with_material_field_circuit_and_link_conductances(
            thermal, source, circuit, bases, fields, additional_links, gate)?;
        Ok(self)
    }

    /// Original contact and material-field circuit with external thermal links.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_contact_field_circuit_and_link_conductances(mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n crate::circuit_enthalpy::contacts::material::ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        fields: Vec<crate::circuit_enthalpy::field_storage::FieldThermalBinding>,
        additional_links: &BTreeMap<u64, f64>, gate: &CancelGate,
    ) -> FluidResult<Self> {
        self.inner = self.inner.with_material_contact_field_circuit_and_link_conductances(
            thermal, source, circuit, fields, additional_links, gate)?;
        Ok(self)
    }

    /// Full material heat transport acts on the same original fluid-loaded
    /// motion. Forces retain the original work-conjugate coordinate transform;
    /// source and finite heat checks finish before any native state publishes.
    /// The existing exterior fluid remains a held transport model and its loss
    /// remains an external sink, not radiation or caloric heat in the solids.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_network(&mut self, thermal: &mut DepositionSimulation<'_>,
        nodal_force: &[f64], heat: FrameHeatInput<'_, '_, '_, '_>,
        electrical: Option<FrameElectricalInput<'_>>, controls: FrameThermalSolve,
        energy_guess_j: &[f64], gate: &CancelGate) -> FluidResult<FrameThermoelasticReport> {
        self.project_forces(nodal_force, gate)?;
        Ok(self.inner.step_heat_network(thermal, &self.native_force, heat, electrical,
            controls, energy_guess_j, gate)?)
    }

    /// Retain complete assembly receipts through the original wet-frame force
    /// transform and joint frame/circuit/thermal publication boundary. The
    /// exterior mechanical fluid remains its original explicit loss owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_assembly<A: HeatAssembly + ?Sized>(
        &mut self, thermal: &mut DepositionSimulation<'_>, nodal_force: &[f64],
        assembly: &A, cx: &fs_exec::Cx<'_>, electrical: Option<FrameElectricalInput<'_>>,
        controls: FrameThermalSolve, energy_guess_j: &[f64], gate: &CancelGate,
    ) -> FluidResult<FrameThermoelasticReport<A::Sample>, A::Error> {
        self.project_forces(nodal_force, gate).map_err(|error| match error {
            FrameFluidThermoelasticError::Thermal(error) => FrameFluidThermoelasticError::Thermal(error),
            FrameFluidThermoelasticError::Fluid(error) => FrameFluidThermoelasticError::Fluid(error),
            FrameFluidThermoelasticError::Heat(FrameHeatError::Frame(error)) =>
                FrameFluidThermoelasticError::Thermal(error),
            FrameFluidThermoelasticError::Heat(FrameHeatError::Heat(error)) =>
                FrameFluidThermoelasticError::Heat(FrameHeatError::Heat(A::Error::from(error))),
        })?;
        Ok(self.inner.step_heat_assembly(thermal, &self.native_force, assembly, cx,
            electrical, controls, energy_guess_j, gate)?)
    }

    pub fn electrical(&self) -> Option<&MaterialCircuit> { self.inner.electrical() }
    /// Electrical work, solid caloric storage, thermoelastic work and fluid
    /// memory are accepted on the original SINGLE transaction boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn step_electrical(&mut self, thermal: &mut DepositionSimulation<'_>,
        nodal_force: &[f64], input: FrameElectricalInput<'_>, controls: FrameThermalSolve,
        energy_guess_j: &[f64], gate: &CancelGate) -> FluidResult<FrameThermoelasticReport> {
        self.project_forces(nodal_force, gate)?;
        Ok(self.inner.step_electrical(thermal, &self.native_force, input,
            controls, energy_guess_j, gate)?)
    }
}

#[cfg(test)]
mod tests;
