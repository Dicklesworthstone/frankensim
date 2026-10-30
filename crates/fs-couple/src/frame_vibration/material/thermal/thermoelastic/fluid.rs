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
pub enum FrameFluidThermoelasticError {
    Thermal(FrameThermalError),
    Fluid(FrameFluidError),
}
impl core::fmt::Display for FrameFluidThermoelasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fluid-loaded frame thermodynamics: {self:?}")
    }
}
impl std::error::Error for FrameFluidThermoelasticError {}
impl From<FrameThermalError> for FrameFluidThermoelasticError {
    fn from(e: FrameThermalError) -> Self { Self::Thermal(e) }
}
impl From<FrameFluidError> for FrameFluidThermoelasticError {
    fn from(e: FrameFluidError) -> Self { Self::Fluid(e) }
}
type FluidResult<T> = core::result::Result<T, FrameFluidThermoelasticError>;
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
    /// layer. Original source pins, reference state, initial displacement and
    /// ALL fluid momenta survive. No eigensolve or Bessel quadrature is repeated.
    /// A warm material-exchange event is not inferred from this constructor.
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
        check(gate)?;
        positive(config.initial_energy_tolerance_j)?;
        let count = model.basis().modes().len();
        // This source factory admits complete mass space and reconstructs its
        // original K from layer laws BEFORE the full wet-coordinate congruence.
        let mut component = model.basis().coupled_component(config.sample_rate_hz,
            &vec![0.0; count], config.modal, config.maximum_projection_terms,
            || gate.is_requested()).map_err(frame_error)?;
        component.restore_states(&model.modal_states(initial_state)?).map_err(frame_error)?;
        let (mut binding, _, mut definitions) = component.into_thermoelastic_parts(
            thermal, source, destinations, config.maximum_absolute_strain, 0, 0,
            config.coupling, gate)?;
        let FrameFluidModalParts { coordinates, structural, memory, connections } =
            model.modal_parts(initial_state, config.sample_rate_hz, config.modal,
                definitions.len(), config.coupling, gate)?;
        for definition in &mut definitions {
            check(gate)?;
            definition.left.shapes = coordinates.project_column(&definition.left.shapes)?;
            definition.right.shapes = coordinates.project_column(&definition.right.shapes)?;
        }
        for member in &mut binding.members {
            check(gate)?;
            for column in &mut member.columns { *column = coordinates.project_column(column)?; }
        }
        // The source laws, reference lengths and initial temperature are the
        // same. Re-evaluate using the actual transformed columns/state rather
        // than transporting rounded free-energy/strain diagnostics blindly.
        let mut states = Vec::with_capacity(binding.members.len());
        for member in &binding.members {
            states.push(resolve_member(member, &temperatures(&member.region, thermal)?,
                structural.states(), structural.sample_period_s())?);
        }
        binding.states = states;
        binding.initial_model = structural.clone();
        binding.definitions = definitions.clone();
        definitions.extend(connections);
        let mut components = vec![structural];
        if let Some(memory) = memory { components.push(memory); }
        let mechanical = match config.prepared {
            Some(budget) => CoupledModalSystem::new_prepared(components, definitions,
                config.coupling, budget, gate),
            None => CoupledModalSystem::new(components, definitions, config.coupling, gate),
        }.map_err(FrameThermalError::from)?;
        let original = model.system().hamiltonian(initial_state);
        let composed = mechanical.total_energy_j().map_err(FrameThermalError::from)?;
        if finite(composed-original)?.abs() > config.initial_energy_tolerance_j {
            return Err(invalid("complete fluid/thermoelastic binding changed original initial storage").into());
        }
        let total = mechanical.mode_count();
        let inner = binding.into_system(mechanical, thermal, gate)?;
        check(gate)?;
        Ok(Self { coordinates, inner, dry_force: vec![0.0; count], native_force: vec![0.0; total] })
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
