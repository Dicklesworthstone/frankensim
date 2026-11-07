//! Native RLC circuits driven by the original saturating solid-conductor field.
//!
//! Every nonlinear interface trial solves the original A–phi electrode problem
//! from the same accepted history. The native circuit owns its original fluxes,
//! charges, mutual inductances and transformers. Actual electrode reactions
//! supply circuit current; field storage, Ohmic heat and Lorentz work retain
//! their original spatial receipts. No affine response is reused across magnetic
//! saturation, and the numerical reference impedance cancels from all physics.

use super::super::{
    MagneticOwner, Prepared, SpatialNonlinearInductionControls, SpatialNonlinearInductionStep,
};
use super::SpatialNonlinearTerminalInduction;
use crate::spatial::nonlinear::SpatialNonlinearMaterial;
use crate::spatial::vector_potential::induction::poll;
pub use crate::spatial::vector_potential::induction::terminals::circuit::SpatialTerminalCircuitError as SpatialNonlinearTerminalCircuitError;
use crate::spatial::vector_potential::induction::terminals::circuit::{
    SpatialTerminalCircuitConfig, SpatialTerminalCircuitInitial,
    SpatialTerminalCircuitNonlinearControls, SpatialTerminalCircuitStep, nonlinear::Circuit,
};
use crate::winding::nonlinear::resistors::GraphResistorMaterial;
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_phs::circuit::{CircuitDae, CircuitGraph, MutualInductance};
use std::sync::Arc;

pub type SpatialNonlinearTerminalCircuitControls =
    SpatialTerminalCircuitNonlinearControls<SpatialNonlinearInductionControls>;
pub type SpatialNonlinearTerminalCircuitStep =
    SpatialTerminalCircuitStep<SpatialNonlinearInductionStep>;
type Result<T> = core::result::Result<T, SpatialNonlinearTerminalCircuitError>;

/// Immutable accepted field and native circuit history. Candidate evaluation
/// and material-bound resistance updates cannot advance the original owner.
#[derive(Clone)]
pub struct SpatialNonlinearTerminalCircuit<'m, C: SpatialConductorMaterial = ResolvedConductorState>
{
    field: SpatialNonlinearTerminalInduction<'m, C>,
    circuit: Circuit,
}

pub struct SpatialNonlinearTerminalCircuitCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    circuit: SpatialNonlinearTerminalCircuit<'m, C>,
    step: SpatialNonlinearTerminalCircuitStep,
}

impl<'m, C: SpatialConductorMaterial> SpatialNonlinearTerminalCircuitCandidate<'m, C> {
    pub fn circuit(&self) -> &SpatialNonlinearTerminalCircuit<'m, C> {
        &self.circuit
    }
    pub fn step(&self) -> &SpatialNonlinearTerminalCircuitStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialNonlinearTerminalCircuit<'m, C>,
        SpatialNonlinearTerminalCircuitStep,
    ) {
        (self.circuit, self.step)
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialNonlinearTerminalCircuit<'m, C> {
    /// Attach complete original electrode nodes without resetting the field or
    /// projecting original capacitor charge and mutually coupled flux linkage.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: SpatialNonlinearTerminalInduction<'m, C>,
        graph: CircuitGraph,
        terminal_nodes: Vec<usize>,
        mutual: Vec<MutualInductance>,
        initial: SpatialTerminalCircuitInitial,
        config: SpatialTerminalCircuitConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let circuit = Circuit::new(
            field.terminals.references(),
            field.problem().mesh().vertex_count,
            graph,
            terminal_nodes,
            mutual,
            initial,
            config,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self { field, circuit })
    }

    pub fn field(&self) -> &SpatialNonlinearTerminalInduction<'m, C> {
        &self.field
    }
    pub fn graph(&self) -> &CircuitGraph {
        self.circuit.graph()
    }
    pub fn dae(&self) -> &CircuitDae {
        self.circuit.dae()
    }
    pub fn mutual_inductance(&self) -> &[MutualInductance] {
        self.circuit.mutual_inductance()
    }
    pub fn terminal_nodes(&self) -> &[usize] {
        self.circuit.terminal_nodes()
    }
    pub fn native_state(&self) -> &[f64] {
        self.circuit.native_state()
    }
    pub fn time_s(&self) -> f64 {
        self.field.time_s()
    }
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        self.circuit.graph_resistor_materials()
    }

    /// Reassemble only authored resistor values from complete original-branch
    /// bindings. Native indices, source units, storage and field history remain.
    pub fn with_graph_resistors(
        &self,
        sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let circuit = self.circuit.with_graph_resistors(sources, cancelled)?;
        let field = self.field.clone();
        poll(cancelled)?;
        Ok(Self { field, circuit })
    }

    /// Native midpoint effort/current and endpoint field reactions share the
    /// same finite-step work. The coupled temporal scheme is first order.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialNonlinearTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalCircuitCandidate<'m, C>> {
        self.step_prepared(
            self.field.field.problem.clone(),
            Arc::clone(&self.field.field.prepared),
            external_inputs,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Prepare original conductivity exactly once before all bounded terminal
    /// Newton trials, preserving the original charge space and source identity.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_conductors_candidate(
        &self,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalCircuitCandidate<'m, C>> {
        let prepared = self.field.field.prepare_conductors(conductors, cancelled)?;
        self.step_prepared(
            self.field.field.problem.clone(),
            prepared,
            external_inputs,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Resolve the original reversible magnetic and electrical sources once.
    /// Parameter work remains explicit; finite thermal consumers must close its
    /// matching magnetic-entropy exchange instead of depositing it as heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_candidate(
        &self,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Vec<SpatialNonlinearMaterial>,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalCircuitCandidate<'m, C>> {
        let (problem, prepared) = self
            .field
            .field
            .prepare_materials(materials, conductors, cancelled)?;
        self.step_prepared(
            problem,
            prepared,
            external_inputs,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &self,
        problem: MagneticOwner<'m>,
        prepared: Arc<Prepared<C>>,
        external_inputs: &[f64],
        impressed: &[[f64; 3]],
        velocity: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialNonlinearTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalCircuitCandidate<'m, C>> {
        let next = self.circuit.step(
            external_inputs,
            duration_s,
            controls,
            |voltage, check| {
                Ok(self
                    .field
                    .step_prepared(
                        problem.clone(),
                        Arc::clone(&prepared),
                        voltage,
                        impressed,
                        velocity,
                        duration_s,
                        controls.field,
                        &mut || check(),
                    )?
                    .into_parts())
            },
            cancelled,
        )?;
        Ok(SpatialNonlinearTerminalCircuitCandidate {
            circuit: Self {
                field: next.field,
                circuit: next.circuit,
            },
            step: next.step,
        })
    }
}
