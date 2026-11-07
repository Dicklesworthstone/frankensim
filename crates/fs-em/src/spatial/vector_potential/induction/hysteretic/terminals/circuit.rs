//! Native RLC circuits driving distributed current and retained vector hysteresis.
//!
//! The existing nonlinear terminal core owns graph preparation, Newton closure,
//! physical branch work and cancellation of numerical reference resistors.
//! This adapter supplies only the original accepted electrode field and its
//! source-prepared candidate advance. Every trial sees the same vector memory.

pub use super::super::super::terminals::circuit::SpatialTerminalCircuitError as SpatialHystereticTerminalCircuitError;
use super::super::super::{
    poll,
    terminals::circuit::{
        SpatialTerminalCircuitConfig, SpatialTerminalCircuitInitial,
        SpatialTerminalCircuitNonlinearControls, SpatialTerminalCircuitStep, nonlinear::Circuit,
    },
};
use super::super::{Prepared, SpatialHystereticInductionControls, SpatialHystereticInductionStep};
use super::SpatialHystereticTerminalInduction;
use crate::winding::nonlinear::resistors::GraphResistorMaterial;
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_phs::circuit::{CircuitDae, CircuitGraph, MutualInductance};
use std::sync::Arc;

pub type SpatialHystereticTerminalCircuitControls =
    SpatialTerminalCircuitNonlinearControls<SpatialHystereticInductionControls>;
pub type SpatialHystereticTerminalCircuitStep =
    SpatialTerminalCircuitStep<SpatialHystereticInductionStep>;
type Result<T> = core::result::Result<T, SpatialHystereticTerminalCircuitError>;

/// Original circuit fluxes/charges, electrode field and every accepted memory.
#[derive(Clone)]
pub struct SpatialHystereticTerminalCircuit<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    field: SpatialHystereticTerminalInduction<'m, C>,
    circuit: Circuit,
}

/// An immutable jointly advanced native circuit and hysteretic electrode field.
pub struct SpatialHystereticTerminalCircuitCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    circuit: SpatialHystereticTerminalCircuit<'m, C>,
    step: SpatialHystereticTerminalCircuitStep,
}

impl<'m, C: SpatialConductorMaterial> SpatialHystereticTerminalCircuitCandidate<'m, C> {
    pub fn circuit(&self) -> &SpatialHystereticTerminalCircuit<'m, C> {
        &self.circuit
    }
    pub fn step(&self) -> &SpatialHystereticTerminalCircuitStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialHystereticTerminalCircuit<'m, C>,
        SpatialHystereticTerminalCircuitStep,
    ) {
        (self.circuit, self.step)
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialHystereticTerminalCircuit<'m, C> {
    /// Connect each original electrode to a native graph node. Mutual storage,
    /// ideal transformers, source-port units and genuine shorts use the shared
    /// graph admission. Initial fluxes and charges retain their authored values.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: SpatialHystereticTerminalInduction<'m, C>,
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

    pub fn field(&self) -> &SpatialHystereticTerminalInduction<'m, C> {
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

    /// Rebind every authored graph resistor at the same accepted time. Only
    /// dissipative coefficients change; native fluxes/charges, electrode maps
    /// and magnetic memories are retained. Finite-thermal callers additionally
    /// pin the original source choices and accepted caloric temperature.
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

    /// Advance the original RLC graph and physical electrode field together.
    /// Native midpoint efforts/flows share the field endpoint terminal work;
    /// the combined nonlinear method remains first order. Every trial restarts
    /// the same accepted vector memories, field and native storage coordinates.
    pub fn step_candidate(
        &self,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialHystereticTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticTerminalCircuitCandidate<'m, C>> {
        self.step_prepared(
            Arc::clone(&self.field.field.prepared),
            external_inputs,
            endpoint_impressed_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Prepare original temperature-refreshed conductor operators once, then
    /// share them across all terminal/Jacobian trials. Fitted magnetic laws,
    /// topology, accepted polarization and circuit storage remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_conductors_candidate(
        &self,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialHystereticTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticTerminalCircuitCandidate<'m, C>> {
        let prepared = self.field.field.prepare_conductors(conductors, cancelled)?;
        self.step_prepared(
            prepared,
            external_inputs,
            endpoint_impressed_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &self,
        prepared: Arc<Prepared<C>>,
        external_inputs: &[f64],
        impressed: &[[f64; 3]],
        dt: f64,
        controls: SpatialHystereticTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticTerminalCircuitCandidate<'m, C>> {
        let candidate = self.circuit.step(
            external_inputs,
            dt,
            controls,
            |voltage, checkpoint| {
                self.field
                    .step_prepared(
                        Arc::clone(&prepared),
                        voltage,
                        impressed,
                        dt,
                        controls.field,
                        &mut || checkpoint(),
                    )
                    .map(|trial| trial.into_parts())
                    .map_err(Into::into)
            },
            cancelled,
        )?;
        Ok(SpatialHystereticTerminalCircuitCandidate {
            circuit: Self {
                field: candidate.field,
                circuit: candidate.circuit,
            },
            step: candidate.step,
        })
    }
}
