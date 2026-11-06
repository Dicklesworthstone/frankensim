//! Native winding circuits with retained three-dimensional vector hysteresis.
//!
//! Complete closed J/I profiles couple the original Whitney A–phi history to
//! native Kirchhoff ports through exact F^T delta(A) work. The shared spatial
//! terminal kernel retains current- and voltage-input causality, graph storage,
//! source-bound copper and ordinary resistors, and passive end self/mutual modes.
//! Every terminal/Jacobian trial starts from the SAME accepted cell memories;
//! remanence, reversals and vector rotation directly affect the winding response.
//!
//! Magnetic pinning, resolved eddy current, stranded copper and ordinary graph
//! losses are separate physical heat receipts. Backward-Euler endpoint losses
//! and end-mode numerical loss remain outside heat. Magnetic laws, geometry,
//! frames and occupancy stay fixed; original conductivity/resistivity sources
//! may refresh at admitted temperatures without restarting field or graph state.
//! This owner supplies no moving-memory transport or temperature-dependent fit.

pub use super::super::nonlinear::winding::{
    SpatialNonlinearWindingConfig, SpatialNonlinearWindings, SpatialWindingCellHeat,
    SpatialWindingConductor,
};
use super::super::{
    coupling::Geometry,
    nonlinear::winding::{SpatialNonlinearWindingError, SpatialNonlinearWindingStep},
    winding_kernel::{FieldAccount, FieldHistory, FieldTrial, Result, Terminal},
};
use super::{
    SpatialHystereticInduction, SpatialHystereticInductionControls, SpatialHystereticInductionStep,
};
use crate::winding::{
    end::{EndWindingInductance, EndWindingMode},
    nonlinear::{
        NonlinearWindingControls, NonlinearWindingInitial, NonlinearWindingPort,
        resistors::GraphResistorMaterial,
    },
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_phs::circuit::{CircuitDae, CircuitGraph};
use std::sync::Arc;

pub type SpatialHystereticWindingControls =
    NonlinearWindingControls<SpatialHystereticInductionControls>;
pub type SpatialHystereticWindingStep = SpatialNonlinearWindingStep<SpatialHystereticInductionStep>;
pub type SpatialHystereticWindingError = SpatialNonlinearWindingError;

/// Accepted original vector memory, original field cochains and graph storage.
/// Failed, interrupted or cancelled candidates leave all accepted state intact.
#[derive(Clone)]
pub struct SpatialHystereticWindingCircuit<'m, C: SpatialConductorMaterial = ResolvedConductorState>
{
    field: SpatialHystereticInduction<'m, C>,
    terminal: Terminal,
}

pub struct SpatialHystereticWindingCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    circuit: SpatialHystereticWindingCircuit<'m, C>,
    step: SpatialHystereticWindingStep,
}
impl<'m, C: SpatialConductorMaterial> SpatialHystereticWindingCandidate<'m, C> {
    pub fn circuit(&self) -> &SpatialHystereticWindingCircuit<'m, C> {
        &self.circuit
    }
    pub fn step(&self) -> &SpatialHystereticWindingStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialHystereticWindingCircuit<'m, C>,
        SpatialHystereticWindingStep,
    ) {
        (self.circuit, self.step)
    }
}

impl<C: SpatialConductorMaterial> FieldHistory for SpatialHystereticInduction<'_, C> {
    fn geometry(&self) -> Geometry<'_> {
        self.problem().geometry()
    }
    fn potential_wb(&self) -> &[f64] {
        self.potential_wb()
    }
    fn time_s(&self) -> f64 {
        self.time_s()
    }
    fn impressed_current_density_a_m2(&self) -> &[[f64; 3]] {
        self.impressed_current_density_a_m2()
    }
    fn conducting_cell(&self, cell: usize) -> bool {
        self.conductors()[cell].is_some()
    }
    fn current_load(
        &self,
        profile: &[[f64; 3]],
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        Ok(self.prepared.coupling.current_load(
            self.problem().geometry(),
            profile,
            tolerance,
            cancelled,
        )?)
    }
    fn stored_energy_j(&self, cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
        Ok(self.field(cancelled)?.stored_energy_j())
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialHystereticWindingCircuit<'m, C> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: SpatialHystereticInduction<'m, C>,
        graph: CircuitGraph,
        windings: SpatialNonlinearWindings,
        sources: Vec<SpatialWindingConductor>,
        initial: NonlinearWindingInitial,
        config: SpatialNonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let terminal_ports =
            vec![NonlinearWindingPort::CurrentInput; windings.branch_indices.len()];
        Self::new_with_ports(
            field,
            graph,
            windings,
            terminal_ports,
            sources,
            initial,
            config,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_with_ports(
        field: SpatialHystereticInduction<'m, C>,
        graph: CircuitGraph,
        windings: SpatialNonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        sources: Vec<SpatialWindingConductor>,
        initial: NonlinearWindingInitial,
        config: SpatialNonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_end_modes(
            field,
            graph,
            windings,
            terminal_ports,
            Vec::new(),
            sources,
            initial,
            config,
            cancelled,
        )
    }

    /// Attach explicit passive end-flux modes outside the resolved field domain.
    /// Signed participation retains winding polarity and permits shared or
    /// dependent modes without inventing leakage or inverting their storage.
    /// Initial currents must reproduce the accepted original impressed current;
    /// constructor admission never re-equilibrates or resets magnetic memory.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_end_modes(
        field: SpatialHystereticInduction<'m, C>,
        graph: CircuitGraph,
        windings: SpatialNonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        end_modes: Vec<EndWindingMode>,
        sources: Vec<SpatialWindingConductor>,
        initial: NonlinearWindingInitial,
        config: SpatialNonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let terminal = Terminal::new(
            &field,
            graph,
            windings,
            terminal_ports,
            end_modes,
            sources,
            initial,
            config,
            |cell, source| {
                field.problem().admit_conductor(cell, source)?;
                Ok(())
            },
            cancelled,
        )?;
        Ok(Self { field, terminal })
    }

    pub fn field(&self) -> &SpatialHystereticInduction<'m, C> {
        &self.field
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.terminal.layout.graph
    }
    pub fn dae(&self) -> &CircuitDae {
        &self.terminal.layout.dae
    }
    pub fn windings(&self) -> &SpatialNonlinearWindings {
        &self.terminal.layout.windings
    }
    pub fn branch_indices(&self) -> &[usize] {
        &self.terminal.layout.windings.branch_indices
    }
    pub fn profiles_per_m2(&self) -> &[Vec<[f64; 3]>] {
        &self.terminal.layout.windings.profiles_per_m2
    }
    /// Independent diagonal end contribution. Shared modes supply additional
    /// self and mutual terms in end_inductance_matrix_h.
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.terminal.layout.windings.end_inductance_h
    }
    pub fn end_inductance(&self) -> &EndWindingInductance {
        &self.terminal.end
    }
    /// Complete row-major end self/mutual matrix in original winding order [H].
    pub fn end_inductance_matrix_h(&self) -> &[f64] {
        self.terminal.end.matrix_h()
    }
    pub fn terminal_ports(&self) -> &[NonlinearWindingPort] {
        &self.terminal.layout.ports
    }
    pub fn conductors(&self) -> &[SpatialWindingConductor] {
        &self.terminal.sources.conductors
    }
    pub fn resistance_ohm(&self) -> &[f64] {
        &self.terminal.sources.resistance_ohm
    }
    pub fn state(&self) -> &[f64] {
        &self.terminal.state
    }
    pub fn currents_a(&self) -> &[f64] {
        &self.terminal.currents_a
    }
    pub fn flux_linkages_wb(&self) -> &[f64] {
        &self.terminal.flux_linkages_wb
    }
    pub fn time_s(&self) -> f64 {
        self.field.time_s()
    }
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        &self.terminal.graph_materials
    }

    /// Bind every ordinary native graph resistor to its actual material source.
    /// Geometry, port order and accepted field/circuit memory are retained.
    /// Thermal consumers must separately pin original sources and admit their
    /// temperature refresh; this method permits explicit source reconfiguration.
    pub fn with_graph_resistors(
        &self,
        sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Ok(Self {
            field: self.field.clone(),
            terminal: self.terminal.with_graph_resistors(sources, cancelled)?,
        })
    }

    pub fn step_to_candidate(
        &self,
        target_s: f64,
        external_inputs: &[f64],
        controls: SpatialHystereticWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticWindingCandidate<'m, C>> {
        self.step_prepared(
            target_s,
            Arc::clone(&self.field.prepared),
            None,
            external_inputs,
            controls,
            cancelled,
        )
    }

    /// Prepare original core conductivity and strand resistivity once, then
    /// share them across all terminal trials from the accepted vector history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_conductors_candidate(
        &self,
        target_s: f64,
        core_conductors: Vec<Option<C>>,
        winding_conductors: Vec<SpatialWindingConductor>,
        external_inputs: &[f64],
        controls: SpatialHystereticWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticWindingCandidate<'m, C>> {
        let prepared = self.field.prepare_conductors(core_conductors, cancelled)?;
        self.step_prepared(
            target_s,
            prepared,
            Some(winding_conductors),
            external_inputs,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &self,
        target_s: f64,
        prepared: Arc<super::Prepared<C>>,
        sources: Option<Vec<SpatialWindingConductor>>,
        external_inputs: &[f64],
        controls: SpatialHystereticWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticWindingCandidate<'m, C>> {
        let sources = match sources {
            None => std::rc::Rc::clone(&self.terminal.sources),
            Some(sources) => self.terminal.bind_sources(
                &self.field,
                sources,
                true,
                |cell, source| {
                    self.field.problem().admit_conductor(cell, source)?;
                    Ok(())
                },
                cancelled,
            )?,
        };
        let candidate = self.terminal.step(
            &self.field,
            target_s,
            sources,
            external_inputs,
            controls,
            |current, dt, check| {
                let (state, step) = self
                    .field
                    .step_prepared_candidate(
                        Arc::clone(&prepared),
                        current,
                        dt,
                        controls.field,
                        &mut || check(),
                    )?
                    .into_parts();
                let account = FieldAccount {
                    stored_before_j: step.stored_before_j,
                    stored_after_j: step.stored_after_j,
                    stored_change_j: step.stored_change_j,
                    heat_j: step.heat_j,
                    impressed_work_j: step.impressed_work_j,
                    numerical_dissipation_j: step.numerical_dissipation_j,
                    mechanical_work_j: 0.0,
                    parameter_work_j: 0.0,
                };
                Ok(FieldTrial {
                    state,
                    step,
                    account,
                })
            },
            cancelled,
        )?;
        Ok(SpatialHystereticWindingCandidate {
            circuit: Self {
                field: candidate.field,
                terminal: candidate.terminal,
            },
            step: candidate.step,
        })
    }
}
