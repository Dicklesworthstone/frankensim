//! Native circuit terminals on the original nonlinear three-dimensional A–phi history.
//!
//! A closed, signed, cellwise J/I profile defines each stranded winding. Its
//! exact Whitney load F is assembled once on the retained mesh. Every terminal
//! trial continues the SAME accepted field and native circuit states; the field
//! is neither re-equilibrated nor replaced by a fitted inductance. Backward
//! Euler uses endpoint winding currents and V = F^T delta(A)/dt + L_end delta(I)/dt
//! + R I. Native power ports cancel this transfer against the original Kirchhoff
//! graph. Reversible saturation, eddy diffusion and prescribed motion in separate
//! conducting regions therefore affect terminal voltage and current directly.
//! Explicit end-region self/mutual storage uses passive shared-flux modes; it
//! adds only caller-declared flux outside the resolved three-dimensional field.
//!
//! Stranded profiles include their physical turns and complete return paths.
//! Their source-bound scalar resistivity and explicit conductor filling give
//! R = integral(rho |J/I|^2 / filling) dV. Distinct insulated bundles can share a
//! cell only within its total filling; they cannot also be resolved eddy cells.
//! Strand geometry is stationary. Moving finite interfaces, solid terminal
//! current injection, hysteresis and homogenization error estimates are outside
//! this owner. Native algebraic/physical gates may refuse a singular topology.
//! Material refresh retains original source claims and accepted cochains;
//! magnetic parameter work and numerical dissipation are never labeled heat.

/// Source-bound ordinary resistors on the same native graph and field history.
pub mod resistors;

use super::super::{
    coupling::Geometry,
    winding_kernel::{FieldAccount, FieldHistory, FieldTrial, Result, Terminal, invalid, poll},
};
use super::{
    SpatialNonlinearInduction, SpatialNonlinearInductionControls, SpatialNonlinearInductionError,
    SpatialNonlinearInductionStep,
};
use crate::{
    spatial::{SpatialMagneticError, nonlinear::SpatialNonlinearMaterial},
    winding::{
        end::{EndWindingError, EndWindingInductance, EndWindingMode},
        nonlinear::{
            NonlinearWindingControls, NonlinearWindingError, NonlinearWindingInitial,
            NonlinearWindingPort, resistors::GraphResistorMaterial,
        },
        ohmic::{OhmicWindingError, WindingGraphResistorHeat},
    },
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_phs::{
    StepRecord,
    circuit::{CircuitDae, CircuitError, CircuitGraph},
};
use fs_solver::{NewtonError, NewtonReport, SolverRunProgress};

#[derive(Clone, Debug)]
pub struct SpatialNonlinearWindings {
    /// Distinct original inductor placeholders, in caller winding order.
    pub branch_indices: Vec<usize>,
    /// Signed complete distributed current divided by terminal current [1/m2].
    /// The values are averages over the represented cells, including filling.
    pub profiles_per_m2: Vec<Vec<[f64; 3]>>,
    /// Independent omitted series self-inductance; zero is admitted [H].
    /// Shared end-flux modes are supplied separately to new_with_end_modes.
    pub end_inductance_h: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct SpatialNonlinearWindingConfig {
    pub max_windings: usize,
    /// Bounds the profiles plus cached Whitney loads: n*(3*cells + edges).
    pub max_profile_values: usize,
    /// Maximum face-current imbalance in each unit-current profile [turns].
    pub current_balance_tolerance_turns: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpatialWindingConductor {
    pub winding: usize,
    pub cell: usize,
    /// Actual stranded conductor volume / represented cell volume in (0,1].
    /// Different entries in a cell are insulated physical bundles.
    pub fill_fraction: f64,
    pub conductor: ResolvedConductorState,
}

#[derive(Clone, Debug)]
pub struct SpatialWindingCellHeat {
    pub source: SpatialWindingConductor,
    pub winding: usize,
    pub cell: usize,
    pub volume_m3: f64,
    pub conductor_volume_m3: f64,
    /// Intrinsic current in the conductor, not its cell-average [A/m2].
    pub current_density_a_m2: [f64; 3],
    pub electric_field_v_m: [f64; 3],
    pub power_w: f64,
    pub heat_j: f64,
}

#[derive(Debug)]
pub enum SpatialNonlinearWindingError {
    Invalid(&'static str),
    Field(SpatialNonlinearInductionError),
    #[cfg(feature = "magnetic-hysteresis")]
    HystereticField(super::super::hysteretic::SpatialHystereticInductionError),
    Spatial(SpatialMagneticError),
    Circuit(CircuitError),
    Interface(NonlinearWindingError),
    Conductor(OhmicWindingError),
    End(EndWindingError),
    Newton(NewtonError),
    TrialBudget,
    Interrupted(SolverRunProgress),
    NotConverged(NewtonReport),
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for SpatialNonlinearWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear spatial winding: {self:?}")
    }
}
impl std::error::Error for SpatialNonlinearWindingError {}
impl From<SpatialMagneticError> for SpatialNonlinearWindingError {
    fn from(value: SpatialMagneticError) -> Self {
        match value {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Spatial(other),
        }
    }
}
impl From<SpatialNonlinearInductionError> for SpatialNonlinearWindingError {
    fn from(value: SpatialNonlinearInductionError) -> Self {
        match value {
            SpatialNonlinearInductionError::Spatial(SpatialMagneticError::Cancelled) => {
                Self::Cancelled
            }
            SpatialNonlinearInductionError::Magnetic(super::NonlinearVectorError::Spatial(
                SpatialMagneticError::Cancelled,
            )) => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<CircuitError> for SpatialNonlinearWindingError {
    fn from(value: CircuitError) -> Self {
        Self::Circuit(value)
    }
}
#[cfg(feature = "magnetic-hysteresis")]
impl From<super::super::hysteretic::SpatialHystereticInductionError>
    for SpatialNonlinearWindingError
{
    fn from(value: super::super::hysteretic::SpatialHystereticInductionError) -> Self {
        use super::super::hysteretic::SpatialHystereticInductionError;
        use fs_material::field_storage::hysteresis::vector::VectorHysteresisError;
        match value {
            SpatialHystereticInductionError::Spatial(SpatialMagneticError::Cancelled)
            | SpatialHystereticInductionError::Material(VectorHysteresisError::Cancelled) => {
                Self::Cancelled
            }
            other => Self::HystereticField(other),
        }
    }
}
impl From<NonlinearWindingError> for SpatialNonlinearWindingError {
    fn from(value: NonlinearWindingError) -> Self {
        match value {
            NonlinearWindingError::Cancelled => Self::Cancelled,
            other => Self::Interface(other),
        }
    }
}
impl From<OhmicWindingError> for SpatialNonlinearWindingError {
    fn from(value: OhmicWindingError) -> Self {
        match value {
            OhmicWindingError::Cancelled => Self::Cancelled,
            other => Self::Conductor(other),
        }
    }
}
impl From<NewtonError> for SpatialNonlinearWindingError {
    fn from(value: NewtonError) -> Self {
        Self::Newton(value)
    }
}
impl From<EndWindingError> for SpatialNonlinearWindingError {
    fn from(value: EndWindingError) -> Self {
        Self::End(value)
    }
}

/// Accepted original field and original graph storage. Failed or cancelled
/// candidates never mutate A, conductor sources, circuit flux/charge or time.
#[derive(Clone)]
pub struct SpatialNonlinearWindingCircuit<'m, C: SpatialConductorMaterial = ResolvedConductorState>
{
    terminal: Terminal,
    field: SpatialNonlinearInduction<'m, C>,
}

pub struct SpatialNonlinearWindingStep<F = SpatialNonlinearInductionStep> {
    /// Includes private terminal transfers. Combined supply is supplied_work_j.
    pub native: StepRecord,
    pub field: F,
    pub cells: Vec<SpatialWindingCellHeat>,
    pub graph_resistors: Vec<WindingGraphResistorHeat>,
    /// Source bindings for ordinary graph losses, in original branch order.
    /// Empty for an explicitly authored numeric-resistor graph.
    pub graph_materials: Vec<GraphResistorMaterial>,
    pub currents_a: Vec<f64>,
    pub terminal_voltage_v: Vec<f64>,
    pub inductive_voltage_v: Vec<f64>,
    pub terminal_residual_v: Vec<f64>,
    pub terminal_current_residual_a: Vec<f64>,
    pub flux_linkages_wb: [Vec<f64>; 2],
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Stable magnetic delta plus native graph and explicit end-storage deltas.
    pub stored_change_j: f64,
    pub heat_j: f64,
    pub supplied_work_j: f64,
    pub mechanical_work_j: f64,
    pub parameter_work_j: f64,
    pub numerical_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub field_trials: usize,
    pub report: NewtonReport,
}
pub struct SpatialNonlinearWindingCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    circuit: SpatialNonlinearWindingCircuit<'m, C>,
    step: SpatialNonlinearWindingStep,
}
impl<'m, C: SpatialConductorMaterial> SpatialNonlinearWindingCandidate<'m, C> {
    pub fn circuit(&self) -> &SpatialNonlinearWindingCircuit<'m, C> {
        &self.circuit
    }
    pub fn step(&self) -> &SpatialNonlinearWindingStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialNonlinearWindingCircuit<'m, C>,
        SpatialNonlinearWindingStep,
    ) {
        (self.circuit, self.step)
    }
}

impl<C: SpatialConductorMaterial> FieldHistory for SpatialNonlinearInduction<'_, C> {
    fn geometry(&self) -> Geometry<'_> {
        super::geometry(self.problem())
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
            super::geometry(self.problem()),
            profile,
            tolerance,
            cancelled,
        )?)
    }
    fn stored_energy_j(&self, cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
        Ok(self.field(cancelled)?.stored_energy_j())
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialNonlinearWindingCircuit<'m, C> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: SpatialNonlinearInduction<'m, C>,
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
        field: SpatialNonlinearInduction<'m, C>,
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

    /// Couple the original field/circuit history with explicit end-flux modes.
    /// Each mode adds ell*n*n^T to the independent end self-inductances.
    /// Signed participation carries winding polarity. Shared or dependent
    /// modes are admitted without inventing leakage or inverting end storage.
    /// These modes must describe flux omitted from the resolved field domain.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_end_modes(
        field: SpatialNonlinearInduction<'m, C>,
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
                field.problem().materials()[cell]
                    .admit_fixed_material_refresh(source.material())
                    .map_err(SpatialMagneticError::Material)?;
                Ok(())
            },
            cancelled,
        )?;
        Ok(Self { terminal, field })
    }

    pub fn field(&self) -> &SpatialNonlinearInduction<'m, C> {
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
    /// Independent diagonal end contribution; shared modes add their own self
    /// and mutual terms. Use end_inductance_matrix_h for the complete block.
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.terminal.layout.windings.end_inductance_h
    }
    /// The admitted original end model, retained across material/thermal refresh.
    pub fn end_inductance(&self) -> &EndWindingInductance {
        &self.terminal.end
    }
    /// Full row-major end self/mutual matrix in original winding order [H].
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

    /// Continue original sources. Velocity is prescribed only in resolved bulk
    /// conductors; all vertices of a stranded winding remain stationary.
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        velocity_m_s: &[[f64; 3]],
        external_inputs: &[f64],
        controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearWindingCandidate<'m, C>> {
        self.step_prepared(
            target_s,
            velocity_m_s,
            self.field.problem.clone(),
            std::sync::Arc::clone(&self.field.prepared),
            None,
            external_inputs,
            controls,
            cancelled,
        )
    }

    /// Refresh original conductor sources at owned temperature with fixed
    /// magnetic laws. No field or circuit history is restarted.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_conductors_candidate(
        &self,
        target_s: f64,
        velocity_m_s: &[[f64; 3]],
        core_conductors: Vec<Option<C>>,
        winding_conductors: Vec<SpatialWindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearWindingCandidate<'m, C>> {
        let prepared = if core_conductors == self.field.conductors() {
            std::sync::Arc::clone(&self.field.prepared)
        } else {
            let coupling =
                self.field
                    .coupling_at(self.field.problem(), core_conductors, cancelled)?;
            std::sync::Arc::new(super::Prepared {
                coupling,
                permeability_bound: self.field.prepared.permeability_bound.clone(),
            })
        };
        self.step_prepared(
            target_s,
            velocity_m_s,
            self.field.problem.clone(),
            prepared,
            Some(winding_conductors),
            external_inputs,
            controls,
            cancelled,
        )
    }

    /// Requery original magnetic and electrical sources. Endpoint operators are
    /// prepared ONCE and shared by every terminal Newton/Jacobian trial. The
    /// field returns explicit magnetic parameter work; a finite thermal caller
    /// must close that with its original magnetic entropy and heat receivers.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_candidate(
        &self,
        target_s: f64,
        velocity_m_s: &[[f64; 3]],
        materials: Vec<SpatialNonlinearMaterial>,
        core_conductors: Vec<Option<C>>,
        winding_conductors: Vec<SpatialWindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearWindingCandidate<'m, C>> {
        let (problem, prepared) =
            self.field
                .prepare_materials(materials, core_conductors, cancelled)?;
        self.step_prepared(
            target_s,
            velocity_m_s,
            problem,
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
        velocity_m_s: &[[f64; 3]],
        problem: super::MagneticOwner<'m>,
        prepared: std::sync::Arc<super::Prepared<C>>,
        sources: Option<Vec<SpatialWindingConductor>>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearWindingCandidate<'m, C>> {
        poll(cancelled)?;
        if velocity_m_s.len() != self.field.problem().mesh().vertex_count
            || velocity_m_s.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("complete finite prescribed velocity required"));
        }
        let sources = match sources {
            None => std::rc::Rc::clone(&self.terminal.sources),
            Some(sources) => self.terminal.bind_sources(
                &self.field,
                sources,
                problem.materials() == self.field.problem().materials(),
                |cell, source| {
                    problem.materials()[cell]
                        .admit_fixed_material_refresh(source.material())
                        .map_err(SpatialMagneticError::Material)?;
                    Ok(())
                },
                cancelled,
            )?,
        };
        for source in &sources.conductors {
            poll(cancelled)?;
            if self.field.problem().mesh().tets[source.cell]
                .iter()
                .any(|&v| velocity_m_s[v as usize].iter().any(|&x| x != 0.0))
            {
                return Err(invalid("stranded winding geometry must remain stationary"));
            }
        }
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
                        problem.clone(),
                        std::sync::Arc::clone(&prepared),
                        current,
                        velocity_m_s,
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
                    mechanical_work_j: step.mechanical_work_j,
                    parameter_work_j: step.parameter_work_j,
                };
                Ok(FieldTrial {
                    state,
                    step,
                    account,
                })
            },
            cancelled,
        )?;
        Ok(SpatialNonlinearWindingCandidate {
            circuit: Self {
                field: candidate.field,
                terminal: candidate.terminal,
            },
            step: candidate.step,
        })
    }
}
