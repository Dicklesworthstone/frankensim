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

use super::{
    SpatialNonlinearInduction, SpatialNonlinearInductionCandidate,
    SpatialNonlinearInductionControls, SpatialNonlinearInductionError,
    SpatialNonlinearInductionStep,
};
use crate::{
    spatial::{SpatialMagneticError, nonlinear::SpatialNonlinearMaterial},
    winding::{
        nonlinear::{
            NonlinearWindingControls, NonlinearWindingError, NonlinearWindingInitial,
            NonlinearWindingPort, native_seed, resistors::GraphResistorMaterial,
        },
        ohmic::{self, OhmicWindingError, WindingGraphResistorHeat},
    },
};
use fs_material::{
    conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial},
    field_storage::admit_source_temperature_refresh,
};
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitError, CircuitGraph, assemble_circuit},
};
use fs_solver::{
    NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress,
    finite_difference::{FallibleFiniteDifference, FiniteDifferenceError},
};
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
};

#[derive(Clone, Debug)]
pub struct SpatialNonlinearWindings {
    /// Distinct original inductor placeholders, in caller winding order.
    pub branch_indices: Vec<usize>,
    /// Signed complete distributed current divided by terminal current [1/m2].
    /// The values are averages over the represented cells, including filling.
    pub profiles_per_m2: Vec<Vec<[f64; 3]>>,
    /// Explicit omitted series self-inductance; zero is admitted [H].
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
    Spatial(SpatialMagneticError),
    Circuit(CircuitError),
    Interface(NonlinearWindingError),
    Conductor(OhmicWindingError),
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
type Result<T> = core::result::Result<T, SpatialNonlinearWindingError>;
fn finite(value: f64) -> Result<f64> {
    Ok(super::finite(value)?)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    Ok(super::poll(cancelled)?)
}
fn invalid(message: &'static str) -> SpatialNonlinearWindingError {
    SpatialNonlinearWindingError::Invalid(message)
}
fn numerical(
    error: FiniteDifferenceError<SpatialNonlinearWindingError>,
) -> SpatialNonlinearWindingError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => invalid(message),
    }
}
fn pairing(a: &[f64], b: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
    if a.len() != b.len() {
        return Err(invalid("complete Whitney source pairing required"));
    }
    a.iter().zip(b).try_fold(0.0, |sum, (&a, &b)| {
        poll(cancelled)?;
        finite(a.mul_add(b, sum))
    })
}
fn impressed(
    profiles: &[Vec<[f64; 3]>],
    currents: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<[f64; 3]>> {
    if currents.len() != profiles.len() || currents.iter().any(|v| !v.is_finite()) {
        return Err(invalid("complete finite winding currents required"));
    }
    let mut values = vec![[0.0; 3]; profiles[0].len()];
    for (profile, &current) in profiles.iter().zip(currents) {
        for (value, density) in values.iter_mut().zip(profile) {
            poll(cancelled)?;
            for j in 0..3 {
                value[j] = finite(density[j].mul_add(current, value[j]))?;
            }
        }
    }
    Ok(values)
}
fn end_energy(inductances: &[f64], currents: &[f64]) -> Result<f64> {
    inductances
        .iter()
        .zip(currents)
        .try_fold(0.0, |sum, (&l, &i)| finite((0.5 * l * i).mul_add(i, sum)))
}
struct Layout {
    graph: CircuitGraph,
    dae: CircuitDae,
    windings: Rc<SpatialNonlinearWindings>,
    ports: Vec<NonlinearWindingPort>,
    external_ports: usize,
    loads_per_a: Rc<[Vec<f64>]>,
}

/// Accepted original field and original graph storage. Failed or cancelled
/// candidates never mutate A, conductor sources, circuit flux/charge or time.
#[derive(Clone)]
pub struct SpatialNonlinearWindingCircuit<'m, C: SpatialConductorMaterial = ResolvedConductorState>
{
    layout: Rc<Layout>,
    graph_materials: Vec<GraphResistorMaterial>,
    field: SpatialNonlinearInduction<'m, C>,
    conductors: Vec<SpatialWindingConductor>,
    cell_resistance_ohm: Vec<f64>,
    resistance_ohm: Vec<f64>,
    state: Vec<f64>,
    currents_a: Vec<f64>,
    flux_linkages_wb: Vec<f64>,
}

pub struct SpatialNonlinearWindingStep {
    /// Includes private terminal transfers. Combined supply is supplied_work_j.
    pub native: StepRecord,
    pub field: SpatialNonlinearInductionStep,
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

fn ports(layout: &Layout, external: &[f64], currents: &[f64], voltage: &[f64]) -> Result<Vec<f64>> {
    if external.len() != layout.external_ports
        || external.iter().any(|x| !x.is_finite())
        || currents.len() != layout.ports.len()
        || voltage.len() != layout.ports.len()
        || currents.iter().chain(voltage).any(|x| !x.is_finite())
    {
        return Err(invalid(
            "complete finite external and winding inputs required",
        ));
    }
    let mut inputs = external.to_vec();
    inputs.extend(layout.ports.iter().enumerate().map(|(w, port)| match port {
        NonlinearWindingPort::CurrentInput => -currents[w],
        NonlinearWindingPort::VoltageInput { .. } => voltage[w],
    }));
    Ok(inputs)
}

fn bind<C: SpatialConductorMaterial>(
    field: &SpatialNonlinearInduction<'_, C>,
    windings: &SpatialNonlinearWindings,
    materials: &[SpatialNonlinearMaterial],
    sources: Vec<SpatialWindingConductor>,
    previous: Option<&[SpatialWindingConductor]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<SpatialWindingConductor>, Vec<f64>, Vec<f64>)> {
    let nc = field.problem().mesh().tets.len();
    let nw = windings.profiles_per_m2.len();
    if materials.len() != nc || sources.len() > nw.checked_mul(nc).unwrap_or(0) {
        return Err(invalid("bounded complete winding material maps required"));
    }
    // Canonical ordering makes source matching independent of caller insertion.
    let mut sources = sources;
    sources.sort_by_key(|s| (s.winding, s.cell));
    if previous.is_some_and(|old| old.len() != sources.len()) {
        return Err(invalid("winding material refresh changed physical support"));
    }
    let mut support = vec![false; nw * nc];
    let mut filling = vec![0.0; nc];
    let mut cell_resistance = Vec::with_capacity(sources.len());
    let mut resistance = vec![0.0; nw];
    for (index, source) in sources.iter().enumerate() {
        poll(cancelled)?;
        if source.winding >= nw
            || source.cell >= nc
            || !positive(source.fill_fraction)
            || source.fill_fraction > 1.0
            || support[source.winding * nc + source.cell]
            || field.conductors()[source.cell].is_some()
        {
            return Err(invalid(
                "distinct stationary stranded cells disjoint from eddy support required",
            ));
        }
        let profile = windings.profiles_per_m2[source.winding][source.cell];
        let squared = finite(super::dot(profile, profile))?;
        if !positive(squared) || !source.conductor.has_admitted_static_sources() {
            return Err(invalid(
                "each stranded source requires an active profile and static resistivity",
            ));
        }
        materials[source.cell]
            .admit_fixed_material_refresh(source.conductor.material())
            .map_err(SpatialMagneticError::Material)?;
        if let Some(old) = previous {
            let old = &old[index];
            if old.winding != source.winding
                || old.cell != source.cell
                || old.fill_fraction.to_bits() != source.fill_fraction.to_bits()
                || !old.conductor.has_same_electrical_source(&source.conductor)
            {
                return Err(invalid(
                    "winding refresh changed geometry, filling or original electrical source",
                ));
            }
            admit_source_temperature_refresh(old.conductor.material(), source.conductor.material())
                .map_err(SpatialMagneticError::Material)?;
            if old
                .conductor
                .with_material(source.conductor.material().clone())
                .map_err(|_| invalid("winding refresh left its original source domain"))?
                != source.conductor
            {
                return Err(invalid(
                    "winding refresh changed its resolved original electrical law",
                ));
            }
        }
        support[source.winding * nc + source.cell] = true;
        filling[source.cell] = finite(filling[source.cell] + source.fill_fraction)?;
        if filling[source.cell] > 1.0 {
            return Err(invalid(
                "insulated winding fractions exceed the represented cell volume",
            ));
        }
        let volume = field.problem().geometry.geometry.vol_signed[source.cell].abs();
        let r = finite(
            volume * (source.conductor.resistivity_ohm_m() / source.fill_fraction) * squared,
        )?;
        if !positive(r) {
            return Err(invalid(
                "positive representable stranded resistance required",
            ));
        }
        cell_resistance.push(r);
        resistance[source.winding] = finite(resistance[source.winding] + r)?;
    }
    for (winding, profile) in windings.profiles_per_m2.iter().enumerate() {
        for (cell, value) in profile.iter().enumerate() {
            poll(cancelled)?;
            if value.iter().any(|&v| v != 0.0) != support[winding * nc + cell] {
                return Err(invalid(
                    "every active winding cell requires its exact conductor source",
                ));
            }
        }
    }
    if resistance.iter().any(|&r| !positive(r)) {
        return Err(invalid(
            "each winding requires represented positive stranded resistance",
        ));
    }
    Ok((sources, cell_resistance, resistance))
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
        poll(cancelled)?;
        let nc = field.problem().mesh().tets.len();
        let ne = field.potential_wb().len();
        let count = windings.profiles_per_m2.len();
        let values = nc
            .checked_mul(3)
            .and_then(|n| n.checked_add(ne))
            .and_then(|n| n.checked_mul(count));
        if count == 0 || count > config.max_windings || config.max_windings > 16
            || values.is_none_or(|n| n > config.max_profile_values)
            || windings.branch_indices.len() != count || windings.end_inductance_h.len() != count
            || terminal_ports.len() != count || graph.node_count > 32 || graph.branches.len() > 64
            || graph.transformers.len() > 16 || graph.branches.iter().any(|&(a,b,_)| a == b)
            || !config.current_balance_tolerance_turns.is_finite() || config.current_balance_tolerance_turns < 0.0
            || windings.end_inductance_h.iter().any(|&l| !l.is_finite() || l < 0.0)
            || windings.profiles_per_m2.iter().any(|p| p.len() != nc || p.iter().flatten().any(|v| !v.is_finite()))
            || terminal_ports.iter().any(|port| matches!(port,
                NonlinearWindingPort::VoltageInput { initial_voltage_v,residual_scale_ohm,current_tolerance_a }
                if !initial_voltage_v.is_finite() || !positive(*residual_scale_ohm) || !positive(*current_tolerance_a))) {
            return Err(invalid("invalid spatial winding geometry, limits, native graph or terminal causality"));
        }
        let mut loads_per_a = Vec::with_capacity(count);
        for profile in &windings.profiles_per_m2 {
            // Closed physical currents cannot disappear at a prescribed A trace.
            for face in &field.problem().geometry.faces {
                poll(cancelled)?;
                let jump = if face.cells.len() == 2 {
                    std::array::from_fn(|j| profile[face.cells[0]][j] - profile[face.cells[1]][j])
                } else {
                    profile[face.cells[0]]
                };
                if finite(face.area_m2 * super::dot(jump, face.outward_normal))?.abs()
                    > config.current_balance_tolerance_turns
                {
                    return Err(invalid(
                        "each winding requires a complete closed current and exterior return path",
                    ));
                }
            }
            loads_per_a.push(field.prepared.coupling.current_load(
                super::geometry(field.problem()),
                profile,
                config.current_balance_tolerance_turns,
                cancelled,
            )?);
        }
        let expected = impressed(&windings.profiles_per_m2, &initial.currents_a, cancelled)?;
        if expected
            .iter()
            .flatten()
            .zip(field.impressed_current_density_a_m2().iter().flatten())
            .any(|(a, b)| (a - b).abs() > 128.0 * f64::EPSILON * a.abs().max(b.abs()))
        {
            return Err(invalid(
                "winding currents disagree with the accepted impressed field history",
            ));
        }
        let mut external_ports = 0;
        for &(_, _, branch) in &graph.branches {
            if let Branch::CurrentSource { port } | Branch::VoltageSource { port } = branch {
                if port >= 32 {
                    return Err(invalid("native source port exceeds capacity"));
                }
                external_ports = external_ports.max(port + 1);
            }
        }
        if external_ports + count > 32 {
            return Err(invalid("combined native port budget exceeded"));
        }
        let mut lowered = graph.clone();
        let mut bound = vec![false; graph.branches.len()];
        for (winding, &branch) in windings.branch_indices.iter().enumerate() {
            poll(cancelled)?;
            if branch >= graph.branches.len()
                || bound[branch]
                || !matches!(graph.branches[branch].2, Branch::Inductor { .. })
            {
                return Err(invalid(
                    "each winding requires a distinct original inductor placeholder",
                ));
            }
            bound[branch] = true;
            lowered.branches[branch].2 = match terminal_ports[winding] {
                NonlinearWindingPort::CurrentInput => Branch::CurrentSource {
                    port: external_ports + winding,
                },
                NonlinearWindingPort::VoltageInput { .. } => Branch::VoltageSource {
                    port: external_ports + winding,
                },
            };
        }
        let dae = assemble_circuit(&lowered)?;
        if dae.system.state_dim() > 128
            || initial.remaining_fluxes_wb.len() != dae.flux_index.len()
            || initial.charges_c.len() != dae.charge_index.len()
            || initial
                .remaining_fluxes_wb
                .iter()
                .chain(&initial.charges_c)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "bounded native state and complete finite initial fluxes/charges required",
            ));
        }
        let (conductors, cell_resistance_ohm, resistance_ohm) = bind(
            &field,
            &windings,
            field.problem().materials(),
            sources,
            None,
            cancelled,
        )?;
        let flux_linkages_wb = loads_per_a
            .iter()
            .enumerate()
            .map(|(w, load)| {
                finite(
                    pairing(load, field.potential_wb(), cancelled)?
                        + windings.end_inductance_h[w] * initial.currents_a[w],
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let layout = Rc::new(Layout {
            graph,
            dae,
            windings: Rc::new(windings),
            ports: terminal_ports,
            external_ports,
            loads_per_a: Rc::from(loads_per_a),
        });
        let initial_voltage = layout
            .ports
            .iter()
            .map(|port| match port {
                NonlinearWindingPort::CurrentInput => 0.0,
                NonlinearWindingPort::VoltageInput {
                    initial_voltage_v, ..
                } => *initial_voltage_v,
            })
            .collect::<Vec<_>>();
        let input = ports(
            &layout,
            &initial.external_inputs,
            &initial.currents_a,
            &initial_voltage,
        )?;
        let mut seed = vec![0.0; layout.dae.system.state_dim()];
        for (&index, &value) in layout
            .dae
            .flux_index
            .iter()
            .zip(&initial.remaining_fluxes_wb)
        {
            seed[index] = value;
        }
        for (&index, &value) in layout.dae.charge_index.iter().zip(&initial.charges_c) {
            seed[index] = value;
        }
        let state = native_seed(&layout.dae, &seed, &input)?;
        let (_, _, g) = layout.dae.system.structure();
        for (winding, port) in layout.ports.iter().enumerate() {
            if let NonlinearWindingPort::VoltageInput {
                current_tolerance_a,
                ..
            } = port
            {
                let p = layout.external_ports + winding;
                let output = (layout.dae.system.differential_dim()..state.len())
                    .try_fold(0.0, |sum, row| {
                        finite(g[row * layout.dae.ports + p].mul_add(0.5 * state[row], sum))
                    })?;
                if finite(-output - initial.currents_a[winding])?.abs() > *current_tolerance_a {
                    return Err(invalid(
                        "native initial current disagrees with accepted winding history",
                    ));
                }
            }
        }
        finite(
            field.field(cancelled)?.stored_energy_j()
                + layout.dae.system.hamiltonian(&state)
                + end_energy(&layout.windings.end_inductance_h, &initial.currents_a)?,
        )?;
        poll(cancelled)?;
        Ok(Self {
            layout,
            graph_materials: Vec::new(),
            field,
            conductors,
            cell_resistance_ohm,
            resistance_ohm,
            state,
            currents_a: initial.currents_a,
            flux_linkages_wb,
        })
    }

    pub fn field(&self) -> &SpatialNonlinearInduction<'m, C> {
        &self.field
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.layout.graph
    }
    pub fn dae(&self) -> &CircuitDae {
        &self.layout.dae
    }
    pub fn windings(&self) -> &SpatialNonlinearWindings {
        &self.layout.windings
    }
    pub fn branch_indices(&self) -> &[usize] {
        &self.layout.windings.branch_indices
    }
    pub fn profiles_per_m2(&self) -> &[Vec<[f64; 3]>] {
        &self.layout.windings.profiles_per_m2
    }
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.layout.windings.end_inductance_h
    }
    pub fn terminal_ports(&self) -> &[NonlinearWindingPort] {
        &self.layout.ports
    }
    pub fn conductors(&self) -> &[SpatialWindingConductor] {
        &self.conductors
    }
    pub fn resistance_ohm(&self) -> &[f64] {
        &self.resistance_ohm
    }
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    pub fn currents_a(&self) -> &[f64] {
        &self.currents_a
    }
    pub fn flux_linkages_wb(&self) -> &[f64] {
        &self.flux_linkages_wb
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
            self.conductors.clone(),
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
            winding_conductors,
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
            winding_conductors,
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
        sources: Vec<SpatialWindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearWindingCandidate<'m, C>> {
        use SpatialNonlinearWindingError::Balance;
        poll(cancelled)?;
        controls.newton.validate()?;
        let dt = target_s - self.time_s();
        if !target_s.is_finite()
            || !positive(dt)
            || controls.maximum_field_trials == 0
            || controls.maximum_newton_iterations == 0
            || !positive(controls.difference_step_a)
            || !positive(controls.terminal_voltage_tolerance_v)
            || !positive(controls.power_tolerance_w)
            || !positive(controls.energy_tolerance_j)
            || velocity_m_s.len() != self.field.problem().mesh().vertex_count
            || velocity_m_s.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "positive finite interval, work budgets and complete velocity required",
            ));
        }
        // Admit external sources before any expensive field trial.
        ports(
            &self.layout,
            external_inputs,
            &self.currents_a,
            &vec![0.0; self.currents_a.len()],
        )?;
        let (conductors, cell_resistance_ohm, resistance_ohm) = if sources == self.conductors
            && problem.materials() == self.field.problem().materials()
        {
            (
                sources,
                self.cell_resistance_ohm.clone(),
                self.resistance_ohm.clone(),
            )
        } else {
            bind(
                &self.field,
                &self.layout.windings,
                problem.materials(),
                sources,
                Some(&self.conductors),
                cancelled,
            )?
        };
        for source in &conductors {
            poll(cancelled)?;
            if self.field.problem().mesh().tets[source.cell]
                .iter()
                .any(|&v| velocity_m_s[v as usize].iter().any(|&x| x != 0.0))
            {
                return Err(invalid("stranded winding geometry must remain stationary"));
            }
        }
        let cancel = RefCell::new(cancelled);
        let trials = Cell::new(0_usize);
        let voltage_input = self
            .terminal_ports()
            .iter()
            .any(|p| matches!(p, NonlinearWindingPort::VoltageInput { .. }));
        let evaluate = |currents: &[f64]| -> Result<Trial<'m, C>> {
            let mut check = || (*cancel.borrow_mut())();
            poll(&mut check)?;
            if trials.get() >= controls.maximum_field_trials {
                return Err(SpatialNonlinearWindingError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            let current = impressed(self.profiles_per_m2(), currents, &mut check)?;
            let field = self.field.step_prepared_candidate(
                problem.clone(),
                std::sync::Arc::clone(&prepared),
                &current,
                velocity_m_s,
                dt,
                controls.field,
                &mut check,
            )?;
            let delta = field
                .state()
                .potential_wb()
                .iter()
                .zip(self.field.potential_wb())
                .map(|(&a, &b)| finite(a - b))
                .collect::<Result<Vec<_>>>()?;
            let field_increments = self
                .layout
                .loads_per_a
                .iter()
                .map(|load| pairing(load, &delta, &mut check))
                .collect::<Result<Vec<_>>>()?;
            let mut linkage = Vec::with_capacity(currents.len());
            let mut inductive_voltage = Vec::with_capacity(currents.len());
            let mut required_voltage = Vec::with_capacity(currents.len());
            for (w, &i) in currents.iter().enumerate() {
                let increment = finite(
                    field_increments[w] + self.end_inductance_h()[w] * (i - self.currents_a[w]),
                )?;
                linkage.push(finite(self.flux_linkages_wb[w] + increment)?);
                let induced = finite(increment / dt)?;
                inductive_voltage.push(induced);
                required_voltage.push(finite(resistance_ohm[w].mul_add(i, induced))?);
            }
            let input = ports(&self.layout, external_inputs, currents, &required_voltage)?;
            let x0 = if voltage_input {
                self.state.clone()
            } else {
                native_seed(self.dae(), &self.state, &input)?
            };
            poll(&mut check)?;
            let (native, _) = self.dae().step_audited(&x0, &input, dt)?;
            poll(&mut check)?;
            let mut voltage = Vec::with_capacity(currents.len());
            let mut voltage_residual = Vec::with_capacity(currents.len());
            let mut current_residual = Vec::with_capacity(currents.len());
            let mut residual = Vec::with_capacity(currents.len());
            for (w, port) in self.terminal_ports().iter().enumerate() {
                let output = finite(native.y[self.layout.external_ports + w])?;
                let v = match port {
                    NonlinearWindingPort::CurrentInput => output,
                    NonlinearWindingPort::VoltageInput { .. } => {
                        let (a, b, _) = self.graph().branches[self.branch_indices()[w]];
                        let potential = |node| {
                            if node == 0 {
                                0.0
                            } else {
                                let index = self.dae().node_potential_index[node - 1];
                                f64::midpoint(x0[index], native.x[index])
                            }
                        };
                        finite(potential(a) - potential(b))?
                    }
                };
                let dv = finite(v - required_voltage[w])?;
                let di = match port {
                    NonlinearWindingPort::CurrentInput => 0.0,
                    NonlinearWindingPort::VoltageInput { .. } => finite(-output - currents[w])?,
                };
                residual.push(match port {
                    NonlinearWindingPort::CurrentInput => dv,
                    NonlinearWindingPort::VoltageInput {
                        residual_scale_ohm, ..
                    } => finite(residual_scale_ohm * di)?,
                });
                voltage.push(v);
                voltage_residual.push(dv);
                current_residual.push(di);
            }
            Ok(Trial {
                field,
                native,
                x0,
                linkage,
                field_increments,
                inductive_voltage,
                voltage,
                voltage_residual,
                current_residual,
                residual,
            })
        };
        let interface = FallibleFiniteDifference::new(
            self.currents_a.len(),
            controls.difference_step_a,
            |x: &[f64]| evaluate(x).map(|trial| trial.residual),
        )
        .map_err(numerical)?;
        let construction =
            NewtonKrylovState::new(&interface, self.currents_a.clone(), controls.newton);
        if let Some(error) = interface.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&interface, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || interface.has_failed()
            });
        if let Some(error) = interface.take_failure() {
            return Err(numerical(error));
        }
        poll(&mut || (*cancel.borrow_mut())())?;
        if result.progress == SolverRunProgress::Paused {
            return Err(SpatialNonlinearWindingError::Cancelled);
        }
        if result.progress != SolverRunProgress::Complete {
            return Err(SpatialNonlinearWindingError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(SpatialNonlinearWindingError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        if trial.voltage_residual.iter().any(|v|v.abs()>controls.terminal_voltage_tolerance_v)
            || self.terminal_ports().iter().zip(&trial.current_residual).any(|(port,di)| {
                matches!(port,NonlinearWindingPort::VoltageInput {current_tolerance_a,..} if di.abs()>*current_tolerance_a)
            }) {
            return Err(Balance("actual native winding current and voltage do not close"));
        }
        let mut check = || (*cancel.borrow_mut())();
        let mut cells = Vec::with_capacity(conductors.len());
        let mut winding_power = vec![0.0; self.currents_a.len()];
        let mut winding_heat = 0.0;
        for (source, &resistance) in conductors.iter().zip(&cell_resistance_ohm) {
            poll(&mut check)?;
            let w = source.winding;
            let volume = self.field.problem().geometry.geometry.vol_signed[source.cell].abs();
            let occupied = finite(volume * source.fill_fraction)?;
            let profile = self.profiles_per_m2()[w][source.cell];
            let density = profile.map(|j| j * solver.x[w] / source.fill_fraction);
            let electric = density.map(|j| source.conductor.resistivity_ohm_m() * j);
            let power = finite(occupied * super::dot(density, electric))?;
            let heat = finite(dt * power)?;
            if density.iter().chain(&electric).any(|x| !x.is_finite())
                || !positive(occupied)
                || power < 0.0
                || heat < 0.0
                || solver.x[w] != 0.0 && (power == 0.0 || heat == 0.0)
                || finite(power - resistance * solver.x[w] * solver.x[w])?.abs()
                    > controls.power_tolerance_w
            {
                return Err(Balance(
                    "intrinsic stranded current, source resistivity and terminal work disagree",
                ));
            }
            winding_power[w] = finite(winding_power[w] + power)?;
            winding_heat = finite(winding_heat + heat)?;
            cells.push(SpatialWindingCellHeat {
                source: source.clone(),
                winding: w,
                cell: source.cell,
                volume_m3: volume,
                conductor_volume_m3: occupied,
                current_density_a_m2: density,
                electric_field_v_m: electric,
                power_w: power,
                heat_j: heat,
            });
        }
        for (w, &power) in winding_power.iter().enumerate() {
            if finite(power - resistance_ohm[w] * solver.x[w] * solver.x[w])?.abs()
                > controls.power_tolerance_w
            {
                return Err(Balance(
                    "total spatial winding heat and terminal series work disagree",
                ));
            }
        }
        let midpoint = trial
            .x0
            .iter()
            .zip(&trial.native.x)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        let (graph_resistors, graph_heat) = ohmic::graph_heat(
            self.graph(),
            self.dae(),
            &midpoint,
            dt,
            controls.power_tolerance_w,
            &mut check,
        )?;
        let external_work = finite(
            dt * external_inputs
                .iter()
                .zip(&trial.native.y)
                .map(|(u, y)| u * y)
                .sum::<f64>(),
        )?;
        let terminal_work = finite(
            dt * solver
                .x
                .iter()
                .zip(&trial.voltage)
                .map(|(i, v)| i * v)
                .sum::<f64>(),
        )?;
        let field_work = pairing(&solver.x, &trial.field_increments, &mut check)?;
        let delta_i = solver
            .x
            .iter()
            .zip(&self.currents_a)
            .map(|(&a, &b)| finite(a - b))
            .collect::<Result<Vec<_>>>()?;
        let end_before = end_energy(self.end_inductance_h(), &self.currents_a)?;
        let end_after = end_energy(self.end_inductance_h(), &solver.x)?;
        let end_numerical = end_energy(self.end_inductance_h(), &delta_i)?;
        let end_work = self
            .end_inductance_h()
            .iter()
            .zip(&delta_i)
            .zip(&solver.x)
            .try_fold(0.0, |sum, ((&l, &di), &i)| finite((l * di).mul_add(i, sum)))?;
        let end_change = self
            .end_inductance_h()
            .iter()
            .zip(&delta_i)
            .zip(self.currents_a.iter().zip(&solver.x))
            .try_fold(0.0, |sum, ((&l, &di), (&old, &new))| {
                finite((l * di).mul_add(f64::midpoint(old, new), sum))
            })?;
        let field_step = trial.field.step();
        let circuit_before = finite(self.dae().system.hamiltonian(&trial.x0))?;
        let circuit_after = finite(self.dae().system.hamiltonian(&trial.native.x))?;
        let stored_before_j = finite(field_step.stored_before_j + circuit_before + end_before)?;
        let stored_after_j = finite(field_step.stored_after_j + circuit_after + end_after)?;
        let stored_change_j =
            finite(field_step.stored_change_j + trial.native.delta_h + end_change)?;
        let heat_j = finite(field_step.heat_j + winding_heat + graph_heat)?;
        let numerical_dissipation_j = finite(field_step.numerical_dissipation_j + end_numerical)?;
        let parameter_work_j = field_step.parameter_work_j;
        let mechanical_work_j = field_step.mechanical_work_j;
        let energy_residual_j = finite(
            stored_change_j + heat_j + numerical_dissipation_j + mechanical_work_j
                - external_work
                - parameter_work_j,
        )?;
        let defects = [
            trial.native.supplied - (external_work - terminal_work),
            trial.native.dissipated - graph_heat,
            trial.native.delta_h - (circuit_after - circuit_before),
            trial.native.balance_residual(),
            field_work - field_step.impressed_work_j,
            end_change + end_numerical - end_work,
            end_after - end_before - end_change,
            terminal_work - field_work - end_work - winding_heat,
            energy_residual_j,
        ];
        if circuit_before < 0.0
            || circuit_after < 0.0
            || heat_j < 0.0
            || numerical_dissipation_j < 0.0
            || trial
                .native
                .x
                .iter()
                .chain(&trial.native.y)
                .any(|x| !x.is_finite())
            || defects
                .iter()
                .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
        {
            return Err(Balance(
                "native circuit/field storage, terminal work, physical heat and numerical loss do not close",
            ));
        }
        let (field, field_step) = trial.field.into_parts();
        let step = SpatialNonlinearWindingStep {
            native: trial.native.clone(),
            field: field_step,
            cells,
            graph_resistors,
            graph_materials: self.graph_materials.clone(),
            currents_a: solver.x.clone(),
            terminal_voltage_v: trial.voltage,
            inductive_voltage_v: trial.inductive_voltage,
            terminal_residual_v: trial.voltage_residual,
            terminal_current_residual_a: trial.current_residual,
            flux_linkages_wb: [self.flux_linkages_wb.clone(), trial.linkage.clone()],
            stored_before_j,
            stored_after_j,
            stored_change_j,
            heat_j,
            supplied_work_j: external_work,
            mechanical_work_j,
            parameter_work_j,
            numerical_dissipation_j,
            energy_residual_j,
            field_trials: trials.get(),
            report: result.report,
        };
        let circuit = Self {
            layout: Rc::clone(&self.layout),
            graph_materials: self.graph_materials.clone(),
            field,
            conductors,
            cell_resistance_ohm,
            resistance_ohm,
            state: trial.native.x,
            currents_a: solver.x,
            flux_linkages_wb: trial.linkage,
        };
        poll(&mut check)?;
        Ok(SpatialNonlinearWindingCandidate { circuit, step })
    }
}

struct Trial<'m, C: SpatialConductorMaterial> {
    field: SpatialNonlinearInductionCandidate<'m, C>,
    native: StepRecord,
    x0: Vec<f64>,
    linkage: Vec<f64>,
    field_increments: Vec<f64>,
    inductive_voltage: Vec<f64>,
    voltage: Vec<f64>,
    voltage_residual: Vec<f64>,
    current_residual: Vec<f64>,
    residual: Vec<f64>,
}
