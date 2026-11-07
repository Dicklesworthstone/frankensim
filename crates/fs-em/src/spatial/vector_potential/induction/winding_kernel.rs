//! Shared native graph terminals on original three-dimensional field histories.
//!
//! Constitutive owners supply only an accepted field and fallible candidate
//! advance. Profile closure, exact Whitney source pairing, terminal causality,
//! end storage, stranded copper and native circuit accounting stay here. No
//! terminal/Jacobian trial can replace or advance the accepted field history.

use super::{
    coupling::Geometry,
    nonlinear::winding::{
        SpatialNonlinearWindingConfig, SpatialNonlinearWindingError, SpatialNonlinearWindingStep,
        SpatialNonlinearWindings, SpatialWindingCellHeat, SpatialWindingConductor,
    },
};
use crate::{
    spatial::{SpatialMagneticError, dot},
    winding::{
        end::{EndWindingInductance, EndWindingMode, EndWindingStep},
        nonlinear::{
            NonlinearWindingControls, NonlinearWindingInitial, NonlinearWindingPort, native_seed,
            resistors::GraphResistorMaterial,
        },
        ohmic,
    },
};
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_material::{
    conductor::ResolvedConductorState, field_storage::admit_source_temperature_refresh,
};
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitGraph, assemble_circuit},
};
use fs_solver::{
    NewtonKrylovState, SolverRunProgress,
    finite_difference::{FallibleFiniteDifference, FiniteDifferenceError},
};
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
};

pub(super) type Result<T> = core::result::Result<T, SpatialNonlinearWindingError>;
pub(super) fn finite(value: f64) -> Result<f64> {
    Ok(super::finite(value)?)
}
pub(super) fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    Ok(super::poll(cancelled)?)
}
pub(super) fn invalid(message: &'static str) -> SpatialNonlinearWindingError {
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

/// Implemented only by the original native field owners inside this module.
pub(super) trait FieldHistory {
    fn geometry(&self) -> Geometry<'_>;
    fn potential_wb(&self) -> &[f64];
    fn time_s(&self) -> f64;
    fn impressed_current_density_a_m2(&self) -> &[[f64; 3]];
    fn conducting_cell(&self, cell: usize) -> bool;
    fn current_load(
        &self,
        profile: &[[f64; 3]],
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>>;
    fn stored_energy_j(&self, cancelled: &mut impl FnMut() -> bool) -> Result<f64>;
}

pub(super) struct Layout {
    pub graph: CircuitGraph,
    pub dae: CircuitDae,
    pub windings: Rc<SpatialNonlinearWindings>,
    pub ports: Vec<NonlinearWindingPort>,
    external_ports: usize,
    loads_per_a: Rc<[Vec<f64>]>,
}
pub(super) struct SourceBinding {
    pub conductors: Vec<SpatialWindingConductor>,
    pub cell_resistance_ohm: Vec<f64>,
    pub resistance_ohm: Vec<f64>,
}

/// Geometry/profile/port data and unchanged source maps are shared across
/// candidates. Only accepted native storage, currents and linkages advance.
#[derive(Clone)]
pub(super) struct Terminal {
    pub layout: Rc<Layout>,
    pub end: Rc<EndWindingInductance>,
    pub graph_materials: Vec<GraphResistorMaterial>,
    pub sources: Rc<SourceBinding>,
    pub state: Vec<f64>,
    pub currents_a: Vec<f64>,
    pub flux_linkages_wb: Vec<f64>,
}

pub(super) struct FieldAccount {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    pub heat_j: f64,
    pub impressed_work_j: f64,
    pub numerical_dissipation_j: f64,
    pub mechanical_work_j: f64,
    pub parameter_work_j: f64,
}
pub(super) struct FieldTrial<S, F> {
    pub state: S,
    pub step: F,
    pub account: FieldAccount,
}
pub(super) struct TerminalCandidate<S, F> {
    pub field: S,
    pub terminal: Terminal,
    pub step: SpatialNonlinearWindingStep<F>,
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

fn bind<S: FieldHistory>(
    field: &S,
    windings: &SpatialNonlinearWindings,
    mut sources: Vec<SpatialWindingConductor>,
    previous: Option<&[SpatialWindingConductor]>,
    mut admit: impl FnMut(usize, &ResolvedConductorState) -> Result<()>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Rc<SourceBinding>> {
    let geometry = field.geometry();
    let nc = geometry.mesh().tets.len();
    let nw = windings.profiles_per_m2.len();
    if sources.len() > nw.checked_mul(nc).unwrap_or(0) {
        return Err(invalid("bounded complete winding material maps required"));
    }
    // Canonical ordering makes source matching independent of caller insertion.
    sources.sort_by_key(|s| (s.winding, s.cell));
    if previous.is_some_and(|old| old.len() != sources.len()) {
        return Err(invalid("winding material refresh changed physical support"));
    }
    let mut support = vec![false; nw * nc];
    let mut filling = vec![0.0; nc];
    let mut cell_resistance_ohm = Vec::with_capacity(sources.len());
    let mut resistance_ohm = vec![0.0; nw];
    for (index, source) in sources.iter().enumerate() {
        poll(cancelled)?;
        if source.winding >= nw
            || source.cell >= nc
            || !positive(source.fill_fraction)
            || source.fill_fraction > 1.0
            || support[source.winding * nc + source.cell]
            || field.conducting_cell(source.cell)
        {
            return Err(invalid(
                "distinct stationary stranded cells disjoint from eddy support required",
            ));
        }
        let profile = windings.profiles_per_m2[source.winding][source.cell];
        let squared = finite(dot(profile, profile))?;
        if !positive(squared) || !source.conductor.has_admitted_static_sources() {
            return Err(invalid(
                "each stranded source requires an active profile and static resistivity",
            ));
        }
        admit(source.cell, &source.conductor)?;
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
        let volume = geometry.geometry.geometry.vol_signed[source.cell].abs();
        let r = finite(
            volume * (source.conductor.resistivity_ohm_m() / source.fill_fraction) * squared,
        )?;
        if !positive(r) {
            return Err(invalid(
                "positive representable stranded resistance required",
            ));
        }
        cell_resistance_ohm.push(r);
        resistance_ohm[source.winding] = finite(resistance_ohm[source.winding] + r)?;
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
    if resistance_ohm.iter().any(|&r| !positive(r)) {
        return Err(invalid(
            "each winding requires represented positive stranded resistance",
        ));
    }
    Ok(Rc::new(SourceBinding {
        conductors: sources,
        cell_resistance_ohm,
        resistance_ohm,
    }))
}

impl Terminal {
    #[allow(clippy::too_many_arguments)]
    pub fn new<S: FieldHistory>(
        field: &S,
        graph: CircuitGraph,
        windings: SpatialNonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        end_modes: Vec<EndWindingMode>,
        sources: Vec<SpatialWindingConductor>,
        initial: NonlinearWindingInitial,
        config: SpatialNonlinearWindingConfig,
        admit: impl FnMut(usize, &ResolvedConductorState) -> Result<()>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let geometry = field.geometry();
        let nc = geometry.mesh().tets.len();
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
        let end = Rc::new(EndWindingInductance::new(
            windings.end_inductance_h.clone(),
            end_modes,
        )?);
        let mut loads_per_a = Vec::with_capacity(count);
        for profile in &windings.profiles_per_m2 {
            // Closed physical currents cannot disappear at a prescribed A trace.
            for face in &geometry.geometry.faces {
                poll(cancelled)?;
                let jump = if face.cells.len() == 2 {
                    std::array::from_fn(|j| profile[face.cells[0]][j] - profile[face.cells[1]][j])
                } else {
                    profile[face.cells[0]]
                };
                if finite(face.area_m2 * dot(jump, face.outward_normal))?.abs()
                    > config.current_balance_tolerance_turns
                {
                    return Err(invalid(
                        "each winding requires a complete closed current and exterior return path",
                    ));
                }
            }
            loads_per_a.push(field.current_load(
                profile,
                config.current_balance_tolerance_turns,
                cancelled,
            )?);
        }
        let expected = impressed(&windings.profiles_per_m2, &initial.currents_a, cancelled)?;
        if expected.len() != field.impressed_current_density_a_m2().len()
            || expected
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
        let sources = bind(field, &windings, sources, None, admit, cancelled)?;
        let end_linkage = end.flux_linkages_wb(&initial.currents_a)?;
        let flux_linkages_wb = loads_per_a
            .iter()
            .enumerate()
            .map(|(w, load)| {
                finite(pairing(load, field.potential_wb(), cancelled)? + end_linkage[w])
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
            field.stored_energy_j(cancelled)?
                + layout.dae.system.hamiltonian(&state)
                + end.stored_energy_j(&initial.currents_a)?,
        )?;
        poll(cancelled)?;
        Ok(Self {
            layout,
            end,
            graph_materials: Vec::new(),
            sources,
            state,
            currents_a: initial.currents_a,
            flux_linkages_wb,
        })
    }

    pub fn bind_sources<S: FieldHistory>(
        &self,
        field: &S,
        sources: Vec<SpatialWindingConductor>,
        unchanged_materials: bool,
        admit: impl FnMut(usize, &ResolvedConductorState) -> Result<()>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Rc<SourceBinding>> {
        poll(cancelled)?;
        if unchanged_materials && sources == self.sources.conductors {
            return Ok(Rc::clone(&self.sources));
        }
        bind(
            field,
            &self.layout.windings,
            sources,
            Some(&self.sources.conductors),
            admit,
            cancelled,
        )
    }

    /// Source changes affect only native dissipation, never accepted storage,
    /// field memory, private-port ordering or the exact Whitney profile loads.
    pub fn with_graph_resistors(
        &self,
        mut sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = self
            .layout
            .graph
            .branches
            .iter()
            .filter(|(_, _, branch)| matches!(branch, Branch::Resistor { .. }))
            .count();
        if sources.len() != count || sources.len() > 64 {
            return Err(invalid(
                "every ordinary graph resistor requires one bounded material source",
            ));
        }
        sources.sort_by_key(|source| source.branch);
        let mut graph = self.layout.graph.clone();
        for (index, source) in sources.iter().enumerate() {
            poll(cancelled)?;
            if index > 0 && sources[index - 1].branch == source.branch
                || !source
                    .conductor
                    .source_property()
                    .is_static_at_fixed_external_conditions()
                || !positive(source.conductor.resistance_ohm())
                || !positive(source.conductor.resistance_ohm().recip())
            {
                return Err(invalid(
                    "unique representable static Ohmic graph sources required",
                ));
            }
            let Some((_, _, branch @ Branch::Resistor { .. })) =
                graph.branches.get_mut(source.branch)
            else {
                return Err(invalid(
                    "material source does not name an original graph resistor",
                ));
            };
            *branch = Branch::Resistor {
                ohms: source.conductor.resistance_ohm(),
            };
        }
        if sources == self.graph_materials {
            poll(cancelled)?;
            return Ok(self.clone());
        }
        let mut lowered = graph.clone();
        for (winding, &branch) in self.layout.windings.branch_indices.iter().enumerate() {
            let port = self.layout.external_ports + winding;
            lowered.branches[branch].2 = match self.layout.ports[winding] {
                NonlinearWindingPort::CurrentInput => Branch::CurrentSource { port },
                NonlinearWindingPort::VoltageInput { .. } => Branch::VoltageSource { port },
            };
        }
        let dae = assemble_circuit(&lowered)?;
        let original = &self.layout.dae;
        if dae.flux_index != original.flux_index
            || dae.charge_index != original.charge_index
            || dae.node_potential_index != original.node_potential_index
            || dae.ports != original.ports
            || dae.system.state_dim() != original.system.state_dim()
            || dae.system.hamiltonian(&self.state).to_bits()
                != original.system.hamiltonian(&self.state).to_bits()
        {
            return Err(invalid(
                "resistance-only reassembly changed original circuit storage or ports",
            ));
        }
        let mut candidate = self.clone();
        candidate.layout = Rc::new(Layout {
            graph,
            dae,
            windings: Rc::clone(&self.layout.windings),
            ports: self.layout.ports.clone(),
            external_ports: self.layout.external_ports,
            loads_per_a: Rc::clone(&self.layout.loads_per_a),
        });
        candidate.graph_materials = sources;
        poll(cancelled)?;
        Ok(candidate)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step<S: FieldHistory, F, FC>(
        &self,
        field: &S,
        target_s: f64,
        sources: Rc<SourceBinding>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls<FC>,
        advance: impl Fn(&[[f64; 3]], f64, &mut dyn FnMut() -> bool) -> Result<FieldTrial<S, F>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TerminalCandidate<S, F>> {
        use SpatialNonlinearWindingError::Balance;
        poll(cancelled)?;
        controls.newton.validate()?;
        let dt = target_s - field.time_s();
        if !target_s.is_finite()
            || !positive(dt)
            || controls.maximum_field_trials == 0
            || controls.maximum_newton_iterations == 0
            || !positive(controls.difference_step_a)
            || !positive(controls.terminal_voltage_tolerance_v)
            || !positive(controls.power_tolerance_w)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(invalid(
                "positive finite interval and terminal work budgets required",
            ));
        }
        ports(
            &self.layout,
            external_inputs,
            &self.currents_a,
            &vec![0.0; self.currents_a.len()],
        )?;
        let geometry = field.geometry();
        let cancel = RefCell::new(cancelled);
        let trials = Cell::new(0_usize);
        let voltage_input = self
            .layout
            .ports
            .iter()
            .any(|p| matches!(p, NonlinearWindingPort::VoltageInput { .. }));
        let evaluate = |currents: &[f64]| -> Result<Trial<S, F>> {
            let mut check = || (*cancel.borrow_mut())();
            poll(&mut check)?;
            if trials.get() >= controls.maximum_field_trials {
                return Err(SpatialNonlinearWindingError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            let current = impressed(&self.layout.windings.profiles_per_m2, currents, &mut check)?;
            let candidate = advance(&current, dt, &mut check)?;
            let delta = candidate
                .state
                .potential_wb()
                .iter()
                .zip(field.potential_wb())
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
            let end = self.end.step(&self.currents_a, currents)?;
            for (w, &i) in currents.iter().enumerate() {
                let increment = finite(field_increments[w] + end.flux_increment_wb[w])?;
                linkage.push(finite(self.flux_linkages_wb[w] + increment)?);
                let induced = finite(increment / dt)?;
                inductive_voltage.push(induced);
                required_voltage.push(finite(sources.resistance_ohm[w].mul_add(i, induced))?);
            }
            let input = ports(&self.layout, external_inputs, currents, &required_voltage)?;
            let x0 = if voltage_input {
                self.state.clone()
            } else {
                native_seed(&self.layout.dae, &self.state, &input)?
            };
            poll(&mut check)?;
            let (native, _) = self.layout.dae.step_audited(&x0, &input, dt)?;
            poll(&mut check)?;
            let mut voltage = Vec::with_capacity(currents.len());
            let mut voltage_residual = Vec::with_capacity(currents.len());
            let mut current_residual = Vec::with_capacity(currents.len());
            let mut residual = Vec::with_capacity(currents.len());
            for (w, port) in self.layout.ports.iter().enumerate() {
                let output = finite(native.y[self.layout.external_ports + w])?;
                let v = match port {
                    NonlinearWindingPort::CurrentInput => output,
                    NonlinearWindingPort::VoltageInput { .. } => {
                        let (a, b, _) =
                            self.layout.graph.branches[self.layout.windings.branch_indices[w]];
                        let potential = |node| {
                            if node == 0 {
                                0.0
                            } else {
                                let index = self.layout.dae.node_potential_index[node - 1];
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
                field: candidate,
                end,
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
        if trial.voltage_residual.iter().any(|v| v.abs() > controls.terminal_voltage_tolerance_v)
            || self.layout.ports.iter().zip(&trial.current_residual).any(|(port, di)| {
                matches!(port, NonlinearWindingPort::VoltageInput { current_tolerance_a, .. } if di.abs() > *current_tolerance_a)
            })
        {
            return Err(Balance("actual native winding current and voltage do not close"));
        }
        let mut check = || (*cancel.borrow_mut())();
        let mut cells = Vec::with_capacity(sources.conductors.len());
        let mut winding_power = vec![0.0; self.currents_a.len()];
        let mut winding_heat = 0.0;
        for (source, &resistance) in sources.conductors.iter().zip(&sources.cell_resistance_ohm) {
            poll(&mut check)?;
            let w = source.winding;
            let volume = geometry.geometry.geometry.vol_signed[source.cell].abs();
            let occupied = finite(volume * source.fill_fraction)?;
            let profile = self.layout.windings.profiles_per_m2[w][source.cell];
            let density = profile.map(|j| j * solver.x[w] / source.fill_fraction);
            let electric = density.map(|j| source.conductor.resistivity_ohm_m() * j);
            let power = finite(occupied * dot(density, electric))?;
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
            if finite(power - sources.resistance_ohm[w] * solver.x[w] * solver.x[w])?.abs()
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
            &self.layout.graph,
            &self.layout.dae,
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
        let end_before = trial.end.stored_before_j;
        let end_after = trial.end.stored_after_j;
        let end_numerical = trial.end.numerical_dissipation_j;
        let end_work = trial.end.endpoint_work_j;
        let end_change = trial.end.stored_change_j;
        let account = &trial.field.account;
        let circuit_before = finite(self.layout.dae.system.hamiltonian(&trial.x0))?;
        let circuit_after = finite(self.layout.dae.system.hamiltonian(&trial.native.x))?;
        let stored_before_j = finite(account.stored_before_j + circuit_before + end_before)?;
        let stored_after_j = finite(account.stored_after_j + circuit_after + end_after)?;
        let stored_change_j = finite(account.stored_change_j + trial.native.delta_h + end_change)?;
        let heat_j = finite(account.heat_j + winding_heat + graph_heat)?;
        let numerical_dissipation_j = finite(account.numerical_dissipation_j + end_numerical)?;
        let parameter_work_j = account.parameter_work_j;
        let mechanical_work_j = account.mechanical_work_j;
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
            field_work - account.impressed_work_j,
            end_change + end_numerical - end_work,
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
        let step = SpatialNonlinearWindingStep {
            native: trial.native.clone(),
            field: trial.field.step,
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
        let terminal = Self {
            layout: Rc::clone(&self.layout),
            end: Rc::clone(&self.end),
            graph_materials: self.graph_materials.clone(),
            sources,
            state: trial.native.x,
            currents_a: solver.x,
            flux_linkages_wb: trial.linkage,
        };
        poll(&mut check)?;
        Ok(TerminalCandidate {
            field: trial.field.state,
            terminal,
            step,
        })
    }
}

struct Trial<S, F> {
    field: FieldTrial<S, F>,
    end: EndWindingStep,
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
